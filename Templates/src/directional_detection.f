!-----------------------------------------------------------------------------------------------------!
        Subroutine directional_detection(dm_response, target)
!-----------------------------------------------------------------------------------------------------!
!       The main subroutine for the calculation of the direct detection differential rate dR/dlogE,   !
!       rate R depending on the earth velocity modulation during the year and total rate R.           !                                                                    !
!       First, the squared ionization amplitude for the DM-electron scattering is computed by         !
!       multipling the dark matter response function with the atomic response functions,              !
!       tabulated in the folder Atomic_responses.                                                     !
!       Than we integrate over the dark matter velocity distribution, considering the motion of the   !
!       Earth respect to the galaxy and the local dark matter distribution.                           !
!                                                                                                     !
!       Please cite:                                                                                  !
!               - arXiv:1912.08204v2 [hep-ph] 19 Aug 2020                                             !
!               - arXiv:1509.01598 [hep-ph] JHEP 05, 046 (2016)                                       !
!                                                                                                     !
!-----------------------------------------------------------------------------------------------------!

        implicit none

        include '../include/maddm.inc'

        double precision dm_response, M_dm, M_e, binding_energy
        integer i, ik, iq, id, tot_num_shell
        character(20) shell_filenames(50), filename
        character(2) target, shell_names(50)
        character(100) maddm_path, atom_resp_path
        double precision diff_rate_shell(gridsize), rate_vs_time_shell(day_bins), tot_rate_shell
        double precision diff_rate(gridsize), rate_vs_time(day_bins), tot_rate
        real(kind=10) atomic_response_matrix(gridsize,gridsize)
        real(kind=10) atomic_response_matrix_tot(gridsize,gridsize)
        real(kind=10) ioniz_amplitude(gridsize,gridsize)
        double precision k_min, k_max, dday
        double precision k_e(gridsize), E_e(gridsize)
        double precision dayvalue(day_bins+1), daymid(day_bins)

        include '../include/coupl.inc'
        include '../include/process_names.inc'
        include '../include/maddm_card.inc'

        write(*,*) '---------------------------------------------------------------'
        write(*,*) '                    directional_detection.f                    '
        write(*,*) '---------------------------------------------------------------'


c       Initialize the parameters and variables
        M_dm = 0.1 * GeV * cs**(-2) ! mdm(1) * GeV * cs**(-2) ! dark matter mass  REMOVE 
        %(electron_mass)s       ! electron mass M_e (GeV)
        write(*,*) 'M_dm           :    ' , M_dm, 'GeV'
        write(*,*) 'M_e            :    ' , M_e, 'GeV'
        %(maddm_path)s          ! MadDM path
        atom_resp_path = trim(maddm_path) // '/Atomic_responses/atomic_response_1/'
        do ik=1,gridsize
            do iq=1,gridsize
                ioniz_amplitude(ik,iq) = 0
                atomic_response_matrix_tot(ik,iq) = 0
            enddo
        enddo
c        dm_response = (4.d0 * 1/sqrt(137.d0) * M_e * 1.0E-5)**2 ! O1 REMOVE
c        dm_response = (3.d0/16.d0)*(16.d0 * 1/sqrt(137.d0) * M_dm * 1.0E-5)**2 ! O4 REMOVE
        dm_response = 16*pi*((M_dm+M_e)**2)*2.56819*1.0E-9

c ------------------------------------------------------------------------------------------------
c Read the atomic response functions and compute the ionization amplitude
c ------------------------------------------------------------------------------------------------

c       Read all the names of the files that contains the tabulated atomic response function
c       of the chosen target. There is one file for each shell.
        call get_shell_filenames(atom_resp_path,target,shell_filenames,tot_num_shell)

c       Get the names of the shell in the order used to read them.
c       The shell names should be in alphabetical order:
c       i=1 -> 1s , i=2 -> 2p , i=3 -> 2s , i=4 -> 3d , i=5 -> 3p , i=6 -> 3s , i=7 -> 4d , 
c       i=8 -> 4p , i=9 -> 4s , i=10 -> 5p , i=11 -> 5s , ...
c       Please check this if you want to select a different combination of shells.
        do i=1,tot_num_shell
            filename = trim(shell_filenames(i))
            shell_names(i) = filename(4:5)
        enddo

c       Read the atomic response function of each file, multiply them for the dark matter
c       response function obtaining the ioniziation amplitude. Integrate it obtaining the rate,
c       and in the end sum over the selected shells.
        do ik=1,gridsize
            diff_rate(ik) = 0
        enddo

        do id=1,day_bins
            rate_vs_time(id) = 0
        enddo

        tot_rate = 0

        do i=1,tot_num_shell ! REMOVE . To sum in every shell: do i=1,tot_num_shell
            call read_atomic_response(atomic_response_matrix,shell_filenames(i),atom_resp_path)
            do ik=1,gridsize
                do iq=1,gridsize
                    ioniz_amplitude(ik,iq) = dm_response * atomic_response_matrix(ik,iq)
                enddo
            enddo

            call Theory_Simulation(M_dm, M_e, ioniz_amplitude,
     &                             binding_energy(target,shell_names(i)),
     &                             diff_rate_shell, rate_vs_time_shell, tot_rate_shell)


            do ik=1,gridsize
                diff_rate(ik) = diff_rate(ik) + diff_rate_shell(ik) 
            enddo

            do id=1,day_bins
                rate_vs_time(id) = rate_vs_time(id) + rate_vs_time_shell(id)
            enddo

            tot_rate = tot_rate + tot_rate_shell

            write(*,*) '|  Reading: ', trim(shell_filenames(i)),
     &                 '  |  shell name: ', shell_names(i), '  |  E_bind: ',
     &                 binding_energy(target,shell_names(i)), 'GeV  |  tot rate (shell):' ,
     &                 tot_rate_shell, '|'

        enddo 


c ------------------------------------------------------------------------------------------------
c       Setup the bins and write the results
c ------------------------------------------------------------------------------------------------

c       Setting up the array of the momentum k_e and the array of the energy E_e of the final electron.
c       The array k_e is written in log scale.
        k_min = 1.0E-7 ! GeV
        k_max = 1.0E-4 ! GeV

        do ik = 1,gridsize
            k_e(ik) = exp(log(k_min) + (log(k_max) - log(k_min))*dble(ik-1)/dble(gridsize-1))
            E_e(ik) = sqrt(M_e**2 + k_e(ik)**2) - M_e
        enddo

c       Setup the array of day values. Each day value is in the centre of the bin
        dday = (day_max - day_min)/dble(day_bins)

        do id = 1, day_bins+1
            dayvalue(id) =  day_min + dble(id-1)*dday
        enddo

        do id = 1, day_bins
             daymid(id) = (dayvalue(id) + dayvalue(id+1))/2.d0
        enddo

c       Write the results
        open(3,file='./output/dRdlogE_sm.dat',status='unknown') 
        open(4,file='./output/rate_vs_time.dat',status='unknown')
        open(5,file='./output/tot_rate.dat',status='unknown')

c       Writing out the differential rate dRdlogE. 
c       ================================================================
c       Differential recoil rate dRdlogE  ./Output/DATA/dRdlogE_sm.dat    
c       E_e(ik), dR/dlogE(ik)
c       ================================================================
        write(3,*) '## Momentum (keV), Energy(keV), dR/dlogE[events/kg/yr]' 
        write(3,*) '## unsmeared  ##'

        do ik = 1, gridsize
            write(3,*) k_e(ik)*1.0E+6, E_e(ik)*1.0E+6, diff_rate(ik)    ! GeV to keV
        enddo

c       Writing out the rate as a function of days to show annual modulation. R = dN/dt(months)
c       ================================================================
c       Recoil Rate R        ./Output/DATA/rate_vs_time.dat
c       day(i), rate_(i)
c       ================================================================
        write(4,*) '## Rates[events/kg/month]'
        write(4,*) '## unsmeared  ##'

        do id = 1, day_bins
            write(4,*) daymid(id), rate_vs_time(id)
        enddo

c       Writing out the total rate. R = dN/dt(years)
c       ================================================================
c       Recoil Rate R        ./Output/DATA/tot_rate.dat
c       rate
c       ================================================================
        write(5,*) '## Rate[events/kg/year]'
        write(5,*) '## unsmeared  ##'

        write(5,'(E11.5)') tot_rate
        

        close(3)
        close(4)
        close(5)


c       Writing out the total rate
        write(*,'(A21,E11.5,A13)') 'Total rate     :    ' , tot_rate , ' events/kg/yr'
        write(*,*) ''

        Return

        End





!-------------------------------------------------------------------------------------------------------!
        subroutine Theory_Simulation(M_dm,M_e,ioniz_amplitude,E_binding,diff_rate,rate_vs_time,tot_rate)
!-------------------------------------------------------------------------------------------------------!
        implicit none
        include '../include/maddm.inc'

        double precision M_dm, M_e, E_binding
        real(kind=10) ioniz_amplitude(gridsize,gridsize)
        double precision dRdlogE, q_min, q_max, k_min, k_max, dday
        double precision Theory(gridsize,day_bins), diff_rate(gridsize), rate_vs_time(day_bins), tot_rate

        integer ik, iq, id
        double precision dayvalue(day_bins+1), daymid(day_bins)
        double precision q_exc(gridsize), k_e(gridsize), E_e(gridsize)
        double precision dlogE(gridsize), dq(gridsize)
        double precision q_upper, k_e_upper, E_e_upper

        include '../include/maddm_card.inc'


c       Initialize the parameters
        tot_rate = 0

        do ik=1,gridsize
            diff_rate(ik)=0
        enddo

        do id=1,day_bins
            rate_vs_time(id)=0
        enddo

c       Range used in the computation of the ioniz_amplitude
        q_min = 1.0E-7 ! GeV
        q_max = 1.0E-3 ! GeV

        k_min = 1.0E-7 ! GeV
        k_max = 1.0E-4 ! GeV

c       Time step
        dday = (day_max - day_min)/dble(day_bins)

c ------------------------------------------------------------------------------------------------
c       Setting up each bin to calculate the double differential distribution
c ------------------------------------------------------------------------------------------------

c       Setting up the array of the momentum k_e and the array of the energy E_e of the final electron.
c       The array k_e is written in log scale.
        do ik = 1,gridsize
            k_e(ik) = exp(log(k_min) + (log(k_max) - log(k_min))*dble(ik-1)/dble(gridsize-1))
            E_e(ik) = sqrt(M_e**2 + k_e(ik)**2) - M_e
        enddo

c       Compute dlogE as dE/E
        do ik=1,gridsize
            if(ik.ne.gridsize) then
                dlogE(ik) = (E_e(ik+1) - E_e(ik))/E_e(ik)
            else
c               Extrapolate the gridsize+1 interval lenght
                k_e_upper = exp(log(k_min) + (log(k_max) - log(k_min))*dble(ik)/dble(gridsize-1))
                E_e_upper = sqrt(M_e**2 + k_e_upper**2) - M_e
                dlogE(ik) = (E_e_upper - E_e(ik))/E_e(ik)
            endif
        enddo

c       Setting up array of the exchanged momentum q_exc. The array q_exc is written in log scale.
        do iq = 1,gridsize
            q_exc(iq) = exp(log(q_min) + (log(q_max) - log(q_min))*dble(iq-1)/dble(gridsize-1))
        enddo

c       Compute dq
        do iq = 1,gridsize
            if(iq.ne.gridsize) then
                dq(iq) = q_exc(iq+1) - q_exc(iq)
            else
c               Extrapolate the gridsize+1 interval lenght
                q_upper = exp(log(q_min) + (log(q_max) - log(q_min))*dble(iq)/dble(gridsize-1))
                dq(iq) = q_upper - q_exc(iq)
            endif
        enddo
        
c       Setting up the day values in the centre of the bin
        do id = 1, day_bins+1
            dayvalue(id) =  day_min + dble(id-1)*dday
        enddo

        do id = 1, day_bins
             daymid(id) = (dayvalue(id) + dayvalue(id+1))/2.d0
        enddo


c ------------------------------------------------------------------------------------------------
c       Calculate the differential rate
c ------------------------------------------------------------------------------------------------

c       Loop over the momenta of the outgoing electron k_e.
        do ik = 1,gridsize

c           Loop over day_bins
            do id = 1, day_bins

c               Compute the differential rate, than use to compute dR/dlogE, total rate and total events
                Theory(ik,id) =
     &             dRdlogE(ik,E_e(ik),E_binding,q_exc,dq,M_dm,M_e,ioniz_amplitude,V_E(daymid(id)))
     &             *dday

                tot_rate         = tot_rate + Theory(ik,id) * dlogE(ik)
                rate_vs_time(id) = rate_vs_time(id) + Theory(ik,id) * dlogE(ik)
                diff_rate(ik)    = diff_rate(ik) + Theory(ik,id)
            enddo
        enddo


        return
        end





!-----------------------------------------------------------------------------------------------!
        Function dRdlogE(ik, ER, Eb, q_exc, dq, M_dm, M_e, ioniz_amplitude, ve)
!-----------------------------------------------------------------------------------------------!
!	This calculates the double differential spectrum for                                        !
!       Dark Matter Directional detection	                                                    !
!	Uses the modulation information of the earth inside function V_E                            !
!	To print out information flag = 1, otherwise flag = 0				                        !
!-----------------------------------------------------------------------------------------------!
        implicit none

        include '../include/maddm.inc'

        integer ik,iq
        double precision dRdlogE, c, ER, Eb, ve, nDM
        double precision M_dm, M_e, v0, vearth, vmin, RhoD
        double precision N_events, eta, r_kin, kNorm, vesc, const_integral
        double precision q_exc(gridsize), dq(gridsize)

        real(kind=10) ioniz_amplitude(gridsize,gridsize)

        include '../include/maddm_card.inc'

c       Parameters and variables
        c       = 3.0E+5 * (km/sec)     ! Speed of light in km/s
c        v0      = vMP*(km/sec)          ! Most Probable velocity of WIMPs in DM Halo
        v0      = 220.d0 * (km/sec)          ! REMOVE
c        vesc    = vescape*(km/sec)      ! Escape velocity of a WIMP from the Galactic Halo
        vesc    = 544.d0 * (km/sec)         ! REMOVE
c        vearth  = ve*(km/sec)           ! Velocity of the Earth, taken from the V_E function (cm/s)
        vearth  = 232.d0 * (km/sec)          ! REMOVE
c        RhoD    = rhoDM*GeV*cs**(-2)*cm**(-3) ! Density of Dark Matter in our local part of the Galaxy
c        nDM = RhoD / M_dm * cm**(-3)         ! Number density of DM
        nDM     = 0.4 / M_dm * cm**(-3)         ! Number density of DM
        const_integral = nDM/(128.d0*pi*(M_dm**2)*(M_e**2)) * c**4 ! Constant in front of the integral

        kNorm   = (v0**3)*pi*(sqrt(pi)*erf(vesc/v0) - 2.d0*(vesc/v0)*exp(-(vesc/v0)**2)) ! Normalization factor for velocity distribution integral
        dRdlogE = 0

c       Loop over the exchanged momenta q_exc
        do iq=1,gridsize

c           Minimum velocity required for recoil.
            vmin = ((ER+Eb)/q_exc(iq) + q_exc(iq)/(2*M_dm))*c

c           Dark matter integrated velocity distibution eta (integrated also in cos(theta))
            if(vmin.le.(vesc-vearth)) then
                eta = (v0**2 * pi)/(2*vearth*kNorm) *
     &                ((-4.d0)*exp(-(vesc/v0)**2)*vearth + sqrt(pi)*v0*(erf((vmin+vearth)/v0) - erf((vmin-vearth)/v0)))
            else if ((vesc-vearth).le.vmin.and.vmin.le.(vesc+vearth)) then
                eta = (v0**2 * pi)/(2*vearth*kNorm) *
     &                ((-2.d0)*exp(-(vesc/v0)**2)*(vesc-vmin+vearth) + sqrt(pi)*v0*(erf(vesc/v0) - erf((vmin-vearth)/v0)))
            else
                eta = 0
            endif

c           Integrate q_exc*eta*ioniz_amplitude over q_exc (second index of ioniz_amplitude).
c           The first index ik is for the momenta of the outgoing electron.
c           Multiply also for the constant in front of the integral
            dRdlogE = dRdlogE + const_integral * dq(iq)/c * q_exc(iq)/c * eta * ioniz_amplitude(ik,iq)

        enddo
      
      return
      end
!-----------------------------------------------------------------------------------------!





!-----------------------------------------------------------------------------------------!
      Real*8 Function V_E(days)
!-----------------------------------------------------------------------------------------!
!     This function calculates the velocity of the earth from the galactic velocity and   !
!     sun proper motion.                                                                  !
!-----------------------------------------------------------------------------------------!
      Implicit none
      
      Real*8	 pi, degtorad, cosbx, cosby, cosbz, g, lambda, Lx, Ly, Lz, L0
      Real*8	 V_E1(3), days, ecc, uearthave, u_E, vuGal(3), vuSun(3), vu_E(3)
      
      vuGal = (/ 0.d0, 230.d0, 0.d0 /) ! velocity of the Galactic rotation
      vuSun = (/ 10.d0, 13.d0, 7.d0 /) ! Sun proper motion relative to nearby stars
      
      pi = 4.d0*Atan(1.d0)
      degtorad = pi/180.d0      ! converting degrees to radians
      cosbx = cos(-5.5303*degtorad) 
      cosby = cos(59.575*degtorad)
      cosbz = cos(29.812*degtorad)
      
!     Following Lewin and Smith Appendix B.
      L0 = 13.d0*degtorad
      Lx = 266.141d0*degtorad
      Ly = -13.3485d0*degtorad
      Lz = 179.3212d0*degtorad
      
      g = degtorad*357.528d0 + 0.9856003d0*days*degtorad 
      lambda = degtorad*280.46d0 + degtorad*0.9856474d0*days +
     &     degtorad*1.915d0*sin(g) + degtorad*0.02d0*sin(2.d0*g) !galactic coordinates of the earth on each day.
      ecc = 0.016722d0          ! Earth eccentricity.
      uearthave = 29.79d0 
      u_E = uearthave*(1.d0 - ecc*sin(lambda - L0)) ! Average velocity of the earth around the Sun.
      vu_E = (/ u_E*cosbx*sin(lambda - Lx), u_E*cosby*sin(lambda - Ly), u_E*cosbz*sin(lambda - Lz) /) 
      
      V_E1(1) = vuGal(1) + vuSun(1) + vu_E(1)
      V_E1(2) = vuGal(2) + vuSun(2) + vu_E(2)
      V_E1(3) = vuGal(3) + vuSun(3) + vu_E(3)
      
      V_E = sqrt(V_E1(1)**2 + V_E1(2)**2 + V_E1(3)**2)
      
      return
      end 





!-----------------------------------------------------------------------------------------!
        subroutine get_shell_filenames(atom_resp_path,target,shell_filenames,num_shell)
!-----------------------------------------------------------------------------------------!
!
!       Read all names of the files that contains the tabulated atomic response
!       function of the chosen target. There is one file for each shell.
!
!-----------------------------------------------------------------------------------------!
            Use, intrinsic :: iso_fortran_env, Only : iostat_end

            character(100) atom_resp_path
            character(2) target
            character(20) filename, shell_filenames(50)
            integer error, num_shell

            include '../include/maddm.inc'
            include '../include/maddm_card.inc'

            num_shell = 0

            call system('ls ' // trim(atom_resp_path) // ' > ' // trim(atom_resp_path) // 'fileNames.txt')
            open(1, FILE = trim(atom_resp_path)//'fileNames.txt', action="read")
                  Do
                        Read(1, *, iostat = error) filename
                        Select Case(error)
                        Case(0)
                              if(filename(1:3).eq.(trim(target)//'_')) then
                                    num_shell = num_shell + 1
                                    shell_filenames(num_shell) = filename
                              endif
                        Case(iostat_end)
                              Exit
                        Case Default
                              Write(*, *) 'Error in reading file'
                              Stop
                        End Select
                  enddo
            close(1)

        end subroutine get_shell_filenames





!---------------------------------------------------------------------------------------------!
        subroutine read_atomic_response(atomic_response_matrix,shell_filename,atom_resp_path)
!---------------------------------------------------------------------------------------------!

                Use, intrinsic :: iso_fortran_env, Only : iostat_end

                include '../include/maddm.inc'                

                real(kind=10) atomic_response_matrix(gridsize,gridsize)
                character(20) shell_filename
                character(100) atom_resp_path
                integer j,k,error

                open(2, FILE = trim(atom_resp_path)//trim(shell_filename), action="read")
c                            write(*,*) 'Opened ', trim(atom_resp_path) // trim(shell_filename)

                Do j=1,gridsize
                Read(2, *, iostat = error) (atomic_response_matrix(j,k), k=1,gridsize)
                Select Case(error)
                Case(0)
c                        do k=1,gridsize
c                                Write(*, *) '(', j , ',' , k , ')' , atomic_response_matrix(j,k)
c                        enddo
c                                Write(*,*) ''
                Case(iostat_end)
                        Exit
                Case Default
                        Write(*, *) 'Error in reading file'
                        Stop
                End Select
                enddo
                close(2)

        end subroutine read_atomic_response





!---------------------------------------------------------------------------------------------!
        function binding_energy(target,shell_name)
!---------------------------------------------------------------------------------------------!
!       Returns the binding energy in GeV of selected target.
!       In case of two distinct binding energies for the different spins inside a shell,
!       we took a simple mean.
!       
!       Please cite: 
!           - https://www.webelements.com/xenon/atoms.html
!---------------------------------------------------------------------------------------------!
            double precision binding_energy
            character(2) target, shell_name

            if(target.eq.'Xe') then

                if(shell_name.eq.'1s') then
                    binding_energy = 33317.56055656222 * 1.0E-9 ! GeV

                else if(shell_name.eq.'2s') then
                    binding_energy = 5149.213639792182 * 1.0E-9 ! GeV

                else if(shell_name.eq.'2p') then
                    binding_energy = 4837.706588171414 * 1.0E-9 ! GeV

                else if(shell_name.eq.'3s') then
                    binding_energy = 1093.2351842564003 * 1.0E-9 ! GeV

                else if(shell_name.eq.'3p') then
                    binding_energy = 958.4299495823894 * 1.0E-9 ! GeV

                else if(shell_name.eq.'3d') then
                    binding_energy = 710.7303605534998 * 1.0E-9 ! GeV

                else if(shell_name.eq.'4s') then
                    binding_energy = 213.7805688618793 * 1.0E-9 ! GeV

                else if(shell_name.eq.'4p') then
                    binding_energy = 163.49493390058458 * 1.0E-9 ! GeV

                else if(shell_name.eq.'4d') then
                    binding_energy = 75.58972072252894 * 1.0E-9 ! GeV

                else if(shell_name.eq.'5s') then
                    binding_energy = 25.69862365041479 * 1.0E-9 ! GeV

                else if(shell_name.eq.'5p') then
                    binding_energy = 12.44330433672413 * 1.0E-9 ! GeV

                else
                    write(*,*) 'Invalid shell_name'

                endif

            endif

            return

        end function binding_energy
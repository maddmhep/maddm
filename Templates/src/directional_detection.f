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

        double precision dm_response, dm_response_c13(gridsize)
        double precision  q_exc(gridsize), E_e(gridsize), k_e(gridsize), v_perp_squared(gridsize) ! REMOVE all
        double precision m_reduced, gamma ! REMOVE all
        real massDM, M_e
        integer i, ik, iq, tot_num_shell
        character(20) shell_filenames(50)
        character(2) target
        character(100) maddm_path, atom_resp_path
        double precision Theory(gridsize,day_bins)
        real(kind=10) atomic_response_matrix(gridsize,gridsize)
        real(kind=10) atomic_response_matrix_tot(gridsize,gridsize)
        real(kind=10) ioniz_amplitude(gridsize,gridsize)
c        type(atomic_response_nl) atom_res_nl(50)

        include '../include/coupl.inc'
        include '../include/process_names.inc'
        include '../include/maddm_card.inc'

        write(*,*) '---------------------------------------------------------------'
        write(*,*) '                    directional_detection.f                    '
        write(*,*) '---------------------------------------------------------------'


c       Initialize the parameters
        massDM = 1 ! mdm(1) REMOVE 
        %(electron_mass)s       ! electron mass M_e (GeV)
        %(maddm_path)s          ! MadDM path
        atom_resp_path = trim(maddm_path) // '/Atomic_responses/'
        do ik=1,gridsize
            do iq=1,gridsize
                ioniz_amplitude(ik,iq) = 0
                atomic_response_matrix_tot(ik,iq) = 0
            enddo
        enddo

c       REMOVE -----------------------------------------------------------------------------------------
        m_reduced = massDM*M_e/(massDM + M_e)
        do iq = 1,gridsize
            q_exc(iq) = exp(log(1.0E-7) + (log(1.0E-3) - log(1.0E-7))*dble(iq-1)/dble(gridsize-1))
            k_e(iq) = exp(log(1.0E-7) + (log(1.0E-3) - log(1.0E-7))*dble(iq-1)/dble(gridsize-1))
            E_e(iq) = sqrt(M_e**2 + k_e(iq)**2) - M_e
            gamma = sqrt(M_e**2 + k_e(iq)**2)/M_e 
            v_perp_squared(iq) = (k_e(iq) * 1.0E+5 / (M_e*gamma) )**2 +
     &                   (q_exc(iq)**2)/(4*m_reduced**2)*(massDM-M_e)*(massDM+M_e) * 1.0E+10 - 
     &                   E_e(iq)/m_reduced * 1.0E+10
            dm_response_c13(iq) = 1.0E-6 * (q_exc(iq) * 1.0E+5 / M_e)**2 * v_perp_squared(iq)
        enddo
c       REMOVE -----------------------------------------------------------------------------------------



c ------------------------------------------------------------------------------------------------
c Read the atomic response functions and compute the ionization amplitude
c ------------------------------------------------------------------------------------------------

c       Read all the names of the files that contains the tabulated atomic response function
c       of the chosen target. There is one file for each shell.
        call get_shell_filenames(atom_resp_path,target,shell_filenames,tot_num_shell)


c       Read the atomic response function of each file, then sum the contributions of each each shell,
c       multiplying for the dark matter response function
        do i=1,(tot_num_shell)
                call read_atomic_response(atomic_response_matrix,shell_filenames(i),atom_resp_path)
                do ik=1,gridsize
                    do iq=1,gridsize
                            ioniz_amplitude(ik,iq) = ioniz_amplitude(ik,iq) +
     &                                               dm_response_c13(iq) * atomic_response_matrix(ik,iq) ! REMOVE dm resp c13
                    enddo
                enddo
        enddo

c ------------------------------------------------------------------------------------------------
c Compute the differential tot_rate and print it
c ------------------------------------------------------------------------------------------------

        call Theory_Simulation(massDM, M_e, ioniz_amplitude, Theory)



        Return

        End





!-----------------------------------------------------------------------------------------------!
        subroutine Theory_Simulation(massDM, M_e, ioniz_amplitude, Theory)
!-----------------------------------------------------------------------------------------------!
        implicit none
        include '../include/maddm.inc'

        real massDM, M_e
        real(kind=10) ioniz_amplitude(gridsize,gridsize)
        double precision dRdlogE, q_min, q_max, k_min, k_max, dday
        double precision Theory(gridsize,day_bins), diff_rate(gridsize), tot_rate(day_bins), tot_events

        Integer ik, iq, id
        double precision dayvalue(day_bins+1)
        double precision q_exc(gridsize), k_e(gridsize), E_e(gridsize)
        double precision dlogE(gridsize)
        double precision k_e_upper, E_e_upper
        double precision daymid(day_bins)

        include '../include/maddm_card.inc'


c       Initialize the parameters
        tot_events = 0

        do ik=1,gridsize
            diff_rate(ik)=0
        enddo

        do id=1,day_bins
            tot_rate(id)=0
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
        enddo

c       Setting up array of the exchanged momentum q_exc. The array q_exc is written in log scale.
        do iq = 1,gridsize
                q_exc(iq) = exp(log(q_min) + (log(q_max) - log(q_min))*dble(iq-1)/dble(gridsize-1))
        enddo
        
c       Setting up the values of "day" in the centre of the bin
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

c           Compute the kinetic energy corresponding to k_e(ik) and the dlogE
            E_e(ik) = sqrt(M_e**2 + k_e(ik)**2) - M_e
            if(ik.ne.gridsize) then
                E_e(ik+1) = sqrt(M_e**2 + k_e(ik+1)**2) - M_e
                dlogE(ik) = log(E_e(ik+1)) - log(E_e(ik))
            else
c               Trying to extrapolate the gridsize+1 interval lenght
                k_e_upper = exp(log(k_min) + (log(k_max) - log(k_min))*dble(ik)/dble(gridsize-1))
                E_e_upper = sqrt(M_e**2 + k_e_upper**2) - M_e
                dlogE(ik) = log(E_e_upper) - log(E_e(ik))
            endif

c           Loop over day_bins
            do id = 1, day_bins
c               Compute the differential rate, than use to compute dR/dlogE, total rate and total events
                Theory(ik,id) =
     &             dRdlogE(ik,E_e(ik),q_exc,massDM,M_e,ioniz_amplitude,V_E(daymid(id)))
     &             *dday

                tot_events   = tot_events   + Theory(ik,id) * dlogE(ik)
                tot_rate(id) = tot_rate(id) + Theory(ik,id) * dlogE(ik)
                diff_rate(ik)  = diff_rate(ik)  + Theory(ik,id)
            enddo
        enddo


c ------------------------------------------------------------------------------------------------
c       Write the results
c ------------------------------------------------------------------------------------------------

        open(3,file='./output/dRdlogE_sm.dat',status='unknown') 
        open(4,file='./output/rate_sm.dat',status='unknown')

c       Writing out the energy distribution dRdlogE. 
c       ================================================================
c       Energy distribution  ./Output/DATA/dRdlogE_sm.dat    
c       E_e(ik), dR/dlogE(ik)
c       ================================================================
        write(3,*) '## Momentum (keV), Energy(keV), dR/dlogE(SI+SD)[events/kg/yr]' 
        write(3,*) '## unsmeared  ##'

        do ik = 1, gridsize
            write(3,*) k_e(ik)*1.0E6, E_e(ik)*1.0E6, diff_rate(ik)    ! GeV to keV
        enddo

c       Writing out the rate as a function of days to show annual modulation. R = dN/dt
c       ================================================================
c       Recoil Rate R        ./Output/DATA/rate_sm.dat
c       day(i), rate_(i)
c       ================================================================
        write(4,*) '## cos\theta, rates(SI+SD)[events/kg/month]'
        write(4,*) '## unsmeared  ##'

        do id = 1, day_bins
            write(4,*) daymid(id), tot_rate(id)
        enddo

        close(3)
        close(4)


c       Writing out the total rate
        write(*,*) 'Total rate: ' , tot_events , 'events/kg/yr'
        write(*,*) ''

        return
        end





!-----------------------------------------------------------------------------------------------!
        Function dRdlogE(ik, ER, q_exc, massDM, M_e, ioniz_amplitude, ve)
!-----------------------------------------------------------------------------------------------!
!	This calculates the double differential spectrum for                                    !
!       Dark Matter Directional detection	                                                !
!	Uses the modulation information of the earth inside function V_E                        !
!	To print out information flag = 1, otherwise flag = 0				        !
!-----------------------------------------------------------------------------------------------!
        implicit none

        include '../include/maddm.inc'

        integer ik,iq
        double precision dRdlogE, c, ER, ve, nDM
        real massDM, M_dm, M_e, v0, vearth, vmin, RhoD
        double precision N_events, eta, r_kin, kNorm, vesc, const_integral
        double precision q_exc(gridsize)

        real(kind=10) ioniz_amplitude(gridsize,gridsize)

        include '../include/maddm_card.inc'


c       Parameters and variables
        c       = 3.d+5                 ! Speed of light in km/s  
        M_dm    = massDM*(GeV/cs**2)    ! dark matter mass (GeV)
        M_e     = M_e*(GeV/cs**2)       ! e mass (GeV)
c        v0      = vMP*(km/sec)          ! Most Probable velocity of WIMPs in DM Halo
        v0      = 220*(km/sec)          ! REMOVE
c        vesc    = vescape*(km/sec)      ! Escape velocity of a WIMP from the Galactic Halo
        vesc    = 544 *(km/sec)         ! REMOVE
c        vearth  = ve*(km/sec)           ! Velocity of the Earth, taken from the V_E function (cm/s)
        vearth  = 244*(km/sec)          ! REMOVE
c        RhoD    = rhoDM*GeV*cs**(-2)*cm**(-3) ! Density of Dark Matter in our local part of the Galaxy
        RhoD = 0.4 ! REMOVE
        nDM     = RhoD / M_dm           ! Number density of DM
        const_integral = RhoD/(128*pi*(M_dm**2)*(M_e**2)) ! Constant in front of the integral

        kNorm   = (v0**3)*pi*(sqrt(pi)*erf(vesc/v0) - 2.d0*(vesc/v0)*exp(-(vesc/v0)**2)) ! Normalisation factor for velocity distribution integral

        dRdlogE = 0

c       Loop over the exchanged momenta q_exc
        do iq=1,gridsize

c           Minimum velocity required for recoil.
            vmin = (ER/q_exc(iq) + q_exc(iq)/(2*M_dm)) * c * (km/sec)

c           Dark matter integrated velocity distibution eta (integrated also in cos(theta))
            if(vmin.le.(vesc-vearth)) then
                eta = (v0**2 * pi)/(2*vearth*kNorm) *
     &                ((-4)*exp(-(vesc/v0)**2)*vearth + sqrt(pi)*v0*(erf((vmin+vearth)/(v0)) - erf((vmin-vearth)/(v0))))
            else if ((vesc-vearth).le.vmin.and.vmin.le.(vesc+vearth)) then
                eta = (v0**2 * pi)/(2*vearth*kNorm) *
     &                ((-2)*exp(-(vesc/v0)**2)*(vesc-vmin+vearth) + sqrt(pi)*v0*(erf(vesc/v0) - erf((vmin-vearth)/(v0))))
            else
                eta = 0
            endif

c           Integrate q_exc*eta*ioniz_amplitude over q_exc (second index of ioniz_amplitude).
c           The first index ik is for the momenta of the outgoing electron.
c           Multiply also for the constant in front of the integral
            dRdlogE = dRdlogE + const_integral * (q_exc(iq) / c) * eta * ioniz_amplitude(ik,iq)

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
c                              do k=1,gridsize
c                                          Write(*, *) '(', j , ',' , k , ')' , this%%atomic_response_matrix(j,k)
c                              enddo
c                                    Write(*,*) ''
                Case(iostat_end)
                        Exit
                Case Default
                        Write(*, *) 'Error in reading file'
                        Stop
                End Select
                enddo
                close(2)

        end subroutine read_atomic_response
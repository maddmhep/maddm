!-----------------------------------------------------------------------------------------------------!
        Subroutine electron_recoil_signal(dm_response,bin_Xenon10,bin_Xenon1T,sig_Xenon10,sig_Xenon1T)
!-----------------------------------------------------------------------------------------------------!
!       The main subroutine for the calculation of the electron recoil signal from DM-electron        !
!       interaction.                                                                                  !
!                                                                                                     !
!       Please cite:                                                                                  !
!               - arXiv:1912.08204v2 [hep-ph] 19 Aug 2020                                             !
!               - arXiv:1509.01598 [hep-ph] JHEP 05, 046 (2016)                                       !
!                                                                                                     !
!-----------------------------------------------------------------------------------------------------!

        implicit none

        include '../include/maddm.inc'

        double precision M_dm, M_e, binding_energy
        integer i, ik, iq, id, iS2, tot_num_shell
        character(20) shell_filenames(50), filename
        character(2) target, shell_names(50)
        character(200) maddm_path, atom_resp_path, Xenon10_efficency_path, Xenon1T_efficency_path
        double precision diff_rate_shell(gridsize_k), rate_vs_time_shell(day_bins), tot_rate_shell
        double precision dRdS2_shell(gridsize_k)
        double precision diff_rate(gridsize_k), rate_vs_time(day_bins), tot_rate
        real(kind=10) atomic_response_matrix(gridsize_k,gridsize_q)
        real(kind=10) atomic_response_matrix_tot(gridsize_k,gridsize_q)
        real(kind=10) ioniz_amplitude(gridsize_k,gridsize_q)
        double precision dday
        double precision k_e(gridsize_k+1), E_e(gridsize_k+1)
        double precision dayvalue(day_bins+1), daymid(day_bins)
        double precision dRdS2_shell_Xenon10(S2_max_Xenon10), dRdS2_shell_Xenon1T(S2_max_Xenon1T)
        double precision dRdS2_Xenon10(S2_max_Xenon10), dRdS2_Xenon1T(S2_max_Xenon1T)
        double precision Xenon10_trigger_efficency(S2_max_Xenon10), Xenon1T_tot_efficency(S2_max_Xenon1T)
        double precision Xenon10_flat_efficency, exposure_Xenon10, exposure_Xenon1T
        parameter(exposure_Xenon10 = 15)        ! Kg days
        parameter(exposure_Xenon1T = 80755.2)     ! Kg days
        double precision tot_sig_Xenon10, tot_sig_Xenon1T

        include '../include/coupl.inc'
        include '../include/process_names.inc'
        include '../include/maddm_card.inc'

c        write(*,*) '---------------------------------------------------------------'
c        write(*,*) '                    electron_recoil_signal.f                   '
c        write(*,*) '---------------------------------------------------------------'


c       Initialize the parameters and variables
        target = 'Xe'
        M_dm = mdm(1) * GeV * cs**(-2) ! dark matter mass
        %(electron_mass)s       ! electron mass M_e (GeV)
        %(maddm_path)s          ! MadDM path
        atom_resp_path = trim(maddm_path) // '/Atomic_responses/atomic_response_1/'
        Xenon10_efficency_path = trim(maddm_path) // '/Efficencies/Xenon10_TriggerEfficiency.txt'
        Xenon1T_efficency_path = trim(maddm_path) // '/Efficencies/Xenon1T_TotalEfficiency.txt'


        do ik=1,gridsize_k
            do iq=1,gridsize_q
                ioniz_amplitude(ik,iq) = 0
                atomic_response_matrix_tot(ik,iq) = 0
            enddo
        enddo

        do iS2 = 1,S2_max_Xenon10
            dRdS2_shell_Xenon10(iS2) = 0
            dRdS2_Xenon10(iS2) = 0
        enddo

        do iS2 = 1,S2_max_Xenon1T
            dRdS2_shell_Xenon1T(iS2) = 0
            dRdS2_Xenon1T(iS2) = 0
        enddo

        do i=1,7
            sig_Xenon10(i) = 0
        enddo

        tot_sig_Xenon10 = 0

        do i=1,4
            sig_Xenon1T(i) = 0
        enddo

        tot_sig_Xenon1T = 0


c       Setting up the array of the momentum k_e and the array of the energy E_e of the final electron.
c       The array k_e is written in log scale.
c       The arrays are of shape gridsize_k+1 because we need (gridsize_k+1)-elements to compute
c       (gridsize_k)-differentials

        open(10000, file='/home/gianmarcolucchetti/thesis/Code/E.dat', status='unknown')
        do ik = 1,gridsize_k+1
            k_e(ik) = exp(log(k_min) + (log(k_max) - log(k_min))*dble(ik-1)/dble(gridsize_k-1))
            E_e(ik) = sqrt(M_e**2 + k_e(ik)**2) - M_e
            write(10000,*) E_e(ik), ','
        enddo
        close(10000)
c ------------------------------------------------------------------------------------------------
c Read the atomic response functions, compute the ionization amplitude and compute the rates
c ------------------------------------------------------------------------------------------------

c       Read all the names of the files that contains the tabulated atomic response function
c       of the chosen target. There is one file for each shell.
        call get_shell_filenames(atom_resp_path,target,shell_filenames,tot_num_shell)

c       Get the names of the shell in the order used to read them.
        do i=1,tot_num_shell
            filename = trim(shell_filenames(i))
            shell_names(i) = filename(4:5)
        enddo

c       Read the atomic response function of each file, multiply them for the dark matter
c       response function obtaining the ioniziation amplitude. Integrate it obtaining the rate,
c       and in the end sum over the selected shells.
        do ik=1,gridsize_k
            diff_rate(ik) = 0
        enddo

        do id=1,day_bins
            rate_vs_time(id) = 0
        enddo

        tot_rate = 0

c       Sum over all shells
        do i=1,tot_num_shell

            tot_rate_shell = 0

            do ik=1,gridsize_k
                diff_rate_shell(ik)=0
            enddo

            do id=1,day_bins
                rate_vs_time_shell(id)=0
            enddo

            call read_atomic_response(atomic_response_matrix,shell_filenames(i),atom_resp_path)

            do ik=1,gridsize_k
                do iq=1,gridsize_q
                    ioniz_amplitude(ik,iq) = atomic_response_matrix(ik,iq) *
     &                                       dm_response * GeVtoKg**2 * gevtopb * picobarn
                enddo
            enddo

            call get_rate(M_dm, M_e, E_e, ioniz_amplitude,
     &                             binding_energy(target,shell_names(i)),
     &                             diff_rate_shell, rate_vs_time_shell, tot_rate_shell)
            
c            open(200+i, file='./output/dRdE_' // trim(shell_names(i)) // '.dat', status='unknown')
c            write(200+i,*) '# shell: ' // shell_names(i)
c            write(200+i,*) '# E (eV)                    dR/dE (Kg^-1 day^-1 KeV^-1)'
c            do ik=1,gridsize_k
c                write(200+i,*) E_e(ik)*1E+9 , diff_rate_shell(ik) / E_e(ik) * 1E-6 / 365.d0 ! E_e (eV), dR/dE (Kg^-1 day^-1 KeV^-1)
c            enddo
c            close(200+i)


c           Avoid to compute the limits for the shells that doesn't contribute to the rate
            if (tot_rate_shell.ne.0) then
                call get_dRdS2_shell(E_e, diff_rate_shell, dRdS2_shell_Xenon10, dRdS2_shell_Xenon1T,
     &                              target, shell_names(i), i)

                do ik=1,gridsize_k
                    diff_rate(ik) = diff_rate(ik) + diff_rate_shell(ik)
                enddo

                do id=1,day_bins
                    rate_vs_time(id) = rate_vs_time(id) + rate_vs_time_shell(id)
                enddo

                do iS2 = 1,S2_max_Xenon10
                    dRdS2_Xenon10(iS2) = dRdS2_Xenon10(iS2) + dRdS2_shell_Xenon10(iS2)
                enddo

                do iS2 = 1,S2_max_Xenon1T
                    dRdS2_Xenon1T(iS2) = dRdS2_Xenon1T(iS2) + dRdS2_shell_Xenon1T(iS2)
                enddo

                tot_rate = tot_rate + tot_rate_shell

            endif

c            write(*,*) '|  Readed: ', trim(shell_filenames(i)),
c     &                 '  |  shell name: ', shell_names(i), '  |  E_bind: ',
c     &                 binding_energy(target,shell_names(i)), 'GeV  |  tot rate (shell):' ,
c     &                 tot_rate_shell, '|'

        enddo

c       Read the efficencies
        call read_efficency(Xenon10_trigger_efficency,S2_max_Xenon10,Xenon10_efficency_path)
        call read_efficency(Xenon1T_tot_efficency,S2_max_Xenon1T,Xenon1T_efficency_path)
        Xenon10_flat_efficency = 0.92

c       Multiply for the exposure and efficency
        do iS2 = 2,S2_max_Xenon10
            dRdS2_Xenon10(iS2) = dRdS2_Xenon10(iS2) * exposure_Xenon10 * Xenon10_trigger_efficency(iS2-1) * Xenon10_flat_efficency
        enddo

        do iS2 = 2,S2_max_Xenon1T
            dRdS2_Xenon1T(iS2) = dRdS2_Xenon1T(iS2) * exposure_Xenon1T * Xenon1T_tot_efficency(iS2-1)
        enddo


c ------------------------------------------------------------------------------------------------
c       Get the number of observed events per bin
c ------------------------------------------------------------------------------------------------

        bin_Xenon10 = (/14,41,68,95,122,149,176,203/)
        bin_Xenon1T = (/150,200,250,300,350/)

c       Xenon10
        do iS2=1,S2_max_Xenon10
            do i=1,7
                if(bin_Xenon10(i).le.iS2.and.iS2.lt.bin_Xenon10(i+1)) then
                    sig_Xenon10(i) = sig_Xenon10(i) + dRdS2_Xenon10(iS2)
                    exit
                endif
            enddo
        enddo

        do i=1,7
            tot_sig_Xenon10 = tot_sig_Xenon10 + sig_Xenon10(i)
        enddo

c        write(*,*) ''
c        write(*,*) 'Xenon10 signal'
c        write(*,*) '------------------------------------'
c        write(*,*) 'bin (S2)   | signal'
c        write(*,*) '------------------------------------'
c        do i=1,7
c            write(*,'(A2,I3,A1,I3,A5,E23.17)') ' [' , bin_Xenon10(i), ',' , bin_Xenon10(i+1), ')  | ' , sig_Xenon10(i)
c        enddo

c       Xenon1T
        do iS2=1,S2_max_Xenon1T
            do i=1,4
                if(bin_Xenon1T(i).le.iS2.and.iS2.lt.bin_Xenon1T(i+1)) then
                    sig_Xenon1T(i) = sig_Xenon1T(i) + dRdS2_Xenon1T(iS2)
                    exit
                endif
            enddo
        enddo

        do i=1,4
            tot_sig_Xenon1T = tot_sig_Xenon1T + sig_Xenon1T(i)
        enddo

c        write(*,*) ''
c        write(*,*) 'Xenon1T signal'
c        write(*,*) '------------------------------------'
c        write(*,*) 'bin (S2)   | signal'
c        write(*,*) '------------------------------------'
c        do i=1,4
c            write(*,'(A2,I3,A1,I3,A5,E23.17)') ' [' , bin_Xenon1T(i), ',' , bin_Xenon1T(i+1), ')  | ' , sig_Xenon1T(i)
c        enddo

c        write(*,*)
c        write(*,*) "Total signal Xenon10:" , tot_sig_Xenon10, "Upper limit for Xenon10: 25.1628"
c        write(*,*) "Total signal Xenon1T:" , tot_sig_Xenon1T, "Upper limit for Xenon1T: 8.6918"

        open(9,file='./output/signal.dat',status='unknown')
        write(9,*) 'Xenon10_bins' , bin_Xenon10
        write(9,*) 'Xenon10_signal' , sig_Xenon10
        write(9,*) 'Xenon1T_bins' , bin_Xenon1T
        write(9,*) 'Xenon1T_signal' , sig_Xenon1T
        close(9)

c ------------------------------------------------------------------------------------------------
c       Setup the bins and write the results
c ------------------------------------------------------------------------------------------------

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
        open(6,file='./output/dRdS2_Xenon10.dat',status='unknown')
        open(7,file='./output/dRdS2_Xenon1T.dat',status='unknown')

c       Writing out the differential rate dRdlogE. 
c       ================================================================
c       Differential recoil rate dRdlogE  ./Output/dRdlogE_sm.dat    
c       E_e(ik), dR/dlogE(ik)
c       ================================================================
        write(3,*) '## Momentum (keV), Energy(keV), dR/dlogE[events/kg/yr]' 
        write(3,*) '## unsmeared  ##'

        do ik = 1, gridsize_k
            write(3,*) k_e(ik)*1.0E+6, E_e(ik)*1.0E+6, diff_rate(ik)    ! GeV to keV
        enddo

c       Writing out the rate as a function of days to show annual modulation. R = dN/dt(months)
c       ================================================================
c       Recoil Rate R        ./Output/rate_vs_time.dat
c       day(i), rate_(i)
c       ================================================================
        write(4,*) '## Rates[events/kg/month]'
        write(4,*) '## unsmeared  ##'

        do id = 1, day_bins
            write(4,*) daymid(id), rate_vs_time(id)
        enddo

c       Writing out the total rate. R = dN/dt(years)
c       ================================================================
c       Recoil Rate R        ./Output/tot_rate.dat
c       rate
c       ================================================================
        write(5,*) '## Rate[events/kg/year]'
        write(5,*) '## unsmeared  ##'

        write(5,'(E11.5)') tot_rate

c       Writing out the differential rate dRdS2 for Xenon10 
c       ================================================================
c       Differential recoil rate dRdS2  ./Output/dRdS2_Xenon10.dat    
c       S2(iS2), dR/dSE(iS2)
c       ================================================================
        write(6,*) '## S2      dR/dSE[events/kg/yr]' 
        write(6,*) '## unsmeared  ##'

        do iS2 = 1, S2_max_Xenon10
            write(6,*) iS2, dRdS2_Xenon10(iS2)
        enddo

c       Writing out the differential rate dRdS2 for Xenon1T 
c       ================================================================
c       Differential recoil rate dRdS2  ./Output/dRdS2_Xenon1T.dat    
c       S2(iS2), dR/dSE(iS2)
c       ================================================================
        write(7,*) '## S2      dR/dSE[events/kg/yr]' 
        write(7,*) '## unsmeared  ##'

        do iS2 = 1, S2_max_Xenon1T
            write(7,*) iS2, dRdS2_Xenon1T(iS2)
        enddo
        

        close(3)
        close(4)
        close(5)
        close(6)
        close(7)

        Return

        End





!-------------------------------------------------------------------------------------------------------!
        subroutine get_rate(M_dm,M_e,E_e,ioniz_amplitude,E_binding,diff_rate,rate_vs_time,tot_rate)
!-------------------------------------------------------------------------------------------------------!
        implicit none
        include '../include/maddm.inc'

        double precision M_dm, M_e, E_binding
        real(kind=10) ioniz_amplitude(gridsize_k,gridsize_q)
        double precision dRdlogE, dday, v_earth
        double precision rate_logE_days(gridsize_k,day_bins), diff_rate(gridsize_k), rate_vs_time(day_bins), tot_rate

        integer ik, iq, id
        double precision dayvalue(day_bins+1), daymid(day_bins)
        double precision q_exc(gridsize_q+1), E_e(gridsize_k+1)
        double precision dlogE(gridsize_k), dq(gridsize_q)

        include '../include/maddm_card.inc'


c       Initialize the parameters
        tot_rate = 0

        do ik=1,gridsize_k
            diff_rate(ik)=0
        enddo

        do id=1,day_bins
            rate_vs_time(id)=0
        enddo

c       Time step
        dday = (day_max - day_min)/dble(day_bins)

c ------------------------------------------------------------------------------------------------
c       Setting up each bin to calculate the double differential distribution
c ------------------------------------------------------------------------------------------------

c       Compute dlogE as dE/E
        do ik= 1,gridsize_k
            dlogE(ik) = (E_e(ik+1) - E_e(ik))/E_e(ik)
        enddo

c       Setting up array of the exchanged momentum q_exc. The array q_exc is written in log scale.
        do iq = 1,gridsize_q+1
            q_exc(iq) = exp(log(q_min) + (log(q_max) - log(q_min))*dble(iq-1)/dble(gridsize_q-1))
        enddo

c       Compute dq
        do iq = 1,gridsize_q
            dq(iq) = q_exc(iq+1) - q_exc(iq)
        enddo
        
c       Setting up the day values in the centre of the bin
        do id = 1,day_bins+1
            dayvalue(id) = day_min + dble(id-1)*dday
        enddo

        do id = 1,day_bins
            daymid(id) = (dayvalue(id) + dayvalue(id+1))/2.d0
        enddo


c ------------------------------------------------------------------------------------------------
c       Calculate the differential rate
c ------------------------------------------------------------------------------------------------

c       Loop over the momenta of the outgoing electron k_e.
        do ik = 1,gridsize_k

c           Loop over day_bins
            do id = 1, day_bins

c               Compute the differential rate, than use to compute dR/dlogE, total rate and total events
                rate_logE_days(ik,id) =
     &             dRdlogE(ik,E_e(ik),E_binding,q_exc,dq,M_dm,M_e,ioniz_amplitude,v_earth(daymid(id)))
     &             *dday/365.0*daytosec

                tot_rate         = tot_rate + rate_logE_days(ik,id) * dlogE(ik)
                rate_vs_time(id) = rate_vs_time(id) + rate_logE_days(ik,id) * dlogE(ik)
                diff_rate(ik)    = diff_rate(ik) + rate_logE_days(ik,id)
            enddo
        enddo


        return
        end





!-----------------------------------------------------------------------------------------------!
        Function dRdlogE(ik, ER, Eb, q_exc, dq, M_dm, M_e, ioniz_amplitude, ve)
!-----------------------------------------------------------------------------------------------!
!	    Calculate the double differential spectrum for dark Matter detection.                   !
!	    Uses the modulation information of the earth inside function v_earth.                   !
!-----------------------------------------------------------------------------------------------!
        implicit none

        include '../include/maddm.inc'

        integer ik,iq
        double precision dRdlogE, ER, Eb, ve, nDM, c
        double precision M_dm, M_e, v0, vearth, vmin, RhoD
        double precision N_events, eta, r_kin, kNorm, vesc, const_integral
        double precision q_exc(gridsize_q+1), dq(gridsize_q)

        real(kind=10) ioniz_amplitude(gridsize_k,gridsize_q)

        include '../include/maddm_card.inc'

c       Parameters and variables
        c       = 29979245800.d0        ! cm sec**(-1)
c        v0      = vMP*(km/sec)          ! Most Probable velocity of WIMPs in DM Halo
        v0      = 220.d0 * (km/sec)          ! REMOVE
c        vesc    = vescape*(km/sec)      ! Escape velocity of a WIMP from the Galactic Halo
        vesc    = 544.d0 * (km/sec)         ! REMOVE
c        vearth  = ve*(km/sec)           ! Velocity of the Earth, taken from the v_earth function (cm/s)
        vearth  = 232.d0 * (km/sec)          ! REMOVE
c        RhoD    = rhoDM*GeV*cs**(-2)*cm**(-3) ! Density of Dark Matter in our local part of the Galaxy
c        nDM = RhoD / M_dm * cm**(-3)         ! Number density of DM
        nDM     = 0.4 / M_dm * cm**(-3)         ! Number density of DM [cm**(-3)]
        const_integral = nDM/(128.d0*pi*((M_dm*GeVtoKg)**2)*((M_e*GeVtoKg)**2)) ! Constant in front of the integral [cm**(-3) Kg**(-4)]

        kNorm   = (v0**3)*pi*(sqrt(pi)*erf(vesc/v0) - 2.d0*(vesc/v0)*exp(-(vesc/v0)**2)) ! Normalization factor for velocity distribution integral [cm**(3) sec**(-3)]

        dRdlogE = 0

c       Loop over the exchanged momenta q_exc
        do iq=1,gridsize_q

c           Minimum velocity required for recoil.
            vmin = ((ER+Eb)/q_exc(iq) + q_exc(iq)/(2*M_dm))*c ! [cm/sec]

c           Dark matter integrated velocity distibution eta [cm**(-1) sec] (integrated also in cos(theta)) 
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
c           Multiply also for the constant in front of the integral, and divide for the Xenon mass.
            dRdlogE = dRdlogE + const_integral * (dq(iq)*GeVtoJ*mtocm**2)/c * (q_exc(iq)*GeVtoJ*mtocm**2)/c
     &                                         * eta * ioniz_amplitude(ik,iq) / (m_Xenon*au)

        enddo
      
      return

      end
!-----------------------------------------------------------------------------------------!





!-----------------------------------------------------------------------------------------------!
        subroutine get_dRdS2_shell(E_e,diff_rate_shell,dRdS2_shell_Xenon10,dRdS2_shell_Xenon1T,target,shell_name,i)
!-----------------------------------------------------------------------------------------------!
!       Get the differential rate over the number of photomultiplied electrons in Xenon10 and
!       Xenon1T.
!-----------------------------------------------------------------------------------------------!
        implicit none

        include '../include/maddm.inc'

        double precision E_e(gridsize_k+1), diff_rate_shell(gridsize_k)
        double precision dRdS2_shell_Xenon10(S2_max_Xenon10), dRdS2_shell_Xenon1T(S2_max_Xenon1T)
        integer iE, i, g2_Xenon10, g2_Xenon1T, n_electr_max
        parameter(n_electr_max = 90)
        double precision dRdne_shell(n_electr_max)
        character(2) target, shell_name
        integer n_1(gridsize_k), n_2(gridsize_k), n_e(gridsize_k), n_electr
        double precision W, f_e, binomial_prob, gauss_prob
        parameter(f_e = 1.d0/1.2)       ! Fraction of quanta observed as electrons
        double precision sigma_Xenon10, sigma_Xenon1T
        integer S2
        

        do n_electr = 0,n_electr_max
            dRdne_shell(n_electr) = 0
        enddo

        do S2 = 1,S2_max_Xenon10
            dRdS2_shell_Xenon10(S2) = 0
        enddo

        do S2 = 1,S2_max_Xenon1T
            dRdS2_shell_Xenon1T(S2) = 0
        enddo

c       Mean energy used to produce a quanta
        W = 13.8E-9 ! GeV

c       Secondary scintillation gain factor and sigma for Xenon10 and Xenon1T
        g2_Xenon10 = 27
        g2_Xenon1T = 33

        sigma_Xenon10 = 6.7
        sigma_Xenon1T = 7.0

c       Write the results in a file
c        open(10+i, file='./output/dRdne' // trim(shell_name) // '.dat', status='unknown')
c        write(10+i,*) '# shell: ' // shell_name
c        write(10+i,*) '#  n_electr              dR/dn_e'

c       From the number of primary (n_1) and secondary (n_2) quanta (photons and electrons)
c       generate the number of electrons (n_e) from a binomial distribution with n_1 + n_2
c       trials and success probability f_e, plus the initial scattered electron     
        do iE=1,gridsize_k

            n_1(iE) = floor(E_e(iE)/W)
    
c           From Table II of arXiv:1703.00910v1
            if (shell_name.eq.'4d') then
                n_2(iE) = 4
            else if (shell_name.eq.'4p') then
                n_2(iE) = 6
            else if (shell_name.eq.'4s') then
                n_2(iE) = 3
            else
                n_2(iE) = 0
            endif

c            write(*,*) shell_name , iE , 'E_e:' , E_e(iE)
c     &                 ,'  |  n_1:' , n_1(iE) , 'n_2:' , n_2(iE)
c     &                 ,'dRdlogE:' , diff_rate_shell(iE)

        enddo
        
c       Compute the differential rate (dRdne) wrt the number of initial scattered electrons (n_electr).
        do n_electr = 1, n_electr_max
            do iE=1,gridsize_k
cc              dR/dn_e = integral(dE * P(n_e|E_e) * dR/dE) = integral(dE * P(n_e|E_e) * dR/dlogE * 1/E)
cc              Where P(n_e|E_e) is given by a binomial distribution with number of events n_e, number of
cc              trials is n_1(E_e)+n_2(E_e) and probability of success is f_e.
cc              N.B. : the binomial prob. is done for n_electr - 1 because we need to start from 0
                dRdne_shell(n_electr) = dRdne_shell(n_electr) +
     &                                  (E_e(iE+1)-E_e(iE)) * diff_rate_shell(iE) / E_e(iE) *
     &                                  binomial_prob(n_electr-1, n_1(iE)+n_2(iE), f_e)
            enddo

c           Output the data
c            write(10+i,*) n_electr, dRdne_shell(n_electr)

        enddo

c        close(10+i)

c       Write the results in a file
c        open(100+i, file='./output/dRdS2_' // trim(shell_name) // '.dat', status='unknown')
c        write(100+i,*) '# shell: ' // shell_name
c        write(100+i,*) '#        S2     dR/dS2 (Xenon10)          dR/dS2 (Xenon1T)'

c       Compute the differential rate (dRdS2) wrt the number of photomultiplied electrons (S2).
        do S2 = 1, S2_max_Xenon10
            do n_electr = 1, n_electr_max
c               dR/dS2 = integral(dn_e * P(S2|n_e*g2,sigma) * dR/dn_e)
c               Where P(S2|n_e*mean,sigma) is the probability to get S2 from a gaussian distribution
c               with mean n_e*g2
                dRdS2_shell_Xenon10(S2) = dRdS2_shell_Xenon10(S2) +
     &                                  dRdne_shell(n_electr) *
     &                                  gauss_prob(S2, n_electr*g2_Xenon10, sqrt(dble(n_electr))*sigma_Xenon10)
            enddo
        enddo

        do S2 = 1, S2_max_Xenon1T
            do n_electr = 1, n_electr_max
c               dR/dS2 = integral(dn_e * P(S2|n_e*g2,sigma) * dR/dn_e)
c               Where P(S2|n_e*mean,sigma) is the probability to get S2 from a gaussian distribution
c               with mean n_e*g2
                dRdS2_shell_Xenon1T(S2) = dRdS2_shell_Xenon1T(S2) +
     &                                  dRdne_shell(n_electr) *
     &                                  gauss_prob(S2, n_electr*g2_Xenon1T, sqrt(dble(n_electr))*sigma_Xenon1T)
            enddo

c           Output the data
c            if(S2.le.S2_max_Xenon10) then
c                write(100+i,*) S2, dRdS2_shell_Xenon10(S2), dRdS2_shell_Xenon1T(S2)
c            else
c                write(100+i,*) S2, 0, dRdS2_shell_Xenon1T(S2)
c            endif

        enddo

c        close(100+i)

        end subroutine get_dRdS2_shell





!-----------------------------------------------------------------------------------------!
      Real*8 Function v_earth(days)
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
      
      v_earth = sqrt(V_E1(1)**2 + V_E1(2)**2 + V_E1(3)**2)
      
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

                real(kind=10) atomic_response_matrix(gridsize_k,gridsize_q)
                character(20) shell_filename
                character(100) atom_resp_path
                integer j,k,error

                open(2, FILE = trim(atom_resp_path)//trim(shell_filename), action="read")
c                write(*,*) 'Opened ', trim(atom_resp_path) // trim(shell_filename)
                
c               Skip the firts 3 rows 
                Read(2,*)
                Read(2,*)
                Read(2,*)

                Do j=1,gridsize_k
                    Read(2, *, iostat = error) (atomic_response_matrix(j,k), k=1,gridsize_q)
                    Select Case(error)
                    Case(0)
c                            do k=1,gridsize_q
c                                    Write(*, *) '(', j , ',' , k , ')' , atomic_response_matrix(j,k)
c                            enddo
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





!---------------------------------------------------------------------------------------------!
        subroutine read_efficency(efficency,S2_max,path)
!---------------------------------------------------------------------------------------------!

                Use, intrinsic :: iso_fortran_env, Only : iostat_end

                integer S2_max
                double precision efficency(S2_max)
                character(200) path
                integer iS2,error

                open(99, FILE = trim(path), action="read")
            
                do iS2 =1,S2_max
                    Read(99, *, iostat = error) efficency(iS2)
            
                    Select Case(error)
                    Case(0)
c                        write(*,*) iS2, efficency(iS2)
                    Case(iostat_end)
                        Exit
                    Case Default
                        write(*, *) 'Error in reading file'
                        Stop
                    End Select
                enddo
                close(99)

        end subroutine read_efficency        
        




!---------------------------------------------------------------------------------------------!
        function binding_energy(target,shell_name)
!---------------------------------------------------------------------------------------------!
!       Returns the binding energy in GeV of selected target.
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





!---------------------------------------------------------------------------------------------!
        function binomial_prob(k, n, p)
!---------------------------------------------------------------------------------------------!
!       Compute the binomial probability given k successes, with n trials and p probability
!       of success.
!---------------------------------------------------------------------------------------------!
            integer k, n
            double precision p
            double precision num, den, n_combinations, binomial_prob

            if (k.le.n) then
                num = 1.d0
                den = 1.d0 
                do i = 1,k
                    num = num * (n - i + 1)
                    den = den * i
                enddo
                n_combinations = num/den
                binomial_prob = n_combinations * p**k * (1.d0-p)**(n-k)
            else
                binomial_prob = 0
            endif

            return

        end function binomial_prob






!---------------------------------------------------------------------------------------------!
        function gauss_prob(val, mean, sigma)
!---------------------------------------------------------------------------------------------!
!       Compute the gaussian probability to obtain a value val
!---------------------------------------------------------------------------------------------!
            double precision gauss_prob, sigma, pi
            integer val, mean

            pi = 3.14159265359

            gauss_prob = 1.d0/(sigma*sqrt(2.d0*pi)) * exp(-1.d0/2.d0*((val-mean)/sigma)**2)

            return

        end function gauss_prob
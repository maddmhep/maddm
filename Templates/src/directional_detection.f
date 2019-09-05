!-----------------------------------------------------------------------------------------------------!
      Subroutine directional_detection(mass,sigma_wnSI,sigma_wpSI,sigma_wnSD,sigma_wpSD,flag,N_events)
!-----------------------------------------------------------------------------------------------------!
!     The main subroutine for the calculation of the direct and directional detection                 !
!     rates R, dN/dE, dN/dcostheta, d2N/dEdcostheta                                                   !
!-----------------------------------------------------------------------------------------------------!
      Implicit none 
      
      Integer  i, j, k, ix, jy, kz, niter, countsm, Enercounter, flag
      parameter(niter = 400)
      double precision mass, test 
      double precision sigma_wnSI, sigma_wpSI  
      double precision sigma_wnSD, sigma_wpSD  
      double precision Th_events(Energy_bins,cos_theta_bins,day_bins) 
      double precision Th_eventsSI(Energy_bins,cos_theta_bins,day_bins) 
      double precision Th_eventsSD(Energy_bins,cos_theta_bins,day_bins) 
      double precision Theory(Energy_bins,cos_theta_bins,day_bins) 
      double precision TheorySI(Energy_bins,cos_theta_bins,day_bins)
      double precision TheorySD(Energy_bins,cos_theta_bins,day_bins)
      double precision sm_Theory(niter,cos_theta_bins,day_bins)      
      double precision sm_TheorySI(niter,cos_theta_bins,day_bins)      
      double precision sm_TheorySD(niter,cos_theta_bins,day_bins)  
    
      double precision sm_Th_events(niter,cos_theta_bins,day_bins)     
      double precision sm_Th_eventsSI(niter,cos_theta_bins,day_bins)     
      double precision sm_Th_eventsSD(niter,cos_theta_bins,day_bins)     
  
      double precision sm_ang_theory(Energy_bins,cos_theta_bins,day_bins)   
      double precision sm_ang_theorySI(Energy_bins,cos_theta_bins,day_bins)   
      double precision sm_ang_theorySD(Energy_bins,cos_theta_bins,day_bins)   
      double precision Evalue(Energy_bins+1), cosvalue(cos_theta_bins+1)
      double precision dayvalue(day_bins+1), Emid(Energy_bins)
      double precision cosmid(cos_theta_bins), daymid(day_bins)
      double precision Evalues(niter+1), Emids(niter)
      double precision dEnergy, dcos, dday, emins, emaxs, dEs   
      double precision dNdE(Energy_bins), dNdcos(cos_theta_bins)  
      double precision dNdESI(Energy_bins), dNdcosSI(cos_theta_bins) 
      double precision dNdESD(Energy_bins), dNdcosSD(cos_theta_bins)   
      double precision Rate(day_bins), Enersm(Energy_bins)  
      double precision RateSI(day_bins), EnersmSI(Energy_bins) 
      double precision RateSD(day_bins), EnersmSD(Energy_bins)  
      double precision dNdE_sm(Energy_bins), dNdcos_sm(cos_theta_bins)  
      double precision dNdE_smSI(Energy_bins), dNdcos_smSI(cos_theta_bins)  
      double precision dNdE_smSD(Energy_bins), dNdcos_smSD(cos_theta_bins)  
      double precision Rate_sm(day_bins), Rate_smSI(day_bins), Rate_smSD(day_bins)
      double precision E_r(Energy_bins), Enerevsm(Energy_bins)
      double precision EnerevsmSI(Energy_bins), EnerevsmSD(Energy_bins)  
      double precision event_counter, Ener_dist(Energy_bins)
      double precision Ener_distSI(Energy_bins), Ener_distSD(Energy_bins)
      double precision efficiency, N_events, count, countSI, countSD

      include '../include/maddm.inc'
      include '../include/maddm_card.inc'

! Energy smearing range.
      emins = -200.d0
      emaxs = 200.d0

      dEs = (emaxs - emins)/dble(niter)
      dEnergy = (En_max - En_min)/dble(Energy_bins)
      dcos    = (cos_max - cos_min)/dble(cos_theta_bins)
      dday    = (day_max - day_min)/dble(day_bins)

!  setting up the values at the center of each bin for calculation.
      do i = 1, niter+1
         Evalues(i) = emins + dble(i-1)*dEs
      enddo
      
      do i = 1, niter
         Emids(i) = (Evalues(i) + Evalues(i+1))/2.d0
      enddo
      
      do ix = 1, Energy_bins+1
         Evalue(ix) = En_min + dble(ix-1)*dEnergy
      enddo

      do ix = 1, Energy_bins
         Emid(ix) = (Evalue(ix) + Evalue(ix+1))/2.d0
      enddo
      
      do jy = 1, cos_theta_bins+1
         cosvalue(jy) = cos_min + dble(jy-1)*dcos
      enddo

      do jy = 1, cos_theta_bins
         cosmid(jy) = (cosvalue(jy) + cosvalue(jy+1))/2.d0 
      enddo

      do kz = 1, day_bins+1
        dayvalue(kz) =  day_min + dble(kz-1)*dday
      enddo

      do kz = 1, day_bins
         daymid(kz) = (dayvalue(kz) + dayvalue(kz+1))/2.d0
      enddo

!---------------------------------------------- SI+SD -----------------------------------!
      call Theory_Simulation(mass, sigma_wnSI, sigma_wpSI,sigma_wnSD, sigma_wpSD, Theory) ! For total rates

!----------------------------------------------- SI -------------------------------------!
      call Theory_Simulation(mass, sigma_wnSI, sigma_wpSI, 0.d0, 0.d0, TheorySI) ! For SI rates only

!----------------------------------------------- SD -------------------------------------!
      call Theory_Simulation(mass, 0.d0, 0.d0, sigma_wnSD, sigma_wpSD, TheorySD) ! For SD rates only

!-------------------------------------------------- -------------------------------------!

! switch to either apply detector resolution or not. We will assume 100% WIMP detection efficiency


      if ( .not. smearing ) then

c==========================================================
c         Not applying any detector resolution
c==========================================================

!     Initialize the distribution arrays
         do i = 1, Energy_bins
            do j = 1, cos_theta_bins
               do k = 1, day_bins
                  Th_events(i,j,k)   = 0.d0
                  Th_eventsSI(i,j,k) = 0.d0
                  Th_eventsSD(i,j,k) = 0.d0
                  dNdE(i)            = 0.d0
                  dNdESI(i)          = 0.d0
                  dNdESD(i)          = 0.d0
                  dNdcos(j)          = 0.d0
                  dNdcosSI(j)        = 0.d0
                  dNdcosSD(j)        = 0.d0
                  rate(k)            = 0.d0
                  rateSI(k)          = 0.d0
                  rateSD(k)          = 0.d0
               enddo
            enddo
         enddo
         
         open(2,file='./output/d2NdEdcos.dat',status='unknown') 
         open(3,file='./output/dNdE.dat',status='unknown') 
         open(4,file='./output/dNdcos.dat',status='unknown') 
         open(5,file='./output/rate.dat',status='unknown') 
 
c        ==================================================
c        Generating Unsmeared Distributions               
c        ==================================================
        
!  writing out the d2NdEdcostheta distribution at after applying an Energy threshold cut.
c        ==================================================
c        double distribution  ./Output/DATA/d2NdEdcos.dat  
c        ==================================================
         
         write(2,*) '## Energy(keV), cos\theta, days, ' 
         write(2,*) '## d2RdEdcos(SI + SD), '
         write(2,*) '## d2RdEdcos(SI), d2RdEdcos(SD) '
         write(2,*) '## unsmeared  ##'

!     Calculating the event rates per bin after applying a thresold cut.
         do i = 1, Energy_bins
            do j = 1, cos_theta_bins
               do k = 1, day_bins
                  if ((Emid(i) .ge. En_threshold) .and. (Emid(i) .le. En_max)) then 

                     Th_events(i,j,k) = Theory(i,j,k)
                     Th_eventsSI(i,j,k) = TheorySI(i,j,k)
                     Th_eventsSD(i,j,k) = TheorySD(i,j,k)
                  
                     write(2,*) Emid(i), cosmid(j), daymid(k), Th_events(i,j,k), 
     &                 Th_eventsSI(i,j,k), Th_eventsSD(i,j,k)

                     dNdE(i)            = dNdE(i)   + Th_events(i,j,k)
                     dNdcos(j)          = dNdcos(j) + Th_events(i,j,k)
                     rate(k)            = rate(k)   + Th_events(i,j,k)   

                     dNdESI(i)          = dNdESI(i)   + Th_eventsSI(i,j,k)
                     dNdcosSI(j)        = dNdcosSI(j) + Th_eventsSI(i,j,k)
                     rateSI(k)          = rateSI(k)   + Th_eventsSI(i,j,k)   

                     dNdESD(i)          = dNdESD(i)   + Th_eventsSD(i,j,k)
                     dNdcosSD(j)        = dNdcosSD(j) + Th_eventsSD(i,j,k)
                     rateSD(k)          = rateSD(k)   + Th_eventsSD(i,j,k) 
  
                  endif
               enddo
            enddo
         enddo
         

         Enercounter = 0
         do i = 1, Energy_bins
            if ((Emid(i) .ge. En_threshold) .and. (Emid(i) .le. En_max)) then 
               Enercounter              = enercounter + 1
               E_r(Enercounter)         = Emid(i)

               Ener_dist(Enercounter)   = dNdE(i)
               Ener_distSI(Enercounter) = dNdESI(i)
               Ener_distSD(Enercounter) = dNdESD(i)

            endif
         enddo
         

!  writing out the energy distribution dN/dE.                          
c        =============================================================
c        Energy distribution  ./Output/DATA/dNdE.dat    
c        =============================================================
         write(3,*) '## Energy(keV), dR/dE(SI+SD)[events/keV/kg/yr],' 
         write(3,*) '## dR/dE(SI)[events/keV/kg/yr],'
         write(3,*) '## dR/dE(SD)[events/keV/kg/yr]'
         write(3,*) '## unsmeared  ##'

         N_events = 0.d0       ! This gives the number of events as given by spin independent interactions

         do i = 1, Energy_bins - int(En_threshold)
            N_events = N_events + Ener_distSI(i)*efficiency(E_r(i),flag) 
            write(3,*) E_r(i), Ener_dist(i)*efficiency(E_r(i),flag),
     &       Ener_distSI(i)*efficiency(E_r(i),flag), Ener_distSD(i)*efficiency(E_r(i),flag)
         enddo

         
!  writing out the angular distribution dN/dcostheta                           
c        ============================================================
c        Angular distribution ./Output/DATA/dNdcos.dat  
c        ============================================================
         write(4,*) '## cos\theta, dN/dcos(SI+SD)[events/kg/year],' 
         write(4,*) '## dN/dcos(SI)[events/kg/year],'
         write(4,*) '## dN/dcos(SD)[events/kg/year]'
         write(4,*) '## unsmeared  ##'

         do j = 1, cos_theta_bins
            write(4,*) cosmid(j), dNdcos(j), dNdcosSI(j), dNdcosSD(j)
         enddo        

!  writing out the total rate as a function of days to show annual modulation. R = dN/dt                  
c        ============================================================
c        Recoil Rate R        ./Output/DATA/rate.dat      
c        ============================================================
         write(5,*) '## cos\theta, rates(SI+SD)[events/kg/year],' 
         write(5,*) '## rates(SI)[events/kg/year],'
         write(5,*) '## rates(SD)[events/kg/year]'
         write(5,*) '## unsmeared  ##'

         do k = 1, day_bins
            write(5,*) daymid(k), rate(k), rateSI(k), rateSD(k)
         enddo

         
         close(2)
         close(3)
         close(4)
         close(5)
         
      else                      ! Switch on the smearing and calculate the distributions
         
c        ==================================================
c        Generating smeared Distributions                 
c        ==================================================
! Start initialization         
         do i = 1, niter
            do j = 1, cos_theta_bins
               do k = 1, day_bins
                  sm_th_events(i,j,k)   = 0.d0
                  sm_th_eventsSI(i,j,k) = 0.d0
                  sm_th_eventsSD(i,j,k) = 0.d0
                  dNdE_sm(i)            = 0.d0
                  dNdE_smSI(i)          = 0.d0
                  dNdE_smSD(i)          = 0.d0
                  dNdcos_sm(j)          = 0.d0
                  dNdcos_smSI(j)        = 0.d0
                  dNdcos_smSD(j)        = 0.d0
                  rate_sm(k)            = 0.d0
                  rate_smSI(k)          = 0.d0
                  rate_smSD(k)          = 0.d0
               enddo
            enddo
         enddo
         
         do i = 1, Energy_bins
            do j = 1, cos_theta_bins
               do k = 1, day_bins
                  sm_ang_theory(i,j,k) = 0.d0
                  sm_ang_theorySI(i,j,k) = 0.d0
                  sm_ang_theorySD(i,j,k) = 0.d0
               enddo
            enddo
         enddo          
! End initialization
         

         count   = 0.d0
         countSI = 0.d0
         countSD = 0.d0
!------------------------------SI + SD--------------------------!       
  
         call GaussmearAngle(Theory, sm_ang_theory)

         call GaussmearEnergy(sm_ang_theory, sm_theory)

!------------------------------- SI ----------------------------!

         call GaussmearAngle(TheorySI, sm_ang_theorySI)

         call GaussmearEnergy(sm_ang_theorySI, sm_theorySI)

!------------------------------- SD ----------------------------!

         call GaussmearAngle(TheorySD, sm_ang_theorySD)

         call GaussmearEnergy(sm_ang_theorySD, sm_theorySD)

!---------------------------------------------------------------!

         do i = 1, niter
            do j = 1, cos_theta_bins
               do k = 1, day_bins
                  count   = count   +  sm_theory(i,j,k) 
                  countSI = countSI +  sm_theorySI(i,j,k) 
                  countSD = countSD +  sm_theorySD(i,j,k) 
               enddo
            enddo
         enddo 


         open(2,file='./output/d2NdEdcos_sm.dat',status='unknown') 
         open(3,file='./output/dNdE_sm.dat',status='unknown') 
         open(4,file='./output/dNdcos_sm.dat',status='unknown') 
         open(5,file='./output/rate_sm.dat',status='unknown') 
         
!  Printing out the d2NdEdcostheta distribution at after applying an Energy threshold cut.
c        =====================================================
c        double distribution  ./Output/DATA/d2NdEdcos_sm.dat  
c        =====================================================

         write(2,*) '## Energy(keV), cos\theta, days, '
         write(2,*) '## d2RdEdcos(SI + SD), d2RdEdcos(SI),'
         write(2,*) '## d2RdEdcos(SD)'
         write(2,*) '## smeared  ##'

         do i = 1, niter
            do j = 1, cos_theta_bins
               do k = 1, day_bins
                  if ((Emids(i) .ge. En_threshold) .and. (Emids(i) .le. En_max))then

                     sm_th_events(i,j,k) = sm_theory(i,j,k)
                     sm_th_eventsSI(i,j,k) = sm_theorySI(i,j,k)
                     sm_th_eventsSD(i,j,k) = sm_theorySD(i,j,k)

                     write(2,*) Emids(i), cosmid(j), daymid(k), sm_th_events(i,j,k),
     &                sm_th_eventsSI(i,j,k), sm_th_eventsSD(i,j,k)

                     dNdE_sm(i)      = dNdE_sm(i) + sm_th_events(i,j,k)
                     dNdE_smSI(i)    = dNdE_smSI(i) + sm_th_eventsSI(i,j,k)
                     dNdE_smSD(i)    = dNdE_smSD(i) + sm_th_eventsSD(i,j,k)

                     dNdcos_sm(j)    = dNdcos_sm(j) + sm_th_events(i,j,k)
                     dNdcos_smSI(j)  = dNdcos_smSI(j) + sm_th_eventsSI(i,j,k)
                     dNdcos_smSD(j)  = dNdcos_smSD(j) + sm_th_eventsSD(i,j,k)

                     rate_sm(k)      = rate_sm(k) + sm_Th_events(i,j,k)
                     rate_smSI(k)    = rate_smSI(k) + sm_Th_eventsSI(i,j,k)
                     rate_smSD(k)    = rate_smSD(k) + sm_Th_eventsSD(i,j,k)

                  endif
               enddo
            enddo
         enddo

         countsm = 0
         do i = 1, niter
            if ((Emids(i) .ge. En_threshold) .and. (Emids(i) .le. En_max)) then
               countsm = countsm + 1
               Enersm(countsm) = Emids(i) 

               Enerevsm(countsm) = dNdE_sm(i)
               EnerevsmSI(countsm) = dNdE_smSI(i)
               EnerevsmSD(countsm) = dNdE_smSD(i)

            endif
         enddo
      
!  writing out the energy distribution dN/dE.      
c        ================================================================
c        Energy distribution  ./Output/DATA/dNdE_sm.dat    
c        Er_i, dN/dE_i, dN/dESI_i, dN/dESD_i
c        ================================================================
         write(3,*) '## Energy(keV), dR/dE(SI+SD)[events/keV/kg/yr],' 
         write(3,*) '## dR/dE(SI)[events/keV/kg/yr],'
         write(3,*) '## dR/dE(SD)[events/keV/kg/yr]  ##'
         write(3,*) '## smeared  ##'

         N_events = 0.d0
         do i = 1, Energy_bins - int(En_threshold)

            N_events = N_events + EnerevsmSI(i)*efficiency(Enersm(i),flag) 

            write(3,*) Enersm(i), Enerevsm(i)*efficiency(Enersm(i),flag), 
     &       EnerevsmSI(i)*efficiency(Enersm(i),flag),
     &       EnerevsmSD(i)*efficiency(Enersm(i),flag)

         enddo

 
!  writing out the Angular Distribution dN/dcostheta.         
c        ================================================================
c        Angular distribution ./Output/DATA/dNdcos_sm.dat  
c        costheta_i, dN/dcos_i, dN/dcosSI_i, dN/dcosSD_i
c        ================================================================
         write(4,*) '## cos\theta, dN/dcos(SI+SD)[events/kg/year],' 
         write(4,*) '## dN/dcos(SI)[events/kg/year],'
         write(4,*) '## dN/dcos(SD)[events/kg/year]'
         write(4,*) '## smeared  ##'

         do j = 1, cos_theta_bins
            write(4,*) cosmid(j), dNdcos_sm(j), dNdcos_smSI(j), dNdcos_smSD(j)
         enddo         

!  writing out the total rate as a function of days to show annual modulation. R = dN/dt         
c        ================================================================
c        Recoil Rate R        ./Output/DATA/rate_sm.dat 
c        day_i, rate_i, rateSI_i, rateSD_i     
c        ================================================================
         write(5,*) '## cos\theta, rates(SI+SD)[events/kg/year],' 
         write(5,*) '## rates(SI)[events/kg/year],'
         write(5,*) '## rates(SD)[events/kg/year]'
         write(5,*) '## smeared  ##'

         do k = 1, day_bins
            write(5,*) daymid(k), rate_sm(k), rate_smSI(k), rate_smSD(k)
         enddo

         close(2)
         close(3)
         close(4)
         close(5)
         
      endif      


      Return
      End
!-----------------------------------------------------------------------------------------------!
      

!-----------------------------------------------------------------------------------------------!
      subroutine Theory_Simulation(mass,sigma_wnSI,sigma_wpSI,sigma_wnSD,sigma_wpSD,Theory)
!-----------------------------------------------------------------------------------------------!
      Implicit none

      Integer ix, jy, kz
      double precision mass, sigma_wnSI, sigma_wpSI
      double precision sigma_wnSD, sigma_wpSD
      double precision Evalue(Energy_bins+1), cosvalue(cos_theta_bins+1)
      double precision dayvalue(day_bins+1), Emid(Energy_bins)
      double precision cosmid(cos_theta_bins), daymid(day_bins)
      double precision dEnergy, dcos, dday
      double precision Theory(Energy_bins, cos_theta_bins, day_bins)

      include '../include/maddm.inc'
      include '../include/maddm_card.inc'

      dEnergy = (En_max - En_min)/dble(Energy_bins)
      dcos    = (cos_max - cos_min)/dble(cos_theta_bins)
      dday    = (day_max - day_min)/dble(day_bins)
      

! Setting up the center of each bin to calculate the double differential distribution
      

      do ix = 1, Energy_bins+1
         Evalue(ix) = En_min + dble(ix-1)*dEnergy
      enddo

      do ix = 1, Energy_bins
         Emid(ix) = (Evalue(ix) + Evalue(ix+1))/2.d0
      enddo
      

      do jy = 1, cos_theta_bins+1
         cosvalue(jy) = cos_min + dble(jy-1)*dcos
      enddo

      do jy = 1, cos_theta_bins
         cosmid(jy) = (cosvalue(jy) + cosvalue(jy+1))/2.d0 
      enddo

      do kz = 1, day_bins+1
        dayvalue(kz) =  day_min + dble(kz-1)*dday
      enddo

      do kz = 1, day_bins
         daymid(kz) = (dayvalue(kz) + dayvalue(kz+1))/2.d0
      enddo


! Loop to calculate the double differential distribution at the center of each bin
      do ix = 1, Energy_bins
         do jy = 1, cos_theta_bins
            do kz = 1, day_bins
               Theory(ix,jy,kz) = 
     &   d2RdEdcostheta(Emid(ix),cosmid(jy),sigma_wnSI,sigma_wpSI,sigma_wnSD,sigma_wpSD,mass,V_E(daymid(kz)),0)*
     &   dEnergy*dcos*dday
            enddo
         enddo
      enddo

      return
      end

!-----------------------------------------------------------------------------------------------!
      Subroutine GaussmearEnergy(uevents, smevents)
!-----------------------------------------------------------------------------------------------!
!     Smears Events in each bin using Gaussian with resolution delta_E = delE*sqrt(ER)	!
!-----------------------------------------------------------------------------------------------!
      Implicit none
      
      Integer	 i, j, k, niter, nenermax, ix, jy, kz, ncosmax, ndatmax
      parameter(niter = 400, nenermax = 100, ncosmax = 20, ndatmax = 10)
      double precision  sigmawn0, MD, enermid(niter)
      double precision  integ, dE, Emins, Emaxs, dEs
      double precision  func(Energy_bins, cos_theta_bins, day_bins)
      double precision  Evalues(niter+1), Emids(niter), counter
      double precision  Evalue(Energy_bins+1), Emid(Energy_bins), ER(Energy_bins)
      double precision  Array(niter,Energy_bins,cos_theta_bins,day_bins), smevents(niter, ncosmax, day_bins)
      double precision  uevents(1:Energy_bins, 1:cos_theta_bins, 1:day_bins), test(niter)
      
      Include '../include/maddm.inc'
      Include '../include/maddm_card.inc'
      
!     Energy range for smearing.
      
      emins = -200.d0
      emaxs = 200.d0
      
!     Calculating the bin widths.
      dEs = (emaxs - emins)/dble(niter)
      dE = (En_max - En_min)/dble(Energy_bins)
      
      
!     Setting up the midpoints of Energy
      do i = 1, niter+1
         Evalues(i) = emins + dble(i-1)*dEs
      enddo
	
      do i = 1, niter
         Emids(i) = (Evalues(i) + Evalues(i+1))/2.d0
      enddo
      
      do i = 1, Energy_bins+1
         Evalue(i) = En_min + dble(i-1)*dE
      enddo
      
      do j = 1, Energy_bins
         ER(j)   = 0.d0
         Emid(j) = (Evalue(j)+Evalue(j+1))/2.d0
         ER(j)   = Emid(j)
      enddo
      
!cccc Initializing arrays
      do i = 1 , niter
         do ix = 1, Energy_bins
            do jy = 1, cos_theta_bins
               do kz = 1, day_bins
                  smevents(i,jy,kz) = 0.d0
                  func(ix,jy,kz)    = 0.d0
               enddo
            enddo
         enddo
      enddo
      
!     Loop that calculates the smearing. This smearing is done according to Mohlabeng et al, arXiv:1503.03937.
      do i = 1, niter 
         do jy = 1, cos_theta_bins
            do kz = 1, day_bins
               do ix = 1, Energy_bins
                  func(ix,jy,kz)    = uevents(ix,jy,kz)
                  Array(i,ix,jy,kz) = func(ix,jy,kz)*Gaussian(ER(ix),Emids(i),lambd*sqrt(ER(ix)))*dE
                  smevents(i,jy,kz) = smevents(i,jy,kz) + Array(i,ix,jy,kz)     
               enddo
            enddo
         enddo
      enddo
      
      return
      end                       ! End of subroutine.

!-------------------------------------------------------------------------------------------------!
      Subroutine GaussmearAngle(uaevents, smaevents)
!-------------------------------------------------------------------------------------------------!
! This subroutine smears the angle with set angular resolution after energy threshold cut is 	  !
! applied. Input is events unsmeared in angle and output is angular smeared events.               !
!-------------------------------------------------------------------------------------------------!
      Implicit none
      
      Integer	 i, j, k, ncosmax, ix, jy, kz
      parameter(ncosmax = 20)
      double precision  func(Energy_bins,cos_theta_bins,day_bins)
      double precision  thetvalue(cos_theta_bins+1), thetmid(cos_theta_bins)
      double precision  thetmids(ncosmax), thetvalues(ncosmax+1)
      double precision  dthet, dthets, thetmin, thetmax, thetmins, thetmaxs
      double precision  uaevents(Energy_bins,cos_theta_bins,day_bins)
      double precision  smaevents(Energy_bins,ncosmax,day_bins)
      double precision  Array(Energy_bins,ncosmax,cos_theta_bins,day_bins)
      double precision  Array1(Energy_bins,ncosmax,cos_theta_bins,day_bins)
      double precision  Array2(Energy_bins,ncosmax,cos_theta_bins,day_bins)
      double precision  Array3(Energy_bins,ncosmax,cos_theta_bins,day_bins)
      double precision  Array4(Energy_bins,ncosmax,cos_theta_bins,day_bins)
      double precision  Array5(Energy_bins,ncosmax,cos_theta_bins,day_bins)
      
      Include '../include/maddm.inc'
      Include '../include/maddm_card.inc'
      
      thetmin = 0.d0
      thetmax = 180.d0
      
      thetmins = 0.d0
      thetmaxs = 180.d0
      
! Setting up the midpoints of the histograms
      
      dthet = (thetmax-thetmin)/dble(cos_theta_bins)
      dthets = (thetmaxs-thetmins)/dble(ncosmax)
      
      
      do i = 1, cos_theta_bins+1
         thetvalue(i) = thetmin + dble(i-1)*dthet
      enddo
      
      do i = 1, cos_theta_bins
         thetmid(i) = (thetvalue(i) + thetvalue(i+1))/2.d0
      enddo
      
      do i = 1, ncosmax+1
         thetvalues(i) = thetmins + dble(i-1)*dthets
      enddo
      
      do i = 1, ncosmax
         thetmids(i) = (thetvalues(i) + thetvalues(i+1))/2.d0
      enddo
      
!cccc  Initialize arrays
      do i = 1, Energy_bins
         do j = 1 , ncosmax
            do jy = 1, cos_theta_bins
               do kz = 1, day_bins
                  Array1(i,j,jy,kz) = 0.d0
                  Array2(i,j,jy,kz) = 0.d0
                  Array3(i,j,jy,kz) = 0.d0
                  Array4(i,j,jy,kz) = 0.d0
                  Array5(i,j,jy,kz) = 0.d0
                  Array(i,j,jy,kz)  = 0.d0
                  func(i,jy,kz)     = 0.d0
                  smaevents(i,j,kz) = 0.d0
               enddo
            enddo
         enddo
      enddo
!cccccccEnd initialization
      
      do i = 1, Energy_bins
         do jy = 1, cos_theta_bins
            do kz = 1, day_bins
               func(i,jy,kz) = uaevents(i,jy,kz) ! setting input array
            enddo
         enddo
      enddo

      
!    Do the angular smearing from (0 - 180) degrees and do folding.

      do i = 1, Energy_bins 
         do j = 1, ncosmax
            do jy = 1, cos_theta_bins
               do kz = 1, day_bins
                  Array1(i,j,jy,kz) = func(i,jy,kz)*Gaussian(thetmid(jy), thetmids(j), sig_theta)*dthets
                  Array2(i,j,jy,kz) = func(i,jy,kz)*Gaussian(thetmid(jy), -thetmids(j), sig_theta)*dthets
                  Array3(i,j,jy,kz) = func(i,jy,kz)*Gaussian(thetmid(jy), 360.d0-thetmids(j), sig_theta)*dthets
                  Array4(i,j,jy,kz) = func(i,jy,kz)*Gaussian(thetmid(jy), thetmids(j)-360.d0, sig_theta)*dthets
                  Array5(i,j,jy,kz) = func(i,jy,kz)*Gaussian(thetmid(jy), 360.d0+thetmids(j), sig_theta)*dthets
! Adding all Gaussian contributions at each smearing angle.
                  Array(i,j,jy,kz)  = Array1(i,j,jy,kz) + Array2(i,j,jy,kz) + Array3(i,j,jy,kz) +
     &                 Array4(i,j,jy,kz) + Array5(i,j,jy,kz) 
                  smaevents(i,j,kz) = smaevents(i,j,kz) + Array(i,j,jy,kz)
               enddo
            enddo
         enddo
      enddo

      
      return
      end


!--------------------------------------------------------------------------------------!
      Function Gaussian(x, y, sigma)
!--------------------------------------------------------------------------------------!
!     Gaussian distribution function for applying detector resolution/smearing         !
!--------------------------------------------------------------------------------------!
      Implicit none
      
      double precision  x, y, sigma
      
      include '../include/maddm.inc'
      include '../include/maddm_card.inc'        
      
      gaussian = (1.d0/sqrt(2.d0*pi))*(1.d0/sigma)*exp(-((x-y)**2)/(2.d0*(sigma)**2))
      
      return 
      end

!--------------------------------------------------------------------------------------!

!----------------------------------------------------------------------------------------------!      
      Function d2RdEdcostheta(ER, costheta,sigmawn0SI,sigmawp0SI,sigmawn0SD,sigmawp0SD, MD, ve,flag)
!----------------------------------------------------------------------------------------------!

      Implicit none
      Integer mater, mater_comp, n, i, flag
      Integer nIo, nsod, nfl, ncarb, nsul  
      double precision ER, costheta, sigmawn0SI, sigmawp0SI
      double precision MD, ve, diff, double_diff, diffcmp1, diffcmp2
      double precision fA, Z, sigmawn0SD, sigmawp0SD
      double precision Zna, zii, Zcarb, Zfl, Zsul
      double precision avgsn, avgsp, Js, Jsna, JsI, Jsc, Jsfl, Jssul
      double precision avgsnna, avgsnI, avgsnc, avgsnfl, avgsns
      double precision avgspna, avgspI, avgspc, avgspfl, avgsps
      double precision Aflour, AIo, Asod, Acarb, Asulp
      double precision Ana(10), Aii(10), Abna(10), AbI(10), Abfl(10)
      double precision Afl(10), Acarbo(10), Abcarbo(10)
      double precision Asulph(10), Absulph(10)
      double precision A(10), Ab(10)

      include '../include/maddm.inc'
      include '../include/maddm_card.inc'

      mater = material

      do i = 1, 10
         A(i) = 0.d0
         Ab(i)= 0.d0
      enddo
c--------------------------------- single component materials ---------------c
      if (material .le. 10 ) then
         call target_material(mater, A, Ab, Z, n)

         call Spin_matrix(mater, Js, avgsp, avgsn)

         call diff_array(ER,costheta,sigmawn0SI,sigmawp0SI,sigmawn0SD,sigmawp0SD,MD, ve, 
     &    mater, A, Ab, Z, n, Js, avgsp, avgsn, diff, flag)

         d2RdEdcostheta = diff  ! differential rates for single materials.
      endif

c---------------------------------- for NaI ---------------------------------c

      if (material .eq. 11) then ! For NaI (Sodium Iodide)

         do i = 1, 1
            Ana(i) = 0.d0
            Abna(i)= 0.d0
            Aii(i)  = 0.d0
            AbI(i) = 0.d0
         enddo
         
         call target_material(6, Ana, Abna, Zna, nsod)
         call Spin_matrix(6, Jsna, avgspna, avgsnna)

         call diff_array(ER,costheta,sigmawn0SI,sigmawp0SI,sigmawn0SD,sigmawp0SD,MD,ve,
     &    6, Ana, Abna, Zna, nsod, Jsna, avgspna, avgsnna, diffcmp1,flag)
         
         Asod = 0.d0
         do i = 1, nsod
            Asod = Asod + Ana(i)*Abna(i) ! Adding all the Na components weighted by their abundances
         enddo
         
         Asod = Asod/nsod       ! Average of all abundance weighted Na contributions
         
         call target_material(7, Aii, AbI, zii, nIo)
         call Spin_matrix(7, JsI, avgspI, avgsnI)
         
         call diff_array(ER,costheta,sigmawn0SI,sigmawp0SI,sigmawn0SD,sigmawp0SD,MD,ve,
     &    7, Aii, AbI, zii, nIo, JsI, avgspI, avgsnI, diffcmp2,flag)


         do i = 1, nIo
            AIo = AIo + Aii(i)*AbI(i) ! Adding all the I components weighted by their abundances
         enddo
         AIo = AIo/nIo          ! Average of all abundance weighted Na contributions
         
         fA = Asod/(Asod + AIo)
                 
         d2RdEdcostheta = fA*diffcmp1 + (1.d0 - fA)*diffcmp2 ! Adding both contributions to the rates.
         

c------------------------------For CF4 -----------------------------------c

      elseif (material .eq. 12) then ! For CF4 (Carbon Tetraflouride)
         
         do i = 1, 2            ! initialize
            Acarbo(i) = 0.d0
            Abcarbo(i)  = 0.d0
         enddo
         
         call target_material(8, Acarbo, Abcarbo, Zcarb, ncarb)
         call Spin_matrix(8, Jsc, avgspc, avgsnc)

         call diff_array(ER,costheta,sigmawn0SI,sigmawp0SI,sigmawn0SD,sigmawp0SD,MD,ve,
     &    8, Acarbo, Abcarbo, Zcarb, ncarb,Jsc,avgspc,avgsnc,diffcmp1,flag)
         
         Acarb = 0.d0
         do i = 1, ncarb
            Acarb = Acarb + Acarbo(i)*Abcarbo(i) ! Adding all C components weighted by their abundances
         enddo
         Acarb = Acarb/ncarb    ! Average of all abundance weighted C contributions
         
         do i = 1, 1            ! initialize
            Afl(i) = 0.d0
            Abfl(i)  = 0.d0
         enddo
         
         call target_material(9, Afl, Abfl, Zfl, nfl)
         call Spin_matrix(9, Jsfl, avgspfl, avgsnfl)

         call diff_array(ER,costheta,sigmawn0SI,sigmawp0SI,sigmawn0SD,sigmawp0SD,MD,ve,
     &    9, Afl, Abfl, Zfl, nfl, Jsfl, avgspfl, avgsnfl, diffcmp2,flag)
         
         Aflour = 0.d0
         do i = 1, nfl
            Aflour = Aflour + Afl(i)*Abfl(i)
         enddo
         Aflour = Aflour/nfl
         
         fA = Aflour/(Aflour + Acarb)
     
         d2RdEdcostheta = (1.d0 - fA)*diffcmp1 + 4.d0*fA*diffcmp2
          
c--------------------------------For CS2--------------------------------------------c

      elseif( material .eq. 13 ) then ! For CS2 (Carbon Disulphide)
         
         call target_material(8, Acarbo, Abcarbo, Zcarb, ncarb)
         call Spin_matrix(8, Jsc, avgspc, avgsnc)

         call diff_array(ER,costheta,sigmawn0SI,sigmawp0SI,sigmawn0SD,sigmawp0SD,MD,ve,
     &    8, Acarbo, Abcarbo, Zcarb, ncarb,Jsc,avgspc,avgsnc,diffcmp1,flag)

         Acarb = 0.d0
         do i = 1, ncarb
            Acarb = Acarb + Acarbo(i)*Abcarbo(i)
         enddo
         Acarb = Acarb/ncarb
         
         call target_material(10, Asulph, Absulph, Zsul, nsul)
         call Spin_matrix(10, Jssul, avgsps, avgsns)

         call diff_array(ER,costheta,sigmawn0SI,sigmawp0SI,sigmawn0SD,sigmawp0SD,MD,ve,
     &   10, Asulph, Absulph,Zsul,nsul,Jssul,avgsps,avgsns,diffcmp2,flag)
         
         Asulp = 0.d0
         do i = 1, nsul
            Asulp = Asulp + Asulph(i)*Absulph(i)
         enddo
         Asulp = Asulp/nsul
         
         fA = Acarb/(Acarb + Asulp)
         
         d2RdEdcostheta = fA*diffcmp1 + 2.d0*(1.d0 - fA)*diffcmp2
         
      endif
      
      return
      end


!-----------------------------------------------------------------------------------------------!
      Subroutine diff_array(ER,costheta,sigmawn0SI,sigmawp0SI,sigmawn0SD,sigmawp0SD, MD, ve, 
     &  mater, A, Ab, Z, n, Js, avgsp, avgsn, diff, flag)
!-----------------------------------------------------------------------------------------------!
!	This calculates the double differential spectrum for                                    !
!       Dark Matter Directional detection	                                                !
!	Uses the modulation information of the earth inside function V_E                        !
!	Uses the Helm Form Factor in Fhelm and contains Target Isotope information	        !
!       as well as spin dependent form factors. The output is an array of isotopes with         !
!       directional information for each.                                                       !
!	To print out information flag = 1, otherwise flag = 0					!
!-----------------------------------------------------------------------------------------------!
      Implicit none
      
      Integer	i, j, k, flag, n, mater
      double precision E0, MD0, rn, muA, v0, vesc, vearth, AMU, rhoD,c, ve
      double precision  kratio
      double precision MD, day, costheta 
      double precision diff, ER
      double precision fn, fp, Z, mu_p, mu_n
      double precision ap, an, sigmawp0SI, sigmawn0SI, sigmawp0SD, sigmawn0SD
      double precision sigmawpSI, sigmawnSI, sigmawpSD, sigmawnSD
      double precision avgsp, avgsn, Js, formsd, mq, lambbda
      double precision IcSD, FWS

      double precision A(10), Ab(10)
      double precision, Allocatable :: MT(:), r(:), muT(:), R0(:), Ic(:)
      double precision, Allocatable :: vmin(:), formfac(:), s(:), ss(:)
      double precision, Allocatable :: st(:), SigmaSI(:), SigmaSD(:)
      double precision, Allocatable :: formfacSD(:), sigmaWN(:)
      
      include '../include/maddm.inc'
      include '../include/maddm_card.inc'

      mater  = material
      
      call target_material(mater, A, Ab, Z, n) ! Subroutine provides the target material according to the input of the user.

      call Spin_matrix(mater, Js, avgsp, avgsn)


! Allocating arrays to the sizes of the different target materials and their isotopes.
      Allocate (MT(n))
      Allocate (r(n))
      Allocate (muT(n))
      Allocate (R0(n))
      Allocate (Ic(n))

      Allocate (vmin(n))
      Allocate (formfac(n))
      Allocate (formfacSD(n))
      Allocate (s(n))
      Allocate (ss(n))
      Allocate (st(n))
      Allocate (SigmaSI(n))
      Allocate (sigmaSD(n))
      Allocate (sigmaWN(n))
      
!	Some variables
            
      MD0       = MD*(GeV/cs**2)
      E0        = (1.d0/2.d0)*MD0*((v0/c)**2)*GeVtoKeV ! most probable kinetic energy of WIMPs
      
      sigmawnSI = sigmawn0SI*picobarn !WIMP-proton spin independent cross-section in cm^2
      sigmawpSI = sigmawp0SI*picobarn !WIMP-neutron spin independent cross-section in cm^2

      sigmawnSD = sigmawn0SD*picobarn !WIMP-neutron spin dependent cross-section in cm^2
      sigmawpSD = sigmawp0SD*picobarn !WIMP-proton spin dependent cross-section in cm^2

      mu_p      = MD0*M_proton/(MD0 + M_proton) ! reduced mass of Proton-WIMP
      mu_n      = MD0*M_neutron/(MD0 + M_neutron) ! reduced mass of Neutron-WIMP

      v0        = vMP*(km/sec) !Most Probable velocity of WIMPs in DM Halo
      vesc      = vescape*(km/sec) ! Escape velocity of a WIMP from the Galactic Halo.
      vearth    = ve*(km/sec) ! Velocity of the Earth, taken from the V_E function.
      AMU       = 0.932d0*(GeV/cs**2) ! Atomic Mass Units
      RhoD      = rhoDM*GeV*cs**(-2)*cm**(-3) !Density of Dark Matter in our local part of the Galaxy.

      c         = 3.d+5 ! Speed of light in km/s

c form factor contribution terms for the spin independent calculation.
      fn  = sqrt((pi/4.d0)*sigmawnSI*(1.d0/mu_n)**2.d0) ! cm/GeV
      fp  = sqrt((pi/4.d0)*sigmawpSI*(1.d0/mu_p)**2.d0) ! cm/GeV
      
c form factor contribution terms for the spin dependent calculation.
      ap = sqrt(sigmawpSD*((1.d0/G_fermi)**2)*(pi/24.d0)*(1.d0/mu_p)**2.d0) !cm.GeV
      an = sqrt(sigmawnSD*((1.d0/G_fermi)**2)*(pi/24.d0)*(1.d0/mu_n)**2.d0) !cm.GeV

! Normalisation factor for velocity distribution integral.
      kratio  = 1.d0/(erf(vesc/v0) - (2.d0/sqrt(pi))*(vesc/v0)*exp(-(vesc/v0)**2))
      diff    = 0.d0

! interaction for spin dependent WIMP-Nucleus cross-section
      IcSD    = (32.d0/pi)*((G_fermi)**2.d0)*(ap*avgsp + an*avgsn)**2.d0 

      do i = 1, n 
         MT(i)      = 0.932d0*A(i)*(GeV/cs**2) ! Mass of target material. 
         r(i)       = (4.d0*MD0*MT(i))/(MD0 + MT(i))**2 ! Kinematic factor relating WIMP and target masses.
         vmin(i)    = sqrt(ER/(E0*r(i)))*v0 ! Minimum velocity required for recoil.
         muT(i)     = MD0*MT(i)/(MD0 + MT(i)) ! reduced mass of target and WIMP.
         muA        = MD0*AMU/(MD0 + AMU)

! interaction for spin independent WIMP-nucleus cross-section
         Ic(i)      = (4.d0/pi)*(Z*fp + (A(i) - Z)*fn)**2.d0
         
         formfac(i) = Fhelm(ER, A(i))**2 ! Helm Form-factor squared.
!         formfac(i) = (FWS(ER, A(i)))**2.d0 ! Wood-Saxon form-factor squared.
         formfacSD(i) = FormSD(ER, A(i), mater) ! Spin dependent form factor.

         SigmaSI(i) = ((muT(i))**2)*Ic(i)*formfac(i)  ! The spin independent WIMP-nucleus cross-section
         SigmaSD(i) = ((muT(i))**2)*IcSD*((Js+1.d0)/Js)*formfacSD(i) !The spin dependent WIMP-nucleus cross-section

         sigmaWN(i) = SigmaSI(i)+ SigmaSD(i) !Adding SI and SD contributions

! Please look at Lewin & Smith Equation 5.8 for the calculation below.
         R0(i)      = (2.d0/sqrt(pi))*(Naa/A(i))*(RhoD/MD0)*SigmaWN(i)*(v0)*daytosec*detector_size
         
         ss(i)      = kratio*(1.d0/(2.d0*E0*r(i)))*(exp(-((vearth*costheta -vmin(i))/v0)**2))
         st(i)      = (1.d0/(2.d0*E0*r(i)))*kratio*exp(-(vesc/v0)**2)
         
         if((vmin(i) .lt. vesc) .and. (exp(-((vearth*costheta -vmin(i))/v0)**2) .gt. exp(-(vesc/v0)**2)))then
            s(i) = ss(i) - st(i) !Makes sure to calculate the distribution in the correct kinematic regime.
         else
            s(i) = 0.d0
         endif
      
      diff = diff + (Ab(i)/100.d0)*R0(i)*s(i)  
      enddo
            
      
      if(flag .eq. 1) then
         write(*,*) '-----------------------------------'
         write(*,*) 'DM Mass is', MD0,'GeV/c^2'
         write(*,*) 'Recoil Energy :', ER,'keV'         
         do i = 1, n
            write(*,*) 'Spin Independent cross-section:', sigmaSI(i)/picobarn ,'pb'
            write(*,*) 'Spin dependent cross-section  :', sigmaSD(i)/picobarn ,'pb'
         enddo
         write(*,*) 'an =', an*cmtoinvGeV
         write(*,*) 'an =', an*cmtoinvGeV
         write(*,*) 'fp =', fp, 'cm/GeV/c^2'
         write(*,*) 'fn =', fn, 'cm/GeV/c^2'
         write(*,*) 'E0 =', E0, 'keV'
         write(*,*) 'Detector Size', detector_size, 'kg'
         write(*,*) '-----------------------------------'         
         do i = 1, n
            write(*,*) 'A =', real(A(i))
            write(*,*) 'Abundance =', real(Ab(i))
         enddo
         write(*,*) '-----------------------------------'         
         write(*,*) 'Velocity of the Earth =', ve/kmtocm, 'km/s'
         write(*,*) 'WIMP escape Velocity = ', vesc/kmtocm, 'km/s' 
         write(*,*) 'v0 =', v0/kmtocm, 'km/s'
         write(*,*) 'local DM density rho_0 =', RhoD, 'GeV/c^2/cm^-3'
         write(*,*) '-----------------------------------'
         do i = 1, n
            write(*,*) 'M_Target, mu_Target:', MT(i), muT(i), 'GeV/c^2' 
            write(*,*) 'r =', r(i)
         enddo
         write(*,*) 'WIMP-Proton reduced mass =', mu_p
         write(*,*) 'WIMP-Neutron reduced mass =', mu_n
         write(*,*) '-----------------------------------'
         do i = 1, n
            write(*,*) 'SI Form_factor^2, q:', formfac(i)
            write(*,*) 'SD Form_factor^2, q,:', formfacSD(i)
         enddo
         write(*,*) '------------------------------------'
         write(*,*) 'Distribution Normalization =', kratio
         do i = 1, n
            write(*,*) 'R0 =', R0(i), 'kg^-1.KeV^-1.yr^-1'
            write(*,*) 'Vmin =', vmin(i)/kmtocm, 'km/s'
         enddo
      endif
      
      
      return
      end
!-----------------------------------------------------------------------------------------!




!-----------------------------------------------------------------------------------------!
      subroutine target_material(material, A, Ab, Z, n)
!-----------------------------------------------------------------------------------------!
!     Contains all the information for target materials                                   !
!-----------------------------------------------------------------------------------------!

      Implicit none

      Integer n, i, material
      Integer nXe, nGe, nSi, nAr, nNa, nI, nNe, nF, nC, nS
      parameter(nxe=9, nGe=5, nSi=3, nAr=2, nNa=1, nI=1, nNe=3, nC=2)
      parameter(nF=1, nS=4)
      double precision Z_Xe, Z_Ge, Z_Si, Z_Ar, Z_Na, Z_I, Z_Ne, Z_C, Z_S, Z_F
      double precision A_Xe(nXe), A_Si(nSi), A_Ge(nGe), A_Ar(nAr), A_Na(nNa)
      double precision A_Ne(nNe), A_I(nI), A_F(nF), A_S(nS), A_C(nC)
      double precision Ab_Xe(nXe), Ab_Ge(nGe), Ab_Si(nSi), Ab_Ar(nAr), Ab_Na(nNa)
      double precision Ab_Ne(nNe), Ab_I(nI), Ab_F(nF), Ab_S(nS), Ab_C(nC)
      double precision A(10), Ab(10), Z


!----------------------------------------Target information---------------------------!

! Xenon
      A_Xe  = (/ 123.9d0, 125.9d0, 127.9d0, 128.9d0, 129.9d0, 130.9d0, 131.9d0, 133.9d0, 135.9d0 /)
      Ab_Xe = (/ 0.09d0, 0.09d0, 1.92d0, 25.44d0, 4.08d0, 21.18d0, 26.89d0, 10.44d0, 8.87d0 /)
      Z_Xe  = 54.d0
! Germanium
      A_Ge  = (/ 69.92d0, 71.92d0, 72.92d0, 73.92d0, 75.92d0 /)
      Ab_Ge = (/ 20.84d0, 27.54d0, 7.73d0, 36.28d0, 7.61d0 /)
      Z_Ge  = 32.d0
! Silicon
      A_si  = (/ 27.98d0, 28.98d0, 29.97d0 /)
      Ab_si = (/ 92.23d0, 4.68d0, 3.087d0 /)
      Z_Si  = 14.d0
! Neon
      A_Ne  = (/ 19.99d0, 20.99d0, 21.99d0 /)
      Ab_Ne = (/ 90.48d0, 0.27d0, 9.25d0 /)
      Z_Ne  = 10.d0
! Argon
      A_Ar  = (/ 37.962, 39.962 /)
      Ab_Ar = (/ 0.0632, 99.6 /)
      Z_Ar  = 18.d0
! Sodium 
      A_Na  = (/ 22.9d0 /)
      Ab_Na = (/ 100.d0 /)
      Z_Na  = 11.d0
! Iodine
      A_I   = (/ 126.9d0 /)
      Ab_I  = (/ 100.d0 /)
      Z_I   = 53.d0
! Carbon
      A_C   = (/ 12.d0, 13.d0 /)
      Ab_C  = (/ 98.89d0, 1.11d0 /)
      Z_C   = 6.d0
! Flourine
      A_F   = ( / 18.998d0 /)
      Ab_F  = ( / 100.d0 /)
      Z_F   = 9.d0
! Sulphur
      A_S   = (/ 31.9d0, 32.97d0, 33.96d0, 35.96d0 /)
      Ab_S  = (/ 94.9d0, 0.76d0, 4.29d0, 0.02d0 /)
      Z_S   = 16.d0
!---------------------------------------------------------------------------------------!

      if (material .eq. 1) then
         n = nXe

         do i = 1, nXe
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nXe
            A(i)  = A_Xe(i)
            Ab(i) = Ab_Xe(i)
         enddo
         Z  = Z_Xe

      elseif(material .eq. 2) then 
         n = nGe

         do i = 1, nGe
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nGe
            A(i)  = A_Ge(i)
            Ab(i) = Ab_Ge(i)
         enddo
         Z  = Z_Ge

      elseif(material .eq. 3) then
         n = nSi

         do i = 1, nSi
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nSi
            A(i)  = A_Si(i)
            Ab(i) = Ab_Si(i)
         enddo
         Z  = Z_Si

      elseif(material .eq. 4) then
         n = nAr

         do i = 1, nAr
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nAr
            A(i)  = A_Ar(i)
            Ab(i) = Ab_Ar(i) 
         enddo
         Z  = Z_Ar

      elseif(material .eq. 5) then
         n = nNe

         do i = 1, nNe
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nNe
            A(i)  = A_Ne(i)
            Ab(i) = Ab_Ne(i)
         enddo
         Z  = Z_Ne         

      elseif(material .eq. 6) then
         n = nNa

         do i = 1, nNa
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nNa
            A(i)  = A_Na(i)
            Ab(i) = Ab_Na(i)
         enddo
         Z  = Z_Na

      elseif(material .eq. 7) then
         n = nI

         do i = 1, nI
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nI
            A(i)  = A_I(i)
            Ab(i) = Ab_I(i)
         enddo
         Z  = Z_I

      elseif(material .eq. 8) then
         n = nC

         do i = 1, nC
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nC
            A(i)  = A_C(i)
            Ab(i) = Ab_C(i)
         enddo
         Z  = Z_C

      elseif(material .eq. 9) then
         n = nF

         do i = 1, nF
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nF
            A(i)  = A_F(i)
            Ab(i) = Ab_F(i)
         enddo
         Z  = Z_F

      elseif(material .eq. 10) then
         n = nS

         do i = 1, nS
            A(i)  = 0.d0
            Ab(i) = 0.d0
         enddo

         do i = 1, nS
            A(i)  = A_S(i)
            Ab(i) = Ab_S(i)
         enddo
         Z  = Z_S

      endif


      return
      end

!-----------------------------------------------------------------------------------------!
      Real*8 Function V_E(days)
!-----------------------------------------------------------------------------------------!
!     This function calculates the velocity of the earth from the galactic velocity and !
!     sun proper motion.                                                                !
!-----------------------------------------------------------------------------------------!
      Implicit none
      
      Real*8	 pi, degtorad, cosbx, cosby, cosbz, g, lambda, Lx, Ly, Lz, L0
      Real*8	 V_E1(3), days, ecc, uearthave, u_E, vuGal(3), vuSun(3), vu_E(3)
      
      vuGal = (/ 0.d0, 230.d0, 0.d0 /) ! velocity of the Galactic rotation
      vuSun = (/ 10.d0, 13.d0, 7.d0 /) ! Sun's proper motion relative to nearby stars
      
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
      
!--------------------------------------------------------------------------------------- -!
      Real*8 Function FHelm(ER, A)
!-----------------------------------------------------------------------------------------!
!     The Helm form factor as a function of recoil energy ER and target atomic mass A   !
!     We use Lewin & Smith as main reference                                            !
!-----------------------------------------------------------------------------------------!
      Implicit none
      
      Integer	 i 
      double precision   ca, cs, pi, A, ER, q, rn, cc, qr, FF
      
      ca  = 0.52d0              ! Eq (4.10)
      cs  = 0.9d0
      pi  = 4.d0*Atan(1.d0)
      cc  = 1.23d0*(A**(1.d0/3.d0)) - 0.6d0
      rn  = sqrt(cc**2 + (7.d0/3.d0)*(pi**2)*ca**2 - 5.d0*cs**2) ! Eq. (4.11) in Lewin & Smith
!     rn  = 1.14d0*(A**(1.d0/3.d0))   ! Used for Figs 6 and 8 in Lewein & Smith
      q   = 6.92d-3*sqrt(A*ER)
      qr  = q*rn
      
      if (q .eq. 0.d0) then
         FF    = 1.d0
         Fhelm = 1.d0
      else
         FF    = (3.d0*(sin(qr) - qr*cos(qr)))/(qr**3)
         Fhelm = FF*exp(-((q*cs)**2)/2.d0)
      endif
      
      return
      end
      
!----------------------------------------------------------------------------------------!
      Real*8 Function FWS(ER, A)
!-----------------------------------------------------------------------------------------!
!     The Wood-saxon form factor as a function of recoil energy ER and                    ! 
!     target atomic mass A                                                                !
!-----------------------------------------------------------------------------------------!
      Implicit none
      
      Integer	 i 
      double precision   ca, cs, pi, A, ER, q, rn, cc, qr,FF

      cs  = 0.9d0
      pi  = 4.d0*Atan(1.d0)
      cc  = 1.27d0*(A**(1.d0/3.d0))
      rn  = sqrt(cc**2 - 5.d0*cs**2)
      q   = 6.92d-3*sqrt(A*ER)
      qr  = q*rn
      
      if (q .eq. 0.d0) then
         FF    = 1.d0
         FWS = 1.d0
      else
         FF    = (3.d0*(sin(qr) - qr*cos(qr)))/(qr**3)
         FWS = FF*exp(-((q*cs)**2)/2.d0)
      endif
      
      return
      end


!--------------------------------------------------------------------------------------- -!
      Real*8 Function FWSmass(m_chi, A)
!-----------------------------------------------------------------------------------------!
!     The Wood-saxon form factor as a function of WIMP mass and target atomic mass A      ! 
!-----------------------------------------------------------------------------------------!
      Implicit none
      
      Integer	 i 
      double precision   ca, cs, pi, A, ER, q, rn, cc, qr,FF
      double precision v0, MT, m_chi, m_chi0

      v0 = 1.d-3
      MT = 0.932d0*A
      ER = (((m_chi**2)*(v0**2)*MT)/(2*(m_chi + MT)**2))*1.d+6

      cs  = 0.9d0
      pi  = 4.d0*Atan(1.d0)
      cc  = 1.27d0*(A**(1.d0/3.d0))
      rn  = sqrt(cc**2 - 5.d0*cs**2)
      q   = 6.92d-3*sqrt(A*ER)
      qr  = q*rn
      
      if (q .eq. 0.d0) then
         FF    = 1.d0
         FWSmass = 1.d0
      else
         FF    = (3.d0*(sin(qr) - qr*cos(qr)))/(qr**3)
         FWSmass = FF*exp(-((q*cs)**2)/2.d0)
      endif
      
      return
      end



!-----------------------------------------------------------------------------------------!
      Function FormSD(ER, A, mater)
!-----------------------------------------------------------------------------------------!
!     Spin dependent formfactor                                                           !
!-----------------------------------------------------------------------------------------!
      implicit none 

      Integer mater
      double precision A, ER, s00n, s01n, s11n, sq, s0
      double precision s000, s010, s110, a0, a1, delu, deld, dels, formsd
      double precision v0, MT, m_chi, m_chi0

      delu = 0.842d0
      deld = -0.427d0
      dels = -0.085d0


      a0 = (11.d0/18.d0)*(delu + deld) + (5.d0/18.d0)*dels
      a1 = (1.d0/3.d0)*(delu - deld)

      call structure_func(ER, A, mater, s00n, s01n, s11n)

      sq = (a0**2.d0)*s00n + a0*a1*s01n + (a1**2.d0)*s11n

      call structure_func(0.d0, A, mater, s000, s010, s110)

      s0 = (a0**2.d0)*s000 + a0*a1*s010 + (a1**2.d0)*s110

      formSD = sq/s0

      return
      end

!-----------------------------------------------------------------------------------------!
      Function FormSDmass(m_chi, A, mater)
!-----------------------------------------------------------------------------------------!
!     Spin dependent formfactor                                                           !
!-----------------------------------------------------------------------------------------!
      implicit none 

      Integer mater
      double precision A, ER, s00n, s01n, s11n, sq, s0
      double precision s000, s010, s110, a0, a1, delu, deld, dels, formsdmass
      double precision v0, MT, m_chi, m_chi0

      delu = 0.842d0
      deld = -0.427d0
      dels = -0.085d0

      v0 = 1.d-3
      MT = 0.932d0*A
      ER = (((m_chi**2)*(v0**2)*MT)/(2*(m_chi + MT)**2))*1.d+6

      a0 = (11.d0/18.d0)*(delu + deld) + (5.d0/18.d0)*dels
      a1 = (1.d0/3.d0)*(delu - deld)

      call structure_func(ER, A, mater, s00n, s01n, s11n)

      sq = (a0**2.d0)*s00n + a0*a1*s01n + (a1**2.d0)*s11n

      call structure_func(0.d0, A, mater, s000, s010, s110)

      s0 = (a0**2.d0)*s000 + a0*a1*s010 + (a1**2.d0)*s110

      formSDmass = sq/s0

      return
      end

!-----------------------------------------------------------------------------------------!
      subroutine Structure_Func(ER, A, mater, s00, s01, s11)
!-----------------------------------------------------------------------------------------!
!     Calculates the structure functions for Spin dependent form factors                  !
!-----------------------------------------------------------------------------------------!
      implicit none 

      Integer mater, i, ncxe, ncge, ncI, ncNa, ncSI
      parameter(ncxe = 9, ncge = 7, ncI = 9, ncNA = 4, ncSI = 1)
      double precision s00, s11, s01, b, ER, A, q, y, rnbig, rnsmall 
      double precision coeff_s00(10), coeff_s01(10), coeff_s11(10)  
      double precision p01_F, p21_F, Q01_F, Q21_F, J_F, pi 


      pi = 4.d0*Atan(1.d0)

      call coefficients(mater,coeff_s00,coeff_s01,coeff_s11)


      rnbig = 1.7d0*A**(1.d0/3.d0) - 0.28d0 
     &  - 0.78d0*(A**(1.d0/3.d0) - 3.8d0 + sqrt((A**(1.d0/3.d0) - 3.8d0)**2 + 0.2d0)) !fm

      rnsmall = 1.5d0*A**(1.d0/3.d0) !fm

      b = A**(1.d0/6.d0) !fm
      q = 6.92d-3*sqrt(A*ER) !fm^-1

      y = ((q*b)/2.d0)**2.d0

      s00 = 0.d0
      s01 = 0.d0
      s11 = 0.d0

      if (mater .eq. 1) then
         do i = 2, ncxe
            s00 = s00 + exp(-2.d0*y)*(coeff_s00(i)*y**(i-1))
            s01 = s01 + exp(-2.d0*y)*(coeff_s01(i)*y**(i-1))
            s11 = s11 + exp(-2.d0*y)*(coeff_s11(i)*y**(i-1))
         enddo
         s00 = s00 + exp(-2.d0*y)*coeff_s00(1)
         s01 = s01 + exp(-2.d0*y)*coeff_s01(1)
         s11 = s11 + exp(-2.d0*y)*coeff_s11(1)

      elseif(mater .eq. 2) then
         do i = 2, ncge
            s00 = s00 + (coeff_s00(i)*y**(i-1))
            s01 = s01 + (coeff_s01(i)*y**(i-1))
            s11 = s11 + (coeff_s11(i)*y**(i-1))
         enddo
         s00 = s00 + coeff_s00(1)
         s01 = s01 + coeff_s01(1)
         s11 = s11 + coeff_s11(1)

      elseif(mater .eq. 3) then

         do i = 1, ncSi
            s00 = s00 + coeff_s00(i)*Exp(-4.428d0*y)
            s01 = s01 + coeff_s01(i)*Exp(-5.413d0*y)
            s11 = s11 + coeff_s11(i)*Exp(-6.264d0*y)
         enddo

      elseif(mater .eq. 4) then
         s00 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0) ! Gaussian Distribution from arXiv: 0803.2360
         s01 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0)
         s11 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0)

      elseif(mater .eq. 5) then
         s00 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0) ! Gaussian Distribution from arXiv: 0803.2360
         s01 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0)
         s11 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0)

      elseif(mater .eq. 6) then
         do i = 2, ncNa
            s00 = s00 + coeff_s00(i)*y**(i-1)
            s01 = s01 + coeff_s01(i)*y**(i-1)
            s11 = s11 + coeff_s11(i)*y**(i-1)
         enddo
         s00 = s00 + coeff_s00(1)
         s01 = s01 + coeff_s01(1)
         s11 = s11 + coeff_s11(1) 

      elseif(mater .eq. 7) then
         do i = 2, ncI
            s00 = s00 + exp(-2.d0*y)*(coeff_s00(i)*y**(i-1))
            s01 = s01 + exp(-2.d0*y)*(coeff_s01(i)*y**(i-1))
            s11 = s11 + exp(-2.d0*y)*(coeff_s11(i)*y**(i-1))
         enddo
         s00 = s00 + exp(-2.d0*y)*coeff_s00(1)
         s01 = s01 + exp(-2.d0*y)*coeff_s01(1)
         s11 = s11 + exp(-2.d0*y)*coeff_s11(1)  

      elseif(mater .eq. 8) then
         s00 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0) ! Gaussian Distribution from arXiv: 0803.2360
         s01 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0)
         s11 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0)

      elseif(mater .eq. 9) then
         p01_F = 0.1145d0*y**2 - 0.6667d0*y + 1.d0
         p21_F = -0.0026d0*y**2 + 0.0100d0*y
         Q01_F = 0.1088d0*y**2 - 0.6667d0*y + 1.d0
         Q21_F = 0.0006d0*y**2 + 0.0041d0*y
         J_F = 1.d0/2.d0

         s00 = ((2.d0*J_F + 1.d0)/16.d0*pi)*(2.610d0)*(p01_F**2.d0 + p21_F**2.d0)*
     &         Exp(-y)
         s01 = ((2.d0*J_F + 1.d0)/16.d0*pi)*(2.707d0)*(Q01_F*p01_F + Q21_F*p21_F)*
     &         Exp(-y)
         s11 = ((2.d0*J_F + 1.d0)/16.d0*pi)*(2.807d0)*(Q01_F**2.d0 + Q21_F**2.d0)*
     &         Exp(-y)

      elseif(mater .eq. 10) then
         s00 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0) ! Gaussian Distribution from arXiv: 0803.2360
         s01 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0)
         s11 = exp(-((q**2.d0)*(rnbig**2.d0))/4.d0)

      endif

      return
      end

!-----------------------------------------------------------------------------------------!
      subroutine Spin_matrix(mater, J, avg_sp, avg_sn)
!-----------------------------------------------------------------------------------------!
!     provides spin matrix coefficients                                                   !
!-----------------------------------------------------------------------------------------!
      implicit none

      integer mater
      double precision J, avg_sp, avg_sn
      double precision J_xe, J_ge, J_I, J_F, J_Na, J_Si, J_Ne, J_Ar
      double precision J_C, J_S
      double precision avg_sp_xe, avg_sp_ge, avg_sp_I, avg_sp_Na, avg_sp_Si
      double precision avg_sn_xe, avg_sn_ge, avg_sn_I, avg_sn_Na, avg_sn_Si
      double precision avg_sp_Ar, avg_sp_Ne, avg_sp_C, avg_sp_F, avg_sp_S
      double precision avg_sn_Ar, avg_sn_Ne, avg_sn_C, avg_sn_F, avg_sn_S

      J_xe = 3.d0/2.d0
      J_ge = 9.d0/2.d0
      J_Si = 1.d0/2.d0
      J_Ar = 1.d0/2.d0
      J_Ne = 1.d0/2.d0
      J_Na = 3.d0/2.d0
      J_I  = 5.d0/2.d0
      J_C  = 1.d0/2.d0
      J_F  = 1.d0/2.d0
      J_S  = 1.d0/2.d0



!from the Nijmegen II fits Arxiv/0406218v1 
      avg_sp_xe = -0.012d0 
      avg_sp_ge = 0.030d0
      avg_sp_Si = -0.002d0
      avg_sp_Ar = 0.d0
      avg_sp_Ne = 0.02d0    !Arxiv/0406218 EOGM(gA/gV = 1)
      avg_sp_Na = 0.2477d0
      avg_sp_I  = 0.354d0
      avg_sp_C  = -0.009d0  !Arxiv/0406218 EOGM(gA/gV = 1)
      avg_sp_F  = 0.415d0   !Arxiv/0406218 EOGM(gA/gV = 1)
      avg_sp_S  = 0.d0

      avg_sn_xe = -0.217d0
      avg_sn_ge = 0.378d0
      avg_sn_Si = 0.13d0
      avg_sn_Ar = 0.d0
      avg_sn_Ne = 0.294d0  !Arxiv/0406218 EOGM(gA/gV = 1)
      avg_sn_Na = 0.0199d0
      avg_sn_I  = 0.064d0
      avg_sn_C  = -0.172d0 !Arxiv/0406218 EOGM(gA/gV = 1) 
      avg_sn_F  = -0.047d0 !Arxiv/0406218 EOGM(gA/gV = 1)
      avg_sn_S  = 0.d0


      if (mater .eq. 1) then
         J = J_xe
         avg_sp = avg_sp_xe
         avg_sn = avg_sn_xe

      elseif(mater .eq. 2) then 
         J = J_ge
         avg_sp = avg_sp_ge
         avg_sn = avg_sn_ge

      elseif(mater .eq. 3) then 
         J = J_Si
         avg_sp = avg_sp_Si
         avg_sn = avg_sn_Si

      elseif(mater .eq. 4) then 
         J = J_Ar
         avg_sp = avg_sp_Ar
         avg_sn = avg_sn_Ar

      elseif(mater .eq. 5) then 
         J = J_Ne
         avg_sp = avg_sp_Ne
         avg_sn = avg_sn_Ne

      elseif(mater .eq. 6) then 
         J = J_Na
         avg_sp = avg_sp_Na
         avg_sn = avg_sn_Na

      elseif(mater .eq. 7) then
         J = J_I
         avg_sp = avg_sp_I
         avg_sn = avg_sn_I

      elseif(mater .eq. 8) then
         J = J_C
         avg_sp = avg_sp_C
         avg_sn = avg_sn_C

      elseif(mater .eq. 9) then
         J = J_F
         avg_sp = avg_sp_F
         avg_sn = avg_sn_F

      elseif(mater .eq. 10) then
         J = J_S
         avg_sp = avg_sp_S
         avg_sn = avg_sn_S

      endif


      return
      end

!-----------------------------------------------------------------------------------------!
      subroutine coefficients(mater,coeff_s00,coeff_s01,coeff_s11)
!-----------------------------------------------------------------------------------------!
!     provides specific coefficients to calculate the nuclear structure functions         !
!     specific materials                                                                  !
!-----------------------------------------------------------------------------------------!
      implicit none 

      integer mater, i, ncxe, ncge, ncI, ncna, ncsi
      integer ncC, ncNe, ncF, ncAr, ncS
      parameter(ncxe = 9, ncge = 7, ncI = 9, ncna = 4, ncSI = 1)
      parameter(ncC = 1, ncNe = 1, ncAr = 1, ncF = 1, ncS = 1)
      double precision coeff_s00(10), coeff_s01(10), coeff_s11(10)
      double precision coeff_Xe_s00(ncxe), coeff_Xe_s11(ncxe), coeff_Xe_s01(ncxe)
      double precision coeff_I_s00(nci), coeff_I_s11(nci), coeff_I_s01(nci)
      double precision coeff_Ge_s00(ncge), coeff_Ge_s11(ncge), coeff_Ge_s01(ncge)
      double precision coeff_Na_s00(ncna), coeff_Na_s11(ncna), coeff_Na_s01(ncna)
      double precision coeff_Si_s00(ncSi), coeff_Si_s11(ncSi), coeff_Si_s01(ncSi)


!  for Xenon 131 from the Nijmegen II experimental fit.
      coeff_Xe_s00 = (/0.0277344, -0.124487, 0.328287, -0.481399, 0.475646, -0.285177, 0.0968193, -0.01709, 0.001237/)
      coeff_Xe_s11 = (/0.0223447, -0.12206, 0.31949, -0.466949, 0.428767, -0.236789, 0.0740837, -0.011966, 0.000780/)
      coeff_Xe_s01 = (/-0.0497844, 0.247247, -0.632306, 0.896416, -0.816445, 0.45235, -0.142686, 0.0233463, -0.001562/)

!  for Iodine 127 from the Nijmegen II experimental fit.      
      coeff_I_s00 = (/0.11663, -0.572149, 1.33797, -1.72517, 1.37742, -0.669986, 0.190522, -0.0291803, 0.0019081/)
      coeff_I_s11 = (/0.05527, -0.303825, 0.794783, -1.17027, 1.06373, -0.571342, 0.172197, -0.0266165, 0.00166238/)
      coeff_I_s01 = (/0.16205, -0.836288, 2.05944, -2.83193, 2.39726, -1.21214, 0.348612, -0.0521813, 0.00320731/)

!  for Germanium 73 from the Nijmegen II experimental fit.      
      coeff_Ge_s00 = (/0.1606, -1.1052, 3.2320, -4.9245, 4.1229, -1.8016, 0.3211/)
      coeff_Ge_s11 = (/0.1164, -0.9228, 2.9753, -4.8709, 4.3099, -1.9661, 0.3624/)
      coeff_Ge_s01 = (/-0.2736, 2.0374, -6.2803, 9.9426, -8.5710, 3.8310, -0.6948/)

!  for Silicon 29 from (Arxiv:0608097v1)
      coeff_Si_s00 = (/0.00818/)
      coeff_Si_s01 = (/-2.06/)
      coeff_Si_s11 = (/1.06/)

!  for Sodium 23 from Nijmegen II fit
      coeff_Na_s00 = (/0.0380, -0.1743, 0.3783, -0.3430/)
      coeff_Na_s01 = (/0.0647, -0.3503, 0.9100, -0.9858/)
      coeff_Na_s11 = (/0.0275, -0.1696, 0.5077, -0.6180/)


      if (mater .eq. 1) then
         do i = 1, ncxe
            coeff_s00(i) = coeff_Xe_s00(i)
            coeff_s01(i) = coeff_Xe_s01(i)
            coeff_s11(i) = coeff_Xe_s11(i)
         enddo

      elseif(mater .eq. 2) then
         do i = 1, ncge
            coeff_s00(i) = coeff_Ge_s00(i)
            coeff_s01(i) = coeff_Ge_s01(i)
            coeff_s11(i) = coeff_Ge_s11(i)
         enddo

      elseif(mater .eq. 3) then
         do i = 1, ncSi
            coeff_s00(i) = coeff_Si_s00(i)
            coeff_s01(i) = coeff_Si_s01(i)
            coeff_s11(i) = coeff_Si_s11(i)
         enddo

      elseif(mater .eq. 4) then
         do i = 1, ncAr
            coeff_s00(i) = 1.d0
            coeff_s01(i) = 1.d0
            coeff_s11(i) = 1.d0
         enddo

      elseif(mater .eq. 5) then
         do i = 1, ncNe
            coeff_s00(i) = 1.d0
            coeff_s01(i) = 1.d0
            coeff_s11(i) = 1.d0
         enddo

      elseif(mater .eq. 6) then
         do i = 1, ncna
            coeff_s00(i) = coeff_na_s00(i)
            coeff_s01(i) = coeff_na_s01(i)
            coeff_s11(i) = coeff_na_s11(i)
         enddo

      elseif(mater .eq. 7) then
         do i = 1, ncI
            coeff_s00(i) = coeff_I_s00(i)
            coeff_s01(i) = coeff_I_s01(i)
            coeff_s11(i) = coeff_I_s11(i)
         enddo

      elseif(mater .eq. 8) then
         do i = 1, ncC
            coeff_s00(i) = 1.d0
            coeff_s01(i) = 1.d0
            coeff_s11(i) = 1.d0
         enddo

      elseif(mater .eq. 9) then
         do i = 1, ncF
            coeff_s00(i) = 1.d0
            coeff_s01(i) = 1.d0
            coeff_s11(i) = 1.d0
         enddo

      elseif(mater .eq. 10) then
         do i = 1, ncS
            coeff_s00(i) = 1.d0
            coeff_s01(i) = 1.d0
            coeff_s11(i) = 1.d0
         enddo

      endif


      return
      end


!-------------------------------------------------------------------------------------------!      
      double precision function efficiency(x, flag)
!-------------------------------------------------------------------------------------------!      
!     Function for the WIMP detection efficiency                                            !
!-------------------------------------------------------------------------------------------!      
      Implicit none

      integer flag
      double precision x, eff 

      efficiency = 0d0
      if(flag .eq. 0) then
         efficiency = 1.d0 ! for 100% detection efficiency

      elseif(flag .eq. 1) then
         efficiency = 0.5d0 ! for 50% detection efficiency

      elseif(flag .eq. 2) then
         efficiency = eff(x)*0.5d0 ! for Lux efficiency with 50% NR acceptance rate.

      endif

      end
!-------------------------------------------------------------------------------------------!




!-------------------------------------------------------------------------------------------!      
      double precision function eff(x)
!-------------------------------------------------------------------------------------------!      
!     Function to do interpolation for the nuclear recoil efficiency for Lux                !
!-------------------------------------------------------------------------------------------!      
      Implicit none

      integer N, i, ii
      double precision x, x1, x2, x3, x4,xmin, xmax, xi, xii
      double precision xn(100),g(100)

      xmin = 0.d0
      xmax = 100.d0
      N = 52 ! the number of points from Dexter


cc read data file
      open(1,file='./Output/DATA/Lux_Efficiency',status='unknown')
      do i=1,N,1
         read(1,*) xn(i), g(i)
      enddo
      close(1)


cc find xi and xii such that  xi <  'x' < xii
      do i=1, N, 1
         
         xi  = xn(i)
         xii = xn(i+1)

         if ( x .ge. xi .and. x .lt. xii ) then
            ii = i
         endif

      enddo
c      print*, x, ii
      
cc check whether 'x' is close to the left side
cc interpolate with 3 points from the left side
      if( ii .lt. 2 .or. ii .eq. 2) then

        x1 = xn(1) 
        x2 = xn(2)
        x3 = xn(3)
        eff = g(1)*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) 
     &         + g(2)*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3))
     &         + g(3)*(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))

cc check whether 'x' is close to the right side
cc interpolate with 3 points from the right side
      else if( ii .gt. N-1 ) then

        x1 = xn(N-2)
        x2 = xn(N-1)
        x3 = xn(N)
        eff = g(N-2)*(x - x2)*(x - x3)/((x1 - x2)*(x1 - x3)) 
     &         + g(N-1)*(x - x1)*(x - x3)/((x2 - x1)*(x2 - x3))
     &         + g(N)  *(x - x1)*(x - x2)/((x3 - x1)*(x3 - x2))
cc otherwise, interpolate with 4 points
      else

           x1 = xn(ii-1)
           x2 = xn(ii)
           x3 = xn(ii+1)
           x4 = xn(ii+2)
       eff =g(ii-1)*(x-x2)*(x-x3)*(x-x4)/((x1-x2)*(x1-x3)*(x1-x4))
     &        +g(ii  )*(x-x1)*(x-x3)*(x-x4)/((x2-x1)*(x2-x3)*(x2-x4)) 
     &        +g(ii+1)*(x-x1)*(x-x2)*(x-x4)/((x3-x1)*(x3-x2)*(x3-x4))
     &        +g(ii+2)*(x-x1)*(x-x2)*(x-x3)/((x4-x1)*(x4-x2)*(x4-x3))

      endif



      end
c ---------------------------------------------------------------------

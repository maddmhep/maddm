!-----------------------------------------------------------------------------------------------------!
        Subroutine directional_detection_electron(mdm,ioniz_amplitude)
!-----------------------------------------------------------------------------------------------------!
!       The main subroutine for the calculation of the direct and directional detection               !
!       rates d^4 R / (dE dq dcostheta dday)                                                          !
!-----------------------------------------------------------------------------------------------------!

        Implicit none

        include '../include/maddm.inc'
        include '../include/maddm_card.inc'
        include '../include/dm_info.inc'

        double precision mdm
        integer, parameter gridsize = 100
        real(kind=10) ioniz_amplitude(gridsize,gridsize)
        double precision Theory(gridsize,gridsize,cos_theta_bins,day_bins)


        call Theory_Simulation(mdm, ioniz_amplitude, gridsize, Theory)

c       HERE PRINT THE DIFFERENTIAL ARRAY INTO /output
c       (do it once evrything works fine)


        Return

        End





!-----------------------------------------------------------------------------------------------!
        subroutine Theory_Simulation(mdm, ioniz_amplitude, gridsize, Theory)
!-----------------------------------------------------------------------------------------------!
        Implicit none

        double precision mdm
        integer gridsize
        real(kind=10) ioniz_amplitude(gridsize,gridsize)
        double precision qMin, qMax, kMin, kMax, dcos, dday
        double precision Theory(Energy_bins,cos_theta_bins,day_bins) 

        Integer iE, iq, ic, id
        double precision cosvalue(cos_theta_bins+1), dayvalue(day_bins+1)
        double precision q(gridsize), k_e(gridsize), E_e(gridsize)
        double precision cosmid(cos_theta_bins), daymid(day_bins)
        double precision dcos, dday
        double precision Theory(Energy_bins, cos_theta_bins, day_bins)

        include '../include/maddm.inc'
        include '../include/maddm_card.inc'

c       Range used in the computation of the ioniz_amplitude
        qMin = 1.0E-7
        qMax = 1.0E-3

        kMin = 1.0E-7
        kMax = 1.0E-4

c       Angle and time range
        dcos = (cos_max - cos_min)/dble(cos_theta_bins)
        dday = (day_max - day_min)/dble(day_bins)
        

c       Setting up each bin to calculate the double differential distribution
        
c       Setting up array of the momentum k_e and the energy E_e of the final electron and
c       the array of the exchanged momentum q. Both k_e and q need to be written in log scale.
        do iE = 1,gridsize
                k_e(iE) = exp(log(kMin) + (log(kMax) - log(kMin))*dble(iE)/dble(gridsize-1)
                E_e(iE) = sqrt(M_e**2 + k_e(iE)**2)
        enddo

        do iq = 1,gridsize
                q(iq) = exp(log(qMin) + (log(qMax) - log(qMin))*dble(iq)/dble(gridsize-1)
        enddo
        
c       Setting up the values of cos(theta) and day in the centre of the bin
        do ic = 1, cos_theta_bins+1
           cosvalue(ic) = cos_min + dble(ic-1)*dcos
        enddo

        do ic = 1, cos_theta_bins
           cosmid(ic) = (cosvalue(ic) + cosvalue(ic+1))/2.d0 
        enddo

        do id = 1, day_bins+1
          dayvalue(id) =  day_min + dble(id-1)*dday
        enddo

        do id = 1, day_bins
           daymid(id) = (dayvalue(id) + dayvalue(id+1))/2.d0
        enddo


c       Loop to calculate the double differential distribution at the center of each bin
        do iE = 1, gridsize
            do iq = 1, gridsize
                do ic = 1, cos_theta_bins
                    do id = 1, day_bins
                        Theory(iE,iq,ic,id) = 
       &                      d2RdEdcostheta(E_e(iE),q(iq),cosmid(ic),mdm,ioniz_amplitude, gridsize,V_E(daymid(id)),0)
       &                      *dEnergy*dcos*dday
                    enddo
                enddo
            enddo
        enddo

        return
        end





!----------------------------------------------------------------------------------------------!      
        Function d2RdEdcostheta(ER, q, costheta, mdm, ioniz_amplitude, gridsize, ve, flag)
!----------------------------------------------------------------------------------------------!

        Implicit none
        Integer n, i, flag
        double precision ER, q, costheta
        double precision mdm, ve, diff, double_diff, diffcmp1, diffcmp2
        integer gridsize
        real(kind=10) ioniz_amplitude(gridsize,gridsize))

        include '../include/maddm.inc'
        include '../include/maddm_card.inc'


        call diff_array(ER,q,costheta,mdm,ioniz_amplitude, gridsize,ve,diff,flag)

        d2RdEdcostheta = diff  ! differential rates for single materials.
      
        return
        end





!-----------------------------------------------------------------------------------------------!
        Subroutine diff_array(ER, q, costheta, mdm, ioniz_amplitude, gridsize, ve, diff, flag)
!-----------------------------------------------------------------------------------------------!
!	This calculates the double differential spectrum for                                    !
!       Dark Matter Directional detection	                                                !
!	Uses the modulation information of the earth inside function V_E                        !                                               !
!	To print out information flag = 1, otherwise flag = 0				        !
!-----------------------------------------------------------------------------------------------!
        Implicit none

        integer i,j,flag
        parameter(niter = 400)
        double precision pi, c
        double precision mdm, M_dm, M_e, ER, q, E0, v0, vesc, vearth, vmin, rhoDM, RhoD
        double precision N_events, FvDM, r_kin, kNorm, vesc, v0, exp1, exp2
        double precision q(gridsize), k_e(gridsize), E_e(gridsize)
        integer gridsize
        real(kind=10) ioniz_amplitude(gridsize,gridsize))
        real(kind=10) diff_rate(gridsize)


c       Parameters and variables   
        pi      = 4.d0*Atan(1.d0)
        c       = 3.d+5                 ! Speed of light in km/s  

        E0      = (1.d0/2.d0)*MD0*((v0/c)**2)*GeVtoKeV ! most probable kinetic energy of WIMPs
        M_dm    = mdm*(GeV/cs**2)       ! dark matter mass (GeV)
        v0      = vMP*(km/sec)          ! Most Probable velocity of WIMPs in DM Halo
        vesc    = vescape*(km/sec)      ! Escape velocity of a WIMP from the Galactic Halo
        vearth  = ve*(km/sec)           ! Velocity of the Earth, taken from the V_E function
        RhoD    = rhoDM*GeV*cs**(-2)*cm**(-3) ! Density of Dark Matter in our local part of the Galaxy

        r_kin   = (4.d0*M_dm*M_e)/(M_dm + M_e)**2 ! Kinematic factor relating WIMP and target masses.
        kNorm   = 1.d0/(erf(vesc/v0) - (2.d0/sqrt(pi))*(vesc/v0)*exp(-(vesc/v0)**2)) ! Normalisation factor for velocity distribution integral

c       Minimum velocity required for recoil.
        vmin = sqrt(ER/(E0*r_kin))*v0

c       Used to simplify the expression
        exp1 = exp(-((vearth*costheta - vmin)/v0)**2)
        exp2 = exp(-(vesc/v0)**2))
                
c       Dark matter integrated velocity distibution
        if ((vmin.lt.vesc).and.(exp1.gt.exp2) then           ! check the boundaries of v
            FvDM = kNorm*(1.d0/(2.d0*E0*r_kin))*(exp1 - exp2)
        else
           FvDM = 0.d0
        endif

        do i=1,(gridsize)
            diff_rate(i) = 0
        enddo

        do i=1,(gridsize)
c               Integrate q*ioniz_amplitude over q
                do j=1,(gridsize)
                        diff_rate(i) = diff_rate(i) + q(i) * ioniz_amplitude(i,j)
                enddo

c               Multiply for the integrated velocity distribution and for the for the constant factor 
c               in front of the integral of eq. 30.
                diff_rate(i) = RhoD/(128*pi*(M_dm**2)*(M_e**2)) * FvDM * diff_rate(i)
        enddo
        

      if(flag .eq. 1) then
         write(*,*) '-----------------------------------'
         write(*,*) 'DM Mass is', M_dm,'GeV/c^2' 
         write(*,*) 'E0 =', E0, 'keV'
         write(*,*) '-----------------------------------'         
         write(*,*) 'Velocity of the Earth =', ve/kmtocm, 'km/s'
         write(*,*) 'WIMP escape Velocity = ', vesc/kmtocm, 'km/s' 
         write(*,*) 'v0 =', v0/kmtocm, 'km/s'
         write(*,*) 'local DM density rho_0 =', RhoD, 'GeV/c^2/cm^-3'
         write(*,*) '-----------------------------------'
         write(*,*) 'M_Target:', M_e, 'GeV/c^2' 
         write(*,*) 'r_kin =', r_kin
         write(*,*) '------------------------------------'
         write(*,*) 'Distribution Normalization =', kratio
         write(*,*) 'Vmin =', vmin(i)/kmtocm, 'km/s'

      endif
      
      
      return
      end
!-----------------------------------------------------------------------------------------!





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
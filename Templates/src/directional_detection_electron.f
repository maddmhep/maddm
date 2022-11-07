!-------------------------------------------------------------------------------------------------------------!
      Subroutine directional_detection_electron(M_DM,sigma_wnSI,sigma_wpSI,sigma_wnSD,sigma_wpSD,flag,N_events)
        !-----------------------------------------------------------------------------------------------------!
        !     The main subroutine for the calculation of the direct and directional detection                 !
        !     rates R, dN/dE, dN/dcostheta, d2N/dEdcostheta                                                   !
        !-----------------------------------------------------------------------------------------------------!

        Implicit none

        include '../include/maddm.inc'
        include '../include/maddm_card.inc'

        Integer flag
        parameter(niter = 400)
        double precision M_DM, M_e
        double precision sigma_wnSI, sigma_wpSI
        double precision sigma_wnSD, sigma_wpSD
        double precision N_events, FvDM, r_kin, kNorm, vesc, v0

cc      Parameters and variables   
        pi      = 4.d0*Atan(1.d0)
        c       = 3.d+5                 ! Speed of light in km/s  
        
        


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        M_DM    = MD*(GeV/cs**2)        ! CHECK THE THING WITH CS (SHOULD BE C) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        E0      = (1.d0/2.d0)*MD0*((v0/c)**2)*GeVtoKeV ! most probable kinetic energy of WIMPs
        sigmawnSI = sigmawn0SI*picobarn ! WIMP-proton spin independent cross-section in cm^2
        sigmawpSI = sigmawp0SI*picobarn ! WIMP-neutron spin independent cross-section in cm^2
        sigmawnSD = sigmawn0SD*picobarn ! WIMP-neutron spin dependent cross-section in cm^2
        sigmawpSD = sigmawp0SD*picobarn ! WIMP-proton spin dependent cross-section in cm^2
        M_e     = 0.0005109989          ! electron mass (GeV)
        v0      = vMP*(km/sec)          ! Most Probable velocity of WIMPs in DM Halo
        vesc    = vescape*(km/sec)      ! Escape velocity of a WIMP from the Galactic Halo
        vearth  = ve*(km/sec)           ! Velocity of the Earth, taken from the V_E function
        RhoD    = rhoDM*GeV*cs**(-2)*cm**(-3) ! Density of Dark Matter in our local part of the Galaxy


cc      Normalisation factor for velocity distribution integral
        kNorm = 1.d0/(erf(vesc/v0) - (2.d0/sqrt(pi))*(vesc/v0)*exp(-(vesc/v0)**2))

cc      Kinematic factor relating WIMP and target masses
        r_kin = (4.d0*M_DM*M_e)/(M_DM + M_e)**2 

cc      Dark matter integrated velocity distibution 
        if ((vmin(i) .lt. vesc) .and.
     &     (exp(-((vearth*costheta -vmin(i))/v0)**2) .gt. exp(-(vesc/v0)**2))) then

            FvDM = kNorm*(1.d0/(2.d0*E0*r_kin))*(exp(-((vearth*costheta -vmin(i))/v0)**2)
     &                                           - exp(-(vesc/v0)**2))
        else
            FvDM = 0.d0

        endif



        Return
        end subroutine 
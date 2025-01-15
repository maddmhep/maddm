!-------------------------------------------------------------------------------------------!
      function Fform(element, object0, m_darkmatter)
!-------------------------------------------------------------------------------------------!
!   Form factor suppression
!   object0 - 1 - Sun, 2 - Earth, 3 - user
!   element - mass number of the element
!   m_darkmatter - mass of dark matter
!   Ftype - which form factor fit to use (FIX THIS! TO BE IMPLEMENTED)
!-------------------------------------------------------------------------------------------!

      implicit none
      include 'maddm.inc'

      double precision m_darkmatter
      integer object0, element

      include 'capturerate_models.inc'

      if (element.eq.1) then
        Fform = 1.0d0
      else
        Fform = Finf(object0, element) - (1.d0 - Finf(object0, element))
     &  *exp(-1.d0*(log(m_darkmatter)/log(mci(object0, element)))**ali(object0, element))
      endif

      return
      end


!-------------------------------------------------------------------------------------------!
      function Akin(x, object0)
!-------------------------------------------------------------------------------------------!
!     For the kinematic suppression factor
!     object0 - celestial body
!-------------------------------------------------------------------------------------------!

      implicit none

      include '../include/maddm.inc'

      double precision x
      integer object0

      include '../include/capturerate_models.inc'

      Akin = 3.d0/2.d0 * x/ (x - 1.d0)**2 * (vesc0(object0)/ vbar(object0))**2

      return
      end
!-------------------------------------------------------------------------------------------!
      function Skin(object0, element, m_darkmatter, nucleon_mass)
!-------------------------------------------------------------------------------------------!
!     Kinematic suppression factor
!     m_darkmatter - mass of dark matter
!     mn - nucleon mass
!     object0 - celestial body (1 - sun, 2 - earth, ...)
!     element - mass number of the scattering target
!-------------------------------------------------------------------------------------------!

      implicit none

      include '../include/maddm.inc'

      double precision m_darkmatter, nucleon_mass, b,x
      integer element, object0

!     the loop over elements goes over all elements. If they are not present, all arrays
!     are set to 0. This needs to be checked, otherwise a division by 0 might occur

      if (Ai(object0, element).eq.0.d0) then
          Skin = 0.d0
          return
      endif

      b = 1.5d0
      x = m_darkmatter/(Ai(object0, element) * nucleon_mass)

      include '../include/capturerate_models.inc'

      Skin = (Akin(x, object0)**b/(1.d0 + Akin(x, object0)**b))**(1.d0/b)

      return
      end
!-------------------------------------------------------------------------------------------!


!-------------------------------------------------------------------------------------------!
      function Ccap(object0, sigmaSI, sigmaSD_p, sigmaSD_n, m_darkmatter)
!-------------------------------------------------------------------------------------------!
!     Returns the capture rate coefficient of the DM by elements in the Sun/Earth/...       !
!     Obtained from Phys. Rep. 267(1996) 195-373                                            !
!     object0 - celestial body for which the capture rate is being calculated                !
!     sigma_X - spin independent/dependent DM nucleon scattering cross section              !
!-------------------------------------------------------------------------------------------!
      implicit none


      include 'maddm.inc'

      double precision sigmaSI, sigmaSD_p, sigmaSD_n, m_darkmatter
      double precision Cscalar, Cspin
      integer object0, i

      include 'capturerate_models.inc'

      Ccap = 0.d0
      Cscalar = 0.d0
      Cspin = 0.d0

c     Calculate the spin independent contribution
      do i= 1, nelements

           if (Ai(object0, i).ne.0.d0) then

                Cscalar = Cscalar + Fform(i, object0, m_darkmatter)*(sigmaSI/1.d-40)
     &               *fi(object0,i)*phii(object0, i)*Skin(object0, i, m_darkmatter, m_nucleon)
     &               /(Ai(object0, i)*m_nucleon)
           endif
      enddo
      Cscalar = Cscalar*c0(object0)*(rho0(object0)/0.3d0)/(m_darkmatter*(vbar(object0)/270.0d0))

c     and the spin dependent contribution
      do i = 1, 2

           if (Ai(object0, i).ne.0.d0) then

                 Cspin = Cspin + Fform(i, object0, m_darkmatter)*(sigmaSD_n/1.d-40)
     &                *fi(object0,i)*phii(object0, i)*Skin(object0, i, m_darkmatter, m_nucleon)
     &                 /(Zi(object0,i)*m_nucleon)

                 Cspin = Cspin + Fform(i, object0, m_darkmatter)*(sigmaSD_p/1.d-40)
     &                *fi(object0,i)*phii(object0, i)*Skin(object0, i, m_darkmatter, m_nucleon)
     &                /(Zi(object0,i)*m_nucleon)
            endif
      enddo
      Cspin = Cspin*c0(object0)*(rho0(object0)/0.3d0)/(m_darkmatter*(vbar(object0)/270.0d0))

      Ccap = Cscalar + Cspin

      return
      end


c This file contains data about the celestial body objects for the purpose of capture rate calculation
c First index represents the celestial body (1 - Sun, 2 - Earth, 3 - ...)
c Second index represents the mass number of the element
c All conventions and default values taken from Jungman et. al. Phys. Reports 267(1996)

c     Initialization

      Ai = 0.d0
      Zi = 0.d0
      fi = 0.d0
      phii = 0.d0
      Finf = 0.d0
      mci = 0.d0
      ali = 0.d0

c -----------------------------
c      1 - SUN
c -----------------------------

      object_name(1) = 'Sun'
      M0(1) = 1.989d+30  ! mass of the sun in kg
      c0(1) = 4.8d24     ! overall capture coefficient in 1/s
      rho0(1) = rhoDM    ! local DM density in GeV/cm^3
      vesc0(1) = 1156.d0 ! sun's escape velocity  in km/s
      vbar(1) = Vmp      ! mean local dm velocity in km/s

c     Solar composition
c     For Hydrogen
      Ai(1,1)   = 1.d0          !mass number
      Zi(1,1)   = 1.d0          !atomic number
      fi(1,1)   = 0.772d0       !f_i
      phii(1,1) = 3.16d0       !phi_i
      Finf(1,1)  = 1.d0         !Finf_i
      mci(1,1) = 1.d0           !mc_i
      ali(1,1)  = 1.d0          !alpha_i

c     For Helium
      Ai(1,4)   = 4.d0
      Zi(1,4)   = 2.d0
      fi(1,4)   = 0.209d0
      phii(1,4) = 3.4d0
      Finf(1,4)  = 0.986d0
      mci(1,4) = 18.2d0
      ali(1,4)  = 1.58d0

c     For Carbon
      Ai(1,12)   = 12.d0
      Zi(1,12)   = 6.d0
      fi(1,12)   = 3.87d-3
      phii(1,12) = 3.23d0
      Finf(1,12)  = 0.788d0
      mci(1,12) = 61.6d0
      ali(1,12)  = 2.69d0

c     For Nitrogen
      Ai(1,14)   = 14.d0
      Zi(1,14)   = 7.d0
      fi(1,14)   = 9.4d-4
      phii(1,14) = 3.23d0
      Finf(1,14)  = 0.613d0
      mci(1,14) = 75.2d0
      ali(1,14)  = 2.69d0

c     For Oxygen
      Ai(1,16)   = 16.d0
      Zi(1,16)   = 8.d0
      fi(1,16)   = 8.55d-3
      phii(1,16) = 3.23d0
      Finf(1,16)  = 0.613d0
      mci(1,16) = 75.2d0
      ali(1,16)  = 2.69d0

c     For Neon
      Ai(1,20)   = 20.d0
      Zi(1,20)   = 10.d0
      fi(1,20)   = 1.51d-3
      phii(1,20) = 3.23d0
      Finf(1,20)  = 0.613d0
      mci(1,20) = 75.2d0
      ali(1,20)  = 2.69d0

c     For Magnesium
      Ai(1,24)   = 24.d0
      Zi(1,24)   = 12.d0
      fi(1,24)   = 7.39d-4
      phii(1,24) = 3.23d0
      Finf(1,24)  = 0.281d0
      mci(1,24) = 71.7d0
      ali(1,24)  = 2.97d0

c     For Silicon
      Ai(1,28)   = 28.d0
      Zi(1,28)   = 14.d0
      fi(1,28)   = 8.13d-4
      phii(1,28) = 3.23d0
      Finf(1,28)  = 0.281d0
      mci(1,28) = 71.7d0
      ali(1,28)  = 2.97d0

c     For Sulphur
      Ai(1,32)   = 32.d0
      Zi(1,32)   = 16.d0
      fi(1,32)   = 4.65d-4
      phii(1,32) = 3.23d0
      Finf(1,32)  = 0.101d0
      mci(1,32) = 57.d0
      ali(1,32)  = 3.1d0

c     For Iron
      Ai(1,56)   = 56.d0
      Zi(1,56)   = 26.d0
      fi(1,56)   = 1.46d-3
      phii(1,56) = 3.23d0
      Finf(1,56)  = 0.00677d0
      mci(1,56) = 29.3d0
      ali(1,56)  = 3.36d0

c -----------------------------
c      2 - EARTH
c -----------------------------


c -----------------------------
c      3 - USER
c -----------------------------
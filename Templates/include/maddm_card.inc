c print out some status updates of the various calculations
      print_out = .false.

c print out the values of the annihilation cross section at x_f (s = 4m^2/(1-2/x_f))
      print_sigmas = .false.

c flag to use the canonical calculation or the coupled relic density calculation
      relic_canonical = .true.

c flags to determine which calculations to perform
c WARNING: Do not turn on if you haven't already generated a specific process in maddm.
      do_relic_density          = .true.
      do_direct_detection       = .false.
      do_directional_detection  = .false.

c flags for calculating an array of taacs values or calculate them on the spot
      calc_taacs_ann_array = .true.
      calc_taacs_dm2dm_array = .true.
      calc_taacs_scattering_array = .true.

c precision of the ODE integrator which solves for relic density.
      eps_ode = 1.0d-2

c use an approximation for the decoupling temperature in the relic density calc?
      xd_approx = .true.

c if xd_approx = .false. use "step back" parameters for the ODE integrator
c (see manual for more details). Slower, but more accurate.
      x_start = 50.0d0
      x_end = 1.0d3
      dx_step = 1.0d0

c integration precision parameters

c Simpson's rule parameters
c     For the calculation of taacs only, initial number of points in the grid (of dm velocity)
      ngrid_init = 20

c     Number of points to add for one width
c     (the code will add points which exponentially increase in distance from the pole)
      nres_points = 20

c precision of the romberg integration for Wij
      eps_wij = 1.0d-2

c minimum number of iterations in the romberg integration algorithm for both Wij
      iter_wij = 2

c Direct Detection nucleon form factors (fPx - proton, fNx - neutron)

c Scalar FF
      SPu = 0.0153d0
      SPd = 0.0191d0
      SPs = 0.0447d0
      SPg = 1.d0 - SPu - SPd - SPs

      SNu = 0.0110d0
      SNd = 0.0273d0
      SNs = 0.0447d0
      SNg = 1.d0 - SNu - SNd - SNs

c Vector FF
      VPu = 2.0d0
      VPd = 1.0d0

      VNu = 1.0d0
      VNd = 2.0d0

c Axial Vector FF
      AVPu = 0.842d0
      AVPd = -0.427d0
      AVPs = -0.085d0

      AVNu = -0.427d0
      AVNd = 0.842d0
      AVNs = -0.085d0

c Sigma(mu,nu) FF
      SigPu = 0.84d0
      SigPd = -0.23d0
      SigPs = -0.046d0

      SigNu = -0.23d0
      SigNd = 0.84d0
      SigNs = -0.046d0

c For Directional detection and direct detection rates

! choose the target material
! Xenon = 1, Germanium = 2, Silicon = 3, Argon = 4, Neon = 5
! Sodium = 6, Iodine = 7, Carbon = 8, Flourine = 9, Sulphur = 10
! For Compounds: NaI = 11, CF4 = 12, CS2 = 13

      material = 1

c Setting up the DM constants
      vMP      = 220.d0  ! The most probable WIMP velocity v_0 (km/s)
      vescape  = 650.d0  ! The escape velocity of WIMPS in the Halo (km/s)
      rhoDM    = 0.3d0   ! The density of WIMPs in our local Universe (GeV/c^2/cm^-3)

c choose the size of the detector
      detector_size = 1000.d0  ! For 1 ton = 1000 kg detector

      En_threshold  = 4.d0     ! in keVnr

! choose the energy resolution coefficient lambda (sigma_E = lambda*sqrt(E))
      lambd = 1.d0

! choose the angular resolution sigma_theta
      sig_theta = 30.d0

! choose the range for the Energy, angle and time simulations
      En_min         = 0.d0            ! minimum recoil Energy
      En_max         = 100.d0          ! maximum recoil energy

      cos_min        = -1.d0           ! Minimum cosine of recoil angle
      cos_max        = 1.d0            ! Maximum cosine of recoil angle

      day_min        = 0.d0            ! 1st day of observation.
      day_max        = 365.d0          ! 1 year later.

      Energy_bins    = 100             ! Number of Energy bins
      cos_theta_bins = 20              ! Number of cosine angle bins
      day_bins       = 10              ! Number of time bins

! Turn on detector smearing or not (both anergy and angle)
      smearing = .false.

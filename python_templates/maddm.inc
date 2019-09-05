c function names
c integrate.f
      double precision simpson, simpson_adaptive
c init.f
      double precision get_Y_eq, get_gstar, get_Wij_ann, Wij_ann_integrand
      double precision get_Wij_ann_nogrid
c interpolate.f
      integer getindex
      double precision get_Wij_value, get_taacs_value, lineint, bessk, bessk1
c phasespace.f
      double precision gmunu, invmass, lambda, p3dot, p4dot, random
c smatrix.f
      double precision smatrix_ann, smatrix_dm2dm, smatrix_scattering
      double precision smatrix_dd, smatrix_dd_eff, smatrix_dd_tot
c relic_canon.f
      double precision relic_canon, taacs_canon, taacs_integrand_canon
c relic_coupled.f
      double precision relic_coupled, Wij_integrand_other, get_taacs_dm2dm, taacs_dm2dm_integrand
      double precision taacs_coupled, taacs_integrand
c relic_density.f
      double precision relic_density, approximate_xd
c tests.f
      double precision crossijk, cross_check_process
c indirect_detection.f
      double precision sigmav_ID, dm_vel_dist_ID, integrand, taacs_ID, vave_glob,vave_indirect
      integer channel_glob
      character(20) halo_profile
      common/idstuff/ vave_glob, channel_glob, vave_indirect, halo_profile

c direct_detection.f
	  double precision sigma_proton_SI, sigma_neutron_SI, sigma_nucleon
	  double precision sigma_proton_SD, sigma_neutron_SD

c directional_detection.f
      double precision d2RdEdcostheta, Fhelm, V_E, Gaussian

c random seed
      integer idum
      common/randint/ idum

c counting parameters
      integer x1, x2, x3, x4

c generic momentum 4-vectors
      double precision p(0:3), q(0:3)

c model parameters
      integer nparticles, isfermion(%(nb_part)s)
      common/modelints/ nparticles, isfermion
      double precision masses(%(nb_part)s), dofs(%(nb_part)s)
      common/modeldbls/ masses, dofs

c relativistic degrees of freedom
      double precision gstar_x(401), gstar_boson(401), gstar_fermion(401)
      common/gstarstuff/ gstar_x, gstar_boson, gstar_fermion

c more variables for the relativistic degrees of freedom
      double precision g_star, g_star_S, temperature

c old gstar and gstarS stuff
      double precision gstar_temperature(100), gstar(100)
      double precision gstarS_temperature(100), gstarS(100)
      common/oldgstarstuff/ gstar_temperature, gstar, gstarS_temperature, gstarS

c variables needed for the phase space routines
      double precision pf1(0:3), pf2(0:3), pf3(0:3), pf4(0:3), pmass(4)

c dm paraeters
      integer maxdms, max_dm2dm
      parameter (maxdms=%(nb_dm)i, max_dm2dm=%(max_dm2dm)i)
      integer ndmparticles, nvar, dm_index(maxdms), dm_antiindex(maxdms)
      common/dmparms_int/ ndmparticles, nvar, dm_index, dm_antiindex
      integer num_xbar, xbar_index(maxdms), num_not_xbar, not_xbar_index(maxdms), rev_index(maxdms)
      common/xbar_int/ num_xbar, xbar_index, num_not_xbar, not_xbar_index, rev_index
      double precision mdm(maxdms), dof_dm(maxdms), dof_total(maxdms), dof_SM(maxdms,maxdms,max_dm2dm)
      double precision dof_SM_total(maxdms,maxdms,max_dm2dm), xf(maxdms), min_xbar_mass
      common/dmparms/ mdm, dof_dm, dof_total, dof_SM, dof_SM_total, xf, min_xbar_mass
      character(10) dm_names(maxdms), dm_antinames(maxdms)
      common/dmnames/ dm_names, dm_antinames
      logical solved_ODE(maxdms), found_xf(maxdms), dm_sm_scattering(maxdms), print_out, print_sigmas, dmsm_warning
      common/dmparms_logical/ solved_ODE, found_xf, dm_sm_scattering, print_out, dmsm_warning

c relic density integrals
      integer fitsize
      parameter (fitsize=2500)
      integer nWij
      common/num_Wijs/ nWij
      double precision betas(fitsize), beta_step_min, beta_step_max
      common/beta_steps/ betas, beta_step_min,beta_step_max
      double precision Wij_ann(maxdms,maxdms,fitsize), Wij_dm2dm(maxdms,maxdms,max_dm2dm,fitsize)
      double precision Wij_scattering(maxdms,maxdms,max_dm2dm,fitsize), x_taacs(fitsize)
      double precision taacs_ann(maxdms,maxdms,fitsize)
      double precision taacs_dm2dm(maxdms,maxdms,max_dm2dm,fitsize)
      double precision taacs_scattering(maxdms,maxdms,max_dm2dm,fitsize)
      common/dmintegrals/ Wij_ann, Wij_dm2dm, Wij_scattering, x_taacs, taacs_ann, taacs_dm2dm, taacs_scattering

c number of annihilation particles and associated diagrams
      integer dd_tot_num_processes, dd_eff_num_processes, dd_num_processes,num_processes,
     .ann_num_processes, dm2dm_num_processes, scattering_num_processes
      integer ann_nprocesses(maxdms,maxdms), dm2dm_nprocesses(maxdms,maxdms)
      integer scattering_nprocesses(maxdms,maxdms), dm2dm_fs(maxdms,maxdms,max_dm2dm,2)
      common/process_int/ num_processes, ann_num_processes, dm2dm_num_processes,
     . scattering_num_processes, ann_nprocesses, dm2dm_nprocesses, scattering_nprocesses, dm2dm_fs
      character(30) process_names(%(nb_me)i)
      character(30) dd_process_names(%(nb_me_dd)i)
      character(30) dd_eff_process_names(%(nb_me_dd_eff)i)
      character(30) dd_tot_process_names(%(nb_me_dd_tot)i)
      common/process_names/ process_names
      common/dd_process_names / dd_process_names
      common/dd_eff_process_names / dd_eff_process_names
      common/dd_tot_process_names / dd_tot_process_names
	  integer dd_process_ids(%(nb_me_dd)i)
	  integer dd_eff_process_ids(%(nb_me_dd_eff)i)
	  integer dd_tot_process_ids(%(nb_me_dd_tot)i)
	  common/ dd_process_ids/ dd_process_ids
	  common/ dd_eff_process_ids/ dd_eff_process_ids
	  common/ dd_tot_process_ids/ dd_tot_process_ids
      logical ann_process_iden_init(maxdms,maxdms, %(max_dm2sm)i), dm2dm_process_iden_init(maxdms,maxdms,%(max_dm2dm)i)
      logical relic_canonical, calc_taacs_ann_array, calc_taacs_dm2dm_array, calc_taacs_scattering_array
      logical do_relic_density, do_direct_detection, do_indirect_detection, xd_approx, only2to2lo
      common/process_logical/ ann_process_iden_init, dm2dm_process_iden_init, relic_canonical,
     . calc_taacs_ann_array, calc_taacs_dm2dm_array, calc_taacs_scattering_array, xd_approx, only2to2lo

c flags for the test subroutines to see if certain quantities have already been calculated
      logical Wij_ann_calc, Wij_dm2dm_calc, Wij_scattering_calc
      logical taacs_ann_calc, taacs_dm2dm_calc, taacs_scattering_calc, dm_sm_scattering_warn
      common/calc_checks/ Wij_ann_calc, Wij_dm2dm_calc, Wij_scattering_calc,
     . taacs_ann_calc, taacs_dm2dm_calc, taacs_scattering_calc, dm_sm_scattering_warn

c passing variables
      double precision x_pass, beta_pass, x_f, sigmav_xf
      common/passing_dble/ x_pass, beta_pass, x_f, sigmav_xf
      integer dmi, dmj, p_k
      common/passing_int/ dmi, dmj, p_k

c precisions and number of iterations for intergration routines
      double precision eps_ode, eps_taacs, eps_wij, x_start, x_end, dx_step
      common/epses/ eps_ode, eps_taacs, eps_wij, x_start, x_end, dx_step
      integer iter_taacs, iter_wij
      common/iterations/ iter_taacs, iter_wij

c path common block for the ODE integrator
      integer kmax, kount
      double precision dxsav, xp(200), yp(50,200)
      common /path/ kmax, kount, dxsav, xp, yp

c constants
      double precision pi, mpl, ZERO, SPu, SPd, SPs, SPg, SNu, SNd, SNs, SNg
      double precision VPu,VPd,VNu,VNd,AVPu,AVPd,AVPs,AVNd,AVNu,AVNs
      double precision SigPu,SigPd,SigNu,SigNd,SigPs,SigNs
      parameter (pi=3.14159265359d0, mpl=1.2209d19, ZERO=0.d0)

c--------------------------Directional detection-------------------------------

c  Directional detection switch
      logical do_directional_detection
      common/ddm/do_directional_detection

!  Declaration for Target Material
      Integer Material
      common/target/Material


!c Declare WIMP constants
      double precision rhoDM, vescape, vMP
      common/const/rhoDM, vescape, vMP

!c Declare detector specifics
      logical smearing
      double precision detector_size, En_threshold, lambd, sig_theta
      common/detector/detector_size, En_threshold, lambd, smearing, sig_theta

!  Declare simulation variables
      Integer Energy_bins, cos_theta_bins, day_bins
      double precision En_min,En_max,cos_min,cos_max,day_min,day_max, m_proton, m_neutron, m_nucleon
      double precision Naa, G_fermi, g1,sec, GeV, cm, cs, fm, km, kmtocm, barn, picobarn, GeVtokev, pbtocm2,
     & GeVtoMeV, invcmtoGeV, cmtoinvGeV,Mevtokev,daytosec, hrtosec,tontokg,hukgtokg, gevtopb, pbtocm3
      common/simu/Energy_bins,cos_theta_bins,day_bins,En_min,En_max,cos_min,cos_max,day_min,day_max


!  Constants
      parameter(m_proton = 0.938d0)
      parameter(m_neutron = 0.940d0)
      parameter(m_nucleon = 0.939d0)
      parameter(Naa       = 6.023d+26) ! Avogadro's number
      parameter(G_fermi   = 1.166d-5) ! Fermi Constant in (GeV**(-2))
      parameter(g1        = 0.357088d0) ! value of the U(1) gauge coupling constant

! base units
      parameter(sec   = 1.d0)
      parameter(GeV   = 1.d0)
      parameter(cm    = 1.d0)
      parameter(cs    = 1.d0)

! Some important constants and conversion factors.

      parameter(fm         = 1.d-13) ! cm
      parameter(km         = 1.d+5)  ! cm
      parameter(kmtocm     = 1.d+5)  !cm
      parameter(barn       = 1.d-24) ! cm^2
      parameter(picobarn    = 1.d-12*1.d-24) ! cm^2
      parameter(GeVtokev   = 1.d+6)
      parameter(GeVtoMeV   = 1.d+3)
      parameter(invcmtoGeV = 1.98d-14) ! converts cm^-1 to GeV
      parameter(cmtoinvGeV = 5.076d+13) ! converts cm to GeV^-1
      parameter(Mevtokev   = 1.d+3)
      parameter(daytosec   = 24.d0*60.d0*60.d0)
      parameter(hrtosec    = 60.d0*60.d0)
      parameter(tontokg    = 1000.d0)
      parameter(hukgtokg   = 100.d0)
      parameter(gevtopb    = 3.89d+8)
      parameter(pbtocm3    = 2.99d-26)
      parameter(pbtocm2    = 1.00d-36)

c------------------------Capture rate---------------------------------
      logical do_capture

      integer nelements, n_celestial_bodies
      parameter (nelements = 294)
      parameter (n_celestial_bodies = 1)
      double precision capture_rate(n_celestial_bodies), c0(n_celestial_bodies), vbar(n_celestial_bodies)
     &  , vesc0(n_celestial_bodies), rho0(n_celestial_bodies), M0(n_celestial_bodies)
      double precision Ccap, Akin, Skin, Fform

      character(30) object_name(n_celestial_bodies)

      double precision, dimension(n_celestial_bodies, nelements) ::  Ai
      double precision, dimension(n_celestial_bodies, nelements) :: Zi
      double precision, dimension(n_celestial_bodies, nelements) :: fi
      double precision, dimension(n_celestial_bodies, nelements) :: phii
      double precision, dimension(n_celestial_bodies, nelements) :: Finf
      double precision, dimension(n_celestial_bodies, nelements) :: mci
      double precision, dimension(n_celestial_bodies, nelements) :: ali
      common/capture_params/ Ai, Zi, fi, phii, Finf, mci, ali,do_capture,
     & capture_rate, c0, M0, rho0, vbar, object_name

c------------------------ Resonance information --------------------------
! s-channel Resonance locations
      Integer nres, ngrid_init, nres_points, grid_npts
      parameter (nres = %(n_resonances)d)
      parameter (grid_npts = 2000)
      double precision grid(grid_npts), grid_ID(grid_npts)
      double precision resonances(nres)
      double precision resonance_widths(nres)
      double precision beta_res(nres)
      double precision beta_res_width(nres)
      common/resonancevars/resonances, resonance_widths, beta_res, beta_res_width, nres_points, ngrid_init
      common/grids/ grid_ID, grid


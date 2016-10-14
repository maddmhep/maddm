c-------------------------------------------------------------------------c
      program maddm
c-------------------------------------------------------------------------c
c                                                                         c
c  This is the main driver for maddm.  It can be used to call all the     c
c  desired subroutines and functions.                                     c
c                                                                         c
c-------------------------------------------------------------------------c
      implicit none
      include 'maddm.inc'
      include 'input.inc'
      include 'coupl.inc'

c parameters used in this routine only
      Integer sm_flag, prnt_tag
      double precision Oh2
      double precision total_events, sigmawnSI, sigmawpSI, gevtopb
      double precision sigmawnSD, sigmawpSD
      character(len=32) outfilename

c include files generated my the python side that contains necessary information
      include 'maddm_card.inc'
      include 'diagrams.inc'

      include 'resonances.inc'

c Sets all the parameters used in MG/ME to the values in the param card
      call setpara('Cards/param_card.dat')

c If the user wants to change any of the parameters from within the fortran code then
c just change the variables given in 'input.inc' and simply run the subroutine coup()
c
c NOTE: If there are any other model parameters that are dependent on the parameters
c being changed then make sure those parameters are also changed as well
c
c Example:
c      del1 = 0.11d0
c      Mx1 = 500.d0
c      call coup()


c subroutine to initialize the necessary parameters for the calculation
      call init_relic()

c ---- RELIC DENSITY CALCULATION -----
c      call to relic_density will automatically set x_f and sigmav_xf.
	  if (do_relic_density) then
            Oh2 = relic_density(relic_canonical)
	  else
	  		Oh2 = -1d0
	  		x_f = -1d0
	  		sigmav_xf=-1d0
      endif

c ---- SIGMA_NUCLEON CALCULATION -----
c     calls are to sigma_nucleon(proton, spin_independent).
c     proton = 1 is for proton, 0 for neutron
c     spin_independent = 1 for SI and 0 for SD.

      gevtopb = 3.9d+8
      prnt_tag = 0

      if (do_direct_detection) then
      	   sigma_proton_SI = sigma_nucleon(1,1)
      	   sigma_neutron_SI = sigma_nucleon(0,1)
      	   sigma_proton_SD = sigma_nucleon(1,0)
      	   sigma_neutron_SD = sigma_nucleon(0,0)

           sigmawnSI  = sigma_neutron_SI*gevtopb
           sigmawpSI  = sigma_proton_SI*gevtopb
           sigmawnSD  = sigma_neutron_SD*gevtopb
           sigmawpSD  = sigma_proton_SD*gevtopb


           if (do_directional_detection) then

              call directional_detection(mdm(1),sigmawnSI,sigmawpSI,sigmawnSD,
     &        sigmawpSD,prnt_tag,total_events)

           endif
      else
      	   sigma_proton_SI = -1d0
      	   sigma_neutron_SI = -1d0
      	   sigma_proton_SD = -1d0
      	   sigma_neutron_SD = -1d0

c      	   write(*,*) 'SI : ', sigma_proton_SI, sigma_neutron_SI !Antony
c      	   write(*,*) 'SD : ', sigma_proton_SD, sigma_neutron_SD !Antony
      endif

      if (.not. smearing ) then
         sm_flag = 0
      else
         sm_flag = 1
      endif

      if (print_sigmas) then
         write(*,*) 'Contributions to the annihilation cross section at E = mdm(1) / x_f'
	 call cross_check_all_process(mdm(1)/x_f,mdm(1)/x_f,1)
      endif

C      Here write the output.
c	  if (iargc() > 0) then
c	  	call getarg(1, outfilename)
c	    open (unit=33, file= outfilename, action="write", status="replace")
c	  else
c	  	open (unit=33,file="output/maddm.out",action="write",status="replace")
c	  endif

	  write(*,*) 'Omegah^2: ', Oh2
	  write(*,*) 'x_f: ', x_f
          write(*,*) 'Wimp_Mass: ', mdm(1), ' GeV'
	  write(*,*) 'sigmav(xf): ', sigmav_xf
	  write(*,*) 'sigmaN_SI_p: ', sigma_proton_SI, ' GeV^-2',
     &                ':', sigma_proton_SI*gevtopb, ' pb'
	  write(*,*) 'sigmaN_SI_n: ', sigma_neutron_SI,' GeV^-2',
     &                ':', sigma_neutron_SI*gevtopb, ' pb'
	  write(*,*) 'sigmaN_SD_p: ', sigma_proton_SD, ' GeV^-2',
     &                ':', sigma_proton_SD*gevtopb, ' pb'
	  write(*,*) 'sigmaN_SD_n: ', sigma_neutron_SD, ' GeV^-2',
     &                ':', sigma_neutron_SD*gevtopb, ' pb'
          write(*,*) 'Nevents: ', Nint(total_events)
          write(*,*) 'smearing: ', sm_flag

c	  close(33)

c-------------------------------------------------------------------------c
c  Other test functions
c-------------------------------------------------------------------------c

C c bessel function check (min_x, max_x, nsteps, logscale)
C       call bessel_check(0.1d0,100.d0,201,1)
C

C c relativisitic degrees of freedom checks (min_temp, max_temp, nsteps, logscale)
C       call dof_check(0.001d0,100.d0,201,1)
C

C c Wij check (prints out all array of Wij's calculated in relicdensity())
C       call Wij_check()
C

C c taacs test (min_x, max_x, nsteps, logscale)
c       call taacs_check(1.00d0,100d0,101,0)
C

C c check all cross sections for a range of inital state momenta (xinit_min, x_init_max, nsteps, p_or_E, logscale)
C       call cross_check_scan_all(1.d0,1.d3,201,1,1)
C

C c cross section for an individual process over range of inital state momenta
C c (dm_i, dm_j, process_k, xinit_min, xinit_max, nsteps, p_or_E, logscale)
C       call cross_check_scan(1,1,1,1.d0,1.0d3,201,1,1)
C

C c cross-section for all processes at a single initial state momentum or energy (x1init, x2init, p_or_E)
C       call cross_check_all_process(1.d3,1.d3,1)
C

C c cross-section check for a single process at a single initial state momentum or energy
C c (dm_i, dm_j, process_k, x1init, x2init, p_or_E)
C       write(*,*) cross_check_process(1,1,1,1.d3,1.d3,1)
C

C c matrix element check for all processes at a single initial state momentum or energy (x1init, x2init, p_or_E)
C       call matrix_element_check_all(1.0d3,1.0d3,0)
C

C c matrix element check for a single process at a single initial state momentum or energy
C c (dm_i, dm_j, process_k, x1init, x2init, p_or_E)
C       call matrix_element_check_process(1,1,1,1.0d3,1.0d3,0)

cc directional detection tests
c Prints out all the ingredients in the directional detection routine.
c      call d2RdEdcos_test(1.d0,1.d0,sigmawnSI,sigmawpSI,sigmawnSD,sigmawpSD,mDM(1),100.d0)

c Prints out the velocity of the earth in the galactic frame on the 100th day of observation
c      call V_E_test(100.d0, 0)
c Prints out array of earth velocities in the galactic frame, written in ./output/VE_test.dat
c      call V_E_test(100.d0, 1)

c Prints out form factors, both spin dependent and spin independent as function of mass
c      call Form_test(1.d0, 0)
c Prints out array of form factors, both spin dependent and spin independent as functions
c of Energy written in ./output/formfac_test.dat
c      call Form_test(1.d0, 1)

c Prints out form factors, both spin dependent and spin independent as function of mass
c      call Form_testmass(100.d0, 0)
c Prints out array of form factors, both spin dependent and spin independent as functions
c of mass written in ./output/formfacmass_test.dat
c      call Form_testmass(100.d0, 1)

c Calculates the Lux_exclusion limit at 90% C.L. The data is found in x,y,z format in
c ./output/Lux_exclusion.dat
c      call Lux_Exclusion(5.d0,7000.d0,3.d-10,7.d-8,25.d0)



      end

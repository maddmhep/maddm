#! /usr/bin/env python

from init import *
from darkmatter import *
from numpy import loadtxt

#Start the user interface
[model_name, project_name,do_relic_density, do_direct_detection, do_directional_detection] =\
			 initialize_MadDM_session(print_banner = True)

# Create an instance of the base darkmatter class and initialize it
# with the user input
# if the project is already existent, do not generate diagrams etc.
dm = darkmatter()
dm.init_from_model(model_name, project_name)


#If it is a new project...
if dm._new_proj == True:
	# Find the DM candidate and coannihilaton particles
	dm.FindDMCandidate()
	dm.FindCoannParticles()
	dm.GetProjectName()


	print "------ Generating Diagrams ------\n",
	if do_relic_density:
		dm.GenerateDiagramsRelicDensity()
	if do_direct_detection:
		dm.GenerateDiagramsDirDetect()
		if do_directional_detection:
			dm._do_directional_detection = do_directional_detection
	print " Done!"
	print "------ Creating the Numerical Session ------\n",
	dm.CreateNumericalSession()
	print " Done!"

	print "Diagnostics:"
	print dm._projectname
	print dm._paramcard
	print dm._projectpath


	print "------ Calculating Dark matter observables ------ ",
	[omegah2, x_freezeout, wimp_mass, sigmav_xf, sigmaN_SI_proton, sigmaN_SI_neutron, \
	  			 sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch] = dm.Calculate()
	print " Done! \n"
	print "\n\n       RESULTS: \n\n"
	print "NOTE: If a value -1 appears, it means MadDM did not calculate the quantity."
	print "----------------------------------------------"
	print "Dark Matter Mass: ", wimp_mass, " GeV"
	if do_relic_density:
	    print "----------------------------------------------"
            print "Relic Density: %.2e" % omegah2
	    print "x_f: %.2f" % x_freezeout
	    print "sigmav(xf): %.2e GeV^-2    =    %.2e pb" % (sigmav_xf, sigmav_xf*GeV2pb)
	    print "----------------------------------------------"
        if do_direct_detection:
	    print "----------------------------------------------"
	    print "Nucleon cross sections: \n"
	    print "sigma_SI(proton): %.2e GeV^-2    =    %.2e pb" % (sigmaN_SI_proton, sigmaN_SI_proton*GeV2pb)
	    print "sigma_SI(neutron): %.2e GeV^-2   =     %.2e pb"% (sigmaN_SI_neutron, sigmaN_SI_neutron*GeV2pb)
	    print "sigma_SD(proton): %.2e GeV^-2    =    %.2e pb" %(sigmaN_SD_proton, sigmaN_SD_proton*GeV2pb)
	    print "sigma_SD(neutron): %.2e GeV^-2    =    %.2e pb"% (sigmaN_SD_neutron, sigmaN_SD_neutron*GeV2pb)
	    print "Total Number of expected events for a 1 ton detector over a year: ", Nevents
	    print "----------------------------------------------"
	    if do_directional_detection:
                if sm_switch == 0:
	            print "----------------------------------------------"
                    print "Generating Unsmeared Distributions"
	            print "----------------------------------------------"
	            print "found in /output"
	            print "d2NdEdcos.dat"
	            print "dNdE.dat"
	            print "dNdcos.dat"
	            print "rate.dat"
	        else:
	            print "----------------------------------------------"
                    print "Generating smeared Distributions"
	            print "----------------------------------------------"
	            print "found in /output"
	            print "d2NdEdcos_sm.dat"
	            print "dNdE_sm.dat"
	            print "dNdcos_sm.dat"
	            print "rate_sm.dat"
	    else:
               print "----------------------------------------------"
               print " "
	       print "Not simulating WIMP-Nucleus scattering events"
               print " "
               print "----------------------------------------------"

	#Determine whether the point is already ruled out
	LUX_data='./ExpData/LUX_data.dat'
	LUX_x, LUX_y = loadtxt(LUX_data, unpack=True, usecols=[0,1])
	print "--------------------------------------------------------"
	print "Running the exclusion analysis on the parameter point..."
	print "Considering relic density and bound on SI cross section from LUX\n"
	excluded = dm.is_excluded(omega_min = 0., omega_max = 0.1, si_limit_x=LUX_x, si_limit_y=LUX_y)
	if excluded[0]:
		print "       The parameter point is "+bcolors.FAIL+"Excluded."+bcolors.ENDC
		print "		  Excluded by relic density: "+str(excluded[1])
		print "		  Excluded by direct detection: "+str(excluded[2])
	else:
		print "       The parameter point is "+bcolors.OKGREEN+"Allowed."+bcolors.ENDC
	print "--------------------------------------------------------"


#If a project already exists skip to the optional parameter scan.
legit_answer = False
while not legit_answer:
	do_param_scan = raw_input('Would you like to perform a parameter scan?[n] (y/n):')
	if do_param_scan == 'y' or do_param_scan == 'Y':
		legit_answer = True
		param_scan_name = raw_input('Enter the name of your parameter scan: ')
		dm.init_param_scan(param_scan_name)

	elif do_param_scan == 'n' or do_param_scan == 'N' or do_param_scan =='':
		legit_answer = True
		print "Finished! Exiting!"
		exit(0)
	else:
		print "Not a legitimate input. Try again."

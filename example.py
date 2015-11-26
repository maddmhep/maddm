#! /usr/bin/env python
from init import *
from darkmatter import *
from numpy import loadtxt

#Create the relic density object. 
dm=darkmatter()
#Initialize it from the rsxSM model in the MadGraph model folder, 
#and store all the results in the Projects/rsxSM subfolder. 
dm.init_from_model('Higgs-portal_UFO', 'Higgs_portal')


# Determine the dark matter candidate...
dm.FindDMCandidate(prompts=False, dm_candidate='')

#...and all the coannihilation partners with the mass splitting 
# defined by |mX0 - mX1| / mX0 < coann_eps.
#dm.FindCoannParticles(prompts = False, coann_eps = 0.1)

#Get the project name with the set of DM particles and see 
#if it already exists.
dm.GetProjectName()

#Generate all 2-2 diagrams DM annihilation diagrams for relic density.        
dm.GenerateDiagramsRelicDensity()    

#Generate the diagrams for direct detection.
print "Generating direct detection diagrams..."
dm.GenerateDiagramsDirDetect()

#Switch to turn on directional detection and recoil rate calculations
dm._do_directional_detection = True

print "Relic density? -"+str(dm._do_relic_density)
print "Direct detection? - "+str(dm._do_direct_detection)
print "Recoil & directional rates? -"+str(dm._do_directional_detection)

#Print some dark matter properties in the mean time.
print "------ Testing the darkmatter object properties ------"
print "Calc. Omega: "+str(dm._do_relic_density)
print "Calc. DD: "+str(dm._do_direct_detection)
print "DM name: "+dm._dm_particles[0].get('name')
print "DM spin: "+str(dm._dm_particles[0].get('spin'))
print "DM mass var: "+dm._dm_particles[0].get('mass')
print "Mass: "+ str(dm.GetMass(dm._dm_particles[0].get('pdg_code')))+"\n"
print "Project: "+dm._projectname

#Output the FORTRAN version of the matrix elements 
#and compile the numerical code.
dm.CreateNumericalSession()

#Calculate relic density, direct detection ...
[omegah2, x_freezeout, wimp_mass, sigmav_xf, sigmaN_SI_proton, sigmaN_SI_neutron, \
                    sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch]\
			= dm.Calculate()

#See if the model point is exluded by relic density or LUX
LUX_data='./ExpData/LUX_data.dat'
LUX_x, LUX_y = loadtxt(LUX_data, unpack=True, usecols=[0,1])
excluded = dm.is_excluded(omega_min = 0., omega_max = 0.1, si_limit_x=LUX_x, si_limit_y=LUX_y)

print " Done! \n"
print "\n\n       RESULTS: \n\n"	
print "NOTE: If a value -1 appears, it means MadDM did not calculate the quantity."
print "----------------------------------------------"
print "Relic Density: %.2e" % omegah2
print "x_f: %.2f" % x_freezeout
print "sigmav(xf): %.2e GeV^-2    =    %.2e pb" % (sigmav_xf, sigmav_xf*GeV2pb)
print "----------------------------------------------"
print "Nucleon cross sections: \n"
print "sigma_SI(proton): %.2e GeV^-2    =    %.2e pb" % (sigmaN_SI_proton, sigmaN_SI_proton*GeV2pb)
print "sigma_SI(neutron): %.2e GeV^-2   =     %.2e pb"% (sigmaN_SI_neutron, sigmaN_SI_neutron*GeV2pb)
print "sigma_SD(proton): %.2e GeV^-2    =    %.2e pb" %(sigmaN_SD_proton, sigmaN_SD_proton*GeV2pb)
print "sigma_SD(neutron): %.2e GeV^-2    =    %.2e pb"% (sigmaN_SD_neutron, sigmaN_SD_neutron*GeV2pb)
print "----------------------------------------------"
if excluded:
	print "       The parameter point is "+bcolors.FAIL+"Excluded."+bcolors.ENDC
else:
	print "       The parameter point is "+bcolors.OKGREEN+"Allowed."+bcolors.ENDC
print"-----------------------------------------------"
print "All results are written in the output directory of the project!"

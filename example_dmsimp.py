#! /usr/bin/env python
from init import *
from darkmatter import *
from numpy import loadtxt


def isok(passed):
    if passed:
        return bcolors.OKGREEN+"OK!"+bcolors.ENDC
    else:
        return bcolors.FAIL+"FAILED!"+bcolors.ENDC


#Create the relic density object.
dm=darkmatter()
#Initialize it from the rsxSM model in the MadGraph model folder,
#and store all the results in the Projects/rsxSM subfolder.
dm.init_from_model('DMsimp_s_spin0_LO_UFO_wgluons', 'DMsimp_TestSuite')


# Determine the dark matter candidate...
dm.FindDMCandidate(prompts=False, dm_candidate='xd', exclude_particles='xc xr')

#print str(dm._bsm_particles)

#...and all the coannihilation partners with the mass splitting
# defined by |mX0 - mX1| / mX0 < coann_eps.
#dm.FindCoannParticles(prompts = True, coann_eps = 0.0)

#Get the project name with the set of DM particles and see
#if it already exists.
dm.GetProjectName()

#Generate all 2-2 diagrams DM annihilation diagrams for relic density.
dm.GenerateDiagramsRelicDensity()

#Generate the diagrams for direct detection.
#print "Generating direct detection diagrams..."
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
#print "Mass: "+ str(dm.GetMass(dm._dm_particles[0].get('pdg_code')))+"\n"
print "Project: "+dm._projectname

#Output the FORTRAN version of the matrix elements
#and compile the numerical code.
dm.CreateNumericalSession()

#Calculate relic density, direct detection ...

gq = 1.0
gdm = 1.0
gg = 0.0

dm.ChangeParameter('gSXd', gdm)
dm.ChangeParameter('gPXd', 0.0)
dm.ChangeParameter('gPu11', 0.0)
dm.ChangeParameter('gPu22', 0.0)
dm.ChangeParameter('gPu33', 0.0)
dm.ChangeParameter('gPd11', 0.0)
dm.ChangeParameter('gPd22', 0.0)
dm.ChangeParameter('gPd33', 0.0)
dm.ChangeParameter('gSu11', gq)
dm.ChangeParameter('gSu22', gq)
dm.ChangeParameter('gSu33', gq)
dm.ChangeParameter('gSd11', gq)
dm.ChangeParameter('gSd22', gq)
dm.ChangeParameter('gSd33', gq)
dm.ChangeParameter('gSg', gg)

#Xd Xd -> SM SM dominated, non-resonant

dm.ChangeParameter('MXd', 100.0)
dm.ChangeParameter('MY0', 300.0)
dm.ChangeParameter('WY0', 4.96)

[omegah2, x_freezeout, wimp_mass, sigmav_xf, sigmaN_SI_proton, sigmaN_SI_neutron, \
                    sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch]\
			= dm.Calculate()

omegah2_nonres = omegah2
sigmaSI_nonres = (sigmaN_SI_proton + sigmaN_SI_neutron)/2.0*GeV2pb

#Xd Xd -> SM SM dominated, resonant

dm.ChangeParameter('MXd', 145.0)
dm.ChangeParameter('MY0', 300.0)
dm.ChangeParameter('WY0', 1.0)

[omegah2, x_freezeout, wimp_mass, sigmav_xf, sigmaN_SI_proton, sigmaN_SI_neutron, \
                    sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch]\
                        = dm.Calculate()

omegah2_res = omegah2
sigmaSI_res = (sigmaN_SI_proton + sigmaN_SI_neutron)/2.0*GeV2pb

#Xd Xd -> SM SM dominated, resonant, narrow width

dm.ChangeParameter('MXd', 40.0)
dm.ChangeParameter('MY0', 100.0)
dm.ChangeParameter('WY0', 1e-5)

[omegah2, x_freezeout, wimp_mass, sigmav_xf, sigmaN_SI_proton, sigmaN_SI_neutron, \
                    sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch]\
                        = dm.Calculate()


omegah2_res2 = omegah2
sigmaSI_res2 = (sigmaN_SI_proton + sigmaN_SI_neutron)/2.0*GeV2pb

#Xd Xd -> Y0 Y0 dominated

dm.ChangeParameter('MXd', 100.0)
dm.ChangeParameter('MY0', 40.0)
dm.ChangeParameter('WY0', 1.72E-03)

[omegah2, x_freezeout, wimp_mass, sigmav_xf, sigmaN_SI_proton, sigmaN_SI_neutron, \
                    sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch]\
                        = dm.Calculate()

omegah2_dark = omegah2
sigmaSI_dark = (sigmaN_SI_proton + sigmaN_SI_neutron)/2.0*GeV2pb


#Check that the results are ok
nonres_oh2_ok = (2.0*(omegah2_nonres - 9.67e+00)/ (omegah2_nonres + 9.67e+00) < 0.05)
res_oh2_ok = (2.0*(omegah2_res - 3.77e-02)/ (omegah2_res + 3.77e-02) < 0.05)
res2_oh2_ok = (2.0*(omegah2_res2 - 3.77e-02)/ (omegah2_res2 + 3.77e-02) < 0.05)
dark_oh2_ok = (2.0*(omegah2_dark - 5.98e-03)/ (omegah2_dark + 5.98e-03) < 0.05)

nonres_SI_ok = (2.0*(sigmaSI_nonres - 1.551E-08)/ (sigmaSI_nonres + 1.551E-08) < 0.05)
res_SI_ok = (2.0*(sigmaSI_res - 2.877E-07)/ (sigmaSI_res + 2.877E-07) < 0.05)
res2_SI_ok = (2.0*(sigmaSI_res2 - 1.22E-6)/ (sigmaSI_res2 + 1.22E-6) < 0.05)
dark_SI_ok = (2.0*(sigmaSI_dark - 4.906E-05)/ (sigmaSI_dark +  4.906E-05 ) < 0.05)

print "-------------- RESULTS FOR MODEL DMSimp ------------------------------------------- [Oh2 OK?]   [DD(SI) OK?] \n"
print "   (MXd, MY0) = (100, 300) GeV:     Oh2 = %.2e     sigSI = %.2e pb           " % (omegah2_nonres, sigmaSI_nonres),
print '{:^15}'.format(isok(nonres_oh2_ok)),
print '{:^25}'.format(isok(nonres_SI_ok))
print "   (MXd, MY0) = (145, 300) GeV:     Oh2 = %.2e     sigSI = %.2e pb           " % (omegah2_res, sigmaSI_res),
print '{:^15}'.format(isok(res_oh2_ok)),
print '{:^25}'.format(isok(res_SI_ok))
print "   (MXd, MY0) = (40, 100)  GeV:     Oh2 = %.2e     sigSI = %.2e pb           " % (omegah2_res2, sigmaSI_res2),
print '{:^15}'.format(isok(res2_oh2_ok)),
print '{:^25}'.format(isok(res2_SI_ok))
print "   (MXd, MY0) = (100, 40)  GeV:     Oh2 = %.2e     sigSI = %.2e pb           " % (omegah2_dark, sigmaSI_dark),
print '{:^15}'.format(isok(dark_oh2_ok)),
print '{:^25}'.format(isok(dark_SI_ok))

print "----------------------------------------------------------------------------------------------------------"

#See if the model point is exluded by relic density or LUX
#LUX_data='./ExpData/LUX_data.dat'
#LUX_x, LUX_y = loadtxt(LUX_data, unpack=True, usecols=[0,1])
#excluded = dm.is_excluded(omega_min = 0., omega_max = 0.1, si_limit_x=LUX_x, si_limit_y=LUX_y)

#print " Done! \n"
#print "\n\n       RESULTS: \n\n"
#print "NOTE: If a value -1 appears, it means MadDM did not calculate the quantity."
#print "----------------------------------------------"
#print "Relic Density: %.2e" % omegah2
#print "x_f: %.2f" % x_freezeout
#print "sigmav(xf): %.2e GeV^-2    =    %.2e pb" % (sigmav_xf, sigmav_xf*GeV2pb)
#print "----------------------------------------------"
#print "Nucleon cross sections: \n"
#print "sigma_SI(proton): %.2e GeV^-2    =    %.2e pb" % (sigmaN_SI_proton, sigmaN_SI_proton*GeV2pb)
#print "sigma_SI(neutron): %.2e GeV^-2   =     %.2e pb"% (sigmaN_SI_neutron, sigmaN_SI_neutron*GeV2pb)
#print "sigma_SD(proton): %.2e GeV^-2    =    %.2e pb" %(sigmaN_SD_proton, sigmaN_SD_proton*GeV2pb)
#print "sigma_SD(neutron): %.2e GeV^-2    =    %.2e pb"% (sigmaN_SD_neutron, sigmaN_SD_neutron*GeV2pb)
#print "----------------------------------------------"
#if excluded:
#	print "       The parameter point is "+bcolors.FAIL+"Excluded."+bcolors.ENDC
#else:
#	print "       The parameter point is "+bcolors.OKGREEN+"Allowed."+bcolors.ENDC
#print"-----------------------------------------------"
#print "All results are written in the output directory of the project!"

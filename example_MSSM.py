#! /usr/bin/env python
from init import *
from darkmatter import *
from numpy import loadtxt

from shutil import copyfile


def isok(passed):
    if passed:
        return bcolors.OKGREEN+"OK!"+bcolors.ENDC
    else:
        return bcolors.FAIL+"FAILED!"+bcolors.ENDC

#Calculate relic density, direct detection ...


#sps_points = ['1a', '1b', '2','3','4','5','9']
#models=['MSSM_UFO','MSSM_UFO_sps1b','MSSM_UFO_sps2','MSSM_UFO_sps3','MSSM_UFO_sps1b','MSSM_UFO_sps5','MSSM_UFO_sps9']
#dmCandidates=['n1','~n1','n1','~n1','~n1','~n1','~n1']
#coannP=['','sl3+ sl6+','','sl1+ sl2+ sl3+ sl4+ sl5+ sl6+','','','x1+ x2+']

sps_points = ['1a', '1b', '2','3','5','9']
models=['MSSM_UFO','MSSM_UFO_sps1b','MSSM_UFO_sps2','MSSM_UFO_sps3','MSSM_UFO_sps5','MSSM_UFO_sps9']
dmCandidates=['n1','~n1','n1','~n1','~n1','~n1']
coannP=['','~sl3+ ~sl6+','','~sl1+ ~sl2+ ~sl3+ ~sl4+ ~sl5+ ~sl6+','','~x1+ ~x2+']

#correct_oh2 = [0.195, 0.374, 7.860, 0.116, 0.0474, 0.338, 0.00117]
#correct_SI_p = [2.18e-10,3.89e-11, 7.39e-11, 6.02e-11,2.07e-10,4.33e-11,5.03e-11]
#correct_SD_p = [6.56e-6,1.21e-6,9.89e-6,1.64e-6,7.58e-6,4.09e-8,1.12e-6]
#correct_SI_n = [2.17e-10,3.87e-11, 7.39e-11, 5.98e-11,2.05e-10,4.31e-11,4.98e-11]
#correct_SD_n = [8.79e-6,1.47e-6,8.05e-6,2.08e-6,7.82e-6,4.77e-7,1.72e-6]

correct_oh2 = [0.195, 0.374, 7.860, 0.116, 0.338, 0.00117]
correct_SI_p = [2.18e-10,3.89e-11, 7.39e-11, 6.02e-11,4.33e-11,5.03e-11]
correct_SD_p = [6.56e-6,1.21e-6,9.89e-6,1.64e-6,4.09e-8,1.12e-6]
correct_SI_n = [2.17e-10,3.87e-11, 7.39e-11, 5.98e-11,4.31e-11,4.98e-11]
correct_SD_n = [8.79e-6,1.47e-6,8.05e-6,2.08e-6,4.77e-7,1.72e-6]

correct_SI=[]
correct_SD=[]
j=0
while j<len(correct_SI_p):
    correct_SI+=[(correct_SI_p[j]+correct_SI_n[j])/2.0]
    correct_SD+=[(correct_SD_p[j]+correct_SD_n[j])/2.0]
    j+=1
print correct_SI
print correct_SD    
omegah2_results = []
SI_results = []
SD_results = []
oh2_ok = []
SI_ok = []
SD_ok = []
i=0
for sps_pt in sps_points:
    #Create the relic density object.
    dm=darkmatter()
    print 'sps_pt :', sps_pt
    print 'models: ', models[i]
    print 'DM: ', dmCandidates[i]
    #Initialize it from the rsxSM model in the MadGraph model folder,
    #and store all the results in the Projects/rsxSM subfolder.
    dm.init_from_model(models[i], 'MSSM_TestSuite_sps'+sps_pt, overwrite=False)
    
    # Determine the dark matter candidate...
    dm.FindDMCandidate(prompts=False, dm_candidate=dmCandidates[i], exclude_particles='all')
    print "DM name: "+dm._dm_particles[0].get('name')
    #print str(dm._bsm_particles)
    
    #...and all the coannihilation partners with the mass splitting
    # defined by |mX0 - mX1| / mX0 < coann_eps.
    if coannP[i]!='':
        dm.FindCoannParticles(coannP[i], prompts = False)
    #Get the project name with the set of DM particles and see
    #if it already exists.
    dm.GetProjectName()
    
    #Generate all 2-2 diagrams DM annihilation diagrams for relic density.
    dm.GenerateDiagramsRelicDensity()
    
    #Generate the diagrams for direct detection.
    #print "Generating direct detection diagrams..."
    dm.GenerateDiagramsDirDetect()
    
    #Switch to turn on directional detection and recoil rate calculations
    dm._do_directional_detection = False
    
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
    
    print "using SPS"+sps_pt+":"
    print "----------------------"
    #copyfile("SPS_cards/sps"+sps_pt+"_slha.out.slha2","Projects/"+dm._projectname+"/Cards/param_card.dat")

    [omegah2, x_freezeout, wimp_mass, sigmav_xf, sigmaN_SI_proton, sigmaN_SI_neutron, \
                    sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch]\
			= dm.Calculate()

    omegah2_results.append(omegah2)
    SI_results.append((sigmaN_SI_proton + sigmaN_SI_neutron)/2.0*GeV2pb)
    SD_results.append((sigmaN_SD_proton + sigmaN_SD_neutron)/2.0*GeV2pb)
    i+=1

print "------- RESULTS FOR MODEL MSSM -------- [Oh2 OK?]   [DD(SI) OK?]   [DD(SD) OK?] \n"
print SI_results
print SD_results
print omegah2_results
print correct_oh2
for i in range(len(sps_points)):
    oh2_ok.append(isok((2.0*abs((omegah2_results[i] - correct_oh2[i])/ (omegah2_results[i] + correct_oh2[i])) < 0.1)))
    SI_ok.append(isok((2.0*abs((SI_results[i] - correct_SI[i])/ (SI_results[i] + correct_SI[i])) < 0.1)))
    SD_ok.append(isok((2.0*abs((SD_results[i] - correct_SD[i])/ (SD_results[i] + correct_SD[i])) < 0.1)))

for i in range(len(sps_points)):
    print '{:^40}'.format('SPS'+sps_points[i]),
    print '{:^20}'.format(oh2_ok[i]),
    print '{:^25}'.format(SI_ok[i]),
    print '{:^15}'.format(SD_ok[i])

print "----------------------------------------------------------------------------------"

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

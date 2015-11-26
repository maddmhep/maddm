#! /usr/bin/env python

try:
    import numpy as np
except Exception, error:
    logging.debug(error)
    print "ERROR: This parameter scanning script requires numpy, which is not installed"
    sys.exit(1)
import fileinput
import sys
import os
import pickle

#Set the proper paths
maddm_path = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
pieces = maddm_path.split('/')
maddm_path = maddm_path.replace(pieces[-1], '')
sys.path.append(maddm_path)

from init import *
from darkmatter import *

#Read in the DM information
with open('dm_object.pik', 'r') as f:
	dm = pickle.load(f)

#Change the directory to MadDM root folder so that code can run properly.
os.chdir(maddm_path)

#Print out some basic information about the dark matter particle.
print "--------------------------------------------"
print "Model: "+dm._modelname
print "Project Name: "+dm._projectname
print "Project Path: "+dm._projectpath
print "Parameter Card: "+dm._paramcard
print "DM particles: "+dm._dm_particles[0].get('name')
print "DM spin (2s+1): "+str(dm._dm_particles[0].get('spin'))
print "--------------------------------------------"

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#CHANGE THIS PART FOR YOUR PARAMETER SCAN
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#Add as many for loops as you wish to scan over parameters
#NOTE:  The default script varies over two parameters (param1 and param2), and outputs the
#		results in the format "param1    param2    omegah^2 ..." on the screen and in the output file.
#		Please make sure to change the number of "for" loops, calls to ChangeParameter() and the
#		output commands to account for the number of parameters you wish to vary. The critical
#		places where changes are usually required are marked with "<-----Change here ----!!!"


#Define the arrays of values your parameters should take.
#If you are using np.arrange, the format is np.arange(init_val, final_val, step_size)
param1_name = 'parameter1' #"<-----Change here ----!!!"
param2_name = 'parameter2' #"<-----Change here ----!!!"
param1_values = np.arange(0.001, 0.501, 0.01) #"<-----Change here ----!!!"
param2_values = np.arange(0.001, 0.501, 0.01) #"<-----Change here ----!!!"

#File in which to write the output.
outputfile_name = dm._projectpath+"/output/parameter_scan.txt" #"<-----Change here ----!!!"

#Check if the file already exists.
if os.path.exists(outputfile_name):
	legit_answer = False
	while (not legit_answer):
	    answer = raw_input("Output file already exists. Overwrite?[n] (y/n):")
            if ((answer == 'y') or (answer == 'Y')):
			legit_answer = True
	    if ((answer == 'n') or (answer == 'N') or (answer == '')):
			legit_answer = True
			print "Nothing to do... exiting."
			sys.exit(0)

outputfile = open(outputfile_name, "w")

for param1 in param1_values:
	for param2 in param2_values: #"<-----Change here ----!!!" (you can add additional loops)
		#Change the parameter with name param1_name in the param_card.dat to the value param1
		dm.ChangeParameter(param1_name, param1)
		#Change the parameter with name param2_name in the param_card.dat to the value param2
		dm.ChangeParameter(param2_name, param2)
		#Calculate relic density, spin independent/dependent scattering cross sections ...
		[omegah2, x_freezeout, wimp_mass, sigmav_xf, sigmaN_SI_proton, sigmaN_SI_neutron, \
	  			 sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch] = dm.Calculate()

		#Output the results on the screen and in the output file
		print param1, param2, omegah2, x_freezeout, sigmav_xf, sigmaN_SI_proton, sigmaN_SI_neutron, \
	  			 sigmaN_SD_proton, sigmaN_SD_neutron, Nevents #"<-----Change here ----!!!"
		outputstring = " %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e  %.3e \n" % \
				(param1, param2, omegah2, x_freezeout, sigmav_xf, sigmaN_SI_proton*GeV2pb, sigmaN_SI_neutron*GeV2pb, \
	  			 sigmaN_SD_proton*GeV2pb, sigmaN_SD_neutron*GeV2pb, Nevents)#"<-----Change here ----!!!"
		outputfile.write(outputstring)

		#If you are generating the recoil rate distributions, angular distributions etc.
		#the following part of the code will change the output file names so that they
		#are not overwritten. For instance, it will move 'rate.dat' to rate_param1_param2.dat
		if dm._do_directional_detection:
			suffix = "_"+str(param1)+"_"+str(param2)+".dat"
			shutil.move(dm._projectpath+"/output/d2NdEdcos.dat", dm._projectpath+"/output/d2NdEdcos"+suffix)
			shutil.move(dm._projectpath+"/output/dNdE.dat", dm._projectpath+"/output/dNdE"+suffix)
			shutil.move(dm._projectpath+"/output/dNdcos.dat", dm._projectpath+"/output/dNdcos"+suffix)
			shutil.move(dm._projectpath+"/output/rate.dat", dm._projectpath+"/output/rate"+suffix)

outputfile.close()
#---------------------------------------------------------------------------
#-------------------------------------------------------------------------

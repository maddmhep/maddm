import os
import shutil
import sys
import subprocess
import fileinput
import glob
import pickle
import logging
import re

# python files from the madgraph source directory
import madgraph.core.base_objects as base_objects
#import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.interface.master_interface as master_interface
import models.model_reader as model_reader
import models.check_param_card as check_param_card
import madgraph.iolibs.export_v4 as export_v4
import madgraph.core.helas_objects as helas_objects
import madgraph.iolibs.file_writers as writers
import madgraph.various.misc as misc

import madgraph.iolibs.drawing_eps as draw


#Set up logging
logging.basicConfig(filename='maddm.log',level=logging.DEBUG, filemode='w')

# Output routines for MadDM
import MGoutput
try:
        from numpy import interp
except Exception, error:
        logging.debug(error)
        print "ERROR: Numpy is not installed on your system. Please install Numpy to proceed!"
        sys.exit(0)

# Root path
rp = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.append(rp)

# Current directory path
pp = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[1]

sys.path.append(rp+'/'+pp)

#Conversion from GeV^-2 to pb etc.
GeV2pb = 3.894E8
pb2cm2  = 1.0E-36
cm22pb  = 1.0E36

#-------------------------------------------------------------------------#
class darkmatter(base_objects.Particle):
#-------------------------------------------------------------------------#
#                                                                          #
#  This class conains all the routines and functions that are needed to          #
#  do the following actions.                                                  #
#                                                                          #
#  1. Determine the Dm candidate and coannihilation particles of a          #
#          given model.                                                    #
#  2. Generate all the relevant annihilation diagrams that are needed to  #
#          calculate the overall relic abundance of DM and direct detection rates.                          #
#  3. Create a fortran project where all the numerical calculations will  #
#          take place.                                                          #
#                                                                          #
#-------------------------------------------------------------------------#

        #-----------------------------------------------------------------------#
        def __init__(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine initializes some of the lists that are used to find        #
        #  the dark matter candidate and coannihilation partners needed for        #
        #  the relic abundance calculation.                                              #
        #                                                                              #
        #-----------------------------------------------------------------------#

          # Initializes useful variables
          self._bsm_particles = []
          self._bsm_masses = []
          self._bsm_final_states = []
          self._dm_particles = []
          self._coann_particles = []
          self._dm_names_list = []
          self._dm_antinames_list = []
          self._dm_thermal_scattering = []
          self._wanted_lorentz = []
          self._wanted_couplings = []
          self._new_proj = True
          self._projectpath = ''
          self._modelname = ''
          self._projectname = ''
          self._paramcard = ''
          self._do_direct_detection = False
          self._do_indirect_detection = False
          self._do_directional_detection = False
          self._do_relic_density = False
          self._DM_all_names = []
          self._eff_operators_SI = {1:'SIEFFS', 2:'SIEFFF', 3:'SIEFFV'}
          self._eff_operators_SD = {1:False, 2:'SDEFFF', 3:'SDEFFV'}
          self._eff_model_dm_names = {1:'~sdm', 2:'~fdm', 3:'~vdm'}
          self._excluded_particles = []
          self._resonances = []
          self._dm_mass = -1.0
          self._allowed_ID_fs=['a a', 'v v', 've ve~', 'vm vm~', 'vt vt~', 'g g']

          #Results of the calculation
          self._omegah2 = -1.0
          self._x_freezeout = -1.0
          self._sigmav_xf = -1.0
          self._sigmaN_SI_proton = -1.0
          self._sigmaN_SI_neutron=-1.0
          self._sigmaN_SD_proton=-1.0
          self._sigmaN_SD_neutron=-1.0
          self._Nevents=-1.0

        #-----------------------------------------------------------------------#
        def init_from_model(self,modelname, projectname, overwrite=False):
        #-----------------------------------------------------------------------#
        #                                                                            #
        #  Given a user input model this routine initializes the project            #
        #  folder and initializes several class variables used in other         #
        #  routines.                                                                   #
        #                                                                            #
        #-----------------------------------------------------------------------#

          # Set up the model and project name
          self._modelname = modelname
          self._projectname = projectname

          #First check if the project already exists. If it does, make sure to ask
          #before deleting it.
          project_list = glob.glob('Projects/'+projectname+'_*')
          if (project_list == []):
                new_project = True
          else:
                new_project = False

          if (new_project == False and (not overwrite)):
                    legit_answer = False
                    while (not legit_answer):
                          answer = raw_input("Project directory already exists. Overwrite?[n] (y/n):")
                          if ((answer == 'y') or (answer == 'Y')):
                                new_project = True
                                legit_answer = True
                                for project in project_list:
                                  shutil.rmtree(project)
                          if ((answer == 'n') or (answer == 'N') or (answer == '')):
                                print "Nothing to do... exiting"
                                sys.exit(0)

          self._new_proj = new_project


          #print "------ Importing Model ------"

          # MadGraph 5 command interface
          self._MG5Cmd = master_interface.MasterCmd()
          self._mgme_dir = self._MG5Cmd._mgme_dir
          self._MG5Cmd.do_import('model %s --modelname' % self._modelname)
          #self._MG5Cmd.do_import('model %s' % self._modelname)
          self._Model = self._MG5Cmd._curr_model

          # Gets a list of all the particles in the user-supplied model
          self._particles = self._Model.get('particles')

          #print "----- Creating Default Param Card ------"

          # Copy over the param_card.dat from v4 model files. If the project already exists
          # do nothing but set the member variables to point to the project files.
          if self._new_proj:

                if self._MG5Cmd._model_v4_path:
                        shutil.copy2(self._mgme_dir+'/Models/'+modelname+'/param_card.dat','Projects/param_card.dat')
                 # If the model is v5, generate the default param card from the python files
                else:
                        write_dir = 'Projects/'
                        model_builder = export_v4.UFO_model_to_mg4(self._Model, write_dir)
                        model_builder.create_param_card()

                # Initialize the paramcard variable
                self._paramcard = 'Projects/param_card.dat'
                self._project_exists = False

                if ((self._modelname == 'mssm') or (self._modelname[0:5] == 'mssm-')):
                        print "Warning: You are using a restricted model. MadDM will automatically change the parameter card format."
                        check_param_card.convert_to_mg5card(self._paramcard, self._paramcard)

          else:
                project_list = glob.glob('Projects/'+self._projectname+'_*')
                self._projectpath = project_list[0]
                self._paramcard = self._projectpath+'/Cards/param_card.dat'

          #print "------ Initialization Complete! ------\n"

        #-----------------------------------------------------------------------#
        def FindDMCandidate(self, prompts = True, dm_candidate='', exclude_particles=''):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine finds the dark matter candidate and assigns it to the        #
        #  self._dm_particles list.         There are two ways in which this is done.        #
        #  Either the user can input the DM candidate or let the following                #
        #  algorithm find the dm candidates.                                                              #
        #  exclude_particles variable signals which mediators/bsm particles     #
        #  should be completely excluded from the calculation.
        #                                                                        #
        #  1. The DM particle must be a BSM particle (pdg_code > 25)                     #
        #  2. The particle should have no charge (electric or color)                   #
        #  3. The particle's width should be 0 or 'ZERO'                        #
        #  4. The assigned DM candidate is the lightest of all the particles        #
        #          that meet the above criteria.                                        #
        #                                                                        #
        #-----------------------------------------------------------------------#

          #If there are particles to exclude, set that string up first
          self._excluded_particles = exclude_particles.split()

          # initialize boolean flags for the multiple ways of inputting the DM particles
          found_DM_particle = False
          self._ask_param_card = False

          # All the particles with pdg_codes > 25 are added to the array of BSM particles
          for i in range(len(self._particles)):
                particle_code = self._particles[i].get('pdg_code')
                if (particle_code > 25):
                  #This stores the LOCATION of the particle in the self._particles list.
                  self._bsm_particles.append(i)


          #Exclude the particles specified by the user
          #for excl_part in self._excluded_particles:
          #    for location in self._bsm_particles:
          #    	if self._particles[location].get('name') == excl_part:
          #    		self._bsm_particles.remove(location)


          # When manually entering in the DM and coannihilation particles, both particle and anti-particle
          # are not both required to be listed. If one is not included then it is added later.
          while (not found_DM_particle):

                if (prompts == True):

                  print "Enter DM candidate (press Enter to automatically find the DM candidate): "
                  print " (Enter 'particles' to see a list of particles)"
                  print " *** Exclude particles via the standard MG \"/\" notation"
                  dm_answer = raw_input()
                else:
                  dm_answer = ''

                # Here is the algorithm for automatically finding the DM particle and coannihilation candidates
                if (dm_answer == '' and dm_candidate==''):

                  # Ask if the param card needs to be changed
                  if (not self._ask_param_card):
                        self.ChangeParamCard(prompts)
                        self._ask_param_card = True

                  # Looping over all the BSM particles we check the criteria described above
                  for i in range(len(self._bsm_particles)):
                        if (self._particles[self._bsm_particles[i]]['charge'] == 0.0):

                          # We separate the particles that have 'ZERO' for the width since we don't need to numerically evaluate it
                          if (self._particles[self._bsm_particles[i]]['width'] == 'ZERO'):

                                # If nothing has been found so far then set the first DM candidate
                                if (len(self._dm_particles) == 0):
                                  self._dm_particles.append(self._particles[self._bsm_particles[i]])
                                  self._dm_mass = abs(self.GetMass(self._particles[self._bsm_particles[i]]['pdg_code']))
                                # If we already found a candidate, comare the masses and keep the one that's lighter
                                elif (self._bsm_masses[i] < self._dm_mass):
                                  self._dm_particles[0] = self._particles[self._bsm_particles[i]]
                                  self._dm_mass = abs(self.GetMass(self._particles[self._bsm_particles[i]]['pdg_code']))

                          # If the width is not set to 'ZERO' then we have to get the width from the param card
                          elif (self.GetWidth(self._particles[self._bsm_particles[i]]['pdg_code']) == 0.0):
                                if (len(self._dm_particles) == 0):
                                  self._dm_particles.append(self._particles[self._bsm_particles[i]])
                                  self._dm_mass = abs(self.GetMass(self._particles[self._bsm_particles[i]]['pdg_code']))
                                elif (self._bsm_masses[i] < self._dm_mass):
                                  self._dm_particles[0] = self._particles[self._bsm_particles[i]]
                                  self._dm_mass = abs(self.GetMass(self._particles[self._bsm_particles[i]]['pdg_code']))

                  # Check to see if we actually found a DM candidate
                  if (self._dm_particles == []):
                        print "ERROR: No dark matter candidates in the model!"
                        sys.exit(1)
                  else:
                        found_DM_particle = True

                elif (dm_answer == '' and dm_candidate!=''):
                  # We loop over all the particles and check to see if the desired DM candidate is indeed in the model
                  for particle in self._particles:
                        if ( (particle['name'] == dm_candidate) or (particle['antiname'] == dm_candidate)):
                          self._dm_particles.append(particle)
                          found_DM_particle = True
                          #self._dm_mass = abs(self.GetMass(self._dm_particles[0]['pdg_code']))
                  if found_DM_particle == False:
                        print "ERROR: Dark Matter candidate "+dm_candidate+" does not exist in this model!"
                        sys.exit(0)
                  else:
                  	       # Now that we have the param card set, we can set up the model reader
          					# This is necessary to get the numerical value for the particle masses
          					self._fullmodel = model_reader.ModelReader(self._Model)
          					self._fullmodel.set_parameters_and_couplings(param_card=self._paramcard)

                elif (dm_answer != ''):

                  if (dm_answer == 'particles'):
                          print ''
                          self._MG5Cmd.do_display('particles')
                          print ''
                  else:
                          #split the string
                          dm_answer_split = dm_answer.split('/')
                          dm_answer = dm_answer_split[0].rstrip()
                          #if there are particles to exclude, exclude them.
                          if (len(dm_answer_split) > 1):
                                        self._excluded_particles = dm_answer_split[1].lstrip().split()
                          # We loop over all the particles and check to see if the desired DM candidate is indeed in the model
                          for particle in self._particles:
                                if ((particle['name'] == dm_answer) or (particle['antiname'] == dm_answer)):
                                  self._dm_particles.append(particle)
                                  found_DM_particle = True


                          # Check to see if we found the desired DM candidate in the model.
                          if (self._dm_particles == []):
                                print "WARNING: Dark Matter candidate not present in the model! Try again."
                          else:
                           		                  	       # Now that we have the param card set, we can set up the model reader
          					# This is necessary to get the numerical value for the particle masses
          					self._fullmodel = model_reader.ModelReader(self._Model)
          					self._fullmodel.set_parameters_and_couplings(param_card=self._paramcard)
          					self._dm_mass = abs(self.GetMass(self._dm_particles[0]['pdg_code']))

                else:
                   print "WARNING: Unknown command! Try again"

          # Print out the DM candidate
          if prompts == True:
                print "-------------------------------"
                print "DARK MATTER CANDIDATE:"
                print "-------------------------------"
                print self._dm_particles[0]



        #-----------------------------------------------------------------------------#
        def FindCoannParticles(self, coann_partners='', prompts = True, coann_eps = 0.1):
        #-----------------------------------------------------------------------------#
        #                                                                                  #
        #  This routine finds the coannihilation particles for the relic                  #
        #  density calculation.         Either the user can manually input the desired #
        #  particles, or the code can search for all the BSM particles that are       #
        #  within an input mass difference with the DM candidate.  All                         #
        #  coannihilation particles are then added to the self._dm_particles                  #
        #  list.                                                                          #
        #                                                                                  #
        #-----------------------------------------------------------------------------#

          # initialize boolean flags for the multiple ways of inputting the Coannihilation particles
          found_DM_particle = False
          legit_coann_answer = False

          # Time to find the coannihilation particles
          while (not legit_coann_answer):
				  if coann_partners=='':
					coann_answer=''
					if prompts == True:
							print "Enter the coannihilation particles:"
							print "(press Enter to automatically find the coannihilation particles)"
							print "(Enter 'particles' to see a list of particles)"
							coann_answer = raw_input()


					# Automatically find the coannihilation candidates
					if (coann_answer == ''):
					  legit_coann_answer = True

					  ratio_answer = ''
					  if prompts == True:
							 print "Enter the mass difference ratio desired for coannihilating particles [0.1]:"
							 print "(Enter 0 for no coannihilating particles)"
							 ratio_answer = raw_input()

					  if ratio_answer == '':
							self._coann_eps = coann_eps
					  else:
							self._coann_eps = float(ratio_answer)

					  print "INFO: Using coann_eps="+str(self._coann_eps)
					  # If the param card hasn't been changed then as if the param card needs to be changed.
					  if (not self._ask_param_card):
							self.ChangeParamCard(prompts)
							self._ask_param_card = True
							self._dm_mass = abs(self.GetMass(self._dm_particles[0]['pdg_code']))

					  # If the user wishes to find coannihilation candidates we simply loop over the rest of the BSM particles
					  # and see which particles have a mass within the input fractional mass difference.
					  if (self._coann_eps > 0.0):

							# Loop over BSM particles
							for i in range(len(self._bsm_particles)):
							  if (self._particles[self._bsm_particles[i]] != self._dm_particles[0]):
									if (abs(self._dm_mass-self._bsm_masses[i])/self._dm_mass <= self._coann_eps):
									  self._coann_particles.append(self._particles[self._bsm_particles[i]])

									# If there are BSM particles that are too small to be included in the coannihilation they
									# are still tabulated to include in the final state particles with the SM particles.
									elif ((self._bsm_masses[i] < (1.0-self._coann_eps)*self._dm_mass) and \
										not ('all' in self._excluded_particles) and \
										not (self._particles[self._bsm_particles[i]] in self._excluded_particles)):
									  		self._bsm_final_states.append(self._particles[self._bsm_particles[i]])


					# This is the case where the user inputs their own set of coannihilation particles
					elif (coann_answer != 'particles'):
					  legit_coann_answer = True
					else:
				  		print ''
				  		self._MG5Cmd.do_display('particles')
				  		print ''
				  else:
				  	   coann_answer = coann_partners
				  	   legit_coann_answer = True

				  # break up the string into a list of particles
				  input_coann_names_list = coann_answer.split(" ")

				  # Loop over the model particles so we can add them to the self._dm_particles list
				  for i in range(len(input_coann_names_list)):
						for particle in self._particles:
						  # Checks to see if either the particle or anti-particle is in the model as well as if the
						  # particle isn't already included in the list of DM particles (to avoid double counting)
						  if ((particle['name'] == input_coann_names_list[i]) and (not (particle in self._dm_particles)) \
								  and (not (particle in self._coann_particles))):
								self._coann_particles.append(particle)
						  elif ((particle['antiname'] == input_coann_names_list[i]) and (not (particle in self._dm_particles)) \
								  and (not (particle in self._coann_particles))):
								self._coann_particles.append(particle)

				  # If the param card hasn't been changed then as if the param card needs to be changed.
				  if (not self._ask_param_card):
						self.ChangeParamCard(prompts)
						self._ask_param_card = True
						if (self._dm_mass < 0.0):
						  self._dm_mass = abs(self.GetMass(self._dm_particles[0]['pdg_code']))

				  # For organizational purposes we put the coannihilation particles by alphabetical order by name
				  coann_ordered = sorted(self._coann_particles, key=lambda k: k['name'])
				  self._dm_particles += coann_ordered

				  # If we found any coannihilation particles then print out the candidates
				  if (len(self._coann_particles) > 0 and prompts == True):
						print "-------------------------------"
						print "COANNIHILATION PARTNERS:"
						print "-------------------------------"
				  		for i in range(len(self._coann_particles)):
				  				print self._coann_particles[i]



        #-----------------------------------------------------------------------#
        def GetProjectName(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  Once we find the DM and Coannihilation particles we can append the        #
        #  the project name and see if the project has already been generated        #
        #                                                                        #
        #-----------------------------------------------------------------------#

          project_suffix = ''

          # For all the DM particles we create lists that contain all the names and append to the project suffix
          for dm_particle in self._dm_particles:
                self._dm_names_list.append(dm_particle['name'])
                self._dm_antinames_list.append(dm_particle['antiname'])
                project_suffix += dm_particle['name']

          # Convert the '+', '-', and '~' to 'p', 'm', and 'x' respectively
          project_suffix = self.Convertname(project_suffix)

          # Now we can set the project name and project path
          self._projectname += '_'+project_suffix
          self._projectpath = os.path.join('Projects', self._projectname)

          if os.path.isdir(self._projectpath):
                self._project_exists = True



        #-----------------------------------------------------------------------#
        def ChangeParamCard(self, prompts = True):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine allows the user to either edit the default param card        #
        #  or enter in the location of the param card that they wish to use.        #
        #  if the param card entered doesn't exist we'll just repeat until a        #
        #  valid paramcard or answer is given.                                        #
        #                                                                        #
        #-----------------------------------------------------------------------#

          # Loop until we have a vaid result to change (or not change) the param card
          change_param_card = False
          while (not change_param_card):

                # Ask if the param_card needs to be changed
                answer = ''
                if prompts == True:
                        print "Would you like to edit the default param_card?[n] (y/n):"
                        print "(or enter the location of param_card to be used)"
                        answer = raw_input()

                # Edit the default param card
                if ((answer == 'y') or (answer == 'Y')):
                  subprocess.call('vi %s' % self._paramcard, shell= True)
                  change_param_card = True
                # Do nothing
                elif ((answer == 'n') or (answer == 'N') or (answer == '')):
                  change_param_card = True
                # Check to see if the entered file exists.        If so, copy it to the param card
                else:
                  if (os.path.isfile(answer)):
                        shutil.copy2(answer,self._paramcard)
                        # If the model is based off of the mssm, then we may need to convert the param card
                        # from the slha1 format to the slha2 format
                        if ((self._modelname == 'mssm') or (self._modelname[0:5] == 'mssm-')):
                          check_param_card.convert_to_mg5card(self._paramcard, self._paramcard)
                        change_param_card = True
                  else:
                        print "\nThe file entered does not exist. Try Again.\n"

          # Now that we have the param card set, we can set up the model reader
          # This is necessary to get the numerical value for the particle masses
          self._fullmodel = model_reader.ModelReader(self._Model)
          self._fullmodel.set_parameters_and_couplings(param_card=self._paramcard)

          # For all the bsm particles we create the array of each particle's mass.
          for i in range(len(self._bsm_particles)):
                self._bsm_masses.append(abs(self.GetMass(self._particles[i]['pdg_code'])))



        #-----------------------------------------------------------------------#
        def GetMass(self, pdg_id):
        #-----------------------------------------------------------------------#
        #                                                                                                                                                #
        #  Finds the mass of a particle from a Model object given a PDG code        #
        #  and returns it. CAUTION: Masses in MadGraph are stored as complex        #
        #  numbers.
        #                                                                                                                                                #
        #-----------------------------------------------------------------------#

          mass = self._fullmodel.get('parameter_dict')[self._Model.get_particle(pdg_id).get('mass')].real

          return mass



        #-----------------------------------------------------------------------#
        def GetWidth(self, pdg_id):
        #-----------------------------------------------------------------------#
        #                                                                                                                                                #
        #  Same as the previous method but for the width.                                                #
        #                                                                                                                                                #
        #-----------------------------------------------------------------------#

          width = self._fullmodel.get('parameter_dict')[self._Model.get_particle(pdg_id).get('width')].real

          return width


        #-----------------------------------------------------------------------#
        def Convertname(self, name):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine takes a character sting, and for every '+', '-', and        #
        #  '~' in the name it is replaced with a 'p', 'm', and 'x'                #
        #  respectively.  This is used for both project and process names.        #
        #                                                                        #
        #-----------------------------------------------------------------------#

          # Gets a list of the individual characters
          name_list = list(name)
          # Loops over the characters and makes appropriate changes
          for i in range(len(name_list)):
                if (name_list[i] == '+'):
                  name_list[i] = 'p'
                elif (name_list[i] == '-'):
                  name_list[i] = 'm'
                elif (name_list[i] == '~'):
                  name_list[i] = 'x'
          # Combines everything to the new name
          name = ''.join(name_list)
          return name

        #-----------------------------------------------------------------------#
        def CheckQuarkMasses(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #        Checks that all the quark masses in the model are set to non-zero        #
        #        values. This is important for direct detection only. If there is        #
        #        a quark with a zero mass, it stops the code.                        #
        #                                                                        #
        #-----------------------------------------------------------------------#

                quarks =[1,2,3,4,5,6]
                stop_code = False
                for quark in quarks:
                        if self.GetMass(quark) == 0:
                                print "ERROR: quark with PDG code "+str(quark)+" has a 0 mass!"
                                stop_code = True

                if stop_code == True:
                        print "ERROR: There are quarks with 0 mass in your model. Direct detection module can not proceed! Exiting."
                        sys.exit(1)


        #-----------------------------------------------------------------------#
        def DiagramsDD(self, i_dm, eff_operators_SI, eff_operators_SD, SI_order, \
                                                SD_order, QED_order):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #    Generates direct detection diagrams. i_dm is the index of DM part. #
        #         Whether spin dependent or spin independent diagrams should be                 #
        #         calculated. XX_order parameters determine maximum order of a                 #
        #         coupling. For instance, if you only want effective vertices for         #
        #         spin independent coupling you would set SI_order = 2 and all others#
        #         to zero. If you want the spin independent full lagrangian + eff.        #
        #         then you need to set SI_order=2 and QED_order=2...                 #
        #         WARNING: This function is a proxy to be used inside                 #
        #         GenerateDiagramsDirDetect() function and no place else!        #
        #                                                                        #
        #-----------------------------------------------------------------------#


                # Generate all the DM particles scattering off of quarks (FOR NOW NO INELASTIC DM)
                #self._MG5Cmd.do_define('nucleon = d u s c b t')
                #self._MG5Cmd.do_define('anti_nucleon = d~ u~ s~ c~ b~ t~')

                quarks = ['d', 'u', 's', 'c', 'b','t']
                antiquarks = ['d~', 'u~', 's~', 'c~', 'b~','t~']


                #loop over quarks
                for i in range(0, 6):
                        try:
                                if self._excluded_particles!=[] and 'all' not in self._excluded_particles:
                                        proc = self._dm_names_list[i_dm]+' '+quarks[i]+' > '+self._dm_names_list[i_dm]+' '+quarks[i]+' /'\
                                                        +' '.join(self._excluded_particles)+' '\
                                                        +str(eff_operators_SD)+'='+str(SD_order)+' '+str(eff_operators_SI)+\
                                                        '='+str(SI_order)+' QED='+str(QED_order)
                                else:
                                        proc = self._dm_names_list[i_dm]+' '+quarks[i]+' > '+self._dm_names_list[i_dm]+' '+quarks[i]+' '\
                                                        +str(eff_operators_SD)+'='+str(SD_order)+' '+str(eff_operators_SI)+\
                                                        '='+str(SI_order)+' QED='+str(QED_order)

                                #proc = '~sc u > ~sc u'
                                print "Trying "+proc

                                #Necessary to generate diagrams for both quarks and anti-quarks in order
                                #to extract the even and odd parts of the amplitude.

                                self._MG5Cmd.do_generate(proc)

                                # Get the matrix elements
                                curr_matrix_elements_DD = helas_objects.HelasMultiProcess(self._MG5Cmd._curr_amps)
                                #print 'Number of processes for quarks: '+str(len(curr_matrix_elements_DD))
                                self._wanted_lorentz += curr_matrix_elements_DD.get_used_lorentz()
                                self._wanted_couplings += curr_matrix_elements_DD.get_used_couplings()

                                #Here pick out the appropriate matrix element array based on coupling order
                                if SI_order ==0 and SD_order==0 and QED_order==2:
                                                self._scatteringDirDetect_me[i_dm] += curr_matrix_elements_DD.get_matrix_elements()
                                elif (SI_order==2 or SD_order==2) and QED_order==0:
                                                self._scatteringDirDetect_eff_me[i_dm] += curr_matrix_elements_DD.get_matrix_elements()
                                elif (SI_order==2 or SD_order==2) and QED_order==2:
                                                self._scatteringDirDetect_tot_me[i_dm] += curr_matrix_elements_DD.get_matrix_elements()

                        except Exception, error:
                                logging.debug(error)
                                print "WARNING: No direct detection diagrams for the chosen coupling order or "+quarks[i]+" quark."

                        #Loop over the antiquarks
                for i in range(0, 6):
                        try:
                                if self._excluded_particles!=[] and 'all' not in self._excluded_particles:
                                        proc2 = self._dm_names_list[i_dm]+' '+antiquarks[i]+' > '+self._dm_names_list[i_dm]+' '+antiquarks[i]+' /'\
                                                        +' '.join(self._excluded_particles)+' '\
                                                        +str(eff_operators_SD)+'='+str(SD_order)+' '+str(eff_operators_SI)+\
                                                        '='+str(SI_order)+' QED='+str(QED_order)
                                else:
                                        proc2 = self._dm_names_list[i_dm]+' '+antiquarks[i]+' > '+self._dm_names_list[i_dm]+' '+antiquarks[i]+' '\
                                                        +str(eff_operators_SD)+'='+str(SD_order)+' '+str(eff_operators_SI)+\
                                                        '='+str(SI_order)+' QED='+str(QED_order)
                                #proc = '~sc u > ~sc u'
                                print "Trying "+proc2

                                #Necessary to generate diagrams for both quarks and anti-quarks in order
                                #to extract the even and odd parts of the amplitude.

                                self._MG5Cmd.do_generate(proc2)

                                # Get the matrix elements
                                curr_matrix_elements_DD = helas_objects.HelasMultiProcess(self._MG5Cmd._curr_amps)
                                #print 'Number of processes for quarks: '+str(len(curr_matrix_elements_DD))
                                self._wanted_lorentz += curr_matrix_elements_DD.get_used_lorentz()
                                self._wanted_couplings += curr_matrix_elements_DD.get_used_couplings()

                                #Here pick out the appropriate matrix element array based on coupling order
                                if SI_order ==0 and SD_order==0 and QED_order==2:
                                                self._scatteringDirDetect_me[i_dm] += curr_matrix_elements_DD.get_matrix_elements()
                                elif (SI_order==2 or SD_order==2) and QED_order==0:
                                                self._scatteringDirDetect_eff_me[i_dm] += curr_matrix_elements_DD.get_matrix_elements()
                                elif (SI_order==2 or SD_order==2) and QED_order==2:
                                                self._scatteringDirDetect_tot_me[i_dm] += curr_matrix_elements_DD.get_matrix_elements()

                        except Exception, error:
                                logging.debug(error)
                                print "WARNING: No direct detection diagrams for the chosen coupling order or "+antiquarks[i]+" quark."


        #-----------------------------------------------------------------------#
        def GenerateDiagramsDirDetect(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  User level function which performs direct detection functions        #
        #  Generates the DM - q,g scattering matrix elements for spin dependent #
        #  and spin independent direct detection cross section                        #
        #  Currently works only with canonical mode and one dm candidate.        #
        #  The function also merges the dark matter model with the effective        #
        #  vertex model.                                                         #
        #                                                                        #
        #-----------------------------------------------------------------------#

          #First check that all the quark masses in the model are set to non-zero values
          self.CheckQuarkMasses()

          self._do_direct_detection = True

          # If the project already exists we don't need to generate diagrams
          if (self._project_exists):
                return

          #Add the effective vertex model to the existing DM model
          #First determine whether DM is a complex or a real field
          if self._dm_particles[0]['name'] == self._dm_particles[0]['antiname']:
                  operators_folder = 'REAL'
          else:
                  operators_folder = 'COMPLEX'
          print "\nAdding effective vertices from folder "+operators_folder+" to the model file...\n"
          print 'model '+rp+'/'+pp+'/EffOperators/'+operators_folder
          try:
                  eff_dm_name = self._eff_model_dm_names[int(self._dm_particles[0]['spin'])]
                  print eff_dm_name
                  mg5_command = 'model '+rp+'/'+pp+'/EffOperators/'+\
                                  operators_folder+' '+eff_dm_name+'='+self._dm_particles[0].get('name')+' --recreate'
                  print mg5_command
                  self._MG5Cmd.do_add(mg5_command)
                  print "Current model: ",
                  print self._MG5Cmd._curr_model.get('modelpath')
                  self._Model = self._MG5Cmd._curr_model
                  #self._MG5Cmd.do_add_model(rp+'/'+pp+'/EffOperators/'+operators_folder)
          except Exception, error:
                logging.debug(error)
                print "ERROR: Failed to add the effective vertices!"\
                          +" Try to add the effective vertices to your model in MadGraph as a check."
                sys.exit(1)
          print "INFO: Finished adding the Effective vertex model."

          #Now figure out the label of the effective vertex to use. The convention is:
          #<SI or SD>EFF<F, S, or V>, i.e. SIEFFV for Spin Independent Effective vertex for Vector Current.
          eff_operators_SI = self._eff_operators_SI[int(self._dm_particles[0]['spin'])]
          eff_operators_SD = self._eff_operators_SD[int(self._dm_particles[0]['spin'])]

          print "INFO: Using SI effective operator coupling order "+str(eff_operators_SI)
          print "INFO: Using SD effective operator coupling order "+str(eff_operators_SD)
          # Arrays which hold the matrix elements for direct detection in the following order:
          # just full lagrangian, just effective lagrangian, full + effective lagrangian.
          # All three are needed to extract the interference terms which give us the scattering
          # amplitude coefficients.
          self._scatteringDirDetect_me = [[] for i in range(len(self._dm_particles))]
          self._scatteringDirDetect_eff_me = [[] for i in range(len(self._dm_particles))]
          self._scatteringDirDetect_tot_me = [[] for i in range(len(self._dm_particles))]


          #print 'dms: '+str(self._dm_names_list)

          for i in [0]:#range(len(self._dm_particles)):
          # Generate the direct detection scattering diagrams in three steps
          # 1. Generate M_tot + M_eff, where tot is the matrix element for the full lagrangian
          # 2. Generate M_tot, M_eff separately so that the cross term can be extracted later.
                #============================
                #SPIN INDEPENDENT PART
                #============================
                print "INFO: Doing the spin independent part..."
                #Generate diagrams

                #ONLY FULL LAGRANGIAN
                #(marked by EFFX=0, where X is some suffix specifying the spin of DM.)
                print "INFO: Generating X Nucleon > X Nucleon diagrams from the full lagrangian..."
                self.DiagramsDD(i, eff_operators_SI, eff_operators_SD, 0, 0, 2)

                #If matrix element generation went well, generate also the matrix elements for the effective
                #operators. Otherwise just continue to the rest of the code.
                if (not self._scatteringDirDetect_me[0]):
                        self._do_direct_detection = False
                        self._do_directional_detection = False
                        continue
                else:
                        #If generation of matrix elements for direct detection using the full model went ok,
                        #generate the effective model matrix elements and the effective+full lagrangian.

                        print "INFO: Generating X Nucleon > X Nucleon diagrams from the effective lagrangian..."
                        #ONLY EFFECTIVE LAGRANGIAN
                        self.DiagramsDD(i, eff_operators_SI, eff_operators_SD, 2, 0, 0)

                        print "INFO: Generating X Nucleon > X Nucleon diagrams from the effective+full lagrangian..."
                        #EFFECTIVE + FULL
                        self.DiagramsDD(i, eff_operators_SI, eff_operators_SD, 2, 0, 2)

                        print "INFO: Done!"

                        #============================
                        #SPIN DEPENDENT PART
                        #============================
                        #Generate diagrams
                        print "INFO: Doing the spin dependent part..."

                        if (eff_operators_SD != False):

                                #ONLY EFFECTIVE LAGRANGIAN
                                print "INFO: Generating X Nucleon > X Nucleon diagrams from the effective lagrangian..."
                                self.DiagramsDD(i, eff_operators_SI, eff_operators_SD, 0, 2, 0)

                                #EFFECTIVE + FULL
                                print "INFO: Generating X Nucleon > X Nucleon diagrams from the effective + full lagrangian..."
                                self.DiagramsDD(i, eff_operators_SI, eff_operators_SD, 0, 2, 2)

                                print "INFO: Done!"

        #-----------------------------------------------------------------------#
        def GenerateDiagramsIDLoopIndiced(self, finalstate='a a', excluded_particles=''):
        #-----------------------------------------------------------------------#
        #   Comment here...
        #-----------------------------------------------------------------------#

          #For neutrinos, define the 'v' symbol to be all neutrinos
          self._MG5Cmd.do_define('v= ve vt vm ve~ vt~ vm~')

          if finalstate not in self._allowed_ID_fs:
          	print "Final state ["+finalstate+"] not supported"
          	sys.exit()

          self._do_indirect_detection = True

          # If the project already exists we don't need to generate diagrams
          if (self._project_exists):
                return

          #If there are any excluded particles in the propagator, set them up
          if (excluded_particles!=''):
               self._excluded_particles = excluded_particles

          #Set up the array for indirect detection matrix elements
#          self._annihilationIndirDet_me = [[] for i in range(len(self._dm_particles))]

          # Set up the initial state multiparticles that contain the particle and antiparticle
          self._DM_all_names = list(set(self._dm_names_list + self._dm_antinames_list))
          self._DM_init_state = []
          for i in range(len(self._dm_particles)):
                if (self._dm_names_list[i] != self._dm_antinames_list[i]):
                  self._DM_init_state.append(self._dm_names_list[i]+' '+self._dm_antinames_list[i])
                else:
                  self._DM_init_state.append(self._dm_names_list[i])
                self._MG5Cmd.do_define('DM_particle'+str(i+1)+' = '+self._DM_init_state[i])

		  # Generate the annihilation diagrams by going through all the combinations of
          # initial state particles
          at_counter = 1
          for i in range(len(self._dm_particles)):

                  # Create the appropriate string of initial to final state particles
                  try:
                        if self._excluded_particles!=[]:
                                proc = self._dm_names_list[i]+' '+self._dm_antinames_list[i]+' > '+finalstate+' /'\
                                +' '.join(self._excluded_particles)\
                                +' [virt=ALL] @'+str(at_counter)\
                                #+' @ID_'+finalstate.replace(' ', '')
                        else:
                                proc = self._dm_names_list[i]+' '+self._dm_antinames_list[i]+' > '+finalstate\
                                +' [virt=ALL] @'+str(at_counter)\
                                #+' @ID_'+finalstate.replace(' ', '')
                        print "Trying "+proc
                        self._MG5Cmd.do_generate(proc)
                        at_counter = at_counter+1

                        # Once we generate the diagrams, we then immediately get the matrix elements for this process
                        # so we can keep track of the number of processes as well as generate the next set
#                        curr_matrix_elements = helas_objects.HelasMultiProcess(self._MG5Cmd._curr_amps)
#                        self._annihilationIndirDet_me[i] = curr_matrix_elements.get_matrix_elements()
#                        self._wanted_lorentz += curr_matrix_elements.get_used_lorentz()
#                        self._wanted_couplings += curr_matrix_elements.get_used_couplings()

                        #check if the project folder has been created. If not create it.
                        if  glob.glob(self._projectpath+'_*')==[]:
                            os.makedirs(self._projectpath)

                        self._MG5Cmd.do_output(self._projectpath+'/ID_loop_induced_'+finalstate.replace(' ', '')+' --standalone')

                  except Exception, error:
                        print "ERROR: Something went wrong with the loop calculation!"
                        print "ERROR: "+str(error)
                        logging.debug(error)
                        continue



        #-----------------------------------------------------------------------#
        def GenerateDiagramsRelicDensity(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine generates all the 2 -> 2 annihilation matrix elements.        #
        #  The SM particles, BSM particles in self._bsm_final_states, and the        #
        #  DM particles are used as final states.  The initial states are        #
        #  looped over all the different combinations of DM particles.                #
        #                                                                        #
        #-----------------------------------------------------------------------#
          #print 'dms: '+str(self._dm_names_list)
          self._do_relic_density = True

          # If the project already exists we don't need to generate diagrams
          if (self._project_exists):
                return

          # Tabulates all the BSM particles so they can be included in the
          # final state particles along with the SM particles.
          print "INFO: Excluding the following BSM particles from the final state: "+str(self._excluded_particles)
          if not ('all' in self._excluded_particles):
	  	       for i in self._bsm_particles:
                          if ((not (self._particles[i] in self._dm_particles)) and\
                          (not (self._particles[i].get('width') == 'ZERO') ) and\
                          (not (self._particles[i].get('name') in self._excluded_particles))):
                          #if (self._bsm_masses[i] < self._dm_mass):
                          #_bsm_particles holds the location in the array of particles
                                self._bsm_final_states.append(self._particles[i])


          # Else generate the matrix elements
          #print "----- Generating Matrix Elements -----\n",

          # Initialize the arrays that will store all the matrix element information
          self._annihilation_me = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._dm2dm_me = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._scattering_me = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

          # Set up the initial state multiparticles that contain the particle and antiparticle
          self._DM_all_names = list(set(self._dm_names_list + self._dm_antinames_list))
          self._DM_init_state = []
          for i in range(len(self._dm_particles)):
                if (self._dm_names_list[i] != self._dm_antinames_list[i]):
                  self._DM_init_state.append(self._dm_names_list[i]+' '+self._dm_antinames_list[i])
                else:
                  self._DM_init_state.append(self._dm_names_list[i])
                self._MG5Cmd.do_define('DM_particle'+str(i+1)+' = '+self._DM_init_state[i])

          # Set up the SM pdg codes
          leptons = range(11, 17)
          quarks = range(1, 7)
          bosons = range(21, 26)
          sm_pdgs = leptons + quarks + bosons

          # Get the names of the SM particles from the pdg codes
          sm_names_list = []
          for particle in self._particles:
                if particle['pdg_code'] in sm_pdgs:
                  sm_names_list.append(particle['name'])
                  if (not particle['self_antipart']):
                        sm_names_list.append(particle['antiname'])

          # With the list of BSM particles tabulated in FindDMCandidate we can get the list of bsm names
          bsm_names_list = []
          for bsm_particle in self._bsm_final_states:
                bsm_names_list.append(bsm_particle.get('name'))
                if (not bsm_particle['self_antipart']):
                  bsm_names_list.append(bsm_particle.get('antiname'))

          # Create a MadGraph multiparticle that conatins all three groups of particles
          sm_names = ' '.join(sm_names_list)
          if bsm_names_list!=[]:
          	bsm_names = ' '.join(bsm_names_list)
          else:
            bsm_names=''
          print "INFO: DM is allowed to annihilate into the following BSM particles"
          print "INFO: (If you want to exclude them, do so using the standard '/' MadGraph notation):"
          print str(bsm_names)
          dm_names = ' '.join(self._DM_all_names)
          self._MG5Cmd.do_define('dm_particles = '+dm_names)
          self._MG5Cmd.do_define('fs_particles = '+sm_names+' '+bsm_names)


          # Generate the annihilation diagrams by going through all the combinations of
          # initial state particles (including coannihilations)
          for i in range(len(self._dm_particles)):
                for j in range(i,len(self._dm_particles)):

                  # Create the appropriate string of initial to final state particles
                  if True:
                        print self._excluded_particles
                  #try:
                        print "Dm_particle: "+str(self._DM_init_state[i])
                        print "dm_names_list: "+str(self._dm_names_list[i])
                        if self._excluded_particles!=[] and 'all' not in self._excluded_particles:
                                proc = 'DM_particle'+str(i+1)+' DM_particle'+str(j+1)+' > fs_particles fs_particles /'\
                                +' '.join(self._excluded_particles)\
                                +' QED=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0'
                        else:
                                proc = 'DM_particle'+str(i+1)+' DM_particle'+str(j+1)+' > fs_particles fs_particles'\
                                +' QED=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0'
                        print "Trying "+proc
                        try:
                           print proc
                           self._MG5Cmd.do_generate(proc)
                        except Exception, error:
                           continue

                        # Once we generate the diagrams, we then immediately get the matrix elements for this process
                        # so we can keep track of the number of processes as well as generate the next set
                        curr_matrix_elements = helas_objects.HelasMultiProcess(self._MG5Cmd._curr_amps)
                        self._annihilation_me[i][j] = curr_matrix_elements.get_matrix_elements()
                        self._wanted_lorentz += curr_matrix_elements.get_used_lorentz()
                        self._wanted_couplings += curr_matrix_elements.get_used_couplings()

						#---------------------------------------------------------------

                        #Find the resonances
                        print "Finding the locations of resonances..."

                        for helas_matrix_element in self._annihilation_me[i][j]:
                        		helas_model = helas_matrix_element.get('processes')[0].get('model')
                        		amp = helas_matrix_element.get_base_amplitude()
                        		diagrams = amp.get('diagrams')

								#Plot the annihilation diagrams
                        		plot = draw.MultiEpsDiagramDrawer(diagrams,
                                          'maddm_plot.eps',
                                          model=helas_model,
                                          amplitude=amp,
                                          legend=amp.get('process').input_string(),
                                          diagram_type='')
                            		plot.draw()

					for diag in helas_matrix_element.get('diagrams'):
										for amp in diag.get('amplitudes'):
											helas_n_initial_states = helas_matrix_element.get('processes')[0].get_ninitial()
											#print "Model name: ",
											#print helas_model.get('name')
											#print helas_model.get('particle_dict').keys()
											s_channels, t_channels = amp.\
												get_s_and_t_channels(helas_n_initial_states, helas_model,0)
											for s_channel in s_channels:
												resonance_pdg = s_channel.get('legs')[-1].get('id')
												#If the resonance pdg code is 0, that's a fake resonance
												#used to emulate a 4 point interaction
												if resonance_pdg ==0:
													continue
												#print "Resonance PDG: ",
												#print resonance_pdg
												#print "s-channel legs: ",
												#print s_channel.get('legs')
												resonance_mass = helas_model.get_particle(resonance_pdg).get('mass')
												resonance_width = helas_model.get_particle(resonance_pdg).get('width')
												self._resonances.append((resonance_pdg, resonance_mass, resonance_width))

					#Select only the unique occurences of resonances
					self._resonances=list(set(self._resonances))
                          #---------------------------------------------------------------

                 # except Exception, error:
                 #       logging.debug(error)

          if len(self._annihilation_me[i][j])==0:
                    print "WARNING: No DM annihilation diagrams. Relic density calculation can not proceed!"
                    self._do_relic_density = False
                    return

          # Generate all the DM -> DM processes
          for i in range(len(self._dm_particles)):
            for j in range(i,len(self._dm_particles)):

                  # Create the appropriate string of initial to final state particles
                  try:
                        if self._excluded_particles!=[]:
                                proc = 'DM_particle'+str(i+1)+' DM_particle'+str(j+1)+' > dm_particles dm_particles /'\
                                +' '.join(self._excluded_particles)+' '\
                                +'QED=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0'
                        else:
                                proc = 'DM_particle'+str(i+1)+' DM_particle'+str(j+1)+' > dm_particles dm_particles'\
                                +' QED=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0'
                        self._MG5Cmd.do_generate(proc)

                        # Get the matrix elements and make sure that we don't have any pure scattering processes
                        curr_matrix_elements = helas_objects.HelasMultiProcess(self._MG5Cmd._curr_amps)
                        self._dm2dm_me[i][j] = curr_matrix_elements.get_matrix_elements()
                        self._wanted_lorentz += curr_matrix_elements.get_used_lorentz()
                        self._wanted_couplings += curr_matrix_elements.get_used_couplings()
                        for me in self._dm2dm_me[i][j]:
                          if (set(me.get('processes')[0].get_initial_ids()) == (set(me.get('processes')[0].get_final_ids()))):
                                self._dm2dm_me[i][j].remove(me)
                  except Exception, error:
                        logging.debug(error)
                        continue

          # Generate all the DM particles scattering off of the thermal background and
          # change to a different DM species (again, no pure scatterings)
          for i in range(len(self._dm_particles)):
                for j in range(len(self._dm_particles)):

                  # We do not want to genereate the purely scattering processes
                  if (i != j):

                        # Create the appropriate strong of initial to final state particles
                        try:
                          if self._excluded_particles!=[]:
                                  proc = 'DM_particle'+str(i+1)+' fs_particles > DM_particle'+str(j+1)+' fs_particles /'
                                  +' '.join(self._excluded_particles)+' '\
                                  +' QED=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0'
                          else:
                                  proc = 'DM_particle'+str(i+1)+' fs_particles > DM_particle'+str(j+1)+' fs_particles'\
                                  +' QED=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0'
                          self._MG5Cmd.do_generate(proc)

                          # Get the matrix elements
                          curr_matrix_elements = helas_objects.HelasMultiProcess(self._MG5Cmd._curr_amps)
                          self._scattering_me[i][j] = curr_matrix_elements.get_matrix_elements()
                          self._wanted_lorentz += curr_matrix_elements.get_used_lorentz()
                          self._wanted_couplings += curr_matrix_elements.get_used_couplings()

                        except Exception, error:
                          logging.debug(error)
                          continue

          # Look at which particles have the thermal scattering diagrams and those that do
          # get flagged so that the fortran side will know how to properly handle the numerical code
          for i in range(len(self._dm_particles)):
                sum_thermal_scattering = 0
                for j in range(len(self._dm_particles)):
                  sum_thermal_scattering += len(self._scattering_me[i][j])
                if (sum_thermal_scattering == 0):
                  self._dm_thermal_scattering.append(False)
                else:
                  self._dm_thermal_scattering.append(True)

          #print "Finished!\n"



        #-----------------------------------------------------------------------#
        def CreateNumericalSession(self, prompts = True):
        #-----------------------------------------------------------------------#
        #                                                                            #
        #  This routine creates several files that are needed to run the        #
        #  FORTRAN side of the code.  The files that are created include:        #
        #                                                                        #
        #  - all the individual matrix element fortran files.                        #
        #  - all the individual pmass include files associated with each        #
        #         matrix element                                                        #
        #  - makefiles                                                                #
        #                                                                          #
        #                                                                        #
        #-----------------------------------------------------------------------#

          # Since we are done generating the matrix elements we convert the names so the
          # particle names can be used in the fortran code.
          for i in range(len(self._dm_particles)):
                self._dm_names_list[i] = self.Convertname(self._dm_names_list[i])
                self._dm_antinames_list[i] = self.Convertname(self._dm_antinames_list[i])
          for i in range(len(self._DM_all_names)):
                self._DM_all_names[i] = self.Convertname(self._DM_all_names[i])

          # If the project already exists just move over the param card and pickle the dm object
          if (self._project_exists):
                # Move over the param_card to the project directory
                #shutil.move('Projects/param_card.dat',self._projectpath+'/Cards/param_card.dat')
                self._paramcard = os.path.join(self._projectpath,'Cards/param_card.dat') #MIHAILO ADDED THIS

                #Dump the darkmatter object so it can be read back in by the param. scan script
                pickle.dump(self, open( self._projectpath+'/dm_object.pik', 'wb' ), -1)

                return

          # Otherwise create a new numerical session

          # Set up the export class
          self._exporter = MGoutput.ProcessExporterFortranMadDM(self._mgme_dir, self._projectpath)
          self._exporter.opt['model'] = self._modelname

          # Copy over the template directory
          self._exporter.copy_template(self._projectpath)

          #Copy the parameter scan default script
          shutil.copy('param_scan_default.py',self._projectpath+'/param_scan_default.py') #MIHAILO ADDED THIS
          #Dump the darkmatter object so it can be read back in by other Python scripts.
          print "PP: "+self._projectpath
          self._paramcard = os.path.join(self._projectpath,'Cards/param_card.dat') #MIHAILO ADDED THIS

          #Finally, pickle the dm object
          pickle.dump(self, open( self._projectpath+'/dm_object.pik', 'wb' ), -1)  #MIHAILO ADDED THIS

          # Arrays that are needed to store the information for the annihilation diagrams
          self._ann_nprocesses = [[0 for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._ann_process_names = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._ann_process_iden_init = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

          # Arrays that are needed to store the information for the DM -> DM diagrams
          self._dm2dm_nprocesses = [[0 for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._dm2dm_process_names = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._dm2dm_process_iden_init = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._dm2dm_final_states = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

          # Arrays that are needed to store the information for the DM SM -> DM SM processes
          self._scattering_nprocesses = [[0 for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._scattering_process_names = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._scattering_initial_dofs = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._scattering_initial_dofs_total = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

          # Arrays needed for Direct Detection
          self._dd_initial_dofs = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._dd_initial_dofs_total = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._dd_process_names = []
          self._dd_process_ids = []
          self._dd_eff_process_names = []
          self._dd_eff_process_ids = []
          self._dd_eff_initial_dofs = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._dd_eff_initial_dofs_total = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._dd_tot_process_names = []
          self._dd_tot_process_ids = []
          self._dd_tot_initial_dofs = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
          self._dd_tot_initial_dofs_total = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

          #Arrays needed for Indirect Detection
#          self._id_ann_nprocesses = [[0 for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
#          self._id_ann_process_names = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]
 #         self._id_ann_process_iden_init = [[[] for i in range(len(self._dm_particles))] for j in range(len(self._dm_particles))]

          path_matrix = os.path.join(self._projectpath, 'matrix_elements')

          #FOR RELIC DENSITY
          #=============================================
          if (self._do_relic_density == True):
          # writing the matrix elements
                  # Annihilation matrix elements
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):

                          # Get the total number of annihilation processes
                          self._ann_nprocesses[i][j] = len(self._annihilation_me[i][j])

                          for me in self._annihilation_me[i][j]:

                                # Get the name of the process (we disregard the first two characters in the name '0_')
                                process_name = me.get('processes')[0].shell_string()
                                process_name = process_name[2:len(process_name)]
                                self._ann_process_names[i][j].append(process_name)

                                # Check to see if the initial state particles are identical
                                initial_state = me.get('processes')[0].get_initial_ids()
                                if (initial_state[0] == initial_state[1]):
                                  self._ann_process_iden_init[i][j].append(True)
                                else:
                                  self._ann_process_iden_init[i][j].append(False)

                                # Using the proess name we create the filename for each process and export the fortran file
                                filename_matrix = os.path.join(path_matrix, 'matrix_' + process_name + '.f')
                                self._exporter.write_matrix_element(\
                                                writers.FortranWriter(filename_matrix),\
                                                me, self._MG5Cmd._curr_fortran_model)

                                # We also need the external mass include file for each process.
                                filename_pmass = os.path.join(self._projectpath, 'include', 'pmass_' + process_name + '.inc')
                                self._exporter.write_pmass_file(\
                                                writers.FortranWriter(filename_pmass), me)


                  # DM -> DM matrix elements
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):

                          # Get the total number of dm2dm processes
                          self._dm2dm_nprocesses[i][j] = len(self._dm2dm_me[i][j])

                          for me in self._dm2dm_me[i][j]:

                                # Get the name of the process (we disregard the first two characters in the name '0_')
                                process_name = me.get('processes')[0].shell_string()
                                process_name = process_name[2:len(process_name)]
                                self._dm2dm_process_names[i][j].append(process_name)

                                # Check to see if the initial state particles are identical
                                initial_state = me.get('processes')[0].get_initial_ids()
                                if (initial_state[0] == initial_state[1]):
                                  self._dm2dm_process_iden_init[i][j].append(True)
                                else:
                                  self._dm2dm_process_iden_init[i][j].append(False)

                                # Get the final states of the DM -> DM processes
                                final_state = me.get('processes')[0].get_final_ids()
                                # Change the PDG codes to the index in the DM particle list
                                for k in range(len(self._dm_particles)):
                                  if abs(final_state[0]) == abs(self._dm_particles[k]['pdg_code']):
                                        final_state[0] = k+1
                                  if abs(final_state[1]) == abs(self._dm_particles[k]['pdg_code']):
                                        final_state[1] = k+1
                                self._dm2dm_final_states[i][j].append(final_state)

                                # Using the proess name we create the filename for each process and export the fortran file
                                filename_matrix = os.path.join(path_matrix, 'matrix_' + process_name + '.f')
                                self._exporter.write_matrix_element(\
                                                writers.FortranWriter(filename_matrix),\
                                                me, self._MG5Cmd._curr_fortran_model)

                                # We also need the external mass include file for each process.
                                filename_pmass = os.path.join(self._projectpath, 'include', 'pmass_' + process_name + '.inc')
                                self._exporter.write_pmass_file(\
                                                writers.FortranWriter(filename_pmass), me)


                  # SM DM -> SM DM matrix elements
                  for i in range(len(self._dm_particles)):
                        for j in range(len(self._dm_particles)):

                          # Add to the total number of dm2dm processes
                          self._scattering_nprocesses[i][j] = len(self._scattering_me[i][j])

                          for me in self._scattering_me[i][j]:

                                # Get the name of the process (we disregard the first two characters in the name '0_')
                                process_name = me.get('processes')[0].shell_string()
                                process_name = process_name[2:len(process_name)]
                                self._scattering_process_names[i][j].append(process_name)

                                # Get the number of degrees of freedom for the SM particle in the initial state
                                initial_state = me.get('processes')[0].get_initial_ids()
                                for k in range(len(self._particles)):
                                  if (abs(initial_state[1]) == self._particles[k]['pdg_code']):
                                        initial_SM_dof = self._particles[k]['spin']*self._particles[k]['color']
                                        self._scattering_initial_dofs[i][j].append(initial_SM_dof)
                                        if (self._particles[k]['self_antipart']):
                                          self._scattering_initial_dofs_total[i][j].append(initial_SM_dof)
                                        else:
                                          self._scattering_initial_dofs_total[i][j].append(2*initial_SM_dof)

                                # Using the proess name we create the filename for each process and export the fortran file
                                filename_matrix = os.path.join(path_matrix, 'matrix_' + process_name + '.f')
                                self._exporter.write_matrix_element(\
                                                writers.FortranWriter(filename_matrix),\
                                                me, self._MG5Cmd._curr_fortran_model)

                                # We also need the external mass include file for each process.
                                filename_pmass = os.path.join(self._projectpath, 'include', 'pmass_' + process_name + '.inc')
                                self._exporter.write_pmass_file(\
                                                writers.FortranWriter(filename_pmass), me)


      #FOR DIRECT DETECTION
      #================================================================================
          self._total_dd_nprocesses = 0
          if (self._do_direct_detection == True):
                self._total_dd_nprocesses = len(self._scatteringDirDetect_me[0])
                print "Number of DD diags: "+str(self._total_dd_nprocesses)
#                        self._scatteringDirDetect_me[i]
#                        Here, we'll have to add the inelastic contributions.
                for j in [0]:#range(len(self._dm_particles)):
                          for me in self._scatteringDirDetect_me[j]:

                                # Get the label of the process (we disregard the first two characters in the name '0_')
                                process_name = me.get('processes')[0].shell_string()
                                #print "\n"
                                #print me.get('processes')[0]
                                #print "\n"
                                process_name = process_name[2:len(process_name)]
                                self._dd_process_names.append(process_name)
                                #print "Process: "+process_name

                                #Also label the processes by the SM particle which DM scatters off.
                                #We use the PDG code of the sm particle. This info is useful in process_names.inc
                                #to identify which smatrix_... call corresponds to which SM particle since we have
                                #to multiply each matrix element by the appropriate form factor.
                                process_ids = me.get('processes')[0].get_initial_ids()
                                self._dd_process_ids.append(process_ids[1])

                                # Get the number of degrees of freedom for the SM particle in the initial state
                                initial_state = me.get('processes')[0].get_initial_ids()
                                for k in range(len(self._particles)):
                                        if (abs(initial_state[1]) == self._particles[k]['pdg_code']):
                                                initial_SM_dof = self._particles[k]['spin']*self._particles[k]['color']
                                                #print "DOF test: "+str(initial_SM_dof)
                                                self._dd_initial_dofs[j][j].append(initial_SM_dof)
                                                if (self._particles[k]['self_antipart']):
                                                           self._dd_initial_dofs_total[j][j].append(initial_SM_dof)
                                                else:
                                                           self._dd_initial_dofs_total[j][j].append(2*initial_SM_dof)

                                # Using the process name we create the filename for each process and export the fortran file
                                filename_matrix = os.path.join(path_matrix, 'DD_matrix_' + process_name + '.f')
                                self._exporter.write_matrix_element(\
                                        writers.FortranWriter(filename_matrix),\
                                        me, self._MG5Cmd._curr_fortran_model)

                                # We also need the external mass include file for each process.
                                filename_pmass = os.path.join(self._projectpath, 'include', 'DD_pmass_' + process_name + '.inc')
                                self._exporter.write_pmass_file(\
                                                writers.FortranWriter(filename_pmass), me)


        #AND FOR THE EFFECTIVE INTERACTIONS
        #===================================================================================
                self._total_dd_eff_nprocesses = len(self._scatteringDirDetect_eff_me[j])
                #print "Number of DD eff. diags: "+str(self._total_dd_eff_nprocesses)
#                        self._scatteringDirDetect_me[i]
#                        Here, we'll have to add the inelastic contributions.
                for j in [0]:#range(len(self._dm_particles)):
                          ii = 1
                          for me in self._scatteringDirDetect_eff_me[j]:

                                # Get the label of the process (we disregard the first two characters in the name '0_')
                                process_name = me.get('processes')[0].shell_string()
                                #print "\n"
                                #print me.get('processes')[0]
                                #print "\n"
                                process_name = process_name[2:len(process_name)]
                                if ii <= self._total_dd_eff_nprocesses /2:
                                        self._dd_eff_process_names.append('SI_'+process_name)
                                else:
                                        self._dd_eff_process_names.append('SD_'+process_name)
                                #print "Process: "+process_name

                                #Also label the processes by the SM particle which DM scatters off.
                                #We use the PDG code of the sm particle. This info is useful in process_names.inc
                                #to identify which smatrix_... call corresponds to which SM particle since we have
                                #to multiply each matrix element by the appropriate form factor.
                                process_ids = me.get('processes')[0].get_initial_ids()
                                self._dd_eff_process_ids.append(process_ids[1])

                                # Get the number of degrees of freedom for the SM particle in the initial state
                                initial_state = me.get('processes')[0].get_initial_ids()
                                for k in range(len(self._particles)):
                                        if (abs(initial_state[1]) == self._particles[k]['pdg_code']):
                                                initial_SM_dof = self._particles[k]['spin']*self._particles[k]['color']
                                                #print "DOF test: "+str(initial_SM_dof)
                                                self._dd_eff_initial_dofs[j][j].append(initial_SM_dof)
                                                if (self._particles[k]['self_antipart']):
                                                           self._dd_eff_initial_dofs_total[j][j].append(initial_SM_dof)
                                                else:
                                                           self._dd_eff_initial_dofs_total[j][j].append(2*initial_SM_dof)

                                # Using the proess name we create the filename for each process and export the fortran file
                                if ii <= self._total_dd_eff_nprocesses /2:
                                        filename_matrix = os.path.join(path_matrix, 'DD_matrix_eff_SI_' + process_name + '.f')
                                else:
                                        filename_matrix = os.path.join(path_matrix, 'DD_matrix_eff_SD_' + process_name + '.f')

                                self._exporter.write_matrix_element(\
                                        writers.FortranWriter(filename_matrix),\
                                        me, self._MG5Cmd._curr_fortran_model)

                                # We also need the external mass include file for each process.
                                filename_pmass = os.path.join(self._projectpath, 'include', 'DD_pmass_eff_' + process_name + '.inc')
                                self._exporter.write_pmass_file(\
                                                writers.FortranWriter(filename_pmass), me)

                                ii=ii+1

        #AND FOR THE EFFECTIVE + TOTAL INTERACTIONS
        #===================================================================================
                self._total_dd_tot_nprocesses = len(self._scatteringDirDetect_tot_me[j])
                #print "Number of DD eff. diags: "+str(self._total_dd_tot_nprocesses)
#                        self._scatteringDirDetect_me[i]
#                        Here, we'll have to add the inelastic contributions.
                for j in [0]:#range(len(self._dm_particles)):
                          ii = 1
                          for me in self._scatteringDirDetect_tot_me[j]:

                                # Get the label of the process (we disregard the first two characters in the name '0_')
                                process_name = me.get('processes')[0].shell_string()
                                #print "\n"
                                #print me.get('processes')[0]
                                #print "\n"
                                process_name = process_name[2:len(process_name)]
                                if ii <=self._total_dd_eff_nprocesses /2:
                                        self._dd_tot_process_names.append('SI_'+process_name)
                                else:
                                        self._dd_tot_process_names.append('SD_'+process_name)
                                #print "Process: "+process_name

                                #Also label the processes by the SM particle which DM scatters off.
                                #We use the PDG code of the sm particle. This info is useful in process_names.inc
                                #to identify which smatrix_... call corresponds to which SM particle since we have
                                #to multiply each matrix element by the appropriate form factor.
                                process_ids = me.get('processes')[0].get_initial_ids()
                                self._dd_tot_process_ids.append(process_ids[1])

                                # Get the number of degrees of freedom for the SM particle in the initial state
                                initial_state = me.get('processes')[0].get_initial_ids()
                                for k in range(len(self._particles)):
                                        if (abs(initial_state[1]) == self._particles[k]['pdg_code']):
                                                initial_SM_dof = self._particles[k]['spin']*self._particles[k]['color']
                                                #print "DOF test: "+str(initial_SM_dof)
                                                self._dd_tot_initial_dofs[j][j].append(initial_SM_dof)
                                                if (self._particles[k]['self_antipart']):
                                                           self._dd_tot_initial_dofs_total[j][j].append(initial_SM_dof)
                                                else:
                                                           self._dd_tot_initial_dofs_total[j][j].append(2*initial_SM_dof)

                                # Using the proess name we create the filename for each process and export the fortran file
                                if ii <= self._total_dd_eff_nprocesses /2:
                                        filename_matrix = os.path.join(path_matrix, 'DD_matrix_tot_SI_' + process_name + '.f')
                                else:
                                        filename_matrix = os.path.join(path_matrix, 'DD_matrix_tot_SD_' + process_name + '.f')

                                self._exporter.write_matrix_element(\
                                        writers.FortranWriter(filename_matrix),\
                                        me, self._MG5Cmd._curr_fortran_model)

                                # We also need the external mass include file for each process.
                                filename_pmass = os.path.join(self._projectpath, 'include', 'DD_pmass_tot_' + process_name + '.inc')
                                self._exporter.write_pmass_file(\
                                                writers.FortranWriter(filename_pmass), me)

                                ii=ii+1

          #FOR INDIRECT DETECTION
#           if (self._do_indirect_detection == True):
#           # writing the matrix elements
#                   # Annihilation matrix elements
#                   for i in range(len(self._dm_particles)):
#
#                           # Get the total number of annihilation processes
#                           self._ann_nprocesses_indirect_det[i] = len(self._annihilationIndirDet_me[i])
#
#                           for me in self._annihilationIndirDet_me[i]:
#
#                                 # Get the name of the process (we disregard the first two characters in the name '0_')
#                                 process_name = me.get('processes')[0].shell_string()
#                                 process_name = process_name[2:len(process_name)]
#                                 self._id_ann_process_names[i].append(process_name)
#
#                                 # Check to see if the initial state particles are identical
#                                 initial_state = me.get('processes')[0].get_initial_ids()
#                                 if (initial_state[0] == initial_state[1]):
#                                   self._id_ann_process_iden_init[i][j].append(True)
#                                 else:
#                                   self._id_ann_process_iden_init[i][j].append(False)
#
#                                 # Using the proess name we create the filename for each process and export the fortran file
#                                 filename_matrix = os.path.join(path_matrix, 'matrix_' + process_name + '.f')
#                                 self._exporter.write_matrix_element(\
#                                                 writers.FortranWriter(filename_matrix),\
#                                                 me, self._MG5Cmd._curr_fortran_model)
#
#                                 # We also need the external mass include file for each process.
#                                 filename_pmass = os.path.join(self._projectpath, 'include', 'pmass_' + process_name + '.inc')
#                                 self._exporter.write_pmass_file(\
#                                                 writers.FortranWriter(filename_pmass), me)


#===================================================================================

          # Create the dm_info.inc file that contains all the DM particles as well as spin and mass information
          self.WriteDMInfo()

          # Create the model_info.txt file that is used to calculate the number of relativistic degrees of freedom.
          self.WriteModelInfo()

          # Create the diagrams.inc file that contains the number of processes for each pair of initial state particles
          self.WriteDiagramInfo()

          # Create the process_names.inc file that contains the names of all the individual processes
          self.WriteProcessNames()

          # Modify the matrix element files for DD so that they return the matrix element value
          # and not the matrix element squared.
          #if self._do_direct_detection == True:
          self.Modify_dd_matrix_files(whichME = 'full') #For the full matrix elements
          self.Modify_dd_matrix_files(whichME = 'eff')        #For the effective vertices
          self.Modify_dd_matrix_files(whichME = 'total')  #For the effective + full vertices

          # Create the smatrix.f file
          self.Write_smatrix()

          # Create the makefile for compiling all the matrix elements
          self.Write_makefile()

          # Write the include file
          self.WriteMadDMinc()

          #Write the include file containing all the resonance locations
          self.Write_resonance_info()

          # This creates the fortran files for the model
          # v4 model
          if self._MG5Cmd._model_v4_path:
                print 'Copy %s model files to directory %s' % (os.path.basename(self._model_v4_path), self._projectpath)
                self._exporter.export_model_files(self._model_v4_path)
                self._exporter.export_helas(pjoin(self._mgme_dir,'HELAS'))
          # v5 model
          else:
                print 'Export UFO model to MG4 format'
                #print self._Model['lorentz']
                self._exporter.convert_model_to_mg4(self._Model, self._wanted_lorentz, self._wanted_couplings)


          #Copy the param_card.dat
          shutil.move('Projects/param_card.dat',self._projectpath+'/Cards/param_card.dat')
          #Make sure that it the parameter card is merged correctly with the effecive vertices
          self.fix_param_card()

          # Copy over the coupl.inc and input.inc files to the project's include directory.
          shutil.copy2(self._projectpath+'/Source/MODEL/coupl.inc', self._projectpath + '/include/coupl.inc')
          shutil.copy2(self._projectpath+'/Source/MODEL/input.inc', self._projectpath + '/include/input.inc')

          #Output the logical flags which signal which calculations should be performed
          self.Modify_maddm_card()

          #Edit the maddm_card.inc file
          if prompts == True:
                 print "Would you like to edit maddm_card.inc?[n] (y/n):"
                 answer = raw_input()
                 if answer == 'y' or answer == 'Y':
                        subprocess.call('vi '+self._projectpath+'/include/maddm_card.inc', shell=True)

          # Now that everything is created we can compile the whole code
          #print "------Compiling the numerical code------"
          #curr_dir = os.getcwd() # get current directory
          makedir = self._projectpath
          if self._do_relic_density and self._do_direct_detection:
                  subprocess.call(['make'],cwd = makedir, shell=True)
          elif self._do_relic_density and not self._do_direct_detection:
            subprocess.call(['make relic_density'],cwd = makedir, shell=True)
          elif self._do_direct_detection and not self._do_relic_density:
             subprocess.call(['make direct_detection'],cwd = makedir, shell=True)
          else:
             print 'ERROR: You can not create a numerical session without specifying'+\
                             ' which DM calculations you would like to perform (i.e. relic density, direct deteciton...)'
             sys.exit(1)

          print "Finished!\n"

        #-----------------------------------------------------------------------#
        def Modify_maddm_card(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This function passes the logical flags of which calculations to        #
        #  perform to the fortran code.                                                #
        #                                                                        #
        #-----------------------------------------------------------------------#
                   __flag1 = 'do_relic_density'
                   __flag2 = 'do_direct_detection'
                   __flag3 = 'do_directional_detection'

                   shutil.move('Projects/'+self._projectname+'/include/maddm_card.inc',\
                          'Projects/'+self._projectname+'/include/maddm_card_temp.inc')
                   fileFrom = open(os.path.join('Projects',self._projectname,'include','maddm_card_temp.inc'),'r')
                   fileTo = open(os.path.join('Projects',self._projectname,'include','maddm_card.inc'),'w')
                   lines = fileFrom.readlines()

                   for line in lines:
                          if __flag1 in line:
                                line_split = line.split('=')
                                line = line_split[0]+' = .'+ str(self._do_relic_density).lower()+'.\n'
                          elif __flag2 in line:
                                line_split = line.split('=')
                                line = line_split[0]+' = .'+ str(self._do_direct_detection).lower()+'.\n'
                          elif __flag3 in line:
                                line_split = line.split('=')
                                line = line_split[0]+' = .'+ str(self._do_directional_detection).lower()+'.\n'

                          fileTo.write(line)

                   os.remove('Projects/'+self._projectname+'/include/maddm_card_temp.inc')
                   fileFrom.close()
                   fileTo.close()

        #-----------------------------------------------------------------------#
        def WriteDMInfo(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This function writes the information about the inital state particle #
        #  dof to a file for use by the FORTRAN part of the code.                #
        #                                                                        #
        #-----------------------------------------------------------------------#

          filename_dof = os.path.join('Projects',self._projectname,'include','dm_info.inc')
          file_writer = open(filename_dof, 'w')

          # Write out some header comments for the file
          file_writer.write("c------------------------------------------------------------------------------c\n")
          file_writer.write("c This file contains the needed information about all the DM particles.\n")

          # The counter is used to write the correct array index needed for the fortran code
          counter = 1
          for dm_particle in self._dm_particles:

                dof = float(dm_particle['spin'])*float(dm_particle['color'])

                # Write out the mass parameter, degrees of freedom, particle name, and the index for the end of the name
                file_writer.write("c------------------------------------------------------------------------------c\n")
                file_writer.write("                 mdm(" + str(counter) + ") = abs(" + dm_particle['mass'] + ")\n")
                file_writer.write("                 dof_dm(" + str(counter) +") = " + str(dof) + "\n")
                if (self._do_relic_density and self._dm_thermal_scattering[counter-1]):
                  file_writer.write("           dm_sm_scattering(" + str(counter) + ") = .true.\n")
                else:
                  file_writer.write("           dm_sm_scattering(" + str(counter) + ") = .false.\n")
                file_writer.write("                 dm_names(" + str(counter) + ") = \'" + self._dm_names_list[counter-1] + "\'\n")
                file_writer.write("                 dm_index(" + str(counter) + ") = " + str(len(dm_particle['name'])) + "\n")
                file_writer.write("                 dm_antinames(" + str(counter) + ") = \'" + self._dm_antinames_list[counter-1] + "\'\n")
                file_writer.write("                 dm_antiindex(" + str(counter) + ") = " + str(len(dm_particle['antiname'])) + "\n")

                counter = counter + 1

          file_writer.write("c------------------------------------------------------------------------------c\n")
          file_writer.close()



        #-----------------------------------------------------------------------#
        def WriteModelInfo(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine writes the mass, degrees of freedom and boson/fermion        #
        #  information for all the relevant particles in the model so that the        #
        #  number of relativistic degrees of freedom can be calculated.                #
        #                                                                        #
        #  The SM particles are in the template by default so this would just        #
        #  add the relevant BSM particles.                                        #
        #                                                                        #
        #-----------------------------------------------------------------------#

          # Flags that are needed to write the model info file
          __flag1 = '#NUM_PARTICLES'
          __flag2 = '#BSM_INFO'

          # Open up the template and output file
          model_info_file = open('Projects/'+self._projectname+'/include/model_info.txt','w')
          model_info_template = open('Projects/'+self._projectname+'/include/model_info_template.txt','r')
          model_info_lines = model_info_template.readlines()

          # Read all the lines from the template file and insert appropriate parts
          # which are flagged with __flag1 and __flag2
          for line in model_info_lines:

                # write out the number of DM paticles
                if __flag1 in line:
                  model_info_file.write(str(17 + len(self._dm_particles) + len(self._bsm_final_states))+'\n')

                # write out the mass, degrees of freedom and fermion/boson info for the bsm particles
                elif __flag2 in line:
                  bsm_particles = self._dm_particles + self._bsm_final_states
                  for bsm_particle in bsm_particles:

                        # Get all the necessary information
                        mass = self.GetMass(bsm_particle['pdg_code'])
                        dof = float(bsm_particle['spin'])*float(bsm_particle['color'])

                        # If the particle has an anti particle double the number of degrees of freedom
                        if (not bsm_particle['self_antipart']):
                          dof *= 2.0

                        # The spin information is written in 2s+1 format (i.e. 1, 3, etc = boosn; 2, 4, etc = fermion)
                        if (int(bsm_particle['spin']) == (int(bsm_particle['spin']/2)*2)):
                          boson_or_fermion = 1
                        else:
                          boson_or_fermion = 0

                        # write the line to the file
                        new_line = str(mass)+'         '+str(dof)+'        '+str(boson_or_fermion)+'\n'
                        model_info_file.write(new_line)

                # If there is no flag then it's a SM particle which we just write to the output file
                else:
                  model_info_file.write(line)

          model_info_file.close()



        #-----------------------------------------------------------------------#
        def WriteDiagramInfo(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine creates the diagrams.inc file which contains the number #
        #  of annihilation diagrams for each pair of DM particles.                #
        #                                                                        #
        #-----------------------------------------------------------------------#

          # open up the file and first print out the number of dm particles
          diagramsfile = open('Projects/'+self._projectname+'/include/diagrams.inc','w')

          #don't include the bsm particles which are not dark matter
          ndm_particles =0
          for dm_part in self._dm_particles:
          	if (dm_part['charge']==0 and \
          		(dm_part['width']=='ZERO' or self.GetWidth(dm_part['pdg_code']) == 0.0 )):
          		ndm_particles = ndm_particles+1

#  	  stringtowrite = 'c Total number of IS particles participating in the coannihilations \n'+\
#  												  '                 ndmparticles = '+str(ndm_particles)+'\n\n' +\
#  												  'c Processes by class'
#  	  diagramsfile.write(stringtowrite)

	  stringtowrite = 'c Total number of IS particles participating in the coannihilations \n'+\
												  '                 ndmparticles = '+str(len(self._dm_particles))+'\n\n' +\
												  'c Processes by class'
	  diagramsfile.write(stringtowrite)


          # Write out how many annihlation processes for each combinations of initial states
          diagramsfile.write('\nc Number of annihilation processes for each DM pair\n')
          for i in range(len(self._dm_particles)):
                for j in range(i,len(self._dm_particles)):
                  diagramsfile.write('                ann_nprocesses('+str(i+1)+','+str(j+1)+') = '+str(self._ann_nprocesses[i][j])+'\n')

          # Write out how many DM -> DM processes for each combinations of initial states
          diagramsfile.write('\nc Number of DM -> DM processes for each DM pair\n')
          for i in range(len(self._dm_particles)):
                for j in range(i,len(self._dm_particles)):
                  diagramsfile.write('                dm2dm_nprocesses('+str(i+1)+','+str(j+1)+') = '+str(self._dm2dm_nprocesses[i][j])+'\n')

          # Write out how many DM/SM scattering processes for each combinations of initial states
          diagramsfile.write('\nc Number of DM/SM scattering processes for each DM pair\n')
          for i in range(len(self._dm_particles)):
                for j in range(len(self._dm_particles)):
                  stringtowrite = '                 scattering_nprocesses('+str(i+1) +','+str(j+1) +') = '+\
                          str(self._scattering_nprocesses[i][j])+'\n'
                  diagramsfile.write(stringtowrite)

          #For direct detection
          diagramsfile.write('\nc Number of DM/SM Direct Detection scattering processes\n')

          diagramsfile.write('\nc Number of DM/SM Effective Direct Detection scattering processes\n')


          diagramsfile.close()



        #-----------------------------------------------------------------------#
        def WriteProcessNames(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine creates the process_names.inc file which contains the        #
        #  names of all the individual processes.  These names are primairly        #
        #  used in the test subroutines on the fortran side.                        #
        #                                                                        #
        #-----------------------------------------------------------------------#

          # This creates the file process_names.inc which will have a list of all the individual process names.
          process_names_file = open('Projects/'+self._projectname+'/include/process_names.inc','w')

          # Write out the list of process names
          process_names_file.write('c List of the process names in order of dmi, dmj, nprocesses(dmi, dmj)\n')
          process_counter = 0

          # Annihilation process names
          if (self._do_relic_density == True):

                  process_names_file.write('c Annihilation process names\n')
                  self._total_annihilation_nprocesses = 0
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):
                          self._total_annihilation_nprocesses += self._ann_nprocesses[i][j]
                          for k in range(self._ann_nprocesses[i][j]):

                                # Increment the process counters and write the process name
                                process_counter += 1
                                process_names_file.write('                process_names('+str(process_counter)+') = \''+self._ann_process_names[i][j][k]+'\'\n')

                  # DM -> DM process names
                  process_names_file.write('\nc DM -> DM process names\n')
                  self._total_dm2dmscattering_nprocesses = 0
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):
                          self._total_dm2dmscattering_nprocesses += self._dm2dm_nprocesses[i][j]
                          for k in range(self._dm2dm_nprocesses[i][j]):

                                # Increment the process counter and write the process name
                                process_counter += 1
                                process_names_file.write('                process_names('+str(process_counter)+') = \''+self._dm2dm_process_names[i][j][k]+'\'\n')

                  # DM/SM scattering process names
                  process_names_file.write('\nc DM/SM scattering process names\n')
                  self._total_scattering_nprocesses = 0
                  for i in range(len(self._dm_particles)):
                        for j in range(len(self._dm_particles)):
                          self._total_scattering_nprocesses += self._scattering_nprocesses[i][j]
                          for k in range(self._scattering_nprocesses[i][j]):

                                # Increment the process counter and write the process name
                                process_counter += 1
                                process_names_file.write('                process_names('+str(process_counter)+') = \''+\
                                          self._scattering_process_names[i][j][k]+'\'\n')

                  # Write out the total number of process names and number of each group of processes
                  process_names_file.write('\nc Total number of processes for each category\n')
                  process_names_file.write('          num_processes = ' + str(process_counter) + '\n')
                  process_names_file.write('          ann_num_processes = ' + str(self._total_annihilation_nprocesses) + '\n')
                  process_names_file.write('          dm2dm_num_processes = ' + str(self._total_dm2dmscattering_nprocesses) + '\n')
                  process_names_file.write('          scattering_num_processes = ' + str(self._total_scattering_nprocesses) + '\n')

                  # Write out the boolean flags for the processes with identical initial state particles
                  process_names_file.write('\nc Boolean operators for identical particles in the initial state\n')
                  process_names_file.write('c Annihialtion diagrams\n')
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):

                          # Boolean operators for the annihilation processes
                          for k in range(self._ann_nprocesses[i][j]):

                                if (self._ann_process_iden_init[i][j][k]):
                                  process_names_file.write('          ann_process_iden_init('+str(i+1)+','+str(j+1)+','+str(k+1)+') = .true.\n')
                                else:
                                  process_names_file.write('          ann_process_iden_init('+str(i+1)+','+str(j+1)+','+str(k+1)+') = .false.\n')

                  process_names_file.write('\nc DM -> DM diagrams\n')
                  for i in range(len(self._dm_particles)):
                        for j in range(i, len(self._dm_particles)):

                          # Boolean operators for the dm2dm processes
                          for k in range(self._dm2dm_nprocesses[i][j]):

                                if (self._dm2dm_process_iden_init[i][j][k]):
                                  process_names_file.write('          dm2dm_process_iden_init('+str(i+1)+','+str(j+1)+','+str(k+1)+') = .true.\n')
                                else:
                                  process_names_file.write('          dm2dm_process_iden_init('+str(i+1)+','+str(j+1)+','+str(k+1)+') = .false.\n')


                  # Write out final state information for all the DM -> DM processes
                  process_names_file.write('\nc Final state information for all the DM -> DM processes\n')
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):
                          for k in range(self._dm2dm_nprocesses[i][j]):

                                process_names_file.write('                dm2dm_fs('+str(i+1)+','+str(j+1)+','+str(k+1)+',1) = '+\
                                          str(self._dm2dm_final_states[i][j][k][0])+'\n')
                                process_names_file.write('                dm2dm_fs('+str(i+1)+','+str(j+1)+','+str(k+1)+',2) = '+\
                                          str(self._dm2dm_final_states[i][j][k][1])+'\n')

                  # Write out the inital state degrees of freedom for all the DM/SM scattering processes
                  process_names_file.write('\nc Initial state degrees of freedom for all the DM/SM scattering processes\n')
                  for i in range(len(self._dm_particles)):
                        for j in range(len(self._dm_particles)):
                          for k in range(self._scattering_nprocesses[i][j]):

                                process_names_file.write('                dof_SM('+str(i+1)+','+str(j+1)+','+str(k+1)+') = '+\
                                          str(self._scattering_initial_dofs[i][j][k])+'\n')

                  # Write out the total inital state degrees of freedom for all the DM/SM scattering processes
                  process_names_file.write('\nc Total initial state degrees of freedom for all the DM/SM scattering processes\n')
                  for i in range(len(self._dm_particles)):
                        for j in range(len(self._dm_particles)):
                          for k in range(self._scattering_nprocesses[i][j]):

                                process_names_file.write('                dof_SM_total('+str(i+1)+','+str(j+1)+','+str(k+1)+') = '+\
                                          str(self._scattering_initial_dofs_total[i][j][k])+'\n')


          #Write process information for direct detection
          if (self._do_direct_detection == True):
                                process_names_file.write('\nc  Total number of Direct Detection processes\n')
                                process_names_file.write('          dd_num_processes = ' + str(self._total_dd_nprocesses) + '\n')

                                process_names_file.write('\nc  Processes relevant for Direct Detection \n')
                                for i in range(len(self._dd_process_names)):
                                        process_names_file.write('                dd_process_names('+str(i+1)+') = \''+\
                                                  self._dd_process_names[i]+'\'\n')

                                process_names_file.write('\n')
                                for i in range(len(self._dd_process_ids)):
                                        process_names_file.write('                dd_process_ids('+str(i+1)+') = '+\
                                                  str(self._dd_process_ids[i])+'\n')

          #Write process information for EFFECTIVE direct detection
                                process_names_file.write('\nc  Processes relevant for Direct Detection (effective vertices) \n')
                                for i in range(len(self._dd_eff_process_names)):
                                        process_names_file.write('                dd_eff_process_names('+str(i+1)+') = \''+\
                                                  self._dd_eff_process_names[i]+'\'\n')

                                process_names_file.write('\n')
                                for i in range(len(self._dd_eff_process_ids)):
                                        process_names_file.write('                dd_eff_process_ids('+str(i+1)+') = '+\
                                                  str(self._dd_eff_process_ids[i])+'\n')

          #Write process information for EFFECTIVE + TOTAL direct detection
                                process_names_file.write('\nc  Processes relevant for Direct Detection (effective + full vertices) \n')
                                for i in range(len(self._dd_tot_process_names)):
                                        process_names_file.write('                dd_tot_process_names('+str(i+1)+') = \''+\
                                                  self._dd_tot_process_names[i]+'\'\n')

                                process_names_file.write('\n')
                                for i in range(len(self._dd_tot_process_ids)):
                                        process_names_file.write('                dd_tot_process_ids('+str(i+1)+') = '+\
                                                  str(self._dd_tot_process_ids[i])+'\n')


          process_names_file.close()

	#-----------------------------------------------------------------------#
	def Write_resonance_info(self):
	#-----------------------------------------------------------------------#
		#Write the locations of resonances
		resonances_inc_file = open('Projects/'+self._projectname+'/include/resonances.inc','w+')
          	for i in range (len(self._resonances)):
          				new_line = '                resonances(%d)     =   %s \n'% (i+1, self._resonances[i][1])
          				resonances_inc_file.write(new_line)
          	for i in range (len(self._resonances)):
          				new_line = '                resonance_widths(%d)     =   %s \n'% (i+1, self._resonances[i][2])
          				resonances_inc_file.write(new_line)

		resonances_inc_file.close()


        #-----------------------------------------------------------------------#
        def Write_smatrix(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine creates the smatrix.f file based on the smatrix                #
        #  template.  This is used to call the appropriate matrix element as        #
        #  well as the appropriate pmass include file for each process.                #
        #                                                                        #
        #-----------------------------------------------------------------------#

          # Flags that are needed to create the matrix element makefile and smatrix.f files
          __flag1 = '#PMASS_MADDM1'
          __flag2 = '#PMASS_MADDM2'
          __flag3 = '#PMASS_MADDM3'
          __flag4 = '#SMATRIX_MADDM1'
          __flag5 = '#SMATRIX_MADDM2'
          __flag6 = '#SMATRIX_MADDM3'
          __flag7 = '#SMATRIX_MADDM_DD'
          __flag8 = '#SMATRIX_MADDM_EFF_DD'
          __flag9 = '#SMATRIX_MADDM_TOT_DD'

          # Edit smatrix.f file to incorporate all different subprocess
          smatrix_file = open('Projects/'+self._projectname+'/matrix_elements/smatrix.f','w')
          smatrix_template = open('Projects/'+self._projectname+'/matrix_elements/smatrix_template.f','r')
          smatrix_lines = smatrix_template.readlines()

          # Read all the lines from the smatrix_template.f and insert appropriate parts
          # which are flagged with __flag1 and __flag2 in the smatrix_template.f file.
          # Write the output into smatrix.f
          for line in smatrix_lines:
                # Include all the combinations of external pmass.inc files for the annihlation processes
                if __flag1 in line:
                  firstentry = True
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):
                          for k in range(self._ann_nprocesses[i][j]):
                                if (firstentry):
                                  firstentry = False
                                  new_line = '                if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  include \'pmass_'+self._ann_process_names[i][j][k]+'.inc\' \n'
                                  smatrix_file.write(new_line)
                                else:
                                  new_line = '                else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  include \'pmass_'+self._ann_process_names[i][j][k]+'.inc\' \n'
                                  smatrix_file.write(new_line)
                  if (self._do_relic_density):
                    if (self._total_annihilation_nprocesses > 0):
                           smatrix_file.write('              endif')

                # Include all the combinations of external pmass.inc files for the DM -> DM processes
                elif __flag2 in line:
                  firstentry = True
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):
                          for k in range(self._dm2dm_nprocesses[i][j]):
                                if (firstentry):
                                  firstentry = False
                                  new_line = '                if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  include \'pmass_'+self._dm2dm_process_names[i][j][k]+'.inc\' \n'
                                  smatrix_file.write(new_line)
                                else:
                                  new_line = '                else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  include \'pmass_'+self._dm2dm_process_names[i][j][k]+'.inc\' \n'
                                  smatrix_file.write(new_line)
                  if (self._do_relic_density):
                    if (self._total_dm2dmscattering_nprocesses > 0):
                          smatrix_file.write('              endif')

                # Include all the combinations of external pmass.inc files for the DM/SM scattering processes
                elif __flag3 in line:
                  firstentry = True
                  for i in range(len(self._dm_particles)):
                        for j in range(len(self._dm_particles)):
                          for k in range(self._scattering_nprocesses[i][j]):
                                if (firstentry):
                                  firstentry = False
                                  new_line = '                if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  include \'pmass_'+self._scattering_process_names[i][j][k]+'.inc\' \n'
                                  smatrix_file.write(new_line)
                                else:
                                  new_line = '                else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  include \'pmass_'+self._scattering_process_names[i][j][k]+'.inc\' \n'
                                  smatrix_file.write(new_line)
                  if (self._do_relic_density):
                    if (self._total_scattering_nprocesses > 0):
                           smatrix_file.write('              endif')

                # Creates all the appropriate calls to the individual annihilaton smatrix.f's
                elif __flag4 in line:
                  firstentry = True
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):
                          for k in range(self._ann_nprocesses[i][j]):
                                if (firstentry):
                                  firstentry = False
                                  new_line = '                if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  call smatrix_'+\
                                           self._ann_process_names[i][j][k]+'(p_ext,smatrix_ann)\n'
                                  smatrix_file.write(new_line)
                                else:
                                  new_line = '                else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  call smatrix_'+\
                                           self._ann_process_names[i][j][k]+'(p_ext,smatrix_ann)\n'
                                  smatrix_file.write(new_line)
                  if (self._do_relic_density):
                          if (self._total_annihilation_nprocesses > 0):
                                smatrix_file.write('          endif')

                # Creates all the appropriate calls to the individual annihilaton smatrix.f's
                elif __flag5 in line:
                  firstentry = True
                  for i in range(len(self._dm_particles)):
                        for j in range(i,len(self._dm_particles)):
                          for k in range(self._dm2dm_nprocesses[i][j]):
                                if (firstentry):
                                  firstentry = False
                                  new_line = '                if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  call smatrix_'+\
                                          self._dm2dm_process_names[i][j][k]+'(p_ext,smatrix_dm2dm)\n'
                                  smatrix_file.write(new_line)
                                else:
                                  new_line = '                else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  call smatrix_'+\
                                          self._dm2dm_process_names[i][j][k]+'(p_ext,smatrix_dm2dm)\n'
                                  smatrix_file.write(new_line)
                  if (self._do_relic_density):
                    if (self._total_dm2dmscattering_nprocesses > 0):
                                smatrix_file.write('          endif')

                # Creates all the appropriate calls to the individual annihilaton smatrix.f's
                elif __flag6 in line:
                  firstentry = True
                  for i in range(len(self._dm_particles)):
                        for j in range(len(self._dm_particles)):
                          for k in range(self._scattering_nprocesses[i][j]):
                                if (firstentry):
                                  firstentry = False
                                  new_line = '                if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  call smatrix_'+\
                                          self._scattering_process_names[i][j][k]+'(p_ext,smatrix_scattering) \n'
                                  smatrix_file.write(new_line)
                                else:
                                  new_line = '                else if ((i.eq.'+str(i+1)+') .and. (j.eq.'+str(j+1)+') .and. ' +\
                                          '(k.eq.'+str(k+1)+')) then\n                  call smatrix_'+\
                                          self._scattering_process_names[i][j][k]+'(p_ext,smatrix_scattering) \n'
                                  smatrix_file.write(new_line)
                  if (self._do_relic_density):
                          if (self._total_scattering_nprocesses > 0):
                                smatrix_file.write('          endif')

                #FOR DIRECT DETECTION - for now only one DM candidate.
                elif __flag7 in line:
                  print "Adding the SMATRIX calls for the direct detection."
                  if self._do_direct_detection == True:
                          firstentry = True
                          print "Number of diagrams: "+str(len(self._scatteringDirDetect_me[0]))
                          for i in range(len(self._scatteringDirDetect_me[0])):
                                        if(firstentry):
                                                  firstentry = False
                                                  new_line = '                if (k.eq.'+str(i+1)+') then\n                 call smatrix_'+\
                                                  self._dd_process_names[i]+'(p_ext,smatrix_dd) \n'
                                                  smatrix_file.write(new_line)
                                        else:
                                                  new_line = '                else if (k.eq.'+str(i+1)+') then\n                  call smatrix_'+\
                                                  self._dd_process_names[i]+'(p_ext, smatrix_dd) \n'
                                                  smatrix_file.write(new_line)
                          if (len(self._scatteringDirDetect_me[0]) > 0):
                                        smatrix_file.write('          endif')

                #FOR EFFECTIVE DIRECT DETECTION - for now only one DM candidate.
                elif __flag8 in line:
                  print "Adding the SMATRIX calls for the effective vertices."
                  if self._do_direct_detection == True:
                          firstentry = True
                          print "Number of diagrams: "+str(len(self._scatteringDirDetect_eff_me[0]))
                          for i in range(len(self._scatteringDirDetect_eff_me[0])):
                                        if(firstentry):
                                                  firstentry = False
                                                  new_line = '                if (k.eq.'+str(i+1)+') then\n                 call smatrix_eff_'+\
                                                                  self._dd_eff_process_names[i]+'(p_ext,smatrix_dd_eff) \n'
                                                  smatrix_file.write(new_line)
                                        else:
                                                  new_line = '                else if (k.eq.'+str(i+1)+') then\n                  call smatrix_eff_'+\
                                                                  self._dd_eff_process_names[i]+'(p_ext, smatrix_dd_eff) \n'
                                                  smatrix_file.write(new_line)
                          if (len(self._scatteringDirDetect_eff_me[0]) > 0):
                                        smatrix_file.write('          endif')

                #FOR EFFECTIVE + FULL DIRECT DETECTION - for now only one DM candidate.
                elif __flag9 in line:
                  print "Adding the SMATRIX calls for the effective + full lagrangian vertices."
                  if self._do_direct_detection == True:
                          firstentry = True
                          print "Number of diagrams: "+str(len(self._scatteringDirDetect_tot_me[0]))
                          for i in range(len(self._scatteringDirDetect_tot_me[0])):
                                        if(firstentry):
                                                  firstentry = False
                                                  new_line = '                if (k.eq.'+str(i+1)+') then\n                 call smatrix_tot_'+\
                                                                self._dd_tot_process_names[i]+'(p_ext,smatrix_dd_tot) \n'
                                                  smatrix_file.write(new_line)
                                        else:
                                                  new_line = '                else if (k.eq.'+str(i+1)+') then\n                  call smatrix_tot_'+\
                                                                  self._dd_tot_process_names[i]+'(p_ext, smatrix_dd_tot) \n'
                                                  smatrix_file.write(new_line)
                          if (len(self._scatteringDirDetect_tot_me[0]) > 0):
                                        smatrix_file.write('          endif')


                # If none of the flags are present then just write the line from the template to the file
                else:
                        smatrix_file.write(line)

          smatrix_template.close()
          smatrix_file.close()

        #-----------------------------------------------------------------------#
        def Modify_dd_matrix_files(self, whichME = ''):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine edits the matrix.f files relevent for direct detec.        #
        #  Since we have three different smatrix functions                        #
        #  (full, eff+full and eff) we have to modify the function names so not #
        #  all calls are to smatrix()                                                #
        #                                                                        #
        #-----------------------------------------------------------------------#
           if self._do_direct_detection:
                   print "Changing the direct detection matrix element files so function names don't conflict: "
                   if (whichME=='full'):
                           processes = self._dd_process_names
                   elif (whichME=='eff'):
                           processes = self._dd_eff_process_names
                   elif (whichME=='total'):
                           processes = self._dd_tot_process_names
                   else:
                           print 'Error: Unknown option for modifying direct detection matrix element files!'
                           sys.exit(1)

                   ii = 1
                   for process in processes:

                           #For effective matrix elements and total we ground even and odd together.
                           #We then need to label them (1-12 Spin independent, 13-24 Spin Dependent)
                           SI_OR_SD = 'SI'
                           if  ii > len(processes) / 2:
                                         SI_OR_SD = 'SD'

                           #print process
                           if (whichME=='full'):
                                        matrixfile_name = 'DD_matrix_'+process+'.f'
                           elif (whichME=='eff'):
                                        matrixfile_name = 'DD_matrix_eff_'+process+'.f'
                           elif (whichME=='total'):
                                        matrixfile_name = 'DD_matrix_tot_'+process+'.f'
                           else:
                                        print 'Error: Unknown option for modifying direct detection matrix element files!'
                                        sys.exit(1)

                           print "Opening: "+matrixfile_name

                           #The name of the file which returns only the matrix element value.
                           shutil.move('Projects/'+self._projectname+'/matrix_elements/'+matrixfile_name,\
                                  'Projects/'+self._projectname+'/matrix_elements/dd_temp.f')
                           matrixfile = open('Projects/'+self._projectname+'/matrix_elements/dd_temp.f','r')
                           matrixfileTo = open('Projects/'+self._projectname+'/matrix_elements/'+matrixfile_name,'w')

                           matrixfile_lines = matrixfile.readlines()

                           __dd_flag1 = 'MATRIX_'

                           for line in matrixfile_lines:
                                if __dd_flag1 in line:
                                  if (whichME=='full'):
                                                line= line
                                  elif (whichME=='eff'):
                                                #print SI_OR_SD
                                                line =line.replace(__dd_flag1, __dd_flag1+'EFF_'+SI_OR_SD+'_')
                                  elif (whichME=='total'):
                                                #print SI_OR_SD
                                                line =line.replace(__dd_flag1, __dd_flag1+'TOT_'+SI_OR_SD+'_')

                                matrixfileTo.write(line)

                           os.remove('Projects/'+self._projectname+'/matrix_elements/dd_temp.f')
                           matrixfile.close()
                           matrixfileTo.close()
                           ii=ii+1

        #-----------------------------------------------------------------------#
        def Write_makefile(self):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  This routine creates the makefile for compiling all the matrix        #
        #  elements generated by madgraph as well as the overall maddm makefile        #
        #                                                                        #
        #-----------------------------------------------------------------------#

          __flag = '#MAKEFILE_MADDM'
          #FOR THE MADDM MAKEFILE
          makefile = open('Projects/'+self._projectname+'/makefile','w')
          makefile_template = open('Projects/'+self._projectname+'/makefile_template','r')

          suffix = ''
          if self._do_relic_density and not self._do_direct_detection:
                  suffix = 'relic_density'
          elif not self._do_relic_density and self._do_direct_detection:
                  suffix = 'direct_detection'

          makefile_lines = makefile_template.readlines()
          for line in makefile_lines:
                if __flag in line:
                        new_line = '\tcd src/ && make '+suffix
                        new_line = new_line + '\n'
                        makefile.write(new_line)
                else:
                        makefile.write(line)

          makefile.close()
          makefile_template.close()

          #FOR THE MATRIX ELEMENTS MAKEFILE
          # Flags that are needed to create the matrix element makefile and smatrix.f files
          # Creates the make file in the matrix_elements folder.
          # includes the object files for all the individual matrix elements
          makefile = open('Projects/'+self._projectname+'/matrix_elements/makefile','w')
          makefile_template = open('Projects/'+self._projectname+'/matrix_elements/makefile_template','r')

          makefile_lines = makefile_template.readlines()
          for line in makefile_lines:
                if __flag in line:
                        new_line = 'objs = smatrix.o'

                        if self._do_relic_density == True:
                        # Add all the annihilation matrix.o files to the makefile
                          for i in range(len(self._dm_particles)):
                                for j in range(i,len(self._dm_particles)):
                                        for k in range(self._ann_nprocesses[i][j]):
                                                new_line = new_line+' matrix_'+self._ann_process_names[i][j][k]+'.o'

                          # Add all the DM -> DM matrix.o files to the makefile
                          for i in range(len(self._dm_particles)):
                                for j in range(i,len(self._dm_particles)):
                                        for k in range(self._dm2dm_nprocesses[i][j]):
                                                new_line = new_line+' matrix_'+self._dm2dm_process_names[i][j][k]+'.o'

                          # Add all the DM/SM scattering matrix.o files to the makefile
                          for i in range(len(self._dm_particles)):
                                for j in range(len(self._dm_particles)):
                                        for k in range(self._scattering_nprocesses[i][j]):
                                                new_line = new_line+' matrix_'+self._scattering_process_names[i][j][k]+'.o'

                          #Add all the Direct Detection matrix.o files to the makefile
                        if (self._do_direct_detection == True):
                                  for i in range(len(self._scatteringDirDetect_me[0])):
                                          new_line = new_line+' DD_matrix_'+self._dd_process_names[i]+'.o'

                                  for i in range(len(self._scatteringDirDetect_eff_me[0])):
                                          new_line = new_line+' DD_matrix_eff_'+self._dd_eff_process_names[i]+'.o'

                                  for i in range(len(self._scatteringDirDetect_tot_me[0])):
                                          new_line = new_line+' DD_matrix_tot_'+self._dd_tot_process_names[i]+'.o'

                        new_line = new_line + '\n'
                        makefile.write(new_line)
                else:
                        makefile.write(line)

          makefile.close()
          makefile_template.close()
          os.remove('Projects/'+self._projectname+'/makefile_template')


        #-----------------------------------------------------------------------#
        def check_param_card(self, path):
                """
                1) Check that no scan parameter are present
                2) Check that all the width are define in the param_card.
                - If a scan parameter is define. create the iterator and recall this fonction
                  on the first element.
                - If some width are set on 'Auto', call the computation tools."""

                pattern_scan = re.compile(r'''^[\s\d]*scan''', re.I+re.M)
                pattern_width = re.compile(r'''decay\s+(\+?\-?\d+)\s+auto(@NLO|)''',re.I)
                text = open(path).read()


                if pattern_scan.search(text):
                        # at least one scan parameter found. create an iterator to go trough the cards
                        main_card = check_param_card.ParamCardIterator(text)
                        self.param_card_iterator = main_card.__iter__()
                        first_card = self.param_card_iterator.next()
                        first_card.write(path)
                        return self.check_param_card(path)

                pdg_info = pattern_width.findall(text)
                if pdg_info:
                            print "Automatically computing the width of particle with PDF code "+str(pdg_info)
                            logging.info('Computing the width set on auto in the param_card.dat')
                            has_nlo = any(nlo.lower()=="@nlo" for _,nlo in pdg_info)
                            pdg = [pdg for pdg,nlo in pdg_info]
                            if not has_nlo:
                            		#--body_decay=2 for only 2 body decay
                                    self._MG5Cmd.do_compute_widths('%s --path=%s' % (' '.join(pdg), path))
                            else:
                                    self._MG5Cmd.do_compute_widths('%s %s --nlo' % (' '.join(pdg), path))

    	#-------------------------------------------------------------------#
	# Routine which writes the maddm.inc file. It takes the existing file and adds
	# the information about the location of s-channel resonances.
	#---------------------------------------------------------------------------------
	def WriteMadDMinc(self):

			print "(PDG, Mass, Width) of Resonances:"
			print self._resonances

			incfile_template = open('Templates/include/maddm.inc','r').read()
			res_dict = {}

# 			chunk_size = 5
 			init_lines=[]
# 			for i,k in enumerate(self._resonances):
# 				init_lines.append("resonances(%d)=%s"%(i+1,k[1]))

			res_dict['resonances_initialization'] = '\n'.join(init_lines)
			res_dict['n_resonances'] = len(self._resonances)

			writer = writers.FortranWriter(os.path.join(self._projectpath, 'include', 'maddm.inc'))
			writer.write(incfile_template % res_dict)
			writer.close()

	#---------------------------------------------------------------------------------


        #-----------------------------------------------------------------------#
        def Calculate(self, output_file='output/maddm.out'):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  Runs the maddm.x code created by CreateNumericalSession()                #
        #  and extracts the output values                                        #
        #  See the list in the 'return' statement for the values which are         #
        #  printed out. All cross sections at this level are in GeV^-2.         #
        #  WARNING: Be careful with x_freezeout. See the MadDM v.1.0 manual        #
        #  for details on how it is calculated. Do not expect it to match         #
        #  with micrOMEGAs, due to a difference in definitions.                        #
        #                                                                        #
        #-----------------------------------------------------------------------#

          cwd = os.getcwd()
          cmd = self._projectpath+'/maddm.x '+output_file
          exec_dir = self._projectpath
          exec_cmd = os.path.join(cwd, cmd)
          self.check_param_card(self._paramcard)
          subprocess.call(exec_cmd, cwd = exec_dir, shell=True)

          #Here we read out the results which the FORTRAN module dumped into a file
          #called 'omega'. The format is such that the first line is always relic density
          # , second line is the nucleon scattering cross section (SI) for proton, third is SI
          # nucleon cross section for the neutron etc. If a quantity was not calculated, we output -1
          result = open(self._projectpath+'/'+output_file, 'r')
          omegah2 = float(result.readline().split()[1])
          x_freezeout = float(result.readline().split()[1])
          wimp_mass = float(result.readline().split()[1])
          sigmav_xf = float(result.readline().split()[1])
          sigmaN_SI_proton = float(result.readline().split()[1])
          sigmaN_SI_neutron = float(result.readline().split()[1])
          sigmaN_SD_proton = float(result.readline().split()[1])
          sigmaN_SD_neutron = float(result.readline().split()[1])
          Nevents  = int(result.readline().split()[1])
          sm_switch = int(result.readline().split()[1])
          result.close()
          #Set up the results in the darkmatter object
          self._omegah2 = omegah2
          self._x_freezeout = x_freezeout
          self._sigmav_xf = sigmav_xf
          self._sigmaN_SI_proton = sigmaN_SI_proton
          self._sigmaN_SI_neutron = sigmaN_SI_neutron
          self._sigmaN_SD_proton = sigmaN_SD_proton
          self._sigmaN_SD_neutron = sigmaN_SD_neutron
          self._Nevents = Nevents
          self._dm_mass = wimp_mass


          return [omegah2, x_freezeout, wimp_mass, sigmav_xf , sigmaN_SI_proton \
                        ,sigmaN_SI_neutron, sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch]


        #-----------------------------------------------------------------------#
        def ChangeParameter(self, flag, value):
        #-----------------------------------------------------------------------#
        #                                                                        #
        #  Changes the parameter in the param card.                                #
        #  flag is a string name of the parameter. Value is the numerical value #
        #  that the parameter should be set to.                                        #
        #                                                                              #
        #-----------------------------------------------------------------------#

          try:
           foundflag = False
           #Open the param card to read in all the lines
           finput = open(self._paramcard, 'r')
           lines = finput.readlines()
           finput.close()
           #Open it again to rewrite it. Necessary in order to
           #not leave junk at the end.
           finput = open(self._paramcard, 'w')
           for line in lines:#fileinput.FileInput(self._paramcard,inplace=1):
            line_list = line.split()
            if flag in line_list:
                #Make sure that the flag is not a part of another variable name
                foundflag = True
                index_of_value = line_list.index(flag)-2
                newline=''
                #Find where the value is in the string
                #Put the line back together
                for j in range(0, index_of_value):
                   if j==0:
                    newline=newline+line_list[j]
                   else:
                    newline=newline+' '+line_list[j]
                newline=newline+' '+str(value)
                for j in range(index_of_value+1, len(line_list)):
                   newline=newline+' '+line_list[j]
                newline=newline+'\n'
                finput.write(newline)
                #print "New Line: "+newline
            else:
            	#print "Old Line: "+line
                finput.write(line)

           finput.close()
           if not foundflag:
                print "ERROR: Parameter "+flag+" not found!"

          except OSError, error:
                logging.debug(error)
                print "ERROR: Could not open the Param Card file!"
                sys.exit(1)

        #-----------------------------------------------------------------------#
        def init_param_scan(self,name_of_script):
        #-----------------------------------------------------------------------#
        #                                                                       #
        #  Starts the param_scan.py program. The function prompts the user to   #
        #  edit the param_scan.py script and set up the parameter scan.         #
        #   Editing the param_scan.py file requires basic knowledge of Python   #
        #-----------------------------------------------------------------------#

                #Edit the parameter scan script
                script_path = self._projectpath+'/'+name_of_script
                if os.path.isfile(script_path):
                    print 'WARNING: The parameter scan under the name '+name_of_script+' already exists. Overwrite? [y] (y/n)'
                    answer = raw_input()
                    if answer =='y' or answer=='':
                        shutil.copy2('param_scan_default.py', script_path)
                        subprocess.call('vi '+self._projectpath+'/%s' % name_of_script, shell= True)
                    elif answer =='n':
                            print 'Continuing to the existing parameter scan.'
                    elif answer !='n' and answer !='y' and answer !='':
                        print 'Not a valid answer... continuing to the existing parameter scan.'
                else:
                    shutil.copy2('param_scan_default.py', script_path)
                    subprocess.call('vi '+self._projectpath+'/%s' % name_of_script, shell= True)

                curr_dir = os.getcwd()
                os.chdir(self._projectpath)

                cmd = os.path.join('./',name_of_script)
                subprocess.call(cmd, shell=True)
                os.chdir(curr_dir)

        #---------------------------------------------------------------------------#
        def is_excluded(self, omega_min = 0., omega_max = 0.1, si_limit_x = [], si_limit_y= []):
        #---------------------------------------------------------------------------#
        #                This function determines whether a model point is excluded or not   #
        #                based on the resulting relic density, spin independent and spin         #
        #                dependent cross sections. The user can supply min/max values for    #
        #                relic density, as well as the exclusion curves for dm nucleon                #
        #                scattering cross section.                                                                                        #
        #                WARNING: The function takes the values for cross sections in pb                #
        #---------------------------------------------------------------------------#
                excluded_by_omegah2 = False
                excluded_by_dir_detect = False

                sigma_si_ave = 0.5*(self._sigmaN_SI_proton + self._sigmaN_SI_neutron)*GeV2pb
                #Check relic density
                if (self._omegah2 < omega_min or self._omegah2 > omega_max):
                        excluded_by_omegah2 = True
                #Check si cross section limit
                if (sigma_si_ave > interp(self._dm_mass, si_limit_x, si_limit_y)):
						excluded_by_dir_detect = True
                #Check sd cross section limit
                #else if (self._sigma_sd_ave > sigma_sd_limit(self._dm_mass)):
                #        return True

                excluded = excluded_by_omegah2 or excluded_by_dir_detect
                return [excluded, excluded_by_omegah2, excluded_by_dir_detect]

        #------------------------------------------------------------------------------#
        def fix_param_card(self):
        #------------------------------------------------------------------------------#
        #  A function which makes sure that the parameter card contains all the params.#
        #  from the effective operator model. This is not an essential function, but   #
        #  it prevents some annoying Warning messages from appearing during the exec-  #
        #  ution of the code.                                                                                                                   #
        #------------------------------------------------------------------------------#
                default_param_card = self._projectpath+'/Cards/param_card_default.dat'
                temp_param_card = self._projectpath+'/Cards/param_card_temp.dat'
                shutil.copy2(default_param_card, temp_param_card)
                print "Parameter card: "+self._paramcard
                try:
                        #open the parameter card file
                        finput = open(self._paramcard, 'r')
                        #open the default parameter card
                        finput_temp = open(temp_param_card, 'r+')
                        lines = finput.readlines()
                        lines_temp = finput_temp.readlines()
                        #Move the temp param card pointer to 0, so we can overwrite things
                        finput_temp.seek(0,0)
                except OSError, error:
                        logging.debug(error)
                        print "Error: Could not open the Param Card file!"
                        sys.exit(1)
                #Loop through all the lines in the param card
                #and for each parameter, make sure that it is the same in the def. param card.
                #If there is a difference, change the parameter in the default_card
                #This works because the default card always has all parameters but the
                #one which is copied may not.
                for line_temp in lines_temp:
                        line_temp_list = line_temp.split()
                        if (line_temp[0]=='#' or len(line_temp_list) < 4 \
                                        or line_temp_list[0] == 'Block' or line_temp_list[0] =='BLOCK'):
                                        finput_temp.write(line_temp)
                                        continue

                        #The parameter name is always after the '#' character
                        parname_index= line_temp_list.index('#')+1
                        param_name = line_temp_list[parname_index]
                        #Now look for this parameter name in the def. param card and see if the value is
                        #the same
                        for line in lines:
                                line_list = line.split()
                                if param_name in line_list:
                                        #print line,
                                        #print ":::"+param_name
                                        line_to_write = line
                                        break
                                else:
                                        line_to_write = line_temp

                        #print "wrote: "+line_to_write,
                        finput_temp.write(line_to_write)

                finput.close()
                finput_temp.close()
                try:
                        shutil.copy2(temp_param_card, self._paramcard)
                except OSError, error:
                        logging.debug(error)
                        print "Error: Could not open the Param Card file!"
                        sys.exit(1)

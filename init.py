import os
import sys
import shutil
import readline
import rlcompleter

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    GRAY = '\033[90m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#Root path
rp = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.append(rp)


model_list2 = os.listdir(os.path.join(rp, 'models'))

def completer(text, state):
    options = [x for x in model_list2 if x.startswith(text)]
    try:
        return options[state]
    except IndexError:
        return None

readline.set_completer(completer)
if 'libedit' in readline.__doc__:
    readline.parse_and_bind("bind ^I rl_complete")
else:
    readline.parse_and_bind("tab: complete")

#readline.parse_and_bind("tab: complete")

#-----------------------------------------------------------------------#
def initialize_MadDM_session(print_banner = False):
#-----------------------------------------------------------------------#
#                                                                       #
#  This routine starts the whole program.  The desired model and        #
#  project name are input by the user.  The model is checked against    #
#  the model directory in the main madgraph folder.                     #
#                                                                       #
#-----------------------------------------------------------------------#

  if (print_banner):
    print "\n"
    print \
  "            ====================================================\n"+\
  "            |                  "+bcolors.OKBLUE+"  MadDM v2.0                     "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "                                                                               \n"+\
  "                #########            Basics tutorial:  susy.phsx.ku.edu/~mihailo \n"+\
  "             ###\\\\####//#####              Launchpad:  launchpad.net/maddm      \n"+\
  "           ######\\\\##//########                                              \n"+\
  "          ########\\\\//###########                                            \n"+\
  "         #########//\\\\############                                           \n"+\
  "        #########//##\\\\############      "+bcolors.FAIL+"             arXiv:1308.4955        \n"+bcolors.ENDC+\
  "       ########//#####\\\\###########                   "+bcolors.FAIL+"arXiv:1505.04190        \n"+bcolors.ENDC+\
  "       ######################### ## ___________________________________________\n"+\
  "       ####################### 0  # "+bcolors.OKGREEN+" _     _               _  _____   _     _  \n"+bcolors.ENDC+\
  "       #############   0  ###    ## "+bcolors.OKGREEN+"| \   / |   ___    ___|| | ___ \ | \   / | \n"+bcolors.ENDC+\
  "       ##############    #########  "+bcolors.OKGREEN+"||\\\\ //|| / __ |  / __ | ||   || ||\\\\ //|| \n"+bcolors.ENDC+\
  "        ##########################  "+bcolors.OKGREEN+"||  V  || ||__||  ||__|| ||___|| ||  V  || \n"+bcolors.ENDC+\
  "         ###################   ##   "+bcolors.OKGREEN+"||     || \_____\ \____| |_____/ ||     || \n"+bcolors.ENDC+\
  "          ############       ###    ___________________________________________\n"+\
  "           ##########    ######                                                 \n"+\
  "             ################                                                   \n"+\
  "                 ########                                                       \n"+\
  "                                                                                    \n"

  # Get the model name and check if it exists in the Models folder
  model_list = os.listdir(os.path.join(rp, 'models'))

  # Loop until the model name is ok
  leg_ans = False
  while not leg_ans:

	  model_name = raw_input('Enter model name (press Enter to list the available models): ')

	  if (model_name == ''):
		  # List the available DM models in the Models directory of MadGraph

		  for mdl in model_list:
			  mdl_dir = os.path.join('models', mdl)
			  if os.path.isdir(os.path.join(rp, mdl_dir)):
				  print mdl+"   "

	  elif model_name in model_list:
			leg_ans = True
	  else:
		  print "The model you entered is not available! Please try again"

  # Loop until the project name is set
  project_name = raw_input("Enter project name ("+model_name+"): ")
  if (project_name == ''):
      project_name = model_name

  #Ask the user which calculations they would like to perform (relic density, direct detection)
  do_direct_detection = False
  do_relic_density = False
  leg_ans = False
  while not leg_ans:

  	answer = raw_input("Calculate relic density? [y] (y/n)")
  	if answer == 'y' or answer=='Y' or answer =='':
  		do_relic_density = True
  		leg_ans = True
  	if answer == 'n' or answer=='N':
  		do_relic_density = False
  		leg_ans=True

  leg_ans = False
  while not leg_ans:

  	answer = raw_input("Calculate direct detection cross sections? [y] (y/n)")
  	if answer == 'y' or answer=='Y' or answer =='':
  		do_direct_detection = True
  		leg_ans = True
  	if answer == 'n' or answer=='N':
  		do_direct_detection = False
  		leg_ans=True

  do_directional_detection = False
  if do_direct_detection == True:
  	leg_ans = False
  	while not leg_ans:
  		answer = raw_input("Simulate DM scattering events off nuclei? [y] (y/n)")
  		if answer == 'y' or answer=='Y' or answer =='':
  			do_directional_detection = True
  			leg_ans = True
  		if answer == 'n' or answer=='N':
  			do_directional_detection = False
  			leg_ans=True

  names = [model_name,project_name, do_relic_density, do_direct_detection, do_directional_detection]
  return names

import collections
import logging
import os
import re
import sys

import MGoutput

MG5MODE = True
import madgraph.various.misc as misc
import madgraph.interface.extended_cmd as cmd
import madgraph.various.banner as banner_mod
import madgraph.interface.common_run_interface as common_run
import madgraph.iolibs.files as files
import models.check_param_card as param_card_mod
        
#import darkmatter as darkmatter

pjoin = os.path.join
logger = logging.getLogger('madgraph.plugin.maddm')

MDMDIR = os.path.dirname(os.path.realpath( __file__ ))

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

#===============================================================================
# CommonRunCmd
#===============================================================================
class MADDMRunCmd(cmd.Cmd):
    
    intro_banner=\
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
  "                                                                                    \n"+\
  "            ====================================================\n"+\
  "            |           "+bcolors.OKBLUE+" A MadGraph5_aMC@NLO plugin.            "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "" 
    
    
    def __init__(self, dir_path, options, *args, **opts):
        """common"""
  
        cmd.Cmd.__init__(self, *args, **opts)
                # Define current MadEvent directory
        if dir_path is None and not MG5MODE:
            dir_path = root_path

        self.dir_path = dir_path
        
        # Define self.proc_characteristics (set of information related to the
        # directory status
        self.get_characteristics()
        
        
    def get_characteristics(self, path=None):
        """reads the proc_characteristics file and initialises the correspondant
        dictionary"""
        
        if not path:
            path = os.path.join(self.dir_path, 'matrix_elements', 'proc_characteristics')
        
        self.proc_characteristics = MGoutput.MADDMProcCharacteristic(path)
        return self.proc_characteristics
    
    
    def do_launch(self, line):
        """run the code"""

        args = line.split()
        self.ask_run_configuration(mode=[])
        self.compile()
        
        nb_output = 1
        output = pjoin(self.dir_path, 'output', 'maddm_%s.out') 
        while os.path.exists(output % nb_output):
            nb_output +=1
        else:
            output = pjoin('output', 'maddm_%s.out') % nb_output
        
        misc.call(['maddm.x', output], cwd =self.dir_path)

        #Here we read out the results which the FORTRAN module dumped into a file
        #called 'omega'. The format is such that the first line is always relic density
        # , second line is the nucleon scattering cross section (SI) for proton, third is SI
        # nucleon cross section for the neutron etc. If a quantity was not calculated, we output -1
        
        result = []
        for line in open(pjoin(self.dir_path, output)):
            misc.sprint(line)
            result.append(float(line.split()[1]))

        (omegah2, x_freezeout, wimp_mass, sigmav_xf , sigmaN_SI_proton \
                        ,sigmaN_SI_neutron, sigmaN_SD_proton, sigmaN_SD_neutron, 
                        Nevents, sm_switch) = result
        Nevents  = int(Nevents)
        sm_switch = int(sm_switch)
        
        misc.sprint(result)
        
        return [omegah2, x_freezeout, wimp_mass, sigmav_xf , sigmaN_SI_proton \
                        ,sigmaN_SI_neutron, sigmaN_SD_proton, sigmaN_SD_neutron, Nevents, sm_switch]

        
        
        
        
        
    
    def ask_run_configuration(self, mode=None):
        """ask the question about card edition / Run mode """
        self.mode = self.ask('', '0', mode=mode, data=self.proc_characteristics, 
                        ask_class=MadDMSelector, timeout=60)
        misc.sprint(self.mode)
        self.maddm_card = MadDMCard(pjoin(self.dir_path, 'Cards', 'maddm_card.dat'))
        self.param_card = param_card_mod.ParamCard(pjoin(self.dir_path, 'Cards', 'param_card.dat'))
        for key, value in self.mode.items():
            if value == 'ON':
                self.mode[key] = True
            else:
                self.mode[key] = False
        
        logger.info("Start computing %s" % ','.join([name for name, value in self.mode.items() if value]))
        return self.mode

    def compile(self):
        """compile the code"""

        if self.mode['relic'] and self.mode['direct']:
            misc.compile(cwd=self.dir_path)
        elif self.mode['relic'] and not self.mode['direct']:
            misc.compile(['relic_density'],cwd=self.dir_path)
        elif self.mode['direct'] and not self.mode['relic_density']:
            misc.compile(['direct_detection'],cwd=self.dir_path)
        else:
            raise Exception, "No computation requested. End the computation"
    
        if os.path.exists(pjoin(self.dir_path, 'src', 'maddm.x')) or os.path.exists(pjoin(self.dir_path, 'maddm.x')):
            misc.sprint("compilation done")
        else:
            raise Exception, 'Compilation of maddm failed'
        
    ############################################################################
    def do_open(self, line):
        """Open a text file/ eps file / html file"""

        args = self.split_arg(line)
        # Check Argument validity and modify argument to be the real path
        self.check_open(args)
        file_path = args[0]

        misc.open_file(file_path)
    
    ############################################################################
    def check_open(self, args):
        """ check the validity of the line """

        if len(args) != 1:
            self.help_open()
            raise self.InvalidCmd('OPEN command requires exactly one argument')

        if args[0].startswith('./'):
            if not os.path.isfile(args[0]):
                raise self.InvalidCmd('%s: not such file' % args[0])
            return True

        # if special : create the path.
        if not self.dir_path:
            if not os.path.isfile(args[0]):
                self.help_open()
                raise self.InvalidCmd('No MadEvent path defined. Unable to associate this name to a file')
            else:
                return True

        path = self.dir_path
        if os.path.isfile(os.path.join(path,args[0])):
            args[0] = os.path.join(path,args[0])
        elif os.path.isfile(os.path.join(path,'Cards',args[0])):
            args[0] = os.path.join(path,'Cards',args[0])
        # special for card with _default define: copy the default and open it
        elif '_card.dat' in args[0]:
            name = args[0].replace('_card.dat','_card_default.dat')
            if os.path.isfile(os.path.join(path,'Cards', name)):
                files.cp(os.path.join(path,'Cards', name), os.path.join(path,'Cards', args[0]))
                args[0] = os.path.join(path,'Cards', args[0])
            else:
                raise self.InvalidCmd('No default path for this file')
        elif not os.path.isfile(args[0]):
            raise self.InvalidCmd('No default path for this file')



class MadDMSelector(common_run.EditParamCard):
    """ """

    def cmdloop(self, intro=None):
        super(MadDMSelector,self).cmdloop(intro)
        return self.run_options

    
    def __init__(self, *args, **opts):

        #0. some default variable
        process_data = opts.pop('data', collections.defaultdict(bool))
        self.run_options = {'relic': 'ON',
                            'direct': 'OFF' if process_data['has_direct_detection'] else 'Not available',
                            'directional': 'OFF' if process_data['has_directional_detection'] else 'Not available'}
        
        #1. Define what to run and create the associated question
        mode = opts.pop('mode', None)  
        if mode:
            for key, value in mode.items():
                if self.run_options[key] in ['ON', 'OFF']:
                    self.run_options[key] = value      
        question = self.create_question()            
                    
        #2. Define the list of allowed argument
        allow_args = [str(0), 'done', str(4), 'param', str(5), 'maddm']
        for i,key in enumerate(['relic', 'direct', 'directional']):
            if self.run_options[key] in ['ON', 'OFF']:
                allow_args.append(str(i+1))
                allow_args.append(key)
        opts['allow_arg'] = allow_args
                
        # 3.initialise the object         
        param_card_path = pjoin(opts['mother_interface'].dir_path, 'Cards', 'param_card.dat')
        super(MadDMSelector,self).__init__(question, [param_card_path], **opts)

        self.me_dir = opts['mother_interface'].dir_path
        self.define_paths(pwd=self.me_dir) 
        
        #4. Initialise the maddm cards
        self.maddm_def = MadDMCard(self.paths['maddm_default'])
        try:
            self.maddm = MadDMCard(self.paths['maddm'])
        except Exception as e:
            logger.error('Current maddm_card is not valid. We are going to use the default one.')
            logger.error('problem detected: %s' % e) 
            files.cp(self.paths['maddm_default'], self.paths['maddm'])
            self.maddm = MadDMCard(self.paths['maddm'])
            
        self.maddm_set = self.maddm_def.keys() + self.maddm_def.hidden_param
        for var in self.pname2block:
            if var in self.maddm_set:
                self.conflict.append(var)  
            
        
    digitoptions = {1: 'relic', 2:'direct', 3:'directional'}  
    def create_question(self):
        """create the new question depending of the status"""
        
        def get_status_str(status):
            if status == 'ON':
                return bcolors.OKGREEN + 'ON' + bcolors.ENDC
            elif status == 'OFF':
                return bcolors.FAIL + 'OFF' + bcolors.ENDC
            else:
                return bcolors.WARNING + status + bcolors.ENDC
        
        question = """\n%(start_green)sHere is the current status of requested run %(stop)s: 
 * Enter the name/number to (de-)activate the corresponding feature
    1. Compute the Relic Density     %(start_underline)srelic%(stop)s       = %(relic)s
    2. Compute Direct Detection      %(start_underline)sdirect%(stop)s      = %(direct)s
    3. Compute Directional Detection %(start_underline)sdirectional%(stop)s = %(directional)s    
%(start_green)s You can also edit the various input card%(stop)s:
 * Enter the name/number to open the editor
 * Enter a path to a file to replace the card
 * Enter %(start_bold)set NAME value%(stop)s to change any parameter to the requested value 
    4. Edit the model parameters    [%(start_underline)sparam%(stop)s]
    5. Edit the MadDM options       [%(start_underline)smaddm%(stop)s]\n""" % \
    {'start_green' : '\033[92m',
     'stop':  '\033[0m',
     'start_underline': '\033[4m',
     'start_bold':'\033[1m', 
     'relic': get_status_str(self.run_options['relic']),
     'direct':get_status_str(self.run_options['direct']),
     'directional':get_status_str(self.run_options['directional'])
     }
        return question

    def define_paths(self, **opt):
        
        super(MadDMSelector, self).define_paths(**opt)
        misc.sprint("pass in define_paths")
        self.paths['maddm'] = pjoin(self.me_dir,'Cards','maddm_card.dat')
        self.paths['maddm_default'] = pjoin(self.me_dir,'Cards','maddm_card_default.dat')
        
    
    
    def default(self, line):
        """Default action if line is not recognized"""
        
        line = line.strip()
        args = line.split()
        if args:
            if args[0].isdigit():
                val = int(args[0])
                if val in [1,2,3]:
                    args[0] = self.digitoptions[val]
                elif val == 4:
                    self.open_file('param')
                    self.value = 'repeat'
                elif val == 5:
                    self.open_file('maddm')
                    self.value = 'repeat'
                elif val !=0:
                    logger.warning("Number not supported: Do Nothing") 
            if '=' in line:
                if '=' in args:
                    args.remove('=')
                if '=' in args[0]:
                    args = args[0].split('=')
                tag, value = args
                if self.run_options[tag] in ['ON', 'OFF']:
                    if value.lower() == 'on':
                        self.run_options[tag] = 'ON'
                    elif value.lower() == 'off':
                        self.run_options[tag] = 'OFF'
                    else:
                        logger.warning('%s is not a valid entry')
                self.value = 'repeat'
            elif args[0] in self.run_options:
                if self.run_options[args[0]] == 'ON':
                    self.run_options[args[0]] = 'OFF'
                elif self.run_options[args[0]] == 'OFF':
                    self.run_options[args[0]] = 'ON'
                else:
                    logger.warning('This entry can not be changed for this running directory.')
                self.value = 'repeat'
        if self.value == 'repeat':
            self.question = self.create_question()
            return line
        else:
            return super(MadDMSelector, self).default(line)
        
    def open_file(self, path):
        
        path = super(MadDMSelector,self).open_file(path)
        
        if path == self.paths['maddm']:
            try:
                self.maddm = MadDMCard(path) 
            except Exception as e:
                logger.error('Current param_card is not valid. We are going to use the default one.')
                logger.error('problem detected: %s' % e)
                logger.error('Please re-open the file and fix the problem.')
                logger.warning('using the \'set\' command without opening the file will discard all your manual change')
    
    def complete_set(self, text, line, begidx, endidx, formatting=True):
        """ Complete the set command"""

        possibilities = super(MadDMSelector,self).complete_set(text, line, begidx, endidx, formatting=False)
        args = self.split_arg(line[0:begidx])
        if len(args)>1 and args[1] == 'maddm_card':
            start = 2
        else:
            start = 1 
            if len(args) ==1:
                possibilities['category of parameter (optional)'] += \
                                self.list_completion(text, ['maddm_card'], line)
                
        if len(args)==start:
            correct = self.maddm_set + ['default']
            possibilities['maddm Card'] = self.list_completion(text, correct, line)   
        elif len(args)==start+1:
            possibilities['maddm Card'] = ['default']
            
        return self.deal_multiple_categories(possibilities, formatting)

    def do_set(self, line):
        """ edit the value of one parameter in the card"""
        
        args = self.split_arg(line)
        # fix some formatting problem
        if '=' in args[-1]:
            arg1, arg2 = args.pop(-1).split('=')
            args += [arg1, arg2]
        if '=' in args:
            args.remove('=')
        args = [ a.lower() for a in args]
        
        misc.sprint(args)
        if args[0] in ['maddm_card']:
            start = 1
            if args[1] == 'default':
                logging.info('replace %s by the default card' % args[0])
                self.maddm = MadDMCard(self.paths['maddm_default'])
        elif args[0] in self.maddm_set and args[0] not in self.conflict:
            start = 0
        else:
            return super(MadDMSelector, self).do_set(line)
        misc.sprint(args, start)
        if args[start+1] == 'default':
            default = self.maddm_def[args[start]]
            self.setDM(args[start], default)
        else:
            self.setDM(args[start], args[start+1])
        # write the new file
        self.maddm.write(self.paths['maddm'])
        
    def setDM(self, name, value):
        logger.info('modify parameter %s of the maddm_card.dat to %s' % (name, value))
        self.maddm.set(name, value, user=True)
        
        
            
        



class InvalidMaddmCard(banner_mod.InvalidRunCard):
    pass
             
class MadDMCard(banner_mod.RunCard):
    """MadDM card use the read/write of the runcard object"""
    
    def __new__(cls, finput=None):
        """Bypass the standard RunCard one"""
        return super(banner_mod.RunCard, cls).__new__(cls, finput)
    
    def default_setup(self):
        """define the default value"""
        
        self.add_param('print_out', False)
        self.add_param('print_sigmas', False)
        
        self.add_param('relic_canonical', True)
        self.add_param('do_relic_density', True, system=True)
        self.add_param('do_direct_detection', False, system=True)
        self.add_param('do_directional_detection', False, system=True)
        
        self.add_param('calc_taacs_ann_array', True)
        self.add_param('calc_taacs_dm2dm_array', True)
        self.add_param('calc_taacs_scattering_array', True)
        
        self.add_param('eps_ode', 0.01)
        self.add_param('xd_approx', True)
        
        self.add_param('x_start', 50.0)
        self.add_param('x_end', 1000.0)
        self.add_param('dx_step', 1.0)
        
        self.add_param('ngrid_init', 20)
        self.add_param('nres_points', 20)
        self.add_param('eps_wij', 0.01)
        self.add_param('iter_wij', 2)
        
        # Direct Detection nucleon form factors (fPx - proton, fNx - neutron)
        #Scalar FF
        self.add_param('SPu', 0.0153)
        self.add_param('SPd',0.0191)
        self.add_param('SPs', 0.0447)
        self.add_param('SPg', 1 - 0.0153 - 0.0191 - 0.0447)
        self.add_param('SNu', 0.0110)
        self.add_param('SNd', 0.0273)
        self.add_param('SNs', 0.0447)
        self.add_param('SNg', 1 - 0.0110 - 0.0273 - 0.0447)
        # Vector FF
        self.add_param('VPu', 2.0)
        self.add_param('VPd', 1.0)
        self.add_param('VNu', 1.0)
        self.add_param('Vnd', 2.0)
        # Axial Vector FF
        self.add_param('AVPu',0.842)
        self.add_param('AVPd',-0.427)
        self.add_param('AVPs',-0.085)
        self.add_param('AVNu',-0.427)
        self.add_param('AVNd',0.842)
        self.add_param('AVNs',-0.085)
        # Sigma(mu,nu) FF
        self.add_param('SigPu',0.84)
        self.add_param('SigPd',-0.23)
        self.add_param('SigPs',-0.046)
        self.add_param('SigNu',-0.23)
        self.add_param('SigNd',0.84)
        self.add_param('SigNs',-0.046)
        #
        # For Directional detection and direct detection rates
        #
        self.add_param('material', 1)
        #Setting up the DM constants
        self.add_param('vMP', 220.0)
        self.add_param('vescape', 650.0)
        self.add_param('rhoDM', 0.3)
        #detector 
        self.add_param('detector_size', 1000.0)
        self.add_param('En_threshold', 4.0)
        self.add_param('lambd', 1.0)
        self.add_param('sig_theta', 3.0)
        
        self.add_param('En_min', 0.0)
        self.add_param('En_max', 100.0)
        self.add_param('cos_min', -1.0)
        self.add_param('cos_max', 1.0)
        self.add_param('day_min', 0.0)
        self.add_param('day_max', 365.0)
        self.add_param('Energy_bins', 100)
        self.add_param('cos_theta_bins', 20)
        self.add_param('day_bins', 10)
        
        self.add_param('smearing', False)
        
    def write(self, output_file, template=None, python_template=False):
        """Write the run_card in output_file according to template 
           (a path to a valid run_card)"""

        if not template:
            template = pjoin(MDMDIR, 'Templates', 'Cards', 'maddm_card.dat')
            python_template = True

        super(MadDMCard, self).write(output_file, template=template,
                                    python_template=python_template) 

    def check_validity(self):
        """ """
        
        if self['SPu'] + self['SPs'] + self['SPd'] + self['SPg'] -1 > 1e-3:
            raise InvalidMaddmCard , 'The sum of SP* parameter should be 1.0 get %s' % (self['SPu'] + self['SPs'] + self['SPd'] + self['SPg'])
        
        if self['SNu'] + self['SNs'] + self['SNd'] - self['SNg'] -1 > 1e-3:
            raise InvalidMaddmCard, 'The sum of SM* parameter should be 1.0 get %s' % (self['SNu'] + self['SNs'] + self['SNd'] + self['SNg'])
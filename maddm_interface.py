from __future__ import absolute_import
from __future__ import print_function
import logging
import os
import collections

from . import maddm_run_interface as maddm_run_interface
from . import dm_tutorial_text as dm_tutorial_text

import madgraph.core.base_objects as base_objects
import madgraph.core.helas_objects as helas_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.interface.master_interface as master_interface
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc
import madgraph.iolibs.files as files
#import darkmatter as darkmatter
import madgraph.interface.common_run_interface as common_run
import models.check_param_card as check_param_card
import models.model_reader as model_reader
from madgraph import MG5DIR

import re
import six
from six.moves import map
from six.moves import range
from six.moves import zip
pjoin = os.path.join

logger = logging.getLogger('madgraph.plugin.maddm')
logger_tuto = logging.getLogger('tutorial_plugin') # -> stdout include instruction in
                                                  #    order to learn maddm
GeV2pb = 3.894E8

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

try:
    import scipy
except ImportError:
    HAS_SCIPY = False 
else:
    HAS_SCIPY = True

try:
    import numpy
except ImportError:
    HAS_NUMPY = False
else:
    HAS_NUMPY = True

if HAS_NUMPY and HAS_SCIPY:
    from . import jfactor_gc as jfactor_gc

class DMError(Exception): pass

# Root path
MDMDIR = os.path.dirname(os.path.realpath( __file__ ))

class MadDM_interface(master_interface.MasterCmd):

    intro_banner=\
  "            ====================================================\n"+\
  "            |                  "+bcolors.OKBLUE+"  MadDM v3.2                     "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "                                                                               \n"+\
  "                #########                                                        \n"+\
  "             ###\\\\####//#####              Launchpad:  launchpad.net/maddm      \n"+\
  "           ######\\\\##//########                                              \n"+\
  "          ########\\\\//###########                     "+bcolors.FAIL+"arXiv:1308.4955        \n"+bcolors.ENDC+\
  "         #########//\\\\############                    "+bcolors.FAIL+"arXiv:1505.04190       \n"+bcolors.ENDC+\
  "        #########//##\\\\############                   "+bcolors.FAIL+"arXiv:1804.00044        \n"+bcolors.ENDC+\
  "       ########//#####\\\\###########                   "+bcolors.FAIL+"arXiv:2107.04598        \n"+bcolors.ENDC+\
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
  "            Need to learn? -> type tutorial                               \n"+\
  "                                                                                    \n"+\
  "            ====================================================\n"+\
  "            |           "+bcolors.OKBLUE+" A MadGraph5_aMC@NLO plugin.            "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "%s" 
    
    _define_options = ['darkmatter', 'coannihilator', 'benchmark']
    # process number to distinguish the different type of matrix element
    process_tag = {'DM2SM': 1999, 
                   'DM2DM': 1998,
                   'DMSM': 1997, # not use so far!
                   'DD': 1996,
                   'ID': 1995}
    
    eff_operators_SI = {1:'SIEFFS', 2:'SIEFFF', 3:'SIEFFV'}
    eff_operators_SD = {1:False, 2:'SDEFFF', 3:'SDEFFV'} 
    
    
    # for install command:
    _install_opts = list(master_interface.MasterCmd._install_opts)
    _install_opts.append('PPPC4DMID')
    _advanced_install_opts = list(master_interface.MasterCmd._advanced_install_opts)
    _advanced_install_opts +=['gsl', 'fitsio','dragon','dragon_data']
    
    install_ad = dict(master_interface.MasterCmd.install_ad)
    install_ad.update({'PPPC4DMID': ['1012.4515'],
                       'dragon':['0807.4730'],
                       'dragon_data':['1712.09755']}
                       )
    install_name =  dict(master_interface.MasterCmd.install_name)
    install_name.update({'PPPC4DMID':'PPPC4DMID'})
    install_name.update({'dragon_data_from_galprop':'dragon_data'})
    # path to PPPC4DMID
    
    options_configuration = dict(master_interface.MasterCmd.options_configuration)
    options_configuration['pppc4dmid_path'] = None
    options_configuration['dragon_path'] = None
    options_configuration['maddm_first_indirect'] = True
    
    _set_options = list(master_interface.MasterCmd._set_options)
    _set_options.append('dragon_data_from_galprop')
    
    
    def post_install_PPPC4DMID(self):
        if os.path.exists(pjoin(MG5DIR, 'PPPC4DMID')):
            self.options['pppc4dmid_path'] = pjoin(MG5DIR, 'PPPC4DMID')
        
        if not maddm_run_interface.HAS_SCIPY:
            logger.critical("Fermi-LAT and HESS limits calculation for indirect detection requires numpy and scipy. Please install the missing modules.")
            logger.info("you can try to use \"pip install scipy\"")

        self.exec_cmd('save options %s' % 'pppc4dmid_path')

        return
    
    def set_configuration(self, config_path=None, final=True, **opts):

        out = master_interface.MasterCmd.set_configuration(self, config_path, final, **opts)
        
        if final:
            if self.options['pppc4dmid_path'] and not os.path.exists(self.options['pppc4dmid_path']):
                self.options['pppc4dmid_path'] = None
            if self.options['dragon_path'] and not os.path.exists(self.options['dragon_path']):
                self.options['dragon_path'] = None                
            #self.do_save('options', to_keep={'pppc4dmid_path':'./PPPC4DMID'})
        
        return out
    
    def preloop(self, *args, **opts):
        super(MadDM_interface, self).preloop(*args, **opts)
        self.prompt = 'MadDM>'
    
    def change_principal_cmd(self, name, *args, **opts):
        out = super(MadDM_interface, self).change_principal_cmd(name,*args, **opts)
        self.prompt = 'MadDM>'
        return out
    
    def __init__(self, *args, **opts):
        
        super(MadDM_interface, self).__init__(*args, **opts)
        # DM candidates and coannihilators (on the basis of Z2 symmetry)
        self._dm_candidate = []
        self._coannihilation = []
        self._z2_odd = []
        self._disconnected_particles = []
        self._param_card = None
        # self.coannihilation_diff = 1
        self.option_2to2lo = False # keep track if the user wants to use this option
        self._relic_off = True
        self._has_indirect = False
        self._has_spectral = False
        self.relic_generated_by_id = False # keep track if relic has been generated by id
        # ID procs, matrix_elements, amps --> lists with 2 elements [tree, loop-induced]
        self._ID_cont_procs           = [base_objects.ProcessDefinitionList() for i in range(2)]
        self._ID_line_procs           = [base_objects.ProcessDefinitionList() for i in range(2)]
        self._ID_cont_matrix_elements = [helas_objects.HelasMultiProcess() for i in range(2)]
        self._ID_line_matrix_elements = [helas_objects.HelasMultiProcess() for i in range(2)]
        self._ID_cont_amps            = [diagram_generation.AmplitudeList() for i in range(2)]
        self._ID_line_amps            = [diagram_generation.AmplitudeList() for i in range(2)]
        # only for display purposes
        self._last_amps = diagram_generation.AmplitudeList()
        #self.dm = darkmatter.darkmatter()
        self._forbid_fast = False
        self._unphysical_particles = []
        self._pdg_particle_map = maddm_run_interface.PDGParticleMap()


################################################################################        
# DISPLAY COMMAND
################################################################################        
    def do_display(self, *args, **opts):
        ''' concerning 'processes', 'diagrams' and 'diagrams_text' it accepts options:
                relic, direct, indirect: displays the generated processes/diagrams of the chosen category;
                all: displays all the processes/diagrams that have been generated;
                last: displays the processes/diagrams generated in the last command.
        '''
        args = list(args) # convert the tuple to list so we can modify it
        arguments = self.split_arg(args[0]) #args[0] is the line of the command

        if all(arg in ['darkmatter', 'coannihilators', 'z2-odd', 'z2-even', 'disconnected'] for arg in arguments):
            to_display = {
                'darkmatter': self._dm_candidate,
                'coannihilators': self._coannihilation,
                'z2-odd': self._z2_odd,
                'z2-even': [p for p in self._curr_model.get('particles') if abs(p.get('pdg_code')) not in [a.get('pdg_code') for a in self._z2_odd+self._disconnected_particles]],
                'disconnected': self._disconnected_particles
            }
            which = [arg for arg in arguments if arg in list(to_display.keys())]
            which = list(set(which))
            for arg in which:
                logger.info(arg + ": " + ", ".join([p.get_name() for p in to_display[arg]]))
            return
        
        if arguments[0] in ['processes', 'diagrams', 'diagrams_text']:
            # check the generated diagrams
            self._relic_amps    = [amp for amp in self._curr_amps if amp.get('process').get('id') in [self.process_tag['DM2SM'], self.process_tag['DMSM']]]
            self._direct_amps   = [amp for amp in self._curr_amps if amp.get('process').get('id') == self.process_tag['DD']]
            self._indirect_amps = self._ID_cont_amps[0] + self._ID_line_amps[0] + self._ID_cont_amps[1] + self._ID_line_amps[1]
            generated_amps = {
                'relic'   : bool(self._relic_amps),
                'direct'  : bool(self._direct_amps),
                'indirect': bool(self._indirect_amps)
            }
            # in the following lines we filter out only the allowed options
            which = [arg for arg in arguments if arg in list(generated_amps.keys()) + ['all', 'last']]
            which = list(set(which))
            for arg in which:
                arguments.remove(arg)
            if len(which) == 0:
                which.append('all')
            if 'all' in which:
                generated_categories = [category for category, generated in generated_amps.items() if generated]
                if len(generated_categories) == 0:
                    logger.error("No %s have been generated yet." % (arguments[0].split('_')[0]))
                    return
                args[0] = " ".join(arguments) + ' ' + ' '.join(generated_categories)
                self.do_display(*args, **opts)
                return
            args[0] = " ".join(arguments)
            N_head = 80
            if 'last' in which:
                header = "==== LAST %s " % arguments[0].upper()
                print(header + "=" * max([N_head - len(header), 0]))
                with misc.TMP_variable(self, ['_curr_amps'], [self._last_amps]):
                    try:
                        super(MadDM_interface, self).do_display(*args, **opts)
                    except AttributeError:
                        logger.error("No processes were generated in the last call.")
                return
            else:
                for category in which:
                    if (category == 'relic') and self._relic_off:
                        continue
                    header = "==== %s %s " % (category.upper(), arguments[0].upper())
                    print(header + "=" * max([N_head - len(header), 0]))
                    with misc.TMP_variable(self, ['_curr_amps'], [getattr(self, "_%s_amps" % category)]):
                        super(MadDM_interface, self).do_display(*args, **opts)
                return
        super(MadDM_interface, self).do_display(*args, **opts)

    def complete_display(self, text, line, begidx, endidx, formatting = True):
        ''' if len(args) > 1 it means we have already written the subcommand, so when we run the superclass completion
            it will not return anything, setting out = None, so that the next method deal_multiple_categories
            would crash, because it tries to get the len of a None object.
            We need to set explicitly out = [] when it is None.
        '''
        args = self.split_arg(line[:begidx])
        
        out = super(MadDM_interface, self).complete_display(text, line, begidx, endidx)

        if out is None:
            out = []
        if not isinstance(out, dict):
            out = {"standard options": out}

        if len(args) > 0 and all(arg in ['darkmatter', 'coannihilators', 'z2-odd', 'z2-even', 'disconnected'] for arg in args[1:]):
            options = ['darkmatter', 'coannihilators', 'z2-odd', 'z2-even', 'disconnected']
            out['maddm options'] = self.list_completion(text, options , line)

        if len(args) > 1 and args[1] in ['processes', 'diagrams', 'diagrams_text']:
            options = ['relic', 'direct', 'indirect', 'all', 'last']
            out['processes|diagrams|diagrams_text options'] = self.list_completion(text, options, line)

        return self.deal_multiple_categories(out, formatting)

    def help_display(self):
        logger.info("**************** MADDM NEW OPTION ***************************")
        logger.info("syntax: display processes|diagrams|diagrams_text [PATH] [OPTIONS]", '$MG:color:BLUE')
        logger.info(" -- displays the generated processes/diagrams, storing the latter")
        logger.info("    in the specified PATH")   
        logger.info(" OPTIONS:", '$MG:BOLD')
        logger.info("    - You can display only certain categories of processes/diagrams,")
        logger.info("      by using 'relic', 'direct', 'indirect' or a combination of them")
        logger.info("     o example: display processes relic indirect",'$MG:color:GREEN')   
        logger.info("    - You can display all of the generated processes/diagrams,")
        logger.info("      by using 'all'. This overwrites the other options.")
        logger.info("      This is the default if no options are given.")
        logger.info("     o example: display processes all",'$MG:color:GREEN')             
        logger.info("    - You can display the processes/diagrams generated in the last call,")
        logger.info("      by using 'last'. This overwrites the other options, except 'all'.")
        logger.info("     o example: display processes last",'$MG:color:GREEN')         
        logger.info("")
        logger.info("syntax: display darkmatter|coannihilators|z2-odd|z2-even|disconnected", '$MG:color:BLUE')
        logger.info(" -- displays the related list of particles") 
        logger.info("**************** MG5AMC OPTION ***************************")
        super(MadDM_interface, self).help_display()

################################################################################        
# J-FACTOR COMMAND
################################################################################        
    def do_jfactor(self, *args, **opts):
        ''' compute the J-factor, leveraging the script jfactor_gc.py
        '''
        if not (HAS_NUMPY and HAS_SCIPY):
            logger.error("numpy and/or scipy are required to run this command")
        parser = jfactor_gc.command_parser()
        arguments = self.split_arg(args[0]) #args[0] is the line of the command
        jfactor_value, density_profile = jfactor_gc.parse_cmd_line(parser=parser, cmd_line=arguments)
        logger.info("J-factor = %.4e GeV^2 cm^{-5}" % jfactor_value)
        logger.info("Density profile: " + str(density_profile))

    def complete_jfactor(self, text, line, begidx, endidx, formatting = True):
        # args = self.split_arg(line[:begidx])
        # 
        # out = super(MadDM_interface, self).complete_display(text, line, begidx, endidx)
        #
        # if out is None:
        #     out = []
        # if not isinstance(out, dict):
        #     out = {"standard options": out}
        #
        # if len(args) > 0 and all(arg in ['darkmatter', 'coannihilators', 'z2-odd', 'z2-even', 'disconnected'] for arg in args[1:]):
        #     options = ['darkmatter', 'coannihilators', 'z2-odd', 'z2-even', 'disconnected']
        #     out['maddm options'] = self.list_completion(text, options , line)
        #
        # if len(args) > 1 and args[1] in ['processes', 'diagrams', 'diagrams_text']:
        #     options = ['relic', 'direct', 'indirect', 'all', 'last']
        #     out['processes|diagrams|diagrams_text options'] = self.list_completion(text, options, line)
        #
        # return self.deal_multiple_categories(out, formatting)
        pass

    def help_jfactor(self):
        parser = jfactor_gc.command_parser()
        logger.info(parser.format_help().replace("usage: %s" % parser.prog, bcolors.BOLD + "Compute the J-Factor for the Galactic Center\n" + bcolors.ENDC + "             jfactor"))

################################################################################        
# DEFINE COMMAND
################################################################################        
    def do_define(self, line, **opts):
        """pass"""
        args = self.split_arg(line)
        if len(args) and args[0] in MadDM_interface._define_options:
            if args[0] == 'darkmatter':
                if len(args) == 1:
                    self.search_dm_candidate([])
                elif args[1].startswith('/'):
                    self.search_dm_candidate([a.replace('/', '') for a in args[1:]])
                elif len(args)>=2:
                    self._dm_candidate = [self._curr_model.get_particle(dm) for dm in args[1:]]
                    #remove None
                    self._dm_candidate = [dm for dm in self._dm_candidate if dm]
                    if not self._dm_candidate:
                        raise DMError('%s is not a valid particle for the model.' % args[1]) 
                    if len(self._dm_candidate) == 1:
                        # No update of the model if 2(or more) DM since DD is not possible
                        self._z2_odd.extend([p for p in self._dm_candidate])
                        self.find_z2_odd_particles()
                        self.update_model_with_EFT()
                    else:
                        txt = self._curr_model.write_param_card()
                        ff = open('/tmp/param_card.dat','w')
                        ff.write(txt)
                        ff.close()
                        self.define_benchmark(answer=True, path='/tmp/param_card.dat')
                        
            elif args[0] == 'coannihilator':
                if not self._dm_candidate:
                    self.search_dm_candidate([])
                for i in range(1,len(args)): # convert pdgs to int
                    try:
                        args[i] = int(args[i])
                    except:
                        pass
                
                if len(args) == 1:
                    self.search_coannihilator()                    
                else:
                    z2_odd_allowed = [p for p in self._z2_odd if p.get_pdg_code() not in [dm.get_pdg_code() for dm in self._dm_candidate]]
                    z2_odd_allowed_names = []
                    for p in z2_odd_allowed:
                        # user may in principle use name or antiname of the particle, or even the pdg (or -pdg for antiparticles)
                        z2_odd_allowed_names.extend([p.get('name'), p.get('antiname'), p.get_pdg_code(), p.get_anti_pdg_code()])
                    for a in args[1:]:
                        if a == '--usepdg':
                            continue
                        elif a not in z2_odd_allowed_names:
                            logger.error("'%s' is not allowed (either does not exist or is not z2-odd or is a dark matter candidate). Particles allowed: " % a + ", ".join([p.get_name() for p in z2_odd_allowed]))
                            continue
                        self._coannihilation.append(self._curr_model.get_particle(a))
                            
                # elif len(args)>2 and isdigit and args[2].startswith('/'):
                #     self.search_coannihilator(gap=args[1], excluded=[a.replace('/', '') for a in args[2:]])
                # elif len(args)>1 and not isdigit and args[1].startswith('/'):
                #     self.search_coannihilator( excluded=[a.replace('/', '') for a in args[1:]])
                # elif len(args)==2 and isdigit:
                #     self.search_coannihilator(gap=args[1])
                # else:
                #     self._coannihilation = [self._curr_model.get_particle(a) for a in args[1:]]

                # avoid duplication
                if None in self._coannihilation:
                    raise self.InvalidCmd('Some of the particle name are invalid. Please retry.')
                all_name = [c.get('name') for c in self._coannihilation]
                self._coannihilation = [c for i,c in enumerate(self._coannihilation) if c.get('name') not in all_name[:i]]
                
            elif args[0] == 'benchmark':
                if len(args)==1:
                    self.define_benchmark()
                elif os.path.exists(args[1]):
                    self.define_benchmark(answer=args[2])
        else:
            return super(MadDM_interface, self).do_define(line,**opts)

    def complete_define(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddm options"""

        out = {}
        args = self.split_arg(line[0:begidx])
        if len(args) == 1:
            out['maddm options'] = self.list_completion(text, self._define_options, line)
        out['particles name']  = self.model_completion(text, line[6:], line)
         
        return self.deal_multiple_categories(out, formatting)

    def help_define(self):
        """ """
        logger.info("**************** MADDM NEW OPTION ***************************")
        logger.info("syntax: define benchmark [PATH]", '$MG:color:BLUE')
        logger.info(" -- define the current benchmark. if no path are specified, a question is asked")
        logger.info("")
        logger.info("syntax: define darkmatter [OPTIONS]", '$MG:color:BLUE')
        logger.info(" -- define the current (set) of darkmatter.")
        logger.info("    If no option specified. It assigned the less massive neutral BSM particle (with zero width).")
        logger.info(" OPTIONS:", '$MG:BOLD')
        logger.info("    - You can specify the darkmatter by specifying their name/pdg code")
        logger.info("     o example: define darkmatter n1 n2",'$MG:color:GREEN')
        logger.info("    - You can remove some particle from the search by prefixing the particle name by '/'.")
        logger.info("     o example: define darkmatter /n1", '$MG:color:GREEN')                    
        logger.info("")        
        logger.info("syntax: define coannihilator [REL_DIFF | PARTICLE(S)] [/ excluded_particles]", '$MG:color:BLUE')
        logger.info(" -- define the current (sets) of coannihilator. ")
        logger.info("    If no option is provided use a  relative difference of 10% is used.")         
        logger.info(" OPTIONS:", '$MG:BOLD')
        logger.info("    - You can specify the coannihilator by specifying their name")
        logger.info("     o example: define coannihilator n2 n3",'$MG:color:GREEN')   
        logger.info("    - You can specify the coannihilator by specifying their name if you use --usepdg")
        logger.info("     o example: define coannihilator 20000011 --usepdg",'$MG:color:GREEN')             
        logger.info("    - When specifying an allowed relative difference, you can exclude some particles")
        logger.info("     o example: define coannihilator 0.05 / n2 n3",'$MG:color:GREEN')         
        logger.info("**************** MG5AMC OPTION ***************************")
        super(MadDM_interface, self).help_define()


    def search_dm_candidate(self, excluded_particles=[]):
        """ automatic search of the dm candidate"""

        if not self._param_card:
            self.define_benchmark()
 
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
        #  5. We then check that such DM candidate has:
        #          - some non zero couplings
        #          - no decay to pair of SM particle                            #
        #                                                                       #
        #-----------------------------------------------------------------------#
        particles = self._curr_model.get('particles')  
        bsm_particles = [p for p in particles 
                         if 25 < p.get('pdg_code') < 99000000 and\
                         (p.get('name') not in excluded_particles and p.get('antiname') not in excluded_particles) and\
                         p.get('charge') == 0]
        
        dm_particles = []
        dm_mass = 0
        # When manually entering in the DM and coannihilation particles, both particle and anti-particle
        # are not both required to be listed. If one is not included then it is added later.
        get_width = lambda p: self._curr_model.get_width(p['pdg_code'])

        # Looping over all the BSM particles we check the criteria described above
        for p in bsm_particles:
            curr_mass = abs(self._curr_model.get_mass(p['pdg_code']))
            if (p['width'] == 'ZERO' or  get_width(p) == 0) and curr_mass>0:
                # If nothing has been found so far then set the first DM candidate
                # If we already found a candidate, compare the masses and keep the one that's lighter
                if not dm_particles or curr_mass < dm_mass:
                    dm_particles = [p]
                    dm_mass = curr_mass
                elif dm_particles and curr_mass == dm_mass:
                    dm_particles.append(p)

        # Check to see if we actually found a DM candidate
        if not dm_particles:
            raise DMError("No dark matter candidates in the model!")
        if dm_mass == 0:
            raise DMError("No dark matter candidate since we found a DarkEnergy candidate: %s" % dm_particles[0]['name'])

        # Validity check
        dm_names = [p.get('name') for p in dm_particles]
        for p in list(dm_particles):
            has_coupling = False
            for vert in self._curr_model.get('interactions'):
                if p in vert['particles']:
                    for coup_name in  vert['couplings'].values():
                        if self._curr_model.get('coupling_dict')[coup_name]:
                            has_coupling = True
                            break
                if has_coupling:
                    break
            else:
                dm_particles.remove(p)
                dm_names.remove(p.get('name'))

        if not dm_particles:
            logger.warning("""found %s but none have them have non zero coupling. Retry by excluding those""", ','.join(dm_names))
            return self.search_dm_candidate(excluded_particles=excluded_particles+dm_names)
        
        if len(dm_particles) >1:
            choice = self.ask("More than one valid candidate found: Please select the one that you want:",
                     default=dm_names[0], choices=dm_names)
            dm_particles = [p for p in dm_particles if p.get('name') == choice]
            assert(dm_particles)             
        # Print out the DM candidate
        logger.info("Found Dark Matter candidate: %s" % dm_particles[0]['name'],  '$MG:BOLD')
        self._dm_candidate = dm_particles
        self._z2_odd.extend([p for p in dm_particles])
        self.find_z2_odd_particles()
        self.update_model_with_EFT()
        
    def define_benchmark(self, path=None, answer=False):
     
     
        self._dm_candidate = []
        self._coannihilation = []
        self._disconnected_particles = []
        question = """Do you want to edit the benchmark (press enter to bypass editing)?\n"""
        question += """ - Press 1 or param to open an editor and edit the file\n"""
        question += """ - You can specify a path to a valid card (and potentially edit it afterwards)\n"""
        question += """ - You can use the 'set' command to edit the card\n """
        possible_answer = ['0', 'done','1', 'param']
#        card = {0:'done', 1:'param'}
        
        if not answer:
            if not path:
                    dirpath = self._curr_model.get('modelpath')
                    path = pjoin(dirpath, 'param_card_orig.dat')
                    if self._param_card:
                        self._param_card.write(path)
                    else:
                        self._curr_model.write_param_card(path)     
            out=''
            while out not in ['0', 'done']:
                timeout = 0
                out = self.ask(question, '0', possible_answer, timeout=int(1.5*timeout),
                                  path_msg='enter path', 
                                  ask_class = EditParamCard,
                                  card=[path],
                                  pwd=os.getcwd(),
                                  param_consistency=False
                                  )
        elif isinstance(answer, str):
            path = answer
            
        if not isinstance(self._curr_model, model_reader.ModelReader):
            self._curr_model = model_reader.ModelReader(self._curr_model)        
        self._curr_model.set_parameters_and_couplings(path)   
        self._param_card = check_param_card.ParamCard(path)
        #self._dm_candidate = darkmatter.darkmatter()
    
      
    # def search_coannihilator(self, gap=0.1, excluded=[]):
    #     """  This routine finds the coannihilation particles for the relic
    #     #  density calculation. The code can search for all the BSM particles that are
    #     #  within an input mass difference with the DM candidate.  All
    #     #  coannihilation particles are then added to the self._dm_particles
    #     #  list.
    #     """

    #     self._coannihilation = []
    #     dm_mass = self._curr_model.get_mass(self._dm_candidate[0])
    #     self.coannihilation_diff = gap 
    #     dm_name = [dm['name'] for dm in self._dm_candidate]
        
    #     bsm_particles = [p for p in self._curr_model.get('particles')                
    #                      if 25 < p.get('pdg_code') < 999000000 and\
    #                       (p.get('name') not in excluded or 
    #                       p.get('antiname') not in excluded or
    #                       str(p.get('pdgcode')) not in excluded)] 

    #     # Loop over BSM particles
    #     for p in bsm_particles:
    #         if p['name'] in dm_name:
    #             continue
            
    #         bsm_mass = self._curr_model.get_mass(p)
    #         if (abs(dm_mass-bsm_mass)/dm_mass <= gap):
    #             self._coannihilation.append(p)
    #         # If there are BSM particles that are too small to be included in the coannihilation they
    #         # are still tabulated to include in the final state particles with the SM particles.
    #         #elif ((self._bsm_masses[i] < (1.0-self._coann_eps)*self._dm_mass) and \
    #         #    not ('all' in self._excluded_particles) and \
    #         #    not (self._particles[self._bsm_particles[i]] in self._excluded_particles)):
    #         #          self._bsm_final_states.append(self._particles[self._bsm_particles[i]])

    #         # For organizational purposes we put the coannihilation particles by
    #         #alphabetical order by name
    #         self._coannihilation.sort(key=lambda p: p['name'])
        
    #     # Validity check
    #     for p in list(self._coannihilation):
    #         has_coupling = False
    #         for vert in self._curr_model.get('interactions'):
    #             if p in vert['particles']:
    #                 for coup_name in  vert['couplings'].values():
    #                     if self._curr_model.get('coupling_dict')[coup_name]:
    #                         has_coupling = True
    #                         break
    #             if has_coupling:
    #                 break
    #         else:
    #             self._coannihilation.remove(p)
        
        
    #     if self._coannihilation:
    #         logger.info("Found coannihilation partners: %s" % ','.join([p['name'] for p in self._coannihilation]),
    #                 '$MG:BOLD')
    #     else:
    #         logger.info("No coannihilation partners found.", '$MG:BOLD')

    def search_coannihilator(self,excluded=None):
        ''' automatically make the coannihilators list equal to the list of Z2-odd particles found, exclude dark matter candidates. '''
        if excluded is None:
            excluded = []
        self._coannihilation[:] = [p for p in self._z2_odd if 
                           (p.get_pdg_code() not in [dm.get_pdg_code() for dm in self._dm_candidate] and 
                           p.get_pdg_code() not in excluded)]
        logger.info("Found coannihilator(s): " + ", ".join([p.get_name() for p in self._coannihilation]))

    def find_z2_odd_particles(self):
        """ Find Z2-odd particles and disconnected particles. """
        # Z2-odd particles: self._dm_candidate + self._coannihilation
        # unknown Z2 parity particles: self._disconnected
        Z2_EVEN = 1
        Z2_ODD = -1
        Z2_UNKN = 0
        # reduction functions
        def remove_any_number_even_particles(z2_parity, particles):
            ''' gets rid of the even particles '''
            to_remove = [i for i, z2 in enumerate(z2_parity) if z2 == Z2_EVEN]
            for i in sorted(to_remove, reverse = True): # reverse indexes upon deletion so smaller indexes won't change
                del z2_parity[i]
                del particles[i]
        def remove_even_number_same_particles(z2_parity, particles):
            ''' gets rid of an even number of particles with the same name '''
            to_remove = []
            unique_particles = set([p.get('name') for p in particles])
            for p_name in unique_particles:
                indexes = [i for i, p in enumerate(particles) if p_name == p.get('name')]
                if len(indexes) < 2:
                    continue
                # now I have a list of indexes labelling particles with the same name
                # keep only an even number of them
                trunc_len = 2 * (len(indexes) // 2)
                to_remove.extend(indexes[0:trunc_len])
            # now I have a list of indexes of the particles to remove
            for i in sorted(to_remove, reverse = True): # reverse indexes upon deletion so smaller indexes won't change
                del z2_parity[i]
                del particles[i]
        def remove_even_number_odd_particles(z2_parity, particles):
            ''' gets rid of an even number of Z2-odd particles '''
            indexes = [i for i, z2 in enumerate(z2_parity) if z2 == Z2_ODD]
            # now I have a list of indexes labelling Z2-odd particles
            # keep only an even number of them
            trunc_len = 2 * (len(indexes) // 2)
            to_remove = indexes[0:trunc_len]
            for i in sorted(to_remove, reverse = True): # reverse indexes upon deletion so smaller indexes won't change
                del z2_parity[i]
                del particles[i]
        # all bsm particles except dark matter have an unknown parity
        z2_odd_particles = [p.get('name') for p in self._z2_odd]
        for p in self._curr_model.get('particles'):
            if p.get('pdg_code') <= 25 or p.get('name') in z2_odd_particles or p.get('pdg_code') in [9000001, 9000002, 9000003, 9000004, 82, 250, 251]: # remove fixed PDG code ghosts
                continue
            self._disconnected_particles.append(p)
        z2_unkn_particles = [p.get('name') for p in self._disconnected_particles]
        def assign_parity(part):
            if part.get('name') in z2_odd_particles:
                return Z2_ODD
            elif part.get('name') in z2_unkn_particles:
                return Z2_UNKN
            else:
                return Z2_EVEN
        # store particle list for each vertex
        vert_particles = [[p for p in vertex['particles']] for vertex in self._curr_model['interactions'].get_type('base')]
        # core of the algorithm: is_modified is a flag which is True when something changes in the particles' list
        is_modified = True
        # ["[%s] -> [%s]" % (",".join([p.get('name') for p in vpart]),",".join(list(map(str,[assign_parity(p) for p in vpart])))) for vpart in vert_particles]
        while(is_modified):
            is_modified = False # if nothing is modified, then the loop will finish
            to_remove = [] # contains indexes of the vertices to remove after the next loop
            for i_vert, this_vertex_particles in enumerate(vert_particles):
                vert_z2 = [assign_parity(p) for p in this_vertex_particles]
                if vert_z2.count(Z2_UNKN) == 0:
                    to_remove.append(i_vert)
                    continue
                remove_any_number_even_particles(vert_z2, this_vertex_particles)
                remove_even_number_same_particles(vert_z2, this_vertex_particles)
                remove_even_number_odd_particles(vert_z2, this_vertex_particles)
                count_unkn = vert_z2.count(Z2_UNKN)
                # only 3 cases are possible at this stage
                if count_unkn == 0: # vert_z2 == []
                    to_remove.append(i_vert)
                elif count_unkn == 1: # vert_z2 == [0] or vert_z2 == [0, -1] or vert_z2 == [-1, 0]
                    is_modified = True
                    if Z2_ODD in vert_z2: # in this case the unknown particle is Z2-odd
                        unkn = this_vertex_particles[vert_z2.index(Z2_UNKN)]
                        self._z2_odd.append(unkn)
                        # update the name list
                        z2_odd_particles.append(unkn.get('name'))
                    else: # in this case the unknown particle is Z2-even
                        unkn = this_vertex_particles[0]
                    # find and remove that particle from the self._disconnected_particles list
                    i_unkn_name = z2_unkn_particles.index(unkn.get('name')) # find the index in the name list
                    del self._disconnected_particles[i_unkn_name] # remove the particle from the list
                    # check: if self._disconnected_particles is empty then we are done: each bsm particle has been catalogued
                    if len(self._disconnected_particles) == 0:
                        is_modified = False # set this to False to break also the outer loop
                        break
                    # update the name list
                    del z2_unkn_particles[i_unkn_name]
                    # add the vertex index to the list of indexes of the vertices to remove at the end of this loop
                    to_remove.append(i_vert)
            for i in sorted(to_remove, reverse = True):
                del vert_particles[i]
        if len(self._z2_odd) != 0:
            self.exec_cmd('display z2-odd', postcmd = False)

    def do_import(self, line,*args, **opts):
        """normal import but perform a cleaning for MadDM  variables"""
        
        lineargs = self.split_arg(line)
        self.check_import(lineargs)
        if lineargs and lineargs[0].startswith('model'):
            #reset DM information
            self._param_card = None
            self._dm_candidate = []
            self._coannihilation = []
            self._disconnected_particles = []
            
        return super(MadDM_interface, self).do_import(line, *args, **opts)

    def clean_process(self):
        """ensure that all processes are cleaned from memory.
        typically called from import model and generate XXX command
        """
        super(MadDM_interface, self).clean_process()
        self._ID_cont_procs           = [base_objects.ProcessDefinitionList() for i in range(2)]
        self._ID_line_procs           = [base_objects.ProcessDefinitionList() for i in range(2)]
        self._ID_cont_matrix_elements = [helas_objects.HelasMultiProcess() for i in range(2)]
        self._ID_line_matrix_elements = [helas_objects.HelasMultiProcess() for i in range(2)]
        self._ID_cont_amps            = [diagram_generation.AmplitudeList() for i in range(2)]
        self._ID_line_amps            = [diagram_generation.AmplitudeList() for i in range(2)]
        self._last_amps               = diagram_generation.AmplitudeList()
        self.option_2to2lo = False
        self._has_indirect = False
        self._has_spectral = False
        self._forbid_fast  = False
        self.relic_generated_by_id = False
        
    def check_save(self, args):
        
        if args[1] == 'maddm_first_indirect':
            args.insert(1, pjoin(MG5DIR, 'input', 'mg5_configuration.txt'))
            print(args)
            return args
        else:
            return super(MadDM_interface, self).check_save(args)
        
    
    def help_generate(self):
        """ """
        logger.info("**************** MADDM NEW OPTION ***************************")
        logger.info("syntax: generate|add relic_density [/ X]", '$MG:color:BLUE')
        logger.info(" -- generate relic density matrix element excluding any particle(s) X") 
        logger.info("    to appear as s/t channel in any diagram")
        logger.info("    - FOR ADVANCED USER", '$MG:BOLD')
        logger.info("      manual definition of the matrix element is possible via the following syntax:")
        logger.info("        add process DM  DM > Z Y @ DM2SM # annihilation" )
        logger.info("        add process DMi  DMj > DMk DMl @ DM2DM # DM-diffusion for co-annihilation")
        logger.info("")
        logger.info("syntax: generate|add direct_detection [/ X]", '$MG:color:BLUE')
        logger.info(" -- generate direct detection matrix element excluding any particle(s) X") 
        logger.info("    to appear as s/t channel in any diagram")        
        logger.info("    - FOR ADVANCED USER", '$MG:BOLD')
        logger.info("      manual definition of the matrix element is possible via the following syntax:")
        logger.info("        add process DM  N > DM N @ DD")
        logger.info("")
        logger.info("syntax: generate|add indirect_detection [2to2lo] [final_states] [/ X]", '$MG:color:BLUE')    
        logger.info(" -- generate indirect detection matrix element")
        logger.info("       X forbids any s/t channel propagator to be present in the Feynman diagram")
        logger.info("    syntax: generate|add indirect_detection 2to2lo / X", '$MG:color:BLUE') 
        logger.info("    -- generate indirect detection matrix element for inclusive method of integration only")
        logger.info("    syntax: generate|add indirect_detection F1 F2 / X", '$MG:color:BLUE') 
        logger.info("    -- generate indirect detection matrix element for a given final state (F1 F2)")
        logger.info("       'three body final states (actually n body finally states) are allowed (those are forbidden when using PPPC4DMID)")
        logger.info("")    
        logger.info("    - FOR ADVANCED USER", '$MG:BOLD')             
        logger.info("      manual definition of the matrix element is possible via the following syntax:")
        logger.info("        note: only for 2to2lo")
        logger.info("        add process DM  N > DM N @ ID")
        logger.info("")
        logger.info("syntax: generate|add indirect_spectral_features [final_states]", '$MG:color:BLUE')    
        logger.info(" -- generate indirect detection matrix elements for final states 'a X'")
        logger.info("    syntax: generate|add indirect_spectral_features F1 F2", '$MG:color:BLUE') 
        logger.info("    -- generate indirect detection matrix element for a given final state (F1 F2) among the allowed ones")
        logger.info("       note: at least one of F1 or F2 must be a photon ('a' or '22')")
        logger.info("")
        logger.info("**************** MG5AMC OPTION ***************************")
        super(MadDM_interface, self).help_generate()

    help_add = help_generate
    
    def do_add(self, line, user=True):
        """ """
        args = self.split_arg(line)    
        if len(args) and args[0] == 'process':
            args.pop(0)
        if len(args) and args[0] in["relic_density", 'relic']:
            # check that relic is empty
            if any(p.get('id') != self.process_tag['DM2SM'] for p in self._curr_proc_defs):
                logger.warning('relic density Feynman diagram already generated (likely due to indirect mode). Avoid to doing it again.')
                return
            if '/' not in line:
                return self.generate_relic([])
            else:
                subline = line.split('/',1)[1]
                excluded = [ a for a in self.split_arg(subline) if not a.startswith('-')]
                return self.generate_relic(excluded)
        elif len(args) and args[0] in ["direct_detection", "direct"]:
            if '/' not in line:
                return self.generate_direct([])
            else:
                subline = line.split('/',1)[1]
                excluded = [ a for a in self.split_arg(subline) if not a.startswith('-')]
                return self.generate_direct(excluded)
        elif len(args) and args[0] in ["indirect_detection", "indirect"]:
            self.generate_indirect(args[1:])
        elif len(args) and args[0] in ["indirect_spectral_features", "spectral"]:
            self.generate_spectral(args[1:])

        else:
            if '@' in line: # this allows to substitute the tag words with the numbers allowed in MadDM
                try:
                    line = re.sub(r'''(?<=@)(%s\b)''' % '\\b|'.join(self.process_tag), 
                                lambda x: repr(self.process_tag[x.group(0)]), line)
                    tag = int(re.search(r'(?<=@)\d+', line).group(0))
                    maddm_can_handle = tag in list(self.process_tag.values())
                except:
                    tag = None
                    maddm_can_handle = False

                if maddm_can_handle and (not self._dm_candidate):
                    self.search_dm_candidate()

                if (tag == self.process_tag['ID']) and user:
                    # do this only if the user has specified ID
                    # in this case we need to handle the various _ID_*_* variables on the basis
                    # of the kind of process the user has input
                    # if we get here from the generate_indirect/spectral commands, then user == False, so we would jump these instructions
                    # if we get here from user-defined relic density or direct detection processes, then '@1995' not in line, so we would jump these instructions
                    
                    if len(self._curr_proc_defs) == 0 or all(p.get('id') != self.process_tag['DM2SM'] for p in self._curr_proc_defs):
                        logger.info('For indirect detection we need to generate relic density matrix-element','$MG:BOLD')
                        self.generate_relic([], user = False)
                        self.relic_generated_by_id = True
                    
                    final_states = re.search(r'(?<=>).*?(?=[@\[])', line).group(0)
                    if not final_states:
                        logger.error("Invalid process.")
                        return
                    line_or_cont = 'line' if (self.indirect_contains_photon(final_states) and ',' not in line) else 'cont' # spectral processes can't contain decay chains
                    id_proc_type_index = 1 if ('noborn' in line) else 0
                    ID_procs           = getattr(self, '_ID_%s_procs' % line_or_cont)[id_proc_type_index]
                    ID_matrix_elements = getattr(self, '_ID_%s_matrix_elements' % line_or_cont)[id_proc_type_index]
                    ID_amps            = getattr(self, '_ID_%s_amps' % line_or_cont)[id_proc_type_index]
                    try:
                        with misc.TMP_variable(self, ['_curr_proc_defs', '_curr_matrix_elements', '_curr_amps'], 
                                            [ID_procs, ID_matrix_elements, ID_amps]):
                            n_actual_amps = len(self._curr_amps)
                            function_return = super(MadDM_interface, self).do_add(line)
                            # if we get here, then the generation has gone well
                            # initialize _last_amps
                            self._last_amps = self._curr_amps[n_actual_amps:]
                            # switch on indirect/spectral
                            self._has_indirect = (self._has_indirect or line_or_cont == 'cont')
                            self._has_spectral = (self._has_spectral or line_or_cont == 'line')
                            return function_return
                    except (self.InvalidCmd, diagram_generation.NoDiagramException) as error:
                        logger.error(error)
                        self.change_principal_cmd('MadGraph')
                        return
                elif maddm_can_handle and user:
                    try:
                        n_actual_amps = len(self._curr_amps)
                        function_return = super(MadDM_interface, self).do_add(line)
                        self._last_amps = self._curr_amps[n_actual_amps:]
                        # switch on relic density in case the input process is related to relic density
                        self._relic_off = (self._relic_off and (tag not in [self.process_tag['DM2SM'], self.process_tag['DMSM']]))
                        return function_return
                    except (self.InvalidCmd, diagram_generation.NoDiagramException) as error:
                        logger.error(error)
                        return
                elif user:
                    logger.error("Process not supported by MadDM, allowed tags are: " + ", ".join(list(self.process_tag.keys())))
                    return
            
            return super(MadDM_interface, self).do_add(line)

    @misc.mute_logger(['madgraph', 'models'], [30,30])    
    def add_model(self,*args,**opts):
        super(MadDM_interface, self).add_model(*args, **opts)
        self.model_nlo = True
        self._pdg_particle_map.set_model_map(self._curr_model)
        possible_added_particles = ['999000006', '999000007', '999000008', '999000009', '999000010'] # from __REAL or __COMPLEX model
        self._unphysical_particles = []
        for p in self._curr_model.get('particles'):
            if p.get('type') != '' or str(p.get('pdg_code')) in possible_added_particles:
                if not p.get('self_antipart'):
                    self._unphysical_particles += [str(p.get('pdg_code')), str(-p.get('pdg_code'))]
                else:
                    self._unphysical_particles += [str(p.get('pdg_code'))]

    def complete_generate(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddm options"""

        args = self.split_arg(line[:begidx])
        out = super(MadDM_interface, self).complete_generate(text, line, begidx, endidx,formatting=False)
        if not isinstance(out, dict):
            out = {"standard options": out}

        if len(args) == 1:
            options = ['relic_density', 'direct_detection', 'indirect_detection', 'indirect_spectral_features']#, 'capture']
            out['maddm options'] = self.list_completion(text, options , line)
        return self.deal_multiple_categories(out, formatting)

    def complete_add(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddm options"""

        args = self.split_arg(line[:begidx])
        out = super(MadDM_interface, self).complete_add(text, line, begidx, endidx,formatting=False)
        if not isinstance(out, dict):
            out = {"standard options": out}
        
        if len(args) == 1:
            options = ['relic_density', 'direct_detection', 'indirect_detection', 'indirect_spectral_features']#, 'capture']
            out['maddm options'] = self.list_completion(text, options , line)
        return self.deal_multiple_categories(out, formatting)


    def do_tutorial(self, line):
        """Activate/deactivate the tutorial mode."""
        
        if line:
            super(MadDM_interface, self).do_tutorial(line)
            
        args = self.split_arg(line)
        #self.check_tutorial(args)
        if not args:
            args = ['maddm']
        tutorials = {'maddm': logger_tuto}
        
        try:
            tutorials[args[0]].setLevel(logging.INFO)
        except KeyError:
            logger_tuto.info("\n\tThanks for using the tutorial!")
            logger_tuto.setLevel(logging.ERROR)
            
    def postcmd(self,stop, line):
        """ finishing a command
        This looks if the command add a special post part.
        This looks if we have to write an additional text for the tutorial."""

        stop = super(MadDM_interface, self).postcmd(stop, line)
        # Print additional information in case of routines fails
        if stop == False:
            return False

        args=line.split()
        # Return for empty line
        if len(args)==0:
            return stop
        
        if args[0] == 'import' and ('__REAL' in line or '__COMPLEX' in line):
            return stop 

        # try to print linked to the first word in command
        #as import_model,... if you don't find then try print with only
        #the first word.
        if len(args)==1:
            command=args[0]
        else:
            command = args[0]+'_'+args[1].split('.')[0]

        try:
            logger_tuto.info(getattr(dm_tutorial_text, command).replace('\n','\n\t'))
        except Exception:
            try:
                logger_tuto.info(getattr(dm_tutorial_text, args[0]).replace('\n','\n\t'))
            except Exception:
                pass
            
        return stop

    def help_output(self):
        """ """
        logger.info("**************** MG5AMC OPTION ***************************")
        super(MadDM_interface, self).help_output()
        logger.info("**************** MADDM NEW OPTION ***************************")
        logger.info("   -- syntax: output [MODE] PATH [OPTIONS]", '$MG:color:BLUE')
        logger.info("")
        logger.info("   MadDM plugin adds two additional MODE for the output of the matrix-element:")
        logger.info("    - indirect:", '$MG:BOLD') 
        logger.info("         is basicaly equivalent to madevent but with a different default run_card.")
        logger.info("         and with some tweak to correctly generate the banner in the context of maddm")
        logger.info("    - maddm:", '$MG:BOLD')
        logger.info("         mode where computation of relic/indirect and direct detection can take place")
        logger.info("         This is the default mode if not explicitly specified")
        logger.info("")        

    def do_output(self, line):
        """ """
        
        if not self._curr_amps:
            self.do_generate('relic_density')
            self.do_add('direct_detection')
            self.do_add('indirect_detection')
            self.do_add('indirect_spectral_features')
            self.history.append('generate relic_density')
            self.history.append('add direct_detection')            
            self.history.append('add indirect_detection')
            self.history.append('add indirect_spectral_features')
        
        args = self.split_arg(line)
        if not args or args[0] not in self._export_formats + ['maddm']:
            if self._curr_amps:
                if any(amp.get('process').get('id') in list(self.process_tag.values())
                        for amp in self._curr_amps):
                    args.insert(0, 'maddm')
            
        if args and args[0] == 'maddm':
            line = ' '.join(args)

        if self._curr_amps:
            super(MadDM_interface, self).do_output(line)
        
        path = self._done_export[0]
        if any(proc for proc in self._ID_cont_procs + self._ID_line_procs):
            self.output_indirect(path, 'Indirect_tree_cont', self._ID_cont_procs[0], self._ID_cont_matrix_elements[0], self._ID_cont_amps[0])
            self.output_indirect(path, 'Indirect_tree_line', self._ID_line_procs[0], self._ID_line_matrix_elements[0], self._ID_line_amps[0])
            self.output_indirect(path, 'Indirect_LI_cont',   self._ID_cont_procs[1], self._ID_cont_matrix_elements[1], self._ID_cont_amps[1])
            self.output_indirect(path, 'Indirect_LI_line',   self._ID_line_procs[1], self._ID_line_matrix_elements[1], self._ID_line_amps[1])

        # find processes names map and write proc_characteristics file
        curr_processes_names_map      = self.extract_processes_names(self._get_processes_from_amplitudes_list(self._curr_amps))
        tree_cont_processes_names_map = self.extract_processes_names(self._get_processes_from_amplitudes_list(self._ID_cont_amps[0]))
        tree_line_processes_names_map = self.extract_processes_names(self._get_processes_from_amplitudes_list(self._ID_line_amps[0]))
        li_cont_processes_names_map   = self.extract_processes_names(self._get_processes_from_amplitudes_list(self._ID_cont_amps[1]))
        li_line_processes_names_map   = self.extract_processes_names(self._get_processes_from_amplitudes_list(self._ID_line_amps[1]))
        # store processes explicitly asked by indirect detection, so that we can consider them only when using fast mode
        indirect_detection_asked_processes = list(tree_cont_processes_names_map.keys()) + list(tree_line_processes_names_map.keys()) + list(li_cont_processes_names_map.keys()) + list(li_line_processes_names_map.keys())
        from . import MGoutput
        proc_path = pjoin(path, 'matrix_elements', 'proc_characteristics')
        proc_charac = MGoutput.MADDMProcCharacteristic(proc_path)
        proc_charac['relic_density_off']      = self._relic_off
        proc_charac['has_indirect_detection'] = self._has_indirect
        proc_charac['has_indirect_spectral']  = self._has_spectral
        proc_charac['forbid_fast']            = self._forbid_fast
        proc_charac['pdg_particle_map']       = self._pdg_particle_map
        proc_charac['processes_names_map']    = { **curr_processes_names_map, **tree_cont_processes_names_map, **tree_line_processes_names_map, **li_cont_processes_names_map, **li_line_processes_names_map }
        proc_charac['indirect_detection_asked_processes'] = indirect_detection_asked_processes
        proc_charac.write(proc_path)

    def _get_processes_from_amplitudes_list(self, amplitude_list):
        return (amplitude["process"] for amplitude in amplitude_list)

    def extract_processes_names(self, processes_iter):
        '''Converts processes names to pdg sequences

        The amplitude name assigned by MadGraph is converted into two forms:
            - human-readable: contains the human-readable particle names;
            - parsing-friendly: contains pdg codes of the particles ready to parse.

        It leverages the methods of the internal Amplitude objects.

        Parameters
        ----------
        processes_iter: iterator over madgraph.core.diagram_generation.Process
            The list of amplitudes to extract the processes names from.

        Returns
        -------
        dict[str,str]
            The processes names with the human-readable string as key and the
            parsing-friendly string as value.

        Notes
        -----
        The name of each process is a sequence of particle strings human-readable
        representing each particle involved in the process. Particles are
        separated using an underscore. E.g. assuming dark matter is called 'n1',
        and has a pdg '52', the following is a possible name for a process:
        n1n1_wpwm_epveemvex
        which represents the process n1 n1 > w+ w-, w+ > e+ ve, w- > e- ve~
        With this method is possible to obtain the parsing-friendly string:
        52.52>6.-6,(23>-11.12),(-23>11.-12)
        where the '>' indicates a process; the '.' separates particles appearing
        in the same initial/final state; the '(...)' contains each decay chain,
        which is a full process using the same rules.
        Moreover, in principle we can also have required s-channels.
        To include them in a parsing-friendly string, we could insert them in 
        '[...]', with the s-channels separated by '|' and the ones formed by 
        multiple particles has them separated with '.'
        Because the required s-channel are only important in the generation and
        they are dropped in the results, they are dropped as well here.

        Example
        -------
        Generating a process with the command:

        n1 n1 > w+|w-|h|z > t t~, (t~ > w- b~, w- > e- ve~), t > u d~ b

        would result in the following human-readable string:

        n1n1_wp_or_wm_or_h_or_z_txt_tx_wmbx_wm_emvex_t_udxb (before computation)
        n1n1_txt_tx_wmbx_wm_emvex_t_udxb (after computation)

        and the parsing-friendly resulting string would be (after computation):

        52.52>6.-6,(-6>-24.-5,(-24>11.-12)),(6>2.-1.5)
        '''

        proc_human_readable = []
        proc_parsing_friendly = []
        for proc in processes_iter:
            str_human_readable = proc.shell_string().split("_", 1)[-1]
            proc_human_readable.append(re.sub(r"_(?:[a-z]+_or_)+[a-z]+", "", str_human_readable, count=0))
            str_parts = [ '.'.join(map(str,proc.get_initial_ids())) + '>' + '.'.join(map(str,proc.get_final_ids())) ]
            decay_chains_dict = self.extract_processes_names(proc["decay_chains"])
            str_parts.extend(map(lambda item: "(" + item + ")", decay_chains_dict.values()))
            proc_parsing_friendly.append(",".join(str_parts))
        return collections.OrderedDict(zip(proc_human_readable, proc_parsing_friendly))

    def output_indirect(self, path, directory, ID_procs, ID_matrix_elements, ID_amps):
        ''' Output commands for indirect_detection or loop-induced '''
        if not ID_amps:
            return

        # force the code to use madevent/reshuffling in case of loop_induced processes
        # if loop_induced:
        #     self._forbid_fast = True

        import aloha.aloha_lib as aloha_lib
        aloha_lib.KERNEL = aloha_lib.Computation()

        with misc.TMP_variable(self,
            ['_curr_proc_defs', '_curr_matrix_elements', '_curr_amps', '_done_export'],
            [ID_procs, ID_matrix_elements, ID_amps, None]):
            super(MadDM_interface, self).do_output('indirect %s/%s' % (path, directory))
            ID_amps[:] = self._curr_amps
            
        #ensure to sync the param_card
        os.remove(pjoin(path, directory, 'Cards', 'param_card.dat'))
        files.ln(pjoin(path, 'Cards', 'param_card.dat'), 
                 pjoin(path, directory, 'Cards'))

        #ensure to sync the param_card
        os.remove(pjoin(path, directory, 'Cards', 'param_card.dat'))
        files.ln(pjoin(path, 'Cards', 'param_card.dat'), 
             pjoin(path, directory, 'Cards'))
            

    def find_output_type(self, path):
        if os.path.exists(pjoin(path,'matrix_elements','proc_characteristics')):
            return 'maddm'
        else:
            return super(MadDM_interface, self).find_output_type(path)

    def check_output(self,*args, **opts):
        
        if 'maddm' in self._export_formats:
            # can happen if running launch directly
            export_formats = list(self._export_formats )
            export_formats.remove('maddm')
            with misc.TMP_variable(self, '_export_formats', export_formats):
                return super(MadDM_interface, self).check_output(*args, **opts)
        else:
            return super(MadDM_interface, self).check_output(*args, **opts)

    def do_launch(self, line):
        
        args = self.split_arg(line)
        (options, args) = madgraph_interface._launch_parser.parse_args(args)
        with misc.TMP_variable(self, '_export_formats', self._export_formats + ['maddm']):
            self.check_launch(args, options)
        options = options.__dict__        
        

        if args[0] != 'maddm':
            return super(MadDM_interface, self).do_launch(line)
        else:
            self._MDM = maddm_run_interface.MADDMRunCmd(dir_path=args[1], options=self.options)
            self._done_export = (args[1], 'plugin')
            
            if options['interactive']:
                return self.define_child_cmd_interface(self._MDM)
            else:
                self.define_child_cmd_interface(self._MDM,  interface=False)
                try:
                    self._MDM.exec_cmd('launch ' + line.replace(args[1], ''))
                except:
                    self._MDM.exec_cmd('quit')
                    raise
                else:
                    self._MDM.exec_cmd('quit')
                self._done_export = (args[1], 'plugin')
                return


    def define_multiparticles(self, label, list_of_particles):
        """define a new multiparticle from a particle list (add both particle and anti-particle)"""
        
        pdg_list = []
        for p in list_of_particles:
            if p.get('name') == p.get('antiname'):
                pdg_list.append(p.get('pdg_code'))
            else:
                pdg_list.append(p.get('pdg_code'))
                pdg_list.append(-1*p.get('pdg_code'))

        self.optimize_order(pdg_list)
        self._multiparticles[label] = pdg_list
        

    def generate_relic(self, excluded_particles=[], user = True):
        """--------------------------------------------------------------------#
        #                                                                      #
        #  This routine generates all the 2 -> 2 annihilation matrix elements. #
        #  The SM particles, BSM particles in self._bsm_final_states, and the  #
        #  DM particles are used as final states.  The initial states are      #
        #  looped over all the different combinations of DM particles.         #
        #  if user = True, then it would mean the user has run this command.   #
        #--------------------------------------------------------------------"""

        #Here we define _bsm_ multiparticle as a special case for exlusion only in the
        #final state. This is different than the standard / notation which will exclude
        #the particles from the propagators as well.
        #Since we use the exluded_particles everywhere, we just use the presence of
        #the _bsm_ in the array to initialize a flag and then remove it from the excluded
        #particles array so as not to conflict with the use in other places.
        self._relic_off = not user # if the user has not asked explicitly for relic density, then its default value in the launch interface is 'OFF'
        if self.relic_generated_by_id:
            logger.info("Relic density computation has been switched on.")
            return

        self.define_multiparticles('_bsm_',[p for p in self._curr_model.get('particles')\
                         if abs(p.get('pdg_code')) > 25 and str(p.get('pdg_code')) not in self._unphysical_particles])
        if 'bsm' in excluded_particles:
            exclude_bsm = True
            while 'bsm' in excluded_particles:
                excluded_particles.remove('bsm')
        else:
            exclude_bsm = False

        if not self._dm_candidate:
            self.search_dm_candidate(excluded_particles)
            if not self._dm_candidate:
                return
            self.search_coannihilator()
            self.history.append('add relic_density %s %s' % ( '/' if excluded_particles else '',
                                                      ' '.join(excluded_particles)))


        # Tabulates all the BSM particles so they can be included in the
        # final state particles along with the SM particles.
        
        #dm_mass = abs(self._curr_model.get_mass(self._dm_candidate[0]))
        #is_lower_mass = lambda p : True #abs(self._curr_model.get_mass(p)) < (1-self.coannihilation_diff)*dm_mass
        
        ids_veto = [p.get('pdg_code') for p in self._dm_candidate + self._coannihilation]     
        bsm_final_states = [p for p in self._curr_model.get('particles') 
                         if abs(p.get('pdg_code')) > 25 and \
                         (p.get('name') not in excluded_particles or p.get('antiname') not in excluded_particles) and\
                         p.get('width') != 'ZERO' and \
                         abs(p.get('pdg_code')) not in ids_veto and \
                         str(p.get('pdg_code')) not in self._unphysical_particles and \
                         abs(p.get('pdg_code')) not in [abs(z.get_pdg_code()) for z in self._z2_odd+self._disconnected_particles]]

        if not exclude_bsm:
            if bsm_final_states:
                logger.info("DM is allowed to annihilate into the following BSM particles: %s",
                         ' '.join([p.get('name') for p in bsm_final_states]))
                logger.info("use generate relic_density / X to forbid the decay the annihilation to X")
                logger.info("if you want to forbid DM to annihilate into BSM particles do relic_density / bsm")
        else:
            bsm_final_states=[]

        # Set up the initial state multiparticles that contain the particle 
        #and antiparticle
        for i,dm in enumerate(self._dm_candidate + self._coannihilation):
            self.define_multiparticles('dm_particle_%s' % i,  [dm])
            #self.do_display('multiparticles')

        self.define_multiparticles('dm_particles', self._dm_candidate+self._coannihilation)

        self.do_display('multiparticles')
        sm_pdgs = list(range(1, 7)) + list(range(11, 17)) + list(range(21, 26)) #quarks/leptons/bosons
        sm_part = [self._curr_model.get_particle(pdg) for pdg in sm_pdgs if self._curr_model.get_particle(pdg)]

        self.define_multiparticles('fs_particles', sm_part +bsm_final_states)

        # generate annihilation diagram
        coupling = "QED<=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0"
        if excluded_particles:
            proc = "dm_particles dm_particles > fs_particles fs_particles / %s %s  @DM2SM" \
                    % (' '.join(excluded_particles), coupling)
        else:
            proc = "dm_particles dm_particles > fs_particles fs_particles %s @DM2SM" \
                   %(coupling)  #changed here

        # reset the last amps
        self._last_amps = diagram_generation.AmplitudeList()
        
        self.do_add('process %s' %proc, user=False)

        # Generate all the DM -> DM processes
        nb_dm = len(self._dm_candidate + self._coannihilation)
        for i in range(nb_dm):
            for j in range(i,nb_dm):
                if  excluded_particles:
                    proc = "DM_particle_%s DM_particle_%s > dm_particles dm_particles / %s %s @DM2DM"\
                       % (i,j,' '.join(excluded_particles), coupling)
                else:
                    proc = "DM_particle_%s DM_particle_%s > dm_particles dm_particles %s @DM2DM"\
                       % (i,j, coupling)
                try:
                    self.do_add('process %s' % proc, user=False)
                except (self.InvalidCmd,diagram_generation.NoDiagramException) :
                    continue
        # Get the matrix elements and make sure that we don't have any pure 
        #scattering processes
        for amp in self._curr_amps[:]:
            if amp.get('process').get('id') != self.process_tag['DM2DM']:
                continue
            if set(amp.get('process').get_initial_ids()) == (set(amp.get('process').get_final_ids())):
                self._curr_amps.remove(amp)


        # Generate all the DM particles scattering off of the thermal background and
        # change to a different DM species (again, no pure scatterings)
        for i in range(nb_dm):
            for j in range(nb_dm):
                if i == j:
                    continue
                if excluded_particles:
                    proc = "DM_particle_%s sm_particles > DM_particle_%s sm_particles / %s %s @DMSM"\
                            % (i,j,' '.join(excluded_particles), coupling)
                else:
                    proc = "DM_particle_%s sm_particles > DM_particle_%s sm_particles %s @DMSM"\
                            % (i,j, coupling)
                try:
                    self.do_add('process %s' % proc, user=False)
                except (self.InvalidCmd,diagram_generation.NoDiagramException) :
                    continue

        self._last_amps = [amp for amp in self._curr_amps if amp.get('process').get('id') in [self.process_tag['DM2SM'], self.process_tag['DMSM']]]

    def generate_direct(self, excluded_particles=[]):
        """User level function which performs direct detection functions        
           Generates the DM - q,g scattering matrix elements for spin dependent 
           and spin independent direct detection cross section                   
           Currently works only with canonical mode and one dm candidate.        
           The function also merges the dark matter model with the effective        
           vertex model.          
        """
        if not self._dm_candidate:
            self.search_dm_candidate(excluded_particles)
            if not self._dm_candidate:
                return
            else:
                self.history.append('add direct_detection %s %s' % ( '/' if excluded_particles else '',
                                                      ' '.join(excluded_particles)))

        if len(self._dm_candidate) > 1:
            logger.warning("More than one DM candidate. Can not run Direct Detection.")
            return 
  
        #generate a special
        
        #Now figure out the label of the effective vertex to use. The convention is:
        #<SI or SD>EFF<F, S, or V>, i.e. SIEFFV for Spin Independent Effective vertex for Vector Current.
        dm_spin = int(self._dm_candidate[0]['spin'])
        eff_operators_SI = self.eff_operators_SI[dm_spin]
        eff_operators_SD = self.eff_operators_SD[dm_spin]
        
        logger.info("Generating X Nucleon > X Nucleon diagrams from the full lagrangian...")
        has_direct = self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'QED', excluded_particles)

        if not has_direct:
            logger.warning("No Direct Detection Feynman Diagram")
            return
        
        logger.info("Generating X Nucleon > X Nucleon diagrams from the effective lagrangian...")
        #ONLY EFFECTIVE LAGRANGIAN
        self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SI',excluded_particles)

        logger.info("INFO: Generating X Nucleon > X Nucleon diagrams from the effective+full lagrangian...")
        #EFFECTIVE + FULL
        self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SI+QED',excluded_particles)
        
        if (eff_operators_SD != False):
            logger.info("Doing the spin dependent part...")
            logger.info("Generating X Nucleon > X Nucleon diagrams from the effective lagrangian...")

            self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SD',excluded_particles)
            #EFFECTIVE + FULL
            logger.info("Generating X Nucleon > X Nucleon diagrams from the effective + full lagrangian...")
            self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SD+QED',excluded_particles)

    #-----------------------------------------------------------------------#
    @misc.mute_logger(['madgraph','aloha','cmdprint'], [30,30,30])
    def DiagramsDD(self, SI_name, SD_name, type, excluded=[]):
        """Generates direct detection diagrams. i_dm is the index of DM part. 
                 Whether spin dependent or spin independent diagrams should be                 
                 calculated. XX_order parameters determine maximum order of  a
                 coupling. For instance, if you only want effective vertices for
                 spin independent coupling you would set SI_order = 2 and all others
                 to zero. If you want the spin independent full lagrangian + eff.
                 then you need to set SI_order=2 and QED_order=2...
                 WARNING: This function is a proxy to be used inside
                 GenerateDiagramsDirDetect() function and no where else!
        """                                                             

        quarks = list(range(1,7)) # ['d', 'u', 's', 'c', 'b','t']
        antiquarks = [-1*pdg for pdg in quarks] # ['d~', 'u~', 's~', 'c~', 'b~','t~']

        def format_order(name, order):
            if not name:
                return ''
            return name+order

        if type == 'SI':
            orders = format_order(SI_name, '==2')
        elif type == 'SD':
            orders = format_order(SD_name, '==2')
        elif type == 'QED':
            orders = format_order(SD_name, '=0') + ' ' + format_order(SI_name, '=0')
        elif type == 'SI+QED':
            orders = format_order(SD_name, '=0') + ' ' + format_order(SI_name, '<=99')
        elif type == 'SD+QED':
            orders = format_order(SD_name, '<=99') + ' ' + format_order(SI_name, '=0')
                
        #loop over quarks
        has_diagram = False

        for i in quarks + antiquarks:
            proc = ' %(DM)s %(P)s > %(DM)s %(P)s %(excluded)s %(orders)s @DD' %\
                    {'DM': self._dm_candidate[0].get('name'),
                     'P': i,
                     'excluded': ('/ %s' % ' '.join(excluded) if excluded else ''),
                     'orders': orders
                     }

            try:
                self.do_add('process %s' %proc, user=False)
            except (self.InvalidCmd, diagram_generation.NoDiagramException) as error:
                logger.debug(error)
            else:
                has_diagram = True

            if self._dm_candidate[0].get('antiname') != self._dm_candidate[0].get('name'):
                proc = ' %(DM)s %(P)s > %(DM)s %(P)s %(excluded)s %(orders)s @DD' %\
                    {'DM': self._dm_candidate[0].get('antiname'),
                     'P': i,
                     'excluded': ('/ %s' % ' '.join(excluded) if excluded else ''),
                     'orders': orders
                     }

                try:
                    self.do_add('process %s' % proc, user=False)
                except (self.InvalidCmd, diagram_generation.NoDiagramException) as error:
                    logger.debug(error)
                else:
                    has_diagram = True
           
        self._last_amps = [amp for amp in self._curr_amps if amp.get('process').get('id') == self.process_tag['DD']]
        return has_diagram

    def check_indirect_and_spectral(self, cmd, argument):
        ''' takes the name of the command in cmd (either 'indirect_detection' or 'indirect_spectral_features')
            check that everything is ok and returns the value to set to _has_<cmd> and a number: if != 0 the generation function should return, otherwise should go ahead
        '''
        if not maddm_run_interface.HAS_NUMPY:
            logger.warning("numpy module not detected on your machine. \n"+
                           "Running indirect detection will not be possible as long as numpy is not installed (try 'pip install numpy')" 
                           )
            ans = self.ask("Do you want to generate the diagrams anyway?", 'n', ['y','n'])
            if ans == 'n':
                return False, 2
        
        if cmd == 'indirect_detection':
            self.install_indirect()

        if not self._dm_candidate:
            self.search_dm_candidate()
            if not self._dm_candidate:
                return
            else:
                self.history.append('add %s %s' % (cmd, ' '.join(argument)))

        if len(self._curr_proc_defs) == 0 or \
           all(p.get('id') != self.process_tag['DM2SM'] for p in self._curr_proc_defs):
            logger.info('For indirect detection we need to generate relic density matrix-element','$MG:BOLD')
            self.generate_relic([], user = False)
            self.relic_generated_by_id = True

        # check the '2to2lo' option: forbid it if loop-induced or 2 to 3 processes have been considered
        if '2to2lo' in argument:
            if not self._forbid_fast:
                self.option_2to2lo = True
            else:
                logger.error("Can't use '2to2lo' option for indirect detection if loop-induced or 2 to 3 processes are considered.")
            return True, 1

        return True, 0

    def is_amplitude_not_empty(self, amplitude):
        ''' tree level amplitudes have only property 'diagrams', while for loop amplitudes have both 'diagrams', 'born_diagrams' and 'loop_diagrams', so we need to check all of them to find out if a process has diagrams. '''
        if isinstance(amplitude, diagram_generation.DecayChainAmplitude):
            return bool(amplitude['decay_chains'])
        
        diagrams = bool(amplitude['diagrams'])
        born_diagrams = bool('born_diagrams' in list(amplitude.keys()) and amplitude['born_diagrams'])
        loop_diagrams = bool('loop_diagrams' in list(amplitude.keys()) and amplitude['loop_diagrams'])
        return any([diagrams, born_diagrams, loop_diagrams])

    def indirect_contains_photon(self, final_states):
        if len(final_states.split()) != 2:
            return False
        for arg in final_states.split():
            try:
                if self._pdg_particle_map.get_pdg(arg) == 22:
                    return True
            except ValueError as err:
                # check if it is a multiparticle
                try:
                    if 22 in self._multiparticles[arg]:
                        return True
                except KeyError:
                    pass
        return False

    def generate_indirect(self, argument):
        """User level function which performs indirect detection functions
           Generates the DM DM > X where X is anything user specified.
           Also works for loop induced processes as well as NLO-QCD.
           Currently works only with canonical mode and one dm candidate.
        related to syntax: generate indirect a g / n3
        """
        self._has_indirect, state = self.check_indirect_and_spectral(cmd = 'indirect_detection', argument = argument)
        if state:
            return
        
        # check if there are photons in final states: if so, redirect the user to use 'indirect_spectral_features'
        if "," not in ' '.join(argument) and self.indirect_contains_photon(' '.join(argument).split('/')[0]) and len(' '.join(argument).split('/')[0]) != 0:
            logger.error("Processes with at least one photon in the final state must be generated through 'indirect_spectral_features' command.")
            self._has_indirect = False
            return
        
        # Generates events for all the DM annihilation channels (also BSM)
        # it uses the '--noloop' option
        backup_amps = [amp for amp in self._ID_cont_amps[0] + self._ID_cont_amps[1]]
        if argument == []:
            # bsm_content = self.indirect_bsm_multiparticle()
            odd_and_disconnected = self.list_odd_disconnected_pdg()
            id_fs_particle = ['all', '/', '22'] + self._unphysical_particles + odd_and_disconnected
            self.exec_cmd('define _id_fs_ = %s' % ' '.join(str(x) for x in id_fs_particle), postcmd=False)
            final_states = ['_id_fs_ _id_fs_']
            # self.exec_cmd('define _q_mdm_ = 1 2 3 4 -1 -2 -3 -4', postcmd=False)
            # set final states: notice that '22 22' has been removed, because now it is handles by 'indirect_spectral_features'
            # final_states = ['_q_mdm_ _q_mdm_', '21 21', '5 -5', '6 -6', '23 23', '24 -24', '25 25', '11 -11', '13 -13', '15 -15', '12 -12', '14 -14', '16 -16']
            # final_states += ['_bsm_ _bsm_'] if bsm_content else []
            for final_state in final_states:
                try: 
                    self.indirect_process_generation([final_state, '--noloop'], self._ID_cont_procs, self._ID_cont_matrix_elements, self._ID_cont_amps)
                except diagram_generation.NoDiagramException:
                    continue
                    # logger.info('no diagram for %s' % final_state)
            self._last_amps = [amp for amp in self._ID_cont_amps[0] + self._ID_cont_amps[1] if amp not in backup_amps]
            self._has_indirect = len(self._last_amps) != 0 or len(backup_amps) != 0
            return

        self.indirect_process_generation(argument, self._ID_cont_procs, self._ID_cont_matrix_elements, self._ID_cont_amps)
        self._last_amps = [amp for amp in self._ID_cont_amps[0] + self._ID_cont_amps[1] if amp not in backup_amps]
        self._has_indirect = len(self._last_amps) != 0 or len(backup_amps) != 0
    
    def generate_spectral(self, argument):
        """ Performs indirect detection in the final states 'a a', 'a z', 'a h' and 'a _bsm_'.
            The processes are loop-induced if the model does not provide effective vertices.
            The user can specify a single final state to analyze.
        """
        print(("-"*99))
        logger.warning("You are using the module indirect_spectral_feature. This module is still in beta release.\nInstabilities at small velocities, vave_indirect_line < 1e-3c, might occur. Consider choosing larger values of vave_indirect_line if possible.")
        print(("-"*99))

        self._has_spectral, state = self.check_indirect_and_spectral(cmd = 'indirect_spectral_features', argument = argument)
        if state:
            return

        # forbid decay syntax
        if "," in ' '.join(argument):
            logger.error("Decay syntax not allowed for line searches.")
            self._has_spectral = False
            return
        # check if at least one photon is present in the final state
        if not self.indirect_contains_photon(' '.join(argument).split('/')[0]) and len(' '.join(argument).split('/')[0]) != 0:
            logger.error("There must be at least one photon in the final state and (for now) only 2-body processes are allowed.")
            self._has_spectral = False
            return

        # Generates events for all the DM annihilation channels (also BSM)
        backup_amps = [amp for amp in self._ID_line_amps[0] + self._ID_line_amps[1]]
        if argument == []:
            # exclude electric and color charged particles because of course they can't generate 'a _bsm_' diagrams
            charged_particles = [pdg_code for pdg_code, particle in six.iteritems(self._curr_model.get('particle_dict')) if particle.get('charge') != 0. or particle.get('color') != 1]
            odd_and_disconnected = self.list_odd_disconnected_pdg()
            bsm_content = self.indirect_bsm_multiparticle(excluded = charged_particles+odd_and_disconnected)
            final_states = ['22 22', '22 23', '22 25']
            final_states += ['22 _bsm_'] if bsm_content else []
            for final_state in final_states:
                try: 
                    self.indirect_process_generation([final_state], self._ID_line_procs, self._ID_line_matrix_elements, self._ID_line_amps)
                except diagram_generation.NoDiagramException:
                    continue
                    # logger.info('no diagram for %s' % final_state)
            self._last_amps = [amp for amp in self._ID_line_amps[0] + self._ID_line_amps[1] if amp not in backup_amps]
            self._has_spectral = len(self._last_amps) != 0 or len(backup_amps) != 0
            return

        self.indirect_process_generation(argument, self._ID_line_procs, self._ID_line_matrix_elements, self._ID_line_amps)
        self._last_amps = [amp for amp in self._ID_line_amps[0] + self._ID_line_amps[1] if amp not in backup_amps]
        self._has_spectral = len(self._last_amps) != 0 or len(backup_amps) != 0

    def indirect_bsm_multiparticle(self, excluded = []):
        ''' it defines the '_bsm_' particle for indirect detection processes generation,
            and it returns its content.
            it is possible to specify a list of excluded particles in addition to the already excluded SM ones.
        '''
        bsm_all  = [ 'all', '/' ] + list(set(['1','2','3','4','5','6','-1','-2','-3','-4','-5','-6','21','22','23','24','-24','25','11','12','13','14','15','16','-11','-12',
                '-13','-14','-15','-16'] + self._unphysical_particles + excluded))
        bsm = " ".join(str(x) for x in bsm_all)
        self.exec_cmd('define _bsm_ = %s' % bsm, postcmd=False)
        return self._multiparticles['_bsm_']

    def list_odd_disconnected_pdg(self):
        odd_and_disconnected = []
        for particle in self._z2_odd + self._disconnected_particles:
            odd_and_disconnected.extend([particle.get_pdg_code(), particle.get_anti_pdg_code()])
        return list(set(odd_and_disconnected))

    def indirect_process_generation(self, argument, ID_procs, ID_matrix_elements, ID_amps):
        ''' generate the processes for indirect detection, either for continuum spectrum (generate_indirect) or line spectrum (generate_spectral).
            'ID_' variables are the lists which will contains the processes/matrix elements/amplitudes definitions: [0]tree level, [1]loop-induced.
        '''
        if not self.options['pythia8_path']:
            logger.warning('''In order to have distribution related to indirect detection:
             Pythia8 needs to be linked to the code. You can install it by typing 'install pythia8'.
             You can compute the rate (no distribution) by typing 'add indirect 2to2lo'. ''')
        #    return

        allow_loops = True
        if '--noloop' in argument:
            argument.remove('--noloop')
            allow_loops = False

        # Check if the argument contains only two particles in the final state
        # if not, force the code to use madevent/reshuffling
        if (len(' '.join(argument).split('/')[0].split())!=2):
            self._forbid_fast = True 

        # First try LO matrix-element
        coupling = "SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0"
        done= []

        # AmplitudeList containing the last amplitudes generated by different dark matter candidates for the same final state
        temp_amps = diagram_generation.AmplitudeList()
        for dm in self._dm_candidate:
            name = dm.get('name')
            antiname = dm.get('antiname')
            if name in done:
                continue
            done += [name, antiname]
            # variable containing the last amplitude generated in this iteration
            last_amp = diagram_generation.Amplitude()
            #We put the coupling order restrictions after the @ID in order to
            #apply it to the entire matrix element.
            proc = '%s %s > %s @ID %s' % (name, antiname, ' '.join(argument), coupling)
            try:
                # flip indirect/standard process definition
                with misc.TMP_variable(self, ['_curr_proc_defs', '_curr_matrix_elements', '_curr_amps'], 
                                     [ID_procs[0], ID_matrix_elements[0], ID_amps[0]]):
                    self.do_add('process %s' % proc, user=False)
                    last_amp = self._curr_amps[-1]
            except diagram_generation.NoDiagramException as error:
                if not self.model_nlo:
                    logger.error(". ".join(error.args))
                if allow_loops and self.model_nlo:
                    with misc.TMP_variable(self, ['_curr_proc_defs', '_curr_matrix_elements', '_curr_amps'], 
                                     [ID_procs[1], ID_matrix_elements[1], ID_amps[1]]):
                        proc = '%s %s > %s %s [noborn=QCD] @ID ' % (name, antiname, ' '.join(argument), coupling)
                        try:                    
                            self.do_add('process %s' % proc, user=False)
                            last_amp = self._curr_amps[-1]
                            self._forbid_fast = True
                        except (self.InvalidCmd, diagram_generation.NoDiagramException) as error:
                            if "cannot handle loop processes" in error.args[0]:
                                logger.error("No amplitudes have been generated at tree-level, and this model does not allow for loop-induced computations.")
                                self.model_nlo = False
                                self.change_principal_cmd('MadGraph')
                                continue
                            proc = '%s %s > %s %s [noborn=QED] @ID ' % (name, antiname, ' '.join(argument), coupling)
                            try:
                                self.do_add('process %s' % proc, user=False)
                                last_amp = self._curr_amps[-1]
                                self._forbid_fast = True
                            except (self.InvalidCmd, diagram_generation.NoDiagramException) as error:
                                # if no loop-induced diagrams have been found then change back the command interface to MadGraph
                                logger.error(error.args[0].replace('noborn = QED', 'noborn = QCD or QED'))
                                self.change_principal_cmd('MadGraph')
            except self.InvalidCmd as error:
                logger.error(". ".join(error.args))
            if self.is_amplitude_not_empty(last_amp):
                temp_amps.append(last_amp)

        # check if the option_2to2lo is True while we have generated a LI process
        if self.option_2to2lo and self._forbid_fast:
            logger.error("Can't use '2to2lo' option for indirect detection if loop-induced or 2 to 3 processes are considered. Please run indirect detection without that option.")
            self.option_2to2lo = False

        return temp_amps # I return the last amps of all the DM candidates, because this function is called for each final state in turn

    def install_indirect(self):
        """Code to install the reduction library if needed"""
        
        opt = self.options
        
        # Check if first time:
        if not opt['maddm_first_indirect']:
            return
        
        logger.info("First time that you asked for indirect detection. Now asking for dependency tool:", '$MG:BOLD')
        to_install = self.ask('install', '0',  ask_class=AskMadDMInstaller, timeout=300, 
                              path_msg=' ')
        
        rename = {'dragon_data_from_galprop': 'dragon_data'}
        key_to_opts = {'PPPC4DMID':'pppc4dmid_path',
                       'dragon': 'dragon_path',
                       'dragon_data_from_galprop': None,
                       'pythia8': 'pythia8_path'}
        
        for key, value in to_install.items():
            if value == 'install':
                if key in rename:
                    key = rename[key]
                self.exec_cmd('install %s' % key)
            # Not install
            elif value == 'off' and key_to_opts[key]:
                key =  key_to_opts[key]
                self.exec_cmd("set %s '' " % key)
                self.exec_cmd('save options %s' % key)
            elif key_to_opts[key]:
                key =  key_to_opts[key]
                self.exec_cmd("set %s %s" % (key,value))
                self.exec_cmd('save options %s' % key) 
        self.exec_cmd('set maddm_first_indirect False')
        self.exec_cmd('save options maddm_first_indirect')
                
            
  
    def update_model_with_EFT(self):
        """ """
        eff_operators_SD = {1:False, 2:'SDEFFF', 3:'SDEFFV'}
        eff_model_dm_names = {1:'~sdm', 2:'~fdm', 3:'~vdm'}

        DM = self._dm_candidate[0]
        
        if self._dm_candidate[0]['self_antipart']: 
            EFT = 'REAL'
        else:
            EFT = 'COMPLEX'

        spin = int(DM['spin'])
        if spin not in eff_model_dm_names:
            logger.warning('This type of spin is NOT supported for direct-detection')
        else:
            eff_dm_name = eff_model_dm_names[int(DM['spin'])]
            mg5_command = 'add model %s %s=%s --recreate --keep_decay' % (pjoin(MDMDIR, 'EffOperators', EFT),
                                                     eff_dm_name,DM.get('name'))
        
            # We want to preserve the following variable while updating the model
            backup_amp = self._curr_amps
            backup_param_card = self._param_card
            backup_dm_candidate = self._dm_candidate
            backup_coannihilation = self._coannihilation
            backup_disconnected_particles = self._disconnected_particles
        
        
            self.exec_cmd(mg5_command, postcmd=False)
            self._curr_amps = backup_amp
            self._param_card = backup_param_card 
            self._dm_candidate = backup_dm_candidate
            self._coannihilation = backup_coannihilation
            self._disconnected_particles = backup_disconnected_particles

        # update the param_card value
        txt = self._curr_model.write_param_card()
        param_card = check_param_card.ParamCard(self._curr_model.write_param_card())
        
        if self._param_card:
            for block in self._param_card:
                if block not in param_card:
                    logger.debug('%s not valid block entry in the card' ,block)
                    continue   
                for param in self._param_card[block]:
                    try:
                        param_card[block].get(param.lhacode).value =  param.value
                    except KeyError:
                        continue #possible in MSSM due to restriction on SLHA2 format
 
        self._param_card = param_card        
        if not isinstance(self._curr_model, model_reader.ModelReader):
            self._curr_model = model_reader.ModelReader(self._curr_model)
        self._curr_model.set_parameters_and_couplings(self._param_card) 
        
import madgraph.interface.common_run_interface as common_run
import madgraph.interface.extended_cmd as cmd

class EditParamCard(common_run.AskforEditCard):
    """a dedicated module for the param"""

    special_shortcut ={}

    def __init__(self, question, card=[], mode='auto', *args, **opt):

        if not card:
            card = ['param']
        return super(EditParamCard,self).__init__(question, card, mode=mode, *args, **opt)

    def do_asperge(self, *args, **opts):
        "Not available"
        logger.warning("asperge not available in this mode")
        
    def do_help(self,*args, **opts):
        """ print standard help """
        
        out =  super(EditParamCard, self).do_help(*args, **opts)

        logger.info("""In order to automatically determine your darkmatter candidate, we need to have a benchmark point.""")
        logger.info("")
        logger.info("  For the darkmatter, we will apply the following algorithm:")
        logger.info("      1. The DM particle must be a BSM particle (pdg_code > 25) ")  
        logger.info("      2. The particle should have no charge (electric or color)")
        logger.info("      3. The particle\'s width should be 0 or 'ZERO'")
        logger.info("      4. The particle should have at least one non vanishing coupling" )
        logger.info("      5. The assigned DM candidate is the lightest of all the particles that meet the above criteria.")  
        logger.info("")
        logger.info("  For coannihilation, we apply the following algorithm:")
        # logger.info("      The code selects all the BSM particles that are within an input mass difference with the DM candidate.")
        # logger.info("      Particles without any non vanishing coupling are discarded")
        logger.info("      We consider a Z2 symmetry under which dark matter candidate is odd, while Standard Model particles are even.")
        logger.info("      The code goes through the interaction vertices and determines the Z2 parity of the BSM particles in the model.")
        logger.info("      The coannnihilators are the Z2-odd particles.")
 
        logger_tuto.info("""
This card is here ONLY to allow to determine automatically darkmatter/coannihilator.
You will have the possibility to edit it later to define another benchmark, perform scan/...
But the diagram generated will depend of the dark matter selected at this stage.

To edit a parameter of the card without having to use a text editor your can type
>set mxd 5

When you are done editing the card, just press enter ( you can also type done or 0)
""")
            
        return out
        
import madgraph.interface.loop_interface as loop_interface
class AskMadDMInstaller(loop_interface.AskLoopInstaller):
    
    local_installer = []
    required = []
    order = ['pythia8', 'PPPC4DMID', 'dragon', 'dragon_data_from_galprop']
    bypassed = []
    
   
    def __init__(self, question, *args, **opts):
        
        self.code = dict([(name, 'install') for name in self.order])
        self.online=True
        #check if some partial installation is already done.  
        if 'mother_interface' in opts:
            mother = opts['mother_interface']
            if  'heptools_install_dir' in mother.options:
                install_dir = mother.options['heptools_install_dir']
            else:
                install_dir = pjoin(MG5DIR, 'HEPTools')


            for c in ['pythia8', 'dragon']:
                if os.path.exists(pjoin(install_dir, c)):
                    self.code[c] =  pjoin(install_dir, c)
                    
            if os.path.exists(pjoin(MG5DIR, 'PPPC4DMID')):
                c= 'PPPC4DMID'
                self.code[c] =  pjoin(MG5DIR, c)
            if os.path.exists(pjoin(install_dir, 'dragon','data')):
                self.code['dragon_data_from_galprop'] =  pjoin(install_dir, 'dragon','data')

        
        # 1. create the question
        question, allowed_answer = self.create_question(first=True)
        
        opts['allow_arg'] = allowed_answer
        
        cmd.OneLinePathCompletion.__init__(self, question, *args, **opts)
    
    def create_question(self, first = False):
        """ """

        question = "For indirect detection, MadDM relies on some external tools." +\
                   "You can decide here which one you want to include." +\
                   "Not installing all the dependencies will limit functionalities."+\
                   "\nWhich one do you want to install? (this needs to be done only once)\n"
        
        allowed_answer = set(['0','done'])

        #order = ['pythia8', 'PPPC4DMID', 'dragon', 'dragon_data_from_galprop']
                
        status = {'off': '%(start_red)sdo not install%(stop)s',
                  'install': '%(start_green)swill be installed %(stop)s',
                  'local': '%(start_green)swill be installed %(stop)s(offline installation from local repository)',
                  }
        
        descript = {'pythia8': ['pythia8', ' shower (precise mode)','[1410.3012]'],
                    'PPPC4DMID': ['PPPC4DMID', 'all (fast mode)' , '[1012.4515]'],
                    'dragon': ['dragon', 'propagation (precise mode)' , '[0807.4730]'],
                    'dragon_data_from_galprop':['dragon_data_from_galprop', 'input for dragon', '[1712.09755]']
                    }
        
        
        for i,key in enumerate(self.order,1):
            if key in self.bypassed and self.code[key] == 'off':
                continue
            if os.path.sep not in self.code[key]:
                question += '%s. %%(start_blue)s%-9s %-5s %-13s%%(stop)s : %s%s\n' % \
                   tuple([i,]+descript[key]+[status[self.code[key]],]+\
                     ['(recommended)' if key in ['ninja','collier'] and self.code[key] in ['install'] else ''])
            else:
                question += '%s. %%(start_blue)s%-9s %-5s %-13s%%(stop)s : %s\n' % tuple([i,]+descript[key]+[self.code[key],])
            if key in self.required:
                continue
            allowed_answer.update([str(i), key])
            if key in self.local_installer:
                allowed_answer.update(['key=local','key=off'])

                
        question += "You can:\n -> hit 'enter' to proceed\n -> type a number to cycle its options\n -> enter the following command:\n"+\
          '    %(start_blue)s{tool_name}%(stop)s [%(start_blue)sinstall%(stop)s|%(start_blue)snoinstall%(stop)s|'+\
          '%(start_blue)s{prefixed_installation_path}%(stop)s]\n'
        if first:
            question += '\n%(start_bold)s%(start_red)sIf you are unsure about what this question means, just type enter to proceed. %(stop)s'

        question = question % {'start_green' : '\033[92m',
                               'start_red' : '\033[91m',
                               'start_blue' : '\033[34m',
                               'stop':  '\033[0m',
                               'start_bold':'\033[1m', 
                               }
        return question, allowed_answer
    
    do_pythia8 = lambda self,line : self.apply_name('pythia8', line)
    do_PPPC4DMID = lambda self,line : self.apply_name('PPPC4DMID', line)
    
    def do_dragon(self, line):
        
        self.apply_name('dragon', line)
        if self.code['dragon'] in ['install']:
            self.code['dragon_data_from_galprop'] = 'install'
        elif self.code['dragon'] != 'off':
            if os.path.exists(pjoin(self.code['dragon'],'data')):
                self.apply_name('dragon_data_from_galprop', pjoin(self.code['dragon'],'data'))
            else:
                self.apply_name('dragon_data_from_galprop', 'install')
        else:
            self.apply_name('dragon_data_from_galprop', 'noinstall')
            #self.code['dragon_data_from_galprop'] = 'off'

    def do_dragon_data(self, line):
        
        self.apply_name('dragon_data_from_galprop', line)
        if self.code['dragon_data_from_galprop'] in ['install']:
            if self.code['dragon'] == 'off':
                self.code['dragon'] = 'install'

    do_dragon_data_from_galprop = do_dragon_data

    complete_pythia8 = loop_interface.AskLoopInstaller.complete_prog
    complete_PPPC4DMID = loop_interface.AskLoopInstaller.complete_prog
    complete_dragon = loop_interface.AskLoopInstaller.complete_prog
    complete_dragon_data = loop_interface.AskLoopInstaller.complete_prog
    complete_dragon_data_from_galprop = loop_interface.AskLoopInstaller.complete_prog

    
    

import logging
import os

import maddm_run_interface as maddm_run_interface
import dm_tutorial_text as dm_tutorial_text

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


class DMError(Exception): pass

# Root path
MDMDIR = os.path.dirname(os.path.realpath( __file__ ))


class MadDM_interface(master_interface.MasterCmd):

    intro_banner=\
  "            ====================================================\n"+\
  "            |                  "+bcolors.OKBLUE+"  MadDM v3.0                     "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "                                                                               \n"+\
  "                #########                                                        \n"+\
  "             ###\\\\####//#####              Launchpad:  launchpad.net/maddm      \n"+\
  "           ######\\\\##//########                                              \n"+\
  "          ########\\\\//###########                                            \n"+\
  "         #########//\\\\############                    "+bcolors.FAIL+"arXiv:1308.4955        \n"+bcolors.ENDC+\
  "        #########//##\\\\############                   "+bcolors.FAIL+"arXiv:1505.04190        \n"+bcolors.ENDC+\
  "       ########//#####\\\\###########                   "+bcolors.FAIL+"arXiv:1804.00444        \n"+bcolors.ENDC+\
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
    options_configuration['pppc4dmid_path'] = './PPPC4DMID'
    options_configuration['dragon_path'] = None
    options_configuration['maddm_first_indirect'] = True
    
    _set_options = list(master_interface.MasterCmd._set_options)
    _set_options.append('dragon_data_from_galprop')
    
    
    def post_install_PPPC4DMID(self):
        if os.path.exists(pjoin(MG5DIR, 'PPPC4DMID')):
            self.options['pppc4dmid_path'] = pjoin(MG5DIR, 'PPPC4DMID')
        
        if not maddm_run_interface.HAS_SCIPY:
            logger.critical("Fermi-LAT limit calculation for indirect detection requires numpy and scipy. Please install the missing modules.")
            logger.info("you can try to use \"pip install scipy\"")
            
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
        self._dm_candidate = []
        self._coannihilation = []
        self._param_card = None
        self.coannihilation_diff = 1
        self._ID_procs = base_objects.ProcessDefinitionList()
        self._ID_matrix_elements = helas_objects.HelasMultiProcess()
        self._ID_amps = diagram_generation.AmplitudeList() 
        #self.dm = darkmatter.darkmatter()
        self._force_madevent_for_ID = False


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
                        raise DMError, '%s is not a valid particle for the model.' % args[1] 
                    if len(self._dm_candidate) == 1:
                        # No update of the model if 2(or more) DM since DD is not possible
                        self.update_model_with_EFT()
            
            elif args[0] == 'coannihilator':
                if not self._dm_candidate:
                    self.search_dm_candidate([])
                try:
                    args[1] = float(args[1])
                    isdigit=True
                except:
                    isdigit=False
                    pass
                
                if len(args) == 1:
                    self.search_coannihilator()                    
                elif '--usepdg' in args:
                    self._coannihilation = [self._curr_model.get_particle(a) for a in args[1:] if a!= '--usepdg']
                elif len(args)>2 and isdigit and args[2].startswith('/'):
                    self.search_coannihilator(gap=args[1], excluded=[a.replace('/', '') for a in args[2:]])
                elif len(args)>1 and not isdigit and args[1].startswith('/'):
                    self.search_coannihilator( excluded=[a.replace('/', '') for a in args[1:]])
                elif len(args)==2 and isdigit:
                    self.search_coannihilator(gap=args[1])
                else:
                    self._coannihilation = [self._curr_model.get_particle(a) for a in args[1:]]
                #avoid duplication
                if None in self._coannihilation:
                    raise self.InvalidCmd('Some of the particle name are invalid. Please retry.')
                all_name = [c.get('name') for c in self._coannihilation]
                self._coannihilation = [c for i,c in enumerate(self._coannihilation) if c.get('name') not in all_name[:i]]
#                self._dm_candidate += self._coannihilation 
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
            raise DMError, "No dark matter candidates in the model!"
        if dm_mass == 0:
            raise DMError, "No dark matter candidate since we found a DarkEnergy candidate: %s" % dm_particles[0]['name']

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

        if not dm_particles:
            logger.warning("""found %s but none have them have non zero coupling. Retry by excluding those""", ','.join(dm_names))
            return self.search_dm_candidate(excluded_particles=excluded_particles+dm_names)
        
        if len(dm_particles) >1:
            choice = self.ask("More than one valid candidate found: Please select the one that you want:",
                     default=dm_names[0], choices=dm_names)
            dm_particles = [p for p in dm_particles if p.get('name') == choice]
             
        # Print out the DM candidate
        logger.info("Found Dark Matter candidate: %s" % dm_particles[0]['name'],  '$MG:BOLD')
        self._dm_candidate = dm_particles
        self.update_model_with_EFT()
        
    def define_benchmark(self, path=None, answer=False):
     
     
        self._dm_candidate = []
        self._coannihilation = []
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
        else:
            path = answer
            
        if not isinstance(self._curr_model, model_reader.ModelReader):
            self._curr_model = model_reader.ModelReader(self._curr_model)        
        self._curr_model.set_parameters_and_couplings(path)   
        self._param_card = check_param_card.ParamCard(path)
        #self._dm_candidate = darkmatter.darkmatter()
    
      
    def search_coannihilator(self, gap=0.1, excluded=[]):
        """  This routine finds the coannihilation particles for the relic
        #  density calculation. The code can search for all the BSM particles that are
        #  within an input mass difference with the DM candidate.  All
        #  coannihilation particles are then added to the self._dm_particles
        #  list.
        """

        self._coannihilation = []
        dm_mass = self._curr_model.get_mass(self._dm_candidate[0])
        self.coannihilation_diff = gap 
        dm_name = [dm['name'] for dm in self._dm_candidate]
        
        bsm_particles = [p for p in self._curr_model.get('particles')                
                         if 25 < p.get('pdg_code') < 999000000 and\
                          (p.get('name') not in excluded or 
                          p.get('antiname') not in excluded or
                          str(p.get('pdgcode')) not in excluded)] 

        # Loop over BSM particles
        for p in bsm_particles:
            if p['name'] in dm_name:
                continue
            
            bsm_mass = self._curr_model.get_mass(p)
            if (abs(dm_mass-bsm_mass)/dm_mass <= gap):
                self._coannihilation.append(p)
            # If there are BSM particles that are too small to be included in the coannihilation they
            # are still tabulated to include in the final state particles with the SM particles.
            #elif ((self._bsm_masses[i] < (1.0-self._coann_eps)*self._dm_mass) and \
            #    not ('all' in self._excluded_particles) and \
            #    not (self._particles[self._bsm_particles[i]] in self._excluded_particles)):
            #          self._bsm_final_states.append(self._particles[self._bsm_particles[i]])

            # For organizational purposes we put the coannihilation particles by
            #alphabetical order by name
            self._coannihilation.sort(key=lambda p: p['name'])
        
        # Validity check
        for p in list(self._coannihilation):
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
                self._coannihilation.remove(p)
        
        
        if self._coannihilation:
            logger.info("Found coannihilation partners: %s" % ','.join([p['name'] for p in self._coannihilation]),
                    '$MG:BOLD')
        else:
            logger.info("No coannihilation partners found.", '$MG:BOLD')
        


    def do_import(self, line,*args, **opts):
        """normal import but perform a cleaning for MadDM  variables"""
        
        lineargs = self.split_arg(line)
        self.check_import(lineargs)
        if lineargs and lineargs[0].startswith('model'):
            #reset DM information
            self._param_card = None
            self._dm_candidate = []
            self._coannihilation = []
            
        return super(MadDM_interface, self).do_import(line, *args, **opts)

    def clean_process(self):
        """ensure that all processes are cleaned from memory.
        typically called from import model and generate XXX command
        """

        super(MadDM_interface, self).clean_process()
        
        # flip indirect/standard process definition
        self._ID_procs = base_objects.ProcessDefinitionList()
        self._ID_matrix_elements = helas_objects.HelasMultiProcess()
        self._ID_amps = diagram_generation.AmplitudeList()
        
    def check_save(self, args):
        
        if args[1] == 'maddm_first_indirect':
            args.insert(1, pjoin(MG5DIR, 'input', 'mg5_configuration.txt'))
            print args
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
        logger.info("**************** MG5AMC OPTION ***************************")
        super(MadDM_interface, self).help_define()

    help_add = help_generate
    
    def do_add(self, line):
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

        else:
            if '@' in line:
                line = re.sub(r'''(?<=@)(%s\b)''' % '\\b|'.join(self.process_tag), 
                              lambda x: `self.process_tag[x.group(0)]`, line)
            return super(MadDM_interface, self).do_add(line)
        
    @misc.mute_logger(['madgraph', 'models'], [30,30])    
    def add_model(self,*args,**opts):
        super(MadDM_interface, self).add_model(*args, **opts)

    def complete_generate(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddm options"""

        args = self.split_arg(line[:begidx])
        out = super(MadDM_interface, self).complete_generate(text, line, begidx, endidx,formatting=False)
        if not isinstance(out, dict):
            out = {"standard options": out}
        
        if len(args) == 1:
            options = ['relic_density', 'direct_detection', 'indirect_detection']#, 'capture']
            out['maddm options'] = self.list_completion(text, options , line)
        return self.deal_multiple_categories(out, formatting)

    def complete_add(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddm options"""

        args = self.split_arg(line[:begidx])
        out = super(MadDM_interface, self).complete_add(text, line, begidx, endidx,formatting=False)
        if not isinstance(out, dict):
            out = {"standard options": out}
        
        if len(args) == 1:
            options = ['relic_density', 'direct_detection', 'indirect_detection']#, 'capture']
            out['maddm options'] = self.list_completion(text, options , line)
        return self.deal_multiple_categories(out, formatting)


    def do_tutorial(self, line):
        """Activate/deactivate the tutorial mode."""

        if line:
            super(MadDM_interface, self).do_tutorial(line)
            if not line.lower().strip() not in ['off', 'quit']:
                return 
            
            
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
            self.history.append('generate relic_density')
            self.history.append('add direct_detection')            
            self.history.append('add indirect_detection')
        
        
        args = self.split_arg(line)
        if not args or args[0] not in self._export_formats + ['maddm']:
            if self._curr_amps:
                if any(amp.get('process').get('id') in self.process_tag.values()
                        for amp in self._curr_amps):
                    args.insert(0, 'maddm')

            
        if args and args[0] == 'maddm':
            line = ' '.join(args)
        
        if self._curr_amps:
            super(MadDM_interface, self).do_output(line)
        
        if self._ID_procs:

            path = self._done_export[0]
            if self._ID_procs!='2to2lo':
                import aloha.aloha_lib as aloha_lib
                aloha_lib.KERNEL = aloha_lib.Computation()

                with misc.TMP_variable(self,
                    ['_curr_proc_defs', '_curr_matrix_elements', '_curr_amps', '_done_export'],
                    [self._ID_procs, self._ID_matrix_elements, self._ID_amps, None]):
                    super(MadDM_interface, self).do_output('indirect %s/Indirect' % path)
                    
                #ensure to sync the param_card
                os.remove(pjoin(path, 'Indirect', 'Cards', 'param_card.dat'))
                files.ln(pjoin(path, 'Cards', 'param_card.dat'), 
                         pjoin(path, 'Indirect', 'Cards'))

                #ensure to sync the param_card
                os.remove(pjoin(path, 'Indirect', 'Cards', 'param_card.dat'))
                files.ln(pjoin(path, 'Cards', 'param_card.dat'), 
                     pjoin(path, 'Indirect', 'Cards'))

            import MGoutput
            proc_path = pjoin(path, 'matrix_elements', 'proc_characteristics')
            proc_charac = MGoutput.MADDMProcCharacteristic(proc_path)
            proc_charac['has_indirect_detection'] = True
            proc_charac.write(proc_path)
            

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
        

    def generate_relic(self, excluded_particles=[]):
        """--------------------------------------------------------------------#
        #                                                                      #
        #  This routine generates all the 2 -> 2 annihilation matrix elements. #
        #  The SM particles, BSM particles in self._bsm_final_states, and the  #
        #  DM particles are used as final states.  The initial states are      #
        #  looped over all the different combinations of DM particles.         #
        #--------------------------------------------------------------------"""

        #Here we define bsm multiparticle as a special case for exlusion only in the
        #final state. This is different than the standard / notation which will exclude
        #the particles from the propagators as well.
        #Since we use the exluded_particles everywhere, we just use the presence of
        #the bsm in the array to initialize a flag and then remove it from the excluded
        #particles array so as not to conflict with the use in other places.
        self.define_multiparticles('bsm',[p for p in self._curr_model.get('particles')\
                         if abs(p.get('pdg_code')) > 25])
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
            self.search_coannihilator(excluded=excluded_particles)
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
                         p.get('width') != 'ZERO' and
                         abs(p.get('pdg_code')) not in ids_veto]
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
        sm_pdgs = range(1, 7) + range(11, 17) + range(21, 26) #quarks/leptons/bosons
        sm_part = [self._curr_model.get_particle(pdg) for pdg in sm_pdgs]

        self.define_multiparticles('fs_particles', sm_part +bsm_final_states)

        # generate annihilation diagram
        coupling = "QED<=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0"
        if excluded_particles:
            proc = "dm_particles dm_particles > fs_particles fs_particles / %s %s  @DM2SM" \
                    % (' '.join(excluded_particles), coupling)
        else:
            proc = "dm_particles dm_particles > fs_particles fs_particles %s @DM2SM" \
                   %(coupling)  #changed here
                   
        self.do_add('process %s' %proc)

        # Generate all the DM -> DM processes
        nb_dm = len(self._dm_candidate + self._coannihilation)
        for i in xrange(nb_dm):
            for j in xrange(i,nb_dm):
                if  excluded_particles:
                    proc = "DM_particle_%s DM_particle_%s > dm_particles dm_particles / %s %s @DM2DM"\
                       % (i,j,' '.join(excluded_particles), coupling)
                else:
                    proc = "DM_particle_%s DM_particle_%s > dm_particles dm_particles %s @DM2DM"\
                       % (i,j, coupling)
                try:
                    self.do_add('process %s' % proc)
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
        for i in xrange(nb_dm):
            for j in xrange(nb_dm):
                if i == j:
                    continue
                if excluded_particles:
                    proc = "DM_particle_%s sm_particles > DM_particle_%s sm_particles / %s %s @DMSM"\
                            % (i,j,' '.join(excluded_particles), coupling)
                else:
                    proc = "DM_particle_%s sm_particles > DM_particle_%s sm_particles %s @DMSM"\
                            % (i,j, coupling)
                try:
                    self.do_add('process %s' % proc)
                except (self.InvalidCmd,diagram_generation.NoDiagramException) :
                    continue

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

        quarks = range(1,7) # ['d', 'u', 's', 'c', 'b','t']
        antiquarks = [-1*pdg for pdg in quarks] # ['d~', 'u~', 's~', 'c~', 'b~','t~']
        
        if type == 'SI':
            orders = '%s==2' % SI_name
        elif type == 'SD':
            orders = '%s==2' % SD_name
        elif type == 'QED':
            orders = '%s=0 %s=0' % (SD_name, SI_name)
        elif type == 'SI+QED':
            orders = '%s=0 %s<=99' % (SD_name, SI_name)
        elif type == 'SD+QED':
            orders = '%s<=99 %s=0' % (SD_name, SI_name)
            
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
                self.do_add('process %s' % proc)
            except (self.InvalidCmd, diagram_generation.NoDiagramException), error:
                logger.debug(error)
                continue # no diagram generated
            has_diagram = True
        return has_diagram

    def generate_indirect(self, argument, user=True):
        """User level function which performs indirect detection functions
           Generates the DM DM > X where X is anything user specified.
           Also works for loop induced processes as well as NLO-QCD.
           Currently works only with canonical mode and one dm candidate.
        related to syntax: generate indirect a g / n3
        """
        
        if not maddm_run_interface.HAS_NUMPY and user:
            logger.warning("numpy module not detected on your machine. \n"+
                           "Running indirect detection will not be possible as long as numpy is not installed (try 'pip install numpy')" 
                           )
            ans = self.ask("Do you want to generate the diagrams anyway?", 'n', ['y','n'])
            if ans == 'n':
                return
        
        self.install_indirect()

        
        if not self._dm_candidate:
            self.search_dm_candidate()
            if not self._dm_candidate:
                return
            else:
                self.history.append('add indirect_detection %s' % ' '.join(argument))

        if len(self._curr_proc_defs) == 0 or \
           all(p.get('id') != self.process_tag['DM2SM'] for p in self._curr_proc_defs):
            logger.info('For indirect detection we need to generate relic density matrix-element','$MG:BOLD')
            self.generate_relic([])

        if '2to2lo' in argument:
            self._ID_procs ='2to2lo'
            return
        
        if not self.options['pythia8_path']:
            logger.warning('''In order to have distribution related to indirect detection:
             Pythia8 needs to be linked to the code. You can install it by typing 'install pythia8'.
             You can compute the rate (no distribution) by typing 'add indirect 2to2lo'. ''')
        #    return

        #Check if the argument contains only two particles in the final state
        #if not, force the code to use madevent
        
        if (len(' '.join(argument).split('/')[0].split())!=2):
            self._force_madevent_for_ID = True
        # Generates events for all the DM annihilation channels (also BSM)
        if argument == []:

            dm_cand = [ dic['pdg_code'] for dic in self._dm_candidate ] 
            DM_scattering = []

            for pdg in dm_cand:
                DM_scattering.append(pdg) , DM_scattering.append(-pdg)

            bsm_all  = [ 'all', '/', '1','2','3','4','5','6','-1','-2','-3','-4','-5','-6','21','22','23','24','-24','25','11','12','13','14','15','16','-11','-12',
                    '-13','-14','-15','-16']

            for m in DM_scattering:
                bsm_all.append(m)
            bsm = " ".join(str(x) for x in bsm_all)

            self.exec_cmd('define q_mdm = 1 2 3 4 -1 -2 -3 -4',postcmd=False)
            self.exec_cmd('define bsm = %s'  % bsm,postcmd=False)

            final_states = ['bsm bsm','q_mdm q_mdm','21 21', '5 -5', '6 -6', '22 22', '23 23', '24 -24','25 25', '11 -11', '13 -13', '15 -15', '12 -12', '14 -14', '16 -16']

            for final_state in final_states:
                try: 
                    self.generate_indirect([final_state, '--noloop'], user=False)
                except diagram_generation.NoDiagramException:
                    continue
                    logger.info('no diagram for %s' % final_state)

            return

    
        allow_loop_induce=True
        if '--noloop' in argument:
            argument.remove('--noloop')
            allow_loop_induce=False

        # flip indirect/standard process definition
        with misc.TMP_variable(self, ['_curr_proc_defs', '_curr_matrix_elements', '_curr_amps'], 
                                     [self._ID_procs, self._ID_matrix_elements, self._ID_amps]):
            # First try LO matrix-element
            coupling = "SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0"
            done= []
            for dm in self._dm_candidate:
                name = dm.get('name')
                antiname = dm.get('antiname')
                if name in done:
                    continue
                done += [name, antiname]
                #We put the coupling order restrictions after the @ID in order to
                #apply it to the entire matrix element.
                proc = '%s %s > %s @ID %s' % (name, antiname, ' '.join(argument), coupling)
                try:
                    self.do_add('process %s' % proc)
                except (self.InvalidCmd, diagram_generation.NoDiagramException), error:
                    if allow_loop_induce:
                        proc = '%s %s > %s %s [noborn=QCD] @ID ' % (name, antiname, ' '.join(argument), coupling)
                        try:                    
                            self.do_add('process %s' % proc)
                        except (self.InvalidCmd, diagram_generation.NoDiagramException), error:
                            proc = '%s %s > %s %s [noborn=QED] @ID ' % (name, antiname, ' '.join(argument), coupling)
                            self.do_add('process %s' % proc)
  
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
        
        
            self.exec_cmd(mg5_command, postcmd=False)
            self._curr_amps = backup_amp
            self._param_card = backup_param_card 
            self._dm_candidate = backup_dm_candidate
            self._coannihilation = backup_coannihilation

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
        logger.info("      The code selects all the BSM particles that are within an input mass difference with the DM candidate.")
        logger.info("      Particles without any non vanishing coupling are discarded")
 
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

    
    

import logging
import os

import maddm_run_interface as maddm_run_interface

import madgraph.core.base_objects as base_objects
import madgraph.core.helas_objects as helas_objects
import madgraph.core.diagram_generation as diagram_generation
import madgraph.interface.master_interface as master_interface
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc
#import darkmatter as darkmatter
import madgraph.interface.common_run_interface as common_run
import models.check_param_card as check_param_card
import models.model_reader as model_reader

#HERE I HARDCODED NUMPY
import numpy as np

import re
pjoin = os.path.join

logger = logging.getLogger('madgraph.plugin.maddm')

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
  "            |                  "+bcolors.OKBLUE+"  MadDM v2.0                     "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "                                                                               \n"+\
  "                #########                                                        \n"+\
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
  "%s" 
    
    _define_options = ['darkmatter', 'coannihilator', 'benchmark']
    # process number to distinguish the different type of matrix element
    process_tag = {'DM2SM': 1999, 
                   'DM2DM': 1998,
                   'DMSM': 1997,
                   'DD': 1996,
                   'ID': 1995}
    
    eff_operators_SI = {1:'SIEFFS', 2:'SIEFFF', 3:'SIEFFV'}
    eff_operators_SD = {1:False, 2:'SDEFFF', 3:'SDEFFV'} 
    
    def preloop(self, *args, **opts):
        super(MadDM_interface, self).preloop(*args, **opts)
        self.prompt = 'MadDM>'
    
    def change_principal_cmd(self, name):
        out = super(MadDM_interface, self).change_principal_cmd(name)
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
                elif len(args)==2:
                    self._dm_candidate = [self._curr_model.get_particle(args[1])]
                    if not self._dm_candidate[0]:
                        raise DMError, '%s is not a valid particle for the model.' % args[1] 
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
        logger.info("    If no option specified. It assigned the less massive neutral BSM particle.")
        logger.info(" OPTIONS:", '$MG:color:BLACK')
        logger.info("    - You can specify the darkmatter by specifying their name/pdg code")
        logger.info("     o example: define darkmatter n1 n2",'$MG:color:GREEN')
        logger.info("    - You can remove some particle from the search by prefixing the particle name by '/'.")
        logger.info("     o example: define darkmatter /n1", '$MG:color:GREEN')                    
        logger.info("")        
        logger.info("syntax: define coannihilator [REL_DIFF | PARTICLE(S)] [/ excluded_particles]", '$MG:color:BLUE')
        logger.info(" -- define the current (sets) of coannihilator. ")
        logger.info("    If no option is provided use a  relative difference of 10% is used.")         
        logger.info(" OPTIONS:", '$MG:color:BLACK')
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
        #                                                                        #
        #-----------------------------------------------------------------------#
        particles = self._curr_model.get('particles')                
        bsm_particles = [p for p in particles 
                         if 25 < p.get('pdg_code') < 99000000 and\
                         (p.get('name') not in excluded_particles or p.get('antiname') not in excluded_particles) and\
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
                if not dm_particles:
                    dm_particles.append(p)
                    dm_mass = curr_mass
                # If we already found a candidate, comare the masses and keep the one that's lighter
                elif (curr_mass < dm_mass):
                    dm_particles[0] = p
                    dm_mass = curr_mass

        # Check to see if we actually found a DM candidate
        if not dm_particles:
            raise DMError, "No dark matter candidates in the model!"
        if dm_mass == 0:
            raise DMError, "No dark matter candidate since we found a DarkEnergy candidate: %s" % dm_particles[0]['name']

        # Print out the DM candidate
        logger.info("Found Dark Matter candidate: %s" % dm_particles[0]['name'],  '$MG:color:BLACK')
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
                                  ask_class = common_run.EditParamCard,
                                  card=[path])
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
        
        if self._coannihilation:
            logger.info("Found coannihilation partners: %s" % ','.join([p['name'] for p in self._coannihilation]),
                    '$MG:color:BLACK')
        else:
            logger.info("No coannihilation partners found.", '$MG:color:BLACK')
        


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
      

    def do_add(self, line):
        """ """
        
        args = self.split_arg(line)    
        if len(args) and args[0] == 'process':
            args.pop(0)
        if len(args) and args[0] in["relic_density", 'relic']:
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
        elif len(args) and args[0] in ["earth_capture", "earth"]:
            self.generate_EarthCapture()
        elif len(args) and args[0] in ["sun_capture", "sun"]:
            self.generate_SolarCapture()

        else:
            if '@' in line:
                line = re.sub(r'''(?<=@)(%s\b)''' % '\\b|'.join(self.process_tag), 
                              lambda x: `self.process_tag[x.group(0)]`, line)
            return super(MadDM_interface, self).do_add(line)
            

    def complete_generate(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddm options"""

        args = self.split_arg(line[:begidx])
        out = super(MadDM_interface, self).complete_generate(text, line, begidx, endidx,formatting=False)
        if not isinstance(out, dict):
            out = {"standard options": out}
        
        if len(args) == 1:
            options = ['relic_density', 'direct_detection', 'indirect_detection', 'sun_capture', 'earth_capture']
            out['maddm options'] = self.list_completion(text, options , line)
        return self.deal_multiple_categories(out, formatting)

    def complete_add(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddm options"""

        args = self.split_arg(line[:begidx])
        out = super(MadDM_interface, self).complete_add(text, line, begidx, endidx,formatting=False)
        if not isinstance(out, dict):
            out = {"standard options": out}
        
        if len(args) == 1:
            options = ['relic_density', 'direct_detection', 'indirect_detection', 'sun_capture', 'earth_capture']
            out['maddm options'] = self.list_completion(text, options , line)
        return self.deal_multiple_categories(out, formatting)



    def do_output(self, line):
        """ """
        
        if not self._curr_amps:
            self.do_generate('relic_density')
            self.do_add('direct_detection')
        
        
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
                    super(MadDM_interface, self).do_output('madevent_maddm %s/Indirect' % path)

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

    def do_launch(self, line):
        
        
        args = self.split_arg(line)
        (options, args) = madgraph_interface._launch_parser.parse_args(args)
        self.check_launch(args, options)
        options = options.__dict__        
        
        
        if args[0] != 'maddm':
            return super(MadDM_interface, self).do_launch(line)
        else:
            self._MDM = maddm_run_interface.MADDMRunCmd(dir_path=args[1], options=self.options)
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
                   
        self.do_generate(proc)

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

        # for i in xrange(nb_dm):
        #     for j in xrange(nb_dm):
        #         if i == j:
        #             continue
        #         if excluded_particles:
        #             proc = "DM_particle_%s sm_particles > DM_particle_%s sm_particles / %s %s @DMSM"\
        #                % (i,j,' '.join(excluded_particles), coupling)
        #         else:
        #             proc = "DM_particle_%s sm_particles > DM_particle_%s sm_particles %s @DMSM"\
        #                % (i,j, coupling)
        #         try:
        #             self.do_add('process %s' % proc)
        #         except (self.InvalidCmd,diagram_generation.NoDiagramException) :
        #             continue

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

        if len(self._dm_candidate) > 1:
            logger.warning("More than one DM candidate. Can not run Direct Detection.")
            return 
  
   
        
        #Now figure out the label of the effective vertex to use. The convention is:
        #<SI or SD>EFF<F, S, or V>, i.e. SIEFFV for Spin Independent Effective vertex for Vector Current.
        dm_spin = int(self._dm_candidate[0]['spin'])
        eff_operators_SI = self.eff_operators_SI[dm_spin]
        eff_operators_SD = self.eff_operators_SD[dm_spin]
        
        logger.info("Generating X Nucleon > X Nucleon diagrams from the full lagrangian...")
        has_direct = self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'QED')

        if not has_direct:
            logger.warning("No Direct Detection Feynman Diagram")
            return
        
        logger.info("Generating X Nucleon > X Nucleon diagrams from the effective lagrangian...")
        #ONLY EFFECTIVE LAGRANGIAN
        self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SI')

        logger.info("INFO: Generating X Nucleon > X Nucleon diagrams from the effective+full lagrangian...")
        #EFFECTIVE + FULL
        self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SI+QED')
        
        if (eff_operators_SD != False):
            logger.info("Doing the spin dependent part...")
            logger.info("Generating X Nucleon > X Nucleon diagrams from the effective lagrangian...")

            self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SD')
            #EFFECTIVE + FULL
            logger.info("Generating X Nucleon > X Nucleon diagrams from the effective + full lagrangian...")
            self.DiagramsDD(eff_operators_SI, eff_operators_SD, 'SD+QED')

    #-----------------------------------------------------------------------#
    def DiagramsDD(self, SI_name, SD_name, type, excluded=[]):
        """Generates direct detection diagrams. i_dm is the index of DM part. 
                 Whether spin dependent or spin independent diagrams should be                 
                 calculated. XX_order parameters determine maximum order of  a
                 coupling. For instance, if you only want effective vertices for
                 spin independent coupling you would set SI_order = 2 and all others
                 to zero. If you want the spin independent full lagrangian + eff.
                 then you need to set SI_order=2 and QED_order=2...
                 WARNING: This function is a proxy to be used inside
                 GenerateDiagramsDirDetect() function and no place else!
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

    def generate_indirect(self, argument):
        """User level function which performs indirect detection functions
           Generates the DM DM > X where X is anything user specified.
           Also works for loop induced processes as well as NLO-QCD.
           Currently works only with canonical mode and one dm candidate.
        related to syntax: generate indirect a g / n3
        """

        if '2to2lo' in argument:
            self._ID_procs ='2to2lo'
            return

        #Check if the argument contains only two particles in the final state
        #if not, force the code to use madevent
        if (len(' '.join(argument).split('/')[0].split())!=2):
            self._force_madevent_for_ID = True

        if not self._dm_candidate:
            self.search_dm_candidate()
            if not self._dm_candidate:
                return
            
        if len(self._curr_proc_defs) == 0:
            self.generate_relic([])

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
                    proc = '%s %s > %s %s [noborn=QCD] @ID ' % (name, antiname, ' '.join(argument), coupling)
                    self.do_add('process %s' % proc)
  
    def update_model_with_EFT(self):
        """ """
        eff_operators_SD = {1:False, 2:'SDEFFF', 3:'SDEFFV'}
        eff_model_dm_names = {1:'~sdm', 2:'~fdm', 3:'~vdm'}

        DM = self._dm_candidate[0]
        
        if self._dm_candidate[0]['self_antipart']: 
            EFT = 'REAL'
        else:
            EFT = 'COMPLEX'

        eff_dm_name = eff_model_dm_names[int(DM['spin'])]
        mg5_command = 'add model %s %s=%s --recreate --keep_decay' % (pjoin(MDMDIR, 'EffOperators', EFT),
                                                     eff_dm_name,DM.get('name'))
        
        # We want to preserve the following variable while updating the model
        backup_amp = self._curr_amps
        backup_param_card = self._param_card
        backup_dm_candidate = self._dm_candidate
        backup_coannihilation = self._coannihilation
        
        self.exec_cmd(mg5_command)
        self._curr_amps = backup_amp
        self._param_card = backup_param_card 
        self._dm_candidate = backup_dm_candidate
        self._coannihilation = backup_coannihilation

        # update the param_card value
        txt = self._curr_model.write_param_card()
        param_card = check_param_card.ParamCard(self._curr_model.write_param_card())

        if self._param_card:
            for block in self._param_card:
                for param in self._param_card[block]:
                    param_card[block].get(param.lhacode).value =  param.value
 
        self._param_card = param_card        
        if not isinstance(self._curr_model, model_reader.ModelReader):
            self._curr_model = model_reader.ModelReader(self._curr_model) 
        self._curr_model.set_parameters_and_couplings(self._param_card) 
        
        
 

        
        
        

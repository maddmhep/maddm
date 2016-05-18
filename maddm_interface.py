import logging
import os

import maddm_run_interface as maddm_run_interface


import madgraph.interface.master_interface as master_interface
import madgraph.interface.madgraph_interface as madgraph_interface
import madgraph.various.misc as misc
#import darkmatter as darkmatter
import madgraph.interface.common_run_interface as common_run
import models.check_param_card as check_param_card
import models.model_reader as model_reader


import re
pjoin = os.path.join

logger = logging.getLogger('madgraph.plugin.maddm')

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



    

class MadDM_interface(master_interface.MasterCmd):

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
  "%s" 
    
    _define_options = ['darkmatter', 'coannihilator', 'benchmark']
    
    # process number to distinguish the different type of matrix element
    process_tag = {'DM2SM': 1999, 
                   'DM2DM': 1998,
                   'DMSM': 1997}
    
    
    def __init__(self, *args, **opts):
        
        super(MadDM_interface, self).__init__(*args, **opts)
        self._dm_candidate = []
        self._coannihilation = []
        self._param_card = None
        self.coannihilation_diff = 1
        #self.dm = darkmatter.darkmatter()
        


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
                    self._dm_candidate = self._curr_model.get_particle(args[1])
            
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
                         if p.get('pdg_code') > 25 and\
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
            if p['width'] == 'ZERO' or  get_width(p) == 0:
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
            raise DMError, "No dark matter candidate since we found a DarkEnergy candidate"

        # Print out the DM candidate
        logger.info("Found Dark Matter candidate: %s" % dm_particles[0]['name'],  '$MG:color:BLACK')
        self._dm_candidate = dm_particles
        
    def define_benchmark(self, path=None, answer=False):
     
     
        self._dm_candidate = []
        self._coannihilation = None
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
        
        
        bsm_particles = [p for p in self._curr_model.get('particles')                
                         if p.get('pdg_code') > 25 and\
                         (p.get('name') not in excluded or 
                          p.get('antiname') not in excluded or
                          str(p.get('pdgcode'))) not in excluded] 

        # Loop over BSM particles
        for p in bsm_particles:
            if p in self._dm_candidate:
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
        misc.sprint(args)
        if len(args) >=2 and args[1] == "relic_density":
            if '/' not in line:
                return self.generate_relic([])
            else:
                subline = line.split('/',1)[1]
                excluded = [ a for a in self.split_arg(subline) if not a.startswith('-')]
                return self.generate_relic(excluded)
        else:
            if '@' in line:
                misc.sprint(line)
                line = re.sub(r'''(?<=@)(%s\b)''' % '\\b|'.join(self.process_tag), 
                              lambda x: `self.process_tag[x.group(0)]`, line)
                misc.sprint(line)
            return super(MadDM_interface, self).do_add(line)
            

    def complete_generate(self, text, line, begidx, endidx, formatting=True):
        """Complete particle information + handle maddm options"""

        args = self.split_arg(line[:begidx])
        out = super(MadDM_interface, self).complete_generate(text, line, begidx, endidx,formatting=False)
        if not isinstance(out, dict):
            out = {"standard options": out}
        
        if len(args) == 1:
            options = ['relic_density']
            out['maddm options'] = self.list_completion(text, options , line)
        return self.deal_multiple_categories(out, formatting)



    def do_output(self, line):
        """ """
        
        if not self._curr_amps:
            self.do_generate('relic_density')
        
        args = self.split_arg(line)
        if not args or args[0] not in self._export_formats + ['maddm']:
            if self._curr_amps:
                if any(amp.get('process').get('id') in self.process_tag.values()
                        for amp in self._curr_amps):
                    args.insert(0, 'maddm')

            
        if args and args[0] == 'maddm':
            line = ' '.join(args)
        
        super(MadDM_interface, self).do_output(line)

    def find_output_type(self, path):
        
        if os.path.exists(pjoin(path,'matrix_elements','proc_characteristics')):
            return 'maddm'
        else:
            return super(MadDM_interface, self).find_output_type(self, path)

    def do_launch(self, line):
        
        
        args = self.split_arg(line)
        (options, args) = madgraph_interface._launch_parser.parse_args(args)
        self.check_launch(args, options)
        options = options.__dict__        
        
        
        if args[0] != 'maddm':
            return super(MadDM_interface, self).do_launch(line)
        else:
            MDM = maddm_run_interface.MADDMRunCmd(dir_path=args[1], options=self.options) 
            if options['interactive']:              
                return self.define_child_cmd_interface(MDM)
            else:
                self.define_child_cmd_interface(MDM,  interface=False)
                MDM.exec_cmd('launch ' + line.replace(args[1], ''))
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


        if not self._dm_candidate:
            self.search_dm_candidate(excluded_particles)
            if not self._dm_candidate:
                return
            self.search_coannihilator(excluded=excluded_particles)

        # Tabulates all the BSM particles so they can be included in the
        # final state particles along with the SM particles.
        
        dm_mass = abs(self._curr_model.get_mass(self._dm_candidate[0]))
        is_lower_mass = lambda p : abs(self._curr_model.get_mass(p)) < (1-self.coannihilation_diff)*dm_mass
                
        bsm_final_states = [p for p in self._curr_model.get('particles') 
                         if abs(p.get('pdg_code')) > 25 and \
                         (p.get('name') not in excluded_particles or p.get('antiname') not in excluded_particles) and\
                         is_lower_mass(p) and
                         p not in self._dm_candidate and 
                         p not in self._coannihilation]
        
        misc.sprint([(p.get('name'), self._curr_model.get_mass(p)) for p in bsm_final_states])
        
        # Set up the initial state multiparticles that contain the particle 
        #and antiparticle
        for i,dm in enumerate(self._dm_candidate + self._coannihilation):
            misc.sprint(i)
            self.define_multiparticles('dm_particle_%s' % i,  [dm])
            self.do_display('multiparticles')
        self.define_multiparticles('dm_particles', self._dm_candidate)
        self.do_display('multiparticles')
        sm_pdgs = range(1, 7) + range(11, 17) + range(21, 26) #quarks/leptons/bosons
        sm_part = [self._curr_model.get_particle(pdg) for pdg in sm_pdgs]
        self.define_multiparticles('fs_particles', sm_part+bsm_final_states)
        
        
        logger.info("DM is allowed to annihilate into the following BSM particles:\n %s", ', '.join(p['name'] for p in bsm_final_states))

        # generate annihilation diagram
        coupling = "QED<=4 SIEFFS=0 SIEFFF=0 SIEFFV=0 SDEFFF=0 SDEFFV=0"
        if excluded_particles:
            proc = "dm_particles dm_particles > fs_particles fs_particles / %s %s  @DM2SM" \
                    % (' '.join(excluded_particles), coupling)
        else:
            proc = "dm_particles dm_particles > fs_particles fs_particles %s @DM2SM" \
                   %(coupling)
                   
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
                misc.sprint(proc)
                self.do_add('process %s' % proc)
        
        # Get the matrix elements and make sure that we don't have any pure 
        #scattering processes
        for amp in self._curr_amps[:]:
            if amp.get('process').get('id') != self.process_tag['DM2DM']:
                continue
            if (set(amp.get('process').get_initial_ids()) == (set(amp.get('process').get_final_ids()))):
                self._curr_amps.remove(amp)


        # Generate all the DM particles scattering off of the thermal background and
        # change to a different DM species (again, no pure scatterings)
        for i in xrange(nb_dm):
            for j in xrange(nb_dm):
                if i == j:
                    continue
                if excluded_particles:
                    proc = "DM_particle_%s fs_particles > DM_particle_%s fs_particles / %s %s @DMSM"\
                       % (i,j,' '.join(excluded_particles), coupling)
                else:
                    proc = "DM_particle_%s fs_particles > DM_particle_%s fs_particles %s @DMSM"\
                       % (i,j, coupling)
                self.do_add('process %s' % proc)



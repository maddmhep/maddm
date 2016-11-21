from __future__ import division
import collections
import math
import logging
import os
import re
import sys
import subprocess


import MGoutput

MG5MODE = True
import madgraph.various.misc as misc
import madgraph.interface.extended_cmd as cmd
import madgraph.various.banner as banner_mod
import madgraph.interface.common_run_interface as common_run
import madgraph.interface.madevent_interface as me5_interface
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

#Conversion from GeV^-2 to pb etc.
GeV2pb = 3.894E8
pb2cm2  = 1.0E-36
cm22pb  = 1.0E36

#===============================================================================
# CommonRunCmd
#===============================================================================
class MADDMRunCmd(cmd.CmdShell):
    
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
  "" 
    
    
    def __init__(self, dir_path, options, *args, **opts):
        """common"""
  
        self.in_scan_mode = False # please modify this value only with "with" statement
        cmd.Cmd.__init__(self, *args, **opts)
        # Define current MadEvent directory
        if dir_path is None and not MG5MODE:
            dir_path = root_path

        self.dir_path = dir_path
        self.param_card_iterator = [] #an placeholder containing a generator of paramcard for scanning
        
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

    def check_param_card(self, path, run=True):
        """
        1) Check that no scan parameter are present
        2) Check that all the width are define in the param_card.
        - If a scan parameter is define. create the iterator and recall this fonction 
          on the first element.
        - If some width are set on 'Auto', call the computation tools."""
        
        pattern_scan = re.compile(r'''^(decay)?[\s\d]*scan''', re.I+re.M)  
        pattern_width = re.compile(r'''decay\s+(\+?\-?\d+)\s+auto(@NLO|)''',re.I)
        text = open(path).read()
               
        if pattern_scan.search(text):
            if not isinstance(self, cmd.CmdShell):
                # we are in web mode => forbid scan due to security risk
                raise Exception, "Scan are not allowed in web mode"
            # at least one scan parameter found. create an iterator to go trough the cards
            main_card = param_card_mod.ParamCardIterator(text)
            self.param_card_iterator = main_card
            first_card = main_card.next(autostart=True)
            first_card.write(path)
            return self.check_param_card(path, run)
        
        pdg_info = pattern_width.findall(text)
        if pdg_info:
            if run:
                if not self.in_scan_mode:
                    logger.info('Computing the width set on auto in the param_card.dat')
                has_nlo = any(nlo.lower()=="@nlo" for _,nlo in pdg_info)
                pdg = [pdg for pdg,nlo in pdg_info]
                if not has_nlo:
                    self.do_compute_widths('%s --path=%s' % (' '.join(pdg), path))
                    #self.run_mg5(['compute_widths %s --path=%s' % (' '.join(pdg), path)])
                else:
                    self.run_mg5(['compute_widths %s --path=%s --nlo' % (' '.join(pdg), path)])
            else:
                logger.info('''Some width are on Auto in the card. 
    Those will be computed as soon as you have finish the edition of the cards.
    If you want to force the computation right now and being able to re-edit
    the cards afterwards, you can type \"compute_wdiths\".''')


    def do_compute_widths(self, line):
                
        if self.maddm_card['only_two_body_decays']:
            line = ' --body_decay=2 ' + line
        #return self.run_mg5(['compute_widths --body_decay=2 ' + line])
        if not hasattr(self, 'mg5'):
            return self.run_mg5(['compute_widths ' + line])
        elif not hasattr(self.mg5, '_curr_model'):
            return self.run_mg5(['compute_widths ' + line])
        else:
            self.mg5.do_compute_widths(line, model=self.mg5._curr_model)#, decaymodel=self.mg5._curr_decaymodel)
    
    def help_compute_widths(self, line):
        
        return self.run_mg5([' help compute_widths ' + line])    


    ############################################################################
    def run_mg5(self, commands, model=None):
        """ """
        
        if not hasattr(self, 'mg5'):
            if self.mother:
                self.mg5 = self.mother
            else:
                import madgraph.interface.master_interface as master_interface
                self.mg5 = master_interface.MasterCmd()

        if not model:
            model = self.proc_characteristics['model']
        if self.mg5._curr_model.get('modelpath+restriction') != model:
            self.mg5.do_import(model)
            
            
        for line in commands:
            self.mg5.exec_cmd(line, errorhandling=False, printcmd=False, 
                              precmd=False, postcmd=False, child=False)
            
    ############################################################################
    def do_launch(self, line):
        """run the code"""

        args = line.split()
        if '-f' in args or '--force' in args:
            force = True
        else:
            force = False
        self.ask_run_configuration(mode=[], force=force)
        self.compile()
        
        nb_output = 1
        output = pjoin(self.dir_path, 'output', 'maddm_%s.out') 
        while os.path.exists(output % nb_output):
            nb_output +=1
        else:
            output = pjoin('output', 'maddm.out')
        
        misc.call(['./maddm.x', output], cwd =self.dir_path)
        #process = subprocess.Popen(['./maddm.x'], cwd =self.dir_path, stdout=subprocess.PIPE)
        #Here we read out the results which the FORTRAN module dumped into a file
        #called 'maddm.out'. The format is such that the first line is always relic density
        # , second line is the nucleon scattering cross section (SI) for proton, third is SI
        # nucleon cross section for the neutron etc. If a quantity was not calculated, we output -1        
        result = []
        #for line in process.stdout:
        for line in open(pjoin(self.dir_path, output)):
            result.append(float(line.split()[1]))

        self.last_results = result
        output_name = ('omegah2', 'x_freezeout', 'wimp_mass', 'sigmav_xf' , 
                       'sigmaN_SI_proton', 'sigmaN_SI_neutron', 'sigmaN_SD_proton',
                        'sigmaN_SD_neutron','Nevents', 'smearing')
        result = dict(zip(output_name, result))
        self.last_results = result

        if self.mode['indirect']:
            with misc.MuteLogger(names=['madevent','madgraph'],levels=[50,50]):
                self.launch_indirect(force)


        if not self.in_scan_mode:
            self.print_results()
        #else:
        #    logger.info("relic density  : %.2e ", self.last_results['omegah2'])
        
        
        if self.param_card_iterator:
            param_card_iterator = self.param_card_iterator
            self.param_card_iterator = []
            param_card_iterator.store_entry(nb_output, result)
            #print results:
            order = []
            if self.mode['relic']:
                order += ['omegah2', 'x_freezeout', 'sigmav_xf']
            if self.mode['direct'] :
                order += ['sigmaN_SI_proton', 'sigmaN_SI_neutron', 'sigmaN_SD_proton',
                        'sigmaN_SD_neutron']                
            if self.mode['directional']:
                order += ['Nevents', 'smearing']
            if self.mode['indirect']:
                for i in range(len(self.maddm_card['halo_dm_velocity'])):
                    order +=['halo_velocity#%s' %i,'indirect#%s' %i, 'indirect_error#%s' %i] 
            to_print = param_card_iterator.write_summary(None, order,nbcol=10)
            for line in to_print.split('\n'):
                if line:
                    logger.info(line)

            #check if the param_card defines a scan.
            with misc.TMP_variable(self, 'in_scan_mode', True):
                with misc.MuteLogger(names=['cmdprint','madevent','madgraph','madgraph.plugin'],levels=[50,50,50,20]):
                    for i,card in enumerate(param_card_iterator):
                        card.write(pjoin(self.dir_path,'Cards','param_card.dat'))
                        self.exec_cmd("launch -f", precmd=True, postcmd=True,
                                                   errorhandling=False)
                        param_card_iterator.store_entry(nb_output+i, self.last_results)
                        logger.info(param_card_iterator.write_summary(None, order, lastline=True,nbcol=5)[:-1])
            param_card_iterator.write(pjoin(self.dir_path,'Cards','param_card.dat'))
            name = misc.get_scan_name('maddm_%s' % (nb_output), 'maddm_%s' % (nb_output+i))
            path = pjoin(self.dir_path, 'output','scan_%s.txt' % name)
            logger.info("write all results in %s" % path ,'$MG:color:BLACK')

            param_card_iterator.write_summary(path, order)
    
        return result

    def launch_indirect(self, force):
        """running the indirect detection"""

        if not self.in_scan_mode: 
            logger.info('Running indirect detection')
        if not hasattr(self, 'me_cmd'):
            self.me_cmd = Indirect_Cmd(pjoin(self.dir_path, 'Indirect'))
        elif self.me_cmd.me_dir != pjoin(self.dir_path, 'Indirect'):
            self.me_cmd.do_quit()
            self.me_cmd = Indirect_Cmd(pjoin(self.dir_path, 'Indirect'))
            
        runcardpath = pjoin(self.dir_path,'Indirect', 'Cards', 'run_card.dat')
        param_path  = pjoin(self.dir_path,'Indirect', 'Cards', 'param_card.dat')
        run_card = banner_mod.RunCard(runcardpath)
        param_card = param_card_mod.ParamCard(param_path)
        mdm = param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        
        for i,v in enumerate(self.maddm_card['halo_dm_velocity']):
            run_card['ebeam1'] = mdm * math.sqrt(1+v**2)
            run_card['ebeam2'] = mdm * math.sqrt(1+v**2)
            run_card.write(runcardpath)
            
            self. me_cmd.do_launch('-f')
            
            self.last_results['halo_velocity#%s' %i] = v
            self.last_results['indirect#%s' %i] = self.me_cmd.results.get_detail('cross')
            self.last_results['indirect_error#%s' %i] = self.me_cmd.results.get_detail('error')
        
    def print_results(self):
        """ print the latest results """
        omega_min = 0
        omega_max = 0.1
        
        if omega_min < self.last_results['omegah2'] < omega_max:
            fail_relic_msg = ''
        else:
            fail_relic_msg = '%s Model excluded (relic not in range [%s,%s])%s' %\
                              (bcolors.FAIL, omega_min, omega_max,bcolors.ENDC)  
        
        
        logger.info("*** RESULTS ***", '$MG:color:BLACK')
        if self.mode['relic']:
            logger.info('   relic density  : %.2e %s', self.last_results['omegah2'],fail_relic_msg)
                        
            logger.info('   x_f            : %.2f', self.last_results['x_freezeout'])
            logger.info('   sigmav(xf)     : %.2e GeV^-2 = %.2e pb', self.last_results['sigmav_xf'],self.last_results['sigmav_xf']*GeV2pb)
            
        if self.mode['direct']:
            logger.info(' sigmaN_SI_p      : %.2e GeV^-2 = %.2e pb',self.last_results['sigmaN_SI_proton'],self.last_results['sigmaN_SI_proton']*GeV2pb)
            logger.info(' sigmaN_SI_n      : %.2e GeV^-2 = %.2e pb',self.last_results['sigmaN_SI_neutron'],self.last_results['sigmaN_SI_neutron']*GeV2pb)
            logger.info(' sigmaN_SD_p      : %.2e GeV^-2 = %.2e pb',self.last_results['sigmaN_SD_proton'],self.last_results['sigmaN_SD_proton']*GeV2pb)
            logger.info(' sigmaN_SD_n      : %.2e GeV^-2 = %.2e pb',self.last_results['sigmaN_SD_neutron'],self.last_results['sigmaN_SD_neutron']*GeV2pb)
        if self.mode['directional']:
            logger.info(' Nevents          : %i', self.last_results['Nevents'])
            logger.info(' smearing         : %.2e', self.last_results['smearing'])
        if self.mode['indirect']:
            for i,v in enumerate(self.maddm_card['halo_dm_velocity']):
                logger.info('Indirect detection cross section at v = %2.e: %.2e+-%2.e', v,
                        self.last_results['indirect#%s' %i],self.last_results['indirect_error#%s' %i])

    
    def is_excluded_relic(self, relic, omega_min = 0., omega_max = 0.1):
        """  This function determines whether a model point is excluded or not
             based on the resulting relic density, spin independent and spin
            dependent cross sections. The user can supply min/max values for
            relic density
        """
        
        


        
    
    def ask_run_configuration(self, mode=None, force=False):
        """ask the question about card edition / Run mode """
        
        if not force:  
            self.mode = self.ask('', '0', mode=mode, data=self.proc_characteristics, 
                            ask_class=MadDMSelector, timeout=60, path_msg=' ')
            if self.mode in ['', '0']:
                process_data = self.proc_characteristics
                self.mode = {'relic': 'ON' if process_data['has_relic_density'] else 'Not available',
                            'direct': 'ON' if process_data['has_direct_detection'] else 'Not available',
                            'directional': 'ON' if process_data['has_directional_detection'] else 'Not available',
                            'indirect': 'ON' if process_data['has_indirect_detection'] else 'Not available'}
                
            self.maddm_card = MadDMCard(pjoin(self.dir_path, 'Cards', 'maddm_card.dat'))
            for key, value in self.mode.items():
                if value == 'ON' or value is True:
                    self.mode[key] = True
    
                else:
                    self.mode[key] = False

            # create the inc file for maddm
            self.maddm_card.set('do_relic_density', self.mode['relic'], user=False)
            self.maddm_card.set('do_direct_detection', self.mode['direct'], user=False)
            self.maddm_card.set('do_directional_detection', self.mode['directional'], user=False)
            self.maddm_card.set('do_indirect_detection', self.mode['indirect'], user=False)
            self.maddm_card.write_include_file(pjoin(self.dir_path, 'include'))
        else:
            if not hasattr(self, 'maddm_card'):
                self.maddm_card = MadDMCard(pjoin(self.dir_path, 'Cards', 'maddm_card.dat'))
            if not hasattr(self, 'mode'):
                process_data = self.proc_characteristics
                self.mode = {'relic': 'ON' if process_data['has_relic_density'] else 'Not available',
                            'direct': 'ON' if process_data['has_direct_detection'] else 'Not available',
                            'directional': 'ON' if process_data['has_directional_detection'] else 'Not available',
                            'indirect': 'ON' if process_data['has_indirect_detection'] else 'Not available'}
            if not os.path.exists(pjoin(self.dir_path, 'include', 'maddm_card.inc')):
                # create the inc file for maddm
                self.maddm_card.set('do_relic_density', self.mode['relic'], user=False)
                self.maddm_card.set('do_direct_detection', self.mode['direct'], user=False)
                self.maddm_card.set('do_directional_detection', self.mode['directional'], user=False)
                self.maddm_card.set('do_indirect_detection', self.mode['indirect'], user=False)
                self.maddm_card.write_include_file(pjoin(self.dir_path, 'include'))                
            
        self.check_param_card(pjoin(self.dir_path, 'Cards', 'param_card.dat'))
        
        if not self.in_scan_mode:    
            logger.info("Start computing %s" % ','.join([name for name, value in self.mode.items() if value]))
        return self.mode


    def compile(self):
        """compile the code"""
        
        if self.in_scan_mode:
            return

        if self.mode['relic'] and self.mode['direct']:
            misc.compile(['all'],cwd=self.dir_path)
        elif self.mode['relic'] and not self.mode['direct']:
            misc.compile(['relic_density'],cwd=self.dir_path)
        elif self.mode['direct'] and not self.mode['relic']:
            misc.compile(['direct_detection'],cwd=self.dir_path)
        else:
            raise Exception, "No computation requested. End the computation"
    
        if os.path.exists(pjoin(self.dir_path, 'src', 'maddm.x')) or os.path.exists(pjoin(self.dir_path, 'maddm.x')):
            logger.info("compilation done")
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
        elif os.path.isfile(os.path.join(path, '%s_card.dat' % args[0])):
            args[0] = os.path.join(path, '%s_card.dat' % args[0])
        elif os.path.isfile(os.path.join(path, 'Cards' ,'%s_card.dat' % args[0])):
            args[0] = os.path.join(path, 'Cards', '%s_card.dat' % args[0])
        elif not os.path.isfile(args[0]):
            raise self.InvalidCmd('No default path for this file')



class MadDMSelector(common_run.EditParamCard):
    """ """

    @property
    def answer(self):
        return self.run_options
    
    def __init__(self, *args, **opts):

        #0. some default variable
        process_data = opts.pop('data', collections.defaultdict(bool))
        self.run_options = {'relic': 'ON' if process_data['has_relic_density'] else 'Not available',
                            'direct': 'ON' if process_data['has_direct_detection'] else 'Not available',
                            'directional': 'ON' if process_data['has_directional_detection'] else 'Not available',
                            'indirect': 'ON' if process_data['has_indirect_detection'] else 'Not available'}
        
        #1. Define what to run and create the associated question
        mode = opts.pop('mode', None)  
        if mode:
            for key, value in mode.items():
                if self.run_options[key] in ['ON', 'OFF']:
                    self.run_options[key] = value      
        question = self.create_question()            
                    
        #2. Define the list of allowed argument
        allow_args = [str(0), 'done', str(4), 'param', str(5), 'maddm']
        for i,key in enumerate(['relic', 'direct', 'directional', 'indirect']):
            if self.run_options[key] in ['ON', 'OFF']:
                allow_args.append(str(i+1))
                allow_args.append(key)
                allow_args.append('%s=on' % key)
                allow_args.append('%s=off' % key)
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
                
        self.has_delphes = False 
     
        
    digitoptions = {1: 'relic', 2:'direct', 3:'directional', 4:'indirect'}
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
    4. Compute Indirect Detection    %(start_underline)sindirect%(stop)s    = %(indirect)s
%(start_green)s You can also edit the various input card%(stop)s:
 * Enter the name/number to open the editor
 * Enter a path to a file to replace the card
 * Enter %(start_bold)sset NAME value%(stop)s to change any parameter to the requested value 
    4. Edit the model parameters    [%(start_underline)sparam%(stop)s]
    5. Edit the MadDM options       [%(start_underline)smaddm%(stop)s]\n""" % \
    {'start_green' : '\033[92m',
     'stop':  '\033[0m',
     'start_underline': '\033[4m',
     'start_bold':'\033[1m', 
     'relic': get_status_str(self.run_options['relic']),
     'direct':get_status_str(self.run_options['direct']),
     'directional':get_status_str(self.run_options['directional']),
     'indirect':get_status_str(self.run_options['indirect'])
     }
        return question

    def define_paths(self, **opt):
        
        super(MadDMSelector, self).define_paths(**opt)
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
                    self.check_coherence(tag)
                self.value = 'repeat'
            elif args[0] in self.run_options:
                if self.run_options[args[0]] == 'ON':
                    self.run_options[args[0]] = 'OFF'
                elif self.run_options[args[0]] == 'OFF':
                    self.run_options[args[0]] = 'ON'
                else:
                    logger.warning('This entry can not be changed for this running directory.')
                self.value = 'repeat'
                self.check_coherence(args[0])
        if hasattr(self, 'value') and self.value == 'repeat':
            self.question = self.create_question()
            return line
        elif os.path.exists(line):
            super(MadDMSelector, self).default(line)
            self.value = 'repeat'
            return line  
        elif not hasattr(self, 'value'):
            return self.run_options
        else:
            return super(MadDMSelector, self).default(line)
        
    
    def do_compute_widths(self, line):
        """normal fct but ensure that self.maddm_card is up-to-date"""
        
        try:
            self.mother_interface.maddm_card = self.maddm
        except Exception,error:
            logger.error("Invalid command: %s " % error)
            return
        return super(MadDMSelector, self).do_compute_widths(line)
            
        
    def do_update(self, line, timer=0):
        """syntax: update dependent: Change the mass/width of particles which are not free parameter for the model.
                    update missing:   add to the current param_card missing blocks/parameters.
           Bypass dependent mode but if request by the user
        """
         
        if timer == 0:
            return super(MadDMSelector, self).do_update(line)
        else: 
            args = self.split_arg(line)
            if args[0] == 'dependent':
                return
            else:
                return super(MadDMSelector, self).do_update(line)
         
    
    def check_coherence(self, last_modified):
        
        if last_modified == 'directional':
            if self.run_options['directional'] == 'ON':
                self.run_options['direct'] = 'ON'
        elif last_modified == 'direct':
            if self.run_options['direct'] == 'OFF':
                self.run_options['directional'] = 'OFF'            
    
    def check_card_consistency(self):
        
        super(MadDMSelector, self).check_card_consistency()
        
        # If direct detection is ON ensure that quark mass are not zero
        if self.run_options['direct'] == 'ON':
            to_change = []
            for i in range(1,7):
                if self.param_card.get_value('mass', i) == 0.:
                    to_change.append(i)
            if to_change:
                logger.warning('For direct detection the quark mass need to be different of zero. Automatically adding such masses (PDG 2014).')
                quark_masses = {1: 4.8e-3, 2: 2.3e-3,  3: 95e-3, 4: 1.275, 5: 4.18, 6: 173.21}
                for i in to_change:
                    self.do_set('param_card mass %s %s' % (i, quark_masses[i]))
        
    def open_file(self, path):
        
        path = super(MadDMSelector,self).open_file(path)
        
        if path == self.paths['maddm']:
            try:
                self.maddm = MadDMCard(path) 
            except Exception as e:
                logger.error('Current maddm_card is not valid. We are going to use the default one.')
                logger.error('problem detected: %s' % e)
                logger.error('Please re-open the file and fix the problem.')
                logger.warning('using the \'set\' command without opening the file will discard all your manual change')
    
    def complete_set(self, text, line, begidx, endidx, formatting=True):
       """ Complete the set command"""
       try:
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
       except Exception, error:
           misc.sprint(error)
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
        
        if args[0] in ['maddm_card']:
            start = 1
            if args[1] == 'default':
                logging.info('replace %s by the default card' % args[0])
                self.maddm = MadDMCard(self.paths['maddm_default'])
        elif args[0] in self.maddm_set and args[0] not in self.conflict:
            start = 0
        else:
            return super(MadDMSelector, self).do_set(line)
        if args[start+1] == 'default':
            default = self.maddm_def[args[start]]
            self.setDM(args[start], default)
        else:
            if args[start] in self.maddm.list_parameter or \
                   args[start] in self.maddm.dict_parameter:
                val = ' '.join(args[start+1:])
                val = val.split('#')[0]
                self.setDM(args[start], val)
            else:
                self.setDM(args[start], args[start+1:])
        # write the new file
        self.maddm.write(self.paths['maddm'])
        
    def setDM(self, name, value):
        logger.info('modify parameter %s of the maddm_card.dat to %s' % (name, value))
        self.maddm.set(name, value, user=True)
        
        
            
        



class InvalidMaddmCard(banner_mod.InvalidRunCard):
    pass
             
class MadDMCard(banner_mod.RunCard):
    """MadDM card use the read/write of the runcard object"""
    
    filename = 'maddm_card'
    default_include_file = 'maddm_card.inc'
    
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
        self.add_param('do_indirect_detection', False, system=True)
        
        self.add_param('calc_taacs_ann_array', True)
        self.add_param('calc_taacs_dm2dm_array', True)
        self.add_param('calc_taacs_scattering_array', True)
        
        self.add_param('eps_ode', 0.01)
        self.add_param('xd_approx', False)
        
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
        self.add_param('only_two_body_decays', True, include=False)
        self.add_param('halo_dm_velocity', [0.001], include=False)
        
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
        
        
class Indirect_Cmd(me5_interface.MadEventCmdShell):
    
    def __init__(self, *args, **opts):
        
        super(Indirect_Cmd, self).__init__(*args, **opts)
        self.history_header = ''
    
    def do_plot(self, line):
        return
    
    def do_madanalysis5_parton(self, line):
        return
    def do_madanalysis5_hadron(self, line):
        return    
        

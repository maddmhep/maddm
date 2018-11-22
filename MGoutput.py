import logging
import math
import os
import shutil
import sys
import collections
import random
from StringIO import StringIO
import re

import maddm_run_interface

# python routines from MadGraph
import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.file_writers as writers
import madgraph.various.misc as misc
import madgraph.iolibs.file_writers as file_writers
import madgraph.various.banner as banner_mod
import madgraph.iolibs.files as files
import aloha
import aloha.create_aloha as create_aloha
from madgraph import MG5DIR
from madgraph.iolibs.files import cp
from madgraph.loop import loop_exporters
from madgraph.core.base_objects import Process
import madgraph.interface.reweight_interface as rwgt_interface
import madgraph.various.banner as bannermod
from madgraph.core import base_objects


class MYStringIO(StringIO):
    """one stringIO behaving like our dedicated writer for writelines"""
    def writelines(self, lines):
        if isinstance(lines, list):
            return StringIO.writelines(self, '\n'.join(lines))
        else:
            return StringIO.writelines(self, lines)
        
# Root path
MDMDIR = os.path.dirname(os.path.realpath( __file__ ))

#usefull shortcut
pjoin = os.path.join
logger = logging.getLogger('madgraph.maddm')

class MADDMProcCharacteristic(banner_mod.ProcCharacteristic):
    
    def default_setup(self):
        """initialize the directory to the default value""" 

        self.add_param('has_relic_density', False)
        self.add_param('has_direct_detection', False)
        self.add_param('has_directional_detection', False)
        self.add_param('has_indirect_detection', False)
        self.add_param('has_capture', False)
        self.add_param('dm_candidate', [0])
        self.add_param('coannihilator', [0])
        self.add_param('model', '')

#-------------------------------------------------------------------------#
class ProcessExporterMadDM(export_v4.ProcessExporterFortranSA):
    """_________________________________________________________________________#
    #                                                                         #
    #  This class is used to export the matrix elements generated from        #
    #  MadGraph into a format that can be used by MadDM                       #
    #                                                                         #
    #  Based off of the ProcessExporterFortranSA class in MadGraph            #
    #                                                                         #
    #_______________________________________________________________________"""

    # forbid the creation of two Matrix element for initial state flipping
    sa_symmetry = True 
    
    # flag to distinguish different type of matrix-element
    DM2SM = 1999  # DM DM > SM SM
    DM2DM = 1998  # DM DM > DM DM
    DMSM  = 1997  # DM SM > DM SM
    DD    = 1996  # DIRECT DETECTION (EFT) DIAGRAM DM QUARK > DM QUARK

    def __init__(self, dir_path = "", opt=None):
        super(ProcessExporterMadDM, self).__init__(dir_path, opt)
        self.resonances = set() 
        self.proc_characteristic = MADDMProcCharacteristic()
        
    def convert_model(self, model, wanted_lorentz = [], wanted_couplings = []):
        """-----------------------------------------------------------------------#  
        #                                                                       #
        #  Create a full valid MG4 model from a MG5 model (coming from UFO)     #
        #                                                                       #
        #  Based off the routine in the ProcessExporterFortran class.  Needed   #
        #  to asjust it so that we can set the self.opt parameter.  By default  #
        #  this parameter is set to 'madevent' and we need to set it to         #
        #  essentially anything other than 'madevent'                           #
        #                                                                       #
        #---------------------------------------------------------------------"""  

        # here we set self.opt so that we get the right makefile in the Model directory
        #self.opt['export_format']='standalone' #modified by antony
        #self.opt['loop_induced']=False #modified by antony
        
        self.proc_characteristic['model'] = model.get('modelpath+restriction')
        
        out =  super(ProcessExporterMadDM, self).convert_model(model, 
                                               wanted_lorentz, wanted_couplings)
        return out

    def write_procdef_mg5(self,*args,**opts):
        """ do nothing """
        return
    
    def make_model_symbolic_link(self):
        """Make the copy/symbolic links"""
        
        model_path = self.dir_path + '/Source/MODEL/'
        if os.path.exists(pjoin(model_path, 'ident_card.dat')):
            files.mv(model_path + '/ident_card.dat', self.dir_path + '/Cards')
        if os.path.exists(pjoin(model_path, 'particles.dat')):
            files.ln(model_path + '/particles.dat', self.dir_path + '/SubProcesses')
            files.ln(model_path + '/interactions.dat', self.dir_path + '/SubProcesses')
        files.cp(model_path + '/param_card.dat', self.dir_path + '/Cards')
        files.mv(model_path + '/param_card.dat', self.dir_path + '/Cards/param_card_default.dat')
 
        for name in ['coupl.inc', 'input.inc']:
            files.ln(pjoin(self.dir_path, 'Source','MODEL', name),
                 pjoin(self.dir_path, 'include'))
         
    def write_procdef_mg5(self,*args):
        return
    
    def pass_information_from_cmd(self, cmd):
        """pass information from the command interface to the exporter.
           Please do not modify any object of the interface from the exporter.
        """
        
        self.model = cmd._curr_model
        self.dm_particles = cmd._dm_candidate + cmd._coannihilation
        
        # update the process python information
        self.proc_characteristic['dm_candidate'] = [p.get('pdg_code') for p in cmd._dm_candidate]
        self.proc_characteristic['coannihilator'] = [p.get('pdg_code') for p in cmd._coannihilation]

    #-----------------------------------------------------------------------#  
    def copy_template(self, model):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine first checks to see if the user supplied project        #
    #  directory is there and asks to overwrite if it is there.  Then it    #
    #  copies over the template files as the basis for the Fortran code.    #
    #                                                                       #
    #-----------------------------------------------------------------------#
        project_path = self.dir_path
        #print 'project path: ' + project_path
        #print self.mgme_dir
        #print self.dir_path
        # Checks to see if projectname directory exists
        logger.info('Initializing project directory: %s', os.path.basename(project_path))

        # Copy over the full template tree to the new project directory
        shutil.copytree(pjoin(MDMDIR, 'Templates'), project_path)

        # Create the directory structure needed for the MG5 files
        os.mkdir(os.path.join(project_path, 'Source'))
        os.mkdir(os.path.join(project_path, 'Source', 'MODEL'))
        os.mkdir(os.path.join(project_path, 'Source', 'DHELAS'))
        os.mkdir(os.path.join(project_path, 'lib'))

        temp_dir = os.path.join(self.mgme_dir, 'Template')        
        # Add make_opts file in Source
        print os.path.join(temp_dir, 'LO/Source', 'make_opts') #modified by antony
        shutil.copy(os.path.join(temp_dir, 'LO/Source', 'make_opts'), #modified by antony
                    os.path.join(self.dir_path, 'Source'))        
 
        # Add the makefile 
        filename = os.path.join(self.dir_path,'Source','makefile')
        self.write_source_makefile(writers.FortranWriter(filename))            

    def get_dd_type(self, process):
        orders = process.get('orders')
        
        efttype = None
        mode = 'bsm'
        
        for coupling, value in orders.items():
            if coupling.startswith('SIEFF') and value > 0:
                if not efttype:
                    efttype = 'SI'
                else:
                    raise Exception
                if value == 99:
                    mode = 'tot'
                elif value ==2:
                    mode = 'eft'
                     
            elif coupling.startswith('SDEFF') and value > 0:
                if not efttype:
                    efttype = 'SD'
                else:
                    raise Exception  
                if value == 99:
                    mode = 'tot'
                elif value ==2:
                    mode = 'eft'       
            #elif value and mode == 'bsm':
            #    mode = 'tot'

#        misc.sprint(orders, mode, efttype)            
        return mode, efttype
        
                
    
    def get_process_name(self, matrix_element, print_id=True):
        """ """
        process = matrix_element.get('processes')[0]
        process_type = process.get('id')
        
        if process_type in [self.DM2SM, self.DM2DM, self.DMSM]:
            return process.shell_string(print_id=print_id)
        if process_type in [self.DD]:
            efttype, siorsd = self.get_dd_type(process)
            if not siorsd:
                return process.shell_string(print_id=print_id)
            else:
                return "%s_%s_%s" % (efttype.upper(),siorsd.upper(), process.shell_string(print_id=print_id))
            
    def generate_subprocess_directory(self, matrix_element, helicity_model, me_number):
        """Routine to generate a subprocess directory.
           For MadDM, this is just a matrix-element in the matrix_elements directory
        """
        
        
        process_type = [p.get('id') for p in matrix_element.get('processes')][0]
        process_name = self.get_process_name(matrix_element) 

        #super(ProcessExporterFortranMaddm,self).generate_subprocess_directory(matrix_element, helicity_model, me):
        path_matrix = pjoin(self.dir_path, 'matrix_elements')        
        
        if process_type in [self.DM2SM, self.DM2DM, self.DD, self.DMSM]:
            # Using the proess name we create the filename for each process and export the fortran file
            filename_matrix = os.path.join(path_matrix, 'matrix_' + process_name + '.f')
        else:
            raise Exception
        
        self.write_matrix_element(writers.FortranWriter(filename_matrix),\
                            matrix_element, helicity_model)

        return 0 # return an integer stating the number of call to helicity routine
    


    #-----------------------------------------------------------------------#
    def write_matrix_element(self, writer, matrix_element, fortran_model):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  Export a matrix element to a matrix.f file in MG4 standalone format  #
    #                                                                       #
    #  This is essentially identical to the write_matrix_element that is    #
    #  found in the standalone output in MadGraph.  The only difference is  #
    #  we include a way to change the name of the smatrix and matrix        #
    #  subroutines so they include the name of the process.                 #
    #                                                                       #
    #-----------------------------------------------------------------------#
        if not matrix_element.get('processes') or \
               not matrix_element.get('diagrams'):
            return 0

        if not isinstance(writer, writers.FortranWriter):
            raise writers.FortranWriter.FortranWriterError(\
                "writer not FortranWriter")


        # first track the S-channel resonances:
        proc_id = matrix_element.get('processes')[0].get('id') 
        if  proc_id == self.DM2SM:
            self.add_resonances(matrix_element, fortran_model)
        
        # track the type of matrix-element and update the process python info
        if proc_id in [self.DM2SM, self.DM2DM, self.DMSM]:
            self.proc_characteristic['has_relic_density'] = True
        if proc_id in [self.DD]:
            self.proc_characteristic['has_direct_detection'] = True
            self.proc_characteristic['has_directional_detection'] = True
            self.proc_characteristic['has_capture'] = True

        # Set lowercase/uppercase Fortran code
        writers.FortranWriter.downcase = False


        replace_dict = super(ProcessExporterMadDM,self).write_matrix_element_v4(
                                            None, matrix_element, fortran_model)


        # Extract the process information to name the subroutine
        process_name = self.get_process_name(matrix_element, print_id=False) 
        replace_dict['proc_prefix'] = process_name + '_'

        #p = pjoin(MDMDIR, 'Templates', 'matrix_elements', 'matrix_template.inc')
        file = open(replace_dict['template_file']).read()
        file = file % replace_dict
        if replace_dict['nSplitOrders'] !='':
            file = file + '\n' + open(replace_dict['template_file2'])\
                                                            .read()%replace_dict
        # Write the file
        writer.writelines(file)

        return replace_dict['return_value']
    
    def add_resonances(self, matrix_element, fortran_model):
        """keep track of all the s-channel resonances in DM DM > SM SM"""
        
        for diag in matrix_element.get('diagrams'):
            for amp in diag.get('amplitudes'):
                init_states = matrix_element.get('processes')[0].get_ninitial()
                #print "Model name: ",
                #print helas_model.get('name')
                #print helas_model.get('particle_dict').keys()
                s_channels, _ = amp.get_s_and_t_channels(init_states, self.model,0)
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
                    resonance_mass = self.model.get_particle(resonance_pdg).get('mass')
                    resonance_width = self.model.get_particle(resonance_pdg).get('width')
                    self.resonances.add((resonance_pdg, resonance_mass, resonance_width))

    
    def finalize(self, matrix_element, cmdhistory, MG5options, outputflag):
        """ """
        
        
        self.global_dict_info = {}
        # Create the dm_info.inc file that contains all the DM particles 
        #as well as spin and mass information
        self.WriteDMInfo(matrix_element)
        # Create the model_info.txt file that is used to calculate the number 
        #of relativistic degrees of freedom.
        self.WriteModelInfo()
        
        # Create the diagrams.inc file that contains the number of processes for each pair of initial state particles
        self.WriteDiagramInfo(matrix_element)

        # Create the process_names.inc file that contains the names of all the individual processes
        self.WriteProcessNames(matrix_element)

        # add quark masses definition in DirectDetection file (needed even for 
        #relic density since the code needs to compile.
        self.WriteMassDirectDetection()
        
        # Create the smatrix.f file
        self.Write_smatrix(matrix_element)
# 
#           # Create the makefile for compiling all the matrix elements
        self.Write_makefile()
# 
#           # Write the include file
        self.WriteMadDMinc()
# 
#           #Write the include file containing all the resonance locations
        self.Write_resonance_info()
        
        # write the details for python
        self.proc_characteristic.write(pjoin(self.dir_path, 'matrix_elements', 'proc_characteristics'))    


        # write maddm_card
        maddm_card = maddm_run_interface.MadDMCard()
        maddm_card.write(pjoin(self.dir_path, 'Cards', 'maddm_card.dat'))
        maddm_card.write(pjoin(self.dir_path, 'Cards', 'maddm_card_default.dat'))


    #-----------------------------------------------------------------------#
    def WriteDMInfo(self, matrix_element):
        """This function writes the information about the inital state particle 
           dof to a file for use by the FORTRAN part of the code."""

        # need to have access to 
        #   self._dm_particles (+coannihilation?)
        #   self._do_relic_density --assume always True here.
        #   self._dm_thermal_scattering[counter-1]
        

        # Look at which particles have the thermal scattering diagrams and those that do
        # get flagged so that the fortran side will know how to properly handle the numerical code
        has_thermalization = set()
        for me in matrix_element:
            proc = me.get('matrix_elements')[0].get('processes')[0]
            if proc.get('id') != 1997: # DM scattering
                continue
            id1, id2 = proc.get_initial_ids()
            has_thermalization.add(id1)
            
        
        #for i in range(len(self._dm_particles)):
        #        sum_thermal_scattering = 0
        #        for j in range(len(self._dm_particles)):
        #          sum_thermal_scattering += len(self._scattering_me[i][j])
        #        if (sum_thermal_scattering == 0):
        #          self._dm_thermal_scattering.append(False)
        #        else:
        #          self._dm_thermal_scattering.append(True)

        path = pjoin(self.dir_path, 'include', 'dm_info.inc')
        writer = file_writers.FortranWriter(path, 'w')
        
        # Write out some header comments for the file
        writer.write_comments("c------------------------------------------------------------------------------c")
        writer.write_comments("c This file contains the needed information about all the DM particles.")

        for i, dm_particle in enumerate(self.dm_particles):
            # Write out the mass parameter, degrees of freedom, particle name, and the index for the end of the name
            writer.write_comments("c------------------------------------------------------------------------------c")
            writer.writelines(""" 
            mdm( %(id)s ) = abs( %(mass)s )
            dof_dm( %(id)s ) = %(dof)s
            dm_sm_scattering( %(id)s ) = %(is_scatter)s
            dm_names( %(id)s ) =\'%(name)s\'
            dm_index( %(id)s ) = %(len_name)s
            dm_antinames( %(id)s ) = \'%(antiname)s\'
            dm_antiindex( %(id)s ) = %(len_antiname)s
            """ % {'id': i+1,
                   'mass': dm_particle['mass'],
                   'dof' : float(dm_particle['spin'])*float(dm_particle['color']),
                   'is_scatter': '.true.' if (dm_particle['pdg_code'] in has_thermalization) else '.false.',
                   'name': dm_particle['name'],
                   'len_name': len(dm_particle['name']),
                   'antiname':dm_particle['antiname'],
                   'len_antiname':len(dm_particle['antiname'])}
                              )
        
        writer.write_comments("c------------------------------------------------------------------------------c")
        writer.close()

    def WriteModelInfo(self):
        """ This routine writes the mass, degrees of freedom and boson/fermion        #
        #  information for all the relevant particles in the model so that the        #
        #  number of relativistic degrees of freedom can be calculated.                #
        #                                                                        #
        #  The SM particles are in the template by default so this would just        #
        #  add the relevant BSM particles."""                                        #


        #need to define 
        #    self._bsm_final_states
        #    
        bsm_particles = [p for p in self.model.get('particles')                
                         if p.get('pdg_code') > 25 and
                         p not in self.dm_particles]

        path = pjoin(self.dir_path, 'include', 'model_info.txt')
        # Open up the template and output file
        model_info_file = open(path, 'w')

        self.global_dict_info['nb_part'] = 17 + len(self.dm_particles) + len(bsm_particles)
        
        # Read all the lines from the template file and insert appropriate parts
        # which are flagged with __flag1 and __flag2
        for line in open(pjoin(self.dir_path, 'include', 'model_info_template.txt')):

            # write out the number of DM paticles
            if "#NUM_PARTICLES" in line:
                nb_particles = 17 + len(self.dm_particles) + len(bsm_particles)
                model_info_file.write("%i\n" % nb_particles)
            # write out the mass, degrees of freedom and fermion/boson info for the bsm particles
            elif "#BSM_INFO" in line:
                for bsm_particle in (self.dm_particles + bsm_particles):
                    # Get all the necessary information
                    mass = str(self.model.get_mass(bsm_particle['pdg_code'])) 
                    dof = float(bsm_particle['spin'])*float(bsm_particle['color'])
                    # If the particle has an anti particle double the number of degrees of freedom
                    if (not bsm_particle['self_antipart']):
                        dof *= 2.0
                    dof = str(dof)
                    # The spin information is written in 2s+1 format 
                    if (int(bsm_particle['spin']) == (int(bsm_particle['spin']/2)*2)):
                        boson_or_fermion = '1'
                    else:
                        boson_or_fermion = '0'
                    # write the line to the file
                    new_line ='         '.join([mass,dof, boson_or_fermion])
                    model_info_file.write(new_line + '\n')
            else:
                # If there is no flag then it's a SM particle which we just write to the output file
                model_info_file.write(line)
        model_info_file.close()
   
    #-----------------------------------------------------------------------#
    def WriteDiagramInfo(self, matrix_element_list):
        """This routine creates the diagrams.inc file which contains the number
           of annihilation diagrams for each pair of DM particles."""
                                                                                

        # open up the file and first print out the number of dm particles
        path  = pjoin(self.dir_path, 'include', 'diagrams.inc')
        fsock = file_writers.FortranWriter(path,'w')

        #don't include the bsm particles which are not dark matter
        ndm_particles =0
        for dm_part in self.dm_particles:
            if (dm_part['charge']==0 and self.model.get_width(dm_part) == 0):
                ndm_particles += 1 

        fsock.write_comments("Total number of IS particles participating in the coannihilations")
        fsock.writelines(' ndmparticles = %s\n\n' %  len(self.dm_particles))
        fsock.write_comments("Processes by class")

        # Compute the number of annihilation processes for each DM pair
        nb_annihilation = collections.defaultdict(int)
        nb_dm2dm = collections.defaultdict(int)
        nb_scattering = collections.defaultdict(int)
        for m in matrix_element_list.get_matrix_elements():
            p = m.get('processes')[0]
            tag = p.get('id')
            ids = self.make_ids(p, tag)
            if tag == self.DM2SM:
                nb_annihilation[ids] +=1
            elif tag == self.DM2DM:
                nb_dm2dm[ids] +=1    
            elif tag == self.DMSM:
                nb_scattering[ids] +=1                                               
            
        fsock.write_comments("Number of annihilation processes for each DM pair")   
        for i,dm1 in enumerate(self.dm_particles):
            for j,dm2 in enumerate(self.dm_particles[i:],i):
                ids = self.make_ids(dm1, dm2, self.DM2SM)
                fsock.writelines(' ann_nprocesses(%i,%i) = %i\n' % (i+1, j+1, nb_annihilation[ids]))

        for i,dm1 in enumerate(self.dm_particles):
            for j,dm2 in enumerate(self.dm_particles[i:],i):
                ids = self.make_ids(dm1, dm2, self.DM2DM)
                fsock.writelines(' dm2dm_nprocesses(%i,%i) = %i\n'% (i+1, j+1, nb_dm2dm[ids]))
        
        fsock.write_comments('\n Number of DM/SM scattering processes for each DM pair\n')

        for i,dm1 in enumerate(self.dm_particles):
            for j,dm2 in enumerate(self.dm_particles):
                ids = self.make_ids(dm1, dm2, self.DMSM)
                fsock.writelines(' scattering_nprocesses(%i,%i) = %i\n' %\
                            (i+1,j+1,nb_scattering[ids]))
   
        fsock.close()
        
    #-----------------------------------------------------------------------#
    def WriteProcessNames(self, matrix_element_list):
        """This routine creates the process_names.inc file which contains the
           names of all the individual processes.  These names are primairly
           used in the test subroutines on the fortran side."""
               
        # This creates the file process_names.inc which will have a list of all the individual process names.
        path = pjoin(self.dir_path, 'include', 'process_names.inc')
        fsock = file_writers.FortranWriter(path,'w')
 
        # Write out the list of process names
        fsock.write_comments("List of the process names in order of dmi, dmj, nprocesses(dmi, dmj)")

         
        # Annihilation process names
        # FOR RELIC DENSITY
 
        #if (self._do_relic_density == True):
        fsock.write_comments("Annihilation process names")
        total_annihilation_nprocesses = 0
        total_dm2dmscattering_nprocesses = 0
        total_scattering_nprocesses = 0
        
        annihilation = collections.defaultdict(list)   # store process name
        annihilation_iden = collections.defaultdict(list) # store if initial particle are identical
        dm2dm = collections.defaultdict(list)          # store process name
        dm2dm_fs = collections.defaultdict(list)       # store the pdg of the final state
        dm2dm_iden = collections.defaultdict(list)    # store if initial particle are identical
        scattering = collections.defaultdict(list)     # store process name
        scattering_sm =  collections.defaultdict(list) # store the pdg of the sm particles
        dd_names = {'bsm':[], 'eft':[],'tot':[]}       # store process name
        dd_initial_state = {'bsm':[], 'eft':[],'tot':[]} # store the pdg of the initial state
        
        for me in matrix_element_list.get_matrix_elements():

            p = me.get('processes')[0]
            tag = p.get('id')
            p1,p2,p3,p4 = p.get('legs')
            name = p.shell_string(print_id=False)
            ids = self.make_ids(p, tag)
            if tag == self.DM2SM:
                total_annihilation_nprocesses += 1
                annihilation[ids].append(name)
                annihilation_iden[ids].append(p1.get('id') == p2.get('id'))
            elif tag == self.DM2DM:
                total_dm2dmscattering_nprocesses += 1 
                
                dm2dm[ids].append(name)
                dm2dm_fs[ids].append([abs(p3.get('id')), abs(p4.get('id'))])
                dm2dm_iden[ids].append(p1.get('id') == p2.get('id'))   
            elif tag == self.DMSM:
                total_scattering_nprocesses +=1
                scattering[ids].append(name)
                scattering_sm[ids].append(abs(p2.get('id')))
            elif tag == self.DD:
                ddtype, si_or_sd = self.get_dd_type(p)
                logger.debug('type: %s' % ddtype)
                dd_names[ddtype].append(self.get_process_name(me, print_id=False))
                dd_initial_state[ddtype].append(p.get_initial_ids()[1])

            #logger.debug('dd_names:')
            #logger.debug(dd_names)
                                                             
        
        # writting the information
        process_counter = 0
        fsock.write_comments("Number of annihilation processes for each DM pair")   
        for i,dm1 in enumerate(self.dm_particles):
            for j,dm2 in enumerate(self.dm_particles[i:],i):
                ids = self.make_ids(dm1, dm2, self.DM2SM)
                for name in annihilation[ids]:
                    process_counter += 1
                    fsock.writelines('process_names(%i) = \'%s\'' % (process_counter, name))


        for i,dm1 in enumerate(self.dm_particles):
            for j,dm2 in enumerate(self.dm_particles[i:],i):
                ids = self.make_ids(dm1, dm2, self.DM2DM)
                for name in dm2dm[ids]:
                    process_counter += 1
                    fsock.writelines('process_names(%i) = \'%s\'\n' % (process_counter, name)) 

        fsock.write_comments('DM/SM scattering process names')

        for i,dm1 in enumerate(self.dm_particles):
            for dm2 in self.dm_particles:
                ids = self.make_ids(dm1, dm2, self.DMSM)
                for name in scattering[ids]:
                    process_counter += 1
                    fsock.writelines('process_names(%i) = \'%s\'\n' % (process_counter, name)) 

        fsock.write_comments('Total number of processes for each category')
        fsock.writelines(" num_processes = %s \n" % process_counter)
        self.global_dict_info['nb_me'] = process_counter
        fsock.writelines(" ann_num_processes = %s \n" % total_annihilation_nprocesses)
        fsock.writelines(" dm2dm_num_processes = %s \n" % total_dm2dmscattering_nprocesses)
        fsock.writelines(" scattering_num_processes = %s \n" %total_scattering_nprocesses)

        # Write out the boolean flags for the processes with identical initial state particles
        fsock.write_comments('Boolean operators for identical particles in the initial state')
        fsock.write_comments('Annihilation diagrams')
        for i,dm1 in enumerate(self.dm_particles):
            for j,dm2 in enumerate(self.dm_particles[i:],i):
                ids = self.make_ids(dm1, dm2, self.DM2SM)              
                for k,iden in enumerate(annihilation_iden[ids]):
                    fsock.writelines(' ann_process_iden_init(%s,%s,%s) = %s' % \
                               (i+1, j+1, k+1, '.true.' if iden else '.false.'))
                  
        fsock.write_comments('DM -> DM diagrams')
        for i,dm1 in enumerate(self.dm_particles):
            for j,dm2 in enumerate(self.dm_particles[i:],i):
                ids = self.make_ids(dm1, dm2,self.DM2DM) 
                #iden = dm1.get('pdg_code') == dm2.get('pdg_code')         
                for k, iden in enumerate(dm2dm_iden[ids]):
                    fsock.writelines(' dm2dm_process_iden_init(%s,%s,%s) = %s' % \
                               (i+1, j+1, k+1, '.true.' if iden else '.false.'))
            
                  
        fsock.write_comments('Final state information for all the DM -> DM processes')
        
        dm_pdg = [abs(p.get('pdg_code')) for p in self.dm_particles]
        for i,dm1 in enumerate(self.dm_particles):
            for j,dm2 in enumerate(self.dm_particles[i:],i):
                ids = self.make_ids(dm1, dm2, self.DM2DM)
                for k,(fs1,fs2) in enumerate(dm2dm_fs[ids]):
                    fs1 = dm_pdg.index(fs1) +1 # get the index in the DM list
                    fs2 = dm_pdg.index(fs2) +1 
                    fsock.writelines(' dm2dm_fs(%s,%s,%s,1) = %s \n' %\
                                                           (i+1, j+1, k+1, fs1))
                    fsock.writelines(' dm2dm_fs(%s,%s,%s,2) = %s \n' %\
                                                           (i+1, j+1, k+1, fs2))

        fsock.write_comments('Initial state degrees of freedom for all the DM/SM scattering processes')
        fsock.write_comments('And Total initial state degrees of freedom for all the DM/SM scattering processes')
        for i,dm1 in enumerate(self.dm_particles):
            for j,dm2 in enumerate(self.dm_particles): 
                ids = self.make_ids(dm1, dm2, self.DMSM)
                for k, pdgsm in enumerate(scattering_sm[ids]):
                    smpart = self.model.get_particle(pdgsm)
                    dof = smpart['spin'] * smpart['color']
                    total_dof = dof if smpart['self_antipart'] else 2*dof
                    
                    fsock.writelines(' dof_SM(%s,%s,%s) = %s\n' %\
                                    (i+1,j+1,k+1, dof))
                    fsock.writelines(' dof_SM_total(%s,%s,%s) = %s\n' %\
                                    (i+1,j+1,k+1, total_dof))

 
        #Write process information for direct detection

        fsock.write_comments('Total number of Direct Detection processes')
        fsock.writelines(' dd_num_processes = %i\n' % len(dd_names['bsm']))
        self.global_dict_info['nb_me_dd'] = len(dd_names['bsm'])
        fsock.write_comments('Processes relevant for Direct Detection')
        for i,name in enumerate(dd_names['bsm']):
            fsock.writelines(' dd_process_names(%i) = \'%s\'\n' % (i+1, name))
            fsock.writelines(' dd_process_ids(%i) = %s' % (i+1,dd_initial_state['bsm'][i]))
        fsock.write_comments('Processes relevant for Direct Detection (effective vertices)')
        for i,name in enumerate(dd_names['eft']):
            fsock.writelines(' dd_eff_process_names(%i) = \'%s\'\n' % (i+1, name))
            fsock.writelines(' dd_eff_process_ids(%i) = %s' % (i+1,dd_initial_state['eft'][i]))
        fsock.write_comments('Processes relevant for Direct Detection (effective vertices + full vertices)')
        self.global_dict_info['nb_me_dd_eff'] = len(dd_names['eft'])
        for i,name in enumerate(dd_names['tot']):
            fsock.writelines(' dd_tot_process_names(%i) = \'%s\'\n' % (i+1, name))
            fsock.writelines(' dd_tot_process_ids(%i) = %s' % (i+1,dd_initial_state['tot'][i]))
        self.global_dict_info['nb_me_dd_tot'] = len(dd_names['tot'])
        fsock.close()
      
    #-----------------------------------------------------------------------#
    def Write_smatrix(self, matrix_element):
        """This routine creates the smatrix.f file based on the smatrix                #
        template.  This is used to call the appropriate matrix element as        #
        well as the appropriate pmass include file for each process."""
        
        replace_dict = {'pmass_ann':None,
                    'pmass_dm2dm': None,
                    'pmass_scattering': None,
                    'smatrix_ann': None,
                    'smatrix_dm2dm': None,
                    'smatrix_scattering': None,
                    'smatrix_dd': None,
                    'smatrix_dd_eff': None,
                    'smatrix_dd_tot': None,
                    }                
        replace_dict['pmass_ann'] = self.get_pmass(matrix_element, self.DM2SM)
        replace_dict['pmass_dm2dm'] = self.get_pmass(matrix_element, self.DM2DM)
        replace_dict['pmass_scattering'] = self.get_pmass(matrix_element, self.DMSM)
    
        replace_dict['smatrix_ann'] = self.get_smatrix(matrix_element, self.DM2SM)
        replace_dict['smatrix_dm2dm'] = self.get_smatrix(matrix_element, self.DM2DM)
        replace_dict['smatrix_scattering'] = self.get_smatrix(matrix_element, self.DMSM)
    
        replace_dict['smatrix_dd'] = self.get_smatrix_dd(matrix_element, 'bsm')
        replace_dict['smatrix_dd_eff'] = self.get_smatrix_dd(matrix_element, 'eft')
        replace_dict['smatrix_dd_tot'] = self.get_smatrix_dd(matrix_element, 'tot')
    
        text = open(pjoin(MDMDIR, 'python_templates','smatrix_template.f')).read()
        to_write = text % replace_dict
        path = pjoin(self.dir_path, 'matrix_elements', 'smatrix.f')
        fsock = file_writers.FortranWriter(path,'w')
        fsock.writelines(to_write)
        fsock.close()
        
        
    def get_pmass(self, matrix_element, flag=None):
        """return the mass of the final state with a if/else if type of entry"""

        pmass = collections.defaultdict(list) 
        output =[]
        
        for me in matrix_element.get_matrix_elements():
            p = me.get('processes')[0]
            tag = p.get('id')
            ids = self.make_ids(p, flag)
            if tag == flag:
                info = MYStringIO()
                self.write_pmass_file(info, me)                  
                pmass[ids].append(info.getvalue())
        
        for i,dm1 in enumerate(self.dm_particles):
            start_second = (i if flag in [self.DM2DM, self.DM2SM] else 0)
            for j,dm2 in enumerate(self.dm_particles[i:], start_second):
                ids = self.make_ids(dm1, dm2, flag)
                for k,info in enumerate(pmass[ids]):
                    output.append('if ((i.eq.%s) .and. (j.eq.%s) .and. (k.eq.%s)) then \n %s \n'\
                                      % (i+1, j+1, k+1, info))
        
        return '%s \n %s' % (' else'.join(output), ' endif' if output else '')
    
    
    def get_smatrix(self, matrix_element, flag):
        """ """

        smatrix = collections.defaultdict(list) 
        output =[]

        for me in matrix_element.get_matrix_elements():
            p = me.get('processes')[0]
            tag = p.get('id')
            ids = self.make_ids(p, flag)
            if tag == flag:
                info = ' call %s_smatrix(p_ext,smatrix)' % p.shell_string(print_id=False)
                smatrix[ids].append(info)

        maxintype = 0 
        for i,dm1 in enumerate(self.dm_particles):
            start_second = (i if flag in [self.DM2DM, self.DM2SM] else 0)
            for j,dm2 in enumerate(self.dm_particles[i:], start_second):
                ids = self.make_ids(dm1, dm2, flag)
                maxintype = max(len(smatrix[ids]), maxintype)
                for k,info in enumerate(smatrix[ids]):
                    output.append('if ((i.eq.%s) .and. (j.eq.%s) .and. (k.eq.%s)) then \n %s \n'\
                                      % (i+1, j+1, k+1, info)) 

        if flag == self.DM2DM:
            self.global_dict_info['max_dm2dm'] = maxintype
        elif flag == self.DM2SM:
            self.global_dict_info['max_dm2sm'] = maxintype
        return '%s \n %s' % (' else'.join(output), ' endif' if output else '')

    def get_smatrix_dd(self, matrix_element, efttype):        
        """ """

        output =[]
        
        for me in matrix_element.get_matrix_elements():
            p = me.get('processes')[0]
            if p.get('id') != self.DD:
                continue
            if self.get_dd_type(p)[0] != efttype:
                continue
            info = ' call %s_smatrix(p_ext, smatrix)' % self.get_process_name(me, print_id=False)
            output.append('if (k.eq.%s) then \n %s \n' % ( len(output)+1, info)) 
        
        return '%s \n %s' % (' else'.join(output), ' endif' if output else '')
                
                           
    
    ## helping method
    @classmethod       
    def make_ids(cls, p, flag, opt=None):
        """get the ids in term of pdg code.
           Two type of input allowed:
             - process and flag
             - DM1 particle, DM2 particle, flag
             """
        if opt is None:
            p1,p2,p3,p4 = p.get('legs')
            if flag in [cls.DM2DM, cls.DM2SM]:
                ids = [abs(p1.get('id')), abs(p2.get('id'))]
                ids.sort()
            else:
                ids = [p1.get('id'), p3.get('id')]               
            return tuple(ids)             
        else:
            # entry are DM1, DM2 , flag
            dm1, dm2, flag = p, flag, opt
            ids =  [abs(dm1.get('pdg_code')), abs(dm2.get('pdg_code'))]  
            if flag in [cls.DM2DM, cls.DM2SM]:
                ids.sort() 
            return tuple(ids)
    
    
    
    #-----------------------------------------------------------------------#
    def Write_makefile(self):
        """This routine creates the makefile for compiling all the matrix        #
           elements generated by madgraph as well as the overall maddm makefile        #
        """
        
        __flag = '#MAKEFILE_MADDM'
        #FOR THE MADDM MAKEFILE
        makefile = open(pjoin(self.dir_path, 'makefile'),'w')
        makefile_template = open(pjoin(self.dir_path,'makefile_template'),'r')
        
        suffix = 'all'
        prop = self.proc_characteristic
        if prop['has_relic_density'] and not prop['has_direct_detection']:
            suffix = 'relic_density'
        elif not prop['has_relic_density'] and prop['has_direct_detection']:
                suffix = 'direct_detection'
        
        makefile_lines = makefile_template.readlines()
        for line in makefile_lines:
            if __flag in line:
                new_line = '\t-cd src/ && make '+suffix
                new_line = new_line + '\n'
                makefile.write(new_line)
            else:
                makefile.write(line)
        
        makefile.close()
        makefile_template.close()
        os.remove(pjoin(self.dir_path,'makefile_template'))

    def WriteMadDMinc(self):
        """ Routine which writes the maddm.inc file. It takes the existing file and adds
         the information about the location of s-channel resonances.
        """
        
        incfile_template = open(pjoin(MDMDIR, 'python_templates', 'maddm.inc')).read()
        res_dict = {}
        #             chunk_size = 5
        init_lines=[]
        #             for i,k in enumerate(self._resonances):
        #                 init_lines.append("resonances(%d)=%s"%(i+1,k[1]))
        res_dict['resonances_initialization'] = '\n'.join(init_lines)
        res_dict['n_resonances'] = len(self.resonances)
        res_dict['nb_part'] = self.global_dict_info['nb_part']
        res_dict['nb_dm'] = len(self.dm_particles)
        res_dict['max_dm2dm'] = 1000#self.global_dict_info['max_dm2dm']
        res_dict['max_dm2sm'] = self.global_dict_info['max_dm2sm']
        res_dict['nb_me'] = self.global_dict_info['nb_me']
        res_dict['nb_me_dd'] = self.global_dict_info['nb_me_dd']
        res_dict['nb_me_dd_eff'] = self.global_dict_info['nb_me_dd_eff']
        res_dict['nb_me_dd_tot'] = self.global_dict_info['nb_me_dd_tot']
        
        writer = writers.FortranWriter(pjoin(self.dir_path, 'include', 'maddm.inc'))
        writer.write(incfile_template % res_dict)
        writer.close()

    #---------------------------------------------------------------------------------

    #-----------------------------------------------------------------------#
    def Write_resonance_info(self):
        """Write the locations of resonances"""
        
        
        resonances_inc_file = open(pjoin(self.dir_path, 'include', 'resonances.inc'), 'w+')
                                   
        for i,res in enumerate(self.resonances):
            new_line = ('                resonances(%d)     =   %s \n' + \
                        '                 resonance_widths(%d)     =   %s \n') %\
                            (i+1, res[1], i+1,res[2])
            resonances_inc_file.write(new_line)

        resonances_inc_file.close()
 

    #-----------------------------------------------------------------------#   
    def WriteMassDirectDetection(self):
        """Adding the light quark mass definition in direct_detection.f"""
        
        writer = open(pjoin(self.dir_path, 'src', 'direct_detection.f'), 'w')
        
        to_replace = {'quark_masses':[]}
        for pdg in range(1,7):
            p = self.model.get_particle(pdg)
            to_replace['quark_masses'].append('M(%s) = %s' % (pdg, p.get('mass')))
        to_replace['quark_masses'] = '\n           '.join(to_replace['quark_masses'])
        
        writer.write(open(pjoin(MDMDIR, 'python_templates', 'direct_detection.f')).read() % to_replace)
        

class Indirect_Reweight(rwgt_interface.ReweightInterface):
    
    def __init__(self, *args, **opts):
        
        self.velocity = 1e-3
        self.shuffling_mode = 'reweight'
        super(Indirect_Reweight, self).__init__(*args, **opts)
        self.output_type = '2.0'
            
    def change_kinematics(self, event):
        
        old_sqrts = event.sqrts
        m = event[0].mass
        new_sqrts = self.maddm_get_sqrts(m, self.velocity)
        jac =  event.change_sqrts(new_sqrts)
        
        def flux(S, M1, M2):
            return S**2 + M1**4 + M2**4 - 2.*S*M1**2 - 20*M1**2*M2**2 - 2.*S*M2**2
    
        jac *= flux(old_sqrts**2, m,m)/flux(new_sqrts**2, m,m)
        
        
        if self.output_type != 'default':
            mode = self.run_card['dynamical_scale_choice']
            if mode == -1:
                if self.dynamical_scale_warning:
                    logger.warning('dynamical_scale is set to -1. New sample will be with HT/2 dynamical scale for renormalisation scale')
                mode = 3
            event.scale = event.get_scale(mode)
            event.aqcd = self.lhe_input.get_alphas(event.scale, lhapdf_config=self.mother.options['lhapdf'])
         
        return jac,event
        
    def calculate_matrix_element(self,*args, **opts):
        if self.shuffling_mode == 'shuffle':
            return 1.0
        else:
            return super(Indirect_Reweight, self).calculate_matrix_element(*args,**opts)
    
    def create_standalone_directory(self, *args, **opts):
        if self.shuffling_mode == 'shuffle':
            return
        return super(Indirect_Reweight, self).create_standalone_directory(*args, **opts)
    
    @staticmethod
    def maddm_get_sqrts(m, ve):

        R = random.random()
        f = lambda v: (-1/ve**2 * math.exp(-v**2/ve**2) * (ve**2+v**2) +1) - R
        df = lambda v: 2. * v**3/ve**4 *math.exp(-v**2/ve**2)
        v = misc.newtonmethod(f, df, ve, error=1e-7, maxiter=250)

        return 2* math.sqrt(m**2 / (1-v/2.)**2)    
        
    def do_change(self, line):    
        
        args = self.split_arg(line)
        if len(args)>1:
            if args[0] == 'velocity':
                self.velocity = bannermod.ConfigFile.format_variable(args[1], float, 'velocity')
                return 
            elif args[0] == 'shuffling_mode':
                if args[1].lower() in ['reweight', 'shuffle']:
                    self.shuffling_mode = args[1]
                return
        return super(Indirect_Reweight, self).do_change(line)


class ProcessExporterIndirectD(object):
    """_________________________________________________________________________#
    #                                                                         #
    #  This class is used to export the matrix elements generated from        #
    #  MadGraph into a format that can be used by MadDM                       #
    #                                                                         #
    #  This is basically a standard MadEvent output with some special trick   #
    #                                                                         #
    #_______________________________________________________________________"""
 
    check = True
    exporter = 'v4'
    output = 'Template'
    sa_symmetry = False
     
    def __new__(cls, path, options):
        if cls is ProcessExporterIndirectD:
            if 'compute_color_flows' in options:
                return loop_exporters.LoopInducedExporterMEGroup.__new__(ProcessExporterIndirectDLI, path, options) 
            else:
                return export_v4.ProcessExporterFortranMEGroup.__new__(ProcessExporterIndirectDLO, path, options)
        
    #===========================================================================
    # link symbolic link
    #===========================================================================
    def finalize(self, matrix_elements, history, mg5options, flaglist):
        """Finalize Standalone MG4 directory"""
        
        # modify history to make an proc_card_mg5 more similar to the real process
        # hidden specific MadDM flag and or replace them by equivalent

        new_history = self.modify_history(history)
        
        super(ProcessExporterIndirectD, self).finalize(matrix_elements, new_history, mg5options, flaglist)
        
        path = pjoin(self.dir_path, 'Cards', 'param_card.dat')
        if os.path.exists(pjoin('..','..', 'Cards', 'param_card.dat')):
            try:
                os.remove(path)
            except Exception:
                pass
            
            files.ln('../../Cards/param_card.dat', starting_dir=os.path.dirname(path),
                     cwd=os.path.dirname(path))   
        
        filename = os.path.join(self.dir_path, 'Cards', 'me5_configuration.txt')
        self.cmd.do_save('options %s' % filename.replace(' ', '\ '), check=False,
                         to_keep={'mg5_path':MG5DIR})
        
        self.write_procdef_mg5( pjoin(self.dir_path, 'SubProcesses', \
                                     'procdef_mg5.dat'),
                                self.cmd._curr_model['name'],
                                self.cmd._generate_info)
        
        self.modify_banner()
        

    def modify_banner(self):
        """enforce that <init> in events.lhe have id 52 for both beam (ensure that py8 accepts it)"""
        
        all_lines = open(pjoin(self.dir_path,'bin','internal','banner.py')).readlines()
        
        for i, line in enumerate(all_lines):
            if 'def get_idbmup(lpp):' in line:
                break

        next_line = all_lines[i+1] 
        nb_space = next_line.find(next_line.strip())
        all_lines.insert(i+2, '%sreturn 52 #enforce DM ID for pythia8\n' % (' '*nb_space))
        
        open(pjoin(self.dir_path,'bin','internal','banner.py'),'w').writelines(all_lines)
        
    def modify_history(self, history):
        """modify history to make an proc_card_mg5 more similar to the real process
           hidden specific MadDM flag and or replace them by equivalent"""       
        
        new_history = []
        next_generate = True #put on True if the next add process need to be switch to generate
        for line in history:
            line = re.sub('\s+', ' ', line)
        
            if line.startswith(('define darkmatter', 'define benchmark','define coannihilator')):
                continue
            
            if line.startswith('generate'):
                next_generate = False
                if line.startswith(('generate relic','generate direct')):
                    next_generate = True
                    continue
                if '@' in line:
                    index = line.find('@')
                    if not line[index:].startswith(('@1995','@ID')):
                        next_generate = True
                        continue  
                    
                if line.startswith('generate indirect'):
                    new_history +=  [p.get('process').nice_string(
                        prefix='add process ' if i>0 else 'generate ') 
                                     for i,p in enumerate(self.cmd._curr_amps)]
                    continue
            
            if line.startswith('add'):
                if line.startswith(('add relic','add direct')):
                    continue
                if '@' in line:
                    index = line.find('@')
                    if not line[index:].startswith(('@1995','@ID')):
                        continue  
                if line.startswith('add indirect'):
                    new_history +=  [p.get('process').nice_string(
                        prefix='generate ' if (i==0 and next_generate) else 'add process ') 
                                     for i,p in enumerate(self.cmd._curr_amps)]
                    next_generate = False
                    continue
                
                if next_generate and line.startswith('add process'):
                    line = line.replace('add process', 'generate')
                    next_generate = False
                
                
            #default behavior propagate
            new_history.append(line)
                  
        return bannermod.ProcCard(new_history)
    #def copy_template(self, model):
    #    
    #    out = super(ProcessExporterIndirectD,self).copy_template(model)        
    #    self.modify_dummy()
    #
    #   return out
    #===========================================================================
    # write_source_makefile
    #===========================================================================
    #def write_source_makefile(self, writer):
    #    """Write the nexternal.inc file for MG4"""
    #
    #    replace_dict = super(ProcessExporterIndirectD, self).write_source_makefile(None)
    #
    #    path = pjoin(MG5DIR,'madgraph','iolibs','template_files','madevent_makefile_source')


        
    #    if writer:
    #        text = open(path).read() % replace_dict
    #        writer.write(text)
    #        
    #    return replace_dict
        
        
    def pass_information_from_cmd(self, cmd):
        """pass information from the command interface to the exporter.
           Please do not modify any object of the interface from the exporter.
        """
        
        self.cmd=cmd
        return super(ProcessExporterIndirectD, self).pass_information_from_cmd(cmd)
             
             
             
    def convert_model(self, model, wanted_lorentz = [],
                             wanted_couplings = []):
        """ Create a full valid MG4 model from a MG5 model (coming from UFO)"""

        # make sure that mass/width external parameter without assoicated particle
        # are removed. This should not be necessary since 2.6.5
        
        def_part = [p['pdg_code'] for p in self.model['particles']]
        for param in self.model['parameters'][('external',)][:]:            
            if param.lhablock in ['MASS','DECAY']:
                if param.lhacode[0] not in def_part:
                    self.model['parameters'][('external',)].remove(param)
                    param = base_objects.ModelVariable(param.name, str(param.value), 'real')
                    self.model['parameters'][tuple()].append(param)
        
        
        super(ProcessExporterIndirectD, self).convert_model(model, 
                                                wanted_lorentz, wanted_couplings)


    #===========================================================================
    # create the run_card
    #=========================================================================== 
    def create_run_card(self, matrix_elements, history):
        """ """
 
        run_card = banner_mod.RunCard()
    
        # pass to new default
        run_card["run_tag"] = "\'DM\'"
        run_card['dynamical_scale_choice'] = 4
        run_card['lpp1'] = 0
        run_card['lpp2'] = 0
        run_card['use_syst'] = False
        run_card.remove_all_cut()
                  
        run_card.write(pjoin(self.dir_path, 'Cards', 'run_card_default.dat'),
                       template=pjoin(MG5DIR, 'Template', 'LO', 'Cards', 'run_card.dat'),
                       python_template=True)
        run_card.write(pjoin(self.dir_path, 'Cards', 'run_card.dat'),
                       template=pjoin(MG5DIR, 'Template', 'LO', 'Cards', 'run_card.dat'),
                       python_template=True)
    
        
         
class ProcessExporterIndirectDLO(ProcessExporterIndirectD,export_v4.ProcessExporterFortranMEGroup):
    pass

class ProcessExporterIndirectDLI(ProcessExporterIndirectD,loop_exporters.LoopInducedExporterMEGroup):
    pass

import os
import shutil
import sys

# python routines from MadGraph
import madgraph.iolibs.export_v4 as export_v4
import madgraph.iolibs.file_writers as writers
import aloha
import aloha.create_aloha as create_aloha
from madgraph import MG5DIR
from madgraph.iolibs.files import cp

# Root path
rp = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[0]
sys.path.append(rp)

# Current directory path
pp = os.path.split(os.path.dirname(os.path.realpath( __file__ )))[1]
sys.path.append(pp)

#-------------------------------------------------------------------------#
class ProcessExporterFortranMadDM(export_v4.ProcessExporterFortran):
#_________________________________________________________________________#
#                                                                         #
#  This class is used to export the matrix elements generated from        #
#  MadGraph into a format that can be used by MadDM                       #
#                                                                         #
#  Based off of the ProcessExporterFortranSA class in MadGraph            #
#                                                                         #
#_________________________________________________________________________#

    #-----------------------------------------------------------------------#  
    def convert_model_to_mg4(self, model, wanted_lorentz = [],
                             wanted_couplings = []):
    #-----------------------------------------------------------------------#  
    #                                                                       #
    #  Create a full valid MG4 model from a MG5 model (coming from UFO)     #
    #                                                                       #
    #  Based off the routine in the ProcessExporterFortran class.  Needed   #
    #  to asjust it so that we can set the self.opt parameter.  By default  #
    #  this parameter is set to 'madevent' and we need to set it to         #
    #  essentially anything other than 'madevent'                           #
    #                                                                       #
    #-----------------------------------------------------------------------#  

        # here we set self.opt so that we get the right makefile in the Model directory
        
        #self.opt['export_format'] = 'not_madevent' #commented by antony
        self.opt['export_format']='standalone' #modified by antony
        self.opt['loop_induced']=False #modified by antony
        
        # create the MODEL
        write_dir=os.path.join(self.dir_path, 'Source', 'MODEL')
        model_builder = export_v4.UFO_model_to_mg4(model, write_dir, self.opt)
        model_builder.build(wanted_couplings)

        # Create and write ALOHA Routine
        aloha_model = create_aloha.AbstractALOHAModel(model.get('name'))
        if wanted_lorentz:
            aloha_model.compute_subset(wanted_lorentz)
        else:
            aloha_model.compute_all(save=False)
        write_dir=os.path.join(self.dir_path, 'Source', 'DHELAS')
        aloha_model.write(write_dir, 'Fortran')

        #copy Helas Template
        cp(MG5DIR + '/aloha/template_files/Makefile_F', write_dir+'/makefile')
        for filename in os.listdir(os.path.join(MG5DIR,'aloha','template_files')):
            if not filename.lower().endswith('.f'):
                continue
            cp((MG5DIR + '/aloha/template_files/' + filename), write_dir)
        create_aloha.write_aloha_file_inc(write_dir, '.f', '.o')

        # Make final link in the Process
        self.make_model_symbolic_link()

    #-----------------------------------------------------------------------#  
    def copy_template(self,project_path):
    #-----------------------------------------------------------------------#
    #                                                                       #
    #  This routine first checks to see if the user supplied project        #
    #  directory is there and asks to overwrite if it is there.  Then it    #
    #  copies over the template files as the basis for the Fortran code.    #
    #                                                                       #
    #-----------------------------------------------------------------------#
        print 'project path: ' + project_path
        print self.mgme_dir
        print self.dir_path
        # Checks to see if projectname directory exists
        print 'Initializing project directory: %s' % os.path.basename(project_path)

        # Copy over the full template tree to the new project directory
        shutil.copytree('Templates', project_path)

        # Create the directory structure needed for the MG5 files
        os.mkdir(os.path.join(project_path, 'Source'))
        os.mkdir(os.path.join(project_path, 'Source', 'MODEL'))
        os.mkdir(os.path.join(project_path, 'Source', 'DHELAS'))
        os.mkdir(os.path.join(project_path, 'lib'))
        os.mkdir(os.path.join(project_path, 'Cards'))

        temp_dir = os.path.join(self.mgme_dir, 'Template')        
        # Add make_opts file in Source
        print os.path.join(temp_dir, 'LO/Source', 'make_opts') #modified by antony
        shutil.copy(os.path.join(temp_dir, 'LO/Source', 'make_opts'), #modified by antony
                    os.path.join(self.dir_path, 'Source'))        
 
        # Add the makefile 
        filename = os.path.join(self.dir_path,'Source','makefile')
        self.write_source_makefile(writers.FortranWriter(filename))            



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

        # Set lowercase/uppercase Fortran code
        writers.FortranWriter.downcase = False

        replace_dict = {}

        # Extrace the process information to name the subroutine
        process_name = matrix_element.get('processes')[0].shell_string()
        replace_dict['process_name'] = process_name[2:len(process_name)]

        # Extract version number and date from VERSION file
        info_lines = self.get_mg5_info_lines()
        replace_dict['info_lines'] = info_lines

        # Extract process info lines
        process_lines = self.get_process_info_lines(matrix_element)
        replace_dict['process_lines'] = process_lines

        # Extract number of external particles
        (nexternal, ninitial) = matrix_element.get_nexternal_ninitial()
        replace_dict['nexternal'] = nexternal

        # Extract ncomb
        ncomb = matrix_element.get_helicity_combinations()
        replace_dict['ncomb'] = ncomb

        # Extract helicity lines
        helicity_lines = self.get_helicity_lines(matrix_element)
        replace_dict['helicity_lines'] = helicity_lines

        # Extract overall denominator
        # Averaging initial state color, spin, and identical FS particles
        den_factor_line = self.get_den_factor_line(matrix_element)
        replace_dict['den_factor_line'] = den_factor_line

        # Extract ngraphs
        ngraphs = matrix_element.get_number_of_amplitudes()
        replace_dict['ngraphs'] = ngraphs

        # Extract nwavefuncs
        nwavefuncs = matrix_element.get_number_of_wavefunctions()
        replace_dict['nwavefuncs'] = nwavefuncs

        # Extract ncolor
        ncolor = max(1, len(matrix_element.get('color_basis')))
        replace_dict['ncolor'] = ncolor

        # Extract color data lines
        color_data_lines = self.get_color_data_lines(matrix_element)
        replace_dict['color_data_lines'] = "\n".join(color_data_lines)

        # Extract helas calls
        helas_calls = fortran_model.get_matrix_element_calls(\
                    matrix_element)
        replace_dict['helas_calls'] = "\n".join(helas_calls)

        # Extract JAMP lines
        jamp_lines = self.get_JAMP_lines(matrix_element)
        replace_dict['jamp_lines'] = '\n'.join(jamp_lines)

        file = open('Templates/matrix_elements/matrix_template.inc').read()
        file = file % replace_dict

        # Write the file
        writer.writelines(file)

        return len(filter(lambda call: call.find('#') != 0, helas_calls))

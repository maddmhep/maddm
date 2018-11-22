## import the required files

import maddm_interface as maddm_interface
import MGoutput as MGoutput
## Define a typical error for the plugin
class MadDMError(Exception): pass

## Does it define a new interface (will be avaible with ./bin/mg5_aMC --mode=maddm
## Put None if no dedicated command are required
new_interface = maddm_interface.MadDM_interface

 
## Does it define a new output mode. Need to be define
new_output = {'maddm': MGoutput.ProcessExporterMadDM,
              'indirect': MGoutput.ProcessExporterIndirectD}
new_reweight = {'indirect': MGoutput.Indirect_Reweight}

## The test/code have been validated up to this version
latest_validated_version = (2,6,4)
minimal_mg5amcnlo_version = (2,6,2)
maximal_mg5amcnlo_version = (1000,1000,1000)

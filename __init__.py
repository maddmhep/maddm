## import the required files

from __future__ import absolute_import
from .data_install import check_python_dependencies

check_python_dependencies(__name__)

from . import maddm_interface as maddm_interface
from . import MGoutput as MGoutput
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
latest_validated_version = (3,8,1)
minimal_mg5amcnlo_version = (2,9,18)
maximal_mg5amcnlo_version = (3,8,1)

## MadGraph7 checks its own version numbering
## (the v0.3.0 release still reports 'version = 0.2.1' in its VERSION file)
latest_validated_mg7_version = (0,3,0)
minimal_mg7_version = (0,2,1)
maximal_mg7_version = (0,3,0)

## MadGraph7 compares its own version against the *_mg5amcnlo_* fields:
## use the MadGraph7 range there when running under MadGraph7
## (versions 0.x, while MG5_aMC is >= 2)
try:
    import madgraph.various.misc as _misc
    _is_mg7 = int(_misc.get_pkg_info()['version'].split('.')[0]) < 2
except Exception:
    _is_mg7 = False
if _is_mg7:
    latest_validated_version = latest_validated_mg7_version
    minimal_mg5amcnlo_version = minimal_mg7_version
    maximal_mg5amcnlo_version = maximal_mg7_version

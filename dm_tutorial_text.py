################################################################################
#
# Copyright (c) 2018 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################

tutorial = """
You have entered tutorial mode. This will introduce you to the main
syntax options of the MadDM plugin.

As compared with MG5aMC, the following functions have additional options
1) define 
2) generate 
3) add
4) output 
You can type help XXXX to see the associated help.

The goal of this tutorial is to learn how to generate processes for dark matter annihilation and scattering
and to use those to compute the relic density, direct and/or indirect detection.

First you need to import your model of interest. If you do not have a favorite model you
can use one of those available online (that can be download automatically). For example:
MadDM> import model DMsimp_s_spin0
"""

help=tutorial

display_modellist="""
The second stage is to define what is your darkmatter candidate
MadDM> define darkmatter xd

If your model has R-parity and that you do not know what is your DM for a given
benchmark, you can ask us to find it for you with:
MadDM> define darkmatter
You will need to define your benchmark point (edit the file) and then
the code will pick the lightest massive neutral BSM particle (with zero width).
"""

import_model = """
You have just imported a new BSM model.
If you want to see the list of models that can be downloaded automatically you can enter:
MadDM> display modellist

""" +display_modellist


define_coannihilator = """Now that the dark matter candidate  and coannihilator(s) are defined.

 You can specify what type of computation you want to do
   - computation of relic density
     MadDM> generate relic_density
   - computation of direct detection (detector on earth for DM interaction with matter)
     MadDM> add direct_detection
   - computation of indirect detection (annihilation of dark matter in halos and detector on or close to the Earth)
     MadDM> add indirect_detection     

   Note: 'generate' remove all previously generated matrix-element.
         'add' keeps those
         for the rest they are identical
"""




define_darkmatter = """
Now that the dark matter candidate is defined. You can either 

1) search/define coannihilator
[Note that we do not have such particles for the proposed model]
MadDM> define coannihilator xr

2) specify what type of computation you want to do
   - computation of relic density
     MadDM> generate relic_density
   - computation of direct detection (detector on earth for DM interaction with matter)
     MadDM> add direct_detection
   - computation of indirect detection (annihilation of dark matter in halos and detector on or close to the Earth)
     MadDM> add indirect_detection     

   Note: 'generate' remove all previously generated matrix-element.
         'add' keeps those
         for the rest they are identical
"""

generate="""

You can add additional type of computation if you want
   - computation of relic density
     MadDM> add relic_density
   - computation of direct detection (detector on earth for DM interaction with matter)
     MadDM> add direct_detection
   - computation of indirect detection (detector on earth for DM interaction with matter)
     MadDM> add indirect_detection 
     
When you have all the type of computation that you want, you can do
     MadDM> output MY_FIRST_MADDM_RUN
   
"""

add=generate

output="""

You have now in your local directory a new directory called "MY_FIRST_MADDM_RUN".
with all the information needed to make the computation requested.

MadDM> launch MY_FIRST_MADDM_RUN
and then follow instruction on the screen to edit your benchmark point/ choose 
what to run/...

"""

launch="""

Here you are, you have finished your first computation!
and you have reached the end of this tutorial.
you can now type 
MadDM> tutorial quit
"""

set_mxd="""
The above is only valid WITHIN the question trigger by the 'launch' command.
At this level of the interface, the 'set' command controls some options of the runs.
To see the list of such options you can do
MadDM> display options
"""

set_my0=set_mxd

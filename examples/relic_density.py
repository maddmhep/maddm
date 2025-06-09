"""

=======================================
Relic Density Computation in MadDM
=======================================

JUST A CHATGPT EXAMPLE

MadDM can compute the relic density of dark matter via numerical solution of the Boltzmann equation, including coannihilation effects and resonances when relevant.

This allows you to constrain models based on the observed dark matter abundance from cosmological data, such as the one measured by Planck.

Example with a simplified scalar mediator model

First, import the model. If not already downloaded, it will be fetched automatically from here <https://feynrules.irmp.ucl.ac.be/wiki/DMsimp>_:

.. code-block:: text

MadDM> import model DMsimp_s_spin0

Then, define the dark matter particle. For scalar DM models, this could be Xr (real scalar), Xc (complex scalar), or Xd (Dirac fermion):

.. code-block:: text

MadDM> define darkmatter Xr
MadDM> generate relic_density

You can now configure your model parameters (in param_card.dat) by pressing 7, or directly editing the file inside the process directory.

Next, generate the process folder. For example:

.. code-block:: text

MadDM> output RD_example
MadDM> launch RD_example
At launch, you can modify the MadDM options (in maddm_card.dat) by pressing 8, or editing the file at RD_example/Cards/maddm_card.dat.

By default, MadDM uses precise thermal averaging and includes coannihilation processes if the relevant particles are close in mass.

Once the setup is complete, press Enter to compute the relic density.

The output will look like:

.. code-block:: text

INFO: compilation done
INFO: MadDM Results
INFO: Omega h^2 (total) = 0.118        ALLOWED     PLANCK 2σ range = [0.114, 0.126]
INFO: 
INFO: Results written in: /Users/yourname/yourprocessfolder/RD_example/output/run_01/MadDM_results.txt

Here, Omega h^2 refers to the dark matter relic abundance predicted by your model.
If it falls within the allowed experimental range (e.g., Planck 2σ), it will be marked as ALLOWED.

You can also inspect additional details of the computation, such as the annihilation channels contributing to freeze-out, in the output folder.

"""
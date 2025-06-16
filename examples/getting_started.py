"""
===============
Getting Started
===============

Get started with the MadDM command line interface.

MadDM tutorial
==============

To begin, you can access a basic tutorial directly in MadDM by typing ``tutorial`` in the prompt shell.
We strongly recommend you to go through this tutorial, as it will guide you through the most important features of MadDM.

Relic Density, Direct Detection, and Indirect Detection
=======================================================

To generate relic density, direct detection, and indirect detection (including loop-induced processes), follow these steps in the command line:

First, import the model, which will be downloaded automatically from `here <https://feynrules.irmp.ucl.ac.be/wiki/DMsimp>`_ if not already present: 

.. code-block:: text

    MadDM> import model DMsimp_s_spin1_MD

Next, define the dark matter (DM) particle manually with:

.. code-block:: text

    MadDM> define darkmatter ~xd

Generally, for DMsimp models, you can choose between Xr (real scalar DM), Xc (complex scalar DM)
and Xd (Dirac spinor DM).
If not specified, MadDM will automatically select the DM candidate using the following rules:

- The DM particle must be a BSM particle (pdg_code > 25).
- The particle must be neutral (no electric or color charge).
- The particle must be stable (its width must be 0 or 'ZERO' in the param_card.dat).
- Among the particles that satisfy the above, the lightest is chosen.

Finally, the candidate must:

- Have some non-zero couplings,
- Not decay into a pair of SM particles.

Then, generate the direct detection (``direct``), indirect detection (``indirect``) or relic density (``relic_rensity``) processes with:

.. code-block:: text

    MadDM> generate direct

If you want to add additional processes, do so by typing:

.. code-block:: text

    MadDM> add indirect
    MadDM> add relic_density

At this point, you can set the parameter values of the model, stored in ``param_card_orig.dat``, by pressing 1. This will open a vim editor where you can modify the parameters.
Next, create the process folder. In this case, we call it ``DMspin1``, but you can choose any name you like. Then launch the process:

.. code-block:: text

    MadDM> output DMspin1
    MadDM> launch DMspin1

Note that once the process folder has been created with the ``output`` command in a previous session, you can always
relaunch the same process in a new session with ``launch myprocessname``, skipping the previous steps.
Once you launch the process, you will be prompted a selection menu to choose which modules you want to run.

.. code-block:: text

        The following switches determine which programs are run:
    /====================== Description =======================|===================== values =====================|========== other options ===========\\
    | 1. Compute the Relic Density                             |                relic = OFF                       |     ON                             |
    | 2. Compute direct detection - nucleon recoil             |               direct = ON                        |     OFF                            |
    | 3. Compute direct detection - electron recoil            |      direct_electron = ON                        |     OFF                            |
    | 4. Compute indirect detection/flux (cont spectrum)       |             indirect = sigmav                    |     flux_source|flux_earth|OFF     |
    | 5. Compute indirect detection in aX (line spectrum)      |             spectral = OFF                       |     Please install module          |
    | 6. Run Multinest scan                                    |             nestscan = OFF                       |     Please install module          |
    \\==================================================================================================================================================/
    You can also edit the various input card:
    * Enter the name/number to open the editor
    * Enter a path to a file to replace the card
    * Enter set NAME value to change any parameter to the requested value
    /=============================================================================\\ 
    |  7. Edit the model parameters    [param]                                    |  
    |  8. Edit the MadDM options       [maddm]                                    |
    \\=============================================================================/
    [60s to answer]

You can select the modules you want to run by entering the corresponding number or name, or by setting the
corresponding variable to ON:

.. code-block:: text

    MadDM> relic = ON

Some modules can be set to specific values, such as:

.. code-block:: text

    MadDM> set indirect = flux_source

You can set the model parameters in the ``param_card.dat`` by pressing 7, or by editing the file directly in the process folder (in this case ``DMspin1/Cards/param_card.dat``).

You can set the MadDM parameters (such as the halo parameters, recoil energy range, and much more) in the ``maddm_card.dat`` by pressing 8,
or by editing the file directly in the process folder (in this case ``DMspin1/Cards/maddm_card.dat``).

Once you are all set, press Enter to launch the process. The results will be displayed in the terminal as shown below:

.. code-block:: text

    INFO: MadDM Results 
    INFO: 
    ****** Relic Density
    OMEGA IS 11.557143 
    INFO: Relic Density       = 1.16e+01       OVERABUNDANT       
    INFO: x_f                 = 1.40e+01             
    INFO: sigmav(xf)          = 8.50e-29 cm^3 s^-1   
    INFO: xsi                 = 1.00e+00             
    INFO:  
    INFO: Channels contributions: 
    INFO: xxdxxdb_ddx         : 20.08 % 
    INFO: xxdxxdb_uux         : 20.08 % 
    INFO: xxdxxdb_ssx         : 20.08 % 
    INFO: xxdxxdb_ccx         : 20.08 % 
    INFO: xxdxxdb_bbx         : 19.69 % 
    INFO: No contribution from processes: xxdxxdb_y1y1, xxdxxdb_vevex, xxdxxdb_vmvmx, xxdxxdb_vtvtx, xxdxxdb_ttx, xxdxxdb_emep, xxdxxdb_mummup, xxdxxdb_tamtap 
    INFO: 
    ****** Direct detection [cm^2]:  
    INFO: SigmaN_SI_p         All DM = 2.05e-40       EXCLUDED  	LZ2024 ul        = 6.48e-47 
    INFO: SigmaN_SI_n         All DM = 2.06e-40       EXCLUDED  	LZ2024 ul        = 6.48e-47 
    INFO: SigmaN_SD_p         All DM = 2.58e-60       ALLOWED   	Pico60 (2019) ul = 4.16e-41 
    INFO: SigmaN_SD_n         All DM = 1.01e-59       ALLOWED   	LZ2024 ul        = 7.59e-42 
    INFO: 
    ****** Indirect detection: 
    INFO: <sigma v> method: inclusive  
    INFO: ====== continuum spectrum final states 
    INFO: DM particle halo velocity: 2e-05 c  
    INFO: *** Print <sigma v> [cm^3 s^-1] with Fermi dSph limits 
    INFO: xxdxxdb_ddx         All DM = 3.48e-29       ALLOWED   	Fermi dSph ul    = 3.30e-27 
    INFO: xxdxxdb_uux         All DM = 3.48e-29       ALLOWED   	Fermi dSph ul    = 3.31e-27 
    INFO: xxdxxdb_ssx         All DM = 3.48e-29       ALLOWED   	Fermi dSph ul    = 3.30e-27 
    INFO: xxdxxdb_ccx         All DM = 3.48e-29       ALLOWED   	Fermi dSph ul    = 4.54e-27 
    INFO: xxdxxdb_bbx         All DM = 3.41e-29       ALLOWED   	Fermi dSph ul    = 5.98e-27 
    INFO: Skipping zero cross section processes for: xxdxxdb_y1y1, xxdxxdb_ttx, xxdxxdb_emep 
    INFO: Using generic Fermi dSph limits for light quarks (u,d,s) 
    INFO: DM DM > SM SM       All DM = 1.73e-28       ALLOWED   	Fermi dSph ul    = 2.96e-27 
    INFO:  
    INFO: *** Fluxes at earth [particle cm^-2 sr^-1]: 
    INFO: gammas Flux         =	5.67e-12             
    INFO: neutrinos_e Flux    =	2.71e-12             
    INFO: neutrinos_mu Flux   =	5.00e-12             
    INFO: neutrinos_tau Flux  =	4.26e-15             
    INFO:  
    INFO:  
    INFO: Results written in: /root/DMspin1/output/run_01/MadDM_results.txt

The results will be stored in the ``MadDM_results.txt`` file located in the process folder (in this case ``DMspin1/output/run_01/MadDM_results.txt``).
In the process folder, you will also find additional files produced by the different MadDM modules. See the tutorials for those modules for additional details.

Exiting MadDM
=============

.. code-block:: text

    MadDM> quit

Install additional tools
========================

.. code-block:: text

    MadDM> install pythia8
    MadDM> install PPPC4DMID
    MadDM> install dragon
    MadDM> install dragon_data_from_galprop
    MadDM> quit

"""

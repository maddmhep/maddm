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

Then, generate the direct detection (``direct`` or ``direct_detection``), indirect detection (``indirect`` or ``indirect_detection``),
gamma-ray line (``spectral`` or ``indirect_spectral_features``) or relic density (``relic_density``) processes with:

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
    | 2. Compute direct detection - nucleon recoil             |       direct_nucleon = ON                        |     OFF                            |
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

    MadDM> indirect = flux_source

Direct detection has two switches, ``direct_nucleon`` (DM-nucleon scattering) and ``direct_electron`` (DM-electron
scattering, available if the model has DM-electron couplings). What they compute is set by ``direct_mode`` in the
``maddm_card.dat``, with two shortcuts:

.. code-block:: text

    MadDM> set limit_and_rates
    MadDM> set limit

``limit_and_rates`` (the default) computes the cross sections, compares them with the exclusion limits and saves the
recoil spectra of the targets in ``direct_nucleon_targets`` (``xe``, ``ge``, ``ar``, ``f``) and
``direct_electron_targets`` (``xe``, ``ar``); ``limit`` computes only the cross sections and the limits, which is faster
(for instance in scans). See the :doc:`direct detection tutorial <plot_direct_detection>` for the details.

You can set the model parameters in the ``param_card.dat`` by pressing 7, or by editing the file directly in the process folder (in this case ``DMspin1/Cards/param_card.dat``).

You can set the MadDM parameters (such as the halo parameters, the direct detection targets, and much more) in the ``maddm_card.dat`` by pressing 8,
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
    INFO: SigmaN_*_eq: DM-nucleon cross sections equivalent to the model in LZ (same expected signal as the reference model of each LZ limit; RAPIDD chiral EFT, p-n interference and two-body currents included)
    INFO: SigmaN_SI_eq        All DM = 5.07e-41       EXCLUDED  	LZ2024 ul        = 6.48e-47
    INFO: SigmaN_SD_n_eq      All DM = 7.77e-36       EXCLUDED  	LZ2024 ul        = 7.59e-42
    INFO: SigmaN_SD_p_eq      All DM = 2.39e-34       EXCLUDED  	LZ2024 ul        = 2.33e-40
    INFO: SigmaN_SD_p         All DM = 6.46e-61       ALLOWED   	Pico60 (2019) ul = 4.16e-41
    INFO: sigma_p and sigma_n (information only: compared separately with the LZ limits they ignore the p-n interference; SD_p is compared with PICO-60 above):
    INFO: SigmaN_SI_p         All DM = 5.13e-41
    INFO: SigmaN_SI_n         All DM = 5.14e-41
    INFO: SigmaN_SD_n         All DM = 2.52e-60
    INFO: Neutrino fog (xenon): reference boundary, NOT an exclusion limit. https://github.com/cajohare/NeutrinoFog; C. A. J. O'Hare, Phys. Rev. Lett. 127, 251802 (2021), arXiv:2109.03116
    INFO: Rescaled from rho = 0.3 GeV/cm^3 (published tables) to rhoDM = 0.3 GeV/cm^3.
    INFO: neutrino_fog_SI = 4.3291631542e-48 # cm^2; -1 means unavailable
    INFO: neutrino_fog_SD_p = 1.9474339879e-40 # cm^2; -1 means unavailable
    INFO: neutrino_fog_SD_n = 1.0584409940e-42 # cm^2; -1 means unavailable
    INFO: Recoil spectra: dRdE_nucleon.dat in /root/DMspin1/output/run_01
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
    INFO: Positron results will come soon.
    INFO:  
    INFO:  
    INFO: Results written in: /root/DMspin1/output/run_01/MadDM_results.txt

For DM-nucleon scattering the verdicts use the cross sections equivalent to the model in LZ (``SigmaN_SI_eq``, ``SigmaN_SD_n_eq``,
``SigmaN_SD_p_eq``, see the direct detection tutorial); this model has only a vector (SI) coupling, and its SI signal seen through the
SD reference models also exceeds the LZ SD limits. The proton and neutron cross sections are given for information.

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
    MadDM> install cosmixs
    MadDM> install RAPIDD
    MadDM> install py_requirements
    MadDM> quit

``pythia8`` and ``cosmixs`` (CosmiXs source spectra) are used by indirect detection, ``RAPIDD`` by direct detection
(recoil spectra and DM-electron scattering), and ``py_requirements`` installs ``numpy``, ``scipy`` and ``matplotlib``
in the Python environment running MadDM. MadDM also offers to install them when a process first needs them.

"""

# sphinx_gallery_thumbnail_path = '_static/thumbnail/getting_started_thumb.png'
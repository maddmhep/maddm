"""
===============
Getting Started
===============

MadDM features its own command line, designed to resemble the familiar MadGraph command line structure. 
To begin, you can access a basic tutorial directly in MadDM by typing ``tutorial`` in the prompt shell.
We strongly recommend you to go through this tutorial, as it will guide you through the most important features of MadDM.

Example
=======

To generate relic density, direct detection, and indirect detection (including loop-induced processes), follow these steps in the command line:

.. code-block:: text

    MadDM> import model DMsimp_s_spin0
    MadDM> define darkmatter xd  # Necessary only if you want to 'force' a DM candidate
    MadDM> generate relic_density  # You can exclude particles with the '/' syntax
    MadDM> add direct_detection
    MadDM> add indirect_detection
    MadDM> add indirect_detection g g  # Computes loop-induced annihilation into gluons
    MadDM> add indirect_spectral_features  # Computes all loop-induced processes involving photons (aa, aZ, aH, etc.)
    MadDM> output MY_MADDM_PROCESS
    MadDM> launch

This will create a folder named ``MY_MADDM_PROCESS`` in your current working directory.
The command ``launch`` will execute the process, compute the relevant quantities and generate events (if applicable).

Install additional tools
========================

.. code-block:: text

    MadDM> install pythia8
    MadDM> install PPPC4DMID
    MadDM> install dragon
    MadDM> install dragon_data_from_galprop
    MadDM> quit

"""

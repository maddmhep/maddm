"""
===========================
Indirect Detection Analysis
===========================

Compute indirect detection signals and limits from dark matter annihilation.

Overview
========

Indirect detection looks for products of dark matter annihilation in astrophysical environment where the dark
matter is denser. For instance, typical benchmarks for gamma-ray searches are the dSphs or the Galactic Center.

Getting Started
===============

To include indirect detection in your MadDM workflow, simply add the appropriate module after generating your process.

Basic Example
=============

.. code-block:: text

    MadDM> import model MyDMmodel
    MadDM> define darkmatter chi
    MadDM> generate indirect_detection
    MadDM> output MyDMmodel
    MadDM> launch MyDMmodel

This setup computes the total annihilation cross sections for all available final states and stores the results in the output directory ``ID_example``.

Adding Loop-Induced Processes
=============================

To compute annihilation channels that are loop-induced (e.g., dark matter annihilation into gluons or photons), you can explicitly specify the final state:

.. code-block:: text

    MadDM> add indirect_detection g g
    MadDM> add indirect_spectral_features

These commands include:
- Annihilation into gluons (via loops)
- All relevant annihilations involving photons: gamma gamma, gamma Z, gamma h, etc.

Computing Spectra
=================

Once loop-induced and tree-level annihilation channels are generated, MadDM can compute the energy spectra of the final state particles:

.. code-block:: text

    MadDM> launch
    MadDM> plot spectrum photon

This will show the photon spectrum and save it as a PDF in the output folder.

Available Final States
----------------------

You can request spectra for:

- photon
- electron
- positron
- neutrino
- antiproton

Example:

.. code-block:: text

    MadDM> plot spectrum positron

Using External Tools
====================

MadDM supports integration with external tools for more realistic propagation models and astrophysical assumptions.

To install them:

.. code-block:: text

    MadDM> install PPPC4DMID
    MadDM> install dragon
    MadDM> install dragon_data_from_galprop

This will allow you to go beyond prompt spectra and compute propagated signals for cosmic rays.

Output Files
============

MadDM will generate:

- ``indirect_detection_results.txt``: total and partial annihilation cross sections.
- ``spectrum_<channel>.dat``: differential spectra for the specified final state.
- ``spectrum_plot_<channel>.pdf``: plot of the spectrum.

Exiting MadDM
=============

.. code-block:: text

    MadDM> quit
"""

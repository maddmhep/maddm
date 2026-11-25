"""
===========================
Direct Detection Tutorial
===========================

Compute direct detection signal and limits for both **nuclear** and **electronic** recoils.

----------------------------------------------------
1. Nuclear Recoils with RAPIDD and LZ/PICO-60 limits
----------------------------------------------------

For nuclear recoil interactions, MadDM compares the model with the upper limits from **LZ 2024** and **PICO-60** (and, with
``xlzd_limits`` in the ``maddm_card.dat``, with the projected XLZD sensitivity), and uses the
`RAPIDD <https://github.com/cheekyparticle/RAPIDD_for_DM>`_ module to evaluate the differential recoil spectra.

The LZ limits are published for reference models (an isoscalar SI coupling, SD couplings to neutrons only or to protons
only), while a model couples to protons and neutrons, whose amplitudes interfere in the nucleus. MadDM therefore compares
with the LZ limits the **equivalent cross sections** ``SigmaN_SI_eq``, ``SigmaN_SD_n_eq`` and ``SigmaN_SD_p_eq``: the cross
section of the reference model of each limit that gives the same expected number of events in LZ as the full model
(SI and SD couplings, protons and neutrons, with the chiral EFT structure factors and their two-body currents).
Comparing the proton and neutron cross sections separately with the limits would ignore the interference: it can exclude
models whose xenon signal cancels, or (for SI couplings to protons only) be too aggressive by :math:`(A/Z)^2 \\simeq 6`,
because the SI limit is isoscalar. See the direct detection section of the documentation for the details and the limits of
validity (the comparison is exact when the spectrum has the same shape as the reference one).


Example with a simplified scalar-mediated DM model
--------------------------------------------------

First, import the model and define the dark matter particle. Generally, for DMsimp models, you can choose between Xr (real scalar DM), Xc (complex scalar DM)
and Xd (Dirac spinor DM):

.. code-block:: text

    MadDM> import model DMsimp_s_spin0
    MadDM> define darkmatter xd

Then, generate the direct detection processes with:

.. code-block:: text

    MadDM> generate direct

At this point, you can set the parameter values of the model, stored in ``param_card_orig.dat``, by pressing 1, but this can also be done when launching the process.
Next, create the process folder. In this case, we call it ``DD_NR_spin0``, but you can choose any name you like. Then launch the process:

.. code-block:: text

    MadDM> output DD_NR_spin0
    MadDM> launch DD_NR_spin0

Select to run the ``direct`` module (but it should be already ON) by entering:

.. code-block:: text

    MadDM> direct = ON

If your model has both DM-nucleon and DM-electron interactions, you can set both ``direct_nucleon`` and ``direct_electron`` modules to ON,
and you will receive results for both types of interactions considered separately.

You can set the model parameters in the ``param_card.dat`` by pressing 7, or by editing the file directly in the process folder (in this case ``DD_NR_spin0/Cards/param_card.dat``).

You can set the MadDM parameters (such as the halo parameters, recoil energy range and much more) in the ``maddm_card.dat`` by pressing 8,
or by editing the file directly in the process folder (in this case ``DD_NR_spin0/Cards/maddm_card.dat``).

Once you are all set, you can run the process by pressing Enter. This module will first evaluate the DM-nucleon cross sections and the
equivalent cross sections and compare them against the experimental limits from LZ 2024 and PICO-60.
Next, with ``direct_mode = limit_and_rates``, it will use RAPIDD to compute the differential nuclear recoil rates.
After running the above commands, the output will look like:

.. code-block:: text

    INFO: MadDM Results
    INFO: Local Galactic halo (SHM): rho = 0.3 GeV/cm^3, vMP = 238 km/s, vescape = 544 km/s
    INFO:
    ****** Direct detection [cm^2]:
    INFO: SigmaN_*_eq: DM-nucleon cross sections equivalent to the model in LZ (same expected signal as the reference model of each LZ limit; RAPIDD chiral EFT, p-n interference and two-body currents included)
    INFO: SigmaN_SI_eq        All DM = 6.41e-48       ALLOWED     LZ2024 ul        = 6.48e-47
    INFO: SigmaN_SD_n_eq      All DM = 9.84e-43       ALLOWED     LZ2024 ul        = 7.59e-42
    INFO: SigmaN_SD_p_eq      All DM = 3.03e-41       ALLOWED     LZ2024 ul        = 2.33e-40
    INFO: SigmaN_SD_p         All DM = 0.00e+00       ALLOWED     Pico60 (2019) ul = 4.16e-41
    INFO: sigma_p and sigma_n (information only: compared separately with the LZ limits they ignore the p-n interference; SD_p is compared with PICO-60 above):
    INFO: SigmaN_SI_p         All DM = 6.16e-48
    INFO: SigmaN_SI_n         All DM = 6.14e-48
    INFO: SigmaN_SD_n         All DM = 0.00e+00
    INFO: Neutrino fog (xenon): reference boundary, NOT an exclusion limit. https://github.com/cajohare/NeutrinoFog; C. A. J. O'Hare, Phys. Rev. Lett. 127, 251802 (2021), arXiv:2109.03116
    INFO: Rescaled from rho = 0.3 GeV/cm^3 (published tables) to rhoDM = 0.3 GeV/cm^3.
    INFO: neutrino_fog_SI = 4.3291631542e-48 # cm^2; -1 means unavailable
    INFO: neutrino_fog_SD_p = 1.9474339879e-40 # cm^2; -1 means unavailable
    INFO: neutrino_fog_SD_n = 1.0584409940e-42 # cm^2; -1 means unavailable
    INFO: Recoil spectra: dRdE_nucleon.dat in /Users/yourname/yourprocessfolder/DD_NR_spin0/output/run_01
    INFO:
    INFO: Results written in: /Users/yourname/yourprocessfolder/DD_NR_spin0/output/run_01/MadDM_results.txt
    quit

where ``All DM`` refers to the DM-nucleon cross sections in case the DM particle contributes to the entire dark matter density of the universe
(default model parameters, with a DM mass of 10 GeV).
In the last column, you can see the upper limits (ul) from the LZ 2024 and PICO-60 experiments, which are compared against the
equivalent cross sections (LZ) and the SD proton cross section (PICO-60), returning the status of the model (ALLOWED or EXCLUDED).
This model has only a scalar (SI) coupling: ``SigmaN_SD_n_eq`` and ``SigmaN_SD_p_eq`` are not zero because they are the SI signal
of the model seen through the SD reference models (LZ sees the sum of the SI and SD signals), and they give approximately the same verdict as
``SigmaN_SI_eq`` (here sigma_eq/limit = 0.10 for SI and 0.13 for SD: different LZ curves and spectral shapes). The proton and neutron cross sections are given for information, without a verdict.
The xenon neutrino fog is a reference boundary, not an exclusion limit. ``DD_limits.pdf`` in the run folder shows the equivalent cross sections
of the run against the LZ limits, the SD proton cross section against the PICO-60 limit, and the neutrino fog as a function of the DM mass.

Results
-------

The results of each run will be stored inside the process folder you created, whose path is shown in the last line of the output.
Each time you launch the same process, MadDM creates a new folder ``output/run_XX`` within the process folder.

In this folder, you will find the cross sections in a file named ``MadDM_results.txt``.

With ``direct_mode = limit_and_rates`` (the default; type ``set limit`` to compute only the exclusion limits) you can also find
the differential recoil rates in the ``dRdE_nucleon.dat`` file, which contains the differential recoil spectra vs energy
for DM-nucleon interactions, one column per target chosen with ``direct_nucleon_targets`` in the ``maddm_card.dat``
(one or more of ``xe``, ``ge``, ``ar``, ``f``; default ``xe``). The spectra are computed with the chiral EFT structure factors
(``dd_nuclear_responses = chiral_eft``, the default; ``nreft`` gives the one-body NREFT responses with the
:math:`\\mathcal O_1` and :math:`\\mathcal O_4` couplings only). The file used below was produced with ``set direct_nucleon_targets xe, ge, ar``.
Below is an example of how to plot the differential rates:
"""
# sphinx_gallery_thumbnail_path = '_static/thumbnail/plot_direct_detection_thumb.png'
# %%
import matplotlib.pyplot as plt
import numpy as np

def plot_dd_rates(dd_file):
    """
    Plots the direct detection rates of all the targets in dRdE_nucleon.dat.

    Parameters:
    dd_file (str): Path to the dRdE_nucleon.dat file.
    """
    # the header gives the targets: E_R [keV] drde_<target> [events/kg/day/keV] ...
    with open(dd_file) as file:
        targets = [name[len('drde_'):] for name in file.readline().split() if name.startswith('drde_')]
    data = np.loadtxt(dd_file)

    # Plotting
    plt.figure(figsize=(10, 6))
    for i, target in enumerate(targets):
        rates = data[:, i + 1]
        mask = rates > 0  # the rate vanishes above the kinematic endpoint of each target
        plt.plot(data[mask, 0], rates[mask], label=target)
    plt.xlabel('Recoil Energy $E_R$ [keV]')
    plt.ylabel('dR/dE [$\\mathrm{events / kg / day / keV}$]')
    plt.title('DD nucleon Differential Rates')
    plt.yscale('log')  # y-axis in log scale
    plt.legend()
    plt.show()

plot_dd_rates('./plot_data/dRdE_nucleon.dat')

# %%
# ---------------------------------------------------------
# 2. Electronic Recoils with XENONnT and XENON10 limits
# ---------------------------------------------------------
#
# MadDM also supports **electronic recoil** constraints: the reference DM-electron cross section is compared with the
# **XENONnT** S2-only upper limit and with the **XENON10** S2-only data. These are particularly relevant for models with sub-GeV DM that interact with electrons.
# With this module you can compute the DM-electron recoil rates and constrain the Lagrangian parameters, or the cross section.
#
# Example with a simplified scalar-mediated DM model
# --------------------------------------------------
#
# First, import the model and define the dark matter particle. Generally, for DMsimp models, you can choose between Xr (real scalar DM), Xc (complex scalar DM)
# and Xd (Dirac spinor DM). Then, generate the direct detection processes with:
#
# .. code-block:: text
#
#     MadDM> import model eDMsimp_s_spin0
#     MadDM> define darkmatter xd
#     MadDM> generate direct
#
# At this point, you can set the parameter values of the model, stored in ``param_card_orig.dat``, by pressing 1, but this can also be done when launching the process.
# Next, create the process folder. In this case, we call it ``DD_ER_spin0``, but you can choose any name you like. Then launch the process:
#
# .. code-block:: text
#
#     MadDM> output DD_ER_spin0
#     MadDM> launch DD_ER_spin0
#
# Select to run the ``direct_electron`` module (but it should be already ON) by entering:
#
# .. code-block:: text
#
#     MadDM> direct_electron = ON
#
# If your model has both DM-nucleon and DM-electron interactions, you can set both ``direct_nucleon`` and ``direct_electron`` modules to ON,
# and you will receive results for both types of interactions considered separately.
#
# You can set the model parameters in the ``param_card.dat`` by pressing 7, or by editing the file directly in the process folder (in this case ``DD_ER_spin0/Cards/param_card.dat``).
#
# You can set the MadDM parameters (such as the halo parameters, recoil energy range and much more) in the ``maddm_card.dat`` by pressing 8,
# or by editing the file directly in the process folder (in this case ``DD_ER_spin0/Cards/maddm_card.dat``).
#
# Note that in the ``direct_electron`` module, if the DM particle has a mass above 1 GeV, it will be ignored by default.
# To change this behavior, set ``direct_electron_mode`` to ``always`` in the ``maddm_card.dat`` file, which you can do by pressing 8 after launching the process.
#
# Once you are all set, you can run the process by pressing Enter. This module will compute the reference cross section :math:`\bar\sigma_e`
# and compare it with the XENONnT limit, both as published and rescaled to the RAPIDD atomic responses; with ``direct_mode = limit_and_rates``
# it will also compute the DM-electron recoil rates and compare them with the XENON10 S2-only data.
#
# After running the above commands, the output will look like:
#
# .. code-block:: text
#
#     INFO: MadDM Results
#     INFO: Local Galactic halo (SHM): rho = 0.3 GeV/cm^3, vMP = 238 km/s, vescape = 544 km/s
#     INFO:
#     ****** Direct detection - electron recoil [cm^2]:
#     INFO: Sigma_e             All DM = 3.23e-43       ALLOWED     XENONnT ul       = 1.54e-41
#     INFO: Sigma_e (rescaled)  All DM = 3.23e-43       ALLOWED     XENONnT ul       = 3.01e-41
#     INFO: Sigma_e             All DM = 3.23e-43       ALLOWED     XENON10 ul       = 4.73e-38
#     INFO: Recoil spectra: dRdE_electron.dat, dRdS2_electron_XENON10.dat, dRdS2_electron_XENON1T.dat in /Users/yourname/yourprocessfolder/DD_ER_spin0/output/run_01
#     INFO:
#     INFO: Results written in: /Users/yourname/yourprocessfolder/DD_ER_spin0/output/run_01/MadDM_results.txt
#
# where ``All DM`` refers to the reference DM-electron cross section in case the DM particle contributes to the entire dark matter density of the universe,
# ``ALLOWED`` indicates that the model is consistent with the corresponding upper limit, shown in the last column. The XENON10 limit is the cross section
# for which the p-value of the most constraining S2 bin falls below 0.1; the p-value itself is written in ``MadDM_results.txt``.
#
# Results
# -------
#
# The results of each run will be stored inside the process folder you created, whose path is shown in the last line of the output.
# Each time you launch the same process, MadDM creates a new folder ``output/run_XX`` within the process folder.
#
# In this folder, you will find the results, such as the cross section and p-values, in a file named ``MadDM_results.txt``.
#
# The raw output of the Fortran module of MadDM, such as the DM mass and the effective DM-electron couplings, is in the ``maddm.out`` file,
# and the XENON10 signal in each S2 bin is in ``signal_e_recoil.dat``.
#
# With ``direct_mode = limit_and_rates`` (the default) you can also find the differential recoil spectra vs energy in ``dRdE_electron.dat``
# (for each target of ``direct_electron_targets``, ``xe`` and/or ``ar``, the rate summed over the shells and then the rate of each shell)
# and, for Xenon, the number of events vs scintillation signal in ``dRdS2_electron_XENON10.dat`` and ``dRdS2_electron_XENON1T.dat``.
# With ``direct_mode = limit`` (``set limit``) only sigma_e and the XENONnT limit are computed: the XENON10 and XENON1T results need the Xenon rates.
# The files used below were produced with ``set direct_electron_targets xe, ar``.
# Below is an example of how to plot the differential rates and number of events vs scintillation signal for electronic recoils:

# %%
import numpy as np
import matplotlib.pyplot as plt

def plot_dRdE_er(file_path):
    """
    Read and plot the dR/dE spectrum (all shells) of each target in dRdE_electron.dat.

    Parameters:
    file_path (str): Path to the dRdE_electron.dat file.
    """
    # header: k [keV] E [keV] drde_<target> [events/kg/day/keV] ... drde_<target>_<shell> ...
    with open(file_path) as file:
        names = [name for name in file.readline().split() if name.startswith('drde_')]
    targets = [name[len('drde_'):] for name in names if '_' not in name[len('drde_'):]]
    data = np.loadtxt(file_path)

    # Plot (points with rate <= 0 are not shown on the log scale)
    plt.figure(figsize=(8, 5))
    for i, target in enumerate(targets):
        rates = data[:, 2 + i]
        mask = rates > 0
        plt.plot(data[mask, 1], rates[mask], linestyle='-', label='dR/dE ' + target)
    plt.yscale('log')
    plt.xlabel("Recoil Energy (keV)")
    plt.ylabel("dR/dE [events/kg/day/keV]")
    plt.title("Differential Recoil Rate")
    plt.legend()
    plt.show()

plot_dRdE_er('./plot_data/dRdE_electron.dat')

# %%
# Plot the number of events vs scintillation signal for Electronic Recoils
import numpy as np
import matplotlib.pyplot as plt

def plot_dRdS2_er(xenon10_file, xenon1t_file):
    """
    Read and plot the dN/dS2 spectra (second column) from the XENON10 and XENON1T files.
    """
    plt.figure(figsize=(8, 5))
    for path, label, color in ((xenon10_file, 'Xenon10', 'red'), (xenon1t_file, 'Xenon1T', 'blue')):
        s2, rates = np.loadtxt(path, usecols=(0, 1)).T
        mask = rates > 0
        plt.plot(s2[mask], rates[mask], linestyle='-', color=color, label=label)
    plt.yscale('log')
    plt.xlabel("S2")
    plt.ylabel("dN/dS2 [events]")
    plt.xlim(14, 271)
    plt.title("Differential Rate vs S2")
    plt.legend()
    plt.show()

plot_dRdS2_er("./plot_data/dRdS2_electron_XENON10.dat", "./plot_data/dRdS2_electron_XENON1T.dat")

# %%
#Exiting MadDM
#=============
#
#.. code-block:: text
#
#    MadDM> quit
==============================
Welcome to MadDM documentation
==============================

MadDM is a numerical tool designed to compute dark matter relic abundance, dark matter nucleus scattering rates, and dark matter indirect detection predictions in a generic model. Based on the existing MadGraph architecture, MadDM is easily integrable into any MadGraph collider study. The tool features a simple Python interface, maintaining the user-friendliness characteristic of MadGraph without sacrificing functionality.

Install MadDM
=============

**MadDM** is a plugin for ``MG5_aMC@NLO`` and for ``MadGraph7`` (**MadGraph**).

Requirements
------------

- One of the following MadGraph versions:

  - `MG5_aMC@NLO <https://launchpad.net/mg5amcnlo/+download>`__ v3.8.x, with Python 3.8 or higher;
  - `MadGraph7 <https://github.com/MadGraphTeam/MadGraph7/releases>`__ v0.2.x, with Python 3.12 or higher.

- Python packages: `numpy <https://numpy.org/>`__, `scipy <https://scipy.org/>`__ and `matplotlib <https://matplotlib.org/>`__.
  Numerical calculations require them; model import, definitions, diagram generation and process export do not.
  Use ``install py_requirements`` at the MadDM prompt to install them, together with Meson and Ninja, in the Python environment running MadGraph.
  MadDM uses Conda only for a Conda-managed interpreter, otherwise pip through the running Python; neither is needed if you install the packages yourself.
  For indirect detection, the installation choice is shown before the Pythia 8 setup, without an answer timeout; other numerical launches are checked at the main prompt, and missing packages give an actionable message.
- A Fortran compiler and ``make`` for the generated processes.
- Optional: `PyMultiNest <https://johannesbuchner.github.io/PyMultiNest/>`__ with the MultiNest library for nested-sampling scans.

Further tools are installed on demand (see `Third-party packages`_ below):

- Pythia 8 needs a C++ compiler, ``make`` and ``rsync``. MadDM supplies a verified ``zlib`` archive when needed.
- `RAPIDD <https://github.com/cheekyparticle/RAPIDD_for_DM>`__, used for nuclear and electron recoils, needs Git, CMake 3.13+, a C/C++ toolchain and GNU GSL.
  MadDM uses an existing GSL that CMake can find, or builds GSL under ``HEPTools`` when it is missing.
  On macOS, build RAPIDD and GSL for the same architecture as the Python running MadDM.
- With Python 3.12 or later, ``sigmav_method = reshuffling`` requires Meson and the Ninja build executable for NumPy F2PY.
  MadDM reuses executables on ``PATH``; they are part of the first numerical indirect-detection dependency check and of ``install py_requirements``.

First-time downloads require internet access.

Installation
------------

1. Download and unpack MadGraph. For MG5_aMC@NLO, take a v3.8.x tarball from the `download page <https://launchpad.net/mg5amcnlo/+download>`__; for MadGraph7, take a v0.2.x tarball from the `releases page <https://github.com/MadGraphTeam/MadGraph7/releases>`__. For example:

   .. code-block:: bash

     tar xzf MG5_aMC_v3.8.1.tar.gz
     cd MG5_aMC_v3_8_1

2. Start MadGraph, with ``./bin/mg5_aMC`` for MG5_aMC@NLO or ``./bin/madgraph`` for MadGraph7.

3. Install MadDM from the MadGraph prompt and then quit MadGraph:

   .. code-block:: text

     MG5_aMC>install maddm
     MG5_aMC>quit

   For a manual installation, for example to use the development branch, move any previous MadDM copy **outside** the MadGraph ``PLUGIN/`` directory and clone MadDM there instead:

   .. code-block:: bash

     git clone -b dev https://github.com/maddmhep/maddm.git PLUGIN/maddm

4. Start MadDM from the MadGraph root directory:

   .. code-block:: bash

     ./bin/mg5_aMC --mode=maddm     # MG5_aMC@NLO
     ./bin/madgraph --mode=maddm    # MadGraph7

5. You are now set, have fun! MadDM offers to install the remaining dependencies when a process first needs them (see below).

When updating MadDM, export old processes again under a new directory name and reuse their input cards (see :ref:`old-process-directories`).

Third-party packages
--------------------

The following packages can be installed with the ``install`` command at the MadDM prompt. Be sure to cite them as indicated during the install phase if you use them.

- ``install pythia8`` installs `Pythia 8 <https://pythia.org/>`__ with the MadGraph installer, together with `LHAPDF <https://www.lhapdf.org/index.html>`__, ``zlib`` and the MadGraph--Pythia interface (``mg5amc_py8_interface``).
- ``install cosmixs`` installs the tabulated `CosmiXs <https://github.com/ajueid/CosmiXs>`__ source spectra; select them with ``indirect_flux_source_method = CosmiXs`` in the run card.
  On the first indirect generation, MadDM offers to install both Pythia 8 and the CosmiXs source spectra.
- ``install RAPIDD`` downloads and compiles RAPIDD. MadDM also offers to do it automatically the first time you generate ``direct_detection``.
- ``install cosmixs_flux`` downloads the `CosmiXsPPPC <https://github.com/CosmiXsPPPC/CosmicRay-flux>`__ propagated fluxes.
  The first ``generate indirect_detection`` or ``add process indirect_detection`` downloads them automatically; use the command to prepare an existing process directory or to retry a failed download.
  The secondary fluxes are included with MadDM.
- ``install neutrino_fog`` downloads the xenon `neutrino-fog tables <https://github.com/cajohare/NeutrinoFog>`__.
  They are downloaded automatically when the direct-detection diagrams are generated (``generate`` or ``add direct_detection``).
- ``install py_requirements`` installs ``numpy``, ``scipy``, ``matplotlib``, Meson and Ninja; ``install reweight_tools`` installs only Meson and Ninja, also in an isolated pip environment under ``HEPTools``.

Launching an analysis does not install any of these dependencies.

Running relic density, direct and indirect detection together
-------------------------------------------------------------

After importing your model and defining its dark matter candidate, generate the processes for all the calculations you want in the same output:

.. code-block:: text

  generate relic_density
  add direct_detection
  add indirect_detection
  output dm_combined
  launch dm_combined

At the launch menu, ``1`` toggles relic density, ``2`` and ``3`` toggle direct detection for DM--nucleon (``direct_nucleon``) and DM--electron (``direct_electron``) scattering, and ``4`` cycles through the continuum indirect modes.
You can also enter ``relic = ON``, ``direct_nucleon = ON``, ``direct_electron = ON`` and ``indirect = sigmav`` for rates, or ``indirect = flux_earth`` for the propagated CosmiXsPPPC antiproton and antideuteron predictions.
``set limit`` computes only the direct detection cross sections and limits, ``set limit_and_rates`` (the default) also the recoil spectra (see the direct detection module below).
For a gamma-ray line analysis, also run ``add indirect_spectral_features`` before ``output``.

If the menu shows **No processes generated**, add the missing process at the main MadDM prompt and export a new output directory.
Use ``add`` to keep the processes already in memory; ``generate`` starts a new set.
Installing RAPIDD enables the nucleon recoil calculation but does not generate relic or indirect matrix elements.
Existing cards can be reused: the obsolete ``indirect_flux_earth_method`` setting is ignored when they are read.

See the :doc:`tutorials and examples <tutorial_and_examples/index>` for complete worked setups.

Overview
========

MadDM can calculate the dark matter relic abundance in models that include a multi-component dark sector, resonance annihilation channels, and co-annihilations.

Direct Detection Module
-----------------------

The direct detection module of MadDM computes spin-independent and spin-dependent cross-sections, as well as differential recoil rates, for both dark matter--nucleon and dark matter--electron scattering.
It then evaluates the statistical significance of the selected model against the most recent experimental data.

For dark matter--nucleon scattering, MadDM compares the model with the LZ 2024 and PICO-60 upper limits and, with ``xlzd_limits = 200ty`` or ``1000ty`` in the ``maddm_card``, with the projected XLZD sensitivity.
RAPIDD calculates the differential recoil rates for xenon, argon, germanium and fluorine targets.

The published limits are given for reference models: an isoscalar spin-independent coupling (:math:`f_p = f_n`), and spin-dependent couplings to neutrons only or to protons only.
A model couples in general to both nucleons, and a xenon detector sees the coherent sum of the proton and neutron amplitudes.
Comparing :math:`\sigma_p` and :math:`\sigma_n` separately with the limits ignores this interference and goes wrong in both directions:

- with destructive interference the signal can vanish while both cross sections are above the limits (spin independent: :math:`f_n/f_p \simeq -Z/(A-Z) \simeq -0.7` in xenon; spin dependent: e.g. an axial coupling to the :math:`u` quark only, for which the naive SD-neutron comparison is about 3 times too aggressive);
- the SI limit is isoscalar: for a coupling to protons only the xenon signal scales as :math:`Z^2\sigma_p` while the limit refers to :math:`A^2\sigma`, so the naive comparison is too aggressive by :math:`(A/Z)^2 \simeq 5.9`;
- in xenon the two-body currents make the SD response to a proton coupling much larger than the one-body estimate.

MadDM therefore uses, for each LZ limit, the **equivalent cross section** computed by RAPIDD (``rapidd.chiral_eft.equivalent_sigma``):

.. math::

   \sigma_{\rm eq} = \sigma_{\rm ref}\,\frac{N_{\rm model}}{N_{\rm ref}(\sigma_{\rm ref})},

where :math:`N` is the expected LZ signal (recoil spectrum times the LZ 2024 nuclear-recoil efficiency), computed with the full model (SI and SD couplings, protons and neutrons, all xenon isotopes, with the chiral EFT structure factors of `Hoferichter et al. 2019 <https://arxiv.org/abs/1812.05617>`__ for SI and `Hoferichter et al. 2020 <https://arxiv.org/abs/2007.08529>`__ for SD, which include the two-body currents) and with the reference model of each curve (Helm form factor and isoscalar coupling for SI, neutron-only or proton-only axial couplings for SD-n and SD-p, as in the LZ analysis), at the same DM mass and halo.
:math:`\sigma_{\rm eq}` does not depend on :math:`\sigma_{\rm ref}` nor on ``rhoDM``, and is compared directly with the published curve: ``SigmaN_SI_eq``, ``SigmaN_SD_n_eq`` and ``SigmaN_SD_p_eq`` give the verdicts of LZ (and of the XLZD projections, which use the same xenon target; their different efficiency is neglected).
A model with both SI and SD couplings has an SI equivalent cross section that includes the SD signal, and vice versa: LZ sees the sum.
The SD-n and SD-p rows are the same LZ data with different reference models, so they give the same verdict up to the relative normalization of the two curves (about 10%).
The comparison is exact when the recoil spectrum of the model has the same shape as the reference one, which holds to a good approximation for standard SI and SD couplings; for momentum- or velocity-dependent interactions and light mediators it is an approximation, still better than the separate comparison.
If the reference models have no signal in LZ (DM mass below threshold) there is no verdict.
The direct-detection likelihood of the nested-sampling scans (``directSI``, ``directSD_p``, ``directSD_n`` in the ``multinest_card``) uses the same equivalent cross sections: the half-Gaussian penalty applies when :math:`\sigma_{\rm eq}` is above the LZ limit (for ``directSD_p`` also when :math:`\sigma_{\rm SD,p}` is above the PICO-60 limit; the most constraining is used).

``SigmaN_SD_p`` is compared directly with the PICO-60 SD-proton limit (fluorine is dominated by the proton spin; PICO-60 uses one-body responses). This comparison is approximate: with the two-body currents and the proton-neutron interference the fluorine signal of an axial coupling is up to 1.5 (coupling to the :math:`u` quark only) or 2.6 (:math:`d` quark only) times smaller than the one implied by :math:`\sigma_p`, so the naive comparison can be too aggressive by these factors (no difference for isoscalar axial or tensor couplings).
``SigmaN_SI_p``, ``SigmaN_SI_n`` and ``SigmaN_SD_n`` are shown as information only, without a verdict; they and their limits (``lim_sigmaN_*``) stay in the scan output, together with ``sigmaN_SI_eq``, ``sigmaN_SD_n_eq``, ``sigmaN_SD_p_eq`` and their limits.
The equivalent cross sections always use the chiral EFT; the ``maddm_card`` parameter ``dd_nuclear_responses`` (``chiral_eft``, default, or ``nreft``: one-body NREFT responses, :math:`\mathcal O_1` and :math:`\mathcal O_4` only) selects only the nuclear responses of the saved recoil spectra (with ``nreft``: DMFormFactor with the shell-model density matrices of `arXiv:2208.07945 <https://arxiv.org/abs/2208.07945>`__ for xenon and fluorine, the fits of `arXiv:1203.3542 <https://arxiv.org/abs/1203.3542>`__ for germanium and argon).

For dark matter--electron scattering, the module evaluates the differential recoil spectrum in xenon and argon and, for xenon, the event rate as a function of the ionization (S2) signal for XENON10 and XENON1T, and compares the reference cross section :math:`\bar\sigma_e` with the XENONnT S2-only limit (heavy mediator, `arXiv:2601.11296 <https://arxiv.org/abs/2601.11296>`__).
That limit is computed by XENONnT with the ionization form factors of Essig et al.; with the atomic responses used by RAPIDD it would be about 2 times weaker.
MadDM shows both: ``Sigma_e`` with the published limit and ``Sigma_e (rescaled)`` with the limit rescaled to the RAPIDD atomic responses (an estimate, ``vendor/RAPIDD_for_DM/lib/experiments/XENONnT/XENONnT_DMe_SI_RAPIDD_FF.dat``; ``lim_sigma_e_rescaled`` in the scan output).
The cross section is also tested against the XENON10 S2-only data as reanalysed in `arXiv:1703.00910 <https://arxiv.org/abs/1703.00910>`__ (7 S2 bins, 15 kg day, no background subtraction): MadDM reports the p-value of the most constraining bin, :math:`P(n \le n_{\rm obs}\,|\,s)`, and the equivalent 90% CL upper limit on :math:`\bar\sigma_e`, ``Sigma_e (XENON10)`` (excluded when the p-value is below 0.1; ``lim_sigma_e_Xenon10`` and ``pvalue_Xenon10`` in the scan output).
The p-value refers to all DM at the local density of the card.

The ``maddm_card`` parameter ``direct_mode`` sets what the two modules compute.
With ``limit`` (``set limit`` at the launch prompt), MadDM computes only the cross sections and compares them with the upper limits, which is faster.
With ``limit_and_rates`` (the default, ``set limit_and_rates``), it also computes and saves the recoil spectra of the targets chosen with ``direct_nucleon_targets`` (one or more of ``xe``, ``ge``, ``ar``, ``f``) and ``direct_electron_targets`` (``xe`` and/or ``ar``), both ``xe`` by default.
The XENON10 and XENON1T results need the xenon DM--electron rates, so they are computed only with ``limit_and_rates`` and ``xe`` among the electron targets.
With ``limit_and_rates``, ``direct_modulation = True`` also computes the annual modulation of the DM--nucleon and DM--electron rates, in ``day_bins`` time bins of the year (see the :doc:`direct detection tutorial <tutorial_and_examples/plot_direct_detection>`).
See :ref:`plotting` for the output files.

MadDM also reports the xenon neutrino fog of `O'Hare (2021) <https://arxiv.org/abs/2109.03116>`__ at the DM mass, as a reference boundary rather than an exclusion limit.
With ``direct_nucleon``, ``DD_limits.pdf`` compares the equivalent cross sections of the run (and :math:`\sigma_{\rm SD,p}` for PICO-60) with the limits, the XLZD projection and the neutrino fog as a function of the DM mass (see :ref:`plotting`).
Above the highest tabulated mass, the DM--nucleon limits are extrapolated linearly in the DM mass up to :math:`10^{15}` GeV.

The upper limits on the cross sections (LZ 2024 SI, SD-neutron and SD-proton, PICO-60 SD-proton, XENONnT DM-electron), the XLZD projections and the neutrino fog assume a local dark matter density of :math:`0.3\,\mathrm{GeV\,cm^{-3}}`: they are rescaled by ``0.3/rhoDM``, with ``rhoDM`` the density of the ``maddm_card`` (see :ref:`galactic-halo`).

Indirect Detection Module
-------------------------

The indirect detection module of MadDM computes the velocity-averaged cross-section for dark matter particles annihilating into multiple final state particles. It further provides the energy spectra of photons, neutrinos, and cosmic rays generated by these final states after decaying, showering, and hadronization. The module automatically computes the flux of prompt neutrinos and gamma rays at detection.

In ``flux_earth`` mode, it computes antiproton and antideuteron fluxes with the CosmiXsPPPC propagated annihilation tables and bundled secondary spectra. Positron Earth fluxes are pending the corresponding inputs. A built-in interface with the **nested sampling** PyMultiNest algorithm enables efficient sampling of the model parameter space. The code also allows testing the model against the Fermi-LAT dwarf spheroidal galaxy likelihood.

MadDM is also able to compute loop-induced processes relevant to indirect detection phenomenology. Specifically, MadDM can automatically compute the annihilation of dark matter into final states like gamma-gamma, gamma-Z, gamma-H, gluon-gluon, and others. In generic dark matter models, these final states are loop-induced processes, where the tree-level amplitude is zero.

For the gamma-ray line final state, MadDM performs an analysis of the signal against Fermi-LAT gamma-ray line data to assess whether the model parameter space is allowed by current exclusion limits at 95% CL.

See :ref:`indirect-detection` for the details.

References
==========

1. **M. Backovic**, **K. Kong**, and **M. McCaskey** (2014), *MadDM v.1.0: Computation of Dark Matter Relic Abundance Using MadGraph5*.
   `arXiv:1308.4955 <https://arxiv.org/abs/1308.4955>`_.
   DOI: `10.1016/j.dark.2014.04.001 <https://doi.org/10.1016/j.dark.2014.04.001>`_.
   *Physics of the Dark Universe*, **5-6**, 18–28.

2. **M. Backović**, **A. Martini**, **O. Mattelaer**, **K. Kong**, and **G. Mohlabeng** (2015), *Direct Detection of Dark Matter with MadDM v.2.0*.
   `arXiv:1505.04190 <https://arxiv.org/abs/1505.04190>`_.
   DOI: `10.1016/j.dark.2015.09.001 <https://doi.org/10.1016/j.dark.2015.09.001>`_.
   *Phys. Dark Univ.*, **9-10**, 37–50.

3. **F. Ambrogi**, **C. Arina**, **M. Backovic**, **J. Heisig**, **F. Maltoni**, **L. Mantani**, **O. Mattelaer**, and **G. Mohlabeng** (2019), *MadDM v.3.0: a Comprehensive Tool for Dark Matter Studies*.
   `arXiv:1804.00044 <https://arxiv.org/abs/1804.00044>`_.
   DOI: `10.1016/j.dark.2018.11.009 <https://doi.org/10.1016/j.dark.2018.11.009>`_.
   *Phys. Dark Univ.*, **24**, 100249.

4. **C. Arina**, **J. Heisig**, **F. Maltoni**, **D. Massaro**, and **O. Mattelaer** (2023), *Indirect dark-matter detection with MadDM v3.2 – Lines and Loops*.
   `arXiv:2107.04598 <https://arxiv.org/abs/2107.04598>`_.
   DOI: `10.1140/epjc/s10052-023-11377-2 <https://doi.org/10.1140/epjc/s10052-023-11377-2>`_.
   *Eur. Phys. J. C*, **83** (3), 241.

Citation
========

If you are using this tool or the results obtained with it in a published work, please cite one or more of the previous references.
A ``CITATION.bib`` file is available in the root of the `MadDM repository <https://github.com/maddmhep/maddm>`__.

.. toctree::
   :hidden:

   tutorial_and_examples/index
   indirect_detection
   galactic_halo
   plotting
   contribute

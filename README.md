# MadDM

MadDM is a numerical tool designed to compute dark matter relic abundance, dark matter nucleus scattering rates, and dark matter indirect detection predictions in a generic model.  
Based on the existing MadGraph 5 architecture, MadDM is easily integrable into any MadGraph collider study. The tool features a simple Python interface, maintaining the user-friendliness characteristic of MadGraph 5 without sacrificing functionality.

## Installation

**MadDM** is a plugin for `MG5_aMC@NLO` (**MadGraph**). For a manual
installation, close MadGraph, move any previous MadDM copies **outside**
the MadGraph `PLUGIN/` directory, and place this `maddm/` folder in
`PLUGIN/`. Keep the new folder separate from older installations. From the
MadGraph root, start the plugin with:

```bash
./bin/mg5_aMC --mode=maddm
```

When updating, export old processes again under a new directory name and
reuse their input cards. See the [indirect detection guide](doc/indirect_detection.rst).

For a standard installation from MadGraph's package manager, enter
`install maddm` at the MadGraph prompt.

### Requirements

- At least MadGraph v3.6 is required to run the latest version of MadDM.
- MadDM is compatible with the LTS version of MadGraph.
- MadDM requires Python 3.8+, `six`, `numpy`, `scipy`, and `matplotlib`.
  At startup, MadDM reports missing packages and the command to install
  them in the Python environment running MadGraph. Interactive launches
  offer to run the command and restart MadDM. Batch launches print the
  command and exit. MadDM uses Conda when it manages the active Python
  and its executable is available, or pip with the active interpreter.
  Generated processes require a Fortran compiler and `make`. Install
  PyMultiNest with the MultiNest library for optional nested scans.
- On first indirect generation, MadDM offers to install Pythia8 and the
  CosmiXs source spectra. Pythia8 uses the MadGraph installer and needs a
  C++ compiler, `make`, and `rsync`. MadDM supplies a verified
  zlib archive when needed. CosmiXsPPPC propagated fluxes are downloaded
  separately, and the secondary flux tables are included with MadDM.
- With Python 3.12 or later, `sigmav_method = reshuffling` requires Meson
  and the Ninja build executable for NumPy F2PY. MadDM reuses executables
  on `PATH`. It installs missing tools with Conda in the active Conda
  environment, or with pip in an isolated environment under MadGraph's
  HEPTools directory. To install them separately, run `install
  reweight_tools` at the MadDM prompt.
- Nuclear and electron recoil calculations use RAPIDD. Its first build
  requires Git, CMake 3.13+, a C/C++ toolchain, and GNU GSL. MadDM uses an
  existing GSL that CMake can find or builds GSL under HEPTools when it is
  missing. First-time downloads require internet access.


### Running relic density, direct and indirect detection together

After importing your model and defining its dark matter candidate, generate
the processes for all the calculations you want in the same output:

```text
generate relic_density
add direct_detection
add indirect_detection
output dm_combined
launch dm_combined
```

At the launch menu, `1` toggles relic density and `4` cycles through the
continuum indirect modes. You can also enter `relic = ON` and
`indirect = sigmav` for rates, or `indirect = flux_earth` for the propagated
CosmiXsPPPC antiproton and antideuteron predictions. For a gamma-ray line
analysis, also run `add indirect_spectral_features` before `output`.

If the menu shows **No processes generated**, add the missing process at
the main MadDM prompt and export a new output directory. Use `add` to keep
the processes already in memory; `generate` starts a new set. Installing
RAPIDD enables the nucleon recoil calculation but does not generate relic
or indirect matrix elements. Existing cards can be reused: the obsolete
`indirect_flux_earth_method` setting is ignored when they are read.

## Overview

MadDM can calculate the dark matter relic abundance in models that include a multi-component dark sector, resonance annihilation channels, and co-annihilations.

### Direct Detection Module

The Direct Detection module of MadDM computes spin-independent and spin-dependent cross-sections, as well as differential recoil rates, for both dark matter–nucleon and dark matter–electron scattering. It then evaluates the statistical significance of the selected model against the most recent experimental data.

For dark matter–nucleon scattering, the RAPIDD code is used to calculate the differential recoil rates for xenon, argon, and germanium targets. Subsequently, it computes the expected signal in xenon and determines the p-value of the input model by employing the likelihood function derived from the 2024 LZ results.

For dark matter–electron scattering, the module evaluates both the differential recoil spectrum and the event rate as a function of the ionization (S2) signal in xenon for XENON10 and XENON1T, and compares the reference cross section $\bar\sigma_e$ with the XENONnT S2-only limit (heavy mediator, [arXiv:2601.11296](https://arxiv.org/abs/2601.11296)). That limit is computed by XENONnT with the ionization form factors of Essig et al.; with the atomic responses used by RAPIDD it would be about 2 times weaker. MadDM shows both: `Sigma_e` with the published limit and `Sigma_e (rescaled)` with the limit rescaled to the RAPIDD atomic responses (an estimate, `vendor/RAPIDD_for_DM/lib/experiments/XENONnT/XENONnT_DMe_SI_RAPIDD_FF.dat`; `lim_sigma_e_rescaled` in the scan output). The cross section is also tested against the XENON10 S2-only data as reanalysed in [arXiv:1703.00910](https://arxiv.org/abs/1703.00910) (7 S2 bins, 15 kg day, no background subtraction): MadDM reports the p-value of the most constraining bin, $P(n \le n_{\rm obs}\,|\,s)$, and the equivalent 90% CL upper limit on $\bar\sigma_e$, `Sigma_e (XENON10)` (excluded when the p-value is below 0.1; `lim_sigma_e_Xenon10` and `pvalue_Xenon10` in the scan output). The p-value refers to all DM at the local density of the card.

The upper limits on the cross sections (LZ 2024 SI and SD-neutron, PICO-60 SD-proton, XENONnT DM-electron) assume a local dark matter density of 0.3 GeV/cm$^3$: they are rescaled by 0.3/`rhoDM`, with `rhoDM` the density of the `maddm_card`.

### Indirect Detection Module

The indirect detection module of MadDM computes the velocity-averaged cross-section for dark matter particles annihilating into multiple final state particles. It further provides the energy spectra of photons, neutrinos, and cosmic rays generated by these final states after decaying, showering, and hadronization. The module automatically computes the flux of prompt neutrinos and gamma rays at detection.

In `flux_earth` mode, it computes antiproton and antideuteron fluxes with the CosmiXsPPPC propagated annihilation tables and bundled secondary spectra. Positron Earth fluxes are pending the corresponding inputs. A built-in interface with the **nested sampling** PyMultiNest algorithm enables efficient sampling of the model parameter space. The code also allows testing the model against the Fermi-LAT dwarf spheroidal galaxy likelihood.

MadDM is also able to compute loop-induced processes relevant to indirect detection phenomenology. Specifically, MadDM can automatically compute the annihilation of dark matter into final states like gamma-gamma, gamma-Z, gamma-H, gluon-gluon, and others. In generic dark matter models, these final states are loop-induced processes, where the tree-level amplitude is zero.

For the gamma-ray line final state, MadDM performs an analysis of the signal against Fermi-LAT gamma-ray line data to assess whether the model parameter space is allowed by current exclusion limits at 95% CL.

### Documentation

The Sphinx documentation covers [indirect detection](doc/indirect_detection.rst),
the [shared Galactic halo](doc/galactic_halo.rst), and
[run outputs and scan plots](doc/plotting.rst). See the
[MadDM documentation site](https://maddmhep.github.io/maddm/dev/) for
the tutorials, or build them with `make -C doc html`.

## References

1.  **M. Backovic**, **K. Kong**, and **M. McCaskey** (2014), *MadDM v.1.0: Computation of Dark Matter Relic Abundance Using MadGraph5*.
    [arXiv:1308.4955](https://arxiv.org/abs/1308.4955).
    DOI: [10.1016/j.dark.2014.04.001](https://doi.org/10.1016/j.dark.2014.04.001).
    *Physics of the Dark Universe*, **5-6**, 18--28.
2.  **M. Backović**, **A. Martini**, **O. Mattelaer**, **K. Kong**, and **G. Mohlabeng** (2015), *Direct Detection of Dark Matter with MadDM v.2.0*.
    [arXiv:1505.04190](https://arxiv.org/abs/1505.04190).
    DOI: [10.1016/j.dark.2015.09.001](https://doi.org/10.1016/j.dark.2015.09.001).
    *Phys. Dark Univ.*, **9-10**, 37--50
3.  **F. Ambrogi**, **C. Arina**, **M. Backovic**, **J. Heisig**, **F. Maltoni**, **L. Mantani**, **O. Mattelaer**, and **G. Mohlabeng** (2019), *MadDM v.3.0: a Comprehensive Tool for Dark Matter Studies*.
    [arXiv:1804.00044](https://arxiv.org/abs/1804.00044).
    DOI: [10.1016/j.dark.2018.11.009](https://doi.org/10.1016/j.dark.2018.11.009).
    *Phys. Dark Univ.*, **24**, 100249.
4.  **C. Arina**, **J. Heisig**, **F. Maltoni**, **D. Massaro**, and **O. Mattelaer** (2023), *Indirect dark-matter detection with MadDM v3.2 -- Lines and Loops*.
    [arXiv:2107.04598](https://arxiv.org/abs/2107.04598).
    DOI: [10.1140/epjc/s10052-023-11377-2](https://doi.org/10.1140/epjc/s10052-023-11377-2).
    *Eur. Phys. J. C*, **83** (3), 241.

Some third-party packages can be installed using the `install` command in the shell. Be sure to cite them as indicated during the install phase if you use them.

Nuclear and electron recoil rates use [RAPIDD](https://github.com/cheekyparticle/RAPIDD_for_DM). MadDM offers to download and compile it when you first generate `direct_detection`. You can also run `install RAPIDD` at the MadDM prompt.
On macOS, build RAPIDD and GSL for the same architecture as the Python running MadDM.

To install the tabulated production spectra separately, run `install cosmixs`
and select `indirect_flux_source_method = CosmiXs` in the run card.
The first `generate indirect_detection` or `add process indirect_detection`
downloads the two annihilation archives from
[CosmiXsPPPC](https://github.com/CosmiXsPPPC/CosmicRay-flux).
MadDM checks the contents and caches 54 compressed tables in
`ExpData/CR_flux`. Run `install cosmixs_flux` to prepare an existing
process directory or to retry a failed download. The release still
includes the secondary fluxes. Launching an analysis does not install
any of these dependencies.

## Citation

If you are using this tool or the results obtained with it in a published work, please cite one or more of the previous references.
A [CITATION.bib](./CITATION.bib) file is also available.

# MadDM

MadDM is a numerical tool designed to compute dark matter relic abundance, dark matter nucleus scattering rates, and dark matter indirect detection predictions in a generic model.
Based on the existing MadGraph architecture, MadDM is easily integrable into any MadGraph collider study. The tool features a simple Python interface, maintaining the user-friendliness characteristic of MadGraph without sacrificing functionality.

**Documentation:** <https://maddmhep.github.io/maddm/dev/> (installation, tutorials, indirect detection, Galactic halo, run outputs and plots).

## Quick installation

MadDM is a plugin for MadGraph. It requires either

- [MG5_aMC@NLO](https://launchpad.net/mg5amcnlo/+download) v3.8.x, with Python 3.8 or higher, or
- [MadGraph7](https://github.com/MadGraphTeam/MadGraph7/releases) v0.2.x, with Python 3.12 or higher,

together with a Fortran compiler and `make`. The numerical calculations need `numpy`, `scipy` and `matplotlib`:
run `install py_requirements` at the MadDM prompt to install them in the Python environment running MadGraph.

At the MadGraph prompt (`./bin/mg5_aMC` or `./bin/madgraph`), run `install maddm`, quit, and start MadDM from the MadGraph root:

```bash
./bin/mg5_aMC --mode=maddm     # MG5_aMC@NLO
./bin/madgraph --mode=maddm    # MadGraph7
```

For a manual installation, move any previous MadDM copy outside the MadGraph `PLUGIN/` directory and place this `maddm/` folder in `PLUGIN/`.
MadDM offers to install Pythia 8 and the CosmiXs spectra, and to compile RAPIDD, when a process first needs them; this can also be done with `install pythia8`, `install cosmixs` and `install RAPIDD` at the MadDM prompt. RAPIDD is shipped with MadDM (`vendor/RAPIDD_for_DM`): use the release archives or `git clone --recursive` (the GitHub "Source code" archives do not contain it).
See the [installation guide](https://maddmhep.github.io/maddm/dev/index.html#install-maddm) for the full list of requirements and options.

## Direct detection

Direct detection uses [RAPIDD](https://github.com/cheekyparticle/RAPIDD_for_DM) (`vendor/RAPIDD_for_DM`):

- DM-nucleon scattering: LZ 2024 (equivalent cross sections with the chiral EFT), PICO-60 and the XLZD projections;
  recoil spectra in xenon, germanium, argon and fluorine (`direct_nucleon_targets`).
- DM-electron scattering: XENONnT (published and rescaled to the RAPIDD atomic responses), DAMIC-M (published and
  rescaled to the RAPIDD silicon response) and DarkSide-50 limits; XENON10 and XENON1T S2-only data, weaker but the
  only tests computed from the recoil spectrum of the model (any operator); recoil spectra in xenon, argon and silicon
  (`direct_electron_targets`; silicon with the QCDark2 RPA dielectric function, also per number of electron-hole pairs).
- The rates use the heavy-mediator (contact) approximation: MadDM identifies the mediators of the direct detection
  diagrams and warns (flags `dd_light_mediator_warning`, `dd_e_light_mediator_warning`) when a t-channel mediator is too
  light compared with the momentum transfer, when there are several t-channel mediators, or when an s/u-channel mediator
  is nearly degenerate with DM plus the SM fermion.

See the [direct detection module](https://maddmhep.github.io/maddm/dev/index.html#direct-detection-module) for details.

## Citation

If you are using this tool or the results obtained with it in a published work, please cite one or more of the [MadDM references](https://maddmhep.github.io/maddm/dev/index.html#references).
A [CITATION.bib](./CITATION.bib) file is also available.

If you use MadDM, please cite also the works of other authors that it uses. For direct detection:

- RAPIDD: D. G. Cerdeño, A. Cheek, E. Reid, H. Schulz, JCAP 08 (2018) 011, [arXiv:1802.03174](https://arxiv.org/abs/1802.03174).
- DM-nucleon scattering: the non-relativistic EFT of A. L. Fitzpatrick et al., JCAP 02 (2013) 004, [arXiv:1203.3542](https://arxiv.org/abs/1203.3542), and N. Anand, A. L. Fitzpatrick, W. C. Haxton, Phys. Rev. C 89, 065501 (2014), [arXiv:1308.6288](https://arxiv.org/abs/1308.6288); the nuclear responses of the chosen nuclear model and the chiral EFT structure factors (references in the [RAPIDD README](vendor/RAPIDD_for_DM/README.md#citation)).
- DM-electron scattering: R. Catena, T. Emken, N. A. Spaldin, W. Tarantino, Phys. Rev. Research 2, 033195 (2020), [arXiv:1912.08204](https://arxiv.org/abs/1912.08204), with the atomic responses of DarkART (T. Emken, [doi:10.5281/zenodo.6046225](https://doi.org/10.5281/zenodo.6046225)).
- Silicon: R. Catena, T. Emken, M. Matas, N. A. Spaldin, E. Urdshals, Phys. Rev. Research 3, 033149 (2021), [arXiv:2105.02233](https://arxiv.org/abs/2105.02233); R. Catena, N. A. Spaldin, Phys. Rev. Research 6, 033230 (2024), [arXiv:2402.06817](https://arxiv.org/abs/2402.06817); QCDark2, C. Dreyer, R. Essig, M. Fernandez-Serra, M. Hott, A. Singal, Phys. Rev. D 114, 035003 (2026), [arXiv:2603.12326](https://arxiv.org/abs/2603.12326); K. Ramanathan, N. Kurinsky, Phys. Rev. D 102, 063026 (2020), [arXiv:2004.10709](https://arxiv.org/abs/2004.10709).
- Halo conventions: D. Baxter et al., Eur. Phys. J. C 81, 907 (2021), [arXiv:2105.00599](https://arxiv.org/abs/2105.00599).
- Experimental limits and projections: LZ ([arXiv:2410.17036](https://arxiv.org/abs/2410.17036)), PICO-60 ([arXiv:1902.04031](https://arxiv.org/abs/1902.04031)), XLZD ([arXiv:2410.17137](https://arxiv.org/abs/2410.17137), [arXiv:2203.02309](https://arxiv.org/abs/2203.02309)), XENONnT ([arXiv:2601.11296](https://arxiv.org/abs/2601.11296)), DAMIC-M ([arXiv:2503.14617](https://arxiv.org/abs/2503.14617)), DarkSide-50 ([arXiv:2207.11968](https://arxiv.org/abs/2207.11968)), XENON10 (R. Essig, T. Volansky, T.-T. Yu, [arXiv:1703.00910](https://arxiv.org/abs/1703.00910)) and XENON1T ([arXiv:1907.11485](https://arxiv.org/abs/1907.11485)). The rescaled DAMIC-M limit also uses QCDark ([arXiv:2306.14944](https://arxiv.org/abs/2306.14944)).

The full list, with the references for each part of RAPIDD, is in the [RAPIDD README](vendor/RAPIDD_for_DM/README.md#citation). The references of the indirect detection inputs are in [CITATION.bib](./CITATION.bib), together with those listed above.

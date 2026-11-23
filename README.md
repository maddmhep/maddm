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
MadDM offers to install Pythia 8, the CosmiXs spectra and RAPIDD when a process first needs them; they can also be installed with `install pythia8`, `install cosmixs` and `install RAPIDD` at the MadDM prompt.
See the [installation guide](https://maddmhep.github.io/maddm/dev/index.html#install-maddm) for the full list of requirements and options.

## Citation

If you are using this tool or the results obtained with it in a published work, please cite one or more of the [MadDM references](https://maddmhep.github.io/maddm/dev/index.html#references).
A [CITATION.bib](./CITATION.bib) file is also available.

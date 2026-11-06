# CosmiXsPPPC Earth fluxes and shared Galactic halo

Extract this archive into the `PLUGIN` directory of MadGraph so that the
plugin directory is `PLUGIN/maddm`. Restart MadGraph. The ZIP includes all
54 CosmiXsPPPC propagated annihilation tables for antiprotons and
antideuterons and the six secondary spectra. There is no separate install
step for those data. Pythia and the optional PPPC4DMID/CosmiXs production
tables can still be used for **source** spectra and the combined dwarf
likelihood. RAPIDD remains the optional nuclear-recoil backend.

Relic density and indirect detection are selectable only if their matrix
elements are in the generated output. The launch menu now identifies
missing processes and prints the corresponding `add` command. You can add
relic density after generating direct detection; see the [combined run
example](README.md#running-relic-density-direct-and-indirect-detection-together).

```text
launch my_model_output
set indirect flux_earth
set fast
set rholocal 0.4
set vMP 238
set vescape 544
set cr_halo_profile NFW
set cr_propagation_scheme QUAINT
set cr_propagation_set MED
done
```

`set precise` is also supported. For existing output directories, the code
updates the inclusive rate routines before rebuilding; rerun `output` with
the updated plugin to also refresh the generated native recoil sources.
With `sigmav_method = reshuffling`, MadDM strips direct detection and other
MadDM generation shortcuts from legacy event banners before reconstructing
the annihilation matrix elements. Fresh indirect exports also exclude them.
The card reader silently ignores the obsolete `indirect_flux_earth_method`
entry in older cards. The propagation selectors are `cr_halo_profile`,
`cr_propagation_scheme` and `cr_propagation_set`.
Each BIG, QUAINT or SLIM secondary is shared across MIN, MED and MAX.

| File in `output/run_XX` | Contents |
|---|---|
| `CosmiXsPPPC_pbar_flux_earth.dat` | IS and top-of-atmosphere antiproton fluxes, per GeV kinetic energy |
| `CosmiXsPPPC_dbar_flux_earth.dat` | IS and top-of-atmosphere antideuteron fluxes, per GeV/n |
| `CosmiXsPPPC_flux_earth.json` | Table names, branching fractions, halo settings, modulation and status |

Both flux files contain separate secondary, DM and total columns for the
all-DM and thermal cases. They describe the model rate, with the pbar
model-point normalization, tilt and Fisk potential when that fit is
available. Without it, the reference secondary is used with `dbar_phi`
(default 0.7 GV). The dedicated fit/limit files and plots still report the
best fit and upper limit where available. The flux files can be generated
with `pbar_limits` and `dbar_events` both off. Integrated photon and neutrino
fluxes remain available internally; the terminal and `MadDM_results.txt`
omit their former summary. Source spectra and differential photon and
neutrino files remain available when the source inputs are configured, with
their existing target J factors.

The code contains a dormant positron charged-flux path with its own particle
mass and Fisk potential. It expects the propagated and secondary positron
tables in the same units and conventions as the charged inputs, under
`ExpData/CR_flux/positron_flux_<scheme>_<set>_<halo>.dat` and
`ExpData/AMS02_positron/secondary_<scheme>_IS.txt`. MadDM skips this path in
runs and prints “Positron results will come soon.”

`rholocal` is the single local density in GeV/cm³. If an older card supplies
only `rhoDM`, it is imported as `rholocal`; if both appear, `rholocal` wins.
MadDM writes only `rholocal` to cards and silently ignores the retired
`gce_rho_sun` and `pbar_rho_sun` entries in older cards.
Recoil rates are proportional to this density; Galactic annihilation fluxes
are proportional to its square. The published direct-detection upper-limit
curves retain the experiments' reference halos. Dwarf J factors and the
dwarf velocity input are independent of the Milky Way.

Select `halo_velocity_profile` for RAPIDD nuclear recoils and Galactic ID.
The card offers `SHM`, `SHM_numeric`, `SHM_beta`, `SHM_wLMC` and `Lisanti`.
See [the shared velocity profiles](docs/HALO_VELOCITIES.md) for the parameters
and the `Halo_velocity.json` diagnostics.

By default, the DD recoil and Galactic ID calculations use the same
truncated Galactic Maxwellian, set by `vMP` and `vescape` in km/s. Recoils
use the Earth's laboratory motion. The pbar, Dbar, GCE and Galactic line
rates average each partial `sigma*v` over the relative speeds of two
independent halo particles. The inclusive and numerical backends use the
same pair kinematics. The update also corrects the earlier inclusive
single-speed calculation, which used the pair speed in place of one
particle's speed. The optional `mw_velocity_model=fixed` uses the legacy
single-speed settings. `mw_shm_npts` controls the default SHM integration
(16 nodes); check convergence for sharp thresholds and narrow resonances.
The model uses the local SHM throughout the Galactic signals and does not
include a position-dependent Galactic velocity model.

| Code and data | Changes |
|---|---|
| `maddm_run_interface.py`, `maddm_interface.py` | Run routing, result summaries, optional source installer and recoil integration |
| `cr_earth.py`, `cr_config.py`, `pbar/core.py`, `pbar_interface.py`, `dbar_interface.py` | Bundled charged-flux construction, solar modulation and normalization |
| `velocity_profiles.py`, `halo_model.py`, `Templates/include/maddm_mw_rates.inc`, `Templates/src/maddm.f`, `Templates/src/indirect_detect_2to2.f` | Shared velocity laws and relative-speed quadrature, inclusive rates and generated-project migration |
| `Templates/src/directional_detection.f`, `python_templates/electron_recoil_signal.f`, `Templates/Cards/maddm_card*.dat` | Common recoil inputs and default card controls |
| `ExpData/` | Existing propagated tables and secondary spectra retained byte for byte |
| `docs/`, `tests/`, `validation/SHM_EARTH/` | Updated notes, regression tests and Linux MadGraph checks |

`SHA256SUMS` records the files shipped in this release. Historical validation
snapshots remain for comparison; their copies omit the obsolete propagation
card entry and the old integrated-flux display block.

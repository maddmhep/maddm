# AMS-02 antiproton limits in MadDM

**Revision 9:** common `rholocal`, automatic BIG/QUAINT/SLIM secondary inputs shared across MIN/MED/MAX, quiet rendering logs, explicit pbar rate/statistic names, revised GCE band/statistical-error plot and Argonne source default. See [README_GCE_CR.md](README_GCE_CR.md), [the detailed PDF note](docs/indirect_detection_notes.pdf), and `CHANGES_SINCE_V8.txt`. Earlier version descriptions below are historical.


**Version 8:** GCE prompt fitting, selectable charged-particle tables and updated run plots are documented in [README_GCE_CR.md](README_GCE_CR.md). This supersedes the fixed QUAINT-MED/NFW scope and plotting defaults described for version 7 below.

**Version 7 addition:** a third results block reports AMS-02/GAPS Dbar event
forecasts from the current model and propagated QUAINT-MED/NFW Dbar table.
Successful pbar flux fits also save `AMS02_pbar_comparison.pdf` automatically,
using the supplied updated script (linear residual panel, x axis to at least
1 TV). See [README_AMS02_DBAR.md](README_AMS02_DBAR.md) for the energy/flux
Jacobian, fitted/fixed solar potential, settings and saved outputs.

In `indirect = sigmav` mode this addition compares each generated two-body channel's rate with the precomputed BR=1 antiproton limit in `ExpData`, alongside the existing Fermi dwarf channel limits. In a flux mode it combines the **current two-body branching fractions** with the bundled **CosmiXsPPPC antiproton QUAINT-MED/NFW propagated fluxes** and profiles a model-specific upper limit. Results are printed and saved in the ordinary run outputs.

The implementation is based on the uploaded `maddm-rc-33-indirect.zip`. Its Git archive comment identifies commit `3b3f8d1141c10c3a80bafb63d47ef47e721c48b8` of the requested `rc/33/indirect` branch. The live GitHub branch could not be fetched, so this package uses the exact uploaded source as its baseline.

## Previous version 6: fast and precise sigmav launches

Install the complete directory as `PLUGIN/maddm`, restart MadDM and relaunch
the existing project. The plugin now synchronizes and verifies the generated
`DO_INDIRECT_DETECTION`, `DO_INDIRECT_SPECTRAL` and `ONLY2TO2LO` switches with
the selected method immediately before compilation. It rebuilds the main
driver, adds missing include dependencies and fixes generated make recipes
that previously ignored compilation failures. An old `maddm.x` cannot now
make a failed build appear successful. Successful builds record their method,
mode and switches in `maddm_build.json` in the generated project.

The supplied fast-run output has no annihilation-rate records. Its archive
alone does not establish which earlier compilation failed or which switches
were used for that run; the stale-binary acceptance and incomplete build
dependencies are confirmed defects repaired here.

For precise runs, MG5 now dispatches the indirect reweight operation to the
active plugin's implementation, even if another MadDM copy advertises the
same plugin key. Its F2PY build uses the Python interpreter running MadDM,
with that interpreter's tools searched first. The environment is restored
after compilation. This also avoids selecting a different environment's
`f2py` executable. Meson and Ninja remain required on Python >=3.12: install
them in the launching environment with the command below if missing. The
diagnostic now exercises the same F2PY wrapper used by real reweight builds.

Full MG5 3.8.0 launches of the supplied Scalar Higgs Portal project were run
on Linux with Python 3.13: precise sigmav followed by fast sigmav in the same
project. Both printed all 12 channel rates and Fermi/AMS-02 channel bounds.
The precise smoke test generated 100 events per velocity calculation; this
checks execution and output, not the statistical precision of a production
event sample. See `validation/AMS02_pbar/sigmav_build_validation_summary.json`
and the accompanying launch logs. The macOS toolchain was not executed.
The existing fast and reshuffling rate implementations are retained. In this
benchmark their totals differ by approximately a factor of two (fast about
1.18e-27 and precise about 2.36e-27 cm^3/s); the repaired workflow does not
establish numerical equivalence of those upstream methods. The saved outputs
make this distinction explicit rather than rescaling either method here.

## Previous update: Python 3.13 and zero-event shower failures

The two new failures have different causes. The attached precise `sigmav`
trace reaches the reweight matrix-element build, then fails with
`FileNotFoundError: No such file or directory: meson`. On Python >=3.12,
NumPy F2PY uses the Meson backend by default. Meson and Ninja must be installed
in the environment that launches MadGraph; this is an environment dependency,
not an antiproton likelihood error. In your active `chatgpt` environment run:

```bash
conda activate chatgpt
python -m pip install meson ninja
# From the MadGraph root, after installing the updated plugin:
python PLUGIN/maddm/check_precise_runtime.py --build-check
```

The last command checks the active interpreter and build tools, compiles a tiny
Fortran extension and imports it. It does not change model parameters or run
outputs. `requirements-precise.txt` lists the added build requirements. The
plugin checks these tools before exporting the precise reweight module and
reports the installation command if they are absent. It keeps the requested
reshuffling method; it does not substitute a different rate calculation.
See the [NumPy F2PY backend documentation](https://numpy.org/doc/2.1/release/1.26.0-notes.html#meson-backend-for-f2py).

The supplied `run_04` LHE file contains 10,000 complete events, but has 9,492
blank lines inside the required event/particle numeric blocks. Pythia 8.317's
LHE reader does not skip these lines; it fails before reading the first event
and reports EOF. This was reproduced with the actual version 8.317 reader.
The earlier patch fixed the Lund setting and filename selection but missed
this event-record formatting problem.

Before showering, MadDM now creates `maddm_shower_events.lhe` in the event
folder. It validates each event header and its NUP particle records, removes
only blank lines inside those numeric records, and preserves every nonblank
input record and its numbers. The supplied compressed/uncompressed file is
untouched. `shower_input.json` records the source filename, event count and
removed-line count. The generated card uses this prepared copy for every LHEF
assignment. The MadDM `Main:NumberOfEvents=-1` sentinel becomes Pythia's allowed
value 0, meaning read until EOF in the custom shower executable.

Both copies of `main101.cc` now check settings and initialization failures,
handle EOF, count successfully showered events rather than the zero-based
attempt index, and refuse to normalize an empty/failed shower. Histograms
use the LHE/Pythia event weights and divide by their sum, so a reweighted sample
keeps its relative event contributions. Constant-weight events give the
ordinary per-annihilation spectrum, including a valid one-event sample.
Previously saved spectra are cleared before the executable runs; invalid
spectra still produce an error rather than being silently replaced by zero.

For an existing generated project, the plugin recognizes the exact previous
shower source, saves `main101.cc.before_shower_fix`, installs the updated
source, and rebuilds the executable. `main101.build.json` records the source
hash and configured Pythia directory, preventing reuse of the old binary.
Unrecognized local source edits are preserved and reported for merging.
Install the complete updated plugin, restart MadDM, and relaunch; a new output
directory is not required. The first flux launch recompiles the shower, so a
working configured Pythia/C++ toolchain is required as usual.

Live validation replayed all 10,000 supplied events through the updated C++
executable with Pythia 8.317: all were accepted and all nine spectrum files
were finite. The resulting gamma spectrum gives a combined 55-dwarf endpoint
of 2.58322887e-25 cm^3/s with profiled J factors and the retained q=3.84145882
threshold. This is a shower replay of the supplied event sample, **not a fresh
MadGraph velocity-reweighted model run**. The actual uploaded reweight Fortran
module also compiled/imported with Python 3.13 and returned a finite positive
matrix element. See `validation/AMS02_pbar/shower_fix_validation_summary.json`
for the exact scope and logs.

## Precise-mode fixes and combined Fermi limit

The precise `sigmav` path keeps the MadDM velocity-reshuffling calculation.
The local `Indirect_Cmd.do_reweight` adapter translates MG5's generated
`--mode=indirect` option to the supported `--plugin=indirect` dispatch before
MG5 looks up the event file. This fixes the reported `No events file
corresponding to --mode=indirect run` exception without changing the rate
method. Beam energies and velocities are written even when no shower card is
present. The spectrum-free channel comparisons from the previous release are
retained for both fast and precise sigmav modes.

For a flux mode, the combined Fermi bound uses the actual gamma spectrum and
the bundled 55-dwarf likelihoods, with profiled J factors as before. It is not
constructed from the individual BR=1 upper limits. The uploaded run contained
only NaN photon yields: its Pythia log reports an unrecognized
`StringZ:deriveBLund` setting and an aborted initialization. The generated
shower card also contained a later subrun assignment to `events.lhe.gz`, even
though that run contained `unweighted_events.lhe.gz`.

The updated shower-card adapter inspects the installed Pythia settings XML.
For current Pythia versions it translates the legacy `deriveBLund=on` flag to
`StringZ:deriveLundPars=1`, which derives bLund from avgZLund; the requested tune
parameters are retained. Legacy installations keep their supported flag.
Every `Beams:LHEF` assignment, including subruns, is set to the actual event
file. The adapter needs a valid configured `pythia8_path` (or `PYTHIA8DATA`
pointing to the settings directory); it never guesses which tune settings an
unidentified installation supports. See the [official Pythia Lund settings](https://pythia.org/latest-manual/ZSelection.html).
A failed shower executable or non-finite spectrum now gives an explicit error
instead of silently becoming a zero signal and a -1 bound.

The combined-limit solver searches in log10(rate) within the requested range,
includes the physical zero-signal point, refines the best fit within bounds,
and solves the upper likelihood crossing above that best fit. The old
unconstrained simplex could leave the requested range. The model-point
likelihood now uses a log-likelihood reference, rather than a cross section
passed as that reference. `Fermi_dSph_limits.json` records the mass, dwarf
count, J-factor treatment, best fit, upper bound, threshold and sampled profile.
An unbracketed limit or invalid spectrum is reported explicitly.

The existing Fermi confidence convention is preserved: chi-square with one
degree of freedom at CDF 0.95, q=3.84145882. The separate antiproton likelihood
continues to use the one-sided threshold q=2.70554345. No BR=1 tables, flux
inputs, or antiproton likelihood parameters change in this update.

Replace the complete plugin, restart MadDM, and relaunch the run. The fixes in
this section are Python changes; no new matrix-element export is required.
Previously saved NaN shower spectra must be regenerated by relaunching the
flux calculation, not reused as valid predictions.

## Plot the saved antiproton fluxes

Use the standalone `examples/plot_ams02_pbar_fluxes.py` (also supplied as a
separate download). It needs NumPy and Matplotlib and imports no MadDM code:

```bash
python plot_ams02_pbar_fluxes.py /path/to/AMS02_pbar_fluxes.dat -o pbar_comparison.pdf
# Optional alternative fitted prediction:
python plot_ams02_pbar_fluxes.py /path/to/AMS02_pbar_fluxes.dat --prediction upper -o pbar_upper.pdf
```

The upper panel shows R^3 times the TOA flux per GV for AMS-02, secondary, DM,
and total, on base-10 logarithmic rigidity and flux axes. The flux-axis range
is chosen around the data and total so a vanishing DM tail does not stretch
the plot over many empty decades. The lower panel uses the exact data nodes
and `(total-data)/total`, with uncertainties `data_error/total`. It does not
refit, interpolate data, or mix secondary components from different fits.
Choose `--prediction model` (default), `thermal`, `best`, or `upper`; each
selects its corresponding secondary/DM/total block from the named header.

The updated default residual axis is linear, following the supplied script.
Use `--residual-scale symlog` for a signed logarithmic alternative, with a
linear region within +/-0.02; adjust `--linthresh` as desired. The rigidity
axis extends to at least 1000 GV. The output is a vector PDF with embedded TrueType fonts and
600-dpi settings for any raster elements. A checked example using the bundled
100 GeV bb fixed-rate benchmark is in
`validation/AMS02_pbar/bb_majorana_100/AMS02_pbar_comparison.pdf`.

## Try it in your MadGraph area

1. Put this complete plugin directory at `MG5_aMC/PLUGIN/maddm` in your test installation, retaining the original plugin separately if needed. Use the ordinary MadDM/MadGraph launcher for your installation.
2. Existing generated outputs such as `test_pbar_v0` can be relaunched after restarting MadDM. Writing the card include automatically repairs the generated declarations and inclusive driver, then the usual compilation rebuilds it. A fresh output can also be generated normally. Do not reuse an executable whose compilation failed.
3. Import your UFO model, define one dark-matter candidate, choose your couplings in the parameter card, and generate the desired indirect processes using your normal workflow. For example, replacing the uppercase placeholders:

```text
import model YOUR_UFO_MODEL
define darkmatter YOUR_DM_PARTICLE
generate indirect_detection
output YOUR_NEW_OUTPUT_DIRECTORY
launch
```

In the launch configuration enable `indirect = sigmav` or a flux mode. For a first tree-level check, use the following MadDM-card choices through the launch card editor:

```text
set sigmav_method inclusive
set pbar_limits True
set mw_velocity_model shm
```

Set the model masses and couplings in the usual parameter card. The pbar computation works in `sigmav` mode: it does not require a gamma-ray injection table or a Pythia shower. In a flux mode the existing gamma/neutrino workflow retains its own dependencies. `pbar_limits` defaults to true; setting it false disables the extra calculation.

The new likelihood requires Python 3.10+ with NumPy 1.24+ and SciPy 1.10+ as a supported environment for this addition. It has been tested here on Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0. No separately installed antiproton package is required. Existing MadGraph, compiler, model, shower and other MadDM requirements are not bundled into the ZIP.

For repeated scans on a system with threaded BLAS, setting `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` before launching MadGraph can reduce overhead from the small fits.

## Fast and precise sigmav mode: channel tables

Select `indirect = sigmav` and then use either `set fast` or `set precise`.
`set fast` selects inclusive Fortran rates. `set precise` selects reshuffling;
`sigmav_method madevent` is also supported. In all three cases the sigmav path
reads the per-channel limits without building a gamma spectrum, showering for
spectra, or fitting an antiproton flux for the branching mixture.
The numerical rate method still needs its ordinary MadGraph runtime.

For each continuum channel the Fermi section reports its rate at
`vave_indirect_cont`, the existing dwarf-table limit, and ALLOWED/EXCLUDED/NO
LIMIT. The antiproton section separately reports its Galactic SHM-averaged rate and the BR=1 AMS-02 bound. The distinction
between the two environments matters for velocity-dependent annihilation. Relic-fraction rescaling is
shown separately when `xsi < 1`. Pure-channel Fermi bounds use the inherited
self-conjugate table normalization and are multiplied by
`norm_Majorana_Dirac()/4` for the caller's particle convention. Pbar source and
density conversions follow the explicit card convention described below.

`ExpData/MadDM_AMS02_pbar_Limit_*.txt` contains all **19 channels and 64 native
masses per channel**, copied unchanged from the completed standalone scan:
1,132 finite limits and 84 absent-signal rows marked -1. The columns are mass
in GeV and the one-sided 95% rate limit in cm^3/s. Exact nodes are read directly;
between adjacent positive rows the pbar lookup interpolates in log mass and
log limit. It never extrapolates or crosses an absent-signal row. The existing
Fermi table interpolation is retained. Additional details are in
`ExpData/AMS02_pbar_channel_limits_README.md`.

These are **individual channel checks**, not a combined exclusion of the
branching mixture. A combined Fermi or pbar bound and likelihood are not
inferred by adding/averaging BR=1 limits. Their global scan fields stay at -1
in sigmav mode. Continuum and SM rate totals are nevertheless calculated from
the rates, so they no longer depend on requesting a spectrum. A flux mode
retains the existing spectrum-based Fermi and mixed-channel pbar likelihoods.

The tabulated pbar analysis fixes A in [0.8,1.2], tilt in [-0.1,0.1] about
10 GeV, phi profiled in [0,2] GV, and all 58 AMS-02 points. Changing the phi
range, pivot, or rigidity cuts cannot be applied to a precomputed table: such
a request is reported unavailable with a reason; use a flux mode to refit it.
Changing local density or the source convention is an exact rescaling and
is supported in either mode.

The result directory now contains `AMS02_pbar_channel_limits.dat` and
`AMS02_pbar_limits.json` in sigmav mode, and per-channel results are also written
to `MadDM_results.txt`. Scan columns named `AMS02_pbar_model#CHANNEL`,
`AMS02_pbar_limit#CHANNEL`, `AMS02_pbar_excluded#CHANNEL`, and
`AMS02_pbar_excluded_thermal#CHANNEL` record rates, limits and flags (0 allowed,
1 excluded, -1 unavailable). `AMS02_pbar_status=tabulated_channels` identifies
this path. `AMS02_pbar_fluxes.dat` belongs to the full likelihood path.

This update fixes fixed-width Fortran rate records, parsing of process names
with underscores, the inclusive card switch, and the numerical results
collector's fixed `[7:]` prefix removal. Uploaded subprocess names such as
`P1995_xsdmxsdm_bbx` now retain the full process name. Numerical rate errors
also receive the same velocity/unit conversion as the rates. Missing inclusive
records or a failed executable produce an explicit error, rather than invented
zero rates or an apparently allowed model. Restart MadDM after installing the
updated plugin; inactive copies should remain outside `PLUGIN`.

## Fix for ConfigFile.get compatibility and plugin warnings

The error `ConfigFile.__getitem__() takes 2 positional arguments but 3 were given` came from calling `maddm_card.get(key, default)`. In the MG5 runtime bundled with the supplied run, `ConfigFile.get` aliases `__getitem__`, which accepts only the key. The antiproton run hook and all adapter settings now use indexed access with a KeyError fallback through `pbar_interface.card_value`. Explicit zero and false values are preserved. No global MadGraph API is monkey-patched.

For this Python-only fix, replace `maddm_run_interface.py` and `pbar_interface.py` in the active `PLUGIN/maddm` directory with the updated versions, restart MadDM and relaunch the existing output. Regenerating or recompiling matrix elements is unnecessary for this fix. The earlier generated-include repair remains included.

The uploaded run includes the actual MadGraph `internal/banner.py`. The reported TypeError was reproduced with that implementation. An additional test executes the actual MadDMCard class on this runtime, reads the uploaded card, writes its Fortran include assignments and runs the actual pbar hook with a controlled partial-rate fixture and real flux fits. It does not rerun MadGraph event generation.

The reported SyntaxWarnings and validation warning name `PLUGIN/maddm-feat-id_spectra`, an additional plugin copy. MadGraph scans directories under PLUGIN when discovering extensions, so an old copy can be imported even while the main run interface comes from PLUGIN/maddm. Move inactive MadDM copies outside PLUGIN rather than renaming them inside that directory, then restart MadGraph. The updated branch's three files named in those SyntaxWarnings compile without invalid-escape warnings.

The warning about `bypass_check=[]` comes from MadGraph's list validation. It was reproduced on the supplied runtime: assigning an empty list emits a warning and leaves the empty default unchanged. It is non-fatal and does not cause the antiproton TypeError. This patch does not alter MadGraph's validation logic or bypass any checks.

## Fix for undeclared antiproton symbols in an existing output

If gfortran reports that `do_pbar_limits` or `vave_indirect_pbar` has no implicit type in `maddm_card.inc`, the generated card assigns new variables while the generated `include/maddm.inc` lacks their declarations. This mismatch can occur when Python/card files and generated Fortran sources come from different versions. The distributed template contains the declarations; the compilation path now also verifies the generated files before writing card assignments.

After installing this updated plugin, restart MadDM and relaunch. The include writer adds missing declarations and independent common blocks without duplicating declarations already present. If the generated main driver also predates the antiproton rate calculation, it adds that block before the output file is closed. Original changed files are retained with suffix `.before_ams02_pbar`.

To repair the existing `test_pbar_v0` immediately, run from the MadGraph root:

```bash
python PLUGIN/maddm/fix_pbar_include.py test_pbar_v0
cd test_pbar_v0
make -j1 all
```

This repair changes only `include/maddm.inc` and, if needed, `src/maddm.f`. It does not alter model parameters or regenerate matrix elements. It addresses the two reported undeclared symbols; a complete local compiler/linker run must still succeed. The helper refuses incompatible existing types or an unrecognizable driver rather than making an unsafe edit.

## Milky Way rate calculation

Antiprotons, antideuterons, the GCE and Galactic gamma lines share the
`vMP` and `vescape` settings of the recoil SHM. MadDM averages each partial
cross section times relative speed over two independent velocities in this
halo. Fast inclusive runs do this in Fortran; precise runs integrate the
generated partonic processes at the quadrature energies. The rates used
for the dwarf likelihood keep their separate `vave_indirect_cont` setting.
The Galactic partonic integrations do not generate a second LHE shower sample.

`mw_velocity_model=fixed` selects the earlier single-speed Galactic
prescription with `vave_indirect_pbar` and `vave_indirect_line`. A missing
Galactic rate is reported instead of silently renormalizing the branching
mixture. See [the shared halo and Earth flux notes](CHANGES_COSMIXSPPPC_SHM.md).

## Channel matching and flux normalization

The mixture construction in this section describes the flux-mode likelihood.
The sigmav path instead reports the independent channel-table checks above.

The adapter uses PDG final states, independent of particle names or their ordering. It supports:

```text
ee mumu tautau nu_e_nu_e nu_mu_nu_mu nu_tau_nu_tau
uu dd ss cc bb tt AA gg WW ZZ hh AZ hZ
```

`AA` means two photons. Partial rates with the same supported final state are summed. Branching fractions are derived from present-day partial annihilation rates, not freeze-out rates or relic-density channel percentages.

At the chosen mass and branching fractions, the DM flux is the sum of the separately interpolated channel fluxes. The upper-limit fit varies their common rate normalization while keeping the current mass and branching fractions fixed. Changing model couplings afterwards can change the branching fractions and requires a new evaluation.

All positive-rate generated channels must be supported. Unsupported BSM channels, gamma-Higgs, non-two-body channels, explicit decay-chain selections, coannihilation or mismatched initial states cause **NOT EVALUATED**, with details in JSON. Such rates are not discarded and the remaining SM subset is not normalized to one. Zero-rate unsupported channels do not affect the spectrum. Only one DM species is accepted; Dirac rates require the particle-antiparticle initial state and a symmetric abundance. The total rate refers to the processes actually generated, so generate all channels whose contributions you want included.

The native table reference is 3e-26 cm^3/s and local NFW density 0.4 GeV/cm^3. Flux units are per GeV per square metre per second per steradian. The code applies no additional source integral, angular factor, or mass-squared factor to these already propagated tables.

The particle/antiparticle source convention is an **explicit input assumption**, controlled by `pbar_table_convention`:

| Setting | Native table assumption | Flux multiplier for MadDM self-conjugate DM | Flux multiplier for symmetric Dirac DM |
| --- | --- | ---: | ---: |
| `self_conjugate`, default | Source proportional to one half of density squared times rate divided by mass squared | 1 | 1/2 |
| `particle_antiparticle` | Source proportional to one quarter of total density squared times particle-antiparticle rate divided by mass squared | 2 | 1 |

MadDM's `norm_Majorana_Dirac` supplies the model convention. The selected native-table convention is stored in every successful JSON result. The supplied table README fixes the rate and density reference but does not fully document the self-conjugacy convention; this could not be independently verified from the live paper here. **Check that the selected card convention matches the intended native table normalization before interpreting the bound.** Choosing the alternative convention is supported and tested; there is no hidden extra factor of two.

`pbar_rho_sun` rescales the DM flux by the square of its ratio to 0.4 GeV/cm^3, keeping the NFW shape fixed. It is independent of the gamma-line profile settings and direct-detection density. If MadDM derives `xsi` below one from the relic density, the interface evaluates both the full-local-density model and a thermal scenario whose flux is multiplied by `xsi` squared. It does not apply that rescaling twice. The reported cross-section upper limit is the full-local-density bound; the two model decisions are separately profiled.

## Likelihood and nuisance parameters

The supplied secondary flux has zero solar potential and is interpreted as an interstellar antiproton flux. Its deformation is

$$
J_{\rm sec}^{\rm IS}[T]
=A\left[\frac{T}{10\,{\rm GeV}}\right]^\delta
J_{\rm sec,0}^{\rm IS}[T],
\qquad 0.8\le A\le1.2,\quad -0.1\le\delta\le0.1.
$$

These are hard bounds without Gaussian penalties. The normalization range applies at the energy pivot; the tilt can produce a larger fractional shift far from that pivot. Positive delta hardens the spectrum. The energy power is a spectral-index deformation, not a detector energy-calibration shift.

A single force-field potential is fitted for both secondary and DM fluxes, with default **0–2 GV**, and is recomputed at every trial rate:

$$
J^{\rm TOA}[T]=
\frac{T\left[T+2m_p\right]}
{\left[T+\phi\right]\left[T+\phi+2m_p\right]}
J^{\rm IS}[T+\phi].
$$

The shift is phi, not phi divided by two. The secondary tilt is evaluated before solar modulation, at the shifted IS energy. The AMS measurements are fluxes per rigidity, so the code converts the model using

$$
T[R]=\sqrt{R^2+m_p^2}-m_p,
\qquad J_R[R]=J_T[T[R]]\frac{R}{\sqrt{R^2+m_p^2}}.
$$

The supplied 58 points at 1.08–427.5 GV and their symmetric exported errors enter a diagonal Gaussian chi-square. Normalization is eliminated analytically within its bounds; tilt and potential are fitted with multiple numerical starts and an independent polish. At fixed mass and mixture, the physical nonnegative signal best fit sets the reference. The upper endpoint uses a one-sided asymptotic 95% threshold of 2.705543454095404 in the profiled chi-square difference.

For each model prediction, `q_upper` is zero below the best-fit rate and is the profiled chi-square difference above it. **EXCLUDED** means `q_upper` exceeds the threshold. This decision uses the model's actual rate, rather than simply testing the secondary-only reference. In the usual connected case it agrees with comparing the model rate against the reported endpoint. If multiple intervals are found the endpoint is the largest detected upper crossing; the actual-point likelihood-ratio decision remains separate. Finite scans cannot guarantee detection of arbitrarily narrow additional minima or intervals.

Both signal and nuisance parameters are reprofiled exactly within this model, with no linearized-background approximation or Asimov calibration factor. Up to 16 exact mass/mixture/configuration results are cached to avoid repeating limits when only the total model rate changes. Each actual model point is still reprofiled.

## Card settings

| Parameter | Default | Purpose |
| --- | ---: | --- |
| `pbar_limits` | `True` | Enable the additional computation when indirect detection is active |
| `vave_indirect_pbar` | `0.00075` | Milky Way `vave` in units of c, used in fixed mode only |
| `pbar_phi_min`, `pbar_phi_max` | `0`, `2` | Force-field potential fit interval in GV |
| `pbar_pivot` | `10` | IS secondary tilt pivot in GeV |
| `pbar_rho_sun` | `0.4` | Local NFW density in GeV/cm^3 |
| `pbar_rigidity_min`, `pbar_rigidity_max` | `0`, `-1` | Data selection in GV; negative maximum uses all points |
| `pbar_table_convention` | `self_conjugate` | Native propagated-table source normalization assumption |

The secondary normalization and tilt bounds remain 0.8–1.2 and -0.1–0.1. Their definitions and default values are shared by all model channels.

## Output

The terminal prints the MW model rate, 95% upper limit, and **ALLOWED / EXCLUDED** decisions for full-local-density and, where relevant, thermal DM. Unavailable/disabled cases print **NOT EVALUATED**.

The following fit outputs describe the full flux-mode likelihood. Sigmav-mode
channel-table outputs are listed above.

Every flux-mode evaluated run has:

- `MadDM_results.txt`: the normal results file gains an `AMS02_pbar` section, partial-channel branching fractions, upper endpoint, model rate, likelihood-ratio values, and decisions.
- `AMS02_pbar_limits.json`: complete input normalization, branching fractions, original process rates and velocity, nuisance fits, sampled likelihood profile, assumptions and status. Failure and unsupported-channel statuses also produce a JSON file.
- `AMS02_pbar_fluxes.dat`: data/errors and fitted secondary, signal and total **TOA flux per GV**, for the full-density model, thermal model, global best fit and upper endpoint. Every prediction uses its own profiled nuisance values. Unavailable results do not leave an old flux file behind.
- Sequential scan summaries: stable numeric columns such as `AMS02_pbar_sigmav`, `AMS02_pbar_model_sigmav`, `AMS02_pbar_q`, `AMS02_pbar_excluded`, and their thermal counterparts. Exclusion flags are 1 for excluded, 0 for allowed, and -1 when unavailable. Missing limits are -1. Per-run JSON retains string status and diagnostics.

The pbar likelihood is **not automatically added to MultiNest's joint likelihood**; it is computed per launched model and the requested sequential-scan outputs are supported. Existing joint indirect likelihood behavior is preserved.

## Validation and an offline check

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest discover -s tests -v
# Optional real-card/runtime test, using an existing generated project:
MADDM_PBAR_TEST_RUN=/path/to/test_pbar_v0 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python validation/AMS02_pbar/run_benchmarks.py
python examples/check_ams02_pbar.py --mass 200 --rates '{"bb":2.1e-26,"WW":0.9e-26}' --xsi 0.5 --output pbar_check
```

The example accepts supplied MW rates; it does not compute a UFO model or matrix elements. A normal MadDM launch obtains those rates from its generated code.

**Version 6: 73 tests passed under both Python 3.12 and 3.13, with the uploaded-runtime and live compiled-shower tests enabled.** These include likelihood/physics checks, all channel tables, generated-card/source compatibility, result writers, stale-build failures, active-plugin reweight dispatch, interpreter/environment restoration and exact-velocity rate reuse. Full MG5 3.8.0 launches of the supplied project now also complete on Linux: precise sigmav, a subsequent fast sigmav launch, and precise flux_source with a live 100-event Pythia shower and both combined bounds. The 100-event run is an execution smoke test, not a production-precision spectrum. Existing fixed-rate and synthetic-spectrum tests remain fixtures; they are identified separately in the validation records. The macOS toolchain remains untested. See `validation/AMS02_pbar/sigmav_build_validation_summary.json` for the current scope, hashes and saved outputs.

For a 100 GeV bb benchmark with supplied rate 3e-26 cm^3/s and `xsi=0.5`:

| DM convention, native table assumed self-conjugate | Upper endpoint in cm^3/s | Full local density | Thermal density |
| --- | ---: | --- | --- |
| Self-conjugate | 1.77608776e-26 | Excluded | Allowed |
| Symmetric Dirac | 3.55217551e-26 | Allowed | Allowed |

For a 200 GeV 70% bb + 30% WW mixture, the self-conjugate endpoint is 1.34552231e-26 cm^3/s. These are fixed-rate adapter checks, not results from a generated physical model. The factor-of-two conversion and thermal decisions are checked by the benchmark script.

## Input limitations retained from the standalone package

- The secondary file is the exact supplied QUAINT zero-modulation prediction. Its header does not provide a transport parameter card, so its exact correspondence to the DM table's QUAINT-MED setup remains unverified. No transport parameters are varied or reconstructed by this code.
- The supplied AMS export has no covariance matrix or original rigidity-bin edges. It is used with diagonal errors and point evaluation, not an exact forward-folded experimental likelihood. A single effective solar potential approximates the time-averaged dataset.
- The source convention of the native propagated tables must match the explicit setting discussed above.
- The confidence threshold is asymptotic; coverage has not been calibrated with mock experiments.
- Mass interpolation is linear in the upstream tabulated flux, with PCHIP energy interpolation in logarithmic coordinates. No smoothing or extrapolation is performed. Near channel thresholds the interpolation remains an approximation. The upstream 1e-100 sentinel is treated as empty; zero tabulated signal returns no bound rather than an invented limit.

All three numerical input files are unchanged from the supplied original project. Their provenance and hashes are in `ExpData/AMS02_pbar/provenance.json`. The full modified/added-file list is in `CHANGED_FILES_AMS02.md`.

## References

- M. Cirelli, M. Di Mauro and A. Kar, *Refined anti-proton and anti-deuteron fluxes from weak-scale Dark Matter*, arXiv:2606.10017. [CosmiXsPPPC propagated tables](https://github.com/CosmiXsPPPC/CosmicRay-flux), [paper](https://arxiv.org/abs/2606.10017).
- AMS Collaboration, *The Alpha Magnetic Spectrometer on the International Space Station: Part II — Results from the First Seven Years*, Physics Reports 894, 1–116, 2021, DOI 10.1016/j.physrep.2020.09.003. [Official publication and Table 6](https://ams02.space/publications/202102). The supplied SSDC export identifies this dataset; the data file is used unchanged.
- Supplied USINE secondary prediction, originally named `local_fluxes_1HBAR_Ekn_Model1DKisoVc_SolMod0DFF_phi0_00GV_QUAINT.out`.

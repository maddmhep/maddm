# Validation scope

For the current Dbar event and automatic plot additions, see
`../AMS02_dbar/README.md` and its validation summary (84 passing tests).

## Previous version 6: complete MG5 3.8.0 launches

`sigmav_build_validation_summary.json` records the current validation.
`sigmav_launches/` contains full fast and precise sigmav launch logs, saved
channel results, generated flags and build evidence for the supplied 1000 GeV
Scalar Higgs Portal project. Precise mode uses the actual MG5 generator and
velocity reweighter at both velocities. Fast mode was launched after precise
mode in the same project, rebuilding the driver with inclusive switches.
Both modes produce 12 channel comparisons with Fermi and AMS-02 tables.
The inherited rate methods are unchanged; in this benchmark the fast total is
about half the precise total. This validation checks successful calculation
and reporting, not numerical equivalence of those upstream methods.
A full equal-velocity precise launch also confirms reuse of all 12 rates
without creating another pbar event run (`precise_same_velocity/`).

A further precise flux_source smoke test used 100 events, the actual compiled
Pythia 8.317 executable, finite spectra and both combined likelihoods. This
small event count verifies the workflow; its numerical Fermi limit is not a
production spectrum prediction. The active-Python build probe also passes.

`sigmav_build_final_tests.log` and `sigmav_build_python313_tests.log` record
73 passing tests under Python 3.12.14 and 3.13.16 respectively, including the
supplied 10,000-event live-shower regression from the previous release.
The new tests also reproduce masked make failure and stale executable reuse,
verify current generated switches, pin active plugin dispatch, restore build
environments, and distinguish equal/different-velocity rate collection.

The validation copy restored source symlinks/executable permissions from the
uploaded ZIP, relocated the model's absolute macOS path to its included UFO
copy, and rebuilt platform-specific objects for Linux using gfortran 13.3.0.
MG5 tag v3.8.0 is used without changes to its physics implementation. These
adaptations affect the private validation copy, not the user's original run
or the bundled model-independent plugin. The macOS toolchain was not run.

## Previous version 5: supplied-event shower and Python 3.13 build checks

- `shower_fix_final_tests.log`: 68 passing tests, including the actual Python
  shower/reader/Fermi dispatch with the compiled C++ executable and supplied
  10,000-event sample. Optional environment variables enable the live test.
- `shower_fix_python313_tests.log`: the same 68 tests also pass under Python 3.13.16.
- `shower_fix_validation_summary.json`: exact current scope and source/input hashes.
- `shower_regression_summary.json` and `shower_replay/`: live executable tests of
  the full sample, one-event EOF, empty input (nonzero exit and no spectra), and
  event-weight normalization. All nine spectra from the full sample are finite.
- `shower_replay_Fermi_dSph_limits.json`: the combined likelihood profile for
  that replayed photon yield at 1000 GeV. This is not a new velocity-reweighted
  model run.
- `python313_build_probe.log`: a genuine F2PY/Meson Fortran compile/import probe.
- `python313_reweight_build.log`, `python313_reweight_import.log` and
  `python313_reweight_module.json`: the uploaded 12-channel reweight Fortran
  module rebuilt/imported under Python 3.13, with a finite positive matrix element.
  Source symlinks from the uploaded ZIP were restored in an isolated build copy;
  Darwin object files were rebuilt for Linux. The uploaded project was untouched.

The live shower cases can be reproduced with a compiled updated main101:

```bash
python validation/AMS02_pbar/run_shower_regressions.py --input-run /path/to/test_pbar_v0 --executable /path/to/test_pbar_v0/bin/internal/main101 --pythia-data /path/to/pythia8/share/Pythia8/xmldoc --output /path/to/shower_check
```

The regular suite also accepts `MADDM_PYTHIA_EXECUTABLE` and
`MADDM_PYTHIA_DATA` alongside `MADDM_PBAR_TEST_RUN` to enable the live dispatch
regression. Without the supplied run and compiled executable, those optional
tests skip; they do not fabricate events. At version 5, no complete MadGraph
import/generation or velocity-reweight launch had been performed. The compiler/shower validation uses
Linux and does not execute the user's macOS toolchain.

## Previous release validation records


The included tests passed on 2026-10-22 with Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0.

- `precise_fermi_final_tests.log`: 60 passing tests including the actual uploaded MG5 event lookup, generated shower card, combined real 55-dwarf likelihood with a specified analytic photon yield, model-point probability at the upper endpoint, and numerical residuals for all plot prediction blocks.
- `precise_fermi_validation_summary.json`: current regression scope, source hashes, and Pythia diagnostics from the supplied run.
- `fermi_combined_fixture.json`: full likelihood profile for a **synthetic** photon yield at 500 GeV (not a UFO result). The profiled-J endpoint is 1.2818974035e-25 cm^3/s at q=3.84145908.
- `bb_majorana_100/AMS02_pbar_comparison.pdf`: checked vector-PDF plot of the existing fixed-rate pbar fixture.
- `sigmav_tests.log`: 47 passing tests with both optional uploaded-runtime tests enabled; the actual card switches and source repair are checked. `tests.log` records the preceding 35-test release.
- `benchmark_summary.json` and `benchmarks.log`: actual profile-fit results for three supplied partial-rate fixtures, including self-conjugate/Dirac normalization and a two-channel mixture.
- The three benchmark directories contain the full sampled likelihood profiles and 58-point flux comparisons, generated from these fixed rates.
- `run_benchmarks.py` reproduces the fixtures and verifies the factor-of-two source conversion and thermal-density decisions.
- `validation_summary.json` records the exact baseline-archive SHA-256 and its Git archive comment, test scope and environment.

These are adapter/likelihood calculations using explicitly supplied rates. They are not results from generated UFO models. The numerical generator calls in the Python hook tests are controlled harness calls; those earlier records did not compile Fortran or run Pythia. The current additional live/build validations are described above; those earlier records did not include full MadGraph generation or madevent/reshuffling launches.

The new sigmav tests exercise all 1,216 native table rows, missing neighboring
rows and mass extrapolation, density/source conversion, actual fast and precise
shortcuts, the actual numerical results collector with different subprocess-ID
lengths, inclusive/reshuffling/madevent hooks without any spectrum access, actual
console and summary writers, thermal decisions, unsupported rows, and rate
aggregation. The full-likelihood regressions remain part of the suite.

The tests compile and execute the actual modified Python method bodies from the main run interface, including the real result and scan writers. Card registration and both template expansions are checked without substituting a different production implementation. The standalone likelihood is checked against an independent full optimizer and an analytic Gaussian limit, and tests verify the rigidity Jacobian, force-field shift, IS secondary tilt, spectral interpolation, branching linearity and treatment of empty tables.

Physical limitations such as missing experimental covariance/bin edges, the secondary transport parameter card, the native-table source convention, and asymptotic confidence coverage are described in the main `README_AMS02_PBAR.md`; software tests do not remove them.

The include compatibility fix is covered by six tests: legacy-header/driver repair with idempotence and backups, preservation of current sources, uppercase continued declarations, rejecting incompatible types, refusing an invalid driver before mutation, and invocation through the actual card include writer. Those earlier include-repair tests are source/ordering checks. The current Fortran compilation evidence concerns the supplied reweight matrix module, not a full MadDM driver rebuild; version 6 adds a full driver rebuild and launches.

The ConfigFile compatibility fix replaces dict-style `.get(key, default)` reads with indexed card access and a KeyError fallback. The tests now use a card object with the actual MG5 one-argument `.get` signature throughout the hook harness. Three additional contract tests and one test of the actual MadDMCard class running on the uploaded project's internal runtime cover this failure. The actual configuration include writer and pbar hook both complete. The partial rates in this test are controlled fixtures, so this is not a full physical model/event-generation validation. Set `MADDM_PBAR_TEST_RUN` to an existing generated run to enable the optional test; without it, that test is skipped.

`config_api_summary.json` records the hash of the exact uploaded runtime source used to reproduce the exception and the non-fatal empty `bypass_check` warning.

## Previous precise-mode regressions

The MG5 reweight regression first reproduces the reported exception with the
uploaded event-lookup methods, then executes the actual local command adapter
and verifies plugin dispatch and the resolved LHE file. It does not generate
or reweight events. The actual uploaded MadGraph PY8Card implementation
produces the shower-card fixture, including its duplicate subrun LHE setting;
normalization to modern and legacy settings is tested with small XML capability
fixtures. Those earlier tests did not use a live Pythia installation; the latest tests above do.

The combined Fermi test executes the actual flux_source dispatch, spectrum
reader, and Fermi likelihood class on all 55 supplied dwarf likelihood tables.
Its explicitly synthetic input is dN/dlog10(x)=10 exp(-10x) on 70 logarithmic
x nodes from 1e-4 to 1, at mDM=500 GeV. The fixture tests numerical and dispatch
correctness; it does not predict the user's Scalar Higgs Portal result. Both
the upper crossing and model-point probability agree with the retained Fermi
threshold q=3.84145882. The pbar one-sided threshold is unchanged.

The standalone plot was rendered and visually inspected. It uses base-10 log
rigidity and flux axes, and a signed symmetric-log residual axis (linear within
+/-0.02). Tests check all four fitted prediction blocks and the exact
(model-data)/model denominator at the 58 measurement nodes.

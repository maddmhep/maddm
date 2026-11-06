# Files changed for AMS-02 antiprotons

**Revision 9:** common `rholocal`, automatic BIG/QUAINT/SLIM secondary inputs shared across MIN/MED/MAX, quiet rendering logs, explicit pbar rate/statistic names, revised GCE band/statistical-error plot and Argonne source default. See [README_GCE_CR.md](README_GCE_CR.md), [the detailed PDF note](docs/indirect_detection_notes.pdf), and `CHANGES_SINCE_V8.txt`. Earlier version descriptions below are historical.


**Version 8:** GCE prompt fitting, selectable charged-particle tables and updated run plots are documented in [README_GCE_CR.md](README_GCE_CR.md). This supersedes the fixed QUAINT-MED/NFW scope and plotting defaults described for version 7 below.

Baseline archive commit `3b3f8d1141c10c3a80bafb63d47ef47e721c48b8`.

Latest update, package version 7, 2026-10-22: antideuteron event forecasts
and automatic antiproton comparison PDFs. Changes relative to version 6:

- `dbar_interface.py` (new): propagated-table mixture, total-energy to
  per-nucleon energy/flux conversion, force-field modulation, integration
  over the two AMS-02 windows and the GAPS window, two-event sensitivity
  normalization, all-DM/thermal forecasts and diagnostic/summary writers.
- `ExpData/AMS02_dbar/dbar_flux_QUAINT_MED_NFW.dat` and `provenance.json`
  (new): unmodified annihilation table from the original supplied archive,
  with exact units, source member and hash.
- `maddm_run_interface.py`: reuse the collected MW partial rates for Dbar;
  add calculation, third console/results block and scan columns; register
  new settings; enable shared MW rate output when only Dbar is requested.
- `pbar_interface.py`: call the supplied plot after a successful flux write;
  remove stale PDFs and allow plotting to be disabled.
- `pbar_plot.py` (new), `examples/plot_ams02_pbar_fluxes.py`: the attached
  updated script, unchanged; linear residual panel, rigidity through 1 TV,
  vector PDF output with 600-dpi settings.
- `Templates/Cards/maddm_card.dat`, `Templates/Cards/maddm_card_full.dat`:
  `dbar_events`, `dbar_phi` and `pbar_plot` settings.
- `requirements-ams02.txt` (new): NumPy, SciPy and Matplotlib requirements.
- `CITATION.bib`, `README.md`, `README_AMS02_PBAR.md`, `README_AMS02_DBAR.md`
  (new), tests and `validation/AMS02_dbar/`: sensitivity reference, usage,
  physical conventions, regression tests and complete MG5 launch evidence.

The pbar/dwarf likelihoods, all 19 pbar bounds, existing flux/data inputs,
Fortran physics sources, C++ shower and precise rate calculations are
unchanged from version 6. The Dbar feature adds signal forecasts, not bounds.

Previous update, package version 6: reliable sigmav builds,
active-plugin precise dispatch and velocity-matched rate reuse.
Production changes relative to version 5:

- `maddm_run_interface.py`: synchronize and verify generated switches before
  compilation; require a fresh driver executable; record build switches;
  explicitly select the active indirect reweighter; reuse complete existing
  rates only when every requested velocity matches the pbar velocity.
- `pbar_interface.py`: repair ignored-error make recipes and include
  dependencies in existing projects; invalidate stale driver objects,
  executables and build metadata; check generated switches.
- `Templates/makefile_template`, `Templates/src/makefile`: make build failures
  propagate and compile all driver objects against current generated includes.
- `MGoutput.py`, `indirect_compat.py`, `f2py_compat.py` (new): bind the real
  reweight compilation to the launching Python and its build tools; prevent
  inherited make jobserver flags from confusing F2PY's Ninja child process.
- `check_precise_runtime.py`: exercise that same interpreter-bound wrapper.
- `tests/test_sigmav_build.py` (new), `tests/test_pbar_integration.py`,
  `tests/test_precise_and_fermi.py`: build-failure/stale-executable regression,
  generated-switch verification, active-plugin dispatch and environment
  restoration, same/different-velocity collection.
- Documentation and `validation/AMS02_pbar/`: complete MG5 3.8.0 launch logs
  and outputs for the supplied model, runtime build probe and regression logs.

All 19 pbar channel tables, secondary/DM/data inputs, likelihood settings,
Fortran physics sources and C++ shower sources are unchanged from version 5.
The complete cumulative list relative to the original branch follows.

Previous update, package version 5: LHE numeric-record formatting, empty-shower
normalization and Python 3.13 reweight build prerequisites. Production changes
relative to package version 4:

- `indirect_compat.py`: make a validated Pythia input copy, remove numeric-block
  blanks, retain all nonblank records; map the all-event sentinel; check
  Meson/Ninja; upgrade the known generated shower source and record build identity.
- `maddm_run_interface.py`: use the prepared LHE and diagnostics, rebuild old
  shower executables, clear old spectra before each shower.
- `MGoutput.py`: check precise F2PY build prerequisites before standalone export.
- `Templates/bin/internal/main101.cc` and `Templates/bin/main101.cc`: check card
  reading/initialization, stop correctly at EOF, count accepted events, reject
  empty/failed normalization, fill/normalize histograms with event weights.
- `check_precise_runtime.py` and `requirements-precise.txt` (new): active-environment
  diagnostic and optional real Fortran compile/import probe; build dependencies.
- `tests/test_shower_input.py` (new), `tests/test_precise_and_fermi.py`: source-upgrade,
  formatting/build-tool regressions and an optional live C++/Python/Fermi test.
- `validation/AMS02_pbar/run_shower_regressions.py` (new): full supplied sample,
  one-event, empty-input and weighted-sample tests of the compiled executable.
- Documentation and `validation/AMS02_pbar/`: current logs, finite replayed
  spectra, combined Fermi endpoint and Python 3.13 matrix-module build evidence.

The earlier fast/precise sigmav channel-table implementation remains included.
All 19 BR=1 pbar tables and the three numerical input flux/data files are
byte-identical to package version 4. This update changes C++ shower sources but
no Fortran model or MadDM driver sources.
The complete cumulative file list relative to the original branch follows.

Install the complete updated plugin at `PLUGIN/maddm`, keeping inactive copies
outside `PLUGIN`, then restart and relaunch. Existing generated Fortran sources
are repaired automatically on include writing and must compile successfully.
The earlier ConfigFile.get compatibility fix remains included.

## Modified existing files

| File | Change |
| --- | --- |
| `MGoutput.py` | Precise F2PY dependency check and interpreter-bound compilation. |
| `Templates/makefile_template`, `Templates/src/makefile` | Propagate build failures and include dependencies. |
| `Templates/bin/internal/main101.cc` and `Templates/bin/main101.cc` | Checked shower initialization/EOF; successful-event and event-weight normalization. |
| `CITATION.bib` | AMS and propagated-flux citations. |
| `README.md` | Feature documentation link. |
| `Templates/Cards/maddm_card.dat` | Antiproton settings. |
| `Templates/Cards/maddm_card_full.dat` | Antiproton settings. |
| `Templates/src/maddm.f` | Separate inclusive Milky Way rates and full-width process records. |
| `maddm_run_interface.py` | Antiproton workflow, indexed ConfigFile access, generated-source repair, sigmav channel tables, corrected rate totals/collector and output, precise reweight dispatch, shower compatibility/validation, combined Fermi profile. |
| `python_templates/maddm.inc` | Antiproton declarations/common blocks. |

## Added files

- `AMS02_PBAR_MANIFEST.sha256`
- `CHANGED_FILES_AMS02.md`
- `ExpData/AMS02_pbar/pbar_ams02.txt`
- `ExpData/AMS02_pbar/pbar_flux_QUAINT_MED_NFW.dat`
- `ExpData/AMS02_pbar/provenance.json`
- `ExpData/AMS02_pbar/secondary_QUAINT_IS.txt`
- `ExpData/MadDM_AMS02_pbar_Limit_*.txt` (19 files)
- `ExpData/AMS02_pbar_channel_limits.json`
- `ExpData/AMS02_pbar_channel_limits_README.md`
- `README_AMS02_PBAR.md`
- `examples/check_ams02_pbar.py`
- `examples/plot_ams02_pbar_fluxes.py`
- `indirect_compat.py`
- `fermi_profile.py`
- `check_precise_runtime.py`
- `requirements-precise.txt`
- `f2py_compat.py`
- `tests/test_sigmav_build.py`
- `validation/AMS02_pbar/sigmav_build_validation_summary.json`
- `validation/AMS02_pbar/sigmav_launches/`
- `validation/AMS02_pbar/sigmav_build_final_tests.log`
- `validation/AMS02_pbar/sigmav_build_python313_tests.log`
- `tests/test_shower_input.py`
- `validation/AMS02_pbar/run_shower_regressions.py`
- `validation/AMS02_pbar/shower_fix_validation_summary.json`
- `validation/AMS02_pbar/shower_fix_final_tests.log`
- `validation/AMS02_pbar/shower_replay/` (live replay logs and finite spectra)
- `fix_pbar_include.py`
- `pbar/__init__.py`
- `pbar/core.py`
- `pbar_interface.py`
- `tests/test_pbar_include_repair.py`
- `tests/test_pbar_integration.py`
- `tests/test_pbar_likelihood.py`
- `tests/test_sigmav_channel_limits.py`
- `tests/test_precise_and_fermi.py`
- `tests/test_pbar_plot.py`
- `validation/AMS02_pbar/README.md`
- `validation/AMS02_pbar/precise_fermi_final_tests.log`
- `validation/AMS02_pbar/precise_fermi_validation_summary.json`
- `validation/AMS02_pbar/fermi_combined_fixture.json`
- `validation/AMS02_pbar/bb_majorana_100/AMS02_pbar_comparison.pdf`
- `validation/AMS02_pbar/bb_WW_majorana_200/AMS02_pbar_fluxes.dat`
- `validation/AMS02_pbar/bb_WW_majorana_200/AMS02_pbar_limits.json`
- `validation/AMS02_pbar/bb_dirac_100/AMS02_pbar_fluxes.dat`
- `validation/AMS02_pbar/bb_dirac_100/AMS02_pbar_limits.json`
- `validation/AMS02_pbar/bb_majorana_100/AMS02_pbar_fluxes.dat`
- `validation/AMS02_pbar/bb_majorana_100/AMS02_pbar_limits.json`
- `validation/AMS02_pbar/benchmark_summary.json`
- `validation/AMS02_pbar/benchmarks.log`
- `validation/AMS02_pbar/config_api_summary.json`
- `validation/AMS02_pbar/run_benchmarks.py`
- `validation/AMS02_pbar/tests.log`
- `validation/AMS02_pbar/validation_summary.json`

The three numerical flux/data inputs and the 19 copied limit tables are
unchanged. The version 6 suite has 73 tests, including live C++ shower checks.
Complete fast/precise sigmav launches of the supplied project now run under
MG5 3.8.0 and Python 3.13 on Linux. This includes generating and compiling the
velocity-reweight matrix module. The macOS toolchain remains untested. Rate
fixtures in other workflow tests remain controlled inputs; the current
validation summary distinguishes those from full launches.

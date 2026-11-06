# Version 7 validation

The Dbar feature is tested both independently and in full MG5 3.8.0 launches
of the supplied 1000 GeV Scalar Higgs Portal project on Linux/Python 3.13.16.
`validation_summary.json` records exact results and hashes.

- `fast_pbar_disabled/`: inclusive sigmav launch with `pbar_limits=False`
  and `dbar_events=True`. The driver still emits MW rates, Dbar events are
  calculated, the pbar likelihood is disabled, and phi=0.7 GV is explicitly
  identified as fixed.
- `precise_sigmav/`: real generation and velocity reshuffling; 12 Fermi and
  pbar channel comparisons, plus Dbar forecasts from the same saved MW rates.
  No additional Dbar event run is launched. Phi=0.7 GV is fixed, and no pbar
  PDF is produced without a fitted flux file.
- `precise_flux_source/`: real generator/reweight/shower run with 100 events
  per velocity calculation. Both combined bounds succeed, Dbar reuses exactly
  the pbar process records and all-DM/thermal model-point potentials, and
  `AMS02_pbar_comparison.pdf` is saved automatically in the output folder.
  The PDF was rendered and visually checked. The small sample verifies
  execution, not production-precision photon spectra or limits.
- `fast_sigmav/`: full fast launch after the precise run, with both pbar and
  Dbar enabled. All three result blocks are saved.

`dbar_final_tests.log` and `dbar_python313_tests.log` record **84 passing
tests with no skips** under Python 3.12.14 and 3.13.16 respectively. They
include the previous live 10,000-event shower regression and 11 new tests:

- Sensitivity flux gives exactly two events, including the window width.
- Independent analytic checks of the K to T/n energy/flux Jacobian and the
  force-field momentum ratio/charge-to-mass shift.
- Adaptive integration agrees with a dense independent trapezoid integral.
- Channel mixtures, rate/density scaling, Dirac convention and thermal scaling.
- Model-point fitted phi selection versus upper-limit/best-fit phi; fixed fallback.
- Unsupported positive rates, zero-rate models and mass-range refusal.
- Saved units, row counts, totals and stale-output removal.
- Actual run hooks reuse rates and work with pbar disabled.
- The supplied plot settings, residual denominator and actual PDF save hook.

The actual MadDMCard regression also checks that Dbar-only operation enables
the shared Fortran rate-output switch while retaining `pbar_limits=False`
in the Python card. No physics or shower Fortran/C++ source changes were
needed for this release.

The validation copy restores ZIP source symlinks/executable permissions,
relocates the absolute macOS model path to its included UFO, and rebuilds
platform-specific objects for Linux using gfortran 13.3.0. Stale model-import
pickle caches in that relocated private UFO copy were cleared before fresh
precise launches. These are validation-environment adaptations; no UFO or
MG5 source changes are included in the plugin. macOS execution is untested.

To reproduce the suite from the plugin directory:

```bash
MADDM_PBAR_TEST_RUN=/path/to/supplied/test_pbar_v0 \
MADDM_PYTHIA_EXECUTABLE=/path/to/updated/main101 \
MADDM_PYTHIA_DATA=/path/to/pythia8/share/Pythia8/xmldoc \
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python -m unittest discover -s tests -v
```

Without the supplied run/executable, optional live checks skip. The remaining
tests use explicitly controlled model rates, not simulated UFO calculations.
Existing pbar/dwarf inputs, likelihoods, all 19 pbar limit tables and rate
implementations are byte-identical to version 6. The previously documented
fast/precise numerical rate difference is not altered by these forecasts.

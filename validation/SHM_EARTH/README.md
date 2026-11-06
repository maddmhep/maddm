# Shared halo and charged Earth-flux validation

The snapshots are from MadGraph5_aMC@NLO 3.8.0 on Linux, using its built-in
`MSSM_SLHA2` UFO and `generate indirect_detection b b~` for neutralino
annihilation. Both `set fast` and `set precise` completed in `flux_earth`
mode with no installed photon/neutrino source inputs. Both wrote pbar and
Dbar model-point spectra, the pbar likelihood products and the GCE fit. The
combined dwarf bound was unavailable because there was no prompt spectrum;
the per-channel dwarf bound and the theoretical rate were reported.

| Result | Inclusive, 16 nodes | MadEvent, 8 nodes |
|---|---:|---:|
| Galactic b bbar rate [cm³/s] | 2.3320383e-28 | 2.331573e-28 |
| AMS-02 pbar 95% upper limit [cm³/s] | 1.093133e-26 | 1.093133e-26 |

The inclusive and numerical rates agree within 0.020%. The Fortran
fixed-speed calculation and the numerical dwarf survey also agree at this
s-wave point after the incoming-particle momentum correction. Neither
summary prints the integrated prompt fluxes. Both summaries show the
positron placeholder.

The `rho_04` and `rho_08` directories record two otherwise identical runs
with both nuclear recoils and charged fluxes enabled. All 3,000 nonzero
Xe, Ar and Ge RAPIDD recoil-rate bins scaled by 2. The pbar and Dbar model
IS and TOA fluxes scaled by 4. Both runs produced the charged spectra with
the pbar likelihood, Dbar forecast and GCE fit disabled, using the reference
secondaries and fixed 0.7 GV modulation. Nuclear cross sections and dwarf
rates did not change. The RAPIDD table has 200 speed rows for this backend.

The offline regression suite ran 106 tests (100 passed, six optional
external-fixture tests skipped). Its new checks cover the exact second
moment of the relative-speed distribution, recoil halo normalization,
charge modulation, pbar/Dbar flux consistency and the dormant positron
output. Native macOS and a live Pythia shower were not exercised here.

# Version 8 validation

Both Python 3.12 and Python 3.13 suites pass all 95 tests without skips; logs are included. These include the uploaded runtime and live shower regressions. Tests cover all 54 table selections, GCE units and asymmetric amplitude fitting, conjugacy factors, missing secondary handling, Dbar counts and pbar plot selection.

`full_fast` records a complete default QUAINT-MED/NFW fast sigmav run. `full_precise_BIG_MIN_Einasto` records a complete precise flux-source run with BIG-MIN and Einasto, using MG5aMC 3.8.0 and the supplied 1000 GeV scalar-Higgs-portal UFO on Linux. Fermi, pbar, Dbar and GCE statuses are all ok in the latter run. Its PDFs were regenerated from the saved products and visually checked. This is a software integration example, not a claim that this model fits the GCE. Absolute paths in snapshots record the validation environment.

The GCE uses 11 detections and reports the 15 upper-limit bins separately. Its fit is prompt-only with diagonal asymmetric systematic errors. The charged secondary MIN/MED/MAX grid is not supplied: QUAINT/BIG reference backgrounds are held fixed within each scheme, and SLIM pbar requires an explicit input. See README_GCE_CR.md for all assumptions. macOS execution has not been tested here.

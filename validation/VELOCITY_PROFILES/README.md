# Validation of the RAPIDD velocity profiles

The joint tests used MadGraph5_aMC@NLO 3.8.0 on Linux, its built-in
MSSM_SLHA2 neutralino, and the RAPIDD C nuclear recoil backend. All
five profiles completed in inclusive `flux_earth` mode with the
same card inputs: vMP = 238 km/s, vescape = 544 km/s, and rholocal =
0.4 GeV/cm³. The default SHM recoil table and spectra match the
previous package byte for byte. SHM_numeric gives the same spectra.

| Profile | RMS relative speed [km/s] | Galactic neutralino-to-neutrino rate [cm³/s] | Xe recoil rate relative to SHM at 20 keV |
|---|---:|---:|---:|
| `SHM` | 401.981 | 3.20296756e-33 | 1.00000000 |
| `SHM_numeric` | 401.981 | 3.20296756e-33 | 1.00000000 |
| `SHM_beta` | 389.353 | 3.00488731e-33 | 1.00626212 |
| `Lisanti` | 376.535 | 2.81030537e-33 | 1.01327778 |
| `SHM_wLMC` | 405.775 | 3.26372325e-33 | 0.99705015 |

The neutrino row is a velocity-suppressed channel in the inclusive
Fortran output. Across the five profiles its rate follows the second
moment of the pair distribution to better than 3e-8 in the ratio.
The recoil column uses the 500th row of RAPIDD's recoil spectrum.

The precise `SHM_wLMC` run used eight Galactic nodes and completed
with RAPIDD nuclear recoils, both charged Earth fluxes, the
antiproton likelihood, the antideuteron forecast and the GCE fit.
The Galactic b bbar rate from precise MadEvent was 2.33200030e-28 cm³/s
(0.00163% difference from the inclusive result). Optional prompt source
spectra were not installed. A separate `flux_earth` run without
RAPIDD used the native default SHM and produced both charged fluxes.

The suite ran 114 tests: 108 passed, and six optional external fixtures
were skipped. The RAPIDD comparison exercised its C integrals for the
hard SHM, numerical SHM, Lisanti and SHM_wLMC at fixed lab velocity.
The beta law was checked independently with the Galactic cutoff.

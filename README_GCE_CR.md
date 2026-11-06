# GCE fit, Dbar sensitivity plots and propagated-table selection

Install the complete updated plugin at `PLUGIN/maddm`, restart MadDM and
relaunch the existing generated project. All new products are written under
`output/run_XX`. NumPy, SciPy and Matplotlib are listed in
`requirements-ams02.txt`. No network download is needed at runtime.

## Shared antiproton and antideuteron settings

```text
set cr_halo_profile NFW
set cr_propagation_scheme QUAINT
set cr_propagation_set MED
set pbar_secondary_file auto
```

The supported DM table choices are:

- `cr_halo_profile`: `NFW`, `Einasto`, `Burkert`.
- `cr_propagation_scheme`: `QUAINT`, `BIG`, `SLIM`.
- `cr_propagation_set`: `MIN`, `MED`, `MAX`.

All 27 combinations for each species are bundled in `ExpData/CR_flux`,
unchanged from the original supplied annihilation archives. The settings
select the actual DM flux table used by the pbar likelihood, pure-channel
bounds, Dbar counts and both charged-particle figures. They appear in result
metadata and plot annotations. The old `profile` setting for gamma-line
searches is separate and retains its existing meaning.

**Secondary inputs:** all three supplied USINE zero-potential spectra are now bundled for each species. `auto` selects the matching BIG, QUAINT or SLIM spectrum. The same scheme-specific secondary is used for MIN, MED and MAX; the DM primary table changes with all three controls. Pbar normalization, tilt and solar potential are still profiled. An explicit `pbar_secondary_file` can override the automatic pbar input. Dbar reference secondary curves use the corresponding scheme and do not enter the DM signal event counts.

**One Galactic density and velocity distribution:** `set rholocal 0.4` sets the density in GeV/cm^3 for recoil rates, pbar, Dbar, the GCE and Galactic gamma lines. Charged DM fluxes scale as `(rholocal/0.4)^2` and recoil rates are linear in the density. The GCE uses `rholocal^2` with the dimensionless Jbar. The shared SHM reads `vMP` and `vescape` from the card: recoil rates use the lab boost and Galactic annihilation uses the distribution of pair relative speeds. Dwarf velocities and J factors are separate. See [shared halo settings](CHANGES_COSMIXSPPPC_SHM.md).


In `sigmav` mode the original QUAINT-MED/NFW defaults retain the bundled BR=1
pbar bounds in `ExpData`. For any other table choice, explicit secondary
override or nuisance/range setting, pure-channel limits are profiled at the
requested mass instead of reusing those default bounds. This takes longer
than a table lookup. The selected numerical annihilation-rate method remains
fast inclusive or precise reshuffling as requested. No combined branching
mixture limit is inferred in sigmav mode. Flux modes continue to profile the
full model mixture. Table selection is part of the cache key.

## Galactic-center excess amplitude fit

New settings:

```text
set gce_fit True
set gce_plot True
set gce_jbar 120.0
set rholocal 0.4
set gce_rsun_kpc 8.122
```

The first standard-analysis block in the attached `GCE_flux_data` is used:
40 by 40 degree ROI, 0.5--1000 GeV and solid angle 0.4775494185 sr.
Its SED values and systematic bounds are converted from ROI-integrated
`MeV cm^-2 s^-1` to `GeV cm^-2 s^-1 sr^-1` by multiplying by 1e-3 and dividing
by that solid angle. The predicted SED is evaluated at the supplied reference
energy of each bin; no different bin averaging is invented.

Prompt gamma yields for all 19 two-body channels and 64 native masses are
bundled, selected without changing values from the official CosmiXs source
table. They are combined with the current model's **Milky Way partial rates
and branching fractions**, not the dwarf-velocity branching fractions. Mass
interpolation is linear at fixed x=E/mDM, without mass extrapolation.
These gamma production yields are distinct from the propagated charged-CR
tables. `ExpData/GCE/provenance.json` records both source and data hashes.

The dimensionless, ROI-averaged Jbar defaults to 120 in the supplied package.
Its convention follows Table I and Eqs. 37--39 of arXiv:2508.02775. The paper's
density reference is 0.345 GeV/cm^3 and Solar
radius 8.122 kpc. Equations 37--39 give the prompt intensity; the implementation
uses denominator 16 pi for particle--antiparticle DM and 8 pi for self-conjugate
DM, following the model's existing MadDM convention. It does not multiply an
averaged intensity by the ROI solid angle a second time. The GCE geometry
remains independent of `cr_halo_profile`: changing a charged-CR table does
not silently replace the explicitly requested Jbar. Set the GCE geometry inputs deliberately if using a different Galactic-center shape; use `rholocal` for its common density normalization.

At fixed mass and branching shape, only a nonnegative total annihilation
amplitude is fitted. Measured bins use the asymmetric errors `data-lower`
and `upper-data` from the **systematic envelope**, with the side chosen by
whether the prediction is below or above the measurement. The statistical
error column is not substituted or added in quadrature. The envelope is
treated as a diagonal 1-sigma error model; no covariance was supplied.

There are 11 measured bins and 15 upper-limit rows with `dnde=-1` in this
block. A `-1` sentinel is never fitted as a negative or invented zero
measurement. The reported chi-square therefore uses the 11 measurements,
with 10 degrees of freedom for the single fitted amplitude. The 15 upper
limits and their systematic upper bands are retained as **separate checks**,
with per-bin predictions and violation flags. No likelihood for those
censored bins can be reconstructed from an upper bound alone. A small
detected-bin chi-square is not a certificate that the upper limits are met.

The error interval is the conditional Delta-chi-square=1 interval around
the minimum, saved as separate lower/upper errors. Geometry uncertainties
are not profiled. No chi-square p-value is asserted for a systematic envelope.

The calculation is **prompt only**: inverse-Compton and bremsstrahlung
components from later equations in the reference paper are not added.
Equation 38 describes multi-body cascade final states. Those can use an
actual MadDM runtime gamma spectrum in a flux mode with
`mw_velocity_model=fixed` when its continuum velocity equals the Milky Way velocity and there is no omitted
separate gamma-line component. Otherwise unsupported nonzero channels are
reported rather than being mapped to a fictitious two-body yield. For ordinary
supported two-body models the bundled CosmiXs mixture is used in all modes.

GCE fitting never changes the model's couplings, its calculated rate or the
pbar/Dbar predictions. It measures the normalization required by that fixed
spectral shape. New numeric scan columns include `GCE_sigmav_best`,
`GCE_sigmav_error_minus`, `GCE_sigmav_error_plus`, `GCE_chi2`, `GCE_ndof` and
`GCE_chi2_ndof`.

## Figures and run products

**GCE:** `GCE_fit.json` and `GCE_fit.txt` give the best-fit rate, asymmetric
errors, chi-square, degrees of freedom, input geometry, branching fractions
and upper-limit checks. `GCE_spectrum.dat` saves all 26 data rows, systematic
bounds, upper limits, model/best-fit predictions and the fitted-bin mask.
`GCE_prompt_spectrum.dat` saves the dense model, best-fit and error-envelope
curves. `GCE_fit.pdf` compares them with the SED on log10 axes.

**Dbar:** `AMS02_dbar_comparison.pdf` compares the current model TOA flux
with the supplied GAPS and two AMS-02 sensitivity windows, annotating the
model rate, Fisk potential and expected DM signal counts in each window.
It uses the same table, branching mixture, solar modulation and integration
as `AMS02_dbar_events.json/dat`. A thermal-density curve/count summary is
included when xsi<1. A supplied secondary reference is shown separately;
it is not included in the DM event forecast. Its plotted values are saved
in `AMS02_dbar_secondary.dat`. `dbar_plot=False` disables this figure.

**Pbar:** `AMS02_pbar_comparison.pdf` now shows secondary, DM and total
fluxes at the **upper-limit nuisance configuration**, plus a separate DM
curve at the **model's coupling-derived rate**, with the latter's model-point
fitted solar potential. Labels distinguish those two configurations.
The residual panel uses `(upper-total - data)/upper-total` and remains linear,
as in the supplied plotting script. The flux axes use log10 scales and the
rigidity axis extends to at least 1 TV. All components are already saved
as the `upper_*` and `model_*` blocks in `AMS02_pbar_fluxes.dat`.

All three figures are vector PDFs with embedded fonts and 600-dpi settings.
Plots can be regenerated from completed run products without new MC events:

```bash
python PLUGIN/maddm/examples/plot_run_products.py /path/to/output/run_XX
```

Default pbar likelihood thresholds and the original 19 pbar bound files are
unchanged. The documented upstream fast/precise rate difference is retained.
See `validation/GCE_CR/` for exact regression scope, example outputs and full
MG5 launch evidence. The macOS toolchain has not been run here.

## Revision 9 outputs and notes

`MadDM_results.txt` prints `AMS02_pbar_sigmav_UL95`, `AMS02_pbar_sigmav_best_fit` and `AMS02_pbar_sigmav_model`, with definitions. The former numerical names remain scan aliases. `q_upper_model` and its thermal version are one-sided profiled Delta chi-square statistics at the theoretical rate and at `xsi^2` times that rate. They are not p-values and coincide when `xsi=1`.

The GCE figure uses a grey systematic envelope, black detected points with statistical errors and black arrows for published upper limits. The fit still uses systematic errors; statistical errors are appended as column 11 of `GCE_spectrum.dat`.

The Pythia default is `Main:methodDbar=4.0000000000e+00` (Argonne). Update an existing explicit card value if changing it is intended. Pythia generates source yields; the new charged signals at Earth use the fixed CosmiXsPPPC propagated tables, which are not regenerated by changing a source card.

Plotting suppresses only Matplotlib/fontTools DEBUG/INFO chatter that conflicts with MG5's coloring formatter. Real warnings/errors and other MadDM logging remain active.

The 22-page implementation note is `docs/indirect_detection_notes.pdf`, with editable LaTeX, separate `references.bib` and example figures. `CHANGES_SINCE_V8.txt` lists this revision's changes. Current validation is under `validation/V9`.

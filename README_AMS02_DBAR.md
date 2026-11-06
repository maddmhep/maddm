# Antideuteron event forecasts and automatic antiproton plots

The forecasts use the model's Galactic SHM-averaged two-body rates and the
bundled CosmiXsPPPC propagated antideuteron tables. The default is
QUAINT-MED/NFW; all three profiles and each BIG, QUAINT or SLIM MIN/MED/MAX
set are included. Antideuteron forecasts work in rate-only or flux modes
and reuse the Galactic rates used for antiprotons and the GCE. See
[the shared halo and Earth flux documentation](CHANGES_COSMIXSPPPC_SHM.md).

## Definition of the event estimate

For each constant-sensitivity window, use exactly

$$
N_{\bar D}=\frac{2}{\Phi_{\mathrm{sens}}\,[T_2-T_1]}
\int_{T_1}^{T_2}\Phi_{\bar D}^{\mathrm{TOA}}[T]\,dT,
\qquad T=K/2.
$$

| Experiment/window | T1 to T2 in GeV/n | Flux sensitivity |
| --- | --- | --- |
| AMS-02 low | 0.18 to 0.70 | 8e-7 |
| AMS-02 high | 2.0 to 4.2 | 2e-7 |
| GAPS | 0.1 to 0.3 | 2e-6 |

Flux sensitivity units are `(GeV/n)^-1 m^-2 s^-1 sr^-1`. Each supplied
sensitivity corresponds to **two events**. A flux equal to the sensitivity
throughout a window therefore yields exactly two expected events. The two
AMS-02 expectations are added because the windows are disjoint. GAPS is
reported separately. These are fractional expected **DM signal counts**, not
integer observations, confidence limits or a detailed detector simulation.
No secondary Dbar background is included. This implements the sensitivity
prescription supplied with the request and used in
[arXiv:2602.15132](https://arxiv.org/abs/2602.15132), Figs. 5--7.

## Energy variable, flux units and solar modulation

The native table tabulates total kinetic energy K in GeV and dPhi/dK, not
energy per nucleon. For Dbar, A=2 and |Z|=1. For each requested TOA T, the code
evaluates the IS table at `K_IS = 2*T + phi` and applies

$$
\frac{d\Phi^{\mathrm{TOA}}}{dT}
=2\,\frac{K_{\mathrm{TOA}}[K_{\mathrm{TOA}}+2m_{\bar D}]}
{K_{\mathrm{IS}}[K_{\mathrm{IS}}+2m_{\bar D}]}
\frac{d\Phi^{\mathrm{IS}}}{dK}\bigg|_{K_{\mathrm{IS}}},
\qquad K_{\mathrm{TOA}}=2T,
\qquad m_{\bar D}=1.875613\,\mathrm{GeV}.
$$

The leading 2 is the differential-flux Jacobian. The equivalent per-nucleon
solar energy shift is `phi/2`, not `phi`. An adaptive integration over each
window includes the transformed native interpolation knots; no endpoint or
mass extrapolation is allowed. The saved plotting grid is not the integration
grid.

When a full pbar fit succeeds, the Dbar calculation uses the fitted **model
point** Fisk potential. It uses the all-DM and thermal model-point fits
separately. It does not use the potential at the upper-limit endpoint or at
the global best-fit rate. No secondary normalization or tilt is applied to
the Dbar DM signal.

In sigmav mode, or when no pbar model-point fit exists, it uses the explicit
`dbar_phi` setting, default **0.7 GV**. This is a fixed input and is identified
as such in the console and JSON; it is not labelled a fitted potential.

## Rates, normalization and table coverage

The DM mass, branching fractions and total rate come from the model point using the
same Galactic SHM-averaged channel rates as the pbar and GCE calculations. The table's reference rate is 3e-26
cm^3/s and reference local density is 0.4 GeV/cm^3. Fluxes scale with the
current rate and density squared. The explicit source convention and local
density settings are shared with pbar: `pbar_table_convention` and
`rholocal`. Thermal predictions use the model rate multiplied by `xsi^2`.
With distinct fitted potentials they also include that change of modulation;
with a common fixed potential the counts scale exactly as `xsi^2`.

The unmodified annihilation tables cover the same 64 native masses and 19
channels as pbar, for all bundled profiles and propagation choices. Decays,
multi-body/BSM final states and mass extrapolation are not
added. An unsupported nonzero channel makes the forecast unavailable rather
than being discarded or renormalized away. Exact zero rates yield zero counts.
Channels populated only by the native 1e-100 sentinel give no tabulated
signal and are identified explicitly; they do not establish an exact
physical zero production rate. `ExpData/AMS02_dbar/provenance.json` records
the source member, units, normalization and hash.

## Settings and outputs

New card settings, enabled by default:

```text
set dbar_events True
set dbar_phi 0.7
set pbar_plot True
```

`dbar_phi` is the fallback described above. Set `dbar_events False` to disable
the forecast; set `pbar_plot False` to disable automatic PDF output. Dbar
forecasts also work with `pbar_limits False`: the driver still collects MW
rates for Dbar without enabling the pbar likelihood. The generated legacy
`DO_PBAR_LIMITS` Fortran switch controls this shared rate output; the user's
Python `pbar_limits` setting remains unchanged.

Every indirect run writes these files in `output/run_XX`:

- `AMS02_dbar_events.json`: status, model rates, mixture, source convention,
  density, fitted/fixed potentials, exact windows, integrated/average fluxes,
  event expectations and numerical integration errors.
- `AMS02_dbar_events.dat`: successful forecasts, one row per case/window.
- `AMS02_dbar_fluxes.dat`: successful IS and TOA Dbar DM fluxes per GeV/n,
  for the all-DM and thermal cases.
- `MadDM_results.txt`: a third block with both AMS windows, their total,
  GAPS and the modulation potentials. Stable numeric event columns are also
  included in parameter-scan summaries.

Unavailable or disabled forecasts have explicit status and -1 numeric
sentinels. Stale forecast/flux files are removed. No Dbar exclusion flag is
assigned from these projected sensitivities.

## Automatic antiproton figure

After a successful pbar flux fit, `AMS02_pbar_comparison.pdf` is saved in the
same run directory. It uses the attached updated plotting script: the upper
panel has R^3 times the measured, secondary, DM and total TOA flux, with
base-10 log axes; the lower panel has a **linear** signed residual axis and
`(model-data)/model` at the measured rigidities. The rigidity axis extends to
at least 1000 GV. No theoretical extrapolation beyond the saved data grid is
invented. The model-point secondary/DM/total block is selected consistently.
The PDF is vector output with embedded fonts and 600-dpi settings.

No automatic plot is made in sigmav mode, which has no fitted flux file.
Stale PDFs are removed before a new attempt. Missing Matplotlib or a plotting
error is reported without cancelling the physics results. Install
`requirements-ams02.txt` in the environment running MadDM to enable plotting.

The standalone script remains available as
`examples/plot_ams02_pbar_fluxes.py`, identical to the attached version; the
integrated implementation is `pbar_plot.py`.

```bash
python -m pip install -r PLUGIN/maddm/requirements-ams02.txt
python PLUGIN/maddm/examples/plot_ams02_pbar_fluxes.py \
    output/run_XX/AMS02_pbar_fluxes.dat -o output/run_XX/AMS02_pbar_comparison.pdf
```

Install the complete updated plugin and restart MadDM. Existing generated
projects can be relaunched. No new model export is needed for this addition.
The previously documented fast/precise rate-method difference is unchanged.

## Validation

84 tests pass with no skips under Python 3.12 and 3.13, including the prior
live 10,000-event shower regression and the new integration/Jacobian/scaling,
rate-reuse, card-switch and plotting tests. Full MG5 3.8.0 runs cover fast
sigmav, precise sigmav and precise flux_source, including Dbar-only operation
with pbar disabled and automatic PDF creation in a flux run. See
`validation/AMS02_dbar/README.md` and `validation_summary.json` for exact
scope. Linux is tested; the macOS toolchain is not executed here.

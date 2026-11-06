# Shared velocity distributions for RAPIDD and Milky Way annihilation

Use `halo_velocity_profile` in `maddm_card.dat` to select the local Galactic
velocity law for both RAPIDD nuclear recoils and Milky Way indirect detection.
For example, in the launch menu:

```text
set halo_velocity_profile SHM_wLMC
set vMP 238
set vescape 544
set halo_lmc_fraction 0.006
```

`vMP` is the Maxwellian scale `v0` used by RAPIDD, in km/s. It is the mode of
the untruncated Maxwellian speed distribution and is not its mean.
`vescape` is the Galactic cutoff for the smooth component, also in km/s.
Both values come from the same MadDM card for DD and ID. The default is the
native hard-cutoff `SHM`, with `vMP = 238` and `vescape = 544`. The indirect
calculation needs no RAPIDD installation; RAPIDD is the backend for the
nuclear recoil calculation. `rholocal` remains the common total local density.

| `halo_velocity_profile` | Galactic velocity law | Additional card inputs |
|---|---|---|
| `SHM` | Native Maxwellian with a hard cutoff | None |
| `SHM_numeric` | Same hard-cutoff Maxwellian, integrated numerically | None |
| `SHM_beta` | Subtracts a constant from the Maxwellian inside the escape sphere, then renormalizes | `halo_beta` (default 1; 0 recovers the hard SHM) |
| `Lisanti` | Smoothly suppressed high-speed tail with a shape exponent | `halo_lisanti_k` (default 1.5; supported from 0.1) |
| `SHM_wLMC` | Mixture of the hard-cutoff SHM and a second, boosted truncated Maxwellian | `halo_lmc_fraction`, `halo_lmc_vb`, `halo_lmc_sigb`, `halo_lmc_vcut`, `halo_lmc_cosb`, `halo_lmc_azimuth_deg` |

The alternatives use the underlying velocity laws defined by
[RAPIDD](https://github.com/cheekyparticle/RAPIDD_for_DM), with inputs
corresponding to `calc_halo` in `rapidd/core.py`. Within the escape sphere,
the unnormalized phase-space factors for `SHM_beta` and `Lisanti` are
respectively

$$
\exp[-v^2/v_0^2] - \beta\exp[-v_{\rm esc}^2/v_0^2],
\qquad
\left[\exp[-v^2/[k v_0^2]]-\exp[-v_{\rm esc}^2/[k v_0^2]]\right]^k.
$$

They vanish outside the sphere. In particular, MadDM enforces this support
for `SHM_beta` before integrating; the checked standalone RAPIDD recoil
integrator omitted the angular escape cut for that option and can return
negative values in the tail. MadDM supplies RAPIDD with a table from the
normalized, cutoff distribution used for annihilation as well. The
standalone `calc_halo` default for `beta` is 0; MadDM sets `halo_beta=1` so
choosing `SHM_beta` selects the smooth cutoff by default.

The LMC defaults are a local density fraction of 0.006, bulk speed 570 km/s,
Maxwellian scale 100 km/s, and a cutoff of 200 km/s around the LMC bulk.
`halo_lmc_cosb` defaults to −0.71 and gives the cosine with the RAPIDD
reference lab direction. That cosine defines a cone; `halo_lmc_azimuth_deg`
chooses a fixed Galactic direction on it. Zero degrees points toward the
projection of the Galactic radial axis onto the plane transverse to the
reference direction. For DD, MadDM boosts this fixed LMC vector into the
lab over the year. The annihilation speed distribution is independent of
the stream orientation because both components are isotropic in their own
rest frames and the observable depends only on the pair speed.

## How the two calculations use the same law

For DD, MadDM computes the annual mean of the laboratory mean inverse
speed, $\eta[v_{\min}]=\int_{v>v_{\min}}f_{\rm lab}[\mathbf v]\,\mathrm d^3v/v$,
and gives that table to RAPIDD for the nuclear recoil spectra. The lab
motion still uses the card's Solar velocity components and `vEarth_mod`.
For ID, two independent draws from the Galactic distribution give the
relative speed $u=|\mathbf v_1-\mathbf v_2|$. MadDM computes
$\langle\sigma u\rangle$ for each annihilation channel using a normalized
quadrature in $u$. The weights for the smooth–smooth, smooth–LMC and
LMC–LMC pairs are $[1-f]^2$, $2f[1-f]$ and $f^2$, where $f$ is the LMC
density fraction. The drift enters the cross pairs; it cancels from
LMC–LMC pairs.

`mw_velocity_model=shm` is the original card spelling for distribution
averaging, including the alternatives selected by `halo_velocity_profile`.
`mw_velocity_model=fixed` retains the legacy single-speed rate prescription
with `SHM` only. `mw_shm_npts` controls the number of relative-speed points
(default 16, range 8–128). Increase it and check convergence for narrow
resonances and channels close to threshold. The setting covers Galactic
antiprotons, antideuterons, the Galactic-center excess, and Galactic gamma
lines. Dwarfs keep their separate velocity input. The model applies the
local law throughout these Galactic calculations and does not change the
cosmic-ray transport tables.

The selector covers RAPIDD nuclear recoils. The native electron and
directional recoil routines still require SHM; MadDM rejects those modes
when a nonstandard profile is selected. The legacy capture module and the
experimental DD upper-limit curves retain their own halo prescriptions.

## Outputs and implementation

| File | Contents |
|---|---|
| `Halo_velocity.json` | Selected law, its active parameters, common density and both output file names |
| `SHM.dat` or `DD_halo.dat` | Annual-mean lab recoil integral passed to RAPIDD (when nuclear recoils are enabled) |
| `MW_relative_speed.dat` | Relative speeds in km/s and normalized pair weights used by Galactic ID |

`velocity_profiles.py` defines the shared radial PDFs and computes the
numerical recoil and pair integrals. `halo_model.py` dispatches the native
SHM, sets the lab boost and writes the tables. `maddm_run_interface.py`
reads the card, passes the recoil table to RAPIDD, and uses the same
relative-speed rule in the inclusive Fortran and precise MadEvent paths.
`tests/test_rapidd_velocity_profiles.py` tests normalization, pair moments,
the LMC cross term and agreement with the RAPIDD C recoil integrals at fixed
lab velocity.

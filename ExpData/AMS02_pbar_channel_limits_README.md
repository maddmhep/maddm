# BR=1 AMS-02 antiproton limits

The 19 `MadDM_AMS02_pbar_Limit_*.txt` files are byte-identical copies of the
standalone AMS-02 scan supplied in this conversation. Each has 64 native
CosmiXs mass nodes (5 GeV to 100 TeV), with two columns: mass [GeV] and the
conditional one-sided 95% upper limit [cm^3/s]. There are 1,132 finite limits
and 84 rows marked -1, denoting no signal in the supplied table/data domain.
An unavailable row never means exclusion or an allowed model.

Fixed inputs: annihilation, BR=1 for one channel, QUAINT-MED, NFW, density
0.4 GeV/cm^3 and the native propagated-table source convention. The secondary
normalization is profiled within [0.8,1.2], its IS tilt within [-0.1,0.1] about
10 GeV, and a shared solar force-field potential within [0,2] GV. All 58
supplied AMS-02 rigidity measurements and their diagonal supplied errors are
used. The endpoint has delta-chi-square 2.705543454095404 relative to the
physical best fit. Complete standalone scan metadata/hashes are retained in
`AMS02_pbar_channel_limits.json`; the original package supplies generation
scripts and all profiles.

Names follow the dwarf-table convention: `gammagamma` = AA, `gammaZ` = AZ,
`nuenue` = nu_e_nu_e, `numunumu` = nu_mu_nu_mu, `nutaunutau` = nu_tau_nu_tau.
Other channel names match the flux tables. No new particle thresholds are
imposed on nonzero off-shell content supplied in the tables.

`ExpConstraints.PBAR_max(mass, final_states, source_factor=1, rho=0.4)` loads
these limits lazily. `final_states` can be a channel name, PDG list, or dotted
PDG string. It returns -1 outside the grid or where a node/neighbor has no
limit. It returns exact stored values at the nodes and uses log-log
interpolation only between adjacent positive rows. Source/density conversion
divides the native bound by `source_factor*(rho/0.4)**2`.

In sigmav mode, MadDM compares each Milky Way partial rate with its channel
bound. Rates for the same supported channel and initial state are added before
comparison. Unsupported initial/final states remain explicit NO LIMIT rows.
No overall mixed-channel confidence bound or likelihood is obtained from
these channel tables. Use the flux-mode likelihood for a model-specific fit.

The parent `README_AMS02_PBAR.md` documents source-convention assumptions,
the missing secondary transport parameter card, absent experimental covariance,
the effective solar modulation prescription, and asymptotic coverage. These
conditional tables do not add propagation/halo marginalization or an official
AMS-02 experimental likelihood.

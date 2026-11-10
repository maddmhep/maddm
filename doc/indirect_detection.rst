.. _indirect-detection:

Indirect detection
==================

MadDM uses the annihilation channels generated for the selected UFO model.
The same launch can also compute relic density and direct detection if the
corresponding processes were generated. See the
:doc:`indirect detection tutorial <tutorial_and_examples/plot_indirect_detection>`
for a complete worked setup.

Run modes and source spectra
----------------------------

In the launch menu, choose ``indirect = sigmav``, ``flux_source``, or
``flux_earth``. The first mode reports partial annihilation rates and
compares supported channels with pure-channel dwarf and antiproton
limits. The comparison tests each theoretical partial rate against its
own limit. A collection of those comparisons does not define a combined
limit on the model's branching mixture. For configurations that differ
from the precomputed antiproton reference, MadDM profiles the selected
pure-channel table instead. The GCE fit and antideuteron forecasts can
use their own tabulated spectra even in rate-only mode.

``set fast`` selects inclusive two-body rates. When source spectra are
requested, MadDM interpolates the CosmiXs yields and combines them with
the partial rates. The published yields include the table's treatment of
off-shell weak-boson pairs. ``set precise`` uses MadEvent to calculate the
rates. With source spectra enabled, Pythia 8 or its VINCIA shower can
shower the generated final states, including explicit multibody or cascade
channels. Pythia's source antideuteron calculation uses the Argonne
wave-function prescription by default (``Main:methodDbar = 4``); it does
not change the antideuteron propagation tables. The VINCIA settings follow the CosmiXs shower and hadronization
setup. The source spectra are particle yields per annihilation before
Galactic transport or solar modulation. Install the source spectra with
``install PPPC4DMID`` and Pythia with ``install pythia8`` in MadDM as needed.
Numerical velocity reweighting uses F2PY with the Python environment that
launches MadDM. On Python 3.12 and later, install Meson and Ninja for this
build.

``flux_earth`` also writes the charged interstellar and top-of-atmosphere
fluxes at the model rate using the CosmiXsPPPC propagation
tables. This does not require a separate transport solver. The charged
tables cover supported two-body annihilations only. Generating a
multibody source spectrum in precise mode does not supply a corresponding
propagated template. Positron yields can be produced at the source; the
Earth flux awaits the positron transport and secondary tables.

Dwarf spheroidal galaxies
-------------------------

When a photon source spectrum is available, MadDM constructs a combined
Fermi-LAT likelihood using the released dwarf likelihoods from
`arXiv:2311.04982 <https://arxiv.org/abs/2311.04982>`__. The combined
limit and the two model likelihoods use the same dwarf sample. In
inclusive mode, ``fermi_fast_dwarfs = top10`` (the default) selects the
ten largest nominal J factors. Set ``fermi_fast_dwarfs = all`` for the
full sample. Precise mode uses the full sample. The selected dwarf names
and result are saved to ``Fermi_dSph_limits.json``. The continuum dwarf
velocity (``vave_indirect_cont``, default :math:`2\times10^{-5}c`) and
the dwarf J factors are independent of the Milky Way halo settings.

Propagated antiprotons and antideuterons
----------------------------------------

Choose ``cr_propagation_scheme`` from ``BIG``, ``SLIM``, and ``QUAINT``,
``cr_propagation_set`` from ``MIN``, ``MED``, and ``MAX``, and
``cr_halo_profile`` from ``NFW``, ``Einasto``, and ``Burkert``. Defaults
are QUAINT, MED, and NFW. On first continuum process generation
(``generate indirect_detection`` or ``add process indirect_detection``),
MadDM downloads the two annihilation archives from
`CosmiXsPPPC <https://github.com/CosmiXsPPPC/CosmicRay-flux>`__
(`arXiv:2606.10017 <https://arxiv.org/abs/2606.10017>`__).
It verifies the decompressed contents against the checksums in
``ExpData/CR_flux_provenance.json`` and caches all 54 tables as
gzip-compressed text in ``ExpData/CR_flux``. The same tables are reused
for subsequent commands and launches. For an existing process directory,
or to retry a failed download, run ``install cosmixs_flux`` at the
main MadDM prompt.

At a fixed mass, MadDM interpolates the table in mass at fixed energy
nodes, forms the sum weighted by model branching fractions :math:`B_f`,
and rescales by :math:`\langle\sigma v\rangle/\langle\sigma v\rangle_{\rm ref}`
and :math:`[\rho_{\rm local}/\rho_{\rm ref}]^2`. The reference rate is
:math:`3\times10^{-26}\,\mathrm{cm^3\,s^{-1}}` and the reference local
density is :math:`0.4\,\mathrm{GeV\,cm^{-3}}`. The native tables
already contain the mass-squared and propagation dependence. The
``pbar_table_convention`` card setting states the reference source
normalization explicitly. For the default self-conjugate convention, a
symmetric particle/antiparticle model has half the flux at fixed total
density and rate. The model fractions include all positive channels;
unsupported channels are reported instead of being dropped and
renormalized. The table does not distinguish polarizations or
extrapolate outside its mass range.

The bundled interstellar secondary spectra are held fixed across
MIN/MED/MAX within each propagation family. The archived BIG and SLIM
secondary inputs were byte-identical for each species. The bundled
configuration therefore uses the same stored secondary for both; the
SLIM secondary provenance has not been independently verified. A
separate ``secondary_SLIM_IS.txt`` in the corresponding species data
directory takes precedence when present. For antiprotons,
``pbar_secondary_file`` can also select an explicit secondary file.
The primary dark matter tables for BIG and SLIM remain distinct.

AMS-02 antiproton likelihood
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fit uses the absolute antiproton flux from
`AMS-02 <https://ams02.space/publications/202102>`__ and evaluates the
model at the tabulated rigidity values. It uses the reported vertical
errors as diagonal Gaussian uncertainties. It does not incorporate
experimental covariance or propagate secondary-production cross-section
uncertainties. For examples of their effect, see
`Heisig, Korsmeier and Winkler <https://arxiv.org/abs/2005.04237>`__ and
`Korsmeier, Donato and Di Mauro <https://arxiv.org/abs/1802.03030>`__.

The secondary interstellar flux has a free multiplicative factor in
:math:`[0.8,1.2]` and a kinetic-energy tilt
:math:`[T/T_{\rm pivot}]^{\delta}` with
:math:`-0.1\leq\delta\leq0.1`. The default pivot is 10 GeV.
A force-field modulation potential is profiled jointly with these
parameters (default range 0--2 GV). The code applies the same potential to
the dark matter and secondary terms and uses the kinetic-energy to
rigidity Jacobian before comparison with the data. At each fixed
nonnegative signal rate it minimizes the chi-square over the nuisance
parameters. The one-sided 95% upper endpoint uses
:math:`\Delta\chi^2=2.7055` relative to the physical best fit. The
one-sided test statistic is zero below the best fit.

``AMS02_pbar_sigmav_UL95`` is the bound for the model's fixed mixture,
``AMS02_pbar_sigmav_best_fit`` is the data-preferred rate, and
``AMS02_pbar_sigmav_model`` is the rate from the UFO parameters. The
``q_upper_model`` and ``q_upper_model_thermal`` fields test the model
rate and the rate scaled by the square of its relic fraction,
respectively. They are test statistics, not p-values. When no fit can be
evaluated, the result carries an explicit status. The default
QUAINT-MED/NFW BR=1 limit curves in ``ExpData`` may be used in
``sigmav`` mode for their matching likelihood settings; masses outside
the tabulated range or between unavailable table rows are not extrapolated.

Antideuteron forecasts
~~~~~~~~~~~~~~~~~~~~~~

MadDM converts the propagated spectrum to kinetic energy per nucleon
and applies force-field modulation. It uses the corresponding
antiproton model-point fit for the solar potential where available,
otherwise ``dbar_phi`` (default 0.7 GV). Following
`arXiv:2602.15132 <https://arxiv.org/abs/2602.15132>`__, the estimated
dark matter signal count in a window is

.. math::

   N = \frac{2}{\Phi_{\rm sens}\,[T_2-T_1]}
       \int_{T_1}^{T_2}\mathrm{d}T\,\Phi_{\bar d}^{\rm DM}[T].

The factor two reflects the two-event reference sensitivity. These are
sensitivity-based estimates; they do not fold the spectrum with a
detector response. The comparison shows the modulated secondary flux
separately; it is not added to the dark matter count.

.. list-table::
   :header-rows: 1

   * - Window
     - Kinetic energy [GeV/n]
     - Two-event flux sensitivity [m⁻² s⁻¹ sr⁻¹ (GeV/n)⁻¹]
   * - GAPS
     - 0.10--0.30
     - :math:`2\times10^{-6}`
   * - AMS-02 low
     - 0.18--0.70
     - :math:`8\times10^{-7}`
   * - AMS-02 high
     - 2.0--4.2
     - :math:`2\times10^{-7}`

The total AMS-02 count is the sum of its two windows.

Galactic Centre excess
----------------------

The prompt fit uses the standard 40° by 40° region and 0.5--1000 GeV
range of `arXiv:2605.22913 <https://arxiv.org/abs/2605.22913>`__. The
ROI-integrated SED is converted from MeV to GeV and divided by the
solid angle (about 0.478 sr) to obtain an ROI-mean intensity. For an
ROI-integrated model prediction, multiply the mean intensity by this
same solid angle. The geometry follows Eqs. 37--39 of
`arXiv:2508.02775 <https://arxiv.org/abs/2508.02775>`__:

.. math::

   \frac{\mathrm{d}\Phi_\gamma}{\mathrm{d}E} =
   \frac{r_\odot\rho_{\rm local}^2\overline J}{2\pi d_{\rm MD} m_\chi^2}
   \langle\sigma v\rangle
   \sum_f B_f\frac{\mathrm{d}N_{\gamma,f}}{\mathrm{d}E}.

The denominator uses :math:`d_{\rm MD}=4` for self-conjugate dark matter
and :math:`d_{\rm MD}=8` for symmetric particle/antiparticle dark matter.
``gce_jbar`` is a dimensionless ROI average and defaults to 120;
``gce_rsun_kpc`` defaults to 8.122 kpc. The charged-particle halo
selector does not change the GCE geometry automatically.

The fit varies a single nonnegative rate at fixed mass and branching
shape. For detected bins, the widths from the lower and upper
systematic envelopes are used as asymmetric diagonal errors.
Statistical errors are plotted separately. Upper-limit bins are shown
and checked in the output, but are not fit as Gaussian detections.
MadDM reports the best-fit rate, asymmetric
:math:`\Delta\chi^2=1` errors, and chi-square. The prompt calculation
does not include inverse-Compton or bremsstrahlung emission. A
generated prompt spectrum outside the table channels is usable only
when the run provides it at the requested fixed Galactic velocity and
it includes the relevant photon contribution.

Implementation and checks
-------------------------

``maddm_run_interface.py`` dispatches the run modes. The
``indirect/`` package contains the dwarf profiling, propagated-flux
selection, antiproton likelihood, antideuteron forecasts and GCE fit.
Run the offline regression suite from the MadDM source directory:

.. code-block:: bash

   python -m unittest discover -s tests -v
   python tests/tools/run_benchmarks.py --output /tmp/maddm-benchmarks

The suite uses small reference files in ``tests/fixtures/``.
Propagation checks that require the full tables are skipped if the cache
is empty. Run ``install cosmixs_flux`` first to include these checks and
to run the benchmark script.
Optional integration checks require a generated MadGraph run or a
locally built Pythia or RAPIDD library. See the test skip messages and
``python tests/tools/run_shower_regressions.py --help`` for the inputs.

Process directories created by older versions
---------------------------------------------

If an old process directory lacks the current Galactic rate driver, repeat
the model and process generation commands in the updated MadDM installation
and use ``output new_process``. Copy your input cards from the old
directory's ``Cards/`` directory before launching the new process. Current
exports include the required Fortran sources, makefiles and shower code;
launch does not rewrite those files in old process directories.

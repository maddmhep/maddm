Run outputs and scan plots
==========================

Each completed run has a ``MadDM_results.txt`` and a short
``README.md`` listing the outputs that were produced. For scans,
``save_output = off`` keeps the main summary and applicable JSON
results. ``spectra`` also keeps the available indirect spectra and
figures. ``all`` also keeps the detailed recoil products and event
samples. The retention mode does not change the numerical scan results.

With Matplotlib installed, MadDM plots available source yields and
direct-detection rates. With ``direct_nucleon``, ``DD_limits.pdf``
compares the SI and SD DM-nucleon cross sections of the run with the
LZ and PICO-60 limits, the XLZD projection (``xlzd_limits = 200ty`` or
``1000ty``, off by default) and the xenon neutrino fog, all rescaled to
``rhoDM``; the curves are saved in ``DD_limits.json``. The neutrino-fog
tables (O'Hare, arXiv:2109.03116) are downloaded when the direct-detection
diagrams are generated, or with ``install neutrino_fog``; ``launch`` never
downloads them. The antiproton comparison shows the secondary,
dark matter and total flux at the upper limit and a separate dark
matter curve at the model rate. The signed residual is
(model minus data) divided by model. The antideuteron comparison
shows the propagated flux, sensitivity windows and signal counts.
The GCE figure shows the systematic band, statistical errors on
detected points, and the best-fit and coupling-predicted curves.

To regenerate the propagated antinucleus and GCE figures or to plot a
single antiproton flux table:

.. code-block:: bash

   python plotting/plot_run_products.py path/to/output/run_01
   python plotting/plot_ams02_pbar_fluxes.py path/to/AMS02_pbar_fluxes.dat -o pbar.pdf

Two standalone scripts are also copied to the parent ``output/``
directory of a MadDM run, beside any scan table. Existing copies are
kept so users can retain local edits:

.. code-block:: bash

   python plot_sigmav_limits.py scan_run_01.txt --outdir figures
   python plot_coupling_mass.py --help

``plot_sigmav_limits.py`` plots available dwarf and antiproton limits,
the GCE best fit, theoretical rates and branching fractions, nuclear
and electron-recoil cross sections and bounds, relic density, freeze-out
rates, and antideuteron counts. It can locate the
lowest chi-square GCE point and integrate the saved source spectra
for multiplicities. For multiplicities, retain spectra during the
scan. Missing observables are omitted.

``plot_coupling_mass.py`` combines available direct-detection,
indirect-detection and relic-density information on a mass versus
coupling plot. For each observable :math:`O\propto g^p`, it projects
the limit from the scan point using
:math:`g_{\rm lim}=g_{\rm scan}[O_{\rm lim}/O_{\rm scan}]^{1/p}`.
The default direct- and indirect-detection powers are both four and
can be changed. The relic-density line uses a reference annihilation
rate of :math:`3\times10^{-26}\,\mathrm{cm^3\,s^{-1}}` by default.
The projection assumes smooth power-law scaling at fixed mass. It can
fail at resonances and thresholds, for significant coannihilation,
or if changing the coupling changes the branching ratios. The
relic-density line does not rerun the Boltzmann calculation.

For several scanned masses or couplings, select the x axis explicitly:

.. code-block:: bash

   python plot_sigmav_limits.py scan.txt --mass-column 'mass#9000005'
   python plot_coupling_mass.py scan.txt --mass-column 'mass#9000005' --coupling-column 'DMINPUTS#1'

Both scripts also accept ``--x-column`` for a different scanned input,
``--select COLUMN=VALUE`` (repeatable) to select rows, and
``--group-by COLUMN [COLUMN ...]`` for separate parameter slices.
If several masses are present, the last one is the reported default;
choose the physical dark-matter mass explicitly when it differs.
Repeated x values are split by the other scanned input parameters.
If repeated rows remain within a slice, they are shown as points.
The plots use base-ten logarithmic axes and vector PDF output.
Uncomputed or out-of-range negative sentinels, zeros and NaNs are not
plotted as positive constraints. The electron legend distinguishes the
published XENONnT limit, its RAPIDD atomic-response rescaling and XENON10.
The two XENONnT conventions do not both enter the coupling envelope;
choose the desired one with ``--electron-limit-column``.
A relic-only scan without a freeze-out rate uses the approximate inverse
abundance scaling instead, marked in the coupling metadata.

Every curve in the scan summaries has a distinct color, supplemented by
markers and dash patterns. The freeze-out diagnostic is labeled as the
reported annihilation rate evaluated at ``x_f``. It is not a
cross section inferred by imposing the measured relic density; the relic
abundance is computed by integrating the Boltzmann equation.

A scan can contain linked external model inputs, for example:

.. code-block:: text

   BLOCK MASS
         54 2*mxr # my0
    5000001 scan1:[1,3,4,5,7,8,10,100] # mxr

The first parameter is recalculated at every point. ``mxr/2``, ``mxr**2``
and chains of relations are also supported. Prefer the exact UFO name in
``ident_card.dat`` or use the explicit ``mass_5000001`` alias if a comment
name is ambiguous. Standard ``scanN`` groups retain MadGraph semantics:
same-number lists advance together; different groups form a Cartesian product.
Linked entries do not add independent scan dimensions. Their evaluated values
are written to the summary. A one-element scan is supported.

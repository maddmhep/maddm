Run outputs and scan plots
==========================

Each completed run has a ``MadDM_results.txt`` and a short
``README.md`` listing the outputs that were produced. For scans,
``save_output = off`` keeps the main summary and applicable JSON
results. ``spectra`` also keeps the available indirect spectra and
figures. ``all`` also keeps the detailed recoil products and event
samples. The retention mode does not change the numerical scan results.

With Matplotlib installed, MadDM plots available source yields and
direct-detection rates. The antiproton comparison shows the secondary,
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
cross sections and bounds, and antideuteron counts. It can locate the
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

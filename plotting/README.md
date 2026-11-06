# MadDM scan plots

Install NumPy and Matplotlib. Keep both scripts in the same directory. The
included plots were generated from the supplied 49-point scan.

## Limits, best-fit profile, and source multiplicities

    python plot_sigmav_limits.py /path/to/output/scan_run_03.txt --outdir figures

The script plots the dwarf, antiproton, and Galactic Center Excess (GCE)
annihilation limits, four direct-detection upper limits, two AMS-02 and one
GAPS antideuteron windows and the AMS-02 sum, and theoretical annihilation
rates and branching fractions. Add --thermal for thermal-fraction event
forecasts. All figures are vector PDFs with base-10 log scales on both axes.

The GCE upper limit is a fixed-mass limit. The script uses an explicit 95%
limit in the scan if present. Otherwise it can invert the split-normal
chi-square of the fitted bins in GCE_spectrum.dat at delta chi-square 2.7055.
The template must reproduce the scan fit. With only the scan's upper
one-sigma error, it uses best fit + sqrt(2.7055) × upper error and labels
the curve approximate. This is not a test of whether DM causes the excess.

The GCE best mass is the grid point with lowest absolute chi-square.
Mass errors come from the nearest crossings of chi-square minimum + 1,
interpolated linearly in log mass; an unbracketed side is not extrapolated.
Rate errors are the asymmetric fit errors at that mass, not a joint
two-parameter interval. If the matching run has its GCE spectrum, an SED
figure is also produced. In the sample, the best row is run_03_06:
32.6 GeV, minimum chi-square 1.26, rate 1.88e-26 cm³/s, mass crossings
25.47 and 34.33 GeV.

The script matches each run name to a subdirectory beside the scan file.
Use --runs-dir to set a different parent. It integrates dN/dlog10(x) over
the saved bins, with x = kinetic energy / DM mass. The CSV records the
integration range and the backend. The neutrino total requires all three
flavors over the same range. Use --spectrum-method to enforce the same
backend across runs and --min-x to impose a common lower cutoff. Missing
spectra are skipped; zero yields are in the CSV but omitted on the log
plot. The attached scan text did not include run folders, so the sample
figures contain neither source multiplicities nor the best-fit SED.

## Coupling versus mass

    python plot_coupling_mass.py /path/to/output/scan_run_03.txt --coupling 0.1 --dd-power 4 --id-power 4 --outdir figures

Replace 0.1 with the coupling used in your scan. Omit it to plot g/g_ref,
as in the sample, or use --coupling-column NAME if it varies across the
scan. Use one row per mass. For prediction P_ref and target T, the script
solves g = g_ref (T/P_ref)^(1/n). Both powers default to four. A zero
benchmark response gives no finite constraint under this scaling. Each DD
bound is interpreted one response at a time, without a mixed likelihood.

The shaded envelope is the strongest individual upper bound at each mass,
not a joint confidence region. The GCE best-fit band is plotted separately.
Use --no-gce-in-envelope to leave its upper limit out of the envelope.
Rescaling assumes the particle makes up all DM and that spectral shapes,
branching fractions, velocity dependence and astrophysics are fixed.
The dwarf rate defaults to tot_SM_xsec for table-based spectra; use
--fermi-rate-column tot_Xsec for full shower spectra if that is the
normalization of your run. The CSV and JSON record the input columns.

The thermal line rescales sigmav(xf) to 3e-26 cm³/s. --relic-power defaults
to the indirect power; --relic-target and --relic-rate-column override the
default target and benchmark column. This is a reference line, not a
Boltzmann calculation. It is reliable only away from resonances, thresholds,
coannihilation and other effects that invalidate the power-law or simple
thermal freeze-out approximation.

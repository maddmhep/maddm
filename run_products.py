"""Optional figures and a file guide for an individual MadDM run."""
import logging
from pathlib import Path
import shutil


logger = logging.getLogger('madgraph.plugin.maddm')
PLOTTERS = ('plot_sigmav_limits.py', 'plot_coupling_mass.py')

FILE_GUIDE = (
    ('MadDM_results.txt', 'Main results: relic abundance, rates, cross sections and limits for the enabled analyses.'),
    ('maddm_card.dat', 'Settings used for this run. `rhoDM` is the local DM density in GeV/cm³, shared by direct detection and the Milky Way signals.'),
    ('maddm.out', 'Raw output from the MadDM calculation, useful for troubleshooting.'),
    ('DDrates.txt', 'Differential nuclear recoil rates for xenon, argon and germanium versus recoil energy in keV. Rates are in events per kg, day and keV.'),
    ('DDrates.npz', 'Numerical arrays for the nuclear recoil calculation.'),
    ('DDrates.pdf', 'Log-scale comparison of the three target spectra in `DDrates.txt`.'),
    ('DDrates_e.txt', 'Differential xenon electron recoil rates from RAPIDD, summed and per shell.'),
    ('DDrates_e.npz', 'Numerical arrays for the RAPIDD electron recoil and detector signals.'),
    ('dRdS2_Xenon*_e_recoil.dat', 'RAPIDD electron recoil signal as a function of S2.'),
    ('dRdE_e_recoil.dat', 'RAPIDD differential electron recoil rate versus electron energy.'),
    ('signal_e_recoil.dat', 'RAPIDD signal per XENON10 S2 bin.'),
    ('SHM.dat', 'SHM velocity table used by RAPIDD, written when vMP, vescape or vSun_* differ from the defaults.'),
    ('rate_vs_time*.dat', 'Annual modulation (direct_modulation = True): RAPIDD recoil rate in each time bin of the year.'),
    ('dRdE_modulation*.dat', 'Annual modulation (direct_modulation = True): RAPIDD differential recoil rate in each time bin.'),
    ('Halo_velocity.json', 'Local density and velocity distribution settings for the Galactic annihilation rates.'),
    ('MW_relative_speed.dat', 'Relative speeds and quadrature weights for the Milky Way annihilation average.'),
    ('*_spectrum_*.dat', 'Source yields per annihilation: log10(x), dN/dlog10(x), and uncertainties when available, with x = kinetic energy / DM mass.'),
    ('source_spectra.pdf', 'Source yields for photons, available neutrino flavors, antiprotons, positrons and antideuterons. These are production spectra, before propagation or solar modulation.'),
    ('*_dphide_*.dat', 'Target-integrated differential photon or neutrino flux at Earth; see the file header for units.'),
    ('CosmiXsPPPC_*_flux_earth.dat', 'Propagated antiproton or antideuteron fluxes; the column names distinguish interstellar and top-of-atmosphere predictions and their components.'),
    ('CosmiXsPPPC_flux_earth.json', 'Propagation settings, table selections and status for charged fluxes at Earth.'),
    ('Fermi_dSph_limits.json', 'Combined Fermi-LAT dwarf likelihood, annihilation limit and the names of the dwarfs used.'),
    ('AMS02_pbar_limits.json', 'Antiproton constraint, fit status and diagnostics. A rate-only run can report precomputed limits by channel.'),
    ('AMS02_pbar_*', 'Antiproton flux comparisons and supporting data, when a flux fit is available.'),
    ('AMS02_dbar_events.*', 'Antideuteron event forecasts and calculation status for the AMS-02 and GAPS energy windows.'),
    ('AMS02_dbar_*', 'Antideuteron flux tables and sensitivity figure, when available.'),
    ('GCE_fit.*', 'Galactic Center Excess fit summary, diagnostics or figure. In JSON, `geometry.rho_sun_GeVcm3` is the density used in the fit.'),
    ('GCE_*spectrum.dat', 'Prompt gamma-ray spectra and the binned data used to fit or check the Galactic Center Excess.'),
    ('Output_Indirect_*', 'Link to retained MadEvent files for the corresponding indirect process group.'),
)


def copy_scan_plotters(output_directory):
    """Put both plotting scripts next to scan tables, preserving local edits."""
    source = Path(__file__).resolve().parent / 'plotting'
    for name in PLOTTERS:
        target = Path(output_directory) / name
        try:
            with (source / name).open('rb') as inp, target.open('xb') as out:
                shutil.copyfileobj(inp, out)
        except FileExistsError:
            pass


def write_run_readme(directory):
    from fnmatch import fnmatchcase

    directory = Path(directory)
    names = {path.name for path in directory.iterdir()} - {'README.md'}
    rows = []
    for pattern, meaning in FILE_GUIDE:
        matches = sorted(name for name in names if fnmatchcase(name, pattern))
        if matches:
            names.difference_update(matches)
            rows.append('| %s | %s |' % (', '.join('`%s`' % name for name in matches), meaning))
    if names:
        rows.append('| %s | Other data or diagnostics from enabled modules. |' %
                    ', '.join('`%s`' % name for name in sorted(names)))
    header = (
        '# MadDM run: %s\n\n'
        'Start with `MadDM_results.txt` for the physics summary and the status of each analysis. '
        '“All DM” assumes the model particle supplies the full dark matter density; '
        '“Thermal” uses the abundance fraction reported by the run. '
        'The files below are the products present for this run.\n\n'
        'A `sigmav` run does not produce source spectra. Flux modes save spectra supplied '
        'by the selected backend when enabled by the output settings. The log-scale figures '
        'show only positive rates or yields; a species with no positive data is omitted.\n\n'
        'For scan plots, use `plot_sigmav_limits.py` and `plot_coupling_mass.py` in the '
        'parent `output` directory. Use `save_output = spectra` or `all` in a flux scan '
        'to retain the source spectra for multiplicity plots.\n\n'
        '| File | Contents |\n| --- | --- |\n' % directory.name
    )
    (directory / 'README.md').write_text(header + '\n'.join(rows) + '\n', encoding='utf-8')


def write_run_products(directory, *, recoil, source_method=None):
    directory = Path(directory)
    if not directory.is_dir():
        return
    plots = []
    if recoil and (directory / 'DDrates.txt').is_file():
        plots.append(('DDrates.pdf', 'recoil_rates_plot', (directory / 'DDrates.txt',)))
    if source_method and any(directory.glob('*_spectrum_' + source_method + '.dat')):
        plots.append(('source_spectra.pdf', 'source_spectra_plot', (directory, source_method)))
    for name, plotter, args in plots:
        try:
            from .plotting import signal_plots
            target = directory / name
            target.unlink(missing_ok=True)
            getattr(signal_plots, plotter)(*args, target)
        except (ImportError, OSError, RuntimeError, ValueError) as error:
            logger.warning('Could not write %s: %s', name, error)
    write_run_readme(directory)

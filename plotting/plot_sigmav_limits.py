#!/usr/bin/env python3
"""Plot MadDM scan results and integrate available per-run source spectra.

Example: python plot_sigmav_limits.py scan_run_03.txt --outdir figures
If run folders are elsewhere, add --runs-dir /path/to/output.
Requires NumPy and Matplotlib.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re

import numpy as np


DD = [('SI', 'n'), ('SI', 'p'), ('SD', 'n'), ('SD', 'p')]
DD_COLORS = ['#176b9a', '#28a09d', '#ad4b22', '#d7902a']
SPECIES = {
    'gammas': r'$\gamma$',
    'neutrinos_e': r'$\nu_e$',
    'neutrinos_mu': r'$\nu_\mu$',
    'neutrinos_tau': r'$\nu_\tau$',
    'neutrinos_total': r'$\nu$ (all flavors)',
    'positrons': r'$e^+$',
    'antiprotons': r'$\bar p$',
    'antideuterons': r'$\overline{\mathrm{D}}$',
}
NU = ('neutrinos_e', 'neutrinos_mu', 'neutrinos_tau')
SOURCE_SPECIES = [s for s in SPECIES if s != 'neutrinos_total']
# 95% one-sided normal quantile squared (profile-likelihood convention).
DELTA95 = 2.705543454095404


def read_scan(path):
    """Read numbered headers, retaining run identifiers and status fields."""
    header, rows = {}, []
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        line = line.strip()
        match = re.match(r'^#\s*\[(\d+)\]\s*:\s*(.+?)\s*$', line)
        if match:
            i, name = int(match[1]) - 1, match[2]
            if i in header or name in header.values():
                raise ValueError(f'Duplicate header at line {number}')
            header[i] = name
        elif line and not line.startswith('#'):
            rows.append((number, line.replace(r'\t', '\t').split()))
    if not header or not rows:
        raise ValueError('Expected numbered MadDM headers and data rows')
    if sorted(header) != list(range(len(header))):
        raise ValueError('Header indices must start at 1 and be consecutive')
    for n, row in rows:
        if len(row) != len(header):
            raise ValueError(f'Line {n}: expected {len(header)} columns, found {len(row)}')
    columns = {}
    for i in sorted(header):
        raw = [row[i] for _, row in rows]
        try:
            columns[header[i]] = np.array([float(v.replace('D', 'e').replace('d', 'e')) for v in raw])
        except ValueError:
            columns[header[i]] = np.array(raw)
    return columns


def numeric(columns, key):
    if key not in columns:
        return None
    try:
        return np.asarray(columns[key], float)
    except (ValueError, TypeError):
        return None


def positive(a):
    a = np.array(a, float, copy=True)
    a[(~np.isfinite(a)) | (a <= 0)] = np.nan
    return a


def choose(columns, candidates):
    for key in candidates:
        value = numeric(columns, key)
        if value is not None and np.any(np.isfinite(positive(value))):
            return key
    return None


def sort_scan(columns, mass_column=None):
    if mass_column is None:
        options = [k for k in columns if k.lower().startswith('mass#')]
        if not options:
            options = [k for k in columns if k.lower() in ('mdm', 'mass', 'mchi')]
        if len(options) != 1:
            raise ValueError('Specify --mass-column when there is not exactly one mass column')
        mass_column = options[0]
    mass = numeric(columns, mass_column)
    if mass is None:
        raise ValueError('Missing numeric x column: ' + mass_column)
    good = np.flatnonzero(np.isfinite(positive(mass)))
    if not len(good):
        raise ValueError('No positive finite x values')
    order = good[np.argsort(mass[good], kind='stable')]
    return mass_column, mass[order], {k: np.asarray(v)[order] for k, v in columns.items()}


def mass_label(column):
    if column.lower().startswith('mass') or column.lower() in ('mdm', 'mchi'):
        return r'Dark matter mass $m_{\rm DM}$ [GeV]'
    return column


def plot_style(show=False):
    import matplotlib
    if not show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11, 'pdf.fonttype': 42, 'savefig.dpi': 600})
    return plt


def style_axes(ax, xlabel, ylabel):
    ax.set(xscale='log', yscale='log', xlabel=xlabel, ylabel=ylabel)
    ax.grid(which='major', alpha=.32)
    ax.grid(which='minor', alpha=.11)
    ax.tick_params(which='both', direction='in', top=True, right=True)


def save_figure(plt, fig, path, show=False):
    fig.savefig(path, dpi=600)
    if show:
        plt.show()
    plt.close(fig)
    print('Saved', path)


def run_folder(root, run):
    if root is None or run is None:
        return None
    run = str(run)
    if run in ('', '.', '..') or Path(run).name != run:
        return None
    path = Path(root) / run
    return path if path.is_dir() else None


def gce_errors(cols, n):
    sym = numeric(cols, 'GCE_sigmav_error')
    sides = []
    for key in ('GCE_sigmav_error_minus', 'GCE_sigmav_error_plus'):
        a = numeric(cols, key)
        out = np.full(n, np.nan) if a is None else np.where(np.isfinite(a) & (a >= 0), a, np.nan)
        if sym is not None:
            out = np.where(np.isfinite(out), out, np.where(sym >= 0, sym, np.nan))
        sides.append(out)
    return tuple(sides)


def gce_upper_from_bins(path, best, model_rate, ndof=None, expected_chi2=None):
    """Invert the saved split-normal objective in the fitted bins only."""
    a = np.loadtxt(path, ndmin=2)
    if a.shape[1] < 10:
        raise ValueError('GCE_spectrum.dat needs at least ten columns')
    a = a[a[:, 9] > 0]
    if ndof is not None and np.isfinite(ndof) and len(a) != ndof + 1:
        raise ValueError('GCE bin count does not match scan ndof')
    y, low, high = a[:, 3:6].T
    if np.isfinite(best) and best > 0 and np.any(a[:, 8] > 0):
        scale, template = best, a[:, 8]
    elif np.isfinite(model_rate) and model_rate > 0:
        scale, template = model_rate, a[:, 7]
    else:
        raise ValueError('No positive normalization for the GCE template')
    lo, hi = y-low, high-y
    if (not np.all(np.isfinite([template, y, lo, hi])) or np.any(lo <= 0) or
            np.any(hi <= 0) or np.any(template < 0) or not np.any(template > 0)):
        raise ValueError('Invalid GCE template or uncertainties')
    def chi(k):
        diff = template*k-y
        return np.sum((diff/np.where(diff >= 0, hi, lo))**2)
    def derivative(k):
        diff = template*k-y
        return np.sum(template*diff/np.where(diff >= 0, hi, lo)**2)
    def bisect(fun, left, right):
        for _ in range(80):
            mid = (left+right)/2
            if fun(mid) < 0:
                left = mid
            else:
                right = mid
        return (left+right)/2
    mask = template > 0
    right = max(1., 2*float(np.max(y[mask]/template[mask])))
    optimum = 0. if derivative(0) >= 0 else bisect(derivative, 0., right)
    minimum = chi(optimum)
    # MadDM's text scan rounds to about three significant figures.
    if np.isfinite(best) and best >= 0 and not np.isclose(optimum*scale, best, rtol=.025, atol=0):
        raise ValueError('GCE bin template does not reproduce the scan best fit')
    if expected_chi2 is not None and np.isfinite(expected_chi2) and not np.isclose(minimum, expected_chi2, rtol=.025, atol=.015):
        raise ValueError('GCE bin template does not reproduce the scan chi-square')
    while chi(right)-minimum < DELTA95:
        right *= 2
    return scale*bisect(lambda k: chi(k)-minimum-DELTA95, optimum, right)


def gce_limits(cols, n, runs_dir=None, best_column=None, upper_column=None):
    """Prefer a scan upper limit, then saved GCE bins, then the scan's +1 sigma error."""
    best_key = (best_column if best_column else
                next((k for k in ('GCE_sigmav_best', 'GCE_sigmav_best_fit') if numeric(cols, k) is not None), None))
    best = numeric(cols, best_key)
    if best is None:
        best = np.full(n, np.nan)
    minus, plus = gce_errors(cols, n)
    upper = np.full(n, np.nan)
    method = np.full(n, '', dtype='U24')
    raw = numeric(cols, choose(cols, [upper_column] if upper_column else
                               ['GCE_sigmav_UL95', 'GCE_sigmav_ul95', 'GCE_sigmav_upper95']))
    if raw is not None:
        valid = np.isfinite(positive(raw))
        upper[valid], method[valid] = raw[valid], 'scan_UL95'
    status = cols.get('GCE_status')
    valid_fit = np.isfinite(best) & (best >= 0)
    if status is not None:
        valid_fit &= (np.char.lower(np.asarray(status, str)) == 'ok')
    best = np.where(valid_fit, best, np.nan)
    models = numeric(cols, 'GCE_model_sigmav')
    models = np.full(n, np.nan) if models is None else models
    ndof, chi2 = numeric(cols, 'GCE_ndof'), numeric(cols, 'GCE_chi2')
    notes = []
    for i, run in enumerate(cols.get('run', [None]*n)):
        if method[i] or not valid_fit[i]:
            continue
        folder = run_folder(runs_dir, run)
        if folder is not None and (folder/'GCE_spectrum.dat').is_file():
            try:
                upper[i] = gce_upper_from_bins(folder/'GCE_spectrum.dat', best[i], models[i],
                                              None if ndof is None else ndof[i],
                                              None if chi2 is None else chi2[i])
                method[i] = 'saved_GCE_bins'
            except (OSError, ValueError) as error:
                notes.append(f'{run}: {error}')
        if not method[i] and best[i] > 0 and np.isfinite(plus[i]) and plus[i] > 0:
            upper[i] = best[i] + math.sqrt(DELTA95)*plus[i]
            method[i] = 'quadratic_approx'
    return best_key, best, minus, plus, upper, method, notes


def gce_profile(cols, mass, best):
    """Profile over repeated masses, compare like fits, and bracket each crossing."""
    chi = numeric(cols, 'GCE_chi2')
    if chi is None:
        return None
    valid = np.isfinite(best) & np.isfinite(chi) & (chi >= 0) & (best >= 0)
    ndof = numeric(cols, 'GCE_ndof')
    if ndof is not None:
        valid &= np.isfinite(ndof) & (ndof > 0)
        if len(np.unique(ndof[valid])) > 1:
            print('GCE mass fit skipped: restrict the scan to fits with the same number of fitted bins.')
            return None
    if not np.any(valid):
        return None
    j = int(np.argmin(np.where(valid, chi, np.inf)))
    compatible = valid.copy()
    if ndof is not None:
        compatible &= ndof == ndof[j]
    # A mass with no compatible fit is a gap; never interpolate across it.
    grid = np.unique(mass)
    values = np.full(len(grid), np.nan)
    for i, m in enumerate(grid):
        idx = np.flatnonzero((mass == m) & compatible)
        if len(idx):
            values[i] = np.min(chi[idx])
    center = int(np.searchsorted(grid, mass[j]))
    target = chi[j] + 1
    def crossing(direction):
        i = center
        while 0 <= i+direction < len(grid):
            k = i + direction
            if not np.isfinite(values[k]):
                return None
            if values[k] >= target:
                t = (target-values[i])/(values[k]-values[i])
                return float(np.exp((1-t)*np.log(grid[i])+t*np.log(grid[k])))
            i = k
        return None
    return j, grid, values, crossing(-1), crossing(1), int(np.sum(valid & ~compatible))


def spectrum_integral(path, min_x=None):
    """Integrate histogram densities over the tabulated log10(x) bins."""
    a = np.loadtxt(path, ndmin=2, usecols=(0, 1))
    t, y = a.T
    if len(t) < 2 or not np.all(np.isfinite(a)) or np.any(np.diff(t) <= 0) or np.any(y < 0):
        raise ValueError('Expected an ordered log10(x) grid with nonnegative yields')
    edges = np.r_[t[0]-(t[1]-t[0])/2, (t[1:]+t[:-1])/2, t[-1]+(t[-1]-t[-2])/2]
    left = edges[0] if min_x is None else max(edges[0], np.log10(min_x))
    if left >= edges[-1]:
        return None, (float(edges[0]), float(edges[-1]))
    widths = np.maximum(0., edges[1:]-np.maximum(edges[:-1], left))
    return float(np.dot(y, widths)), (float(left), float(edges[-1]))


def multiplicities(cols, mass, runs_dir, output, method=None, min_x=None):
    if min_x is not None and (not np.isfinite(min_x) or not 0 < min_x < 1):
        raise ValueError('--min-x must be between 0 and 1')
    results = {s: np.full(len(mass), np.nan) for s in SPECIES}
    rows, notes, methods = [], [], set()
    matched_folders = 0
    for i, run in enumerate(cols.get('run', [None]*len(mass))):
        folder = run_folder(runs_dir, run)
        if folder is None:
            continue
        matched_folders += 1
        candidates = {}
        for s in SOURCE_SPECIES:
            prefix = s+'_spectrum_'
            for path in folder.glob(prefix+'*.dat'):
                key = path.stem[len(prefix):]
                candidates.setdefault(key, {})[s] = path
        if method:
            use = next((k for k in sorted(candidates) if k.lower() == method.lower()), None)
        else:
            priority = {'pythia8': 0, 'vincia': 1, 'cosmixs': 2, 'pppc4dmid': 3}
            use = min(candidates, key=lambda k: (-len(candidates[k]), priority.get(k.lower(), 99), k), default=None)
        if not use:
            continue
        methods.add(use)
        bounds = {}
        for s, path in candidates[use].items():
            try:
                value, span = spectrum_integral(path, min_x)
            except (OSError, ValueError) as error:
                notes.append(f'{run}/{path.name}: {error}')
                continue
            if value is None:
                continue
            results[s][i], bounds[s] = value, span
            rows.append((run, mass[i], use, s, value, *span))
        # Report a total only if all three flavors cover the same range.
        if all(s in bounds for s in NU) and np.allclose([bounds[s] for s in NU], bounds[NU[0]], atol=1e-7, rtol=0):
            total = sum(results[s][i] for s in NU)
            results['neutrinos_total'][i] = total
            rows.append((run, mass[i], use, 'neutrinos_total', total, *bounds[NU[0]]))
    if rows:
        with output.open('w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['run', 'mass_GeV', 'method', 'species', 'N_per_annihilation', 'log10_x_min', 'log10_x_max'])
            writer.writerows(rows)
    if len(methods) > 1:
        notes.append('Multiple spectrum backends used; select --spectrum-method for a uniform comparison.')
    return results, {'matched_run_folders': matched_folders, 'integrals': len(rows), 'methods': sorted(methods),
                     'min_x': min_x, 'notes': notes, 'table': output.name if rows else None}


def channels(columns, explicit=None):
    if explicit:
        return [k for k in explicit if numeric(columns, k) is not None]
    keys = [k for k in columns if k.startswith('taacsID#')]
    if keys:
        return keys
    return [k for k in columns if 'lim_'+k in columns and
            not k.startswith(('lim_', 'AMS02_', 'GCE_', 'sigmaN_', 'sigmaP_', 'sigma_'))]


def generate(columns, outdir, stem='scan', mass_column=None, fermi_column=None, pbar_column=None,
             gce_column=None, gce_upper_column=None, total_column=None, channel_columns=None,
             thermal=False, show=False, runs_dir=None, spectrum_method=None, min_x=None):
    plt = plot_style(show)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mass_column, x, cols = sort_scan(columns, mass_column)
    repeated = len(np.unique(x)) < len(x)
    manifest = {'x_column': mass_column, 'points': len(x), 'repeated_x_scatter': repeated, 'plots': {}, 'skipped': []}
    xlabel = mass_label(mass_column)

    def plot(kind, curves, ylabel, title, band=None):
        entries = []
        for key, label, props in curves:
            y = numeric(cols, key)
            if y is not None and np.any(np.isfinite(positive(y))):
                entries.append((key, positive(y), label, props))
        if not entries:
            manifest['skipped'].append(kind)
            return
        fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
        if len(entries) > 10:
            ax.set_prop_cycle(color=plt.get_cmap('tab20').colors)
        for key, y, label, props in entries:
            fmt = dict(marker='o', ms=2.7, lw=1.65, ls='-')
            fmt.update(props)
            if repeated:
                fmt['ls'] = 'None'
            line, = ax.plot(x, y, label=label, **fmt)
            if band is not None and key == band[0] and not repeated:
                lo, hi = band[1:]
                ax.fill_between(x, positive(y-lo), positive(y+hi), color=line.get_color(), alpha=.17,
                                label=r'GCE fixed-mass $\Delta\chi^2=1$')
        style_axes(ax, xlabel, ylabel)
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=9, ncol=2 if len(entries) > 5 else 1)
        path = outdir/(stem+'_'+kind+'.pdf')
        save_figure(plt, fig, path, show)
        manifest['plots'][kind] = {'file': path.name, 'columns': [k for k, *_ in entries]}

    f = choose(cols, [fermi_column or 'Fermi_sigmav'])
    p = choose(cols, [pbar_column] if pbar_column else ['AMS02_pbar_sigmav_UL95', 'AMS02_pbar_sigmav'])
    g, best, minus, plus, upper, method, notes = gce_limits(cols, len(x), runs_dir, gce_column, gce_upper_column)
    cols['GCE_fit_curve'] = best
    cols['GCE_UL_exact'] = np.where((method != '') & (method != 'quadratic_approx'), upper, np.nan)
    cols['GCE_UL_approx'] = np.where(method == 'quadratic_approx', upper, np.nan)
    plot('sigmav_limits', [(f, 'Fermi-LAT dwarfs: upper limit', {'color':'#176b9a'}),
                           (p, 'AMS-02 antiprotons: upper limit', {'color':'#ad4b22'}),
                           ('GCE_fit_curve', 'GCE: best fit', {'color':'#3e8a65'}),
                           ('GCE_UL_exact', 'GCE: 95% upper limit', {'color':'#3e8a65','ls':'--'}),
                           ('GCE_UL_approx', 'GCE: approximate 95% upper limit', {'color':'#3e8a65','ls':':'})],
         r'$\langle\sigma v\rangle$ [cm$^3$ s$^{-1}$]', 'Annihilation limits and Galactic Center Excess fit',
         ('GCE_fit_curve', minus, plus))
    manifest['GCE_upper_limit'] = {'methods': {m:int(np.sum(method == m)) for m in sorted(set(method)-{''})},
                                   'delta_chi2': DELTA95, 'notes': notes}
    profile = gce_profile(cols, x, best)
    if profile is not None:
        j, grid, chi, lower, higher, incompatible = profile
        run = str(cols.get('run', ['']*len(x))[j])
        ndof = numeric(cols, 'GCE_ndof')
        profile_info = {'run': run, 'mass_column':mass_column, 'mass_best_scan':float(x[j]), 'mass_interval_delta_chi2_1':[lower,higher],
                        'sigmav_best_cm3_s':float(best[j]),
                        'sigmav_error_minus_cm3_s':float(minus[j]) if np.isfinite(minus[j]) else None,
                        'sigmav_error_plus_cm3_s':float(plus[j]) if np.isfinite(plus[j]) else None,
                        'chi2':float(numeric(cols,'GCE_chi2')[j]),
                        'ndof_at_fixed_mass':int(ndof[j]) if ndof is not None else None,
                        'incompatible_ndof_points':incompatible,
                        'mass_interval_method':'nearest crossings; piecewise linear interpolation in log mass'}
        manifest['GCE_best_fit'] = profile_info
        manifest['plots']['gce_best_fit'] = plot_gce_profile(plt, outdir, stem, x, best, minus, plus, profile_info,
                                                           grid, chi, xlabel, show)
        folder = run_folder(runs_dir, run)
        if folder is not None and (folder/'GCE_spectrum.dat').is_file():
            try:
                rates = numeric(cols, 'GCE_model_sigmav')
                gce_upper_from_bins(folder/'GCE_spectrum.dat', best[j], np.nan if rates is None else rates[j],
                                    profile_info['ndof_at_fixed_mass'], profile_info['chi2'])
                manifest['plots']['gce_best_fit_sed'] = plot_gce_sed(plt, folder, outdir, stem, profile_info, show)
            except (OSError, ValueError) as error:
                notes.append(f'{run}: best-fit SED not plotted: {error}')
    else:
        manifest['skipped'].append('gce_best_fit')

    dd_curves = [('lim_sigmaN_'+a+'_'+n, a+' '+n, {'color':c,'ls':'-' if n=='n' else '--'})
                 for (a,n),c in zip(DD,DD_COLORS)]
    plot('dd_limits', dd_curves, r'Nucleon cross-section upper limit [cm$^2$]', 'Direct-detection upper limits')
    events = [('GAPS_dbar_events','GAPS (0.1–0.3 GeV/n)'), ('AMS02_dbar_events_low','AMS-02 low (0.18–0.70 GeV/n)'),
              ('AMS02_dbar_events_high','AMS-02 high (2.0–4.2 GeV/n)'), ('AMS02_dbar_events_total','AMS-02 total')]
    for suffix in ('', '_thermal') if thermal else ('',):
        key = 'AMS02_dbar_events_total'+suffix
        low, high = (numeric(cols, 'AMS02_dbar_events_'+w+suffix) for w in ('low','high'))
        if low is not None and high is not None:
            existing = numeric(cols, key)
            existing = np.full(len(x), np.nan) if existing is None else existing
            cols[key] = np.where((existing >= 0) & np.isfinite(existing), existing,
                                 np.where((low >= 0)&(high >= 0), low+high, np.nan))
    event_curves = [(k,l,{}) for k,l in events]
    if thermal:
        event_curves += [(k+'_thermal',l+'; thermal',{'ls':':'}) for k,l in events]
    plot('dbar_events', event_curves, 'Expected signal events', 'Antideuteron signal forecasts')
    total = choose(cols, [total_column] if total_column else ['tot_Xsec','taacsID','tot_SM_xsec'])
    mw = choose(cols, ['AMS02_pbar_sigmav_model','AMS02_pbar_model_sigmav','AMS02_dbar_model_sigmav','GCE_model_sigmav'])
    plot('sigmav_total', [(total,'Theory total (continuum scan)',{}),(mw,'Theory total (Milky Way)',{'ls':'--'})],
         r'$\langle\sigma v\rangle_{\mathrm{tot}}$ [cm$^3$ s$^{-1}$]', 'Theoretical total annihilation rate')
    brcurves = []
    if total:
        denominator = positive(numeric(cols, total))
        for k in channels(cols, channel_columns):
            values = numeric(cols, k)
            if values is None:
                continue
            cols['BR#'+k] = np.where(values >= 0, values/denominator, np.nan)
            brcurves.append(('BR#'+k, k.split('#',1)[-1], {}))
    plot('branching_ratios', brcurves, r'$\langle\sigma v\rangle_i/\langle\sigma v\rangle_{\rm tot}$', 'Channel branching fractions')
    manifest['branching_denominator'] = total

    mult, info = multiplicities(cols, x, runs_dir, outdir/(stem+'_source_multiplicities.csv'), spectrum_method, min_x)
    cols.update({'N#'+s: v for s,v in mult.items()})
    main_species = ['gammas','neutrinos_total','positrons','antiprotons','antideuterons']
    if any(np.any(np.isfinite(mult[s]) & ~np.isfinite(mult['neutrinos_total'])) for s in NU):
        for s in NU:
            cols['N#'+s] = np.where(np.isfinite(mult['neutrinos_total']), np.nan, mult[s])
        main_species += list(NU)
    plot('source_multiplicities', [('N#'+s,SPECIES[s],{}) for s in main_species],
         'Particles per annihilation', 'Integrated source multiplicities')
    manifest['source_multiplicities'] = info
    (outdir/(stem+'_plot_manifest.json')).write_text(json.dumps(manifest, indent=2, allow_nan=False)+'\n')
    if not info['matched_run_folders']:
        print('Source multiplicities: no matching run folders; supply --runs-dir if they are elsewhere.')
    for kind in manifest['skipped']:
        print('Skipped', kind, '(no usable data)')
    return manifest


def fit_text(info, with_chi=False):
    m, lo, hi = info['mass_best_scan'], *info['mass_interval_delta_chi2_1']
    if lo is not None and hi is not None:
        mass = r'$m_{\rm DM}=%.1f^{+%.1f}_{-%.1f}$ GeV' % (m, hi-m, m-lo)
    else:
        mass = r'$m_{\rm DM}=%.3g$ GeV; ' % m
        mass += ('interval not bracketed' if lo is None and hi is None else
                 ('lower crossing %.3g GeV' % lo if hi is None else 'upper crossing %.3g GeV' % hi))
    r, a, b = (info['sigmav_best_cm3_s'], info['sigmav_error_minus_cm3_s'], info['sigmav_error_plus_cm3_s'])
    if r > 0 and a is not None and b is not None:
        exponent = math.floor(math.log10(r))
        q = 10.**exponent
        rate = r'$\langle\sigma v\rangle=(%.3g^{+%.3g}_{-%.3g})\times10^{%d}$ cm$^3$ s$^{-1}$' % (r/q,b/q,a/q,exponent)
    else:
        rate = r'$\langle\sigma v\rangle=%.3g$ cm$^3$ s$^{-1}$' % r
    result = mass+'\n'+rate
    if with_chi:
        result += '\n'+r'$\chi^2_{\min}=%.3g$' % info['chi2']
    return result


def plot_gce_profile(plt, outdir, stem, x, best, minus, plus, info, grid, chi, xlabel, show):
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(8, 7.7), sharex=True, constrained_layout=True,
                                     gridspec_kw={'height_ratios':[1.2,1]})
    top.set_title('Galactic Center Excess: best fit on the scanned mass grid')
    repeated = len(np.unique(x)) < len(x)
    top.plot(x, positive(best), 'o' if repeated else 'o-', ms=3, color='#287d61', label='Best fit at each mass')
    if not repeated:
        top.fill_between(x, positive(best-minus), positive(best+plus), color='#287d61', alpha=.16)
    # Draw the same linear-in-log-mass profile used to locate the crossings.
    # Connecting raw points directly on two log axes instead interpolates log(chi2).
    dense = np.geomspace(grid[0], grid[-1], max(200, 25*len(grid)))
    brackets = [edge for edge in info['mass_interval_delta_chi2_1'] if edge is not None]
    dense = np.unique(np.r_[dense, grid, brackets])
    bottom.plot(dense, positive(np.interp(np.log(dense), np.log(grid), chi)), color='#176b9a')
    bottom.plot(grid, positive(chi), 'o', ms=3, color='#176b9a')
    bottom.axhline(info['chi2']+1, color='#79535e', ls='--', label=r'$\chi^2_{\rm min}+1$')
    center = info['mass_best_scan']
    lower, upper = info['mass_interval_delta_chi2_1']
    for ax in (top, bottom):
        ax.axvline(center, color='#222222', ls=':', lw=1)
        for edge in info['mass_interval_delta_chi2_1']:
            if edge is not None:
                ax.axvline(edge, color='#787878', ls='--', lw=.8)
    if info['sigmav_best_cm3_s'] > 0:
        top.scatter(center, info['sigmav_best_cm3_s'], color='#222222', s=34, zorder=4)
    style_axes(top, '', r'Best-fit $\langle\sigma v\rangle$ [cm$^3$ s$^{-1}$]')
    style_axes(bottom, xlabel, r'Profile $\chi^2$')
    if lower is not None and upper is not None:
        for ax in (top, bottom):
            ax.axvspan(lower, upper, color='#6e7e88', alpha=.10, zorder=0)
        xmin, xmax = max(grid[0], lower/1.65), min(grid[-1], upper*1.65)
        bottom.set_xlim(xmin, xmax)
        # The full mass range remains visible on the annihilation limits plot.
        near = (x >= xmin) & (x <= xmax)
        local = positive(np.r_[best[near], (best-minus)[near], (best+plus)[near]])
        if np.isfinite(local).any():
            top.set_ylim(np.nanmin(local)/1.25, np.nanmax(local)*1.35)
        local_chi = positive(chi[(grid >= xmin) & (grid <= xmax)])
        if np.isfinite(local_chi).any():
            bottom.set_ylim(np.nanmin(local_chi)/1.25, max(np.nanmax(local_chi),info['chi2']+1)*1.25)
    top.text(.97,.06,fit_text(info),ha='right',va='bottom',transform=top.transAxes,
             bbox=dict(facecolor='white',edgecolor='#dddddd',alpha=.94),fontsize=10)
    bottom.legend(frameon=False, loc='lower right')
    path = outdir/(stem+'_gce_best_fit.pdf')
    save_figure(plt, fig, path, show)
    return {'file':path.name}


def plot_gce_sed(plt, folder, outdir, stem, info, show):
    a = np.loadtxt(folder/'GCE_spectrum.dat', ndmin=2)
    if a.shape[1] < 10:
        raise ValueError('GCE_spectrum.dat needs at least ten columns')
    y = a[:,3]
    used = (a[:,9]>0) & (y>0)
    fig, ax = plt.subplots(figsize=(8,5.5),constrained_layout=True)
    ax.fill_between(a[:,0],positive(a[:,4]),positive(a[:,5]),color='#999999',alpha=.2,label='Systematic band')
    err = a[used,10] if a.shape[1] > 10 and np.all(a[used,10] >= 0) else None
    ax.errorbar(a[used,0],y[used],yerr=err,fmt='o',ms=3,color='#444444',label='Fitted data')
    censored = (a[:,9]<=0) & (a[:,6]>0)
    if np.any(censored):
        ax.scatter(a[censored,0],a[censored,6],marker='v',facecolors='none',edgecolors='#555555',label='Published upper limits')
    if (folder/'GCE_prompt_spectrum.dat').is_file():
        e, f = np.loadtxt(folder/'GCE_prompt_spectrum.dat', usecols=(0,2), ndmin=2).T
    else:
        e, f = a[:,0], a[:,8]
    ax.plot(e,positive(f),color='#287d61',lw=2,label='Best-fit prompt DM')
    style_axes(ax,r'Photon energy [GeV]',r'$E^2 d\Phi/dE$ [GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
    ax.set_title('Galactic Center Excess at the best scanned mass')
    ax.text(.96,.95,fit_text(info,True),ha='right',va='top',transform=ax.transAxes,
            bbox=dict(facecolor='white',edgecolor='#dddddd',alpha=.94),fontsize=10)
    ax.legend(frameon=False,loc='lower left')
    path = outdir/(stem+'_gce_best_fit_sed.pdf')
    save_figure(plt,fig,path,show)
    return {'file':path.name}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--outdir', type=Path)
    parser.add_argument('--mass-column','--x-column',dest='mass_column')
    for name in ('fermi','pbar','gce','gce-upper','total'):
        parser.add_argument('--'+name+'-column')
    parser.add_argument('--channel-columns',nargs='+',help='Override auto-detected partial-rate columns')
    parser.add_argument('--thermal',action='store_true',help='Also show relic-rescaled event forecasts')
    parser.add_argument('--runs-dir',type=Path,help='Directory containing the per-point run folders (default: scan directory)')
    parser.add_argument('--spectrum-method',help='Restrict source integration to this spectrum backend')
    parser.add_argument('--min-x',type=float,help='Lower integration cut on kinetic energy / DM mass')
    parser.add_argument('--show',action='store_true')
    args = parser.parse_args()
    kwargs = {k:v for k,v in vars(args).items() if k not in ('input','outdir')}
    kwargs['runs_dir'] = args.runs_dir or args.input.parent
    try:
        generate(read_scan(args.input),args.outdir or args.input.parent,args.input.stem,**kwargs)
    except (ValueError,OSError,ImportError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()

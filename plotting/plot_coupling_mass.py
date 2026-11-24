#!/usr/bin/env python3
"""Rescale MadDM bounds to the coupling–mass plane for a one-coupling model.

Example: python plot_coupling_mass.py scan_run_03.txt --coupling 0.1 --dd-power 4 --id-power 4
Use the actual coupling of the scan as --coupling. Omit it for g/g_ref.
Standalone script. Requires NumPy and Matplotlib.
With multiple mass columns, the last one is selected by default and reported.
Use --mass-column 'mass#9000005' to choose explicitly.
The CSV retains all scanned masses. Missing analyses are skipped.
Electron cross sections use --electron-power (defaults to --dd-power).
The rescaled electron limit is shown as an alternative to the published limit.
Use --electron-limit-column lim_sigma_e_rescaled to include it in the envelope.
Use --mass-column 'mass#9000005' (or --x-column) to choose the axis.
Use --select 'BLOCK#1=0.1' to fix parameters; repeat for additional selections.
Repeated x values are split by the other scanned inputs automatically.
Use --group-by COLUMN [COLUMN ...] to choose slices explicitly.
Only available analyses are plotted; nonpositive/sentinel values are omitted on log axes.
"""
import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np

DD = [('SI', 'n'), ('SI', 'p'), ('SD', 'n'), ('SD', 'p')]
DD_COLORS = ['#176b9a', '#28a09d', '#ad4b22', '#d7902a']
DELTA95 = 2.705543454095404

# Deliberately separated hues; continuous curves use solid/dashed line styles.
CURVE_COLORS = ('#0072B2', '#D55E00', '#009E73', '#CC79A7', '#6F3C9F',
                '#E69F00', '#56B4E9', '#8C564B', '#444444', '#B2A000',
                '#E7298A', '#00A6A6', '#A65628', '#5F8A30', '#A23B72', '#7B83CC')


CURVE_LINESTYLES = ('-', '--', '-.', ':')


def curve_style(index):
    if index < len(CURVE_COLORS):
        color = CURVE_COLORS[index]
    else:
        import colorsys
        rgb = colorsys.hsv_to_rgb(((index-len(CURVE_COLORS))*.61803398875+.07)%1., .72, .72)
        color = '#'+''.join('%02x' % round(255*channel) for channel in rgb)
    return color, None


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


def plot_style(show=False):
    import matplotlib
    if not show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 11, 'pdf.fonttype': 42, 'savefig.dpi': 600})
    return plt


def style_axes(ax, xlabel, ylabel):
    ax.set_xscale('log', base=10)
    ax.set_yscale('log', base=10)
    ax.set(xlabel=xlabel, ylabel=ylabel)
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


def parameter_columns(columns):
    """MadDM writes scanned inputs before xsi and calculated observables."""
    result = []
    for key in columns:
        low = key.lower()
        if low in ('xsi', 'omegah2', 'omegah^2', 'omega', 'x_f', 'sigmav(xf)') or low.startswith(
                ('sigma', 'lim_', 'taacs', 'err_', 'tot_', 'fermi_', 'ams02_', 'gce_', 'gaps_', 'pvalue', 'like_')):
            break
        if low != 'run' and numeric(columns, key) is not None:
            result.append(key)
    return result


def sort_scan(columns, mass_column=None):
    """Use a reported default, with an explicit override for any scan axis."""
    if mass_column is None:
        options = [k for k in columns if k.lower().startswith('mass#')]
        if not options:
            options = [k for k in columns if k.lower() in ('mdm', 'mass', 'mchi')]
        if not options:
            options = parameter_columns(columns)
        if not options:
            raise ValueError('No scan parameter found. Specify --x-column COLUMN.')
        mass_column = options[-1]
        if len(options) > 1:
            print('Available axes: '+', '.join(options)+'. Using '+mass_column+
                  '; override with --mass-column or --x-column.')
    mass = numeric(columns, mass_column)
    if mass is None:
        raise ValueError('Missing numeric x column: '+mass_column)
    good = np.flatnonzero(np.isfinite(mass) & (mass > 0))
    if not len(good):
        raise ValueError('No positive finite x values in '+mass_column)
    order = good[np.argsort(mass[good], kind='stable')]
    return mass_column, mass[order], {k: np.asarray(v)[order] for k,v in columns.items()}


def mass_label(column):
    if column.lower().startswith('mass#'):
        return 'Mass [GeV] ('+column+')'
    if column.lower() in ('mass', 'mdm', 'mchi'):
        return r'Dark matter mass $m_{\rm DM}$ [GeV]'
    return column


def scan_slices(columns, axis=None, select=None, group_by=None):
    """Keep unrelated parameter slices separate instead of connecting a grid."""
    mask = np.ones(len(next(iter(columns.values()))), dtype=bool)
    for expression in select or ():
        key, sep, value = expression.partition('=')
        if not sep or key not in columns:
            raise ValueError('--select requires an existing COLUMN=VALUE: '+expression)
        values = numeric(columns, key)
        if values is None:
            mask &= np.asarray(columns[key], str) == value
        else:
            try:
                target = float(value)
            except ValueError:
                raise ValueError('Expected a numeric selection for '+key)
            mask &= np.isclose(values, target, rtol=1e-7, atol=0.)
    if not mask.any():
        raise ValueError('No rows match --select.')
    axis, x, cols = sort_scan({k: np.asarray(v)[mask] for k,v in columns.items()}, axis)
    if group_by is None:
        group_by = []
        if len(np.unique(x)) < len(x):
            for key in parameter_columns(cols):
                if key != axis and any(len(np.unique(cols[key][x == m])) > 1 for m in np.unique(x)):
                    group_by.append(key)
    for key in group_by:
        if key not in cols or key == axis:
            raise ValueError('--group-by must name existing parameters other than the x axis: '+key)
    groups = {}
    for i in range(len(x)):
        signature = tuple(str(cols[key][i]) for key in group_by)
        groups.setdefault(signature, []).append(i)
    result = []
    for number, (signature, indices) in enumerate(groups.items(), 1):
        context = ', '.join(k+'='+v for k,v in zip(group_by, signature))
        suffix = ('_slice%02d_' % number)+re.sub(r'[^a-zA-Z0-9_.-]+', '_', context)[:120] if len(groups)>1 else ''
        result.append((suffix, context, {k: np.asarray(v)[indices] for k,v in cols.items()}))
    return axis, result


def dd_label(kind, nucleon, suffix=''):
    experiment = ('PICO-60 (2019)' if kind == 'SD' and nucleon == 'p' else 'LZ (2024)')
    if suffix:
        # _XLZD: projected XLZD sensitivity (maddm_card xlzd_limits = 200ty or 1000ty)
        experiment = {'_LZ': 'LZ (2024)', '_XLZD': 'XLZD (projected)'}.get(
            suffix.upper(), suffix.strip('_').replace('_', ' '))
    target = 'proton' if nucleon == 'p' else 'neutron'
    return experiment+'; '+kind+' '+target


# DM-nucleon channels of the plots: (prediction column, legend label, [(limit column, label)], color). With the
# equivalent cross sections of MadDM (sigmaN_<SI|SD_n|SD_p>_eq: the cross section of the reference model of each LZ
# limit with the same expected LZ signal as the model) the LZ and XLZD limits are compared with them, and the PICO-60
# SD-proton limit with sigmaN_SD_p. Scans written before (no sigmaN_*_eq column) compare sigma_p and sigma_n.
DD_EQ = (('SI_eq', 'SI'), ('SD_n_eq', 'SD neutron'), ('SD_p_eq', 'SD proton'))
DD_EQ_COLORS = ['#176b9a', '#28a09d', '#ad4b22', '#d7902a']


def dd_channels(cols):
    if any('sigmaN_' + key in cols for key, _ in DD_EQ):
        out = []
        for (key, what), color in zip(DD_EQ, DD_EQ_COLORS):
            rate = 'sigmaN_' + key
            limits = [(lim, experiment + '; ' + what) for lim, experiment in
                      (('lim_' + rate, 'LZ (2024)'), ('lim_' + rate + '_XLZD', 'XLZD (projected)')) if lim in cols]
            out.append((rate, r'Prediction; $\sigma_{\rm eq}$ ' + what, limits, color))
        out.append(('sigmaN_SD_p', r'Prediction; $\sigma_{\rm SD}$ proton',
                    [(lim, 'PICO-60 (2019); SD proton') for lim in ('lim_sigmaN_SD_p',) if lim in cols], DD_EQ_COLORS[3]))
        return out
    out = []
    for (kind, nucleon), color in zip(DD, DD_COLORS):
        rate = 'sigmaN_' + kind + '_' + nucleon
        limit = 'lim_' + rate
        limits = [(key, dd_label(kind, nucleon, key[len(limit):])) for key in cols
                  if key == limit or key.startswith(limit + '_')]
        out.append((rate, 'Prediction; ' + kind + ' ' + ('proton' if nucleon == 'p' else 'neutron'), limits, color))
    return out


def fog_predictions(cols):
    """(neutrino-fog column, prediction compared with it, label): the equivalent cross sections when available"""
    eq = any('sigmaN_' + key in cols for key, _ in DD_EQ)
    return [('neutrino_fog_SI', 'sigmaN_SI_eq' if eq else 'sigmaN_SI_p', 'Xe neutrino fog; SI' + ('' if eq else ' (equal p/n)')),
            ('neutrino_fog_SD_p', 'sigmaN_SD_p_eq' if eq else 'sigmaN_SD_p', 'Xe neutrino fog; SD proton'),
            ('neutrino_fog_SD_n', 'sigmaN_SD_n_eq' if eq else 'sigmaN_SD_n', 'Xe neutrino fog; SD neutron')]


def electron_label(key):
    return {'lim_sigma_e': 'XENONnT; published convention',
            'lim_sigma_e_rescaled': 'XENONnT; RAPIDD atomic responses',
            'lim_sigma_e_Xenon10': 'XENON10; S2-only'}.get(key, key.replace('lim_sigma_e_', '').replace('_', ' '))


def rescale(reference, prediction, target, power):
    """Solve target = prediction (g/g_ref)^power; skip zero/unavailable response."""
    if prediction is None or target is None:
        return np.full(len(reference), np.nan)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore', under='ignore'):
        return positive(reference * np.exp((np.log(positive(target))-np.log(positive(prediction)))/power))


def _generate(input_path, *, outdir=None, mass_column=None, coupling=None, coupling_column=None,
             dd_power=4., id_power=4., relic_power=None, relic_target=3e-26, relic_rate_column='sigmav(xf)',
             fermi_rate_column=None, pbar_rate_column=None, gce_rate_column=None, runs_dir=None,
             gce_upper_column=None, no_gce_in_envelope=False, electron_power=None,
             electron_limit_column=None, xlabel=None, show=False, columns=None, stem=None, context=''):
    input_path = Path(input_path)
    mass_column, x, cols = sort_scan(read_scan(input_path) if columns is None else columns, mass_column)
    repeated = len(np.unique(x)) != len(x)
    if coupling is not None and coupling_column is not None:
        raise ValueError('Use either --coupling or --coupling-column.')
    if coupling_column:
        ref = numeric(cols, coupling_column)
        if ref is None:
            raise ValueError('Missing numeric benchmark coupling column: '+coupling_column)
    else:
        ref = np.full(len(x), 1. if coupling is None else coupling)
    if np.any(~np.isfinite(ref) | (ref <= 0)):
        raise ValueError('Benchmark couplings must be positive and finite.')
    normalized = coupling is None and coupling_column is None
    relic_power = id_power if relic_power is None else relic_power
    electron_power = dd_power if electron_power is None else electron_power
    if any(not np.isfinite(p) or p <= 0 for p in (dd_power, electron_power, id_power, relic_power, relic_target)):
        raise ValueError('Powers and thermal target must be positive and finite.')
    runs_dir = Path(runs_dir) if runs_dir is not None else input_path.parent
    root = Path(outdir) if outdir is not None else input_path.parent
    root.mkdir(parents=True,exist_ok=True)

    curves, upper_keys, dd_lines = {}, [], []
    references = {}
    targets = {}
    def add(key, rate_key, target, power, is_upper=False):
        y = rescale(ref, numeric(cols, rate_key), target, power)
        curves[key] = y
        references[key] = rate_key
        if is_upper and np.isfinite(y).any():
            upper_keys.append(key)
        return y
    # coupling at which each prediction reaches its limit: the equivalent cross sections for LZ / XLZD,
    # sigma_SD_p for PICO-60 (sigma_p, sigma_n for the scans written before the equivalent cross sections)
    for rate, _, limits, color in dd_channels(cols):
        for limit_key, label in limits:
            key = limit_key[len('lim_sigmaN_'):]
            add(key, rate, numeric(cols, limit_key), dd_power, 'XLZD' not in limit_key)
            targets[key] = limit_key
            dd_lines.append((key, label, color, ':' if 'XLZD' in limit_key else '--'))

    for column, rate, label in fog_predictions(cols):
        add(column, rate, numeric(cols, column), dd_power, False)
        targets[column] = column
        dd_lines.append((column, label, '#777777', ':'))

    electron_limits = [key for key in cols if key == 'lim_sigma_e' or key.startswith('lim_sigma_e_')]
    if electron_limit_column is not None and numeric(cols, electron_limit_column) is None:
        raise ValueError('Missing numeric electron limit column: ' + electron_limit_column)
    primary_electron = electron_limit_column or choose(cols, ['lim_sigma_e', 'lim_sigma_e_rescaled'])
    if electron_limit_column and electron_limit_column not in electron_limits:
        electron_limits.append(electron_limit_column)
    for index, limit_key in enumerate(electron_limits):
        suffix = limit_key[len('lim_sigma_e'):].strip('_')
        key = 'electron' + ('_' + suffix if suffix else '')
        # The published and rescaled curves are two conventions for the same
        # experiment. The selected convention enters the exclusion envelope.
        alternative = limit_key in ('lim_sigma_e', 'lim_sigma_e_rescaled')
        include = limit_key == primary_electron or not alternative
        add(key, 'sigma_e', numeric(cols, limit_key), electron_power, include)
        targets[key] = limit_key
        dd_lines.append((key, electron_label(limit_key),
                         ['#a12b89', '#cc72b2', '#7353a6'][index % 3],
                         ':' if limit_key == 'lim_sigma_e_rescaled' else '-.'))
    fermi_rate = fermi_rate_column or choose(cols, ['Fermi_model_sigmav','tot_SM_xsec','tot_Xsec','taacsID'])
    pbar_rate = pbar_rate_column or choose(cols, ['AMS02_pbar_sigmav_model','AMS02_pbar_model_sigmav'])
    gce_rate = gce_rate_column or choose(cols, ['GCE_model_sigmav'])
    pbar_limit = choose(cols, ['AMS02_pbar_sigmav_UL95','AMS02_pbar_sigmav'])
    add('dwarfs', fermi_rate, numeric(cols, 'Fermi_sigmav'), id_power, True)
    add('pbar', pbar_rate, numeric(cols, pbar_limit), id_power, True)
    _, best, minus, plus, upper, methods, notes = gce_limits(cols,len(x),runs_dir,upper_column=gce_upper_column)
    add('gce_best', gce_rate, best, id_power)
    add('gce_lo68', gce_rate, best-minus, id_power)
    add('gce_hi68', gce_rate, best+plus, id_power)
    add('gce_ul95', gce_rate, upper, id_power, not no_gce_in_envelope)
    if choose(cols, [relic_rate_column]):
        add('relic_target', relic_rate_column, relic_target, relic_power)
        relic_method = 'freezeout_rate'
    else:
        omega = choose(cols, ['Omegah^2','Omegah2','omegah2','Omega','omega'])
        abundance = numeric(cols, omega)
        inverse = None if abundance is None else 1./positive(abundance)
        cols['inverse_relic_abundance'] = inverse if inverse is not None else np.full(len(x), np.nan)
        add('relic_target', 'inverse_relic_abundance', 1./.12, relic_power)
        relic_method = 'inverse_abundance_approximation'


    envelope = np.full(len(x),np.nan)
    sections = [np.arange(len(x))]
    if upper_keys:
        stack = np.array([curves[k] for k in upper_keys])
        available = np.any(np.isfinite(stack),axis=0)
        envelope[available] = np.nanmin(stack[:,available],axis=0)
        # A limit that ends at the edge of its tabulated coverage must not
        # be interpolated towards the bound from a different experiment.
        owner = np.argmin(np.where(np.isfinite(stack),stack,np.inf),axis=0)
        left, right = np.arange(len(x)-1), np.arange(1,len(x))
        gap = (owner[left] != owner[right]) & (~np.isfinite(stack[owner[left],right]) |
                                              ~np.isfinite(stack[owner[right],left]))
        sections = np.split(np.arange(len(x)), np.flatnonzero(gap)+1)
    usable = any(np.isfinite(y).any() for y in curves.values())
    if not usable:
        print('No usable coupling constraints'+(' for '+context if context else '')+
              '; writing the table and metadata without a PDF.')
    curves['upper_envelope'] = envelope
    prefix = root/((stem or input_path.stem)+'_coupling_mass')
    table = Path(str(prefix)+'.csv')
    mass_columns = [key for key in cols if key != mass_column and key.lower().startswith('mass#')]
    input_parameters = [key for key in parameter_columns(cols) if key != mass_column]
    with table.open('w',newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['run',mass_column,*input_parameters,'benchmark_coupling',*curves,'gce_ul95_method'])
        for i,m in enumerate(x):
            writer.writerow([cols.get('run',['']*len(x))[i],m,*[cols[key][i] for key in input_parameters],
                             '' if normalized else ref[i],
                             *[y[i] if np.isfinite(y[i]) else '' for y in curves.values()],methods[i]])

    if usable:
        plt = plot_style(show)
        fig, ax = plt.subplots(figsize=(9.4,6.4),constrained_layout=True)
        plotted_styles = []
        def draw(key,label,color,ls='-',lw=1.6):
            y = curves[key]
            if np.isfinite(y).any():
                color, marker = curve_style(len(plotted_styles))
                ls = CURVE_LINESTYLES[len(plotted_styles) % len(CURVE_LINESTYLES)]
                plotted_styles.append(dict(curve=key,color=color,marker=marker,linestyle=ls))
                ax.plot(x,y,marker=None,lw=lw,label=label,color=color,ls=ls)
        for key, label, color, linestyle in dd_lines:
            draw(key, label, color, linestyle)
        draw('dwarfs','Fermi-LAT dwarfs','#654c9c','--')
        draw('pbar','AMS-02 antiprotons','#bd412f','-.')
        draw('gce_best','GCE best fit','#278559','-')
        valid = np.isfinite(curves['gce_lo68']) & np.isfinite(curves['gce_hi68'])
        if valid.any() and not repeated:
            band_color = next((style['color'] for style in plotted_styles if style['curve']=='gce_best'), '#009E73')
            ax.fill_between(x,curves['gce_lo68'],curves['gce_hi68'],where=valid,alpha=.14,color=band_color,
                            label=r'GCE fixed-mass $\Delta\chi^2=1$')
        if np.isfinite(curves['gce_ul95']).any():
            approximate = methods=='quadratic_approx'
            curves['gce_exact_plot'] = np.where(~approximate,curves['gce_ul95'],np.nan)
            curves['gce_approx_plot'] = np.where(approximate,curves['gce_ul95'],np.nan)
            draw('gce_exact_plot','GCE upper limit','#278559','--')
            draw('gce_approx_plot','GCE upper limit (approx.)','#278559',':')
            del curves['gce_exact_plot'],curves['gce_approx_plot']
        draw('relic_target',r'Reference from $\langle\sigma v\rangle(x_f)$','#65574e','--',1.9)
        if np.isfinite(envelope).any() and not repeated:
            first = True
            for segment in sections:
                if np.isfinite(envelope[segment]).any():
                    ax.plot(x[segment], envelope[segment], '-', lw=2.3, color='#181b20',
                            label='Strongest available upper bound' if first else None)
                    first = False
            style_axes(ax,mass_label(mass_column),'')
            ax.autoscale_view()
            ymin,ymax = ax.get_ylim()
            for segment in sections:
                ax.fill_between(x[segment],envelope[segment],ymax,where=np.isfinite(envelope[segment]),
                                color='#445264',alpha=.08,zorder=0)
            ax.set_ylim(ymin,ymax)
        axis_label = xlabel or mass_label(mass_column)
        if np.isfinite(curves['relic_target']).any():
            axis_label += '\nFreeze-out target approximate near resonances and thresholds.'
        style_axes(ax,axis_label,r'Coupling ratio $g/g_{\rm ref}$' if normalized else r'Coupling $g$')
        if any(key.startswith('neutrino_fog_') for key in curves if np.isfinite(curves[key]).any()):
            fig.supxlabel('NeutrinoFog: github.com/cajohare/NeutrinoFog; O’Hare (2021), arXiv:2109.03116\n'
                          'Xe reference boundary; not an exclusion', fontsize=7)
        ax.set_title('Coupling constraints from a fixed-shape rescaling'+('\n'+context if context else ''))
        ax.legend(frameon=False,fontsize=8.5,ncol=2,loc='best')
        save_figure(plt,fig,Path(str(prefix)+'.pdf'),show)
    meta = {'reference_coupling':coupling,'coupling_column':coupling_column,'axis_is_g_over_gref':normalized,
            'mass_column':mass_column,'other_mass_columns':mass_columns, 'input_parameters':input_parameters,
            'plot_created':usable, 'plot_styles':plotted_styles if usable else [],
            'slice':context,'relic_target_method':relic_method, 'repeated_x_scatter':repeated,
            'powers':{'direct':dd_power,'electron':electron_power,'indirect':id_power,'freezeout':relic_power},
            'rate_columns':references,'limit_columns':targets,
            'electron_limit_column':primary_electron,
            'direct_reference_sources': {'XLZD': 'RAPIDD lib/experiments/XLZD: SI arXiv:2410.17137, SD arXiv:2203.02309; projected median 90% CL, rescaled to rhoDM', 'neutrino_fog': 'https://github.com/cajohare/NeutrinoFog; arXiv:2109.03116; xenon, rescaled to rhoDM; reference, not exclusion'},
            'envelope_inputs':upper_keys,'relic_target_cm3_s':relic_target,
            'GCE_upper_limit_methods':{k:int(np.sum(methods==k)) for k in sorted(set(methods)-{''})},
            'notes':notes+['Envelope = strongest available individual constraint, no joint likelihood.',
                           'Full DM density assumed; spectra, branching fractions and halo settings fixed at each mass.',
                           'Thermal target is an approximate freeze-out criterion, unreliable near thresholds, resonances or with coannihilation.']}
    Path(str(prefix)+'.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
    return curves,meta


def generate(input_path, *, select=None, group_by=None, **kwargs):
    input_path = Path(input_path)
    axis, slices = scan_slices(read_scan(input_path), kwargs.get('mass_column'), select, group_by)
    kwargs['mass_column'] = axis
    results = {}
    for suffix, context, cols in slices:
        results[context or 'all'] = _generate(input_path, columns=cols, stem=input_path.stem+suffix,
                                              context=context, **kwargs)
    if len(slices) == 1:
        return next(iter(results.values()))
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path)
    p.add_argument('--outdir',type=Path)
    p.add_argument('--mass-column','--x-column',dest='mass_column')
    p.add_argument('--xlabel', help='Override the mass-axis label')
    p.add_argument('--electron-power', type=float, help='Electron cross-section scaling exponent (default: --dd-power)')
    p.add_argument('--electron-limit-column', help='Electron convention included in the envelope (default: published)')
    group = p.add_mutually_exclusive_group()
    group.add_argument('--coupling','--g-ref',type=float,help='Benchmark coupling used to produce the scan')
    group.add_argument('--coupling-column',help='Use a per-mass benchmark coupling column instead')
    for name in ('dd','id','relic'):
        p.add_argument('--'+name+'-power',type=float,default=None if name=='relic' else 4.)
    p.add_argument('--relic-target',type=float,default=3e-26,help='Thermal target [cm^3/s]')
    p.add_argument('--relic-rate-column',default='sigmav(xf)')
    for name in ('fermi','pbar','gce'):
        p.add_argument('--'+name+'-rate-column')
    p.add_argument('--gce-upper-column')
    p.add_argument('--runs-dir',type=Path)
    p.add_argument('--no-gce-in-envelope',action='store_true',help='Show GCE separately from the exclusion envelope')
    p.add_argument('--select', action='append', help='Fix COLUMN=VALUE; repeat as needed')
    p.add_argument('--group-by', nargs='+', help='Separate curves by these scanned parameters')
    p.add_argument('--show',action='store_true')
    try:
        args = vars(p.parse_args())
        generate(args.pop('input'),**args)
    except (ValueError,OSError,ImportError) as error:
        p.error(str(error))


if __name__ == '__main__':
    main()

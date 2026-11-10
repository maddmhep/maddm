"""MadDM-to-AMS-02 adapter. Imports MadGraph neither directly nor indirectly.

Rates are present-day partial sigma*v in cm^3/s, evaluated at the separate
Milky Way velocity. Never use relic fractions or normalize away a missing
channel. Only one DM species and inclusive two-body SM channels are supported.
"""
import json
import math
import re
from . import cr_config
from pathlib import Path


def card_value(card, name, default=None):
    """Read dicts and MadGraph ConfigFile without assuming dict.get semantics.

    MG5 ConfigFile aliases get to __getitem__, so get(key, default) raises
    TypeError. Indexing preserves its case/auto handling and explicit zeroes.
    """
    try:
        return card[name]
    except KeyError:
        return default


def check_generated_switches(project_dir, indirect, spectral, inclusive):
    text = (Path(project_dir)/'include/maddm_card.inc').read_text()
    expected = dict(do_indirect_detection=bool(indirect), do_indirect_spectral=bool(spectral),
                    only2to2lo=bool(inclusive))
    for key, value in expected.items():
        found = re.findall(r'(?im)^\s*'+key+r'\s*=\s*\.(true|false)\.', text)
        if len(found) != 1 or (found[0].lower() == 'true') != value:
            raise RuntimeError('Generated '+key+' disagrees with the selected MadDM mode')
    return expected


PDG_CHANNELS = {
    (-11, 11): 'ee', (-13, 13): 'mumu', (-15, 15): 'tautau',
    (-12, 12): 'nu_e_nu_e', (-14, 14): 'nu_mu_nu_mu', (-16, 16): 'nu_tau_nu_tau',
    (-1, 1): 'dd', (-2, 2): 'uu', (-3, 3): 'ss', (-4, 4): 'cc',
    (-5, 5): 'bb', (-6, 6): 'tt', (21, 21): 'gg', (22, 22): 'AA',
    (-24, 24): 'WW', (23, 23): 'ZZ', (25, 25): 'hh', (22, 23): 'AZ', (23, 25): 'hZ',
}
# Stable numerical observables, including unavailable sentinels, for scans.
NUMERIC_KEYS = (
    'AMS02_pbar_sigmav_UL95', 'AMS02_pbar_sigmav_best_fit', 'AMS02_pbar_sigmav_model',
    'AMS02_pbar_q_upper_model', 'AMS02_pbar_q_upper_model_thermal',
    'AMS02_pbar_sigmav', 'AMS02_pbar_sigmav_best', 'AMS02_pbar_model_sigmav',
    'AMS02_pbar_q', 'AMS02_pbar_q_thermal', 'AMS02_pbar_excluded',
    'AMS02_pbar_excluded_thermal', 'AMS02_pbar_chi2_background', 'AMS02_pbar_chi2_best',
    'AMS02_pbar_norm_upper', 'AMS02_pbar_slope_upper', 'AMS02_pbar_phi_upper',
    'AMS02_pbar_velocity',
)


def parse_fortran_record(line, processes_names_map):
    """Read full-width/legacy indirect records without splitting UFO names."""
    match = re.match(r'^\s*(?:mw_shm_sv|pbar_sv|sigma\*v):\s*(\S+)\s+(\S+)', line)
    if not match:
        raise ValueError('Malformed indirect rate record: '+line.rstrip())
    raw, value = match.groups()
    name = resolve_process_name(raw, processes_names_map)
    rate = float(value.replace('D', 'e').replace('d', 'e'))
    if not math.isfinite(rate) or rate < 0:
        raise ValueError('Invalid annihilation rate for '+name)
    return name, rate


def resolve_process_name(raw, processes_names_map):
    if raw in processes_names_map:
        return raw
    unprefixed = re.sub(r'^(?:P)?\d+_', '', raw)
    if unprefixed in processes_names_map:
        return unprefixed
    matches = [name for name in processes_names_map
               if raw.startswith(name+'_') or raw.startswith(name+'/')]
    if matches:
        return max(matches, key=len)
    # Read an old A12 field only if its truncated prefix identifies one process.
    if len(raw) == 12:
        matches = [name for name in processes_names_map if name.startswith(raw)]
        if len(matches) == 1:
            return matches[0]
        if matches:
            raise ValueError('Ambiguous truncated process name '+raw+'; recompile the updated driver')
    return raw


CHANNEL_FILES = {ch: 'MadDM_AMS02_pbar_Limit_'+label+'.txt' for ch, label in (
    ('ee', 'ee'), ('mumu', 'mumu'), ('tautau', 'tautau'),
    ('nu_e_nu_e', 'nuenue'), ('nu_mu_nu_mu', 'numunumu'), ('nu_tau_nu_tau', 'nutaunutau'),
    ('uu', 'uu'), ('dd', 'dd'), ('ss', 'ss'), ('cc', 'cc'), ('bb', 'bb'), ('tt', 'tt'),
    ('AA', 'gammagamma'), ('gg', 'gg'), ('WW', 'WW'), ('ZZ', 'ZZ'),
    ('hh', 'hh'), ('AZ', 'gammaZ'), ('hZ', 'hZ'))}


def source_scaling(card, majorana_dirac_factor):
    """Flux scaling shared by the precomputed tables and the full likelihood."""
    if majorana_dirac_factor not in (4, 8):
        raise ValueError('Unknown MadDM particle/antiparticle normalization')
    convention = card_value(card, 'pbar_table_convention', 'self_conjugate')
    if convention not in ('self_conjugate', 'particle_antiparticle'):
        raise ValueError('Unknown pbar_table_convention')
    source_factor = (1. if convention == 'self_conjugate' else 2.)/(majorana_dirac_factor/4.)
    rho = cr_config.local_density(card)
    if not math.isfinite(rho) or rho <= 0:
        raise ValueError('Require positive finite rhoDM')
    return convention, source_factor, rho


class SingleChannelLimits:
    """BR=1 table lookup, with no interpolation across absent rows/extrapolation.

    Exact native nodes are returned unchanged; between two adjacent positive
    rows use linear interpolation in log(mass), log(sigmav95).
    """
    def __init__(self, directory=None):
        self.directory = Path(directory) if directory else Path(__file__).resolve().parents[1]/'ExpData'
        self.tables = {}

    def upper_limit(self, mass, channel, source_factor=1., rho=.4):
        import numpy as np
        if not isinstance(channel, str) or channel not in CHANNEL_FILES:
            try:
                states = channel.split('.') if isinstance(channel, str) else channel
                channel = PDG_CHANNELS.get(tuple(sorted(int(x) for x in states)))
            except (ValueError, TypeError):
                return -1.
        if channel not in CHANNEL_FILES or not math.isfinite(float(mass)):
            return -1.
        if not math.isfinite(source_factor) or source_factor <= 0 or not math.isfinite(rho) or rho <= 0:
            raise ValueError('Invalid antiproton table normalization')
        if channel not in self.tables:
            data = np.loadtxt(self.directory/CHANNEL_FILES[channel], usecols=(0, 1))
            if (data.ndim != 2 or data.shape[0] < 2 or not np.all(np.isfinite(data))
                    or np.any(data[:, 0] <= 0) or np.any(np.diff(data[:, 0]) <= 0)
                    or np.any((data[:, 1] <= 0) & (data[:, 1] != -1.))):
                raise ValueError('Invalid BR=1 antiproton table for '+channel)
            self.tables[channel] = data
        data = self.tables[channel]
        i = int(np.searchsorted(data[:, 0], mass))
        if i < len(data) and data[i, 0] == mass:
            value = float(data[i, 1])
        elif i == 0 or i == len(data) or np.any(data[i-1:i+1, 1] <= 0):
            return -1.
        else:
            value = float(np.exp(np.interp(np.log(mass), np.log(data[i-1:i+1, 0]),
                                          np.log(data[i-1:i+1, 1]))))
        return value/(source_factor*(rho/.4)**2) if value > 0 else -1.


def evaluate_channel_limits(mass, records, dm_pdg, majorana_dirac_factor, xsi, card, lookup, *, profiled=False):
    """Compare MW partial rates individually; never construct a mixture bound."""
    for name, expected in [('pbar_phi_min', 0.), ('pbar_phi_max', 2.),
                           ('pbar_pivot', 10.), ('pbar_rigidity_min', 0.)]:
        if not profiled and float(card_value(card, name, expected)) != expected:
            raise ValueError('BR=1 tables require the default '+name+'='+str(expected)+
                             '; use a flux mode for a new profiled likelihood')
    if not profiled and float(card_value(card, 'pbar_rigidity_max', -1.)) > 0:
        raise ValueError('BR=1 tables use all AMS-02 points; use a flux mode for a rigidity cut')
    if not math.isfinite(float(xsi)) or not 0 < float(xsi) <= 1:
        raise ValueError('Require 0 < xsi <= 1')
    convention, factor, rho = source_scaling(card, majorana_dirac_factor)
    expected = sorted([abs(int(dm_pdg))]*2 if majorana_dirac_factor == 4
                      else [-abs(int(dm_pdg)), abs(int(dm_pdg))])
    groups = {}
    for r in records:
        rate = float(r['sigmav_cm3s'])
        if not math.isfinite(rate) or rate < 0:
            raise ValueError('Invalid partial annihilation rate')
        final = tuple(sorted(int(x) for x in r['final_states']))
        supported = (sorted(int(x) for x in r['initial_states']) == expected
                     and not r.get('decay_chains') and final in PDG_CHANNELS)
        key = final if supported else ('unsupported', r['process'])
        if key not in groups:
            groups[key] = dict(channel=PDG_CHANNELS.get(final), final_states=list(final),
                               processes=[], sigmav_cm3s=0., supported=supported)
        groups[key]['processes'].append(r['process'])
        groups[key]['sigmav_cm3s'] += rate
    rows = []
    for row in groups.values():
        row['process'] = '+'.join(row['processes'])
        upper = lookup(mass, row['final_states'], source_factor=factor, rho=rho) if row['supported'] else -1.
        row.update(sigmav95_cm3s=upper, thermal_sigmav_cm3s=row['sigmav_cm3s']*float(xsi)**2,
                   excluded=None if upper <= 0 else row['sigmav_cm3s'] > upper,
                   excluded_thermal=None if upper <= 0 else row['sigmav_cm3s']*float(xsi)**2 > upper,
                   status='ok' if upper > 0 else 'no_limit')
        if upper <= 0:
            row['reason'] = ('unsupported initial/final state or decay chain' if not row['supported'] else
                             'outside table range or no signal in adjacent native table rows')
        rows.append(row)
    return dict(status='tabulated_channels', evaluation='single_channel_tables', mass_GeV=float(mass),
                channels=rows, rates=records, total_sigmav_cm3s=sum(r['sigmav_cm3s'] for r in rows),
                sigmav95_cm3s=None, table_convention=convention, source_factor=factor,
                rho_local_GeVcm3=rho, xsi=float(xsi), propagation=cr_config.selection(card)['propagation'], halo=cr_config.selection(card)['halo'],
                bound_source='runtime pure-channel profiles' if profiled else 'bundled BR=1 tables',
                combined_limit_evaluated=False,
                warnings=['Individual channel comparisons; no combined mixed-channel exclusion or likelihood.',
                          'BR=1 bounds fix secondary normalization to [0.8,1.2], tilt to [-0.1,0.1], and profile phi in [0,2] GV.',
                          'The source convention is an explicit card assumption.'])


def channel_mixture(records, dm_pdg, majorana_dirac_factor):
    """Aggregate independent partial rates without losing unsupported states."""
    if majorana_dirac_factor not in (4, 8):
        raise ValueError('Unknown MadDM particle/antiparticle normalization')
    expected = sorted([abs(int(dm_pdg))]*2 if majorana_dirac_factor == 4
                      else [-abs(int(dm_pdg)), abs(int(dm_pdg))])
    partial = {}
    errors = []
    total = 0.
    for record in records:
        rate = float(record['sigmav_cm3s'])
        if not math.isfinite(rate) or rate < 0:
            raise ValueError('Invalid partial rate for ' + record.get('process', '?'))
        if rate == 0:
            continue
        total += rate
        ini = sorted(int(x) for x in record['initial_states'])
        final = tuple(sorted(int(x) for x in record['final_states']))
        channel = PDG_CHANNELS.get(final)
        if (ini != expected or record.get('decay_chains') or channel is None):
            errors.append({'process': record.get('process', '?'), 'initial_states': ini,
                           'final_states': list(final), 'sigmav_cm3s': rate,
                           'reason': 'initial state, decay chain, or final state is unsupported'})
            continue
        partial[channel] = partial.get(channel, 0.) + rate
    if errors:
        return {'status': 'unsupported_channels', 'unsupported': errors, 'total_sigmav_cm3s': total}
    if total <= 0:
        return {'status': 'zero_rate', 'total_sigmav_cm3s': total}
    return {'status': 'ok', 'partial_rates_cm3s': partial, 'total_sigmav_cm3s': total,
            'branching': {ch: rate/total for ch, rate in partial.items()}}


def empty_observables(velocity):
    values = dict.fromkeys(NUMERIC_KEYS, -1.)
    values['AMS02_pbar_velocity'] = float(velocity)
    values['AMS02_pbar_status'] = 'not_evaluated'
    return values


class AntiprotonLimits:
    """Lazy, reusable table reader. No dependence on gamma-ray spectrum settings."""
    def __init__(self):
        self.table = None
        self.likelihoods = {}
        self.limits = {}  # only exact mass/mixture/config keys; no rounded reuse

    def configure(self, card):
        from .pbar.core import FitConfig, FluxTable, Likelihood
        path = cr_config.table_path('pbar', card)
        secondary = cr_config.secondary_path(card)
        if self.table is None or self.table.path != path:
            self.table = FluxTable(path)
        lo, hi = float(card_value(card,'pbar_phi_min',0.)), float(card_value(card,'pbar_phi_max',2.))
        rmin, rmax = float(card_value(card,'pbar_rigidity_min',0.)),float(card_value(card,'pbar_rigidity_max',-1.))
        pivot = float(card_value(card,'pbar_pivot',10.))
        import hashlib
        key = (lo, hi, rmin, rmax, pivot, str(secondary), hashlib.sha256(secondary.read_bytes()).hexdigest())
        if key not in self.likelihoods:
            self.likelihoods[key] = Likelihood(FitConfig(phi_bounds_GV=(lo,hi),
                rigidity_min_GV=rmin, rigidity_max_GV=rmax if rmax>0 else math.inf,pivot_GeV=pivot),secondary_path=secondary)
        return self.likelihoods[key], (str(path), key)

    def channel_limit(self, mass, final_states, *, source_factor, rho, card):
        likelihood, config_key = self.configure(card)
        channel = PDG_CHANNELS.get(tuple(sorted(map(int,final_states))))
        if channel is None: return -1.
        template = self.table.template(mass,{channel:1.},rho_local_GeVcm3=rho,source_factor=source_factor)
        key = (float(mass), ((channel,1.),), rho, source_factor, config_key)
        if key not in self.limits:
            self.limits[key] = likelihood.upper_limit(template)
        result = self.limits[key]
        upper = result.get('sigmav95_cm3s')
        return float(upper) if upper is not None and upper>0 else -1.

    def evaluate(self, mass, records, dm_pdg, majorana_dirac_factor, xsi, card):
        mixture = channel_mixture(records, dm_pdg, majorana_dirac_factor)
        if mixture['status'] != 'ok':
            return mixture
        likelihood, config_key = self.configure(card)
        convention, source_factor, rho = source_scaling(card, majorana_dirac_factor)
        if not math.isfinite(float(xsi)) or not 0. < float(xsi) <= 1.:
            raise ValueError('Require 0 < xsi <= 1 for the thermal density rescaling')
        template = self.table.template(mass, mixture['branching'], rho_local_GeVcm3=rho,
                                       source_factor=source_factor)
        key = (float(mass), tuple(sorted(mixture['branching'].items())), rho,
               source_factor, config_key)
        if key not in self.limits:
            if len(self.limits) >= 16:
                self.limits.pop(next(iter(self.limits)))
            self.limits[key] = likelihood.upper_limit(template)
        # Copy before attaching point-specific rates/profiles to a reusable limit.
        import copy
        result = copy.deepcopy(self.limits[key])
        result.update(mixture)
        # The upper_limit status must survive an empty tabulated signal.
        result['status'] = self.limits[key]['status']
        result['rates'] = records
        result['secondary_path'] = str(likelihood.secondary_path)
        result['propagation'] = cr_config.selection(card)['propagation']
        result['halo'] = cr_config.selection(card)['halo']
        result['warnings'].append('Secondary spectrum follows BIG/QUAINT/SLIM and is shared across MIN/MED/MAX; only the DM propagated table changes with MIN/MED/MAX.')
        result['table_convention'] = convention
        result['source_factor'] = source_factor
        result['xsi'] = float(xsi)
        result['velocity'] = float(card_value(card, 'vave_indirect_pbar', .00075))
        result['warnings'].append('Table source convention is an explicit card assumption, not inferred from the native table format.')
        if result['sigmav95_cm3s'] is None:
            return result
        rate = mixture['total_sigmav_cm3s']
        best = result['best_fit']
        points = {}
        for label, value in [('all_dm', rate), ('thermal', rate*float(xsi)**2)]:
            fit = likelihood.profile(template, value)
            if fit.chi2 < best['chi2']-1.e-4:
                raise RuntimeError('Model point improves the scanned global minimum; the limit needs a denser signal scan')
            q = max(0., fit.chi2-best['chi2']) if value >= best['sigmav_cm3s'] else 0.
            points[label] = dict(fit.to_dict(), q_upper=q,
                       excluded=bool(q > likelihood.config.delta_chi2))
        result['model_points'] = points
        return result


def observables(result, velocity):
    values = empty_observables(velocity)
    values['AMS02_pbar_status'] = result.get('status', 'not_evaluated')
    values['AMS02_pbar_model_sigmav'] = result.get('total_sigmav_cm3s', -1.)
    values['AMS02_pbar_sigmav_model'] = values['AMS02_pbar_model_sigmav']
    if result.get('status') == 'tabulated_channels':
        for row in result['channels']:
            name = row['channel'] if row['supported'] else row['process']
            values['AMS02_pbar_model#'+name] = row['sigmav_cm3s']
            values['AMS02_pbar_limit#'+name] = row['sigmav95_cm3s']
            for flag in ('excluded', 'excluded_thermal'):
                values['AMS02_pbar_'+flag+'#'+name] = -1. if row[flag] is None else float(row[flag])
        return values
    if result.get('status') != 'ok' or result.get('sigmav95_cm3s') is None:
        return values
    b, u = result['best_fit'], result['upper_limit']
    a, t = result['model_points']['all_dm'], result['model_points']['thermal']
    values.update(AMS02_pbar_sigmav=result['sigmav95_cm3s'], AMS02_pbar_sigmav_best=b['sigmav_cm3s'],
                  AMS02_pbar_q=a['q_upper'], AMS02_pbar_q_thermal=t['q_upper'],
                  AMS02_pbar_excluded=float(a['excluded']), AMS02_pbar_excluded_thermal=float(t['excluded']),
                  AMS02_pbar_chi2_background=result['background']['chi2'], AMS02_pbar_chi2_best=b['chi2'],
                  AMS02_pbar_norm_upper=u['norm'], AMS02_pbar_slope_upper=u['slope'],
                  AMS02_pbar_phi_upper=u['phi_GV'])
    values.update(AMS02_pbar_sigmav_UL95=result['sigmav95_cm3s'],
                  AMS02_pbar_sigmav_best_fit=b['sigmav_cm3s'],AMS02_pbar_sigmav_model=values['AMS02_pbar_model_sigmav'],
                  AMS02_pbar_q_upper_model=a['q_upper'],AMS02_pbar_q_upper_model_thermal=t['q_upper'])
    return values


def save_plot(out_dir, enabled=True):
    """Use the antiproton plotting script after a successful flux-file write."""
    out = Path(out_dir)
    pdf = out/'AMS02_pbar_comparison.pdf'
    pdf.unlink(missing_ok=True)
    data = out/'AMS02_pbar_fluxes.dat'
    if not enabled or not data.is_file():
        return None
    if "." in __package__:
        from ..plotting import plot_ams02_pbar_fluxes as pbar_plot
    else:
        from plotting import plot_ams02_pbar_fluxes as pbar_plot
    values = pbar_plot.read_fluxes(data, 'upper')
    model = pbar_plot.read_fluxes(data, 'model')
    label=None
    metadata=out/'AMS02_pbar_limits.json'
    if metadata.is_file():
        meta=json.loads(metadata.read_text())
        conf=meta.get('template',{})
        label=meta.get('propagation',conf.get('propagation',''))+' / '+meta.get('halo',conf.get('halo',''))
        if meta.get('secondary_path'):
            label+='\nSecondary: '+Path(meta['secondary_path']).stem.replace('_IS','').replace('secondary_','')+' reference'
    pbar_plot.make_plot(values, pdf, prediction='upper', residual_scale='linear',model_dm=model['dm'],setup=label)
    return str(pdf)


def write_outputs(out_dir, result, engine=None, *, detailed=True):
    """Always write a diagnostic JSON, including unavailable/error status.

    In detailed mode a successful fit also writes flux components per GV.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir/'AMS02_pbar_comparison.pdf').unlink(missing_ok=True)
    (out_dir/'AMS02_pbar_limits.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    channel_file = out_dir/'AMS02_pbar_channel_limits.dat'
    if detailed and result.get('status') == 'tabulated_channels':
        lines = ['# channel mass_GeV sigmav_MW_cm3s sigmav95_cm3s excluded_AllDM excluded_Thermal',
                 '# -1 means no limit/decision; pure-channel comparisons, no combined mixture bound']
        for row in result['channels']:
            flags = [-1 if row[key] is None else int(row[key]) for key in ('excluded', 'excluded_thermal')]
            lines.append('%s %.12e %.12e %.12e %d %d' %
                         (row['channel'] if row['supported'] else row['process'], result['mass_GeV'],
                          row['sigmav_cm3s'], row['sigmav95_cm3s'], flags[0], flags[1]))
        channel_file.write_text('\n'.join(lines)+'\n')
    elif channel_file.exists():
        channel_file.unlink()
    if not detailed or result.get('status') != 'ok' or not result.get('sigmav95_cm3s') or engine is None:
        stale = out_dir/'AMS02_pbar_fluxes.dat'
        if stale.exists():
            stale.unlink()
        return
    from .pbar.core import FitConfig
    import numpy as np
    c = dict(result['config'])
    if c['rigidity_max_GV'] is None:
        c['rigidity_max_GV'] = math.inf
    likelihood = next(l for l in engine.likelihoods.values() if l.config == FitConfig(**c)
                      and (not result.get('secondary_path') or str(l.secondary_path)==result['secondary_path']))
    template = engine.table.template(result['template']['mass_GeV'], result['branching'],
                         rho_local_GeVcm3=result['template']['rho_local_GeVcm3'],
                         source_factor=result['source_factor'])
    blocks = [likelihood.rigidity, likelihood.observed, likelihood.errors]
    for point in [result['model_points']['all_dm'], result['model_points']['thermal'],
                  result['best_fit'], result['upper_limit']]:
        sec, dm = likelihood.components(template, point['sigmav_cm3s'],
                               point['norm'], point['slope'], point['phi_GV'])
        blocks.extend([sec, dm, sec+dm])
    np.savetxt(out_dir/'AMS02_pbar_fluxes.dat', np.column_stack(blocks), fmt='%.12e',
          header='R_GV data_per_GV error_per_GV '
                 'model_sec model_dm model_total thermal_sec thermal_dm thermal_total '
                 'best_sec best_dm best_total upper_sec upper_dm upper_total\n'
                 'Flux units: GV^-1 m^-2 s^-1 sr^-1; TOA; each prediction uses its own profiled nuisance parameters')


# Documentation shared by the result writer and regression checks.
LEGACY_KEYS = {'AMS02_pbar_sigmav','AMS02_pbar_sigmav_best','AMS02_pbar_model_sigmav',
               'AMS02_pbar_q','AMS02_pbar_q_thermal'}
RESULT_COMMENTS = {
 'AMS02_pbar_sigmav_UL95': 'cm^3/s; one-sided 95% CL upper limit, fixed model branching fractions',
 'AMS02_pbar_sigmav_best_fit': 'cm^3/s; pbar-data best fit after profiling secondary and solar parameters',
 'AMS02_pbar_sigmav_model': 'cm^3/s; theoretical total at the model couplings and Milky Way velocity',
 'AMS02_pbar_q_upper_model': 'one-sided profiled Delta chi2 at theoretical rate; exclusion threshold 2.7055434541',
 'AMS02_pbar_q_upper_model_thermal': 'same statistic at xsi^2 times theoretical rate; equals AllDM when xsi=1',
}

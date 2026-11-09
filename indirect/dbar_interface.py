"""Two-event sensitivity forecasts from propagated CosmiXs Dbar fluxes.

The native table uses TOTAL kinetic energy K and flux per dK. Experimental
windows use T=K/A and flux per dT, requiring the Jacobian dK/dT=A=2.
These are DM signal forecasts, not an exclusion test or a detector simulation.
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from . import pbar_interface
from . import cr_config
from .pbar.core import FluxTable, SIGMAV_REF, force_field


MASS_GEV = 1.875613
A = 2.
DATA_DIR = Path(__file__).resolve().parents[1]/'ExpData/AMS02_dbar'
WINDOWS = {'AMS02_low': (0.18, 0.70, 8.e-7),
           'AMS02_high': (2.0, 4.2, 2.e-7),
           'GAPS': (0.1, 0.3, 2.e-6)}
NUMERIC_KEYS = ('AMS02_dbar_events_low', 'AMS02_dbar_events_high', 'AMS02_dbar_events_total',
                'AMS02_dbar_events_low_thermal', 'AMS02_dbar_events_high_thermal',
                'AMS02_dbar_events_total_thermal', 'GAPS_dbar_events', 'GAPS_dbar_events_thermal',
                'AMS02_dbar_phi', 'AMS02_dbar_phi_thermal', 'AMS02_dbar_model_sigmav')


def flux_per_nucleon(template, t_GeV_per_n, sigmav_cm3s, phi_GV):
    """TOA dPhi/dT: K_IS=2*T+|Z|*phi, p_TOA^2/p_IS^2, and dK/dT=2."""
    t = np.asarray(t_GeV_per_n, float)
    if (np.any(~np.isfinite(t)) or np.any(t <= 0) or not math.isfinite(phi_GV)
            or phi_GV < 0 or not math.isfinite(sigmav_cm3s) or sigmav_cm3s < 0):
        raise ValueError('Require positive T/n, nonnegative finite rate and phi')
    if sigmav_cm3s == 0:
        return np.zeros_like(t)
    return A*sigmav_cm3s/SIGMAV_REF*force_field(template, A*t, phi_GV, mass_GeV=MASS_GEV)


def window_forecast(flux, low, high, sensitivity, points=()):
    if not 0 < low < high or not math.isfinite(sensitivity) or sensitivity <= 0:
        raise ValueError('Invalid forecast window or sensitivity')
    knots = [float(t) for t in points if low < t < high]
    integral, error = quad(lambda t: float(flux(t)), low, high, points=knots,
                           epsabs=1.e-30, epsrel=1.e-7, limit=max(100, len(knots)+10))
    if not math.isfinite(integral) or integral < 0:
        raise ValueError('Invalid integrated Dbar signal flux')
    average = integral/(high-low)
    return dict(T_min_GeV_per_n=low, T_max_GeV_per_n=high,
                sensitivity_per_GeV_per_n=sensitivity, reference_events=2.,
                integrated_flux=integral, integration_absolute_error=error,
                average_flux_per_GeV_per_n=average, expected_events=2.*average/sensitivity)


def modulation(card, pbar_result):
    fallback = float(pbar_interface.card_value(card, 'dbar_phi', .7))
    if not math.isfinite(fallback) or fallback < 0:
        raise ValueError('dbar_phi must be finite and nonnegative [GV]')
    if pbar_result and pbar_result.get('status') == 'ok':
        points = pbar_result.get('model_points', {})
        if all(name in points for name in ('all_dm', 'thermal')):
            values = {name: float(points[name]['phi_GV']) for name in ('all_dm', 'thermal')}
            if any(not math.isfinite(v) or v < 0 for v in values.values()):
                raise ValueError('Invalid fitted pbar solar potential')
            return values, 'pbar model-point nuisance fits (all-DM and thermal separately)'
    return dict(all_dm=fallback, thermal=fallback), 'fixed dbar_phi; no pbar model-point flux fit available'


class AntideuteronEvents:
    def __init__(self):
        self.table = None

    def evaluate(self, mass, records, dm_pdg, majorana_dirac_factor, xsi, card, pbar_result=None):
        mixture = pbar_interface.channel_mixture(records, dm_pdg, majorana_dirac_factor)
        if mixture['status'] not in ('ok', 'zero_rate'):
            return mixture
        xsi = float(xsi)
        if not math.isfinite(xsi) or not 0 < xsi <= 1:
            raise ValueError('Require 0 < xsi <= 1')
        convention, factor, rho = pbar_interface.source_scaling(card, majorana_dirac_factor)
        path = cr_config.table_path('dbar', card)
        if self.table is None or self.table.path != path:
            self.table = FluxTable(path)
        # A zero-rate model needs no arbitrary branching mixture.
        template = self.table.template(mass, mixture.get('branching', {'bb': 1.}),
                                        rho_local_GeVcm3=rho, source_factor=factor)
        phi, phi_source = modulation(card, pbar_result)
        rate = mixture['total_sigmav_cm3s']
        cases = {}
        for name, value in [('all_dm', rate), ('thermal', rate*xsi**2)]:
            flux = lambda t: flux_per_nucleon(template, t, value, phi[name])
            knots = (self.table.energy-phi[name])/A
            windows = {key: window_forecast(flux, *window, points=knots)
                       for key, window in WINDOWS.items()}
            cases[name] = dict(sigmav_cm3s=value, phi_GV=phi[name], windows=windows,
                               AMS02_total_events=windows['AMS02_low']['expected_events']
                                                  +windows['AMS02_high']['expected_events'])
        grid = np.unique(np.r_[np.geomspace(.1, 100., 360),
                              [v for low, high, _ in WINDOWS.values() for v in (low, high)]])
        columns = [grid]
        for name in ('all_dm', 'thermal'):
            value = cases[name]['sigmav_cm3s']
            columns += [flux_per_nucleon(template, grid, value, 0.),
                        flux_per_nucleon(template, grid, value, phi[name])]
        secondary = cr_config.secondary_path(card, species='dbar', required=False)
        secondary_flux = None
        if secondary:
            from .pbar.core import Spectrum
            a = np.loadtxt(secondary)
            spectrum = Spectrum(a[:,0],a[:,1])  # native secondary is already per GeV/n
            u = grid+phi['all_dm']/2.
            secondary_flux = spectrum(u)*grid*(grid+MASS_GEV)/(u*(u+MASS_GEV))
        return dict(status='ok', mass_GeV=float(mass), propagation=cr_config.selection(card)['propagation'], halo=cr_config.selection(card)['halo'],
                    total_sigmav_cm3s=rate, branching=mixture.get('branching', {}), rates=records,
                    table=template.metadata, species='antideuteron', mass_number=2, absolute_charge=1,
                    mass_particle_GeV=MASS_GEV, table_convention=convention, source_factor=factor,
                    rho_local_GeVcm3=rho, xsi=xsi, solar_modulation_source=phi_source, cases=cases,
                    method='N=2*integral_T1^T2(Phi_TOA(T)dT)/((T2-T1)*sensitivity)',
                    flux_units='(GeV/n)^-1 m^-2 s^-1 sr^-1',
                    zero_model_rate=rate == 0.,
                    no_tabulated_signal=all(np.all(component.flux == 0.) for _, component in template.components),
                    notes=['DM signal only; no secondary Dbar background added.',
                           'Two-event sensitivity forecast, not a detector response or a confidence bound.',
                           'The table uses total kinetic energy; both energy and flux Jacobian are converted.',
                           'All-DM uses full density; thermal uses xsi^2 times the rate.',
                           'Native source convention is the same explicit assumption as for pbar.'],
                    secondary_path=str(secondary) if secondary else None,
                    _secondary_flux=secondary_flux,
                    _flux_columns=np.column_stack(columns))


def empty_observables():
    return dict(dict.fromkeys(NUMERIC_KEYS, -1.), AMS02_dbar_status='not_evaluated')


def observables(result):
    values = empty_observables()
    values['AMS02_dbar_status'] = result.get('status', 'not_evaluated')
    if result.get('status') != 'ok':
        return values
    for case, suffix in [('all_dm', ''), ('thermal', '_thermal')]:
        row = result['cases'][case]
        for window in ('low', 'high'):
            values['AMS02_dbar_events_'+window+suffix] = row['windows']['AMS02_'+window]['expected_events']
        values['AMS02_dbar_events_total'+suffix] = row['AMS02_total_events']
        values['GAPS_dbar_events'+suffix] = row['windows']['GAPS']['expected_events']
        values['AMS02_dbar_phi'+suffix] = row['phi_GV']
    values['AMS02_dbar_model_sigmav'] = result['total_sigmav_cm3s']
    return values


def write_outputs(out_dir, result, *, detailed=True):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    public = {key: value for key, value in result.items() if not key.startswith('_')}
    (out/'AMS02_dbar_events.json').write_text(json.dumps(public, indent=2, allow_nan=False)+'\n')
    for name in ('AMS02_dbar_events.dat', 'AMS02_dbar_fluxes.dat', 'AMS02_dbar_secondary.dat', 'AMS02_dbar_comparison.pdf'):
        (out/name).unlink(missing_ok=True)
    if not detailed or result.get('status') != 'ok':
        return
    header = '# case experiment T_min_GeV_per_n T_max_GeV_per_n sensitivity average_flux expected_events\n'
    lines = [header.rstrip()]
    for case in ('all_dm', 'thermal'):
        for name, row in result['cases'][case]['windows'].items():
            lines.append('%s %s %.12e %.12e %.12e %.12e %.12e' %
                         (case, name, row['T_min_GeV_per_n'], row['T_max_GeV_per_n'],
                          row['sensitivity_per_GeV_per_n'], row['average_flux_per_GeV_per_n'],
                          row['expected_events']))
    (out/'AMS02_dbar_events.dat').write_text('\n'.join(lines)+'\n')
    np.savetxt(out/'AMS02_dbar_fluxes.dat', result['_flux_columns'], fmt='%.12e',
               header='T_GeV_per_n model_IS_per_T model_TOA_per_T thermal_IS_per_T thermal_TOA_per_T\n'
                      'Flux units: (GeV/n)^-1 m^-2 s^-1 sr^-1; see AMS02_dbar_events.json for phi and rates')
    (out/'AMS02_dbar_secondary.dat').unlink(missing_ok=True)
    if result.get('_secondary_flux') is not None:
        np.savetxt(out/'AMS02_dbar_secondary.dat', np.c_[result['_flux_columns'][:,0],result['_secondary_flux']],
                   header='T_GeV_per_n secondary_TOA_per_T; reference secondary, not included in signal counts')


def save_plot(out_dir, result):
    if result.get('status')!='ok':return
    if "." in __package__:
        from ..plotting import signal_plots
    else:
        from plotting import signal_plots
    signal_plots.dbar_plot(result,Path(out_dir)/'AMS02_dbar_comparison.pdf')


def print_results(result, emit):
    emit('====== Antideuteron event forecasts: %s / %s' % (result.get('propagation','unknown'),result.get('halo','unknown')))
    if result.get('status') != 'ok':
        emit('Dbar events: NOT EVALUATED (%s)' % result.get('status', 'not_evaluated'))
        return
    emit('Two-event sensitivity estimate; DM signal only')
    emit('Solar modulation: %s' % result['solar_modulation_source'])
    for case, label in [('all_dm', 'All DM'), ('thermal', 'Thermal')]:
        if case == 'thermal' and result['xsi'] >= 1.: continue
        row = result['cases'][case]
        emit('%s: phi = %.6g GV' % (label, row['phi_GV']))
        for name, window in row['windows'].items():
            emit('%s %s [%.3g, %.3g] GeV/n: %.6e expected events' %
                 (name, label, window['T_min_GeV_per_n'], window['T_max_GeV_per_n'], window['expected_events']))
        emit('AMS-02 %s total (two windows): %.6e expected events' % (label, row['AMS02_total_events']))
    if result['no_tabulated_signal']:
        emit('No populated Dbar signal in the selected native table channels.')


def save_summary(out, result):
    out.write('\n# Antideuteron event forecasts: %s / %s; two-event sensitivities\n' % (result.get('propagation','unknown'),result.get('halo','unknown')))
    out.write('AMS02_dbar_status = %s\n' % result.get('status', 'not_evaluated'))
    for key, value in observables(result).items():
        if key != 'AMS02_dbar_status': out.write('%s = %.10e\n' % (key, value))
    if result.get('status') == 'ok':
        out.write('# Solar modulation: '+result['solar_modulation_source']+'\n')
    out.write('# Forecast diagnostics: AMS02_dbar_events.json; flux tables are saved with detailed output.\n')

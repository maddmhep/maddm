"""CosmiXsPPPC charged fluxes for the model point, independent of likelihoods."""
from dataclasses import dataclass
from pathlib import Path
import json
import math

import numpy as np

from . import cr_config, dbar_interface, pbar_interface
from .pbar.core import FluxTable, Spectrum, SIGMAV_REF, MP, force_field


@dataclass(frozen=True)
class Particle:
    mass: float
    energy_scale: int  # dK/dT, equal to A for a nucleus and 1 otherwise


PARTICLES = {'pbar': Particle(MP, 1), 'dbar': Particle(dbar_interface.MASS_GEV, 2),
             'positron': Particle(0.00051099895, 1)}
ACTIVE = ('pbar', 'dbar')
POSITRON_NOTICE = 'Positron results will come soon.'
DIAGNOSTICS = 'CosmiXsPPPC_flux_earth.json'


class EarthFluxes:
    """Both table inputs are IS. Secondaries are already per displayed energy."""

    def __init__(self):
        self.tables = {}

    def evaluate(self, species, mass, mixture, majorana_dirac_factor, xsi, card, pbar_result=None):
        if mixture['status'] not in ('ok', 'zero_rate'):
            return mixture
        particle = PARTICLES[species]
        path = cr_config.table_path(species, card)
        table = self.tables.get(species)
        if table is None or table.path != path:
            table = self.tables[species] = FluxTable(path)
        convention, factor, rho = pbar_interface.source_scaling(card, majorana_dirac_factor)
        rate = mixture['total_sigmav_cm3s']
        if not 0. < xsi <= 1.:
            raise ValueError('Require 0 < xsi <= 1')
        template = (table.template(mass, mixture['branching'], rho_local_GeVcm3=rho,
                                   source_factor=factor) if rate else None)
        secondary_path = cr_config.secondary_path(card, species=species)
        data = np.loadtxt(secondary_path, ndmin=2)
        if data.shape[1] != 2:
            raise ValueError('Secondary table must contain energy and differential flux')
        secondary = Spectrum(*data.T)
        if species == 'positron':
            positron_phi = float(cr_config.value(card, 'positron_phi', .7))
            if not math.isfinite(positron_phi) or positron_phi < 0:
                raise ValueError('Require finite nonnegative positron_phi [GV]')
            phi = dict(all_dm=positron_phi, thermal=positron_phi)
            phi_source = 'fixed positron_phi'
        else:
            phi, phi_source = dbar_interface.modulation(card, pbar_result)
        a = particle.energy_scale
        pivot = float(cr_config.value(card, 'pbar_pivot', 10.))
        if species == 'pbar' and (not math.isfinite(pivot) or pivot <= 0):
            raise ValueError('Require positive finite pbar_pivot')

        # Write only the common domain of the inputs, including the IS points
        # shifted by either Fisk potential. No flux is extrapolated.
        lower = secondary.energy[0]
        upper = secondary.energy[-1] - max(phi.values())/a
        if template is not None:
            lower = max(lower, table.energy[0]/a)
            upper = min(upper, (table.energy[-1]-max(phi.values()))/a)
        if not 0 < lower < upper:
            raise ValueError('No common kinetic-energy range after solar modulation')
        grid = np.geomspace(lower, upper, 500)
        if template is not None:
            grid = np.unique(np.r_[grid, table.energy[(table.energy > lower*a) & (table.energy < upper*a)]/a])
        columns = [grid]
        cases = {}
        for case, scale in [('all_dm', 1.), ('thermal', xsi*xsi)]:
            nuisance = ((pbar_result or {}).get('model_points', {}).get(case, {}) if species == 'pbar' else {})
            norm, slope = float(nuisance.get('norm', 1.)), float(nuisance.get('slope', 0.))

            def secondary_per_k(k):
                tilt = (k/pivot)**slope if slope else 1.
                return norm*secondary(k/a)*tilt/a

            sec_is = secondary_per_k(a*grid)*a
            sec_toa = a*force_field(secondary_per_k, a*grid, phi[case], mass_GeV=particle.mass)
            if template is None:
                dm_is = dm_toa = np.zeros_like(grid)
            else:
                dm_scale = a*rate*scale/SIGMAV_REF
                dm_is = dm_scale*template(a*grid)
                dm_toa = dm_scale*force_field(template, a*grid, phi[case], mass_GeV=particle.mass)
            columns.extend((sec_is, dm_is, sec_is+dm_is, sec_toa, dm_toa, sec_toa+dm_toa))
            cases[case] = dict(sigmav_cm3s=rate*scale, phi_GV=phi[case],
                               secondary_norm=norm, secondary_slope=slope)

        energy = 'T_GeV_per_n' if a > 1 else 'K_GeV'
        units = '(GeV/n)^-1' if a > 1 else 'GeV^-1'
        units += ' m^-2 s^-1 sr^-1'
        names = ['%s_%s_%s' % (case, frame, component) for case in cases for frame in ('IS', 'TOA')
                 for component in ('secondary', 'DM', 'total')]
        return dict(status='ok', **cr_config.selection(card), model_mass_GeV=float(mass),
                    flux_units=units, energy_column=energy, source_factor=factor,
                    table_convention=convention, rho_local_GeVcm3=rho, xsi=xsi,
                    dm_table=str(path), secondary_table=str(secondary_path), cases=cases,
                    branching=mixture.get('branching', {}), solar_modulation_source=phi_source,
                    no_tabulated_signal=template is None or
                    all(np.all(component.flux == 0.) for _, component in template.components),
                    _columns=np.column_stack(columns), _header=' '.join([energy]+names)+'\nFlux units: '+units)


def write_outputs(directory, results, *, detailed=True):
    directory = Path(directory)
    summary = dict(positron=dict(status='pending_tables', message=POSITRON_NOTICE))
    for species in ACTIVE:
        result = results[species]
        target = directory/('CosmiXsPPPC_%s_flux_earth.dat' % species)
        target.unlink(missing_ok=True)
        summary[species] = {k: v for k, v in result.items() if not k.startswith('_')}
        if detailed and result['status'] == 'ok':
            np.savetxt(target, result['_columns'], fmt='%.12e', header=result['_header'])
            summary[species]['file'] = target.name
    (directory/DIAGNOSTICS).write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    return summary

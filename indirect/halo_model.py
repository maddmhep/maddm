"""Shared Galactic velocity settings for recoil and annihilation rates.

The card selection (halo_velocity_profile, vMP, vescape and the profile parameters) defines both
the Galactic pair average of the annihilation rates (galactic_rule) and the RAPIDD direct detection
halo tables (rapidd.core.halo_kwargs_from_card), with the same velocity laws."""
from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
import math

import numpy as np
from scipy.special import erf

from .cr_config import value, local_density
from .velocity_profiles import PROFILE_NAMES, PROFILE_INPUTS, RadialLaw, pair_rule


C_KMS = 299792.458
RULE_FILE = 'mw_shm_quadrature.dat'
VELOCITY_OUTPUT = 'Halo_velocity.json'
_PROFILE_DEFAULTS = {name: default for name, default, _ in PROFILE_INPUTS}


@dataclass(frozen=True)
class VelocityProfile:
    name: str
    primary: RadialLaw
    stream: RadialLaw = None
    fraction: float = 0.
    bulk_speed: float = 0.
    cosb: float = 0.

    @property
    def is_shm(self):
        return self.name in ('SHM', 'SHM_numeric')

    def metadata(self):
        result = dict(halo_velocity_profile=self.name, vmp_kms=self.primary.scale,
                      vescape_kms=self.primary.cutoff)
        if self.name == 'SHM_beta':
            result['halo_beta'] = self.primary.parameter
        elif self.name == 'Lisanti':
            result['halo_lisanti_k'] = self.primary.parameter
        elif self.stream is not None:
            result.update(halo_lmc_fraction=self.fraction, halo_lmc_vb=self.bulk_speed,
                          halo_lmc_sigb=self.stream.scale, halo_lmc_vcut=self.stream.cutoff,
                          halo_lmc_cosb=self.cosb)
        return result


def parameters(card):
    v0 = float(value(card, 'vmp', 238.))
    vesc = float(value(card, 'vescape', 544.))
    if not (0. < v0 < .1*C_KMS and 0. < vesc < .1*C_KMS and vesc/v0 >= .05):
        raise ValueError('Require positive Galactic vMP and vescape below 0.1 c, with vescape/vMP >= 0.05')
    return v0, vesc


def velocity_profile(card):
    """Read the single card selection consumed by RAPIDD and Galactic ID (with mw_velocity_model = fixed
    the annihilation rates use fixed speeds and the selection enters only the RAPIDD recoil rates)."""
    names = {p.casefold(): p for p in PROFILE_NAMES}
    raw = str(value(card, 'halo_velocity_profile', 'SHM')).strip()
    if raw.casefold() not in names:
        raise ValueError('halo_velocity_profile must be one of '+', '.join(PROFILE_NAMES))
    name = names[raw.casefold()]
    v0, vesc = parameters(card)

    def get(key):
        x = float(value(card, key, _PROFILE_DEFAULTS[key]))
        if not math.isfinite(x):
            raise ValueError(key+' must be finite')
        return x

    if name == 'SHM_beta':
        beta = get('halo_beta')
        if not 0. <= beta <= 1.:
            raise ValueError('Require 0 <= halo_beta <= 1')
        return VelocityProfile(name, RadialLaw(v0, vesc, 'beta', beta))
    if name == 'Lisanti':
        k = get('halo_lisanti_k')
        if k < .1:
            raise ValueError('Require halo_lisanti_k >= 0.1 for the RAPIDD Lisanti law')
        return VelocityProfile(name, RadialLaw(v0, vesc, 'lisanti', k))
    primary = RadialLaw(v0, vesc)
    if name == 'SHM_wLMC':
        fraction, vb, width, cut, cosb = [get(k) for k in (
            'halo_lmc_fraction', 'halo_lmc_vb', 'halo_lmc_sigb', 'halo_lmc_vcut', 'halo_lmc_cosb')]
        if not (0. <= fraction <= 1. and -1. <= cosb <= 1.):
            raise ValueError('Require 0 <= halo_lmc_fraction <= 1 and -1 <= halo_lmc_cosb <= 1')
        if not (0. <= vb < .1*C_KMS and 0. < width < .1*C_KMS and
                0. < cut and vb+cut < .1*C_KMS and cut/width >= .05):
            raise ValueError('Require positive LMC width and cutoff, vcut/sigb >= 0.05, and vb+vcut < 0.1 c')
        return VelocityProfile(name, primary, RadialLaw(width, cut), fraction, vb, cosb)
    return VelocityProfile(name, primary)


@lru_cache(maxsize=16)
def relative_speed_rule(v0, vesc, n=16):
    """Normalized quadrature for u/c from two independent Galactic velocities.

    Both particles have f(v) proportional to exp(-|v|²/v0²) inside vesc.
    For q=u/v0, z=vesc/v0, integrate the pair's average velocity over the
    intersection of the two escape spheres, then normalize the speed PDF.
    The observer's boost cancels from their relative velocity.
    """
    parameters({'vmp': v0, 'vescape': vesc})
    if not 8 <= n <= 128:
        raise ValueError('Require 8 <= mw_shm_npts <= 128')
    nodes, weights = np.polynomial.legendre.leggauss(n)
    inner_nodes, inner_weights = np.polynomial.legendre.leggauss(40)
    z = vesc/v0
    # Beyond q=12 the omitted Gaussian probability is below 1e-29.
    qmax = min(2*z, 12.)
    q = (nodes+1)*qmax/2
    y_max = np.minimum(z-q/2, 8.)
    y = (inner_nodes[None, :]+1)*y_max[:, None]/2
    yp = y+q[:, None]/2
    transverse = -np.expm1(-2*(z-yp)*(z+yp))
    area = y_max/2*np.sum(inner_weights*np.exp(-2*y*y)*transverse, axis=1)
    weights = weights*q*q*np.exp(-q*q/2)*area
    weights /= weights.sum()
    return tuple((float(speed*v0/C_KMS), float(weight)) for speed, weight in zip(q, weights))


def rms_speed(rule):
    return math.sqrt(sum(w*u*u for u, w in rule))


def galactic_rule(profile, n=16):
    """Return normalized quadrature at relative speeds in units of c."""
    a = profile.primary
    if profile.name == 'SHM' or (profile.stream is not None and profile.fraction == 0.):
        return relative_speed_rule(a.scale, a.cutoff, n)
    if profile.fraction == 1.:
        # A common bulk velocity cancels from a pair of LMC particles.
        a = profile.stream
        return relative_speed_rule(a.scale, a.cutoff, n)
    return tuple((u/C_KMS, w) for u, w in
                 pair_rule(a, n, profile.stream, profile.fraction, profile.bulk_speed))


def _lab_motion_parameters(card):
    v0, _ = parameters(card)
    solar = np.array([value(card, k, default) for k, default in
                      [('vSun_r', 11.1), ('vSun_phi', 12.2), ('vSun_theta', 7.3)]], dtype=float)
    orbit = float(value(card, 'vEarth_mod', 29.8))
    if not np.all(np.isfinite(solar)) or not math.isfinite(orbit) or orbit < 0.:
        raise ValueError('Require finite Sun velocity components and nonnegative Earth orbital speed')
    return solar+np.array([0., v0, 0.]), orbit


def earth_velocity(card, days):
    """Laboratory velocity in Galactic coordinates, in km/s."""
    sun, orbit = _lab_motion_parameters(card)
    phase = np.atleast_1d(np.asarray(days, float))*0.0172
    return sun[:, None] + orbit*(
        np.array([.9941, .1088, .0042])[:, None]*np.cos(phase)
        + np.array([-.0504, .4946, -.8677])[:, None]*np.sin(phase))


def earth_speed(card, days):
    return np.linalg.norm(earth_velocity(card, days), axis=0)


def mean_inverse_speed(vmin, v0, vesc, boost):
    """Analytic lab-frame SHM integral in (km/s)^-1; includes b≥vesc."""
    x = np.asarray(vmin, float)/v0
    y, z = boost/v0, vesc/v0
    ez = math.exp(-z*z)
    norm = erf(z)-2*z*ez/math.sqrt(math.pi)
    if y < 1.e-8:
        return np.where(x < z, 2*(np.exp(-x*x)-ez)/(math.sqrt(math.pi)*norm*v0), 0.)
    prefactor = 1./(2*y*v0*norm)
    whole = erf(x+y)-erf(x-y)-4*y*ez/math.sqrt(math.pi)
    tail = erf(z)-erf(x-y)-2*(z+y-x)*ez/math.sqrt(math.pi)
    return np.maximum(0., np.where(x <= y-z, 1./boost,
           np.where(x <= z-y, prefactor*whole, np.where(x < z+y, prefactor*tail, 0.))))


def write_velocity_diagnostics(directory, card, profile, rule, *, write_grid=True):
    directory = Path(directory)
    data = profile.metadata()
    data.update(rhoDM_gev_cm3=local_density(card),
                mw_velocity_model=value(card, 'mw_velocity_model', 'shm'))
    if rule is not None:
        rule_file = 'MW_relative_speed.dat'
        if write_grid:
            np.savetxt(directory/rule_file, np.array(rule)*[C_KMS, 1.], fmt='%.16e',
                       header='relative_speed_km_s normalized_weight')
            data['mw_quadrature_file'] = rule_file
        data.update(mw_shm_npts=len(rule),
                    mw_rms_relative_speed_kms=rms_speed(rule)*C_KMS)
    (directory/VELOCITY_OUTPUT).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def write_rule(project, rule):
    path = Path(project)/'Cards'/RULE_FILE
    if rule is None:
        path.unlink(missing_ok=True)
    else:
        np.savetxt(path, rule, fmt='%.17e', header=str(len(rule)), comments='')

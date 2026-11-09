"""RAPIDD velocity laws and deterministic speed integrals in km/s.

Each radial law is isotropic in its own rest frame. Its PDF is used both
for the lab-frame recoil integral and for the Galactic pair average.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh_tridiagonal


PROFILE_NAMES = ('SHM', 'SHM_numeric', 'SHM_beta', 'SHM_wLMC', 'Lisanti')
# RAPIDD calc_halo inputs, with beta=1 making SHM_beta the smooth-cutoff choice.
PROFILE_INPUTS = (
    ('halo_beta', 1., 'SHM_beta: smooth-cutoff parameter in [0, 1]'),
    ('halo_lisanti_k', 1.5, 'Lisanti: dimensionless tail-shape parameter'),
    ('halo_lmc_fraction', .006, 'SHM_wLMC: local LMC density fraction in [0, 1]'),
    ('halo_lmc_vb', 570., 'SHM_wLMC: Galactic bulk speed [km/s]'),
    ('halo_lmc_sigb', 100., 'SHM_wLMC: Maxwellian scale about the bulk [km/s]'),
    ('halo_lmc_vcut', 200., 'SHM_wLMC: cutoff about the bulk [km/s]'),
    ('halo_lmc_cosb', -.71, 'SHM_wLMC: cosine with the RAPIDD reference lab direction'),
    ('halo_lmc_azimuth_deg', 0., 'SHM_wLMC: azimuth from the projected Galactic radial axis [deg]'),
)


@lru_cache(maxsize=12)
def _legendre(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return (x+1.)/2., w/2.


def _integrate_nodes(lo, hi, n, power=1):
    """Gauss nodes on intervals; a power map smooths a Lisanti endpoint."""
    x, w = _legendre(n)
    dx = np.asarray(hi-lo)[..., None]
    if power == 1:
        return np.asarray(lo)[..., None]+dx*x, dx*w
    return np.asarray(hi)[..., None]-dx*x**power, dx*(w*power*x**(power-1))


@dataclass(frozen=True)
class RadialLaw:
    scale: float
    cutoff: float
    shape: str = 'maxwell'
    parameter: float = 0.

    @property
    def bound(self):
        # At >12 widths the Maxwellian tail is negligible, also for the
        # subtracted beta and Lisanti laws. Keeps quadrature resolved for
        # an arbitrarily cold component. The physical cutoff is retained.
        return min(self.cutoff, 12*self.scale)

    @property
    def power(self):
        return 4 if self.shape == 'lisanti' else 1

    def _phase(self, x):
        """Unnormalized f_3D at x=v/scale, inside the escape sphere."""
        z = self.cutoff/self.scale
        d = np.maximum(0., (z-x)*(z+x))
        g = np.exp(-x*x)
        if self.shape == 'beta':
            return g*((1-self.parameter)-self.parameter*np.expm1(-d))
        if self.shape == 'lisanti':
            return g*(-np.expm1(-d/self.parameter))**self.parameter
        return g

    @lru_cache(maxsize=32)
    def radial_norm(self):
        """Integral of x² f_3D(x) dx; the angular factor cancels."""
        return quad(lambda x: x*x*self._phase(x), 0., self.bound/self.scale,
                    epsabs=1.e-16, epsrel=1.e-11)[0]

    def pdf_over_r(self, r):
        x = r/self.scale
        return x*self._phase(x)/(self.scale**2*self.radial_norm())

    def inverse_speed_interval(self, lo, hi):
        """Integral of p(r)/r over the intersection with the radial support."""
        hi = np.minimum(hi, self.bound)
        lo = np.minimum(lo, hi)
        if self.shape == 'maxwell':
            x, y = lo/self.scale, hi/self.scale
            return np.exp(-x*x)*(-np.expm1(-(y-x)*(y+x)))/(2*self.scale*self.radial_norm())
        r, w = _integrate_nodes(lo, hi, 36, self.power)
        return np.sum(w*self.pdf_over_r(r), axis=-1)

    def eta(self, vmin, boost):
        """Lab-frame mean inverse speed for a boost of given magnitude."""
        vmin = np.asarray(vmin, dtype=float)
        if boost < 1.e-9*self.scale:
            return self.inverse_speed_interval(vmin, self.bound)
        end = np.full_like(vmin, self.bound)
        # |r-boost| and max(vmin,|r-boost|) change branch at these radii.
        edges = np.sort(np.stack((np.zeros_like(vmin), end, np.minimum(boost, end),
                                  np.minimum(abs(vmin-boost), end), np.minimum(vmin+boost, end)), axis=-1), axis=-1)
        r, w = _integrate_nodes(edges[..., :-1], edges[..., 1:], 48, self.power)
        span = np.maximum(0., r+boost-np.maximum(vmin[..., None, None], abs(r-boost)))
        return np.sum(w*self.pdf_over_r(r)*span, axis=(-2, -1))/(2*boost)


def _centered_pair_kernel(speed, a, b):
    """p(speed)/speed for the difference of two independent isotropic draws."""
    lo, hi = np.clip(speed-b.bound, 0., a.bound), np.minimum(speed+b.bound, a.bound)
    edges = np.sort(np.stack((lo, hi, np.clip(speed, lo, hi),
                              np.clip(b.bound-speed, lo, hi)), axis=-1), axis=-1)
    r, w = _integrate_nodes(edges[..., :-1], edges[..., 1:], 48, a.power)
    v = speed[..., None, None]
    inner = b.inverse_speed_interval(abs(v-r), v+r)
    return .5*np.sum(w*a.pdf_over_r(r)*inner, axis=(-2, -1))


def _pair_measure(a, b, drift, order):
    """Resolve pair density on a positive reference quadrature."""
    breaks = np.unique([0., abs(a.bound-b.bound), a.bound, b.bound, a.bound+b.bound])
    if drift <= 1.e-9:
        u, du = _integrate_nodes(breaks[:-1], breaks[1:], order)
        u, du = u.ravel(), du.ravel()
        pdf = u*_centered_pair_kernel(u, a, b)
    else:
        # The centered pair vector is isotropic. Add a fixed drift with a
        # second angular integral and split at its kinematic branch points.
        edges = np.unique(np.r_[drift+breaks, abs(drift-breaks)])
        u, du = _integrate_nodes(edges[:-1], edges[1:], order)
        u, du = u.ravel(), du.ravel()
        lo = np.minimum(abs(u-drift), breaks[-1])
        hi = np.minimum(u+drift, breaks[-1])
        pdf = np.empty_like(u)
        for start in range(0, len(u), 24):
            batch = slice(start, start+24)
            bounds = np.clip(breaks, lo[batch, None], hi[batch, None])
            s, ds = _integrate_nodes(bounds[:, :-1], bounds[:, 1:], 32)
            pdf[batch] = u[batch]/(2*drift)*np.sum(ds*_centered_pair_kernel(s, a, b), axis=(1, 2))
    weights = du*pdf
    total = float(weights.sum())
    if not np.isclose(total, 1., atol=3.e-5, rtol=0.):
        raise ValueError('Could not normalize the selected relative-speed distribution (integral %.8g)' % total)
    return u, weights/total


def _gauss_from_measure(speeds, weights, n):
    """Compress a resolved measure to n weighted speeds via its Jacobi matrix."""
    keep = weights > 0
    u, w = speeds[keep], weights[keep]
    scale = u.max()/2
    x = u/scale-1
    basis = np.zeros((n, len(u)))
    basis[0] = np.sqrt(w/w.sum())
    diagonal, off = np.empty(n), np.empty(n-1)
    for j in range(n):
        r = x*basis[j]
        diagonal[j] = np.dot(basis[j], r)
        if j < n-1:
            # Reorthogonalize to maintain accuracy at large mw_shm_npts.
            vectors = basis[:j+1]
            for _ in range(2):
                r -= np.einsum('i,ij->j', np.einsum('ij,j->i', vectors, r), vectors)
            off[j] = np.linalg.norm(r)
            if off[j] <= np.finfo(float).eps:
                raise ValueError('Relative-speed quadrature lost rank; reduce mw_shm_npts')
            basis[j+1] = r/off[j]
    values, vectors = eigh_tridiagonal(diagonal, off, lapack_driver='stev')
    w = vectors[0]**2
    return (values+1)*scale, w/w.sum()


@lru_cache(maxsize=24)
def pair_rule(primary, n, stream=None, fraction=0., bulk_speed=0.):
    """n speed/weight pairs in km/s for the selected Galactic distribution."""
    if not 8 <= n <= 128:
        raise ValueError('Require 8 <= mw_shm_npts <= 128')
    components = [(primary, primary, 0., (1-fraction)**2)]
    if stream is not None and fraction:
        components.extend(((primary, stream, bulk_speed, 2*fraction*(1-fraction)),
                           (stream, stream, 0., fraction*fraction)))
    parts = []
    for a, b, drift, weight in components:
        if weight:
            u, w = _pair_measure(a, b, drift, max(64, n))
            parts.append((u, weight*w))
    u, w = _gauss_from_measure(np.concatenate([a for a, _ in parts]),
                                np.concatenate([b for _, b in parts]), n)
    return tuple(zip(map(float, u), map(float, w)))

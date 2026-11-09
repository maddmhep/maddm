"""Bounded one-dimensional profile solver for the existing dwarf likelihood."""
import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def sum_bin_loglikes(grid, loglikes, flux):
    """Vectorized equivalent of the existing per-bin interp1d likelihoods."""
    grid, loglikes, flux = np.asarray(grid), np.asarray(loglikes), np.asarray(flux)
    index = np.sum(grid <= flux[:, None], axis=1)-1
    index = np.clip(index, 0, grid.shape[1]-2)
    row = np.arange(len(flux))
    lo, hi = grid[row, index], grid[row, index+1]
    weight = (flux-lo)/(hi-lo)
    values = (1-weight)*loglikes[row, index]+weight*loglikes[row, index+1]
    values = np.where((flux < grid[:, 0]) | (flux > grid[:, -1]), -1.e5, values)
    return float(np.sum(values))


def upper_limit(loglike, ll_null, minimum, maximum, q_threshold, step_size_scaling=1.):
    """Find the upper crossing above the physical best fit, including rate=0.

    Work in log10(rate), keep all optimization inside the requested interval,
    and return the largest detected upper crossing. No p-value tolerance band
    or unconstrained simplex refinement is used.
    """
    if not (0 < minimum < maximum and q_threshold > 0 and step_size_scaling > 0):
        raise ValueError('Invalid Fermi rate range, threshold, or scan resolution')
    low, high = math.log10(minimum), math.log10(maximum)
    cache = {}

    def at(x):
        x = float(x)
        if x not in cache:
            value = float(loglike(10.**x))
            if not math.isfinite(value):
                raise ValueError('Non-finite combined dwarf likelihood')
            cache[x] = value
        return cache[x]

    grid = np.linspace(low, high, max(41, int(math.ceil(4*(high-low)*step_size_scaling))+1))
    values = np.array([at(x) for x in grid])
    candidates = [(float(ll_null), None)] + list(zip(values.tolist(), grid.tolist()))
    for i in range(1, len(grid)-1):
        if values[i] > values[i-1]+1.e-10 and values[i] >= values[i+1]:
            fit = minimize_scalar(lambda x: -at(x), bounds=(grid[i-1], grid[i+1]),
                                  method='bounded', options={'xatol': 1.e-7})
            if not fit.success:
                raise RuntimeError('Fermi best-fit refinement failed')
            candidates.append((at(fit.x), float(fit.x)))
    ll_best, x_best = max(candidates, key=lambda p: p[0])
    if x_best == high:
        raise RuntimeError('Fermi best fit reaches sigmavmax; enlarge the search interval')
    search = np.unique(np.r_[grid if x_best is None else grid[grid > x_best],
                             low if x_best is None else x_best])

    def q(x):
        return 2*(ll_best-at(x))-q_threshold

    q_values = [q(x) for x in search]
    brackets = [(a, b) for a, b, qa, qb in zip(search[:-1], search[1:], q_values[:-1], q_values[1:])
                if qa <= 0 < qb]
    if not brackets or q_values[-1] <= 0:
        raise RuntimeError('Combined Fermi upper crossing not bracketed; adjust sigmavmin/sigmavmax')
    root = brentq(q, *brackets[-1], xtol=1.e-7, rtol=1.e-12)
    residual = q(root)
    if abs(residual) > 2.e-3:
        raise RuntimeError('Combined Fermi upper crossing failed its likelihood tolerance')
    if max(cache.values()) > ll_best+1.e-4:
        raise RuntimeError('Dwarf scan found a better minimum during crossing refinement; increase scan resolution')
    return dict(status='ok', sigmav95_cm3s=10.**root,
                sigmav_best_cm3s=0. if x_best is None else 10.**x_best,
                loglike_best=ll_best, loglike_null=float(ll_null),
                loglike_upper=at(root), q_threshold=float(q_threshold),
                q_at_limit=float(2*(ll_best-at(root))), crossing_residual=float(residual),
                crossing_count=len(brackets), rate_range_cm3s=[minimum, maximum],
                profile=[dict(sigmav_cm3s=10.**x, loglike=ll) for x, ll in sorted(cache.items())])

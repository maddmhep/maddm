"""AMS-02 antiproton likelihood; fixed QUAINT-MED/NFW annihilation templates.

Fluxes entering the likelihood are per GV, per m^2, per s, per sr.
The only runtime dependencies are NumPy and SciPy. No import-time fits or I/O.
"""
from dataclasses import dataclass, asdict
from itertools import product
from pathlib import Path
import hashlib
import math

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize, minimize_scalar, brentq

MP = 0.9382720813  # GeV; |Z|=A=1 for antiprotons
SIGMAV_REF = 3.e-26  # cm^3 s^-1, native table normalization
DATA_DIR = Path(__file__).resolve().parents[1] / 'ExpData' / 'AMS02_pbar'
CHANNELS = tuple('ee mumu tautau nu_e_nu_e nu_mu_nu_mu nu_tau_nu_tau uu dd ss cc bb tt AA gg WW ZZ hh AZ hZ'.split())


def kinetic_energy(rigidity):
    """TOA kinetic energy [GeV] from positive absolute rigidity [GV]."""
    r = np.asarray(rigidity, dtype=float)
    if np.any(~np.isfinite(r)) or np.any(r <= 0):
        raise ValueError('Absolute rigidity must be finite and positive')
    return r*r / (np.hypot(r, MP) + MP)


def dtdr(rigidity):
    r = np.asarray(rigidity, dtype=float)
    return r / np.hypot(r, MP)


def force_field(flux_is, t_toa, phi_GV, *, mass_GeV=MP, charge=1):
    """TOA flux per total kinetic energy, with shift |Z|*phi and p_TOA²/p_IS²."""
    t = np.asarray(t_toa, dtype=float)
    if not np.isfinite(phi_GV) or phi_GV < 0 or np.any(t <= 0):
        raise ValueError('Require positive kinetic energy and nonnegative finite phi')
    u = t + abs(charge)*phi_GV
    return flux_is(u) * t*(t+2*mass_GeV) / (u*(u+2*mass_GeV))


class Spectrum:
    """Shape-preserving log-log interpolation, no extrapolation or smoothing.

    Zero-to-positive intervals use flux-linear interpolation in log energy;
    all-zero intervals stay exactly zero. This avoids logarithms of zero.
    """
    def __init__(self, energy_GeV, flux_per_GeV):
        e, f = np.asarray(energy_GeV, float), np.asarray(flux_per_GeV, float)
        if (e.ndim != 1 or f.shape != e.shape or e.size < 2
                or not np.all(np.isfinite(e)) or not np.all(np.isfinite(f))
                or np.any(e <= 0) or np.any(np.diff(e) <= 0) or np.any(f < 0)):
            raise ValueError('Need ordered positive energies and finite nonnegative flux')
        self.energy = e.copy()
        self.flux = f.copy()
        self.x = np.log(e)
        self._ip = PchipInterpolator(self.x, np.log(np.maximum(f, 1.e-300)), extrapolate=False)
        self._dip = self._ip.derivative()

    def derivative(self, energy_GeV):
        """Derivative per GeV of the same interpolant used by __call__."""
        e = np.asarray(energy_GeV, float)
        value = self(e)  # validates the interpolation domain
        x = np.log(np.clip(e, self.energy[0], self.energy[-1]))
        j = np.clip(np.searchsorted(self.x, x, side='right')-1, 0, self.x.size-2)
        zero = (self.flux[j] == 0) | (self.flux[j+1] == 0)
        linear = (self.flux[j+1]-self.flux[j]) / ((self.x[j+1]-self.x[j])*e)
        return np.where(zero, linear, value*self._dip(x)/e)

    def minimum(self, energy_low, energy_high):
        """Rigorous interval lower bound: PCHIP is monotone on each segment."""
        lo, hi = np.broadcast_arrays(energy_low, energy_high)
        self(lo)
        self(hi)
        if np.any(hi < lo):
            raise ValueError('Reversed energy interval')
        out = np.minimum(self(lo), self(hi))
        for e, f in zip(self.energy, self.flux):
            out = np.where((lo <= e) & (e <= hi), np.minimum(out, f), out)
        return out

    def __call__(self, energy_GeV):
        e = np.asarray(energy_GeV, float)
        if (np.any(~np.isfinite(e)) or np.any(e < self.energy[0]*(1-1e-13))
                or np.any(e > self.energy[-1]*(1+1e-13))):
            raise ValueError(f'Energy outside spectrum [{self.energy[0]:g}, {self.energy[-1]:g}] GeV')
        x = np.log(np.clip(e, self.energy[0], self.energy[-1]))
        j = np.clip(np.searchsorted(self.x, x, side='right')-1, 0, self.x.size-2)
        y = np.exp(self._ip(x))
        zero = (self.flux[j] == 0) | (self.flux[j+1] == 0)
        w = (x-self.x[j]) / (self.x[j+1]-self.x[j])
        linear = (1-w)*self.flux[j] + w*self.flux[j+1]
        return np.where(zero, linear, y)


class DMTemplate:
    """Sum of separately interpolated channel fluxes at SIGMAV_REF."""
    def __init__(self, components, metadata=None):
        self.components = tuple(components)
        self.metadata = dict(metadata or {})

    def __call__(self, energy_GeV):
        return sum(w*f(energy_GeV) for w, f in self.components)

    def derivative(self, energy_GeV):
        return sum(w*f.derivative(energy_GeV) for w, f in self.components)

    @classmethod
    def from_propagated_flux(cls, energy_GeV, flux_per_GeV, *, sigmav_ref_cm3s=SIGMAV_REF):
        """Already propagated local IS pbar flux; NOT an injection dN/dE.

        Input units: GeV^-1 m^-2 s^-1 sr^-1. The source convention/density
        belongs to the caller; no further 1/(4*pi), 1/mass^2 or 1/2 is applied.
        """
        if not np.isfinite(sigmav_ref_cm3s) or sigmav_ref_cm3s <= 0:
            raise ValueError('Reference annihilation rate must be finite and positive')
        spectrum = Spectrum(energy_GeV, np.asarray(flux_per_GeV)*SIGMAV_REF/sigmav_ref_cm3s)
        return cls([(1., spectrum)], {'source': 'user-provided propagated IS flux',
                    'input_sigmav_ref_cm3s': float(sigmav_ref_cm3s)})


class FluxTable:
    def __init__(self, path=DATA_DIR / 'pbar_flux_QUAINT_MED_NFW.dat'):
        self.path = Path(path)
        fields = self.path.stem.split('_')
        self.propagation = '-'.join(fields[-3:-1]) if len(fields) >= 5 else 'QUAINT-MED'
        self.halo = fields[-1] if len(fields) >= 5 else 'NFW'
        with self.path.open() as stream:
            header = stream.readline().split()[3:]
        if tuple(header) != CHANNELS:
            raise ValueError('Unexpected DM table channel order')
        a = np.loadtxt(self.path)
        self.masses = np.unique(a[:, 0])
        energy = np.unique(a[:, 1])
        if a.shape != (len(self.masses)*len(energy), len(CHANNELS)+2):
            raise ValueError('Unexpected DM table shape')
        cube = a.reshape(len(self.masses), len(energy), len(CHANNELS)+2)
        if not np.all(cube[:, :, 0] == self.masses[:, None]) or not np.all(cube[:, :, 1] == energy):
            raise ValueError('DM table must contain a common ordered energy grid')
        self.energy = 10.**energy
        self.flux = cube[:, :, 2:].copy()
        if np.any(~np.isfinite(self.flux)) or np.any(self.flux < 0):
            raise ValueError('Invalid DM flux')
        # Upstream 1e-100 marks unpopulated channels; do not turn it into a bound.
        self.flux[self.flux <= 1.e-90] = 0.

    def template(self, mass_GeV, branching=None, *, rho_local_GeVcm3=0.4, source_factor=1.):
        """Fixed spectral mixture. Missing channels have zero branching fraction.

        Linear interpolation in mass of the tabulated flux at fixed energies,
        matching the supplied upstream script. No mass extrapolation.
        source_factor is a multiplier relative to the native table convention.
        """
        m = float(mass_GeV)
        if not np.isfinite(m) or not self.masses[0] <= m <= self.masses[-1]:
            raise ValueError(f'Mass outside [{self.masses[0]:g}, {self.masses[-1]:g}] GeV')
        if not np.isfinite(rho_local_GeVcm3) or rho_local_GeVcm3 <= 0:
            raise ValueError('Local density must be finite and positive')
        if not np.isfinite(source_factor) or source_factor <= 0:
            raise ValueError('Source factor must be finite and positive')
        br = {'bb': 1.} if branching is None else dict(branching)
        if not br or set(br)-set(CHANNELS):
            raise ValueError(f'Branching keys must belong to {CHANNELS}')
        if any(not np.isfinite(v) or v < 0 for v in br.values()) or not np.isclose(sum(br.values()), 1., atol=1.e-8, rtol=0):
            raise ValueError('Branching fractions must be finite, nonnegative, and sum to one')
        j = int(np.clip(np.searchsorted(self.masses, m, side='right')-1, 0, len(self.masses)-2))
        w = (m-self.masses[j])/(self.masses[j+1]-self.masses[j])
        flux = (1-w)*self.flux[j] + w*self.flux[j+1]
        scale = (rho_local_GeVcm3/0.4)**2 * source_factor
        components = [(float(v)*scale, Spectrum(self.energy, flux[:, CHANNELS.index(ch)]))
                      for ch, v in br.items() if v > 0]
        return DMTemplate(components, {'mass_GeV': m, 'branching': br,
                          'propagation': self.propagation, 'halo': self.halo,
                          'rho_local_GeVcm3': rho_local_GeVcm3, 'source_factor': source_factor,
                          'normalization': 'native CosmiXsPPPC annihilation tables',
                          'table_sha256': hashlib.sha256(self.path.read_bytes()).hexdigest()})


@dataclass(frozen=True)
class FitConfig:
    norm_bounds: tuple = (0.8, 1.2)
    slope_bounds: tuple = (-0.1, 0.1)
    phi_bounds_GV: tuple = (0., 2.)
    pivot_GeV: float = 10.
    rigidity_min_GV: float = 0.
    rigidity_max_GV: float = math.inf
    delta_chi2: float = 2.705543454095404  # one-sided asymptotic 95%

    def __post_init__(self):
        for b in (self.norm_bounds, self.slope_bounds, self.phi_bounds_GV):
            if len(b) != 2 or not np.all(np.isfinite(b)) or b[0] >= b[1]:
                raise ValueError('Each nuisance bound must be finite and increasing')
        if self.norm_bounds[0] <= 0 or self.phi_bounds_GV[0] < 0:
            raise ValueError('Require positive normalization and nonnegative phi')
        if (not np.isfinite(self.pivot_GeV) or self.pivot_GeV <= 0
                or not np.isfinite(self.delta_chi2) or self.delta_chi2 <= 0
                or not np.isfinite(self.rigidity_min_GV) or self.rigidity_min_GV < 0
                or np.isnan(self.rigidity_max_GV)
                or self.rigidity_max_GV <= self.rigidity_min_GV):
            raise ValueError('Invalid pivot, fit range, or confidence threshold')


@dataclass
class ProfilePoint:
    sigmav_cm3s: float
    chi2: float
    norm: float
    slope: float
    phi_GV: float
    nuisance_at_bound: list

    def to_dict(self):
        return asdict(self)


class Likelihood:
    """Gaussian likelihood for the supplied rigidity flux measurements.

    The file contains symmetric supplied errors and zero horizontal errors.
    Compare at the supplied abscissae; do not invent bin edges or covariance.
    """
    def __init__(self, config=None, *, data_path=DATA_DIR/'pbar_ams02.txt',
                 secondary_path=DATA_DIR/'secondary_QUAINT_IS.txt',
                 observed_flux=None, errors=None):
        self.config = config or FitConfig()
        self.data_path, self.secondary_path = Path(data_path), Path(secondary_path)
        a = np.loadtxt(self.data_path)
        if (a.ndim != 2 or a.shape[1] != 6 or np.any(~np.isfinite(a))
                or np.any(a[:, 0] <= 0) or np.any(np.diff(a[:, 0]) <= 0)
                or np.any(a[:, 1:3] != 0) or np.any(a[:, 4:] <= 0)
                or not np.allclose(a[:, 4], a[:, 5], rtol=1.e-10, atol=0)):
            raise ValueError('Expected ordered rigidity data, zero x errors, symmetric positive errors')
        mask = ((a[:, 0] >= self.config.rigidity_min_GV)
                & (a[:, 0] <= self.config.rigidity_max_GV))
        if np.count_nonzero(mask) < 5:
            raise ValueError('Too few data points in selected fit range')
        self.rigidity = a[mask, 0]
        self.observed = a[mask, 3] if observed_flux is None else np.asarray(observed_flux, float)
        self.errors = a[mask, 4] if errors is None else np.asarray(errors, float)
        if (self.observed.shape != self.rigidity.shape or self.errors.shape != self.rigidity.shape
                or np.any(~np.isfinite(self.observed)) or np.any(~np.isfinite(self.errors))
                or np.any(self.errors <= 0)):
            raise ValueError('Observation/error override must match selected data and be finite')
        self.t = kinetic_energy(self.rigidity)
        self.jacobian = dtdr(self.rigidity)
        s = np.loadtxt(self.secondary_path)
        if s.ndim != 2 or s.shape[1] != 2:
            raise ValueError('Secondary input must have exactly two columns: IS kinetic energy [GeV], flux [GeV^-1 m^-2 s^-1 sr^-1]')
        self.secondary = Spectrum(s[:, 0], s[:, 1])
        # Check the full modulation domain before any fit. Never extrapolate.
        self.secondary(self.t + self.config.phi_bounds_GV[0])
        self.secondary(self.t + self.config.phi_bounds_GV[1])

    def components(self, template, sigmav_cm3s, norm, slope, phi_GV):
        c = self.config
        t_is = self.t + phi_GV
        # Deform the IS secondary, then modulate; evaluate tilt at T_IS.
        sec = norm * force_field(self.secondary, self.t, phi_GV) * (t_is/c.pivot_GeV)**slope
        dm = np.zeros_like(sec) if sigmav_cm3s == 0 else sigmav_cm3s/SIGMAV_REF * force_field(template, self.t, phi_GV)
        return sec*self.jacobian, dm*self.jacobian

    def profile(self, template, sigmav_cm3s, *, extra_starts=()):
        """Re-fit the nuisance parameters at one fixed physical signal rate.

        A is eliminated exactly by bounded weighted linear least squares.
        Only delta and phi require numerical minimization. All fits use
        deterministic multiple starts; no tangent approximation is used.
        """
        if not np.isfinite(sigmav_cm3s) or sigmav_cm3s < 0:
            raise ValueError('Annihilation rate must be finite and nonnegative')
        c = self.config
        bounds = [c.slope_bounds, c.phi_bounds_GV]
        for phi in c.phi_bounds_GV:
            if sigmav_cm3s:
                values = template(self.t + phi)
                if np.any(~np.isfinite(values)) or np.any(values < 0):
                    raise ValueError('Template must return finite nonnegative local IS flux')
        y = self.observed/self.errors

        def evaluate(x, details=False, gradient=False):
            slope, phi = x
            sec, dm = self.components(template, sigmav_cm3s, 1., slope, phi)
            b, d = sec/self.errors, dm/self.errors
            norm = float(np.clip(np.dot(b, y-d)/np.dot(b, b), *c.norm_bounds))
            residual = y - norm*b - d
            chi2 = float(np.dot(residual, residual))
            if gradient:
                shifted = self.t+phi
                pre = self.t*(self.t+2*MP)/(shifted*(shifted+2*MP))
                loss = 1/shifted+1/(shifted+2*MP)
                tilt = (shifted/c.pivot_GeV)**slope
                sec_prime = (self.secondary.derivative(shifted)*pre*tilt*self.jacobian
                             + sec*(slope/shifted-loss))
                dm_prime = (np.zeros_like(sec) if sigmav_cm3s == 0 else
                            sigmav_cm3s/SIGMAV_REF*pre*template.derivative(shifted)*self.jacobian-dm*loss)
                grad = -2*np.array([np.dot(residual, norm*b*np.log(shifted/c.pivot_GeV)),
                                   np.dot(residual, (norm*sec_prime+dm_prime)/self.errors)])
                return chi2, grad
            return (chi2, norm) if details else chi2

        # Interior and boundary starts resolve secondary/phi degeneracy.
        starts = list(product([c.slope_bounds[0], 0., c.slope_bounds[1]],
                              [c.phi_bounds_GV[0], np.mean(c.phi_bounds_GV), c.phi_bounds_GV[1]]))
        starts += list(extra_starts)
        fits = []
        for start in starts:
            x = np.clip(start, [b[0] for b in bounds], [b[1] for b in bounds])
            r = minimize(lambda x: evaluate(x, gradient=True), x, jac=True, method='L-BFGS-B', bounds=bounds,
                         options={'ftol': 2.e-13, 'gtol': 2.e-6, 'maxiter': 350, 'maxls': 40})
            if r.success and np.isfinite(r.fun):
                fits.append(r)
        if not fits:
            raise RuntimeError('Nuisance minimization failed at every start')
        best = min(fits, key=lambda r: r.fun)
        # Independent direction-set polish; retain only converged results.
        check = minimize(evaluate, best.x, method='Powell', bounds=bounds,
                         options={'ftol': 1.e-11, 'xtol': 1.e-9, 'maxiter': 250})
        if check.success and check.fun < best.fun:
            best = check
        value, norm = evaluate(best.x, True)
        slope, phi = map(float, best.x)
        at_bound = [name for name, v, b in zip(['norm', 'slope', 'phi_GV'],
                    [norm, slope, phi], [c.norm_bounds, c.slope_bounds, c.phi_bounds_GV])
                    if min(abs(v-b[0]), abs(v-b[1])) < 1.e-5*(b[1]-b[0])]
        return ProfilePoint(float(sigmav_cm3s), value, norm, slope, phi, at_bound)

    def upper_limit(self, template):
        """Physical mu>=0 best fit and asymptotic one-sided profile upper limit.

        Adaptive rate scale, a coarse global scan, then local refinements;
        find the uppermost crossing detected on a finite scan up to a rigorous
        exclusion ceiling. Arbitrarily narrow disconnected regions are not
        guaranteed to be resolved by the scan.
        """
        c = self.config
        background = self.profile(template, 0.)
        # Adapt signal scaling, avoiding optimizers acting directly on 1e-26.
        phis = np.linspace(*c.phi_bounds_GV, 65)
        dm_grid = np.array([force_field(template, self.t, p)*self.jacobian for p in phis])
        if np.any(~np.isfinite(dm_grid)) or np.any(dm_grid < 0):
            raise ValueError('Template returned invalid fluxes')
        if np.max(dm_grid) == 0:
            return self._result(template, background, background, None, [],
                                status='no_signal_in_table', crossing_count=0)
        norm_ref = np.linalg.norm(dm_grid[len(phis)//2]/self.errors)
        if norm_ref == 0:
            norm_ref = np.max(np.linalg.norm(dm_grid/self.errors, axis=1))
        alpha_scale = 1./norm_ref
        cache = {0.: background}

        def point(x):
            x = float(x)
            if x not in cache:
                cache[x] = self.profile(template, SIGMAV_REF*alpha_scale*x,
                            extra_starts=[(background.slope, background.phi_GV)])
            return cache[x]

        # Bound the DM flux from below over the entire continuous phi interval.
        # Every Spectrum is shape preserving between knots. Minima of the
        # separate positive components give a lower envelope for the mixture.
        # The force-field prefactor decreases with phi. Once this lower
        # envelope alone overshoots one datum by sqrt(chi0+q) sigma, all
        # larger rates are excluded regardless of nuisance parameters.
        lower_is = sum(w*f.minimum(self.t+c.phi_bounds_GV[0],
                        self.t+c.phi_bounds_GV[1]) for w, f in template.components)
        shifted = self.t+c.phi_bounds_GV[1]
        lower_dm = lower_is*self.t*(self.t+2*MP)/(shifted*(shifted+2*MP))*self.jacobian
        guaranteed_ratio = np.max(lower_dm/self.errors)
        if guaranteed_ratio <= 0:
            raise RuntimeError('No positive continuous-phi signal lower bound; narrow the phi interval')
        ceiling = max(16., 2.*(np.max(abs(self.observed/self.errors))
                          + np.sqrt(background.chi2+c.delta_chi2)+1.) / (alpha_scale*guaranteed_ratio))
        grid = np.r_[0., np.geomspace(1.e-3, ceiling, 48)]
        values = np.array([point(x).chi2 for x in grid])
        candidates = [background]
        for i in range(1, len(grid)-1):
            if values[i] <= values[i-1] and values[i] <= values[i+1]:
                r = minimize_scalar(lambda x: point(x).chi2,
                        bounds=(grid[i-1], grid[i+1]), method='bounded',
                        options={'xatol': 2.e-7, 'maxiter': 120})
                if not r.success:
                    raise RuntimeError('Signal best-fit minimization failed')
                candidates.append(point(r.x))
        best = min(candidates, key=lambda p: p.chi2)
        xbest = best.sigmav_cm3s/(SIGMAV_REF*alpha_scale)
        threshold = best.chi2+c.delta_chi2
        # Refine the high-rate scan; this checks for extra crossing intervals.
        high_grid = np.unique(np.r_[grid[grid > xbest],
                              np.geomspace(max(xbest, 1.e-3), ceiling, 40), xbest])
        q = np.array([point(x).chi2-threshold for x in high_grid])
        if q[-1] <= 0:
            raise RuntimeError('Failed to bracket the upper limit')
        brackets = [(a, b) for a, b, qa, qb in zip(high_grid[:-1], high_grid[1:], q[:-1], q[1:])
                    if qa <= 0 < qb]
        if not brackets:
            raise RuntimeError('No upper-limit crossing found')
        roots = [brentq(lambda x: point(x).chi2-threshold, a, b,
                        xtol=2.e-7, rtol=2.e-8) for a, b in brackets]
        upper = point(max(roots))
        if abs(upper.chi2-threshold) > 2.e-4:
            raise RuntimeError('Upper-limit crossing did not meet chi2 tolerance')
        if min(p.chi2 for p in cache.values()) < best.chi2-1.e-4:
            raise RuntimeError('Profile scan found a better minimum; increase signal scan resolution')
        return self._result(template, background, best, upper,
                    sorted(cache.values(), key=lambda p: p.sigmav_cm3s),
                    status='ok', crossing_count=len(roots))

    def _result(self, template, background, best, upper, profile, **kwargs):
        config = asdict(self.config)
        if math.isinf(config['rigidity_max_GV']):
            config['rigidity_max_GV'] = None
        warnings = [
            'Diagonal supplied flux errors; no experimental covariance was supplied.',
            'Point evaluation at supplied rigidity coordinates; original bin edges are absent.',
            'One effective force-field potential for the time-averaged AMS exposure.',
            'The scheme-specific secondary reference is shared across MIN/MED/MAX by the adopted analysis prescription.',
            'Native DM table source convention retained; verify any MadDM particle/antiparticle factor.',
            'Asymptotic profile threshold, not Monte Carlo calibrated coverage.',
        ]
        if best.nuisance_at_bound or (upper and upper.nuisance_at_bound):
            warnings.append('One or more nuisance parameters reach hard bounds; inspect fit records.')
        if kwargs.get('crossing_count', 0) > 1:
            warnings.append('Disconnected profile intervals found; reported endpoint is the largest upper crossing.')
        return dict(**kwargs, template=template.metadata, config=config,
                    n_data=len(self.rigidity), rigidity_range_GV=[float(self.rigidity[0]), float(self.rigidity[-1])],
                    background=background.to_dict(), best_fit=best.to_dict(),
                    upper_limit=upper.to_dict() if upper else None,
                    sigmav95_cm3s=upper.sigmav_cm3s if upper else None,
                    delta_chi2_at_limit=upper.chi2-best.chi2 if upper else None,
                    improvement_over_background=background.chi2-best.chi2,
                    profile=[p.to_dict() for p in profile], warnings=warnings,
                    input_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in [self.data_path, self.secondary_path]},
                    versions={'numpy': np.__version__, 'scipy': __import__('scipy').__version__})

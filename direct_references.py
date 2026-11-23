"""Neutrino-fog references for the DM-nucleon cross sections (xenon).

The tables of C. A. J. O'Hare (github.com/cajohare/NeutrinoFog, arXiv:2109.03116) are downloaded
by `generate`/`add direct_detection` (or `install neutrino_fog`), never by `launch`.
They are interpolated in log mass/log cross section, without extrapolation, and rescaled by
FOG_RHO/rhoDM like the experimental limits: the fog is the cross section at which the DM signal
reaches a given rate in a xenon detector, and the DM rate scales as rhoDM * sigma
(standard halo model, rho_0 = 0.3 GeV/cm^3 in WIMPFuncs.dRdE of NeutrinoFog).
The fog is a reference boundary, not an exclusion limit.
"""
import bisect
import functools
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FOG_SOURCE = ('https://github.com/cajohare/NeutrinoFog; C. A. J. O\'Hare, '
              'Phys. Rev. Lett. 127, 251802 (2021), arXiv:2109.03116')
FOG_KEYS = ('neutrino_fog_SI', 'neutrino_fog_SD_p', 'neutrino_fog_SD_n')
FOG_FILES = dict(zip(FOG_KEYS, ('Xe_SI.txt', 'Xe_SDp.txt', 'Xe_SDn.txt')))
FOG_RHO = 0.3  # local DM density [GeV/cm^3] of the published tables


def fog_path(key, root=ROOT):
    return Path(root)/'NeutrinoFog'/'data'/'floors'/FOG_FILES[key]


@functools.lru_cache(maxsize=32)
def _table(path, mtime_ns):
    rows = []
    for line in Path(path).read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        row = tuple(map(float, line.split()))
        if len(row) != 2 or any(not math.isfinite(x) or x <= 0 for x in row):
            raise ValueError('Invalid reference table: %s' % path)
        rows.append(row)
    if len(rows) < 2 or any(b[0] <= a[0] for a, b in zip(rows, rows[1:])):
        raise ValueError('Reference masses must be strictly increasing: %s' % path)
    return tuple(math.log10(r[0]) for r in rows), tuple(math.log10(r[1]) for r in rows)


def table(path):
    """(log10 masses, log10 cross sections) of a two-column table."""
    path = Path(path)
    return _table(str(path.resolve()), path.stat().st_mtime_ns)


def interpolate(path, mass):
    x, y = table(path)
    mass = float(mass)
    if not math.isfinite(mass) or mass <= 0:
        return -1.
    v = math.log10(mass)
    if v < x[0] or v > x[-1]:
        return -1.
    i = max(0, min(bisect.bisect_right(x, v)-1, len(x)-2))
    return 10.**(y[i]+(y[i+1]-y[i])*(v-x[i])/(x[i+1]-x[i]))


def rho_factor(rhoDM):
    rhoDM = float(rhoDM)
    if not math.isfinite(rhoDM) or rhoDM <= 0:
        raise ValueError('rhoDM must be positive, got %s' % rhoDM)
    return FOG_RHO/rhoDM


def evaluate(mass, rhoDM=FOG_RHO, root=ROOT):
    """Neutrino fog at the DM mass [cm^2], rescaled to rhoDM; -1 if unavailable or outside the tables."""
    values, warnings = {}, []
    factor = rho_factor(rhoDM)
    folder = fog_path(FOG_KEYS[0], root).parent
    if not folder.is_dir():
        values.update(dict.fromkeys(FOG_KEYS, -1.))
        values['DD_reference_warnings'] = ['Neutrino fog unavailable: %s is missing. Run `install neutrino_fog` '
                                           'at the MadDM prompt.' % folder]
        return values
    for key in FOG_KEYS:
        path = fog_path(key, root)
        if not path.is_file():
            values[key] = -1.
            warnings.append('%s unavailable: %s is missing. Run `install neutrino_fog` at the MadDM prompt.'
                            % (key, path))
            continue
        try:
            value = interpolate(path, mass)
        except (OSError, ValueError) as error:
            values[key] = -1.
            warnings.append('%s unavailable: %s' % (key, error))
            continue
        if value == -1.:
            warnings.append('DM mass %s GeV outside the %s table; no extrapolation.' % (mass, key))
        values[key] = value*factor if value > 0 else -1.
    values['DD_reference_warnings'] = warnings
    return values


def curves(rhoDM=FOG_RHO, root=ROOT):
    """{key: (masses, cross sections rescaled to rhoDM)} of the available fog tables."""
    factor = rho_factor(rhoDM)
    result = {}
    for key in FOG_KEYS:
        path = fog_path(key, root)
        try:
            x, y = table(path)
        except (OSError, ValueError):
            continue
        result[key] = ([10.**v for v in x], [10.**v*factor for v in y])
    return result


def report(values, rhoDM=FOG_RHO):
    """Lines describing the fog at the current point (screen and MadDM_results.txt)."""
    if not any(k in values for k in FOG_KEYS):
        return []
    lines = ['Neutrino fog (xenon): reference boundary, NOT an exclusion limit. ' + FOG_SOURCE,
             'Rescaled from rho = %.2g GeV/cm^3 (published tables) to rhoDM = %.4g GeV/cm^3.' % (FOG_RHO, rhoDM)]
    lines += ['%s = %.10e # cm^2; -1 means unavailable' % (k, values[k]) for k in FOG_KEYS if k in values]
    lines += ['WARNING: '+s for s in values.get('DD_reference_warnings', [])]
    return lines

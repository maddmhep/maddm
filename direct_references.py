"""Published nuclear-recoil projections and neutrino-fog references.

Tables are interpolated in log mass/log cross section, without extrapolation.
They retain their published halo assumptions: no detector likelihood is rerun.
"""
import bisect
import functools
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FOG_SOURCE = ('https://github.com/cajohare/NeutrinoFog; C. A. J. O\'Hare, '
              'Phys. Rev. Lett. 127, 251802 (2021), arXiv:2109.03116')
XLZD_SOURCE = 'XLZD Design Book, arXiv:2410.17137, Fig. 2; user-supplied smoothed digitization'
FOG_KEYS = ('neutrino_fog_SI', 'neutrino_fog_SD_p', 'neutrino_fog_SD_n')
XLZD_KEYS = ('lim_sigmaN_SI_XLZD_200ty', 'lim_sigmaN_SI_XLZD_1000ty')


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


def interpolate(path, mass):
    path = Path(path)
    x, y = _table(str(path.resolve()), path.stat().st_mtime_ns)
    mass = float(mass)
    if not math.isfinite(mass) or mass <= 0:
        return -1.
    v = math.log10(mass)
    if v < x[0] or v > x[-1]:
        return -1.
    i = max(0, min(bisect.bisect_right(x, v)-1, len(x)-2))
    return 10.**(y[i]+(y[i+1]-y[i])*(v-x[i])/(x[i+1]-x[i]))


def evaluate(mass, xlzd=False, root=ROOT):
    root = Path(root)
    values, warnings = {}, []
    for key, filename in zip(FOG_KEYS, ('Xe_SI.txt', 'Xe_SDp.txt', 'Xe_SDn.txt')):
        path = root/'NeutrinoFog'/'data'/'floors'/filename
        try:
            values[key] = interpolate(path, mass)
        except (OSError, ValueError) as error:
            values[key] = -1.
            warnings.append('%s unavailable: %s' % (key, error))
        if values[key] == -1. and path.is_file():
            warnings.append('DM mass %s GeV outside the %s table; no extrapolation.' % (mass, key))
    if xlzd:
        for exposure, key in zip((200, 1000), XLZD_KEYS):
            path = root/'ExpData'/'XLZD'/('XLZD_SI_90CL_%sty_digitized_smoothed.txt' % exposure)
            try:
                values[key] = interpolate(path, mass)
            except (OSError, ValueError) as error:
                values[key] = -1.
                warnings.append('%s unavailable: %s' % (key, error))
            if values[key] == -1. and path.is_file():
                warnings.append('DM mass %s GeV outside the XLZD %s t yr table; no extrapolation.' % (mass, exposure))
    values['DD_reference_warnings'] = warnings
    values['neutrino_fog_source'] = FOG_SOURCE
    values['neutrino_fog_target'] = 'natural xenon; native published halo (rho=0.3 GeV/cm^3)'
    if xlzd:
        values['XLZD_source'] = XLZD_SOURCE
    return values


def report(values):
    lines = []
    if any(k in values for k in FOG_KEYS):
        lines += ['Neutrino fog: xenon reference boundary, NOT an exclusion limit.', FOG_SOURCE,
                  'Published native halo assumptions; not recalculated for the run-card halo.']
        lines += ['%s = %.10e # cm^2; -1 means unavailable' % (k, values[k]) for k in FOG_KEYS if k in values]
    if any(k in values for k in XLZD_KEYS):
        lines += ['XLZD: projected 90% CL SI sensitivity, NOT a current exclusion.', XLZD_SOURCE,
                  'Native digitized reference curves; no run-card halo rescaling. SI assumes equal proton/neutron couplings.']
        lines += ['%s = %.10e # cm^2; -1 means unavailable' % (k, values[k]) for k in XLZD_KEYS if k in values]
    for exposure, key in zip((200, 1000), XLZD_KEYS):
        upper = values.get(key, -1.)
        if upper <= 0:
            continue
        for nucleon in ('p', 'n'):
            prediction = values.get('sigmaN_SI_'+nucleon, -1.)
            if prediction < 0 or not math.isfinite(prediction):
                continue
            ratio = prediction/upper
            lines.append('XLZD %s t yr SI %s: sigma/reference = %.6e; %s projected curve (All DM)' %
                         (exposure, nucleon, ratio, 'above' if ratio > 1. else 'below or on'))
            xsi = values.get('xsi', 1.)
            if 0. < xsi < 1.:
                ratio *= xsi
                lines.append('XLZD %s t yr SI %s: sigma/reference = %.6e; %s projected curve (thermal fraction)' %
                             (exposure, nucleon, ratio, 'above' if ratio > 1. else 'below or on'))
    lines += ['WARNING: '+s for s in values.get('DD_reference_warnings', [])]
    return lines

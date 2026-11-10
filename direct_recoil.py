"""Native SI/SD nucleon cross sections from the MadDM quark projections."""

import math

# Same zero-momentum matching and mass convention as the RAPIDD electron-port
# revision, so a missing RAPIDD installation does not change the nucleon limits.
QUARKS = ('d', 'u', 's', 'c', 'b', 't')
QUARK_MASSES = (0.00504, 0.00255, 0.101, 1.27, 4.7, 172.)
NUCLEON_MASSES = {'p': 0.93827208943, 'n': 0.93956542194}
GEV2_TO_CM2 = 3.893793721e-28


def has_projection(values):
    """The Fortran returns four -1 values when no projection is available."""
    return len(values) == 4 and all(math.isfinite(v) for v in values) and any(v != -1. for v in values)


def native_nucleon_cross_sections(card, alphas, mass):
    """Return the spin independent and dependent p/n cross sections in cm²."""
    projections = [(i, alphas.get('alpha_' + q, ())) for i, q in enumerate(QUARKS)]
    projections = [(i, values) for i, values in projections if has_projection(values)]
    results = {}
    for target, mn in NUCLEON_MASSES.items():
        tag = target.upper()
        si = sd = -1.
        if projections:
            c1 = axial = 0.
            for i, (si_even, si_odd, sd_even, sd_odd) in projections:
                q = QUARKS[i]
                scalar = card['S' + tag + q] if i < 3 else (2./27.) * card['S' + tag + 'g']
                c1 += si_even * mn * scalar / QUARK_MASSES[i]
                if i < 2:
                    c1 += si_odd * card['V' + tag + q]
                if i < 3:
                    axial += sd_even * card['AV' + tag + q] + sd_odd * card['Sig' + tag + q]
            reduced_mass = mass * mn / (mass + mn)
            factor = GEV2_TO_CM2 * reduced_mass**2 / math.pi
            si, sd = factor * c1**2, 3. * factor * axial**2
        results['sigmaN_SI_' + target] = si
        results['sigmaN_SD_' + target] = sd
    return results

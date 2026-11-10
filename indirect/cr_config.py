"""Shared selection for propagated charged cosmic ray tables."""
from pathlib import Path
import math

ROOT = Path(__file__).resolve().parents[1]
SPECIES = ('pbar', 'dbar', 'positron')
# Local DM density [GeV/cm^3] of the CosmiXs/PPPC propagated tables and of the
# precomputed AMS-02 antiproton limits; they are rescaled as (rhoDM/TABLE_RHO_GEVCM3)^2.
TABLE_RHO_GEVCM3 = 0.4


def value(card, name, default):
    try: return card[name]
    except KeyError: return default


def selection(card):
    profile = str(value(card, 'cr_halo_profile', 'NFW'))
    scheme = str(value(card, 'cr_propagation_scheme', 'QUAINT')).upper()
    setup = str(value(card, 'cr_propagation_set', 'MED')).upper()
    profiles = {p.lower(): p for p in ('NFW', 'Einasto', 'Burkert')}
    if profile.lower() not in profiles or scheme not in ('QUAINT', 'BIG', 'SLIM') or setup not in ('MIN', 'MED', 'MAX'):
        raise ValueError('Require cr_halo_profile=NFW/Einasto/Burkert, cr_propagation_scheme=QUAINT/BIG/SLIM, cr_propagation_set=MIN/MED/MAX')
    return dict(halo=profiles[profile.lower()], scheme=scheme, setup=setup,
                propagation=scheme+'-'+setup)


def table_path(species, card):
    if species not in SPECIES: raise ValueError('Unknown propagated species')
    s = selection(card)
    path = ROOT/'ExpData/CR_flux'/(species+'_flux_'+s['scheme']+'_'+s['setup']+'_'+s['halo']+'.dat.gz')
    if not path.is_file():
        hint = (' Run `install cosmixs_flux` at the main MadDM prompt.'
                if species in ('pbar', 'dbar') else '')
        raise ValueError('Missing propagated table: ' + str(path) + '.' + hint)
    return path


def secondary_path(card, species='pbar', required=True):
    if species not in SPECIES: raise ValueError('Unknown propagated species')
    s = selection(card)
    override = str(value(card, 'pbar_secondary_file', 'auto')) if species == 'pbar' else 'auto'
    if override.lower() != 'auto':
        path = Path(override).expanduser().resolve()
    else:
        directory = ROOT/('ExpData/AMS02_'+species)
        path = directory/('secondary_'+s['scheme']+'_IS.txt')
        if species in ('pbar', 'dbar') and s['scheme'] == 'SLIM' and not path.is_file():
            # The bundled BIG and SLIM secondary inputs were byte-identical.
            # A separately provided SLIM file takes precedence.
            path = directory/'secondary_BIG_IS.txt'
    if not path.is_file():
        if required: raise ValueError('Missing secondary spectrum: '+str(path))
        return None
    return path


def native_channel_limits(card):
    s = selection(card)
    defaults = [('pbar_phi_min',0.),('pbar_phi_max',2.),('pbar_pivot',10.),('pbar_rigidity_min',0.)]
    return (s['halo']=='NFW' and s['propagation']=='QUAINT-MED'
            and str(value(card,'pbar_secondary_file','auto')).lower()=='auto'
            and all(float(value(card,k,v))==v for k,v in defaults)
            and float(value(card,'pbar_rigidity_max',-1.))<=0)


def local_density(card):
    """Local DM density rhoDM [GeV/cm^3], shared by direct detection and the
    Milky Way signals; dwarf J factors are independent inputs.

    Legacy pbar_rho_sun/gce_rho_sun/rho_sun card entries are accepted by
    MadDM for old cards but do not override this common normalization.
    """
    rho = float(value(card, 'rhoDM', .3))
    if not math.isfinite(rho) or rho <= 0:
        raise ValueError('Require positive finite rhoDM in GeV/cm^3')
    return rho

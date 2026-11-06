"""Shared selection for propagated charged cosmic ray tables."""
from pathlib import Path
import math

ROOT = Path(__file__).resolve().parent
SPECIES = ('pbar', 'dbar', 'positron')


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
    path = ROOT/'ExpData/CR_flux'/(species+'_flux_'+s['scheme']+'_'+s['setup']+'_'+s['halo']+'.dat')
    if not path.is_file(): raise ValueError('Missing propagated table: '+str(path))
    return path


def secondary_path(card, species='pbar', required=True):
    if species not in SPECIES: raise ValueError('Unknown propagated species')
    s = selection(card)
    override = str(value(card, 'pbar_secondary_file', 'auto')) if species == 'pbar' else 'auto'
    if override.lower() != 'auto':
        path = Path(override).expanduser().resolve()
    else:
        path = ROOT/('ExpData/AMS02_'+species)/('secondary_'+s['scheme']+'_IS.txt')
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
    """Single Milky Way density; dwarf J factors are independent inputs.

    Legacy pbar_rho_sun/gce_rho_sun/rho_sun card entries are accepted by
    MadDM for old cards but do not override this common normalization.
    """
    rho = float(value(card, 'rholocal', .4))
    if not math.isfinite(rho) or rho <= 0:
        raise ValueError('Require positive finite rholocal in GeV/cm^3')
    return rho

"""Mark physics tests that need the external CosmiXsPPPC table cache."""

import unittest

from indirect import cr_config


def requires_flux(*species):
    available = True
    for name in species:
        try:
            cr_config.table_path(name, {})
        except ValueError:
            available = False
            break
    return unittest.skipUnless(available, 'run install cosmixs_flux at the MadDM prompt')


requires_pbar = requires_flux('pbar')
requires_dbar = requires_flux('dbar')
requires_both = requires_flux('pbar', 'dbar')

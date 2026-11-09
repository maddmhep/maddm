#!/usr/bin/env python3
"""
Antiproton likelihood from partial rates
========================================

Evaluate model partial rates without importing MadGraph.

This is a likelihood/adapter demonstration. It does not generate a model or
compute any matrix elements. A normal MadDM launch computes its own rates.
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from indirect.pbar_interface import AntiprotonLimits, PDG_CHANNELS, write_outputs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mass', type=float, required=True)
    p.add_argument('--rates', type=json.loads, required=True, help='JSON channel:partial sigma*v [cm^3/s]')
    p.add_argument('--dirac', action='store_true')
    p.add_argument('--xsi', type=float, default=1.)
    p.add_argument('--output', type=Path, default=Path('pbar_check_output'))
    a = p.parse_args()
    reverse = {v: k for k, v in PDG_CHANNELS.items()}
    if not isinstance(a.rates, dict) or set(a.rates)-set(reverse):
        p.error('Unknown channel; choices: ' + ', '.join(sorted(reverse)))
    ini = [-1000022, 1000022] if a.dirac else [1000022, 1000022]
    records = [dict(process=ch, initial_states=ini, final_states=list(reverse[ch]),
                    sigmav_cm3s=rate, decay_chains=False, source='user-supplied demonstration rates')
               for ch, rate in a.rates.items()]
    engine = AntiprotonLimits()
    result = engine.evaluate(a.mass, records, 1000022, 8 if a.dirac else 4, a.xsi, {})
    result['validation_scope'] = 'User-supplied rates; this example computes no matrix elements'
    write_outputs(a.output, result, engine)
    print('Status:', result['status'])
    print('95% upper limit [cm^3/s]:', result.get('sigmav95_cm3s'))
    for name, value in result.get('model_points', {}).items():
        print(name + ':', 'EXCLUDED' if value['excluded'] else 'ALLOWED', 'q_upper=', value['q_upper'])
    print('Details:', a.output/'AMS02_pbar_limits.json')


if __name__ == '__main__':
    main()

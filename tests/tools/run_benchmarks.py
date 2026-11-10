"""Reproduce fixed partial-rate adapter benchmarks, not MG model generation.

From the plugin root: python tests/tools/run_benchmarks.py --output /tmp/maddm-benchmarks
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from indirect.pbar_interface import AntiprotonLimits, write_outputs


def record(channel, rate, dirac=False):
    return dict(process=channel, initial_states=[-1000022 if dirac else 1000022, 1000022],
                final_states=[5, -5] if channel == 'bb' else [24, -24],
                sigmav_cm3s=rate, decay_chains=False, source='fixed input for adapter validation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    destination = parser.parse_args().output
    destination.mkdir(parents=True, exist_ok=True)
    engine = AntiprotonLimits()
    maj = engine.evaluate(100., [record('bb', 3.e-26)], 1000022, 4, .5, {})
    dirac = engine.evaluate(100., [record('bb', 3.e-26, True)], 1000022, 8, .5, {})
    mixed = engine.evaluate(200., [record('bb', 2.1e-26), record('WW', .9e-26)], 1000022, 4, .5, {})
    assert abs(dirac['sigmav95_cm3s']/maj['sigmav95_cm3s']-2.) < 2.e-5
    assert maj['model_points']['all_dm']['excluded'] and not maj['model_points']['thermal']['excluded']
    assert not dirac['model_points']['all_dm']['excluded']
    # Choosing native particle/antiparticle tables changes the explicit factor;
    # with a Dirac model it recovers factor one and the cached Majorana limit.
    alternative = engine.evaluate(100., [record('bb', 3.e-26, True)], 1000022, 8, .5,
                                    {'pbar_table_convention': 'particle_antiparticle'})
    assert alternative['source_factor'] == 1.
    assert alternative['sigmav95_cm3s'] == maj['sigmav95_cm3s']
    summary = []
    for label, result in [('bb_majorana_100', maj), ('bb_dirac_100', dirac), ('bb_WW_majorana_200', mixed)]:
        result['validation_scope'] = 'Fixed Milky Way partial rates; no MG matrix elements generated'
        write_outputs(destination/label, result, engine)
        summary.append(dict(name=label, mass_GeV=result['template']['mass_GeV'],
                    branching=result['branching'], source_factor=result['source_factor'],
                    sigmav95_cm3s=result['sigmav95_cm3s'],
                    excluded_all_dm=result['model_points']['all_dm']['excluded'],
                    excluded_thermal=result['model_points']['thermal']['excluded'],
                    delta_chi2_at_limit=result['delta_chi2_at_limit']))
    (destination/'benchmark_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

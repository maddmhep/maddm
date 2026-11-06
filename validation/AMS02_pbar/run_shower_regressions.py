#!/usr/bin/env python3
"""Replay the supplied LHE input and check the actual C++ shower executable.

Requires a compiled updated main101 and matching Pythia XML data directory.
No MadGraph import, model regeneration, or alteration of the supplied run.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import indirect_compat


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-run', type=Path, required=True)
    parser.add_argument('--event-run', default='run_04')
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--pythia-data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run, output = args.input_run.resolve(), args.output.resolve()
    source = indirect_compat.event_file(run/'Indirect_tree_cont/Events'/args.event_run)
    basecard = (run/'Indirect_tree_cont/Source/spectrum.cmnd').read_text()
    environment = dict(os.environ, PYTHIA8DATA=str(args.pythia_data.resolve()))
    result = {}
    for name in ('full_sample', 'one_event', 'empty', 'weighted'):
        directory = output/name/'Events/replay'
        directory.mkdir(parents=True, exist_ok=True)
        settings = directory.parent.parent/'Source/spectrum.cmnd'
        settings.parent.mkdir(exist_ok=True)
        lhe = directory/'maddm_shower_events.lhe'
        if name == 'full_sample':
            result['input_preparation'] = indirect_compat.prepare_lhe(source, lhe)
            text = lhe.read_text()
            # The optional records and model header are preserved in every case.
            first = re.search(r'<event(?:\s[^>]*)?>.*?</event>', text, re.S)
            if not first: raise ValueError('No event in the prepared sample')
            header = text[:first.start()]
            one = first.group(0)
            settings.write_text(basecard)
        else:
            settings.write_text(re.sub(r'(?im)^\s*Main:NumberOfEvents\s*=.*$',
                                       'Main:NumberOfEvents = -1', basecard))
            body = one if name == 'one_event' else ''
            if name == 'weighted':
                # Four 500 GeV neutrinos: two flavor samples with weights 1:3.
                # Incoming total energy is 2 TeV. Each event has two neutrinos
                # and two antineutrinos; the existing output counts neutrinos.
                initial = ('300 -1 0 0 0 0 0 0 0 1000 1000 0 9\n'*2)
                for flavor, weight in ((12, 1), (14, 3)):
                    particles = ''.join(
                        f'{pid} 1 1 2 0 0 {px} {py} 0 500 0 0 9\n'
                        for pid, px, py in [(flavor, 500, 0), (-flavor, -500, 0),
                                            (flavor, 0, 500), (-flavor, 0, -500)])
                    body += f'<event>\n6 1995 {weight} 1000 0.007 0.08\n'+initial+particles+'</event>\n'
                value = settings.read_text()
                value = re.sub(r'(?im)^\s*TimeShower:weakShower\s*=.*$', 'TimeShower:weakShower = off', value)
                value = re.sub(r'(?im)^\s*Main:mDM\s*=.*$', 'Main:mDM = 1000', value)
                settings.write_text(value)
            lhe.write_text(header+body+'\n</LesHouchesEvents>\n')
        indirect_compat.prepare_shower_card(settings, lhe, '', str(args.pythia_data.resolve()))
        for p in directory.glob('*_spectrum_pythia8.dat'): p.unlink()
        log = directory/'pythia8.log'
        with log.open('w') as out:
            process = subprocess.run([str(args.executable.resolve())], cwd=directory,
                                     env=environment, stdout=out, stderr=subprocess.STDOUT)
        match = re.search(r'MADDM_SHOWER_ACCEPTED\s+(\d+)\s+SUM_WEIGHTS\s+([^\s]+)', log.read_text())
        accepted = int(match[1]) if match else 0
        spectra = list(directory.glob('*_spectrum_pythia8.dat'))
        if name == 'empty':
            assert process.returncode != 0 and accepted == 0 and not spectra
        else:
            assert process.returncode == 0, log
            expected = {'full_sample': result['input_preparation']['events'], 'one_event': 1, 'weighted': 2}[name]
            assert accepted == expected, (name, accepted, expected)
            for p in spectra:
                a = np.loadtxt(p)
                indirect_compat.validate_spectrum(10.**a[:, 0], a[:, 1], p.name)
            if name == 'full_sample':
                assert np.max(np.loadtxt(directory/'gammas_spectrum_pythia8.dat')[:, 1]) > 0
            if name == 'weighted':
                electron = np.loadtxt(directory/'neutrinos_e_spectrum_pythia8.dat')
                muon = np.loadtxt(directory/'neutrinos_mu_spectrum_pythia8.dat')
                delta = electron[1, 0]-electron[0, 0]
                np.testing.assert_allclose(np.sum(electron[:, 1])*delta, .5, atol=1.e-7)
                np.testing.assert_allclose(np.sum(muon[:, 1])*delta, 1.5, atol=1.e-7)
        result[name] = dict(returncode=process.returncode, accepted_events=accepted,
                            sum_weights=float(match[2]) if match else None,
                            spectra_written=len(spectra), status='passed')
    (output/'shower_regression_summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

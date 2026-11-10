"""Offline integration checks of actual MadDM method bodies and real pbar fits.

The offline checks do not require MadGraph or Fortran. AST extraction supplies the actual
modified methods with a small run-state harness. No generated matrix element
or MadGraph shower run is simulated or claimed to have been tested.
Run from the plugin root: python -m unittest discover -s tests -v
"""
import ast
import io
import json
import logging
import math
from pathlib import Path
import os
import re
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from indirect import pbar_interface as adapter
from indirect import dbar_interface
from indirect import gce_interface
from indirect import cr_earth
from indirect import halo_model
from indirect.pbar.core import FitConfig, FluxTable, Likelihood

SOURCE = ROOT/'maddm_run_interface.py'
TREE = ast.parse(SOURCE.read_text())
RUN = next(n for n in TREE.body if isinstance(n, ast.ClassDef) and n.name == 'MADDMRunCmd')
STR = next(n for n in TREE.body if isinstance(n, ast.ClassDef) and n.name == 'StrProcess')
ENV = dict(halo_model=halo_model, cr_earth=cr_earth, MadGraph5Error=RuntimeError, pbar_interface=adapter, dbar_interface=dbar_interface, gce_interface=gce_interface, os=os, re=re, math=math, pjoin=os.path.join, logger=logging.getLogger('pbar_test'))
exec(compile(ast.Module(body=[STR], type_ignores=[]), str(SOURCE), 'exec'), ENV)

METHODS = ['is_spectral_finalstate', '_pbar_process_record', '_collect_pbar_rates', '_get_milky_way_rates', '_mw_velocity_info', 'calculate_pbar_limits', 'calculate_dbar_events',
           'print_pbar_results', 'save_pbar_summary', 'save_summary_single',
           'write_scan_output', 'form_n', 'form_s']
for node in RUN.body:
    if isinstance(node, ast.FunctionDef) and node.name in METHODS:
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(SOURCE), 'exec'), ENV)


from flux_test_data import requires_pbar


class ConfigFileContract(dict):
    """MG5 ConfigFile get is an alias of __getitem__, not dict.get."""
    def __getitem__(self, name):
        return dict.__getitem__(self, name.lower())
    get = __getitem__


class Harness:
    def __init__(self, path, fast=True):
        self.dir_path, self.run_name = str(path), 'run_01'
        self._mw_shm_rule = self._cr_earth_summary = None
        self._two2twoLO = fast
        self.maddm_card = ConfigFileContract(pbar_limits=True, pbar_plot=False, vave_indirect_pbar=.00075,
                              vave_indirect_cont=2.e-5, vave_indirect_line=.00075,
                              pbar_phi_min=0., pbar_phi_max=2., pbar_pivot=10.,
                              pbar_table_convention='self_conjugate', rhoDM=.4)
        self.proc_characteristics = {'dm_candidate': [1000022],
                                    'indirect_detection_asked_processes': ['xx_bbx']}
        self.processes_names_map = {'xx_bbx': '1000022 1000022 > 5 -5'}
        self._pbar_fortran_rates = {'xx_bbx': 3.e-26}
        self.last_results = dict(xsi=.5, Fermi_sigmav=7.e-26, taacsID=8.e-30, run='run_01')
        self.str_processes = {'original_marker': object()}
        self.indirect_directories = {'Indirect_tree_cont': True, 'Indirect_tree_line': False}
        self.indirect_directories_cross_section_contribution = {'Indirect_tree_cont': 8.e-30}
        self.calls = []

    def norm_Majorana_Dirac(self):
        return 4

    def launch_indirect(self, force, directory, velocity):
        self.calls.append((self.run_name, directory, velocity))
        self.last_results['taacsID#xx_bbx'] = 3.e-26
        path = Path(self.dir_path)/directory/'Cards'
        path.mkdir(parents=True, exist_ok=True)
        (path/'run_card.dat').write_text('Milky Way run card\n')
        (path/'reweight_card.dat').write_text('Milky Way reweight card\n')


for name in METHODS:
    setattr(Harness, name, ENV[name])


def record(final, rate=1.e-26, initial=('1000022', '1000022'), name='test', decay=False):
    return dict(process=name, initial_states=list(initial), final_states=list(map(str, final)),
                sigmav_cm3s=rate, decay_chains=decay)


class CardAccessTests(unittest.TestCase):
    def test_madgraph_get_signature_and_indexed_fallback(self):
        card=ConfigFileContract(vave_indirect_pbar=.00075)
        with self.assertRaises(TypeError): card.get('vave_indirect_pbar', .1)
        self.assertEqual(adapter.card_value(card,'vave_indirect_pbar',.1),.00075)
        self.assertEqual(adapter.card_value(card,'absent',42),42)

    def test_false_zero_and_case_are_preserved(self):
        card=ConfigFileContract(pbar_limits=False,pbar_phi_min=0.)
        self.assertFalse(adapter.card_value(card,'PBAR_LIMITS',True))
        self.assertEqual(adapter.card_value(card,'PBAR_PHI_MIN',1.),0.)
        self.assertEqual(adapter.card_value({'plain':0.},'plain',1.),0.)

    def test_no_dict_style_card_get_in_antiproton_code(self):
        nodes=ast.parse((ROOT/'indirect/pbar_interface.py').read_text())
        bad=[n for n in ast.walk(nodes) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
             and isinstance(n.func.value,ast.Name) and n.func.value.id=='card' and n.func.attr=='get']
        self.assertEqual(bad,[])
        method=next(n for n in RUN.body if isinstance(n,ast.FunctionDef) and n.name=='calculate_pbar_limits')
        bad=[n for n in ast.walk(method) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
             and n.func.attr=='get' and isinstance(n.func.value,ast.Attribute) and n.func.value.attr=='maddm_card']
        self.assertEqual(bad,[])


class MappingTests(unittest.TestCase):
    def test_all_nineteen_channels_and_order(self):
        self.assertEqual(len(adapter.PDG_CHANNELS), 19)
        for final, channel in adapter.PDG_CHANNELS.items():
            a = adapter.channel_mixture([record(final[::-1])], 1000022, 4)
            self.assertEqual(a['branching'], {channel: 1.})

    def test_aggregation_and_partial_rates(self):
        result = adapter.channel_mixture([record([5, -5], 1.e-26), record([-5, 5], 2.e-26),
                                          record([-24, 24], 1.e-26)], 1000022, 4)
        self.assertAlmostEqual(result['branching']['bb'], .75)
        self.assertAlmostEqual(result['branching']['WW'], .25)
        self.assertAlmostEqual(result['total_sigmav_cm3s']/1.e-26, 4.)

    def test_unknown_channel_not_discarded_or_renormalized(self):
        result = adapter.channel_mixture([record([5, -5]), record([9000001, 9000001])], 1000022, 4)
        self.assertEqual(result['status'], 'unsupported_channels')
        self.assertNotIn('branching', result)
        for bad in [record([5, -5], decay=True), record([5, -5, 21]),
                    record([5, -5], initial=('1000022', '1000023')),
                    record([22, 25])]:
            self.assertEqual(adapter.channel_mixture([bad], 1000022, 4)['status'], 'unsupported_channels')
        with self.assertRaises(ValueError):
            adapter.channel_mixture([record([5, -5], -1.)], 1000022, 4)
        with self.assertRaises(ValueError):
            adapter.channel_mixture([record([5, -5], float('nan'))], 1000022, 4)
        self.assertEqual(adapter.channel_mixture([record([5, -5], 0.)], 1000022, 4)['status'], 'zero_rate')

    def test_dirac_initial_state(self):
        self.assertEqual(adapter.channel_mixture([record([5, -5])], 1000022, 8)['status'], 'unsupported_channels')
        r = record([5, -5], initial=('-1000022', '1000022'))
        self.assertEqual(adapter.channel_mixture([r], 1000022, 8)['status'], 'ok')

    def test_fortran_parser_untruncated_names_and_exponents(self):
        names = {'dmx1dmx1_bbx': '1000022 1000022 > 5 -5'}
        name, value = adapter.parse_fortran_record('pbar_sv:dmx1dmx1_bbx 1.234567890123456D-26 cm^3/s', names)
        self.assertEqual(name, 'dmx1dmx1_bbx')
        self.assertEqual(value, 1.234567890123456e-26)


class RunHookTests(unittest.TestCase):
    def test_numerical_rate_collection_restores_dwarf_state_and_cards(self):
        with tempfile.TemporaryDirectory() as path:
            h = Harness(path, fast=False)
            cards = Path(path)/'Indirect_tree_cont'/'Cards'; cards.mkdir(parents=True)
            (cards/'run_card.dat').write_text('original dwarf run card\n')
            results, processes = h.last_results, h.str_processes
            contributions = dict(h.indirect_directories_cross_section_contribution)
            got = h._collect_pbar_rates()
            self.assertEqual(got[0]['sigmav_cm3s'], 3.e-26)
            self.assertEqual(got[0]['initial_states'], ['1000022', '1000022'])
            self.assertEqual(h.calls, [('run_01_pbar', 'Indirect_tree_cont', .00075)])
            self.assertIs(h.last_results, results)
            self.assertIs(h.str_processes, processes)
            self.assertEqual(h.run_name, 'run_01')
            self.assertEqual(h.last_results['Fermi_sigmav'], 7.e-26)
            self.assertEqual(h.indirect_directories_cross_section_contribution, contributions)
            self.assertEqual((cards/'run_card.dat').read_text(), 'original dwarf run card\n')
            self.assertFalse((cards/'reweight_card.dat').exists())

    def test_rate_collection_restores_after_failure(self):
        with tempfile.TemporaryDirectory() as path:
            h = Harness(path, fast=False)
            results = h.last_results
            def failure(*args):
                raise RuntimeError('Generator failed')
            h.launch_indirect = failure
            with self.assertRaises(RuntimeError):
                h._collect_pbar_rates()
            self.assertIs(h.last_results, results)
            self.assertEqual(h.run_name, 'run_01')

    def test_fast_and_generated_rates_missing_not_used(self):
        with tempfile.TemporaryDirectory() as path:
            h = Harness(path)
            h._pbar_fortran_rates = {}
            with self.assertRaises(RuntimeError):
                h._collect_pbar_rates()
            h._pbar_fortran_rates = {'xx_bbx': 1.e-26}
            h.proc_characteristics['indirect_detection_asked_processes'].append('xx_ww')
            with self.assertRaises(RuntimeError):
                h._collect_pbar_rates()

    def test_disabled_and_failure_never_report_allowed(self):
        with tempfile.TemporaryDirectory() as path:
            h = Harness(path)
            h.maddm_card['pbar_limits'] = False
            h.calculate_pbar_limits(100.)
            self.assertEqual(h.last_results['AMS02_pbar_status'], 'disabled')
            self.assertEqual(h.last_results['AMS02_pbar_excluded'], -1.)
            h.maddm_card['pbar_limits'] = True
            h._pbar_fortran_rates = {}
            h.calculate_pbar_limits(100.)
            self.assertEqual(h.last_results['AMS02_pbar_status'], 'calculation_error')
            buf = io.StringIO(); h.save_pbar_summary(buf)
            self.assertIn('NOT EVALUATED', buf.getvalue())
            self.assertNotIn('= ALLOWED', buf.getvalue())
            payload = json.loads((Path(path)/'output/run_01/AMS02_pbar_limits.json').read_text())
            self.assertEqual(payload['status'], 'calculation_error')

    @requires_pbar
    def test_real_fit_through_actual_run_hook_and_output(self):
        with tempfile.TemporaryDirectory() as path:
            h = Harness(path)
            h.calculate_pbar_limits(100.)
            r = h.last_results
            self.assertEqual(r['AMS02_pbar_status'], 'ok')
            self.assertLess(abs(r['AMS02_pbar_sigmav']/1.77608776e-26-1.), 3.e-6)
            self.assertEqual(r['AMS02_pbar_excluded'], 1.)
            self.assertEqual(r['AMS02_pbar_excluded_thermal'], 0.)
            self.assertEqual(r['Fermi_sigmav'], 7.e-26)
            self.assertEqual(r['taacsID'], 8.e-30)
            output = Path(path)/'output/run_01'
            flux = np.loadtxt(output/'AMS02_pbar_fluxes.dat')
            self.assertEqual(flux.shape, (58, 15))
            np.testing.assert_allclose(flux[:, 3]+flux[:, 4], flux[:, 5], rtol=3.e-12)
            buf = io.StringIO(); h.save_pbar_summary(buf)
            self.assertIn('AMS02_pbar_AllDM = EXCLUDED', buf.getvalue())
            self.assertIn('AMS02_pbar_Thermal = ALLOWED', buf.getvalue())
            with self.assertLogs('pbar_test', level='INFO') as logs:
                h.print_pbar_results()
            self.assertTrue(any('All DM: EXCLUDED' in line for line in logs.output))
            # Exercise the actual complete summary and numerical scan writers.
            h.param_card = types.SimpleNamespace(get_value=lambda *args: 100.)
            h.MadDM_version = '3.3'
            h.line_experiments = []
            h.is_spectral_finalstate = lambda proc: '22' in proc.final_states
            h.maddm_card.update(sigmav_method='inclusive', indirect_flux_source_method='CosmiXs')
            h.last_results.update(tot_SM_xsec=8.e-30, pvalue_th=.8, like_th=0.,
                                  pvalue_nonth=.5, like_nonth=1.)
            for species in ('gammas', 'neutrinos_e', 'neutrinos_mu', 'neutrinos_tau'):
                h.last_results['flux_'+species] = 4.e-12
            h.save_summary_single(indirect=True, fluxes_source=True, fluxes_earth=True)
            summary = (output/'MadDM_results.txt').read_text()
            self.assertNotIn('Flux_gammas', summary)
            self.assertNotIn('Flux_neutrinos', summary)
            self.assertIn('Positron results will come soon.', summary)
            self.assertIn('AMS02_pbar_AllDM = EXCLUDED', summary)
            self.assertIn('AMS02_pbar_Thermal = ALLOWED', summary)
            h.save_MadDM_card = lambda **kwargs: None
            scan = str(output/'scan.txt')
            keys = ['run'] + list(adapter.NUMERIC_KEYS)
            h.write_scan_output(scan, keys, header=True)
            h.write_scan_output(scan, keys)
            rows = [line for line in Path(scan).read_text().splitlines() if line and not line.startswith('#')]
            self.assertEqual(len(rows[0].split()), len(keys))
            self.assertEqual(float(rows[0].split()[keys.index('AMS02_pbar_excluded')]), 1.)
            # Same mass and mixture, changed rate: cached bound, reprofiled model.
            h._pbar_fortran_rates = {'xx_bbx': 1.e-27}
            h._milky_way_rates = None
            h.calculate_pbar_limits(100.)
            self.assertEqual(h.last_results['AMS02_pbar_excluded'], 0.)
            self.assertEqual(len(h._pbar_engine.limits), 1)


class GeneratedRuntimeTests(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('MADDM_PBAR_TEST_RUN'), 'set MADDM_PBAR_TEST_RUN to a generated run with bundled internal/banner.py')
    @requires_pbar
    def test_real_card_class_and_run_hook(self):
        run=Path(os.environ['MADDM_PBAR_TEST_RUN']).resolve()
        sys.path.insert(0,str(run/'Indirect_tree_cont/bin'))
        from internal import banner
        source=ast.parse(SOURCE.read_text())
        classes=[n for n in source.body if isinstance(n,ast.ClassDef) and n.name in ('InvalidMaddmCard','MadDMCard')]
        import math
        env=dict(banner_mod=banner,MDMDIR=str(ROOT),pjoin=os.path.join,os=os,re=re,
                 logger=logging.getLogger('real_card_test'),math=math,pbar_interface=adapter,halo_model=halo_model,
                 secure_float_f77=lambda s:float(s.replace('D','e').replace('d','e')))
        exec(compile(ast.Module(body=classes,type_ignores=[]),str(SOURCE),'exec'),env)
        card=env['MadDMCard'](str(run/'Cards/maddm_card.dat'))
        with self.assertRaises(TypeError): card.get('vave_indirect_pbar',.1)
        self.assertEqual(adapter.card_value(card,'vave_indirect_pbar',0.),.00075)
        self.assertEqual(adapter.card_value(card,'missing_test_setting',42),42)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'include').mkdir();(root/'src').mkdir()
            (root/'include/maddm.inc').write_bytes((run/'include/maddm.inc').read_bytes())
            (root/'src/maddm.f').write_bytes((run/'src/maddm.f').read_bytes())
            card.write_include_file(str(root/'include'))
            assignments=(root/'include/maddm_card.inc').read_text().lower()
            self.assertIn('do_pbar_limits',assignments)
            self.assertIn('vave_indirect_pbar',assignments)
            h=Harness(directory)
            h.maddm_card=card
            h.calculate_pbar_limits(100.)
            self.assertEqual(h.last_results['AMS02_pbar_status'],'ok')
            self.assertLess(abs(h.last_results['AMS02_pbar_sigmav']/1.77608776e-26-1.),3.e-6)
            self.assertEqual(h.last_results['AMS02_pbar_excluded'],1.)
            self.assertTrue((root/'output/run_01/AMS02_pbar_limits.json').is_file())


class SourceWiringTests(unittest.TestCase):
    def test_launch_hook_after_spectral_and_before_output(self):
        launch = next(n for n in RUN.body if isinstance(n, ast.FunctionDef) and n.name == 'do_launch')
        calls = [(n.lineno, n.func.attr) for n in ast.walk(launch)
                 if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)]
        by_name = {name: line for line, name in calls}
        self.assertLess(by_name['launch_spectral_computations'], by_name['calculate_pbar_limits'])
        self.assertLess(by_name['calculate_pbar_limits'], by_name['print_results'])
        self.assertIn('order.extend(pbar_interface.NUMERIC_KEYS)', SOURCE.read_text())

    def test_card_defaults_and_both_template_expansions(self):
        card_class = next(n for n in TREE.body if isinstance(n, ast.ClassDef) and n.name == 'MadDMCard')
        setup = next(n for n in card_class.body if isinstance(n, ast.FunctionDef) and n.name == 'default_setup')
        scope = {"halo_model": halo_model}
        exec(compile(ast.Module(body=[setup], type_ignores=[]), str(SOURCE), 'exec'), scope)
        class Card(dict):
            def __init__(self):
                self.parameters = {}
            def __getitem__(self, name):
                return super().__getitem__(name.lower())
            def add_param(self, name, value, **metadata):
                self[name.lower()] = value
                self.parameters[name.lower()] = metadata
            def fill_jfactors(self):
                pass
        card = Card()
        scope['default_setup'](card)
        self.assertTrue(card['pbar_limits'])
        self.assertEqual(card.parameters['pbar_limits']['fortran_name'], 'do_pbar_limits')
        self.assertTrue(card.parameters['vave_indirect_pbar']['include'])
        for name in ['maddm_card.dat', 'maddm_card_full.dat']:
            rendered = (ROOT/'Templates/Cards'/name).read_text() % card
            self.assertIn('0.00075 = vave_indirect_pbar', rendered)
            self.assertIn('True = pbar_limits', rendered)

    def test_cards_and_fortran_use_independent_velocity(self):
        for name in ['maddm_card.dat', 'maddm_card_full.dat']:
            text = (ROOT/'Templates/Cards'/name).read_text()
            self.assertIn('%(vave_indirect_pbar)s', text)
            self.assertIn('%(pbar_limits)s', text)
        text = (ROOT/'Templates/src/maddm.f').read_text()
        block = text[text.index('c Independent antiproton rates'):text.index('      close(33)')]
        self.assertIn('vID_natural = sqrt(3.d0) * vave_indirect_pbar', block)
        self.assertNotIn('IS_LINE_PROCESS', block)
        self.assertIn("fmt='(A,A,A,ES24.16,A)'", block)
        include = (ROOT/'python_templates/maddm.inc').read_text()
        self.assertIn('common/pbar_velocity/ vave_indirect_pbar', include)


if __name__ == '__main__':
    unittest.main()

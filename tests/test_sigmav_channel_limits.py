"""Execute the actual rate-reader, sigmav hooks and output writers offline.

Matrix-element rates are controlled fixtures, not calculated UFO predictions.
No MadGraph installation or Fortran compiler is available in this environment.
"""
import ast
import collections
from contextlib import nullcontext
import io
import json
import logging
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import six

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pbar_interface as adapter
import dbar_interface
import gce_interface
import cr_earth
import halo_model

SOURCE = ROOT/'maddm_run_interface.py'
TREE = ast.parse(SOURCE.read_text())
RUN = next(x for x in TREE.body if isinstance(x, ast.ClassDef) and x.name == 'MADDMRunCmd')


class Logger:
    """Support MG's optional $MG:BOLD argument with the standard test logger."""
    def __getattr__(self, name):
        method = getattr(logging.getLogger('sigmav_test'), name)
        return lambda message, *args, **kw: method(message, *(() if args == ('$MG:BOLD',) else args), **kw)


COLORS = types.SimpleNamespace(**dict.fromkeys(
    ['ENDC', 'FAIL', 'OKGREEN', 'GRAY', 'BOLD', 'WARNING', 'OKBLUE'], ''))


class RunCard(dict):
    def __init__(self, path):
        super().__init__()
        self.path = path
    def remove_all_cut(self):
        pass
    def write(self, path):
        Path(path).write_text(json.dumps(self))


class EventBackend:
    """One-pb controlled channel fixtures, converted by real launch_indirect."""
    def __init__(self, directory, **kw):
        self.me_dir = directory
        self.Presults = {'xsec_xx_bbx': 1., 'xerr_xx_bbx': 0.,
                         'xsec_xx_wpwm': 2., 'xerr_xx_wpwm': 0.}
        self.commands = []
    def RunWebHandling(self, path):
        return nullcontext()
    def import_command_file(self, commands):
        self.commands.append(commands)
    def do_launch(self, command):
        self.commands.append(command)


ENV = dict(halo_model=halo_model, cr_earth=cr_earth, MadGraph5Error=RuntimeError, pbar_interface=adapter, dbar_interface=dbar_interface, gce_interface=gce_interface, pjoin=os.path.join, MDMDIR=str(ROOT), os=os,
           np=np, math=math, re=re, six=six, collections=collections, HAS_NUMPY=True,
           logger=Logger(), bcolors=COLORS, pb2cm3=2.99792458e-26,
           banner_mod=types.SimpleNamespace(RunCard=RunCard), Indirect_Cmd=EventBackend,
           Indirect_PY8Card=lambda path: {'Main:NumberOfEvents': -1},
           misc=types.SimpleNamespace(TMP_variable=lambda *a: nullcontext(),
                                      MuteLogger=lambda *a, **k: nullcontext()))
for name in ('ExpConstraints', 'StrProcess'):
    node = next(x for x in TREE.body if isinstance(x, ast.ClassDef) and x.name == name)
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(SOURCE), 'exec'), ENV)

METHODS = ('_read_fortran_indirect_record', '_check_inclusive_rates', '_update_indirect_totals',
           'launch_indirect', 'launch_indirect_computations', 'calculate_fermi_limits',
           '_pbar_process_record', '_collect_pbar_rates', '_get_milky_way_rates', '_mw_velocity_info', 'calculate_pbar_limits',
           'print_pbar_results', 'save_pbar_summary', 'print_results', 'print_ind',
           'save_summary_single', 'form_n', 'form_s', 'det_message_screen', 'is_spectral_finalstate')
for node in RUN.body:
    if isinstance(node, ast.FunctionDef) and node.name in METHODS:
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(SOURCE), 'exec'), ENV)


class ForbiddenSpectra:
    def __getattr__(self, name):
        raise AssertionError('Sigmav mode must not access spectra: '+name)


class Harness:
    def __init__(self, directory, method='inclusive'):
        self.dir_path, self.run_name = str(directory), 'run_01'
        self.in_scan_mode, self.param_card_iterator = False, []
        self._mw_shm_rule = self._cr_earth_summary = None
        self._two2twoLO = method == 'inclusive'
        self.mode = dict(indirect='sigmav', relic=False, direct=False,
                         direct_electron=False, capture=False, spectral=False)
        self.maddm_card = dict(sigmav_method=method, indirect_flux_source_method='pythia8',
            vave_indirect_cont=2.e-5, vave_indirect_line=.00075, vave_indirect_pbar=.00075,
            pbar_limits=True, pbar_plot=False, pbar_phi_min=0., pbar_phi_max=2., pbar_pivot=10.,
            rholocal=.4, pbar_table_convention='self_conjugate')
        self.processes_names_map = {'xx_bbx': '1000022.1000022>5.-5',
                                    'xx_wpwm': '1000022.1000022>24.-24'}
        self.proc_characteristics = dict(dm_candidate=[1000022],
            indirect_detection_asked_processes=list(self.processes_names_map))
        self.str_processes, self._pbar_fortran_rates = {}, {}
        self.last_results = dict(xsi=1., **{'Omegah^2': -1., 'x_f': -1., 'sigmav(xf)': -1.})
        self.last_results.update(run=self.run_name, tot_SM_xsec=-1.)
        self.param_card = types.SimpleNamespace(get_value=lambda *args: 100.)
        self.MadDM_version, self.line_experiments = '3.3', []
        self.limits = ENV['ExpConstraints']()
        self.vave_indirect_cont_range = (3.e-6, 1.4e-4)
        self.vave_indirect_line_range = (5.e-4, .001)
        self.indirect_directories = dict(Indirect_tree_cont=True, Indirect_tree_line=False)
        self.indirect_directories_cross_section_contribution = dict.fromkeys(self.indirect_directories, 0.)
        self.Spectra = ForbiddenSpectra()
        self.Fermi = types.SimpleNamespace(Fermi_sigmav_lim=lambda *a, **kw:
                    (_ for _ in ()).throw(AssertionError('Spectrum-based Fermi likelihood called')))
        cards = Path(directory)/'Indirect_tree_cont/Cards'
        cards.mkdir(parents=True)
        (cards/'run_card.dat').write_text('original dwarf card')
        (Path(directory)/'Cards').mkdir()
        (Path(directory)/'Cards/pythia8_card.dat').write_text('unused shower card')
        (Path(directory)/'output'/self.run_name).mkdir(parents=True)

    def norm_Majorana_Dirac(self):
        return 4


for name in METHODS:
    setattr(Harness, name, ENV[name])


class LookupTests(unittest.TestCase):
    def test_all_native_rows_and_no_extrapolation(self):
        limits = adapter.SingleChannelLimits()
        count, empty = 0, 0
        for channel, filename in adapter.CHANNEL_FILES.items():
            for mass, upper in np.loadtxt(ROOT/'ExpData'/filename):
                self.assertEqual(limits.upper_limit(mass, channel), upper)
                count += 1
                empty += upper == -1
            self.assertEqual(limits.upper_limit(4.99, channel), -1.)
            self.assertEqual(limits.upper_limit(100001., channel), -1.)
        self.assertEqual(count, 1216)
        self.assertEqual(empty, 84)
        self.assertEqual(limits.upper_limit(100., [5, -5]), limits.upper_limit(100., 'bb'))

    def test_log_interpolation_and_absent_neighbor(self):
        limits = adapter.SingleChannelLimits()
        mass = math.sqrt(100.*110.)
        expected = math.sqrt(limits.upper_limit(100., 'bb')*limits.upper_limit(110., 'bb'))
        self.assertAlmostEqual(limits.upper_limit(mass, 'bb')/expected, 1., places=12)
        self.assertEqual(limits.upper_limit(45., 'nu_e_nu_e'), -1.)
        self.assertEqual(limits.upper_limit(50., 'nu_e_nu_e') > 0, True)

    def test_density_dirac_and_expdata_api(self):
        constraints = ENV['ExpConstraints']()
        raw = constraints.PBAR_max(100., '5.-5')
        self.assertAlmostEqual(raw/1.776087755903e-26, 1., places=12)
        self.assertEqual(constraints.PBAR_max(100., [-5, 5], source_factor=.5), 2.*raw)
        self.assertEqual(constraints.PBAR_max(100., '5.-5', rho=.8), raw/4.)
        self.assertEqual(constraints.PBAR_max(100., '22.25'), -1.)


class RunTests(unittest.TestCase):
    def test_actual_numerical_collector_preserves_subprocess_names(self):
        cls = next(x for x in TREE.body if isinstance(x, ast.ClassDef) and x.name == 'Indirect_Cmd')
        method = next(x for x in cls.body if isinstance(x, ast.FunctionDef)
                      and x.name == 'make_make_all_html_results')
        class Combined(list):
            xsec, xerru = 3., .1
        for prefix in ('P1_', 'P1995_', 'P10000_'):
            with self.subTest(prefix=prefix):
                combined = Combined([types.SimpleNamespace(name=prefix+'long_dm_name_bbx',
                                    xsec=3., xerru=.1, compute_values=lambda: None)])
                scope = dict(ENV, sum_html=types.SimpleNamespace(collect_result=lambda *a, **kw: combined))
                exec(compile(ast.Module(body=[method], type_ignores=[]), str(SOURCE), 'exec'), scope)
                backend = types.SimpleNamespace(results=types.SimpleNamespace(current={'run_name': 'run_01'}))
                self.assertEqual(scope[method.name](backend), (3., .1))
                self.assertEqual(backend.Presults, {'xsec_long_dm_name_bbx': 3., 'xerr_long_dm_name_bbx': .1})

    def test_actual_fast_and_precise_shortcuts(self):
        selector = next(x for x in TREE.body if isinstance(x, ast.ClassDef) and x.name == 'MadDMSelector')
        methods = {x.name: x for x in selector.body if isinstance(x, ast.FunctionDef)
                   and x.name in ('pass_to_fast_mode', 'pass_to_precise_mode')}
        class Selector:
            availmode = {'forbid_fast': False}
            answer = {'indirect': 'sigmav', 'spectral': 'OFF'}
            def __init__(self): self.commands = []
            def do_set(self, command): self.commands.append(command)
        for name, method in (('pass_to_fast_mode', 'inclusive'), ('pass_to_precise_mode', 'reshuffling')):
            scope = {'logger': Logger()}
            exec(compile(ast.Module(body=[methods[name]], type_ignores=[]), str(SOURCE), 'exec'), scope)
            instance = Selector()
            scope[name](instance)
            self.assertIn('sigmav_method '+method, instance.commands)

    def complete(self, h):
        if h._two2twoLO:
            for name, dwarf, mw in [('xx_bbx', 1.e-27, 3.e-26), ('xx_wpwm', 2.e-27, 6.e-26)]:
                self.assertTrue(h._read_fortran_indirect_record('sigma*v:'+name+' '+str(dwarf), 100., h.last_results))
                self.assertTrue(h._read_fortran_indirect_record('pbar_sv:'+name+' '+str(mw), 100., h.last_results))
            h._check_inclusive_rates()
        else:
            h.launch_indirect(False, 'Indirect_tree_cont', h.maddm_card['vave_indirect_cont'])
        h._update_indirect_totals()
        h.launch_indirect_computations(100.)
        h.calculate_pbar_limits(100.)

    def test_fast_precise_and_madevent_without_any_spectra(self):
        for method in ('inclusive', 'reshuffling', 'madevent'):
            with self.subTest(method=method), tempfile.TemporaryDirectory() as directory:
                h = Harness(directory, method)
                with patch.object(adapter.AntiprotonLimits, 'evaluate', side_effect=AssertionError('Flux fit called')):
                    self.complete(h)
                self.assertGreater(h.last_results['tot_SM_xsec'], 0.)
                self.assertGreater(h.last_results['lim_taacsID#xx_bbx'], 0.)
                self.assertGreater(h.last_results['lim_taacsID#xx_wpwm'], 0.)
                self.assertEqual(h.last_results['Fermi_sigmav'], -1.)
                self.assertEqual(h.last_results['AMS02_pbar_status'], 'tabulated_channels')
                self.assertEqual(h.last_results['AMS02_pbar_excluded'], -1.)
                self.assertFalse(hasattr(h, '_pbar_engine'))
                pbar = h._pbar_result['channels'][0]
                dwarf = h.last_results['taacsID#xx_bbx']
                self.assertGreater(pbar['sigmav_cm3s'], dwarf)
                self.assertEqual(pbar['sigmav95_cm3s'], 1.776087755903e-26)
                if method != 'inclusive':
                    self.assertAlmostEqual(pbar['sigmav_cm3s']/dwarf, .00075/2.e-5, places=10)
                with self.assertLogs('sigmav_test', level='INFO') as logs:
                    h.print_results()
                output = '\n'.join(logs.output)
                self.assertIn('xx_bbx', output)
                self.assertIn('Fermi dSph ul', output)
                self.assertIn('AMS-02 pbar ul', output)
                self.assertNotIn('gamma spectrum', output)
                self.assertNotIn('NOT EVALUATED (calculation_error)', output)
                h.save_summary_single(indirect=True)
                saved = Path(directory)/'output/run_01'
                text = (saved/'MadDM_results.txt').read_text()
                self.assertIn('AMS02_pbar_status = tabulated_channels', text)
                self.assertIn('AllDM=', text)
                self.assertTrue((saved/'AMS02_pbar_channel_limits.dat').is_file())
                self.assertFalse((saved/'AMS02_pbar_fluxes.dat').exists())
                payload = json.loads((saved/'AMS02_pbar_limits.json').read_text())
                self.assertFalse(payload['combined_limit_evaluated'])

    def test_thermal_allowed_excluded_and_no_limit(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory)
            h.last_results['xsi'] = .5
            self.complete(h)
            bb = next(x for x in h._pbar_result['channels'] if x['channel'] == 'bb')
            self.assertTrue(bb['excluded'])
            self.assertFalse(bb['excluded_thermal'])
            self.assertEqual(h.last_results['AMS02_pbar_excluded#bb'], 1.)
            self.assertEqual(h.last_results['AMS02_pbar_excluded_thermal#bb'], 0.)
            self.assertIn('ALLOWED', h.det_message_screen(bb['sigmav95_cm3s'], bb['sigmav95_cm3s']))
            record = dict(process='nu', initial_states=['1000022']*2, final_states=[12, -12], sigmav_cm3s=1.e-26)
            result = adapter.evaluate_channel_limits(10., [record], 1000022, 4, 1., {}, h.limits.PBAR_max)
            self.assertIsNone(result['channels'][0]['excluded'])
            self.assertEqual(adapter.observables(result, .00075)['AMS02_pbar_excluded#nu_e_nu_e'], -1.)

    def test_custom_nuisance_range_profiles_instead_of_reusing_fixed_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory)
            h.maddm_card['pbar_phi_max'] = 1.2
            self.complete(h)
            self.assertEqual(h.last_results['AMS02_pbar_status'], 'tabulated_channels')
            self.assertEqual(h.last_results['AMS02_pbar_excluded'], -1.)
            self.assertEqual(h._pbar_result['bound_source'], 'runtime pure-channel profiles')
            self.assertTrue(all(row['sigmav95_cm3s']>0 for row in h._pbar_result['channels']))

    def test_full_names_legacy_names_and_missing_rates(self):
        names = {'long_dm_name_long_dm_name_mummup': '', 'long_dm_name_long_dm_name_tamtap': ''}
        name = next(iter(names))
        self.assertEqual(adapter.parse_fortran_record('sigma*v: '+name+' 1.2D-26', names), (name, 1.2e-26))
        with self.assertRaises(ValueError):
            adapter.parse_fortran_record('sigma*v:'+name[:12]+' 1.D-26', names)
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory)
            with self.assertRaisesRegex(RuntimeError, 'Missing inclusive'):
                h._check_inclusive_rates()
            self.assertFalse(h._read_fortran_indirect_record('Wimp_Mass: 100.', 100., h.last_results))

    def test_unavailable_dwarf_velocity_does_not_disable_mw_check(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory)
            h.maddm_card['vave_indirect_cont'] = .002
            self.complete(h)
            self.assertEqual(h.last_results['lim_taacsID#xx_bbx'], -1.)
            self.assertGreater(h.last_results['AMS02_pbar_limit#bb'], 0.)

    def test_unsupported_rows_and_channel_aggregation_have_no_global_decision(self):
        constraints = ENV['ExpConstraints']()
        def record(name, final, rate):
            return dict(process=name, initial_states=[1000022]*2, final_states=final, sigmav_cm3s=rate)
        result = adapter.evaluate_channel_limits(100., [record('one', [5, -5], 1.e-26),
            record('two', [-5, 5], 2.e-26), record('BSM', [9999]*2, 1.e-20)],
            1000022, 4, 1., {}, constraints.PBAR_max)
        self.assertEqual(len(result['channels']), 2)
        self.assertAlmostEqual(result['channels'][0]['sigmav_cm3s']/1.e-26, 3.)
        self.assertTrue(result['channels'][0]['excluded'])
        self.assertIsNone(result['channels'][1]['excluded'])
        self.assertFalse(result['combined_limit_evaluated'])
        self.assertEqual(adapter.observables(result, .00075)['AMS02_pbar_excluded'], -1.)


class UploadedCardTests(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('MADDM_PBAR_TEST_RUN'), 'requires uploaded generated run')
    def test_actual_mg5_card_switches_inclusive_flag_and_repairs_driver(self):
        run = Path(os.environ['MADDM_PBAR_TEST_RUN']).resolve()
        sys.path.insert(0, str(run/'Indirect_tree_cont/bin'))
        from internal import banner
        nodes = [x for x in TREE.body if isinstance(x, ast.ClassDef)
                 and x.name in ('InvalidMaddmCard', 'MadDMCard')]
        scope = dict(ENV, banner_mod=banner,
                     secure_float_f77=lambda x: float(x.replace('D', 'e').replace('d', 'e')))
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), 'exec'), scope)
        card = scope['MadDMCard'](str(run/'Cards/maddm_card.dat'))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'include').mkdir(); (root/'src').mkdir()
            (root/'include/maddm.inc').write_bytes((run/'include/maddm.inc').read_bytes())
            (root/'src/maddm.f').write_bytes((run/'src/maddm.f').read_bytes())
            card.set('sigmav_method', 'inclusive')
            card.set('only2to2lo', False, user=False)
            card.write_include_file(str(root/'include'))
            text = (root/'include/maddm_card.inc').read_text()
            self.assertRegex(text, r'ONLY2TO2LO\s*=\s*\.TRUE\.')
            rate_lines = [line for line in (root/'src/maddm.f').read_text().splitlines()
                          if "'sigma*v:'" in line and line[:1].lower() != 'c']
            self.assertEqual(len(rate_lines), 2)
            self.assertTrue(all('A12' not in line for line in rate_lines))
            self.assertEqual((root/'src/maddm.f').read_text().count("'sigma*v:'"), 2)
            card.set('sigmav_method', 'reshuffling')
            card.write_include_file(str(root/'include'))
            text = (root/'include/maddm_card.inc').read_text()
            self.assertRegex(text, r'ONLY2TO2LO\s*=\s*\.FALSE\.')
            card.set('pbar_limits', False)
            card.set('dbar_events', True)
            card.write_include_file(str(root/'include'))
            text = (root/'include/maddm_card.inc').read_text()
            self.assertRegex(text, r'DO_PBAR_LIMITS\s*=\s*\.TRUE\.')
            self.assertFalse(card['pbar_limits'])
            card.set('dbar_events', False)
            card.set('gce_fit', False)
            card.write_include_file(str(root/'include'))
            self.assertRegex((root/'include/maddm_card.inc').read_text(), r'DO_PBAR_LIMITS\s*=\s*\.FALSE\.')


if __name__ == '__main__':
    unittest.main()

"""Regressions for MG5 plugin reweighting, shower cards and combined dwarfs.

Most checks are offline and execute the generated MG5 parser or the bundled
55-dwarf likelihoods with an analytic photon fixture. An optional live test
showers the supplied existing events with the compiled C++ executable and
executes the actual Python shower/reader/Fermi dispatch. It generates no new
MadGraph hard events or matrix-element predictions.
"""
import ast
import json
import math
import os
from pathlib import Path
import shlex
import sys
import tempfile
import types
import unittest
from unittest.mock import patch
import shutil
import stat
import subprocess

import numpy as np
from scipy import stats
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
from scipy.special import gammainc

from test_sigmav_channel_limits import ROOT, TREE, RUN, ENV, Harness
from indirect import indirect_compat
from indirect import fermi_profile
ENV.update(indirect_compat=indirect_compat, fermi_profile=fermi_profile)


class ReweightTests(unittest.TestCase):
    @unittest.skipUnless(os.environ.get('MADDM_PBAR_TEST_RUN'), 'requires a generated MG5 parser')
    def test_reported_failure_and_plugin_command_fix_on_actual_parser(self):
        run = Path(os.environ['MADDM_PBAR_TEST_RUN'])
        source = run/'Indirect_tree_cont/bin/internal/common_run_interface.py'
        tree = ast.parse(source.read_text())
        scope = {'pjoin': os.path.join, 'os': os}
        for name in ('check_decay_events', 'get_events_path'):
            method = next(x for x in ast.walk(tree) if isinstance(x, ast.FunctionDef) and x.name == name)
            exec(compile(ast.Module(body=[method], type_ignores=[]), str(source), 'exec'), scope)
        class Base:
            InvalidCmd = RuntimeError
            mode = 'madevent'
            run_name = 'run_01'
            split_arg = staticmethod(shlex.split)
            def set_run_name(self, name): self.run_name = name
            def do_reweight(self, line):
                self.command = line
                args = self.split_arg(line)
                self.check_decay_events(args)
                self.factory_seen = scope['misc'].from_plugin_import([], 'new_reweight', 'indirect')
                self.selected = args[0]
                self.plugin = next((a.split('=', 1)[1] for a in args if a.startswith('--plugin=')), None)
                return self.selected
        Base.check_decay_events, Base.get_events_path = scope['check_decay_events'], scope['get_events_path']
        cls = next(x for x in TREE.body if isinstance(x, ast.ClassDef) and x.name == 'Indirect_Cmd')
        method = next(x for x in cls.body if isinstance(x, ast.FunctionDef) and x.name == 'do_reweight')
        mini = ast.ClassDef(name='Indirect_Cmd', bases=[ast.Name(id='Base', ctx=ast.Load())],
                            keywords=[], body=[method], decorator_list=[])
        ast.fix_missing_locations(mini)
        fake_misc = types.SimpleNamespace(from_plugin_import=lambda *a, **k: 'old_plugin',
                       TMP_variable=lambda obj, key, value: patch.object(obj, key, value))
        scope = dict(Base=Base, misc=fake_misc,
                     indirect_compat=types.SimpleNamespace(require_reweight_build_tools=lambda: None),
                     MGoutput=types.SimpleNamespace(Indirect_Reweight='active_plugin'))
        exec(compile(ast.Module(body=[mini], type_ignores=[]), 'maddm_run_interface.py', 'exec'), scope)
        with tempfile.TemporaryDirectory() as directory:
            event = Path(directory)/'Events/run_01/unweighted_events.lhe'
            event.parent.mkdir(parents=True)
            event.write_text('<LesHouchesEvents></LesHouchesEvents>')
            old = Base(); old.me_dir = directory
            with self.assertRaisesRegex(RuntimeError, 'No events file corresponding to --mode=indirect'):
                old.do_reweight('--mode=indirect -from_cards')
            for command in ('--mode=indirect -from_cards', 'run_01 --mode=indirect -from_cards',
                            'run_01 --plugin=indirect -from_cards'):
                new = scope['Indirect_Cmd'](); new.me_dir = directory
                self.assertEqual(new.do_reweight(command), str(event))
                self.assertEqual(new.plugin, 'indirect')
                self.assertEqual(new.factory_seen, 'active_plugin')
                self.assertEqual(fake_misc.from_plugin_import([], 'new_reweight', 'indirect'), 'old_plugin')
                self.assertNotIn('--mode=indirect', new.command)

    def test_beam_card_written_without_shower_card(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory, 'reshuffling')
            (Path(directory)/'Cards/pythia8_card.dat').unlink()
            with patch.object(indirect_compat, 'require_reweight_build_tools'):
                h.launch_indirect(False, 'Indirect_tree_cont', .00075)
            card = json.loads((Path(directory)/'Indirect_tree_cont/Cards/run_card.dat').read_text())
            self.assertGreater(card['ebeam1'], 100.)
            self.assertEqual(card['ebeam1'], card['ebeam2'])


class ShowerTests(unittest.TestCase):
    def fixture(self, directory, modern=True):
        root = Path(directory)
        xml = root/'share/Pythia8/xmldoc'
        xml.mkdir(parents=True)
        setting = '<mode name="StringZ:deriveLundPars"/>' if modern else '<flag name="StringZ:deriveBLund"/>'
        (xml/'ZSelection.xml').write_text('<settings>'+setting+'</settings>')
        return root

    def test_modern_and_legacy_tunes_and_duplicate_event_filename(self):
        for modern in (True, False):
            with self.subTest(modern=modern), tempfile.TemporaryDirectory() as directory:
                root = self.fixture(directory, modern)
                card = root/'spectrum.cmnd'
                card.write_text('Main:NumberOfEvents = -1\nStringZ:deriveBLund = on\nStringZ:aLund = 0.601\n'
                                'StringZ:avgZLund = 0.540\nBeams:LHEF=unweighted_events.lhe.gz\n'
                                'Main:subrun=0\nBeams:LHEF=events.lhe.gz\n')
                event = root/'unweighted_events.lhe'
                event.write_text('fixture')
                indirect_compat.prepare_shower_card(card, event, root)
                text = card.read_text()
                self.assertIn('Main:NumberOfEvents = 0', text)
                self.assertEqual(text.count('Beams:LHEF = unweighted_events.lhe\n'), 2)
                self.assertIn('StringZ:avgZLund = 0.540', text)
                self.assertIn('StringZ:aLund = 0.601', text)
                self.assertIn('StringZ:deriveLundPars = 1' if modern else 'StringZ:deriveBLund = on', text)
                if modern: self.assertNotIn('deriveBLund', text)

    def test_event_file_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(OSError): indirect_compat.event_file(root)
            (root/'unweighted_events.lhe').write_text('original')
            (root/'events.lhe.gz').write_bytes(b'reweighted')
            self.assertEqual(indirect_compat.event_file(root).name, 'unweighted_events.lhe')
            self.assertEqual(indirect_compat.event_file(root, True).name, 'events.lhe.gz')

    def test_nan_spectrum_is_rejected(self):
        indirect_compat.validate_spectrum([.1, .2], [0., 1.])
        with self.assertRaisesRegex(ValueError, 'finite nonnegative'):
            indirect_compat.validate_spectrum([.1, .2], [np.nan, np.nan])

    @unittest.skipUnless(os.environ.get('MADDM_PBAR_TEST_RUN'), 'requires a generated MG5 card runtime')
    def test_actual_generated_shower_card_preserves_tune_and_event_file(self):
        run = Path(os.environ['MADDM_PBAR_TEST_RUN']).resolve()
        sys.path.insert(0, str(run/'Indirect_tree_cont/bin'))
        from internal import banner
        cls = next(x for x in TREE.body if isinstance(x, ast.ClassDef) and x.name == 'Indirect_PY8Card')
        scope = {'banner_mod': banner}
        exec(compile(ast.Module(body=[cls], type_ignores=[]), 'maddm_run_interface.py', 'exec'), scope)
        default = str(run/'Cards/pythia8_card_default.dat')
        card = scope['Indirect_PY8Card'](default)
        card.read(str(run/'Cards/pythia8_card.dat'), setter='user')
        with tempfile.TemporaryDirectory() as directory:
            root = self.fixture(directory)
            command = root/'spectrum.cmnd'
            card.write(str(command), default, direct_pythia_input=True)
            before = command.read_text()
            self.assertIn('StringZ:deriveBLund', before)
            indirect_compat.prepare_shower_card(command, root/'unweighted_events.lhe.gz', root)
            after = command.read_text()
            self.assertNotIn('deriveBLund', after)
            self.assertIn('StringZ:deriveLundPars = 1', after)
            lhe_lines = [line for line in after.splitlines() if line.lower().startswith('beams:lhef')]
            self.assertGreaterEqual(len(lhe_lines), 2)
            self.assertTrue(all('unweighted_events.lhe.gz' in line for line in lhe_lines))


class ProfileSolverTests(unittest.TestCase):
    def test_gaussian_upper_branch_and_physical_boundary(self):
        q = stats.chi2.ppf(.95, 1)
        for mean in (0., 5.):
            width = .5
            ll = lambda rate: -.5*((rate-mean)/width)**2
            result = fermi_profile.upper_limit(ll, ll(0.), 1.e-5, 100., q)
            expected = mean+math.sqrt(q)*width
            self.assertAlmostEqual(result['sigmav95_cm3s']/expected, 1., places=6)
            self.assertGreaterEqual(result['sigmav95_cm3s'], result['sigmav_best_cm3s'])
            self.assertAlmostEqual(result['q_at_limit'], q, places=5)

    def test_unbracketed_interval_is_explicit(self):
        with self.assertRaisesRegex(RuntimeError, 'not bracketed'):
            fermi_profile.upper_limit(lambda rate: -.5*(rate/10.)**2, 0., 1.e-5, 1., 3.84)

    def test_vectorized_interpolation_matches_original(self):
        grid = np.array([[0., .1, .4, 1.], [0., .2, .5, 2.], [0., 1., 3., 5.]])
        ll = np.array([[-3., 0., -2., -10.], [-1., 0., -5., -30.], [-2., 0., -1., -5.]])
        old = [interp1d(x, y, bounds_error=False, fill_value=-1.e5) for x, y in zip(grid, ll)]
        for flux in ([.15, .4, 1.1], [0., 0., 0.], [1., 2., 5.], [2., -.1, 3.]):
            self.assertAlmostEqual(fermi_profile.sum_bin_loglikes(grid, ll, np.array(flux)),
                                   sum(float(f(x)) for f, x in zip(old, flux)), places=10)

    def test_photon_zero_intervals_preserved(self):
        spectrum = indirect_compat.photon_spectrum(100., [.01, .1, 1.], [0., 1., 0.])
        self.assertEqual(spectrum(1.), 0.)
        self.assertGreater(spectrum(5.), 0.)
        self.assertEqual(spectrum(101.), 0.)


class CombinedDwarfTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        scope = dict(ENV, HAS_SCIPY=True, stats=stats, interp1d=interp1d,
                     quad=quad, minimize_scalar=minimize_scalar, gammainc=gammainc)
        node = next(x for x in TREE.body if isinstance(x, ast.ClassDef) and x.name == 'Fermi_bounds')
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(ROOT/'maddm_run_interface.py'), 'exec'), scope)
        cls.fermi_class = scope['Fermi_bounds']
        node = next(x for x in RUN.body if isinstance(x, ast.FunctionDef) and x.name == 'read_py8spectra')
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(ROOT/'maddm_run_interface.py'), 'exec'), ENV)

    def test_actual_flux_source_hook_returns_combined_55_dwarf_limit(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory, 'reshuffling')
            h.mode['indirect'] = 'flux_source'
            h.read_py8spectra = types.MethodType(ENV['read_py8spectra'], h)
            h.Spectra = types.SimpleNamespace(spectra={'x': [], 'log10x': [], 'gammas': []}, errors={})
            h.Fermi = self.fermi_class()
            h.last_results.update(taacsID=3.e-27, tot_SM_xsec=3.e-27)
            h.indirect_directories_cross_section_contribution['Indirect_tree_cont'] = 3.e-27
            h.me_cmd = types.SimpleNamespace(run_name='run_01')
            x = np.geomspace(1.e-4, 1., 70)
            # Analytic positive gamma fixture; this is not a model prediction.
            y = 10*np.exp(-10*x)
            def shower(unused):
                out = Path(directory)/'Indirect_tree_cont/Events/run_01'
                out.mkdir(parents=True)
                np.savetxt(out/'gammas_spectrum_pythia8.dat', np.column_stack((np.log10(x), y, .01*y)))
            h.run_pythia8_for_flux = shower
            h.neu_oscillations = lambda: None
            h.calculate_fluxes = lambda: None
            h.launch_indirect_computations(500.)
            self.assertEqual(len(h.Fermi.dw_in), 55)
            self.assertEqual(h.last_results['Fermi_status'], 'ok')
            self.assertGreater(h.last_results['Fermi_sigmav'], 1.e-35)
            self.assertLess(h.last_results['Fermi_sigmav'], 1.e-15)
            self.assertAlmostEqual(h.last_results['Fermi_sigmav']/1.281897403532432e-25, 1., places=5)
            output = Path(directory)/'output/run_01/Fermi_dSph_limits.json'
            result = json.loads(output.read_text())
            self.assertAlmostEqual(result['q_at_limit'], stats.chi2.ppf(.95, 1), places=5)
            self.assertGreater(h.last_results['Fermi_sigmav'], result['sigmav_best_cm3s'])
            pvalue, _ = h.Fermi.Fermi_sigmav_lim(500., x, y, maj_dirac=h.norm_Majorana_Dirac(),
                                               sigmav_th=result['sigmav95_cm3s'])
            self.assertAlmostEqual(pvalue, .95, places=5)

    def test_fast_dwarf_selection_and_model_likelihood_use_same_sample(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory, 'inclusive')
            h.mode['indirect'] = 'flux_source'
            h.maddm_card['indirect_flux_source_method'] = 'cosmixs'
            h.Fermi = self.fermi_class()
            h.last_results.update(taacsID=3.e-27, tot_SM_xsec=3.e-27, xsi=.5)
            x = np.geomspace(1.e-4, 1., 70)
            y = 10*np.exp(-10*x)
            h.Spectra = types.SimpleNamespace(spectra={'x': x, 'gammas': y})
            top = sorted(h.Fermi.dw_in, key=lambda name: (-h.Fermi.dw_in[name]['Jfac'], name))[:10]
            result_file = Path(directory)/'output/run_01/Fermi_dSph_limits.json'
            limits = {}
            for selection, names in [('top10', top), ('all', list(h.Fermi.dw_in))]:
                h.maddm_card['fermi_fast_dwarfs'] = selection
                h.calculate_fermi_limits(500.)
                result = json.loads(result_file.read_text())
                self.assertEqual(result['dwarf_names'], names)
                self.assertEqual(h.last_results['Fermi_n_dwarfs'], len(names))
                limits[selection] = result['sigmav95_cm3s']
                for rate, field in [(3.e-27, 'like_nonth'), (3.e-27*.5**2, 'like_th')]:
                    _, like = h.Fermi.Fermi_sigmav_lim(500., x, y, maj_dirac=h.norm_Majorana_Dirac(),
                                                       sigmav_th=rate, dwarf_selection=selection)
                    self.assertAlmostEqual(h.last_results[field], like)
            self.assertNotAlmostEqual(limits['top10']/limits['all'], 1., places=5)

    @unittest.skipUnless(os.environ.get('MADDM_PYTHIA_EXECUTABLE') and
                         os.environ.get('MADDM_PYTHIA_DATA') and os.environ.get('MADDM_PBAR_TEST_RUN'),
                         'requires compiled main101, matching Pythia data and a generated run')
    def test_actual_shower_reader_and_combined_fermi_dispatch(self):
        """Run the flux hooks and C++ shower on an external event sample."""
        run = Path(os.environ['MADDM_PBAR_TEST_RUN']).resolve()
        event_name = os.environ.get('MADDM_PBAR_TEST_EVENT_RUN', 'run_01')
        try:
            original = indirect_compat.event_file(run/'Indirect_tree_cont/Events'/event_name)
        except (OSError, ValueError) as error:
            self.skipTest(str(error))
        sys.path.insert(0, str(run/'Indirect_tree_cont/bin'))
        from internal import banner
        scope = dict(ENV, stat=stat, subprocess=subprocess,
                     misc=types.SimpleNamespace(get_shell_type=lambda: 'bash',
                                                call=subprocess.call),
                     MadGraph5Error=RuntimeError, ModuleMissing=RuntimeError)
        py8class = next(x for x in TREE.body if isinstance(x, ast.ClassDef) and x.name == 'Indirect_PY8Card')
        py8scope = {'banner_mod': banner}
        exec(compile(ast.Module(body=[py8class], type_ignores=[]), str(ROOT/'maddm_run_interface.py'), 'exec'), py8scope)
        scope['Indirect_PY8Card'] = py8scope['Indirect_PY8Card']
        method = next(x for x in RUN.body if isinstance(x, ast.FunctionDef) and x.name == 'run_pythia8_for_flux')
        exec(compile(ast.Module(body=[method], type_ignores=[]), str(ROOT/'maddm_run_interface.py'), 'exec'), scope)
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory, 'reshuffling')
            h.mode['indirect'] = 'flux_source'
            h.param_card = types.SimpleNamespace(get_value=lambda *args: 1000.)
            h.read_py8spectra = types.MethodType(ENV['read_py8spectra'], h)
            h.run_pythia8_for_flux = types.MethodType(scope['run_pythia8_for_flux'], h)
            h.Spectra = types.SimpleNamespace(spectra={s: [] for s in (
                'x', 'log10x', 'antideuterons', 'gammas', 'antiprotons', 'positrons',
                'neutrinos_e', 'neutrinos_mu', 'neutrinos_tau')}, errors={})
            h.Fermi = self.fermi_class()
            h.last_results.update(taacsID=3.e-27, tot_SM_xsec=3.e-27)
            h.indirect_directories_cross_section_contribution['Indirect_tree_cont'] = 3.e-27
            h.me_cmd = types.SimpleNamespace(run_name='run_01', options={'run_mode': 0})
            xml = Path(os.environ['MADDM_PYTHIA_DATA']).resolve()
            h.mg5 = types.SimpleNamespace(options={'pythia8_path': str(xml.parents[2])})
            event = Path(directory)/'Indirect_tree_cont/Events/run_01'
            event.mkdir(parents=True)
            (Path(directory)/'Indirect_tree_cont/Source').mkdir()
            event_filename = 'unweighted_events.lhe.gz' if original.suffix == '.gz' else 'unweighted_events.lhe'
            shutil.copyfile(original, event/event_filename)
            expected_events = indirect_compat.prepare_lhe(original, event/'count_sample.lhe')['events']
            (event/'count_sample.lhe').unlink()
            for name in ('pythia8_card.dat', 'pythia8_card_default.dat'):
                shutil.copyfile(run/'Cards'/name, Path(directory)/'Cards'/name)
            binary_dir = Path(directory)/'bin/internal'; binary_dir.mkdir(parents=True)
            source = binary_dir/'main101.cc'
            source.write_bytes((ROOT/'Templates/bin/internal/main101.cc').read_bytes())
            shutil.copyfile(os.environ['MADDM_PYTHIA_EXECUTABLE'], binary_dir/'main101')
            (binary_dir/'main101').chmod(0o755)
            _, stamp, identity, _ = indirect_compat.shower_build_state(
                directory, h.mg5.options['pythia8_path'])
            stamp.write_text(json.dumps(identity))
            h.neu_oscillations = lambda: None
            h.calculate_fluxes = lambda: None
            from unittest.mock import patch
            with patch.dict(os.environ, PYTHIA8DATA=str(xml)):
                h.launch_indirect_computations(1000.)
            self.assertEqual(h.last_results['Fermi_status'], 'ok')
            self.assertGreater(h.last_results['Fermi_sigmav'], 0.)
            self.assertTrue(np.isfinite(h.Spectra.spectra['gammas']).all())
            self.assertIn(f'MADDM_SHOWER_ACCEPTED {expected_events}', (event/'pythia8.log').read_text())
            self.assertEqual(json.loads((event/'shower_input.json').read_text())['events'], expected_events)
            self.assertIn('maddm_shower_events.lhe', (Path(directory)/'Indirect_tree_cont/Source/spectrum.cmnd').read_text())


if __name__ == '__main__':
    unittest.main()

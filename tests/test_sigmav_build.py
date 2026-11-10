"""Regression checks for fresh builds, interpreter binding and safe rate reuse."""
import ast
from contextlib import contextmanager
import os
import types
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from indirect import indirect_compat
from indirect import pbar_interface
from test_sigmav_channel_limits import Harness


class BuildTests(unittest.TestCase):
    def test_maddm_output_uses_active_package_exporters(self):
        source = Path(__file__).resolve().parents[1] / 'maddm_interface.py'
        interface = next(node for node in ast.parse(source.read_text()).body
                         if isinstance(node, ast.ClassDef) and node.name == 'MadDM_interface')
        method = next(node for node in interface.body
                      if isinstance(node, ast.FunctionDef) and node.name == 'check_output')

        class MG5Output:
            _export_formats = ['madevent']

            def check_output(self, args):
                # As in MG5, a recognized format bypasses plugin discovery.
                if args[0] not in self._export_formats:
                    raise AssertionError('An older plugin could be selected')
                return args.pop(0), args

        @contextmanager
        def temporary_value(obj, name, value):
            previous = getattr(obj, name)
            setattr(obj, name, value)
            try:
                yield
            finally:
                setattr(obj, name, previous)

        package = types.ModuleType('_maddm_output_test')
        package.MGoutput = types.SimpleNamespace(
            ProcessExporterMadDM=object(), ProcessExporterIndirectD=object())
        test_class = ast.ClassDef(name='MadDM_interface',
                                  bases=[ast.Name(id='MG5Output', ctx=ast.Load())],
                                  keywords=[], body=[method], decorator_list=[])
        module = ast.fix_missing_locations(ast.Module(body=[test_class], type_ignores=[]))
        scope = dict(__package__=package.__name__, MG5Output=MG5Output,
                     misc=types.SimpleNamespace(TMP_variable=temporary_value))
        with patch.dict(sys.modules, {package.__name__: package}):
            exec(compile(module, str(source), 'exec'), scope)
            command = scope['MadDM_interface']()
            for mode, exporter in [('maddm', package.MGoutput.ProcessExporterMadDM),
                                   ('indirect', package.MGoutput.ProcessExporterIndirectD)]:
                self.assertEqual(command.check_output([mode, 'run']), ('plugin', ['run']))
                self.assertIs(command._export_plugin, exporter)
                self.assertEqual(command._export_formats, ['madevent'])
            self.assertEqual(command.check_output(['madevent', 'run']), ('madevent', ['run']))

    @unittest.skipUnless(shutil.which('make'), 'requires make')
    def test_exported_makefile_propagates_failed_source_build(self):
        root_path = Path(__file__).resolve().parents[1]
        tree = ast.parse((root_path/'MGoutput.py').read_text())
        method = next(n for n in ast.walk(tree)
                      if isinstance(n, ast.FunctionDef) and n.name == 'Write_makefile')
        scope = {'pjoin': os.path.join, 'os': os}
        exec(compile(ast.Module(body=[method], type_ignores=[]), 'MGoutput.py', 'exec'), scope)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shutil.copy2(root_path/'Templates/makefile_template', root/'makefile_template')
            export = types.SimpleNamespace(dir_path=directory, proc_characteristic={
                'has_relic_density': False, 'has_direct_detection': False})
            scope['Write_makefile'](export)
            for name in ('Source/DHELAS', 'Source/MODEL', 'matrix_elements', 'lib', 'src'):
                (root/name).mkdir(parents=True, exist_ok=True)
            for name in ('Source/DHELAS', 'Source/MODEL'):
                (root/name/'makefile').write_text('all:\n\ttrue\n')
            # BSD ar (macOS) rejects an archive without members.
            (root/'matrix_elements/makefile').write_text(
                'all:\n\techo stub > stub.txt\n\tar cr libmatrix_elements.a stub.txt\n')
            (root/'src/makefile').write_text('all:\n\tfalse\n')
            (root/'src/maddm.x').write_text('old binary')
            result = subprocess.run(['make', 'all'], cwd=root, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('false', result.stdout)
            self.assertFalse((root/'maddm.x').exists())

    def test_generated_flags_must_match_the_selected_method(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root/'include').mkdir()
            p = root/'include/maddm_card.inc'
            p.write_text('      DO_INDIRECT_DETECTION=.TRUE.\n'
                         '      DO_INDIRECT_SPECTRAL=.FALSE.\n      ONLY2TO2LO=.FALSE.\n')
            with self.assertRaisesRegex(RuntimeError, 'only2to2lo'):
                pbar_interface.check_generated_switches(root, True, False, True)
            pbar_interface.check_generated_switches(root, True, False, False)

    def test_reweight_uses_active_python_and_restores_environment_after_failure(self):
        original = dict(PATH=os.environ.get('PATH'), F2PY='wrong-python-f2py')
        with patch.dict(os.environ, F2PY=original['F2PY']), \
             patch.object(indirect_compat, 'require_reweight_build_tools'):
            try:
                with indirect_compat.reweight_build_environment():
                    command = shlex.split(os.environ['F2PY'])
                    self.assertEqual(command[0], sys.executable)
                    self.assertEqual(Path(command[1]).name, 'f2py_compat.py')
                    self.assertEqual(os.environ['PATH'].split(os.pathsep)[0], str(Path(sys.executable).parent))
                    raise RuntimeError('compile failed')
            except RuntimeError:
                pass
            self.assertEqual(os.environ['F2PY'], original['F2PY'])
            self.assertEqual(os.environ.get('PATH'), original['PATH'])


class ReuseTests(unittest.TestCase):
    def test_equal_velocities_reuse_complete_channel_rates(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory, 'reshuffling')
            h.maddm_card['vave_indirect_pbar'] = h.maddm_card['vave_indirect_cont']
            h.last_results.update({'taacsID#xx_bbx': 1.e-27, 'taacsID#xx_wpwm': 2.e-27})
            def forbidden(*args): raise AssertionError('duplicate numerical launch')
            h.launch_indirect = forbidden
            rows = h._collect_pbar_rates()
            self.assertEqual([r['sigmav_cm3s'] for r in rows], [1.e-27, 2.e-27])
            self.assertTrue(all(r['source'] == 'same_velocity_reuse' for r in rows))
            self.assertEqual(h.run_name, 'run_01')

    def test_different_velocities_require_independent_rates(self):
        with tempfile.TemporaryDirectory() as directory:
            h = Harness(directory, 'reshuffling')
            h.last_results.update({'taacsID#xx_bbx': 1.e-27, 'taacsID#xx_wpwm': 2.e-27})
            calls = []
            def numerical(force, folder, velocity):
                calls.append((h.run_name, velocity))
                h.last_results.update({'taacsID#xx_bbx': 3.e-27, 'taacsID#xx_wpwm': 6.e-27})
            h.launch_indirect = numerical
            rows = h._collect_pbar_rates()
            self.assertEqual(calls, [('run_01_pbar', .00075)])
            self.assertEqual([r['sigmav_cm3s'] for r in rows], [3.e-27, 6.e-27])
            self.assertEqual(h.last_results['taacsID#xx_bbx'], 1.e-27)


if __name__ == '__main__':
    unittest.main()

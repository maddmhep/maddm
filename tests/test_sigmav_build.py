"""Regression checks for fresh builds, interpreter binding and safe rate reuse."""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import indirect_compat
import pbar_interface
from test_sigmav_channel_limits import Harness


class BuildTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('make'), 'requires make')
    def test_old_make_masks_failure_and_new_build_never_reuses_old_binary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'src').mkdir(); (root/'include').mkdir()
            bad = root/'Source/DHELAS'; bad.mkdir(parents=True)
            (bad/'makefile').write_text('all:\n\tfalse\n')
            (root/'makefile').write_text('all:\n\t-cd Source/DHELAS/ && make\n\t-mv src/maddm.x .\n')
            (root/'src/makefile').write_text('objs_all=maddm.o\nincdir=../include\nall:\n\ttrue\n')
            (root/'maddm.x').write_text('old executable')
            (root/'src/maddm.o').write_text('old object')
            (root/'maddm_build.json').write_text('{"switches": "old"}')
            before = subprocess.run(['make', 'all'], cwd=root, capture_output=True)
            self.assertEqual(before.returncode, 0)
            self.assertEqual((root/'maddm.x').read_text(), 'old executable')
            pbar_interface.prepare_generated_build(root)
            after = subprocess.run(['make', 'all'], cwd=root, capture_output=True)
            self.assertNotEqual(after.returncode, 0)
            self.assertFalse((root/'maddm.x').exists())
            self.assertFalse((root/'src/maddm.o').exists())
            self.assertFalse((root/'maddm_build.json').exists())
            pbar_interface.prepare_generated_build(root)
            self.assertEqual((root/'src/makefile').read_text().count('# MadDM: recompile source objects'), 1)
            self.assertTrue((root/'makefile.before_maddm_build_fix').exists())

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

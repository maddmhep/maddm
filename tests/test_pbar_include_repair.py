"""Regression for generated-card variables absent from a stale maddm.inc."""
import ast
import logging
import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import pbar_interface as adapter


class IncludeRepairTests(unittest.TestCase):
    def fixture(self, root, header='      double precision existing\n', driver=None):
        (root/'include').mkdir()
        (root/'src').mkdir()
        (root/'include/maddm.inc').write_text(header)
        (root/'src/maddm.f').write_text(driver or '      program maddm\n      implicit none\n      close(33)\n      end\n')

    def test_reported_legacy_header_and_driver_are_repaired_once(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.fixture(root)
            before = (root/'include/maddm.inc').read_bytes()
            self.assertEqual(adapter.repair_generated_pbar_sources(root), ['include/maddm.inc','src/maddm.f'])
            header = (root/'include/maddm.inc').read_text()
            self.assertIn('double precision vave_indirect_pbar',header)
            self.assertIn('logical do_pbar_limits',header)
            self.assertIn('common/pbar_velocity/ vave_indirect_pbar',header)
            self.assertEqual((root/'include/maddm.inc.before_ams02_pbar').read_bytes(),before)
            self.assertIn("'pbar_sv:'",(root/'src/maddm.f').read_text())
            self.assertEqual(adapter.repair_generated_pbar_sources(root),[])
            self.assertEqual((root/'include/maddm.inc').read_text(),header)

    def test_current_generated_sources_are_not_changed(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            self.fixture(root, (ROOT/'python_templates/maddm.inc').read_text(),
                         (ROOT/'Templates/src/maddm.f').read_text())
            self.assertEqual(adapter.repair_generated_pbar_sources(root),[])
            self.assertEqual(list(root.rglob('*.before_ams02_pbar')),[])

    def test_partial_uppercase_continuation_declarations_not_duplicated(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            header='      DOUBLE PRECISION EXISTING,\n     & VAVE_INDIRECT_PBAR\n'
            self.fixture(root, header, adapter.PBAR_DRIVER_BLOCK+'      close(33)\n')
            self.assertEqual(adapter.repair_generated_pbar_sources(root),['include/maddm.inc'])
            changed=(root/'include/maddm.inc').read_text().lower()
            self.assertEqual(changed.count('double precision'),1)
            self.assertIn('logical do_pbar_limits',changed)
            self.assertEqual(adapter.repair_generated_pbar_sources(root),[])

    def test_incompatible_type_is_rejected_without_mutation(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            self.fixture(root,'      integer do_pbar_limits\n')
            before=(root/'include/maddm.inc').read_bytes()
            with self.assertRaises(ValueError): adapter.repair_generated_pbar_sources(root)
            self.assertEqual((root/'include/maddm.inc').read_bytes(),before)
            self.assertEqual(list(root.rglob('*.before_ams02_pbar')),[])

    def test_invalid_driver_is_rejected_before_header_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            self.fixture(root,driver='      program maddm\n      end\n')
            before=(root/'include/maddm.inc').read_bytes()
            with self.assertRaises(ValueError): adapter.repair_generated_pbar_sources(root)
            self.assertEqual((root/'include/maddm.inc').read_bytes(),before)

    def test_actual_card_include_writer_repairs_before_assignments(self):
        tree=ast.parse((ROOT/'maddm_run_interface.py').read_text())
        cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='MadDMCard')
        method=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='write_include_file')
        env=dict(pbar_interface=adapter,os=os,logger=logging.getLogger('repair_test'))
        exec(compile(ast.Module(body=[method],type_ignores=[]),'maddm_run_interface.py','exec'),env)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            self.fixture(root)
            class Card(dict):
                includepath=['maddm_card.inc']
                def __init__(self): super().__init__(pbar_limits=True, dbar_events=True, gce_fit=False)
                def set(self, key, value, **options): self[key] = value
                def check_validity(self): pass
                def update_system_parameter_for_include(self): pass
                def write_one_include_file(self, output, name, supplied):
                    header=(Path(output)/'maddm.inc').read_text()
                    assert 'logical do_pbar_limits' in header
                    assert 'double precision vave_indirect_pbar' in header
                    (Path(output)/name).write_text('      do_pbar_limits = .true.\n      vave_indirect_pbar = .00075d0\n')
            env['write_include_file'](Card(),str(root/'include'))
            self.assertTrue((root/'include/maddm_card.inc').is_file())


if __name__=='__main__': unittest.main()

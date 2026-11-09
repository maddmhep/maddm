"""Actual failure regressions: LHE numeric blanks, old binaries and build tools."""
import ast
import gzip
import hashlib
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

from indirect import indirect_compat

ROOT = Path(__file__).resolve().parents[1]
PARTICLE = '300 -1 0 0 0 0 0 0 0 1000 1000 0 9\n'
LHE = ('<LesHouchesEvents version="3.0">\n<header>keep header\n\n</header>\n'
       '<init>\n52 52 1000 1000 0 0 0 0 -4 1\n1 0 1 1\n</init>\n'
       '<event>\n\n2 1 1 1000 0.007 0.1\n\n'+PARTICLE+'\n'+PARTICLE+
       '# keep metadata\n<rwgt><wgt id="a">1</wgt></rwgt>\n</event>\n</LesHouchesEvents>\n')


class ShowerInputTests(unittest.TestCase):
    def test_plain_and_gzip_copy_preserve_all_nonblank_records(self):
        for compressed in (False, True):
            with self.subTest(compressed=compressed), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                src, dest = root/('original.lhe.gz' if compressed else 'original.lhe'), root/'shower.lhe'
                if compressed:
                    with gzip.open(src, 'wt') as out: out.write(LHE)
                else: src.write_text(LHE)
                before = src.read_bytes()
                result = indirect_compat.prepare_lhe(src, dest)
                self.assertEqual(result['events'], 1)
                self.assertEqual(result['removed_numeric_blank_lines'], 3)
                self.assertEqual(src.read_bytes(), before)
                self.assertEqual([line for line in LHE.splitlines() if line.strip()],
                                 [line for line in dest.read_text().splitlines() if line.strip()])
                self.assertIn('keep header\n\n</header>', dest.read_text())

    def test_empty_truncated_and_nonfinite_events_are_rejected(self):
        for text in ('<LesHouchesEvents>\n</LesHouchesEvents>', LHE.replace('</event>', ''),
                     LHE.replace('1000 1000 0 9', 'nan 1000 0 9')):
            with self.subTest(text=text[:30]), tempfile.TemporaryDirectory() as directory:
                src, dest = Path(directory)/'original.lhe', Path(directory)/'shower.lhe'
                src.write_text(text)
                with self.assertRaises(ValueError): indirect_compat.prepare_lhe(src, dest)
                self.assertFalse(dest.exists())

    def test_original_is_never_replaced(self):
        with tempfile.TemporaryDirectory() as directory:
            src = Path(directory)/'events.lhe'; src.write_text(LHE)
            with self.assertRaises(ValueError): indirect_compat.prepare_lhe(src, src)
            self.assertEqual(src.read_text(), LHE)

class BuildTests(unittest.TestCase):
    def test_missing_meson_has_actionable_active_environment_diagnostic(self):
        with patch.object(indirect_compat.sys, 'version_info', (3, 13, 0)), \
             patch.object(indirect_compat.shutil, 'which', return_value=None):
            with self.assertRaisesRegex(RuntimeError, 'pip install meson ninja'):
                indirect_compat.require_reweight_build_tools()
        with patch.object(indirect_compat.sys, 'version_info', (3, 13, 0)), \
             patch.object(indirect_compat.shutil, 'which', return_value='/tools/bin/tool'):
            indirect_compat.require_reweight_build_tools()

    def test_actual_reweight_hook_checks_dependencies_before_export(self):
        tree = ast.parse((ROOT/'MGoutput.py').read_text())
        cls = next(x for x in tree.body if isinstance(x, ast.ClassDef) and x.name == 'Indirect_Reweight')
        method = next(x for x in cls.body if isinstance(x, ast.FunctionDef) and x.name == 'create_standalone_directory')
        called = []
        class Base:
            def create_standalone_directory(self, *args, **opts): called.append('export'); return 'exported'
        helper = types.SimpleNamespace(require_reweight_build_tools=lambda: called.append('tools'))
        mini = ast.ClassDef(name='Indirect_Reweight', bases=[ast.Name(id='Base', ctx=ast.Load())],
                            keywords=[], body=[method], decorator_list=[])
        ast.fix_missing_locations(mini)
        scope = dict(Base=Base, indirect_compat=helper)
        exec(compile(ast.Module(body=[mini], type_ignores=[]), 'MGoutput.py', 'exec'), scope)
        instance = scope['Indirect_Reweight'](); instance.shuffling_mode = 'reweight'
        self.assertEqual(instance.create_standalone_directory(), 'exported')
        self.assertEqual(called, ['tools', 'export'])

    def test_existing_binary_requires_a_build_stamp_and_custom_source_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root/'bin/internal/main101.cc'
            source.parent.mkdir(parents=True)
            template = ROOT/'Templates/bin/internal/main101.cc'
            source.write_bytes(template.read_bytes())
            source.with_suffix('').write_bytes(b'old binary')
            path, stamp, identity, build = indirect_compat.shower_build_state(root, '/pythia')
            self.assertTrue(build)
            import json
            stamp.write_text(json.dumps(identity))
            self.assertFalse(indirect_compat.shower_build_state(root, '/pythia')[3])
            source.write_text('// local customization\n')
            self.assertTrue(indirect_compat.shower_build_state(root, '/pythia')[3])
            self.assertEqual(source.read_text(), '// local customization\n')


if __name__ == '__main__':
    unittest.main()

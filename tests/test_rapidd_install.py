from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import rapidd_install


class RAPIDDInstallTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.lib = rapidd_install.library_path(self.root)
        self.lib.parent.mkdir(parents=True)
        source = rapidd_install.ROOT
        shutil.copytree(source / 'rapidd', self.root / 'rapidd',
                        ignore=shutil.ignore_patterns('__pycache__', 'testoutput'))
        shutil.copytree(source / 'lib' / 'efficiency_tables',
                        self.root / 'lib' / 'efficiency_tables')
        shutil.copyfile(source / 'lib' / 'CMakeLists.txt', self.root / 'lib' / 'CMakeLists.txt')

    def test_existing_library_must_load_or_it_triggers_rebuild(self):
        self.lib.write_bytes(b'not a native library')
        with self.assertRaises(rapidd_install.UnusableLibrary):
            rapidd_install.verify(self.root)
        with patch.object(rapidd_install, 'build') as build:
            rapidd_install.ensure(self.root)
        build.assert_called_once_with(self.root)

    def test_unloadable_replacement_does_not_overwrite_existing_library(self):
        original = b'original library'
        self.lib.write_bytes(original)
        cache = self.lib.parent / 'CMakeCache.txt'
        cache.write_text('old machine')

        def build_unloadable(command, cwd):
            self.assertNotEqual(cwd, cache.parent)
            self.assertFalse((cwd / cache.name).exists())
            if '--build' in command:
                (cwd / self.lib.name).write_bytes(b'unloadable replacement')

        with patch.object(rapidd_install.shutil, 'which', return_value='cmake'), \
                patch.object(rapidd_install, '_run_build', side_effect=build_unloadable):
            with self.assertRaisesRegex(rapidd_install.RAPIDDError, 'newly built RAPIDD'):
                rapidd_install.build(self.root)
        self.assertEqual(self.lib.read_bytes(), original)
        self.assertEqual(cache.read_text(), 'old machine')


if __name__ == '__main__':
    unittest.main()

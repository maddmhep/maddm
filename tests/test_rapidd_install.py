from pathlib import Path
import os
import shutil
import subprocess
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
        # Build the minimal layout locally; the test suite must work without
        # an initialized submodule or a network connection.
        for folder in ('efficiency_tables', 'experiments', 'halo_table'):
            (self.root / 'lib' / folder).mkdir()
        (self.root / 'lib' / 'CMakeLists.txt').touch()
        package = self.root / 'rapidd'
        package.mkdir()
        (package / '__init__.py').touch()
        (package / 'electron.py').touch()
        for index in range(1, 5):
            path = self.root / 'lib' / 'atomic_responses' / ('atomic_response_%d' % index)
            path.mkdir(parents=True)
            for shell in range(11):
                (path / ('Xe_%d.txt' % shell)).touch()
        (package / 'alpha_map.py').write_text('from . import core\n')
        (package / 'core.py').write_text(
            'from ctypes import CDLL\nfrom pathlib import Path\n'
            'CDLL(str(Path(__file__).resolve().parents[1] / "lib" / "build" / %r))\n'
            % self.lib.name)

    def test_launch_loader_does_not_install_missing_library(self):
        with patch.object(rapidd_install, 'build') as build:
            with self.assertRaisesRegex(rapidd_install.UnusableLibrary, 'library is missing'):
                rapidd_install.load(self.root)
        build.assert_not_called()

    def test_working_library_is_not_rebuilt_during_preparation(self):
        with patch.object(rapidd_install, 'verify') as verify, \
                patch.object(rapidd_install, 'build') as build:
            rapidd_install.ensure(self.root)
        verify.assert_called_once_with(self.root)
        build.assert_not_called()

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


@unittest.skipUnless(shutil.which('git'), 'Git required for checkout tests')
class RAPIDDFetchTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.base = Path(tmp.name)
        self.upstream = self.base / 'upstream'
        self.upstream.mkdir()
        self.git(self.upstream, 'init', '--quiet')
        for file in ('lib/CMakeLists.txt', 'rapidd/core.py'):
            path = self.upstream / file
            path.parent.mkdir(exist_ok=True)
            path.write_text('reference revision')
        self.git(self.upstream, 'add', '.')
        self.git(self.upstream, 'commit', '--quiet', '-m', 'old nuclear-only revision')
        self.old_revision = self.git(self.upstream, 'rev-parse', 'HEAD').strip()
        (self.upstream / 'rapidd' / 'electron.py').write_text('electron port')
        self.git(self.upstream, 'add', 'rapidd/electron.py')
        self.git(self.upstream, 'commit', '--quiet', '-m', 'reference')
        self.revision = self.git(self.upstream, 'rev-parse', 'HEAD').strip()
        (self.upstream / 'rapidd/core.py').write_text('newer revision')
        self.git(self.upstream, 'commit', '--quiet', '-am', 'newer')

    @staticmethod
    def git(path, *args):
        return subprocess.check_output(
            ['git', '-c', 'user.name=MadDM tests', '-c', 'user.email=test@example.invalid',
             *args], cwd=path, text=True)

    def test_first_install_fetches_pinned_zip_revision_and_reuses_it(self):
        root = self.base / 'zip' / 'vendor' / 'RAPIDD_for_DM'
        root.mkdir(parents=True)
        with patch.multiple(rapidd_install, SOURCE_URL=self.upstream.as_uri(),
                            SOURCE_REVISION=self.revision):
            rapidd_install.fetch_source(root)
        self.assertEqual(self.git(root, 'rev-parse', 'HEAD').strip(), self.revision)
        with patch.object(rapidd_install, '_git') as git:
            rapidd_install.fetch_source(root)
            git.assert_not_called()

    def test_git_checkout_initializes_submodule_at_index_revision(self):
        project = self.base / 'project'
        project.mkdir()
        self.git(project, 'init', '--quiet')
        root = project / 'vendor' / 'RAPIDD_for_DM'
        with patch.dict(os.environ, {'GIT_ALLOW_PROTOCOL': 'file'}):
            self.git(project, 'submodule', 'add', '--quiet', self.upstream.as_uri(),
                     'vendor/RAPIDD_for_DM')
            self.git(root, 'checkout', '--quiet', '--detach', self.revision)
            self.git(project, 'add', 'vendor/RAPIDD_for_DM')
            self.git(project, 'submodule', 'deinit', '--force', '--all')
            rapidd_install.fetch_source(root)
        self.assertEqual(self.git(root, 'rev-parse', 'HEAD').strip(), self.revision)

    def test_old_installation_is_upgraded_to_the_pinned_electron_port(self):
        root = self.base / 'zip' / 'vendor' / 'RAPIDD_for_DM'
        root.parent.mkdir(parents=True)
        self.git(self.base, 'clone', '--quiet', str(self.upstream), str(root))
        self.git(root, 'checkout', '--quiet', '--detach', self.old_revision)
        self.assertFalse(rapidd_install.has_electron_source(root))
        with patch.multiple(rapidd_install, SOURCE_URL=self.upstream.as_uri(),
                            SOURCE_REVISION=self.revision):
            rapidd_install.fetch_source(root)
        self.assertEqual(self.git(root, 'rev-parse', 'HEAD').strip(), self.revision)
        self.assertTrue(rapidd_install.has_electron_source(root))

    def test_failed_zip_download_leaves_destination_empty(self):
        root = self.base / 'zip' / 'vendor' / 'RAPIDD_for_DM'
        root.mkdir(parents=True)
        with patch.object(rapidd_install, 'SOURCE_URL', (self.base / 'absent').as_uri()):
            with self.assertRaisesRegex(rapidd_install.RAPIDDError, 'download failed'):
                rapidd_install.fetch_source(root)
        self.assertEqual(list(root.iterdir()), [])
        self.assertEqual(list(root.parent.iterdir()), [root])


if __name__ == '__main__':
    unittest.main()

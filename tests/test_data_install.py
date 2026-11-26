"""Check data installation with local upstream archives and no network access."""

import gzip
import hashlib
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import data_install


class PythonDependencyTests(unittest.TestCase):
    def test_conda_targets_the_python_prefix(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'env with spaces'
            (prefix / 'conda-meta').mkdir(parents=True)
            with patch.object(data_install.sys, 'prefix', str(prefix)), \
                 patch.object(data_install.shutil, 'which', return_value='/tools/conda'):
                self.assertEqual(data_install.conda_install_command(['matplotlib']),
                                 ['/tools/conda', 'install', '--prefix', str(prefix),
                                  '--yes', 'matplotlib'])
            with patch.object(data_install.sys, 'prefix', directory), \
                 patch.dict(data_install.os.environ, CONDA_PREFIX=str(prefix)):
                self.assertIsNone(data_install.conda_install_command(['matplotlib']))

    def test_manual_choice_does_not_install_or_exit_the_program(self):
        question = []
        def ask(text, default, choices, **kwargs):
            question.append((text, default, choices, kwargs))
            return 'n'
        with patch.object(data_install, 'missing_python_requirements', return_value=['numpy','scipy','matplotlib']), \
             patch.object(data_install, 'python_install_command', return_value=('pip', ['python','-m','pip'])), \
             patch.object(data_install.subprocess, 'run') as run:
            with self.assertRaisesRegex(data_install.DataInstallError, 'install py_requirements'):
                data_install.ensure_python_requirements(ask=ask)
        run.assert_not_called()
        self.assertIn('through pip [y]', question[0][0])
        self.assertIn('yourself [n]', question[0][0])
        self.assertEqual(question[0][1:3], ('n', ['y','n']))
        self.assertIs(question[0][3]['timeout'], False)

    def test_explicit_install_uses_selected_manager_and_checks_result(self):
        for manager, command in [('conda',['conda','install','numpy','meson']),
                                 ('pip',['python','-m','pip','install','numpy','meson'])]:
            with self.subTest(manager=manager), \
                 patch.object(data_install, 'missing_python_requirements', side_effect=[['numpy','meson'], []]), \
                 patch.object(data_install, 'python_install_command', return_value=(manager, command)), \
                 patch.object(data_install.subprocess, 'run') as run:
                data_install.install_python_requirements()
                run.assert_called_once_with(command, check=True)

    def test_acceptance_installs_and_no_manager_preserves_manual_route(self):
        with patch.object(data_install, 'missing_python_requirements', return_value=['scipy']), \
             patch.object(data_install, 'python_install_command', return_value=('conda', ['conda'])), \
             patch.object(data_install, 'install_python_requirements') as install:
            data_install.ensure_python_requirements(ask=lambda *args, **kwargs: 'y', build_tools=True)
            install.assert_called_once_with(True, tools_bin=None)
        with patch.object(data_install, 'missing_python_requirements', return_value=['scipy']), \
             patch.object(data_install, 'python_install_command', return_value=(None, None)), \
             patch.object(data_install.subprocess, 'run') as run:
            with self.assertRaisesRegex(data_install.DataInstallError, 'yourself'):
                data_install.ensure_python_requirements()
            run.assert_not_called()

    def test_pip_is_bound_to_the_running_python_and_is_not_assumed(self):
        with patch.object(data_install, 'conda_install_command', return_value=None), \
             patch.object(data_install, 'find_spec', return_value=object()):
            self.assertEqual(data_install.python_install_command(['matplotlib']),
                             ('pip', [data_install.sys.executable, '-m', 'pip', 'install', 'matplotlib']))
        with patch.object(data_install, 'conda_install_command', return_value=None), \
             patch.object(data_install, 'find_spec', return_value=None):
            self.assertEqual(data_install.python_install_command(['matplotlib']), (None, None))


def record(name, content, **details):
    return dict(file=name, bytes=len(content), sha256=hashlib.sha256(content).hexdigest(),
                **details)


class DataInstallTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.manifest = dict(repository='https://github.com/test/data', revision='pinned',
                             sources=[])

    def flux_archive(self, species, correct=True):
        name = species + '.dat'
        contents = (species + '\n1 2 3\n').encode()
        archive = self.root / (species + '.zip')
        with zipfile.ZipFile(archive, 'w') as handle:
            handle.writestr(species + '/' + name, contents if correct else b'changed')
        self.manifest['sources'].append(record(name + '.gz', contents,
                                              source_archive=archive.name, upstream_file=name))
        return contents

    def test_flux_install_reuses_tables_and_repairs_a_missing_one(self):
        expected = {species: self.flux_archive(species) for species in ('pbar', 'dbar')}
        cache = self.root / 'cache'

        def download(url, path):
            self.assertIn('/test/data/pinned/', url)
            shutil.copyfile(self.root / Path(url).name, path)

        with patch.object(data_install, '_manifest', return_value=self.manifest), \
                patch.object(data_install, '_download', side_effect=download) as downloader:
            data_install.install_cosmixs_flux(cache)
            self.assertEqual(downloader.call_count, 2)
            self.assertEqual(set(cache.iterdir()),
                             {cache / (species + '.dat.gz') for species in expected})
            for species, contents in expected.items():
                self.assertEqual(gzip.decompress((cache / (species + '.dat.gz')).read_bytes()),
                                 contents)
            data_install.install_cosmixs_flux(cache)
            self.assertEqual(downloader.call_count, 2)
            (cache / 'dbar.dat.gz').unlink()
            data_install.install_cosmixs_flux(cache)
            self.assertEqual(downloader.call_count, 3)
            self.assertTrue(downloader.call_args.args[0].endswith('/dbar.zip'))

    def test_bad_table_does_not_publish_a_partial_installation(self):
        self.flux_archive('dbar')
        self.flux_archive('pbar', correct=False)
        cache = self.root / 'cache'

        def download(url, path):
            shutil.copyfile(self.root / Path(url).name, path)

        with patch.object(data_install, '_manifest', return_value=self.manifest), \
                patch.object(data_install, '_download', side_effect=download):
            with self.assertRaisesRegex(data_install.DataInstallError, 'SHA-256 mismatch'):
                data_install.install_cosmixs_flux(cache)
        self.assertEqual(list(cache.iterdir()), [])

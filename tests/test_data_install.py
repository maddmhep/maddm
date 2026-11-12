"""Check data installation with local upstream archives and no network access."""

import gzip
import hashlib
from io import StringIO
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
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

    def test_batch_reports_every_missing_dependency_without_prompting(self):
        import __main__
        output = StringIO()
        options = SimpleNamespace(plugin='maddm', file='commands.txt')
        with patch.object(__main__, 'options', options, create=True), \
             patch.object(data_install, 'find_spec', side_effect=lambda name: object() if name == 'six' else None), \
             patch.object(data_install, 'conda_install_command', return_value=None), \
             patch.object(data_install.sys, 'stderr', output), \
             patch('builtins.input', side_effect=AssertionError('batch startup prompted')):
            with self.assertRaises(SystemExit) as exit_code:
                data_install.check_python_dependencies('PLUGIN.maddm')
        self.assertEqual(exit_code.exception.code, 1)
        self.assertIn('Missing Python packages: numpy, scipy, matplotlib', output.getvalue())
        self.assertIn('-m pip install numpy scipy matplotlib', output.getvalue())


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

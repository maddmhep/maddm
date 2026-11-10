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

"""Install the pinned propagation tables at generation time."""

import gzip
import hashlib
from http.client import HTTPException
import json
import logging
from pathlib import Path
import shutil
import tempfile
from urllib.parse import quote
from urllib.request import urlopen
import zipfile
import zlib


logger = logging.getLogger('madgraph.plugin.maddm')
ROOT = Path(__file__).resolve().parent
FLUX_DIR = ROOT / 'ExpData' / 'CR_flux'


class DataInstallError(RuntimeError):
    """A required dataset could not be installed or found."""


def _manifest(name):
    with (ROOT / 'ExpData' / name).open(encoding='utf-8') as stream:
        return json.load(stream)


def _raw_url(manifest, path):
    repo = manifest['repository'].split('github.com/', 1)[1]
    return 'https://raw.githubusercontent.com/%s/%s/%s' % (
        repo, manifest['revision'], quote(path))


def _download(url, path):
    with urlopen(url, timeout=60) as source, path.open('wb') as target:
        shutil.copyfileobj(source, target)


def _copy_verified(source, target, entry):
    digest = hashlib.sha256()
    length = 0
    for chunk in iter(lambda: source.read(256 * 1024), b''):
        length += len(chunk)
        if length > entry['bytes']:
            break
        digest.update(chunk)
        target.write(chunk)
    if length != entry['bytes'] or digest.hexdigest() != entry['sha256']:
        raise DataInstallError('Size or SHA-256 mismatch for %s' % entry['file'])


def _available(path, size=None):
    try:
        return path.is_file() and (path.stat().st_size == size if size is not None
                                   else path.stat().st_size > 0)
    except OSError:
        return False


def install_cosmixs_flux(destination=FLUX_DIR):
    """Download only the two annihilation ZIPs, verify each table, and cache gzip files."""
    destination = Path(destination)
    manifest = _manifest('CR_flux_provenance.json')
    missing = [entry for entry in manifest['sources']
               if not _available(destination / entry['file'])]
    if not missing:
        return

    logger.info('Downloading the pbar and Dbar fluxes from CosmiXsPPPC: %s',
                manifest['repository'])
    try:
        destination.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.cosmixs-', dir=destination) as work:
            stage = Path(work)
            for name in sorted({entry['source_archive'] for entry in missing}):
                archive = stage / name
                _download(_raw_url(manifest, name), archive)
                with zipfile.ZipFile(archive) as upstream:
                    for entry in missing:
                        if entry['source_archive'] != name:
                            continue
                        member = Path(name).stem + '/' + entry['upstream_file']
                        with upstream.open(member) as source, (stage / entry['file']).open('wb') as raw:
                            with gzip.GzipFile(filename='', mode='wb', compresslevel=6,
                                               fileobj=raw, mtime=0) as compressed:
                                _copy_verified(source, compressed, entry)
                archive.unlink()
            # Do not publish any table until all requested inputs have passed validation.
            for entry in missing:
                filename = entry['file']
                (stage / filename).replace(destination / filename)
    except (DataInstallError, OSError, KeyError, zipfile.BadZipFile, EOFError, HTTPException, zlib.error) as error:
        raise DataInstallError('Could not install CosmiXsPPPC flux tables in %s: %s'
                               % (destination, error)) from error
    logger.info('Installed %d CosmiXsPPPC flux tables in %s.', len(missing), destination)

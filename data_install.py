"""Check and install MadDM's runtime dependencies, build tools and data."""

import gzip
import hashlib
from http.client import HTTPException
import importlib
from importlib.util import find_spec
import json
import logging
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
from urllib.parse import quote
from urllib.request import urlopen
import zipfile
import zlib


logger = logging.getLogger('madgraph.plugin.maddm')
ROOT = Path(__file__).resolve().parent
FLUX_DIR = ROOT / 'ExpData' / 'CR_flux'
PYTHON_REQUIREMENTS = ('six', 'numpy', 'scipy', 'matplotlib')
PYTHON_INSTALL_RESTART = 'MADDM_PYTHON_INSTALL_RESTART'
ZLIB_ARCHIVE = (
    'https://zlib.net/fossils/zlib-1.3.2.tar.gz',
    'bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16',
)
GSL_ARCHIVE = (
    'https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz',
    '6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190',
)
COSMIXS_SOURCE_ARCHIVE = (
    'https://madgraph.mi.infn.it/Downloads/maddm/PPPC4DMID_V1.16.tar.gz',
    'd7bbeadb1a8656e18e33c197f6e44b62d6a6fb094bc2183d0eeba6dea3af20ba',
)


class DataInstallError(RuntimeError):
    """A MadDM dataset or dependency could not be installed."""


def conda_install_command(packages):
    """Use Conda only when it manages the interpreter's actual environment."""
    prefix = Path(sys.prefix)
    if not (prefix / 'conda-meta').is_dir():
        return None
    conda = shutil.which(os.environ.get('CONDA_EXE') or 'conda') or shutil.which('conda')
    if conda:
        return [conda, 'install', '--prefix', str(prefix), '--yes', *packages]
    return None


def check_python_dependencies(plugin_name):
    """Run before importing modules that depend on the scientific Python stack."""
    import __main__

    options = getattr(__main__, 'options', None)
    cli = getattr(options, 'plugin', None) == plugin_name.rsplit('.', 1)[-1]
    restarted = cli and os.environ.pop(PYTHON_INSTALL_RESTART, None)
    missing = [name for name in PYTHON_REQUIREMENTS if find_spec(name) is None]
    if not missing:
        for name in PYTHON_REQUIREMENTS:
            try:
                importlib.import_module(name)
            except (ImportError, OSError, RuntimeError) as error:
                message = 'MadDM could not import %s with %s.\n%s' % (name, sys.executable, error)
                if cli:
                    raise SystemExit(message) from None
                raise ImportError(message) from error
        return

    command = conda_install_command(missing)
    manager = 'Conda' if command else 'pip'
    if command is None:
        command = [sys.executable, '-m', 'pip', 'install', *missing]
    message = ('MadDM requires %s.\n'
               'Missing Python packages: %s\n'
               'Python: %s\nEnvironment: %s\n'
               'Install them with:\n  %s' %
               (', '.join(PYTHON_REQUIREMENTS), ', '.join(missing), sys.executable,
                sys.prefix, shlex.join(command)))
    if not cli:
        raise ImportError(message)
    print(message, file=sys.stderr)
    interactive = (not getattr(options, 'file', None) and not getattr(options, 'web', False)
                   and not getattr(__main__, 'args', ()) and sys.stdin is not None
                   and sys.stdin.isatty())
    if interactive and not restarted:
        try:
            answer = input('Install the missing packages into this environment using %s? [y/N] ' % manager)
            if answer.strip().lower() in ('y', 'yes'):
                subprocess.run(command, check=True)
                original = getattr(sys, 'orig_argv', None)
                argv = [sys.executable, *(original[1:] if original else sys.argv)]
                env = dict(os.environ, **{PYTHON_INSTALL_RESTART: '1'})
                os.execve(sys.executable, argv, env)
        except (EOFError, KeyboardInterrupt):
            pass
        except (OSError, subprocess.CalledProcessError) as error:
            print('MadDM dependency setup failed: %s' % error, file=sys.stderr)
    raise SystemExit(1)


def heptools_dir(mg5_dir, configured=None):
    return Path(mg5_dir) / Path(configured or 'HEPTools').expanduser()


def verified_archive(cache, url, sha256):
    """Reuse a verified download and publish replacements atomically."""
    def valid(path):
        if not path.is_file():
            return False
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(256 * 1024), b''):
                digest.update(chunk)
        return digest.hexdigest() == sha256

    archive = cache / url.rsplit('/', 1)[-1]
    try:
        if valid(archive):
            return archive
        cache.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.download-', dir=cache) as work:
            temporary = Path(work) / archive.name
            _download(url, temporary)
            if not valid(temporary):
                raise DataInstallError('SHA-256 mismatch for %s' % archive.name)
            temporary.replace(archive)
    except (OSError, HTTPException) as error:
        raise DataInstallError('Could not download %s: %s' % (url, error)) from error
    return archive


def install_cosmixs_source(destination, heptools):
    """Install the MadDM-format CosmiXs table without replacing legacy PPPC tables."""
    destination = Path(destination)
    if _available(destination / 'CosmiXs.npy'):
        return
    logger.info('Downloading CosmiXs source spectra [arXiv:2312.01153].')
    archive = verified_archive(heptools / 'maddm_cache', *COSMIXS_SOURCE_ARCHIVE)
    try:
        destination.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.cosmixs-', dir=destination) as work:
            temporary = Path(work) / 'CosmiXs.npy'
            with tarfile.open(archive, 'r:gz') as tables:
                with tables.extractfile('PPPC4DMID/CosmiXs.npy') as source, temporary.open('wb') as target:
                    shutil.copyfileobj(source, target)
            temporary.replace(destination / temporary.name)
    except (OSError, KeyError, tarfile.TarError) as error:
        raise DataInstallError('Could not install CosmiXs source spectra: %s' % error) from error


def reweight_tools_bin(heptools):
    environment = heptools / ('maddm_build_tools_py%d%d' % sys.version_info[:2])
    return environment / ('Scripts' if os.name == 'nt' else 'bin')


def install_reweight_tools(heptools):
    """Reuse build tools or install them with Conda or in a HEPTools virtualenv."""
    from .indirect.indirect_compat import missing_reweight_build_tools

    bin_dir = reweight_tools_bin(heptools)
    missing = missing_reweight_build_tools(bin_dir)
    if not missing:
        return
    command = conda_install_command(missing)
    try:
        if command is None:
            logger.info('Installing %s with pip in %s.', ', '.join(missing), bin_dir.parent)
            python = bin_dir / ('python.exe' if os.name == 'nt' else 'python')
            if not python.is_file():
                subprocess.run([sys.executable, '-I', '-m', 'venv', str(bin_dir.parent)], check=True)
            command = [str(python), '-I', '-m', 'pip', 'install', *missing]
        else:
            logger.info('Installing %s with Conda in %s.', ', '.join(missing), sys.prefix)
        subprocess.run(command, check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise DataInstallError('Could not install Meson and Ninja: %s' % error) from error
    remaining = missing_reweight_build_tools(bin_dir)
    if remaining:
        raise DataInstallError('Build executables are still missing: %s' % ', '.join(remaining))


def install_gsl(heptools):
    """Build GNU GSL in the same HEPTools prefix used by MadGraph."""
    prefix = heptools / 'GSL'
    if (prefix / 'bin' / 'gsl-config').is_file() and (prefix / 'include/gsl/gsl_version.h').is_file():
        return prefix
    logger.info('Building GNU GSL in %s.', prefix)
    archive = verified_archive(heptools / 'maddm_cache', *GSL_ARCHIVE)
    try:
        with tempfile.TemporaryDirectory(prefix='.gsl-', dir=heptools) as work:
            subprocess.run(['tar', '-xzf', str(archive)], cwd=work, check=True)
            source = Path(work) / 'gsl-2.8'
            subprocess.run(['sh', 'configure', '--prefix=' + str(prefix)], cwd=source, check=True)
            jobs = min(os.cpu_count() or 1, 8)
            subprocess.run(['make', '-j', str(jobs)], cwd=source, check=True)
            subprocess.run(['make', 'install'], cwd=source, check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise DataInstallError('Could not build GNU GSL: %s' % error) from error
    return prefix


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

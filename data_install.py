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
PYTHON_REQUIREMENTS = ('numpy', 'scipy', 'matplotlib')
BUILD_REQUIREMENTS = ('ninja', 'meson')
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


def missing_python_requirements(build_tools=False, tools_bin=None):
    """Return missing or broken imports, respecting externally installed tools."""
    missing = []
    for name in PYTHON_REQUIREMENTS:
        try:
            importlib.import_module(name)
        except (ImportError, OSError, RuntimeError, ValueError, AttributeError) as error:
            logger.debug('Cannot import %s: %s', name, error)
            missing.append(name)
    if build_tools:
        paths = [str(Path(sys.executable).parent), str(Path(sys.prefix) / 'bin')]
        if tools_bin is not None:
            paths.append(str(tools_bin))
        paths.append(os.environ.get('PATH', ''))
        search = os.pathsep.join(paths)
        missing.extend(name for name in BUILD_REQUIREMENTS if shutil.which(name, path=search) is None)
    return missing


def python_install_command(packages):
    """Select the manager for this interpreter; never assume pip is installed."""
    command = conda_install_command(packages)
    if command:
        return 'conda', command
    if find_spec('pip') is not None:
        return 'pip', [sys.executable, '-m', 'pip', 'install', *packages]
    return None, None


def install_python_requirements(build_tools=True, tools_bin=None):
    """Explicit `install py_requirements` target, also used after consent."""
    missing = missing_python_requirements(build_tools, tools_bin=tools_bin)
    if not missing:
        logger.info('All requested MadDM Python requirements are available.')
        return
    manager, command = python_install_command(missing)
    if command is None:
        raise DataInstallError('Missing packages: %s. No conda or pip is available for %s. '
                               'Install these packages yourself in this Python environment, then retry.' %
                               (', '.join(missing), sys.executable))
    logger.info('Installing %s through %s in %s: %s', ', '.join(missing), manager,
                sys.prefix, shlex.join(command))
    try:
        subprocess.run(command, check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise DataInstallError('Dependency installation failed: %s. Install the packages yourself or '
                               'retry `install py_requirements`.' % error) from error
    importlib.invalidate_caches()
    remaining = missing_python_requirements(build_tools, tools_bin=tools_bin)
    if remaining:
        raise DataInstallError('Packages still unavailable: %s. Check the installation in %s, then retry.' %
                               (', '.join(remaining), sys.executable))


def ensure_python_requirements(ask=None, build_tools=False, tools_bin=None):
    """Offer installation only when the user requests a numerical calculation."""
    missing = missing_python_requirements(build_tools, tools_bin=tools_bin)
    if not missing:
        return
    required = PYTHON_REQUIREMENTS + (BUILD_REQUIREMENTS if build_tools else ())
    logger.warning('MadDM numerical calculations require: %s. Missing or unusable: %s. '
                   'Python: %s. You can use `install py_requirements` at the main MadDM prompt '
                   'or install the packages yourself.', ', '.join(required), ', '.join(missing), sys.executable)
    manager, command = python_install_command(missing)
    if ask is not None and command is not None:
        question = ('Do you want to install the required packages through %s [y] '
                    'or do you want to install them yourself [n]?' % manager)
        # Declining leaves model import and diagram generation available.
        answer = ask(question, 'n', ['y', 'n'], timeout=False)
        if str(answer).strip().lower() == 'y':
            install_python_requirements(build_tools, tools_bin=tools_bin)
            return
    raise DataInstallError('Missing packages: %s. Install them in %s yourself, or run '
                           '`install py_requirements`, then launch the calculation again. '
                           'Model import and diagram generation remain available.' %
                           (', '.join(missing), sys.executable))


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


NEUTRINO_FOG_REVISION = '0df3d0cd2322602dc147a66d94055f7f64879d80'
NEUTRINO_FOG_ARCHIVE = (
    'https://codeload.github.com/cajohare/NeutrinoFog/tar.gz/' + NEUTRINO_FOG_REVISION,
    '43c53d675ef844919090e11cc1ba87d4bc6d55f1edf432725298aeb3752d44aa',
)


def install_neutrino_fog(root=ROOT):
    """Fetch the pinned upstream project when the DM-nucleon diagrams are generated, never import its code."""
    root = Path(root)
    destination = root/'NeutrinoFog'
    required = ['data/floors/Xe_SI.txt', 'data/floors/Xe_SDp.txt', 'data/floors/Xe_SDn.txt', 'LICENSE']
    marker = destination/'maddm_revision.txt'
    if marker.is_file() and marker.read_text().strip() == NEUTRINO_FOG_REVISION and all(
            (destination/p).is_file() for p in required):
        return destination
    if destination.exists():
        raise DataInstallError('Existing %s is incomplete or from another revision; preserve or move it before retrying.' % destination)
    root.mkdir(parents=True, exist_ok=True)
    logger.info('Downloading NeutrinoFog for DD references: https://github.com/cajohare/NeutrinoFog (arXiv:2109.03116).')
    with tempfile.TemporaryDirectory(prefix='.neutrino-fog-', dir=root) as work:
        stage = Path(work)
        archive = verified_archive(stage, *NEUTRINO_FOG_ARCHIVE)
        with tarfile.open(archive, 'r:gz') as source:
            prefix = 'NeutrinoFog-'+NEUTRINO_FOG_REVISION
            # Explicit copying avoids links and archive path traversal on all supported Python versions.
            unpack = stage/'project'
            for member in source.getmembers():
                parts = Path(member.name).parts
                if not parts or parts[0] != prefix or '..' in parts or member.name.startswith('/'):
                    raise DataInstallError('Unsafe NeutrinoFog archive path')
                target = unpack.joinpath(*parts[1:])
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                elif member.isfile():
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with source.extractfile(member) as inp, target.open('wb') as out:
                        shutil.copyfileobj(inp, out)
                else:
                    raise DataInstallError('Unexpected link in NeutrinoFog archive')
        if not all((unpack/p).is_file() for p in required):
            raise DataInstallError('NeutrinoFog archive is missing required reference tables')
        (unpack/'maddm_revision.txt').write_text(NEUTRINO_FOG_REVISION+'\n')
        unpack.rename(destination)
    return destination

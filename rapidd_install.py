"""Fetch, build and check RAPIDD against the Python running MadDM."""

import importlib
import logging
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

logger = logging.getLogger('madgraph.plugin.maddm')
ROOT = Path(__file__).resolve().parent / 'vendor' / 'RAPIDD_for_DM'
SOURCE_URL = 'https://github.com/cheekyparticle/RAPIDD_for_DM.git'
# Match the electron recoil port pinned by MadDM dev.
SOURCE_REVISION = 'fa7046f8d4d0a9433fc5cb0d770fd2e8ec012a42'


class RAPIDDError(RuntimeError):
    pass


class UnusableLibrary(RAPIDDError):
    pass


def library_path(root=ROOT):
    name = {'Linux': 'libRAPIDD.so', 'Darwin': 'libRAPIDD.dylib',
            'Windows': 'RAPIDD.dll'}.get(platform.system())
    if not name:
        raise RAPIDDError('RAPIDD is not supported on %s.' % platform.system())
    return Path(root).resolve() / 'lib' / 'build' / name


def has_source(root=ROOT):
    root = Path(root)
    return (root / 'lib' / 'CMakeLists.txt').is_file() and (root / 'rapidd' / 'core.py').is_file()


def has_electron_source(root=ROOT):
    root = Path(root)
    return has_source(root) and (root / 'rapidd' / 'electron.py').is_file()


def _git(git, args, cwd):
    result = subprocess.run([git, *args], cwd=cwd, capture_output=True, text=True)
    if result.returncode:
        detail = (result.stderr.strip() or result.stdout.strip()).splitlines()[-10:]
        raise RAPIDDError('RAPIDD download failed (git %s): %s'
                          % (' '.join(args), '\n'.join(detail)))
    return result.stdout


def fetch_source(root=ROOT):
    """Initialize the recorded submodule or fetch the ZIP release's revision."""
    root = Path(root).absolute()
    if has_electron_source(root):
        return
    git = shutil.which('git')
    if not git:
        raise RAPIDDError('Git is required to download RAPIDD. Install Git and retry.')

    project = root.parent.parent
    submodule = 'vendor/RAPIDD_for_DM'
    logger.info('Downloading RAPIDD into %s.', root)
    try:
        index = ''
        if root == project / submodule and (project / '.git').exists():
            index = _git(git, ['ls-files', '--stage', '--', submodule], project)
        if not has_source(root) and any(line.startswith('160000 ') and line.endswith('\t' + submodule)
               for line in index.splitlines()):
            _git(git, ['submodule', 'update', '--init', '--recursive', '--', submodule], project)
        elif not has_source(root):
            # A ZIP has no Git index. Fetch the release revision into a temporary
            # checkout, so an interrupted download never looks installed.
            if root.is_symlink() or (root.exists() and
                                    (not root.is_dir() or any(root.iterdir()))):
                raise RAPIDDError('Incomplete RAPIDD source at %s. Move it aside and retry.'
                                  % root)
            root.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(prefix='.rapidd-download-', dir=root.parent) as tmp:
                checkout = Path(tmp) / 'source'
                checkout.mkdir()
                _git(git, ['init', '--quiet'], checkout)
                _git(git, ['remote', 'add', 'origin', SOURCE_URL], checkout)
                _git(git, ['fetch', '--depth=1', '--no-tags', 'origin', SOURCE_REVISION], checkout)
                _git(git, ['checkout', '--quiet', '--detach', 'FETCH_HEAD'], checkout)
                _git(git, ['submodule', 'update', '--init', '--recursive'], checkout)
                if not has_electron_source(checkout):
                    raise RAPIDDError('The downloaded RAPIDD revision lacks the electron recoil port.')
                if root.is_dir():
                    root.rmdir()  # Fails rather than overwriting an unexpected nonempty folder.
                checkout.rename(root)
        if has_source(root) and not has_electron_source(root):
            # Older installations can be updated in place without deleting local changes.
            if not (root / '.git').exists() or _git(
                    git, ['status', '--porcelain', '--untracked-files=no'], root).strip():
                raise RAPIDDError('Move the existing RAPIDD folder aside or preserve its local '
                                  'changes, then run `install RAPIDD` again: %s' % root)
            _git(git, ['fetch', '--depth=1', '--no-tags', SOURCE_URL, SOURCE_REVISION], root)
            _git(git, ['checkout', '--quiet', '--detach', 'FETCH_HEAD'], root)
    except OSError as error:
        raise RAPIDDError('Could not install RAPIDD in %s: %s' % (root, error)) from error
    if not has_electron_source(root):
        raise RAPIDDError('RAPIDD electron recoil source is missing after initializing %s.' % root)


_PROBE = """
import sys
sys.path.insert(0, sys.argv[1])
try:
    from rapidd.alpha_map import sigmaSI_nucleon_mdm, sigmaSD_nucleon_mdm
    from rapidd.calc_dRdE import DDrate_save, DDrate_electron_save
    from rapidd.core import halo_path
    from rapidd import electron
except (OSError, AttributeError, ImportError) as exc:
    print(str(exc), file=sys.stderr)
    sys.exit(2)
"""


def verify(root=ROOT):
    """Check the wrapper and its native symbols in a fresh Python process."""
    lib = library_path(root)
    if not lib.is_file():
        raise UnusableLibrary('RAPIDD library is missing: %s' % lib)
    if not has_electron_source(root):
        raise UnusableLibrary('Installed RAPIDD source lacks the electron recoil port.')
    data = Path(root) / 'lib' / 'atomic_responses'
    if any(len(list((data / ('atomic_response_%d' % index)).glob('Xe_*.txt'))) != 11
           for index in range(1, 5)):
        raise UnusableLibrary('The RAPIDD xenon atomic responses are incomplete: %s' % data)
    try:
        process = subprocess.run([sys.executable, '-c', _PROBE, str(Path(root).resolve())],
                                 capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise RAPIDDError('Could not check RAPIDD with %s: %s' % (sys.executable, error)) from error
    if process.returncode:
        detail = process.stderr.strip() or 'exit status %s' % process.returncode
        error_type = UnusableLibrary if process.returncode == 2 or process.returncode < 0 else RAPIDDError
        raise error_type(detail)


def _run_build(command, directory):
    result = subprocess.run(command, cwd=directory, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        detail = '\n'.join(result.stdout.strip().splitlines()[-20:])
        raise RAPIDDError('RAPIDD build failed (%s):\n%s' % (' '.join(command), detail))


def _verify_new_build(root):
    try:
        verify(root)
    except UnusableLibrary as error:
        hint = (' If macOS reports Library Validation, use a Python environment that '
                'allows locally built native libraries.' if platform.system() == 'Darwin' else '')
        raise RAPIDDError('The newly built RAPIDD could not be loaded by %s.%s\n%s'
                          % (sys.executable, hint, error)) from error


def build(root=ROOT):
    """Configure afresh, then install only a library verified by this Python."""
    root = Path(root).resolve()
    target = library_path(root)
    cmake = shutil.which('cmake')
    if not cmake:
        raise RAPIDDError('CMake is required to build RAPIDD. Install CMake, GSL and a C compiler '
                          'in your environment, then run `install RAPIDD` in the MadDM shell.')
    fetch_source(root)
    source = root / 'lib'
    logger.info('Building RAPIDD from source for %s (%s).', sys.executable, platform.machine())
    try:
        # The temporary build lives on the same filesystem as the installed library.
        # An old CMake cache or a binary from another machine is never reused.
        with tempfile.TemporaryDirectory(prefix='.maddm-rapidd-', dir=source) as temp:
            stage = Path(temp)
            work = stage / 'lib' / 'build'
            work.mkdir(parents=True)
            configure = [cmake, '-DCMAKE_BUILD_TYPE=Release']
            if platform.system() == 'Darwin':
                configure.append('-DCMAKE_OSX_ARCHITECTURES=%s' % platform.machine())
            prefix = Path(sys.prefix)
            if (prefix / 'include' / 'gsl' / 'gsl_math.h').is_file() and not os.environ.get('GSL_ROOT_DIR'):
                configure.append('-DGSL_ROOT_DIR=%s' % prefix)
            _run_build(configure + [str(source)], work)
            _run_build([cmake, '--build', '.', '--config', 'Release', '--target', 'RAPIDD'], work)

            candidate = work / target.name
            if not candidate.is_file() and (work / 'Release' / target.name).is_file():
                (work / 'Release' / target.name).replace(candidate)
            shutil.copytree(root / 'rapidd', stage / 'rapidd',
                            ignore=shutil.ignore_patterns('__pycache__', 'testoutput'))
            for name in ('efficiency_tables', 'atomic_responses', 'experiments', 'halo_table'):
                shutil.copytree(source / name, stage / 'lib' / name)
            shutil.copy2(source / 'CMakeLists.txt', stage / 'lib' / 'CMakeLists.txt')
            _verify_new_build(stage)
            target.parent.mkdir(parents=True, exist_ok=True)
            candidate.replace(target)
            _verify_new_build(root)
    except OSError as error:
        raise RAPIDDError('Could not build or install RAPIDD in %s: %s' % (root, error)) from error
    logger.info('RAPIDD installed and verified with the active Python.')


def ensure(root=ROOT):
    if not has_electron_source(root):
        build(root)
        return
    try:
        verify(root)
    except UnusableLibrary as error:
        if library_path(root).is_file():
            logger.info('The existing RAPIDD library could not be loaded; rebuilding it locally.')
            logger.debug('RAPIDD loader: %s', error)
        try:
            build(root)
        except RAPIDDError as failure:
            raise RAPIDDError('%s\nOriginal RAPIDD load check: %s' % (failure, error)) from failure


def load(root=ROOT):
    """Load an installed RAPIDD without downloading or building anything."""
    root = Path(root).resolve()
    package = sys.modules.get('rapidd')
    if package is None:
        verify(root)
        if not sys.path or sys.path[0] != str(root):
            sys.path.insert(0, str(root))
    elif Path(getattr(package, '__file__', None) or '.').resolve().parent != root / 'rapidd':
        raise RAPIDDError('A different RAPIDD package is already imported. Restart MadDM to use '
                          'the RAPIDD installation managed by MadDM.')
    try:
        return tuple(importlib.import_module('rapidd.' + name)
                     for name in ('alpha_map', 'calc_dRdE', 'core'))
    except (ImportError, OSError, AttributeError) as error:
        raise RAPIDDError('Could not import RAPIDD with %s: %s' % (sys.executable, error)) from error

"""Build and check the bundled RAPIDD against the Python running MadDM."""

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


_PROBE = """
import sys
sys.path.insert(0, sys.argv[1])
try:
    from rapidd.alpha_map import sigmaSI_nucleon_mdm, sigmaSD_nucleon_mdm
    from rapidd.calc_dRdE import DDrate_save
    from rapidd.core import halo_path
except (OSError, AttributeError) as exc:
    print(str(exc), file=sys.stderr)
    sys.exit(2)
"""


def verify(root=ROOT):
    """Check the wrapper and its native symbols in a fresh Python process."""
    lib = library_path(root)
    if not lib.is_file():
        raise UnusableLibrary('RAPIDD library is missing: %s' % lib)
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
    if not has_source(root):
        raise RAPIDDError('RAPIDD source is missing. Run `install RAPIDD` in the MadDM shell.')
    cmake = shutil.which('cmake')
    if not cmake:
        raise RAPIDDError('CMake is required to build RAPIDD. Install CMake, GSL and a C compiler '
                          'in your environment, then run `install RAPIDD` in the MadDM shell.')
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
            shutil.copytree(source / 'efficiency_tables', stage / 'lib' / 'efficiency_tables')
            _verify_new_build(stage)
            target.parent.mkdir(parents=True, exist_ok=True)
            candidate.replace(target)
            _verify_new_build(root)
    except OSError as error:
        raise RAPIDDError('Could not build or install RAPIDD in %s: %s' % (root, error)) from error
    logger.info('RAPIDD installed and verified with the active Python.')


def ensure(root=ROOT):
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
    """Prepare RAPIDD for an existing generated run and return its modules."""
    root = Path(root).resolve()
    package = sys.modules.get('rapidd')
    if package is None:
        ensure(root)
        if not sys.path or sys.path[0] != str(root):
            sys.path.insert(0, str(root))
    elif Path(getattr(package, '__file__', None) or '.').resolve().parent != root / 'rapidd':
        raise RAPIDDError('A different RAPIDD package is already imported. Restart MadDM to use '
                          'its bundled RAPIDD.')
    try:
        return tuple(importlib.import_module('rapidd.' + name)
                     for name in ('alpha_map', 'calc_dRdE', 'core'))
    except (ImportError, OSError, AttributeError) as error:
        raise RAPIDDError('Could not import RAPIDD with %s: %s' % (sys.executable, error)) from error

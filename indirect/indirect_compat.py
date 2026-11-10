"""Small MG5/Pythia compatibility helpers; no MadGraph import required."""
from functools import lru_cache
from contextlib import contextmanager
from pathlib import Path
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import shlex
import sys


# MadDM shortcuts may also be entered as `add process <shortcut>`.
MADDM_PROCESS_SHORTCUTS = frozenset((
    'relic', 'relic_density', 'direct', 'direct_detection', 'indirect',
    'indirect_detection', 'spectral', 'indirect_spectral_features',
))


def is_maddm_process_shortcut(process):
    """Recognize MadDM shortcuts without matching physical process syntax."""
    tokens = process.split()
    return bool(tokens and tokens[0].lower() in MADDM_PROCESS_SHORTCUTS and '>' not in process)


def require_reweight_build_tools():
    """MG5 uses F2PY's Meson backend by default on Python >= 3.12."""
    if sys.version_info[:2] < (3, 12):
        return
    search_path = str(Path(sys.executable).parent)+os.pathsep+os.environ.get('PATH', '')
    missing = [name for name in ('meson', 'ninja') if shutil.which(name, path=search_path) is None]
    if missing:
        raise RuntimeError('Precise velocity reweighting needs '+', '.join(missing)+
                           ' on PATH for the F2PY build. In the active MadGraph Python environment run:\n'
                           +sys.executable+' -m pip install meson ninja\n'
                           'Then restart MadDM and relaunch. No rate-method fallback was applied.')


@contextmanager
def reweight_build_environment():
    """Bind F2PY and build tools to the Python actually running MadDM."""
    require_reweight_build_tools()
    previous = {key: os.environ.get(key) for key in ('PATH', 'F2PY')}
    try:
        os.environ['PATH'] = str(Path(sys.executable).parent)+os.pathsep+os.environ.get('PATH', '')
        os.environ['F2PY'] = shlex.join([sys.executable, str(Path(__file__).with_name('f2py_compat.py'))])
        yield
    finally:
        for key, value in previous.items():
            if value is None: os.environ.pop(key, None)
            else: os.environ[key] = value


def prepare_lhe(source, destination):
    """Make a Pythia-readable copy; preserve every nonblank event record.

    Pythia's LHE reader consumes exactly NUP particle lines and does not skip
    blank lines in that numeric block. MG5 reweight output can insert one
    immediately after the event header. Never alter the original event file.
    """
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError('Shower LHE copy must differ from the original')
    opener = gzip.open if source.suffix == '.gz' else open
    temporary = destination.with_name(destination.name+'.tmp')
    events = blanks = remaining = 0
    header = inside = False
    try:
        with opener(source, 'rt') as inp, temporary.open('w') as out:
            for number, line in enumerate(inp, 1):
                stripped = line.strip()
                if re.match(r'^<event(?:\s|>)', stripped):
                    if inside:
                        raise ValueError('Nested or unclosed LHE event at line '+str(number))
                    inside, header = True, True
                elif inside and (header or remaining):
                    if not stripped:
                        blanks += 1
                        continue
                    tokens = stripped.split()
                    try:
                        if header:
                            if len(tokens) != 6:
                                raise ValueError('expected 6 event-header fields')
                            remaining = int(tokens[0])
                            int(tokens[1])
                            if remaining < 2 or not all(math.isfinite(float(v)) for v in tokens[2:]):
                                raise ValueError('invalid event-header fields')
                            header = False
                        else:
                            if len(tokens) != 13:
                                raise ValueError('expected 13 particle fields')
                            for token in tokens[:6]: int(token)
                            if not all(math.isfinite(float(v)) for v in tokens[6:]):
                                raise ValueError('non-finite particle fields')
                            remaining -= 1
                    except ValueError as error:
                        raise ValueError('Malformed LHE numeric record at line '+str(number)+': '+str(error)) from error
                elif stripped == '</event>':
                    if not inside or header or remaining:
                        raise ValueError('Incomplete LHE event at line '+str(number))
                    events += 1
                    inside = False
                out.write(line)
        if inside or header or remaining or not events:
            raise ValueError('Shower LHE input has no complete events or is truncated: '+str(source))
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return dict(source=source.name, shower_input=destination.name,
                events=events, removed_numeric_blank_lines=blanks)


def shower_build_state(project, pythia_path):
    """Check the build stamp without modifying the exported C++ source."""
    source = Path(project)/'bin/internal/main101.cc'
    if not source.is_file():
        raise RuntimeError('Missing Pythia spectrum source. Export a new process directory with MadDM.')
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    stamp = source.with_name('main101.build.json')
    desired = dict(source_sha256=source_hash, pythia_path=str(Path(pythia_path).resolve()))
    try:
        current = json.loads(stamp.read_text())
    except (OSError, ValueError):
        current = None
    return source, stamp, desired, current != desired or not source.with_suffix('').is_file()


@lru_cache(maxsize=8)
def lund_settings(pythia_path, data_path=''):
    roots = [Path(data_path)] if data_path else []
    if pythia_path:
        base = Path(pythia_path)
        roots += [base/'share/Pythia8/xmldoc', base/'share/pythia8/xmldoc',
                  base/'xmldoc', base]
    for root in roots:
        if not root.is_dir():
            continue
        # Installed and source distributions both put this setting in XML.
        files = list(root.glob('*Z*.xml')) + list(root.glob('*Fragmentation*.xml'))
        if not files and root.name == 'xmldoc':
            files = list(root.glob('*.xml'))
        text = '\n'.join(p.read_text(errors='replace') for p in files)
        legacy = bool(re.search(r'name\s*=\s*[\"\']StringZ:deriveBLund[\"\']', text))
        modern = bool(re.search(r'name\s*=\s*[\"\']StringZ:deriveLundPars[\"\']', text))
        if legacy or modern:
            return legacy, modern
    raise RuntimeError('Cannot find the installed Pythia Lund settings XML. Check pythia8_path or PYTHIA8DATA.')


def event_file(directory, reshuffling=False):
    directory = Path(directory)
    names = ['unweighted_events.lhe.gz', 'unweighted_events.lhe', 'events.lhe.gz', 'events.lhe']
    if reshuffling:
        names = names[2:]+names[:2]
    for name in names:
        p = directory/name
        if p.is_file() and p.stat().st_size:
            return p
    raise OSError('No nonempty LHE input for the shower in '+str(directory))


def prepare_shower_card(path, lhe_file, pythia_path, data_path=''):
    """Keep the requested Lund parametrization and use the actual event file.

    Modern deriveLundPars=1 is the replacement for deriveBLund=on. No tune
    parameters are dropped and no physics default is substituted.
    """
    path = Path(path)
    lines = path.read_text().splitlines()
    has_legacy, has_modern = lund_settings(str(pythia_path or ''), str(data_path or ''))
    modern_explicit = any(re.match(r'^\s*StringZ:deriveLundPars\s*=', s, re.I) for s in lines)
    out, changes = [], []
    for line in lines:
        match = re.match(r'^\s*([^!=#]+?)\s*=\s*([^!#]+)', line)
        if not match:
            out.append(line)
            continue
        key, value = match.group(1).strip().lower(), match.group(2).strip()
        if key == 'beams:lhef':
            # main101 runs in the event directory; a basename is space-safe.
            out.append('Beams:LHEF = '+Path(lhe_file).name)
        elif key == 'main:numberofevents' and int(value) < 0:
            # MadDM's -1 means all input events; Pythia permits 0, not -1.
            # main101 interprets 0 as read until EOF.
            out.append('Main:NumberOfEvents = 0')
            changes.append('Main:NumberOfEvents='+value+' -> 0 (all LHE events)')
        elif key == 'stringz:deriveblund' and not has_legacy:
            if not has_modern:
                raise RuntimeError('Pythia does not support the requested Lund parametrization')
            if not modern_explicit:
                if value.lower() not in ('on', 'off', 'true', 'false', '1', '0'):
                    raise ValueError('Invalid StringZ:deriveBLund value '+value)
                mode = 1 if value.lower() in ('on', 'true', '1') else 0
                out.append('StringZ:deriveLundPars = '+str(mode))
                changes.append('StringZ:deriveBLund='+value+' -> StringZ:deriveLundPars='+str(mode))
        elif key == 'stringz:derivelundpars' and not has_modern:
            mode = int(value)
            if not has_legacy or mode not in (0, 1):
                raise ValueError('This Pythia cannot represent StringZ:deriveLundPars='+value)
            out.append('StringZ:deriveBLund = '+('on' if mode else 'off'))
        else:
            out.append(line)
    if not any(re.match(r'^\s*Beams:LHEF\s*=', s, re.I) for s in out):
        out.append('Beams:LHEF = '+Path(lhe_file).name)
    path.write_text('\n'.join(out)+'\n')
    return changes


def validate_spectrum(x, values, name='spectrum'):
    import numpy as np
    x, values = np.asarray(x, float), np.asarray(values, float)
    if (x.ndim != 1 or x.size < 2 or values.shape != x.shape
            or np.any(~np.isfinite(x)) or np.any(x <= 0) or np.any(np.diff(x) <= 0)
            or np.any(~np.isfinite(values)) or np.any(values < 0)):
        raise ValueError(name+' must have an ordered positive x grid and finite nonnegative yields; '
                         'check pythia8.log for an aborted shower')
    return x, values


def photon_spectrum(mass, x, dndlogx):
    """dN/dE with log-log interpolation and flux-linear zero transitions."""
    import math
    import numpy as np
    x, values = validate_spectrum(x, dndlogx, 'Gamma spectrum')
    if not math.isfinite(float(mass)) or mass <= 0:
        raise ValueError('Require a finite positive DM mass')
    energy = mass*x
    flux = values/(energy*np.log(10.))
    loge = np.log(energy)

    def spectrum(e):
        if e < energy[0] or e > energy[-1]:
            return 0.
        j = int(np.clip(np.searchsorted(energy, e, side='right')-1, 0, len(energy)-2))
        weight = (np.log(e)-loge[j])/(loge[j+1]-loge[j])
        if flux[j] == 0 or flux[j+1] == 0:
            return float((1-weight)*flux[j]+weight*flux[j+1])
        return float(np.exp((1-weight)*np.log(flux[j])+weight*np.log(flux[j+1])))
    return spectrum

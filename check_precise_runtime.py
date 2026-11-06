#!/usr/bin/env python3
"""Check the active Python/F2PY build tools before a precise MadDM run.

Run with the same Python used to launch MadGraph. --build-check also compiles
and imports a tiny independent Fortran extension, without changing MadDM runs.
"""
import argparse
import importlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import indirect_compat


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-check', action='store_true')
    args = parser.parse_args()
    print('Python:', sys.executable, sys.version.split()[0])
    print('NumPy:', np.__version__)
    for tool in ('f2py', 'meson', 'ninja', 'gfortran'):
        print(tool+':', shutil.which(tool) or 'not found')
    print('MadDM F2PY interpreter:', sys.executable)
    try:
        indirect_compat.require_reweight_build_tools()
        if args.build_check:
            with tempfile.TemporaryDirectory(prefix='maddm_precise_probe_') as directory:
                path = Path(directory)
                (path/'probe.f').write_text('      subroutine double_value(x,y)\n'
                                          '      double precision x,y\n'
                                          'cf2py intent(in) x\n'
                                          'cf2py intent(out) y\n'
                                          '      y=2.d0*x\n      end\n')
                command = [sys.executable, str(Path(__file__).with_name('f2py_compat.py').resolve()), '-c', 'probe.f',
                           '-m', 'maddm_precise_probe']
                with indirect_compat.reweight_build_environment():
                    result = subprocess.run(command, cwd=path, text=True,
                                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                if result.returncode:
                    raise RuntimeError('F2PY build probe failed:\n'+result.stdout[-6000:])
                sys.path.insert(0, directory)
                try:
                    module = importlib.import_module('maddm_precise_probe')
                    if module.double_value(3.) != 6.:
                        raise RuntimeError('F2PY probe returned an unexpected value')
                finally:
                    sys.path.remove(directory)
                    sys.modules.pop('maddm_precise_probe', None)
                print('F2PY Fortran compile/import probe: passed')
        print('Precise runtime prerequisites: passed')
    except (RuntimeError, ImportError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

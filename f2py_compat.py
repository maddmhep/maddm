#!/usr/bin/env python3
"""Run NumPy F2PY in MadDM's Python, without stale make jobserver handles."""
import os


def main():
    # Meson's subprocesses do not inherit make's jobserver file descriptors.
    # An inherited MAKEFLAGS can therefore break recent Ninja versions.
    os.environ.pop('MAKEFLAGS', None)
    os.environ.pop('MFLAGS', None)
    from numpy.f2py.f2py2e import main as f2py_main
    return f2py_main()


if __name__ == '__main__':
    main()

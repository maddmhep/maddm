#!/usr/bin/env python3
"""Repair an existing generated MadDM output; no model parameters are changed.

Run: python /path/to/PLUGIN/maddm/fix_pbar_include.py /path/to/test_pbar_v0
Original changed files are preserved with suffix .before_ams02_pbar.
"""
import argparse
from pbar_interface import repair_generated_pbar_sources

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project', help='Generated MadDM project directory')
    args = parser.parse_args()
    changed = repair_generated_pbar_sources(args.project)
    print('Updated: '+', '.join(changed) if changed else 'Generated antiproton declarations and driver are already current.')
    print('Next, run make -j1 all in the generated project directory.')

#!/usr/bin/env python3
"""Regenerate available pbar, Dbar and GCE PDFs from a completed run directory."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from indirect import pbar_interface
from plotting import signal_plots


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory',type=Path)
    out=parser.parse_args().run_directory
    if (out/'AMS02_pbar_fluxes.dat').is_file():
        print(pbar_interface.save_plot(out))
    if (out/'AMS02_dbar_fluxes.dat').is_file():
        result=json.loads((out/'AMS02_dbar_events.json').read_text())
        result['_flux_columns']=np.loadtxt(out/'AMS02_dbar_fluxes.dat')
        sec=out/'AMS02_dbar_secondary.dat'
        result['_secondary_flux']=np.loadtxt(sec)[:,1] if sec.is_file() else None
        signal_plots.dbar_plot(result,out/'AMS02_dbar_comparison.pdf')
        print(out/'AMS02_dbar_comparison.pdf')
    if (out/'GCE_spectrum.dat').is_file():
        result=json.loads((out/'GCE_fit.json').read_text())
        result['_bins']=np.loadtxt(out/'GCE_spectrum.dat')
        result['_spectrum']=np.loadtxt(out/'GCE_prompt_spectrum.dat')
        signal_plots.gce_plot(result,out/'GCE_fit.pdf')
        print(out/'GCE_fit.pdf')


if __name__=='__main__':main()

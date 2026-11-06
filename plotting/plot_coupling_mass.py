#!/usr/bin/env python3
"""Rescale MadDM bounds to the coupling–mass plane for a one-coupling model.

Example: python plot_coupling_mass.py scan_run_03.txt --coupling 0.1 --dd-power 4 --id-power 4
Use the actual coupling of the scan as --coupling. Omit it for g/g_ref.
Keep this script alongside plot_sigmav_limits.py.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from plot_sigmav_limits import (DD, DD_COLORS, choose, gce_limits, mass_label, numeric,
                                plot_style, positive, read_scan, save_figure, sort_scan, style_axes)


def rescale(reference, prediction, target, power):
    """Solve target = prediction (g/g_ref)^power; skip zero/unavailable response."""
    if prediction is None or target is None:
        return np.full(len(reference), np.nan)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore', under='ignore'):
        return positive(reference * np.exp((np.log(positive(target))-np.log(positive(prediction)))/power))


def generate(input_path, *, outdir=None, mass_column=None, coupling=None, coupling_column=None,
             dd_power=4., id_power=4., relic_power=None, relic_target=3e-26, relic_rate_column='sigmav(xf)',
             fermi_rate_column=None, pbar_rate_column=None, gce_rate_column=None, runs_dir=None,
             gce_upper_column=None, no_gce_in_envelope=False, show=False):
    input_path = Path(input_path)
    mass_column, x, cols = sort_scan(read_scan(input_path), mass_column)
    if len(np.unique(x)) != len(x):
        raise ValueError('Use a scan slice with one row per mass before rescaling to a coupling curve.')
    if coupling is not None and coupling_column is not None:
        raise ValueError('Use either --coupling or --coupling-column.')
    if coupling_column:
        ref = numeric(cols, coupling_column)
        if ref is None:
            raise ValueError('Missing numeric benchmark coupling column: '+coupling_column)
    else:
        ref = np.full(len(x), 1. if coupling is None else coupling)
    if np.any(~np.isfinite(ref) | (ref <= 0)):
        raise ValueError('Benchmark couplings must be positive and finite.')
    normalized = coupling is None and coupling_column is None
    relic_power = id_power if relic_power is None else relic_power
    if any(not np.isfinite(p) or p <= 0 for p in (dd_power, id_power, relic_power, relic_target)):
        raise ValueError('Powers and thermal target must be positive and finite.')
    runs_dir = Path(runs_dir) if runs_dir is not None else input_path.parent
    root = Path(outdir) if outdir is not None else input_path.parent
    root.mkdir(parents=True,exist_ok=True)

    curves, upper_keys = {}, []
    references = {}
    def add(key, rate_key, target, power, is_upper=False):
        y = rescale(ref, numeric(cols, rate_key), target, power)
        curves[key] = y
        references[key] = rate_key
        if is_upper and np.isfinite(y).any():
            upper_keys.append(key)
        return y
    for kind, nucleon in DD:
        rate = 'sigmaN_'+kind+'_'+nucleon
        add(kind+'_'+nucleon, rate, numeric(cols, 'lim_'+rate), dd_power, True)
    fermi_rate = fermi_rate_column or choose(cols, ['tot_SM_xsec','tot_Xsec','taacsID'])
    pbar_rate = pbar_rate_column or choose(cols, ['AMS02_pbar_sigmav_model','AMS02_pbar_model_sigmav'])
    gce_rate = gce_rate_column or choose(cols, ['GCE_model_sigmav'])
    pbar_limit = choose(cols, ['AMS02_pbar_sigmav_UL95','AMS02_pbar_sigmav'])
    add('dwarfs', fermi_rate, numeric(cols, 'Fermi_sigmav'), id_power, True)
    add('pbar', pbar_rate, numeric(cols, pbar_limit), id_power, True)
    _, best, minus, plus, upper, methods, notes = gce_limits(cols,len(x),runs_dir,upper_column=gce_upper_column)
    add('gce_best', gce_rate, best, id_power)
    add('gce_lo68', gce_rate, best-minus, id_power)
    add('gce_hi68', gce_rate, best+plus, id_power)
    add('gce_ul95', gce_rate, upper, id_power, not no_gce_in_envelope)
    add('relic_target', relic_rate_column, relic_target, relic_power)

    envelope = np.full(len(x),np.nan)
    sections = [np.arange(len(x))]
    if upper_keys:
        stack = np.array([curves[k] for k in upper_keys])
        available = np.any(np.isfinite(stack),axis=0)
        envelope[available] = np.nanmin(stack[:,available],axis=0)
        # A limit that ends at the edge of its tabulated coverage must not
        # be interpolated towards the bound from a different experiment.
        owner = np.argmin(np.where(np.isfinite(stack),stack,np.inf),axis=0)
        left, right = np.arange(len(x)-1), np.arange(1,len(x))
        gap = (owner[left] != owner[right]) & (~np.isfinite(stack[owner[left],right]) |
                                              ~np.isfinite(stack[owner[right],left]))
        sections = np.split(np.arange(len(x)), np.flatnonzero(gap)+1)
    if not any(np.isfinite(y).any() for y in curves.values()):
        raise ValueError('No usable limits, fit amplitudes, or freeze-out rates found.')
    curves['upper_envelope'] = envelope
    prefix = root/(input_path.stem+'_coupling_mass')
    table = prefix.with_suffix('.csv')
    with table.open('w',newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['run',mass_column,'benchmark_coupling',*curves,'gce_ul95_method'])
        for i,m in enumerate(x):
            writer.writerow([cols.get('run',['']*len(x))[i],m,'' if normalized else ref[i],
                             *[y[i] if np.isfinite(y[i]) else '' for y in curves.values()],methods[i]])

    plt = plot_style(show)
    fig, ax = plt.subplots(figsize=(9.4,6.4),constrained_layout=True)
    def draw(key,label,color,ls='-',lw=1.6):
        y = curves[key]
        if np.isfinite(y).any():
            ax.plot(x,y,marker='o',ms=2.4,lw=lw,label=label,color=color,ls=ls)
    for (kind,nucleon),color in zip(DD,DD_COLORS):
        draw(kind+'_'+nucleon,f'DD: {kind} {nucleon}',color,'-' if nucleon=='n' else '--')
    draw('dwarfs','Fermi-LAT dwarfs','#654c9c','--')
    draw('pbar','AMS-02 antiprotons','#bd412f','-.')
    draw('gce_best','GCE best fit','#278559','-')
    valid = np.isfinite(curves['gce_lo68']) & np.isfinite(curves['gce_hi68'])
    if valid.any():
        ax.fill_between(x,curves['gce_lo68'],curves['gce_hi68'],where=valid,alpha=.14,color='#278559',
                        label=r'GCE fixed-mass $\Delta\chi^2=1$')
    if np.isfinite(curves['gce_ul95']).any():
        approximate = methods=='quadratic_approx'
        curves['gce_exact_plot'] = np.where(~approximate,curves['gce_ul95'],np.nan)
        curves['gce_approx_plot'] = np.where(approximate,curves['gce_ul95'],np.nan)
        draw('gce_exact_plot','GCE upper limit','#278559','--')
        draw('gce_approx_plot','GCE upper limit (approx.)','#278559',':')
        del curves['gce_exact_plot'],curves['gce_approx_plot']
    draw('relic_target',r'Thermal target at freeze-out','#65574e','--',1.9)
    if np.isfinite(envelope).any():
        first = True
        for segment in sections:
            if np.isfinite(envelope[segment]).any():
                ax.plot(x[segment], envelope[segment], 'o-', ms=2.4, lw=2.3, color='#181b20',
                        label='Strongest available upper bound' if first else None)
                first = False
        style_axes(ax,mass_label(mass_column),'')
        ax.autoscale_view()
        ymin,ymax = ax.get_ylim()
        for segment in sections:
            ax.fill_between(x[segment],envelope[segment],ymax,where=np.isfinite(envelope[segment]),
                            color='#445264',alpha=.08,zorder=0)
        ax.set_ylim(ymin,ymax)
    xlabel = mass_label(mass_column) + '\nIndividual bounds; freeze-out target approximate near resonances and thresholds.'
    style_axes(ax,xlabel,r'Coupling ratio $g/g_{\rm ref}$' if normalized else r'Coupling $g$')
    ax.set_title('Coupling constraints from a fixed-shape rescaling')
    ax.legend(frameon=False,fontsize=8.5,ncol=2,loc='best')
    save_figure(plt,fig,prefix.with_suffix('.pdf'),show)
    meta = {'reference_coupling':coupling,'coupling_column':coupling_column,'axis_is_g_over_gref':normalized,
            'powers':{'direct':dd_power,'indirect':id_power,'freezeout':relic_power},
            'rate_columns':references,'envelope_inputs':upper_keys,'relic_target_cm3_s':relic_target,
            'GCE_upper_limit_methods':{k:int(np.sum(methods==k)) for k in sorted(set(methods)-{''})},
            'notes':notes+['Envelope = strongest available individual constraint, no joint likelihood.',
                           'Full DM density assumed; spectra, branching fractions and halo settings fixed at each mass.',
                           'Thermal target is an approximate freeze-out criterion, unreliable near thresholds, resonances or with coannihilation.']}
    prefix.with_suffix('.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
    return curves,meta


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path)
    p.add_argument('--outdir',type=Path)
    p.add_argument('--mass-column')
    group = p.add_mutually_exclusive_group()
    group.add_argument('--coupling','--g-ref',type=float,help='Benchmark coupling used to produce the scan')
    group.add_argument('--coupling-column',help='Use a per-mass benchmark coupling column instead')
    for name in ('dd','id','relic'):
        p.add_argument('--'+name+'-power',type=float,default=None if name=='relic' else 4.)
    p.add_argument('--relic-target',type=float,default=3e-26,help='Thermal target [cm^3/s]')
    p.add_argument('--relic-rate-column',default='sigmav(xf)')
    for name in ('fermi','pbar','gce'):
        p.add_argument('--'+name+'-rate-column')
    p.add_argument('--gce-upper-column')
    p.add_argument('--runs-dir',type=Path)
    p.add_argument('--no-gce-in-envelope',action='store_true',help='Show GCE separately from the exclusion envelope')
    p.add_argument('--show',action='store_true')
    try:
        args = vars(p.parse_args())
        generate(args.pop('input'),**args)
    except (ValueError,OSError,ImportError) as error:
        p.error(str(error))


if __name__ == '__main__':
    main()

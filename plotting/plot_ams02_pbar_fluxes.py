#!/usr/bin/env python3
"""Plot the named columns of MadDM's AMS02_pbar_fluxes.dat as a vector PDF.

Usage: python plot_ams02_pbar_fluxes.py path/to/AMS02_pbar_fluxes.dat -o pbar.pdf
Dependencies: numpy, matplotlib. No MadDM import, fitting, or interpolation.
The lower panel uses the signed fractional residual (model-data)/model.
"""
import argparse
from pathlib import Path
import numpy as np


PREDICTIONS = {'model': 'model', 'thermal': 'thermal', 'best': 'best', 'upper': 'upper'}


def quiet_plot_logging():
    """Keep third-party warnings/errors; omit rendering DEBUG/INFO chatter.

    MG5's coloring formatter can mutate %-format DEBUG records, causing
    logging tracebacks while fontTools subsets PDF fonts.
    """
    import logging
    for root in ('matplotlib', 'fontTools'):
        logging.getLogger(root).setLevel(logging.WARNING)
        for name, logger in list(logging.Logger.manager.loggerDict.items()):
            if name.startswith(root+'.') and isinstance(logger,logging.Logger):
                logger.setLevel(logging.WARNING)


def read_fluxes(path, prediction='model'):
    path = Path(path)
    names = None
    for line in path.read_text().splitlines():
        if line.lstrip().startswith('#') and 'R_GV' in line and 'data_per_GV' in line:
            names = line.lstrip().lstrip('#').split()
            break
    if names is None:
        raise ValueError('Missing named column header (R_GV, data_per_GV, error_per_GV, ...)')
    a = np.loadtxt(path, ndmin=2)
    if a.shape[1] != len(names):
        raise ValueError('Numeric columns do not match the header')
    prefix = PREDICTIONS[prediction]
    needed = ['R_GV', 'data_per_GV', 'error_per_GV']+[prefix+'_'+x for x in ('sec', 'dm', 'total')]
    if any(name not in names for name in needed):
        raise ValueError('Missing columns for prediction '+prediction)
    r, data, error, sec, dm, total = [a[:, names.index(name)] for name in needed]
    if (np.any(~np.isfinite(np.array([r, data, error, sec, dm, total])))
            or len(r) < 2 or np.any(r <= 0) or np.any(np.diff(r) <= 0)
            or np.any(data <= 0) or np.any(error <= 0) or np.any(sec < 0)
            or np.any(dm < 0) or np.any(total <= 0)):
        raise ValueError('Require ordered positive rigidities, positive data/errors/total, and finite nonnegative components')
    if not np.allclose(sec+dm, total, rtol=5.e-8, atol=0.):
        raise ValueError('The selected secondary plus DM flux does not equal its total')
    return dict(R=r, data=data, error=error, secondary=sec, dm=dm, total=total,
                residual=(total-data)/total, residual_error=error/total)


def make_plot(flux, output, prediction='model', title=None, residual_scale='linear', linthresh=.02,
              model_dm=None, setup=None):
    quiet_plot_logging()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator, NullFormatter

    settings = {'font.family': 'DejaVu Sans', 'font.size': 11,
                'axes.labelsize': 12, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
                'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 600}
    with plt.rc_context(settings):
        fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(7.2, 7.4),
                  gridspec_kw={'height_ratios': [3., 1.25], 'hspace': .06}, layout='constrained')
        r, scale = flux['R'], flux['R']**3
        top.errorbar(r, scale*flux['data'], yerr=scale*flux['error'],
                     fmt='o', ms=3.6, color='#20242b', ecolor='#6e747c', elinewidth=.85,
                     capsize=1.8, label='AMS-02', zorder=5)
        tag=' (UL fit)' if prediction=='upper' else ''
        for key, color, style, label in [('secondary', '#dd8a19', '--', 'Secondary'+tag),
                                         ('dm', '#a54884', ':', 'Dark matter'+tag),
                                         ('total', '#166eae', '-', 'Total'+tag)]:
            y = scale*flux[key]
            # Exact zero DM bins stay absent on a logarithmic flux axis.
            top.plot(r, np.where(y > 0, y, np.nan), color=color, ls=style, lw=2., label=label)
        if model_dm is not None:
            values=scale*np.asarray(model_dm)
            top.plot(r,np.where(values>0,values,np.nan),color='#278455',ls='-.',lw=1.9,
                     label='DM at model rate (model-fit modulation)')
        top.set_xscale('log', base=10)
        top.set_yscale('log', base=10)
        # Show the measured flux range, rather than letting a vanishing DM
        # tail stretch the axis over tens of decades. All data are retained.
        envelope = scale*np.minimum(flux['data'], flux['total'])
        ceiling = scale*np.maximum(flux['data']+flux['error'], flux['total'])
        top.set_ylim(.2*np.min(envelope), 3.*np.max(ceiling))
        top.set_ylabel(r'$R^3 J_{\bar p}^{\mathrm{TOA}}$ [GV$^2$ m$^{-2}$ s$^{-1}$ sr$^{-1}$]')
        top.legend(frameon=False, loc='upper left', ncol=2, fontsize=8.5)
        if title:
            top.set_title(title, fontsize=11, pad=12)
        top.text(.98, .03, 'Prediction: '+prediction+('\n'+setup if setup else ''), transform=top.transAxes,
                 ha='right', va='bottom', fontsize=9, color='#505660')
        bottom.errorbar(r, flux['residual'], yerr=flux['residual_error'],
                        fmt='o', ms=3.5, color='#166eae', ecolor='#6193b7', elinewidth=.9,
                        capsize=1.8, zorder=4)
        bottom.axhline(0., lw=1., ls='--', color='#606974')
        if residual_scale == 'symlog':
            bottom.set_yscale('symlog', base=10, linthresh=linthresh, linscale=.7)
            bottom.text(.98, .04, 'Symmetric log; linear within '+f'{linthresh:g}',
                        ha='right', transform=bottom.transAxes, fontsize=8, color='#606974')
        else:
            bottom.set_yscale('linear')
        bottom.set_ylabel(r'$(J_{\mathrm{model}}-J_{\mathrm{data}})/J_{\mathrm{model}}$', fontsize=10)
        bottom.set_xlabel(r'Rigidity $R$ [GV]')
        extent = max(np.max(np.abs(flux['residual'])+flux['residual_error'])*1.35, 2*linthresh)
        bottom.set_ylim(-extent, extent)
        # 1 TV = 1000 GV (about 1 TeV for highly relativistic antiprotons).
        # Retain all data should the supplied spectrum already exceed 1 TV.
        bottom.set_xlim(r[0]/1.08, max(1000., r[-1]*1.08))
        for ax in (top, bottom):
            ax.grid(True, which='major', color='#d6dce3', lw=.6, alpha=.75)
            ax.grid(True, which='minor', axis='x', color='#e6e9ed', lw=.4, alpha=.5)
            ax.tick_params(which='both', direction='in', top=True, right=True)
            ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
            ax.xaxis.set_minor_formatter(NullFormatter())
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, format='pdf', dpi=600, metadata={'Title': 'AMS-02 antiproton flux comparison'})
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('-o', '--output', type=Path, help='PDF path; default next to input')
    parser.add_argument('--prediction', choices=PREDICTIONS, default='upper',
                        help='Choose the matching secondary/DM/total columns; default: upper, plus model DM')
    parser.add_argument('--residual-scale', choices=['symlog', 'linear'], default='linear')
    parser.add_argument('--linthresh', type=float, default=.02,
                        help='Linear half-width of the signed symmetric-log residual axis')
    parser.add_argument('--title')
    args = parser.parse_args()
    if args.linthresh <= 0:
        parser.error('--linthresh must be positive')
    output = args.output or args.input.with_name('AMS02_pbar_comparison.pdf')
    if output.suffix.lower() != '.pdf':
        parser.error('--output must end with .pdf')
    try:
        values = read_fluxes(args.input, args.prediction)
        model_dm = read_fluxes(args.input,'model')['dm'] if args.prediction=='upper' else None
        make_plot(values, output, args.prediction, args.title, args.residual_scale, args.linthresh,model_dm=model_dm)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    print(output)


if __name__ == '__main__':
    main()

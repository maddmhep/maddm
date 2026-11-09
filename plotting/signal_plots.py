"""Noninteractive, vector PDF figures for MadDM run products."""
from pathlib import Path
import numpy as np


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


def axes():
    quiet_plot_logging()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'savefig.dpi':600})
    fig,ax=plt.subplots(figsize=(8.,6.),constrained_layout=True)
    ax.set_xscale('log');ax.set_yscale('log')
    ax.grid(which='major',alpha=.35);ax.grid(which='minor',alpha=.12)
    ax.tick_params(which='both',direction='in',top=True,right=True)
    return plt,fig,ax


def _plot_yields(curves, output, *, xlabel, ylabel, title):
    """Use gaps for zero yields rather than drawing across them on log axes."""
    valid = []
    for x, y, label, color, style in curves:
        mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if np.any(mask):
            isolated = np.flatnonzero(mask & ~np.r_[False, mask[:-1]] & ~np.r_[mask[1:], False]).tolist()
            valid.append((np.where(mask, x, np.nan), np.where(mask, y, np.nan),
                          label, color, style, isolated))
    if not valid:
        return False
    plt, fig, ax = axes()
    try:
        for x, y, label, color, style, isolated in valid:
            ax.plot(x, y, label=label, color=color, ls=style, lw=1.9,
                    marker='o' if isolated else None, markersize=3.5, markevery=isolated)
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        ax.legend(loc='best', frameon=False)
        fig.savefig(output, format='pdf', dpi=600)
    finally:
        plt.close(fig)
    return True


def recoil_rates_plot(data, output):
    """Plot the RAPIDD xenon, argon and germanium rates as written to disk."""
    energy, xe, ar, ge = np.loadtxt(data, ndmin=2, usecols=(0, 1, 2, 3)).T
    return _plot_yields([(energy, xe, 'Xenon', '#1c6ba0', '-'),
                         (energy, ar, 'Argon', '#c2542e', '--'),
                         (energy, ge, 'Germanium', '#27824d', '-.')], output,
                        xlabel=r'Recoil energy $E_R$ [keV]',
                        ylabel=r'$dR/dE_R$ [events kg$^{-1}$ day$^{-1}$ keV$^{-1}$]',
                        title='Nuclear recoil spectra')


def source_spectra_plot(directory, method, output):
    """Plot only the source yields supplied by this run and spectrum backend."""
    species = [
        ('gammas', r'$\gamma$', '#1c6ba0', '-'),
        ('neutrinos_e', r'$\nu_e$', '#7652a2', '--'),
        ('neutrinos_mu', r'$\nu_\mu$', '#7652a2', '-.'),
        ('neutrinos_tau', r'$\nu_\tau$', '#7652a2', ':'),
        ('antiprotons', r'$\bar p$', '#c2542e', '-'),
        ('positrons', r'$e^+$', '#27824d', '-'),
        ('antideuterons', r'$\bar D$', '#b1368a', '-'),
    ]
    curves = []
    for key, label, color, style in species:
        path = Path(directory) / (key + '_spectrum_' + method + '.dat')
        if not path.is_file():
            continue
        logx, yield_per_logx = np.loadtxt(path, usecols=(0, 1), ndmin=2).T
        curves.append((10.**logx, yield_per_logx, label, color, style))
    return _plot_yields(curves, output, xlabel=r'$x = E_{\rm kin}/m_{\rm DM}$',
                        ylabel=r'$dN/d\log_{10}x$ per annihilation',
                        title='Particle yields at production (%s)' % method)


def dbar_plot(result, output):
    plt,fig,ax=axes()
    a=result['_flux_columns'];t=a[:,0]
    line=np.where(a[:,2]>0,a[:,2],np.nan)
    ax.plot(t,line,color='#b1368a',lw=2.3,label='DM at model rate')
    if result['xsi']<1:
        ax.plot(t,np.where(a[:,4]>0,a[:,4],np.nan),color='#b1368a',lw=1.8,ls=':',label='DM, thermal density')
    sec=result.get('_secondary_flux')
    if sec is not None:ax.plot(t,sec,'--',color='#444444',lw=2,label='Secondary reference')
    scales=[2.e-6,8.e-7,2.e-7,float(np.nanmax(line)) if np.any(np.isfinite(line)) else 1.e-12]
    if sec is not None:scales.append(float(np.max(sec)))
    ymax=max(scales)*10.
    signal_peak=float(np.nanmax(line)) if np.any(np.isfinite(line)) else 1.e-12
    ymin=max(1.e-20,min(2.e-9,signal_peak*1.e-3))
    for name,color,window in [('GAPS','#9b82d7',(.1,.3,2.e-6)),('AMS-02 low','#79a6ec',(.18,.7,8.e-7)),('AMS-02 high','#79a6ec',(2.,4.2,2.e-7))]:
        lo,hi,sens=window
        ax.fill_between([lo,hi],sens,ymax,color=color,alpha=.32,zorder=0)
        ax.plot([lo,hi],[sens,sens],color=color,lw=1.5)
        ax.text(np.sqrt(lo*hi),sens*1.3,name,ha='center',fontsize=9)
    events=result['cases']['all_dm']
    text=('Model DM signal events\nAMS-02 low: %.3g\nAMS-02 high: %.3g\nAMS-02 sum: %.3g\nGAPS: %.3g' %
          (events['windows']['AMS02_low']['expected_events'],events['windows']['AMS02_high']['expected_events'],events['AMS02_total_events'],events['windows']['GAPS']['expected_events']))
    text+='\nRate: %.3g cm$^3$/s\nFisk potential: %.3g GV'%(result['total_sigmav_cm3s'],events['phi_GV'])
    if result['xsi']<1:
        therm=result['cases']['thermal']
        text+='\nThermal AMS-02 sum: %.3g\nThermal GAPS: %.3g'%(therm['AMS02_total_events'],therm['windows']['GAPS']['expected_events'])
    ax.text(.97,.61,text,transform=ax.transAxes,ha='right',va='top',fontsize=10,
            bbox=dict(facecolor='white',alpha=.9,edgecolor='#dddddd'))
    ax.set_xlim(.1,100.);ax.set_ylim(ymin,ymax)
    ax.set_xlabel(r'Kinetic energy per nucleon $T$ [GeV/n]')
    ax.set_ylabel(r'$\Phi_{\bar D}^{\mathrm{TOA}}$ [m$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV/n)$^{-1}$]')
    ax.set_title(r'$m_{\mathrm{DM}}=%.4g$ GeV; '%result['mass_GeV']+result['propagation']+' / '+result['halo'])
    ax.legend(loc='lower left',frameon=False)
    fig.savefig(output,format='pdf',dpi=600);plt.close(fig)


def gce_plot(result, output):
    plt,fig,ax=axes()
    a=result['_bins'];used=a[:,9]>0
    # Publish the full systematic envelope, including censored bins, independently
    # of the statistical error bars used only for plotting the detections.
    ax.fill_between(a[:,0],np.maximum(a[:,4],1.e-30),np.maximum(a[:,5],1.e-30),
                    color='#777777',alpha=.38,label='GCE systematic envelope',zorder=1)
    stat=a[:,10] if a.shape[1]>10 else None
    ax.errorbar(a[used,0],a[used,3],yerr=stat[used] if stat is not None else None,
                xerr=np.array([a[used,0]-a[used,1],a[used,2]-a[used,0]]),
                fmt='s',ms=3.8,color='black',elinewidth=1.1,label='GCE data, statistical errors',zorder=5)
    if np.any(~used):
        # Use the actual published statistical upper limit; the grey envelope
        # remains the separate systematic information used by the checks.
        upper=a[~used,6]
        ax.errorbar(a[~used,0],upper,yerr=upper*.28,uplims=True,fmt='v',ms=3,
                    color='black',label='Published upper limits (not fitted)',zorder=4)
    s=result['_spectrum']
    ax.plot(s[:,0],np.where(s[:,2]>0,s[:,2],np.nan),color='#1c6ba0',lw=2.2,label='Best-fit prompt DM')
    ax.fill_between(s[:,0],np.maximum(s[:,3],1.e-30),np.maximum(s[:,4],1.e-30),color='#1c6ba0',alpha=.18,label=r'$\Delta\chi^2=1$ amplitude interval')
    ax.plot(s[:,0],np.where(s[:,1]>0,s[:,1],np.nan),color='#c2542e',ls='--',lw=1.8,label='DM at model rate')
    ax.text(.97,.94,'Best-fit rate: %.3g cm$^3$/s\nError: -%.2g / +%.2g\n$\\chi^2$/dof = %.3g / %d\n$\\overline{J}$ = %.3g'%
            (result['sigmav_best_cm3s'],result['error_minus_cm3s'],result['error_plus_cm3s'],result['chi2'],result['ndof'],result['geometry']['jbar']),
            transform=ax.transAxes,ha='right',va='top',fontsize=10)
    peaks=max(float(np.max(s[:,1:3])),float(np.max(a[used,5])))
    ax.set_xlim(.5,1000.);ax.set_ylim(1.e-10,peaks*5)
    ax.set_xlabel(r'Photon energy $E_\gamma$ [GeV]')
    ax.set_ylabel(r'$E_\gamma^2 d\Phi/dE_\gamma$ [GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]')
    ax.set_title('Galactic-center excess: prompt spectral fit')
    ax.legend(loc='lower left',frameon=False,fontsize=9)
    fig.savefig(output,format='pdf',dpi=600);plt.close(fig)

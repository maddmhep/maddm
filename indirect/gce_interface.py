"""Prompt GCE amplitude fit using the supplied SED and dimensionless Jbar.

Measured bins use asymmetric systematic-band errors. Sentinel -1 rows are
upper limits, never negative measurements: they are checked but not assigned
a fabricated Gaussian likelihood. No ICS or bremsstrahlung is implied.
"""
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

from . import pbar_interface, indirect_compat, cr_config

ROOT = Path(__file__).resolve().parents[1]/'ExpData/GCE'
KPC_CM = 3.085677581491367e21
REF = 1.e-26
NUMERIC_KEYS = ('GCE_sigmav_best','GCE_sigmav_error_minus','GCE_sigmav_error_plus',
                'GCE_sigmav_error','GCE_chi2','GCE_ndof','GCE_chi2_ndof','GCE_model_sigmav')


def read_data(path=ROOT/'GCE_flux_data.txt'):
    lines = Path(path).read_text().splitlines()
    selected, rows = False, []
    for line in lines:
        if '###### 40x40deg^2 ROI' in line:
            if selected: break
            selected = True
            continue
        if selected:
            try: values = list(map(float,line.split()))
            except ValueError: continue
            if len(values)==8: rows.append(values)
    a=np.array(rows,float)
    if a.shape!=(26,8) or np.any(np.diff(a[:,0])<=0):
        raise ValueError('Expected 26 ordered bins in the first standard 40x40 GCE block')
    # SED file: MeV cm^-2 s^-1 integrated over the ROI; convert to GeV/sr.
    energy=a[:,:3]*1.e-3
    flux=a[:,3:]*1.e-3/.4775494185
    used=a[:,3]>=0
    if np.any(flux[used,0]-flux[used,3]<=0) or np.any(flux[used,4]-flux[used,0]<=0):
        raise ValueError('Systematic bounds must bracket every measured flux')
    return energy, flux, used


def fit_amplitude(pred_ref, observed, lower_error, upper_error):
    """Convex split-normal least squares; physical amplitude >=0, Delta chi2=1."""
    p,y,lo,hi=map(lambda a:np.asarray(a,float),(pred_ref,observed,lower_error,upper_error))
    if p.shape!=y.shape or lo.shape!=y.shape or hi.shape!=y.shape or not np.all(np.isfinite([p,y,lo,hi])) or np.any(p<0) or np.any(lo<=0) or np.any(hi<=0) or not np.any(p>0):
        raise ValueError('Invalid GCE template or errors, or no photons in measured bins')
    def chi(a):
        diff=a*p-y
        return float(np.sum((diff/np.where(diff>=0,hi,lo))**2))
    def derivative(a):
        diff=a*p-y
        return float(np.sum(p*diff/np.where(diff>=0,hi,lo)**2))
    right=max(1.,2*float(np.max(y[p>0]/p[p>0])))
    best=0. if derivative(0.)>=0 else brentq(derivative,0.,right,xtol=1.e-12)
    minimum=chi(best)
    low=0. if chi(0.)<=minimum+1. else brentq(lambda a:chi(a)-minimum-1.,0.,best,xtol=1.e-12)
    right=max(right,2*best,1.)
    while chi(right)<minimum+1.: right*=2
    high=brentq(lambda a:chi(a)-minimum-1.,best,right,xtol=1.e-12)
    return dict(amplitude=best,amplitude_low=low,amplitude_high=high,chi2=minimum,
                ndof=int(len(y)-1),lower_interval_at_boundary=low==0.)


class GCEFit:
    def __init__(self):
        with np.load(ROOT/'CosmiXs_gamma.npz',allow_pickle=False) as data:
            self.masses=data['masses']; self.x=10.**data['log10x']
            self.channels=data['channels'].tolist(); self.yields=data['dndlog10x']

    def spectrum(self, mass, branching):
        if not self.masses[0]<=mass<=self.masses[-1]: raise ValueError('Mass outside CosmiXs gamma grid')
        j=int(np.clip(np.searchsorted(self.masses,mass,side='right')-1,0,len(self.masses)-2))
        w=(mass-self.masses[j])/(self.masses[j+1]-self.masses[j])
        grid=(1-w)*self.yields[j]+w*self.yields[j+1]
        values=sum(float(v)*grid[:,self.channels.index(ch)] for ch,v in branching.items())
        f=indirect_compat.photon_spectrum(mass,self.x,values)
        return lambda e: np.array([f(float(v)) for v in np.atleast_1d(e)])

    def evaluate(self, mass, records, dm_pdg, md_factor, card, runtime_spectrum=None):
        mixture=pbar_interface.channel_mixture(records,dm_pdg,md_factor)
        if mixture['status']=='ok':
            spectrum=self.spectrum(mass,mixture['branching'])
            source='CosmiXs prompt gamma yields, mixture from Milky Way rates'
            rate=mixture['total_sigmav_cm3s']
        elif runtime_spectrum is not None:
            spectrum=runtime_spectrum; rate=sum(r['sigmav_cm3s'] for r in records)
            source='runtime MadDM photon spectrum at the same requested velocity'
        else:
            return dict(status=mixture['status'],reason='No supported two-body mixture; an actual runtime photon spectrum at the GC velocity is required for other final states',rates=records)
        get=pbar_interface.card_value
        jbar=float(get(card,'gce_jbar',120.));rho=cr_config.local_density(card);rsun=float(get(card,'gce_rsun_kpc',8.122))
        if md_factor not in (4,8) or any(not math.isfinite(v) or v<=0 for v in [jbar,rho,rsun,mass]): raise ValueError('Invalid GCE geometry or DM convention')
        # Eq.37 is Dirac (16 pi); self-conjugate models require 8 pi.
        prefactor=rsun*KPC_CM*rho**2*jbar/(2*math.pi*md_factor*mass**2)
        energies,flux,used=read_data()
        eref=energies[:,0]
        prediction=prefactor*REF*eref**2*spectrum(eref)
        y=flux[:,0]
        fit=fit_amplitude(prediction[used],y[used],y[used]-flux[used,3],flux[used,4]-y[used])
        best,low,high=[fit[k]*REF for k in ('amplitude','amplitude_low','amplitude_high')]
        dense=np.geomspace(.5,1000.,350)
        dense_ref=prefactor*REF*dense**2*spectrum(dense)
        compatibility=[dict(E_GeV=float(e),sys_upper=float(upper),quoted_upper=float(ul),
                             best_fit_SED=float(pred*fit['amplitude']),exceeds_sys_upper=bool(pred*fit['amplitude']>upper))
                       for e,pred,upper,ul in zip(eref[~used],prediction[~used],flux[~used,4],flux[~used,2])]
        return dict(status='ok',mass_GeV=float(mass),sigmav_best_cm3s=best,error_minus_cm3s=best-low,
                    error_plus_cm3s=high-best,interval68_cm3s=[low,high],chi2=fit['chi2'],ndof=fit['ndof'],
                    chi2_ndof=fit['chi2']/fit['ndof'],measurements=int(used.sum()),upper_limit_bins=int((~used).sum()),
                    geometry=dict(jbar=jbar,rho_sun_GeVcm3=rho,rsun_kpc=rsun,solid_angle_sr=.4775494185,
                                  reference_rho_sun_GeVcm3=.345,density_rescaling_from_reference=(rho/.345)**2,
                                  denominator_pi=2*md_factor,reference='arXiv:2508.02775 Eqs.37-39/Table I'),
                    source=source,branching=mixture.get('branching'),rates=records,model_sigmav_cm3s=rate,
                    upper_limit_checks=compatibility,lower_interval_at_boundary=fit['lower_interval_at_boundary'],
                    assumptions=['Prompt photons only; no ICS or bremsstrahlung.',
                                 'Systematic envelope treated as asymmetric diagonal 1-sigma errors; no covariance supplied.',
                                 'Upper-limit rows checked separately, not fitted as negative/zero measurements.',
                                 'Mass and branching shape fixed; only nonnegative annihilation amplitude fitted.',
                                 'Jbar=%g is selected independently of the charged-CR halo table.' % jbar],
                    _bins=np.c_[energies,y,flux[:,3],flux[:,4],flux[:,2],prediction*rate/REF,prediction*best/REF,used.astype(int),flux[:,1]],
                    _spectrum=np.c_[dense,dense_ref*rate/REF,dense_ref*best/REF,dense_ref*low/REF,dense_ref*high/REF])


def observables(result):
    values=dict.fromkeys(NUMERIC_KEYS,-1.);values['GCE_status']=result.get('status','not_evaluated')
    if result.get('status')=='ok':
        values.update(GCE_sigmav_best=result['sigmav_best_cm3s'],GCE_sigmav_error_minus=result['error_minus_cm3s'],
            GCE_sigmav_error_plus=result['error_plus_cm3s'],GCE_sigmav_error=.5*(result['error_minus_cm3s']+result['error_plus_cm3s']),
            GCE_chi2=result['chi2'],GCE_ndof=result['ndof'],GCE_chi2_ndof=result['chi2_ndof'],GCE_model_sigmav=result['model_sigmav_cm3s'])
    return values


def summary(out,result):
    out.write('\n# Galactic-center excess: prompt spectrum, systematic-band fit\n')
    for k,v in observables(result).items(): out.write('%s = %s\n' % (k,v if isinstance(v,str) else '%.12e'%v))
    if result.get('status')=='ok': out.write('# Upper-limit rows checked separately; fitted detections=%d\n'%result['measurements'])


def write_outputs(out_dir,result,plot=True,*,detailed=True):
    out=Path(out_dir);out.mkdir(exist_ok=True,parents=True)
    for name in ['GCE_spectrum.dat','GCE_prompt_spectrum.dat','GCE_fit.pdf']:(out/name).unlink(missing_ok=True)
    (out/'GCE_fit.json').write_text(json.dumps({k:v for k,v in result.items() if not k.startswith('_')},indent=2,allow_nan=False)+'\n')
    if not detailed:
        (out/'GCE_fit.txt').unlink(missing_ok=True)
        return
    with (out/'GCE_fit.txt').open('w') as stream:summary(stream,result)
    if result.get('status')!='ok':return
    np.savetxt(out/'GCE_spectrum.dat',result['_bins'],header='Eref_GeV Emin_GeV Emax_GeV data_SED sys_low sys_high quoted_UL model_SED best_fit_SED used_in_fit stat_error_SED; SED GeV cm^-2 s^-1 sr^-1')
    np.savetxt(out/'GCE_prompt_spectrum.dat',result['_spectrum'],header='E_GeV model_SED best_fit_SED low68_SED high68_SED; SED GeV cm^-2 s^-1 sr^-1')
    if plot:
        if "." in __package__:
            from ..plotting import signal_plots
        else:
            from plotting import signal_plots
        signal_plots.gce_plot(result,out/'GCE_fit.pdf')

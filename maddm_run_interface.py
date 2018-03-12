from __future__ import division
import collections
import math
import logging
import os
import re
import sys
import subprocess
import auxiliary as aux
import timeit
import stat
import shutil
import threading, subprocess
import json

# FF check if these are already imported somewhere else?                                                                                                                            
import scipy
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.optimize import brute
from scipy.optimize import fmin
from scipy.special import gammainc

import shutil
import MGoutput
from madgraph import MadGraph5Error
from models import check_param_card

MG5MODE = True
import madgraph
import madgraph.madevent.sum_html as sum_html
import madgraph.various.misc as misc
import madgraph.interface.extended_cmd as cmd
import madgraph.various.banner as banner_mod
import madgraph.interface.common_run_interface as common_run
import madgraph.interface.madevent_interface as me5_interface
import madgraph.iolibs.files as files
import models.check_param_card as param_card_mod
        
#import darkmatter as darkmatter

import numpy as np

try:
    import pymultinest
except:
    print('WARNING: Multinest module not found! All multinest parameter scanning features will be disabled.')


#import types

pjoin = os.path.join
logger = logging.getLogger('madgraph.plugin.maddm')
#logger.setLevel(10) #level 20 = INFO

MDMDIR = os.path.dirname(os.path.realpath( __file__ ))
PPPCDIR = os.getcwd()+'/PPPC4DMID/'



os.system('rm /home/ucl/cp3/fambrogi/NEW_MADDM/MadDM_MARCH/maddm_dev2/TUTTO/Indirect/RunWeb')

#Is there a better definition of infinity?
__infty__ = float('inf')
__mnestlog0__ = -1.0E90

class ExpConstraints:

    def __init__(self):

        self._allowed_final_states = ['qqx', 'ccx', 'gg', 'bbx', 'ttx', 'e+e-', 'mu+mu-', 'ta+ta-', 'w+w-', 'zz', 'hh', 
                                      'hess2013','hess2016', 'aaER16','aaIR90','aaNFWcR3','aaNFWR41']

        self._oh2_planck = 0.1198
        self._oh2_planck_width = 0.0015

        self._dd_si_limit_file = pjoin(MDMDIR, 'ExpData', 'Xenont1T_data_2017.dat')
        self._dd_sd_proton_limit_file = pjoin(MDMDIR, 'ExpData', 'Pico60_sd_proton.dat') # <---------CHANGE THE FILE!!!
        self._dd_sd_neutron_limit_file = pjoin(MDMDIR, 'ExpData', 'Lux_2017_sd_neutron.dat')
        self._id_limit_file = {'qqx'   :pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_qq.dat'),
                               'ccx'   :pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_cc.dat'),
                               'gg'    :pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_gg.dat'),
                               'bbx'   :pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_bb.dat'),
                               'ttx'   :pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_tt.dat'),
                               'e+e-'  :pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_ee.dat'),
                               'mu+mu-':pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_mumu.dat'),
                               'ta+ta-':pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_tautau.dat'),
                               'w+w-'  :pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_ww.dat'),
                               'zz'    :pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_zz.dat'),
                               'hh'    :pjoin(MDMDIR, 'ExpData', 'MadDM_FermiLim_hh.dat'),

                               'hess2013': pjoin(MDMDIR,'ExpData', 'hess_I_2013_einasto.dat'),
                               'hess2016': pjoin(MDMDIR,'ExpData', 'hess_2016_einasto.dat'),

                               'aaER16':pjoin(MDMDIR,'ExpData', 'Fermi_lines_2015_Einasto_R16.dat'),
                               'aaIR90':pjoin(MDMDIR,'ExpData', 'Fermi_lines_2015_Isothermal_R90.dat'),
                               'aaNFWcR3':pjoin(MDMDIR,'ExpData', 'Fermi_lines_2015_NFWcontracted_R3.dat'),
                               'aaNFWR41':pjoin(MDMDIR,'ExpData', 'Fermi_lines_2015_NFW_R41.dat')}

        self._id_limit_vel = {'qqx':2.0E-5,   'ccx':2.0E-5, 'gg':2.0E-5, 'bbx':2.0E-5,'ttx':2.0E-5,'e+e-':2.0E-5,'mu+mu-':2.0E-5,'ta+ta-':2.0E-5,
                              'w+w-':2.0E-5,  'zz':2.0E-5,  'hh':2.0E-5,
                              'aaER16':1.0E-3,'aaIR90':1.0E-3,'aaNFWcR3':1.0E-3,'aaNFWR41':1.0E-3 ,
                              'hess2013': 1.0E-3 , 'hess2016': 1.0E-3 } 

        self._id_limit_mdm = dict()
        self._id_limit_sigv = dict()

        #In case there is a measurement of the cross seciton
        self._sigma_SI = -1.0
        self._sigma_SI_width = -1.0
        self._sigma_SDp = -1.0
        self._sigma_SDn = -1.0
        self._sigma_SDp_width = -1.0
        self._sigma_SDn_width = -1.0

        self._sigma_ID = dict()
        self._sigma_ID_width = dict()
        for item in self._allowed_final_states:
            self._sigma_ID[item] = -1.0
            self._sigma_ID_width[item] = -1.0

        self.load_constraints()

        logger.info('Loaded experimental constraints. To change, use the set command')
        logger.info('Omega h^2 = %.4e +- %.4e' %(self._oh2_planck, self._oh2_planck_width))
        logger.info('Spin Independent cross section: %s' % self._dd_si_limit_file)
        logger.info('Spin Dependent cross section (p): %s' % self._dd_sd_proton_limit_file)
        logger.info('Spin Dependent cross section (n): %s' % self._dd_sd_neutron_limit_file)

#        for chan in self._allowed_final_states:
#            logger.info('Indirect Detection cross section for final state %s at velocity %.2e: %s'\
#                        % (chan, self._id_limit_vel[chan] ,self._id_limit_file[chan]))


    def load_constraints(self):
        #Load in direct detection constraints
        if self._dd_si_limit_file!='':
            self._dd_si_limit_mDM, self._dd_si_limit_sigma = np.loadtxt(self._dd_si_limit_file, unpack=True, comments='#')
        if self._dd_sd_proton_limit_file!='':
            self._dd_sd_p_limit_mDM, self._dd_sd_p_limit_sigma = np.loadtxt(self._dd_sd_proton_limit_file, unpack=True, comments='#')
        if self._dd_sd_neutron_limit_file!='':
            self._dd_sd_n_limit_mDM, self._dd_sd_n_limit_sigma = np.loadtxt(self._dd_sd_neutron_limit_file, unpack=True, comments='#')
                                                                                                            
        for channel, limit_file in self._id_limit_file.iteritems():
            if limit_file != '': # FF : need to redo the limit files as two columns                                                                             
             if 'MadDM_FermiLim' in limit_file:
                self._id_limit_mdm[channel]  = np.loadtxt(limit_file, unpack=True)[0]
                self._id_limit_sigv[channel] = np.loadtxt(limit_file, unpack=True)[3]
             else:  self._id_limit_mdm[channel] ,  self._id_limit_sigv[channel] = np.loadtxt(limit_file, unpack=True , comments = '#')

            else:
                self._id_limit_mdm[channel] = False
                self._id_limit_sigv[channel] = False

    #Returns a value in cm^2
    def SI_max(self, mdm):
        if (mdm < np.min(self._dd_si_limit_mDM) or mdm > np.max(self._dd_si_limit_mDM)):
            logger.warning('Dark matter mass value '+str(mdm)+' is outside the range of SI limit')
            return __infty__
        else:
            return np.interp(mdm, self._dd_si_limit_mDM, self._dd_si_limit_sigma)

    #Returns a value in cm^2
    def SD_max(self,mdm, nucleon):
        if nucleon not in ['n','p']:
            logger.error('nucleon can only be p or n')
            return __infty__
        else:
            if nucleon=='p':
                if (mdm < np.min(self._dd_sd_p_limit_mDM) or mdm > np.max(self._dd_sd_n_limit_mDM)):
                    logger.warning('Dark matter mass value '+str(mdm)+' is outside the range of SD limit')
                    return __infty__
                else:
                    return np.interp(mdm, self._dd_sd_p_limit_mDM, self._dd_sd_p_limit_sigma)
            elif nucleon=='n':
                if (mdm < np.min(self._dd_sd_n_limit_mDM) or mdm > np.max(self._dd_sd_n_limit_mDM)):
                    logger.warning('Dark matter mass value '+str(mdm)+' is outside the range of SD limit')
                    return __infty__
                else:
                    return np.interp(mdm, self._dd_sd_n_limit_mDM, self._dd_sd_n_limit_sigma)

    #Returns a value in cm^3/s
    def ID_max(self,mdm, channel):
        if (mdm < np.min(self._id_limit_mdm[channel]) or mdm > np.max(self._id_limit_mdm[channel])):
            logger.warning('Dark matter mass value %.2e for channel %s is outside the range of ID limit' % (mdm, channel))
            return __infty__
        else:
            return np.interp(mdm, self._id_limit_mdm[channel], self._id_limit_sigv[channel])

# This class holds some functionalities to load spectra from the PPPC4DMID files or the output from pythia8,
# interpolates the PPPC4DMID spectra for arbitrary values of DM etc.
class Spectra:

    def __init__(self):
        
        self.spectra_id    = {'px':'antiprotons', 'gx':'gammas', 'nuex':'neutrinos_e', 'numux':'neutrinos_mu', 'nutaux':'neutrinos_tau', 'ex':'positrons'}

        # this will be calculated when if needed
        self.spectra       = {'x':[] , 'antiprotons':[], 'gammas':[], 'neutrinos_e':[], 'neutrinos_mu':[], 'neutrinos_tau':[], 'positrons':[] }     
        self.spectra_earth = {'x':[] , 'positrons':[] , 'antiprotons':[], 'neutrinos_tau':[],'neutrinos_mu':[], 'neutrinos_e':[] }

        self.flux_source   = {'e':[] , 'positrons': {'dNdE':[],'dPhidE':[] }, 'antiprotons': {'dNdE':[] , 'dPhidE':[] } , 
                              'neutrinos_tau':{'dNdE':[] , 'dPhidE':[] },'neutrinos_mu':{'dNdE':[],'dPhidE':[] }, 
                              'neutrinos_e':{'dNdE':[] , 'dPhidE':[] } , 'gammas':{'dNdE':[] , 'dPhidE':[] } }


        # fluxes_earth only filled if using PPPC4DMID method (no reading from Dragon output)
        self.flux_earth_positrons    = {'e':[] , 'positrons': { 'dPhidlogE':[] } } 
        self.flux_earth_neutrinos    = {'e':[] , 'neutrinos_tau':{'dNdE':[]}, 'neutrinos_mu':{'dNdE':[]}, 'neutrinos_e':{'dNdE':[] } }

        self.channels      = ['ee', 'mumu', 'tautau', 'qq', 'cc', 'bb', 'tt', 'ZZ', 'WW', 'hh', 'gammagamma', 'gg']

        self.map_allowed_final_state_PPPC = {'qqx':'qq', 'ccx':'cc', 'gg':'gg', 'bbx':'bb', 'ttx':'tt',
                                             'e+e-':'ee', 'mu+mu-':'mumu', 'ta+ta-':'tautau', 'w+w-':'WW', 'zz':'ZZ', 'hh':'hh' }

    def check_mass(self,mdm):
        if (mdm < 5.0 or mdm > 100000):
            logger.error('DM mass outside the range of available PPPC spectra. Please run spectra generation with pythia8') # Break and ask the user to download the Tables
            return 
        else: return True 

    def load_PPPC_source(self,corr = '',load = True):
        if load: 
          if (not os.path.isfile(PPPCDIR+'/PPPC_Tables_EW.npy') and not os.path.isfile(PPPCDIR+'/PPPC_Tables_noEW.npy')):
            logger.error('PPPC4DMID Spectra at source not found! Please install by typing install PPPC4DMID') # Break and ask the user to download the Tables       
            return

          if not corr:
             dic =  np.load(PPPCDIR+'/PPPC_Tables_noEW.npy').item()
             misc.sprint('PPPC4DMID Spectra at source loaded')
             return dic

          elif corr == 'ew':
             dic =  np.load(PPPCDIR+'/PPPC_Tables_noEW.npy').item()
             misc.sprint('PPPC4DMID Spectra at source (with EW corrections) loaded')
             return dic

    # not yet implemented (tables missing) 
    def load_PPPC_earth(self, prof = 'Ein'):
         if (not os.path.isfile(PPPCDIR+'/PPPC_Tables_epEarth_'+prof+'.npy') ):
            logger.error('PPPC4DMID Spectra at Earth not found! Please install by typing install PPPC4DMID') # Break and ask the user to download the Tables                     
            return
        
         else:
            #print 'FF dic_path', '/PPPC_Tables_epEarth_'+prof+'.npy' 

            dic = np.load(PPPCDIR+'/PPPC_Tables_epEarth_'+prof+'.npy').item()
            #print "FF dic['x']" , dic['x']
            return dic 


    # FF this function extracts the values of the spectra interpolated linearly between two values mdm_1 and mdm_2                                              
    # mdm is the DM candidate mass, spectrum is gammas, positron etc, channel is the SM annihilation e.g. bbar, hh etc.                        
    # FF remember to CHECK if it works when min and max are the same values, i.e. exactly for a value in the Masses lists!!!                              
    def interpolate_spectra(self, sp_dic, mdm = '', spectrum = '' , channel = '' , earth = False , prof = 'Ein' , prop = 'MED' , halo_func = 'MF1'):
        M   = sp_dic['Masses']
        dm_min =  max([m for m in M if m <= mdm])  # extracting lower mass limit to interpolate from                                                                  
        dm_max =  min([m for m in M if m >= mdm])  # extracting upper mass limit to interpolate from                                                                    

        interpolated = []

        # FF if the MD requested is exactly one of the values in the table, returns the values with no interpolation                                                   
        if not earth: # two spectra: source and earth. First method is for source, secodn is for earth (when earth == True)
            key = channel
        elif earth:
            key = channel+'_'+prof+'_'+prop+'_'+halo_func

        if dm_min == dm_max :
               #print 'FF spec_1 ' , sp_dic[spectrum][ str(dm_min) ][channel] , ' channel ', channel
               return sp_dic[spectrum][ str(dm_min) ][key]
        spec_1 = sp_dic[spectrum][ str(dm_min) ][key]
        spec_2 = sp_dic[spectrum][ str(dm_max) ][key]
     
        #print 'FF spec_1 spec_2 ', spec_1 , spec_2 , key 
        for x in range(len(spec_1)): # the spectrum values are ordered as 'x' vaues extracted from the Log[10,x] in the PPPC Tables                       
               interp_function = interp1d([dm_min, dm_max], [spec_1[x],spec_2[x]] )
               value =  interp_function(mdm)
               interpolated.append(value)

        return interpolated

      

class Fermi_bounds:

    def __init__(self):

        self.nBin = 24
        self.j0 = 3.086e21 # convention spectra    
        self.dSph_jfac_file =  pjoin(MDMDIR, 'Fermi_Data', 'Jfactors.dat')
        self.dSph_ll_files_path = pjoin(MDMDIR, 'Fermi_Data', 'likelihoods')
        self.dwarves_list = ['coma_berenices', 'draco', 'segue_1', 'ursa_major_II', 'ursa_minor', 'reticulum_II' ] # 6 highest J
        self.dwarveslist_all = self.extract_dwarveslist() 
        self.dw_in = self.dw_dic()
        self.ll_tot = ''

# This function reads the list of dwarves from the Jfactor.dat file
    def extract_dwarveslist(self):
        dwarveslist_all = []
        op = open(self.dSph_jfac_file,'r')
        for line in [ line for line in op.readlines() if 'Row' in line]:
            dwarveslist_all.append (line.split(' ')[3].replace('\n',''))
        return dwarveslist_all

    def dw_dic(self):
        
        #dSph_jfac_file, dSph_ll_files_path, nBin = self.dSph_jfac_file , self.dSph_ll_files_path , self.nBin
        #dwarflist_all = self.extract_dwarveslist(dSph_jfac_file) # create the list of dwarves reading the Jfactor dat file                                                          
        dwarves = self.dwarves_list
        dwarves_all = self.dwarveslist_all  

        dSph_ll_files = [ 'like_' + dwarf + '.txt' for dwarf in dwarves_all ]

        jbar_dSph       = np.loadtxt(self.dSph_jfac_file, unpack=True)[3]
        jbar_dSph_error = np.loadtxt(self.dSph_jfac_file, unpack=True)[4]

        dw_dic_coll = {}

        for dwarf in dwarves_all:
            for D in dwarves:
               if D == dwarf:
                  elem = 'like_'+str(dwarf)+'.txt'

                  if elem in dSph_ll_files:
                     dwindex = dSph_ll_files.index(elem)
                     likefile = dSph_ll_files[dwindex]
                     data = np.loadtxt(self.dSph_ll_files_path+'/' + likefile, unpack=True)
                     efluxes = (1e-3*data[2]).reshape(self.nBin,-1) # convert flux list to 2D list of fluxes for all ebins, convert to GeV                                               
                     logLikes = data[3].reshape(self.nBin,-1)

                     likes = [ interp1d(f,l,bounds_error=False, fill_value=-1e5) for f,l in zip(efluxes,logLikes) ]

                     ll_null = 0.0
                     for like in likes:
                         ll_null+=like(0)

                     jfact = jbar_dSph[dwindex]
                     jfacterr = jbar_dSph_error[dwindex]

                     dict = { 'Jfac': jfact,
                              'Jfac_err': jfacterr,
                              'll_null': ll_null,
                              'likelihood': likes }

                     dw_dic_coll[dwarf] = dict

        #print "dwarfs included:", dw_dic_coll.keys()
        return dw_dic_coll

    def eflux(self,spectrum, emin=1e2, emax=1e5, quiet=False):
       """ Integrate a generic spectrum, multiplied by E, to get the energy flux.                                                                                             
       """
       espectrum = lambda e: spectrum(e)*e
       tol = min(espectrum(emin),espectrum(emax))*1e-10
       try:
          return quad(espectrum,emin,emax,epsabs=tol,full_output=True)[0]
       except Exception, msg:
          logger.info('Numerical error "%s" when calculating integral flux.' % msg)
       return np.nan

    def marg_like_dw(self,dw_in_i,pred,marginalize):
      
        j0, nBin = self.j0 , self.nBin

        j,jerr,like_inter = dw_in_i['Jfac'], dw_in_i['Jfac_err'], dw_in_i['likelihood']
       
        # function to be minimized:                                                                                                                                           
        def chi2min(x):

            flux_like = 0.0

            for i in range(nBin):

                like_i = like_inter[i]
                pred_i = pred[i]
                flux = pred_i*10**(j+x*jerr)/j0
                flux_like += like_i(flux)
            jfac_like = - 0.5*x**2
            return -2.0*(flux_like+jfac_like)

        if marginalize:
           res = minimize_scalar(chi2min,method='bounded',bounds=(-5.0,5.0))

           jsigma = res.x
           ll_max = -0.5*res.fun
        else:
           ll_max = -0.5*chi2min(0.0)
           jsigma = 0.0

        return ll_max, jsigma


    def res_tot_dw(self,pred,marginalize):

        dw_in = self.dw_in
        ll_tot = 0.0
        ll_null = 0.0
        j_factors = {}

        for k,v in dw_in.iteritems():
            marg_like = self.marg_like_dw(v,pred,marginalize)
            ll_tot += marg_like[0]
            ll_null += v['ll_null']
            j_factors[k] = marg_like[1]

        pval = self.compute_pvalue(ll_tot,ll_null)
        self.ll_tot = ll_tot
        return ll_tot, ll_null, pval, j_factors
       

    def compute_pvalue(self,ll_tot,ll_null):
    
        pvalue1 = lambda x: 1-gammainc(1/2.0,x/2.0)
        ts = -2*ll_tot+2*ll_null
        oneortwosidedfactor = 1
    
        if ts>=0.0:
            pval = pvalue1(ts)/oneortwosidedfactor
        else:
            if oneortwosidedfactor < 2:
                pval=1
            else:
                pval = 1-pvalue1(-ts)/oneortwosidedfactor

        return 1-pval

    def Fermi_sigmav_lim(self, mDM, x = '' , dndlogx = '' , marginalize = True, sigmavmin=1e-35, sigmavmax=1e-15, step_size_scaling=1.0, cl_val = 0.95, maj_dirac=''):
        np.seterr(divide='ignore', invalid='ignore')   # Keep numpy from complaining about dN/dE = 0...                                                                       
        j0 , nBin = self.j0 , self.nBin # convention spectra                                                                                                                     
        dw_in = self.dw_in
        sigmav0 = 1e-26

        emins = [0.5, 0.666760716, 0.8891397050000001, 1.1856868500000002, 1.5811388300000002, 2.10848252, 2.81170663, 3.7494710500000004, 5.0, 6.667607159999999, 
                 8.89139705 , 11.856868500000001, 15.8113883, 21.0848252, 28.117066299999998, 37.494710500000004, 50.0, 66.6760716, 88.91397049999999, 
                 118.568685, 158.11388300000002, 210.848252, 281.170663, 374.94710499999997]
        emaxs = [0.666760716, 0.8891397050000001, 1.1856868500000002, 1.5811388300000002, 2.10848252, 2.81170663, 3.7494710500000004, 5.0, 6.667607159999999, 
                 8.89139705, 11.856868500000001, 15.8113883, 21.0848252, 28.117066299999998, 37.494710500000004, 50.0, 66.6760716, 88.91397049999999, 
                 118.568685, 158.11388300000002, 210.848252, 281.170663, 374.94710499999997, 500.0]

    # Load the spectrum E, dN/dE                                                                                                                                              
    #energy,dnde = np.loadtxt(specfile,unpack=True)                                                                                                                           
        logx = np.log10(x)
        
        energy = mDM*10**logx

        dnde = (dndlogx/(mDM*10**logx*2.30259))
        log_energy = np.log10(energy)
        log_dnde = np.log10(dnde)
        log_interp = interp1d(log_energy,log_dnde)
        spectrum = lambda e: np.nan_to_num(10**( log_interp(np.log10(e)) )) if (e <= energy.max() and e >= energy.min()) else 0.0

        pred = np.array([self.eflux(spectrum,e1,e2)/mDM**2*sigmav0*j0/(2.*maj_dirac*math.pi) for e1,e2 in zip(emins,emaxs)])

        def find_sigmav(x,pred,dw_in,marginalize):

            pred_sigma = pred*10**(-x)/sigmav0
            pvalue = self.res_tot_dw(pred_sigma,marginalize)[2]

            atanhpval = np.arctanh(pvalue - 1e-9)
            atanclval = np.arctanh(cl_val - 1e-9)
            return (atanhpval-atanclval)**2

        find_sig = lambda x: find_sigmav(x,pred,dw_in,marginalize)
       
        # brute methode:
        brute_range_min = -np.log10(sigmavmax)
        brute_range_max = -np.log10(sigmavmin)
        if marginalize:
            num_steps = int(step_size_scaling*2.0*(brute_range_max-brute_range_min))
        else:
            num_steps = int(step_size_scaling*5.0*(brute_range_max-brute_range_min))

        res = brute(find_sig,[(brute_range_min,brute_range_max)], Ns=num_steps, full_output=False, finish=fmin)
        sigmav_ul = 10**(-res[0])
    
        p_value = -1

        pred_sigma = pred*sigmav_ul/sigmav0
        result  = self.res_tot_dw(pred_sigma,marginalize)
        p_value = result[2]     

        if p_value <= cl_val*0.98 or p_value >= cl_val*1.02:
           sigmav_ul= -1
           print " WARNING: increase range (sigmavmin,sigmavmax) and/or step_size_scaling!"
        
        return [sigmav_ul , p_value , self.ll_tot]    




class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    GRAY = '\033[90m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#Conversion from GeV^-2 to pb etc.
GeV2pb = 3.894E8
pb2cm2  = 1.00E-36
cm22pb  = 1.00E36
pb2cm3 = 2.99E-26



#===============================================================================
# CommonRunCmd
#===============================================================================
class MADDMRunCmd(cmd.CmdShell):
    
    intro_banner=\
  "            ====================================================\n"+\
  "            |                  "+bcolors.OKBLUE+"  MadDM v3.0                     "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "                                                                               \n"+\
  "                #########                                                        \n"+\
  "             ###\\\\####//#####              Launchpad:  launchpad.net/maddm      \n"+\
  "           ######\\\\##//########                                              \n"+\
  "          ########\\\\//###########                                            \n"+\
  "         #########//\\\\############                                           \n"+\
  "        #########//##\\\\############      "+bcolors.FAIL+"             arXiv:1308.4955        \n"+bcolors.ENDC+\
  "       ########//#####\\\\###########                   "+bcolors.FAIL+"arXiv:1505.04190        \n"+bcolors.ENDC+\
  "       ######################### ## ___________________________________________\n"+\
  "       ####################### 0  # "+bcolors.OKGREEN+" _     _               _  _____   _     _  \n"+bcolors.ENDC+\
  "       #############   0  ###    ## "+bcolors.OKGREEN+"| \   / |   ___    ___|| | ___ \ | \   / | \n"+bcolors.ENDC+\
  "       ##############    #########  "+bcolors.OKGREEN+"||\\\\ //|| / __ |  / __ | ||   || ||\\\\ //|| \n"+bcolors.ENDC+\
  "        ##########################  "+bcolors.OKGREEN+"||  V  || ||__||  ||__|| ||___|| ||  V  || \n"+bcolors.ENDC+\
  "         ###################   ##   "+bcolors.OKGREEN+"||     || \_____\ \____| |_____/ ||     || \n"+bcolors.ENDC+\
  "          ############       ###    ___________________________________________\n"+\
  "           ##########    ######                                                 \n"+\
  "             ################                                                   \n"+\
  "                 ########                                                       \n"+\
  "                                                                                    \n"+\
  "            ====================================================\n"+\
  "            |           "+bcolors.OKBLUE+" A MadGraph5_aMC@NLO plugin.            "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "" 
    
    
    def __init__(self, dir_path, options, *args, **opts):
        """common"""
  
        self.in_scan_mode = False # please modify this value only with "with" statement
        cmd.Cmd.__init__(self, *args, **opts)
        # Define current MadEvent directory
        if dir_path is None and not MG5MODE:
            dir_path = root_path

        self.dir_path = dir_path
        self.param_card_iterator = [] #an placeholder containing a generator of paramcard for scanning
        
        # Define self.proc_characteristics (set of information related to the
        # directory status
        self.get_characteristics()

        self._two2twoLO = False
        self._dNdE_setup = False #If true, this flag allows the code to skip loading and caldulating dNdE
        #self._fit_parameters= []

        self._run_couplings = False

        #A shared param card object to the interface.
        #All parts of the code should read this card.
        #Set at the beginning of launch()
        self.param_card = None

        self.multinest_running = False
        self.auto_width = set() # keep track of width set on Auto

        self.limits = ExpConstraints()
        self.options = {}

        self.Spectra = Spectra()
        self.Fermi   = Fermi_bounds()
        self.MadDM_version = '3.0'

    def preloop(self,*args,**opts):
        super(Indirect_Cmd,self).preloop(*args,**opts)
        self.prompt = 'Maddm:%s' % self.prompt 
    
    def get_characteristics(self, path=None):
        """reads the proc_characteristics file and initialises the correspondant
        dictionary"""
        
        if not path:
            path = os.path.join(self.dir_path, 'matrix_elements', 'proc_characteristics')
        
        self.proc_characteristics = MGoutput.MADDMProcCharacteristic(path)
        return self.proc_characteristics

    def check_param_card(self, path, run=True):
        """
        1) Check that no scan parameters are present
        2) Check that all the width are defined in the param_card.
        - If a scan parameter is defined. create the iterator and recall this function 
          on the first element.
        - If some widths are set on 'Auto', call the computation tools.
        - keep track of the width of auto if some present"""
        
        pattern_scan = re.compile(r'''^(decay)?[\s\d]*scan''', re.I+re.M)  
        pattern_width = re.compile(r'''decay\s+(\+?\-?\d+)\s+auto(@NLO|)''',re.I)
        text = open(path).read()
               
        if pattern_scan.search(text):
            if not isinstance(self, cmd.CmdShell):
                # we are in web mode => forbid scan due to security risk
                raise Exception, "Scans are not allowed in web mode"
            # at least one scan parameter found. create an iterator to go trough the cards
            main_card = param_card_mod.ParamCardIterator(text)
            self.param_card_iterator = main_card
            first_card = main_card.next(autostart=True)
            first_card.write(path)
            return self.check_param_card(path, run)
        
        pdg_info = pattern_width.findall(text)
        if pdg_info:
            pdg = [pdg for pdg,nlo in pdg_info]
            self.auto_width = set(pdg)
            if run:
                if not self.in_scan_mode:
                    logger.info('Computing the width set on auto in the param_card.dat')
                has_nlo = any(nlo.lower()=="@nlo" for _,nlo in pdg_info)
                if not has_nlo:
                    self.do_compute_widths('%s --path=%s' % (' '.join(pdg), path))
                    #self.run_mg5(['compute_widths %s --path=%s' % (' '.join(pdg), path)])
                else:
                    self.run_mg5(['compute_widths %s --path=%s --nlo' % (' '.join(pdg), path)])
            else:
                logger.info('''Some widths are on Auto in the card.
    Those will be computed as soon as you have finish editing the cards.
    If you want to force the computation right now and re-edit
    the cards afterwards, you can type \"compute_wdiths\".''')


    def do_compute_widths(self, line):
                
        if self.maddm_card['only_two_body_decays']:
            line = ' --body_decay=2 ' + line
        #return self.run_mg5(['compute_widths --body_decay=2 ' + line])
        if not hasattr(self, 'mg5'):
            return self.run_mg5(['compute_widths ' + line])
        elif not hasattr(self.mg5, '_curr_model'):
            return self.run_mg5(['compute_widths ' + line])
        else:
            self.mg5.do_compute_widths(line, model=self.mg5._curr_model)#, decaymodel=self.mg5._curr_decaymodel)
    
    def help_compute_widths(self, line):
        
        return self.run_mg5([' help compute_widths ' + line])    


    ############################################################################
    def run_mg5(self, commands, model=None):
        """ """
        
        if not hasattr(self, 'mg5'):
            if self.mother:
                self.mg5 = self.mother
            else:
                import madgraph.interface.master_interface as master_interface
                self.mg5 = master_interface.MasterCmd()

        if not model:
            model = self.proc_characteristics['model']
        if self.mg5._curr_model.get('modelpath+restriction') != model:
            self.mg5.do_import(model)
            
            
        for line in commands:
            self.mg5.exec_cmd(line, errorhandling=False, printcmd=False, 
                              precmd=False, postcmd=False, child=False)
            
    ############################################################################
    def do_launch(self, line):
        """run the code"""

        param_path  = pjoin(self.dir_path,'Cards', 'param_card.dat')
        self.param_card = param_card_mod.ParamCard(param_path)

        args = line.split()
        if '-f' in args or '--force' in args:
            force = True
        else:
            force = False
        self.ask_run_configuration(mode=[], force=force)


        if not self.multinest_running:
            self.compile()

        nb_output = 1
 
        # FF Does this do anything ?
        output = pjoin(self.dir_path, 'output', 'maddm_%s.out') 
        while os.path.exists(output % nb_output):
            nb_output +=1
        else:
            output = pjoin('output', 'maddm.out')
        misc.sprint(output)


        misc.call(['./maddm.x', output], cwd =self.dir_path)
        #process = subprocess.Popen(['./maddm.x'], cwd =self.dir_path, stdout=subprocess.PIPE)
        #Here we read out the results which the FORTRAN module dumped into a file
        #called 'maddm.out'. The format is such that the first line is always relic density
        # , second line is the nucleon scattering cross section (SI) for proton, third is SI
        # nucleon cross section for the neutron etc. If a quantity was not calculated, we output -1
        result = {}
        result['GeV2pb*pb2cm2']   = GeV2pb*pb2cm2 # conversion factor                                                                                                           

        sigv_indirect = 0.
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])


        result['tot_SM_xsec'] = -1
        sigv_indirect = 0.

        for line in open(pjoin(self.dir_path, output)):
      
                splitline = line.split()
                #If capture rate is calculated.
                if 'ccap' in line:
                    oname = splitline[0].strip(':')+'_'+splitline[1]
                    val = splitline[2]
                    result[oname.split(':')[0] ] = val

                else:
                    if self._two2twoLO:
                        if 'sigma*v' in line: 
                            sigv_temp = float(splitline[1])
                            oname = splitline[0].split(':',1)[1] .split('_')
                            #oname2 = oname.split('_')
                            oname = oname[0]+'_'+oname[1] #To eliminate the annoying suffix coming from the '/' notation
                            result['taacsID#'     + oname] = sigv_temp 
                            result['err_taacsID#' + oname] = 0 
                            if oname.split('_')[1] in ['uux','ddx','ssx']:
                                result['lim_taacsID#'+oname] = self.limits.ID_max(mdm, 'qqx')
                            elif oname.split('_')[1] in self.limits._allowed_final_states:
                                result['lim_taacsID#'+oname] = self.limits.ID_max(mdm, oname.split('_')[1]) 
                            elif oname.split('_')[1] not in self.limits._allowed_final_states:
                                result['lim_taacsID#'+oname] = -1
                            sigv_indirect += sigv_temp
                    result[splitline[0].split(':')[0]] = float(splitline[1])

        result['taacsID'] = sigv_indirect
                            
        np_names = ['g','nue','numu','nutau']

        if str(self.mode['indirect']).startswith('flux'):
            for chan in np_names:
                result['flux_%s' % chan] = -1.0

        # Calculating the xsi factor
        if result['Omegah^2'] < 0:
            result['xsi'] = 1.0   
        elif result['Omegah^2'] < self.limits._oh2_planck and result['Omegah^2'] > 0:
            result['xsi'] = result['Omegah^2'] / self.limits._oh2_planck
        else: result['xsi'] = 1.0

        if self.mode['direct']:
            result['lim_sigmaN_SI_n'] = self.limits.SI_max(mdm)
            result['lim_sigmaN_SI_p'] = self.limits.SI_max(mdm)
            result['lim_sigmaN_SD_p'] = self.limits.SD_max(mdm, 'p')
            result['lim_sigmaN_SD_n'] = self.limits.SD_max(mdm, 'n')
            result['GeV2pb*pb2cm2']   = GeV2pb*pb2cm2 # conversion factor 

        self.last_results = result
        self.last_results['run'] = nb_output
        
        
        #print 'FF last result', result 
                  
#        if self.mode['indirect'] and not self._two2twoLO:
#            with misc.MuteLogger(names=['madevent','madgraph'],levels=[50,50]):
#                self.launch_indirect(force)

        if self.mode['indirect']:
            ### FF TO_DO ADD HERE the loading of the PPPC tables !!!
            #with misc.MuteLogger(names=['madevent','madgraph'],levels=[50,50]):
                self.launch_indirect(force)

        #if sigv_indirect:
        #    result['taacsID'] = sigv_indirect
        #    result['err_taacsID'] = math.sqrt(sigv_indirect_error)


        
        if not self.in_scan_mode and not self.multinest_running:
            self.print_results()

        # Saving or removing output    
        # self.save_remove_indirect_output(scan = False, point_number= nb_output )


        

        # --------------------------------------------------------------------#
        #   THIS PART IS FOR MULTINEST SCANS
        # --------------------------------------------------------------------#

        #multinest can be launched only after one launch has been executed
        if self.mode['nestscan'] and not self.multinest_running:
            self.multinest_running = True
            self.launch_multinest()
            self.multinest_running = False



        # --------------------------------------------------------------------#
        #   THIS PART IS FOR SEQUENTIAL SCANS
        # --------------------------------------------------------------------#

        #    logger.info("relic density  : %.2e ", self.last_results['omegah2'])
        if self.param_card_iterator:

            param_card_iterator = self.param_card_iterator

            parameters, values =  param_card_iterator.param_order , param_card_iterator.itertag
            self.param_card_iterator = []

            self.save_remove_indirect_output(scan = True, point_number= nb_output ) ## this is to remove or save spectra, not the scan summary file!

            # *** Initialize a list containing the desired variables in the summary output of the scan
            order = ['run']

            print 'FF prima' , parameters, values 
            # *** saving the name of the iterated parameters, and the values in the results dictionary
            for par,val in zip(parameters, values):
                order.append(par)
                self.last_results[par] = val
                #print 'FF ' , self.last_results


            # *** Relic density
            if self.mode['relic']:
                # OLD order += ['omegah2', 'x_freezeout', 'sigmav_xf']
                order += ['Omegah^2','x_f', 'sigmav(xf)']
            order.append('xsi')

            # *** Direct Detection
            if self.mode['direct'] :
                # OLD order += ['sigmaN_SI_proton', 'sigmaN_SI_neutron', 'sigmaN_SD_proton',
                #        'sigmaN_SD_neutron']
                order += ['sigmaN_SI_p', 'lim_sigmaN_SI_p', 
                          'sigmaN_SI_n', 'lim_sigmaN_SI_n',
                          'sigmaN_SD_p', 'lim_sigmaN_SD_p',
                          'sigmaN_SD_n', 'lim_sigmaN_SD_n']

           
            if self.mode['direct'] == 'directional':
                order += ['Nevents', 'smearing']

            if self.mode['capture']:
                detailled_keys = [k for k in self.last_results if k.startswith('ccap_') and '#' not in k]
                for key in detailled_keys:
                    order += [key]
 
            # *** Indirect detection
            if self.mode['indirect']:
    
                #if not self._two2twoLO:
                #order +=['halo_velocity']#,'indirect', 'indirect_error']
                detailled_keys = [k for k in self.last_results if k.startswith('taacsID#') ]
                if len(detailled_keys)>1:
                    for key in detailled_keys:
                        #Useful is there is "/" syntax, because "no" suffix is added.                       
                        clean_key_list = key.split("_")
                        clean_key = clean_key_list[0]+"_"+clean_key_list[1]

                        order +=[clean_key]
                        order +=['lim_'+clean_key]

                order.append('taacsID')
                order.append('tot_SM_xsec')
                order.append('Fermi_sigmav')
                order.append('pvalue')
                order.append('like_tot')

                if self.mode['indirect'].startswith('flux'):
                    for channel in np_names:
                        order.append('flux_%s' % channel)

                #if self.maddm_card['sigmav_method'] == 'inclusive':
                #        order.append('tot_SM_xsec')
            #logger.info(order)            
            # FF to be fixed!
            #<=-------------- Mihailo commented out max_col = 10

#            print 'FF the order is', order
#            print 'FF last results', self.last_results 
            
            # removing extra entries if not in last_results
            #print 'FF the order is ', order 
            #order.append('taacsID')
            for elem in order:
                if elem not in self.last_results.keys():
                   order.remove(elem)


            #to_print = param_card_iterator.write_summary(None, order, nbcol=10)#, max_col=10)
            '''
            for line in to_print.split('\n'):
                if line:
                    logger.info(line)
                    ## added by chiara to check the width (next three lines)
                    # self._param_card = param_card_mod.ParamCard('/Users/arina/Documents/physics/software/maddm_dev2/test_width/Cards/param_card.dat')
                    # width = self.param_card.get_value('width', 5000000)
                    # logger.warning('--> WY0: %.2e' % width)
            '''
            #print 'FF testing writer'
            summary_file = pjoin(self.dir_path, 'output','PROVA_SUMMARY.txt')
            self.write_scan_output(out_path = summary_file , keys = order, header = True )
            self.write_scan_output(out_path = summary_file , keys = order )
            #print 'FF nb_output' , nb_output
            #check if the param_card defines a scan.
            with misc.TMP_variable(self, 'in_scan_mode', True):
                with misc.MuteLogger(names=['cmdprint','madevent','madgraph','madgraph.plugin'],levels=[50,50,50,20]):
                    
                    for i,card in enumerate(param_card_iterator):
                        print 'FF parameters, values ' ,parameters, values 
                        point_number = nb_output+i+1
                        self.last_results['run'] = point_number
                        
                        card.write(pjoin(self.dir_path,'Cards','param_card.dat'))
                        self.exec_cmd("launch -f", precmd=True, postcmd=True, errorhandling=False)

                        for par,val in zip(param_card_iterator.param_order, param_card_iterator.itertag):
                            self.last_results[par] = val
                        #print 'FF last results in lop' , self.last_results
                        #param_card_iterator.store_entry(nb_output+i+1, self.last_results)
                        #print 'FF ' , self.last_results
                        self.write_scan_output(out_path = summary_file , keys = order, header = False)

                        ### the following three lines are added by chiara to check the widht = auto function 
                        # self._param_card = param_card_mod.ParamCard('/Users/arina/Documents/physics/software/maddm_dev2/test_width/Cards/param_card.dat')
                        # width = self.param_card.get_value('width', 5000000)
                        # logger.warning('--> try again WY0: %.2e' % width)
                        #<=-------------- Mihailo commented out max_col = 10
                        #logger.info('Results for the point \n' + param_card_iterator.write_summary(None, order, lastline=True,nbcol=10)[:-1])#, max_col=10)[:-1])
                        self.save_remove_indirect_output(scan = True, point_number= nb_output+1+i )



            param_card_iterator.write(pjoin(self.dir_path,'Cards','param_card.dat'))
            #name = misc.get_scan_name('maddm_%s' % (nb_output), 'maddm_%s' % (nb_output+i))
            #name = misc.get_scan_name('maddm_%s' % (nb_output), 'maddm_%s' % (nb_output+1))                                                                                         
            #path = pjoin(self.dir_path, 'output','scan_%s.txt' % name)
            #logger.info("write all results in %s" % path ,'$MG:color:BLACK')

            #print 'FF the order is', order
            #print 'FF last results', self.last_results
            #print 'FF param card iterator ' , param_card_iterator 
            #param_card_iterator.write_summary(path, order)




    def launch_multinest(self):

        if self.in_scan_mode:
            logger.error('You can not use scan syntax in the param_card.dat and run multinest!')
            return

        if not os.path.exists(pjoin(self.dir_path, 'multinest_chains')):
            os.mkdir(pjoin(self.dir_path, 'multinest_chains'))


        mnest = Multinest(self)
        mnest.load_parameters(pjoin(self.dir_path, 'Cards', 'multinest_card.dat'))
        #the following two must be set for the code to work

        resume_chain = False
        prefix_found = False
        filelist = os.listdir(pjoin(self.dir_path,'multinest_chains'))
        for f in filelist:
            if f.startswith(mnest.options['prefix']):
                prefix_found = True

        if prefix_found:
            if self.ask('A multinest chain with prefix %s already exists. Overwrite [y] or resume [n]? [n]'\
                                % mnest.options['prefix'], 'n', ['y','n'],timeout=60):
                for file in filelist:
                    if file.startswith(mnest.options['prefix']):
                        os.remove(pjoin(self.dir_path,'multinest_chains',file))
            else:
                resume_chain = True

        mnest.launch(resume = resume_chain)
        mnest.write_log()


    def launch_indirect(self, force):
        """running the indirect detection"""

        print "FF self.mode['indirect'] ", self.mode['indirect']

        if not os.path.exists(pjoin(self.dir_path, 'Indirect')):
            self._two2twoLO = True
            #return
        elif self.maddm_card['sigmav_method'] == 'inclusive':
            self._two2twoLO = True
            #return 
        
        if not self.in_scan_mode: 
            logger.info('Running indirect detection')
        if not hasattr(self, 'me_cmd'):
            try:
                os.remove(pjoin(self.dir_path,'indirect', 'RunWeb'))
            except Exception:
                pass
            self.me_cmd = Indirect_Cmd(pjoin(self.dir_path, 'Indirect'))
        elif self.me_cmd.me_dir != pjoin(self.dir_path, 'Indirect'):
            self.me_cmd.do_quit()
            self.me_cmd = Indirect_Cmd(pjoin(self.dir_path, 'Indirect'))

        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        # FF madevent or Reshuffling methods fro sigmav
        print 'FF sigmav method *** ', self.maddm_card['sigmav_method']

        if self.maddm_card['sigmav_method'] != 'inclusive':
           runcardpath = pjoin(self.dir_path,'Indirect', 'Cards', 'run_card.dat')
           run_card = banner_mod.RunCard(runcardpath)

           vave_temp = self.maddm_card['vave_indirect']

           # ensure that VPM is the central one for the printout (so far)
           self.last_results['taacsID'] = 0.0
           #for i,v in enumerate(scan_v):

           run_card['ebeam1'] = mdm * math.sqrt(1+vave_temp**2)
           run_card['ebeam2'] = mdm * math.sqrt(1+vave_temp**2)
           run_card['use_syst'] = False
           run_card.remove_all_cut()
        
           misc.sprint("check how to set nevents!")
           if os.path.exists(pjoin(self.dir_path, 'Cards', 'pythia8_card.dat')):
               py8card = Indirect_PY8Card(pjoin(self.dir_path, 'Cards', 'pythia8_card.dat'))
               if py8card['Main:NumberOfEvents'] != -1:
                   run_card['nevents'] = py8card['Main:NumberOfEvents']
        
               run_card.write(runcardpath)
        
           if self.maddm_card['sigmav_method'] == 'madevent':
               self.me_cmd.do_launch('-f')
           elif self.maddm_card['sigmav_method'] == 'reshuffling': 
               cmd = ['launch',
                   'reweight=indirect',
                   'edit reweight --before_line="launch" change velocity %s' % vave_temp]
               misc.sprint("using reshuffling")
               self.me_cmd.import_command_file(cmd)
                        
           for key, value in self.me_cmd.Presults.iteritems():
               clean_key_list = key.split("/")
               clean_key =clean_key_list[len(clean_key_list)-1].split('_')[1] +'_'+  clean_key_list[len(clean_key_list)-1].split('_')[2] 
 
               if key.startswith('xsec'):
                #<------- FIX THIS. GET RID OF VAVE_TEMP. THIS WILL JUST GET INTEGRATED BY MADEVENT
                  self.last_results['taacsID#%s' %(clean_key)] = value* pb2cm3
                  self.last_results['taacsID'] += value* pb2cm3
               elif key.startswith('xerr'):
                  self.last_results['err_taacsID#%s' %(clean_key)] = value * pb2cm3


           self.run_pythia8_for_flux() 
              
           if self.maddm_card['indirect_flux_source_method'] == 'pythia8':
                logger.info('Calculating Fermi limit using using pythia8 gamma rays spectrum')
           elif 'pythia' not in self.maddm_card['indirect_flux_source_method']:
                logger.warning('Since pyhtia8 is run, using pythia8 gamma rays spectrum (not PPPC4DMID Tables)')

           self.read_py8spectra()


#        else self.maddm_card['sigmav_method'] == 'inclusive':
        else:
            logger.info('Calculating Fermi limit using PPPC4DMID spectra')
            self.read_PPPCspectra()



        # ****** Calculating Fermi Limits

        self.last_results['Fermi_sigmav'] = -1
        x, gammas = self.Spectra.spectra['x'] , self.Spectra.spectra['gammas']
        #print ' FF spectra', x , '' , gammas                                                                                                                               
        if not gammas:
            logger.error('The gamma spectrum is empty!Will not calculate Fermi limit')
            self.last_results['Fermi_sigmav'] = -1
            self.last_results['pvalue']       = -1
            self.last_results['like_tot']     = -1
        elif gammas: 
            #print 'FF majorana dirac' , self.norm_Majorana_Dirac() 
            sigmav = self.Fermi.Fermi_sigmav_lim(mdm, x , gammas ,maj_dirac= self.norm_Majorana_Dirac() ) # Returns a 3d array with [sigmav,pvalue,likel.]
            self.last_results['Fermi_sigmav'] = sigmav[0] # FF store the results from the Fermi limit calculation                                    
            self.last_results['pvalue']       = sigmav[1]
            self.last_results['like_tot']     = sigmav[2]
            #print 'FF sigmav , value , likeli ' , sigmav[0] , sigmav[1] , sigmav[2]
        #print "FF self.mode['indirect'] ***************  " + self.mode['indirect']

        # ****** Calculating Fluxes Source

        if self.mode['indirect'].startswith('flux') :
            if 'pythia' in self.maddm_card['indirect_flux_source_method'] or self.maddm_card['sigmav_method'] != 'inclusive'  :
                logger.info('Calculating cosmic rays fluxes using pythia8 gamma rays spectrum')
            elif 'PPPC' in self.maddm_card['indirect_flux_source_method'] and self.maddm_card['sigmav_method'] == 'inclusive':
                logger.info('Calculating cosmic rays fluxes using the spectra from the PPPC4DMID Tables')
            self.calculate_fluxes_source_n()
        # FF saving PPPC spectra in the output folder (py8 spectra are copied from the Indirect fodler already)                                                  
        # if 'PPPC' in self.maddm_card['indirect_flux_source_method']:
        #    self.save_PPPC_spectra()
            

        #print 'FF save?' , self.maddm_card['save_output']
        #print 'FF new: profile, method, function', self.maddm_card['dm_profile'] , self.maddm_card['prop_method'] , self.maddm_card['halo_funct']       


        # ****** Calculating Fluxes Earth                                                                                                                                        
            if 'earth' in self.mode['indirect']: 
                 self.neu_oscillations() # neutrinos oscillations
                  
                 if 'PPPC' in self.maddm_card['indirect_flux_earth_method']:
                     self.read_PPPC_positrons_earth()
                 elif 'dragon' in self.maddm_card['indirect_flux_earth_method']:
                     logger.info('Calling DRAGON for positrons and antiprotons propagation')
                     self.run_Dragon()                                                                                                  

              
                    
    def run_pythia8_for_flux(self):
        """ compile and run pythia8 for the flux"""
        
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        
        # Spectra produced by Pythia8
        # FF maybe not used? check
        self.spectra = ['ex','gx','nuex','numux','nutaux','px','restx']  

        #compile the pythia8 script
        if not os.path.exists(pjoin(self.dir_path,'bin','internal','main101')):
            if not hasattr(self, 'mg5'):
                self.run_mg5('')
            py8 = self.mg5.options['pythia8_path']
            files.cp(pjoin(py8, 'share','Pythia8','examples','Makefile'), 
               pjoin(self.dir_path,'bin','internal'))
            files.cp(pjoin(py8, 'share','Pythia8','examples','Makefile.inc'), 
               pjoin(self.dir_path,'bin','internal'))
            try:
                misc.compile(['main101'], cwd=pjoin(self.dir_path,'bin','internal'))
            except MadGraph5Error,e:
                print e
                logger.critical('Indirect detection, py8 script can not be compiled. Skip flux calculation')
                return

        
        # Now write the card.
        if not self.in_scan_mode: 
            pythia_cmd_card = pjoin(self.dir_path, 'Indirect' ,'Source', "spectrum.cmnd")
                                            
            # Now write Pythia8 card
            # Start by reading, starting from the default one so that the 'user_set'
            # tag are correctly set.
            PY8_Card = Indirect_PY8Card(pjoin(self.dir_path, 'Cards', 
                                                    'pythia8_card_default.dat'))
            PY8_Card['Main:spareParm1'] = mdm
            PY8_Card.read(pjoin(self.dir_path, 'Cards', 'pythia8_card.dat'),
                                                                  setter='user')
            PY8_Card.write(pythia_cmd_card, 
                           pjoin(self.dir_path, 'Cards', 'pythia8_card_default.dat'),
                            direct_pythia_input=True)

            
        run_name = self.me_cmd.run_name
        # launch pythia8
        pythia_log = pjoin(self.dir_path , 'Indirect', 'Events', run_name, 'pythia8.log')

        # Write a bash wrapper to run the shower with custom environment variables
        wrapper_path = pjoin(self.dir_path,'Indirect', 'Events',run_name,'run_shower.sh')
        wrapper = open(wrapper_path,'w')
        shell = 'bash' if misc.get_shell_type() in ['bash',None] else 'tcsh'
        shell_exe = None
        if os.path.exists('/usr/bin/env'):
            shell_exe = '/usr/bin/env %s'%shell
        else: 
            shell_exe = misc.which(shell)
        if not shell_exe:
            raise self.InvalidCmd('No s hell could be found in your environment.\n'+
              "Make sure that either '%s' is in your path or that the"%shell+\
              " command '/usr/bin/env %s' exists and returns a valid path."%shell)
        
        pythia_main = pjoin(self.dir_path,'bin','internal','main101')
        exe_cmd = "#!%s\n%s" % (shell_exe, pythia_main)
        wrapper.write(exe_cmd)
        wrapper.close()

        # Set it as executable
        st = os.stat(wrapper_path)
        os.chmod(wrapper_path, st.st_mode | stat.S_IEXEC)
        # No need for parallelization ?
        logger.info('Follow Pythia8 shower by running the '+
                'following command (in a separate terminal):\n    tail -f %s' % pythia_log)
        
        options = self.me_cmd.options
        if options['run_mode']==1:
            cluster = self.me_cmd.cluster    
            ret_code = cluster.launch_and_wait(wrapper_path, 
                    argument= [], stdout= pythia_log, stderr=subprocess.STDOUT,
                                  cwd=pjoin(self.dir_path,'Indirect','Events',run_name))
        else:                
            ret_code = misc.call(wrapper_path, stdout=open(pythia_log,'w'), stderr=subprocess.STDOUT,
                                  cwd=pjoin(self.dir_path,'Indirect','Events',run_name))

        #### WORK ON CLUSTER !!!

        ### FF Fix to make it work on cluster!
        #if ret_code != 0:
        #    raise self.InvalidCmd, 'Pythia8 shower interrupted with return code %d.\n'%ret_code+ ' You can find more information in this log file:\n%s' % pythia_log

        '''                                           
        for sp,k in self.Spectra.spectra_id.iteritems() :
            sp_name = sp + '_lhe.dat'
        
            sp_out = pjoin(self.dir_path , 'Indirect', 'Events', run_name, sp_name )
            #raw_input('FF :')
            out_dir = pjoin(self.dir_path,'output', sp_name.replace('lhe','pythia8') )
            print 'FF copying pythia spectra'
            shutil.copyfile(sp_out , out_dir )
        '''


        ### FF I move the spectra at source produced by Pythia8 to the output directory
        # Save the values of x (same for all spectra) and the values of the spectra from the Pytha file 

    # reading the spectra from the pythia8 output
    def read_py8spectra(self):
        run_name = self.me_cmd.run_name
        for sp,k in self.Spectra.spectra_id.iteritems() :
            
            sp_name = sp + '_lhe.dat'
            out_dir = pjoin(self.dir_path,'Indirect', 'Events', run_name, sp_name )
            if sp == 'gx': # FF the x values are the same for any spectra
                  x = np.loadtxt(out_dir , unpack = True )[0]
                  self.Spectra.spectra['x'] = [ np.power(10,num) for num in x]     # from log[10,x] to x
            self.Spectra.spectra[k] = np.loadtxt(out_dir , unpack = True )[1].tolist()                                  
       


    # FF this function reads the spectra from the PPPC tables from each annihilation channel, and adds them up according to the BR 
    def read_PPPCspectra(self):
         mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

         if 'PPPC4DMID' in self.maddm_card['indirect_flux_source_method'] or 'inclusive' in self.maddm_card['sigmav_method']:
            if self.Spectra.check_mass(mdm):
               if '_ew' in self.maddm_card['indirect_flux_source_method']:
                    PPPC_source = self.Spectra.load_PPPC_source(corr = 'ew')
               else:
                    PPPC_source = self.Spectra.load_PPPC_source(corr = '')

         # FF here I combine the spectra: multiply each channel by the BR                                                                                                  
         # I create a temporary list that at each loop in the channel, retains the partial sum of each x value multiplied by the channel BR                                      
         # At the end of the loop over the channels, this is the combine spectrum (do it for all the spectra - gammas, nue, positron etc. )

         self.Spectra.spectra['x'] = PPPC_source['x']

         # print 'FF last results:', self.last_results
         available_channels = {}
         for x in self.last_results.keys():
              if 'err' not in x and 'taacsID#' in x and 'lim_' not in x:
                  available_channels[x] = x.split('_')[1] # list of available SM channel xsections

    
         self.last_results['available_channels'] = available_channels.keys()
         #print 'FF the available channels are' , available_channels.keys()
         # {'taacsID#xxdxxdb_ccx': 'ccx', 'taacsID#xxdxxdb_y0y0': 'y0y0', 'taacsID#xxdxxdb_ttx': 'ttx', 'taacsID#xxdxxdb_ssx': 'ssx', 
         # 'taacsID#xxdxxdb_uux': 'uux', 'taacsID#xxdxxdb_ddx': 'ddx', 'taacsID#xxdxxdb_bbx': 'bbx'}

         # self.Spectra.spectra_id = {'px':'antiprotons', 'gx':'gammas', 'nuex':'neutrinos_e', 'numux':'neutrinos_mu', 'nutaux':'neutrinos_tau', 'ex':'positrons'}

         # FF Check that at lest one SM channel is available
         if not any(i in self.Spectra.map_allowed_final_state_PPPC.keys() for i in available_channels.values()):
                  logger.error('No SM annihilation channel available, cannot use PPPC4DMID Tables!')
                  return
                  
         for sp, sp_t in self.Spectra.spectra_id.iteritems():
              #spec_temp = [ 0 for i in range(len(x)) ] 
              self.last_results['tot_SM_xsec'] = 0
              # self.Spectra.map_allowed_final_state_PPPC = {'qqx':'qq', 'ccx':'cc', 'gg':'gg', 'bbx':'bb', 'ttx':'tt',
              #                                             'e+e-':'ee', 'mu+mu-':'mumu', 'ta+ta-':'tautau', 'w+w-':'WW', 'zz':'ZZ', 'hh':'hh' }              
              temp_dic = {}

              '''
              for ch in self.Spectra.map_allowed_final_state_PPPC.keys() :
                  #print 'FF available_channels.values()', available_channels.values() # ['ccx', 'y0y0', 'ttx', 'ssx', 'uux', 'ddx', 'bbx']
                  #quarks = ['ddx','uux','ssx']
                  #if ch not in available_channels.values() and ch not in quarks: continue
                  if ch not in available_channels.values() : continue                                                                                         
                  key = [ k for k in  available_channels.keys() if available_channels[k] == ch ][0]
                  #print 'FF key-ch-ch_t ', key , ' ' , ch , ' ' , ch_t 
                  #print 'FF lista' , [ self.last_results[i] for i in available_channels if ch in i  ] 
                  
                  ch_BR       =  self.last_results[key]
                  if ch_BR > 0:
                      temp_dic[ch] = {}
                      #print 'FF non zero XSec ' , ch   
                      self.last_results['tot_SM_xsec'] +=  ch_BR
                      print 'FF the total SM Xsec is ' , self.last_results['tot_SM_xsec']

                      interp_spec = self.Spectra.interpolate_spectra(PPPC_source , mdm = mdm, spectrum = sp_t , channel = self.Spectra.map_allowed_final_state_PPPC[ch] )
                      temp_dic[ch]['spec'] = interp_spec
                      temp_dic[ch]['xsec'] = ch_BR 

              '''
              for CH_k in available_channels.keys():
                  CH = available_channels[CH_k]
                  #print 'FF available ' , CH 
                  if CH in self.Spectra.map_allowed_final_state_PPPC.keys():
                     ch = self.Spectra.map_allowed_final_state_PPPC[CH]  # ch is the name of the channels in the Tables

                  
                  # Mapping liht quarks to qq in the Tables
                  elif CH == 'ssx' : ch = 'qq'
                  elif CH == 'uux' : ch = 'qq' 
                  elif CH == 'ddx' : ch = 'qq'
                  else: continue
                  ch_BR =  self.last_results[CH_k]                                                                                                           
                  #print 'FF CH , ch ' , CH , ch 
                  if ch_BR > 0:                                                                                                                                              
                      temp_dic[CH] = {}                                                                                                                                      
                      #print 'FF non zero XSec ' , ch                                                                                                                      
                      self.last_results['tot_SM_xsec'] +=  ch_BR                                                                                                          
                      #print 'FF the total SM Xsec is ' , self.last_results['tot_SM_xsec']                                                                              
                      interp_spec = self.Spectra.interpolate_spectra(PPPC_source , mdm = mdm, spectrum = sp_t , channel = ch , earth = False )        
                      temp_dic[CH]['spec'] = interp_spec                                                                                                                
                      temp_dic[CH]['xsec'] = ch_BR
                  #print 'FF interp_spec ', interp_spec , ' BR  ' , ch_BR , ' ch ' , key 
                  #spec_temp = [ elem + ch_BR * interp for elem,interp in zip(spec_temp,interp_spec) ]
                  #print 'FF spec_temp' , spec_temp 
              # if there is at least one SM channel, the XSec cannot be zero
              #print 'FF check dic squarks' , temp_dic
              # FF here, adding up the specta weighthed by the BR (i.e. cross section of each channel / total SM cross section )
              sum_spec = []
              #print 'FF temp_dic', temp_dic
              if not bool(temp_dic): 
                 logger.error('There is no annihilation process into SM with cross section > 0!')
                 return 
              for num in range(0,len(temp_dic[temp_dic.keys()[0]]['spec']) ):
                  #print 'FF len' , len(temp_dic[temp_dic.keys()[0]]['spec']) , num 
                  val = 0
                  for k in temp_dic.keys():
                      #print 'FF dic[k][spec]',  temp_dic[k]['spec']
                      val = val + temp_dic[k]['spec'][num] * ( temp_dic[k]['xsec'] / self.last_results['tot_SM_xsec'] )
                      #print 'FF val ', val , ' ' , temp_dic[k]['xsec'] / self.last_results['tot_SM_xsec']  
                  sum_spec.append(val)

              #self.Spectra.spectra[sp_t] = [ val/ self.last_results['tot_SM_xsec'] for val in spec_temp ]   
              self.Spectra.spectra[sp_t] = sum_spec
              #print 'FF self.Spectra.spectra[sp_t]' , self.Spectra.spectra[sp_t] 

         #print self.Spectra.spectra

    
    def read_PPPC_positrons_earth(self):

         PROF = self.maddm_card['dm_profile']
         PROP = self.maddm_card['prop_method']
         HALO = self.maddm_card['halo_funct']
         logger.debug('FF profile,prop,halo ***  ' + PROF + ' ' + PROP + ' ' + HALO ) 

         mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

         if 'PPPC4DMID' in self.maddm_card['indirect_flux_earth_method'] or self.maddm_card['sigmav_method'] == 'inclusive':
            if self.Spectra.check_mass(mdm):
               PPPC_earth = self.Spectra.load_PPPC_earth(prof = PROF) 

         E = 10**np.array(PPPC_earth['x']) # converting Log10(E) in E ; will save fluxes as E vs dN/dE
         self.Spectra.flux_earth_positrons['e'] = E

         print 'FF arrays' , PPPC_earth['x'] ,E

         channels = self.last_results['available_channels']
         #print 'FF channels earth' , channels

         # filling a dicitonary with all non -zero cross section positrons fluxes from SM annihilations channels
         temp_dic = {}
         for CH in channels: # e.g.  ccx, ttx ecc.
            C = CH.split('_')[1]
            if C in self.Spectra.map_allowed_final_state_PPPC.keys():
               ch = self.Spectra.map_allowed_final_state_PPPC[C]  # ch is the name of the channels in the Tables                                                           
                  # Mapping light quarks to qq in the Tables                                                                                                                      
            elif C == 'ssx' : ch = 'qq'
            elif C == 'uux' : ch = 'qq'
            elif C == 'ddx' : ch = 'qq'
            else: continue
            #print 'FF CH , ch ' , CH , ch 
            ch_BR = self.last_results[CH]
            if ch_BR > 0:
                temp_dic[C] = {}
                # FIX interpolation !!!
                #interp_spec = self.Spectra.interpolate_spectra(PPPC_earth , mdm = mdm, spectrum = 'positrons' , channel = ch )
                interp_spec = self.Spectra.interpolate_spectra(PPPC_earth, mdm = mdm, spectrum = 'positrons' , channel = ch , \
                                    earth = True , prof = PROF , prop = PROP , halo_func = HALO)
                
                temp_dic[C]['spec'] = interp_spec
                temp_dic[C]['xsec'] = ch_BR

         sum_spec = []
 
         #print 'FF temp dic', temp_dic
         for num in range(0,len(temp_dic[temp_dic.keys()[0]]['spec']) ):
            val = 0
            for k in temp_dic.keys():
                val = val + temp_dic[k]['spec'][num] * ( temp_dic[k]['xsec'] / self.last_results['tot_SM_xsec'] )
            sum_spec.append(val)

         self.Spectra.flux_earth_positrons['positrons']['dPhidlogE'] = sum_spec



    def neu_oscillations(self):
        nue   = self.Spectra.flux_source['neutrinos_e']['dNdE']
        numu  = self.Spectra.flux_source['neutrinos_mu']['dNdE']
        nutau = self.Spectra.flux_source['neutrinos_tau']['dNdE']
        P_e_mu   = 0.199
        P_e_tau  = 0.221
        P_mu_tau = 0.368

        nue_earth , numu_earth , nutau_earth = [] , [] , []

        for x_e, x_mu, x_tau in zip(nue,numu,nutau):
            nue_earth  .append( (1-P_e_mu-P_e_tau)  * x_e   + (P_e_mu) * x_mu + (P_e_tau)  * x_tau )
            numu_earth .append( (1-P_e_mu-P_mu_tau) * x_mu  + (P_e_mu) * x_e  + (P_mu_tau) * x_tau )
            nutau_earth.append( (1-P_e_tau-P_mu_tau)* x_tau + (P_e_tau)* x_e  + (P_mu_tau) * x_mu  )

        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        self.Spectra.flux_earth_neutrinos['e'  ]           =  self.Spectra.flux_source['e']
        self.Spectra.flux_earth_neutrinos['neutrinos_e'  ]['dNdE'] = nue_earth
        self.Spectra.flux_earth_neutrinos['neutrinos_mu' ]['dNdE'] = numu_earth
        self.Spectra.flux_earth_neutrinos['neutrinos_tau']['dNdE'] = nutau_earth
        
        

        #print 'FF these are the oscillated neutrinos: e , mu , tau'
        #for e,mu,tau in zip(nue_earth, numu_earth, nutau_earth):
        #    print 'FF e , mu , tau ' , e , mu , tau 
        # 'neutrinos_mu', 'x', 'neutrinos_e', 'neutrinos_tau'
        #print 'FF neutrinos', self.Spectra.spectra.keys()
 

    '''
    def save_PPPC_spectra(self, PPPC = False, out_dir = '' ):

      if not out_dir: out_dir = pjoin(self.dir_path,'output')                                                                                                    
      if PPPC:                                                                                                                                                                  
        #        ['positrons', 'gammas', 'neutrinos_mu', 'antiprotons', 'x', 'neutrinos_e', 'neutrinos_tau']                                                               
        for spec in self.Spectra.spectra.keys():                                                                                                                                
            name = spec.replace('antiprotons','p').replace('gammas','g').replace('neutrinos_','nu')                                                                              
            out_spec = open( pjoin(out_dir, name+'_PPPC4DMID.dat') , 'w' )                                                                                                     
            for x,s in zip(self.Spectra.spectra['x'], self.Spectra.spectra[spec] ):                                                                                          
                #print '%.3g  %.3e' % ( np.log10(x) , s )                                                                                                                     
                out_spec.write( ' %.4e  %.4e \n' % ( np.log10(x) , s ) )                                                                                                       
            out_spec.close() 
    ''' 


    # out dir must be the Indirect/Events/run_xx folder

    '''
    def save_spectra_fluxes(self , earth = False, what = 'spectrum_source' , out_dir = ''):

        if not out_dir: out_dir = pjoin(self.dir_path,'output')        

        for spec in self.Spectra.spectra.keys():
                 # saving spectra
                 name = spec.replace('antiprotons','p').replace('gammas','g').replace('neutrinos_','nu')                                                                     
                 out_spec = open( pjoin(out_dir, name+'_' + what + '.dat') , 'w' )
                 if 'PPPC' in  self.maddm_card['indirect_flux_source_method']: method =  self.maddm_card['indirect_flux_source_method']
                 else: method = 'PPPC4DMID' 
                 out_spec.write ('# Spectra from ' + method + ' Tables \n')
                 out_spec.write ('# x=Ekin/mDM      dn/dlogx   ')
                 for x,s in zip(self.Spectra.spectra['x'], self.Spectra.spectra[spec] )
                      #print '%.3g  %.3e' % ( np.log10(x) , s )                                                                                                             
                     out_spec.write( ' %.4e     %.4e \n' % ( np.log10(x) , s ) )
                 out_spec.close()
        
         
        
        if earth:
     ''' 
    def calculate_fluxes_source_n(self):
        #print 'FF the spectra are', self.Spectra.spectra                                                                                                                     
        # self.Spectra.spectra = = {'x':[] , 'antiprotons':[], 'gammas':[], 'neutrinos_e':[], 'neutrinos_mu':[], 'neutrinos_tau':[], 'positrons':[] }
        cr_spectra = self.Spectra.spectra.keys()
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        if self.mode['indirect'].startswith('flux'):

             x = np.array (self.Spectra.spectra['x'] )
             E = x*mdm
             self.Spectra.flux_source['e'] = E

             for spec in cr_spectra:
                   if spec == 'x': continue
                   dNdlogx = self.Spectra.spectra[spec] # as it is saved from PPPC or Pythia
                   dNdE = dNdlogx / (E*2.30259)          
                   self.Spectra.flux_source[spec]['dNdE'] = dNdE 
                   Phi = self.Phi(chan= spec)

                   '''
                   Only needed for dPhidE
                   for energy in energies:
                       self.dNdE_dPhidE( energy,  x=x, dndlogx=dndlogx)
                       dNdE   = self.dNdE_int
                       dPhidE = self.dPhidE
                       #print 'FF res dNdE', dNdE , dPhidE
                       dNdE_a  .append(dNdE)
                       dPhidE_a.append(dPhidE)
                   self.Spectra.flux_source[chan]['dNdE']   = dNdE_a
                   self.Spectra.flux_source[chan]['dPhidE'] = dPhidE_a
                   '''

                   self.last_results['flux_%s' % spec] = Phi 




    '''
    def calculate_fluxes_source(self):
        
        #print 'FF the spectra are', self.Spectra.spectra 
        CR = {'gammas':'g'      , 'neutrinos_e':'nue' , 'neutrinos_mu':'numu' , 'neutrinos_tau':'nutau' , 'antiprotons':'p' , 'positrons':'e' }

        self.CR_fluxes_source = { 'Energies':'', 'dPhidE_g':'', 'dPhidE_nue':'', 'dPhidE_numu':'' , 'dPhidE_nutau':'', 'dPhidE_e':'', 'dPhidE_p':'' }

        if self.mode['indirect'].startswith('flux'):           
           #Here we need to skip this part if the scan is being conducted because                                                                                              
           #the value of dark matter mass could be 'scan: ...'                                                                                                            
             mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
               #run_name = self.me_cmd.run_name
             npts = self.maddm_card['npts_for_flux']
             energies = [1.0*(mdm)/ npts * n for n in range(npts)] ##  logscale is better, to be changed...                                                                      
             self.CR_fluxes['Energies'] = energies
             # here: add dn/de spectra!

             for chan,chan_n in CR.iteritems():
                   phis= []
                   for energy in energies:
                       phis.append(self.dPhidE(energy, channel=chan))
                                 
                   self.CR_fluxes_source['dPhidE_'+chan_n] = phis
                   # saving energy-integrated flux
                   self.last_results['flux_%s' % chan_n] = self.Phi(chan=chan)
                   

                   #if chan_n == 'g': chan_n = 'gamma'
                   #flux_filename = pjoin(self.dir_path,'output', 'dPhidE_%s.txt' %chan_n)
                     #flux_filename = pjoin(self.dir_path, 'output', 'dPhidE_%s.txt' %chan)                                                                                

                   # FFFFFFF

                   #aux.write_data_to_file(energies, phis, filename=flux_filename, header='# Energy [GeV]    Differential Flux [GeV^-1 cm^-2 s^-1 sr^-1]')
                     ## FF This is the total integrated flux for neutral particles
                     #output_name.append('flux_%s' % chan)                                                                                                                
                     #for chan,chan_n in cr_names.iteritems():
                     #   Espectrum= []
                     #   for energy in energies:
                     #     Espectrum.append(self.dNdx(energy, channel=chan))

                         # FF replace with the correct way of saving outputs!
#                        dNde_filename = pjoin(self.dir_path,'output', 'dNdE_%s.txt' %chan_n)
#                        aux.write_data_to_file(energies, Espectrum, filename=dNde_filename, header='# E_kin [GeV]   dNdE [GeV^-1]')
    '''

    # dndlogx for each differen channel
    def dNdE_dPhidE(self, energy, channel='' ):
        # input: values of x=Ekin/mDM and dn/dlogx (as PPPC tables)
        # - transforms into E and dN/dE
        # - gives the differential flux dPhi/dE averaged over the opening angle
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        E = self.Spectra.flux_source['e']
        dndlogx = np.array( self.Spectra.spectra[channel] )
        dNdE    = dndlogx /(E*2.30259)
 
        # up to here only conversion from dN/dlogx to dN/dE

        sigv = self.last_results['taacsID']
        halo_profile = self.maddm_card['halo_profile']
        jfact = self.maddm_card['jfactors'][halo_profile]

        dNdE_int = np.interp(energy, E, dNdE ) # interpolated value for an energy of value = interp.

        # *** generic differential flux dPhi/dE
        #dPhidE = 1.0/(self.norm_Majorana_Dirac() *2* math.pi*mdm*mdm)*sigv*jfact * dNdE_int # diff.flux for the energy == interp. 
        
        return dNdE_int


    def Phi(self, chan= ''):

        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        sigv = self.last_results['taacsID']
        halo_profile = self.maddm_card['halo_profile']
        jfact = self.maddm_card['jfactors'][halo_profile]

        npts =  self.maddm_card['npts_for_flux']
        grid = [de*2*mdm/npts for de in range(0, npts+1)]
            #CHECK HERE THAT THE JACOBIAN FROM X TO E IS CORRECT                                                                                                                
        #dNdE, dPhidE = self.dNdE_dPhidE(interp, x=x, dndlogx=dndlogx)                                                                                                         
        integrate_dNdE = aux.integrate( self.dNdE_dPhidE  , grid , channel = chan ) # the last are the arguments of the function                                         
        #print 'FF integrate_dNdE' , integrate_dNdE                                                                                                                             
        phi = 1.0/(self.norm_Majorana_Dirac()*2*math.pi*mdm*mdm)*sigv*jfact*integrate_dNdE                                                                                     
        #print 'FF phi is', phi                                                                                                                                               
        return phi







    
    '''
    def dPhidE_n(self, energies,  x, dndlogx ):
        """generic differential flux from dSPhs, channel can be photons or neutrinos"""
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        sigv = self.last_results['taacsID']
        halo_profile = self.maddm_card['halo_profile']
        jfact = self.maddm_card['jfactors'][halo_profile]

        dphi = 1.0/(self.norm_Majorana_Dirac() *2* math.pi*mdm*mdm)*sigv*self.dNdx(energy, channel)*jfact          ### factor 1/4 for majorana and 1/8 for dirac                
        # expression below for dphi is just to check                                                                                                                             
        #dphi = self.dNdx(energy, channel)                                                                                                                                      
        return dphi
    '''




    '''
    # OLD VERSION reading the spectra for each channel
    def dNdx(self, x, channel=''):
        ## changes by chiara below (i.e. changed name of variables and corrected a typo dNdx instead of dNdE and changed path to files)
        X, dndlogx = np.array(self.Spectra.spectra['x']) , np.array(self.Spectra.spectra[channel])
        #print 'FF x , dndlogx', X , dndlogx             
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        # THIS IS TO TEST
        #for i in range(len(dNdx_x)):
        #    print 10**(dNdx_x[i])*mdm, dNdx_y[i]/(10**(dNdx_x[i])*mdm*2.30259)  ## 2.30259 = Log10; 

        # FF changed this since the spectra are stored as function of "x" and not "Log10[x]" anymore
        # return np.interp(x, mdm*np.power(10,dNdx_x),dNdx_y/(np.power(10,dNdx_x)*mdm*2.30259))   ## dNdx_x = E /mDM; dNdlogx = dNdx_y
        return np.interp(x, mdm* X, dndlogx /(X *mdm*2.30259) )  ## dNdx_x = E /mDM; dNdlogx = dNdx_y                                                
    '''
        
    def norm_Majorana_Dirac(self):
        # (0=own anti) -> this works when there is the factor 2pi in the denominator
        dirac_maj = self.param_card.get_value('qnumbers ' + str (self.proc_characteristics['dm_candidate'][0] ), 4) 
        if   dirac_maj == 0:
            logger.debug('FF DM is a Majorana particle')
            return 4
        elif dirac_maj == 1: 
            logger.debug('FF DM is a Dirac particle')
            return 8 
    '''
    OLD version
    def dPhidE(self,  energy, channel=''):
        """generic differential flux from dSPhs, channel can be photons or neutrinos"""

        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        sigv = self.last_results['taacsID']
        halo_profile = self.maddm_card['halo_profile']
        jfact = self.maddm_card['jfactors'][halo_profile]
        
        dphi = 1.0/(self.norm_Majorana_Dirac() *2* math.pi*mdm*mdm)*sigv*self.dNdx(energy, channel)*jfact          ### factor 1/4 for majorana and 1/8 for dirac
        # expression below for dphi is just to check
        #dphi = self.dNdx(energy, channel)
         
        return dphi
    '''
    '''
    def Phi(self, x, dndlogx):                                                                                                                     
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])                                                                            
        sigv = self.last_results['taacsID']                                                                                                                                  
        halo_profile = self.maddm_card['halo_profile']                                                                                                                      
        jfact = self.maddm_card['jfactors'][halo_profile]                                                                                                                    
        
        npts =  self.maddm_card['npts_for_flux']                                                                                                                               
        grid = [de*2*mdm/npts for de in range(0, npts+1)]                                                                                                                     
            #CHECK HERE THAT THE JACOBIAN FROM X TO E IS CORRECT           
        #dNdE, dPhidE = self.dNdE_dPhidE(interp, x=x, dndlogx=dndlogx)                                                                                                       
        integrate_dNdE = aux.integrate( self.dNdE_dPhidE  , grid , x=x, dndlogx=dndlogx ) # the last are the arguments of the function                                       
                           
        #print 'FF integrate_dNdE' , integrate_dNdE 
        phi = 1.0/(self.norm_Majorana_Dirac() *math.pi*mdm*mdm)*sigv*jfact*integrate_dNdE                                                                                         
     
        #print 'FF phi is', phi 
        return phi 
    '''


    '''
    def Phi(self,chan='' ):
        if not self.last_results['taacsID']:
            logger.error('You can not calculate the flux before calculating <sigmav>!')
            return -1.0
        else:
            #param_card = param_card_mod.ParamCard(self.dir_path, 'Cards', 'param_card.dat')
            mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
            sigv = self.last_results['taacsID']
            halo_profile = self.maddm_card['halo_profile']
            jfact = self.maddm_card['jfactors'][halo_profile]
            
            npts =  self.maddm_card['npts_for_flux']
            grid = [de*2*mdm/npts for de in range(0, npts+1)]
            
            #CHECK HERE THAT THE JACOBIAN FROM X TO E IS CORRECT
            integrate_dNdE = aux.integrate(self.dNdx, grid, channel=chan)  # 1/mdm

            integrate_dNdE = aux.integrate(self.dNdx, grid)  # 1/mdm                                                                                            
            #logger.debug('sigmav: %.5e' % sigv)
            #logger.debug('mdm: %.5e' % mdm)
            #logger.debug(jfact)
            
            # FF is it ok to put twice the Dirac-Majorana factor ??? It is already included in the dNdX caculation
            phi = 1.0/(self.norm_Majorana_Dirac *math.pi*mdm*mdm)*sigv*jfact*integrate_dNdE
            
            return phi
    '''


    '''
    # the channels is not needed anymore
    def Phi(self,chan = ''):
        if not self.last_results['taacsID']:
            logger.error('You can not calculate the flux before calculating <sigmav>!')
            return -1.0
        else:
            #param_card = param_card_mod.ParamCard(self.dir_path, 'Cards', 'param_card.dat')                                                                                       
            mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
            sigv = self.last_results['taacsID']
            halo_profile = self.maddm_card['halo_profile']
            jfact = self.maddm_card['jfactors'][halo_profile]

            npts =  self.maddm_card['npts_for_flux']
            grid = [de*2*mdm/npts for de in range(0, npts+1)]
            #CHECK HERE THAT THE JACOBIAN FROM X TO E IS CORRECT                                                                                                                 
            integrate_dNdE = aux.integrate(self.dNdx, grid, channel = chan)  # 1/mdm                                                            
            #logger.debug('sigmav: %.5e' % sigv)                                                                                                                                 
            #logger.debug('mdm: %.5e' % mdm)                                                                                                                                     
            #logger.debug(jfact)                                                                                                                                                 
            phi = 1.0/(self.norm_Majorana_Dirac()*2*math.pi*mdm*mdm)*sigv*jfact*integrate_dNdE
            return phi
    '''
    # add output directory as input!

    def run_Dragon(self, out_dir = ''):
 
        point_name = self.last_results['point_number']
        
        dragon_dir = '/scratch/federico/MadDM/DRAGON/DRAGON-master/'

        if not os.path.exists(pjoin(dragon_dir,'DRAGON') ):
           logger.error('The DRAGON executable cannot be find! Will not run positrons and antiprotons propagation!')
           return
       
        template_card = pjoin(MDMDIR,'Templates','Cards','dragon_card.xml')
        dr_path = 'path to Dragon installation - read from where?'
        mDM = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        if self.maddm_card['sigmav_method'] =='inclusive': sigv = self.last_results['tot_SM_xsec']
        else:                                              sigv = self.last_results['taacsID']

        # !!! Fix this when calling the function in the loop (scan mode) !!!
        if not out_dir:
           out_dir = pjoin(self.dir_path, 'output')

        dragon_input = pjoin( self.dir_path, 'output' , str(point_name) + '_DRAGON.xml')

        def write_dragon_input(template= template_card , mdm = mDM , sigmav = sigv , dragon_input = dragon_input ):
            # convert spectra in dndx in dnde
            #print 'FF self.Spectra.spectra keys', self.Spectra.spectra.keys()
            energy , positrons, antiprotons = self.Spectra.flux_source['e'] , self.Spectra.flux_source['positrons']['dNdE'] , \
                                                                              self.Spectra.flux_source['antiprotons']['dNdE']
                  
            out_pos  = open(pjoin(self.dir_path, 'output','positrons_dndne.txt') , 'w')
            out_anti = open(pjoin(self.dir_path, 'output','antiprotons_dndne.txt') ,'w')
            print 'FF directories ' , pjoin(self.dir_path, 'output')

            for en, pos, anti in zip (energy, positrons, antiprotons):
                out_pos .write(str(en)+ ' ' + str(pos)  +'\n')
                out_anti.write(str(en)+ ' ' + str(anti) +'\n')
            out_pos .close()
            out_anti.close()

            xml_in = open(dragon_input, 'w') 
            inp = open(template,'r') # read template cardq
            for line in inp.readlines():
                line = line.replace('MadDM_dm_mass'    , str(mdm)    )
                line = line.replace('MadDM_sigmav'     , str(sigmav) )
                line = line.replace('MadDM_Positrons'  , pjoin(self.dir_path, 'output', 'positrons_dndne.txt'  ) )
                line = line.replace('MadDM_Antiprotons', pjoin(self.dir_path, 'output', 'antiprotons_dndne.txt') )
                xml_in.write(line)
            print 'FF written Dragon input card in ' , dragon_input
            xml_in.close()
            

        write_dragon_input(template= template_card , mdm = mDM , sigmav = sigv )
        os.chdir(dragon_dir)
        os.system('./DRAGON ' + dragon_input )
        os.chdir(self.dir_path)
        os.remove(pjoin(self.dir_path, 'output','positrons_dndne.txt') )
        os.remove(pjoin(self.dir_path, 'output','antiprotons_dndne.txt') )


    def write_scan_output(self, out_path = '', keys = '', header = False):


        # writing the parameters                                                                                                                                                 
        if out_path and header:
           nice_keys = []
           for k in keys:
               k = k.replace('taacsID#','')
               k = k.replace('taacsID','tot_Xsec')
               nice_keys.append(k)

           summary = open(out_path, 'w')                                                                                                                            

           for k in nice_keys:
                ind = nice_keys.index(k) + 1
                if ind <=9: ind = '0'+str(ind)
                summary.write( '# [' + str(ind) + ']' + ' : ' + k + '\n' )

           summary.write('\n\n\n')
           summary.close()

        elif (out_path and not header):
            s = '\t'
            summary = open(out_path, 'a+')
            summary.write('{:9d}'.format(int(self.last_results['run'])) + s)

            for k in keys:
                num = self.form_n( self.last_results[k] )
                summary.write(num + s)
            summary.write('\n')
            summary.close()


    
    def print_results(self):
    
        """ print the latest results """
        #omega_min = 0
        #omega_max = 0.12

        # Determine xsi for thermal scenarios with DM < DM tot.                
        dm_scen = 0
        if self.last_results['xsi'] > 0 and self.last_results['xsi'] < 1:
             xsi  = self.last_results['xsi']
             dm_scen = True
        else: 
             xsi = 1
             dm_scen = False

        xsi2 = xsi**2

        mdm= self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        pass_message  = '%s ALLOWED  %s' % (bcolors.OKGREEN, bcolors.ENDC)
        fail_message  = '%s EXCLUDED %s' % (bcolors.FAIL, bcolors.ENDC)
        nolim_message = '%s NO LIMIT %s' % (bcolors.GRAY, bcolors.ENDC)

        #skip this if there is a sequential scan going on.

        '''
        if not self.param_card_iterator:
            pass_relic = pass_message if self.last_results['Omegah^2'] < self.limits._oh2_planck else fail_message
        else:
            pass_relic         = ''
        '''
        pass_relic = pass_message if self.last_results['Omegah^2'] < self.limits._oh2_planck else fail_message                                                                 

        #if omega_min < self.last_results['omegah2'] < omega_max:
        #    fail_relic_msg = ''
        #else:
        #    fail_relic_msg = '%s Model excluded (relic not in range [%s,%s])%s' %\
        #                      (bcolors.FAIL, omega_min, omega_max,bcolors.ENDC)

        omega    = self.last_results['Omegah^2']
        x_f      = self.last_results['x_f']
        sigma_xf = self.last_results['sigmav(xf)']
        xsi      = self.last_results['xsi']        

        logger.info("MadDM Results", '$MG:color:BLACK')
        if dm_scen: 
              logger.info("Define xsi = Relic density/Planck measurement for thermal scenarios\n")
              logger.info("Rescaling theory prediction for xsi(direct det.) and xsi^2(indirect det.)\n")

        if self.mode['relic']:
            logger.info("\n***** Relic Density")

            logger.info( self.form_s('Relic Density') + '=\t' + self.form_s(self.form_n(omega    ) ) )
            logger.info( self.form_s('x_f'          ) + '=\t' + self.form_s(self.form_n(x_f      ) ) )
            logger.info( self.form_s('sigmav(xf)'   ) + '=\t' + self.form_s(self.form_n(sigma_xf ) ) )
            logger.info( self.form_s('xsi'          ) + '=\t' + self.form_s(self.form_n(xsi      ) ) )


#            logger.info(self.form_s('Relic density') +'= '+ '\t' + self.form_n(omega) + '\t'+ pass_relic)
#            logger.info(self.form_s('x_f')           +'= '+ '\t' + self.form_n(x_f)                     )
#            logger.info(self.form_s('sigma(xf)')     +'= '+ '\t' + '%.2e GeV^-2 = %.2e cm^3/s', self.last_results['sigmav(xf)'],self.last_results['sigmav(xf)']*GeV2pb*pb2cm3 )
#            logger.info(self.form_s('sigma(xf)')     +'= '+ '\t' + ' = %.2e cm^3/s', self.last_results['sigmav(xf)']*GeV2pb*pb2cm3 )

            #logger.info('   x_f            : %.2f', self.last_results['x_f'])
            #logger.info('   sigmav(xf)     : %.2e GeV^-2 = %.2e cm^3/s', self.last_results['sigmav(xf)'],self.last_results['sigmav(xf)']*GeV2pb*pb2cm3)

        if self.mode['direct']:
            sigN_SI_p , sigN_SI_n = self.last_results['sigmaN_SI_p'] , self.last_results['sigmaN_SI_n']
            sigN_SD_p , sigN_SD_n = self.last_results['sigmaN_SD_p'] , self.last_results['sigmaN_SD_n']

            units = self.last_results['GeV2pb*pb2cm2']
            direct_names = [ { 'n':'SigmaN_SI_p','sig':sigN_SI_p*units , 'lim':self.limits.SI_max(mdm)     , 'exp':'Xenon1ton' },
                             { 'n':'SigmaN_SI_n','sig':sigN_SI_n*units , 'lim':self.limits.SI_max(mdm)     , 'exp':'Xenon1ton' },
                             { 'n':'SigmaN_SD_p','sig':sigN_SD_p*units , 'lim':self.limits.SD_max(mdm,'p') , 'exp':'Pico60' },
                             { 'n':'SigmaN_SD_n','sig':sigN_SD_n*units , 'lim':self.limits.SD_max(mdm,'n') , 'exp':'Lux2017' } ]

            self.last_results['direct_results'] = direct_names

            logger.info('\n***** Direct detection [cm^2]: ')  
            for D in direct_names:
                    th_cross = D['sig']* xsi 
                    ul = D['lim']
                    mess_th       = self.det_message_screen(th_cross, D['lim'] )
                    mess_alldm    = self.det_message_screen(D['sig'], D['lim'] )
                    self.print_ind(D['n'],th_cross , D['sig'], ul,  thermal=dm_scen ,direc = False , exp=D['exp']) # leave False to avout GeV^-2

                    #if dm_scen: # thermal scenario 
                    #      logger.info(self.form_s(D['n']+' (Thermal)') +'= '+'\t'+ self.form_n(th_cross)+' GeV^-2 = '+self.form_n(th_cross*GeV2pb*pb2cm2)+ \
                    #                 ' cm^2'+'\t'+ mess_th    + '\t'      + self.form_s(' (All DM)')        +'= '+'\t'+ self.form_n(D['sig'])+' GeV^-2 = '+self.form_n(D['sig']                  #                   *GeV2pb*pb2cm2)+' cm^2'+'\t'+ mess_alldm + '\t'      +'ul ='+ self.form_n(D['lim']) + ' cm^2'+ '\t'+ self.form_s(D['exp']) )
                    #else: 
                    #      logger.info(self.form_s(D['n']+' All DM') +'= '+'\t'+ self.form_n(D['sig'])+' GeV^-2 = '+self.form_n(D['sig']*GeV2pb*pb2cm2)+' cm^2'+'\t'+ mess_alldm  +\
                    #      '\t'      +'ul= '+ self.form_n(D['lim']) + ' cm^2'+ '\t'+ self.form_s(D['exp']) )
     
            '''
            logger.info(' sigmaN_SI_p    th: %.2e GeV^-2 = %.2e cm^2   \t ul: %.2e cm^2 %s',xsi*sigN_SI_p, xsi*sigN_SI_p*GeV2pb*pb2cm2,self.limits.SI_max(mdm), pass_dd_si_proton)
            logger.info(' sigmaN_SI_n    th: %.2e GeV^-2 = %.2e cm^2   \t ul: %.2e cm^2 %s',xsi*sigN_SI_n, xsi*sigN_SI_n*GeV2pb*pb2cm2,self.limits.SI_max(mdm), pass_dd_si_neutron)
            logger.info(' sigmaN_SD_p    th: %.2e GeV^-2 = %.2e cm^2   \t ul: %.2e cm^2 %s',xsi*sigN_SD_p, xsi*sigN_SD_p*GeV2pb*pb2cm2,self.limits.SD_max(mdm,'p'),pass_dd_sd_proton
            logger.info(' sigmaN_SD_n    th: %.2e GeV^-2 = %.2e cm^2   \t ul: %.2e cm^2 %s',xsi*sigN_SD_n, xsi*sigN_SD_n*GeV2pb*pb2cm2,self.limits.SD_max(mdm,'n'),pass_dd_sd_neutron)          '''

         
#        if self.mode['direct'] == 'directional':
#            logger.info(' Nevents          : %i', self.last_results['Nevents'])
#            logger.info(' smearing         : %.2e', self.last_results['smearing'])
        
        if self.mode['capture']:
            logger.info('\n capture coefficients: ')
            detailled_keys = [k for k in self.last_results.keys() if k.startswith('ccap')]
            for key in detailled_keys:
                logger.info(' %s            : %.2e 1/s' % (key, self.last_results[key]))


        if self.mode['indirect']:
            logger.info('\n****** Indirect detection [cm^3/s]: ')

            detailled_keys = [k for k in self.last_results.keys() if k.startswith('taacsID#')]
            v = self.maddm_card['vave_indirect']
            logger.info('<sigma v> method: %s ' % self.maddm_card['sigmav_method'] ) 
            logger.info('DM particle halo velocity: %s ' % self.limits._id_limit_vel['ttx']) # the velocity should be the same
            if len(detailled_keys)>0:
                logger.info('FF  use the generic qq limits for light quarks (u,d,s)' )

                for key in detailled_keys:
                    #logger.info(key)
                    clean_key_list = key.split("#")
                    #logger.info(clean_key_list)
                    clean_key = clean_key_list[1] #clean_key_list[0]+"_"+clean_key_list[1]
                    #print 'FF for the method ', self.mode['indirect'] , ' the clean_key is: ',  clean_key 
                    finalstate = clean_key.split("_")[1]
                    f_original = clean_key.split("_")[1]
                    if 'ss' in finalstate or 'uu' in finalstate or 'dd' in finalstate:
                        finalstate = 'qqx'
                    ## The inclusive method does not give the sigma for the BSM states production
                    s_alldm = self.last_results[key]
                    s_thermal  = s_alldm*xsi2

                    if s_alldm <= 10**(-100): 
                       logger.info('Skipping zero cross section processes for ' + f_original )
                       continue 
                    if finalstate in self.limits._allowed_final_states :
                        s_ul   = self.limits.ID_max(mdm, finalstate)

                        self.last_results['Fermi_lim_'+key] = self.limits.ID_max(mdm, finalstate)
                    else:  s_ul   = '-1'
                    #here check if the cross section is allowed. Do this channel by channel
                    # Make sure that the velocity at which the limit is evaluated matches the dm velocity in the calculation.
                    if not self.param_card_iterator:

                      if finalstate not in self.limits._allowed_final_states :
                          #message = nolim_message
                          self.print_ind(clean_key, s_thermal, s_alldm, 0,  thermal= dm_scen)

                          #logger.info('     %s \t sigmav(th): %.3e \t sigmav(ul): %s \t [cm^3/s] \t  %s' % (clean_key, s_theo, s_ul, message))

                      else:
                          if (v == self.limits._id_limit_vel[finalstate]):
                             #s_thermal  = s_alldm*xsi2
                             #mess_ther  = self.det_message_screen(s_therm,ul)
                             #mess_alldm = self.det_message_screen(s_alldm,ul)
                             #self. print_ind(self,what, sig_th, sig_alldm, ul,  thermal=False)
#                                if self.last_results[key] < ( s_ul ):
#                                    message = pass_message 
#                                elif self.last_results[key] > ( s_ul ):
#                                    message = fail_message
#                                else: message = nolim_message
                             self.print_ind(clean_key, s_thermal, s_alldm, s_ul,  thermal= dm_scen)

                          else:
                                message = '%s Sigmav/limit velocity mismatch %s' %(bcolors.GRAY, bcolors.ENDC)

                        #self.print_ind(self, clean_key, s_thermal, sig_alldm, s_ul,  thermal= dm_scen)
                        #logger.info('     %s \t  sigmav(th): %.2e \t sigmav(ul): %.3g \t [cm^3/s] \t  %s' % (clean_key, s_theo, s_ul, message))

            #print 'FF last results ' , self.last_results 

            sigtot_alldm     = self.last_results['taacsID']
            sigtot_SM_alldm  = self.last_results['tot_SM_xsec']
            sigtot_th        = sigtot_alldm * xsi2
            sigtot_SM_th     = sigtot_SM_alldm * xsi2
            fermi_ul         = self.last_results['Fermi_sigmav']

            if self.maddm_card['sigmav_method']!= 'inclusive':     
                 self.print_ind('DM DM > all',sigtot_th , sigtot_alldm, fermi_ul,  thermal= dm_scen)
            else: 
                 self.print_ind('DM DM > SM SM', sigtot_SM_th , sigtot_SM_alldm, fermi_ul,  thermal= dm_scen)

            '''
            if self.maddm_card['indirect_flux_source_method'] == 'pythia8':
                if self.last_results['taacsID']*xsi2 > self.last_results['Fermi_sigmav'] and self.last_results['Fermi_sigmav'] > 0:   message = fail_message
                elif self.last_results['taacsID']*xsi2 < self.last_results['Fermi_sigmav'] : message = pass_message
                elif self.last_results['Fermi_sigmav'] < 0: message = nolim_message
                logger.info('     %s\t sigmav(th): %.2e  \t sigmav(ul): %.3g \t [cm^3/s] \t %s ' % ('DM DM > all',  self.last_results['taacsID']*xsi2 ,  self.last_results['Fe                rmi_sigmav'], message))
            '''


            if str(self.mode['indirect']).startswith('flux'):
                logger.info('*** Fluxes at source:')
                np_names = {'gammas':'g'      , 'neutrinos_e':'nue' , 'neutrinos_mu':'numu' , 'neutrinos_tau':'nutau'}
                for chan in np_names.keys():
                    logger.info( self.form_s(chan + ' Flux') + '=\t' + self.form_s(self.form_n (self.last_results['flux_%s' % chan]) ))

                    
                    #logger.info('%10s : %.3e \t particles/[cm^2 s sr]' %(chan, self.last_results['flux_%s' % chan] ))
                print '\n'


        # Internal: save the results as a numpy dictionary
        # np.save(pjoin(self.dir_path, 'output','Results'), self.last_results)


        if not self.param_card_iterator:
            self.save_summary_single(relic = True, direct = self.mode['direct'], \
                                indirect = self.mode['indirect'], 
                                fluxes_source= self.mode['indirect'].startswith('flux') if isinstance(self.mode['indirect'],str) else self.mode['indirect'], 
                                fluxes_earth = False )                  
            logger.info('Results written in ' +pjoin(self.dir_path, 'output', 'MadDM_results.txt') )

        self.save_remove_indirect_output(scan = False)


    # point number must be e.g. run_xx for madevent , or 1,2, ecc. for PPPC 
    def save_remove_indirect_output(self, scan = False, point_number= -1):

        if point_number < 10: # to uniform incusive vs madevent/reshuffling naming convention
           point_number = '0'+ str(point_number)

        ind_mode = self.mode['indirect']      
        save_switch = self.maddm_card['save_output']

        spectrum_method = self.maddm_card['indirect_flux_source_method']

        # Setting the various paths
        source_indirect = pjoin(self.dir_path,    'Indirect'      )
        events          = pjoin(source_indirect,  'Events'        )
        dir_point       = pjoin(events, 'run_'+ str(point_number) )

        # If indirect is called, then spectra_source == True by def. The other two options depends on the user choice (sigmav, flux_source, flux_earth)
        spec_source, flux_source , flux_earth = '','',''
   
        if ind_mode:
           spec_source = True
           if 'source' in ind_mode:
              flux_source = True
           elif 'earth' in ind_mode  and 'inclusive' in  self.maddm_card['sigmav_method']:
              flux_source , flux_earth = True , True

        if save_switch == 'off' and 'dragon' not in self.mode['indirect']:
           spec_source, flux_source , flux_earth  = False, False, False
            
        #print 'FF ind_mode , spec_source , flux_source , flux_earth ' , ind_mode , spec_source , flux_source , flux_earth 

        # Saving the output in the general output directory for a single point        
        if not scan: 
               out_dir = pjoin(self.dir_path, 'output')
               self.save_spec_flux(out_dir = out_dir, spec_source = spec_source, flux_source = flux_source, flux_earth = flux_earth)    
               logger.info('Output files saved in ' + pjoin(self.dir_path, 'output') )

        elif   scan:
             if not os.path.exists(source_indirect) : os.makedirs(source_indirect)
             if not os.path.exists(events)          : os.makedirs(events)
             if not os.path.exists(pjoin(self.dir_path, 'output' , 'Output_Indirect')) :
                    os.symlink(events, pjoin(self.dir_path, 'output' , 'Output_Indirect') ) # Creating the symbolic link to the Events directory (inside Indirect)
             if not os.path.exists(dir_point)       : os.makedirs(dir_point)

             out_dir = dir_point # all the various output must be saved here             


             print 'FF out_dir to be removed **************** ', out_dir

   
             if save_switch == 'off':
                 if self.maddm_card['sigmav_method'] != 'inclusive': # here, need to remove everyhting i.e. the whole run_xx folder
                    shutil.rmtree(out_dir)
                                      
             elif save_switch == 'all':
                 if self.maddm_card['sigmav_method'] == 'inclusive': # saving spectra onyl in inclusive mode since madevent/resh have their own pythia8 spectra             
                    self.save_spec_flux(out_dir = out_dir, spec_source = spec_source, flux_source = flux_source, flux_earth = flux_earth)
                 
                 elif self.maddm_card['sigmav_method'] != 'inclusive':
                    self.save_spec_flux(out_dir = out_dir, spec_source = False, flux_source = flux_source, flux_earth = flux_earth)

             elif save_switch == 'spectra':
                  if self.maddm_card['sigmav_method'] == 'inclusive':
                       self.save_spec_flux(out_dir = out_dir, spec_source = spec_source, flux_source = flux_source, flux_earth = flux_earth)
                  else:
                       print 'FF removing extra pythia files'
                       if os.path.isfile( pjoin(out_dir,'unweighted_events.lhe.gz') ):
                          os.remove( pjoin(out_dir,'unweighted_events.lhe.gz') )
                          os.remove( pjoin(out_dir,'run_shower.sh') )
       


              

    # This aux function saves the spectra/fluxes in the given directory out_dir
    def save_spec_flux(self, out_dir = '' , spec_source = False, flux_source = False, flux_earth = False):
        if not out_dir: out_dir = pjoin(self.dir_path, 'output') # default output directory if not defined

        spec_method = self.maddm_card['indirect_flux_source_method']

        if spec_source:
           x = self.Spectra.spectra['x']
           for spec in self.Spectra.spectra.keys():
               header = '# x=Ekin/mDM     dn/dlogx   ' + spec + '\t' + self.maddm_card['indirect_flux_source_method'] + ' spectra at source'
               if 'x' not in spec:
                   dndlogx = self.Spectra.spectra[spec]
                   aux.write_data_to_file(x , dndlogx  , filename = out_dir + '/' + spec + '_spectrum_source.txt' , header = header )
                     
        if flux_source:
           logger.debug('FF saving flux source')
           e = self.Spectra.flux_source['e']
           for flux in self.Spectra.flux_source.keys():
               header = '# E[GeV]   dN/dE ' + flux + '\t' + self.maddm_card['indirect_flux_source_method'] +' flux at source'
               if flux != 'e':
                   dNdE = self.Spectra.flux_source[flux]['dNdE']

                   aux.write_data_to_file(e , dNdE  , filename = out_dir + '/' + flux + '_flux_source.txt' , header = header )
                   print 'FF the flux is ', out_dir + '/' + flux + '_flux_source.txt' 

        if flux_earth:
           logger.debug('FF saving flux earth')
           e = self.Spectra.flux_earth_neutrinos['e']
           for flux in self.Spectra.flux_earth_neutrinos.keys():
               header = '# E[GeV]   dN/dE ' + flux + '\t' + self.maddm_card['indirect_flux_source_method'] +' flux at earth'
               if flux != 'e':
                   dNdE = self.Spectra.flux_earth_neutrinos[flux]['dNdE']
                   aux.write_data_to_file(e , dNdE  , filename = out_dir + '/' + flux + '_flux_earth.txt' , header = header )
                   print 'FF the neutrino earth flux is ', out_dir + '/' + flux + '_flux_source.txt'


           if 'PPPC' in self.maddm_card['indirect_flux_earth_method']: # for PPPC4DMID, I save also the positron at earth
               e = self.Spectra.flux_earth_positrons['e']
               header = '# E[GeV]   dPhi/dlogE [particles/(cm2 s sr)]  positrons \t PPPC4DMID flux at earth'
               dPhidlogE = self.Spectra.flux_earth_positrons['positrons']['dPhidlogE']
               print 'FF e , dPhidloge', e , dPhidlogE 
               aux.write_data_to_file(e , dPhidlogE  , filename = out_dir + '/positrons_flux_earth.txt' , header = header )

        misc.sprint("CODE NEED CLEANING HERE!!!")
        '''
        else: out_dir = pjoin(self.dir_path, 'output')

        if save_switch == 'off':
            if self.maddm_card['sigmav_method'] != 'inclusive': # here, need to remove everyhting i.e. the whole run_xx folder
                shutil.rmtree(dir_point)
                    
        elif save_switch == 'all':
            if self.maddm_card['sigmav_method'] == 'inclusive': # saving spectra onyl in inclusive mode since madevent/resh have their own pythia8 spectra             
                for spec in self.Spectra.spectra.keys():
                    if 'earth' in spec or 'x' in spec: 
                        continue
                    out_spec = pjoin(out_dir, spec + '_spectrum_source_'+spectrum_method+'.txt')
                    x = self.Spectra.spectra['x']
                    aux.write_data_to_file(x , self.Spectra.spectra[spec] , filename = out_spec , header = '# x=Ekin/mDM      dn/dlogx' )
                    
                    if 'flux' in self.maddm_card['sigmav_method']:
            #if self.madd_card['sigmav_method'] == 'inclusive':                                   
                        print 'FF saving spectra at source PPPC'
                    # self.Spectra.spectra
                        for spec in self.Spectra.spectra.keys():
                            x = self.CR_fluxes['Energies']
                            flux = self.CR_fluxes['dPhidE_'+  spec]
                            out_flux = pjoin(out_dir, spec + '_flux_source_'+spectrum_method+'.txt')
                            aux.write_data_to_file(x , flux  , filename = out_flux , header= h_flux)

                    if 'earth' in self.maddm_card['sigmav_method']:
                        for spec in self.Spectra.spectra_earth.keys():
                            if 'x' in spec: 
                                continue
                            
                            if 'neutrino' in spec:
                                out_spec = pjoin(out_dir, spec + '_flux_earth_'+spectrum_method+'.txt')
                                x = self.Spectra.spectra['x'] # !!! the x for oscillated neutrinos is the same for the spectra but not the same for positrons!!!
                                aux.write_data_to_file(x , self.Spectra.spectra_earth[spec] , filename = out_spec , header = h_flux )
    '''
    
    def print_ind(self,what, sig_th, sig_alldm, ul,  thermal=False ,direc = False , exp='Fermi'):
        alldm_mess = self.det_message_screen(sig_alldm , ul)

        if not direc:
        #if thermal:
        #   th_mess    = self.det_message_screen(sig_th    , ul)
        #   line = self.form_s( what +' (Thermal)= ') + '\t' + self.form_n(sig_th) + '\t' + self.form_s(th_mess) 
        #   line = line + '\t' + self.form_s('(All DM)= ') + self.form_n(sig_alldm) + '\t' + self.form_s(alldm_mess)
        #   line = line + '\t' + self.form_s('Fermi ul= ') + self.form_n(ul)

            if thermal:
                th_mess    = self.det_message_screen(sig_th    , ul)
                line = self.form_s(what+ '(Thermal)') + '=\t'+ self.form_n(sig_th)  + '   ' + self.form_s(th_mess)
                line = line +              '\t' + self.form_s('(All DM)= '  + self.form_n(sig_alldm) ) + '   ' + self.form_s(alldm_mess)
                line = line + '\t' + self.form_s(exp+' ul = '+ self.form_n(ul) )
    
    #           else:
    #              line = self.form_s( what+'(All DM)=') + '\t' + self.form_n(sig_alldm) + '\t' + self.form_s(alldm_mess) 
    #              line = line + '\t' + self.form_s(exp+'= ') + self.form_n(ul)
            else:                                                                                                                                                           
                line = self.form_s(what+ '(All DM)')+'=\t' + self.form_n(sig_alldm) + '   ' + self.form_s(alldm_mess)                                                           
                line = line + '\t' + self.form_s(exp+' ul = ' + self.form_n(ul) )    
    
            #{ 'n':'SigmaN_SD_n','sig':sigN_SD_n*units , 'lim':self.limits.SD_max(mdm,'n') , 'exp':'Lux2017' } ]

        elif direc:
            un = self.last_results['GeV2pb*pb2cm2']
            if thermal:
                th_mess    = self.det_message_screen(sig_th    , ul)
                line = self.form_s(what) + '\t' + self.form_s('(Thermal)= '+ self.form_n(sig_th)    +'GeV^-2')+' = '+self.form_s(self.form_n(sig_th*un)+'cm^2')
                line = line + '  ' + self.form_s(th_mess)

                line = line +              '\t' + self.form_s('(All DM)= '  + self.form_n(sig_alldm)+'GeV^-2')+' = '+ self.form_s(self.form_n(sig_alldm*un) +'cm^2')
                line = line + '  ' + self.form_s(alldm_mess) + '\t' + self.form_s(exp+' ul= ') +  self.form_n(ul)
            else:
                line = self.form_s(what) + '\t' + self.form_s('(All DM)= '  + self.form_n(sig_alldm)+'GeV^-2')+' = '+ self.form_s(self.form_n(sig_alldm*un) +'cm^2')
                line = line + '  ' + self.form_s(alldm_mess) + '\t'+  self.form_s(exp+' ul= ') +  self.form_n(ul)

        logger.info(line)
          




    def save_summary_single(self, relic = True , direct = False , indirect = False , fluxes_source = False , fluxes_earth = False):
        
        def form_s(stringa):
            formatted = '{:22}'.format(stringa)
            return  formatted
        
        def form_n(num):
            formatted = '{:3.2e}'.format(num)
            return formatted

        out = open(pjoin(self.dir_path, 'output','MadDM_results.txt'),'w')
        out.write('#############################################\n')
        out.write('#                MadDM v. ' + str(self.MadDM_version) +'                #\n' )
        out.write('#############################################\n\n\n')

        if relic:
            out.write('#############################################\n')
            out.write('# Relic Density                             #\n')
            out.write('#############################################\n\n')

            relic, planck , message = self.last_results['Omegah^2'] , self.limits._oh2_planck , self.det_message(self.last_results['Omegah^2'], self.limits._oh2_planck) 
    
            out.write(form_s('Omegah2')             + '= ' + form_n(relic)  + '\n')
            out.write(form_s('Omegah_Planck')       + '= ' + form_n(planck)  + '\n')
            if self.last_results['xsi'] > 0: 
                out.write(form_s('xsi') + '= ' + form_n(self.last_results['xsi']) +' \t # xsi = (Omega/Omega_Planck)\n' )
            out.write(form_s('x_f')                  + '= ' + form_n(self.last_results['x_f'])        + '\n' ) 
            out.write(form_s('sigmav_xf')           + '= ' + form_n(self.last_results['sigmav(xf)']) + '\n' ) 

        if direct:

            out.write('\n#############################################\n')
            out.write('# Direct Detection [cm^2]                   #\n')
            out.write('#############################################\n\n')


            for D in self.last_results['direct_results']:
                cross = D['sig']
                ul = D['lim']
                exp = D['exp']
                out.write(form_s(D['n']) + '= ' + form_s('['+ form_n(cross) + ',' + form_n(ul) + ']' ) + '# '+exp + '\n')


        if indirect:      
            sigmav_meth = self.maddm_card['sigmav_method']

            out.write('\n#############################################\n')
            out.write('# Indirect Detection [cm^3/s]               #\n')
            out.write('#############################################\n\n')

            out.write('# Annihilation cross section computed with the method: ' + sigmav_meth +' \n')

            tot_th , tot_ul , tot_sm = self.last_results['taacsID'] , self.last_results['Fermi_sigmav'] , self.last_results['tot_SM_xsec']
            #pvalue     = self.last_results['pvalue']
            #likelihood = self.last_results['likelihood'] 
            pvalue , likelihood = 1, 2 
            if sigmav_meth !='inclusive':
                out.write('# Fermi Limit for DM annihilation computed with Pythia8 spectra  \n\n')
                out.write(form_s('Total_xsec') +'= ' + form_s('['+ form_n(tot_th) +','+ form_n(tot_ul) +']') + '\n') 

            else:
                if 'PPPC' in  self.maddm_card['indirect_flux_source_method']: 
                    method =  self.maddm_card['indirect_flux_source_method'] 
                else: 
                    method = 'PPPC4DMID'
                out.write('# Fermi Limit computed with ' + method + ' spectra\n\n')
                lista = [ proc for proc in self.last_results.keys() if 'taacsID#' in proc and 'lim_' not in proc and 'err' not in proc ]
                for name in lista:
                    proc = name.replace('taacsID#','')
                    proc_th , proc_ul = self.last_results[name] , self.last_results['lim_'+name]
                    out.write(form_s(proc)      + '= '+ form_s('['+ form_n(proc_th)+',' + form_n(proc_ul)+']') + '\n')
           
                out.write(form_s('TotalSM_xsec')+ '= '+ form_s('['+ form_n(tot_sm)+ ',' + form_n(tot_ul) +']') + '\n')

            out.write( form_s('Fermi_Likelihood')+ '= '+ form_s(form_n(likelihood))  +'\n' )
            out.write( form_s('Fermi_pvalue'    )+ '= '+ form_s(form_n(pvalue    ))  +'\n ')

        if fluxes_source:
           out.write('\n#############################################\n')
           out.write('# CR Flux at source [particles/(cm^2 s sr)] #\n')
           out.write('#############################################\n\n')
           out.write('# Fluxes calculated using the spectra from ' + self.maddm_card['indirect_flux_source_method'] + '\n\n' )
           for name in ['neutrinos_e','neutrinos_mu','neutrinos_tau','gammas']:
                out.write(form_s('Flux_'+name)+ '= '+ form_s(form_n(self.last_results['flux_'+name] )) + '\n' )
 
        if fluxes_earth:
            out.write('############################################\n')
            out.write('# CR Flux at Earth [particles/(cm^2 s sr)] #\n')
            out.write('############################################\n')


        #'''
        #pass_message  = '%s ALLOWED  %s' % (bcolors.OKGREEN, bcolors.ENDC)
        #fail_message  = '%s EXCLUDED %s' % (bcolors.FAIL, bcolors.ENDC)
        #nolim_message = '%s NO LIMIT %s' % (bcolors.GRAY, bcolors.ENDC)
        #'''

    def det_message_screen(self,n1,n2):
        if n2 <= 0 :                 return '%s NO LIMIT %s' % (bcolors.GRAY, bcolors.ENDC)
        elif   n1 > n2 and n2 > 0 : return '%s EXCLUDED %s' % (bcolors.FAIL, bcolors.ENDC)
        elif   n2 > n1            : return '%s ALLOWED %s'  % (bcolors.OKGREEN, bcolors.ENDC) 
        elif   n1 <= 0            : return 'No Theory Prediction'

    def det_message(self,n1,n2):
        if n2 < 0 :                 return 'NO LIMIT'
        elif   n1 > n2 and n2 > 0 : return 'EXCLUDED'
        elif   n2 > n1            : return 'ALLOWED'  
        elif   n1 <= 0            : return 'No Theory Prediction'        
    
    def form_s(self,stringa):
        formatted = '{:25}'.format(stringa)
        return  formatted
    def form_n(self,num):
        formatted = '{:3.2e}'.format(num)
        return formatted

    
    def ask_run_configuration(self, mode=None, force=False):
        """ask the question about card edition / Run mode """
        
        process_data = self.proc_characteristics
        self.mode, cmd_quest = self.ask('', '0', mode=mode, 
                        data=self.proc_characteristics,
                        ask_class=MadDMSelector, timeout=60, path_msg=' ',
                        return_instance=True, force=force)
        
        # automatically switch to keep_wgt option
        #edit the maddm_card to be consistent with self.mode
        cmd_quest.get_cardcmd()

        self.maddm_card = cmd_quest.maddm
        for key, value in self.mode.items():
            if value == 'ON' or value is True:
                self.mode[key] = True

            elif value == 'OFF':
                self.mode[key] = False
        self.mode['capture'] = False
        # create the inc file for maddm
        logger.debug('2to2 in ask_run_configuration: %s' % self._two2twoLO)

        self.auto_width = set() #ensure to reset auto_width! at the 
        self.check_param_card(pjoin(self.dir_path, 'Cards', 'param_card.dat'))
        self.param_card = check_param_card.ParamCard(pjoin(self.dir_path, 'Cards', 'param_card.dat'))
        
        
        #set self._two2twoLO if we need to use inclusive method
        #If the Indirect subfolder is not created, that means that the code is
        #using the 2to2 at LO which is handled by maddm.f.
        if self.maddm_card['sigmav_method'] in ['inclusive']:
            self._two2twoLO = True
        elif not self.mode['indirect']:
            self._two2twoLO = False
        elif os.path.exists(pjoin(self.dir_path, 'Indirect')):
            self._two2twoLO = False
        else:
            misc.sprint(os.listdir(self.dir_path))
            raise Exception, 'Madevent output for indirect detection not available'
            
        logger.debug('2to2: %s' % self._two2twoLO)

        #set fortran switch and write include file
        if not self.in_scan_mode:
            # create the inc file for maddm
            self.maddm_card.set('do_relic_density', self.mode['relic'], user=False)
            self.maddm_card.set('do_direct_detection', True if self.mode['direct'] else False, user=False)
            self.maddm_card.set('do_directional_detection', self.mode['direct'] == 'directional', user=False)
            self.maddm_card.set('do_capture', self.mode['capture'], user=False)
            self.maddm_card.set('do_indirect_detection', True if self.mode['indirect'] else False, user=False)
            self.maddm_card.set('do_flux', True if (self.mode['indirect'] and self.mode['indirect'] != 'sigmav') else False, user=False)
            self.maddm_card.set('only2to2lo', self._two2twoLO, user=False)
            #self.maddm_card.set('run_multinest', self.mode['run_multinest'], user=False)
            self.maddm_card.write_include_file(pjoin(self.dir_path, 'include'))        
        
        if not self.in_scan_mode and not self.mode['nestscan']:
            logger.info("Start computing %s" % ','.join([name for name, value in self.mode.items() if value]))
        return self.mode


    def compile(self):
        """compile the code"""

        #logger.info(self.mode)

        if self.in_scan_mode:
            return

        self.maddm_card.write_include_file(pjoin(self.dir_path,'include'))
        misc.compile(['all'],cwd=self.dir_path)

        # if self.mode['relic'] and self.mode['direct']:
        #     misc.compile(['all'],cwd=self.dir_path)
        # elif self.mode['relic'] and not self.mode['direct']:
        #     misc.compile(['relic_density'],cwd=self.dir_path)
        # elif self.mode['direct'] and not self.mode['relic']:
        #     misc.compile(['direct_detection'],cwd=self.dir_path)
        # elif self.mode['indirect'] and not self.mode['relic']:
        #     misc.compile(['relic_density'],cwd=self.dir_path)
        # else:
        #     raise Exception, "No computation requested. End the computation"
    
        if os.path.exists(pjoin(self.dir_path, 'src', 'maddm.x')) or os.path.exists(pjoin(self.dir_path, 'maddm.x')):
            logger.info("compilation done")
        else:
            raise Exception, 'Compilation of maddm failed'
        
    ############################################################################
    def do_open(self, line):
        """Open a text file/ eps file / html file"""

        args = self.split_arg(line)
        # Check Argument validity and modify argument to be the real path
        self.check_open(args)
        file_path = args[0]

        misc.open_file(file_path)
    
    ############################################################################
    def check_open(self, args):
        """ check the validity of the line """

        if len(args) != 1:
            self.help_open()
            raise self.InvalidCmd('OPEN command requires exactly one argument')

        if args[0].startswith('./'):
            if not os.path.isfile(args[0]):
                raise self.InvalidCmd('%s: not such file' % args[0])
            return True

        # if special : create the path.
        if not self.dir_path:
            if not os.path.isfile(args[0]):
                self.help_open()
                raise self.InvalidCmd('No MadEvent path defined. Unable to associate this name to a file')
            else:
                return True

        path = self.dir_path
        if os.path.isfile(os.path.join(path,args[0])):
            args[0] = os.path.join(path,args[0])
        elif os.path.isfile(os.path.join(path,'Cards',args[0])):
            args[0] = os.path.join(path,'Cards',args[0])
        # special for card with _default define: copy the default and open it
        elif '_card.dat' in args[0]:
            name = args[0].replace('_card.dat','_card_default.dat')
            if os.path.isfile(os.path.join(path,'Cards', name)):
                files.cp(os.path.join(path,'Cards', name), os.path.join(path,'Cards', args[0]))
                args[0] = os.path.join(path,'Cards', args[0])
            else:
                raise self.InvalidCmd('No default path for this file')
        elif os.path.isfile(os.path.join(path, '%s_card.dat' % args[0])):
            args[0] = os.path.join(path, '%s_card.dat' % args[0])
        elif os.path.isfile(os.path.join(path, 'Cards' ,'%s_card.dat' % args[0])):
            args[0] = os.path.join(path, 'Cards', '%s_card.dat' % args[0])
        elif not os.path.isfile(args[0]):
            raise self.InvalidCmd('No default path for this file')

class Indirect_PY8Card(banner_mod.PY8Card):
    
    def default_setup(self):
        """ Sets up the list of available PY8 parameters."""
        
        self.add_param("Main:timesAllowErrors", 10, hidden=True, comment="allow a few failures before quitting")
        self.add_param("Main:NumberOfEvents", -1,  comment="number of events to go trough")
        self.add_param("Main:spareParm1", 1000, hidden=True, comment=" mass of the Dark matter")
        self.add_param("Main:spareWord1" , './', hidden=True, comment="specify output dir") 
        # Init
        self.add_param("Init:showChangedSettings", True, hidden=True, comment="list changed settingspython")
        self.add_param("Init:showChangedParticleData", True, hidden=True, comment="print changed particle and decay data")
        # Next
        self.add_param("Next:numberCount", 1000, hidden=True, comment="print message every n events")
        self.add_param("Next:numberShowInfo", 1, hidden=True, comment="print event information n times")
        self.add_param("Next:numberShowProcess", 1, hidden=True, comment="print process record n times")
        self.add_param("Next:numberShowEvent", 1, hidden=True, comment="print event record n times")
        #Beam
        self.add_param("Beams:frameType", 4, hidden=True,comment='Tell Pythia8 that an LHEF input is used.')
        self.add_param("Beams:LHEF", "unweighted_events.lhe.gz", hidden=True)
        self.add_param("PDF:lepton", False, hidden=True, comment="radiation from initial particles")
        # Parton level
        self.add_param("PartonLevel:MPI", False, hidden=True, comment="multiparton interactions")
        self.add_param("PartonLevel:ISR", False, hidden=True, comment="initial-state radiation")
        self.add_param("PartonLevel:FSR", True, hidden=True, comment="final-state radiation")
        # Weakshower <- allow the user to switch this ON
        self.add_param("TimeShower:weakShower", False, comment="Run weak-shower for FSR")
        self.add_param("TimeShower:weakShowerMode", 0, comment="Determine which branchings are allowed (0 -> W and Z)")
        self.add_param("TimeShower:pTminWeak", 0.1)
        self.add_param("WeakShower:vetoWeakJets", True)
        self.add_param("WeakShower:singleEmission", True)
        self.add_param("WeakShower:enhancement",1.0,comment="enhanced weak shower rate")
        # decay all particles
        self.add_param("13:mayDecay", True, hidden=True, comment="decays muons")
        self.add_param("211:mayDecay", True, hidden=True, comment="decays pions")
        self.add_param("321:mayDecay", True, hidden=True, comment="decays Kaons")
        self.add_param("130:mayDecay", True, hidden=True, comment="decays K_LO")
        self.add_param("310:mayDecay", True, hidden=True, comment="decays K_SO")
        self.add_param("2112:mayDecay", True, hidden=True, comment="decays neutron")


class MadDMSelector(cmd.ControlSwitch, common_run.AskforEditCard):
    """ """

    to_control= [('relic', 'Compute the Relic Density'),
                      ('direct', 'Compute direct(ional) detection'),
                      ('indirect', 'Compute indirect detection/flux'),
                      ('nestscan', 'Run Multinest scan'),
                ]
    to_init_card = ['param', 'maddm','pythia8']
    PY8Card_class = Indirect_PY8Card
    
    integer_bias = len(to_control) + 1 # integer corresponding to the first entry in self.cards
    
    ####################################################################
    # everything related to relic option
    ####################################################################    
    def set_default_relic(self):
        """set the default value for relic="""
        
        if self.availmode['has_relic_density']:
            self.switch['relic'] = 'ON'
        else:
            self.switch['relic'] = 'Not Avail.'

    def get_allowed_relic(self):
        """Specify which parameter are allowed for relic="""
        
        
        if hasattr(self, 'allowed_relic'):
            return getattr(self, 'allowed_relic')
        
        if self.availmode['has_relic_density']:
            self.allowed_relic = ['ON', 'OFF']
        else:
            self.allowed_relic = []
        
        return self.allowed_relic

    ####################################################################
    # everything related to direct option
    ####################################################################    
    def set_default_direct(self):
        """set the default value for relic="""
        
        if self.availmode['has_directional_detection']:
            self.switch['direct'] = 'directional'
        elif self.availmode['has_direct_detection']:
            self.switch['direct'] = 'direct'        
        else:
            self.switch['direct'] = 'Not Avail.'

    def get_allowed_direct(self):
        """Specify which parameter are allowed for relic="""
        
        
        if hasattr(self, 'allowed_direct'):
            return getattr(self, 'allowed_direct')

        if self.availmode['has_directional_detection']:
            self.allowed_direct =  ['directional', 'direct','OFF']
        elif self.availmode['has_direct_detection']:
            self.allowed_direct =  ['direct','OFF']
        else:
            return []

    def check_value_direct(self, value):
        """ allow diret=ON in top of standard mode """
        
        if value in self.get_allowed('direct'):
            return True
        
        value =value.lower()
        if value in ['on'] and self.availmode['has_direct_detection' ]:
            return 'direct'
        
    def get_cardcmd_direct(self, value):
        """return the command to set the maddm_card consistent with the switch"""
        
        cmd =[]
        if value == 'directional':
            cmd.append('set do_directional_detection True')
            value = 'direct'
        else:
            cmd.append('set do_directional_detection False')
        
        if value in ['ON', 'direct']:
            cmd.append('set do_direct_detection True')
        else:
            cmd.append('set do_direct_detection False')

        return cmd


    ####################################################################
    # everything related to indirect option
    ####################################################################    
    
    # TODO -> add check if PY8/Dragon are available for the switch 
    
    def set_default_indirect(self):
        """set the default value for relic="""
        
        if self.availmode['has_indirect_detection']:
            self.switch['indirect'] = 'sigmav'     
        else:
            self.switch['indirect'] = 'Not Avail.'

    def get_allowed_indirect(self):
        """Specify which parameter are allowed for relic="""
        
        
        if hasattr(self, 'allowed_indirect'):
            return getattr(self, 'allowed_indirect')

        if self.availmode['has_indirect_detection']:
            self.allowed_indirect =  ['OFF', 'sigmav', 'flux_source', 'flux_earth']
        else:
            return []

    
    def check_value_indirect(self, value):
        """ allow diret=ON in top of standard mode """
        
        other_valid = ['source_PPPC4DMID', 'source_py8', 
                       'earth_PPPC4DMID+dragon',
                     'earth_PPPC4DMID', 'earth_py8+dragon'] 
        
        if value in self.get_allowed('indirect'):
            return True
        
        value = value.lower()
        # valid but hidden options
        if value in other_valid:
            return True
        if value.startswith('flux_') and value[5:] in other_valid:
            return value[5:] 

        # handle other alias        
        if value in ['on'] and self.availmode['has_indirect_detection' ]:
            return 'sigmav'

        if value in ['dragon']:
            return 'earth_py8+dragon'
        
        if value in ['py8', 'py', 'pythia', 'pythia8']:
            return 'source_py8'
                
    ####################################################################
    # everything related to multinset option
    ####################################################################    
    def set_default_nestscan(self):
        """set the default value for nestscan="""
        
        if aux.module_exists('pymultinest'):
            self.switch['nestscan'] = 'OFF'
        else:
            self.switch['nestscan'] = 'Not Avail.'

    def get_allowed_nestscan(self):
        """Specify which parameter are allowed for relic="""
            
        if hasattr(self, 'allowed_nestscan'):
            return getattr(self, 'allowed_nestscan')
        
        if aux.module_exists('pymultinest'):
            self.allowed_nestscan = ['ON', 'OFF']
        else:
            self.allowed_nestscan = []
        
        return self.allowed_nestscan

    
    def __init__(self, question, *args, **opts):

        self.me_dir = opts['mother_interface'].dir_path
        self.availmode = opts.pop('data', collections.defaultdict(bool))
        # call the initialisation of the ControlSwith part

        cmd.ControlSwitch.__init__(self, self.to_control, opts['mother_interface'], *args, **opts)
        # initialise the various card to control
        question = ''
        
        param_card_path = pjoin(opts['mother_interface'].dir_path, 'Cards', 'param_card.dat')
        pythia8_card_path = pjoin(opts['mother_interface'].dir_path, 'Cards', 'pythia8_card.dat')

        cards = [param_card_path, 'maddm', 'multinest', pythia8_card_path]
        common_run.AskforEditCard.__init__(self, question, cards,
                                            *args, **opts)
        self.question = self.create_question()
                
    def init_maddm(self, path):
        """ initialize cards for the reading/writing of maddm"""
        
        self.maddm_def = MadDMCard(self.paths['maddm_default'])
        try:
            self.maddm = MadDMCard(self.paths['maddm'])
        except Exception as e:
            logger.error('Current maddm_card is not valid. We are going to use the default one.')
            logger.error('problem detected: %s' % e) 
            files.cp(self.paths['maddm_default'], self.paths['maddm'])
            self.maddm = MadDMCard(self.paths['maddm'])
        
        self.special_shortcut.update({'nevents': ([int],['Main:numberOfEvents %(0)s'])})
        self.special_shortcut_help.update({'nevents': 'number of events to generate for indirect detection if sigmav_method is madevent/reshuffling'})
        
        self.maddm_set = list(set(self.maddm_def.keys() + self.maddm_def.hidden_param))
        return self.maddm.keys() 

    def get_cardcmd(self):
        """ return the list of command that need to be run to have a consistent 
            set of cards with the switch value choosen """
        
        cmd = super(MadDMSelector,self).get_cardcmd()
        for c in cmd:
            self.exec_cmd(c)        

        return cmd 

    def create_question(self):
        """create the new question depending of the status"""
        
        # Technical note: 
        # Note that some special trigger happens for 
        #    7/8 trigger via self.trigger_7 and self.trigger_8
        #    If you change the numbering, please change the name of the 
        #    trigger function accordingly.
        
        question = cmd.ControlSwitch.create_question(self, help_text=False)
        question +="""\n%(start_green)s You can also edit the various input card%(stop)s:
 * Enter the name/number to open the editor
 * Enter a path to a file to replace the card
 * Enter %(start_bold)sset NAME value%(stop)s to change any parameter to the requested value
 /=============================================================================\ 
 |  5. Edit the model parameters    [%(start_underline)sparam%(stop)s]                                    |  
 |  6. Edit the MadDM options       [%(start_underline)smaddm%(stop)s]                                    |
"""

        current_val  = self.answer # use that to be secure with conflict -> always propose card
        if current_val['nestscan'] == "ON" or self.switch["nestscan"] ==  "ON":
            question += """ |  7. Edit the Multinest options  [%(start_underline)smultinest%(stop)s]                                 |\n"""
    
        if current_val['indirect'].startswith('flux') or self.switch["indirect"].startswith('flux'):
            question += """ |  8. Edit the Showering Card for flux  [%(start_underline)sflux%(stop)s]                                |\n"""
        
        question+=""" \=============================================================================/\n"""
        self.question =  question % {'start_green' : '\033[92m',
                         'stop':  '\033[0m',
                         'start_underline': '\033[4m',
                         'start_bold':'\033[1m', 
                        }
        return self.question
    
    def define_paths(self, **opt):
        
        super(MadDMSelector, self).define_paths(**opt)
        self.paths['maddm'] = pjoin(self.me_dir,'Cards','maddm_card.dat')
        self.paths['maddm_default'] = pjoin(self.me_dir,'Cards','maddm_card_default.dat')
        self.paths['multinest'] = pjoin(self.me_dir,'Cards','multinest_card.dat')
        self.paths['multinest_default'] = pjoin(self.me_dir,'Cards','multinest_card_default.dat')
        self.paths['flux'] = pjoin(self.me_dir,'Cards','pythia8_card.dat')
        self.paths['flux_default'] = pjoin(self.me_dir,'Cards','pythia8_card_default.dat')
    
    # TODO HERE!
    def default(self, line):
        """Default action if line is not recognized"""
        
        try:
            return cmd.ControlSwitch.default(self, line, raise_error=True)
        except cmd.NotValidInput:
            return common_run.AskforEditCard.default(self, line)     
        
    def trigger_7(self, line):
        """ trigger function for default function:
            allows to modify the line/ trigger action.
            
            check that nestscan can be ON
               if it is, set it ON,
               otherswise, print a warning and set the line to reask the 
               question (return False)
        """
        
        #  1) do like the user type "multinest=ON"
        #  2) go to the line edition
        self.set_switch('nestscan', "ON", user=True) 
        if self.switch['nestscan'] == "ON":
            return '7 %s' % line
        # not valid nestscan - > reask question
        else:
            return 

    def trigger_8(self, line):
        """ trigger function for default function:
            allows to modify the line/ trigger action.
            
            check that indirect can be ON
               if it is, set it ON,
               otherswise, print a warning and set the line to reask the 
               question (return False)
        """
        
        #  1) do like the user type "indirect=flux_source" (if not in flux)
        #  2) forbid sigmav_method = inclusive
        #  3) forbid PPPC4IDMD
        #  4) go to the line edition
        #1.
        if not self.switch['indirect'].startswith('flux_'):
            self.set_switch('indirect', "flux_source", user=True) 
            if not self.switch['indirect'] == 'flux_source':
                # not valid indirect - > reask question
                return 
            logger.warning('switching indirect mode to flux_source since sigmav does not support flux')
           
        #2. check that sigmv_method is not inclusive
        if self.maddm['sigmav_method'] == 'inclusive':
            self.setDM('sigmav_method', 'reshuffling',loglevel=30)
            
        #3. ensure pythia8 is on
        if self.maddm['indirect_flux_source_method'] != 'pythia8':
            self.setDM('indirect_flux_source_method', 'pythia8',loglevel=30)
        
        return '8 %s' % line
    
    trigger_flux = trigger_8
    
    def do_compute_widths(self, line):
        """normal fct but ensure that self.maddm_card is up-to-date"""
        
        try:
            self.mother_interface.maddm_card = self.maddm
        except Exception,error:
            logger.error("Invalid command: %s " % error)
            return
        return super(MadDMSelector, self).do_compute_widths(line)
            
    def do_help(self, line,conflict_raise=False, banner=True):
        """proxy for do_help"""
        
        if banner:                      
            logger.info('*** HELP MESSAGE ***', '$MG:color:BLACK')
        
        if line:
            card = common_run.AskforEditCard.do_help(self, line,conflict_raise=conflict_raise, banner=False)
        args = self.split_arg(line)
        
        start = 0
        if args[start] in ['maddm', 'maddm_card']:
            card = 'maddm'
            start += 1
            
        #### MADDM CARD 
        if args[start] in [l.lower() for l in self.maddm.keys()] and card in ['', 'maddm']:
            if args[start] not in self.maddm_set:
                args[start] = [l for l in self.maddm_set if l.lower() == args[start]][0]

            if args[start] in self.conflict and not conflict_raise:
                conflict_raise = True
                logger.info('**   AMBIGUOUS NAME: %s **', args[start], '$MG:color:BLACK')
                if card == '':
                    logger.info('**   If not explicitely speficy this parameter  will modif the maddm_card file', '$MG:color:BLACK')

            self.maddm.do_help(args[start])
        if banner:                      
            logger.info('*** END HELP ***', '$MG:color:BLACK')  
        return card    
    
    
    def do_update(self, line, timer=0):
        """syntax: update dependent: Change the mass/width of particles which are not free parameter for the model.
                    update missing:   add to the current param_card missing blocks/parameters.
           Bypass dependent mode but if request by the user
        """
         
        if timer == 0:
            return super(MadDMSelector, self).do_update(line)
        else: 
            args = self.split_arg(line)
            if args[0] == 'dependent':
                return
            else:
                return super(MadDMSelector, self).do_update(line)
         
    def update_to_full(self, line):
        """ trigger via update to_full LINE"""
        
        logger.info("update the maddm_card by including all the hidden parameter")
        self.maddm.use_full_template = True
        self.maddm.write(self.paths['maddm'], write_hidden=True)


    def check_card_consistency(self):
        
        super(MadDMSelector, self).check_card_consistency()

        #if there are any new jfactors, make sure to write them in
        #logger.debug('Updating the Jfactor file')
        #self.write_jfactors()

        # If direct detection is ON ensure that quark mass are not zero
        if self.run_options['direct'] != 'OFF':
            to_change = []
            for i in range(1,7):
                if self.param_card.get_value('mass', i) == 0.:
                    to_change.append(i)
            if to_change:
                logger.warning('For direct detection the quark mass need to be different of zero. Automatically adding such masses (PDG 2014).')
                quark_masses = {1: 4.8e-3, 2: 2.3e-3,  3: 95e-3, 4: 1.275, 5: 4.18, 6: 173.21}
                for i in to_change:
                    self.do_set('param_card mass %s %s' % (i, quark_masses[i]))
        
        
    def reload_card(self, path):
        """ensure that maddm object are kept in sync"""
        
        if path == self.paths['maddm']:
            try:
                self.maddm = MadDMCard(path) 
            except Exception as e:
                logger.error('Current maddm_card is not valid. We are going to use the default one.')
                logger.error('problem detected: %s' % e)
                logger.error('Please re-open the file and fix the problem.')
                logger.warning('using the \'set\' command without opening the file will discard all your manual change')
        else:
            return super(MadDMSelector,self).reload_card(path)
        
        
    def complete_set(self, text, line, begidx, endidx, formatting=True):
        """ Complete the set command"""
        possibilities = super(MadDMSelector,self).complete_set(text, line, begidx, endidx, formatting=False)
        args = self.split_arg(line[0:begidx])
        if len(args)>1 and args[1] == 'maddm_card':
            start = 2
        else:
            start = 1 
            if len(args) ==1:
                possibilities['category of parameter (optional)'] += \
                                self.list_completion(text, ['maddm_card'], line)
                
        if len(args)==start:
            correct = self.maddm_set + ['default']
            possibilities['maddm Card'] = self.list_completion(text, correct, line)   
        elif len(args)==start+1:
            allowed_for_run = []
            if args[-1].lower() in self.maddm.allowed_value:
                allowed_for_run = self.maddm.allowed_value[args[-1].lower()]
                if '*' in allowed_for_run: 
                    allowed_for_run.remove('*')
            elif isinstance(self.maddm[args[-1]], bool):
                allowed_for_run = ['True', 'False']
            opts = [str(i) for i in  allowed_for_run]
            possibilities['maddm Card'] = self.list_completion(text, opts)
        
        return self.deal_multiple_categories(possibilities, formatting)


    def do_set(self, line):
        """ edit the value of one parameter in the card"""
        
        args = self.split_arg(line)
        # fix some formatting problem
        if '=' in args[-1]:
            arg1, arg2 = args.pop(-1).split('=')
            args += [arg1, arg2]
        if '=' in args:
            args.remove('=')
        args = [ a.lower() for a in args]
        
        if args[0] in ['maddm_card']:
            start = 1
            if args[1] == 'default':
                logging.info('replace %s by the default card' % args[0])
                self.maddm = MadDMCard(self.paths['maddm_default'])
        elif args[0] in self.maddm_set and args[0] not in self.conflict:
            start = 0

        #For setting the exp. constraints
        elif args[0] == 'dd_si_limits':
            if len(args) >2:
                self.limits._sigma_SI = float(args[1])
                self.limits._sigma_SI_width = float(args[2])
            else:
                logger.info('The file you use for direct detection limits will be interpreted as:')
                logger.info('Column 1 - dark matter mass in GeV')
                logger.info('Column 2 - upper limit on direct detection cross section in pb')
                self.limits._dd_si_limit_file = args[1]
                self.limits.load_constraints()

        elif args[0] == 'dd_sd_limits':

            if len(args) > 3:
                if args[1] == 'p':
                    self.limits._sigma_SDp = args[2]
                    self.limits._sigma_SDp_width = args[3]
                elif args[1] == 'n':
                    self.limits._sigma_SDn = args[2]
                    self.limits._sigma_SDn_width = args[3]

                else:
                    logger.error('set dd_sd_limits needs to be formatted as:')
                    logger.error('set dd_sd_limits <p or n> <observed val> <uncertainty>, or')
                    logger.error('set dd_sd_limits <p or n> <limit filename>, or')
            else:
                logger.info('The file you use for direct detection limits will be interpreted as:')
                logger.info('Column 1 - dark matter mass in GeV')
                logger.info('Column 2 - upper limit on direct detection cross section in pb')
                if args[1] == 'p':
                    self.limits._dd_sd_proton_limit_file = args[2]
                elif args[2] == 'n':
                    self.limits._dd_sd_neutron_limit_file = args[2]
                else:
                    logger.error('set dd_sd_limits needs to be formatted as:')
                    logger.error('set dd_sd_limits <p or n> <observed val> <uncertainty>, or')
                    logger.error('set dd_sd_limits <p or n> <limit filename>, or')

                self.ExpConstraints.load_constraints()

        elif args[0] == 'relic_limits':
            logger.info('The info is interpreted as: <Oh^2>, CL width')
            self.limits._oh2_planck = float(args[1])
            self.limits._oh2_planck_width = float(args[2])

        elif args[0] == 'id_limits':
             if len(args)!= 4:
                 logger.warning('You need to provide the following <ave. velocity> <ann. channel> <file path>')
                 logger.warning('or  <ave. velocity> <ann. channel> <obs. cross section> <obs. uncertainty>')
                 logger.warning('Annihilation channel can be: '+ str(self._allowed_final_states))
             logger.info('The file you use for indirect detection limits will be interpreted as:')
             logger.info('Column 1 - dark matter mass in GeV')
             logger.info('Column 2 - upper limit on the total annihilation cross section in cm^3/s at specified average velocity')

             if len(args) > 4:
                 if args[2] in self.limits._allowed_final_states:
                    self.limits._id_limit_vel[args[2]] = float(args[1])
                    self.limits._sigma_ID[args[2]] = float(args[3])
                    self.limits._sigma_ID_width[args[2]] = float(args[4])
                 else:
                     logger.error('Final state not allowed for ID limits!')
             else:
                 vel = float(args[1])
                 channel = args[2]
                 id_file = args[3]
                 self.limits._id_limit_vel[channel] = vel
                 self.limits._id_limit_file[channel] = id_file

                 self.limits.load_constraints()

        else:
            return super(MadDMSelector, self).do_set(line)
        if args[start+1] == 'default':
            default = self.maddm_def[args[start]]
            self.setDM(args[start], default)
        else:
            if args[start] in self.maddm.list_parameter or \
                   args[start] in self.maddm.dict_parameter:
                val = ' '.join(args[start+1:])
                val = val.split('#')[0]
                self.setDM(args[start], val)
            else:
                self.setDM(args[start], args[start+1:][0])
        # write the new file
        self.maddm.write(self.paths['maddm'])
        
    def setDM(self, name, value, loglevel=20):
        logger.log(loglevel,'modify parameter %s of the maddm_card.dat to %s' % (name, value))
        self.maddm.set(name, value, user=True)


class InvalidMaddmCard(banner_mod.InvalidRunCard):
    pass
             
class MadDMCard(banner_mod.RunCard):
    """MadDM card use the read/write of the runcard object"""
    
    filename = 'maddm_card'
    default_include_file = 'maddm_card.inc'
    initial_jfactors = {}
    initial_distances = {}
    full_template = False
    
    def __new__(cls, finput=None):
        """Bypass the standard RunCard one"""
        return super(banner_mod.RunCard, cls).__new__(cls, finput)

    def fill_jfactors(self, filename=pjoin(MDMDIR,'Jfactors','jfactors.dat')):

        if MadDMCard.initial_jfactors:
            self['jfactors'] = dict(MadDMCard.initial_jfactors)
            self['distances'] = dict(MadDMCard.initial_distances)
            return
        try:
            infile = open(filename,'r')
            lines = infile.readlines()
            infile.close()
            temp = dict()
            temp2=dict()
            for line in lines:
                if line.startswith('#') or line.strip() == '':
                    continue
                spline = line.split()
                jfact = spline[0]
                dist = float(spline[1])
                val = float(spline[2].rstrip())
                temp[jfact] = val
                temp2[jfact]= dist
            self['jfactors'] = temp
            self['distances'] = temp2
            MadDMCard.initial_jfactors = dict(self['jfactors'])
            MadDMCard.initial_distances = dict(self['distances'])
            
        except OSError:
            logger.error('could not open Jfactor file %s ' % filename)
            return False

        #logger.info('Loaded the following Jfactors:')
        #logger.info('            Object      Distance [kpc]   Jfactor [GeV^2/cm^5]')
        #for jf, val in MadDMCard.initial_jfactors.iteritems():
        #    dist = self['distances'][jf]
        #    logger.info('%20s     %.2e           %.2e' %(jf, dist, val))

    def write_jfactors(self):

        if self['jfactors'] == MadDMCard.initial_jfactors:
            return
        if not madgraph.ReadWrite:
            return
        MadDMCard.initial_jfactors = dict(self['jfactors'])
        MadDMCard.initial_distances = dict(self['distances'])
        fsock = open(pjoin(MDMDIR,'Jfactors','jfactors.dat'),'w')
        for data, value in self['jfactors'].items():
            if data == '__type__':
                continue
            dist = self['distances'][data]
            fsock.write('%s %s %s\n' % (data, dist, value))
        fsock.close()
            

    def default_setup(self):
        """define the default value"""

        self.add_param('print_out', False, hidden=True, comment="status update" )
        self.add_param('print_sigmas', False, comment="print out the values of the annihilation cross section at x_f (s = 4m^2/(1-2/x_f))")
        
        self.add_param('relic_canonical', True)
        self.add_param('do_relic_density', True, system=True)
        self.add_param('do_direct_detection', False, system=True)
        self.add_param('do_directional_detection', False, system=True)
        self.add_param('do_capture', False, system=True)
        self.add_param('do_flux', False, system=True, include=False)
        

        self.add_param('do_indirect_detection', False, system=True)
        self.add_param('only2to2lo', False, system=True)
        self.add_param('run_multinest', False, system=True, include=False)

        
        self.add_param('calc_taacs_ann_array', True)
        self.add_param('calc_taacs_dm2dm_array', True)
        self.add_param('calc_taacs_scattering_array', True)
        
        self.add_param('eps_ode', 0.01)
        self.add_param('xd_approx', False)
        
        self.add_param('x_start', 50.0, hidden=True)
        self.add_param('x_end', 1000.0, hidden=True)
        self.add_param('dx_step', 1.0, hidden=True)
        
        self.add_param('ngrid_init', 50, hidden=True, comment="Initial number of points in the grid to integrate over (of dm velocity)")
        self.add_param('nres_points', 50, hidden=True, comment="Number of points to add for one width around each resonance peak. (the code will add points which exponentially increase in distance from the pole)")
        self.add_param('eps_wij', 0.01, hidden=True, comment="precision of the romberg integration for Wij (irrelevant if simpson's rule used)")
        self.add_param('iter_wij', 2, hidden=True, comment="minimum number of iterations in the romberg integration algorithm for both Wij (irrelevant if simpson's rule used)")
        
        # Direct Detection nucleon form factors (fPx - proton, fNx - neutron)
        #Scalar FF
        self.add_param('SPu', 0.0153, hidden=True)
        self.add_param('SPd',0.0191, hidden=True)
        self.add_param('SPs', 0.0447, hidden=True)
        self.add_param('SPg', 1 - 0.0153 - 0.0191 - 0.0447, hidden=True)
        self.add_param('SNu', 0.0110, hidden=True)
        self.add_param('SNd', 0.0273, hidden=True)
        self.add_param('SNs', 0.0447, hidden=True)
        self.add_param('SNg', 1 - 0.0110 - 0.0273 - 0.0447, hidden=True)
        # Vector FF
        self.add_param('VPu', 2.0, hidden=True)
        self.add_param('VPd', 1.0, hidden=True)
        self.add_param('VNu', 1.0, hidden=True)
        self.add_param('Vnd', 2.0, hidden=True)
        # Axial Vector FF
        self.add_param('AVPu',0.842, hidden=True)
        self.add_param('AVPd',-0.427, hidden=True)
        self.add_param('AVPs',-0.085, hidden=True)
        self.add_param('AVNu',-0.427, hidden=True)
        self.add_param('AVNd',0.842, hidden=True)
        self.add_param('AVNs',-0.085, hidden=True)
        # Sigma(mu,nu) FF
        self.add_param('SigPu',0.84, hidden=True)
        self.add_param('SigPd',-0.23, hidden=True)
        self.add_param('SigPs',-0.046, hidden=True)
        self.add_param('SigNu',-0.23, hidden=True)
        self.add_param('SigNd',0.84, hidden=True)
        self.add_param('SigNs',-0.046, hidden=True)
        #
        # For Directional detection and direct detection rates
        #
        self.add_param('material', 1, allowed = range(1,14),
                        comment="""Choose the target material
    - 1: Xenon
    - 2: Germanium 
    - 3: Silicon
    - 4: Argon
    - 5: Neon
    - 6: Sodium
    - 7: Iodine
    - 8: Carbon
    - 9: Flourine
    - 10: Sulphur
    For Compounds:
    - 11: NaI 
    - 12: CF4
    - 13: CS2""")
        
        #Setting up the DM constants
        self.add_param('vMP', 220.0)
        self.add_param('vescape', 650.0)
        self.add_param('rhoDM', 0.3)
        #detector 
        self.add_param('detector_size', 1000.0)
        self.add_param('En_threshold', 4.0)
        self.add_param('lambd', 1.0)
        self.add_param('sig_theta', 3.0)
        
        self.add_param('En_min', 0.0)
        self.add_param('En_max', 100.0)
        self.add_param('cos_min', -1.0)
        self.add_param('cos_max', 1.0)
        self.add_param('day_min', 0.0)
        self.add_param('day_max', 365.0)
        self.add_param('Energy_bins', 100)
        self.add_param('cos_theta_bins', 20)
        self.add_param('day_bins', 10)
        self.add_param('smearing', False)
        
        #For the solar/earth capture rate

        #indirect detection
        self.add_param('vave_indirect', 0.00002, include=True)
        self.add_param('halo_profile', 'Draco', include=True, legacy=True)   ##### this naming doesn't make sense fix it!!!!!!!!

        self.add_param('jfactors', {'__type__':1.0}, include=False, hidden=True)
        self.add_param('distances', {'__type__':1.0}, include=False, hidden=True)
        self.add_param('npts_for_flux', 200, include=False, hidden=True) #number of points for the flux diff. distribution

        self.add_param('only_two_body_decays', True, include=False)
        #self.add_param('num_of_elements', 60)
        
        self.fill_jfactors()

        self.add_param('indirect_flux_source_method', 'pythia8', comment='choose between pythia8,PPPC4DMID and PPPC4DMID_ew', include=False,
                       allowed=['pythia8','PPPC4DMID','PPPC4DMID_ew'])
        self.add_param('indirect_flux_earth_method', 'dragon', comment='choose between dragon and PPPC4DMID_ep', include=False,
                       allowed=['dragon', 'PPPC4DMID_ep'])
        self.add_param('sigmav_method', 'reshuffling', comment='choose between inclusive, madevent, reshuffling', include=False,
                       allowed=['inclusive', 'madevent', 'reshuffling'])


        self.add_param('save_output', 'off', comment='choose between off, spectra or all to save scan output', include=False,
                       allowed=['spectra', 'all', 'off'])

        self.add_param('dm_profile', 'Ein', comment='choose the halo density profile: Ein, Moo, NFW, Iso', include = False, \
                           hidden = True, allowed = ['Ein','Iso','NFW','Moo'])
        self.add_param('prop_method', 'MED', comment='choose the propagation method: MIN, MED, MAX'        , include = False, \
                           hidden = True, allowed = ['MIN','MED','MAX'])
        self.add_param('halo_funct' , 'MF1', comment='choose the halo function: MF1, MF2, MF3'             , include = False, \
                           hidden = True , allowed = ['MF1','MF2','MF3'])





    def write(self, output_file, template=None, python_template=False,
              write_hidden=False):
        """Write the run_card in output_file according to template 
           (a path to a valid run_card)"""

        #self.write_jfactors() 
        

        if not template:
            if self.full_template:
                template = pjoin(MDMDIR, 'Templates', 'Cards', 'maddm_card_full.dat')
            elif 'spd' in self.user_set and 'print_out' in self.user_set:
                self.full_template = True
                template = pjoin(MDMDIR, 'Templates', 'Cards', 'maddm_card_full.dat')
            else:
                template = pjoin(MDMDIR, 'Templates', 'Cards', 'maddm_card.dat')
            python_template = True

        super(MadDMCard, self).write(output_file, template=template,
                                    python_template=python_template,
                                    write_hidden=write_hidden) 

    def check_validity(self):
        """ """
        
        super(MadDMCard, self).check_validity()
        
        if self['SPu'] + self['SPs'] + self['SPd'] + self['SPg'] -1 > 1e-3:
            raise InvalidMaddmCard , 'The sum of SP* parameter should be 1.0 get %s' % (self['SPu'] + self['SPs'] + self['SPd'] + self['SPg'])
        
        if self['SNu'] + self['SNs'] + self['SNd'] - self['SNg'] -1 > 1e-3:
            raise InvalidMaddmCard, 'The sum of SM* parameter should be 1.0 get %s' % (self['SNu'] + self['SNs'] + self['SNd'] + self['SNg'])
        
        if self['sigmav_method'] == 'inclusive':
            if self['indirect_flux_source_method'] == 'pythia8':
                if self['do_indirect_detection']:
                    logger.warning('since sigmav_method is on inclusive, indirect_flux_source_method has been switch to PPPC4MID')
                #if self['do_flux']:
                    #logger.warning('since sigmav_method is on inclusive, indirect_flux_source_method has been switch to PPPC4MID')
                self['indirect_flux_source_method'] == 'PPPC4DMID'
            if self['indirect_flux_earth_method'] != 'PPPC4DMID':
                if self['do_flux']:
                    logger.warning('since sigmav_method is on inclusive, indirect_flux_earth_method has been switch to PPPC4MID')
                self['indirect_flux_earth_method'] == 'PPPC4DMID' 
        ## FF must be fixed       

        elif self['indirect_flux_earth_method'] == 'PPPC4DMID+dragon':
            if self['indirect_flux_source_method'] != 'PPPC4DMID':
                logger.warning('since indirect_flux_earth_method is on PPPC4DMID+dragon, indirect_flux_source_method has been switch to PPPC4DMID')
                self['indirect_flux_source_method'] == 'PPPC4DMID'
        elif self['indirect_flux_earth_method'] == 'PPPC4DMID':
            if self['indirect_flux_source_method'].lower() not in ['PPPC4DMID', 'none']:
                logger.warning('since indirect_flux_earth_method is on PPPC4DMID, indirect_flux_source_method has been switch to none')
            self['indirect_flux_source_method'] = 'none'



                
                
class Indirect_Cmd(me5_interface.MadEventCmdShell):
    
    def __init__(self, *args, **opts):
        
        super(Indirect_Cmd, self).__init__(*args, **opts)
        self.history_header = ''
        self.options['automatic_html_opening'] = False
    
    def preloop(self,*args,**opts):
        super(Indirect_Cmd,self).preloop(*args,**opts)
        self.prompt = 'Indirect:%s' % self.prompt
    
    def do_plot(self, line):
        return
    
    def do_madanalysis5_parton(self, line):
        return
    def do_madanalysis5_hadron(self, line):
        return

    def make_make_all_html_results(self, folder_names = []):
        """keep track of the P results via an instance variable.
           Do not do the html output
        """
        run = self.results.current['run_name']
        Presults = sum_html.collect_result(self, folder_names=folder_names)

        self.Presults = {}
        for i, P_comb in enumerate(Presults):
            P_comb.compute_values()
            self.Presults['xsec_%s' % P_comb.name[7:]] = P_comb.xsec
            self.Presults['xerr_%s' % P_comb.name[7:]] = P_comb.xerru

        return Presults.xsec, Presults.xerru

        


class Priors:
    priors = ['uniform', 'loguniform', 'user']

class Likelihoods:
    likelihoods = ['gaussian', 'half_gauss', 'user', 'off']
    observables = ['relic', 'directSI','directSD_p', 'directSD_n', 'indirect', 'capture']


class Multinest(object):

    def __init__(self, run_interface):

        self.options = {
            'prior':'loguniform',
            'loglikelihood':{'relic':'gaussian', 'directSI':'half_gauss', 'directSD_p':'half_gauss','directSD_n':'half_gauss', 'indirect':'half_gauss'},
            'livepts':50000,
            'sampling_efficiency':'model',
            'parameters':[],
            'prefix':'mnest_',
            'half_gauss_width':{'spinSI':0.01, 'spinSD_p':0.01, 'spinSD_n':0.01} #log of the width
        }

        self.maddm_run = run_interface

        self.param_blocks, _ = self.maddm_run.param_card.analyze_param_card()
        #self.parameter_vars = [] #names of parameters to scan over
        self.output_observables = []

        self.counter = 0

    #Starts the multinest run.
    def launch(self, resume=True):

        self.param_card_orig = param_card_mod.ParamCard(self.maddm_run.param_card)

        for pdg in self.maddm_run.auto_width:
            self.param_card_orig.get('decay', pdg).value = 'Auto'
        
        if self.options['loglikelihood'] == {} or self.options['prior'] =='':
            logger.error("You have to set the priors and likelihoods before launching!")
            return False

        if len(self.options['parameters'])==0:
            logger.error("Multinest needs you need to set up parameters to scan over before launching! [multinest_card]")
            return False

        #Here print out some stuff about which parameters you're scanning over etc.
        parameters = self.options['parameters']
        logger.info("Scanning over the following [parameter, min, max]:")
        for i in range(len(parameters)):
            logger.info(str(parameters[i]))


        #if output_observables not set, automatically set it according to observables which are calculated
        #this is needed because maddm.out file contains too much information
        if self.output_observables == [] or 'default' in self.output_observables:
            if 'default' in self.output_observables:
                self.output_observables.remove('default')
            
            if self.maddm_run.mode['relic']:
                self.output_observables.append('omegah2')
                self.output_observables.append('sigmav_xf')
                self.output_observables.append('x_freezeout')

            if self.maddm_run.mode['direct']:
                self.output_observables.append('sigmaN_SI_neutron')
                self.output_observables.append('sigmaN_SI_proton')
                self.output_observables.append('sigmaN_SD_neutron')
                self.output_observables.append('sigmaN_SD_proton')
                self.output_observables.append('Nevents')

            if self.maddm_run.mode['indirect']:
                self.output_observables.append('taacsID')
                detailled_keys = [k for k in self.maddm_run.last_results.keys() if k.startswith('taacsID#')]
                for key in detailled_keys:
                    self.output_observables.append(key)

            if self.maddm_run.mode['capture']:
                detailled_keys = [k for k in self.maddm_run.last_results.keys() if k.startswith('ccap_')]
                for key in detailled_keys:
                    self.output_observables.append(key)

            #FIX THIS! HERE ADD FLUXES        
            #if self.maddm_run.mode['CR_flux']:
            #    detailled_keys = [k for k in self.maddm_run.last_results.keys() if k.startswith('??????????')]
            #    for key in detailled_keys:
            #        self.output_observables.append(key)
               

        # number of parameters to output
        # this includes the parameters which are scanned over (they go in first)
        # and the output parameters like relic density, dd cross section etc ...
        n_parameters = len(parameters)
        n_dimensions = n_parameters
        #mnest.parameter_vars = [parameters[i][0] for i in range(n_parameters)]

        n_parameters=n_parameters+len(self.output_observables)


        logger.info("Multinest will run with the following parameters: " )
        for key, item in self.options.iteritems():
            logger.info("%20s :  %s" %(key, item))

        pymultinest.run(self.myloglike, self.myprior,
                        n_dims=n_dimensions,
                        n_params=n_parameters,
                        importance_nested_sampling = False,
                        resume = resume,
                        verbose = False,
                        sampling_efficiency = self.options['sampling_efficiency'],
                        n_live_points = self.options['livepts'],
                        outputfiles_basename = pjoin(self.maddm_run.dir_path,'multinest_chains', self.options['prefix']))

        logger.info('Output written in %s' % pjoin(self.maddm_run.dir_path,'multinest_chains'))
        logger.info('Output  of .txt file formatted as:')
        logger.info('column[0] : Sample Probability (Weight)')
        logger.info('column[1] : -2 log(Likelihood)')
        for i, var in enumerate(self.options['parameters']):
            logger.info('column[%i] : %s' % (i+2, var[0]))
        for i, obs in enumerate(self.output_observables):
            logger.info('column[%i] : %s' % (i+2+len(self.options['parameters']), obs))

        return True


    def write_log(self, file =''):
        if file =='':
            file =  pjoin(self.maddm_run.dir_path,'multinest_chains', self.options['prefix']+'info.log')
            with open(file, 'w+') as f:
                f.write('# Please refer to the Multinest README file for more info about output.\n')
                f.write('#  options \n')
                for option in self.options:
                    f.write('%s  :  %s \n' % (option, self.options[option]))
                f.write('#  output format\n')
                f.write('column[0] : Sample Probability (Weight)\n')
                f.write('column[1] : -2 log(Likelihood)\n')
                for i,var in enumerate(self.options['parameters']):
                    f.write('column[%i] : %s\n' % (i+2, var[0]))
                for i, obs in enumerate(self.output_observables):
                    f.write('column[%i] : %s\n' % (i+2+len(self.options['parameters']), obs))



    def load_parameters(self, multinest_card):
        with open(multinest_card) as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('#'):
                    continue
                elif 'half_gauss' in line:
                    spline = line.split('#')[0].split()
                    type = spline[1]
                    half_gauss = spline[2]
                    self.options['loglikelihood'][type] = half_gauss
                    if len(spline)==4:
                        width_val = float(spline[3])
                        self.options['half_gauss_width'][type] = width_val
                else:
                    spline = line.split('#')[0].split()
                    if len(spline) ==0:
                        continue
                    #logger.debug(spline)
                    opt1 = spline[0]
                    if opt1 == 'scan_parameter':
                        var_name = spline[1]
                        min = spline[2]
                        max = spline[3]
                        self.options['parameters'].append([var_name, float(min), float(max)])
                    elif opt1 =='output_variables':
                        for j in range(1, len(spline)):
                            if spline[j] not in self.output_observables:
                                self.output_observables.append(spline[j])
                    elif len(spline)==3:
                        opt2 = spline[1]
                        if aux.isfloat(spline[2]):
                            value = int(spline[2])
                        else:
                            value = spline[2]

                        self.options.get(opt1)[opt2] = value
                    elif len(spline) ==2:
                        if aux.isfloat(spline[1]):
                            value = int(spline[1])
                        else:
                            value = spline[1]

                        self.options[opt1] = value

        logger.debug(self.options['parameters'])


    def change_parameter(self, name, val):

        if not self.param_blocks:
            logger.error('Parameter blocks have not been read in properly! Can not change the param. value!')
        else:
            try:
                block, lhaid = self.param_blocks[name][0]
                lhaid2 = lhaid[0]
                self.maddm_run.param_card[block].get(lhaid2).value = val
            except:
                logger.error('change_parameter() can not find parameter %s in the param_card.dat file' % name)
                sys.exit()

    def myprior(self, cube, ndim, nparams):


        if self.options['prior']=='uniform':
            for i in range(ndim):
                param_min = self.options['parameters'][i][1]
                param_max = self.options['parameters'][i][2]

                cube[i] = param_min + cube[i] * (param_max - param_min)
        elif self.options['prior']=='loguniform':
            for i in range(ndim):
                param_min = self.options['parameters'][i][1]
                param_max = self.options['parameters'][i][2]
                cube[i] = np.power(10., np.log10(param_min) + cube[i]*(np.log10(param_max) - np.log10(param_min)))
        elif self.options['prior']=='user': #USER DEFINED
            for i in range(ndim):
                param_min = self.options['parameters'][i][1]
                param_max = self.options['parameters'][i][2]
                cube[i] = cube[i]
        else:
            logger.error('Not a valid prior choice! Only allowed priors are %s '% str(Priors.priors.keys()))


    def myloglike(self,cube, ndim, nparams):

        chi = 1.0

        #Change the parameters and write them into an appropriate param_card.dat file.
        # reset param_card to original one (important for auto width)
        self.maddm_run.param_card = param_card_mod.ParamCard(self.param_card_orig)
        for i in range(len(self.options['parameters'])):
            logger.debug('Changing parameter %s to %.3e' %( self.options['parameters'][i][0], cube[i]))
            self.change_parameter(self.options['parameters'][i][0].lower(), cube[i])

        self.maddm_run.param_card.write(str(pjoin(self.maddm_run.dir_path, 'Cards', 'param_card.dat')))

        #Execute MadDM with the new parameters
        self.maddm_run.do_launch('-f')

        results = self.maddm_run.last_results
        results = collections.OrderedDict(sorted(results.items()))

        try:
#            if self.output_observables != []:
#            logger.debug('output obs: %s ' % self.output_observables)
            for i, observable in enumerate(self.output_observables):
                #if it's dm-nucleon cross section convert the units to cm2
                if observable.startswith('sigmaN'):
                    cube[ndim+i] = results[observable] * GeV2pb *pb2cm2
                elif observable in results:
                    cube[ndim+i] = results[observable]
                elif observable in self.param_blocks:
                    block, lhaid = self.param_blocks[observable][0]
                    lhaid2 = lhaid[0]
                    cube[ndim+i] = self.maddm_run.param_card.get_value(block, lhaid2)
                    logger.warning('LHAID: %.3e  VALUE: %.3e' %(lhaid2, self.maddm_run.param_card.get_value(block, lhaid2)))
                    
#            else:
#                for i, observable in enumerate(results):
#                    cube[ndim+i] = results[observable]
        except:
            logger.error('Observable %s does not exist' % observable)
            return

        if self.maddm_run.mode['relic']:
            omegah2 = results['omegah2']
        if self.maddm_run.mode['direct']:
            spinSI = 0.5*(results['sigmaN_SI_proton'] + results['sigmaN_SI_neutron']) * GeV2pb * pb2cm2
            spinSDp = results['sigmaN_SD_proton']  * GeV2pb * pb2cm2
            spinSDn = results['sigmaN_SD_neutron'] * GeV2pb * pb2cm2
        #<=========== for ID we will need each channel separately.
        if self.maddm_run.mode['indirect']:
            sigmavID = {'tot': results['taacsID']}
            id_vel = self.maddm_run.maddm_card['vave_indirect']

            detailled_keys = [k.split("#")[1] for k in results.keys() if k.startswith('taacsID#')]
            for key in detailled_keys:
                logger.debug('taacsID key: %s', key)
                sigmavID[key] = results['taacsID#%s' %(key)]

#            logger.debug('sigmavID : %s' % sigmavID)


        mdm = self.maddm_run.param_card.get_value('mass', self.maddm_run.proc_characteristics['dm_candidate'][0])
        #logger.debug('MDM: %.3e', mdm)
        logger.debug(self.maddm_run.mode['relic'])
        logger.debug(self.maddm_run.mode['direct'])

        for obs, likelihood in self.options.get('loglikelihood').iteritems():

            #relic density
            if obs == 'relic' and self.maddm_run.mode['relic']:

                #if relic density evaluates to -1 (too early freezeout for example)
                #then discard this point by assigning some high log likelihood.
                if omegah2 < 0:
                    chi+= __mnestlog0__
                else:
                    if likelihood == 'gaussian':
                        chi += -0.5*pow(omegah2 - self.maddm_run.limits._oh2_planck,2)/pow(self.maddm_run.limits._oh2_planck_width,2)
                    elif likelihood == 'half_gauss':
                        if omegah2 > 0:
                            #chi += np.log(0.5*(np.tanh((self.maddm_run.limits._oh2_planck - omegah2)\
                            #                           /self.maddm_run.limits._oh2_planck)+1.000001))
                            if omegah2 > self.maddm_run.limits._oh2_planck:
                                chi+= -0.5*pow(np.log10(self.maddm_run.limits._oh2_planck/omegah2),2)\
                                      /pow(np.log10(self.maddm_run.limits._oh2_planck_width),2)
                    elif likelihood =='user':
                        chi+=0
                        #
                        # HERE ADD YOUR OWN LOG LIKELIHOOD FOR RELIC DENSITY
                        #
                    elif likelihood == '':
                        chi+=0
                    else:
                        logger.error('You are not using a valid likelihood function for relic density. Omitting the contribution!')

            #direct detection (SI)
            if obs == 'directSI' and self.maddm_run.mode['direct']:

                    if likelihood=='half_gauss':

                        if spinSI > self.maddm_run.limits.SI_max(mdm):
                            chi+= -0.5*pow(np.log10(self.maddm_run.limits.SI_max(mdm)/spinSI),2)\
                                      /pow(self.options['half_gauss_width']['spinSI'],2)

                    elif likelihood == 'gaussian':
                        if self.maddm_run.limits._sigma_SI > 0:
                            chi+=  -0.5*pow(spinSI - self.maddm_run.limits._sigma_SI,2)/pow(self.maddm_run.limits._sigma_SI_width,2)
                        else:
                            logger.error('You have to set up the sigma_SI(_width) to a positive value to use gaussian likelihood!')
                    elif likelihood == 'user':
                        #
                        # HERE ADD YOUR OWN LOG LIKELIHOOD FOR SI
                        #
                        chi+=0

                    elif likelihood == '':
                        chi+=0
                    else:
                        logger.error('You are not using a valid likelihood function for SI direct detection. Omitting the contribution!')

            #direct detection (SD) proton and neutron
            if obs.startswith('directSD') and self.maddm_run.mode['direct']:
                nucleon = obs.split('_')[1]
                if likelihood=='half_gauss':
                    if nucleon =='p':
                        if spinSDp > self.maddm_run.limits.SD_max(mdm, 'p'):
                            chi+= -0.5*pow(np.log10(self.maddm_run.limits.SD_max(mdm, 'p')/spinSDp),2)\
                                  /pow(self.options['half_gauss_width']['spinSDp'],2)
                    elif nucleon == 'n':
                        if spinSDp > self.maddm_run.limits.SD_max(mdm,'n'):
                            chi+= -0.5*pow(np.log10(self.maddm_run.limits.SD_max(mdm, 'n')/spinSDn),2)\
                                  /pow(self.options['half_gauss_width']['spinSDn'],2)
                elif likelihood == 'gaussian':
                    if nucleon == 'p' and self.maddm_run.limits._sigma_SDp > 0:
                        chi+=  -0.5*pow(spinSDp - self.maddm_run.limits._sigma_SDp,2)/pow(self.maddm_run.limits._sigma_SDp_width,2)
                    else:
                        logger.error('You have to set up the sigma_SDp(_width) to a positive value to use gaussian likelihood!')
                    if nucleon == 'n' and self.maddm_run.limits.sigma_SDn > 0:
                        chi+=  -0.5*pow(spinSDn - self.maddm_run.limits.sigma_SDn,2)/pow(self.maddm_run.limits._sigma_SDn_width,2)
                    else:
                        logger.error('You have to set up the sigma_SDp(_width) to a positive value to use gaussian likelihood!')
                elif likelihood == 'user':
                    #
                    # HERE ADD YOUR OWN LOG LIKELIHOOD FOR SD
                    #
                    chi+=0
                elif likelihood == '':
                    chi+=0
                else:
                    logger.error('You are not using a valid likelihood function for SD direct detection. Omitting the contribution!')

            #indirect detection
            if obs == 'indirect' and self.maddm_run.mode['indirect']:

                id_vel = self.maddm_run.maddm_card['vave_indirect']

                for channel in [k for k in results.keys() if k.startswith('taacsID#')]:

                        channel = channel.split('#')[1]
                        finalstate = channel.split('_')[1]

                        BR = sigmavID[channel]/sigmavID['tot']
                        if finalstate in self.maddm_run.limits._allowed_final_states:
                            if self.maddm_run.limits._id_limit_vel[finalstate] == id_vel\
                                    and self.maddm_run.limits._id_limit_sigv:

                                if likelihood=='half_gauss':
                                    if sigmavID[channel] > self.maddm_run.limits.ID_max(mdm,finalstate)/max(BR, 1e-10):
                                        chi += -0.5*pow(np.log10(self.maddm_run.limits.ID_max(mdm,finalstate)\
                                                             /(max(BR, 1e-10)*sigmavID[channel])),2)/pow(0.01,2)
                                elif likelihood=='gaussian':
                                    chi +=  -0.5*pow(self.maddm_run.limits._sigma_ID[finalstate] - self.maddm_run.limits.sigma_SDn,2)\
                                                /pow(self.maddm_run.limits._sigma_ID_width[finalstate],2)
                                elif likelihood =='user':
                                    #
                                    # HERE ADD YOUR OWN LOG LIKELIHOOD FOR ID
                                    #
                                    chi +=0
                                elif likelihood == '':
                                    chi+=0
                            else:
                                logger.warning('Omitting the likelihood contribution for channel %s.' % finalstate)
                                logger.warning('The limit is not set or velocities mismatch!')
                        else:
                            logger.warning('No limit for channel %s. Omitting the likelihood contribution.' % finalstate)

        #Print the counter on the screen

        if self.counter % 100 ==0:
            toprint = 'Scanned over %i points.' % self.counter
            logger.info(toprint)
        self.counter += 1

        return chi



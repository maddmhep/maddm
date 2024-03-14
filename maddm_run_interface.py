from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
import collections
import math
import logging
from multiprocessing import process
import os
import re
import sys
import subprocess

from numpy.core.fromnumeric import var
from . import auxiliary as aux
import stat
import shutil
from copy import copy, deepcopy

from . import MGoutput
from madgraph import MadGraph5Error
from models import check_param_card
import six
from six.moves import map
from six.moves import range
from six.moves import zip

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

logger = logging.getLogger('madgraph.plugin.maddm')

try:
    import pymultinest
except ImportError:
    pass

try:
    from scipy.interpolate import interp1d
    from scipy.integrate import quad
    from scipy.optimize import brute, fmin, minimize_scalar, bisect
    from scipy.special import gammainc
except ImportError as error:
    print(error)
    logger.warning('scipy module not found! Some Indirect detection features will be disabled.')
    HAS_SCIPY = False 
else:
    HAS_SCIPY = True

try:
    import numpy as np 
except ImportError:
    logger.warning('numpy module not found! Indirect detection features will be disabled.')
    HAS_NUMPY = False
else:
    HAS_NUMPY = True

if HAS_NUMPY and HAS_SCIPY:
    from .jfactor_gc import PROFILES, jfactor_gc
        
class ModuleMissing(Exception): pass

pjoin = os.path.join

logger_tuto = logging.getLogger('tutorial_plugin')
#logger.setLevel(10) #level 20 = INFO

MDMDIR = os.path.dirname(os.path.realpath( __file__ ))


#Is there a better definition of infinity?
__infty__ = float('inf')
__mnestlog0__ = -1.0E90

class ExpConstraints:

    def __init__(self):

        self._settable = ['1.-1', '2.-2', '3.-3', '4.-4', '21.21', '5.-5', '6.-6', '11.-11', '13.-13', '15.-15', '24.-24', '23.23', '25.25']

        self._oh2_planck = 0.1200 # from 1807.06209 Tab. 2 (TT,TE,EE+lowE+lensing) quoted also in abstract
        self._oh2_planck_width = 0.0012

        # self._dd_si_limit_file = pjoin(MDMDIR, 'ExpData', 'Xenon1T_data_2018.dat')
        # self._dd_sd_proton_limit_file = pjoin(MDMDIR, 'ExpData', 'Pico60_sd_proton.dat') # <---------CHANGE THE FILE!!!
        # self._dd_sd_neutron_limit_file = pjoin(MDMDIR, 'ExpData', 'Lux_2017_sd_neutron.dat')

        self._dd_limit_file = {
            'si'        : pjoin(MDMDIR, 'ExpData', 'Xenon1T_data_2018.dat'),
            'sd_proton' : pjoin(MDMDIR, 'ExpData', 'Pico60_sd_proton.dat'),
            'sd_neutron': pjoin(MDMDIR, 'ExpData', 'Lux_2017_sd_neutron.dat')
        }

        self._dd_limit_mdm = dict()
        self._dd_limit_sigma = dict()

        self._id_limit_file = {
            # cont final states
            '1.-1'  : pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_qq.dat'), # same qqx
            '2.-2'  : pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_qq.dat'), # same qqx
            '3.-3'  : pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_qq.dat'), # same qqx
            '4.-4'  : pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_cc.dat'),
            '5.-5'  : pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_bb.dat'),
            '6.-6'  : pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_tt.dat'),
            '11.-11': pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_ee.dat'),
            '13.-13': pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_mumu.dat'),
            '15.-15': pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_tautau.dat'),
            '21.21' : pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_gg.dat'),
            '23.23' : pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_ZZ.dat'),
            '24.-24': pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_WW.dat'),
            '25.25' : pjoin(MDMDIR, 'ExpData', 'MadDM_Fermi_Limit_hh.dat'),
            # line final states = <final state>_<experiment label><year>R<roi>
            '22.22_hess2013R1'  : pjoin(MDMDIR,'ExpData', 'hess_I_2013_einasto.dat'),
            '22.22_hess2016R1'  : pjoin(MDMDIR,'ExpData', 'hess_2016_einasto.dat'),
            '22.22_hess2018R1'  : pjoin(MDMDIR,'ExpData', 'HESS_2018_lines_R1_Einasto.dat'), # <sigma v> limits
            '22.22_fermi2015R3' : pjoin(MDMDIR,'ExpData', 'Fermi_2015_lines_R3_NFWcontracted.dat'), # <sigma v> limits
            '22.22_fermi2015R16': pjoin(MDMDIR,'ExpData', 'Fermi_2015_lines_R16_Einasto.dat'), # <sigma v> limits
            '22.22_fermi2015R41': pjoin(MDMDIR,'ExpData', 'Fermi_2015_lines_R41_NFW.dat'), # <sigma v> limits
            '22.22_fermi2015R90': pjoin(MDMDIR,'ExpData', 'Fermi_2015_lines_R90_Isothermal.dat')} # <sigma v> limits

        self._id_limit_file_flux = {
            # line final states
            '22.22_hess2018R1'    : pjoin(MDMDIR,'ExpData', 'HESS_2018_lines_R1_Einasto.dat'), # Flux limits
            '22.22_fermi2015R3'   : pjoin(MDMDIR,'ExpData', 'Fermi_2015_lines_R3_NFWcontracted.dat'), # Flux limits
            '22.22_fermi2015R16'  : pjoin(MDMDIR,'ExpData', 'Fermi_2015_lines_R16_Einasto.dat'), # Flux limits
            '22.22_fermi2015R41'  : pjoin(MDMDIR,'ExpData', 'Fermi_2015_lines_R41_NFW.dat'), # Flux limits
            '22.22_fermi2015R90'  : pjoin(MDMDIR,'ExpData', 'Fermi_2015_lines_R90_Isothermal.dat') # Flux limits
        }
        
        self._id_limit_vel = {
            # cont velocity
            '1.-1':2.0E-5,
            '2.-2':2.0E-5,
            '3.-3':2.0E-5,
            '4.-4':2.0E-5,
            '21.21':2.0E-5,
            '5.-5':2.0E-5,
            '6.-6':2.0E-5,
            '11.-11':2.0E-5,
            '13.-13':2.0E-5,
            '15.-15':2.0E-5,
            '24.-24':2.0E-5,
            '23.23':2.0E-5,
            '25.25':2.0E-5,
            # line velocity
            '22.22_hess2013'    : 1e-3,
            '22.22_hess2016'    : 1e-3,
            '22.22_hess2018R1'  : 1e-3,
            '22.22_fermi2015R3' : 1e-3,
            '22.22_fermi2015R16': 1e-3,
            '22.22_fermi2015R41': 1e-3,
            '22.22_fermi2015R90': 1e-3
        }

        self._id_limit_mdm = dict()
        self._id_limit_sigv = dict()
        self._id_limit_flux = dict()

        #In case there is a measurement of the cross seciton
        self._sigma_SI = -1.0
        self._sigma_SI_width = -1.0
        self._sigma_SDp = -1.0
        self._sigma_SDn = -1.0
        self._sigma_SDp_width = -1.0
        self._sigma_SDn_width = -1.0

        self._sigma_ID = dict()
        self._sigma_ID_width = dict()
        for item in self._id_limit_file.keys():
            self._sigma_ID[item] = -1.0
            self._sigma_ID_width[item] = -1.0

        if HAS_NUMPY:
            self.load_constraints()


        logger.info('Loaded experimental constraints. To change, use the set command')
        logger.info('Omega h^2 = %.4e +- %.4e' %(self._oh2_planck, self._oh2_planck_width))
        if HAS_NUMPY:
            logger.info('Spin Independent cross section: %s' % self._dd_limit_file['si'])
            logger.info('Spin Dependent cross section (p): %s' % self._dd_limit_file['sd_proton'])
            logger.info('Spin Dependent cross section (n): %s' % self._dd_limit_file['sd_neutron'])
        else:
            logger.info('Spin (in)dependent not available due to the missing python module: numpy')
#        for chan in self._allowed_final_states:
#            logger.info('Indirect Detection cross section for final state %s at velocity %.2e: %s'\
#                        % (chan, self._id_limit_vel[chan] ,self._id_limit_file[chan]))

    def safe_load_file(self, filename, variables, key, unpack = True, usecols = None, comments = '#'):
        try:
            temp = np.loadtxt(filename, unpack = unpack, usecols = usecols, comments = comments)
            for i, var in enumerate(variables):
                selfvar = getattr(self, var)
                selfvar[key] = temp[i]
        except Exception: # generic error from reading
            logger.warning("Can't read file '%s': missing or corrupted file." % filename)

    def load_constraints(self):
        #Load in direct detection constraints
        for spin_type, limit_file in six.iteritems(self._dd_limit_file):
            self.safe_load_file(filename = limit_file, variables = ['_dd_limit_mdm', '_dd_limit_sigma'], key = spin_type)
                                                                                                            
        for channel, limit_file in six.iteritems(self._id_limit_file):
            self.safe_load_file(filename = limit_file, variables = ['_id_limit_mdm', '_id_limit_sigv'], key = channel, usecols = (0,1))

        for channel, limit_file in six.iteritems(self._id_limit_file_flux):
            self.safe_load_file(filename = limit_file, variables = ['_id_limit_mdm', '_id_limit_flux'], key = channel, usecols = (0,2))

    #Returns a value in cm^2
    def SI_max(self, mdm):
        if not HAS_NUMPY:
            logger.warning("missing numpy module for SI limit")
            return -1
        
        if (mdm < np.min(self._dd_limit_mdm['si']) or mdm > np.max(self._dd_limit_mdm['si'])):
            logger.warning('Dark matter mass value '+str(mdm)+' is outside the range of SI limit')
            return -1
        else:
            return np.interp(mdm, self._dd_limit_mdm['si'], self._dd_limit_sigma['si'])

    #Returns a value in cm^2
    def SD_max(self, mdm, nucleon):
        if not HAS_NUMPY:
            logger.warning("missing numpy module for SD limit")
            return -1
            
        if nucleon not in ['n','p']:
            logger.error('nucleon can only be p or n')
            return -1
        else:
            if nucleon=='p':
                if (mdm < np.min(self._dd_limit_mdm['sd_proton']) or mdm > np.max(self._dd_limit_mdm['sd_proton'])):
                    logger.warning('Dark matter mass value '+str(mdm)+' is outside the range of SD limit')
                    return -1
                else:
                    return np.interp(mdm, self._dd_limit_mdm['sd_proton'], self._dd_limit_sigma['sd_proton'])
            elif nucleon=='n':
                if (mdm < np.min(self._dd_limit_mdm['sd_proton']) or mdm > np.max(self._dd_limit_mdm['sd_proton'])):
                    logger.warning('Dark matter mass value '+str(mdm)+' is outside the range of SD limit')
                    return -1
                else:
                    return np.interp(mdm, self._dd_limit_mdm['sd_neutron'], self._dd_limit_sigma['sd_neutron'])

    # Returns a value in cm^3/s
    def ID_max(self, mdm, channel, sigmav = True):
        ''' get the upper limit either for sigmav (sigmav = True) or flux (sigmav = False). '''
        if not HAS_NUMPY:
            logger.warning("missing numpy module for ID limit")
            return -1
        id_limit_dict = self._id_limit_sigv if sigmav else self._id_limit_flux
        try:
            if (mdm < np.min(self._id_limit_mdm[channel]) or mdm > np.max(self._id_limit_mdm[channel])):
                logger.warning('Dark matter mass value %.2e for channel %s is outside the range of ID limit' % (mdm, channel))
                return -1
            else:
                return np.interp(mdm, self._id_limit_mdm[channel], id_limit_dict[channel])
        except KeyError:
            logger.debug("Channel '%s' not found." % channel)
            return -1

################################################################################
##    Spectra
################################################################################
class Spectra:
    """This class holds some functionalities to load spectra from the PPPC4DMID 
    files or the output from pythia8, interpolates the PPPC4DMID spectra for 
    arbitrary values of DM etc.
    """
    
    PPPCdata = {}
    PPPC_type = None
    

    def __init__(self):
        
        self.spectra_id    = {'px':'antiprotons', 'gx':'gammas', 'nuex':'neutrinos_e', 'numux':'neutrinos_mu', 'nutaux':'neutrinos_tau', 'ex':'positrons'}

        # this contains the values in x and dndlogx as extracted from PPPC or generated by pythia8
        self.spectra       = {'log10x': [], 'x':[], 'antiprotons':[], 'gammas':[], 'neutrinos_e':[], 'neutrinos_mu':[], 'neutrinos_tau':[], 'positrons':[] }
        self.errors       = {'antiprotons':[], 'gammas':[], 'neutrinos_e':[], 'neutrinos_mu':[], 'neutrinos_tau':[], 'positrons':[] }

        # this contains the spectrua converted to E and dN/dE (needed by DRAGON)
        self.flux_source   = {'e':[], 'positrons': {'dNdE':[]}, 'antiprotons': {'dNdE':[] } , 
                              'neutrinos_tau':{'dNdE':[]},'neutrinos_mu':{'dNdE':[] }, 
                              'neutrinos_e':{'dNdE':[]} , 'gammas':{'dNdE':[]} }

        # fluxes_earth only filled if using PPPC4DMID method (no reading from Dragon output)
        self.flux_earth_positrons    = {'e':[] , 'positrons': { 'dPhidlogE':[] } }
        
        # this contains the fluxes at earth, and the values of the dndlogx after oscillations for neutrinos
        self.flux_earth              = {'e':self.flux_source['e'] , 'neutrinos_tau': {'dPhidE':[],'dndlogx':[] }, 'neutrinos_mu': {'dPhidE':[],'dndlogx':[] },
                                                 'neutrinos_e'  : {'dPhidE':[],'dndlogx':[] }, 'gammas'      : {'dPhidE':[],'dndlogx':[] } }

        self.channels      = ['ee', 'mumu', 'tautau', 'qq', 'cc', 'bb', 'tt', 'ZZ', 'WW', 'hh', 'gammagamma', 'gg']

        # self.map_allowed_final_state_PPPC = {'qqx':'qq', 'ccx':'cc', 'gg':'gg', 'bbx':'bb', 'ttx':'tt',
                                             # 'emep':'ee', 'mummup':'mumu', 'tamtap':'tautau', 'wpwm':'WW', 'zz':'ZZ', 'hh':'hh' }
        self.map_allowed_final_state_PPPC = {'1.-1':'qq', '2.-2':'qq', '3.-3':'qq', '4.-4':'cc', '21.21':'gg', '5.-5':'bb', '6.-6':'tt',
                                             '11.-11':'ee', '13.-13':'mumu', '15.-15':'tautau', '24.-24':'WW', '23.23':'ZZ', '25.25':'hh' }

    def check_mass(self,mdm):
        if (mdm < 5.0 or mdm > 100000):
            logger.error('DM mass outside the range of available PPPC spectra. Please run spectra generation with pythia8') # Break and ask the user to download the Tables
            return 
        else: return True 

    def load_tabulated_sourcespectra(self, PPPCDIR, kind = '',load = True):
        """routine to load the PPPC table if installed already"""
        
        if not load:
            return
        
        if Spectra.PPPC_type == ('source', kind):
            return Spectra.PPPCdata
        
        if (not os.path.isfile(PPPCDIR+'/PPPC_Tables_EW.npy') and not os.path.isfile(PPPCDIR+'/PPPC_Tables_noEW.npy')):
            logger.error('PPPC4DMID Spectra at source not found! Please install by typing install PPPC4DMID') # Break and ask the user to download the Tables       
            return

        if (os.path.isfile(PPPCDIR+'/PPPC_Tables_EW.npy') and os.path.isfile(PPPCDIR+'/PPPC_Tables_noEW.npy') and not os.path.isfile(PPPCDIR+'/CosmiXs.npy')):
            logger.error("PPPC4DMID is not updated, please update it with 'install PPPC4DMID'")
            return

        if kind == 'ew':
            dic =  np.load(PPPCDIR+'/CosmiXs.npy', allow_pickle= True, encoding="latin1").item()
            logger.info('CosmiXs spectra at source (updated version of PPPC) loaded')
        elif kind == 'pppc_ew':
            dic =  np.load(PPPCDIR+'/PPPC_Tables_EW.npy', allow_pickle= True, encoding="latin1").item()
            logger.info('PPPC4DMID Spectra at source (with EW corrections) loaded')
        else:
            dic =  np.load(PPPCDIR+'/PPPC_Tables_noEW.npy', allow_pickle = True, encoding="latin1").item()
            kind = "pppc"
            logger.info('PPPC4DMID Spectra at source loaded')

        Spectra.PPPCdata = dic
        Spectra.PPPC_type = ('source', kind)
        return dic

    def load_PPPC_earth(self, PPPCDIR, prof = 'Ein'):
        
        if Spectra.PPPC_type == ('earth', prof):
            return Spectra.PPPCdata
        
        if (not os.path.isfile(PPPCDIR+'/PPPC_Tables_epEarth_'+prof+'.npy') ):
            logger.error('PPPC4DMID Spectra at Earth not found! Please install by typing install PPPC4DMID') # Break and ask the user to download the Tables                     
            return
        
        else:
            dic = np.load(PPPCDIR+'/PPPC_Tables_epEarth_'+prof+'.npy' , allow_pickle= True, encoding="latin1").item()
            Spectra.PPPCdata = dic
            Spectra.PPPC_type = ('earth', prof)
            
            return dic 

    # this function extracts the values of the spectra interpolated linearly between two values mdm_1 and mdm_2                                              
    # mdm is the DM candidate mass, spectrum is gammas, positron etc, channel is the SM annihilation e.g. bbar, hh etc.                        
    def interpolate_spectra(self, sp_dic, mdm = '', spectrum = '' , channel = '' , earth = False , prof = 'Ein' , prop = 'MED' , halo_func = 'MF1'):
        
        if not HAS_SCIPY:
            raise ModuleMissing('scipy module is required for this functionality.')
        
        M   = sp_dic['Masses']
        dm_min =  max([m for m in M if m <= mdm])  # extracting lower mass limit to interpolate from                                                                  
        dm_max =  min([m for m in M if m >= mdm])  # extracting upper mass limit to interpolate from                                                                    

        interpolated = []
        if not earth: # two spectra: source and earth. First method is for source, second is for spectra at earth (when earth == True)
            key = channel
        elif earth:
            key = channel+'_'+prof+'_'+prop+'_'+halo_func

        if dm_min == dm_max :
            return sp_dic[spectrum][ str(dm_min) ][key]
        spec_1 = sp_dic[spectrum][ str(dm_min) ][key]
        spec_2 = sp_dic[spectrum][ str(dm_max) ][key]
     
        for x in range(len(spec_1)): # the spectrum values are ordered as 'x' vaues extracted from the Log[10,x] in the PPPC Tables                       
            interp_function = interp1d([dm_min, dm_max], [spec_1[x],spec_2[x]] )
            value =  interp_function(mdm)
            interpolated.append(value)

        return interpolated

    def initialize_spectra(self):
        for k in self.spectra.keys():
            self.spectra[k] = []
      
################################################################################                                                                                            
##    Fermi                                                                                                                                                               
################################################################################ 
class Fermi_bounds:
    """This class holds all the functionalities needed to evaluate the Fermi LAT ul on the annihilation cross section,
    as well as the likelihood and pvalue of the model tested.
    """
    
    def __init__(self):

        self.nBin = 24
        self.j0 = 3.086e21 # convention spectra    
        self.dSph_jfac_file     =  pjoin(MDMDIR, 'Fermi_Data', 'Jfactors.dat')
        self.dSph_ll_files_path =  pjoin(MDMDIR, 'Fermi_Data', 'likelihoods')
        self.dwarves_list = ['coma_berenices', 'draco', 'segue_1', 'ursa_major_II', 'ursa_minor', 'reticulum_II' ] # dSphs with the 6 highest Jfactors
        self.dwarveslist_all = self.extract_dwarveslist() 
        if HAS_NUMPY and HAS_SCIPY:
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
                    
        return dw_dic_coll

    
    def eflux(self,spectrum, emin=1e2, emax=1e5, quiet=False):
        """ Integrate a generic spectrum, multiplied by E, to get the energy flux.                                                                                             
        """
        
        if not HAS_SCIPY:
            raise ModuleMissing('scipy module is required for this functionality.')
        
        espectrum = lambda e: spectrum(e)*e
        tol = min(espectrum(emin),espectrum(emax))*1e-10
        try:
            return quad(espectrum,emin,emax,epsabs=tol,full_output=True)[0]
        except Exception as msg:
            logger.info('Numerical error "%s" when calculating integral flux.' % msg)
        return np.nan

    def marg_like_dw(self,dw_in_i,pred,marginalize):
      
        if not HAS_SCIPY:
            raise ModuleMissing('scipy module is required for this functionality.')
        
        j0, nBin = self.j0 , self.nBin

        j,jerr,like_inter = dw_in_i['Jfac'], dw_in_i['Jfac_err'], dw_in_i['likelihood']
       
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

        if not HAS_SCIPY:
            raise ModuleMissing('scipy module is required for this functionality.')

        dw_in = self.dw_in
        ll_tot = 0.0
        ll_null = 0.0
        j_factors = {}
        for k,v in six.iteritems(dw_in):
            marg_like = self.marg_like_dw(v,pred,marginalize)
            ll_tot += marg_like[0]
            ll_null += v['ll_null']
            j_factors[k] = marg_like[1]

        pval = self.compute_pvalue(ll_tot,ll_null)
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

    # this function return the Fermi experimental UL on sigmav if the argument sigmav_th= False;
    # otherwise it calculates the likelihood and p-value for the given point 
    def Fermi_sigmav_lim(self, mDM, x = '' , dndlogx = '' , marginalize = True, sigmav_th = False , maj_dirac='', \
                               sigmavmin=1e-35, sigmavmax=1e-15, step_size_scaling=1.0, cl_val = 0.95):
        stop = False
        if not HAS_NUMPY:
            logger.warning("Fermi limit ignored due to missing numpy module")
            stop = True
        if not HAS_SCIPY:
            logger.warning("Fermi limit ignored due to missing scipy module")
            stop = True
            
        if stop:
            if sigmav_th:
                return -1, -1
            else:
                return -1      
        
        np.seterr(divide='ignore', invalid='ignore')   # Keep numpy from complaining about dN/dE = 0...                                                                
        j0 , nBin = self.j0 , self.nBin
        
        dw_in = self.dw_in
        sigmav0 = 1e-26

        emins = [0.5, 0.666760716, 0.8891397050000001, 1.1856868500000002, 1.5811388300000002, 2.10848252, 2.81170663, 3.7494710500000004, 5.0, 6.667607159999999, 
                 8.89139705 , 11.856868500000001, 15.8113883, 21.0848252, 28.117066299999998, 37.494710500000004, 50.0, 66.6760716, 88.91397049999999, 
                 118.568685, 158.11388300000002, 210.848252, 281.170663, 374.94710499999997]
        emaxs = [0.666760716, 0.8891397050000001, 1.1856868500000002, 1.5811388300000002, 2.10848252, 2.81170663, 3.7494710500000004, 5.0, 6.667607159999999, 
                 8.89139705, 11.856868500000001, 15.8113883, 21.0848252, 28.117066299999998, 37.494710500000004, 50.0, 66.6760716, 88.91397049999999, 
                 118.568685, 158.11388300000002, 210.848252, 281.170663, 374.94710499999997, 500.0]

        logx = np.log10(x) 
        energy = mDM*10**logx
        dnde = (dndlogx/(mDM*10**logx*2.30259))
        log_energy = np.log10(energy)
        log_dnde = np.log10(dnde)
        log_interp = interp1d(log_energy,log_dnde)
        spectrum = lambda e: np.nan_to_num(10**( log_interp(np.log10(e)) )) if (e <= energy.max() and e >= energy.min()) else 0.0

        pred = np.array([self.eflux(spectrum,e1,e2)/mDM**2*sigmav0*j0/(2.*maj_dirac*math.pi) for e1,e2 in zip(emins,emaxs)])

        if not sigmav_th:
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

             pred_sigma = pred*sigmav_ul/sigmav0
             result  = self.res_tot_dw(pred_sigma,marginalize)
             p_value = result[2]     

             if p_value <= cl_val*0.98 or p_value >= cl_val*1.02:
                 sigmav_ul= -1
                 print(" WARNING: increase range (sigmavmin,sigmavmax) and/or step_size_scaling!")
        
             return sigmav_ul    

        elif sigmav_th:
             pred_sigma = pred*sigmav_th/sigmav0
             result = self.res_tot_dw(pred_sigma,marginalize)
             return result[2] , result[0]

class RegionOfInterest(object):
    ''' Definition of a region of interest (ROI) with default profile, masks and
        instructions on J-factor computation.

        Parameters
        ----------
        amplitude : float
            The angle of the ROI amplitude in deg, named :math:`\alpha_1`.
        default_profile : jfactor_gc.DensityProfile subclass
            The default density profile of the ROI, the one it is optimized to.
        upper_limit_label : str
            The label associated to the upper limit file to be used for this ROI.
        jfactor_strategy : str, default="normal"
            The strategy computation for the J-factor, choose between
            `"normal"` or `"inverted"`.
        mask_latitude : float, default=0.
            The angle related to the mask over the latitude in deg, named
            :math:`\lambda`.
        mask_longitude : float, default=180.
            The angle related to the mask over the longitude in deg, named
            :math:`\beta`.
        mask_inner_amplitude : float, default=0.
            The angle related to the circle-shaped mask over the galactic center,
            named :math:`\alpha_1`.
        cached_jfactors : dict[jfactor_gc.DensityProfile subclass, float], optional
            The J-factor already knows for certain profiles related to this ROI.
        r_max : float, optional
            The maximum line of sight value to integrate over.
        likelihood : np.ndarray
            The likelihood data: a multidimensional array containing the following
            values: "Energy(GeV), UL_flux(cm^{-2} s^{-1}), x0 || A, sigma || 0" as
            columns.
    '''
    def __init__(self, amplitude, default_profile, upper_limit_label, jfactor_strategy="normal", mask_latitude=0., mask_longitude=180., mask_inner_amplitude=0., cached_jfactors={}, r_max=np.inf, likelihood=None):
        self._amplitude = amplitude
        self._mask_latitude = mask_latitude
        self._mask_longitude = mask_longitude
        self._mask_inner_amplitude = mask_inner_amplitude
        self.default_profile = default_profile
        self.cached_jfactors = cached_jfactors
        self.jfactor_strategy = jfactor_strategy
        self._r_max = r_max
        self.upper_limit_label = upper_limit_label
        self.likelihood = likelihood

    def to_rad(self, ang):
        return ang * np.pi/180.

    def jfactor(self, profile):
        ''' Get the J-factor for the roi specified.
            First, it checks whether the profile has already a J-factor computed, in case, it directly returns that one.
            Second, it checks if the new profile and the default profile are equal excepts for the normalization (it also keep into account r_sun).
            in case it is True and the default J-factor is cached, it rescales the latter according to the new normalization,
            in case it is False, it computes the new J-factor and caches it for the profile specified.
            In case the new profile is equal but normalisation to the default profile, but the default J-factor is not cached, it first
            computes and cache the default J-factor.
        '''
        if profile in self.cached_jfactors: # search in the cached_j dict
            return self.cached_jfactors[profile]
        if profile.eq_but_norm(self.default_profile):
            if self.default_profile not in self.cached_jfactors:
                # compute the default J-factor
                self.cached_jfactors[self.default_profile] = jfactor_gc(
                                        ang_2    = self.to_rad(self._amplitude),
                                        ang_lat  = self.to_rad(self._mask_latitude),
                                        ang_long = self.to_rad(self._mask_longitude),
                                        ang_1    = self.to_rad(self._mask_inner_amplitude),
                                        profile  = self.default_profile,
                                        strategy = self.jfactor_strategy,
                                        r_max    = self._r_max,
                                        mask     = True
                                    )
            return np.power(profile.rho_s/self.default_profile.rho_s, 2) * self.cached_jfactors[self.default_profile]
        # otherwise compute and add to cache
        self.cached_jfactors[profile] = jfactor_gc(
                                ang_2    = self.to_rad(self._amplitude),
                                ang_lat  = self.to_rad(self._mask_latitude),
                                ang_long = self.to_rad(self._mask_longitude),
                                ang_1    = self.to_rad(self._mask_inner_amplitude),
                                profile  = self.default_profile,
                                strategy = self.jfactor_strategy,
                                r_max    = self._r_max,
                                mask     = True
                            )
        return self.cached_jfactors[profile]

class RegionOfInterestList(object):
    ''' Handle the RegionOfInterest lists
    '''
    def __init__(self, start_list):
        self._rois = []
        self._rois_amplitudes = {}
        for element in start_list:
            self.append(element)

    def append(self, roi):
        if not isinstance(roi, RegionOfInterest):
            raise ValueError("Only objects of type RegionOfInterest can be appended to RegionOfInterestList.")
        self._rois.append(roi)
        self._rois_amplitudes[roi._amplitude] = len(self._rois) - 1

    def __iter__(self):
        return iter(self._rois)

    def __getitem__(self, index):
        if isinstance(index, (float, int)):
            return self._rois[self._rois_amplitudes[index]]
        raise NotImplementedError("Can not accept index of " + str(type(index)).replace("<","").replace(">",""))

    def __len__(self):
        return len(self._rois)

    def __copy__(self):
        cls_ = self.__class__
        newobj = cls_.__new__(cls_)
        newobj.__dict__.update(self.__dict__)
        return newobj

    def __deepcopy__(self, memo):
        cls_ = self.__class__
        newobj = cls_.__new__(cls_)
        memo[id(self)] = newobj
        for k, v in self.__dict__.items():
            setattr(newobj, k, deepcopy(v, memo))
        return newobj

line_experiments_regions_of_interest = {
    "Fermi-LAT_2015": RegionOfInterestList([
        RegionOfInterest(
            amplitude=3.,
            default_profile=PROFILES.NFW(r_s = 20.0, gamma = 1.3, rho_sun = 0.4, r_sun = 8.5),
            upper_limit_label="(22)(22)_fermi2015R3",
            jfactor_strategy="normal",
            mask_latitude=5.,
            mask_longitude=6.,
            mask_inner_amplitude=0.,
            cached_jfactors={PROFILES.NFW(r_s = 20.0, gamma = 1.3, rho_sun = 0.4, r_sun = 8.5) : 1.497e+23},
            likelihood=np.loadtxt(pjoin(MDMDIR, 'Fermi_line_likelihoods', 'R3_gamma_lines_ULflux_like.dat'), unpack = True),
            r_max=np.inf
        ),
        RegionOfInterest(
            amplitude=16.,
            default_profile=PROFILES.Einasto(r_s = 20.0, alpha = 0.17, rho_sun = 0.4, r_sun = 8.5),
            upper_limit_label="(22)(22)_fermi2015R16",
            jfactor_strategy="normal",
            mask_latitude=5.,
            mask_longitude=6.,
            mask_inner_amplitude=0.,
            cached_jfactors={PROFILES.Einasto(r_s = 20.0, alpha = 0.17, rho_sun = 0.4, r_sun = 8.5) : 9.39e+22},
            likelihood=np.loadtxt(pjoin(MDMDIR, 'Fermi_line_likelihoods', 'R16_gamma_lines_ULflux_like.dat'), unpack = True),
            r_max=np.inf
        ),
        RegionOfInterest(
            amplitude=41.,
            default_profile=PROFILES.NFW(r_s = 20.0, gamma = 1.0, rho_sun = 0.4, r_sun = 8.5),
            upper_limit_label="(22)(22)_fermi2015R41",
            jfactor_strategy="normal",
            mask_latitude=5.,
            mask_longitude=6.,
            mask_inner_amplitude=0.,
            cached_jfactors={PROFILES.NFW(r_s = 20.0, gamma = 1.0, rho_sun = 0.4, r_sun = 8.5) : 9.16e+22},
            r_max=np.inf
        ),
        RegionOfInterest(
            amplitude=90.,
            default_profile=PROFILES.Isothermal(r_s = 5.0, rho_sun = 0.4, r_sun = 8.5),
            upper_limit_label="(22)(22)_fermi2015R90",
            jfactor_strategy="normal",
            mask_latitude=5.,
            mask_longitude=6.,
            mask_inner_amplitude=0.,
            cached_jfactors={PROFILES.Isothermal(r_s = 5.0, rho_sun = 0.4, r_sun = 8.5) : 6.94e+22},
            r_max=np.inf
        )
    ]),
    "HESS_2018": RegionOfInterestList([
        RegionOfInterest(
            amplitude=1.,
            default_profile=PROFILES.Einasto(r_s = 20.0, alpha = 0.17, rho_sun = 0.39, r_sun = 8.5),
            upper_limit_label="(22)(22)_hess2018R1",
            jfactor_strategy="normal",
            mask_latitude=0.3,
            mask_longitude=0.,
            mask_inner_amplitude=0.,
            cached_jfactors={PROFILES.Einasto(r_s = 20.0, alpha = 0.17, rho_sun = 0.39, r_sun = 8.5) : 4.66e21},
            r_max=np.inf
        ),
    ]),
}

class GammaLineExperiment(object):
    ''' this class holds all the generic functionalities regarding the upper limits on gamma ray line searches '''

    def __init__(self, name, energy_resolution, detection_range, roi_list, majorana_dirac_factor, check_profile_message = lambda is_optimized, roi: {True: "", False: ""}.get(is_optimized), arxiv_number = None, hidden = False):
        # experiment name
        self.name = name
        # energy resolution function
        self.energy_resolution = energy_resolution
        # detection range
        self.detection_range = detection_range
        # roi_list: a list of RegionOfInterest objects
        self.roi_list = roi_list
        # Majorana-Dirac factor: if the particle is auto-conjugate then this factor -> 1 (Majorana or neutral scalar), else -> 2 (Dirac or charged scalar)
        self.majorana_dirac_factor = majorana_dirac_factor
        # check profile output: function which gives the message to display when the profile not optimized for the ROI
        self.check_profile_message = check_profile_message
        # ArXiv number
        self.arxiv_number = arxiv_number
        # hidden?
        self.hidden = hidden

    def __str__(self):
        return self.name + " [arXiv:%s]" % self.arxiv_number

    def check_profile(self, roi, profile):
        ''' Check the profile with the default one for the ROI, returning True or False, which is then passed to the function check_profile_message. The result of is_optimized is returned as well. '''
        is_optimized = profile.is_optimized(self.roi_list[roi].default_profile)
        return is_optimized, self.check_profile_message(is_optimized, roi)

    def get_J(self, roi, profile):
        return self.roi_list[roi].jfactor(profile)

    def flux(self, mdm, sigmav, jfact, is_aa):
        ''' returns the integrated flux for the process dm dm > a X. If X != a, then takes half of the flux. '''
        coeff = 2. if is_aa else 1.
        return coeff * sigmav * jfact / (8 * np.pi * mdm * mdm * self.majorana_dirac_factor)

    def get_roi(self):
        return self.roi_list._rois_amplitude.keys()

    def get_name(self):
        return self.name

    def get_flux_ul(self, e_peak, roi, id_constraints):
        ''' Returns the upper limit on flux for a given ROI (expressed in degrees).
            The limits are taken from the ExpConstraint class, the interpolation function is evaluated at the energy of the peak.
        '''
        ul_label = self.roi_list[roi].upper_limit_label
        flux_ul = id_constraints.ID_max(mdm = e_peak, channel = ul_label, sigmav = False)
        return flux_ul

    def get_sigmav_ul(self, e_peak, roi, profile, id_constraints, is_aa):
        ''' Returns the upper limit on flux for a given ROI (expressed in degrees).
            The limits are taken from the ExpConstraint class, the interpolation function is evaluated at the energy of the peak.
            It returns np.inf if the profile is not compatible with default one for the ROI and if there are problems in ID_max,
            because that would be used for comparison between the other experiments. Eventually it would become -1 if also other
            constraints are invalid.
        '''
        default_profile = self.roi_list[roi].default_profile
        ul_label = self.roi_list[roi].upper_limit_label
        if profile.eq_but_norm(default_profile):
            coeff = 1. if is_aa else 2.
            sigmav_ul = id_constraints.ID_max(mdm = e_peak, channel = ul_label, sigmav = True) * np.power(default_profile.rho_s / profile.rho_s, 2) * coeff
            return sigmav_ul if sigmav_ul > 0 else np.inf
        else:
            return np.inf

class GammaLineExperimentsList(object):
    ''' Handle the GammaLineExperiments lists: automatically excludes the hidden experiments when iterating.
        Contains also a list with the ROI to analyse during this run.
    '''
    def __init__(self):
        self._experiments = []
        self._rois_to_analyse = []

    def append(self, experiment, roi_to_analyse):
        if not isinstance(experiment, GammaLineExperiment):
            raise ValueError("Only objects of type GammaLineExperiment can be appended to GammaLineExperimentsList.")
        self._experiments.append(experiment)
        self._rois_to_analyse.append(roi_to_analyse)

    def discard_hidden(self):
        return [exp for exp in self._experiments if not exp.hidden]

    def __iter__(self):
        return iter(self.discard_hidden())

    def iternames(self):
        return iter([exp.get_name() for exp in self])

    def iterrois(self):
        return iter([roi for roi, exp in zip(self._rois_to_analyse, self._experiments) if not exp.hidden])

    def __str__(self):
        return '[' + ', '.join([str(exp) for exp in self]) + ']'

    def __getitem__(self, index):
        return self.discard_hidden()[index]

    def __len__(self):
        return len(self.discard_hidden())

    def __copy__(self):
        cls_ = self.__class__
        newobj = cls_.__new__(cls_)
        newobj.__dict__.update(self.__dict__)
        return newobj

    def __deepcopy__(self, memo):
        cls_ = self.__class__
        newobj = cls_.__new__(cls_)
        memo[id(self)] = newobj
        for k, v in self.__dict__.items():
            setattr(newobj, k, deepcopy(v, memo))
        return newobj

class Fermi2015GammaLineLikelihood(object):
    ''' this class holds the functionality to compute the likelihood and the p-value for some ROI of the Fermi-LAT gamma-line 2015 analysis (arXiv: 1506.00013) 
        likelihoods were kindly provided by Alessandro Cuoco, who obtain them for R3 and R16 in an approximate way (see 1603.08228 for reference)
    '''
    def __init__(self, experiment):
        self.experiment = experiment

    def get_energy_bin(self, mdm, energy_means):
        ''' find the energy bin, by minimising abs(E_i - m) and returns the index '''
        if mdm < energy_means[0] or mdm > energy_means[-1]:
            logger.warning("Dark matter mass is out of the range allowed")
            return None # something
        this_bin = -1
        delta = -1
        for i_bin, e in enumerate(energy_means):
            if np.abs(e - mdm) < delta or this_bin == -1:
                delta = np.abs(e - mdm)
                this_bin = i_bin
        return this_bin

    def loglike_linear(self, flux, A):
        return A * flux

    def loglike_parabolic(self, flux, x_zero, sigma):
        ''' set loglike = 0 for flux < x_zero in order to avoid non-zero loglike for flux = 0 (and likelihood depending on mass in that case) '''
        return np.power(flux - x_zero, 2) / np.power(sigma, 2) if flux >= x_zero else 0.

    def loglike(self, peak, roi, profile):
        ''' returns the -2*log(like) evaluated at the flux, for a certain ROI (in degrees) '''
        try:
            like_file = self.experiment.roi_list[roi].likelihood
        except KeyError:
            raise ValueError("ROI of amplitude %f is not allowed, ignore likelihood computation" % roi)
        if like_file is None:
            raise IOError("No likelihood file available for ROI of amplitude %f." % roi)
        energy, _, x_zero, sigma = like_file[:4]
        this_bin = self.get_energy_bin(mdm = peak.e_peak, energy_means = energy)
        if not this_bin:
            # out of mass range
            return -1.
        if sigma[this_bin] == 0.:
            logger.debug("loglike ROI %d: Linear: m = %.3e, x_zero = %.3e, flux = %.3e" % (roi, peak.e_peak, x_zero[this_bin], peak.flux))
            return self.loglike_linear(flux = peak.flux, A = x_zero[this_bin])
        else:
            logger.debug("loglike ROI %d: Parabolic: m = %.3e, x_zero = %.3e, sigma = %.3e, flux = %.3e" % (roi, peak.e_peak, x_zero[this_bin], sigma[this_bin], peak.flux))
            return self.loglike_parabolic(flux = peak.flux, x_zero = x_zero[this_bin], sigma = sigma[this_bin])

    def compute_pvalue(self, peak, roi, profile):
        ''' compute the p-value assuming as null hypothesis the absence of signal (flux = 0), it assumes the chi-square distribution with one degree of freedom '''
        try:
            loglike_h1 = self.loglike(peak, roi, profile)
            loglike_h0 = self.loglike(GammaLineSpectrum.Peak(label = peak.label + '_H0', e_peak = peak.e_peak, full_width_half_max = peak.fwhm, shape = lambda e: 0., flux = 0.), roi, profile)
        except (ValueError, IOError):
            logger.warning("Error during likelihood computation, will not compute p-value.")
            return -1
        return 1 - gammainc(0.5, (loglike_h1 - loglike_h0)/2.) # loglike already have -2 factor


class GammaLineSpectrum(object):
    ''' class to handle gamma line spectra, peaks, peaks merging, fluxes 
        the spectrum is the flux spectrum, not gamma only, because experiments measure fluxes
        so we need to multiply the energy peak (approximately gaussian) for <sigma v>
    '''
    FWHM_MERGED_FRAC = 1.5 # coefficient to multiply to the FWHM after having done all the possible merging: it expresses the maximum fraction of FWHM (with respect to the experiment's one) allowed for the merged line to be considered a line
    ERROR_TEST = collections.OrderedDict([ # order the error function in the order of testing, errors are concatenated as strings
        ('1', 'not_a_line'),
        ('2', 'not_enough_separated')
    ])

    class Peak(object):
        ''' class to handle gamma peaks in the spectrum '''
        def __init__(self, label, e_peak, full_width_half_max, shape, flux):
            self.label   = label # the label associated to the peak i.e. the final state particles in pdg
            self.e_peak  = e_peak
            self.fwhm    = full_width_half_max
            self.flux    = flux
            self.flux_UL = __infty__
            self._shape  = shape
            self.error   = []

        def set_error(self, code_error):
            self.error.append(code_error)
        
        def get_error_num(self):
            return int("".join(self.error)) if len(self.error) != 0 else 0

        def __call__(self, energy):
            return self._shape(energy)

        def __eq__(self, other):
            ''' equality is guaranteed if the final state (aka the label) is the same '''
            assert isinstance(other, GammaLineSpectrum.Peak)
            return self.label == other.label

        def __ne__(self, other):
            return not self.__eq__(other)

        def __add__(self, other):
            assert isinstance(other, GammaLineSpectrum.Peak)
            label      = "%s+%s" % (self.label, other.label) # new label made of the sum of the various final states that have been merged.
            shape      = lambda e: self(e) + other(e)
            min_result = minimize_scalar(lambda e: - shape(e), method = 'bounded', bounds = tuple(sorted([self.e_peak, other.e_peak])))
            e_peak     = min_result.x
            max_value  = - min_result.fun
            max_fwhm   = self.fwhm + other.fwhm # in order to be sure that f(a) and f(b) have different sign, use the sum of the fwhms to find b from a
            fwhm       = np.abs(bisect(lambda e: shape(e) - max_value/2., a = min_result.x, b = min_result.x - max_fwhm) - bisect(lambda e: shape(e) - max_value/2., a = min_result.x, b = min_result.x + max_fwhm))
            return GammaLineSpectrum.Peak(
                label               = label, 
                e_peak              = e_peak, 
                full_width_half_max = fwhm, 
                shape               = shape,
                flux                = self.flux + other.flux
            )

        def __sub__(self, other):
            assert isinstance(other, GammaLineSpectrum.Peak)
            return self.e_peak - other.e_peak

        def __lt__(self, other):
            assert isinstance(other, GammaLineSpectrum.Peak) or isinstance(other, float) or isinstance(other, int)
            if isinstance(other, GammaLineSpectrum.Peak):
                return self.e_peak < other.e_peak
            if isinstance(other, float) or isinstance(other, int):
                return self.e_peak < other

        def __gt__(self, other):
            assert isinstance(other, GammaLineSpectrum.Peak) or isinstance(other, float) or isinstance(other, int)
            if isinstance(other, GammaLineSpectrum.Peak):
                return self.e_peak > other.e_peak
            elif isinstance(other, float) or isinstance(other, int):
                return self.e_peak > other

        def __le__(self, other):
            return self.__gt__(other)

        def __ge__(self, other):
            return self.__lt__(other)

        def __str__(self):
            return "Peak %(label)s (e_peak = %(e_peak).3e GeV, fwhm = %(fwhm).3e GeV" % vars(self) + ", flux = %.3e cm^-2 s^-1)" % self.flux

        def __copy__(self):
            cls_ = self.__class__
            newobj = cls_.__new__(cls_)
            newobj.__dict__.update(self.__dict__)
            return newobj

        def __deepcopy__(self, memo):
            cls_ = self.__class__
            newobj = cls_.__new__(cls_)
            memo[id(self)] = newobj
            for k, v in self.__dict__.items():
                setattr(newobj, k, deepcopy(v, memo))
            return newobj
    #######################################################

    def __init__(self, line_exp, pdg_particle_map, param_card, n_fwhm_separation, peak_height_factor):
        self.pdg_particle_map = pdg_particle_map
        self.param_card       = param_card
        self.mass_dict        = {}
        for p in pdg_particle_map.keys():
            try:
                self.mass_dict["%s" % p] = param_card.get_value("mass", abs(int(p)))
            except KeyError:
                continue
        self.line_exp           = line_exp
        self.spectrum           = []
        self.n_fwhm_separation  = n_fwhm_separation
        self.peak_height_factor = peak_height_factor

    def __iter__(self):
        return iter(self.spectrum)

    def __str__(self):
        return '[' + ', '.join([str(peak) for peak in self.spectrum]) + ']'

    def __getitem__(self, index):
        return self.spectrum[index]

    def __len__(self):
        return len(self.spectrum)

    def __copy__(self):
        cls_ = self.__class__
        newobj = cls_.__new__(cls_)
        newobj.__dict__.update(self.__dict__)
        return newobj

    def __deepcopy__(self, memo):
        cls_ = self.__class__
        newobj = cls_.__new__(cls_)
        memo[id(self)] = newobj
        for k, v in self.__dict__.items():
            setattr(newobj, k, deepcopy(v, memo))
        return newobj

    def append(self, element):
        assert isinstance(element, GammaLineSpectrum.Peak)
        return self.spectrum.append(element)

    def index(self, element):
        assert isinstance(element, GammaLineSpectrum.Peak)
        return self.spectrum.index(element)

    def sort(self, reverse = False, key = lambda item: item.e_peak):
        self.spectrum.sort(reverse = reverse, key = key)

    def e_peak_two_body(self, mdm, mx):
        ''' find energy of the peak for final state ax, with x = a, z, h '''
        return mdm * (1. - np.power(mx/2./mdm, 2))

    def e_peak_three_body(self, mdm, mx):
        pass # for future implementations

    def sigma(self, full_width_half_max):
        ''' standard deviation of a gaussian from a full width at half maximum '''
        return full_width_half_max/2/np.sqrt(2*np.log(2))

    def gaussian_spectra(self, e_peak, full_width_half_max, flux):
        ''' line spectrum with energy resolution (gaussian shape), the coeff depends on the final state (2 for aa, 1 for others).
            the normalization factor of the gaussian = 1, so that the shape computed at the e_peak is the flux.
        '''
        return lambda energy: flux * np.exp(-0.5 * np.power((energy-e_peak)/self.sigma(full_width_half_max), 2)) / (np.sqrt(2 * np.pi) * self.sigma(full_width_half_max))

    def find_lines(self, mdm, computed_sigmav, jfact):
        ''' find lines peaks, FWHM and shape of the spectrum. Append Peak objects with key final_state and _peak, _FWHM, _shape '''
        for fs, sigmav in computed_sigmav.items():
            particles = self.pdg_particle_map.find_particles(fs)
            # it's a gamma spectrum, so remove one 22 and get the other particles
            particles.remove('22')
            particles = list(set(particles))
            if len(particles) == 1:
                mx        = self.mass_dict[particles[0]]
                e_peak_fs = self.e_peak_two_body(mdm, mx)
            elif len(particles) == 2:
                # for future implementations
                raise Exception("Only 2-body final states at the moment, sorry, not '%s'" % self.pdg_particle_map.format_particles(fs))
            else:
                raise Exception("Can analyze only 2-body and 3-body final states, not '%s'" % self.pdg_particle_map.format_particles(fs))
            # check if the line is in the detection range of the experiment
            if e_peak_fs <= 0.:
                continue # exclude negative energy peaks, meaning of a forbidden process
            fwhm_fs  = self.line_exp.energy_resolution(e_peak_fs) * e_peak_fs
            flux     = self.line_exp.flux(mdm, sigmav, jfact, is_aa = fs == '22.22')
            shape_fs = self.gaussian_spectra(e_peak = e_peak_fs, full_width_half_max = fwhm_fs, flux = flux)
            self.append(GammaLineSpectrum.Peak(
                label               = fs,
                e_peak              = e_peak_fs,
                full_width_half_max = fwhm_fs,
                shape               = shape_fs,
                flux                = flux
            ))
        logger.debug(self)
        return self

    def merge_lines(self):
        ''' returns the final spectrum, made of merged lines '''
        comparisons = {}
        for peak_1 in self:
            for peak_2 in self[self.index(peak_1)+1:]:
                diff_energies = np.abs(peak_1 - peak_2)
                # this is a tuple: for each couple the first term is a tuple containing the peaks to sum, the second is a boolean to indicate if they can be merged and the third is a quantity to be used to decide among more couples which one should be merged first.
                # comparisons["%s--%s" % (peak_1.label, peak_2.label)] = ((peak_1, peak_2), diff_energies < peak_1.fwhm and diff_energies < peak_2.fwhm, np.amin([diff_energies - peak_1.fwhm, diff_energies - peak_2.fwhm]))# check with difference, but probably it's better to have the ratio
                discriminant = diff_energies / np.amin([peak_1.fwhm, peak_2.fwhm])
                comparisons["%s--%s" % (peak_1.label, peak_2.label)] = ((peak_1, peak_2), discriminant < 1, discriminant)
        # check if the couples which can be merged haven't got any common final state: in case, take the case for which the third term of the tuple is minimum.
        lines_to_merge = {} # dict containing a label of the peaks to be merged 'label_1--label_2' a tuple with the actual peaks
        # order the dictionary according to the differences between the signals peaks and each FWHM:
        # in this way we start our analysis from the very minimum and we can drop the other final states which may be merged because this value is certainly higher.
        comparisons = collections.OrderedDict(sorted(list(comparisons.items()), key = lambda item: item[1][2]))
        for comp_lines, (peaks_to_sum, merging_condition, discriminant) in comparisons.items():
            if not merging_condition or any(l in lm for l in comp_lines.split('--') for lm in lines_to_merge.keys()):
                continue
            lines_to_merge[comp_lines] = peaks_to_sum
        # the lines_to_merge list contains the labels of the lines to be merged.
        new_spectrum = [peaks_to_sum[0] + peaks_to_sum[1] for peaks_to_sum in lines_to_merge.values()]
        for peak in self.spectrum:
            # update the final spectrum with the missing lines which have not been merged
            if any(peak.label in lm for lm in lines_to_merge.keys()):
                continue
            new_spectrum.append(peak)
        if len(lines_to_merge) == 0:
            # nothing has been merged, so check if the peaks are inside one FWHM from the detection range, otherwise drop them.
            self.spectrum = [peak for peak in new_spectrum if peak > self.line_exp.detection_range[0] - peak.fwhm/2. and peak < self.line_exp.detection_range[1] + peak.fwhm/2.]
            logger.debug(self)
            return self
        else:
            # something has been merged, so check if we can merge something more in the new spectrum.
            self.spectrum[:] = new_spectrum
            return self.merge_lines()

    def test_for_errors(self):
        ''' check if there are errors related to the peaks '''
        for code, error_name in self.ERROR_TEST.items():
            for peak in self.spectrum:
                if getattr(self, "error__" + error_name).__call__(peak):
                    peak.set_error(code)

    def error__not_a_line(self, peak):
        ''' check if the FWHM of each peak is compatible with the experimental one at that energy peak '''
        return peak.fwhm > self.line_exp.energy_resolution(peak.e_peak) * peak.e_peak * self.FWHM_MERGED_FRAC

    def error__not_enough_separated(self, peak):
        # check for incompatible errors (peaks with those errors should not be tested)
        if any(code in peak.error for code in ['1']):
            return False
        min_separation = self.n_fwhm_separation * peak.fwhm
        test_spectrum  = np.array([p for p in self if p != peak])
        # find neighbourhood of the peak
        get_e_peak     = lambda a_peak: a_peak.e_peak
        energy_peaks   = np.array(list(map(get_e_peak, test_spectrum)))
        condition      = np.logical_and(energy_peaks < peak.e_peak + min_separation, energy_peaks > peak.e_peak - min_separation)
        neighbourhood  = test_spectrum[condition]
        # check: if it is alone, then return False
        if len(neighbourhood) == 0:
            return False # no peaks too close
        # check the ratio flux/flux_UL of the peak with respect each other peak in the neighbourhood
        get_height = lambda a_peak: a_peak.flux/a_peak.flux_UL
        height_peaks = np.array(list(map(get_height, neighbourhood)))
        # height_peaks elements are of the order of 1e-10 -- 1e-30, so they can be considered zero if less than 1e-200 let's say
        if any(height_peaks < 1e-200): # happens if one peak has flux_UL == __infty__ (it is out of range), in this case can't do comparison, so assume they are too close
            logger.warning("One peak has upper limit equal to __infty__, I can't compare its height with others, I assume they can't be neglected.")
            return True # because that means peaks are too close to each other, but none of them is clearly higher
        this_peak_height = peak.flux/peak.flux_UL
        # what we do is to compare the ratios between the flux and the upper limit on the flux
        if all(height_peaks < this_peak_height/self.peak_height_factor):
            return False # in case the peak is clearly higher than the others, then no error because it can be considered separated!
        return True


# PDG-particle map
class PDGParticleMap(dict):
    ''' builds a map {pdg : particle label} and allows lookup of pdg/particle '''
    def __init__(self, *args, **kwargs):
        super(PDGParticleMap, self).__init__(*args)
        model = kwargs.get('model', None)
        if model:
            self.set_model_map(model)

    def set_model_map(self, model):
        for pdg_code, particle in six.iteritems(model.get('particle_dict')):
            self[str(pdg_code)] = particle.get_name() 

    def get_pdg(self, particle):
        ''' if 'particle' is a pdg_code itself inside the self dictionary, then return it,
            otherwise obtain the pdg_code from the particle label.
        '''
        if particle in list(self.keys()):
            return int(particle)
        this_pdg = [pdg for pdg, name in six.iteritems(self) if name == particle]
        if len(this_pdg) == 0:
            # check if it is a multiparticle
            raise ValueError("No particle '%s' in the model." % particle)
        elif len(this_pdg) > 1:
            logger.warning("Particle '%s' has multiple pdg codes in this model, please check the model." % particle)
        return int(this_pdg[0])

    def __getitem__(self, pdg):
        try:
            return super(PDGParticleMap, self).__getitem__(str(pdg))
        except KeyError:
            logger.error("No PDG code '%d' in the model, please check the model." % pdg)
            return '?'

    def find_particles(self, string):
        ''' find all particles pdg of the form (pdg) in a string and return a list of them '''
        matches = re.findall(r"-?\d+", string)
        return matches

    def format_particles(self, string):
        ''' given a string with pdg codes, substitutes them with the particle labels
            finds pdg codes with a regex expression, builds a set of unique codes, substitutes
        '''
        pdg_list = set(self.find_particles(string))
        for pdg in pdg_list:
            string = string.replace(str(pdg), self[pdg])
        return string

    def format_process(self, string):
        ''' From the string written with pdgs and symbols to a nice string with particles:
            52.52>6.-6,(-6>-24.-5,(-24>11.-12)),(6>2.-1.5)
            becomes: n1 n1 > t t~, (t~ > w- b~, (w- > e- ve~)), (t > u d~ b)
        '''
        return re.sub(r"-?\d+", lambda s: self[s.group(0)], string).replace(">", " > ").replace(".", " ").replace(",", ", ")

    def format_process_readable(self, string):
        ''' From the string written with pdgs and symbols to a raw readable string:
            52.52>6.-6,(-6>-24.-5,(-24>11.-12)),(6>2.-1.5)
            becomes: n1n1_ttx_tx_wm_bx_w_em_vex_t_udxb
        '''
        return re.sub(r"-?\d+", lambda s: self[s.group(0)], string)\
            .replace(">", "_")\
            .replace(",", "_")\
            .replace(".", "")\
            .replace("(", "")\
            .replace(")", "")\
            .replace("+", "p")\
            .replace("-", "m")\
            .replace("~", "x")

class StrProcess(object):
    ''' Class to handle processes names
        In principle it could be used as key is a dictionary, becuase it has
        both the __hash__ and __eq__ methods implemented.
    '''
    def __init__(self, string_human_friendly, processes_names_map):
        self.processes_names_map = processes_names_map
        self.human_friendly = string_human_friendly
        self._particles_pdg, self.initial_states, self.final_states = self._process()

    def __hash__(self):
        return self.human_friendly

    def __str__(self):
        return str(self._particles_pdg)

    def __eq__(self, other):
        if not isinstance(other, StrProcess):
            raise NotImplementedError
        return self.human_friendly == other.human_friendly

    def _process(self):
        particles_pdg = {}
        initial_states = []
        final_states = []
        try:
            parsing_friendly = self.processes_names_map[self.human_friendly]
        except KeyError:
            logger.warning("The process %s does not have a descrition with pdg final states" % self.human_friendly)
            return particles_pdg, initial_states, final_states
        particles_pdg.update(self._extract_particles(parsing_friendly))
        initial_states.extend(particles_pdg["initial_states"])
        final_states.extend(self._extract_final_states(particles_pdg))
        return particles_pdg, initial_states, final_states

    def _match_nested_parentheses(self, text):
        ''' To be substituted with a recursive regex from Python 3 on '''
        stack = []
        matches = []
        for i, char in enumerate(text):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    start = stack.pop()
                    matches.append(text[start:i+1])
        return matches

    def _extract_particles(self, parsing_friendly_string):
        process_pdg = {}
        hard_process = parsing_friendly_string.split(",", 1)[0]
        # decay_chains = re.findall(r"\((?:[^()]|(?R))+\)", hard_process)
        decay_chains = self._match_nested_parentheses(hard_process)
        ini_s, fin_s = map(lambda part_s: re.findall(r"-?\d+", part_s), hard_process.split(">"))
        process_pdg["initial_states"] = ini_s
        process_pdg["final_states"] = fin_s
        process_pdg["decay_chains"] = list(map(
            lambda item: self._extract_particles(item[1:-1]),
            decay_chains
        ))
        return process_pdg

    def in_initial_states(self, pdg):
        return str(pdg) in self.initial_states

    def in_final_states(self, pdg):
        return str(pdg) in self.final_states

    def _extract_final_states(self, processes_list):
        final_states = processes_list["final_states"]
        if processes_list["decay_chains"]:
            to_substitute = []
            for decay_chain in processes_list["decay_chains"]:
                ini_s = decay_chain["initial_states"][0]
                to_substitute.append(final_states.count(ini_s))
                final_states = list(filter(lambda item: item != ini_s, final_states))
            # Do 2 different loop because in one we simply collect the number of occurrences
            # of a certain particle, removing it from the final state, while in the other
            # we add the new final states
            # If we would do only one loop, then it could happen that in the first loop we
            # add new final states that could by coincidence be particles to substitute in the
            # next decay chains, which is not what we want to obtain.
            for n_to_sub, decay_chain in zip(to_substitute, processes_list["decay_chains"]):
                new_final_states = self._extract_final_states(processes_list["decay_chains"])
                final_states.extend(n_to_sub * new_final_states)
        return final_states

    #


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

def secure_float_f77(d):
    try:
        return float(d)
    except ValueError:
        m=re.search(r'''([+-]?[\d.]*)([+-]\d*)''', d)
        if m:
            return float(m.group(1))*10**(float(m.group(2)))
        raise


#===============================================================================
# CommonRunCmd
#===============================================================================
class MADDMRunCmd(cmd.CmdShell):
    
    intro_banner=\
  "            ====================================================\n"+\
  "            |                  "+bcolors.OKBLUE+"  MadDM v3.2                     "+bcolors.ENDC+"|\n"\
  "            ====================================================\n"+\
  "                                                                               \n"+\
  "                #########                                                        \n"+\
  "             ###\\\\####//#####              Launchpad:  launchpad.net/maddm      \n"+\
  "           ######\\\\##//########                                              \n"+\
  "          ########\\\\//###########                     "+bcolors.FAIL+"arXiv:1308.4955        \n"+bcolors.ENDC+\
  "         #########//\\\\############                    "+bcolors.FAIL+"arXiv:1505.04190       \n"+bcolors.ENDC+\
  "        #########//##\\\\############                   "+bcolors.FAIL+"arXiv:1804.00044        \n"+bcolors.ENDC+\
  "       ########//#####\\\\###########                   "+bcolors.FAIL+"arXiv:2107.04598        \n"+bcolors.ENDC+\
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
  "            Need to learn? -> type tutorial                               \n"+\
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
        self.indirect_directories = dict(list(zip(['Indirect_tree_cont', 'Indirect_tree_line', 'Indirect_LI_cont', 'Indirect_LI_line'], list(map(os.path.isdir, [pjoin(dir_path, 'Indirect_tree_cont'), pjoin(dir_path, 'Indirect_tree_line'), pjoin(dir_path, 'Indirect_LI_cont'), pjoin(dir_path, 'Indirect_LI_line')])))))
        self.indirect_directories_cross_section_contribution = dict(list(zip(['Indirect_tree_cont', 'Indirect_tree_line', 'Indirect_LI_cont', 'Indirect_LI_line'], [0., 0., 0., 0.])))
        self.param_card_iterator = [] #a placeholder containing a generator of paramcard for scanning
        
        # Define self.proc_characteristics (set of information related to the
        # directory status
        self.get_characteristics()

        self._two2twoLO = False # is set in ask_run_configuration
        self._dNdE_setup = False #If true, this flag allows the code to skip loading and caldulating dNdE
        #self._fit_parameters= []

        self._run_couplings = False

        #A shared param card object to the interface.
        #All parts of the code should read this card.
        #Set at the beginning of launch()
        self.param_card = None

        self.multinest_running = False
        self.auto_width = set() # keep track of width set on Auto

        self.pdg_particle_map = PDGParticleMap(self.proc_characteristics['pdg_particle_map'].items())
        self.str_processes = {}

        self.limits = ExpConstraints()
        self.vave_indirect_cont_range = (3e-6, 1.4e-4)
        self.vave_indirect_line_range = (5e-4, 1e-3) # about 180--300 km/s
        self.options = options

        self.Spectra = Spectra()
        self.Fermi   = Fermi_bounds()
        self.MadDM_version = '3.2'

        self.processes_names_map = self.proc_characteristics['processes_names_map']

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
                raise Exception("Scans are not allowed in web mode")
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
    the cards afterwards, you can type \"compute_widths\".''')

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
                self.mg5  = master_interface.MasterCmd()

        if not commands:
            return

        if not model:
            model = self.proc_characteristics['model']
        if  not self.mg5._curr_model or\
                     self.mg5._curr_model.get('modelpath+restriction') != model:
            self.mg5.do_import('model %s ' % model)

        
        for line in commands:
            self.mg5.exec_cmd(line, errorhandling=False, printcmd=False, 
                              precmd=False, postcmd=False, child=False)
            
    ############################################################################
    def do_launch(self, line):
        """run the code"""

        args = line.split()
        if '-f' in args or '--force' in args:
            force = True
        else:
            force = False
        

        # determine run_name for name in the output directory:
        if '-n' in args:
            self.run_name = args[args.index('-n')+1]
            if os.path.exists(pjoin(self.dir_path, 'output', self.run_name)):
                shutil.rmtree(pjoin(self.dir_path, 'output', self.run_name))
                for indirect_directory in self.indirect_directories.keys():
                    try:
                        shutil.rmtree(pjoin(self.dir_path, indirect_directory, 'Events', self.run_name))
                    except Exception:
                        pass
        else:
            i = 1
            while os.path.exists(pjoin(self.dir_path, 'output', 'run_%02d' %i)) or\
                  os.path.exists(pjoin(self.dir_path, 'output', 'run_%02d_01' %i)):
                i += 1
            self.run_name = 'run_%02d' % i

        self.ask_run_configuration(mode=[], force=force)

        if not self.multinest_running:
            self.compile()

        if self.param_card_iterator:
            self.run_name += '_01'

        if self.in_scan_mode or self.param_card_iterator:
            logger.info("Running %s ..." % self.run_name, "$MG:BOLD")

        # create output directory.
        os.mkdir(pjoin(self.dir_path, 'output', self.run_name))
        

        output = pjoin(self.dir_path, 'output', self.run_name, 'maddm.out') 
        misc.call(['./maddm.x', pjoin('output', self.run_name, 'maddm.out')], cwd =self.dir_path)
        #Here we read out the results which the FORTRAN module dumped into a file
        #called 'maddm.out'. The format is such that the first line is always relic density
        # , second line is the nucleon scattering cross section (SI) for proton, third is SI
        # nucleon cross section for the neutron etc. If a quantity was not calculated, we output -1

        # it computes also the percentages for each annihilation channel for relic density
        # %_<process>

        # Define a dictionary holding the results
        result = {}
        # result['GeV2pb*pb2cm2']   = GeV2pb*pb2cm2 # conversion factor                                                                                                           

        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        result['tot_SM_xsec'] = -1

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
                            sigv_temp = secure_float_f77(splitline[1])
                            oname = splitline[0].split(':',1)[1].split('_')
                            oname = oname[0]+'_'+oname[1] #To eliminate the annoying suffix coming from the '/' notation
                            if oname not in self.proc_characteristics['indirect_detection_asked_processes']:
                                continue
                            str_proc = StrProcess(oname, self.processes_names_map)
                            self.str_processes[oname] = str_proc
                            result['taacsID#'     + oname] = sigv_temp 
                            result['err_taacsID#' + oname] = 0 
                            result['lim_taacsID#' + oname] = self.limits.ID_max(mdm, ".".join(str_proc.final_states))
                    if '%' in line:
                        oname = splitline[0].split(':',1)[1].split('_')
                        oname = oname[0] + '_' + oname[1]
                        str_proc = StrProcess(oname, self.processes_names_map)
                        self.str_processes[oname] = str_proc
                        result["%%_relic_%s" % oname] = secure_float_f77(splitline[1])
                    else:
                        result[splitline[0].split(':')[0]] = secure_float_f77(splitline[1])
                            
        np_names = ['g','nue','numu','nutau']
        result['sigmav(xf)'] *= GeV2pb*pb2cm3

        if str(self.mode['indirect']).startswith('flux'):
            for chan in np_names + ['gammas','neutrinos_e', 'neutrinos_mu' , 'neutrinos_tau']: # set -1 to the possible cases
                result['flux_%s' % chan] = -1.0

        # Calculating the xsi factor. Set == 1 if relic is not evaluated
        if self.mode['relic']:
           if result['Omegah^2'] < 0:
               result['xsi'] = 1.0   
           elif result['Omegah^2'] < self.limits._oh2_planck and result['Omegah^2'] > 0:
               result['xsi'] = result['Omegah^2'] / self.limits._oh2_planck
           else: result['xsi'] = 1.0

        else: result['xsi'] = 1.0

        if self.mode['direct']:
            result['sigmaN_SI_n']    *= GeV2pb*pb2cm2
            result['sigmaN_SI_p']    *= GeV2pb*pb2cm2
            result['sigmaN_SD_p']    *= GeV2pb*pb2cm2
            result['sigmaN_SD_n']    *= GeV2pb*pb2cm2
            result['lim_sigmaN_SI_n'] = self.limits.SI_max(mdm)
            result['lim_sigmaN_SI_p'] = self.limits.SI_max(mdm)
            result['lim_sigmaN_SD_p'] = self.limits.SD_max(mdm, 'p')
            result['lim_sigmaN_SD_n'] = self.limits.SD_max(mdm, 'n')

        self.last_results = result
        self.last_results['run'] = self.run_name

                          
#        if self.mode['indirect'] and not self._two2twoLO:
#            with misc.MuteLogger(names=['madevent','madgraph'],levels=[50,50]):
#                self.launch_indirect(force)

        # initialize important variables at zero
        ## cross section contributions from different directories
        self.indirect_directories_cross_section_contribution = dict(zip(['Indirect_tree_cont', 'Indirect_tree_line', 'Indirect_LI_cont', 'Indirect_LI_line'], [0., 0., 0., 0.]))
        ## Spectra object
        self.Spectra.initialize_spectra()

        if self.mode['indirect']:
            for directory in [d for d, v in self.indirect_directories.items() if 'cont' in d and v]:
                self.launch_indirect(force, directory, self.maddm_card['vave_indirect_cont'])
            # compute total cross section only for continuum final states
            self.last_results['taacsID'] = 0.
            for process, sigmav in six.iteritems({k.replace('taacsID#', ''): v for k, v in six.iteritems(self.last_results) if k.startswith('taacsID#')}):
                if self.is_spectral_finalstate(self.str_processes[process]):
                    continue
                self.last_results['taacsID'] += sigmav
            self.launch_indirect_computations(mdm)

        if self.mode['spectral']:
            for directory in [d for d, v in self.indirect_directories.items() if 'line' in d and v]:
                self.launch_indirect(force, directory, self.maddm_card['vave_indirect_line'])
            self.launch_spectral_computations(mdm)

        if not self.in_scan_mode and not self.multinest_running:
            self.print_results()

        # Saving the output for single point or scan
        self.save_remove_output(scan = self.in_scan_mode or self.param_card_iterator)

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
        if self.param_card_iterator:
            
            param_card_iterator = self.param_card_iterator

            parameters, values =  param_card_iterator.param_order , param_card_iterator.itertag
            self.param_card_iterator = []

            ## this is to remove or save spectra, not the scan summary file!
            #self.save_remove_output(scan = True)

            # *** Initialize a list containing the desired variables in the summary output of the scan
            order = ['run']

            # *** saving the name of the iterated parameters, and the values in the results dictionary
            for par,val in zip(parameters, values):
                order.append(par)
                self.last_results[par] = val

            # *** Relic density
            if self.mode['relic']:
                order += ['Omegah^2','x_f', 'sigmav(xf)']
                for percent in [k for k in self.last_results if k.startswith('%_relic')]:
                    order += [percent]
            order.append('xsi')

            # *** Direct Detection
            if self.mode['direct'] :
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
            if self.mode['indirect'] or self.mode['spectral']:
                halo_vel_cont = self.maddm_card['vave_indirect_cont']
                halo_vel_line = self.maddm_card['vave_indirect_line']
             
                #if not self._two2twoLO:
                detailled_keys = [k for k in self.last_results if k.startswith('taacsID#') ]
                # halo velocity condition
                fermi_dsph_vel = halo_vel_cont > self.vave_indirect_cont_range[0] and halo_vel_cont < self.vave_indirect_cont_range[1] # range of validity of Fermi limits
                line_gc_vel = halo_vel_line > self.vave_indirect_line_range[0] and halo_vel_line < self.vave_indirect_line_range[1] # range of validity of line limits
                for key in detailled_keys:
                    order +=[key]

                    # add the lim key only if it has been computed
                    # if fermi_dsph_vel and line_gc_vel:
                    #     order +=['lim_'+key]
                    # elif fermi_dsph_vel and not line_gc_vel:
                    #     if self.is_spectral_finalstate(key):
                    #         continue
                    #     order +=['lim_'+key]
                    # elif not fermi_dsph_vel and line_gc_vel:
                    #     if not self.is_spectral_finalstate(key):
                    #         continue
                    #     order +=['lim_'+key]

                    # or give line limits anyway
                    if fermi_dsph_vel or self.is_spectral_finalstate(self.str_processes[key.replace("taacsID#")]):
                        order +=['lim_'+key]

                if self.mode['indirect']:
                    order.append('taacsID')
                    order.append('tot_SM_xsec')
                    order.append('Fermi_sigmav')

                if self.last_results['xsi'] >0 and self.last_results['xsi'] <1: # thermal and non thermal case 
                   order = order + ['pvalue_th','like_th','pvalue_nonth','like_nonth']
                else:
                   order = order + ['pvalue_nonth','like_nonth']

                if str(self.mode['indirect']).startswith('flux'):
                   #for channel in self.Spectra.spectra.keys(): #['neutrinos_e', 'neutrinos_mu' , 'neutrinos_tau']
                   for channel in ['gammas','neutrinos_e', 'neutrinos_mu' , 'neutrinos_tau']:
                       if 'antip' in channel or 'pos' in channel: continue
                       order.append('flux_%s' % channel)

                # add keys related to lines
                # density profile (no 'profile_name' because it is a string)
                # order += ['profile_r_s', 'profile_rho_s', 'profile_gamma', 'profile_alpha']
                # observables for each experiment
                for line_exp in self.line_experiments:
                    str_part = "line_%s_" % line_exp.get_name()
                    order.append(str_part + "Jfactor")
                    order.append(str_part + "roi")
                    for i in range(len([k for k in detailled_keys if self.is_spectral_finalstate(self.str_processes[key.replace("taacsID#")])])):
                        order.append(str_part + "peak_%d"        % (i+1))
                        order.append(str_part + "flux_%d"        % (i+1))
                        order.append(str_part + "flux_UL_%d"     % (i+1))
                        if str_part + "like_%d" % (i+1) in list(self.last_results.keys()):
                            order.append(str_part + "like_%d"        % (i+1))
                            order.append(str_part + "pvalue_%d"      % (i+1))
                        order.append(str_part + "peak_%d_error"  % (i+1))

            # remove elements which have not been computed
            order[:] = [elem for elem in order if elem in list(self.last_results.keys())]

            run_name = str(self.run_name).rsplit('_',1)[0]
            summary_file = pjoin(self.dir_path, 'output','scan_%s.txt' % run_name)
            line_spectrum_logging_file = pjoin(self.dir_path, 'output','line_spectrum_scan_%s.log' % run_name)

            # scan summary for observables
            self.write_scan_output(out_path = summary_file , keys = order, header = True )
            self.write_scan_output(out_path = summary_file , keys = order )

            with misc.TMP_variable(self, 'in_scan_mode', True):
                with misc.MuteLogger(names=['cmdprint','madevent','madgraph','madgraph.plugin'],levels=[50,50,50,20]):

                    for i,card in enumerate(param_card_iterator):
                        card.write(pjoin(self.dir_path,'Cards','param_card.dat'))
                        self.exec_cmd("launch -f -n %s_%02d" % (run_name, i+2),
                                       precmd=True, postcmd=True, errorhandling=False)
                        for par,val in zip(param_card_iterator.param_order, param_card_iterator.itertag):
                            self.last_results[par] = val
                        self.write_scan_output(out_path = summary_file , keys = order, header = False)

                        ### the following three lines are added by chiara to check the widht = auto function 
                        # self._param_card = param_card_mod.ParamCard('/Users/arina/Documents/physics/software/maddm_dev2/test_width/Cards/param_card.dat')
                        # width = self.param_card.get_value('width', 5000000)
                        # logger.warning('--> try again WY0: %.2e' % width)
                        #<=-------------- Mihailo commented out max_col = 10
                        #logger.info('Results for the point \n' + param_card_iterator.write_summary(None, order, lastline=True,nbcol=10)[:-1])#, max_col=10)[:-1])
                        
                        #self.save_remove_output(scan = True)




            param_card_iterator.write(pjoin(self.dir_path,'Cards','param_card.dat'))

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

    def is_spectral_finalstate(self, str_proc):
        return str_proc.in_final_states("22")

    def launch_indirect(self, force, indirect_directory, halo_vel):
        """running the indirect detection"""

        # if not os.path.exists(pjoin(self.dir_path, 'Indirect')):
        #     self._two2twoLO = True
            #return
        # if self.maddm_card['sigmav_method'] == 'inclusive':
        #     self._two2twoLO = True
            #return 
        
        if not self.in_scan_mode: 
            logger.info('Running indirect detection \'%s\'' % indirect_directory)
        
        fast_mode = self._two2twoLO and 'tree' in indirect_directory

        if not fast_mode:   
            if not hasattr(self, 'me_cmd'):
                self.me_cmd = Indirect_Cmd(pjoin(self.dir_path, indirect_directory), force_run=True)
                #force_run = True means no crash associated with RunWeb -> we check this later
            elif self.me_cmd.me_dir != pjoin(self.dir_path, indirect_directory):
                self.me_cmd.do_quit('') # put an argument otherwise it crashes
                self.me_cmd = Indirect_Cmd(pjoin(self.dir_path, indirect_directory), force_run=True)
                #force_run = True means no crash associated with RunWeb -> we check this later
                
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        if not fast_mode: #self.maddm_card['sigmav_method'] != 'inclusive':
            runcardpath = pjoin(self.dir_path, indirect_directory, 'Cards', 'run_card.dat')
            run_card = banner_mod.RunCard(runcardpath)

            vave_temp = math.sqrt(3)/2.0 * halo_vel

            run_card['ebeam1'] = mdm * math.sqrt(1+vave_temp**2)
            run_card['ebeam2'] = mdm * math.sqrt(1+vave_temp**2)
            run_card['use_syst'] = False
            run_card.remove_all_cut()
        

            if os.path.exists(pjoin(self.dir_path, 'Cards', 'pythia8_card.dat')):
                py8card = Indirect_PY8Card(pjoin(self.dir_path, 'Cards', 'pythia8_card.dat'))
                if py8card['Main:NumberOfEvents'] != -1:
                    run_card['nevents'] = py8card['Main:NumberOfEvents']
                    run_card['iseed'] = int(np.random.uniform(0, 10000))
                run_card.write(runcardpath)
        
            if __debug__:
                set_level = 10
            else:
                set_level = 50
            logger.info("Computing sigmav with method: %s" % self.maddm_card['sigmav_method'])
            with self.me_cmd.RunWebHandling(pjoin(self.dir_path, indirect_directory)):
                with misc.TMP_variable(banner_mod.RunCard, 'get_idbmup', lambda *args: 52):
                    #this with statement ensure that only one process is running within
                    #that directory.
                    with misc.MuteLogger(['madgraph','madevent','cmdprint'], [set_level]*3):
                        #mute logger          
                        if self.maddm_card['sigmav_method'] == 'madevent':
                            if os.path.exists(pjoin(self.dir_path,indirect_directory,'Cards','reweight_card.dat')):
                                os.remove(pjoin(self.dir_path,indirect_directory,'Cards','reweight_card.dat'))
                            self.me_cmd.do_launch('%s -f' % self.run_name)
                        elif self.maddm_card['sigmav_method'] == 'reshuffling':
                            cmd = ['launch %s' % self.run_name,
                               'reweight=indirect',
                               'edit reweight --replace_line="change velocity" --before_line="launch" change velocity %s' % vave_temp]
                            self.me_cmd.import_command_file(cmd)
            
            fermi_dsph_vel = self.maddm_card['vave_indirect_cont'] > self.vave_indirect_cont_range[0] and self.maddm_card['vave_indirect_cont'] < self.vave_indirect_cont_range[1] # range of validity of Fermi limits
            line_gc_vel = self.maddm_card['vave_indirect_line'] > self.vave_indirect_line_range[0] and self.maddm_card['vave_indirect_line'] < self.vave_indirect_line_range[1] # range of validity of line limits
            velocity_in_range = {'cont': fermi_dsph_vel, 'line': line_gc_vel}.get(indirect_directory.split('_')[-1], True)
            
            for key, value in six.iteritems(self.me_cmd.Presults):
                clean_key_list = key.split("/")
                # clean_key = clean_key_list[-1].split('_')[1] +'_'+  clean_key_list[-1].split('_')[2]
                clean_key = clean_key_list[-1].split("_", 1)[-1]
                if key.startswith('xsec'):
                    value = halo_vel * math.sqrt(3)/2 * 2 * value
                    self.last_results['taacsID#%s' %(clean_key)] = value* pb2cm3
                    self.indirect_directories_cross_section_contribution[indirect_directory] += self.last_results['taacsID#%s' %(clean_key)]
                    str_proc = StrProcess(clean_key, self.processes_names_map)
                    self.str_processes[clean_key] = str_proc
                    self.last_results['lim_taacsID#'+clean_key] = self.limits.ID_max(mdm, ".".join(str_proc.final_states)) if velocity_in_range else -1
                    
                elif key.startswith('xerr'):
                    self.last_results['err_taacsID#%s' %(clean_key)] = value * pb2cm3

        # the following is done above
        # *** Multiply all the calculated indirect cross section by sqrt(3)/2 * 2 *(vave_indirect) value (since is relative velocity)
        # for key,value in (self.last_results).iteritems():
        #     if 'taacs' in key and 'lim' not in key and 'err' not in key or 'xsec' in key:
        #         new_value = halo_vel * math.sqrt(3)/2 * 2 * value
        #         self.last_results[key] = new_value

    def launch_indirect_computations(self, mdm):
        ''' do indirect detection calculations for continuum spectra '''
        # check that whether we are in fast mode at least for one of the directories computed
        cont_spectra_directories = [d for d, v in self.indirect_directories.items() if ('cont' in d) and v]
        fast_mode = any(self._two2twoLO and 'tree' in directory for directory in cont_spectra_directories)
        # compute spectra with Pythia or PPPC only in the continuum case
        if not fast_mode:
            if self.mode['indirect'] != 'sigmav':
                if self.maddm_card['indirect_flux_source_method'].startswith('PPPC') or self.maddm_card['indirect_flux_source_method'].lower().startswith('cosmixs'):
                    if self.read_tabulatedspectra():
                        logger.info(f"Calculating Fermi dSph limit using {self.maddm_card['indirect_flux_source_method']} spectra")
                    else:
                        logger.info('no PPPC4DMID')
                else:
                    for id_dir in cont_spectra_directories:
                        self.run_pythia8_for_flux(id_dir)
                    if self.maddm_card['indirect_flux_source_method'] == 'pythia8' or self.maddm_card['indirect_flux_source_method'] == 'vincia':
                        logger.info('Calculating Fermi dSph limit using pythia8 ' + '(vincia) '*(self.maddm_card['indirect_flux_source_method'] == 'vincia') + 'gamma rays spectrum')
                    elif 'pythia' not in self.maddm_card['indirect_flux_source_method'] and 'vincia' not in self.maddm_card['indirect_flux_source_method']:
                        logger.warning('Since pythia8 is run, using pythia8 gamma rays spectrum (not PPPC4DMID Tables)')
                    self.read_py8spectra(cont_spectra_directories)
            elif self.mode['indirect'] == 'sigmav':
                logger.warning('no gamma spectrum since in sigmav mode')

        elif self.mode['indirect'] == 'sigmav':
            logger.warning('no gamma spectrum since in sigmav mode')
        elif self.read_tabulatedspectra():   # if in fast mode, use PPPC. return False if PPPC4DMID not installed!
            logger.info(f"Calculating Fermi dSph limit using {self.maddm_card['indirect_flux_source_method']} spectra")

        # ****** Calculating Fermi dSph Limits
        self.calculate_fermi_limits(mdm)
        
        # ****** Calculating Fluxes at detection
        if str(self.mode['indirect']).startswith('flux'):
            self.neu_oscillations() # neutrinos oscillations                                                                                                                   

        # ****** Calculating Fluxes Earth                                                                                                                       
            if 'earth' in self.mode['indirect']:
                if 'PPPC' in self.maddm_card['indirect_flux_earth_method']:
                    self.read_PPPC_positrons_earth()
                elif 'dragon' in self.maddm_card['indirect_flux_earth_method']:
                    logger.info('Calling DRAGON for positrons and antiprotons propagation.')
                    self.run_Dragon()                                                                                                                

            self.calculate_fluxes() # calculating dPhidE  
                                                                                                                    

    def calculate_fermi_limits(self, mdm):
        """setup the computation for the fermi dSph limits"""

        # check that the vave is in the range (1-50 km/s) = (3*10^-6 - 1.5*10^-4)/c                                                                                                 
        self.last_results['Fermi_sigmav'] = -1
        self.last_results['pvalue_th']    = -1
        self.last_results['like_th']      = -1
        self.last_results['pvalue_nonth'] = -1
        self.last_results['like_nonth']   = -1

        halo_vel = self.maddm_card['vave_indirect_cont']
        
        if halo_vel < self.vave_indirect_cont_range[0] or halo_vel > self.vave_indirect_cont_range[1]:
           logger.error('The DM velocity in the dwarfs halo is not in the [3*10^-6 - 1.5*10^-4]/c range - will not calculate Fermi limits!')        
           logger.error('Please rerun with the correct velocity to re-calculate the correct sigmav for Fermi limits.')
           return 0
    
        x, gammas = self.Spectra.spectra['x'] , self.Spectra.spectra['gammas']

        if not gammas:
            if not self.options['pppc4dmid_path'] and self.maddm_card['indirect_flux_source_method'].startswith('PPP'):
                pass
            else:
                logger.error('The gamma spectrum is empty! Will not calculate Fermi limit')

        elif gammas:
            
            sigmav = self.Fermi.Fermi_sigmav_lim(mdm, x , gammas ,maj_dirac= self.norm_Majorana_Dirac() )
            self.last_results['Fermi_sigmav'] = sigmav # returns Fermi exp UL
                    
            if 'PPPC' in self.maddm_card['indirect_flux_source_method']:
                sigmav_th = self.last_results['tot_SM_xsec']
            else :
                sigmav_th = self.last_results['taacsID']

            pvalue_nonth , like_nonth = self.Fermi.Fermi_sigmav_lim(mdm, x , gammas ,maj_dirac= self.norm_Majorana_Dirac() , sigmav_th = sigmav_th )
            self.last_results['pvalue_nonth'] = pvalue_nonth 
            self.last_results['like_nonth']   = like_nonth
            logger.debug('sigmav = %s , pvalue-nonth = %s , like-nonth = %s '  %(sigmav , pvalue_nonth , like_nonth ) )

            # Evaluates also the thermal scenarios likelihood and p-value when 0 < xsi < 1
            if self.last_results['xsi'] > 0 and self.last_results['xsi'] < 1:                                                                                  
               sigmav_th = sigmav_th * self.last_results['xsi']**2
               pvalue_th , like_th  = self.Fermi.Fermi_sigmav_lim(mdm, x , gammas ,maj_dirac= self.norm_Majorana_Dirac(), sigmav_th = sigmav_th )
               self.last_results['pvalue_th'] = pvalue_th
               self.last_results['like_th']   = like_th
            
               logger.debug('sigmav = %s , p-th=%s , like-th=%s , p-nonth=%s , like-nonth=%s '  %(sigmav, pvalue_th , like_th , pvalue_nonth , like_nonth ) )


    def launch_spectral_computations(self, mdm):
        """running the indirect detection for spectral features"""
        profiles = {
            'gnfw'      : PROFILES.NFW,
            'nfw'       : lambda **kwargs: (kwargs.pop('gamma', None), PROFILES.NFW(gamma = 1.0, **kwargs))[1], # to make canonical nfw keeping the freedom of kwargs, we need to pop "gamma" from kwargs and set it explicitly: make a lambda and built a tuple inside it: kwargs.pop modifies the dictionary and PROFILES.NFW(gamma=1.0, **kwargs) uses the modified kwargs because the tuple is instantiated and then processed, the [1] element is then returned from the lambda
            'einasto'   : PROFILES.Einasto,
            'isothermal': PROFILES.Isothermal,
            'burkert'   : PROFILES.Burkert
        }
        # set the default Fermi-LAT 2015 ROI for each profile (for HESS only one ROI)
        fermilat_2015_rois_default = {
            'gnfw'      : 3.,
            'nfw'       : 41.,
            'einasto'   : 16.,
            'isothermal': 90.,
            'burkert'   : 90. 
        }
        if self.maddm_card['roi_fermi_2015'] == 'default':
            roi_fermi_2015 = fermilat_2015_rois_default[self.maddm_card['profile']]
        else:
            roi_fermi_2015 = int(self.maddm_card['roi_fermi_2015'].strip('rR'))

        ### Fill list with GammaLineExperiment
        if not self.in_scan_mode:
            self.line_experiments = GammaLineExperimentsList()
            r_sun = self.maddm_card['r_sun']
            rho_sun = self.maddm_card['rho_sun']
            ### Fermi-LAT 2015
            log_masses = [-0.939, -0.818, -0.69, -0.562, -0.438, -0.313, -0.188, -0.0638, 0.0608, 0.188, 0.313, 0.438, 0.562, 0.687, 0.812, 0.936, 1.06, 1.19, 1.31, 1.44, 1.56, 1.69, 1.81, 1.94, 2.06, 2.19, 2.31, 2.44, 2.56, 2.69]
            e_resolution = [0.14, 0.122, 0.109, 0.0929, 0.0852, 0.0748, 0.0686, 0.0629, 0.059, 0.0581, 0.0543, 0.0505, 0.0457, 0.0424, 0.04, 0.0381, 0.0371, 0.0352, 0.0343, 0.0333, 0.0324, 0.0319, 0.0319, 0.0324, 0.0324, 0.0329, 0.0338, 0.0348, 0.0348, 0.0343]
            if HAS_SCIPY: # so that it does nothing in case scipy is missing: later the program would not use this function
                resolution_interp_func = interp1d(log_masses, e_resolution, kind = 'linear')
                def energy_resolution_fermi_line_2015(m):
                    if np.log10(m) <= -0.939:
                        return resolution_interp_func(-0.939)
                    elif np.log10(m) >= 2.69:
                        return resolution_interp_func(2.69)
                    else:
                        return resolution_interp_func(np.log10(m))
            else:
                def energy_resolution_fermi_line_2015(m):
                    return -1.

            self.line_experiments.append(
                experiment = GammaLineExperiment(
                    name                  = "Fermi-LAT_2015",
                    energy_resolution     = energy_resolution_fermi_line_2015,
                    detection_range       = [0.214, 462.],
                    roi_list              = line_experiments_regions_of_interest["Fermi-LAT_2015"],
                    majorana_dirac_factor = self.norm_Majorana_Dirac() / 4., # divide by 4, because that methods return 4 if Majorana, 8 if Dirac
                    check_profile_message = lambda is_optimized, roi: {True: "", False: "ROI %d is not optimized for this profile!" % roi}.get(is_optimized),
                    arxiv_number          = "1506.00013",
                    hidden                = (self.maddm_card['toggle_fermi_2015'] == 'off')
                ),
                roi_to_analyse = roi_fermi_2015
            )
            ### HESS 2018
            self.line_experiments.append(
                experiment = GammaLineExperiment(
                    name                  = "HESS_2018",
                    energy_resolution     = lambda m: 0.1,
                    detection_range       = [306.9, 63850.],
                    roi_list              = line_experiments_regions_of_interest["HESS_2018"],
                    majorana_dirac_factor = self.norm_Majorana_Dirac() / 4., # divide by 4, because that methods return 4 if Majorana, 8 if Dirac
                    check_profile_message = lambda is_optimized, roi: {True: "", False: "The chosen profile is not the default for this ROI"}.get(is_optimized),
                    arxiv_number          = "1805.05741",
                    hidden                = (self.maddm_card['toggle_hess_2018'] == 'off')
                ),
                roi_to_analyse = 1.
            )
            ### TEMPLATE EXPERIMENT
            template_roi = self.maddm_card["template_line_experiment_roi"]
            template_profile = profiles[self.maddm_card['template_line_experiment_profile']](
                r_s     = self.maddm_card['template_line_experiment_r_s'],
                gamma   = self.maddm_card['template_line_experiment_gamma'],
                alpha   = self.maddm_card['template_line_experiment_alpha'],
                rho_sun = rho_sun,
                r_sun   = r_sun
            )
            self.line_experiments.append(
                experiment = GammaLineExperiment(
                    name                  = self.maddm_card["template_line_experiment_name"],
                    energy_resolution     = lambda m: self.maddm_card["template_line_experiment_energy_resolution"],
                    detection_range       = [self.maddm_card["template_line_experiment_detection_range_min"], self.maddm_card["template_line_experiment_detection_range_max"]],
                    roi_list              = RegionOfInterestList([
                        RegionOfInterest(
                            amplitude=template_roi,
                            default_profile=template_profile,
                            upper_limit_label="(22)(22)_template",
                            jfactor_strategy=self.maddm_card["template_line_experiment_jfactor_strategy"],
                            mask_latitude = self.maddm_card["template_line_experiment_mask_latitude"], # deg
                            mask_longitude = self.maddm_card["template_line_experiment_mask_longitude"], # deg
                            mask_inner_amplitude = self.maddm_card["template_line_experiment_mask_inner_angle"], # deg
                            cached_jfactors={},
                            r_max=np.inf if self.maddm_card['template_line_experiment_r_max'] == "inf" else self.maddm_card['template_line_experiment_r_max'],
                        ),
                    ]),
                    majorana_dirac_factor = self.norm_Majorana_Dirac() / 4., # divide by 4, because that methods return 4 if Majorana, 8 if Dirac
                    check_profile_message = lambda is_optimized, roi: {True: "", False: "The chosen profile is not the default for this ROI"}.get(is_optimized),
                    arxiv_number          = self.maddm_card["template_line_experiment_arxiv"],
                    hidden                = (self.maddm_card['toggle_template_line_experiment'] == 'off')
                ),
                roi_to_analyse = template_roi
            )
            if (self.maddm_card['template_line_experiment_constraints_file'] != 'None') and (self.maddm_card['toggle_template_line_experiment'] != 'off'):
                self.limits._id_limit_file['22.22_template'] = pjoin(MDMDIR, 'ExpData', self.maddm_card['template_line_experiment_constraints_file'])
                self.limits._id_limit_file_flux['22.22_template'] = pjoin(MDMDIR, 'ExpData', self.maddm_card['template_line_experiment_constraints_file'])
                self.limits._id_limit_vel['22.22_template'] = self.vave_indirect_line_range[-1]
                logger.info("Reload constraints to include new gamma ray line experiment...")
                self.limits.load_constraints()
              
        # find profile
        density_profile = profiles[self.maddm_card['profile']](
            r_s     = self.maddm_card['r_s'],
            gamma   = self.maddm_card['gamma'],
            alpha   = self.maddm_card['alpha'],
            rho_sun = self.maddm_card["rho_sun"],
            r_sun   = self.maddm_card['r_sun']
        )
        # fill last_results with the profile parameters
        self.last_results["density_profile"] = density_profile
        # check scipy
        if HAS_SCIPY is False:
            logger.warning("scipy module not available, line limits computation is disabled.")
            return
        # check halo velocity
        halo_vel = self.maddm_card['vave_indirect_line']
        # if not (halo_vel > self.vave_indirect_line_range[0] and halo_vel < self.vave_indirect_line_range[1]):
        #     logger.warning('The DM velocity in the galactic center is not in the [%.3e - %.3e] c range - line limits application is this case is not guaranteed.' % (self.vave_indirect_line_range[0], self.vave_indirect_line_range[1]))        
        #     for line_exp in self.line_experiments:
        #         str_part = "line_%s_" % line_exp.get_name()    
        #         self.last_results[str_part + "Jfactor"] = -1
        #     return 0   
        logger.info("Calculating line limits from " + ', '.join([name.replace('_', ' ') for name in self.line_experiments.iternames()]))
        # <sigma v> of the various final states
        sigmavs = {".".join(self.str_processes[k.replace("taacsID#","")].final_states): v for k, v in six.iteritems(self.last_results) if k.startswith('taacsID#') and self.is_spectral_finalstate(self.str_processes[k.replace("taacsID#","")])}
        # dict for <sigma v> ul
        sigmav_ul = {k: collections.OrderedDict([(name, np.inf) for name in self.line_experiments.iternames()]) for k in self.last_results.keys() if "lim_taacsID#" in k and self.is_spectral_finalstate(self.str_processes[k.replace("lim_taacsID#","")])} # if there is at least one '22' then we treat it as a spectral final state
        # compute the main results
        for line_exp, line_exp_roi in zip(self.line_experiments, self.line_experiments.iterrois()):
            gamma_line_spectrum = GammaLineSpectrum(
                line_exp = line_exp, 
                pdg_particle_map = self.pdg_particle_map, 
                param_card = self.param_card, 
                n_fwhm_separation = self.maddm_card['n_fwhm_separation'], 
                peak_height_factor = self.maddm_card['peak_height_factor']
            )
            self.last_results["line_%s_roi" % line_exp.get_name()] = line_exp_roi
            str_part = "line_%s_" % line_exp.get_name()
            # likelihood object
            if line_exp.get_name() == "Fermi-LAT_2015":
                exp_likelihood = Fermi2015GammaLineLikelihood(line_exp)
            else:
                exp_likelihood = None
            # fill last_results with all the keys (with value -1, which means 'not computed') of line analysis
            for i in range(len(sigmavs)):
                self.last_results[str_part + "peak_%d"        % (i+1)] = -1
                self.last_results[str_part + "peak_%d_states" % (i+1)] = '/'
                self.last_results[str_part + "peak_%d_error"  % (i+1)] = -1 # -1 means not measured
                self.last_results[str_part + "flux_%d"        % (i+1)] = -1
                self.last_results[str_part + "flux_UL_%d"     % (i+1)] = -1
                # consider likelihood only when it is present
                if exp_likelihood:
                    self.last_results[str_part + "like_%d"        % (i+1)] = -1
                    self.last_results[str_part + "pvalue_%d"      % (i+1)] = -1
            # fill last_results dict with the new values
            self.last_results[str_part + "Jfactor"    ] = line_exp.get_J(roi = line_exp_roi, profile = density_profile)
            self.last_results[str_part + "roi_warning"] = line_exp.check_profile(roi = line_exp_roi, profile = density_profile)[1]          
            gamma_line_spectrum.find_lines(mdm, sigmavs, self.last_results[str_part + "Jfactor"])
            original_gamma_spectrum = deepcopy(gamma_line_spectrum)
            gamma_line_spectrum.merge_lines()
            gamma_line_spectrum.sort()
            # the following dicts are for the testing errors
            for i, peak in enumerate(gamma_line_spectrum):
                self.last_results[str_part + "peak_%d"        % (i+1)] = peak.e_peak
                self.last_results[str_part + "peak_%d_states" % (i+1)] = self.pdg_particle_map.format_particles(peak.label).replace(".", "")
                # flux upper limits
                peak.flux_UL = line_exp.get_flux_ul(e_peak = peak.e_peak, roi = line_exp_roi, id_constraints = self.limits)
                self.last_results[str_part + "flux_%d"    % (i+1)] = peak.flux
                self.last_results[str_part + "flux_UL_%d" % (i+1)] = peak.flux_UL
            # first: merging algorithm
            # second: check for errors after having computed the flux upper limits
            # third: if there are errors, then set flux and flux_UL to -1
            gamma_line_spectrum.test_for_errors()
            for i, peak in enumerate(gamma_line_spectrum):
                self.last_results[str_part + "peak_%d_error"  % (i+1)] = peak.get_error_num()
                # if peak.error:
                #     self.last_results[str_part + "flux_%d"    % (i+1)] = -1
                #     self.last_results[str_part + "flux_UL_%d" % (i+1)] = -1
                if exp_likelihood: # compute likelihoods only if there are no errors and if exists
                    try:
                        self.last_results[str_part + "like_%d"   % (i+1)] = exp_likelihood.loglike(peak, line_exp_roi, density_profile)
                        self.last_results[str_part + "pvalue_%d" % (i+1)] = exp_likelihood.compute_pvalue(peak, line_exp_roi, density_profile)
                    except IOError as err:
                        logger.warning(err.args[0])
                        continue
            # compute <sigma v> ul to display only if the profile chosen is equal to the default one for the ROIs
            # fill a list with the limits from each experiment and then take the minimum
            possible_final_states = { ".".join(self.str_processes[sv_k.replace("lim_taacsID#","")].final_states) : sv_k for sv_k in sigmav_ul.keys() }
            for original_peak in original_gamma_spectrum: # fill value only for peaks which have been detected (are in the detection range), otherwise keep -1
                # rescale it wirth respect to gamma gamma, according to (2)*(mdm^2/E_peak^2) : the 2 factor is comprised in the method get_sigmav_ul
                sigmav_ul[possible_final_states[original_peak.label]][line_exp.get_name()] = line_exp.get_sigmav_ul(e_peak = original_peak.e_peak, roi = line_exp_roi, profile = density_profile, id_constraints = self.limits, is_aa = original_peak.label == '22.22') * np.power(mdm/original_peak.e_peak, 2)
        # get the <sigma v> ul: drop the -1 and get the minimum, in case everything is -1 then the list is empty, so return -1
        logger.debug(sigmav_ul)
        for k, v in sigmav_ul.items():
            limits = np.array(list(v.values()))
            index_of_min = np.argmin(limits)
            strongest_limit = limits[index_of_min]
            self.last_results[k] = strongest_limit if not np.isinf(strongest_limit) else -1
            self.last_results[k + '_experiment'] = list(v.keys())[index_of_min] if not np.isinf(strongest_limit) else ''
        # in case all of the line peaks are out of range of detection for a certain experiment:
        # gamma_line_spectrum would be an empty list
        # self.last_results["line_%s_peak_%d" % (exp_name, i+1)], self.last_results["line_%s_flux_%d" % (exp_name, i+1)], self.last_results["line_%s_flux_UL_%d" % (exp_name, i+1)] will be all -1

    def run_pythia8_for_flux(self, indirect_directory):
        """ compile and run pythia8 for the flux"""
        
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        #compile the pythia8 script
        if not os.path.exists(pjoin(self.dir_path,'bin','internal','main101')):
            
            if not hasattr(self, 'mg5'):
                self.run_mg5([]) # just initialize self.mg5
            
            
            py8 = self.mg5.options['pythia8_path']
            if not py8:
                logger.critical('Indirect detection: py8 code is not installed (or linked). Skip flux calculation')
                logger.info('you can install pythia8 via the command: \'install pythia8\'')
                return
                
            files.cp(pjoin(py8, 'share','Pythia8','examples','Makefile'), 
               pjoin(self.dir_path,'bin','internal'))
            files.cp(pjoin(py8, 'share','Pythia8','examples','Makefile.inc'), 
               pjoin(self.dir_path,'bin','internal'))
            try:
                misc.compile(['main101'], cwd=pjoin(self.dir_path,'bin','internal'))
            except MadGraph5Error as e:
                logger.debug(str(e))
                logger.critical('Indirect detection, py8 script can not be compiled. Skip flux calculation')
                return

        
        # Now write the card.
        pythia_cmd_card = pjoin(self.dir_path, indirect_directory ,'Source', "spectrum.cmnd")
        # Start by reading, starting from the default one so that the 'user_set'
        # tag are correctly set.
        PY8_Card = Indirect_PY8Card(pjoin(self.dir_path, 'Cards', 'pythia8_card_default.dat'), mode_vincia=(self.maddm_card['indirect_flux_source_method'] == "vincia"))
        PY8_Card['Main:mDM'] = mdm
        PY8_Card.read(pjoin(self.dir_path, 'Cards', 'pythia8_card.dat'), setter='user')

        PY8_Card.write(pythia_cmd_card, 
                       pjoin(self.dir_path, 'Cards', 'pythia8_card_default.dat'),
                        direct_pythia_input=True)
            
        run_name = self.me_cmd.run_name
        # launch pythia8
        pythia_log = pjoin(self.dir_path, indirect_directory, 'Events', run_name, 'pythia8.log')

        # Write a bash wrapper to run the shower with custom environment variables
        wrapper_path = pjoin(self.dir_path,indirect_directory, 'Events',run_name,'run_shower.sh')
        wrapper = open(wrapper_path,'w')
        shell = 'bash' if misc.get_shell_type() in ['bash',None] else 'tcsh'
        shell_exe = None
        if os.path.exists('/usr/bin/env'):
            shell_exe = '/usr/bin/env %s'%shell
        else: 
            shell_exe = misc.which(shell)
        if not shell_exe:
            raise self.InvalidCmd('No shell could be found in your environment.\n'+
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
                                  cwd=pjoin(self.dir_path,indirect_directory,'Events',run_name))
        else:                
            ret_code = misc.call(wrapper_path, stdout=open(pythia_log,'w'), stderr=subprocess.STDOUT,
                                  cwd=pjoin(self.dir_path,indirect_directory,'Events',run_name))

        #### WORK ON CLUSTER !!!

        ### FF Fix to make it work on cluster!
        # if ret_code != 0:
        #    raise self.InvalidCmd, 'Pythia8 shower interrupted with return code %d.\n'%ret_code+ ' You can find more information in this log file:\n%s' % pythia_log


    # reading the spectra from the pythia8 output
    def read_py8spectra(self, directories):
        run_name = self.me_cmd.run_name
        temp_spectra = {}
        temp_errors = {}
        for sp in self.Spectra.spectra.keys(): 
            if sp in ['x', 'log10x']: continue
            sp_name = sp + '_spectrum_pythia8.dat'
            out_dir = lambda indirect_directory: pjoin(self.dir_path, indirect_directory, 'Events', run_name, sp_name)
            # load x values: they are the same for each file
            if sp == 'gammas': # x values are the same for all the spectra
                self.Spectra.spectra['log10x'] = np.loadtxt( out_dir(directories[0]) , unpack = True )[0]
                self.Spectra.spectra['x'] = np.power(10,self.Spectra.spectra['log10x'])     # from log[10,x] to x
            temp_spectra[sp] = []
            temp_errors[sp] = []
            for id_dir in directories: # loop over the indirect directories which have run pythia
                file_rows = np.loadtxt( out_dir(id_dir) , unpack = True )
                br_factor = self.indirect_directories_cross_section_contribution[id_dir]/self.last_results['taacsID']
                if len(file_rows) > 2: # errors are present
                    temp_errors[sp].append( np.power(file_rows[2] * br_factor, 2) )
                temp_spectra[sp].append( file_rows[1] * br_factor )
            self.Spectra.spectra[sp] = np.sum(temp_spectra[sp], axis = 0).tolist() # conversion to list is needed, because in the program we use bool(a_list) to see if a_list is not empty (this doesn't work for numpy.array objects)
            self.Spectra.errors[sp] = np.sqrt(np.sum(temp_errors[sp], axis = 0)).tolist()

    def cont_spectra_available_channels(self):
        available_channels = {}
        for x in self.last_results.keys():
            if 'err' not in x and 'taacsID#' in x and 'lim_' not in x:
                available_channels[x] = self.str_processes[x.replace("taacsID#","")]
                if self.is_spectral_finalstate(available_channels[x]): # remove spectral finalstates
                    del available_channels[x]
        return available_channels

    # This function reads the spectra from the PPPC tables from each annihilation channel, and adds them up according to the BR 
    def read_tabulatedspectra(self):
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        if not self.options['pppc4dmid_path']:
            logger.error('PPPC4DMID not installed. Please install by typing "install PPPC4DMID".')
            return

        if not HAS_SCIPY:
            logger.error('using PPPC4DMID requires scipy module. Please install it (for example with "pip install scipy")')
            return

        if self.maddm_card['indirect_flux_source_method'].lower() in ['pppc4dmid', 'pppc4dmid_ew', 'cosmixs'] or 'inclusive' in self.maddm_card['sigmav_method']:
            if self.Spectra.check_mass(mdm):
                if self.maddm_card['indirect_flux_source_method'].lower() == 'cosmixs':
                    kind = 'ew'
                elif self.maddm_card['indirect_flux_source_method'].lower() == 'pppc4dmid_ew':
                    kind = 'pppc_ew'
                else:
                    kind = 'pppc'
                PPPC_source = self.Spectra.load_tabulated_sourcespectra(self.options['pppc4dmid_path'], kind = kind)

        # Combine the spectra: multiply each channel by the BR                                                                                                  
        # I create a temporary list that at each loop in the channel, retains the partial sum of each x value multiplied by the channel BR                                      
        # At the end of the loop over the channels, this is the combined spectrum (do it for all the spectra - gammas, nue, positron etc. )

        self.Spectra.spectra['x'] = PPPC_source['x']
        self.Spectra.spectra['log10x'] = np.log10(self.Spectra.spectra['x'])

        available_channels = self.cont_spectra_available_channels()
        available_final_states_strings = list(map(lambda item: ".".join(item.final_states), available_channels.values()))
        # {'taacsID#xxdxxdb_ccx': 'ccx', 'taacsID#xxdxxdb_y0y0': 'y0y0', 'taacsID#xxdxxdb_ttx': 'ttx', 'taacsID#xxdxxdb_ssx': 'ssx', 
        # 'taacsID#xxdxxdb_uux': 'uux', 'taacsID#xxdxxdb_ddx': 'ddx', 'taacsID#xxdxxdb_bbx': 'bbx'}

        # self.Spectra.spectra_id = {'px':'antiprotons', 'gx':'gammas', 'nuex':'neutrinos_e', 'numux':'neutrinos_mu', 'nutaux':'neutrinos_tau', 'ex':'positrons'}

        # Check that at lest one SM channel is available
        if not any(i in self.Spectra.map_allowed_final_state_PPPC.keys() for i in available_final_states_strings):
            logger.error('No SM annihilation channel available, cannot use PPPC4DMID Tables!')
            return
                  
        for sp, sp_t in six.iteritems(self.Spectra.spectra_id):
            self.last_results['tot_SM_xsec'] = 0
            # self.Spectra.map_allowed_final_state_PPPC = {'qqx':'qq', 'ccx':'cc', 'gg':'gg', 'bbx':'bb', 'ttx':'tt',
            #                                             'e+e-':'ee', 'mu+mu-':'mumu', 'ta+ta-':'tautau', 'w+w-':'WW', 'zz':'ZZ', 'hh':'hh' }              
            temp_dic = {}

            for CH_k in available_channels.keys():
                CH = ".".join(available_channels[CH_k].final_states)
                if CH in self.Spectra.map_allowed_final_state_PPPC.keys():
                    ch = self.Spectra.map_allowed_final_state_PPPC[CH]  # ch is the name of the channels in the Tables
                    
                # Mapping liht quarks to qq in the Tables
                # elif CH == 'ssx' : ch = 'qq'
                # elif CH == 'uux' : ch = 'qq' 
                # elif CH == 'ddx' : ch = 'qq'
                else: continue
                ch_BR =  self.last_results[CH_k]                                                                                                           
                if ch_BR > 0:                                                                                                                                              
                    temp_dic[CH] = {}                                                                                                                                      
                    self.last_results['tot_SM_xsec'] +=  ch_BR                                                                                                          
                    interp_spec = self.Spectra.interpolate_spectra(PPPC_source , mdm = mdm, spectrum = sp_t , channel = ch , earth = False )        
                    temp_dic[CH]['spec'] = interp_spec                                                                                                                
                    temp_dic[CH]['xsec'] = ch_BR
            sum_spec = []
            if not bool(temp_dic): 
                logger.error('There is no annihilation process into SM with cross section > 0!')
                return 
            for num in range(0,len(temp_dic[list(temp_dic.keys())[0]]['spec']) ):
                val = 0
                for k in temp_dic.keys():
                    val = val + temp_dic[k]['spec'][num] * ( temp_dic[k]['xsec'] / self.last_results['tot_SM_xsec'] )
                sum_spec.append(val)

            self.Spectra.spectra[sp_t] = sum_spec
        return True

    
    def read_PPPC_positrons_earth(self):
        # This funcion reads the positron flux at earth from the PPPC Tables
         
        if not self.options['pppc4dmid_path']:
            logger.error('PPPC4DMID not installed. Please install by typing "install PPPC4DMID".')
            return 
                 
        PROF = self.maddm_card['dm_profile']
        PROP = self.maddm_card['prop_method']
        HALO = self.maddm_card['halo_funct']
        logger.debug('FF profile,prop,halo ***  ' + PROF + ' ' + PROP + ' ' + HALO ) 

        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        if 'PPPC4DMID' in self.maddm_card['indirect_flux_earth_method'] or self.maddm_card['sigmav_method'] == 'inclusive':
            if self.Spectra.check_mass(mdm):
                PPPC_earth = self.Spectra.load_PPPC_earth(self.options['pppc4dmid_path'], prof = PROF)
               
        E = 10**np.array(PPPC_earth['x']) # converting Log10(E) in E ; will save fluxes as E vs dPhi/dE
        self.Spectra.flux_earth_positrons['e'] = E

        channels = self.cont_spectra_available_channels()
        available_final_states_strings = list(map(lambda item: ".".join(item.final_states), channels.values()))

        # Check that at lest one SM channel is available
        if not any(i in self.Spectra.map_allowed_final_state_PPPC.keys() for i in available_final_states_strings):
            logger.error('No SM annihilation channel available, cannot use PPPC4DMID Tables!')
            return

        # filling a dictionary with all non -zero cross section positrons fluxes from SM annihilations channels
        temp_dic = {}
        for CH, proc in channels.items(): # e.g.  ccx, ttx ecc.
            C = ".".join(proc.final_states)
            if C in self.Spectra.map_allowed_final_state_PPPC.keys():
               ch = self.Spectra.map_allowed_final_state_PPPC[C]  # ch is the name of the channels in the Tables                                                           
            #       # Mapping light quarks to qq in the Tables                                                                                                                      
            # elif C == 'ssx' : ch = 'qq'
            # elif C == 'uux' : ch = 'qq'
            # elif C == 'ddx' : ch = 'qq'
            else: continue

            ch_BR = self.last_results[CH]
            if ch_BR > 0:
                temp_dic[C] = {}
                interp_spec = self.Spectra.interpolate_spectra(PPPC_earth, mdm = mdm, spectrum = 'positrons' , channel = ch , \
                                    earth = True , prof = PROF , prop = PROP , halo_func = HALO)
                
                temp_dic[C]['spec'] = interp_spec
                temp_dic[C]['xsec'] = ch_BR

        sum_spec = []
        for num in range(0,len(temp_dic[list(temp_dic.keys())[0]]['spec']) ):
            val = 0
            for k in temp_dic.keys():
                val = val + temp_dic[k]['spec'][num] * ( temp_dic[k]['xsec'] / self.last_results['tot_SM_xsec'] )
            sum_spec.append(val)

        self.Spectra.flux_earth_positrons['positrons']['dPhidlogE'] = sum_spec
         
    def neu_oscillations(self):
        # This function calculates neutrino osccillations and gives the fluxes at earth
        nue   = self.Spectra.spectra['neutrinos_e'  ]
        numu  = self.Spectra.spectra['neutrinos_mu' ]
        nutau = self.Spectra.spectra['neutrinos_tau']
        P_e_mu   = 0.199
        P_e_tau  = 0.221
        P_mu_tau = 0.368

        nue_earth , numu_earth , nutau_earth = [] , [] , []

        for x_e, x_mu, x_tau in zip(nue,numu,nutau):
            nue_earth  .append( (1-P_e_mu-P_e_tau)  * x_e   + (P_e_mu) * x_mu + (P_e_tau)  * x_tau )
            numu_earth .append( (1-P_e_mu-P_mu_tau) * x_mu  + (P_e_mu) * x_e  + (P_mu_tau) * x_tau )
            nutau_earth.append( (1-P_e_tau-P_mu_tau)* x_tau + (P_e_tau)* x_e  + (P_mu_tau) * x_mu  )

        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        self.Spectra.flux_earth['e'  ]           =  self.Spectra.flux_source['e']
        self.Spectra.flux_earth['neutrinos_e'  ]['dndlogx'] = nue_earth
        self.Spectra.flux_earth['neutrinos_mu' ]['dndlogx'] = numu_earth
        self.Spectra.flux_earth['neutrinos_tau']['dndlogx'] = nutau_earth


    def calculate_fluxes(self):
        # declaring what to do
        if 'PPPC' in self.maddm_card['indirect_flux_source_method']  or 'cosmixs' in self.maddm_card['indirect_flux_source_method'].lower():
            logger.info(f"Calculating cosmic rays fluxes using gammas and neutrinos spectra from the {self.maddm_card['indirect_flux_source_method']} tables.")
            if not self.options['pppc4dmid_path']:
                logger.error("PPPC4DMID not installed, will not calculate fluxes.")
                return
        elif 'pythia' in self.maddm_card['indirect_flux_source_method'] or 'vincia' in self.maddm_card['indirect_flux_source_method']:
            logger.info('Calculating cosmic rays fluxes using pythia8 ' + '(vincia) '*(self.maddm_card['indirect_flux_source_method'] == 'vincia') + 'gammas and neutrinos spectra.')
        else:
            return
             
        # self.Spectra.spectra = = {'x':[] , 'antiprotons':[], 'gammas':[], 'neutrinos_e':[], 'neutrinos_mu':[], 'neutrinos_tau':[], 'positrons':[] }
        cr_spectra = list(self.Spectra.spectra.keys())
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        if str(self.mode['indirect']).startswith('flux'):

             x = np.array (self.Spectra.spectra['x'] )
             E = x*mdm
             self.Spectra.flux_source['e'] = E
             for spec in cr_spectra:
                   if spec in ['x', 'log10x']: continue
                   if ('positrons' in spec) or ('antiprotons' in spec) or ('antideuterons' in spec): continue
                   self.dNdE_dPhidE(channel = spec)           
                   Phi = self.Phi(chan= spec) # Phi takes care of the J factor so it is calculate at earth!
                   
                   #Only needed for dPhidE
                   #for energy in energies:
                   #    self.dNdE_dPhidE( energy,  x=x, dndlogx=dndlogx)
                   #    dNdE   = self.dNdE_int
                   #    dPhidE = self.dPhidE
                   #    dNdE_a  .append(dNdE)
                   #    dPhidE_a.append(dPhidE)
                   #self.Spectra.flux_source[chan]['dNdE']   = dNdE_a
                   #self.Spectra.flux_source[chan]['dPhidE'] = dPhidE_a
                   
                   self.last_results['flux_%s' % spec] = Phi 

    # dndlogx for each different channel
    def dNdE_dPhidE(self, channel=''):
        # input: values of x=Ekin/mDM and dn/dlogx (as PPPC tables)
        # - transforms into E and dN/dE (source: this is needed by DRAGON )
        # - gives the differential flux dPhi/dE (calculated for the same values as dE)
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])

        # spectra input
        E = self.Spectra.flux_source['e']
        # astrophysical input
        sigv = self.last_results['taacsID']
        halo_profile = self.maddm_card['halo_profile']
        jfact = self.maddm_card['jfactors'][halo_profile]

        # *** generic differential flux dPhi/dE at earth        
        if 'gammas' in channel:
           dndlogx = np.array( self.Spectra.spectra[channel] )
           dNdE = dndlogx / (E*2.30259) 
           self.Spectra.flux_source[channel]['dNdE'] = dNdE  # simply convertion from dNdlogx to dN/dE
           dPhidE = 1.0/(self.norm_Majorana_Dirac() *2* math.pi*mdm*mdm)*sigv*jfact * dNdE  # diff.flux for the energy == interp.
                   
        elif 'neutrinos' in channel:
           dndlogx = np.array( self.Spectra.flux_earth[channel]['dndlogx'] )
           dNdE = dndlogx / (E*2.30259)
           self.Spectra.flux_source[channel]['dNdE'] = dNdE  # simply convertion from dNdlogx to dN/dE                                                                             
           dPhidE = 1.0/(self.norm_Majorana_Dirac() *2* math.pi*mdm*mdm)*sigv*jfact * dNdE  # diff.flux for the energy == interp.  
            
        self.Spectra.flux_earth[channel]['dPhidE'] = dPhidE

    
    def dNdE_int(self,energy, channel= ''):
        if channel != 'positrons':
           # returns the interpolated value for the Phi integral calculation
           E    = self.Spectra.flux_source['e']
           dNdE = self.Spectra.flux_source[channel]['dNdE']
           return np.interp(energy, E, dNdE )

        else:
            E = self.Spectra.flux_earth_positrons['e']
            dPhidlogE = self.Spectra.flux_earth_positrons['positrons']['dPhidlogE']
            return np.interp(energy, E, dPhidlogE )

    def dNdlog10x_int(self, log10x, channel= ''):
        return np.interp(log10x, self.Spectra.spectra['log10x'], self.Spectra.spectra[channel])

    def Phi(self, chan= '' , positrons = False):
        # this function integrates over the energy the differential spectrum to give the total flux 
        mdm = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        sigv = self.last_results['taacsID']
        halo_profile = self.maddm_card['halo_profile']
        jfact = self.maddm_card['jfactors'][halo_profile]

        npts =  self.maddm_card['npts_for_flux']

        if not positrons:
            x_min, x_max = self.Spectra.spectra['log10x'][[0,-1]]
            grid = np.linspace(x_min, x_max, npts+1)
            integrate_dNdE = aux.integrate(self.dNdlog10x_int, grid, channel = chan) # the last are the arguments of the function                                         
            phi = 1.0/(self.norm_Majorana_Dirac()*2*math.pi*mdm*mdm)*sigv*jfact*integrate_dNdE                                                                                     

        else:
            grid = [de*mdm/npts for de in range(0, npts+1)]
            integrate_dPhidE = aux.integrate(self.dNdE_int, grid, channel = 'positrons')
            phi = integrate_dPhidE # flux at Earth already calculated in tables
            
        return phi

        
    def norm_Majorana_Dirac(self):
        # (0=own anti) -> this works when there is the factor 2pi in the denominator
        dirac_maj = self.param_card.get_value('qnumbers ' + str (self.proc_characteristics['dm_candidate'][0] ), 4) 
        if   dirac_maj == 0:
            return 4
        elif dirac_maj == 1: 
            return 8 
    
    def run_Dragon(self):
        """ calling the software DRAGON for CR propagation (only if the DM halo velocity is compatible with Milky Way)"""

        # check that the vave is in the range (150-350 km/s) = (5*10^-4 - 1*10^-3)/c 
        halo_vel = self.maddm_card['vave_indirect_cont']  
        if halo_vel < (5*10**(-4)) or halo_vel > ( 10**(-3) ):
           logger.error('The DM velocity in the Milky Way is not in the [5*10^-4 - 1*10^-3]/c range - will not run Dragon!')
           logger.error('Please rerun with the correct velocity to re-calculate the correct sigmav for CR propagation.')
           return 

        run_name = self.last_results['run']        
        dragon_dir = self.options['dragon_path']
        dragon_bin = pjoin(dragon_dir,'bin','DRAGON')

        if not dragon_dir:
            exe = misc.which('DRAGON')
            if exe:
                dragon_dir = os.path.dirname(exe)
        
        if not dragon_dir or not os.path.exists(dragon_bin):
            logger.error('The DRAGON executable cannot be find! Will not run positrons and antiprotons propagation!')
            logger.info("To specify the path to DRAGON_PATH, please use command 'set dragon_path DRAGON_DIRECTORY'.")
            return
       
        template_card = pjoin(MDMDIR,'Templates','Cards','dragon_card.xml')
        mDM = self.param_card.get_value('mass', self.proc_characteristics['dm_candidate'][0])
        if self.maddm_card['sigmav_method'] =='inclusive': sigv = self.last_results['tot_SM_xsec']
        else:                                              sigv = self.last_results['taacsID']

        out_dir = pjoin(self.dir_path, 'output', run_name)
        dragon_input = pjoin( self.dir_path, 'output' , run_name, run_name+'_DRAGON.xml')

        def write_dragon_input(template= template_card , mdm = mDM , sigmav = sigv , dragon_input = dragon_input ):
            energy , positrons, antiprotons = self.Spectra.flux_source['e'] , self.Spectra.flux_source['positrons']['dNdE'] , \
                                                                              self.Spectra.flux_source['antiprotons']['dNdE']
            
            out_pos  = open(pjoin(out_dir, 'positrons_dndne.txt') , 'w')
            out_anti = open(pjoin(out_dir, 'antiprotons_dndne.txt') ,'w')

            for en, pos, anti in zip (energy, positrons, antiprotons):
                out_pos .write(str(en)+ ' ' + str(pos)  +'\n')
                out_anti.write(str(en)+ ' ' + str(anti) +'\n')
            out_pos .close()
            out_anti.close()

            xml_in = open( dragon_input, 'w')
            xml_in.write(open(template,'r').read() %
                         {'dm_mass': mdm,
                          'sigmav': sigmav,
                          'Positrons': pjoin(out_dir, 'positrons_dndne.txt'),
                          'Antiprotons' : pjoin(out_dir, 'antiprotons_dndne.txt')
                          })
            xml_in.close()

        write_dragon_input(template= template_card , mdm = mDM , sigmav = sigv , dragon_input = dragon_input )
        
        misc.call([dragon_bin, dragon_input], cwd=dragon_dir)

    #        if not __debug__:
    #            logger.debug('keep dragon positrons/antiprotons files due to debug mode.')
    #            os.remove(pjoin(out_dir,'positrons_dndne.txt') )
    #            os.remove(pjoin(out_dir,'antiprotons_dndne.txt') )

    
    def save_MadDM_card(self,scan=False):
        # saves the maddm run card inside the output directory
        run_name = self.run_name
        if not scan:
            if not os.path.isfile(pjoin(self.dir_path, 'output', run_name, 'maddm_card.dat') ):
                shutil.copyfile (pjoin(self.dir_path, 'Cards', 'maddm_card.dat') , pjoin(self.dir_path, 'output', run_name, 'maddm_card.dat') )
        else:
            run_name = run_name.rsplit('_',1)[0]
            if not os.path.isfile(pjoin(self.dir_path, 'output', run_name+'_maddm_card.dat') ):
                shutil.copyfile (pjoin(self.dir_path, 'Cards', 'maddm_card.dat') , pjoin(self.dir_path, 'output', run_name+'_maddm_card.dat') )
        
    
    def write_scan_output(self, out_path = '', keys = '', header = False, write_strings = False):
        if out_path and header:
            nice_keys = []
            for k in keys:
                if 'taacsID#' in k:
                    k = k.replace('taacsID#','')
                # if '%_relic_' in k:
                #     k = k.replace('%_relic_','')
                #     k = '%_relic_' + self.pdg_particle_map.format_particles(k)
                k = k.replace('taacsID','tot_Xsec')
                nice_keys.append(k)
            
            summary = open(out_path, 'w')                                                                                                                            
            
            for k in nice_keys:
                ind = nice_keys.index(k) + 1
                if ind <=9: ind = '0'+str(ind)
                summary.write( '# [' + str(ind) + ']' + ' : ' + k + '\n' )

            summary.write('\n\n\n')
            summary.close()

            self.save_MadDM_card(scan=True) 
        elif (out_path and not header):
            s = '\t'
            summary = open(out_path, 'a+')
            summary.write('%s\t' % self.last_results['run'])

            formatter = self.form_n if not write_strings else self.form_s

            for k in keys:
                if 'run' in k: continue 
                num = formatter( self.last_results[k] )
                summary.write(num + s)
            summary.write('\n')
            summary.close()
    
    def print_results(self):    
        """ print the latest results on screen """
        
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

        under_message  = '%s UNDERABUNDANT     %s' % ('\033[33m', bcolors.ENDC)
        above_message  = '%s OVERABUNDANT      %s' % (bcolors.FAIL, bcolors.ENDC)
        within_message = '%s WITHIN EXP ERROR  %s' % (bcolors.OKGREEN, bcolors.ENDC)

        # skip this if there is a sequential scan going on.
        if self.last_results['Omegah^2'] < 0.:
            pass_relic = ""
        if self.last_results['Omegah^2'] < self.limits._oh2_planck - 2*self.limits._oh2_planck_width:
            pass_relic = under_message
        elif self.last_results['Omegah^2'] > self.limits._oh2_planck + 2*self.limits._oh2_planck_width:
            pass_relic = above_message
        else:
            pass_relic = within_message                                                                

        omega    = self.last_results['Omegah^2']
        x_f      = self.last_results['x_f']
        sigma_xf = self.last_results['sigmav(xf)']
        xsi      = self.last_results['xsi']        

        logger.info("MadDM Results", '$MG:BOLD')
        if dm_scen: 
              logger.info("Define xsi = Relic density/Planck measurement for thermal scenarios.")
              logger.info("Rescaling theory prediction for xsi(direct det.) and xsi^2(indirect det.) for thermal scenarios.\n")

        if self.mode['relic']:
            logger.info("\n****** Relic Density\nOMEGA IS %f" % omega)
            logger.info( self.form_s('Relic Density') + '= ' + self.form_n(omega)   + '      '  +  pass_relic )
            logger.info( self.form_s('x_f'          ) + '= ' + self.form_s(self.form_n(x_f      ) ) )
            logger.info( self.form_s('sigmav(xf)'   ) + '= ' + self.form_s(self.form_n(sigma_xf ) + ' cm^3 s^-1') )
            logger.info( self.form_s('xsi'          ) + '= ' + self.form_s(self.form_n(xsi      ) ) )
            logger.info('')
            logger.info('Channels contributions:')
            skip = []
            for proc in [k for k in self.last_results.keys() if k.startswith('%_relic_')]:
                if self.last_results[proc] == 0.:
                    skip.append(proc.replace(r"%_relic_", ""))
                    continue
                logger.info( self.form_s(proc.replace(r"%_relic_", "")) + ': %.2f %%' % self.last_results[proc] )
            if len(skip) != 0:
                logger.info('No contribution from processes: %s', ', '.join(skip))

        if self.mode['direct']:
            # units = self.last_results['GeV2pb*pb2cm2']
            direct_names = [ { 'n': 'SigmaN_SI_p', 'sig': self.last_results['sigmaN_SI_p'], 'lim': self.last_results['lim_sigmaN_SI_p'], 'exp': 'Xenon1ton' },
                             { 'n': 'SigmaN_SI_n', 'sig': self.last_results['sigmaN_SI_n'], 'lim': self.last_results['lim_sigmaN_SI_n'], 'exp': 'Xenon1ton' },
                             { 'n': 'SigmaN_SD_p', 'sig': self.last_results['sigmaN_SD_p'], 'lim': self.last_results['lim_sigmaN_SD_p'], 'exp': 'Pico60'    },
                             { 'n': 'SigmaN_SD_n', 'sig': self.last_results['sigmaN_SD_n'], 'lim': self.last_results['lim_sigmaN_SD_n'], 'exp': 'Lux2017'   } ]

            self.last_results['direct_results'] = direct_names

            logger.info('\n****** Direct detection [cm^2]: ')  
            for D in direct_names:
                    th_cross = D['sig']* xsi 
                    ul = D['lim']
                    mess_th       = self.det_message_screen(th_cross, D['lim'] )
                    mess_alldm    = self.det_message_screen(D['sig'], D['lim'] )
                    self.print_ind(D['n'],th_cross , D['sig'], ul,  thermal=dm_scen ,direc = False , exp=D['exp']) # leave False to avout GeV^-2
         
#        if self.mode['direct'] == 'directional':
#            logger.info(' Nevents          : %i', self.last_results['Nevents'])
#            logger.info(' smearing         : %.2e', self.last_results['smearing'])
        
        if self.mode['capture']:
            logger.info('\n capture coefficients: ')
            detailled_keys = [k for k in self.last_results.keys() if k.startswith('ccap')]
            for key in detailled_keys:
                logger.info(' %s            : %.2e s^-1' % (key, self.last_results[key]))

        # INDIRECT DETECTION
        detailled_keys = [k for k in self.last_results.keys() if k.startswith('taacsID#')]

        def print_sigmav_with_limits(detailled_keys, filter_, exp_label, no_lim):
            ''' filter is a lambda for the processes to display 
                this method returns True if there are light quarks in the final states
            '''
            skip = []
            light_s = False
            if len(detailled_keys) > 0:
                for key in detailled_keys:
                    process = key.split('#')[-1]
                    if not filter_(process):
                        continue
                    finalstate = ".".join(self.str_processes[process].final_states)
                    if '3.-3' in finalstate or '2.-2' in finalstate or '1.-1' in finalstate:
                        light_s = True
                    s_alldm   = self.last_results[key]
                    s_thermal = s_alldm * xsi2
                    if s_alldm <= 10**(-100): 
                        skip.append(process)
                        continue
                    s_ul = self.last_results['lim_' + key]
                    self.print_ind(process, s_thermal, s_alldm, s_ul, exp = exp_label, thermal= dm_scen, no_lim = no_lim) 

            if len(skip) >= 1:        
                logger.info('Skipping zero cross section processes for: %s', ', '.join(skip))
            return light_s
        
        if self.mode['indirect'] or self.mode['spectral']:
            logger.info('\n****** Indirect detection:')
            logger.info('<sigma v> method: %s ' % self.maddm_card['sigmav_method'])

        if self.mode['indirect']:
            logger.info('====== continuum spectrum final states')
            halo_vel = self.maddm_card['vave_indirect_cont']
            logger.info('DM particle halo velocity: %s c ' % halo_vel)
            logger.info('*** Print <sigma v> [cm^3 s^-1] with Fermi dSph limits')
            sigtot_alldm     = self.last_results['taacsID']
            sigtot_SM_alldm  = self.last_results['tot_SM_xsec']
            sigtot_th        = sigtot_alldm * xsi2
            sigtot_SM_th     = sigtot_SM_alldm * xsi2

            if halo_vel > self.vave_indirect_cont_range[0] and halo_vel < self.vave_indirect_cont_range[1]:
                light_s = print_sigmav_with_limits(detailled_keys, filter_ = lambda process: not self.is_spectral_finalstate(self.str_processes[process.replace("taacsID#","")]), exp_label = 'Fermi dSph', no_lim = False)
                if light_s:
                    logger.info('Using generic Fermi dSph limits for light quarks (u,d,s)')

                if self.maddm_card['sigmav_method']!= 'inclusive':     
                    self.print_ind('DM DM > all',sigtot_th , sigtot_alldm, self.last_results['Fermi_sigmav'],  thermal= dm_scen)
                else: 
                    self.print_ind('DM DM > SM SM', sigtot_SM_th , sigtot_SM_alldm, self.last_results['Fermi_sigmav'],  thermal= dm_scen)
                logger.info('')
            
                if str(self.mode['indirect']).startswith('flux'):
                    logger.info('*** Fluxes at earth [particle cm^-2 sr^-1]:')
                    #np_names = {'gammas':'g'      , 'neutrinos_e':'nue' , 'neutrinos_mu':'numu' , 'neutrinos_tau':'nutau'}
                    for chan in ['gammas','neutrinos_e', 'neutrinos_mu' , 'neutrinos_tau']:
                        logger.info( self.form_s(chan + ' Flux') + '=\t' + self.form_s(self.form_n (self.last_results['flux_%s' % chan]) ))
                    logger.info('')
            
            else:
                light_s = print_sigmav_with_limits(detailled_keys, filter_ = lambda process: not self.is_spectral_finalstate(self.str_processes[process.replace("taacsID#","")]), exp_label = 'Fermi dSph', no_lim = True)
                if light_s:
                    logger.info('Using generic Fermi dSph limits for light quarks (u,d,s)')
                
                if self.maddm_card['sigmav_method']!= 'inclusive':     
                    self.print_ind('DM DM > all',sigtot_th , sigtot_alldm, self.last_results['Fermi_sigmav'],  thermal= dm_scen, no_lim = True)
                else: 
                    self.print_ind('DM DM > SM SM', sigtot_SM_th , sigtot_SM_alldm, self.last_results['Fermi_sigmav'],  thermal= dm_scen, no_lim = True)
            
                logger.info('')
                logger.warning('Fermi dSph limits cannot be calculated since DM halo velocity not compatible with dSph.')
   
        if self.mode['spectral']:
            logger.info('====== line spectrum final states')
            halo_vel = self.maddm_card['vave_indirect_line']
            logger.info('DM particle halo velocity: %s c ' % halo_vel)
            logger.info('*** Print <sigma v> [cm^3 s^-1] with %s line limits' % ', '.join([name.replace('_', ' ') for name in self.line_experiments.iternames()]))

            # if halo_vel > self.vave_indirect_line_range[0] and halo_vel < self.vave_indirect_line_range[1]:
            print_sigmav_with_limits(detailled_keys, filter_ = lambda process: self.is_spectral_finalstate(self.str_processes[process.replace("taacsID#","")]), exp_label = 'Line GC', no_lim = False)

            if not (halo_vel > self.vave_indirect_line_range[0] and halo_vel < self.vave_indirect_line_range[1]):
                logger.warning('DM halo velocity not compatible with GC, line limits application is this case is not guaranteed.')
            logger.info('')
            logger.info('*** Line limits from ' + ', '.join([name.replace('_', ' ') for name in self.line_experiments.iternames()]))
            self.print_line_results(xsi2)

            # else:
            #     print_sigmav_with_limits(detailled_keys, filter_ = lambda process: self.is_spectral_finalstate(process.split('_')[-1]), exp_label = 'Line GC', no_lim = True)
            #     logger.info('')
            #     logger.warning('Line limits cannot be calculated since DM halo velocity not compatible with GC.')

        # Internal: save the results as a numpy dictionary
        # np.save(pjoin(self.dir_path, 'output','Results'), self.last_results)

        logger.info('')

        if not self.param_card_iterator:
            self.save_summary_single(relic = self.mode['relic'], direct = self.mode['direct'], \
                                indirect = self.mode['indirect'], spectral = self.mode['spectral'],
                                fluxes_source= self.mode['indirect'].startswith('flux') if isinstance(self.mode['indirect'],str) else self.mode['indirect'], 
                                fluxes_earth = False )                  
            logger.info('Results written in: ' + pjoin(self.dir_path, 'output', self.run_name, 'MadDM_results.txt') )

            # self.save_remove_output(scan = False)

    def save_direct_detection(self, output_dir, output_run):
        ''' direct detection output is created in the output directory, move to run directory '''
        for file_name in ['d2NdEdcos.dat','dNdcos.dat','dNdE.dat','rate.dat']:
            f_path = pjoin(output_dir, file_name)
            if os.path.isfile(f_path):
                shutil.move( f_path, pjoin(output_run, file_name) )

    def remove_direct_detection(self, output_dir):
        ''' direct detection output is created in the output directory, remove those files '''
        for file_name in ['d2NdEdcos.dat','dNdcos.dat','dNdE.dat','rate.dat']:
            f_path = pjoin(output_dir, file_name)
            if os.path.isfile(f_path):
                os.remove(f_path)

    def create_symlink_to_events(self, indirect_directory, events_run, output_run):
        ''' when sigmav_method != inclusive, then we need to create a symlink to the events_run directory in the output_dir '''
        if not os.path.islink( pjoin(output_run, 'Output_' + indirect_directory) ):
            os.symlink( pjoin(events_run), pjoin(output_run, 'Output_' + indirect_directory) )  

    def remove_lhe_file(self, events_run):
        ''' removes the .lhe and .sh files from the events_run directory '''
        if os.path.isfile( pjoin(events_run,'unweighted_events.lhe.gz') ):
            os.remove( pjoin(events_run,'unweighted_events.lhe.gz') )
        if os.path.isfile( pjoin(events_run,'run_shower.sh') ):
            os.remove( pjoin(events_run,'run_shower.sh') )

    def remove_events_run(self, events_run):
        ''' remove the entire events_run directory '''
        if os.path.isdir(events_run):
            shutil.rmtree(events_run)

    def output_single_directory(self, indirect_directory, scan = False):
        ''' saves/removes the output related to a single directory '''
        # run name
        run_name = self.last_results['run']

        # sigmav method
        sigmav_method = self.maddm_card['sigmav_method']
        # save_output switch
        save_switch = self.maddm_card['save_output']

        # Setting the various paths
        source_indirect = pjoin(self.dir_path  , indirect_directory)
        events_dir      = pjoin(source_indirect, 'Events'          )
        events_run      = pjoin(events_dir     , run_name          )
        output_dir      = pjoin(self.dir_path  , 'output'          )
        output_run      = pjoin(output_dir     , run_name          )

        # check if the indirect_directory exists: if it does not exist, then do nothing
        if not os.path.isdir(source_indirect):
            logger.debug("'%s' directory does not exist. Nothing to do here." % indirect_directory)
            return

        operations = [
            lambda **kwargs: self.create_symlink_to_events(indirect_directory = kwargs['indirect_directory'], events_run = kwargs['events_run'], output_run = kwargs['output_run']),
            lambda **kwargs: self.remove_events_run(events_run = kwargs['events_run']),
            lambda **kwargs: self.remove_lhe_file(events_run = kwargs['events_run'])
        ]

        if not scan:
            # output for single point run
            # save_switch has no effect in this case
            do_operation = [
                sigmav_method != 'inclusive',
                False,
                False
            ]
        else:
            # output for scan
            # everything depends on save_switch
            do_operation = [
                sigmav_method != 'inclusive' and save_switch != 'off',
                save_switch == 'off', # remove the entire events_run directory if 'off'
                save_switch == 'spectra' # remove lhe file only in this case, for 'all' keep it, for 'off' it already removes the entire events_run directory
            ]

        for do_it, operation in zip(do_operation, operations):
            if not do_it:
                continue
            operation.__call__(
                indirect_directory = indirect_directory,
                events_run         = events_run,
                output_run         = output_run,
            )

    def output_entire_run(self, scan = False):
        ''' saves/removes the output related to the entire run '''
        # indirect mode
        ind_mode = self.mode['indirect']
        # sigmav method
        sigmav_method = self.maddm_card['sigmav_method']
        # save_output switch
        save_switch = self.maddm_card['save_output']
        # spectrum method
        source_method = self.maddm_card['indirect_flux_source_method']
        earth_method  = self.maddm_card['indirect_flux_earth_method']

        # Setting the various paths
        output_dir      = pjoin(self.dir_path, 'output'                )
        output_run      = pjoin(output_dir   , self.last_results['run'])

        # define operations
        operations = [
            lambda **kwargs: self.save_MadDM_card(),
            lambda **kwargs: self.save_direct_detection(output_dir = kwargs['output_dir'], output_run = kwargs['output_run']),
            lambda **kwargs: self.save_spec_flux(out_run = kwargs['output_run'], spec_source = kwargs['spec_source'], flux_source = kwargs['flux_source'], flux_earth = kwargs['flux_earth']),
            lambda **kwargs: self.remove_direct_detection(output_dir = kwargs['output_dir'])
        ]

        # PPPC4DMID files
        # during launch_indirect:
        # if ind_mode != sigmav, then 'source' spectra are always computed, while 'earth' is done only if it is selected
        # that means it computes 'source' spectra even if it is selected 'earth'
        # of course this is valid only for ind_mode != OFF
        spec_source = ind_mode and ind_mode != 'sigmav'
        flux_source = ind_mode and 'source' in ind_mode
        flux_earth  = ind_mode and 'earth' in ind_mode and 'PPPC4DMID' in earth_method # positron flux computed inside the variables related to save_spec_flux only if the choice is PPPC4DMID, because DRAGON already writes its own output files.

        if not scan:
            # output for single point run
            # save_switch has no effect in this case
            do_operation = [
                True,
                True,
                True,
                False
            ]
        else:
            # output for scan
            # everything depends on save_switch
            do_operation = [
                False,
                save_switch == 'all',
                save_switch != 'off',
                save_switch != 'all'
            ]

        for do_it, operation in zip(do_operation, operations):
            if not do_it:
                continue
            operation.__call__(
                output_dir  = output_dir,
                output_run  = output_run,
                spec_source = spec_source,
                flux_source = flux_source,
                flux_earth  = flux_earth
            )

    def save_remove_output(self, scan = False):
        ''' It saves the various output files in case of single point.
            It saves the spectra files in each directory.
            It saves/removes outputs according to the 'save_output' switch in scan mode.
        '''
        # run name
        run_name = self.last_results['run']
        assert run_name == self.run_name
        # manage outputs
        self.output_entire_run(scan)
        for indirect_directory in self.indirect_directories.keys():
            self.output_single_directory(indirect_directory, scan)

    # This aux function saves the spectra/fluxes in the given directory out_dir
    def save_spec_flux(self, out_run, spec_source = False, flux_source = False, flux_earth = False):

        spec_method = self.maddm_card['indirect_flux_source_method']
       
        if spec_source:
            x = self.Spectra.spectra['log10x']
            for spec in self.Spectra.spectra.keys():
                header = '# ' + spec + '\t' + self.maddm_card['indirect_flux_source_method'] + ' spectra at source\n' + '# Log10(x=Ekin/mDM)    dn/dlogx'
                if 'x' not in spec:
                    dndlogx = self.Spectra.spectra[spec]
                    errors = self.Spectra.errors[spec]
                    aux.write_data_to_file(x, dndlogx, errors = errors, filename = pjoin(out_run, spec + '_spectrum_' + spec_method +'.dat'), header = header)

        if flux_earth:
            e = self.Spectra.flux_source['e']
            for flux in self.Spectra.flux_earth.keys():
                header = '# ' + flux + '\t' + self.maddm_card['indirect_flux_earth_method'] +' flux at earth\n' + '# E[GeV]    dPhi/dE [particles/(cm2 s sr)]'
                if flux != 'e':
                    dPhidE = self.Spectra.flux_earth[flux]['dPhidE']
                    aux.write_data_to_file(e, dPhidE, errors = False, filename = pjoin(out_run, flux + '_dphide_' + spec_method + '.dat'), header = header)


            if 'PPPC' in self.maddm_card['indirect_flux_earth_method']: # for PPPC4DMID, I save also the positron at earth
                e = self.Spectra.flux_earth_positrons['e']
                header = '# positrons \t PPPC4DMID flux at earth\n' + '# E[GeV]    dPhi/dlogE [particles/(cm2 s sr)]'
                dPhidlogE = self.Spectra.flux_earth_positrons['positrons']['dPhidlogE']
                aux.write_data_to_file(e, dPhidlogE, errors = False, filename = pjoin(out_run, 'positrons_dphide_' + spec_method + '.dat'), header = header)

    
    def print_ind(self,what, sig_th, sig_alldm, ul,  thermal=False ,direc = False , exp='Fermi dSph', no_lim = False):
        ''' to print the output on screen. Set no_lim = True to avoid printing the experimental ul (e.g. Fermi dSph ul for wrong dm velocity) '''
        
        alldm_mess = self.det_message_screen(sig_alldm , ul)

        if not direc and not no_lim:

            if thermal:
                th_mess    = self.det_message_screen(sig_th    , ul)
                line = self.form_s(what) + self.form_s('Thermal = '+ self.form_n(sig_th))  + '   ' + self.form_s(th_mess)
                line = line +              '\t' + self.form_s('All DM = '  + self.form_n(sig_alldm) ) + '   ' + self.form_s(alldm_mess)
                line = line + '\t' + '{0:16}'.format(exp+' ul') + ' = ' +  self.form_n(ul)

            else:                                                                                                                                                           
                line = self.form_s(what)  + self.form_s('All DM = ' + self.form_n(sig_alldm) ) + '   ' + self.form_s(alldm_mess)                               
                line = line + '\t' + '{0:16}'.format(exp+' ul') + ' = ' +  self.form_n(ul)

        elif not direc and no_lim:

            if thermal:
                th_mess    = self.det_message_screen(sig_th    , -1)
                line = self.form_s(what) + self.form_s('Thermal = ' +  self.form_n(sig_th))  + '   ' + self.form_s(th_mess)
                line = line +              '\t' + self.form_s('All DM = ' + self.form_n(sig_alldm) ) + '   ' + self.form_s(alldm_mess)
            else:
                line = self.form_s(what)+ self.form_s('All DM = ' + self.form_n(sig_alldm)) + '   ' + self.form_s(alldm_mess)

        logger.info(line)
    
    def print_line_results(self, xsi2):
        ''' method to print and format the line analysis results '''
        def colored_message(n1, n2):
            ''' don't use det_message_screen because you can't ensure alignment:
                there the colors are already inside the string, so they are kept into account
                while filling the available space;
                so return the color and the message separately to handle alignment.
            '''
            if n2 <= 0:
                return "NO LIMIT", bcolors.GRAY
            elif n1 > n2 and n2 > 0:
                return "EXCLUDED", bcolors.FAIL
            elif n2 > n1:
                return "ALLOWED", bcolors.OKGREEN
        logger.info("Density profile: %s" % self.last_results["density_profile"])
        for line_exp in self.line_experiments:
            str_part = "line_%s_" % line_exp.get_name()
            logger.info("==== %s %s====" % (line_exp.get_name().replace('_', ' '), "="*max([0, 100 - len(line_exp.get_name())])))
            logger.info("ROI: %.1f" % self.last_results[str_part + "roi"])
            str_part_peak = str_part + 'peak'
            energy_peaks = collections.OrderedDict(sorted([(k, v) for k, v in self.last_results.items() if str_part_peak in k and not any(['states' in k, 'error' in k, 'like' in k, 'pvalue' in k]) and v != -1], key = lambda item: item[1])) # key = "line_<exp_name>_peak_<num>", value = energy peak
            roi_warning = self.last_results[str_part + "roi_warning"]
            if roi_warning != "":
                logger.warning(roi_warning)
            logger.info("J = %.6e GeV^2 cm^-5" % self.last_results[str_part + "Jfactor"])
            logger.info("detection range: %.4e -- %.4e GeV" % (line_exp.detection_range[0], line_exp.detection_range[1]))
            if len(energy_peaks) is 0:
                logger.info(bcolors.BOLD + "No peaks found: out of detection range." + bcolors.ENDC)
            else:
                # find first column maximum length for nice table format
                # for each peak, the first column has a length = len("peak" + "<peak_number>" + "(<label_merged>)")
                first_col_len = max([len("peak_%d(%s)" % (int(k.split('_')[-1]), self.last_results[k + '_states'])) for k in energy_peaks.keys()])
                row_format = r"{" + "0:%ds" % first_col_len + r"}"
                if self.mode['relic']:
                    mid_line_len = 10+1+8+3+10+1+8+3+10
                    headers_low = "{0:s}   {2:^12s}   {3:^10s} {1:8s}   {4:^10s} {1:8s}   {5:^10s}".format(" "*first_col_len, "", "", "All DM", "Thermal", "UL")
                else:
                    mid_line_len = 10+1+8+3+10
                    headers_low = "{0:s}   {2:^12s}   {3:^10s} {1:8s}   {4:^10s}".format(" "*first_col_len, "", "", "All DM", "UL")
                len_headers = len(headers_low)
                headers_up_format = r"{" + "0:%ds" % (len_headers - mid_line_len) + r"}{" + "1:^%ds" % (mid_line_len) + r"}"
                headers_mid_format = r"{" + "0:%ds" % first_col_len + r"}{" + "1:^%ds" % (len_headers - mid_line_len - first_col_len) + r"}{" + "2:-^%ds" % (mid_line_len) + r"}"
                logger.info("-"*len_headers)
                logger.info(headers_up_format.format(" ", "Flux [cm^-2 s^-1]"))
                logger.info(headers_mid_format.format(" ", "Energy [GeV]", ""))
                logger.info(headers_low)
                logger.info("-"*len_headers)
                for k, peak in energy_peaks.items():
                    num = int(k.split('_')[-1])
                    error_code = self.last_results[str_part + "peak_%d_error" % num]
                    flux    = self.last_results[str_part + "flux_%d"    % num]
                    flux_UL = self.last_results[str_part + "flux_UL_%d" % num]
                    allowed_or_excluded, color       = colored_message(flux, flux_UL)
                    allowed_or_excluded_th, color_th = colored_message(flux*xsi2, flux_UL)
                    error_string = "   %s# error: %s%s" % ('\033[33m', error_code, bcolors.ENDC) if error_code != 0 else ''
                    if self.mode['relic']:
                        logger.info(row_format.format("peak_%d(%s)" % (num, self.last_results[k + "_states"])) + "   {0:^12.4e}   {1:^10.4e} {6:s}{2:<8s}{8:s}   {3:^10.4e} {7:s}{4:<8s}{8:s}   {5:^10.4e}".format(peak, flux, allowed_or_excluded, flux*xsi2, allowed_or_excluded_th, flux_UL, color, color_th, bcolors.ENDC) + error_string)
                    else:
                        logger.info(row_format.format("peak_%d(%s)" % (num, self.last_results[k + "_states"])) + "   {0:^12.4e}   {1:^10.4e} {4:s}{2:<8s}{5:s}   {3:^10.4e}".format(peak, flux, allowed_or_excluded, flux_UL, color, bcolors.ENDC) + error_string)
                    try:
                        like    = self.last_results[str_part + "like_%d"    % num]
                        pvalue  = self.last_results[str_part + "pvalue_%d"  % num]
                        logger.debug("peak_%d: loglike = %.4e, p-value = %f" % (num, like, pvalue))
                    except KeyError:
                        continue
                logger.info("-"*len_headers)

    def save_summary_single(self, relic = False, direct = False , indirect = False , spectral = False , fluxes_source = False , fluxes_earth = False):

        point = self.last_results['run']
        # creting symlink to Events folder in Inidrect output directory
        #if not os.path.exists(pjoin(self.dir_path, 'output' , point)) :
        #    os.symlink(pjoin(self.dir_path,'Indirect','Events'), pjoin(self.dir_path, 'output' , 'Output_Indirect') )
                    
        def form_s(stringa):
            formatted = '{0:30}'.format(stringa)
            return  formatted
        
        def form_n(num):
            formatted = '{0:3.2e}'.format(num)
            return formatted

        out = open(pjoin(self.dir_path, 'output', point, 'MadDM_results.txt'),'w')
 
        out.write('################################################\n')
        out.write('#                 MadDM v. ' + str(self.MadDM_version) +'                 #\n' )
        out.write('################################################\n\n\n')

        if relic:
            out.write('################################################\n')
            out.write('# Relic Density                                #\n')
            out.write('################################################\n\n')

            relic, planck , message = self.last_results['Omegah^2'] , self.limits._oh2_planck , self.det_message(self.last_results['Omegah^2'], self.limits._oh2_planck) 
    
            out.write(form_s('Omegah2')       + '= ' + form_n(relic)  + '\n')
            out.write(form_s('Omegah_Planck') + '= ' + form_n(planck) + '\n')
            if self.last_results['xsi'] > 0: 
                out.write(form_s('xsi')   + '= ' + form_n(self.last_results['xsi'])        + ' \t # xsi = (Omega/Omega_Planck)\n' )
            out.write(form_s('x_f')       + '= ' + form_n(self.last_results['x_f'])        + '\n' ) 
            out.write(form_s('sigmav_xf') + '= ' + form_n(self.last_results['sigmav(xf)']) + ' \t # cm^3 s^-1\n' ) 
            out.write("# % of the relic density channels\n")
            for proc in [k for k in self.last_results.keys() if k.startswith('%_relic_')]:
                out.write( form_s("%_" + proc.replace('%_relic_','')) + '= %.2f %%\n' % self.last_results[proc] )


        if direct:

            for name in ['d2NdEdcos.dat','dNdE.dat','dNdcos.dat','rate.dat']:
                if os.path.isfile(pjoin(self.dir_path, 'output',name)):
                   shutil.move(pjoin(self.dir_path, 'output',name) ,  pjoin(self.dir_path, 'output', point , name))
            
            out.write('\n################################################\n')
            out.write('# Direct Detection [cm^2]                      #\n')
            out.write('################################################\n\n')

            for D in self.last_results['direct_results']:
                cross = D['sig']
                ul = D['lim']
                exp = D['exp']
                out.write(form_s(D['n']) + '= ' + form_s('['+ form_n(cross) + ',' + form_n(ul) + ']' ) + '# '+exp + '\n')

        if indirect or spectral:      

            sigmav_meth = self.maddm_card['sigmav_method'] # method actually used
            method = self.maddm_card['indirect_flux_source_method']
            
            out.write('\n################################################\n')
            out.write('# Indirect Detection                           #\n')
            out.write('################################################\n')

            out.write('# Results in brackets display [prediction, upper limit]\n\n')
            out.write('# Annihilation cross section computed with the method: ' + sigmav_meth)

            def collect_processes(processes, filter_):
                proc_list = []
                for name in processes:                                                                                                                                     
                    proc = name.replace('taacsID#','')
                    if not filter_(proc):
                        continue
                    proc_list.append(proc)
                return proc_list

            def print_sigmav(proc_list, fileout):
                ''' writes the line proc_name = sigmav, proc is in PDG format. '''
                for proc in proc_list:
                    #proc = [proc_names for proc_names, proc_pdg in self.processes_names_map.iteritems() if proc == proc_pdg][0] # this allows to convert back from PDG
                    proc_th, proc_ul = self.last_results['taacsID#' + proc] , self.last_results['lim_taacsID#' + proc]
                    fileout.write(form_s(proc) + '= '+ form_s('['+ form_n(proc_th)+',' + form_n(proc_ul)+']') + '\n')

            detailled_keys = [key for key in self.last_results.keys() if key.startswith('taacsID#')]

            cont_procs = collect_processes(detailled_keys, filter_ = lambda process: not self.is_spectral_finalstate(self.str_processes[process.replace("taacsID#","")]))
            if cont_procs:
                out.write('\n# <sigma v>[cm^3 s^-1] of continuum spectrum final states and Fermi dSph limits (if available, else -1)\n')
                print_sigmav(cont_procs, out)
            line_procs = collect_processes(detailled_keys, filter_ = lambda process: self.is_spectral_finalstate(self.str_processes[process.replace("taacsID#","")]))
            if line_procs:
                out.write('\n# <sigma v>[cm^3 s^-1] of line spectrum final states and %s line limits (if available, else -1)]' % ', '.join([name.replace('_', ' ') for name in self.line_experiments.iternames()]) + '\n')
                print_sigmav(line_procs, out)

            if indirect:
                tot_th , tot_ul , tot_sm = self.last_results['taacsID'] , self.last_results['Fermi_sigmav'] , self.last_results['tot_SM_xsec']
                pvalue_th     = self.last_results['pvalue_th']
                likelihood_th = self.last_results['like_th'] 
                pvalue_nonth  = self.last_results['pvalue_nonth']
                likelihood_nonth = self.last_results['like_nonth']
                """
                if sigmav_meth !='inclusive':
                    out.write('# Fermi dSph Limit for DM annihilation computed with Pythia8 spectra  \n\n')
                    out.write(form_s('Total_xsec') +'= ' + form_s('['+ form_n(tot_th) +','+ form_n(tot_ul) +']') + '\n') 

                else:
                    if 'PPPC' in  self.maddm_card['indirect_flux_source_method']: 
                        method =  self.maddm_card['indirect_flux_source_method'] 
                    else: 
                        method = 'PPPC4DMID'
                    out.write('# Fermi dSph Limit computed with ' + method + ' spectra\n\n')
                    lista = [ proc for proc in self.last_results.keys() if 'taacsID#' in proc and 'lim_' not in proc and 'err' not in proc ]
                    for name in lista:
                        proc = name.replace('taacsID#','')
                        proc_th , proc_ul = self.last_results[name] , self.last_results['lim_'+name]
                        out.write(form_s(proc)      + '= '+ form_s('['+ form_n(proc_th)+',' + form_n(proc_ul)+']') + '\n')
               
                    out.write(form_s('TotalSM_xsec')+ '= '+ form_s('['+ form_n(tot_sm)+ ',' + form_n(tot_ul) +']') + '\n')
                """

                out.write('\n# Global Fermi dSph Limit computed with ' + method + ' spectra\n')  
                if sigmav_meth != 'inclusive':                                                                            
                    out.write(form_s('Total_xsec')+ '= '+ form_s('['+ form_n(tot_th)+ ',' + form_n(tot_ul) +']') + '\n')
                else:
                    out.write(form_s('TotalSM_xsec')+ '= '+ form_s('['+ form_n(tot_sm)+ ',' + form_n(tot_ul) +']') + '\n')
                out.write( form_s('Fermi_Likelihood')+ '= '+ form_s(form_n(likelihood_nonth))  +'\n' )
                out.write( form_s('Fermi_pvalue'    )+ '= '+ form_s(form_n(pvalue_nonth    ))  +'\n')

                if self.last_results['xsi'] < 1:
                    out.write( form_s('Fermi_Likelihood(Thermal)')+ '= '+ form_s(form_n(likelihood_th))  +'\n' )
                    out.write( form_s('Fermi_pvalue(Thermal)'    )+ '= '+ form_s(form_n(pvalue_th    ))  +'\n ')

            # line limits
            if spectral:
                out.write('\n# Gamma-line spectrum and line limits\n')
                out.write('# peak[GeV], flux[cm^-2 s^-1], J-factor[GeV^2 cm^-5]\n') #, ROI[deg], rho_s[GeV], r_s[kpc]')
                # out.write('\n' + form_s("Density_profile") + '= %s\n' % self.last_results["density_profile"].get_name())
                # for k, v in self.last_results["density_profile"].get_parameters_items(): # print astrophysical parameters
                #     out.write(form_s(k.replace('profile_', '')) + '= ' + form_n(v) + '\n')
                for line_exp in self.line_experiments:
                    str_part = "line_%s_" % line_exp.get_name()
                    out.write('# %s (ArXiv:%s)\n' % (line_exp.get_name().replace('_', ' '), line_exp.arxiv_number))
                    out.write(form_s("ROI") + '= %1.1f\n' % self.last_results[str_part + "roi"])
                    out.write(form_s("J-factor") + '= ' + form_n(self.last_results[str_part + "Jfactor"]) + '\n')
                    str_part_peak = str_part + 'peak'
                    energy_peaks = collections.OrderedDict(sorted([(k, v) for k, v in six.iteritems(self.last_results) if str_part_peak in k and '_states' not in k and '_error' not in k and v != -1], key = lambda item: item[1])) # key = "line_<exp_name>_peak_<num>", value = energy peak
                    if len(energy_peaks) is 0:
                        # this happens when all the peaks are -1, so either if peaks are out of detection range or halo velocity is not compatible with galactic center
                        # if velocity is in the range, print out that peaks are not in the detection range, otherwise print out all -1
                        if (self.maddm_card['vave_indirect_line'] > self.vave_indirect_line_range[0] and self.maddm_card['vave_indirect_line'] < self.vave_indirect_line_range[1]):
                            out.write('# No peaks found: out of detection range.\n')
                            continue
                        else: # recompute energy_peaks without assuming values != -1
                            energy_peaks = collections.OrderedDict(sorted([(k, v) for k, v in six.iteritems(self.last_results) if str_part_peak in k and '_states' not in k and '_error' not in k], key = lambda item: item[1])) # key = "line_<exp_name>_peak_<num>", value = energy peak
                    peaks_string  = ''
                    fluxes_string = ''
                    like_string = ''
                    for k, peak in six.iteritems(energy_peaks):
                        num        = int(k.split('_')[-1])
                        error_code = self.last_results[str_part + "peak_%d_error"  % num]
                        error_str  = "# error: %s" % error_code if error_code != 0 else ''
                        flux       = self.last_results[str_part + "flux_%d"        % num]
                        flux_UL    = self.last_results[str_part + "flux_UL_%d"     % num]
                        peaks_string  += form_s("peak_%d(%s)" % (num, self.last_results[k + "_states"])) + '= ' + form_n(peak) + ' \t ' + error_str + '\n'
                        fluxes_string += form_s("flux_%d" % num) + '= ['+ form_n(flux) + ',' + form_n(flux_UL) + ']\n'
                        try:
                            like       = self.last_results[str_part + "like_%d"        % num]
                            pvalue     = self.last_results[str_part + "pvalue_%d"      % num]
                            like_string   += form_s("-2*log(Likelihood)_%d" % num) + '= '+ form_n(like) + '\n' + form_s("p-value_%d" % num) + '= '+ form_n(pvalue) + '\n'
                        except KeyError:
                            continue
                    out.write(peaks_string)
                    out.write(fluxes_string)
                    out.write(like_string)


        if fluxes_source:
           out.write('\n#################################################\n')
           out.write('# CR Flux at Earth [particles cm^-2 s^-1 sr^-1) #\n')
           out.write('#################################################\n\n')
           out.write('# Fluxes calculated using the spectra from ' + self.maddm_card['indirect_flux_earth_method'] + '\n\n' )
           for name in ['neutrinos_e','neutrinos_mu','neutrinos_tau','gammas']:
                out.write(form_s('Flux_'+name)+ '= '+ form_s(form_n(self.last_results['flux_'+name] )) + '\n' )
#           if self.last_results['flux_positrons']:
#               out.write(form_s('Flux_positrons')+ '= '+ form_s(form_n(self.last_results['flux_positrons'] )) + '\n' )


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
        formatted = '{0:20}'.format(stringa)
        return  formatted
    def form_n(self,num):
        formatted = '{0:8.2e}'.format(num)
        return formatted
    
    def ask_run_configuration(self, mode=None, force=False):
        """ask the question about card edition / Run mode """
        process_data = self.proc_characteristics
        if not hasattr(self, 'mode') or not force:
            self.mode, cmd_quest = self.ask('', '0', mode=mode, 
                        data=self.proc_characteristics,
                        ask_class=MadDMSelector, timeout=60, path_msg=' ',
                        return_instance=True, force=force)
            # automatically switch to keep_wgt option
            #edit the maddm_card to be consistent with self.mode
            cmd_quest.get_cardcmd()
            # write Cards/.lastmode to recover 
            cmd_quest.write_switch()
    
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
        
        # if 'forbid_fast' is True, then we can't use 'inclusive' --> switch to 'reshuffling'
        # if self.proc_characteristics['forbid_fast'] and self.maddm_card['sigmav_method'] in ['inclusive']:
        #     logger.error("Can't use '%s' method, because it is forbidden in this case. Switched automatically to 'reshuffling'." % self.maddm_card['sigmav_method'])
        #     cmd_quest.setDM('sigmav_method','reshuffling',loglevel=30)

        def check_mode_and_dirs(mode, one_dir_exists):
            ''' check:
                if MODE 'OFF': OK
                if MODE 'ON' : check if at least one dirs for that MODE exists: if so: OK, if not: PROBLEM!
                we need to check each mode with this method and see if both independently are OK.
                Indeed: if 2to2lo is not set by sigmav_method being inclusive, then we need to be sure the
                user has not used the option during generation, while then setting a different sigmav_method from inclusive
                So the condition is that if the MODE is OFF then no problem at all, but if it's ON, we need to check
                that there is at least one directory for that MODE (tree or LI doesn't care), otherwise that would
                mean they have not been created as a consequence for the choice of the user.
                mode and one_dir_exists are both bool.
            '''
            if not mode:
                return True
            elif one_dir_exists:
                return True
            else:
                return False

        # set self._two2twoLO if we need to use inclusive method
        if self.maddm_card['sigmav_method'] in ['inclusive']:
            self._two2twoLO = True
        # if the 'if' above is not executed then we have that sigmav_method can be either 'reshuffling' or 'madevent'
        # so we should check if there is at least one directory for 'cont' or 'line' if respectively the modes 'indirect' or 'spectral' are ON
        # in this case we can go ahead (see the __doc__ of the function above), otherwise we can't use Madevent
        elif check_mode_and_dirs(self.mode['indirect'], any(exist for dir_, exist in self.indirect_directories.items() if 'cont' in dir_)) and\
             check_mode_and_dirs(self.mode['spectral'], any(exist for dir_, exist in self.indirect_directories.items() if 'line' in dir_)):
             self._two2twoLO = False
        else:
            misc.sprint(os.listdir(self.dir_path))
            raise Exception('Madevent output for indirect detection not available')
        
        logger.debug('2to2: %s' % self._two2twoLO)

        #set fortran switch and write include file
        if not self.in_scan_mode:
            # create the inc file for maddm
            self.maddm_card.set('do_relic_density', self.mode['relic'], user=False)
            self.maddm_card.set('do_direct_detection', True if self.mode['direct'] else False, user=False)
            self.maddm_card.set('do_directional_detection', self.mode['direct'] == 'directional', user=False)
            self.maddm_card.set('do_capture', self.mode['capture'], user=False)
            self.maddm_card.set('do_indirect_detection', True if self.mode['indirect'] else False, user=False)
            self.maddm_card.set('do_indirect_spectral', self.mode['spectral'], user=False)
            self.maddm_card.set('do_flux', True if (self.mode['indirect'] and self.mode['indirect'] != 'sigmav') else False, user=False)
            self.maddm_card.set('only2to2lo', self._two2twoLO, user=False)
            #self.maddm_card.set('do_indirect_spectral', self.mode['spectral'], user=False)
            #self.maddm_card.set('run_multinest', self.mode['run_multinest'], user=False)

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
            raise Exception('Compilation of maddm failed')
        
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
    def __init__(self, *args, **kwargs):
        self.mode_vincia = kwargs.pop("mode_vincia", False)
        super(Indirect_PY8Card, self).__init__(*args, **kwargs)
    
    def default_setup(self):
        """ Sets up the list of available PY8 parameters."""
        
        self.add_param("Main:timesAllowErrors", 10, hidden=True, comment="allow a few failures before quitting")
        self.add_param("Main:NumberOfEvents", -1,  comment="number of events to go trough")
        self.add_param("Main:mDM", 1000.0, hidden=True, comment=" mass of the Dark matter")
        self.add_param("Main:outdir" , './', hidden=True, comment="specify output dir") 
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
        self.add_param("PartonShowers:model", 1)
        self.add_param("TimeShower:weakShower", False, comment="Run weak-shower for FSR")
        self.add_param("TimeShower:weakShowerMode", 0, comment="Determine which branchings are allowed (0 -> W and Z)")
        self.add_param("TimeShower:pTminWeak", 0.1)
        self.add_param("WeakShower:vetoWeakJets", True)
        self.add_param("WeakShower:singleEmission", False)
        self.add_param("WeakShower:enhancement",1.0,comment="enhanced weak shower rate")
        # decay all particles
        self.add_param("13:mayDecay", True, hidden=True, comment="decays muons")
        self.add_param("211:mayDecay", True, hidden=True, comment="decays pions")
        self.add_param("321:mayDecay", True, hidden=True, comment="decays Kaons")
        self.add_param("130:mayDecay", True, hidden=True, comment="decays K_LO")
        self.add_param("310:mayDecay", True, hidden=True, comment="decays K_SO")
        self.add_param("2112:mayDecay", True, hidden=True, comment="decays neutron")
        # disabling some events checks for bug below 100 GeV
        self.add_param("Check:event", False, hidden=True, comment="avoid pythia crushing")
        self.add_param("Check:beams", False, hidden=True, comment="avoid pythia crushing")
        # Tuning for the following parameters taken from Tab 2 of https://arxiv.org/pdf/2303.11363.pdf
        self.add_param("StringZ:deriveBLund", True)
        self.add_param("StringZ:aLund", 0.601)
        self.add_param("StringZ:bLund", 0.897)
        self.add_param("StringZ:avgZLund", 0.540)
        self.add_param("StringPT:Sigma", 0.307)
        self.add_param("StringZ:aExtraDiquark", 1.671)
        # Coalescence parameters
        self.add_param("Main:p_coalescence", 0.180)
        self.add_param("Main:d_coalescence", 2.0)
        self.add_param("Main:sigma_coalescence", 1.0)
        # Vincia?
        if self.mode_vincia:
            self.add_param("Vincia:ewMode", 3, comment="Vincia EW mode")

        # self.add_param("PartonShowers:model", 2)
        # self.add_param("Vincia:ewMode", 3, comment="Vincia EW mode")
        # self.add_param("StringZ:deriveBLund", False)
        # self.add_param("StringZ:aLund", 0.337409)
        # self.add_param("StringZ:bLund", 0.784682)
        # self.add_param("StringZ:avgZLund", 0.55)
        # self.add_param("StringPT:Sigma", 0.296569)
        # self.add_param("StringZ:aExtraDiquark", 1.246986)

class MadDMSelector(cmd.ControlSwitch, common_run.AskforEditCard):
    """ """

    to_control= [('relic', 'Compute the Relic Density'),
                 ('direct', 'Compute direct(ional) detection'),
                 ('indirect', 'Compute indirect detection/flux (cont spectrum)'),
                 ('spectral', 'Compute indirect detection in aX (line spectrum)'),
                 ('nestscan', 'Run Multinest scan'),
                ]
    to_init_card = ['param', 'maddm', 'pythia8']
    PY8Card_class = Indirect_PY8Card
    
    integer_bias = len(to_control) + 1 # integer corresponding to the first entry in self.cards
    
    ####################################################################
    # everything related to relic option
    ####################################################################    
    def set_default_relic(self):
        """set the default value for relic=
           if relic has been generated when calling indirect detection, then it is set as 'OFF' by default.
           the 'Not Avail.' case happens when relic has not been generated neither explicitly nor through indirect detection"""
        
        if self.availmode['relic_density_off']: # this can be True only if self.availmode['has_relic_density'] is True as well, that would correspond to generate relic_density during ID
            self.switch['relic'] = 'OFF'
        elif self.availmode['has_relic_density']: # otherwise that would mean that relic density has/hasn't been asked explicitly
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
        """set the default value for direct="""
        
        if self.availmode['has_directional_detection']:
            self.switch['direct'] = 'directional'
        elif self.availmode['has_direct_detection']:
            self.switch['direct'] = 'direct'        
        else:
            self.switch['direct'] = 'Not Avail.'

    def get_allowed_direct(self):
        """Specify which parameter are allowed for direct="""
        
        
        if hasattr(self, 'allowed_direct'):
            return getattr(self, 'allowed_direct')

        if self.availmode['has_directional_detection']:
            self.allowed_direct =  ['directional', 'direct','OFF']
        elif self.availmode['has_direct_detection']:
            self.allowed_direct =  ['direct','OFF']
        else:
            return []

    def check_value_direct(self, value):
        """ allow direct=ON in top of standard mode """
        
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
        """set the default value for indirect="""
        
        if not HAS_NUMPY:
            self.switch['indirect'] = 'Not Avail. (numpy missing)'
        elif self.availmode['has_indirect_detection']:
            self.switch['indirect'] = 'sigmav'     
        else:
            self.switch['indirect'] = 'Not Avail.'

    def print_options_indirect(self):
        """print statement for the options"""
        
        if not self.availmode['has_indirect_detection']:
            return "Please install module"
        elif not HAS_NUMPY:
            return "numpy not available"
        else:
            return self.print_options('indirect', keep_default=True)

    def get_allowed_indirect(self):
        """Specify which parameter are allowed for indirect="""
        
        
        if hasattr(self, 'allowed_indirect'):
            return getattr(self, 'allowed_indirect')

        if not HAS_NUMPY:
            self.allowed_indirect =  ['OFF']
        elif self.availmode['has_indirect_detection']:
            self.allowed_indirect =  ['OFF', 'sigmav', 'flux_source', 'flux_earth']
        else:
            self.allowed_indirect = []
        return self.allowed_indirect

    
    def check_value_indirect(self, value):
        """ allow indirect=ON in top of standard mode """
     
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
    # everything related to spectral option
    ####################################################################    
    def set_default_spectral(self):
        """set the default value for spectral="""
        
        if not HAS_NUMPY:
            self.switch['spectral'] = 'Not Avail. (numpy missing)'
        elif not HAS_SCIPY:
            self.switch['spectral'] = 'Not Avail. (scipy missing)'
        elif self.availmode['has_indirect_spectral']:
            self.switch['spectral'] = 'ON'     
        else:
            self.switch['spectral'] = 'Not Avail.'

    def print_options_spectral(self):
        """print statement for the options"""
        if not self.availmode['has_indirect_spectral']:
            return "Please install module"
        elif not HAS_NUMPY:
            return "numpy not available"
        elif not HAS_SCIPY:
            return "scipy not available"
        else:
            return self.print_options('spectral', keep_default=True)

    def get_allowed_spectral(self):
        """Specify which parameter are allowed for spectral="""
        if hasattr(self, 'allowed_indirect_spectral'):
            return getattr(self, 'allowed_indirect_spectral')

        if not HAS_NUMPY:
            self.allowed_indirect_spectral =  ['OFF']
        if not HAS_SCIPY:
            self.allowed_indirect_spectral =  ['OFF']
        elif self.availmode['has_indirect_spectral']:
            self.allowed_indirect_spectral =  ['OFF', 'ON']
        else:
            self.allowed_indirect_spectral = []
        return self.allowed_indirect_spectral
                
    ####################################################################
    # everything related to multinest option
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
        
        # if os.path.exists(pjoin(self.me_dir, 'Cards', '.lastrun')):
        #     self.default_switch = {}
        #     for line in  open(pjoin(self.me_dir, 'Cards', '.lastrun')):
        #         if not line: continue
        #         key, value = line.split()
        #         self.default_switch[key] = value
        

        cmd.ControlSwitch.__init__(self, self.to_control, opts['mother_interface'], *args, **opts)
        # initialise the various card to control
        question = ''
        
        param_card_path = pjoin(opts['mother_interface'].dir_path, 'Cards', 'param_card.dat')
        pythia8_card_path = pjoin(opts['mother_interface'].dir_path, 'Cards', 'pythia8_card.dat')

        cards = [param_card_path, 'maddm', 'multinest', pythia8_card_path]
        
        tmp_allow_arg = list(self.allow_arg)
        opts = dict(opts)
        opts['allow_arg'] = []
        
        common_run.AskforEditCard.__init__(self, question, cards,
                                            *args, **opts )

        self.allow_arg += tmp_allow_arg
        self.question = self.create_question()

               
    def write_switch(self, path=None):
        """store value of the switch for the default at next run"""
        
        if not path:
            path = pjoin(self.me_dir, 'Cards', '.lastrun')
        fsock = open(path,'w')
        for key, value in self.answer.items():
            fsock.write('%s %s\n' % (key,value))
        fsock.close()
           
           
    def detect_card_type(self, path):
        """detect card type -> add multinest"""
        
        output = super(MadDMSelector, self).detect_card_type(path)
        
        if output == 'unknown':
            text = open(path).read()
            if 'livepoints' in text: return 'multinest'
            if 'prefix' in text: return 'multinest'
            
        return output
     
    def init_maddm(self, path):
        """ initialize cards for the reading/writing of maddm"""
        
        self.maddm_def = MadDMCard(self.paths['maddm_default'], consistency=False)
        try:
            self.maddm = MadDMCard(self.paths['maddm'], consistency=False)
        except Exception as e:
            logger.error('Current maddm_card is not valid. We are going to use the default one.')
            logger.error('problem detected: %s' % e) 
            files.cp(self.paths['maddm_default'], self.paths['maddm'])
            self.maddm = MadDMCard(self.paths['maddm'])

        self.already_warned = set([])
        self.old_vave_indirect = self.maddm['vave_indirect']
        
        self.special_shortcut.update({'nevents': ([int],['Main:numberOfEvents %(0)s']),
                                      'fast':([],[lambda self: self.pass_to_fast_mode]),
                                      'precise':([],[lambda self: self.pass_to_precise_mode])
                                      })
        self.special_shortcut_help.update({'nevents': 'number of events to generate for indirect detection if sigmav_method is madevent/reshuffling',
                                           'fast': 'modify the maddm_card to favor fast computation (but less precise) of the indirect mode',
                                           'precise': 'modify the maddm_card to favor accuracy over speed computation for the indirect mode'
                                          }
                                          )
    
        self.maddm_set = list(set(list(self.maddm_def.keys()) + self.maddm_def.hidden_param))
        return list(self.maddm.keys()) 

    
    def pass_to_fast_mode(self):
        """pass to fast mode according to the paper"""
        
        if self.availmode['forbid_fast']:
            logger.error("setting fast mode is only valid when loop-induced processes are not considered.")
            return
        indirect = self.answer['indirect']
        if indirect in ['OFF', None]:
            logger.error("setting fast mode is only valid when indirect mode is getting called.")
            return 
        self.do_set("sigmav_method inclusive")
        self.do_set("indirect_flux_source_method CosmiXs")
        self.do_set("indirect_flux_earth_method PPPC4DMID_ep")
        
    def pass_to_precise_mode(self):
        """pass to precise mode according to the paper"""
        
        indirect = self.answer['indirect']
        spectral = self.answer['spectral']
        if indirect in ['OFF', None] and spectral in ['OFF', None]:
            logger.error("setting precise mode is only valid when indirect or spectral mode is getting called.")
            return 
        self.do_set("sigmav_method reshuffling")
        self.do_set("indirect_flux_source_method pythia8")
        self.do_set("indirect_flux_earth_method dragon")
        self.do_set("Main:numberOfEvents 1000000")
        self.do_set("TimeShower:weakShower = on")

    # def toggle_vincia(self):
    #     """Toggle Vincia"""
    #     
    #     indirect = self.answer['indirect']
    #     if 'flux' not in indirect:
    #         logger.error("setting Vincia is only valid when computation of the spectra is activated")
    #         return
    #
    #     old_value = self.maddm['indirect_flux_source_method']
    #     if old_value == "pythia8":
    #         self.do_set("indirect_flux_source_method vincia")
    #         self.PY8Card_class = Indirect_VinciaCard
    #     else:
    #         self.do_set("indirect_flux_source_method pythia8")
    #         self.PY8Card_class = Indirect_PY8Card
    #     self.do_set("TimeShower:weakShower = on")
    #
    #     pythia8_card_path = pjoin(self.me_dir, 'Cards', 'pythia8_card.dat')
    #     default_card_path = pjoin(self.me_dir, 'Cards', self.maddm['indirect_flux_source_method'] + '_card_default.dat')
    #     shutil.copy(default_card_path, pythia8_card_path)
    #     logger.info("Set new pythia8_card.dat on " + self.maddm['indirect_flux_source_method'] + " template", "$MG:BOLD")

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
        #    8/9 trigger via self.trigger_8 and self.trigger_9
        #    If you change the numbering, please change the name of the 
        #    trigger function accordingly.
        
        question = cmd.ControlSwitch.create_question(self, help_text=False)
        question +="""\n%(start_green)s You can also edit the various input card%(stop)s:
 * Enter the name/number to open the editor
 * Enter a path to a file to replace the card
 * Enter %(start_bold)sset NAME value%(stop)s to change any parameter to the requested value
 /=============================================================================\ 
 |  6. Edit the model parameters    [%(start_underline)sparam%(stop)s]                                    |  
 |  7. Edit the MadDM options       [%(start_underline)smaddm%(stop)s]                                    |
"""

        current_val  = self.answer # use that to be secure with conflict -> always propose card
        if current_val['nestscan'] == "ON" or self.switch["nestscan"] ==  "ON":
            question += """ |  8. Edit the Multinest options  [%(start_underline)smultinest%(stop)s]                                 |\n"""
    
        if current_val['indirect'].startswith('flux') or self.switch["indirect"].startswith('flux'):
            question += """ |  9. Edit the Showering Card for flux  [%(start_underline)sflux%(stop)s]                                |\n"""
        
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
    

    def default(self, line):
        """Default action if line is not recognized"""
        
        try:
            return cmd.ControlSwitch.default(self, line, raise_error=True)
        except cmd.NotValidInput as error:
            return common_run.AskforEditCard.default(self, line)     
        
    def trigger_8(self, line):
        """ trigger function for default function:
            allows to modify the line/ trigger action.
            
            check that nestscan can be ON
               if it is, set it ON,
               otherwise, print a warning and set the line to reask the 
               question (return False)
        """
        
        #  1) do like the user type "multinest=ON"
        #  2) go to the line edition
        self.set_switch('nestscan', "ON", user=True) 
        if self.switch['nestscan'] == "ON":
            return '8 %s' % line
        # not valid nestscan - > reask question
        else:
            return 

    def trigger_9(self, line):
        """ trigger function for default function:
            allows to modify the line/ trigger action.
            
            check that indirect can be ON
               if it is, set it ON,
               otherwise, print a warning and set the line to reask the 
               question (return False)
        """
        
        #  1) do like the user type "indirect=flux_source" (if not in flux)
        #  2) forbid sigmav_method = inclusive
        #  3) forbid PPPC4DMID
        #  4) go to the line edition
        #1.
        if not self.switch['indirect'].startswith('flux_'):
            self.set_switch('indirect', "flux_source", user=True) 
            if not self.switch['indirect'] == 'flux_source':
                # not valid indirect - > reask question
                return 
            logger.warning('switching indirect mode to flux_source since sigmav does not support flux')
           
        #2. check that sigmav_method is not inclusive
        if self.maddm['sigmav_method'] == 'inclusive':
            self.setDM('sigmav_method', 'reshuffling',loglevel=30)
            
        #3. ensure pythia8 is on
        if self.maddm['indirect_flux_source_method'] != 'pythia8' and self.maddm['indirect_flux_source_method'] != 'vincia':
            self.setDM('indirect_flux_source_method', 'pythia8',loglevel=30)
        
        return '9 %s' % line
    
    trigger_flux = trigger_9
    
    def do_compute_widths(self, line):
        """normal fct but ensure that self.maddm_card is up-to-date"""
        
        try:
            self.mother_interface.maddm_card = self.maddm
        except Exception as error:
            logger.error("Invalid command: %s " % error)
            return
        return super(MadDMSelector, self).do_compute_widths(line)
            
    def do_help(self, line, conflict_raise=False, banner=True):
        """proxy for do_help"""
             
        if line:
            if banner:                      
                logger.info('*** HELP MESSAGE ***', '$MG:BOLD')
            card = common_run.AskforEditCard.do_help(self, line, conflict_raise=conflict_raise, banner=False)
        else:
            if banner:
                logger.info('*** HELP MESSAGE FOR CARD EDITION ***', '$MG:BOLD')
            card = common_run.AskforEditCard.do_help(self, line, conflict_raise=conflict_raise, banner=False)
            logger.info('*** HELP MESSAGE FOR CHANGING SWITCH ***', '$MG:BOLD')
            card = cmd.ControlSwitch.do_help(self, line, list_command=False)
            if banner:                      
                logger.info('*** END HELP ***', '$MG:BOLD')
            
            logger_tuto.info("""
This question allows you BOTH to define what you are going to run (via the value at the top).
But also to edit the content of the various file defining the  run/benchmark (bottom).

To bypass the computation of relic density you can do
> relic=OFF            
to make the computation of the directional detection
> direct=directional
                
You can also edit the card referenced.
Note that you can 
   1) edit any parameter like this:
       > set mxd 10
       [use auto completion if you need to search a name]
   2) run a scan over parameter space
        > set mxd scan:[10, 20, 40]
        > set mxd scan:[10**i for i in range(5)]
   3) ask to compute the width automatically for a particle
        > set wy0 Auto  
        
When you are done with such edition, just press enter (or write 'done' or '0')          
""") 

            return
            
        args = self.split_arg(line)
        
        if not args:
            args =['']
        
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
                logger.info('**   AMBIGUOUS NAME: %s **', args[start], '$MG:BOLD')
                if card == '':
                    logger.info('**   If not explicitely speficy this parameter  will modif the maddm_card file', '$MG:BOLD')

            self.maddm.do_help(args[start])
            
        ### CHECK if a help exist (only for those wihtout a do_xxx)
        if len(args) == 1:
            if not hasattr(self, 'do_%s' % args[0]) and hasattr(self, 'help_%s' % args[0]):
                getattr(self, 'help_%s' %args[0])()
            
        if banner:                      
            logger.info('*** END HELP ***', '$MG:BOLD')  
        return card    

    def help_spectral(self):
        logger.info("spectral flag can take two values: ON/OFF")
        logger.info("  It controls if you are going to compute the indirect detection into the final state(s) a X (X = a, h, z)")

    def help_indirect(self):
        
        logger.info("indirect flag can take three values: sigmav, flux_source, flux_earth")
        logger.info("   for each of those flag, we propose two mode of computation fast or precise")
        logger.info("   this can be setup (and tweak) within the maddm_card.") 
        logger.info("    A quick way to setup such mode is to type: 'set fast' or 'set precise'")
        logger.info()
        logger.info("  sigmav: ", "$MG:BOLD")
        logger.info("     computation of the velocity-weighted annihilation cross section at present time.")
        logger.info('     ')
        logger.info("  flux_source", "$MG:BOLD")
        logger.info("")
        logger.info("     computation of the energy spectra before any propagation.")
        logger.info("")
        logger.info("  flux_earth:", "$MG:BOLD")
        logger.info("     This produces the flux at Earth of e+ and anti proton")
                
    def help_direct(self):
        
        logger.info("direct flag can take three values: direct, directional")
        logger.info('     ')
        logger.info("  direct: ", "$MG:BOLD")
        logger.info("     Theoretical elastic spin-independent and spin-dependent cross section dark matter off nucleons")
        logger.info('     ')
        logger.info("  directional", "$MG:BOLD")
        logger.info("     Directional event rate (double differential event rate)")
        
    def help_relic(self):
        logger.info("relic flag can take two values: ON/OFF")
        logger.info("  It controls if you are going to compute the relic density for the current benchmark")
      
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
        self.maddm.full_template = True
        self.maddm.write(self.paths['maddm'], write_hidden=True)


    def postcmd(self, stop, line):
        """temporary to support 2.6.3"""
        
        stop = super(MadDMSelector, self).postcmd(stop, line)
        curr_version = misc.get_pkg_info()['version'] 
        if misc.get_older_version(curr_version, '2.6.4') != '2.6.4':
            common_run.AskforEditCard.postcmd(self, stop, line)
        return stop
    
    def check_card_consistency(self):

        
        super(MadDMSelector, self).check_card_consistency()

        #if there are any new jfactors, make sure to write them in
        #logger.debug('Updating the Jfactor file')
        #self.write_jfactors()

        # check if user uses old variable 'vave_indirect':
        # check if in user_set: warn the user and put it in self.maddm.already_warned set
        # however, it could happen the user uses the old name again: check if the new value is different from the older one,
        # if this is the case, then warn the user again.
        if ('vave_indirect' in self.maddm.user_set) and ('vave_indirect' not in self.already_warned) and (self.maddm['vave_indirect'] != self.old_vave_indirect):
            logger.warning("You are using an old parameter name 'vave_indirect'. Fall back to the updated name 'vave_indirect_cont'.")
            self.already_warned.add('vave_indirect')
            self.old_vave_indirect = self.maddm['vave_indirect']
            self.setDM('vave_indirect_cont', self.maddm.get('vave_indirect'), loglevel = 30)

        # If direct detection is ON ensure that quark mass are not zero
        if self.switch['direct'] != 'OFF':
            to_change = []
            for i in range(1,7):
                if self.param_card.get_value('mass', i) == 0.:
                    to_change.append(i)
            if to_change:
                logger.warning('For direct detection the quark mass need to be different of zero. Automatically adding such masses (PDG 2014).')
                quark_masses = {1: 4.8e-3, 2: 2.3e-3,  3: 95e-3, 4: 1.275, 5: 4.18, 6: 173.21}
                for i in to_change:
                    self.do_set('param_card mass %s %s' % (i, quark_masses[i]))
                    
        def_key = [k.lhacode for k in self.param_card['mass']]
        for key in self.param_card_default['mass']:
            if key.lhacode not in def_key:
                to_add = check_param_card.Parameter(block='mass', lhacode=key.lhacode, value=1e10, comment='for DD')
                self.param_card['mass'].append(to_add)
                self.do_set('param_card mass %s %s' % (key.lhacode, 1e10))
                self.do_set('param_card decay %s %s' % (key.lhacode, 0))
                self.param_card.write(self.paths['param'])

        # check consistency of sigmav_method, on the basis of forbid_fast option
        if self.availmode['forbid_fast'] and self.maddm['sigmav_method'] == 'inclusive':
            logger.warning("setting 'sigmav_method' to 'inclusive' is only valid when loop-induced processes are not considered. Switched to 'reshuffling'.")
            self.setDM('sigmav_method','reshuffling',loglevel=30)

        # check whether the user selects PPPC4DMID
        if self.maddm['indirect_flux_source_method'].startswith("PPPC"):
            answer = ''
            while(answer.lower() not in ['y', 'n']):
                answer = input(bcolors.WARNING + "You selected " + self.maddm['indirect_flux_source_method'] + " as the method to compute indirect detection spectra.\nHowever, MadDM ships an improved version of the spectra, that can be selected with 'CosmiXs', do you want to use these? (Y/n) " + bcolors.ENDC).lower()
                if answer == "" or answer == "y":
                    logger.warning("Switching to 'CosmiXs'")
                    self.setDM('indirect_flux_source_method','CosmiXs',loglevel=30)
                    break
                if answer == "n":
                    logger.warning("Keep 'indirect_flux_source_method' to " + self.maddm['indirect_flux_source_method'])
                    break


        # check pythia card consistency
        card_parameters = {}
        if self.maddm['indirect_flux_source_method'] == 'pythia8':
            card_parameters = {
                "PartonShowers:model": 1,
                "StringZ:deriveBLund": True,
                "StringZ:aLund": 0.601,
                "StringZ:bLund": 0.897,
                "StringZ:avgZLund": 0.540,
                "StringPT:Sigma": 0.307,
                "StringZ:aExtraDiquark": 1.671
            }
        elif self.maddm['indirect_flux_source_method'] == 'vincia':
            card_parameters = {
                "PartonShowers:model": 2,
                "StringZ:deriveBLund": False,
                "StringZ:aLund": 0.337409,
                "StringZ:bLund": 0.784682,
                "StringZ:avgZLund": 0.55,
                "StringPT:Sigma": 0.296569,
                "StringZ:aExtraDiquark": 1.246986
            }
        for par, value in card_parameters.items():
            self.do_set('pythia8_card %s %s' % (par, value))
        
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
                self.limits._sigma_SI = secure_float_f77(args[1])
                self.limits._sigma_SI_width = secure_float_f77(args[2])
            else:
                logger.info('The file you use for direct detection limits will be interpreted as:')
                logger.info('Column 1 - dark matter mass in GeV')
                logger.info('Column 2 - upper limit on direct detection cross section in pb')
                self.limits._dd_limit_file['si'] = args[1]
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
                    self.limits._dd_limit_file['sd_proton'] = args[2]
                elif args[2] == 'n':
                    self.limits._dd_limit_file['sd_neutron'] = args[2]
                else:
                    logger.error('set dd_sd_limits needs to be formatted as:')
                    logger.error('set dd_sd_limits <p or n> <observed val> <uncertainty>, or')
                    logger.error('set dd_sd_limits <p or n> <limit filename>, or')

                self.ExpConstraints.load_constraints()

        elif args[0] == 'relic_limits':
            logger.info('The info is interpreted as: <Oh^2>, CL width')
            self.limits._oh2_planck = secure_float_f77(args[1])
            self.limits._oh2_planck_width = secure_float_f77(args[2])

        elif args[0] == 'id_limits':
            if len(args)!= 4:
                logger.warning('You need to provide the following <ave. velocity> <ann. channel> <file path>')
                logger.warning('or  <ave. velocity> <ann. channel> <obs. cross section> <obs. uncertainty>')
                logger.warning('Annihilation channel can be (in pdg numbers): '+ str(self.limits._settable))
            logger.info('The file you use for indirect detection limits will be interpreted as:')
            logger.info('Column 1 - dark matter mass in GeV')
            logger.info('Column 2 - upper limit on the total annihilation cross section in cm^3/s at specified average velocity')

            if len(args) > 4:
                if args[2] in self.limits._settable:
                    self.limits._id_limit_vel[args[2]] = secure_float_f77(args[1])
                    self.limits._sigma_ID[args[2]] = secure_float_f77(args[3])
                    self.limits._sigma_ID_width[args[2]] = secure_float_f77(args[4])
                else:
                    logger.error('Final state not allowed for ID limits!')
            else:
                vel = secure_float_f77(args[1])
                channel = args[2]
                id_file = args[3]
                self.limits._id_limit_vel[channel] = vel
                self.limits._id_limit_file[channel] = id_file

                self.limits.load_constraints()
        elif args[0].lower() in self.switch:
            return self.default(line.strip())
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
        logger.log(loglevel,'modify parameter %s of the maddm_card.dat to %s' % (name, value), '$MG:BOLD')
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
    
    def __new__(cls, finput=None, **opt):
        """Bypass the standard RunCard one"""
        return super(banner_mod.RunCard, cls).__new__(cls, finput, **opt)

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
                dist = secure_float_f77(spline[1])
                val = secure_float_f77(spline[2].rstrip())
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
        self.add_param('do_indirect_spectral', False, system=True)
        self.add_param('only2to2lo', False, system=True)
        self.add_param('run_multinest', False, system=True, include=False)

        self.add_param('calc_taacs_ann_array', True)
        self.add_param('calc_taacs_dm2dm_array', True)
        self.add_param('calc_taacs_scattering_array', True)
        
        self.add_param('eps_ode', 0.01)
        self.add_param('xd_approx', False)
        self.add_param('running_as', True, comment='choose whether to run the strong coupling up to 2*m_chi energy scale', include=True)
        
        self.add_param('x_start', 50.0, hidden=True)
        self.add_param('x_end', 1000.0, hidden=True)
        self.add_param('dx_step', 1.0, hidden=True)
        
        self.add_param('ngrid_init', 50, hidden=True, comment="Initial number of points in the grid to integrate over (of dm velocity)")
        self.add_param('nres_points', 100, hidden=True, comment="Number of points to add for one width around each resonance peak. (the code will add points which exponentially increase in distance from the pole)")
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
        self.add_param('material', 1, allowed = list(range(1,14)),
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

        # velocities for indirect detection
        self.add_param('vave_indirect_cont', 0.00002, include=True, fortran_name='vave_indirect')
        self.add_param('vave_indirect', -1.0, include=False, hidden=True)
        self.add_param('vave_indirect_line', 0.00075, hidden = False, include = True)

        #indirect detection
        self.add_param('halo_profile', 'Draco', include=True, legacy=True)   ##### this naming doesn't make sense fix it!!!!!!!!

        self.add_param('jfactors', {'__type__':1.0}, include=False, hidden=True, system=True)
        self.add_param('distances', {'__type__':1.0}, include=False, hidden=True, system=True)
        self.add_param('npts_for_flux', 200, include=False, hidden=True) #number of points for the flux diff. distribution

        self.add_param('only_two_body_decays', True, include=False)
        #self.add_param('num_of_elements', 60)
        
        self.fill_jfactors()

        self.add_param('indirect_flux_source_method', 'pythia8', comment='choose between pythia8, vincia, PPPC4DMID and PPPC4DMID_ew', include=False,
                       allowed=['pythia8','vincia','PPPC4DMID','PPPC4DMID_ew','CosmiXs'])
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

        # line analysis parameters
        self.add_param('profile', 'einasto', comment='halo density profile for gamma-line searches, choose between gnfw, einasto, nfw, isothermal, burkert', include = False, \
                           hidden = False, allowed = ['gnfw', 'einasto', 'nfw', 'isothermal', 'burkert'])
        self.add_param('r_s', 20.0, comment='scale radius (in kpc)', include = False, \
                           hidden = False)
        self.add_param('gamma', 1.3, comment='gamma parameter, relevant for gnfw profile', include = False, \
                           hidden = False)
        self.add_param('alpha', 0.17, comment='alpha parameter, relevant for einasto profile', include = False, \
                           hidden = False)
        self.add_param('r_sun', 8.5, comment='distance of the galactic center from the Sun', include = False, \
                           hidden = True)
        self.add_param('rho_sun', 0.4, comment='dark matter density at the position of the Sun', include = False, \
                           hidden = True)
        self.add_param('roi_fermi_2015', 'default', comment='ROI relevant for the Fermi-LAT 2015 results, set default to pick the ROI optimized for the chosen profile', include = False, \
                           hidden = False, allowed = ['default', 'r3', 'r16', 'r41', 'r90'])
        self.add_param('n_fwhm_separation', 5., comment='number of FWHM to take as the minimum energy separation between peaks to consider them well separated', include = False, \
                           hidden = True)
        self.add_param('peak_height_factor', 10., comment="minimum ratio between peaks' heights to have one of them significantly higher than the other", include = False, \
                           hidden = True)
        self.add_param('toggle_fermi_2015', 'on', comment="toggle experiment Fermi-LAT 2015 related to gamma-line searches (arxiv: 1506.00013)", include = False, \
                           hidden = True, allowed = ['on', 'off'])
        self.add_param('toggle_hess_2018', 'on', comment="toggle experiment HESS 2018 related to gamma-line searches (arxiv: 1805.05741)", include = False, \
                           hidden = True, allowed = ['on', 'off'])
        # template experiment
        self.add_param('toggle_template_line_experiment', 'off', comment="toggle experiment HESS 2018 related to gamma-line searches (arxiv: 1805.05741)", include = False, \
                           hidden = True, allowed = ['on', 'off'])
        self.add_param('template_line_experiment_name', 'TEMPLATE_name', comment="name of the template experiment for gamma-line searches", include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_arxiv', 'TEMPLATE_arxiv', comment="arxiv number of the template experiment gamma-line searches", include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_roi', 1., comment="default ROI of the template experiment gamma-line searches", include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_jfactor_strategy', 'normal', comment='strategy for J-factor computation', include = False, \
                           hidden = True, allowed = ['normal', 'inverted'])
        self.add_param('template_line_experiment_profile', 'einasto', comment='default halo density profile for the template experiment gamma-line searches, choose between gnfw, einasto, nfw, isothermal, burkert', include = False, \
                           hidden = True, allowed = ['gnfw', 'einasto', 'nfw', 'isothermal', 'burkert'])
        self.add_param('template_line_experiment_r_max', 'inf', comment='maximum value of line of sight (in kpc) to integrate over', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_r_s', 20.0, comment='scale radius (in kpc) of the default profile related to the template line experiment gamma-line searches', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_gamma', 1.3, comment='gamma parameter, relevant for gnfw profile related to the template line experiment gamma-line searches', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_alpha', 0.17, comment='alpha parameter, relevant for einasto profile related to the template line experiment gamma-line searches', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_mask_latitude', 0., comment='angle beta of the mask (in deg): observing the galactic centre from the position of the Sun, mask P, if abs(latitude(P)) < beta (refer to MadDM 3.2 documentation)', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_mask_longitude', 180., comment='angle lambda of the mask (in deg): observing the galactic centre from the position of the Sun, mask P, if abs(longitude(P)) > lambda (refer to MadDM 3.2 documentation)', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_mask_inner_angle', 0., comment='angle alpha_1 of the mask (in deg): observing the galactic centre from the position of the Sun, mask P, if abs(arctan(P.y/P.x)) < alpha_1/2 (refer to MadDM 3.2 documentation)', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_energy_resolution', 10., comment='percent value to be multiplied by the energy of the peak to find the resolution at that energy', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_detection_range_min', 0., comment='detection range minimum (in GeV)', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_detection_range_max', 1.e+6, comment='detection range maximum (in GeV)', include = False, \
                           hidden = True)
        self.add_param('template_line_experiment_constraints_file', 'None', comment='file containing 3 columns [ DM mass (GeV), <sigmav> (cm^3 s^-1), flux (cm^-2 s^-1) ] related to the constraints on gamma-line searches for the template experiment; comments must be prepended with \'#\'; this file must be placed in $MADDM_PATH/ExpData/', include = False, \
                           hidden = True)

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
            raise InvalidMaddmCard('The sum of SP* parameter should be 1.0 get %s' % (self['SPu'] + self['SPs'] + self['SPd'] + self['SPg']))
        
        if self['SNu'] + self['SNs'] + self['SNd'] - self['SNg'] -1 > 1e-3:
            raise InvalidMaddmCard('The sum of SM* parameter should be 1.0 get %s' % (self['SNu'] + self['SNs'] + self['SNd'] + self['SNg']))
        
        if self['sigmav_method'] == 'inclusive':
            if self['indirect_flux_source_method'] == 'pythia8' or self['indirect_flux_source_method'] == 'vincia':
                if self['do_indirect_detection']:
                    logger.warning("since sigmav_method is on inclusive, indirect_flux_source_method has been switched to 'CosmiXs'")
                #if self['do_flux']:
                    #logger.warning('since sigmav_method is on inclusive, indirect_flux_source_method has been switch to PPPC4DMID')
                self['indirect_flux_source_method'] = 'CosmiXs'
            #if self['indirect_flux_earth_method'] != 'PPPC4DMID':
            #    if self['do_flux']:
            #        logger.warning('since sigmav_method is on inclusive, indirect_flux_earth_method has been switch to PPPC4DMID_ep')
            #    self['indirect_flux_earth_method'] = 'PPPC4DMID_ep' 
                
        elif self['indirect_flux_earth_method'] == 'PPPC4DMID_ep':
            logger.warning('since pythia8 is used to generate spectra at source, indirect_flux_earth_method has been switched to DRAGON')
            self['indirect_flux_earth_method'] = 'dragon' 
                
                
class Indirect_Cmd(me5_interface.MadEventCmdShell):
    
    options_configuration = dict(me5_interface.MadEventCmdShell.options_configuration)
    options_configuration.update({'pppc4dmid_path': './PPPC4DMID',
                                  'dragon_path': './DRAGON'})
    
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
    def create_root_file(self,*args, **opts):
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

    def check_nb_events(self):
        """Redefine maximum number of events for MadDM"""
        
        nb_event = int(self.run_card['nevents'])
        max_nb_events_allowed = int(1e8)
        if nb_event > max_nb_events_allowed:
            logger.warning("Attempting to generate more than 1M events")
            logger.warning("Limiting number to 1M. Use multi_run for larger statistics.")
            path = pjoin(self.me_dir, 'Cards', 'run_card.dat')
            os.system(r"""perl -p -i.bak -e "s/\d+\s*=\s*nevents/%d = nevents/" %s""" \
                                                                         % (max_nb_events_allowed, path))
            self.run_card['nevents'] = max_nb_events_allowed

        return
        


class Priors:
    priors = ['uniform', 'loguniform', 'user']

class Likelihoods:
    likelihoods = ['gaussian', 'half_gauss', 'user', 'off']
    observables = ['relic', 'directSI','directSD_p', 'directSD_n', 'indirect', 'spectral', 'capture']


class Multinest(object):

    def __init__(self, run_interface):

        self.options = {
            'prior':'loguniform',
            'loglikelihood':{'relic':'gaussian', 'directSI':'half_gauss', 'directSD_p':'half_gauss', 'directSD_n':'half_gauss', 'indirect':'', 'spectral':''},
            'livepts':10000,
            'sampling_efficiency':'model',
            'parameters':[],
            'prefix':'mnest_',
            'half_gauss_width':{'spinSI':0.01, 'spinSD_p':0.01, 'spinSD_n':0.01} #log of the width
        }

        self.maddm_run = run_interface
        self.line_experiment = "line_Fermi-LAT_2015_"

        self.param_blocks, _ = self.maddm_run.param_card.analyze_param_card()
        #self.parameter_vars = [] #names of parameters to scan over
        self.output_observables = []
        self.likelihood_parts = []

        self.counter = 0
        #self.Fermi   = Fermi_bounds()
        
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
                self.output_observables.append('Omegah^2')
                self.output_observables.append('sigmav(xf)')
                self.output_observables.append('x_f')
                self.output_observables.append('xsi')
                
            if self.maddm_run.mode['direct']:
                self.output_observables.append('sigmaN_SI_n')
                self.output_observables.append('sigmaN_SI_p')
                self.output_observables.append('sigmaN_SD_n')
                self.output_observables.append('sigmaN_SD_p')
                self.output_observables.append('Nevents')

            if self.maddm_run.mode['indirect']:
                if 'inclusive' in self.maddm_run.maddm_card['sigmav_method']:
                    self.output_observables.append('tot_SM_xsec')
                else:
                    self.output_observables.append('taacsID')

                #self.output_observables.append('Fermi_sigmav')
                #self.output_observables.append('like_nonth')
                
                #detailled_keys = [k for k in self.maddm_run.last_results.keys() if k.startswith('taacsID#')]
                #for key in detailled_keys:
                #    self.output_observables.append(key)

            if self.maddm_run.mode['spectral']:
                xsecs = [k for k in self.maddm_run.last_results.keys() if k.startswith("taacsID#") and self.maddm_run.is_spectral_finalstate(self.maddm_run.str_processes[k.replace("taacsID#","")])]
                for proc in xsecs:
                    self.output_observables.append(proc)

            if self.maddm_run.mode['capture']:
                detailled_keys = [k for k in self.maddm_run.last_results.keys() if k.startswith('ccap_')]
                for key in detailled_keys:
                    self.output_observables.append(key)

            #FIX THIS! HERE ADD FLUXES        
            #if self.maddm_run.mode['CR_flux']:
            #    detailled_keys = [k for k in self.maddm_run.last_results.keys() if k.startswith('??????????')]
            #    for key in detailled_keys:
            #        self.output_observables.append(key)
        
        # add likelihoods parts for each category
        self.likelihood_parts = [mode for mode in ['relic', 'direct', 'indirect', 'spectral'] if self.maddm_run.mode[mode]]

        # number of parameters to output
        # this includes the parameters which are scanned over (they go in first)
        # and the output parameters like relic density, dd cross section etc ...
        n_parameters = len(parameters)
        n_dimensions = n_parameters
        #mnest.parameter_vars = [parameters[i][0] for i in range(n_parameters)]

        n_parameters=n_parameters+len(self.output_observables)+len(self.likelihood_parts)


        logger.info("Multinest will run with the following parameters: " )
        for key, item in six.iteritems(self.options):
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
        for i, obs in enumerate(self.likelihood_parts):
            logger.info('column[%i] : -2*log(L_%s)' % (i+2+len(self.options['parameters'])+len(self.output_observables), obs))


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
                for i, obs in enumerate(self.likelihood_parts):
                    f.write('column[%i] : -2*log(L_%s)\n' % (i+2+len(self.options['parameters'])+len(self.output_observables), obs))



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
                        width_val = secure_float_f77(spline[3])
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
                        self.options['parameters'].append([var_name, secure_float_f77(min), secure_float_f77(max)])
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
                self.maddm_run.param_card[block].get(lhaid).value = val
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
            logger.error('Not a valid prior choice! Only allowed priors are %s '% str(list(Priors.priors.keys())))


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
            omegah2 = results['Omegah^2']
        if self.maddm_run.mode['direct']:
            spinSI = 0.5*(results['sigmaN_SI_p'] + results['sigmaN_SI_n'])# * GeV2pb * pb2cm2
            spinSDp = results['sigmaN_SD_p']#  * GeV2pb * pb2cm2
            spinSDn = results['sigmaN_SD_n']# * GeV2pb * pb2cm2
        #<=========== for ID we will need each channel separately.
        
        #if self.maddm_run.mode['indirect']:   
            # sigmavID = {'tot': results['taacsID']}
            # id_vel = self.maddm_run.maddm_card['vave_indirect']
            #detailled_keys = [k.split("#")[1] for k in results.keys() if k.startswith('taacsID#')]
            #for key in detailled_keys:
            #    logger.debug('taacsID key: %s', key)
            #    sigmavID[key] = results['taacsID#%s' %(key)]
#            logger.debug('sigmavID : %s' % sigmavID)

        mdm = self.maddm_run.param_card.get_value('mass', self.maddm_run.proc_characteristics['dm_candidate'][0])
        #logger.debug('MDM: %.3e', mdm)
        logger.debug(self.maddm_run.mode['relic'])
        logger.debug(self.maddm_run.mode['direct'])

        ndim_obs = ndim + len(self.output_observables)

        for obs, likelihood in six.iteritems(self.options.get('loglikelihood')):

            #relic density
            if obs == 'relic' and self.maddm_run.mode['relic']:

                #if relic density evaluates to -1 (too early freezeout for example)
                #then discard this point by assigning some high log likelihood.
                if omegah2 < 0:
                    relic_contrib= __mnestlog0__
                else:
                    if likelihood == 'gaussian':
                        relic_contrib = -0.5*pow(omegah2 - self.maddm_run.limits._oh2_planck,2)/pow(self.maddm_run.limits._oh2_planck_width,2)
                    elif likelihood == 'half_gauss':
                        if omegah2 > 0:
                            #relic_contrib = np.log(0.5*(np.tanh((self.maddm_run.limits._oh2_planck - omegah2)\
                            #                           /self.maddm_run.limits._oh2_planck)+1.000001))
                            if omegah2 > self.maddm_run.limits._oh2_planck:
                                relic_contrib+= -0.5*pow(np.log10(self.maddm_run.limits._oh2_planck/omegah2),2)\
                                      /pow(np.log10(self.maddm_run.limits._oh2_planck_width),2)
                    elif likelihood =='user':
                        relic_contrib=0
                        #
                        # HERE ADD YOUR OWN LOG LIKELIHOOD FOR RELIC DENSITY
                        #
                    elif likelihood == '':
                        relic_contrib=0
                    else:
                        logger.error('You are not using a valid likelihood function for relic density. Omitting the contribution!')
                        relic_contrib += 0
                # add to cube
                cube[ndim_obs+self.likelihood_parts.index('relic')] = relic_contrib
                chi += relic_contrib

            #direct detection (SI)
            if self.maddm_run.mode['direct']:
                cube[ndim_obs+self.likelihood_parts.index('direct')] = 0
            if obs == 'directSI' and self.maddm_run.mode['direct']:
                if likelihood=='half_gauss':

                    if spinSI > self.maddm_run.limits.SI_max(mdm):
                        direct_contrib_1= -0.5*pow(np.log10(self.maddm_run.limits.SI_max(mdm)/spinSI),2)\
                                    /pow(self.options['half_gauss_width']['spinSI'],2)

                elif likelihood == 'gaussian':
                    if self.maddm_run.limits._sigma_SI > 0:
                        direct_contrib_1=  -0.5*pow(spinSI - self.maddm_run.limits._sigma_SI,2)/pow(self.maddm_run.limits._sigma_SI_width,2)
                    else:
                        logger.error('You have to set up the sigma_SI(_width) to a positive value to use gaussian likelihood!')
                elif likelihood == 'user':
                    #
                    # HERE ADD YOUR OWN LOG LIKELIHOOD FOR SI
                    #
                    direct_contrib_1=0

                elif likelihood == '':
                    direct_contrib_1=0
                else:
                    logger.error('You are not using a valid likelihood function for SI direct detection. Omitting the contribution!')
                    direct_contrib_1 = 0

                # add to cube
                cube[ndim_obs+self.likelihood_parts.index('direct')] += direct_contrib_1
                chi += direct_contrib_1

            #direct detection (SD) proton and neutron
            if obs.startswith('directSD') and self.maddm_run.mode['direct']:
                nucleon = obs.split('_')[1]
                if likelihood=='half_gauss':
                    if nucleon =='p':
                        if spinSDp > self.maddm_run.limits.SD_max(mdm, 'p'):
                            direct_contrib_2= -0.5*pow(np.log10(self.maddm_run.limits.SD_max(mdm, 'p')/spinSDp),2)\
                                  /pow(self.options['half_gauss_width']['spinSDp'],2)
                    elif nucleon == 'n':
                        if spinSDp > self.maddm_run.limits.SD_max(mdm,'n'):
                            direct_contrib_2= -0.5*pow(np.log10(self.maddm_run.limits.SD_max(mdm, 'n')/spinSDn),2)\
                                  /pow(self.options['half_gauss_width']['spinSDn'],2)
                elif likelihood == 'gaussian':
                    if nucleon == 'p' and self.maddm_run.limits._sigma_SDp > 0:
                        direct_contrib_2=  -0.5*pow(spinSDp - self.maddm_run.limits._sigma_SDp,2)/pow(self.maddm_run.limits._sigma_SDp_width,2)
                    else:
                        logger.error('You have to set up the sigma_SDp(_width) to a positive value to use gaussian likelihood!')
                    if nucleon == 'n' and self.maddm_run.limits.sigma_SDn > 0:
                        direct_contrib_2=  -0.5*pow(spinSDn - self.maddm_run.limits.sigma_SDn,2)/pow(self.maddm_run.limits._sigma_SDn_width,2)
                    else:
                        logger.error('You have to set up the sigma_SDp(_width) to a positive value to use gaussian likelihood!')
                elif likelihood == 'user':
                    #
                    # HERE ADD YOUR OWN LOG LIKELIHOOD FOR SD
                    #
                    direct_contrib_2=0
                elif likelihood == '':
                    direct_contrib_2=0
                else:
                    logger.error('You are not using a valid likelihood function for SD direct detection. Omitting the contribution!')
                    direct_contrib_2=0
                # add to cube
                cube[ndim_obs+self.likelihood_parts.index('direct')] += direct_contrib_2
                chi += direct_contrib_2

            # indirect detection:
            # using Fermi Likelihood for combined gamma spectrum (only SM channels and relative SM xsection for PPPC sepctra)
            if obs == 'indirect' and self.maddm_run.mode['indirect']:
                if likelihood != '':
                    logger.warning("Likelihood for indirect detection (continuum spectrum) can't be specified, because it is taken from the Fermi-LAT searches in dSph.")
                indirect_contrib = results['like_nonth']
                # add to cube
                cube[ndim_obs+self.likelihood_parts.index('indirect')] = indirect_contrib
                chi += indirect_contrib

                
                '''
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
                  '''

            if obs == 'spectral' and self.maddm_run.mode['spectral']:
                if likelihood != '':
                    logger.warning("Likelihood for indirect detection (line spectrum) can't be specified, because it is taken from the Fermi-LAT 2015 searches in the Galactic Centre.")
                # find the most constraining peak: -2log(L) higher, but exclude -1 (not computed or out of range peaks)
                like_peaks = [v for k, v in six.iteritems(results) if k.startswith(self.line_experiment + 'like') and v != -1.]
                spectral_contrib = max(like_peaks) if len(like_peaks) != 0 else 0.
                # add to cube
                cube[ndim_obs+self.likelihood_parts.index('spectral')] = spectral_contrib
                chi += spectral_contrib

                
        #Print the counter on the screen

        if self.counter % 100 ==0:
            toprint = 'Scanned over %i points.' % self.counter
            logger.info(toprint)
        self.counter += 1

        return chi



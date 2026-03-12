from unittest import case
import sys
sys.path.append('/Users/davidimig/projects/squirel/plotting')
from utils import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams


#!/usr/bin/env python
# coding: utf-8

# The goal of this script is to plot the total CMB temperature spectrum C_l^TT,
# as well as its decomposition in different contributions:
# - T + SW: intrinsic temperature plus Sachs-Wolfe correction
# - early-ISW: early integrated Sachs-Wolfe
# - late-ISW: late integrated Sachs-Wolfe
# - Doppler contribution
# - total unlensed spectrum
# - total lensed spectrum


from math import pi


####################################################
#
# Cosmological parameters and other CLASS parameters
#
####################################################
theta_s100 = 1.051783             # Angular size of the sound horizon, exactly 100(ds_dec/da_dec)
omega_m_LCDM = 0.1431354439            # Reduced total matter density (Omega*h^2) (Exactly, omega_m = omega_b + omega_cdm + omega_mnu with Mnu=0.06 eV)
omega_cdm_LCDM = 0.1201075             # Reduced cold dark matter density in absence of LiMRs (Omega*h^2)
omega_b = omega_m_LCDM - omega_cdm_LCDM
N_mnu = 1
M_mnu = 0.06
omega_gamma = 2.47298e-05
def_a_eq = (1+7/8*pow(4/11,4/3)*3.044)*omega_gamma/(omega_cdm_LCDM+omega_b)

nu_ur_array = [
    3.044,
    2.0308,
    1.0176,
    0.00441
]
common = {
    'output':'tCl,pCl,lCl',
    'l_max_scalars':2500,
    'lensing': 'yes',
    'background_verbose': 3,
    #'thermodynamics_verbose': 1,
    'perturbations_verbose': 1,
    'input_verbose': 1,
    #'transfer_verbose': 1,
    #'primordial_verbose': 1,
    #'harmonic_verbose': 1,
    #'fourier_verbose': 1,
    #'lensing_verbose': 1,
    'output_verbose': 1,
    'N_ur':nu_ur_array[N_mnu],
    'omega_b':omega_b,
    'A_s':2.100549e-09,
    'n_s':0.9660499,
    'tau_reio':0.05430842,
    '100*theta_s':theta_s100,
}

###############
#
# Initiate a CLASS instance
#
M = Class()
#
###############
#
# call CLASS for each contribution
#
###############
#
# 
contribs = ['tsw','eisw','lisw','dop']
fig = plt.figure(figsize=(16, 16))
cls = {}
cases = ['LCDM', 'FD', 'BE'] # not acutally going to use BE distribution, just for looping over znrs and coloring 
z_NRs = {
    'FD': 1e3,
    'BE': 1e4,
}
Delta_Neff = 0.25

for contrib in contribs:
    params = {}
    ax = plt.subplot(2,2,contribs.index(contrib)+1)
    ax.set_title(contrib)
    for case in cases:
        if N_mnu>0:
            m_ncdm_str = ','.join([str(M_mnu/N_mnu) for x in range(N_mnu)]) # Baseline for LCDM and DR cosmologies, to be appended for LiMR cosmologies
            T_ncdm_str = ','.join(['0.71611' for x in range(N_mnu)])  # Standard neutrino temperature today in units of T_CMB0
            params['N_ncdm'] = N_mnu
            deg_ncdm_str = ','.join(['1' for x in range(N_mnu)]) # g = 2 for massive neutrinos
        
        if case == 'LCDM':
            params['omega_cdm'] = omega_cdm_LCDM
        else:
            sigma = 0.5 #doesn't effect anything
            z_NR = z_NRs[case]
            head_dir='../data/cosmos/'
            qc_dict, m_dict = LiMR_parameters(Delta_Neff,z_NR,[sigma],head_dir)
            omega_chi0 = omega_chi(Delta_Neff,z_NR)
            A_func, B_func, C_func = chi_scaling(Delta_Neff,z_NR,'FD',sigma)
            omega_c = def_a_eq**(-1)*(1+7/8*pow(4/11,4/3)*3.044)*omega_gamma + (B_func(def_a_eq)-C_func(def_a_eq))*omega_chi0-omega_b
            params['omega_cdm'] = omega_c
            
            if N_mnu>0:
                m_ncdm_str = str(m_dict['FD'])+','+m_ncdm_str # Append LiMR mass to baseline neutrino masses if any
                T_ncdm_str = str(qc_dict['FD'])+','+T_ncdm_str # Append LiMR characteristic momentum to baseline neutrino temperatures if any	
                deg_ncdm_str = str(float(g_dict['FD']/2))+','+deg_ncdm_str # Potentially non-standard deg_ncdm for LiMR, default for massive nus
            else:
                m_ncdm_str = str(m_dict['FD'])
                T_ncdm_str = str(qc_dict['FD'])
                deg_ncdm_str = str(float(g_dict['FD']/2))

            params['ncdm_psd_parameters'] = ncdm_flags('FD',sigma)
            params['ncdm_quadrature_strategy'] = ncdm_strategy('FD',sigma,N_mnu)
            params['ncdm_N_momentum_bins'] = ncdm_bins('FD',N_mnu)
            params['ncdm_a'] = ncdm_scaler('FD',N_mnu)
            params['ncdm_maximum_q'] = ncdm_qmax('FD',N_mnu)
            params['N_ncdm'] = N_mnu+1

        params['m_ncdm'] = m_ncdm_str
        params['T_ncdm'] = T_ncdm_str
        params['deg_ncdm'] = deg_ncdm_str

        M.empty() # reset input
        M.set(common) # new input
        M.set(params) # new input

        M.set({'temperature contributions':contrib})

        cls[contrib] = M.raw_cl(3000)

        ell = cls[contrib]['ell']
        factor = 1.e10*ell*(ell+1.)/2./pi
        plt.plot(ell,factor*cls[contrib]['tt'],color=cosmo_color(case),label=z_NRs[case] if case != 'LCDM' else 'LCDM')
    ax.set_xlim(2., 3000)
    ax.set_ylim(0,5)
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel(r"$\ell (\ell+1) C_l^{TT} / 2 \pi \,\,\, [\times 10^{10}]$")
    ax.grid()

#
plt.legend(loc='lower left',bbox_to_anchor=(1.4, 0.5))
figname = 'cltt_terms.pdf'
plot_dir = 'extra_plots/'
log(f'Saving: {plot_dir + figname}')
plt.savefig(plot_dir+figname,bbox_inches='tight')

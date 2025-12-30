from utils import *
import numpy as np
import matplotlib.pyplot as plt
from classy import Class

# this scripts illustrates how to adjust omega_c to fix a_eq in LiMR cosmologies
# this builds on the varying_neff.py CLASS script and various studies which fix a_eq with extra ur species

Delta_Neffs = [0.1,0.05,0.02]
Delta_Neffs = [0.1,0.1,0.1]
z_NRs =[1e3,1e4,1e5]
# Delta_Neffs = [.3,.6,.9]
# z_NRs = np.ones(3)*10
omega_chis = []

omega_b = 2.255065e-02
omega_cdm_standard = 1.193524e-01
omega_gamma = 2.47298e-05
omega_m_standard = omega_b + omega_cdm_standard

common_settings = {# fixed LambdaCDM parameters (Planck 18 + lensing + BAO bestfit)
                   'omega_b':omega_b,
                   'A_s':2.100549e-09,
                   'n_s':0.9660499,
                   'tau_reio':0.05430842,
                   # output and precision parameters
                   'output':'tCl,pCl,lCl,mPk',
                   'P_k_max_1/Mpc':10,
                   'l_max_scalars':2500,
                   'lensing': 'yes',
                   'N_ur':3.044,
                   '100*theta_s':1.041783,
}

M = {}
M['LCDM'] = Class()
M['LCDM'].set(common_settings)
M['LCDM'].compute()

clTT = {}
clTT['LCDM'] = M['LCDM'].lensed_cl()['tt'][2:]
ell = M['LCDM'].lensed_cl()['ell'][2:]
#
def_a_eq = (1+7/8*pow(4/11,4/3)*3.044)*omega_gamma/(omega_cdm_standard+omega_b)
log(f' Default a_eq = {def_a_eq}, corresponding to z_eq = {1/def_a_eq-1}')
for i, Delta_Neff in enumerate(Delta_Neffs):
    omega_chis.append(omega_chi(Delta_Neff,z_NRs[i]))
    a_NR = 1/(z_NRs[i]+1)
    print()
    log(f'For Delta_Neff = {Delta_Neff}, z_NR = {z_NRs[i]}, omega_chi = {omega_chis[i]}')
    # test_CLASS(case='FD', Delta_Neff=Delta_Neff, z_NR=z_NRs[i])

    # this is the fundamentally new component
    A_func, B_func, C_func = chi_scaling(Delta_Neff,z_NRs[i],case='FD',sigma=None)

    omega_c = def_a_eq**(-1)*(1+7/8*pow(4/11,4/3)*3.044)*omega_gamma + (B_func(def_a_eq)-C_func(def_a_eq))*omega_chis[i]-omega_b

    log(f'For Delta_Neff = {Delta_Neff}, z_NR = {z_NRs[i]}, omega_c = {omega_c}')
    
    T0_dict, m_dict = LiMR_parameters(Delta_Neff,z_NRs[i])
    T0, m = T0_dict['FD'], m_dict['FD']

    M[i] = Class()
    # # set input parameters
    M[i].set(common_settings)
    M[i].set({
        'N_ncdm':1,
        'm_ncdm':m,
        'T_ncdm':T0,
        'omega_cdm':omega_c,
        # 'omega_m':omega_m_standard
    })

    M[i].compute()
    
    a_eq = M[i].get_current_derived_parameters(['a_eq'])['a_eq']
    log(f'Output a_eq = {a_eq}, corresponding to z_eq = {1/a_eq-1}')
    log(f'Percent different from def a_eq = {100*(a_eq-def_a_eq)/def_a_eq}')

    clTT[i] = M[i].lensed_cl()['tt'][2:]

    plt.plot(ell,(clTT[i]-clTT['LCDM'])/clTT['LCDM'],label=f'Delta N_eff = {Delta_Neff}, z_NR = {z_NRs[i]}')
plt.legend()
plt.show()
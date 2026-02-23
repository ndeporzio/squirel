from utils import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

h = 0.67810                       # Dimensionless reduced Hubble parameter (H_0 / (100km/s/Mpc))
theta_s100 = 1.051783             # Angular size of the sound horizon, exactly 100(ds_dec/da_dec)
omega_m_LCDM = 0.1431354439
omega_cdm_LCDM = 0.1201075             # Reduced cold dark matter density in absence of LiMRs (Omega*h^2)
omega_b = omega_m_LCDM - omega_cdm_LCDM
m_nu_on_omega_nu = 93.13858

def plot_m_and_f(fixed1='h', fixed2='omega_m', head_dir='../data/cosmos/'):
    m_ncdms = [1, 0.25, 0.30, 0.33, 0.35, 0.37]
    f_ncdms = [0.00001, 0.0023, 0.005, 0.007, 0.009, 0.011]
    # m_ncdms = 30*np.ones(6)
    # f_ncdms = [0.00001, 0.01, 0.02, 0.03, 0.04, 0.05]
    # f_ncdms = [0.00001, 0.1, 0.3, 0.5, 0.8, 0.999]
    ells = []
    lensed_clTTs = []
    for i in range(len(m_ncdms)):
        m_ncdm = m_ncdms[i]
        f_ncdm = f_ncdms[i]
        omega_ncdm = f_ncdm*omega_cdm_LCDM
        T_ncdm = pow(omega_ncdm * m_nu_on_omega_nu / m_ncdm, 1/3) * nu_factor
        cosmo = Class()
        cosmo.set({
			'output':'tCl,pCl,lCl,mPk',
			'P_k_max_1/Mpc':10,
			'l_max_scalars':2500,
			'lensing': 'yes',
			'omega_m':omega_m_LCDM,
			'omega_b':omega_b,
			'N_ncdm':1,
			'm_ncdm':m_ncdm,
			'T_ncdm':T_ncdm,
			# 'omega_ncdm':omega_ncdm,
			'output_verbose': 1,
			'background_verbose': 3,
            '100*theta_s':theta_s100,
        })
        cosmo.compute()
        f_ncdm = cosmo.get_current_derived_parameters(['f_ncdm1'])['f_ncdm1']
        print(f_ncdm)
        ells.append(cosmo.lensed_cl()['ell'][2:])
        lensed_clTTs.append(cosmo.lensed_cl()['tt'][2:])
        plt.plot(ells[i], (lensed_clTTs[i]-lensed_clTTs[0])/lensed_clTTs[0], label=f'f_ncdm={f_ncdm:.1f}')
        
    errbars_f = "../data/Planck_errbars.txt"
    data = np.loadtxt(
        errbars_f,
        delimiter=",",
        comments="#",
        skiprows=3,
        usecols=(2, 4, 5, 6)
    )

    planck_ell     = data[:, 0]
    planck_dl_tt   = data[:, 1]/muK2 
    sigma_minus    = data[:, 2]/muK2
    sigma_plus     = data[:, 3]/muK2
    plt.fill_between(planck_ell, -sigma_minus/planck_dl_tt, sigma_plus/planck_dl_tt, color='gray', alpha=0.3, zorder=-1)
    # plt.errorbar(planck_ell, np.zeros(len(planck_ell)), yerr=(sigma_minus/planck_dl_tt, sigma_plus/planck_dl_tt), color='gray', alpha=0.3, zorder=-1)

    set_xy_lims(xmin=2., xmax=2500, ymin=-0.05, ymax=0.05)
    plt.legend(loc='lower center')
		
if __name__ == '__main__':
	fig = plt.figure(figsize=(8, 8))
	ax = plt.subplot(1, 1, 1)
	plot_m_and_f()
	figname = 'm_and_f.pdf'
	plot_dir = 'extra_plots/'
	save_fig(figname, plot_dir)

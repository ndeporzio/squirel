import numpy as np
import matplotlib.pyplot as plt
from scipy.special import zeta
import matplotlib as mpl
from scipy.integrate import quad
from scipy import constants
from classy import Class
# modules for data I/O
import pickle
import os 
T0CMB = 2.7255

#
# Cosmology functions
#
	
# defining distributions
def f_FD(xi):
	return 1./(np.exp(xi)+1)
def f_BE(xi):
	return 1./(np.exp(xi)-1)
def f_RD(xi):
	return 2.19*pow(xi,-5./2.)*np.exp(-0.74*pow(xi,2.))
def f_LN(xi,sigma):
	return 1./(xi*sigma*np.sqrt(2*np.pi))*np.exp(-pow(np.log(xi),2)/(2*pow(sigma,2)))

# n-th moment of f(xi)
def Qn(f,n,sigma=1):
	if f in [f_FD,f_BE,f_RD]:
		return quad(lambda xi: xi**(2+n)*f(xi), 0, np.inf)[0]
	elif f == f_LN:
		return quad(lambda xi: xi**(2+n)*f(xi,sigma), 0, np.inf)[0]
	else: 
		print('[utils.py] (ERROR) Distrib not recognised, assuming FD.')
		return quad(lambda xi: xi**(2+n)*f_FD(xi), 0, np.inf)[0]

# characteristic momentum from 1st moment and Delta Neff (in units of the photon temp)
def T0chi(Q1,Delta_Neff,g):
	T0nu = 1.95
	return pow(Delta_Neff/(g/2)/(Q1/(7*np.pi**4/120)),1/4)*T0nu/T0CMB

# mass in eV
def m_chi(T0chi,z_NR,Q0,Q1):
	K_to_eV = constants.physical_constants['Boltzmann constant in eV/K'][0]
	return T0chi*T0CMB*K_to_eV*(z_NR+1)*Q1/Q0
	
def LiMR_parameters(Delta_Neff,z_NR):
	Q0_dict = {
		'FD':Qn(f_FD,0),
		'BE':Qn(f_BE,0),
		'RD':Qn(f_RD,0),
		'LN':Qn(f_LN,0,0.5),
	}

	Q1_dict = {
		'FD':Qn(f_FD,1),
		'BE':Qn(f_BE,1),
		'RD':Qn(f_RD,1),
		'LN':Qn(f_LN,1,0.5),
	}

	g_dict = {
		'g':1,
		'g':2,
		'g':2,
		'g':2,
	}

	T0_dict = {
		'FD': T0chi(Q1_dict['FD'],Delta_Neff,2),
		'BE': T0chi(Q1_dict['BE'],Delta_Neff,1),
		'RD': T0chi(Q1_dict['RD'],Delta_Neff,2),
		'LN': T0chi(Q1_dict['LN'],Delta_Neff,2),
	}

	m_dict = {
		'FD': m_chi(T0_dict['FD'],z_NR,Q0_dict['FD'],Q1_dict['FD']),
		'BE': m_chi(T0_dict['BE'],z_NR,Q0_dict['BE'],Q1_dict['BE']),
		'RD': m_chi(T0_dict['RD'],z_NR,Q0_dict['RD'],Q1_dict['RD']),
		'LN': m_chi(T0_dict['LN'],z_NR,Q0_dict['LN'],Q1_dict['LN']),
	}

	return T0_dict, m_dict
	
def fill_LiMR_parameters(Delta_Neff=0.3, z_NR=1000, output_dir='../data/distribution_data/'):
	filename = output_dir + f'LiMR_parameters_DNeff={Delta_Neff}_zNR={z_NR}.pkl'
	if os.path.exists(filename):
		with open(filename, 'rb') as f:
			params = pickle.load(f)
		print('[utils.py] Loaded LiMR parameters from:', filename)
	else:
		params = LiMR_parameters(Delta_Neff, z_NR)
		with open(filename, 'wb') as f:
			pickle.dump(params, f)
		print('[utils.py] Saved LiMR parameters to:', filename)
	return params

def_T0_dict, def_m_dict = fill_LiMR_parameters()
	
def run_CLASS_and_save(case, Delta_Neff=0.3, z_NR=1e3, T0_dict=def_T0_dict, m_dict=def_m_dict, fixed='h', output_dir='../data/distribution_data/'):
	h = 0.67810                       # Dimensionless reduced Hubble parameter (H_0 / (100km/s/Mpc))
	theta_s100 = 1.041783             # Angular size of the sound horizon, exactly 100(ds_dec/da_dec)
	omega_m = 0.1431354439            # Reduced total matter density (Omega*h^2) (Exactly, omega_m = omega_b + omega_cdm + omega_mnu with Mnu=0.06 eV)
	omega_cdm = 0.1201075             # Reduced cold dark matter density in absence of LiMRs (Omega*h^2)
	
    # Initialize CLASS
	cosmo = Class()

    # Use default parameters
	cosmo.set({'output':'tCl,pCl,lCl,mPk',
			   'P_k_max_1/Mpc':10,
			   'l_max_scalars':2500,
			   'lensing': 'yes',
			   'write_background':'yes',            # Write background parameter table
			   'background_verbose': 3,
			   'thermodynamics_verbose': 1,
			   'perturbations_verbose': 1,
			   'input_verbose': 1,
			   'transfer_verbose': 1,
			   'primordial_verbose': 1,
			   'harmonic_verbose': 1,
			   'fourier_verbose': 1,
			   'lensing_verbose': 1,
			   'output_verbose': 1,
			})
	
	print('[utils.py] Computing CLASS with', case, 'parameters.')
    # Modify parameters according to case
	if case == 'LCDM':
		cosmo.set({'N_ncdm':1,
                   'm_ncdm':0.06,
                   'T_ncdm':0.71611,
                   'N_ur':2.0308,
                   'omega_m':omega_m,
                   })
		root = output_dir+case+'_fixed='+fixed
	elif case == 'DR':
		cosmo.set({'N_ncdm':1,
                   'm_ncdm':0.06,
                   'T_ncdm':0.71611,
                   'N_ur':2.0308+Delta_Neff,
                   'omega_m':omega_m,
                   })
		root = output_dir+case+'_DNeff='+Delta_Neff+'_fixed='+fixed
	elif case == 'FD' or case == 'BE' or case == 'RD' or case == 'LN':
		cosmo.set({'N_ncdm':2,
                   'm_ncdm':[0.06, m_dict[case]],
                   'T_ncdm':[0.71611, T0_dict[case]],
                   'omega_m':omega_m,
                   'N_ur':2.0308,
                   })
		root = output_dir+case+'_DNeff='+Delta_Neff+'_zNR='+z_NR+'_fixed='+fixed
	else:
		print('[utils.py] (ERROR) Case not recognised, using LCDM parameters.')
		cosmo.set({'N_ncdm':1,
                   'm_ncdm':0.06,
                   'T_ncdm':0.71611,
                   'N_ur':2.0308,
                   'omega_m':omega_m,
                   })
		root = output_dir+'LCDM_fixed='+fixed
	if fixed == 'theta_s100':
		print('[utils.py] Fixed 100*theta_s requested, H0 will instead vary.')
		cosmo.set({'100*theta_s':theta_s100})
	
	# pip installed class allows writing parameters but not output files oddly, those have to be saved manually
	cosmo.set({
		'write parameters': 'yes',
		'root':root+'_',
	})

	cosmo.compute()

	# save output to file
	filename = root + '_output.pkl'
	bg = cosmo.get_background()
	unlensed_cls = cosmo.raw_cl()
	lensed_cls = cosmo.lensed_cl()
	kvec = np.logspace(-4,1,2500)
	pk = [cosmo.pk(kk, 0.0) for kk in kvec]
	derived = cosmo.get_current_derived_parameters(['rs_rec', 'H0', 'z_reio', 'a_eq', '100*theta_s', 'z_rec'])
	output_data = {
		'background': bg,
		'unlensed_cls': unlensed_cls,
		'lensed_cls': lensed_cls,
		'kvec': kvec,
		'pk': pk,
		'derived': derived,
	}
	with open(filename, 'wb') as f:
		pickle.dump(output_data, f)
	print('[utils.py] CLASS output saved to:', filename)

	cosmo.struct_cleanup()
	cosmo.empty()

	return output_data
	
def fill_cosmos(Delta_Neff=0.3, z_NR=1e3, fixed='h', output_dir='../data/distribution_data/'):
	for case in ['LCDM', 'DR', 'FD', 'BE', 'RD', 'LN']:
		filename = ''
		if case == 'LCDM':
			filename = output_dir+case+'_fixed='+fixed+'_output.pkl'
		elif case == 'DR':
			filename = output_dir+case+'_DNeff='+str(Delta_Neff)+'_fixed='+fixed+'_output.pkl'
		else:
			filename = output_dir+case+'_DNeff='+str(Delta_Neff)+'_zNR='+str(z_NR)+'_fixed='+fixed+'_output.pkl'
		if os.path.exists(filename):
			with open(filename, 'rb') as f:
				output_data = pickle.load(f)
			print('[utils.py] Loaded CLASS output for', case, 'from:', filename)
		else:
			T0_dict, m_dict = fill_LiMR_parameters(Delta_Neff, z_NR, output_dir)
			output_data = run_CLASS_and_save(case, Delta_Neff, z_NR, T0_dict, m_dict, fixed, output_dir)
		globals()[f'cosmo_{case}'] = output_data

#
# Plotting functions
#
def plot_distributions(Delta_Neff=0.3,z_NR=1e3):
	fill_cosmos(Delta_Neff=0.3,z_NR=1e3)
	qarr = np.geomspace(1e-2, 100., 20000)

	Tnu0 = 1.95
	K_TO_CM = 4.366

def plot_distributions():
	temp_Gauss_func = lambda Amp,ystar,sigma,y : Amp * np.exp(-(y-ystar)**2/(2*sigma**2))
	temp_FD_func    = lambda y : 1./(np.exp(y)+1)
	qarr = np.geomspace(1e-2, 50., 10000)

	Tnu0 = 1.95
	K_TO_CM = 4.366

	plt.plot(qarr, 0.04 * qarr * Tnu0**3 * K_TO_CM**3 * temp_FD_func(qarr) * 6 * qarr**2 / (2 * np.pi**2),
			c=cosmo_color('LCDM'), lw=2.2)
	plt.plot(qarr, 0.04 * qarr * Tnu0**3 * K_TO_CM**3 * temp_Gauss_func(33.8595, 0.1, 0.294218, qarr) * 6 * qarr**2 / (2 * np.pi**2), 
		c=cosmo_color('LEDR'), lw=2.2)
	plt.plot(qarr, 0.40 * qarr * Tnu0**3 * K_TO_CM**3 * temp_Gauss_func(0.0000161608, 30, 4.82113, qarr) * 6 * qarr**2 / (2 * np.pi**2), 
		c=cosmo_color('HE'),lw=2.2)
	plt.plot(qarr, 0.40 * qarr * Tnu0**3 * K_TO_CM**3 * temp_Gauss_func(0.000125406, 3.0, 8.82654, qarr) * 6 * qarr**2 / (2 * np.pi**2), 
		c=cosmo_color('HEDR'), lw=2.2)
	plt.plot(qarr, 0.04 * qarr * Tnu0**3 * K_TO_CM**3 * (temp_Gauss_func(0.0743352, 3.5, 0.508274, qarr) + temp_FD_func(qarr / 0.7003)) * 6 * qarr**2 / (2 * np.pi**2),
		c=cosmo_color('LTM'), lw=2.2)


def add_cosmo_cases():
	xoff, yoff, offset = 0.15, -0.0, 0.03
	plt.text(0.2 - xoff, 0.90 - yoff,r"\boldmath{$\Lambda$}\textbf{CDM}",transform=plt.gca().transAxes,
		fontsize=14,
		color=cosmo_color('LCDM'),
		rotation=0)
	plt.text(0.35, 0.09, r"\boldmath{$\Lambda$}\textbf{CDM}",
		transform=plt.gca().transAxes,
		fontsize=10,
		color=cosmo_color('LCDM'),
		rotation=42)
	plt.text(0.38 - xoff, 0.89 - yoff,r"$\sum m_\nu = 0.12\,\mathrm{eV}$, $N_{\rm eff}^\nu = 3.044$",
		transform=plt.gca().transAxes,
		fontsize=14,
		color='k',
		rotation=0)

	plt.text(0.2 - xoff, 0.82 + offset - yoff,r"\textbf{L}\boldmath{$\nu$}\textbf{-DR}",
		transform=plt.gca().transAxes,
		fontsize=14,
		color=cosmo_color('LEDR'),
		rotation=0)
	plt.text(0.11, 0.09,r"\textbf{L}\boldmath{$\nu$}\textbf{-DR}",
			transform=plt.gca().transAxes,
		fontsize=10,
		color=cosmo_color('LEDR'),
		rotation=62)
	plt.text(0.38 - xoff, 0.81 + offset - yoff,r"$\sum m_\nu = 0.12\,\mathrm{eV}$, $N_{\rm eff}^\nu = 0.5$, $N_{\rm eff}^{\rm DR} = 2.544$",
		transform=plt.gca().transAxes,
		fontsize=14,
		color='k',
		rotation=0)

	plt.text(0.2 - xoff, 0.74 + 2 * offset - yoff, r"\textbf{H}\boldmath{$\nu$}",
		transform=plt.gca().transAxes,
		fontsize=14,
		color=cosmo_color('HE'),
		rotation=0)
	plt.text(0.77, 0.32, r"\textbf{H}\boldmath{$\nu$}",
		transform=plt.gca().transAxes,
		fontsize=12,
		color=cosmo_color('HE'),
		rotation=80)
	plt.text(0.38 - xoff, 0.73 + 2 * offset - yoff,r"$\sum m_\nu = 1.20\,\mathrm{eV}$, $N_{\rm eff}^\nu = 3.044$",
		transform=plt.gca().transAxes,
		fontsize=14,
		color='k',
		rotation=0)

	plt.text(0.2 - xoff, 0.66 + 3 * offset - yoff,r"\textbf{H}\boldmath{$\nu$}\textbf{-DR}",
		transform=plt.gca().transAxes,
		fontsize=14,
		color=cosmo_color('HEDR'),
		rotation=0)
	plt.text(0.61, 0.09,r"\textbf{H}\boldmath{$\nu$}\textbf{-DR}",
		transform=plt.gca().transAxes,
		fontsize=10,
		color=cosmo_color('HEDR'),
		rotation=63)
	plt.text(0.38 - xoff, 0.65 + 3 * offset - yoff,r"$\sum m_\nu = 1.20\,\mathrm{eV}$, $N_{\rm eff}^\nu = 1.5$, $N_{\rm eff}^{\rm DR} = 1.544$",
		transform=plt.gca().transAxes,
		fontsize=14,
		color='k',
		rotation=0)

	plt.text(0.2 - xoff, 0.58 + 4 * offset - yoff,r"\textbf{LT+Mid}",transform=plt.gca().transAxes,
		fontsize=14,
		color=cosmo_color('LTM'),
		rotation=0)
	plt.text(0.455,0.24,r"\textbf{LT+Mid}",transform=plt.gca().transAxes,
		fontsize=10,
		color=cosmo_color('LTM'),
		rotation=83)
	plt.text(0.38 - xoff, 0.57 + 4 * offset - yoff,r"$\sum  m_\nu = 0.12\,\mathrm{eV}$, $N_{\rm eff}^{\nu} = 3.044$",
		transform=plt.gca().transAxes,
		fontsize=14,
		color='k',
		rotation=0)

	plt.text(0.2 - xoff, 0.55 + 2 * offset - yoff,r"$\Omega_{\nu,0} / \Omega_{\mathrm{m},0} = 0.009$, $N_{\rm eff} = 3.044$",
		transform=plt.gca().transAxes,
		fontsize=16,
		color='k',
		rotation=0)
	
def set_xy_lims(xmin=10, xmax=1000, ymin=1, ymax=1000):
	plt.xlim(xmin, xmax)
	plt.ylim(ymin, ymax)

def add_xy_labels(xlabel=r'$m_\mathrm{lightest}\,\mathrm{[meV]}$', ylabel=r'$n_\nu^\mathrm{loc.}\,\mathrm{[cm}^{-3}\mathrm{]}$', fontsize=20):
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)

def set_xy_scales(xscale='linear', yscale='log'):
	plt.xscale(xscale)
	plt.yscale(yscale)

HIGHP_LABEL = r'$\mathrm{High-}p_\nu$'
LOWT_LABEL = r'$\mathrm{Low-}T_\nu\mathrm{+DR}$'
LCDM_LABEL = r'$\Lambda\mathrm{CDM}$'

def cosmo_color(case='LCDM'):
    # colors_dict = {
    #     'FD': '#dc267f',
    #     'BE': '#785ef0',
    #     'RD': '#648fff',
    #     'LN': '#ffb000',
    #     # '??': '#fe6100',
    # }
	colors_dict = {
	 'LCDM': '#003f5c',
	 'LEDR': '#58508d',
	 'HE': '#bc5090',
	 'HEDR': '#ff6361',
	 'LTM': '#ffa600'
	 }

	try:
		return colors_dict[case]
	except:
		# print('[utils.py] (ERROR) Case not recognised, available cosmo scenarios are:')
		# print('[utils.py] \t(1) Key: FD - FD distribution')
		# print('[utils.py] \t(2) Key: BE - BE distribution')
		# print('[utils.py] \t(3) Key: RD - Out-of-equilibrium relativistic decay product distribution')
		# print('[utils.py] \t(4) Key: LN - Log-normal proxy distribution')
		print('[utils.py] (ERROR) Case not recognised, available cosmo scenarios are:')
		print('[utils.py] \t(1) Key: LCDM - LCDM with FD distribution')
		print('[utils.py] \t(2) Key: LEDR - Low energy neutrinos with additional DR')
		print('[utils.py] \t(3) Key: HE - High Energy neutrinos')
		print('[utils.py] \t(4) Key: HEDR - High Energy neutrinos with additional DR')
		print('[utils.py] \t(5) Key: LTM - Low temperature, but with additional Gaussian component')
		return 'k'
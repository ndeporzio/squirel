from utils import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({
	"text.usetex": True,
	"font.family": "serif",
	"xtick.direction": "in",
	"ytick.direction": "in",
})

def cosmo_densities(plot_dir='plots/', figname='cosmo_densities.pdf'):

	z_NRs = np.logspace(3,5,3)
	Delta_Neffs = [0.05, 0.1, 0.15][::-1]

	fig = plt.figure(figsize=(8, 5))
	ax = plt.subplot(1, 1, 1)

	# Delta_Neffs = [0.1,0.1,0.1]
	# # Delta_Neffs = [0.1,0.05,0.02]
	# z_NRs =[1e3,1e4,1e5]
	# # z_NRs = np.ones(3)*1e4
	# Delta_Neffs = np.linspace(0.1,0.02,5)
	# z_NRs = np.logspace(3,5,5)
	# Delta_Neffs = [0.1,0.05,0.02,0.06]
	# z_NRs =[1e3,1e4,1e5,38602.36]

	plot_cosmo_densities(Delta_Neffs,z_NRs,N_mnu=1,M_mnu=0.06,head_dir='../data/cosmos/')

	plt.xlabel(r'$\mathrm{Scale\,\,Factor}\,\,a/a_0 = (1 + z)^{-1}$', fontsize=14)
	plt.ylabel(r'$\mathrm{Density\,\,Parameter}\,\,\Omega_{i}(a)$', fontsize=14)
	plt.xlim(1e-6, 1e0)
	plt.ylim(1e-4, 2e0)
	fig.tight_layout()

	save_fig(figname, plot_dir)

def Hubble_changes(plot_dir='extra_plots/', figname='Hubble_changes.pdf'):
	fig = plt.figure(figsize=(8, 5))
	ax = plt.subplot(1, 1, 1)

	Delta_Neffs = [0.1,0.1,0.1]
	# Delta_Neffs = [0.1,0.05,0.02]
	z_NRs =[1e3,1e4,1e5]
	# z_NRs = np.ones(3)*1e4
	fixed2 = 'omega_m'
	# fixed2 = 'a_eq'
	plot_Hubbles(Delta_Neffs,z_NRs,N_mnu=1,M_mnu=0.06,fixed2=fixed2,head_dir='../data/cosmos/')

	plt.xlabel(r'$\mathrm{Scale\,\,Factor}\,\,a/a_0 = (1 + z)^{-1}$', fontsize=14)
	# # plt.ylabel(r'$\mathrm{Density\,\,Parameter}\,\,\Omega_{i}(a)$', fontsize=16)
	plt.ylabel(r'$\Delta H(a)/H(a)|_{\Lambda \rm{CDM}}$', fontsize=14)
	plt.xlim(1e-6, 1e0)
	plt.ylim(0, 0.02)

	figname = name_fig(figname, plot_dir)
	plt.savefig(plot_dir+figname)
	log(f'Saving: {plot_dir + figname}')

def vs_LCDM_Cl(plot_dir='plots/', figname='cosmos_vs_LCDM.pdf'):
	fig = plt.figure(figsize=(8, 4))
	ax = plt.subplot(1, 1, 1)

	# Delta_Neffs = [0.1,0.1,0.1]
	# Delta_Neffs = [0.1,0.05,0.02,0.06]
	# z_NRs = np.ones(3)*1e4
	# Delta_Neffs = np.linspace(0.1,0.02,5)
	# Delta_Neffs = [0.06, 0.14, 0.28, 0.52, 0.87]
	# z_NRs = [38602.36, 31059.37, 25982.36, 22331.94, 19580.91]


	lensed = True
	spectras = ['tt','ee']
	# fixed1 = 'h'
	fixed2 = 'omega_m'
	fixed1 = 'theta_s100'
	# fixed2 = 'a_eq'
	for spectra in spectras:
		plot_vs_LCDM_Cl(Delta_Neffs,z_NRs,N_mnu=1,M_mnu=0.06,lensed=lensed,spectra=spectra, fixed1=fixed1,fixed2=fixed2,head_dir='../data/cosmos/')

		set_xy_scales(xscale='linear', yscale='linear')

		ax.set_xlabel(r"$\ell$", fontsize=14)
		if spectra == 'tt':
			ax.set_ylabel(r"$(C^{\rm TT}_\ell-C^{\rm TT}_\ell|_{\Lambda {\rm CDM}})/C^{\rm TT}_\ell|_{\Lambda {\rm CDM}}$", fontsize=14)
			set_xy_lims(xmin=2., xmax=2500, ymin=-0.05, ymax=0.05)
		elif spectra == 'ee':
			ax.set_ylabel(r"$(C^{\rm EE}_\ell-C^{\rm EE}_\ell|_{\Lambda {\rm CDM}})/C^{\rm EE}_\ell|_{\Lambda {\rm CDM}}$", fontsize=14)
			set_xy_lims(xmin=2., xmax=2500, ymin=-0.05, ymax=0.05)
			plot_dir = 'extra_plots/'

		text = False
		if text:
			cs = IBM_cscheme()
			for i, Delta_Neff in enumerate(Delta_Neffs):
				ax.text(.5,8*Delta_Neff,r'$\Delta N_{\mathrm{eff}}='+str(Delta_Neff)+', z_{\mathrm{NR}}=10^{'+str(np.log10(z_NRs[i]))+'}$',fontsize=13,transform=ax.transAxes,color=cs[2*i%5])

		spectra_name = spectra + '_'
		fixed_str = 'fixed='+fixed1+','+fixed2+'_'
		figname2 = fixed_str + spectra_name +figname
		save_fig(figname2, plot_dir)

		plt.cla() 

def subplots_vs_LCDM_Cl(plot_dir='plots/', figname='cosmos_vs_LCDM.pdf'):
	Delta_Neffs = np.linspace(0.05,0.25,5)
	# Delta_Neffs = np.linspace(0.04,0.20,5)
	z_NRs = np.logspace(3,5,3)[0:2]
	# z_NRs = np.logspace(3,5,5)[0:3]
	# z_NRs = np.logspace(3,6,7)[0:6]
	# Delta_Neffs = np.linspace(0.01,0.05,5)

	fig = plt.figure(figsize=(10, 4))
	# fig = plt.figure(figsize=(12, 4))

	lensed = True
	spectras = ['tt']
	# fixed1 = 'h'
	fixed2 = 'omega_m'
	fixed1 = 'theta_s100'
	fixed2 = 'a_eq'
	fixed_str = 'fixed='+fixed1+','+fixed2+'_'
	for spectra in spectras:
		for i, z_NR in enumerate(z_NRs):
			ax = plt.subplot(1, len(z_NRs), i+1, sharex = ax if i > 0 else None, sharey = ax if i > 0 else None)
			viable_idx = len(Delta_Neffs)
			if z_NR == 1e5 and Delta_Neffs[-1] > 0.21:
				viable_idx -= 1
			fixed_zNRs = np.ones(len(Delta_Neffs))*z_NR

			f_LiMRs = plot_vs_LCDM_Cl(Delta_Neffs[:viable_idx],fixed_zNRs[:viable_idx],N_mnu=1,M_mnu=0.06,lensed=lensed,spectra=spectra, fixed1=fixed1,fixed2=fixed2,head_dir='../data/cosmos/')
			set_xy_scales(xscale='linear', yscale='linear')

			ax.set_xlabel(r"$\ell$", fontsize=14)
			if i == 0:
				if spectra == 'tt':
					ax.set_ylabel(r"$(C^{\rm TT}_\ell-C^{\rm TT}_\ell|_{\Lambda {\rm CDM}})/C^{\rm TT}_\ell|_{\Lambda {\rm CDM}}$", fontsize=14)
					set_xy_lims(xmin=2., xmax=2500, ymin=-0.05, ymax=0.05)
				elif spectra == 'ee':
					ax.set_ylabel(r"$(C^{\rm EE}_\ell-C^{\rm EE}_\ell|_{\Lambda {\rm CDM}})/C^{\rm EE}_\ell|_{\Lambda {\rm CDM}}$", fontsize=14)
					set_xy_lims(xmin=2., xmax=2500, ymin=-0.05, ymax=0.05)
					plot_dir = 'extra_plots/'
			else:
				ax.tick_params(labelleft=False)
			
			ax.text(0.5,0.9,r'$z_{\mathrm{NR}}=10^{'+str(round(np.log10(z_NRs[i]),2))+'}$',horizontalalignment='center',fontsize=14,transform=ax.transAxes,color='k')
			
			for i in range(len(f_LiMRs)):
				if len(f_LiMRs) == 4:
					ylocs = [0.15,0.15,0.05,0.05]
					xlocs = [0.325,0.675,0.325,0.675]
				elif len(f_LiMRs) == 5:
					ylocs = [0.25,0.15,0.15,0.05,0.05]
					xlocs = [0.325,0.325,0.675,0.325,0.675]
				ax.text(xlocs[i],ylocs[i],r'$f_\chi='+f'{f_LiMRs[i]:.1}'+r'$',horizontalalignment='center',fontsize=14,transform=ax.transAxes,color=IBM_cscheme()[i])

		fig.tight_layout(rect=[0, 0, 0.92, 1])
		fig.subplots_adjust(wspace=0.0)
		
		cbar_ax = fig.add_axes([plt.gca().get_position().x1 + 0.02, plt.gca().get_position().y0, 0.02, plt.gca().get_position().height])
		cmap = mpl.colors.ListedColormap(IBM_cscheme()[0:len(Delta_Neffs)])
		spacing = (Delta_Neffs[-1]-Delta_Neffs[0])/(len(Delta_Neffs))
		bounds = [Delta_Neffs[0]-spacing/2, Delta_Neffs[-1]+spacing/2]
		norm = mpl.colors.Normalize(vmin=bounds[0], vmax=bounds[-1])
		cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')
		cb1.set_label(r'$\Delta N_{\mathrm{eff}}$', fontsize=14)
		cb1.set_ticks(Delta_Neffs)

		# fig.tight_layout()

		spectra_name = spectra + '_'
		figname2 = 'subplots_' + fixed_str + spectra_name +figname
		plt.savefig(plot_dir+figname2)
		log(f'Saving: {plot_dir + figname2}')

		plt.cla() 

def distribs_comparison(Delta_Neff=0.3,z_NR=1e3,sigma_array=[0.04,1.5],ins=True, plot_dir='plots/', figname='distribs_comparison.pdf'):

	fig = plt.figure(figsize=(8, 7))
	
	ax = plt.subplot(1, 1, 1)

	ymax = plot_distributions(ax,Delta_Neff,z_NR,sigma_array=sigma_array)
	# Plot settings

	set_xy_lims(xmin=1e-4, xmax=1e4, ymin=1e-1, ymax=1e2)
	set_xy_scales(xscale='log', yscale='log')

	# set_xy_lims(xmin=1e-4, xmax=1e4, ymin=0, ymax=15)
	# set_xy_scales(xscale='log', yscale='linear')

	#ax.set_xlabel(r"$q_\nu \equiv p_\nu/T_{\nu, 0}$")
	ax.set_xlabel(r"$\xi_\chi \equiv q_\chi/q^c_\chi$",fontsize=16)
	ax.set_ylabel(r"$\frac{\mathrm{d}\rho^\mathrm{NR}_\chi}{\mathrm{d} \log \xi_\chi} \, \mathrm{[eV\ cm}^{-3}\mathrm{]}$",fontsize=16)	#ax.set_ylabel(r"$\frac{\mathrm{d}\rho^\mathrm{NR}_\nu}{\mathrm{d} \log q_\nu} \, \mathrm{[eV\ cm}^{-3}\mathrm{]}$")
	# add text saying what DNeff and z_NR are
	ax.text(0.7,0.75,r'$\Delta N_{\mathrm{eff}}='+str(Delta_Neff)+'$',fontsize=16,transform=ax.transAxes,color='k')
	ax.text(0.7,0.7,r'$z_{\mathrm{NR}}=10^{'+str(np.log10(z_NR))+'}$',fontsize=16,transform=ax.transAxes,color='k')

	
	ax.legend(fontsize=14)
	# add_cosmo_cases()

	# and add an inset box zoomed in on the y < 20 region
	if ins:
		# axins = fig.add_axes([0.55, 0.25, 0.32, 0.32])
		axins = fig.add_axes([0.17, 0.5, 0.3, 0.35])
		plot_distributions(axins,Delta_Neff,z_NR,sigma_array=sigma_array,ins=True)
		set_xy_lims(xmin=1e-4, xmax=1e4, ymin=0, ymax=5)
		set_xy_scales(xscale='log', yscale='linear')
	
	save_fig(figname, plot_dir)

def evols_comparison(Delta_Neff=0.3,z_NR=1e3,sigma_array=[0.04,1.5],plot_dir='plots/', figname='evols_comparison.pdf'):
	fig = plt.figure(figsize=(8, 7))

	ax = plt.subplot(1, 1, 1)

	plot_evolution(Delta_Neff,z_NR,sigma_array=sigma_array)

	set_xy_lims(xmin=0.1*z_NR, xmax=10*z_NR, ymin=1, ymax=1.7)
	set_xy_scales(xscale='log', yscale='linear')

	ax.set_xlabel(r"$z$",fontsize=16)
	ax.set_ylabel(r"$\rho_\chi/\rho_\mathrm{asympt.}$",fontsize=16)
	plt.text(0.85*z_NR, 1.05, r'$z_{\mathrm{NR}}$', rotation=90, fontsize=16, color='gray')

	save_fig(figname, plot_dir)

def distrib_cls(Delta_Neff=0.2,z_NR=1e3,sigma_array=[0.04,1.5],plot_dir='plots/', figname='distrib_cls.pdf'):
	fig = plt.figure(figsize=(8, 7))

	ax = plt.subplot(1, 1, 1)

	denom = 'FD'
	fixed1 = 'h'
	fixed2 = 'omega_m'
	fixed1 = 'theta_s100'
	fixed2 = 'a_eq'
	plot_distrib_cls(Delta_Neff,z_NR,sigma_array=sigma_array,denom=denom,lensed=True,fixed1=fixed1,fixed2=fixed2)

	# set_xy_lims(xmin=2., xmax=2500, ymin=-0.015, ymax=0.01)
	set_xy_lims(xmin=2., xmax=2500, ymin=-0.002, ymax=0.002)
	set_xy_scales(xscale='linear', yscale='linear')

	ax.set_xlabel(r"$\ell$")
	ax.set_ylabel(r"$(C^{\rm TT}_\ell-C^{\rm TT}_\ell|_{{\rm FD}})/C^{\rm TT}_\ell|_{ {\rm FD}}$")
	fixed_str = 'fixed='+fixed1+','+fixed2+'_'
	figname = fixed_str + figname
	save_fig(figname, plot_dir)

def chi_scaling(plot_dir='extra_plots/', figname='chi_scaling.pdf'):
	fig = plt.figure(figsize=(8, 5))
	ax = plt.subplot(1, 1, 1)

	Delta_Neff = 0.1
	z_NR = 1e3
	sigma = None

	plot_scaling(Delta_Neff,z_NR,case='FD',sigma=sigma)
	
	ax.set_xlabel(r'$a$')
	ax.set_ylabel(r'$R(a)$')

	save_fig(figname, plot_dir)

def vs_LCDM_contours(plot_dir='plots/', figname='LCDM_vs_DR_vs_LiMR_triangle.pdf', chain_dir='/Users/davidimig/projects/squirel/data/chains'):
	
	log('Generating comparison contours for LCDM vs DR vs FD LiMRs.')

	params = [
    'H0',
    'raw_omega_b',
	# 'omega_b',
    'omega_cdm',
    'ln10^{10}A_s',
    'tau_reio',
	# 'z_reio',
    'n_s',

    'delta_Neff',
    'log10z_tr',

	'100theta_s',
	'rs_rec',
	'ra_rec',

    # 'Omega_Lambda',
    'omega_m',
	'a_eq',
    # 'log10omega0_ncdm1',
    # 'log10m_ncdm_in_eV',
    # 'T_ncdm1',
    'log10f_ncdm1',
    ]

	probes = [
		'PR3',
	]
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
	]
	cases = [
		'LCDM', 
		'DR', 
		'FD',
		# 'zup_FD',
		]
	
	# cases = ['FD', 'logzNR4_FD', 'fixed_zNR_FD'][::-1]
	# cases = ['logzNR4_FD', 'LCDM', 'DR', 'FD', 'fixed_zNR_FD']
	# figname = 'fixed_zNR_FD_triangle.pdf'
	# mcmc_plot_lims['omega_cdm'] = [0.05, None]

	roots = [
		case + '_' + probes[0] + '_' + nu_hierachs[0] + '_' for case in cases
	][::-1] # reverse to plot FD first since it has largest contours

	g = plots.get_subplot_plotter(chain_dir=chain_dir)
	samples = get_samples(g, roots=roots)
	roots = [root.replace('fixed_zNR_FD','logzNR5_FD') for root in roots]
	roots = [root.replace('zup_FD','FD') for root in roots]

	g.triangle_plot(
		samples, 
		params, 
		filled=True,
		param_limits=mcmc_plot_lims,
		contour_args = [
			{
				# 'alpha':1,
				# 'alpha':np.linspace(1,0.1,10)[i],
				'color': cosmo_color(root.split('_')[0]),
				'lw': 2,
			} for root in roots
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cosmo_color(root.split('_')[0]),
				'lw': 1,
			} for root in roots
		],
	)

	# figname = name_fig(figname, plot_dir)
	g.export(plot_dir+figname)
	log(f'Saving: {plot_dir + figname}')
	
def LiMR_contours(plot_dir='plots/', figname='PR3_LiMR_triangle.pdf', chain_dir='/Users/davidimig/projects/squirel/data/chains'):
	
	log('Generating PR3 LiMR distribution comparison contours.')

	params = [
    # 'H0',
    # 'raw_omega_b',
    # 'omega_cdm',
    # 'ln10^{10}A_s',
    # 'tau_reio',
    # 'n_s',

    'delta_Neff',
    'log10z_tr',
	'sigma',
	# 'log10omega0_ncdm1',
    ]
	
	probe = 'PR3'
	nu_hierach = 'N_mnu=1_Mmnu=0.06'
	cases = [
		'FD',
		'BE',
		'LN',
	]
	roots = [case + '_' + probe + '_' + nu_hierach + '_' for case in cases]

	g = plots.get_subplot_plotter(chain_dir=chain_dir)
	samples = get_samples(g, roots=roots)

	g.triangle_plot(
		samples, 
		params, 
		filled=True,
		param_limits=mcmc_plot_lims,
		contour_args = [
			{
				'alpha':1,
				# 'alpha':np.linspace(1,0.1,10)[i],
				'color': cosmo_color(case),
				'lw': 2,
			} for i, case in enumerate(cases)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cosmo_color(case),
				'lw': 1,
			} for case in cases
		],
	)

	# figname = name_fig(figname, plot_dir)
	g.export(plot_dir+figname)
	log(f'Saving: {plot_dir + figname}')

def S4_contours(plot_dir='plots/', figname='S4_triangle.pdf', chain_dir='../data/chains'):
	
	log('Generating S4 LiMR distribution comparison contours.')

	params = [
		'delta_Neff',
		'log10z_tr',
		'sigma',
	]

	probe = 'S4'
	nu_hierach = 'N_mnu=1_Mmnu=0.06'
	cases = [
		# 'FD',
		'LN',
		'LN_FD_fid',
	]
	roots = [case + '_' + probe + '_' + nu_hierach + '_' for case in cases]

	g = plots.get_subplot_plotter(chain_dir=chain_dir)
	samples = get_samples(g, roots=roots)

	cs = IBM_cscheme()
	g.triangle_plot(
		samples, 
		params, 
		filled=True,
		param_limits=mcmc_plot_lims,
		contour_args = [
			{
				# 'alpha':1,
				'alpha':np.linspace(1,0.1,10)[i],
				'color': cosmo_color(case) if case != 'LN_FD_fid' else cs[2],
				'lw': 2,
			} for i, case in enumerate(cases)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cosmo_color(case) if case != 'LN_FD_fid' else cs[2],
				'lw': 1,
			} for case in cases
		],
	)

	# figname = name_fig(figname, plot_dir)
	g.export(plot_dir+figname)
	log(f'Saving: {plot_dir + figname}')

def dataset_contours(plot_dir='plots/', figname='FD_datasets_triangle.pdf', chain_dir='/Users/davidimig/projects/squirel/data/chains'):
	
	log('Generating comparison contours for FD LiMRs in different datasets.')

	params = [
    'H0',
    'raw_omega_b',
	# 'omega_b',
    'omega_cdm',
    'ln10^{10}A_s',
    'tau_reio',
	# 'z_reio',
    'n_s',

    'delta_Neff',
    'log10z_tr',

	# '100theta_s',
	# 'rs_rec',
	# 'ra_rec',

    # 'Omega_Lambda',
    # 'omega_m',
	# 'a_eq',
    # 'log10omega0_ncdm1',
    # 'log10m_ncdm_in_eV',
    # 'T_ncdm1',
    # 'log10f_ncdm1',
    ]

	probes = [
		'PR3',
		'PR3+lens',
		'PR3+eBOSS',
		'PR4',
		# 'priors',
	]
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
		# 'N_mnu=2_Mmnu=0.11',
	]
	case = 'FD'
	roots = [
		case + '_' + probe + '_' + nu_hierachs[0] + '_' for probe in probes
	]
	# roots = [roots[0]]

	g = plots.get_subplot_plotter(chain_dir=chain_dir)
	samples = get_samples(g, roots=roots)

	cs = IBM_cscheme()
	g.triangle_plot(
		samples, 
		params, 
		filled=True,
		param_limits=mcmc_plot_lims,
		contour_args = [
			{
				# 'alpha':1,
				# 'alpha':np.linspace(1,0.1,10)[i],
				'color': cs[i],
				'lw': 2,
			} for i in range(len(roots))
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cs[i],
				'lw': 1,
			} for i in range(len(roots))
		],
	)

	# figname = name_fig(figname, plot_dir)
	g.export(plot_dir+figname)
	log(f'Saving: {plot_dir + figname}')

def nus_contours(plot_dir='plots/', figname='nu_prescriptions_triangle.pdf', chain_dir='../data/chains'):
	
	log('Generating comparison contours for FD LiMRs indifferent neutrino prescriptions.')

	params = [
    'H0',
    'raw_omega_b',
    'omega_cdm',
    'ln10^{10}A_s',
    'tau_reio',
    'n_s',

    'delta_Neff',
    'log10z_tr',
	# 'sigma',
    ]

	probe = 'PR3'
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
		# 'N_mnu=3_Mmnu=0.06',
		'N_mnu=2_Mmnu=0.11',
	]
	cases = [
		'FD',
		# 'LN',
	]
	roots = [cases[0] + '_' + probe + '_' + hierarch + '_' for hierarch in nu_hierachs]
	# roots.append(cases[1] + '_' + probe + '_' + nu_hierachs[1])

	g = plots.get_subplot_plotter(chain_dir=chain_dir)
	samples = get_samples(g, roots=roots)

	cs = IBM_cscheme()
	g.triangle_plot(
		samples, 
		params, 
		filled=True,
		param_limits=mcmc_plot_lims,
		contour_args = [
			{
				# 'alpha':1,
				'alpha':np.linspace(1,0.1,10)[i],
				'color': cosmo_color(root[0:2]) if i == 0 else cs[2*i],
				'lw': 2,
			} for i, root in enumerate(roots)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cosmo_color(root[0:2]) if i == 0 else cs[2*i],
				'lw': 1,
			} for i, root in enumerate(roots)
		],
	)

	# figname = name_fig(figname, plot_dir)
	g.export(plot_dir+figname)
	log(f'Saving: {plot_dir + figname}')

def triangle_3d(plot_dir='plots/', figname='3D_LiMR_triangle.pdf', chain_dir='/Users/davidimig/projects/squirel/data/chains'):
	params = [
    'delta_Neff',
    'log10z_tr',
	'log10f_ncdm1',
    ]
	
	probe = 'PR3'
	nu_hierach = 'N_mnu=1_Mmnu=0.06'
	case = 'FD'
	root = case + '_' + probe + '_' + nu_hierach + '_'

	g = plots.get_subplot_plotter(chain_dir=chain_dir)
	samples = get_samples(g, roots=[root])

	g.plot_3d(
		samples,
		params,
		color_bar=True,
		alpha_samples=True,
	)

	# g.triangle_plot(
	# 	samples, 
	# 	params, 
	# 	filled=True,
	# 	param_limits=mcmc_plot_lims,
	# 	contour_args = [
	# 		{
	# 			'alpha':1,
	# 			# 'alpha':np.linspace(1,0.1,10)[i],
	# 			'color': cosmo_color(case),
	# 			'lw': 2,
	# 		} for i, case in enumerate(cases)
	# 	],
	# 	line_args = [
	# 		{
	# 			'alpha': 1,
	# 			'zorder': -1,
	# 			'color': cosmo_color(case),
	# 			'lw': 1,
	# 		} for case in cases
	# 	],
	# )

	# figname = name_fig(figname, plot_dir)
	g.export(plot_dir+figname)
	log(f'Saving: {plot_dir + figname}')

def priors_contours(plot_dir='extra_plots/', figname='FD_triangle.pdf', chain_dir='/Users/davidimig/projects/squirel/data/chains'):
	
	# log('Generating prior only (no data) comparison contours for FD LiMRs of different parametrizations.')
	log('Generating FD comparison contours for different parametrizations of relic abundance.')

	params = [
    'H0',
    'raw_omega_b',
	# 'omega_b',
    'omega_cdm',
    'ln10^{10}A_s',
    # 'tau_reio',
	'z_reio',
    'n_s',

    'delta_Neff',
    'log10z_tr',
	'log10omega_LiMR',

	'100theta_s',
	'rs_rec',
	'ra_rec',

    'Omega_Lambda',
    'omega_m',
	'a_eq',
    'log10omega0_ncdm1',
    'log10m_ncdm_in_eV',
    'T_ncdm1',
    'log10f_ncdm1',
	'f_ncdm1'
    ]

	probes = [
		'PR3',
		# 'PR3+lens',
		# 'PR3+eBOSS',
		# 'PR4',
		'priors',
	]
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
		# 'N_mnu=2_Mmnu=0.11',
	]
	cases = ['FD', 'om_L_FD']
	roots = [
		case + '_' + probes[0] + '_' + nu_hierachs[0] + '_' for case in cases
	]
	# roots.append( cases[0] + '_' + probes[0] + '_' + nu_hierachs[0] + '_' )
	# roots = [roots[0]]

	g = plots.get_subplot_plotter(chain_dir=chain_dir)
	samples = get_samples(g, roots=roots)

	cs = IBM_cscheme()
	g.triangle_plot(
		samples, 
		params, 
		filled=True,
		param_limits=mcmc_plot_lims,
		contour_args = [
			{
				# 'alpha':1,
				# 'alpha':np.linspace(1,0.1,10)[i],
				'color': cs[i],
				'lw': 2,
			} for i in range(len(roots))
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cs[i],
				'lw': 1,
			} for i in range(len(roots))
		],
	)

	# figname = name_fig(figname, plot_dir)
	g.export(plot_dir+figname)
	log(f'Saving: {plot_dir + figname}')

def loglike_contours(plot_dir='extra_plots/', figname='loglikes_triangle.pdf', chain_dir='/Users/davidimig/projects/squirel/data/chains', Dneffbound=True):
	
	log('Generating FD LiMR contours color coded by log-likelihood.')

	params = [
    'H0',
    'raw_omega_b',
    'omega_cdm',
    'ln10^{10}A_s',
    'tau_reio',
    'n_s',

    'delta_Neff',
    'log10z_tr',
	# 'sigma',
	# 'log10omega0_ncdm1',
	# 'omega_m',
    ]
	
	probe = 'PR3'
	nu_hierach = 'N_mnu=1_Mmnu=0.06'
	cases = [
		'FD',
		# 'fixed_zNR_FD',
		# 'logzNR4_FD',
		# 'Tup_FD',
	]
	# figname = 'Tup_FD_triangle.pdf'

	roots = [case + '_' + probe + '_' + nu_hierach + '_' for case in cases]

	g = plots.get_subplot_plotter(chain_dir=chain_dir)
	samples = get_samples(g, roots=roots)
	samples[0].addDerived(
		samples[0].loglikes,
		name="minus_log_P",
		label=r"$-\log(P)$",
	)
	samples[0].updateSettings({'max_scatter_points': 3000})
	samples.append(samples[0]) # necessary to add contour lines to scatter plot
	

	g.triangle_plot(
		samples, 
		params, 
		plot_3d_with_param="minus_log_P",
		filled=False,
		param_limits=mcmc_plot_lims,
		contour_args = [
			{
				# 'alpha':1,
				# 'alpha':np.linspace(1,0.1,10)[i],
				'color': cosmo_color('FD'),
				'lw': 1,
			} for i in range(len(samples))
		],
		line_args = [
			{
				# 'alpha':1,
				'zorder': -1,
				'color': cosmo_color('FD'),
				'lw': 1,
			} for i in range(len(samples))
		],
		# max_scatter_points=50000,
	)

	if Dneffbound:
		# Add vertical lines at each DNeff bound
		x_param, y_param = 'delta_Neff', 'log10z_tr'

		bounds = {
			'DR_bound': 2.9533e-01, # from my DR contours, consistent with 2207.13133
			'logzNR4_bound': 1.0345e-01,
			'logzNR5_bound': 7.9728e-02,
		}

		ax = g.get_axes_for_params(x_param, y_param)
		for bound in bounds:
			ax.axvline(bounds[bound], c=cosmo_color(bound.split('_')[0]), ls='--', lw=2)

	figname = name_fig(figname, plot_dir)
	g.export(plot_dir+figname)
	log(f'Saving: {plot_dir + figname}')

if __name__ == '__main__':
	log('Starting main.py')
	# loglike_contours()
	# test_CLASS(case='FD',fixed1='theta_s100',fixed2='a_eq',save=False)
	# scaling = input('[main.py] Scaling (y/n): ')
	density = input('[main.py] Densities (y/n): ')
	gen_cls = input('[main.py] LCDM vs DR vs LiMR Cls (y/n): ')
	dists = input('[main.py] Distributions (y/n): ')
	evols = input('[main.py] Evolutions (y/n): ')
	dist_cls = input('[main.py] distribution comparison Cls (y/n): ')
	LiMR_triangle = input('[main.py] LCDM vs DR vs LiMR triangle (y/n): ')
	distribs_triangle = input('[main.py] LiMR distribs triangle (y/n): ')
	S4distribs_triangle = input('[main.py] S4 triangle (y/n): ') 
	datasets_triangle = input('[main.py] FD datasets triangle (y/n): ')
	nus_triangle = input('[main.py] Neutrino prescriptions triangle (y/n): ')
	Hubbles = input('[main.py] Hubble changes (y/n): ')
	priors = input('[main.py] Parametrization contours (y/n): ')
	
	# 3d_triangle = input('[main.py] 3D LiMR triangle (y/n): ')


	# if scaling == 'y':
	# 	chi_scaling()

	if density == 'y':
		cosmo_densities()
	if gen_cls == 'y':
		# vs_LCDM_Cl()
		subplots_vs_LCDM_Cl()
	if Hubbles == 'y':
		Hubble_changes()
		
	Delta_Neff = 0.1
	z_NR = 1e3
	sigma_min = 0.04
	sigma_max = 1.5
	n_sigma = 10
	sigma_array = np.linspace(sigma_min,sigma_max,n_sigma)
	if dists == 'y':
		distribs_comparison(Delta_Neff=Delta_Neff,z_NR=z_NR,sigma_array=np.linspace(sigma_min,sigma_max,200))
	if evols == 'y':
		evols_comparison(Delta_Neff=Delta_Neff,z_NR=z_NR,sigma_array=sigma_array)
	if dist_cls == 'y':
		distrib_cls(Delta_Neff=Delta_Neff,z_NR=z_NR,sigma_array=sigma_array)

	if LiMR_triangle == 'y':
		vs_LCDM_contours()
	if datasets_triangle == 'y':
		dataset_contours()
	if distribs_triangle == 'y':
		LiMR_contours()
	if S4distribs_triangle == 'y':
		S4_contours()
	if nus_triangle == 'y':
		nus_contours()
	if priors == 'y':
		priors_contours()
	# if 3d_triangle == 'y':
	# 	triangle_3d()
from utils import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({
	"text.usetex": True,
	"font.family": "serif",
})

def cosmo_densities(plot_dir='plots/', figname='cosmo_densities.pdf'):
	fig = plt.figure(figsize=(8, 5))
	ax = plt.subplot(1, 1, 1)

	Delta_Neffs = [0.1,0.1,0.1]
	Delta_Neffs = [0.1,0.05,0.02]
	z_NRs =[1e3,1e4,1e5]

	plot_cosmo_densities(Delta_Neffs,z_NRs,N_mnu=1,M_mnu=0.06,head_dir='../data/')

	plt.xlabel(r'$\mathrm{Scale\,\,Factor}\,\,a/a_0 = (1 + z)^{-1}$', fontsize=16)
	plt.ylabel(r'$\mathrm{Energy\,\,Density}\,\,\Omega_{i}(a)$', fontsize=16)
	plt.xlim(1e-6, 1e0)
	plt.ylim(1e-4, 2e0)

	log(f'Saving: {plot_dir + figname}')
	plt.savefig(plot_dir + figname)

def vs_LCDM_Cl(plot_dir='plots/', figname='cosmos_vs_LCDM_.pdf'):
	fig = plt.figure(figsize=(6, 5))
	ax = plt.subplot(1, 1, 1)

	Delta_Neffs = [0.1,0.1,0.1]
	Delta_Neffs = [0.1,0.05,0.02]
	z_NRs =[1e3,1e4,1e5]
	# z_NRs = np.ones(5)*1e5

	lensed = True
	spectras = ['tt','ee']
	fixed1 = 'h'
	fixed2 = 'omega_m'
	fixed1 = 'theta_s100'
	fixed2 = 'a_eq'
	for spectra in spectras:
		plot_vs_LCDM_Cl(Delta_Neffs,z_NRs,N_mnu=1,M_mnu=0.06,lensed=lensed,spectra=spectra, fixed1=fixed1,fixed2=fixed2,head_dir='../data/')

		set_xy_scales(xscale='linear', yscale='linear')

		ax.set_xlabel(r"$\ell$")
		if spectra == 'tt':
			ax.set_ylabel(r"$(C^{\rm TT}_\ell-C^{\rm TT}_\ell|_{\Lambda {\rm CDM}})/C^{\rm TT}_\ell|_{\Lambda {\rm CDM}}$")
			set_xy_lims(xmin=2., xmax=2500, ymin=-0.025, ymax=0.02)
		elif spectra == 'ee':
			ax.set_ylabel(r"$(C^{\rm EE}_\ell-C^{\rm EE}_\ell|_{\Lambda {\rm CDM}})/C^{\rm EE}_\ell|_{\Lambda {\rm CDM}}$")
			set_xy_lims(xmin=2., xmax=2500, ymin=-0.05, ymax=0.05)

		text = False
		if text:
			cs = IBM_cscheme()
			for i, Delta_Neff in enumerate(Delta_Neffs):
				ax.text(.5,8*Delta_Neff,r'$\Delta N_{\mathrm{eff}}='+str(Delta_Neff)+', z_{\mathrm{NR}}=10^{'+str(np.log10(z_NRs[i]))+'}$',fontsize=13,transform=ax.transAxes,color=cs[2*i%5])

		spectra = spectra + '_'
		fixed_str = 'fixed='+fixed1+','+fixed2+'_'
		log(f'Saving: {plot_dir + fixed_str + spectra +figname}')
		plt.savefig(plot_dir + fixed_str +spectra + figname)
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
	ax.set_xlabel(r"$\xi_\chi \equiv q_\chi/T_\chi^0$",fontsize=16)
	ax.set_ylabel(r"$\frac{\mathrm{d}\rho^\mathrm{NR}_\chi}{\mathrm{d} \log \xi_\chi} \, \mathrm{[eV\ cm}^{-3}\mathrm{]}$",fontsize=16)	#ax.set_ylabel(r"$\frac{\mathrm{d}\rho^\mathrm{NR}_\nu}{\mathrm{d} \log q_\nu} \, \mathrm{[eV\ cm}^{-3}\mathrm{]}$")
	
	ax.legend(fontsize=14)
	# add_cosmo_cases()

	# and add an inset box zoomed in on the y < 20 region
	if ins:
		# axins = fig.add_axes([0.55, 0.25, 0.32, 0.32])
		axins = fig.add_axes([0.17, 0.5, 0.3, 0.35])
		plot_distributions(axins,Delta_Neff,z_NR,sigma_array=sigma_array,ins=True)
		set_xy_lims(xmin=1e-4, xmax=1e4, ymin=0, ymax=5)
		set_xy_scales(xscale='log', yscale='linear')
	
	log(f'Saving: {plot_dir + figname}')
	plt.savefig(plot_dir + figname)

def evols_comparison(Delta_Neff=0.3,z_NR=1e3,sigma_array=[0.04,1.5],plot_dir='plots/', figname='evols_comparison.pdf'):
	fig = plt.figure(figsize=(8, 7))

	ax = plt.subplot(1, 1, 1)

	plot_evolution(Delta_Neff,z_NR,sigma_array=sigma_array)

	set_xy_lims(xmin=0.1*z_NR, xmax=10*z_NR, ymin=1, ymax=1.7)
	set_xy_scales(xscale='log', yscale='linear')

	ax.set_xlabel(r"$z$",fontsize=16)
	ax.set_ylabel(r"$\rho_\chi/\rho_\mathrm{analytic}$",fontsize=16)
	plt.text(0.85*z_NR, 1.05, r'$z_{\mathrm{NR}}$', rotation=90, fontsize=16, color='gray')

	log(f'Saving: {plot_dir + figname}')
	plt.savefig(plot_dir + figname)

def distrib_cls(Delta_Neff=0.2,z_NR=1e3,sigma_array=[0.04,1.5],plot_dir='plots/', figname='distrib_cls.pdf'):
	fig = plt.figure(figsize=(8, 7))

	ax = plt.subplot(1, 1, 1)
	denom = 'FD'
	fixed1 = 'h'
	fixed2 = 'omega_m'
	fixed1 = 'theta_s100'
	# fixed2 = 'a_eq'
	plot_distrib_cls(Delta_Neff,z_NR,sigma_array=sigma_array,denom=denom,lensed=True,fixed1=fixed1,fixed2=fixed2)

	# set_xy_lims(xmin=2., xmax=2500, ymin=-0.015, ymax=0.01)
	set_xy_lims(xmin=2., xmax=2500, ymin=-0.002, ymax=0.002)
	set_xy_scales(xscale='linear', yscale='linear')

	ax.set_xlabel(r"$\ell$")
	ax.set_ylabel(r"$(C^{\rm TT}_\ell-C^{\rm TT}_\ell|_{{\rm FD}})/C^{\rm TT}_\ell|_{ {\rm FD}}$")
	fixed_str = 'fixed='+fixed1+','+fixed2+'_'
	log(f'Saving: {plot_dir + fixed_str + figname}')
	plt.savefig(plot_dir + fixed_str + figname)

def chi_scaling(plot_dir='extra_plots/', figname='chi_scaling.pdf'):
	fig = plt.figure(figsize=(8, 5))
	ax = plt.subplot(1, 1, 1)

	Delta_Neff = 0.1
	z_NR = 1e3
	sigma = None

	plot_scaling(Delta_Neff,z_NR,case='FD',sigma=sigma)
	
	ax.set_xlabel(r'$a$')
	ax.set_ylabel(r'$R(a)$')

	log(f'Saving: {plot_dir + figname}')
	plt.savefig(plot_dir + figname)
	

if __name__ == '__main__':
	log('Starting main.py')
	
	# test_CLASS(case='FD',fixed1='theta_s100',fixed2='a_eq',save=False)
	# scaling = input('[main.py] Scaling (y/n): ')
	density = input('[main.py] Densities (y/n): ')
	gen_cls = input('[main.py] LCDM vs DR vs LiMR Cls (y/n): ')
	dists = input('[main.py] Distributions (y/n): ')
	evols = input('[main.py] Evolutions (y/n): ')
	dist_cls = input('[main.py] distribution comparison Cls (y/n): ')

	# if scaling == 'y':
	# 	chi_scaling()
	if density == 'y':
		cosmo_densities()
	if gen_cls == 'y':
		vs_LCDM_Cl()
		
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

	
	# usedef = input('[main.py] Use default LiMR parameters? (y/n): ')
	# if usedef == 'n':
	# 	Delta_Neff = float(input('[main.py] Desired Delta_Neff: '))
	# 	z_NR = float(input('[main.py] Desired z_NR: '))
	# else:
	# 	Delta_Neff = 0.3
	# 	z_NR = 1000

	# print('[main.py] Using LiMR parameters: Delta_Neff =', Delta_Neff, ', z_NR =', z_NR)
	# distribs_comparison()


	# if dists == 'y':
	# 	distribs_comparison()
	# usedef = input('[main.py] Use default LiMR parameters? (y/n): ')
	# if usedef == 'n':
	# 	Delta_Neff = float(input('[main.py] Desired Delta_Neff: '))
	# 	z_NR = float(input('[main.py] Desired z_NR: '))
	# else:
	# 	Delta_Neff = 0.3
	# 	z_NR = 1000

	# print('[main.py] Using LiMR parameters: Delta_Neff =', Delta_Neff, ', z_NR =', z_NR)
	# distribs_comparison()
	# cosmo_distributions()

	# density = input('[main.py] Densities (y/n): ')
	# dists = input('[main.py] Distributions (y/n): ')
	# cmb = input('[main.py] CMB Sensitivity (y/n): ')
	# contour = input('[main.py] Countour Plot (y/n): ')
	# mpk = input('[main.py] Matter power spectrum (y/n): ')
	# appcontour = input('[main.py] Appendix contour (y/n): ')
	# nuonly = input('[main.py] Nu Appendix (y/n): ')
	# if density.lower()[0] == 'y':
	# 	cosmo_densities()
	# if dists.lower()[0] == 'y':
	# 	cosmo_distributions()
	# if cmb.lower()[0] == 'y':
	# 	cosmo_cmb_sensitivity()
	# if contour.lower()[0] == 'y':
	# 	cosmo_contour()
	# if mpk.lower()[0] == 'y':
	# 	cosmo_mpk()
	# if appcontour.lower()[0] == 'y':
	# 	cosmo_appendix_contour()
	# if nuonly.lower()[0] == 'y':
	# 	cosmo_appendix_nu()
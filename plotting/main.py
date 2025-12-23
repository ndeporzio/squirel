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

	Delta_Neffs = [0.2,0.2,0.2]
	z_NRs =[1e3,1e4,1e5]

	plot_cosmo_densities(Delta_Neffs,z_NRs,N_mnu=1,M_mnu=0.06,head_dir='../data/')

	plt.xlabel(r'$\mathrm{Scale\,\,Factor}\,\,a/a_0 = (1 + z)^{-1}$', fontsize=16)
	plt.ylabel(r'$\mathrm{Energy\,\,Density}\,\,\Omega_{i}(a)$', fontsize=16)
	plt.xlim(1e-6, 1e0)
	plt.ylim(1e-4, 2e0)

	log(f'Saving: {plot_dir + figname}')
	plt.savefig(plot_dir + figname)

def vs_LCDM_Cl(plot_dir='plots/', figname='cosmos_vs_LCDM.pdf'):
	fig = plt.figure(figsize=(6, 5))
	ax = plt.subplot(1, 1, 1)

	Delta_Neffs = [0.2,0.2,0.2]
	z_NRs =[1e3,1e4,1e5]
	lensed = True
	plot_vs_LCDM_Cl(Delta_Neffs,z_NRs,N_mnu=1,M_mnu=0.06,lensed=lensed,fixed='h',head_dir='../data/')

	set_xy_lims(xmin=2., xmax=2500, ymin=-0.08, ymax=0.06)
	set_xy_scales(xscale='linear', yscale='linear')

	ax.set_xlabel(r"$\ell$")
	ax.set_ylabel(r"$(C^{\rm TT}_\ell-C^{\rm TT}_\ell|_{\Lambda {\rm CDM}})/C^{\rm TT}_\ell|_{\Lambda {\rm CDM}}$")

	log(f'Saving: {plot_dir + figname}')
	plt.savefig(plot_dir + figname)

def distribs_comparison(Delta_Neff=0.3,z_NR=1e3,plot_dir='plots/', figname='distribs_comparison.pdf'):

	fig = plt.figure(figsize=(14, 12))
	
	ax = plt.subplot(1, 1, 1)

	ymax = plot_distributions_lin(Delta_Neff,z_NR,sigma_min=0.04,sigma_max=1.5,n_sigma=200)
	# Plot settings

	set_xy_lims(xmin=1e-3, xmax=1e4, ymin=0, ymax=50)
	set_xy_scales(xscale='log', yscale='linear')

	#ax.set_xlabel(r"$q_\nu \equiv p_\nu/T_{\nu, 0}$")
	ax.set_xlabel(r"$\xi_\chi \equiv q_\chi/T_\chi^0$")
	ax.set_ylabel(r"$\frac{\mathrm{d}\rho^\mathrm{NR}_\chi}{\mathrm{d} \log \xi_\chi} \, \mathrm{[eV\ cm}^{-3}\mathrm{]}$")	#ax.set_ylabel(r"$\frac{\mathrm{d}\rho^\mathrm{NR}_\nu}{\mathrm{d} \log q_\nu} \, \mathrm{[eV\ cm}^{-3}\mathrm{]}$")

	# add_cosmo_cases()

	log(f'Saving: {plot_dir + figname}')
	plt.savefig(plot_dir + figname)

def chi_scaling(plot_dir='plots/', figname='chi_scaling.pdf'):
	fig = plt.figure(figsize=(8, 5))
	ax = plt.subplot(1, 1, 1)

	Delta_Neff = 0.3
	z_NR = 1e3
	sigma = None

	plot_scaling(Delta_Neff,z_NR,case='FD',sigma=sigma)
	
	ax.set_xlabel(r'$a$')
	ax.set_ylabel(r'$R(a)$')

	log(f'Saving: {plot_dir + figname}')
	plt.savefig(plot_dir + figname)
	

if __name__ == '__main__':
	log('Starting main.py')
	
	# test_CLASS(case='FD',Delta_Neff=0.3,z_NR=1e3)
	# chi_scaling()
	# scaling = input('[main.py] Scaling (y/n): ')
	density = input('[main.py] Densities (y/n): ')
	cls = input('[main.py] Cls (y/n): ')
	dists = input('[main.py] Distributions (y/n): ')
	evols = input('[main.py] Evolutions (y/n): ')
	# if scaling == 'y':
	# 	chi_scaling()
	if density == 'y':
		cosmo_densities()
	if cls == 'y':
		vs_LCDM_Cl()
		
	# Delta_Neff = 0.3
	# z_NR = 1e3
	# sigma_min = 0.04
	# sigma_max = 1.5
	# n_sigma = 2
	# if dists == 'y':
	# 	distribs_comparison()
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
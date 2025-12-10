from utils import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({
	"text.usetex": True,
	"font.family": "serif",
})


def cosmo_comparison(plot_dir='plots/', figname='cosmo_comparison.pdf'):
	usedef = input('[main.py] Use default LiMR parameters? (y/n): ')
	if usedef == 'n':
		Delta_Neff = float(input('[main.py] Desired Delta_Neff: '))
		z_NR = float(input('[main.py] Desired z_NR: '))
	else:
		Delta_Neff = 0.3
		z_NR = 1000

	print('[main.py] Using LiMR parameters: Delta_Neff =', Delta_Neff, ', z_NR =', z_NR)

	fig = plt.figure(figsize=(14, 12))
	
	ax = plt.subplot(2, 2, 1)

	ymax = plot_distributions(Delta_Neff,z_NR)
	# Plot settings

	set_xy_lims(xmin=1e-3, xmax=1e4, ymin=1e-4, ymax=1e3)
	set_xy_scales(xscale='log', yscale='log')

	#ax.set_xlabel(r"$q_\nu \equiv p_\nu/T_{\nu, 0}$")
	ax.set_xlabel(r"$\xi_\chi \equiv q_\chi/T_\chi^0$")
	ax.set_ylabel(r"$\frac{\mathrm{d}\rho^\mathrm{NR}_\chi}{\mathrm{d} \log \xi_\chi} \, \mathrm{[eV\ cm}^{-3}\mathrm{]}$")
	#ax.set_ylabel(r"$\frac{\mathrm{d}\rho^\mathrm{NR}_\nu}{\mathrm{d} \log q_\nu} \, \mathrm{[eV\ cm}^{-3}\mathrm{]}$")

	# add_cosmo_cases()

	# ax = plt.subplot(2, 2, 2)

	# plot_energy_evolution()
	# set_xy_lims(xmin=1e0, xmax=1e4, ymin=0.95, ymax=1.30)
	# set_xy_scales(xscale='log', yscale='linear')

	# ax.set_xlabel(r"$z$")
	# ax.set_ylabel(r"$(\rho_\nu+\rho_{\rm DR})/\rho_\nu|_{\Lambda {\rm CDM}}$")

	# ax = plt.subplot(2, 2, 3)

	# plot_distribution_Cl()
	# set_xy_lims(xmin=2., xmax=2500, ymin=-0.004, ymax=0.004)
	# set_xy_scales(xscale='log', yscale='linear')

	# ax.set_xlabel(r"$\ell$")
	# ax.set_ylabel(r"$(C^{\rm TT}_\ell-C^{\rm TT}_\ell|_{\Lambda {\rm CDM}})/C^{\rm TT}_\ell|_{\Lambda {\rm CDM}}$")

	# ax = plt.subplot(2, 2, 4)

	# plot_rho_hist(case='LCDM')
	# plot_rho_hist(case='LEDR')
	# plot_rho_hist(case='HE', nratio=0.1)
	# plot_rho_hist(case='HEDR', nratio=0.1)
	# plot_rho_hist(case='LTM')
	# set_xy_lims(xmin=0., xmax=50., ymin=0., ymax=1.)
	# set_xy_scales(xscale='linear', yscale='linear')

	# ax.set_xlabel(r"$\rho_{\nu,0}^\mathrm{NR}\,\mathrm{[eV\ cm}^{-3}\mathrm{]}$")
	# ax.set_ylabel(r"$\mathrm{Posterior\ Density}$")

	# plt.text(0.13, 0.92, r"$\textbf{Data:\ }\mathrm{Planck\ TTTEEE+lowE\ +\ BAO}$",
	# 	transform=plt.gca().transAxes,
	# 	fontsize=14,
	# 	color='k',
	# 	rotation=0)
	# plt.text(0.13, 0.83, r'$\rho_{\nu,0}^\mathrm{NR} = (\sum m_\nu) \times 2 T_{\nu,0}^3 \int{\frac{\mathrm{d} q_\nu}{2\pi^2} \, q_\nu^2 f_\nu(q_\nu)}$',
	# 	transform=plt.gca().transAxes,
	# 	fontsize=14,
	# 	color='k')
	rcParams.update({
		"text.usetex": True,
	})
	print('[main.py] Saving:', plot_dir + figname)
	plt.savefig(plot_dir + figname)

if __name__ == '__main__':
	print('[main.py] Starting main.py')

	# fill_LiMR_parameters(Delta_Neff, z_NR)
	cosmo_comparison()
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
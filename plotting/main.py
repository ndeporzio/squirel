from utils import *
from matplotlib import rcParams

# global rcParams
rcParams.update({
	"text.usetex": True,
	"font.family": "serif",
	"xtick.direction": "in",
	"ytick.direction": "in",
	"legend.frameon": False,
	"axes.linewidth": 1,
	"axes.labelsize": 16,
	"xtick.labelsize": 12,
	"ytick.labelsize": 12,
    "axes.axisbelow": False,
})

# specfically for non-getdist plots
rcparams_nongetdist = {
	"xtick.major.size": 6,
	"ytick.major.size": 6,
	"xtick.minor.size": 3,
	"ytick.minor.size": 3,
	"xtick.minor.visible": True,
	"ytick.minor.visible": True,
    "axes.axisbelow": False,
}

# Figure 1
@mpl.rc_context(rcparams_nongetdist)
def cosmo_densities(plot_dir='plots/', figname='cosmo_densities.pdf'):

	z_NRs = np.logspace(3,5,3)
	Delta_Neffs = [0.05, 0.1, 0.15][::-1]

	fig = plt.figure(figsize=(8, 5))
	ax = plt.subplot(1, 1, 1)

	plot_cosmo_densities(Delta_Neffs,z_NRs,N_mnu=1,M_mnu=0.06,head_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/cosmos/'))

	plt.xlabel(r'$\mathrm{Scale\,\,Factor}\,\,a/a_0 = (1 + z)^{-1}$', fontsize=14)
	plt.ylabel(r'$\mathrm{Density\,\,Parameter}\,\,\Omega_{i}(a)$', fontsize=14)
	plt.xlim(1e-6, 1e0)
	plt.ylim(1e-4, 2e0)
	fig.tight_layout()

	save_fig(figname, plot_dir)

# Figure 2
@mpl.rc_context(rcparams_nongetdist)
def subplots_vs_LCDM_residuals(plot_dir='plots/', figname='residuals_vs_LCDM.pdf',showpeaks=False):
	Delta_Neffs = np.linspace(0.05,0.25,5)
	z_NRs = np.logspace(3,5,3)[0:2]

	fig = plt.figure(figsize=(10, 4))
	# showpeaks = True

	lensed = True
	spectras = ['tt']
	fixed1 = 'h'
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

			f_LiMRs = plot_vs_LCDM_residuals(Delta_Neffs[:viable_idx],fixed_zNRs[:viable_idx],N_mnu=1,M_mnu=0.06,lensed=lensed,spectra=spectra, fixed1=fixed1,fixed2=fixed2,head_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/cosmos/'),showpeaks=showpeaks)
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

# Figure 3
@mpl.rc_context(rcparams_nongetdist)
def full_distribs_comparison(Delta_Neff=0.2,z_NR=1e3, plot_dir='plots/', figname='distribs_comparison.pdf'):
	fig = plt.figure(figsize=(6, 12))
	
	ax = plt.subplot(3, 1, 1)
	plot_distributions(ax,Delta_Neff,z_NR,sigma_array=sigma_array_analytic)
	
	set_xy_lims(xmin=1e-4, xmax=1e4, ymin=0, ymax=50)
	set_xy_scales(xscale='log', yscale='linear')

	ax.set_xlabel(r"$\xi \equiv q/q_c$",fontsize=16)
	ax.set_ylabel(r"$\mathrm{d}\rho^\mathrm{NR}_\chi/\mathrm{d} \log \xi \; \mathrm{[eV\ cm}^{-3}\mathrm{]}$",fontsize=16)	
	ax.text(0.075,0.425,r'$\Delta N_{\mathrm{eff}}={'+f'{Delta_Neff:.{2}f}'+'}$',fontsize=14,transform=ax.transAxes,color='k')
	ax.text(0.075,0.35,r'$z_{\mathrm{NR}}=10^{'+str(np.log10(z_NR))+'}$',fontsize=14,transform=ax.transAxes,color='k')

	# Add text labels for each line at top-left, descending
	line_labels = [
		(r'$\mathrm{FD}$', cosmo_color('FD')),
		(r'$\mathrm{BE}$', cosmo_color('BE')),
		(r'$\mathrm{NT}$', cosmo_color('RD')),
		(r'$\mathrm{LN}\;(\sigma_{\mathrm{LN}}=0.04)$', cosmo_color('LNsharp')),
		(r'$\mathrm{LN}\;(\sigma_{\mathrm{LN}}=1.5)$', cosmo_color('LNwide')),
	]
	for i, (label, color) in enumerate(line_labels):
		ax.text(0.05, 0.90 - i*0.08, label, fontsize=13, transform=ax.transAxes, color=color, weight='bold')

	cbar_ax = inset_axes(ax, width="100%", height="100%",
                     bbox_to_anchor=(0.825, 0.35, 0.04, 0.6),
                     bbox_transform=ax.transAxes)
	sigma_array = sigma_array_analytic

	cmap = mpl.colors.LinearSegmentedColormap.from_list('sigma_cmap', [cosmo_color('LNsharp'), cosmo_color('LNwide')], N=len(sigma_array))
	bounds = [sigma_array[0], sigma_array[-1]]
	norm = mpl.colors.Normalize(vmin=bounds[0], vmax=bounds[-1])
	cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')
	cb1.set_label(r'$\sigma_\mathrm{LN}$', fontsize=16)
	cb1.set_ticks(np.concatenate(([sigma_array[0]], np.linspace(0,sigma_array[-1],6)[1:])))
	cb1.set_ticks(np.linspace(0,sigma_array[-1],31)[1:],minor=True)

	ax = plt.subplot(3, 1, 2)
	plot_evolution(Delta_Neff,z_NR,sigma_array=sigma_array_CLASS)

	set_xy_lims(xmin=0.1*z_NR, xmax=10*z_NR, ymin=1, ymax=1.7)
	set_xy_scales(xscale='log', yscale='linear')
	ax.legend(fontsize=14, frameon=False, handlelength=1.2, handleheight=1.2)

	# Add text labels for each line at top-left, descending
	for i, (label, color) in enumerate(line_labels):
		ax.text(0.05, 0.90 - i*0.08, label, fontsize=13, transform=ax.transAxes, color=color, weight='bold')

	ax.set_xlabel(r"$z$",fontsize=16)
	ax.set_ylabel(r"$\rho_\chi/\rho_\mathrm{asmpt}$",fontsize=16)
	plt.text(0.75*z_NR, 1.07, r'$z_{\mathrm{NR}}$', rotation=90, fontsize=16, color='k')


	ax = plt.subplot(3, 1, 3)
	plot_distrib_cls(Delta_Neff,z_NR,sigma_array=sigma_array_CLASS)
	
	zoom = False
	if zoom:
		set_xy_lims(xmin=2., xmax=2500, ymin=-0.008, ymax=0.0025)
	else:
		set_xy_lims(xmin=2., xmax=2500, ymin=-0.01, ymax=0.01)

	for i, (label, color) in enumerate(line_labels):
		if i<3:
			ax.text(0.05, 0.90 - i*0.08, label, fontsize=13, transform=ax.transAxes, color=color, weight='bold')
		else:
			ax.text(0.15, 0.90 - (i-3)*0.08, label, fontsize=13, transform=ax.transAxes, color=color, weight='bold')

	ax.set_xlabel(r"$\ell$")
	ax.set_ylabel(r"$(C^{\rm TT}_\ell-C^{\rm TT}_\ell|_{{\rm FD}})/C^{\rm TT}_\ell|_{ {\rm FD}}$")

	# change y ticks to scientific notation 
	fmt = ScalarFormatter(useMathText=True)
	fmt.set_powerlimits((0, 0)) 
	ax.yaxis.set_major_formatter(fmt)
	
	# fig.tight_layout(rect=[0, 0, 0.93, 1])
	fig.tight_layout()

	save_fig(figname, plot_dir)

# Figure 4
def vs_LCDM_contours(plot_dir='plots/', figname='LCDM_vs_DR_vs_LiMR_triangle.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):
	
	log('Generating comparison contours for LCDM vs DR vs FD LiMRs.')

	params = [
    'H0',
    'omega_cdm',
    'delta_Neff',
    'log10z_tr',
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
		]

	roots = [
		case + '_' + probes[0] + '_' + nu_hierachs[0] + '_' for case in cases
	][::-1] # reverse to plot FD first since it has largest contours

	vs_LCDM_labels = [
		r'$\mathrm{FD\;LiMR}$',
		r'$\mathrm{DR}$',
		r'$\Lambda\mathrm{CDM}$',
	]

	width_inch = 8
	height_inch = 8
	g = plots.get_subplot_plotter(chain_dir=chain_dir,width_inch=width_inch, subplot_size_ratio=height_inch/width_inch)
	g.settings.legend_fontsize = 16
	samples = get_samples(g, roots=roots)

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
		legend_labels = vs_LCDM_labels,
		legend_loc = 'upper right',
		label_order = -1,
	)

	for i, param in enumerate(params):
		param_ax = g.get_axes_for_params(param)
		param_ax.set_title('')
		fig_pos = param_ax.get_position()
		x_fig = fig_pos.x0 + fig_pos.width / 2
		y_fig = fig_pos.y1
		n_written = 0
		for j, s in enumerate(samples):
			stats = s.getMargeStats()
			par = stats.parWithName(param)
			if par is None:
				continue
			limit_idx = 2 if par.limits[0].limitType() == 'one tail upper limit' else 1 # print 95% upper limits but 1 sigma symmetric limits
			title_str = s.getInlineLatex(param, limit=limit_idx)

			# if symmetric, replace ± with explicit ^{+err}_{-err}
			if r'\pm' in title_str:
				parts = title_str.split(r'\pm')
				title_str = parts[0] + r'^{+' + parts[1].strip() + r'}_{-' + parts[1].strip() + r'}'
			g.fig.text(
				x_fig, y_fig + 0.005 + n_written * 0.035,
				f'${title_str}$',
				ha='center', va='bottom',
				fontsize=13,
				color=cosmo_color(roots[j].split('_')[0]),
			)
			n_written += 1	
	
	# g.export(plot_dir+figname)
	g.fig.savefig(plot_dir + figname, bbox_inches='tight', pad_inches=0.1)
	log(f'Saving: {plot_dir + figname}')
	
# Figure 5
@mpl.rc_context(rcparams_nongetdist)
def LiMR3d_posterior(plot_dir='plots/', figname='3d_posterior.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):

	log('Generating scatter plot with posterior contours.')

	params = [
		'delta_Neff',
		'log10z_tr',
		'log10f_ncdm1',
	]

	probe = 'PR3'
	nu_hierach = 'N_mnu=1_Mmnu=0.06'
	case = 'FD'
	root = case + '_' + probe + '_' + nu_hierach + '_'

	width_inch = 4
	height_inch = 3.5
	g = plots.get_single_plotter(chain_dir=chain_dir, width_inch=width_inch, ratio=height_inch/width_inch)
	g.settings.colormap_scatter = mcolors.LinearSegmentedColormap.from_list('plasma_trunc', mpl.colormaps['plasma'](np.linspace(0, 0.95, 256)))
	g.settings.colorbar_label_rotation = 90
	samples = get_samples(g, roots=[root])
	# samples[0].updateSettings({'contours': [0.95]}) # uncomment this line to plot only 95% contour
	samples.append(samples[0]) # necessary to add contour lines to scatter plot

	xmin = 0
	xmax = 0.5
	ymin = mcmc_plot_lims['log10z_tr'][0]
	ymax = mcmc_plot_lims['log10z_tr'][1]
	# Plot contour lines
	g.plot_3d(
		samples,
		params,
		color_bar=True,
		lims=[xmin, xmax, ymin, ymax],
		lws=[1],
		ls=['-'],
		colors=['k'],
	)

	cb_ax = g.last_colorbar.ax
	current_min, current_max = cb_ax.get_ylim()
	fmax = current_max
	fmin = -6

	# Plot annotations for the 1D bounds.
	x_param, y_param = 'delta_Neff', 'log10z_tr'

	ax = g.get_axes_for_params(x_param, y_param)
	bayesian_color = IBM_cscheme()[2]
	frequentist_color = IBM_cscheme()[0]
	draw_bayesian_1d_bounds = True
	draw_frequentist_1d_bounds = True
	markers = {'Bayes':'<', 'Freq':'X'}

	def draw_bounds(bound_dict, color, marker, draw_DR=False):
		for bound_name, bound_values in bound_dict.items():
			if bound_name == 'DR_bound':
				if draw_DR:
					ax.axvline(bound_values[1], c='gray', ls='--', lw=1.0, alpha=1)
				continue
			y = bound_values[0]
			x_tail = bound_values[1]
			ax.scatter([x_tail], [y], marker=marker, s=50, color=color, edgecolors='k', linewidths=0.5, zorder=5)

	if draw_bayesian_1d_bounds:
		bounds = {
			'DR_bound': [min(mcmc_plot_lims['log10z_tr']), 2.9533e-01],
			'logzNR3_bound': [3, 2.3081e-01],
			'logzNR35_bound': [3.5, 2.0805e-01],
			'logzNR4_bound': [4, 1.0345e-01],
			'logzNR45_bound': [4.5, 5.6053e-02],
			'logzNR5_bound': [5, 7.9728e-02],
		}
		draw_bounds(bounds, bayesian_color, markers['Bayes'], draw_DR=True)

	if draw_frequentist_1d_bounds:
		freq_bounds = {
			'DR_bound': [min(mcmc_plot_lims['log10z_tr']), 2.3e-01],
			'logzNR3_bound': [3, 1.5e-01],
			'logzNR35_bound': [3.5, 1.3e-01],
			'logzNR4_bound': [4, 7.0e-02],
			'logzNR45_bound': [4.5, 4.0e-02],
			'logzNR5_bound': [5, 5.0e-02],
		}
		draw_bounds(freq_bounds, frequentist_color, markers['Freq'])

	legend_elements = [Line2D([0], [0], color='k', lw=1, ls='-', label='FD posterior')]
	legend_labels = [r'$\mathrm{Posterior\; (FD)}$']

	if draw_bayesian_1d_bounds:
		legend_elements.append(Line2D([0], [0], color=bayesian_color, marker=markers['Bayes'], linestyle='None', markersize=7, markerfacecolor=bayesian_color, markeredgecolor='k', markeredgewidth=0.5))
		# legend_labels.append(r'$<95\%\;\mathrm{U.L.\;(FD_\mathrm{1D})}$')
		legend_labels.append(r'$\mathrm{Bayes.\;U.L.\;(FD_\mathrm{1D})}$')
		ax.text(0.305,2,r'$\mathrm{Bayes.\;U.L.\;(DR)}$',fontsize=11.5,color='gray',rotation=270)
	
	if draw_frequentist_1d_bounds:
		legend_elements.append(Line2D([0], [0], color=frequentist_color, marker=markers['Freq'], linestyle='None', markersize=7, markerfacecolor=frequentist_color, markeredgecolor='k', markeredgewidth=0.5))
		# legend_labels.append(r'$<95\%\;\mathrm{U.L.\;(FD_\mathrm{1D})}$')
		legend_labels.append(r'$\mathrm{Freq.\;U.L.\;(FD_\mathrm{1D})}$')

	ax.legend(
		handles=legend_elements,
		handlelength=1,
		loc='upper right',
		frameon='true',
		edgecolor='black',
		facecolor='white',
		framealpha=1,
		labels=legend_labels,
		fancybox=False,
	)

	# Fix colorbar size
	old_pos = cb_ax.get_position()
	cb_ax.remove()
	new_cb_ax = g.fig.add_axes([old_pos.x0, old_pos.y0, old_pos.width * 0.35, old_pos.height])
	g.last_colorbar.ax = new_cb_ax
	g.last_colorbar._draw_all()
	new_cb_ax.tick_params(labelsize=g.settings.axes_fontsize)
	g.last_colorbar.set_label(r'$\log_{10}\left(f_\chi\right)$', fontsize=14, rotation=270, labelpad=16)
	new_cb_ax.set_ylim(fmin, fmax)

	g.fig.savefig(plot_dir+figname, bbox_inches='tight', pad_inches=0.1)
	log(f'Saving: {plot_dir + figname}')

# Figure 6
def SquiRel_constraints(plot_dir='plots/', figname='SquiRel_constraints.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):

	log('Generating SquiRel PR3 constraint contours')

	params = [
    'delta_Neff',
    'log10z_tr',
	'sigma',
    ]

	probes = [
		'PR3',
	]
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
	]
	cases = [
		'FD',
		'LN',
		]

	roots = [
		case + '_' + probes[0] + '_' + nu_hierachs[0] + '_' for case in cases
	] 

	SquiRel_labels = [
		r'$\mathrm{FD\;LiMR}$',
		r'$\mathrm{SquiRel}$',
	]

	width_inch = 4
	height_inch = 4
	g = plots.get_subplot_plotter(chain_dir=chain_dir,width_inch=width_inch, subplot_size_ratio=height_inch/width_inch)
	g.settings.scaling_factor = 1
	g.settings.legend_fontsize = 13
	samples = get_samples(g, roots=roots)

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
			} for i, root in enumerate(roots)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cosmo_color(root.split('_')[0]),
				'lw': 1,
			} for i, root in enumerate(roots)
		],
		legend_labels = SquiRel_labels,
		legend_loc = 'upper right',
	)

	for i, param in enumerate(params):
		param_ax = g.get_axes_for_params(param)
		param_ax.set_title('')
		fig_pos = param_ax.get_position()
		x_fig = fig_pos.x0 + 0.005 if i == 1 else fig_pos.x0 + fig_pos.width / 2
		y_fig = fig_pos.y1
		n_written = 0
		for j, s in enumerate(samples[::-1]):
			stats = s.getMargeStats()
			par = stats.parWithName(param)
			if par is not None and par.limits[0].limitType() != 'none':
				limit_idx = 2 if par.limits[0].limitType() == 'one tail upper limit' else 1 # print 95% upper limits but 1 sigma symmetric limits
				title_str = s.getInlineLatex(param, limit=limit_idx)

				# if symmetric, replace ± with explicit ^{+err}_{-err}
				if r'\pm' in title_str:
					parts = title_str.split(r'\pm')
					title_str = parts[0] + r'^{+' + parts[1].strip() + r'}_{-' + parts[1].strip() + r'}'
				g.fig.text(
					x_fig, y_fig + 0.005 + n_written * 0.05,
					f'${title_str}$',
					ha='left' if i == 1 else 'center', va='bottom',
					fontsize=11,
					color=IBM_cscheme()[2*j],
				)
				n_written += 1	
	
	# g.export(plot_dir+figname)
	g.fig.savefig(plot_dir + figname, bbox_inches='tight', pad_inches=0.1)
	log(f'Saving: {plot_dir + figname}')

# Figure 7
def SquiRel_forecasts(plot_dir='plots/', figname='SquiRel_forecasts.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):

	log('Generating SquiRel forecast contours')

	cs = [
		IBM_cscheme()[2],
		# IBM_cscheme()[0],
		IBM_cscheme()[1],
		IBM_cscheme()[3],
		IBM_cscheme()[4],
	]

	params = [
    'delta_Neff',
    'log10z_tr',
	'sigma',
    ]

	probes = [
		# 'PR3',
		'S4',
	]
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
	]
	cases = [
		'LN',
		'LN_logzNR3FDfid',
		'LN_logzNR4FDfid',
		'LN_logzNR4NTfid',
		]

	roots = [case + '_' + probes[0] + '_' + nu_hierachs[0] + '_' for case in cases]

	SquiRel_labels = [
		# r'$\textit{Planck}\mathrm{\;PR3}$',
		r'$\mathrm{S4\mbox{-}like}\;(\Lambda\mathrm{CDM\;Fid.})$',
		r'$\mathrm{S4\mbox{-}like\;(LiMR1\;Fid.)}$',
		r'$\mathrm{S4\mbox{-}like\;(LiMR2\;Fid.)}$',
		r'$\mathrm{S4\mbox{-}like\;(LiMR3\;Fid.)}$',
	]

	scaling_factor = 0.75
	width_inch = 8*scaling_factor
	height_inch = 8*scaling_factor
	g = plots.get_subplot_plotter(chain_dir=chain_dir,width_inch=width_inch, subplot_size_ratio=height_inch/width_inch)
	g.settings.legend_fontsize = 16
	samples = get_samples(g, roots=roots)

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
			} for i, root in enumerate(roots)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cs[i],
				'lw': 1,
			} for i, root in enumerate(roots)
		],
		legend_labels = SquiRel_labels,
		legend_loc = 'upper right',
	)
	g.legend.set_bbox_to_anchor((1.0, 1.1))


	fiducial_points = [
		{'delta_Neff': 0.20, 'log10z_tr': 3.0, 'sigma': 0.55},
		{'delta_Neff': 0.10, 'log10z_tr': 4.0, 'sigma': 0.55},
		{'delta_Neff': 0.10, 'log10z_tr': 4.0, 'sigma': 1.12},
	]

	pairs = [
		('delta_Neff', 'log10z_tr'),
		('delta_Neff', 'sigma'),
		('log10z_tr', 'sigma'),
	]

	for x_param, y_param in pairs:
		ax = g.get_axes_for_params(x_param, y_param)
		for fid in fiducial_points:
			ax.plot(fid[x_param], fid[y_param], 
				marker='*', markersize=9, color='k', zorder=10)

	for i, param in enumerate(params):
		param_ax = g.get_axes_for_params(param)
		param_ax.set_title('')
		fig_pos = param_ax.get_position()
		x_fig = fig_pos.x0 + 0.005 if i == 1 else fig_pos.x0 + fig_pos.width / 2
		y_fig = fig_pos.y1
		n_written = 0
		for j, s in enumerate(samples[::-1]):
			stats = s.getMargeStats()
			par = stats.parWithName(param)
			if par is not None and par.limits[0].limitType() != 'none':
				limit_idx = 2 if par.limits[0].limitType() == 'one tail upper limit' else 1 # print 95% upper limits but 1 sigma symmetric limits
				title_str = s.getInlineLatex(param, limit=limit_idx)

				# if symmetric, replace ± with explicit ^{+err}_{-err}
				if r'\pm' in title_str:
					parts = title_str.split(r'\pm')
					title_str = parts[0] + r'^{+' + parts[1].strip() + r'}_{-' + parts[1].strip() + r'}'
				g.fig.text(
					x_fig, y_fig + 0.005 + n_written * 0.04,
					f'${title_str}$',
					ha='left' if i == 1 else 'center', va='bottom',
					fontsize=13,
					color=cs[::-1][j],
				)
				n_written += 1	
	
	g.fig.savefig(plot_dir + figname, bbox_inches='tight', pad_inches=0.1)
	log(f'Saving: {plot_dir + figname}')

# Figure 8
def SO_vs_S4_forecast(plot_dir='plots/', figname='SO_vs_S4_forecast.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):

	log('Generating SO vs S4 forecast contours')

	cs = [
		IBM_cscheme()[0],
		IBM_cscheme()[3],
	]

	params = [
    'delta_Neff',
    'log10z_tr',
	'sigma',
    ]

	probes = [
		'SO',
		'S4',
	]
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
	]
	cases = [
		'LN_logzNR4FDfid',
		]

	roots = [
		cases[0] + '_' + probe + '_' + nu_hierachs[0] + '_' for probe in probes
	] 

	SquiRel_labels = [
		r'$\mathrm{SO\mbox{-}like\;(LiMR2\;Fid.)}$',
		r'$\mathrm{S4\mbox{-}like\;(LiMR2\;Fid.)}$',
	]

	width_inch = 4
	height_inch = 4
	g = plots.get_subplot_plotter(chain_dir=chain_dir,width_inch=width_inch, subplot_size_ratio=height_inch/width_inch)
	g.settings.scaling_factor = 1
	g.settings.legend_fontsize = 13
	samples = get_samples(g, roots=roots)

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
			} for i, root in enumerate(roots)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cs[i],
				'lw': 1,
			} for i, root in enumerate(roots)
		],
		legend_labels = SquiRel_labels,
		legend_loc = 'upper right',
	)

	fiducial_params = {
    'delta_Neff': 0.10,
    'log10z_tr': 4.0,
    'sigma': 0.55,
	}

	pairs = [
		('delta_Neff', 'log10z_tr'),
		('delta_Neff', 'sigma'),
		('log10z_tr', 'sigma'),
	]

	for x_param, y_param in pairs:
		ax = g.get_axes_for_params(x_param, y_param)
		ax.plot(fiducial_params[x_param], fiducial_params[y_param], 
				marker='*', markersize=9, color='k', zorder=10)

	for i, param in enumerate(params):
		param_ax = g.get_axes_for_params(param)
		param_ax.set_title('')
		fig_pos = param_ax.get_position()
		x_fig = fig_pos.x0 + 0.005 if i == 1 else fig_pos.x0 + fig_pos.width / 2
		y_fig = fig_pos.y1
		n_written = 0
		for j, s in enumerate(samples[::-1]):
			stats = s.getMargeStats()
			par = stats.parWithName(param)
			if par is not None and par.limits[0].limitType() != 'none':
				limit_idx = 2 if par.limits[0].limitType() == 'one tail upper limit' else 1 # print 95% upper limits but 1 sigma symmetric limits
				title_str = s.getInlineLatex(param, limit=limit_idx)
				print(param, par.limits[0].limitType(), title_str)

				# if symmetric, replace ± with explicit ^{+err}_{-err}
				if r'\pm' in title_str:
					parts = title_str.split(r'\pm')
					title_str = parts[0] + r'^{+' + parts[1].strip() + r'}_{-' + parts[1].strip() + r'}'
				g.fig.text(
					x_fig, y_fig + 0.005 + n_written * 0.05,
					f'${title_str}$',
					ha='left' if i == 1 else 'center', va='bottom',
					fontsize=11,
					color=cs[::-1][j],
				)
				n_written += 1	
	
	# g.export(plot_dir+figname)
	g.fig.savefig(plot_dir + figname, bbox_inches='tight', pad_inches=0.1)
	log(f'Saving: {plot_dir + figname}')

# Figure 9
def pessimistic_S4_forecast(plot_dir='plots/', figname='S4_pessimistic_forecast.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):

	log('Generating best case LiMR S4 forecast contour in case SO measures LCDM')

	cs = [
		IBM_cscheme()[1],
		IBM_cscheme()[3],
	]

	params = [
    'delta_Neff',
    'log10z_tr',
	'sigma',
    ]

	probes = [
		'S4',
	]
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
	]
	cases = [
		'LN_logzNR4FDfid_SO_bestcase',
		'LN_logzNR4FDfid',
		]

	roots = [
		case + '_' + probes[0] + '_' + nu_hierachs[0] + '_' for case in cases
	] 

	SquiRel_labels = [
		r'$\mathrm{S4\mbox{-}like\;(LiMR4\;Fid.)}$',
		r'$\mathrm{S4\mbox{-}like\;(LiMR2\;Fid.)}$',
	]

	width_inch = 4
	height_inch = 4
	g = plots.get_subplot_plotter(chain_dir=chain_dir,width_inch=width_inch, subplot_size_ratio=height_inch/width_inch)
	g.settings.scaling_factor = 1
	g.settings.legend_fontsize = 13
	samples = get_samples(g, roots=roots)

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
			} for i, root in enumerate(roots)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cs[i],
				'lw': 1,
			} for i, root in enumerate(roots)
		],
		legend_labels = SquiRel_labels,
		legend_loc = 'upper right',
	)

	fiducial_points = [
		{'delta_Neff': 0.10, 'log10z_tr': 4.0, 'sigma': 0.55},
		{'delta_Neff': 0.03, 'log10z_tr': 4.0, 'sigma': 0.55},
	]

	pairs = [
		('delta_Neff', 'log10z_tr'),
		('delta_Neff', 'sigma'),
		('log10z_tr', 'sigma'),
	]

	for x_param, y_param in pairs:
		ax = g.get_axes_for_params(x_param, y_param)
		for fid in fiducial_points:
			ax.plot(fid[x_param], fid[y_param], 
				marker='*', markersize=9, color='k', zorder=10)

	for i, param in enumerate(params):
		param_ax = g.get_axes_for_params(param)
		param_ax.set_title('')
		fig_pos = param_ax.get_position()
		x_fig = fig_pos.x0 + 0.005 if i == 1 else fig_pos.x0 + fig_pos.width / 2
		y_fig = fig_pos.y1
		n_written = 0
		for j, s in enumerate(samples[::-1]):
			stats = s.getMargeStats()
			par = stats.parWithName(param)
			if par is not None and par.limits[0].limitType() != 'none':
				limit_idx = 2 if par.limits[0].limitType() == 'one tail upper limit' else 1 # print 95% upper limits but 1 sigma symmetric limits
				if param != 'sigma' or j == 0:
					title_str = s.getInlineLatex(param, limit=limit_idx)
				else:
					title_str = r'\sigma_\mathrm{LN} < 1.39' # had to hardcode in as 2sigma< instead of 95%< since montepython won't compute latter

				# if symmetric, replace ± with explicit ^{+err}_{-err}
				if r'\pm' in title_str:
					parts = title_str.split(r'\pm')
					title_str = parts[0] + r'^{+' + parts[1].strip() + r'}_{-' + parts[1].strip() + r'}'
				g.fig.text(
					x_fig, y_fig + 0.005 + n_written * 0.05,
					f'${title_str}$',
					ha='left' if i == 1 else 'center', va='bottom',
					fontsize=11,
					color=cs[::-1][j],
				)
				n_written += 1	
	
	# g.export(plot_dir+figname)
	g.fig.savefig(plot_dir + figname, bbox_inches='tight', pad_inches=0.1)
	log(f'Saving: {plot_dir + figname}')

# Figure 10
@mpl.rc_context(rcparams_nongetdist)
def DNeff_profiles(plot_dir='plots/', figname='DNeff_profiles.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):
	
	log('Comparing profile likelihoods for DR and FD1D LiMRs.')
	
	fig = plt.figure(figsize=(10, 4))

	ax = plt.subplot(1, 2, 1)

	ymin = -4
	plot_extrap_profs(ymin, chain_dir)

	set_xy_lims(xmin=-0.4, xmax=0.25, ymin=ymin, ymax=6)

	ax.set_xlabel(r'$\Delta N_{\rm eff}$')
	ax.set_ylabel(r'$\Delta \chi^2$')

	cmap = mpl.colors.ListedColormap(IBM_cscheme())
	class _FDHandle:
		pass

	legend_elements = [
		Line2D([0], [0], color='k', lw=1.5, label=r'$\rm DR$'),
		(_FDHandle(), r'$\rm FD_{1D}$'),
	]

	ax.legend(
		[Line2D([0], [0], color='k', lw=1.5), _FDHandle()],
		[r'$\rm DR$', r'$\rm FD_{1D}$'],
		handler_map={_FDHandle: HandlerColorbar(cmap)},
		fontsize=14, handlelength=1.2, loc='lower right'
	)

	ax = plt.subplot(1, 2, 2)
	plot_fc_corrected(chain_dir)

	set_xy_lims(xmin=0, xmax=0.25, ymin=0, ymax=7)

	ax.set_xlabel(r'$\Delta N_{\rm eff}$')
	# ax.set_ylabel(r'$\Delta \chi^2$')

	ax.legend(
		[Line2D([0], [0], color='k', lw=1.5), _FDHandle()],
		[r'$\rm DR$', r'$\rm FD_{1D}$'],
		handler_map={_FDHandle: HandlerColorbar(cmap)},
		fontsize=14, handlelength=1.2, loc='lower right'
	)

	fig.tight_layout(rect=[0, 0, 0.92, 1])
	# fig.subplots_adjust(wspace=0.0)

	log10zNRs = np.linspace(3, 5, 5)
	cbar_ax = fig.add_axes([plt.gca().get_position().x1 + 0.02, plt.gca().get_position().y0, 0.02, plt.gca().get_position().height])
	spacing = (log10zNRs[-1]-log10zNRs[0])/(len(log10zNRs)-1)
	bounds = [log10zNRs[0]-spacing/2, log10zNRs[-1]+spacing/2]
	norm = mpl.colors.Normalize(vmin=bounds[0], vmax=bounds[-1])
	cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')
	cb1.set_label(r'$\log_{10}\left(z_\mathrm{NR}\right)$', fontsize=16)
	cb1.set_ticks(log10zNRs)

	save_fig(figname, plot_dir)

# Figure 11
def dataset_triangle(plot_dir='plots/', figname='full_datasets_triangle.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):
	
	log('Generating full comparison contours for different datasets and analysis choices.')

	cs = [
		'#dc267f',
		# '#785ef0',
		'#648fff',
		# '#ffb000',
		'#fe6100',
	]

	params = [
    'H0',
    'raw_omega_b',
    'omega_cdm',
    'ln10^{10}A_s',
    'tau_reio',
    'n_s',
    'delta_Neff',
    'log10z_tr',
    ]

	probes = [
		'PR3',
		# 'PR4',
		# 'PR3+eBOSS',
		# 'PR3+lens',
		'PR4',

	]
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
		'N_mnu=2_Mmnu=0.11',
	]
	case = 'FD'

	roots = [
		case + '_' + probe + '_' + nu_hierachs[0] + '_' for probe in probes
	]
	roots.append(case + '_' + probes[0] + '_' + nu_hierachs[1] + '_' )

	dataset_labels = [
		'PR3',
		'PR4',
		# 'PR3+BAO',
		# 'PR3+lens.',
		r'$\mathrm{PR3}\;(\Sigma m_\nu = 0.11\,\mathrm{eV})$',
	]

	width_inch = 8
	height_inch = 8
	g = plots.get_subplot_plotter(chain_dir=chain_dir,width_inch=width_inch, subplot_size_ratio=height_inch/width_inch)
	g.settings.legend_fontsize = 16
	samples = get_samples(g, roots=roots)

	g.triangle_plot(
		samples,
		params,
		filled=True,
		param_limits=mcmc_plot_lims,
		contour_args = [
			{
				# 'alpha':0.1,
				# 'alpha':np.linspace(1,0.1,10)[i],
				'color': cs[i],
				'lw': 2,
			} for i, root in enumerate(roots)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': cs[i],
				'lw': 1,
			} for i, root in enumerate(roots)
		],
		legend_labels = dataset_labels,
		legend_loc = 'upper right',
	)

	g.export(plot_dir + figname)
	log(f'Saving: {plot_dir + figname}')

if __name__ == '__main__':
	log('Starting main.py')

	density = input('[main.py] Fig. 1 (y/n): ')
	gen_cls = input('[main.py] Fig. 2 (y/n): ')
	dists = input('[main.py] Fig. 3 (y/n): ')
	LiMR_triangle = input('[main.py] Fig. 4 (y/n): ')
	color_posterior = input('[main.py] Fig. 5 (y/n): ')
	squirel_constraints = input('[main.py] Fig. 6 (y/n): ')
	squirel_forecasts = input('[main.py] Fig. 7 (y/n): ')
	so_forecast = input('[main.py] Fig. 8 (y/n): ')
	pessimistic_forecast = input('[main.py] Fig. 9 (y/n): ')
	profiles = input('[main.py] Fig. 10 (y/n): ')
	compare_datasets = input('[main.py] Fig. 11 (y/n): ')


	if density == 'y':
		cosmo_densities()
	if gen_cls == 'y':
		subplots_vs_LCDM_residuals()
		
	Delta_Neff = 0.10
	z_NR = pow(10,4)
	sigma_min = 0.04
	sigma_max = 1.5
	sigma_array_CLASS = np.concatenate((np.array([0.04]), np.linspace(0.3,1.5,9)))
	n_sigma_CLASS = len(sigma_array_CLASS)
	n_sigma_analytic = 293 # so that steps size is 0.005
	sigma_array_analytic = np.linspace(sigma_min,sigma_max,n_sigma_analytic)
	if dists == 'y':
		full_distribs_comparison(Delta_Neff=Delta_Neff,z_NR=z_NR)
	if LiMR_triangle == 'y':
		vs_LCDM_contours()
	if color_posterior == 'y':
		LiMR3d_posterior()
	if profiles == 'y':
		DNeff_profiles()
	if squirel_constraints == 'y':
		SquiRel_constraints()
	if squirel_forecasts == 'y':
		SquiRel_forecasts()
	if so_forecast == 'y':
		SO_vs_S4_forecast()
	if pessimistic_forecast == 'y':
		pessimistic_S4_forecast()
	if profiles == 'y':
		DNeff_profiles()
	if compare_datasets == 'y':
		dataset_triangle()
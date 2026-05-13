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
})

# specfically for non-getdist plots
rcparams_nongetdist = {
	"xtick.major.size": 6,
	"ytick.major.size": 6,
	"xtick.minor.size": 3,
	"ytick.minor.size": 3,
	"xtick.minor.visible": True,
	"ytick.minor.visible": True,
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

	legend_elements = [
		Line2D([0], [0], color=cosmo_color('FD'), lw=1.5, label=r'$\mathrm{FD}$'),
		Line2D([0], [0], color=cosmo_color('BE'), lw=1.5, label=r'$\mathrm{BE}$'),
		Line2D([0], [0], color=cosmo_color('RD'), lw=1.5, label=r'$\mathrm{NT}$'),
		Line2D([0], [0], color=cosmo_color('LNsharp'), lw=1.5, label=r'$\mathrm{LN}\;(\sigma_{\mathrm{LN}}=0.04)$'),
		Line2D([0], [0], color=cosmo_color('LNwide'), lw=1.5, label=r'$\mathrm{LN}\;(\sigma_{\mathrm{LN}}=1.5)$'),		
	]
	ax.legend(handles=legend_elements, fontsize=14, frameon=False, ncols=1, handlelength=1.2, loc='upper left')

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
	# cbar_ax.tick_params(zorder=10)
	# sigma_markers = np.linspace(0,sigma_array[-1],6)[1:-1]
	# sigma_markers = [0.55]
	# for sigma_marker in sigma_markers:
	# 	# color = cmap(norm(sigma_marker))
	# 	color = 'gray'
	# 	cbar_ax.axhline(sigma_marker, color=color, lw=1, alpha=1)
	# cb1.solids.set_alpha(0.7)

	ax = plt.subplot(3, 1, 2)
	plot_evolution(Delta_Neff,z_NR,sigma_array=sigma_array_CLASS)

	set_xy_lims(xmin=0.1*z_NR, xmax=10*z_NR, ymin=1, ymax=1.7)
	set_xy_scales(xscale='log', yscale='linear')
	ax.legend(fontsize=14, frameon=False, handlelength=1.2, handleheight=1.2)

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
	# samples[0].updateSettings({'contours': [0.95]}) # to plot only 95% contour
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

	Dneffbound = True
	if Dneffbound:
		x_param, y_param = 'delta_Neff', 'log10z_tr'

		bounds = {
			'DR_bound': [min(mcmc_plot_lims['log10z_tr']),2.9533e-01], # from my DR contours, consistent with 2207.13133
			'logzNR3_bound': [3,2.3081e-01],
			'logzNR35_bound': [3.5,2.0805e-01],
			'logzNR4_bound': [4,1.0345e-01],
			'logzNR45_bound': [4.5,5.6053e-02],
			'logzNR5_bound': [5,7.9728e-02],
		}

		ax = g.get_axes_for_params(x_param, y_param)
		arrow_lw=2
		arrow_head_width=0.1
		arrow_base_height=0.2
		arrow_head_length=0.01
		arrow_dx = 0.03 

		for bound in bounds:
			if bound == 'DR_bound':
				ax.axvline(bounds[bound][1], c='k', ls='--', lw=1.0, alpha=0.5)
			else:
				y = bounds[bound][0]
				x_tail = bounds[bound][1]
				x_head = x_tail - arrow_dx

				# Shaft
				ax.plot([x_tail, x_head + arrow_head_length], [y, y],
						color='k', lw=arrow_lw, alpha=1, solid_capstyle='butt')

				# Base
				ax.plot([x_tail, x_tail], [y - arrow_base_height/2, y + arrow_base_height/2],
						color='k', lw=arrow_lw, alpha=1, solid_capstyle='butt')

				# Arrowhead
				head = Polygon(
					[[x_head, y],
					[x_head + arrow_head_length, y + arrow_head_width/2],
					[x_head + arrow_head_length, y - arrow_head_width/2]],
					closed=True, color='k', alpha=1, 
				)
				head.set_clip_box(ax.bbox)
				head.set_clip_on(True)
				ax.add_patch(head)
		
	class ArrowHandler(HandlerBase):
		def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
			y_mid = height / 2
			head_length = width * 0.3
			base_height = height 

			# Shaft
			shaft = Line2D(
				[width, head_length], [y_mid, y_mid],
				color='k', lw=1.5, solid_capstyle='butt', transform=trans
			)
			# Base
			base = Line2D(
				[width, width], [y_mid - base_height/2, y_mid + base_height/2],
				color='k', lw=1.5, solid_capstyle='butt', transform=trans
			)
			# Arrowhead
			head = Polygon(
				[[0, y_mid],
				[head_length, y_mid + height*0.3],
				[head_length, y_mid - height*0.3]],
				closed=True, color='k', transform=trans
			)
			return [shaft, base, head]

	arrow_proxy = Line2D([0], [0], color='k')  # dummy handle

	legend_elements = [
		Line2D([0], [0], color='k', lw=1, ls='-', label='FD posterior'),
		Line2D([0], [0], color='k', lw=1, ls='--', alpha=0.5, label='DR bound'),
		arrow_proxy,
	]

	leg_labels = [
		r'$95\%\;\mathrm{Posterior\; (FD)}$',
		r'$<95\%\;\mathrm{U.L.\;(DR)}$',
		r'$<95\%\;\mathrm{U.L.\;(FD_\mathrm{1D})}$',
	]

	ax.legend(
		handles=legend_elements,
		handlelength=1,
		loc='upper right', 
		frameon='true', 
		edgecolor='black', 
		facecolor='white', 
		framealpha=1,
		handler_map={arrow_proxy: ArrowHandler()},
		labels = leg_labels,
		fancybox=False,
		)

	# Fix colorbar size
	old_pos = cb_ax.get_position()
	cb_ax.remove()
	new_cb_ax = g.fig.add_axes([old_pos.x0, old_pos.y0, old_pos.width * 0.35, old_pos.height])
	g.last_colorbar.ax = new_cb_ax
	g.last_colorbar._draw_all()
	new_cb_ax.tick_params(labelsize=g.settings.axes_fontsize)
	g.last_colorbar.set_label(r'$\log_{10}\left(f_\chi\right)$', fontsize=14, rotation=90, labelpad=3)
	new_cb_ax.set_ylim(fmin, fmax)

	g.fig.savefig(plot_dir+figname, bbox_inches='tight', pad_inches=0.1)
	log(f'Saving: {plot_dir + figname}')


# Figure 6
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

# Figure 7
def SquiRel_contours(plot_dir='plots/', figname='SquiRel_contours.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):

	log('Generating SquiRel contours')

	params = [
    'delta_Neff',
    'log10z_tr',
	'sigma',
    ]

	probes = [
		'PR3',
		'S4',
	]
	nu_hierachs = [
		'N_mnu=1_Mmnu=0.06',
	]
	cases = [
		'LN',
		'LN_logzNR4FDfid',
		]

	roots = [
		cases[0] + '_' + probe + '_' + nu_hierachs[0] + '_' for probe in probes
	] 
	roots.append(cases[1] + '_' + probes[1] + '_' + nu_hierachs[0] + '_' )

	SquiRel_labels = [
		r'$\textit{Planck}\mathrm{\;PR3}$',
		r'$\mathrm{S4\mbox{-}like}\;(\Lambda\mathrm{CDM\;Fid.})$',
		r'$\mathrm{S4\mbox{-}like\;(LiMR\;Fid.)}$',
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
				'color': IBM_cscheme()[::-1][2*i], # invert coloring to mimic vs LCDM contours
				'lw': 2,
			} for i, root in enumerate(roots)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': IBM_cscheme()[::-1][2*i],
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
		x_fig = fig_pos.x0 + fig_pos.width / 2
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
					x_fig, y_fig + 0.005 + n_written * 0.035,
					f'${title_str}$',
					ha='center', va='bottom',
					fontsize=13,
					color=IBM_cscheme()[2*j],
				)
				n_written += 1	
	
	# g.export(plot_dir+figname)
	g.fig.savefig(plot_dir + figname, bbox_inches='tight', pad_inches=0.1)
	log(f'Saving: {plot_dir + figname}')

# Figure 8
def dataset_triangle(plot_dir='plots/', figname='full_datasets_triangle.pdf', chain_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data/chains')):
	
	log('Generating full comparison contours for different datasets and analysis choices.')

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
		'PR4',
		'PR3+eBOSS',
		'PR3+lens',
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
		'PR3+BAO',
		'PR3+lens.',
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
				'color': IBM_cscheme()[i],
				'lw': 2,
			} for i, root in enumerate(roots)
		],
		line_args = [
			{
				'alpha': 1,
				'zorder': -1,
				'color': IBM_cscheme()[i],
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

	density = input('[main.py] Densities (y/n): ')
	gen_cls = input('[main.py] LCDM vs DR vs LiMR Cls (y/n): ')
	dists = input('[main.py] Distributions (y/n): ')
	LiMR_triangle = input('[main.py] LCDM vs DR vs LiMR triangle (y/n): ')
	color_posterior = input('[main.py] Scatter plot comparison (y/n): ')
	profiles = input('[main.py] Profiles (y/n): ')
	SquiRel_triangle = input('[main.py] SquiRel contours (y/n): ')
	compare_datasets = input('[main.py] Dataset contours (y/n): ')

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
	if SquiRel_triangle == 'y':
		SquiRel_contours()
	if compare_datasets == 'y':
		dataset_triangle()
#! /usr/bin/env python3

from csv import DictReader
from york import YorkReg
from ancova import ancova
from pylab import *
from matplotlib import ticker
from scipy.interpolate import interp1d
from scipy.stats import norm, chi2

from matplotlib.patches import Ellipse
from matplotlib.legend_handler import HandlerTuple

style.use('../mydefault.mplstyle')

from matplotlib.markers import MarkerStyle
hdiamond = MarkerStyle("d")
hdiamond._transform.rotate_deg(90)

style.use('../mydefault.mplstyle')

f95 = norm.ppf(0.975)

ANDERSON_COLOR = (.45, .65, 1.)
FIEBIG_COLOR = (.95, .5, .1)
DVHLGB_COLOR = (.6,.9,0)

def cov_ellipse(cov, q = .95):
	"""
	Parameters
	----------
	cov : (2, 2) array
		Covariance matrix.
	q : float
		Confidence level, should be in (0, 1)

	Returns
	-------
	width, height, rotation :
		 The lengths of two axises and the rotation angle in degree
	for the ellipse.
	"""

	r2 = chi2.ppf(q, 2)
	val, vec = linalg.eigh(cov)
	width, height = 2 * sqrt(val[:, None] * r2)
	rotation = degrees(arctan2(*vec[::-1, 0]))

	return width, height, rotation

def reformulate_calibration(A, B, CM):
	Xbary = -CM[0,1] / CM[0,0]
	Tbary = Xbary**-.5 - 273.15
	Bbary = A * Xbary + B
	sBbary = (CM[0,0] * Xbary + 2 * CM[0,1] * Xbary + CM[1,1])**.5
	return Tbary, A, Bbary, sBbary

with open('../1_compile_D47_data/piasecki_sites.csv') as f:
	piasecki_sites = [{k: r[k] for k in r} for r in DictReader(f)]
piasecki_sites = {r['Site']: r for r in piasecki_sites}

with open('benthic_micro_habitats.csv') as f:
	habitats = [{k: r[k] for k in r} for r in DictReader(f)]
habitats = {_['Genus']: _['Micro-habitat'] for _ in habitats}

with open('../7_assign_Tcalcif/foram_D47_calibration_data.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in data:
	tobepopped = []
	for k in r:
		if k not in ['Sample', 'Species', 'Type', 'Site', 'Ref', 'Twoa23_vs_Tiso_species', 'Twoa23_500m_vs_Tiso_species_500m', 'Twoa23_1500m_vs_Tiso_species_1500m', 'Tiso_species_offset']:
			if r[k] == '':
				tobepopped.append(k)
			else:
				if k == 'Depth' and '-' in r[k]:
					r['Depth'] = tuple(float(_) for _ in r[k][1:-1].split('-'))
				else:
					r[k] = float(r[k])
	for k in tobepopped:
		r.pop(k)


with open('../8_seawater_chemistry/bottom_seawater_chemistry.csv') as f:
	bswchem = list(DictReader(f))
for r in bswchem:
	for k in r:
		if k not in ['Site']:
			r[k] = float(r[k])
bswchem = {r['Site']: r for r in bswchem}


with open('equilibrium_data.csv') as f:
	eqdata = list(DictReader(f))
	for r in eqdata:
		for k in r:
			if k not in ['Sample', 'Lab', 'Ref']:
				r[k] = float(r[k])

# MARK Foram regressions


kw_errorbar_foram = dict(
	ls = 'None',
	marker = 'None',
	ecolor = 'k',
	alpha = 0.15,
	capsize = 0,
	elinewidth = 1,
	zorder = 100,
	)

kw_plot_foram = dict(
	ls = 'None',
	mew = .8,
	mec = 'k',
	)

kw_lsce = dict(
	mfc = [.7]*3,
	zorder = 200,
	)

kw_bergen = dict(
	mfc = [1]*3,
	zorder = 100,
	)

kw_plot_benthic = kw_plot_foram | dict(
	marker = 's',
	ms = 5,
	)

kw_plot_planktic = kw_plot_foram | dict(
	marker = 'o',
	ms = 5,
	)

kw_yorkreg = dict(
	lw = 0.8,
	color = 'k',
	zorder = 15,
	dashes = (6,3,2,3),
	)

kw_lgbdvh = dict(
	ls = 'None',
	mfc = DVHLGB_COLOR,
	mec = 'k',
	mew = 0.8,
	zorder = 200,
	ms = 7,
	)

regressions = {}

for r in data:
	r['nan'] = nan
	r['SE_nan'] = nan

for figtitle, TWOA_FOR_PLANKTICS_BELOW_ZERO, TWOA_FOR_PLANKTICS_WITHOUT_OFFSET in [
	('1_Regression_concordant_planktics', False, False),
	('2_Regression_planktics', False, False),
	('4_Regression_KON97', False, False),
	('3_Regression_full', False, False),
	('3_Regression_full', 'Twoa23_1500m', False),
	('3_Regression_full', 'Twoa23_1500m', 'nan'),
	]:

	print(f'Preparing {figtitle}...')

	fig = figure(figsize = (5,5))
	subplots_adjust(.15, .1, .95, .95)
	
	show_benthics = 'planktic' not in figtitle
	show_discordant = 'concordant' not in figtitle
	Tiso_field = 'Tiso_KON97' if 'KON97' in figtitle else 'Tiso_species'
	
	if figtitle == '1_Regression_concordant_planktics':
		leg_title_pforams, = plot([], [], ls = 'None', marker = 'None',
		label = r'$\hspace{-3.3}\mathbf{Concordant\ planktic\ foraminifera\ only:}$',
		)
	elif figtitle == '2_Regression_planktics':
		leg_title_pforams, = plot([], [], ls = 'None', marker = 'None',
		label = r'$\hspace{-3.3}\mathbf{All\ planktic\ foraminifera:}$',
		)
	else:
		leg_title_pforams, = plot([], [], ls = 'None', marker = 'None',
		label = r'$\hspace{-2.5}\mathbf{Planktic\ foraminifera:}$',
		)

	legs = [leg_title_pforams]

	'''PERAL PLANKTICS'''
	G = [
		r for r in data 
		if 'Peral' in r['Ref']
		and r['Type'] == 'planktic'
		and (show_discordant or r['Twoa23_vs_Tiso_species'] == '')
		]

	X = [(r[Tiso_field]+273.15)**-2 for r in G]
	eX = [f95 * 2 * r['SE_'+Tiso_field] * (r[Tiso_field]+273.15)**-3 for r in G]
	Y = [r['D47_ICDES'] for r in G]
	eY = [r['SE_D47_ICDES']*f95 for r in G]

	if TWOA_FOR_PLANKTICS_BELOW_ZERO:
		for k,r in enumerate(G):
			if r['Tiso_species'] < 0:
				X[k] = (r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-2
				eX[k] = f95 * 2 * r['SE_'+TWOA_FOR_PLANKTICS_BELOW_ZERO] * (r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-3

	if TWOA_FOR_PLANKTICS_WITHOUT_OFFSET:
		for k,r in enumerate(G):
			if r['Tiso_species_offset'] == 'all planktics':
				X[k] = (r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-2
				eX[k] = f95 * 2 * r['SE_'+TWOA_FOR_PLANKTICS_WITHOUT_OFFSET] * (r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-3

	peral_planktic_data = (X[:], Y[:], [_/f95 for _ in eX], [_/f95 for _ in eY])

	errorbar(X, Y, eY, eX, **kw_errorbar_foram)
	leg_peralp, = plot(X, Y, label = 'Peral et al. (2018)', **kw_plot_planktic, **kw_lsce)
	legs += [leg_peralp]
	
	'''MEINICKE'''
	G = [
		r for r in data 
		if 'Meinicke' in r['Ref']
		and r['Type'] == 'planktic'
		and (show_discordant or r['Twoa23_vs_Tiso_species'] == '')
		]

	X = [(r[Tiso_field]+273.15)**-2 for r in G]
	eX = [f95 * 2 * r['SE_'+Tiso_field] * (r[Tiso_field]+273.15)**-3 for r in G]
	Y = [r['D47_ICDES'] for r in G]
	eY = [r['SE_D47_ICDES']*f95 for r in G]	

	if TWOA_FOR_PLANKTICS_BELOW_ZERO:
		for k,r in enumerate(G):
			if r['Tiso_species'] < 0:
				X[k] = (r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-2
				eX[k] = f95 * 2 * r['SE_'+TWOA_FOR_PLANKTICS_BELOW_ZERO] * (r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-3

	if TWOA_FOR_PLANKTICS_WITHOUT_OFFSET:
		for k,r in enumerate(G):
			if r['Tiso_species_offset'] == 'all planktics':
				X[k] = (r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-2
				eX[k] = f95 * 2 * r['SE_'+TWOA_FOR_PLANKTICS_WITHOUT_OFFSET] * (r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-3

	meinicke_data = (X[:], Y[:], [_/f95 for _ in eX], [_/f95 for _ in eY])
	errorbar(X, Y, eY, eX, **kw_errorbar_foram)
	leg_meinicke, = plot(X, Y, label = 'Meinicke et al. (2020)', **kw_plot_planktic, **kw_bergen)
	legs += [leg_meinicke]
	
	if figtitle == '1_Regression_concordant_planktics':
		concordant_planktic_data = list(zip(*(
			list(zip(*peral_planktic_data))
			+ list(zip(*meinicke_data))
			)))

		ps, pi = ancova(*peral_planktic_data, *zip(*[_ for _ in zip(*meinicke_data) if not isnan(_[0])]))

		with open('../output/saved_values.py', 'a') as fid:
			fid.write(f'''
concordant_planktic_ancova_p_slope = {ps}
concordant_planktic_ancova_p_intercept = {pi}
''')

		if min(ps, pi) >= 0.05:
			print(f'  The two concordant planktic data sets are statistically indistiguishable (p_slope = {ps:.3f}, p_intercept = {pi:.2f})')
		else:
			print(f'  The two concordant planktic data sets are statistically distiguishable (p_slope = {ps:.3f}, p_intercept = {pi:.2f})')

	if show_benthics:

		leg_title_bforams, = plot([], [], ls = 'None', marker = 'None',
			label = r'$\hspace{-2.5}\mathbf{Benthic\ foraminifera:}$',
			)
		legs += [leg_title_bforams]

		'''PERAL BENTHICS''' # uses Twoa23
		G = [r for r in data if 'Peral' in r['Ref'] and r['Type'] == 'benthic']
		X = [(r['Twoa23']+273.15)**-2 for r in G]
		eX = [f95*2*r['SE_Twoa23']*(r['Twoa23']+273.15)**-3 for r in G]
		Y = [r['D47_ICDES'] for r in G]
		eY = [r['SE_D47_ICDES']*f95 for r in G]	
		peral_benthic_data = (X[:], Y[:], [_/f95 for _ in eX], [_/f95 for _ in eY])
		errorbar(X, Y, eY, eX, **kw_errorbar_foram)
		leg_peralb, = plot(X, Y, label = 'Peral et al. (2018)', **(kw_plot_benthic | kw_lsce | dict(mfc = [.5]*3)))
		legs += [leg_peralb]

		'''PIASECKI'''
		G = [r for r in data if 'Piasecki' in r['Ref'] and r['Type'] == 'benthic']

		X, Y, eX, eY = [], [], [], []
		for site in sorted({r['Site'] for r in G}):

			if (
				piasecki_sites[site]['Ref'] not in ['O.A.', 'World Ocean Atlas 2009 volume 1: Temperature (2010)']
				):
				Tfield = 'Tpub'
			else:
				Tfield = 'Twoa23'

			H = [r for r in G if r['Site'] == site]
			X.append((H[0][Tfield]+273.15)**-2)
			eX.append(f95*2*H[0][f'SE_{Tfield}']*(H[0][Tfield]+273.15)**-3)
			D47 = sum([r['D47_ICDES']/r['SE_D47_ICDES']**2 for r in H]) / sum([1/r['SE_D47_ICDES']**2 for r in H])
			sD47 = sum([1/r['SE_D47_ICDES']**2 for r in H])**-.5
			Y.append(D47)
			eY.append(sD47*f95)

		piasecki_data = (X[:], Y[:], [_/f95 for _ in eX], [_/f95 for _ in eY])
		errorbar(X, Y, eY, eX, **kw_errorbar_foram)
		leg_piasecki, = plot(X, Y, label = 'Piasecki et al. (2019)', **(kw_plot_benthic | kw_bergen | dict(mfc = 'k')))
		legs += [leg_piasecki]
	else:
		piasecki_data, peral_benthic_data = [], []



	'''OTHER SAMPLES'''
	G = [_ for _ in eqdata if 'Anderson' in _['Ref']]
	X = [(r['T']+273.15)**-2 for r in G]
	eX = [f95 * 2 * r['SE_T'] * (r['T']+273.15)**-3 for r in G]
	Y = [r['D47_ICDES'] for r in G]
	eY = [r['SE_D47_ICDES']*f95 for r in G]	
	errorbar(X, Y, eY, eX, **kw_errorbar_foram)
	leg_eq_lsce, = plot(X, Y, label = 'LSCE (Anderson et al., 2021)', **(kw_lgbdvh | dict(marker = (3,0,0))))

	G = [_ for _ in eqdata if 'Fiebig' in _['Ref']]
	X = [(r['T']+273.15)**-2 for r in G]
	eX = [f95 * 2 * r['SE_T'] * (r['T']+273.15)**-3 for r in G]
	Y = [r['D47_ICDES'] for r in G]
	eY = [r['SE_D47_ICDES']*f95 for r in G]	
	errorbar(X, Y, eY, eX, **kw_errorbar_foram)
	leg_eq_gu, = plot(X, Y, label = 'Goethe Univ. (Fiebig et al., 2021)', **(kw_lgbdvh | dict(marker = (3,0,180))))

	G = [_ for _ in eqdata if 'Huyghe' in _['Ref']]
	X = [(r['T']+273.15)**-2 for r in G]
	eX = [f95 * 2 * r['SE_T'] * (r['T']+273.15)**-3 for r in G]
	Y = [r['D47_ICDES'] for r in G]
	eY = [r['SE_D47_ICDES']*f95 for r in G]	
	errorbar(X, Y, eY, eX, **kw_errorbar_foram)
	leg_eq_pecten, = plot(X, Y, label = '$\mathit{A.}$$\mathit{Colbecki}$ (Huyghe et al., 2022)', **(kw_lgbdvh | dict(marker = (4,0,45), ms = 6)))


	'''REGRESSIONS'''
	planktic_data = zip(*(
		list(zip(*peral_planktic_data))
		+ list(zip(*meinicke_data))
		))

	A, B, CM, rchisq = YorkReg(*planktic_data)
	
	if figtitle.startswith('1_'):
		print(f'  RMSWD for concordant planktics is {rchisq**.5:.1f}')
	if figtitle.startswith('3_'):
		print(f'  RMSWD for all planktics is {rchisq**.5:.1f}')
	
	xmin, xmax = 308**-2, 269**-2

	xi = linspace(xmin, xmax)
	yi = A * xi + B
	syi = (CM[0,0]*xi**2 + 2*CM[0,1]*xi + CM[1,1])**.5

	plot(xi,yi+f95*syi,'-', **kw_yorkreg)
	if figtitle.startswith('3_'):
# 		T0, A, B0, sB0 = reformulate_calibration(A, B, CM)
		label = f'\nRegression for all planktics (95 % CL)\nΔ$_{{47}}$ = {A/1000:.2f} $\\times$ $10^3$ / T$^2$ + {B:.4f}'
	else:
		label = f'York regression of planktic data (95 % CL)'
	leg_plreg, = plot(xi,yi-f95*syi,'-', label = label, **kw_yorkreg)

# 	legs += [leg_plreg]

	leg_calibs, = plot([], [], ls = 'None', marker = 'None', label = '$\mathbf{I-CDES\ calibrations:}\hspace{-3.5}$')
	
	with open('anderson_data.csv') as fid: # Anderson et al. (2021)
		anderson_data = list(DictReader(fid))

	for _ in anderson_data:
		_['D47'] = float(_['D47 I-CDES90'])
		_['SE_D47'] = float(_['SE'])
		_['T'] = float(_['Temperature'])
		_['SE_T'] = float(_['Temperature error'])
		
	_samples = [
		'LF2012-9_7-A',
		'LF2012-D1-A',
		'LJ2010-12A-Z1A',
		'LJ2010-12A-Z2A',
		'LJ2010-5B-A',
		'LV26NOV10-2A',
		]
	G = [_ for _ in anderson_data if _['Sample'] in _samples]
	X = [(r['T']+273.15)**-2 for r in G]
	eX = [f95 * 2 * r['SE_T'] * (r['T']+273.15)**-3 for r in G]
	Y = [r['D47'] for r in G]
	eY = [r['SE_D47']*f95 for r in G]	
	errorbar(X, Y, eY, eX, **(kw_errorbar_foram | dict(zorder = -200)))
	leg_eq_lakes, = plot(X, Y, label = 'Lacustrine carbonates (Anderson et al., 2021)', **(kw_lgbdvh | dict(marker = (4,0,0), ms = 6, zorder = -100)))

	c_anderson = dict(
		T = array([_['T'] for _ in anderson_data]),
		sT = array([_['SE_T'] for _ in anderson_data]),
		D47 = array([_['D47'] for _ in anderson_data]),
		sD47 = array([_['SE_D47'] for _ in anderson_data]),
		)

	c_anderson['X'] = [(t+273.15)**-2 for t in c_anderson['T']]
	c_anderson['sX'] = [2*st*(t+273.15)**-3 for t,st in zip(c_anderson['T'], c_anderson['sT'])]

	Aa, Ba, CMa, _ = YorkReg(c_anderson['X'], c_anderson['D47'], c_anderson['sX'], c_anderson['sD47'])

	yi = Aa * xi + Ba
# 	yi = 39100 * xi + 0.154
	leg_anderson, = plot(xi, yi, '-', color = ANDERSON_COLOR, alpha = 1, lw = 1.5, zorder = 5, label = 'MIT calibration (eq. 2)')

# 	plot(c_anderson['X'], c_anderson['D47'], 'r+')

	T = xi**-.5
	yi = 1.038*(-5.897/T - 3.521e3/T**2 + 2.391e7/T**3 - 3.541e9/T**4) + 0.1856
	leg_fiebig, = plot(xi, yi, '-', color = FIEBIG_COLOR, alpha = 1, lw = 1.5, zorder = 5, label = 'Fiebig et al. (2021)')

	leg_dvhlgb_title, = plot([], [], ls = 'None', marker = 'None', label = '\n$\mathbf{Devils\ Hole\ &\ Laghetto\ Basso:}\hspace{-3.5}$')
	leg_other_title, = plot([], [], ls = 'None', marker = 'None', label = '\n$\mathbf{Additional\ “cold”\ natural\ carbonates:}\hspace{-3.5}$')

	if figtitle.startswith('3_'):
		G = [
			r for r in data 
			if r['Type'] == 'planktic'
			and r['Twoa23_vs_Tiso_species'] == 'discordant'
			and not (TWOA_FOR_PLANKTICS_WITHOUT_OFFSET == 'nan' and r['Tiso_species_offset'] == 'all planktics')
			]
		X = [
			(r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-2
			if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
			(r[Tiso_field]+273.15)**-2
			for r in G]
		Y = [r['D47_ICDES'] for r in G]
		
		leg_discordant, = plot(X, Y, 'ko', mew = 0, ms = 1.5, zorder = 1000, label = 'Discordant planktics')
		legs = legs[:3] + [(leg_peralp, leg_meinicke, leg_discordant, leg_discordant)] + legs[3:]
		
	legend1 = legend(
		legs,
		[_.get_label() if 'get_label' in dir(_) else _[-1].get_label() for _ in legs],
		loc = 2,
		fontsize = 9,
		labelspacing = 0.3,
		frameon = False,
		ncols = 2,
		handlelength = 1.5,
		handler_map = {tuple: HandlerTuple(ndivide = 2)},
		)

	legend2 = legend(
		[leg_plreg],
		[leg_plreg.get_label()],
		loc = 2,
		fontsize = 9,
		labelspacing = 0.3,
		frameon = False,
		bbox_to_anchor = (0.01, .87),
		)

	for t in legend2.get_texts():
		t.set_y(t.get_position()[1]-1.6)
	
	legs = [
		leg_calibs,
		leg_anderson,
		leg_fiebig,
		leg_dvhlgb_title,
		leg_eq_lsce,
		leg_eq_gu,
		leg_other_title,
		leg_eq_pecten,
		leg_eq_lakes,
		]
	legend(
		legs,
		[_.get_label() for _ in legs],
		loc = 4,
		fontsize = 9,
		labelspacing = 0.4,
		frameon = False,
		markerfirst = False,
		)
		
	gca().add_artist(legend1)
	gca().add_artist(legend2)

	Ti = [30,20,10,0]
	xticks([(t+273.15)**-2 for t in Ti])
	gca().set_xticklabels([f"${t}\\,$°C" for t in Ti])
	ylabel(f'Δ$_{{47}}$ (I-CDES, ‰)')
	text(
		((10+273.15)**-2 + (20+273.15)**-2)/2,
		0.548,
		f'$1\\,/\\,T^2$',
		va = 'top',
		ha = 'center',
		size = 10,
		)

	axis([xmin, xmax, 0.558, 0.712])
	grid(False)

	savefig(f'plots/{figtitle}{"_WOA_for_cold_planktics" if TWOA_FOR_PLANKTICS_BELOW_ZERO else ""}{"_toggle_no_offset" if TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else ""}.pdf')
	close(fig)
	

	if figtitle in ['3_Regression_full']:

		'''BENTHIC RESIDUALS'''	
		# MARK Benthic residuals vs ref or species

		# Define planktic regression baseline
		planktic_data = zip(*(
			list(zip(*peral_planktic_data))
			+ list(zip(*meinicke_data))
			))
		A, B, CM, rchisq = YorkReg(*planktic_data)

		X0 = -CM[0,1]/CM[0,0]
		T0 = X0**-.5 - 273.15

		_str_ = '' if TWOA_FOR_PLANKTICS_BELOW_ZERO else '_with_cold_Tiso'
		with open('../output/saved_values.py', 'a') as fid:
			fid.write(f'''
planktic_reg_slope{_str_} = {A:.0f}
planktic_reg_slope_se{_str_} = {CM[0,0]**.5:.0f}
planktic_reg_intercept{_str_} = {B:.4f}
planktic_reg_intercept_se{_str_} = {CM[1,1]**.5:.4f}
planktic_reg_Tzero{_str_} = {T0:.2f}
planktic_reg_Tzero_intercept{_str_} = {B+A*X0:.4f}
planktic_reg_Tzero_intercept_se{_str_} = {(CM[1,1] - CM[0,1]**2 / CM[0,0])**.5:.4f}
''')

		fig = figure(figsize = (7,5))

		x1, x2, x3, x4 = .13, .5, .63, .88
		y1, y2, y3, y4 = .03, .57, .6, .97
		ax_site_residuals = fig.add_axes((x1, y3, x2-x1, y4-y3))
		ax_site_zscores = fig.add_axes((x3, y3, x4-x3, y4-y3))
		ax_species_residuals = fig.add_axes((x1, y1, x2-x1, y2-y1))
		ax_species_zscores = fig.add_axes((x3, y1, x4-x3, y2-y1))

		for plottype, axr, axz in [
			('by_site', ax_site_residuals, ax_site_zscores),
			('by_species', ax_species_residuals, ax_species_zscores),
			]:
			
			abd = []

			if plottype == 'by_species':
				G = [r for r in data if r['Type'] == 'benthic']
				species = sorted({r['Species'] for r in G})

				for r in G:

					if (
						r['Ref'] == 'Piasecki et al. (2019)'
						and piasecki_sites[r['Site']]['Ref'] not in ['O.A.', 'World Ocean Atlas 2009 volume 1: Temperature (2010)']
						):
						Tfield = 'Tpub'
					else:
						Tfield = 'Twoa23'

					x = (r[Tfield]+273.15)**-2
					y = r['D47_ICDES']
					sx = 2*r[f'SE_{Tfield}']*(r[Tfield]+273.15)**-3
					sy = r['SE_D47_ICDES']
					abd += [{
						'species': r['Species'],
						'genus': r['Species'].split(' ')[0],
						'Ref': r['Ref'],
						'r': ((y-B)/A)**-.5 - x**-.5,
						'sr': ((0.5 / A * sy * ((y-B)/A)**-1.5)**2 + (0.5 * sx * x**-1.5)**2)**.5,
# 						'r': (y - A * x - B),
# 						'sr': (sy**2 + A**2 * sx**2)**.5,
						'omega': bswchem[r['Site']]['Omega_cc'],
						's_omega': bswchem[r['Site']]['SE_Omega_cc'],
						'salinity': bswchem[r['Site']]['Salinity'],
						's_salinity': bswchem[r['Site']]['SE_Salinity'],
						'habitat': habitats[r['Species'].split(' ')[0]]
						}]

					colorlist = [
						(0,0,0),
						(0,.5,1),
						(1,0,0),
						(1,.5,.75),
						(1,1,1),
						(1,1,0),
						(1,.75,0),
						(1,.5,.75),
						(0,1,0),
						(0,1,1),
						]

# 				print(sorted([_['r'] for _ in abd if _['genus'] == 'Cibicidoides']))


			elif plottype == 'by_site':
				G = [r for r in data if r['Type'] == 'benthic']
				sites = sorted({r['Site'] for r in G})

				for s in sites:

					if (
						s in piasecki_sites
						and piasecki_sites[s]['Ref'] not in ['O.A.', 'World Ocean Atlas 2009 volume 1: Temperature (2010)']
						):
						Tfield = 'Tpub'
					else:
						Tfield = 'Twoa23'

					H = [r for r in G if r['Site'] == s]

					x = (H[0][Tfield]+273.15)**-2
					sx = 2*H[0][f'SE_{Tfield}']*(H[0][Tfield]+273.15)**-3
					y = sum([r['D47_ICDES']/r['SE_D47_ICDES']**2 for r in H]) / sum([1/r['SE_D47_ICDES']**2 for r in H])
					sy = sum([1/r['SE_D47_ICDES']**2 for r in H])**-.5
					
					abd += [{
						'site': s,
						'ref': H[0]['Ref'],
						'r': ((y-B)/A)**-.5 - x**-.5,
						'sr': ((0.5 / A * sy * ((y-B)/A)**-1.5)**2 + (0.5 * sx * x**-1.5)**2)**.5,
# 						'r': (y - A * x - B),
# 						'sr': (sy**2 + A**2 * sx**2)**.5,
						'omega': bswchem[s]['Omega_cc'],
						's_omega': bswchem[s]['SE_Omega_cc'],
						'salinity': bswchem[s]['Salinity'],
						's_salinity': bswchem[s]['SE_Salinity'],
						}]

			abd = sorted(abd, key = lambda r: r['r'])
			
			_str_ = '' if TWOA_FOR_PLANKTICS_BELOW_ZERO else '_with_cold_Tiso'
			if plottype == 'by_site':
				with open('../output/saved_values.py', 'a') as fid:
					fid.write(f'''
min_benthic_residuals{_str_} = {abd[0]['r']:.1f}
max_benthic_residuals{_str_} = {abd[-1]['r']:.1f}
''')
					
			kw = dict(
				ls = 'None',
				mec = 'k',
				mew = 1,
				ms = 4,
				)

			kweb = dict(
				ls = 'None',
				marker = 'None',
				ecolor = 'k',
				elinewidth = 1,
				capthick = 1,
				capsize = 1.5,
				)

			sca(axr)
			Y = array([_['r'] for _ in abd])
			eY = array([_['sr']*f95 for _ in abd])
			X = arange(len(Y))
			errorbar(X, Y, eY, **kweb)

			if plottype == 'by_species':
				for genus, color in zip(
					sorted({r['genus'] for r in abd}, key = lambda _: habitats[_]), colorlist):

					k = [r['genus'] == genus for r in abd]
					mrk = {
						'infaunal': 'd',
						'epifaunal': hdiamond,
						}[habitats[genus]]
					plot(X[k], Y[k], marker = mrk, mfc = color, label = genus, **kw)
			elif plottype == 'by_site':
				for ref, color in [
					('Peral et al. (2018)', (.5, .5, .5)),
					('Piasecki et al. (2019)', (0, 0, 0)),
					]:
					k = [r['ref'] == ref for r in abd]
					plot(X[k], Y[k], marker = 's', mfc = color, label = ref, **kw)

			axhline(0, color = 'k', alpha = 1/3, lw = 1)
			legend(
				ncols = 2 if plottype == 'by_species' else 1,
				labelspacing = 0.2,
				handlelength = 0.,
				borderaxespad = 0.,
				fancybox = False,
				framealpha = 1,
				edgecolor = 'k',
				borderpad = 0.75,
				prop = dict(style = 'italic', size = 8) if plottype == 'by_species' else dict(size = 8),
				)
	
			xticks([])
			ylabel(f'T residuals (°C) relative to\nplanktic foraminifer regression')
			axr.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.0f}$' if x else '$0$'))
			axr.yaxis.set_major_locator(ticker.MultipleLocator(3 if plottype == 'by_site' else 10))
			
			if plottype == 'by_species':
				axis([None, None, axis()[2]+1, axis()[3]+11])

			sca(axz)

			abd = sorted(abd, key = lambda r: r['r']/r['sr'])

			_str_ = '' if TWOA_FOR_PLANKTICS_BELOW_ZERO else '_with_cold_Tiso'
			if plottype == 'by_site':
				with open('../output/saved_values.py', 'a') as fid:
					fid.write(f'''
min_benthic_zscore{_str_} = {abd[0]['r']/abd[0]['sr']:.1f}
max_benthic_zscore{_str_} = {abd[-1]['r']/abd[-1]['sr']:.1f}
''')

			Y = array([_['r']/_['sr'] for _ in abd])
			eY = array([f95 for _ in abd])
			X = arange(len(Y))
	# 		errorbar(X, Y, eY, **kweb)

			if plottype == 'by_species':
				for genus, color in zip(
					sorted({r['genus'] for r in abd}, key = lambda _: habitats[_]), colorlist):

					k = [r['genus'] == genus for r in abd]
					mrk = {
						'infaunal': 'd',
						'epifaunal': hdiamond,
						}[habitats[genus]]
					plot(X[k], Y[k], marker = mrk, mfc = color, label = genus, **kw)
			elif plottype == 'by_site':
				for ref, color in [
					('Peral et al. (2018)', (.5, .5, .5)),
					('Piasecki et al. (2019)', (0, 0, 0)),
					]:
					k = [r['ref'] == ref for r in abd]
					plot(X[k], Y[k], marker = 's', mfc = color, label = ref, **kw)

# 			axhline(0, color = 'k', alpha = 1/3, lw = 1)
			axhspan(-f95, f95, color = 'k', alpha = 0.1, lw = 0)
	# 		legend(fontsize = 9, ncols = 2)
	
			xticks([])
			ylabel(f'Z-score residuals relative to\nplanktic foraminifer regression')
			axz.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.0f}$' if x else '$0$'))
			axz.yaxis.set_ticks_position('both')

			x = axis()[1]
			_ymin, _ymax = axis()[-2:]
			for tl in axz.yaxis.get_ticklabels():
				y = tl.get_unitless_position()[1]
				if y:
					p = norm.sf(abs(y))
					p = float(f'{p:.2g}')
				else:
					p = 1
				if _ymin < y < _ymax:
					text(x*1.05, y, str(p).replace('e-0','e-'), va = 'center', ha = 'left', size = 9)

			text(1.35, .5, 'one-sided p-value', rotation = -90, va = 'center', ha = 'center', transform = axz.transAxes)
		
			annotate('',
				xy = (max(X)*0.05, f95),
				xytext = (max(X)*0.05, -f95),
				arrowprops = dict(arrowstyle = '<->'),
				)
			text(max(X)*0.1, 0,
				'Expected 95 %\nconfidence\ninterval',
				va = 'center', ha = 'left', size = 8
				)


			if TWOA_FOR_PLANKTICS_BELOW_ZERO and TWOA_FOR_PLANKTICS_WITHOUT_OFFSET:
				_str_ = '' if TWOA_FOR_PLANKTICS_BELOW_ZERO else '_with_cold_Tiso'
				if plottype == 'by_site':
					with open('../output/saved_values.py', 'a') as fid:
						fid.write(f'''
zscore_outliers_bysite{_str_} = {len([_ for _ in abd if abs(_['r']/_['sr']) <= f95])}
zscore_total_bysite{_str_} = {len(abd)}
''')
				elif plottype == 'by_species':
					with open('../output/saved_values.py', 'a') as fid:
						fid.write(f'''
zscore_outliers_byspecies{_str_} = {len([_ for _ in abd if abs(_['r']/_['sr']) <= f95])}
zscore_total_byspecies{_str_} = {len(abd)}
''')
			
		savefig(f'plots/benthic_residuals{"_WOA_for_cold_planktics" if TWOA_FOR_PLANKTICS_BELOW_ZERO else ""}{"_toggle_no_offset" if TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else ""}.pdf')
		close(fig)



		'''BENTHIC RESIDUALS AS A FUNCTION OF SEAWATER CHEMISTRY'''
		# MARK Benthic residuals vs sewater chemistry

		abd = []

		G = [r for r in data if 'Peral' in r['Ref'] and r['Type'] == 'benthic']

		for r in G:
			Tfield = 'Twoa23'
			x = (r[Tfield]+273.15)**-2
			y = r['D47_ICDES']
			sx = 2*r[f'SE_{Tfield}']*(r[Tfield]+273.15)**-3
			sy = r['SE_D47_ICDES']
			abd += [{
				'Ref': "Peral et al. (2018)",
				'r': ((y-B)/A)**-.5 - x**-.5,
				'sr': ((0.5 / A * sy * ((y-B)/A)**-1.5)**2 + (0.5 * sx * x**-1.5)**2)**.5,
# 				'r': (y - A * x - B),
# 				'sr': (sy**2 + A**2 * sx**2)**.5,
				'omega': bswchem[r['Site']]['Omega_cc'],
				's_omega': bswchem[r['Site']]['SE_Omega_cc'],
				'salinity': bswchem[r['Site']]['Salinity'],
				's_salinity': bswchem[r['Site']]['SE_Salinity'],
				}]

		G = [r for r in data if 'Piasecki' in r['Ref']]
		for site in sorted({r['Site'] for r in G}):
			if (
				piasecki_sites[site]['Ref'] not in ['O.A.', 'World Ocean Atlas 2009 volume 1: Temperature (2010)']
				):
				Tfield = 'Tpub'
			else:
				Tfield = 'Twoa23'
			H = [r for r in G if r['Site'] == site]
			x = (H[0][Tfield]+273.15)**-2
			sx = 2*H[0][f'SE_{Tfield}']*(H[0][Tfield]+273.15)**-3
			y = sum([r['D47_ICDES']/r['SE_D47_ICDES']**2 for r in H]) / sum([1/r['SE_D47_ICDES']**2 for r in H])
			sy = sum([1/r['SE_D47_ICDES']**2 for r in H])**-.5
			abd += [{
				'Ref': "Piasecki et al. (2019)",
				'r': ((y-B)/A)**-.5 - x**-.5,
				'sr': ((0.5 / A * sy * ((y-B)/A)**-1.5)**2 + (0.5 * sx * x**-1.5)**2)**.5,
# 				'r': (y - A * x - B),
# 				'sr': (sy**2 + A**2 * sx**2)**.5,
				'omega': bswchem[site]['Omega_cc'],
				's_omega': bswchem[site]['SE_Omega_cc'],
				'salinity': bswchem[site]['Salinity'],
				's_salinity': bswchem[site]['SE_Salinity'],
				}]

		abd = sorted(abd, key = lambda x: x['r'])

		kw = dict(
			ls = 'None',
			marker = 's',
			mec = 'k',
			mew = 1,
			ms = 4,
			)

		kweb = dict(
			ls = 'None',
			marker = 'None',
			ecolor = [.75]*3,
			elinewidth = 1,
			capthick = 1,
			capsize = 1.5,
			)

		fig = figure(figsize = (5,2.6))
		subplots_adjust(.16,.18,.975,.95,.05)

# 		ax1 = subplot(1,3,1)	
# 		Y = array([_['r'] for _ in abd])*1000
# 		eY = array([_['sr']*f95 for _ in abd])*1000
# 		isperal = array([True if 'Peral' in _['Ref'] else False for _ in abd])
# 		X = arange(len(Y))
# 		errorbar(X, Y, eY, **kweb)
# 		plot(X[isperal], Y[isperal], mfc = [.5]*3, label = 'Peral et al. (2018)', **kw)
# 		plot(X[~isperal], Y[~isperal], mfc = 'k', label = 'Piasecki et al. (2019)', **kw)
# 		axhline(0, color = 'k', alpha = 1/3, lw = 1)
# 	
# 		legend(fontsize = 7.5, loc = 'lower right')
# 		xticks([])
# 		ylabel(f'Δ$_{{47}}$ residuals (ppm) relative to\nplanktic foraminifer regression\n')
# 		ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.0f}$' if x else '$0$'))

# 		ax2 = subplot(1,3,2)	
		ax2 = subplot(1,2,1)
		Y = array([_['r'] for _ in abd])
		eY = array([_['sr']*f95 for _ in abd])
		isperal = array([True if 'Peral' in _['Ref'] else False for _ in abd])
		X = array([_['omega'] for _ in abd])
		eX = array([_['s_omega']*f95 for _ in abd])
		errorbar(X, Y, eY, eX, **kweb)
		plot(X[isperal], Y[isperal], mfc = [.5]*3, label = 'Peral et al. (2018)', **kw)
		plot(X[~isperal], Y[~isperal], mfc = 'k', label = 'Piasecki et al. (2019)', **kw)
		ylabel(f'T residuals (°C) relative to\nplanktic foraminifer regression')
	
		xlabel(f'Calcite saturation')
		ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.0f}$' if x else '$0$'))
		ax2.xaxis.set_major_locator(ticker.MultipleLocator(1))
		ax2.yaxis.set_major_locator(ticker.MultipleLocator(3))
		ax2.axis([None, None, None, 4])

# 		ax3 = subplot(1,3,3)
		ax3 = subplot(1,2,2)
		Y = array([_['r'] for _ in abd])
		eY = array([_['sr']*f95 for _ in abd])
		isperal = array([True if 'Peral' in _['Ref'] else False for _ in abd])
		X = array([_['salinity'] for _ in abd])
		eX = array([_['s_salinity']*f95 for _ in abd])
		errorbar(X, Y, eY, eX, **kweb)
		plot(X[isperal], Y[isperal], mfc = [.5]*3, label = 'Peral et al. (2018)', **kw)
		plot(X[~isperal], Y[~isperal], mfc = 'k', label = 'Piasecki et al. (2019)', **kw)
	
		xlabel(f'Salinity')
		ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ''))
		ax3.yaxis.set_ticks_position('none')
		ax3.xaxis.set_major_locator(ticker.MultipleLocator(1))
		ax3.yaxis.set_major_locator(ticker.MultipleLocator(3))
		ax3.axis([None, None, None, 4])

# 		ylabel(f'Δ$_{{47}}$ residuals (ppm) relative to\nplanktic foraminifer regression\n')

		savefig(f'plots/benthic_residuals_vs_seawater_chemistry{"_WOA_for_cold_planktics" if TWOA_FOR_PLANKTICS_BELOW_ZERO else ""}{"_toggle_no_offset" if TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else ""}.pdf')
		close(fig)



	if figtitle == '3_Regression_full':

		'''PLANKTIC RESIDUALS'''
		# MARK Planktic residuals vs ref or Tiso_offset

		for calib in ['planktic', 'fiebig', 'anderson']:
			for whichplot in ['by_ref', 'by_Tiso_offset']:

				fig = figure(figsize = (3.5,4.5))
				ax1, ax2 = subplot(211),  subplot(212)
				subplots_adjust(.22,.03,.95,.98,.35, .1)

				abd = []

				G = [r for r in data if r['Type'] == 'planktic']

				for r in G:
					if TWOA_FOR_PLANKTICS_BELOW_ZERO and r['Tiso_species'] < 0:
						Tfield = TWOA_FOR_PLANKTICS_BELOW_ZERO
					else:
						Tfield = 'Tiso_species'

					T = r[Tfield] + 273.15
					sT = r[f'SE_{Tfield}']

					x = T**-2
					y = r['D47_ICDES']
					sx = 2*sT*T**-3
					sy = r['SE_D47_ICDES']
					abd += [{
						'sample': r['Sample'],
						'species': 'Trilobatus sacculifer' if r['Species'] == 'Trilobatus trilobus' else r['Species'],
						'offset': r['Tiso_species_offset'],
						'genus': r['Species'].split(' ')[0],
						'ref': r['Ref'],
						'r': {
							'planktic': y - A * x - B,
							'fiebig': y - (1.038*(-5.897/T - 3.521e3/T**2 + 2.391e7/T**3 - 3.541e9/T**4) + 0.1856),
							'anderson': y - Aa * x - Ba,
							}[calib],
						'sr': {
							'planktic': (sy**2 + A**2 * sx**2)**.5,
							'fiebig': (sy**2 + (1.038*sT*(5.897/T**2 + 3.521e3*2/T**3 - 2.391e7*3/T**4 + 3.541e9*4/T**5))**2)**.5,
							'anderson': (sy**2 + Aa**2 * sx**2)**.5,
							}[calib],
						'omega': bswchem[r['Site']]['Omega_cc'],
						's_omega': bswchem[r['Site']]['SE_Omega_cc'],
						'salinity': bswchem[r['Site']]['Salinity'],
						's_salinity': bswchem[r['Site']]['SE_Salinity'],
						}]

				abd = sorted(abd, key = lambda r: r['r'])
				
				kw = dict(
					ls = 'None',
					marker = 'd',
					mec = 'k',
					mew = 0.7,
					ms = 4,
					)

				kweb = dict(
					ls = 'None',
					marker = 'None',
					ecolor = 'k',
					elinewidth = 0.7,
					capthick = 0.7,
					capsize = 1,
					)

				sca(ax1)
				Y = array([_['r'] for _ in abd])*1000
				eY = array([_['sr']*f95 for _ in abd])*1000
				X = arange(len(Y))
				errorbar(X, Y, eY, **kweb)

				if whichplot == 'by_ref':
					for ref, color in [
						('Peral et al. (2018)', (.5, .5, .5)),
						('Meinicke et al. (2020)', (1, 1, 1)),
						]:
						k = [r['ref'] == ref for r in abd]
						plot(X[k], Y[k], mfc = color, label = ref, **kw)
				elif whichplot == 'by_Tiso_offset':

						label = '$^{18}$α based on species'
						color = (1,1,1)
						k = [r['species'] == r['offset'] for r in abd]
						plot(X[k], Y[k], mfc = color, label = label, **kw)

						label = '$^{18}$α based on genus'
						color = (.5,.5,.5)
						k = [r['species'] != r['offset'] and r['offset'].endswith('spp.') for r in abd]
						plot(X[k], Y[k], mfc = color, label = label, **kw)

						label = 'P. obliquiloculata'
						color = (0,0,0)
						k = [r['offset'] == 'all planktics' for r in abd]
						plot(X[k], Y[k], mfc = color, label = label, **kw)

				axhline(0, color = 'k', alpha = 1/3, lw = 1)
				legend(fontsize = 8, labelspacing = 0.2)

				xticks([])
				if calib == 'planktic':
					ylabel(f'Δ$_{{47}}$ residuals (ppm) relative to\nplanktic foraminifer regression')
				elif calib == 'fiebig':
					ylabel(f'Δ$_{{47}}$ residuals (ppm) relative to\nFiebig et al. (2021) calibration')
				elif calib == 'anderson':
					ylabel(f'Δ$_{{47}}$ residuals (ppm) relative to\nMIT calibration')
				gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.0f}$' if x else '$0$'))
		# 		axis([None, None, None, 200])

				sca(ax2)

				abd = sorted(abd, key = lambda r: r['r']/r['sr'])

				Y = array([_['r']/_['sr'] for _ in abd])
				eY = array([f95 for _ in abd])
				X = arange(len(Y))
		# 		errorbar(X, Y, eY, **kweb)

				if whichplot == 'by_ref':

					for ref, color in [
						('Peral et al. (2018)', (.5, .5, .5)),
						('Meinicke et al. (2020)', (1, 1, 1)),
						]:
						k = [r['ref'] == ref for r in abd]
						plot(X[k], Y[k], mfc = color, label = ref, **kw)

				elif whichplot == 'by_Tiso_offset':

					label = '$^{18}$α based on species'
					color = (1,1,1)
					k = [r['species'] == r['offset'] for r in abd]
					plot(X[k], Y[k], mfc = color, label = label, **kw)

					label = '$^{18}$α based on genus'
					color = (.5,.5,.5)
					k = [r['species'] != r['offset'] and r['offset'].endswith('spp.') for r in abd]
					plot(X[k], Y[k], mfc = color, label = label, **kw)

					label = 'P. obliquiloculata'
					color = (0,0,0)
					k = [r['offset'] == 'all planktics' for r in abd]
					plot(X[k], Y[k], mfc = color, label = label, **kw)


				annotate('',
					xy = (max(X)*1.01, f95),
					xytext = (max(X)*1.01, -f95),
					arrowprops = dict(arrowstyle = '<->'),
					)
				text(max(X)*0.95, -1,
					'Expected 95 %\nconfidence\ninterval',
					va = 'center', ha = 'right', size = 8
					)
				
# 				ax2.axis([None, X.max()*1.1, None, None])

				if calib == 'fiebig' and whichplot == 'by_ref':
					with open('../output/saved_values.py', 'a') as fid:
						fid.write(f'''
N_planktic_samples = {len(Y)}
N_planktic_samples_with_fiebig_absZ_greater_than_2 = {len(Y[abs(Y)>2])}
''')

				axhline(0, color = 'k', alpha = 1/3, lw = 1)
				axhspan(-f95, f95, color = 'k', alpha = 1/10, lw = 1)
		# 		legend(fontsize = 9, ncols = 2)

				xticks([])
				if calib == 'planktic':
					ylabel(f'Z-score residuals relative to\nplanktic foraminifer regression')
				elif calib == 'fiebig':
					ylabel(f'Z-score residuals relative to\nFiebig et al. (2021) calibration')
				elif calib == 'anderson':
					ylabel(f'Z-score residuals relative to\nMIT calibration')
				ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.0f}$' if x else '$0$'))

	# 			ax2.yaxis.set_ticks_position('both')

	# 			x = axis()[1]
	# 			for tl in ax2.yaxis.get_ticklabels():
	# 				y = tl.get_unitless_position()[1]
	# 				if y:
	# 					p = norm.sf(abs(y))
	# 					p = float(f'{p:.2g}')
	# 				else:
	# 					p = 1
	# 				text(x*1.05, y, str(p), va = 'center', ha = 'left', size = 9)
	# 
	# 			text(1.2, 0.5, 'one-sided p-value', rotation = -90, va = 'center', ha = 'left', transform = ax2.transAxes)
	
				savefig(f'plots/Planktic_residuals{"" if calib == "planktic" else "_vs_"+calib}{"_WOA_for_cold_planktics" if TWOA_FOR_PLANKTICS_BELOW_ZERO else ""}{"_toggle_no_offset" if TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else ""}_{whichplot}.pdf')
				close(fig)

	if figtitle == '3_Regression_full':


		'''REGRESSION_ELLIPSE_PLOT'''
		# MARK Regression ellipses

		T0 = 15
		DRAW_PERAL = False
		DRAW_MEINICKE = False
		DRAW_ALL_PLANKTICS = True
		DRAW_CONCORDANT_PLANKTICS = True
		DRAW_DISCORDANT_PLANKTICS = True
		DRAW_FIEBIG = True
		DRAW_ANDERSON = True
		DRAW_DVHLGB = True

		G = [r for r in data if r['Type'] == 'planktic']
		_peral_planktic_data = [r for r in G if 'Peral' in r['Ref']]
		_meinicke_planktic_data = [r for r in G if 'Meinicke' in r['Ref']]
		_concordant_planktic_data = [r for r in G if r['Twoa23_vs_Tiso_species'] == '']
		_discordant_planktic_data = [r for r in G if r['Twoa23_vs_Tiso_species'] == 'discordant']

		peral_planktic_XsXYsY = [
			(
				(
					(r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-2 - (T0+273.15)**-2
					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					(
						(r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-2 - (T0+273.15)**-2
						if r['Tiso_species_offset'] == 'all planktics' and TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else
						(r['Tiso_species']+273.15)**-2 - (T0+273.15)**-2
						)
					),
				r['D47_ICDES'],
				(
					2*r['SE_'+TWOA_FOR_PLANKTICS_BELOW_ZERO]*(r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-3
					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					(
						2*r['SE_'+TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]*(r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-3
						if r['Tiso_species_offset'] == 'all planktics' and TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else
						2*r['SE_Tiso_species']*(r['Tiso_species']+273.15)**-3
						)
					),
				r['SE_D47_ICDES'],
				)
			for r in _peral_planktic_data
			]

		meinicke_planktic_XsXYsY = [
			(
				(
					(r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-2 - (T0+273.15)**-2
					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					(
						(r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-2 - (T0+273.15)**-2
						if r['Tiso_species_offset'] == 'all planktics' and TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else
						(r['Tiso_species']+273.15)**-2 - (T0+273.15)**-2
						)
					),
				r['D47_ICDES'],
				(
					2*r['SE_'+TWOA_FOR_PLANKTICS_BELOW_ZERO]*(r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-3
					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					(
						2*r['SE_'+TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]*(r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-3
						if r['Tiso_species_offset'] == 'all planktics' and TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else
						2*r['SE_Tiso_species']*(r['Tiso_species']+273.15)**-3
						)
					),
				r['SE_D47_ICDES'],
				)
			for r in _meinicke_planktic_data
			]

		all_planktic_XsXYsY = [
			(
				(
					(r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-2 - (T0+273.15)**-2
					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					(
						(r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-2 - (T0+273.15)**-2
						if r['Tiso_species_offset'] == 'all planktics' and TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else
						(r['Tiso_species']+273.15)**-2 - (T0+273.15)**-2
						)
					),
				r['D47_ICDES'],
				(
					2*r['SE_'+TWOA_FOR_PLANKTICS_BELOW_ZERO]*(r[TWOA_FOR_PLANKTICS_BELOW_ZERO]+273.15)**-3
					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					(
						2*r['SE_'+TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]*(r[TWOA_FOR_PLANKTICS_WITHOUT_OFFSET]+273.15)**-3
						if r['Tiso_species_offset'] == 'all planktics' and TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else
						2*r['SE_Tiso_species']*(r['Tiso_species']+273.15)**-3
						)
					),
				r['SE_D47_ICDES'],
				)
			for r in G
			]

		concordant_planktic_XsXYsY = [
			(
				(
# 					(r['Twoa23']+273.15)**-2 - (T0+273.15)**-2
# 					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					(r['Tiso_species']+273.15)**-2 - (T0+273.15)**-2
					),
				r['D47_ICDES'],
				(
# 					2*r['SE_Twoa23']*(r['Twoa23']+273.15)**-3
# 					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					2*r['SE_Tiso_species']*(r['Tiso_species']+273.15)**-3
					),
				r['SE_D47_ICDES'],
				)
			for r in _concordant_planktic_data
			]

		discordant_planktic_XsXYsY = [
			(
				(
# 					(r['Twoa23']+273.15)**-2 - (T0+273.15)**-2
# 					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					(r['Tiso_species']+273.15)**-2 - (T0+273.15)**-2
					),
				r['D47_ICDES'],
				(
# 					2*r['SE_Twoa23']*(r['Twoa23']+273.15)**-3
# 					if r['Tiso_species'] < 0 and TWOA_FOR_PLANKTICS_BELOW_ZERO else
					2*r['SE_Tiso_species']*(r['Tiso_species']+273.15)**-3
					),
				r['SE_D47_ICDES'],
				)
			for r in _discordant_planktic_data
			]

		Aperal, Bperal, CMperal, _ = YorkReg(*zip(*peral_planktic_XsXYsY))
		Ameinicke, Bmeinicke, CMmeinicke, _ = YorkReg(*zip(*meinicke_planktic_XsXYsY))
		Ap, Bp, CMp, _ = YorkReg(*zip(*all_planktic_XsXYsY))
		Ap, Bp, CMp, _ = YorkReg(*zip(*all_planktic_XsXYsY))
		Ac, Bc, CMc, _ = YorkReg(*zip(*concordant_planktic_XsXYsY))
		Ad, Bd, CMd, _ = YorkReg(*zip(*discordant_planktic_XsXYsY))

		fig = figure(figsize = (3.6,4))
		fig.subplots_adjust(0.19, 0.14, .95, .95)
		ax = subplot(111)

		F = array([[.001, 1]]) # change scale of slope
		kw = dict(ls = '-', marker = 'None', lw = 1, color = 'k', dashes = (None, None))

		kw['dashes'] = (6,2,2,2)
		if DRAW_CONCORDANT_PLANKTICS:
			A, B, CM = Ac, Bc, CMc
			w,h,r = cov_ellipse(F.T * CM * F, 0.95)
			plot([], [], label = 'Concordant planktics', **kw)
			ax.add_artist(
				Ellipse(
					xy = (A / 1000, B), width = w, height = h, angle = r,
					fc = 'None', ec = kw['color'], lw = kw['lw'], ls = (0, kw['dashes']), zorder = 100
					)
				)


		kw['color'] = [0.5]*3
		if DRAW_DISCORDANT_PLANKTICS:
			A, B, CM = Ad, Bd, CMd
			w,h,r = cov_ellipse(F.T * CM * F, 0.95)
			plot([], [], label = 'Discordant planktics', **kw)
			ax.add_artist(
				Ellipse(
					xy = (A / 1000, B), width = w, height = h, angle = r,
					fc = 'None', ec = kw['color'], lw = kw['lw'], ls = (0, kw['dashes']), zorder = 100
					)
				)

		kw['color'] = 'k'
		kw['dashes'] = (None, None)
		if DRAW_ALL_PLANKTICS:
			A, B, CM = Ap, Bp, CMp
			w,h,r = cov_ellipse(F.T * CM * F, 0.95)
			plot([], [], label = 'All planktics (eqs. 5-6)', **kw)
			ax.add_artist(
				Ellipse(
					xy = (A / 1000, B), width = w, height = h, angle = r,
					fc = 'None', ec = kw['color'], lw = kw['lw'], zorder = 100
					)
				)


		kw['dashes'] = (1,1)
		kw['color'] = 'k'
		if DRAW_PERAL:
			A, B, CM = Aperal, Bperal, CMperal
			w,h,r = cov_ellipse(F.T * CM * F, 0.95)
			plot([], [], label = 'Reprocessed Peral et al.', **kw)
			ax.add_artist(
				Ellipse(
					xy = (A / 1000, B), width = w, height = h, angle = r,
					fc = 'None', ec = kw['color'], lw = kw['lw'], ls = (0, kw['dashes']), zorder = 100
					)
				)

		kw['color'] = [.5]*3
		if DRAW_MEINICKE:
			A, B, CM = Ameinicke, Bmeinicke, CMmeinicke
			w,h,r = cov_ellipse(F.T * CM * F, 0.95)
			plot([], [], label = 'Reprocesssed Meinicke et al.', **kw)
			ax.add_artist(
				Ellipse(
					xy = (A / 1000, B), width = w, height = h, angle = r,
					fc = 'None', ec = kw['color'], lw = kw['lw'], ls = (0, kw['dashes']), zorder = 100
					)
				)

		kw['color'] = [0.5]*3

		if DRAW_FIEBIG:
			Bf = 1.038 * (
				- 5.897/(T0+273.15)
				- 3.521e3/(T0+273.15)**2
				+ 2.391e7/(T0+273.15)**3
				- 3.541e9/(T0+273.15)**4
				) + 0.1856
			X0 = (T0+273.15)**-2
			Af = 1.038 * (
				- 0.5* 5.897 * X0**-.5
				- 3.521e3
				+ 1.5* 2.391e7 * X0**.5
				- 2* 3.541e9 * X0
				) / 1000

			ax.errorbar(Af, Bf, 0.005, elinewidth = 1.2, ecolor = FIEBIG_COLOR, capthick = 1.2, capsize = 3, ls = 'None', marker = 'None', zorder = 1)
			ax.plot(Af, Bf, 'wo', ms = 6, mec = FIEBIG_COLOR, mew = 1.2, label = 'Fiebig et al. (2021)')

		if DRAW_ANDERSON:
			with open('anderson_data.csv') as fid: # Anderson et al. (2021)
				anderson_data = list(DictReader(fid))

			for _ in anderson_data:
				_['D47'] = float(_['D47 I-CDES90'])
				_['SE_D47'] = float(_['SE'])
				_['T'] = float(_['Temperature'])
				_['SE_T'] = float(_['Temperature error'])
		
			c_anderson = dict(
				T = array([_['T'] for _ in anderson_data]),
				sT = array([_['SE_T'] for _ in anderson_data]),
				D47 = array([_['D47'] for _ in anderson_data]),
				sD47 = array([_['SE_D47'] for _ in anderson_data]),
				)

			c_anderson['X'] = [(t+273.15)**-2 - (T0+273.15)**-2 for t in c_anderson['T']]
			c_anderson['sX'] = [2*st*(t+273.15)**-3 for t,st in zip(c_anderson['T'], c_anderson['sT'])]

			Aa, Ba, CMa, _ = YorkReg(c_anderson['X'], c_anderson['D47'], c_anderson['sX'], c_anderson['sD47'])


			w,h,r = cov_ellipse(F.T * CMa * F, 0.95)
			kw = dict(ls = '-', marker = 'None', lw = 1, color = ANDERSON_COLOR)
			plot([], [], label = 'MIT calibration (eq. 2)', **kw)
			ax.add_artist(
				Ellipse(
					xy = (Aa / 1000, Ba), width = w, height = h, angle = r,
					fc = 'None', ec = kw['color'], lw = kw['lw']*1.2,
					)
				)

		if DRAW_DVHLGB:
			T = array([[d['T'] for d in eqdata if d['Sample'] == s][0] for s in ['DVH-2', 'LGB-2']])
			sT = array([[d['SE_T'] for d in eqdata if d['Sample'] == s][0] for s in ['DVH-2', 'LGB-2']])
			X = (T+273.15)**-2 - (T0+273.15)**-2
			sX = 2*sT*(T+273.15)**-3


			D47 = [[d['D47_ICDES'] for d in eqdata if d['Sample'] in s] for s in [['DVH-2', 'DHC2-8'], ['LGB-2']]]
			sD47 = [[d['SE_D47_ICDES'] for d in eqdata if d['Sample'] in s] for s in [['DVH-2', 'DHC2-8'], ['LGB-2']]]

			for k in range(2):
				w = array(sD47[k])**-2
				w /= w.sum()
				D47[k] = (D47[k]*w).sum()
				sD47[k] = ((sD47[k]*w)**2).sum()**.5

			Ae, Be, CMe, _ = YorkReg(X, D47, sX, sD47)
			
			_color_ = (.8,0,1)

			w,h,r = cov_ellipse(F.T * CMe * F, 0.95)
			kw = dict(ls = '-', marker = 'None', lw = 1, color = _color_)
			plot([], [], label = 'Devils Laghetto calibration (eq. 4)', **kw)
			ax.add_artist(
				Ellipse(
					xy = (Ae / 1000, Be), width = w, height = h, angle = r,
					fc = 'None', ec = kw['color'], lw = kw['lw']*1.2,
					)
				)





		ax.legend(fontsize = 8, labelspacing = 0.3, loc = 'lower left', bbox_to_anchor = (0, 0), frameon = False)
		ax.axis([29, 48, .6085, .6325])

		ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
		ax.yaxis.set_major_locator(ticker.MultipleLocator(0.01))

		ax.set_xlabel(f'Regression slope / $10^3$ at {T0:.0f} °C')
		ax.set_ylabel(f'Δ$_{{47}}$ (‰, I-CDES) at {T0:.0f} °C')
		fig.savefig(f'plots/planktic_regression_ellipses{"_WOA_for_cold_planktics" if TWOA_FOR_PLANKTICS_BELOW_ZERO else ""}{"_toggle_no_offset" if TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else ""}.pdf')
		close(fig)



		'''COMPARE_REGRESSIONS'''
		# MARK Compare regressions

		fig = figure(figsize = (3.6,4))
		subplots_adjust(.23, .15, .95, .95)

		planktic_data = list(zip(*(
			list(zip(*peral_planktic_data))
			+ list(zip(*meinicke_data))
			)))

		A, B, CM, rchisq = YorkReg(*planktic_data)

		Ti = linspace(-2,37)
		xi = (Ti+273.15)**-2
		D47i_planktic = A*xi+B
		sD47i_planktic = (CM[0,0]*xi**2 + 2*CM[0,1]*xi + CM[1,1])**.5
		sTi_planktic = 0.5 * sD47i_planktic / A * ((D47i_planktic-B)/A)**-1.5

		# 		obsx, obsy = array(planktic_data[:2])
		# 		plot(obsx**-.5 - 273.15, obsx**-.5 - ((obsy-B)/A)**-.5, 'ko', mew = 0, alpha = .2, zorder = 2)

		with open('anderson_data.csv') as fid: # Anderson et al. (2021)
			anderson_data = list(DictReader(fid))

		for _ in anderson_data:
			_['D47'] = float(_['D47 I-CDES90'])
			_['SE_D47'] = float(_['SE'])
			_['T'] = float(_['Temperature'])
			_['SE_T'] = float(_['Temperature error'])
	
		c_anderson = dict(
			T = array([_['T'] for _ in anderson_data]),
			sT = array([_['SE_T'] for _ in anderson_data]),
			D47 = array([_['D47'] for _ in anderson_data]),
			sD47 = array([_['SE_D47'] for _ in anderson_data]),
			)

		c_anderson['X'] = [(t+273.15)**-2 for t in c_anderson['T']]
		c_anderson['sX'] = [2*st*(t+273.15)**-3 for t,st in zip(c_anderson['T'], c_anderson['sT'])]

		Aa, Ba, CMa, _ = YorkReg(c_anderson['X'], c_anderson['D47'], c_anderson['sX'], c_anderson['sD47'])
		_X0 = -CMa[0,1]/CMa[0,0]
		_T0 = _X0**-.5 - 273.15
		

		if TWOA_FOR_PLANKTICS_BELOW_ZERO:

			with open('../output/saved_values.py', 'a') as fid:
				fid.write(f'''
mit_reg_slope{_str_} = {Aa:.0f}
mit_reg_slope_se{_str_} = {CMa[0,0]**.5:.0f}
mit_reg_intercept{_str_} = {Ba:.4f}
mit_reg_intercept_se{_str_} = {CMa[1,1]**.5:.4f}
mit_reg_Tzero{_str_} = {_T0:.2f}
mit_reg_Tzero_intercept{_str_} = {Ba+Aa*_X0:.4f}
mit_reg_Tzero_intercept_se{_str_} = {(CMa[1,1] - CMa[0,1]**2 / CMa[0,0])**.5:.4f}
''')

		# T = ((D47i_planktic - 0.154)/39100)**-.5 - 273.15
		# legs_anderson, = plot(Ti, T-Ti, '-', color = ANDERSON_COLOR, alpha = 1, lw = 1.5, zorder = 2, label = 'Anderson et al. (2021)')

		T = ((D47i_planktic - Ba)/Aa)**-.5 - 273.15
		sD47 = (CMa[0,0]*xi**2 + 2*CMa[0,1]*xi + CMa[1,1])**.5
		sT = 0.5 * sD47 / Aa * ((D47i_planktic-Ba)/Aa)**-1.5
		print(f'Avg 95% CL T for Anderson et al. regression is ± {sT.mean()*f95:.1f} C')

		T_mit = T*1

		legs_anderson, = plot(Ti, T-Ti, '-', color = ANDERSON_COLOR, alpha = 1, lw = 1.5, zorder = 20, label = 'MIT calibration (eq. 2)')
# 		legs_anderson_cl, = plot(Ti, T-Ti-f95*sT, ':', color = ANDERSON_COLOR, alpha = 1, lw = .8, zorder = 2, dashes = (6,2), label = '95 % CL for Anderson et al. (2021)')
# 		plot(Ti, T-Ti+f95*sT, ':', color = ANDERSON_COLOR, alpha = 1, lw = .8, zorder = 2, dashes = (6,2))

		# Meinicke 2020
		T = ((D47i_planktic - 0.1518)/39700)**-.5 - 273.15
		legs_meinicke, = plot(Ti, T-Ti, '-', color = 'k', alpha = .5, lw = 1.2, zorder = 2, dashes = (6,2), label = 'Meinicke et al. (2021)')

		# Peral 2018
		T = ((D47i_planktic - 0.1815)/36995)**-.5 - 273.15
		legs_peral, = plot(Ti, T-Ti, '-', color = 'k', alpha = .5, lw = 1.2, zorder = 2, dashes = (8,2,2,2), label = 'Peral et al. (2022)')

		# Fiebig
		_T = linspace(-10,40, 4001) + 273.15
		_D47 = 1.038*(-5.897/_T - 3.521e3/_T**2 + 2.391e7/_T**3 - 3.541e9/_T**4) + 0.1856
		_f = interp1d(_D47, _T)
		T = _f(D47i_planktic) - 273.15
		sD47i_planktic = (0.006 - T/30*0.001)/f95
		sT = -(_f(D47i_planktic+0.1*sD47i_planktic) - _f(D47i_planktic-0.1*sD47i_planktic))/0.2
		print(f'Avg 95% CL T for Fiebig et al. regression is ± {sT.mean()*f95:.1f} C')

		T_fiebig = T*1

		legs_fiebig, = plot(Ti, T-Ti, '-', color = FIEBIG_COLOR, alpha = 1, lw = 1.5, zorder = 2, label = 'Fiebig et al. (2021)')
# 		legs_fiebig_cl, = plot(Ti, T-Ti-f95*sT, ':', color = FIEBIG_COLOR, alpha = 1, lw = .8, zorder = 2, dashes = (6,2), label = '95 % CL for Fiebig et al. (2021)')
# 		plot(Ti, T-Ti+f95*sT, ':', color = FIEBIG_COLOR, alpha = 1, lw = .8, zorder = 2, dashes = (6,2))

		# 		# Concordant planktic data
		# 		Acp, Bcp = YorkReg(*concordant_planktic_data)[:2]
		# 		T = ((D47i_planktic - Bcp)/Acp)**-.5 - 273.15
		# 		plot(Ti, T-Ti, '-', color = (.7, .7, .7), alpha = 1, lw = 1.5, zorder = 2, label = 'Planktic regression excluding discordant samples')

		DRAW_DVHLGB = True
		if DRAW_DVHLGB:
			T = array([[d['T'] for d in eqdata if d['Sample'] == s][0] for s in ['DVH-2', 'LGB-2']])
			sT = array([[d['SE_T'] for d in eqdata if d['Sample'] == s][0] for s in ['DVH-2', 'LGB-2']])
			X = (T+273.15)**-2
			sX = 2*sT*(T+273.15)**-3

			D47 = [[d['D47_ICDES'] for d in eqdata if d['Sample'] in s] for s in [['DVH-2', 'DHC2-8'], ['LGB-2']]]
			sD47 = [[d['SE_D47_ICDES'] for d in eqdata if d['Sample'] in s] for s in [['DVH-2', 'DHC2-8'], ['LGB-2']]]

			for k in range(2):
				w = array(sD47[k])**-2
				w /= w.sum()
				D47[k] = (D47[k]*w).sum()
				sD47[k] = ((sD47[k]*w)**2).sum()**.5

			Ae, Be, CMe, _ = YorkReg(X, D47, sX, sD47)

			dvhlgbreg_rmse = (
				(array([
					d['D47_ICDES'] - Ae / (d['T']+273.15)**2 - Be for d in eqdata
					])**2).sum()
				/ (len(eqdata) - 2)
				)**.5

			_X0 = -CMe[0,1]/CMe[0,0]
			_T0 = _X0**-.5 - 273.15

			if TWOA_FOR_PLANKTICS_BELOW_ZERO:
				with open('../output/saved_values.py', 'a') as fid:
					fid.write(f'''
dvhlgbreg_rmse = {dvhlgbreg_rmse:.4f}

dvhlgb_reg_slope{_str_} = {Ae:.0f}
dvhlgb_reg_slope_se{_str_} = {CMe[0,0]**.5:.0f}
dvhlgb_reg_intercept{_str_} = {Be:.4f}
dvhlgb_reg_intercept_se{_str_} = {CMe[1,1]**.5:.4f}
dvhlgb_reg_Tzero{_str_} = {_T0:.2f}
dvhlgb_reg_Tzero_intercept{_str_} = {Be+Ae*_X0:.4f}
dvhlgb_reg_Tzero_intercept_se{_str_} = {(CMe[1,1] - CMe[0,1]**2 / CMe[0,0])**.5:.4f}
''')

			T = ((D47i_planktic - Be)/Ae)**-.5 - 273.15
# 			_color_ = [_*.8 for _ in DVHLGB_COLOR]
			_color_ = (.8,0,1)
			legs_dvhlgb, = plot(Ti, T-Ti, '-', color = _color_, alpha = 1, lw = 1.5, zorder = 10, label = 'Devils Laghetto calibration (eq. 4)')

			sD47 = (CMe[0,0]*xi**2 + 2*CMe[0,1]*xi + CMe[1,1])**.5
			sT = 0.5 * sD47 / Ae * ((D47i_planktic-Be)/Ae)**-1.5
			print(f'Avg 95% CL T for Devils Laghetto regression is ± {sT.mean()*f95:.1f} C')

			T_dvhlgb = T*1

			print('\nThree I-CDES calibration comparison:')
			for a,b,c,d in zip(Ti[::5], T_mit[::5], T_fiebig[::5], T_dvhlgb[::5]):
				print(f'{a:.1f} C: calibration spread of ± {(max([b,c,d]) - min([b,c,d]))/2:.2f} C')
			print()

			_f = interp1d(Ti, T-Ti)
			
			legs_dvhlgb_dot, = plot(X**-.5-273.15, _f(X**-.5-273.15), 'o', mew = 1.25, mec = _color_, mfc = 'w', ms = 5, zorder = 11)

		legs_ts = fill_between(
			Ti, f95*sTi_planktic, -f95*sTi_planktic,
			facecolor = [0,0,0,0.08], edgecolor = [.7]*3, label = '95 % CL for planktic foram. (eqs. 5-6)', lw = 1)

# 		plot(Ti, f95*sTi_planktic, '-', lw = 0.75, color = [0.7]*3)
# 		plot(Ti, -f95*sTi_planktic, '-', lw = 0.75, color = [0.7]*3)

		# axis([Ti[0], Ti[-1], -3.2, 4.2])
		axis([Ti[0], Ti[-1], -2.6, 3.3])

		legs = [
			legs_anderson,
			legs_fiebig,
			] + ([(legs_dvhlgb, legs_dvhlgb_dot)] if DRAW_DVHLGB else []) + [
			legs_ts,
			]
		legend1 = legend(
			legs,
			[_.get_label() if 'get_label' in dir(_) else _[0].get_label() for _ in legs],
			loc = 'lower left',
			fontsize = 8,
			frameon = False,
			bbox_to_anchor = (0., 0.),
			)

		legs = [
			legs_meinicke,
			legs_peral,
			]
# 			legs_anderson,
# 			legs_fiebig,
# 			] + ([legs_dvhlgb] if DRAW_DVHLGB else [])
		legend(
			legs,
			[_.get_label() for _ in legs],
			loc = 'upper right',
			fontsize = 8,
			labelspacing = 0.2,
			frameon = False,
			bbox_to_anchor = (0.99, 0.99),
		# 			handlelength = 1.53,
			handlelength = 3.3,
			ncols = 1,
			markerfirst = False,
			)

		gca().add_artist(legend1)

		xlabel('Temperature (°C)')
		ylabel('Temperature offset between regressions (°C)\nrelative to planktic foraminifer regression', labelpad = 10)
		xticks([0,10,20,30])
		gca().yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.0f}$' if x else '$0$'))

		savefig(f'plots/compare_regressions{"_WOA_for_cold_planktics" if TWOA_FOR_PLANKTICS_BELOW_ZERO else ""}{"_toggle_no_offset" if TWOA_FOR_PLANKTICS_WITHOUT_OFFSET else ""}.pdf')
		close(fig)


		# MARK Save planktic CSV

		with open('planktic_foram_summary.csv', 'w') as fid:
			fid.write('Sample,Species,T18,SE_T18,Twoa,SE_Twoa,D47_ICDES,SE_D47_ICDES,Ref')
			for r in sorted(data, key = lambda _: (_['Tiso_species'] if _['Type'] == 'planktic' else 100)):
				if r['Type'] != 'planktic' or r['Tiso_species_offset'] == 'all planktics':
					continue
				fid.write(f"\n{r['Sample']},{r['Species']},{r['Tiso_species']:.2f},{r['SE_Tiso_species']:.2f},{r['Twoa23']:.2f},{r['SE_Twoa23']:.2f},{r['D47_ICDES']:.4f},{r['SE_D47_ICDES']:.4f},{r['Ref']}")


ps, pi = ancova(*peral_planktic_data, *zip(*[_ for _ in zip(*meinicke_data) if not isnan(_[0])]))

with open('../output/saved_values.py', 'a') as fid:
	fid.write(f'''
planktic_ancova_p_slope = {ps}
planktic_ancova_p_intercept = {pi}
''')

if min(ps, pi) >= 0.05:
	print(f'  The two planktic data sets are statistically indistiguishable (p_slope = {ps:.3f}, p_intercept = {pi:.2f})')
else:
	print(f'  The two planktic data sets are statistically distiguishable (p_slope = {ps:.3f}, p_intercept = {pi:.2f})')

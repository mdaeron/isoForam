#! /usr/bin/env python3

import sys
from csv import DictReader
from pylab import *
from matplotlib.patches import Polygon
import netCDF4 as nc
from scipy.interpolate import interp1d
from matplotlib.scale import FuncScale
from tqdm import tqdm
from matplotlib import ticker
from scipy.stats import norm
from scipy.stats.qmc import MultivariateNormalQMC
from statsmodels.nonparametric.smoothers_lowess import lowess

sys.path.append('../output')
import saved_values

sys.path.append('../6_assign_atlas_T')
from woa import get_T_from_woa, get_fuzzy_T_from_woa

sys.path.append('../2_compile_d18O_data')
from KON97_intercepts import KON97_intercepts

style.use('../mydefault.mplstyle')

RUN_LOWESS = False

if RUN_LOWESS:
	Ni = 2**13 # 2**12 takes 5 mn on 2021 MacBook Pro
	seed = 47

CRAMER_COLOR_T = (1, .5, .25)
CRAMER_COLOR_SW = (.25, .5, 1)

# IODP SITE 1410
lat = 41+19.6987/60
lon = -49-10.1995/60
depth = 3387.5

GET_BOTTOM_δ18Osw_FOR_IODP_SITE_1410 = False

if GET_BOTTOM_δ18Osw_FOR_IODP_SITE_1410:

	d18Osw_model_sigma = 0.1

	ds = nc.Dataset('../5_assign_d18Osw/breitkreutz_2018/D18O_Breitkreuz_et_al_2018.nc')

	lats = ds['lat_1deg_center'][:,0]
	lons = ds['lon_1deg_center'][0,:]
	zc = ds['depth_center'][:]
	ze = ds['depth_edge'][:]
	depths = -concatenate((ze[:1], zc[:]))

	d18o = ds['D18O_1deg'][:,:,:,:] # month, depth, lat, lon
	d18o = concatenate((d18o[:,:1,:,:], d18o[:,:,:,:]), axis = 1)

	d18o = d18o.filled(fill_value = nan)
	d18o_mask = (isnan(d18o)).astype(int)

	glon, glat = meshgrid(lons, lats)
	gx = cos(glon * pi / 180) * cos(glat * pi / 180)
	gy = sin(glon * pi / 180) * cos(glat * pi / 180)
	gz = sin(glat * pi / 180)

	print('Extracting seawater d18O...')

	x = cos(lon * pi / 180) * cos(lat * pi / 180)
	y = sin(lon * pi / 180) * cos(lat * pi / 180)
	z = sin(lat * pi / 180)
	sqdistance = (gx-x)**2 + (gy-y)**2 + (gz-z)**2

	i = [i for i, _ in enumerate(depths) if _ >= depth][0]

	sqdistance += d18o_mask[0,i,:,:] * 10
	j,k = [int(_) for _ in where(sqdistance == sqdistance.min())]

	fig = figure(figsize = (8,4))
	ax1, ax2 = subplot(121), subplot(122)
	subplots_adjust(.15, .15, .95, .9, .25)

	X, Y, M = depths[:], d18o[:,:,j,k], d18o_mask[0,:,j,k]
	X, Y = X[M<1], Y[:,M<1]

	maxdepth = X[-1]

	d18values, d18values_500m = [], []
	for y in Y:
		sca(ax1)
		plot(y, -X, 'b-', alpha = .1)
		f = interp1d(X,y)
		d18values += [f(depth)]

	kw = dict(elinewidth = 1.5, alpha = 1, capsize = 5, marker = 'None', ls = 'None', capthick = 1.5)

	d18values = array(d18values)
	d18, sd18 = d18values.mean(), (d18Osw_model_sigma**2 + d18values.std(ddof = 1)**2)**.5
	sca(ax1)
	errorbar(d18, -depth, None, 1.96*sd18, ecolor = 'b', **kw)
	xlabel('d18Osw')
	ylabel('depth')
	text(.5, .97, f'IODP-1410', va = 'top', ha = 'center', transform = fig.transFigure)

	savefig(f'd18Osw.pdf')
	close(fig)

	d18Osw_modern = d18 * 1.
	SE_d18Osw_modern = sd18 * 1.
	print(f'd18Osw_modern = {d18Osw_modern:+.2f} ‰')
else:
	d18Osw_modern = 0.16


GET_BOTTOM_TEMPERATURE_FOR_IODP_SITE_1410 = False

if GET_BOTTOM_TEMPERATURE_FOR_IODP_SITE_1410:
	T_modern, SE_T_modern, depth, method = get_fuzzy_T_from_woa(lat, lon, depth, 'IODP-1410', plotdir = '.')

# PLANKTIC Δ47 REGRESSION
_A_   = saved_values.planktic_reg_slope
_sA_  = saved_values.planktic_reg_slope_se
_B_   = saved_values.planktic_reg_intercept
_T0_  = saved_values.planktic_reg_Tzero
_B0_  = saved_values.planktic_reg_Tzero_intercept
_sB0_ = saved_values.planktic_reg_Tzero_intercept_se

# DVHLGB Δ47 REGRESSION
_Ae_   = saved_values.dvhlgb_reg_slope
_sAe_  = 0.
_Be_   = saved_values.dvhlgb_reg_intercept
_T0e_  = 0.
_B0e_  = 0.
_sB0e_ = 0.

# ANDERSON Δ47 REGRESSION
_Aa_   = saved_values.mit_reg_slope
_sAa_  = 0.
_Ba_   = saved_values.mit_reg_intercept
_T0a_  = 0.
_B0a_  = 0.
_sB0a_ = 0.

# 18α LAW FOR CIBICIDOIDES
_18alpha_intercept_ = KON97_intercepts['Cibicidoides spp.']['avg']

# DATA FROM MECKLER ET AL. (2022) SUPPLEMENT
with open('meckler_2022_data.csv') as fid:
	data = list(DictReader(fid))

for r in data:
	r['D47_ICDES'] = float(r['D47 avg I-CDES (‰)'])
	r['SE_D47_ICDES'] = float(r['D47 SE (‰)'])
	r['Age'] = float(r['Age (Ma CENOGRID)'])
	r['T47_pub'] = float(r['Temp (°C)'])
	r['95CL_T47_pub'] = float(r['Temp 95% CI (°C)'])
	if r['d18O Cibicidoides (‰ VPDB)']:
		r['d18Ocibi'] = float(r['d18O Cibicidoides (‰ VPDB)'])
	else:
		r['d18Ocibi'] = nan

d18Oc = array([r['d18Ocibi'] for r in data])

X = array([r['Age'] for r in data])
D47 = array([r['D47_ICDES'] for r in data])
sD47 = array([r['SE_D47_ICDES'] for r in data])
T_pub = array([r['T47_pub'] for r in data])
eT_pub = array([r['95CL_T47_pub'] for r in data])



if RUN_LOWESS:
	dist = MultivariateNormalQMC(
		mean = X*0,
		cov = diag(X*0+1),
		seed = seed,
		)
	errors = dist.random(Ni)



# RECENT TEMPERATURE RECONSTRUCTION FROM CRAMER ET AL. (2011) SUPPLEMENT
with open('cramer_2011_S7.csv') as fid:
	cramerdata = list(DictReader(fid))

for r in cramerdata:
	r['T'] = float(r['Temperature'])
	r['SL'] = float(r['Sea level'])
	r['Age'] = float(r['Age'])

X_cramer_recent = array([r['Age'] for r in cramerdata])
T_cramer_recent = array([r['T'] for r in cramerdata])
d18Osw_cramer_recent = array([r['SL'] * -0.011 + d18Osw_modern for r in cramerdata])

recalculated_d18Oc_cramer_recent = -(T_cramer_recent - 16.1)/4.76 + d18Osw_cramer_recent - 0.27
recalculated_alpha18_cramer_recent = (1000+recalculated_d18Oc_cramer_recent) / (1000+d18Osw_cramer_recent) * 1.03092
recalculated_T_cramer_recent = 18030 / (1000*log(recalculated_alpha18_cramer_recent) - _18alpha_intercept_) - 273.15

# ANCIENT TEMPERATURE & δ18Osw RECONSTRUCTION FROM CRAMER ET AL. (2011) SUPPLEMENT
with open('cramer_2011_S4.csv') as fid:
	cramerdata = list(DictReader(fid))

for r in cramerdata:
	r['Age'] = float(r['Age'])
	if r['Temperature (long)']:
		r['T'] = float(r['Temperature (long)'])
	else:
		r['T'] = nan
	if r['Temperature min (long)']:
		r['T_min'] = float(r['Temperature min (long)'])
	if r['Temperature max (long)']:
		r['T_max'] = float(r['Temperature max (long)'])
# 	if r['SLice (long)']:
# 		r['d18Osw'] = float(r['SLice (long)']) * -0.011 - 0.3
	if r['SLice (long)']:
		r['d18Osw_long'] = float(r['SLice (long)']) * -0.011 - 0.19
	if r['SLice max (long)']:
		r['d18Osw_long_max'] = float(r['SLice max (long)']) * -0.011 - 0.19
	if r['SLice min (long)']:
		r['d18Osw_long_min'] = float(r['SLice min (long)']) * -0.011 - 0.19
	if r['SLice']:
		r['d18Osw'] = float(r['SLice']) * -0.011 - 0.19
	else:
		r['d18Osw'] = nan
# 	r['d18Obf_recalc'] = -(r['T'] - 16.1) / 4.76 + (r['d18Osw']-0.27)
# 	r['d18Osw_recalc'] = (
# 		(1000 + r['d18Obf_recalc'])
# 		* 1.03092
# 		/ exp((18030/(r['T']+273.15) + _18alpha_intercept_)/1000)
# 		) - 1000

cramerdata = [r for r in cramerdata if r['Age'] <= 70]

X_cramer_ancient = array([r['Age'] for r in cramerdata])
T_cramer_ancient = array([r['T'] for r in cramerdata])
T_cramer_ancient_min = array([r['T_min'] for r in cramerdata])
T_cramer_ancient_max = array([r['T_max'] for r in cramerdata])
d18Osw_cramer_ancient = array([r['d18Osw_long'] for r in cramerdata])
d18Osw_cramer_ancient_max = array([r['d18Osw_long_max'] for r in cramerdata])
d18Osw_cramer_ancient_min = array([r['d18Osw_long_min'] for r in cramerdata])
d18Osw_cramer_ancient_short = array([r['d18Osw'] for r in cramerdata])
# d18Osw_cramer_ancient_recalc = array([r['d18Osw_recalc'] for r in cramerdata])

recalculated_d18Oc_cramer_ancient = -(T_cramer_ancient - 16.1)/4.76 + d18Osw_cramer_ancient - 0.27
recalculated_alpha18_cramer_ancient = (1000+recalculated_d18Oc_cramer_ancient) / (1000+d18Osw_cramer_ancient) * 1.03092
recalculated_T_cramer_ancient = 18030 / (1000*log(recalculated_alpha18_cramer_ancient) - _18alpha_intercept_) - 273.15

X_cramer = array([x for x in X_cramer_recent] + [x for x in X_cramer_ancient])
T_cramer = array([t for t in T_cramer_recent] + [t for t in T_cramer_ancient])
recalculated_T_cramer = array([t for t in recalculated_T_cramer_recent] + [t for t in recalculated_T_cramer_ancient])
recalculated_d18Oc_cramer = array([t for t in recalculated_d18Oc_cramer_recent] + [t for t in recalculated_d18Oc_cramer_ancient])

fT_cramer = interp1d(X_cramer, T_cramer)
recalculated_fT_cramer = interp1d(X_cramer, recalculated_T_cramer)
recalculated_fd18Oc_cramer = interp1d(X_cramer, recalculated_d18Oc_cramer)

fig = figure(figsize = (5,5))
Dd18c = d18Oc - recalculated_fd18Oc_cramer(X)
xDd18c = array([age for _, age in zip(Dd18c, X) if not isnan(_)])
Dd18c = array([_ for _ in Dd18c if not isnan(_)])

# ax = subplot(121)

plot(xDd18c, Dd18c, 'r+')
xlabel('Age')
ylabel('δ$^{18}O_c$ difference between Meckler et al. (2022) and Cramer at al. (2011)')
title(f'avg = {Dd18c.mean():+0.2f} ± {Dd18c.std(ddof = 1) / sqrt(len(Dd18c)):0.2f} ‰')

# ax2 = subplot(122, sharey = ax)
# hist(Dd18c, orientation = 'horizontal')
# 
savefig('d18Oc_meckler_vs_cramer.pdf')
close(fig)

for figtitle, label, calib in (
	('planktic_regression', 'new planktic calibration (this study)', (_A_, _sA_, _B_, _T0_, _B0_, _sB0_)),
	('devils_laghetto_regression', '“Devils Laghetto” calibration\n(Anderson et al., 2021; Fiebig et al., 2021)', (_Ae_, _sAe_, _Be_, _T0e_, _B0e_, _sB0e_)),
	('anderson_regression', 'MIT calibration (Anderson et al., 2021)', (_Aa_, _sAa_, _Ba_, _T0a_, _B0a_, _sB0a_)),
	('fiebig_regression', 'Fiebig et al. (2021) calibration', 'fiebig'),
	):

	if isinstance(calib, tuple):
		A, sA, B, T0, B0, sB0 = calib

		# COMPUTE T47 USING NEW CALIBRATION
		T47_new = ((D47 - B) / A)**-0.5 - 273.15
		sT47_naive = sD47 * 0.5 / A * ((D47 - B) / A)**-1.5
		if sA:
			sD47_fromcalib = (sA**2 * ((T47_new+273.15)**-2 - (T0+273.15)**-2)**2 + sB0**2)**0.5
			sT47_fromcalib = sD47_fromcalib * 0.5 / A * ((D47 - B) / A)**-1.5
			sT47_new_full = (sT47_naive**2 + sT47_fromcalib**2)**0.5
		else:
			sT47_new_full = sT47_naive
	else:
		_ = linspace(-10,40,51) + 273.15

		_fiebig_ = interp1d(
		1.038*(-5.897/_ - 3.521e3/_**2 + 2.391e7/_**3 - 3.541e9/_**4) + 0.1856,
		_ - 273.15,
		)
		T47_new = _fiebig_(D47)
		sT47_new_full = sD47 * (_fiebig_(D47-0.01) - _fiebig_(D47+0.01))/0.02		
	
	fig = figure(figsize = (7,6))

	x1, x2, x3, x4 = .08, .51, .54, .97
	y1, y2, y3, y4 = .1, .49, .51, .88

	ax_ul = fig.add_axes((x1, y3, x2-x1, y4-y3))
	ax_ur = fig.add_axes((x3, y3, x4-x3, y4-y3), sharex = ax_ul)
	ax_ll = fig.add_axes((x1, y1, x2-x1, y2-y1), sharex = ax_ul)
	ax_lr = fig.add_axes((x3, y1, x4-x3, y2-y1), sharex = ax_ul)

	T47_minus_T18 = {}
	eT47_minus_T18 = {}
	for ax1, ax2, _T_, _eT_, oldnew in (
		(ax_ul, ax_ll, T_pub, eT_pub, 'old'),
		(ax_ur, ax_lr, T47_new, 1.96*sT47_new_full, 'new'),
		):

# 		f = lambda t: exp((18030 / (t + 273.15) + _18alpha_intercept_)/1000)
# 		alpha18 = f(_T_)
# 		deriv_alpha18 = f(_T_+0.5) - f(_T_-0.5)
# # 	
# 		d18sw_with_our_alpha = (1000 + d18Oc) / alpha18 * 1.03092 - 1000
# 		ed18sw_with_our_alpha = - _eT_ * deriv_alpha18 * 1e3

		d18sw = 0.27 + d18Oc + (_T_ - 16.1) / 4.76 # Lynch-Stielglitz (1999) from Cramer (2011)
		ed18sw = _eT_ / 4.76

		kw = dict(
			ls = 'None',
			marker = 'o',
			mfc = 'w',
			mec = 'k',
			ms = 5,
			mew = 1,
			zorder = 101,
			)

		kweb = dict(
			ls = 'None',
			marker = 'None',
			ecolor = 'k',
			elinewidth = 1,
			capthick = 0,
			capsize = 0,
			alpha = 0.3,
			zorder = 100,
			)

		ax1.errorbar(X, _T_, _eT_, **kweb)
		ax1.plot(X, _T_, **kw)


		lowess_color = (0,0,0,0.15)
		bw = 20
		bandwidth_for_lowess_comparison = bw
		it = 0
		frac = bw/(X.max()-X.min())
		if RUN_LOWESS:
			Xi = linspace(0, ceil(X.max()), 2*int(ceil(X.max()))+1)

			Tqmc = _T_ + errors * _eT_/1.96
			Ti = array([lowess(T, X, frac = frac, it = it, xvals = Xi) for T in tqdm(Tqmc)])
			Tmin, Tmax = quantile(Ti, [0.025, 0.975], axis = 0)

			savetxt(f'lowess_results/Tmin_{figtitle}_{oldnew}.csv', Tmin, delimiter = ',')
			savetxt(f'lowess_results/Tmax_{figtitle}_{oldnew}.csv', Tmax, delimiter = ',')
			savetxt('lowess_results/Xi.csv', Xi, delimiter = ',')
		else:
			Tmin = loadtxt(f'lowess_results/Tmin_{figtitle}_{oldnew}.csv', delimiter = ',')
			Tmax = loadtxt(f'lowess_results/Tmax_{figtitle}_{oldnew}.csv', delimiter = ',')
			Xi = loadtxt(f'lowess_results/Xi.csv', delimiter = ',')

# 		ax1.fill_between(Xi, Tmin, Tmax, color = lowess_color, lw = 0)		

		T47_minus_T18[oldnew] = _T_ - fT_cramer(X)
		eT47_minus_T18[oldnew] = _eT_
		
		if figtitle == 'planktic_regression' and oldnew == 'new':
			savetxt('lowess_results/T47_minus_T18.csv', T47_minus_T18['new'], delimiter = ',')
			savetxt('lowess_results/eT47_minus_T18.csv', eT47_minus_T18['new'], delimiter = ',')

		ax2.errorbar(X, d18sw, ed18sw, **kweb)
		ax2.plot(X, d18sw, **kw)
# 		ax2.plot(X, d18sw_with_our_alpha, 'c+', zorder = 200)

		if RUN_LOWESS:
			d18swqmc = d18sw + errors * ed18sw/1.96
			d18swi = array([lowess(_, X, frac = frac, it = it, xvals = Xi) for _ in tqdm(d18swqmc)])
			d18swmin, d18swmax = quantile(d18swi, [0.025, 0.975], axis = 0)

			savetxt(f'lowess_results/d18swmin_{figtitle}_{oldnew}.csv', d18swmin, delimiter = ',')
			savetxt(f'lowess_results/d18swmax_{figtitle}_{oldnew}.csv', d18swmax, delimiter = ',')
		else:
			d18swmin = loadtxt(f'lowess_results/d18swmin_{figtitle}_{oldnew}.csv', delimiter = ',')
			d18swmax = loadtxt(f'lowess_results/d18swmax_{figtitle}_{oldnew}.csv', delimiter = ',')

# 		ax2.fill_between(Xi, d18swmin, d18swmax, color = lowess_color, lw = 0)

		if figtitle == 'planktic_regression':
			post45Ma_d18sw, SE_post45Ma_d18sw = zip(*[(d, ed/1.96) for a, d, ed in zip(X, d18sw, ed18sw) if a > 45 and not isnan(d)])
			w = array([_**-2 for _ in SE_post45Ma_d18sw])
			w /= w.sum()
			avg_post45Ma_d18sw = (array(post45Ma_d18sw) * w).sum()
			SE_avg_post45Ma_d18sw = ((array(SE_post45Ma_d18sw) * w)**2).sum()**.5
# 			print(f'{figtitle}: ({oldnew} calib): {avg_post45Ma_d18sw:.2f} ± {SE_avg_post45Ma_d18sw:.2f} (1SE)')
		
			post45Ma_d18sw_cramer = array([d for d,a in zip(d18Osw_cramer_ancient, X_cramer_ancient) if a > 45])
			avg_post45Ma_d18sw_cramer = post45Ma_d18sw_cramer.mean()
			SE_avg_post45Ma_d18sw_cramer = post45Ma_d18sw_cramer.std(ddof = 1) / len(post45Ma_d18sw_cramer)**.5
# 			print(f'Avg d18Osw > 45 Ma of Cramer et al. (2011) = {avg_post45Ma_d18sw_cramer:+.3f} ± {SE_avg_post45Ma_d18sw_cramer:.3f} (1SE)')
	
# 		d18sw_old_avg = mean([sw for sw, _ in zip(d18sw, X) if _>45 and ~isnan(sw)])
# 		print(d18sw_old_avg)
# 		ax2.plot(65, d18sw_old_avg, 'k>', ms= 9, mew = 1, mfc = (.5,.75,1))
	
		kwline = dict(
			ls = '-',
			marker = 'None',
			lw = 1,
			color = CRAMER_COLOR_T,
			)

		kwfill = dict(
			ls = '-',
			lw = .5,
			color = [(.75 + .25 * _) for _ in CRAMER_COLOR_T],
			zorder = 0,
			)

# 		kw_recalc = dict(
# 			ls = '-',
# 			marker = 'None',
# 			lw = 1,
# 			color = CRAMER_COLOR_T,
# 			)

		ax1.plot(X_cramer_recent, T_cramer_recent, **kwline)
		ax1.fill_between(
			X_cramer_ancient,
			T_cramer_ancient_max,
			T_cramer_ancient_min,
			**kwfill,
			)
		ax1.plot(X_cramer_ancient, T_cramer_ancient, label = 'Cramer et al. (2011)', **kwline)

# 		ax1.plot(X_cramer_recent, recalculated_T_cramer_recent, 'c-', alpha = .4)
# 		ax1.plot(X_cramer_ancient, recalculated_T_cramer_ancient, 'c-', alpha = .4)

		kw['color'] = CRAMER_COLOR_SW

# 		ax2.plot(X_cramer_recent, d18Osw_cramer_recent, **kw)

# 		ax2.plot(X_cramer_ancient, d18Osw_cramer_ancient_short, alpha = 0.35, **kwline)

		ax2.fill_between(
			X_cramer_ancient,
			d18Osw_cramer_ancient_max,
			d18Osw_cramer_ancient_min,
			**kwfill,
			)
		ax2.plot(X_cramer_ancient, d18Osw_cramer_ancient, label = 'Cramer et al. (2011)', **kwline)
	
# 		ax2.plot(X_cramer_ancient, d18Osw_cramer_ancient_recalc, label = 'Cramer et al. (2011)', **kw_recalc)

		ax1.axis([0, 66, -6, 26])
		ax2.axis([0, 66, 2.3, -3.3])

		ax2.set_xlabel('Age (Ma)')
	
# 		ax1.set_xscale(FuncScale(ax1, (lambda _: _**.5, lambda _: _**2)))
	
		ax2.axhline(0, xmax = 0.14, color = (.4,.8,0), lw = 1.5, dashes = (3,1,), zorder = 0, label = 'Modern')
		ax2.axhline(1, xmax = 0.14, color = (0,.5,1), lw = 1.5, dashes = (3,1), zorder = 0, label = 'LGM')
# 		ax2.axhline(-0.9, xmin = 0.65, color = 'k', lw = 1, dashes = (8,2,2,2))

		ax2.legend(loc = 'upper left', fontsize = 9, frameon = False, labelspacing = 0.2, handlelength = 2.6)
		ax1.legend(loc = 'upper left', fontsize = 9, frameon = False, labelspacing = 0.2, handlelength = 2.6)


	ax_ul.set_ylabel('Reconstructed deep ocean T (°C)')
	ax_ll.set_ylabel('$δ^{18}O_{sw}$ (‰ VSMOW)')
	ax_ur.tick_params(length = 0, labelleft = False)
	ax_lr.tick_params(axis = 'y', length = 0, labelleft = False)
	ax_ur.tick_params(length = 0, labelbottom = False)
	ax_ul.tick_params(axis = 'x', length = 0, labelbottom = False)


	ax_ul.set_title('Meckler et al. (2022) data with original\ncalibration of Meinicke et al. (2021)', size = 10, pad = 12)
	ax_ur.set_title(f'Meckler et al. (2022) data, reprocessed\nwith {label}', size = 10, pad = 12)

	fig.savefig(f'meckler_{figtitle}.pdf')
	close(fig)

	fig = figure(figsize = (5,5))

	x1, x2, x3, x4 = .1, .49, .51, .95
	y1, y2, y3, y4 = .1, .6, .66, .9

	ax_ul = fig.add_axes((x1, y3, x2-x1, y4-y3))
	ax_ur = fig.add_axes((x3, y3, x4-x3, y4-y3), sharex = ax_ul)
	ax_lr = fig.add_axes((x3, y1, x4-x3, y2-y1), sharex = ax_ul)
	ax_ll = fig.add_axes((x1, y1, x2-x1, y2-y1), sharex = ax_ul, sharey = ax_lr)

	sca(ax_ll)
	ax_ll.errorbar(T47_minus_T18['old'], X, None, eT47_minus_T18['old'], **kweb)
	ax_ll.plot(T47_minus_T18['old'], X, **kw)
	ax_ll.axvline(0, **kwline)
	ax_ll.set_xlim([-15,15])
	ax_ll.set_ylim(ax_ll.get_ylim()[::-1])
	xlabel('$T_{47} - T_{18}$ (°C)')
	ylabel('Age (Ma)')

	sca(ax_lr)
	ax_lr.errorbar(T47_minus_T18['new'], X, None, eT47_minus_T18['new'], **kweb)
	ax_lr.plot(T47_minus_T18['new'], X, **kw)
	ax_lr.axvline(0, **kwline)
	ax_lr.set_xlim([-15,15])
	xlabel('$T_{47} - T_{18}$ (°C)')
	
	sca(ax_ul)
	xlabel('$T_{47} - T_{18}$ (°C)')
	ylabel('KDE')
	yticks([])
	ti = linspace(*ax_ll.get_xlim(), 1000)
	pdf = ti*0
	for T,sT in zip(T47_minus_T18['old'], eT47_minus_T18['new']/1.96):
		pdf += norm(loc = T, scale = sT).pdf(ti)
	axvline(0, **kwline)
	plot(ti, pdf, 'k-')
	title('Using\nMeinicke et al. (2021) calibration', size = 8)

	sca(ax_ur)
	xlabel('$T_{47} - T_{18}$ (°C)')
	yticks([])
	ti = linspace(*ax_ll.get_xlim(), 1000)
	pdf = ti*0
	for T,sT in zip(T47_minus_T18['new'], eT47_minus_T18['new']/1.96):
		pdf += norm(loc = T, scale = sT).pdf(ti)
	axvline(0, **kwline)
	plot(ti, pdf, 'k-')
	title(f'Using\n{label}', size = 8)

	savefig(f'Tdiff_{figtitle}.pdf')
	close(fig)


if RUN_LOWESS:

	Xi = linspace(0, ceil(X.max()), 2*int(ceil(X.max()))+1)
	bandwidths = linspace(5,25,81)
	lowbound_T47_minus_T18 = zeros((bandwidths.size, Xi.size))
	highbound_T47_minus_T18 = zeros((bandwidths.size, Xi.size))

	it = 0
	for k,bw in tqdm(enumerate(bandwidths), total = bandwidths.size):
		frac = bw/(X.max()-X.min())
		Yqmc = T47_minus_T18['new'] + errors * eT47_minus_T18['new']/1.96
		Yi = array([lowess(Y, X, frac = frac, it = it, xvals = Xi) for Y in Yqmc])
		Ymin, Ymax = quantile(Yi, [0.025, 0.975], axis = 0)
		lowbound_T47_minus_T18[k,:] = Ymin
		highbound_T47_minus_T18[k,:] = Ymax

	savetxt('lowess_results/lowbound_T47_minus_T18.csv', lowbound_T47_minus_T18, delimiter = ',')
	savetxt('lowess_results/highbound_T47_minus_T18.csv', highbound_T47_minus_T18, delimiter = ',')
	savetxt('lowess_results/Xi.csv', Xi, delimiter = ',')
	savetxt('lowess_results/bandwidths.csv', bandwidths, delimiter = ',')

else:
	lowbound_T47_minus_T18 = loadtxt('lowess_results/lowbound_T47_minus_T18.csv', delimiter = ',')
	highbound_T47_minus_T18 = loadtxt('lowess_results/highbound_T47_minus_T18.csv', delimiter = ',')
	Xi = loadtxt('lowess_results/Xi.csv', delimiter = ',')
	bandwidths = loadtxt('lowess_results/bandwidths.csv', delimiter = ',')

T47_minus_T18 = loadtxt('lowess_results/T47_minus_T18.csv', delimiter = ',')
eT47_minus_T18 = loadtxt('lowess_results/eT47_minus_T18.csv', delimiter = ',')

with open('westerhold_2020.csv') as fid:
	whdata = list(DictReader(fid))

for _ in whdata:
	_['Age'] = float(_['Age'])
	_['d18Oc'] = float(_['d18Oc'])

fig = figure(figsize = (5.6,5))
show_bw = 10.
k_show_bw = where(bandwidths == show_bw)[0][0]

x1, x2, x3, x4 = .11, .84, .86, .98
y00, y0, y1, y2, y3, y4 = .09, .31, .35, .60, .64, .98

ax_ul = fig.add_axes((x1, y3, x2-x1, y4-y3))
ax_ur = fig.add_axes((x3, y3, x4-x3, y4-y3), sharey = ax_ul)
ax_ll = fig.add_axes((x1, y1, x2-x1, y2-y1), sharex = ax_ul)
ax_lr = fig.add_axes((x3+.02, y1, (x4-x3)/15, y2-y1))
ax_wh = fig.add_axes((x1, y00, x2-x1, y0-y00), sharex = ax_ul)

setp(ax_ul.get_xticklabels(), visible = False)
setp(ax_ll.get_xticklabels(), visible = False)
setp(ax_ur.get_yticklabels(), visible = False)
ax_ur.tick_params(length = 0)

age, bw = np.meshgrid(Xi, bandwidths)
z = lowbound_T47_minus_T18

sca(ax_ll)
hotzone = lowbound_T47_minus_T18.copy()
hotzone[hotzone<0] = 0
coldzone = highbound_T47_minus_T18.copy()
coldzone[coldzone>0] = 0
bothzones = hotzone + coldzone

from matplotlib.colors import CenteredNorm, ListedColormap, LinearSegmentedColormap

colors1 = cm.BuPu_r(linspace(0., 1, 256))
colors2 = cm.YlOrBr(linspace(0, 1, 256))
colors = vstack((colors1, colors2))
cmap = LinearSegmentedColormap.from_list('my_cmap', colors)
my_cmap = cmap(arange(cmap.N))
my_cmap[cmap.N//2-4:cmap.N//2+4, -1] = [.75,.5,.25,0,0,.25,.5,.75]
my_cmap = ListedColormap(my_cmap)

cf = contourf(age, bw, bothzones, levels = linspace(-3.3,3.3,64), cmap = my_cmap, norm = CenteredNorm(), extend = 'neither', algorithm = 'serial', antialiased = False)
# contour(age, bw, lowbound_T47_minus_T18, levels = [0], colors = 'k', linewidths = 0.7)
# cf_cold = contourf(age, bw, coldzone, levels = linspace(-3,0,10), cmap = 'PuOr_r', vmin = -4, vmax = 4, extend = 'neither')
ylabel('LOWESS bandwidth (Ma)', labelpad = 10, size = 8)
ax_ll.yaxis.set_major_locator(ticker.MultipleLocator(5))
text(.45, .95,
	'Lower bound (2.5 % quantile, in red) and\nupper bound (97.5 % quantile, in blue) of\nthe LOWESS regression of $ΔT_{47-18}$',
	transform = gca().transAxes,
	size = 8,
	va = 'top',
	ha = 'center',
	linespacing = 1.3,
	)

cbarticks = [-3,0,3]
cbar = colorbar(cf, cax = ax_lr, ticks = cbarticks, extend = 'neither')
cbar.ax.set_yticklabels([f'{t:+.0f} °C' if t else '0 °C' for t in cbarticks])
ax_lr.grid(False)


# cbarticks = [-3,-2,-1,0]
# cbar_cold = colorbar(cf_cold, cax = ax_urr, ticks = cbarticks, extend = 'both')
# cbar_cold.ax.set_yticklabels([f'{"+" if t>0 else ""}{t:.0f} °C' for t in cbarticks])

kw['ms'] = 4
kweb['ecolor'] = 'k'
kweb['alpha'] = 0.5
kweb['capsize'] = 2
kweb['capthick'] = 1
kweb['elinewidth'] = 1

sca(ax_ul)
errorbar(X, T47_minus_T18, eT47_minus_T18, **kweb)
plot(X, T47_minus_T18, **kw)
axhline(0, **{**kwline, **dict(color = 'k', dashes = (3,2))})

Ymin = lowbound_T47_minus_T18[k_show_bw,:]
Ymax = highbound_T47_minus_T18[k_show_bw,:]
fill_between(Xi, Ymax, Ymin, color = [.9]*3, lw = 0.8, label = f'95 % CL of LOWESS regression (bandwidth = {show_bw:.0f} Ma)', ec = [.7]*3)
# plot(Xi, Ymax, '-', color = , lw = 0.8)
# plot(Xi, Ymin, '-', color = [.6,.9,.6], lw = 0.8)

ax_ul.set_ylim([-16,12])
ylabel('$ΔT_{47-18}$ (°C)')
legend(fontsize = 8,loc = 'lower center', frameon = False)

max_age = axis()[1]

sca(ax_ur)
# xlabel('$T_{47} - T_{18}$ (°C)')
# ylabel('KDE')
xticks([])
ti = linspace(*ax_ul.get_ylim(), 1000)
pdf = ti*0
for T,sT in zip(T47_minus_T18, eT47_minus_T18/1.96):
	pdf += norm(loc = T, scale = sT).pdf(ti)
axhline(0, **{**kwline, **dict(color = 'k', dashes = (3,2))})
fill_betweenx(ti, pdf*0, pdf, color = 'k', lw = 0, alpha = 0.1)
plot(pdf, ti, 'k-', lw = 0.75)
axis([0, None, None, None])
text(.5, .05, 'Kernel\nDensity\nEstimation', ha = 'center', va = 'bottom', transform = gca().transAxes, size = 8, linespacing = 1.5)
gca().yaxis.set_major_formatter(lambda x, pos: f'{x:+.0f}' if x else '0')

sca(ax_wh)
whcolor = [.5]*3
X,Y = zip(*[(_['Age'], _['d18Oc']) for _ in whdata])
plot(X, Y, '-', lw  = 0.5, color = whcolor)
text(50, 3,
	'Westerhold et al. (2020)',
	size = 8,
	va = 'center',
	ha = 'center',
	linespacing = 1.5,
	color = [x*0.8 for x in whcolor],
	)
xlabel('Age (Ma)')
ylabel('$δ^{18}O_{VPDB}$ (‰)', labelpad = 6)
gca().yaxis.set_major_locator(ticker.MultipleLocator(2))
_ax_ = axis()
dims = [_ax_[0], max_age, _ax_[3]+1, _ax_[2]]
axis(dims)

geoltime = [
	0,
	'Pleistocene',
	2.58,
	'Plio',
	5.333,
	'Miocene',
	23.03,
	'Oligocene',
	33.9,
	'Eocene',
	56,
	'Paleocene',
	max_age
	]
h = 0.12
for t in geoltime[2:-2:2]:
	axvline(t, color = 'k', lw = 1, ymax = h)
# 	axvline(t, color = [.8]*3, lw = 1, zorder = -100)
axhline(dims[2] + (dims[3]-dims[2])*h, color = 'k', lw = 1)
for k in range(6):
	if geoltime[2*k+1].startswith('Plei'):
		text(
			(geoltime[2*k] + geoltime[2*k+2])/2 + 2.4,
			dims[2] + (dims[3]-dims[2])*h*0.47 + 2.2,
			geoltime[2*k+1],
			va = 'center',
			ha = 'center',
			size = 6,
			rotation = -35,
			)
	else:
		text(
			(geoltime[2*k] + geoltime[2*k+2])/2,
			dims[2] + (dims[3]-dims[2])*h*0.47,
			geoltime[2*k+1],
			va = 'center',
			ha = 'center',
			size = 6,
			)
savefig('bwplot.pdf')
close(fig)

fig = figure(figsize = (3.5,3.5))
subplots_adjust(.15, .15, .95, .95)
plot(X_cramer_recent, T_cramer_recent, '-', lw = 0.8, color = [0.7]*3, zorder = 1)
# plot(X_cramer_ancient, T_cramer_ancient, 'k-', lw = 0.8)

Tmin_old = loadtxt(f'lowess_results/Tmin_planktic_regression_old.csv', delimiter = ',')
Tmax_old = loadtxt(f'lowess_results/Tmax_planktic_regression_old.csv', delimiter = ',')

Tmin_new = loadtxt(f'lowess_results/Tmin_planktic_regression_new.csv', delimiter = ',')
Tmax_new = loadtxt(f'lowess_results/Tmax_planktic_regression_new.csv', delimiter = ',')

fill_between(
	Xi,
	Tmin_old,
	Tmax_old,
	label = 'Meckler et al. (2022) with original calibration',
	color = [1, 0.7, 0, 0.4],
	ec = [1, 0.7, 0, 1],
	lw = 0.8,
	zorder = 100,
	)

fill_between(
	Xi,
	Tmin_new,
	Tmax_new,
	label = 'Meckler et al. (2022) with new calibration',
	color = [0.05, .5, 0.86, 0.3],
	ec = [0.05, .5, 0.86, 1],
	lw = 0.8,
	zorder = 200,
	)

fill_between(
	X_cramer_ancient,
	T_cramer_ancient_max,
	T_cramer_ancient_min,
	label = 'Cramer et al. (2011)',
	color = (0, 0, 0, 0.15),
	ec = (.5, .5, .5, 1),
	lw = 0.8,
	zorder = 0.9
	)

text(.95, .05, f'LOWESS bandwidth = {bandwidth_for_lowess_comparison} Ma', ha = 'right', va = 'bottom', transform = gca().transAxes, size = 9)
legend(fontsize = 8, loc = 'upper center')
xlabel('Age (Ma)')
ylabel('Reconstructed T (°C)')
axis([Xi[0], Xi[-1], -4, 27])
savefig('compare_lowess.pdf')
close(fig)

#! /usr/bin/env python3
'''
Compile foram peral_data from original sources.
'''

import numpy as np
from csv import DictReader
from pylab import *
import matplotlib.ticker as ticker

style.use('../mydefault.mplstyle')


def deltas_to_ratios(d13C_VPDB, d18O_VSMOW,
	R13_VPDB = 0.01118,  # (Chang & Li, 1990)
	R18_VSMOW = 0.0020052,  # (Baertschi, 1976)
	R17_VSMOW = 0.00038475,  # (Assonov & Brenninkmeijer, 2003, rescaled to R13_VPDB)
	LAMBDA_17 = 0.528,  # (Barkan & Luz, 2005)
	D17O = 0 # in permil
	):

	d13C_VPDB = np.asarray(d13C_VPDB)
	d18O_VSMOW = np.asarray(d18O_VSMOW)
	if d13C_VPDB.shape != d18O_VSMOW.shape:
		raise ValueError('d13C_VPDB and d18O_VSMOW must both be floats or both be arrays of the same shape.')

	R13 = R13_VPDB * (1 + d13C_VPDB/1e3)
	R18 = R18_VSMOW * (1 + d18O_VSMOW/1e3)
	R17 = np.exp(D17O/1e3) * R17_VSMOW * (1 + d18O_VSMOW/1e3)**LAMBDA_17
	
	R45 = 2 * R17 + R13
	R46 = 2 * R18 + 2 * R17 * R13 + R17**2
	R47 = 2 * R13 * R18 + 2 * R17 * R18 + R13 * R17**2
	
	return(R45, R46, R47)

def old_D47_to_D47_ICDES(d13, d18c, D47, sD47 = None):
	R47_of_VPDBCO2 = deltas_to_ratios(0, (1.03092*1.01025-1)*1000)[2]
	d18g = (1000 + d18c) * 1.01025 * 1.03092 - 1000 # vs VSMOW
	R47 = deltas_to_ratios(d13, d18g)[2]
	d47 = (R47 / R47_of_VPDBCO2 - 1) * 1000
	D47new = -0.038039 - 0.000183 * d47 + 0.942603 * D47
	if sD47 is None:
		return D47new
	sD47new = 0.942603 * sD47
	return D47new, sD47new

def recompute_pooled_repeatability(sX, N):

	sX = np.asarray(sX)
	N  = np.asarray(N)

	internal_SD = sX * N**.5
	sum_of_sq_residuals = (N-1) * internal_SD**2
	pooled_SD = (sum_of_sq_residuals.sum() / (N.sum()-len(N)))**.5

	return pooled_SD


### PERAL

species = {
	'bullo'      : 'G. bulloides',
	'CWuel'      : 'Cibicides wuellerstorfi',
	'inflata200' : 'G. inflata',
	'inflata250' : 'G. inflata',
	'inflata315' : 'G. inflata',
	'inflata355' : 'G. inflata',
	'inflata400' : 'G. inflata',
	'inflata450' : 'G. inflata',
	'menardi'    : 'G. menardii menardi',
	'orbulina'   : 'Orbulina universa',
	'pachyD'     : 'N. pachyderma (d.)',
	'pachyS'     : 'N. pachyderma (s.)',
	'ruber'      : 'G. ruber',
	'truncaS'    : 'G. truncatulinoides (s.)',
	'truncaD'    : 'G. truncatulinoides (d.)',
	'UviMed'     : 'Uvigerina mediterranea',
	}

lalo = {
	'2FPA1'         : (+43.67,   -2.00,  664),
	'MD00-2360'     : (-20.08, +112.67),
	'MD02-2577'     : (+28.84,  -86.67),
	'MD03-2680'     : (+61.06,  -24.55),
	'MD04-2720'     : (-49.13,  +71.36),
	'MD08-3179'     : (+37.86,  -30.30),
	'MD08-3182'     : (+52.71,  -35.94),
	'MD12-3401'     : (-44.69,  +80.40),
	'MD12-3426'     : (+19.73, +114.61),
	'MD95-2014'     : (+60.59,  -22.08),
	'MOCOSED'       : (+73.04,  -11.93, 1839),
	'SU90-03'       : (+40.05,  -30.00),
	### Piasecki et al. (2019)
	'13MC-G'        : (+24.37,  -83.24,  348),
	'19MC-G'        : (+24.42,  -83.21,  173),
	'50MC-G'        : (+24.41,  -83.22,  198),
	'53MC-G'        : (+24.38,  -83.23,  302),
	'89MC-G'        : (+24.56,  -79.24,  353),
	'94MC-G'        : (+24.57,  -79.23,  259),
	'GS06-144-19'   : (+63.83,   +5.27,  830),
	'GS07-150-17-2' : ( -4.47,  -37.21, 1000),
	'GS07-150-22-1' : ( -4.33,  -37.16,  598),
	'MP43-BC'       : (+39.72,  +16.97,  246),
	'MP46-MC'       : (+39.54,  +17.25,  582),
	'SO213-54-4'    : (-43.72, -120.67, 3840),
	'SO213-71-2'    : (-45.58, -157.90,  689),
	}

with open('../0_reprocess_peral/peral_samples.csv') as f:
	peral_data = [{k: float(r[k]) if k != 'Sample' else r[k] for k in r} for r in DictReader(f)]

for r in peral_data:
	r['Species'] = species[r['Sample'].split('_')[1]]
	r['Type'] = 'benthic' if r['Species'] in ['Cibicides wuellerstorfi', 'Uvigerina mediterranea'] else 'planktic'
	r['Site'] = r['Sample'].split('_')[0]
	if r['Type'] == 'benthic':
		r['Lat'], r['Lon'], r['Depth'] = lalo[r['Site']]
	else:
		r['Lat'], r['Lon'] = lalo[r['Site']][:2]
	r['Ref'] = 'Peral et al. (2018)'

print(f'Processed {len(peral_data)} samples from Peral et al. (2018).')


### PIASECKI

with open('../1_compile_D47_data/piasecki_sites.csv') as f:
	piasecki_sites = [r for r in DictReader(f)]

with open('../1_compile_D47_data/piasecki_samples.csv') as f:
	piasecki_data = [r for r in DictReader(f)]

for r in piasecki_data:
# 	r['Sample'] = f"{r['Site']}_{r['Species'].replace(' ','_')}"
	r['Type'] = 'benthic'
	r['Lat'], r['Lon'], r['Depth'] = lalo[r['Site']] if r['Site'] in lalo else ('', '', '')
	r['Tatlas'] = [_ for _ in piasecki_sites if _['Site'] == r['Site']][0]['T']
	r['sTatlas'] = 0.5
	r['Ref'] = 'Piasecki et al. (2019)'
	r['N'] = int(r['N'])
# 	r['sD47pub'] = float(r['sD47pub'])
	r['D47_ICDES'], r['SE_D47_ICDES'] = float(r['D47']), float(r['sD47'])
	
sigma = recompute_pooled_repeatability(
	[r['SE_D47_ICDES'] for r in piasecki_data],
	[r['N']       for r in piasecki_data],
	)

with open('../output/saved_values.py', 'a') as fid:
	fid.write(f'\next_D47_sigma_piasecki = {sigma:.4f}')

for r in piasecki_data:
	r['SE_D47_ICDES'] = sigma/r['N']**.5

print(f'Processed {len(piasecki_data)} samples from Piasecki et al. (2019).')


### MEINICKE

with open('../1_compile_D47_data/meinicke_samples.csv') as f:
	meinicke_data = [r for r in DictReader(f)]


for r in meinicke_data:
	r['d18Oc'] = f"{(1000 + float(r['d18Oc']))/1.03092 - 1000:.2f}"
	r['sd18Oc'] = f"{float(r['sd18Oc'])/1.03092:.2f}"
	r['Type'] = 'planktic'
	r['N'] = int(r['N'])
	r['Ref'] = 'Meinicke et al. (2020)'
	r['D47_ICDES'], r['SE_D47_ICDES'] = old_D47_to_D47_ICDES(float(r['d13C']), float(r['d18Oc']), float(r['pubD47']),float(r['SE_pubD47']))
	

R = 1000 * array([(float(_['reprocessedD47']) - _['D47_ICDES']) for _ in meinicke_data])

rmse = (R**2).mean()**.5
with open('../output/saved_values.py', 'a') as fid:
	fid.write(f'\nEqA7_rmse = {rmse/1e3:.4f}')


fig = figure(figsize = (4,3))
subplots_adjust(.05, .25, .95, .95)
hist(R, bins = linspace(-3,6,18+1), histtype = 'stepfilled', color = [.75]*3, lw = 0, zorder = 2)
hist(R, bins = linspace(-3,6,18+1), histtype = 'step', color = 'k', lw = 0.8)
xlabel('Difference (ppm) between reprocessed and converted\nΔ$_{47}$ (I-CDES) values for the Meinicke et al. (2020) data set', labelpad = 10)
yticks([])
axis([-2.5, 5.5, 0, None])
gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.0f}$' if x else '$0$'))
savefig('EqA7_performance.pdf')
close(fig)


for r in meinicke_data:
	r['D47_ICDES'] = float(r['reprocessedD47'])

sigma = recompute_pooled_repeatability(
	[r['SE_D47_ICDES']  for r in meinicke_data],
	[r['N']             for r in meinicke_data],
	)

for r in meinicke_data:
	r['SE_D47_ICDES'] = sigma/r['N']**.5

with open('../output/saved_values.py', 'a') as fid:
	fid.write(f'\next_D47_sigma_meinicke = {sigma:.4f}')

print(f'Processed {len(meinicke_data)} samples from Meinicke et al. (2020).')


### SAVE COMPILATION

with open('../1_compile_D47_data/foram_D47_compilation.csv', 'w') as f:
	f.write('Sample,N,Species,Type,Site,Depth,Lat,Lon,Tatlas,sTatlas,d13C_VPDB,SE_d13C_VPDB,d18O_VPDB,SE_d18O_VPDB,D47_ICDES,SE_D47_ICDES,Ref')
	for r in peral_data:
		f.write(f"\n{r['Sample']},{int(r['N'])},{r['Species']},{r['Type']},{r['Site']},{r['Depth'] if 'Depth' in r else ''},{r['Lat']},{r['Lon']},{r['Tatlas']},{r['sTatlas']}")
		f.write(f",{r['d13C_VPDB']},{r['SE_d13C_VPDB']},{r['d18O_VPDB']},{r['SE_d18O_VPDB']},{r['D47_ICDES']},{r['SE_D47_ICDES']},{r['Ref']}")
	for r in piasecki_data:
		f.write(f"\n{r['Sample']},{int(r['N'])},{r['Species']},{r['Type']},{r['Site']},{r['Depth'] if 'Depth' in r else ''},{r['Lat']},{r['Lon']},{r['Tatlas']},{r['sTatlas']}")
		f.write(f",{r['d13C']},{r['sd13C']},{r['d18Oc']},{r['sd18Oc']},{r['D47_ICDES']:.4f},{r['SE_D47_ICDES']:.4f},{r['Ref']}")
	for r in meinicke_data:
		f.write(f"\n{r['Sample']},{int(r['N'])},{r['Species']},{r['Type']},{r['Site']},{r['Depth'] if 'Depth' in r else ''},{r['Lat']},{r['Lon']},,")
		f.write(f",{r['d13C']},{r['sd13C']},{r['d18Oc']},{r['sd18Oc']},{r['D47_ICDES']:.4f},{r['SE_D47_ICDES']:.4f},{r['Ref']}")

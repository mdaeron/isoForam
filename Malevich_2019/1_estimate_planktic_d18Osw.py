#! /usr/bin/env python3
'''
Assign calcification temperatures.
'''

from scipy.interpolate import interp1d
from csv import DictReader
import netCDF4 as nc
from matplotlib import ticker
from pylab import *
from tqdm import tqdm

style.use('../mydefault.mplstyle')

import sys

sys.path.append('../4_species_depth')
from species_depth import species_depth

sys.path.append('../3_species_lookup')
from species_lookup import species_lookup

d18Osw_model_sigma = 0.1
Tsw_model_sigma = 1.0
GENERATE_PLOTS = False

### SEAWATER d18O GRIDDED MODEL

ds = nc.Dataset('../5_assign_d18Osw/breitkreutz_2018/D18O_Breitkreuz_et_al_2018.nc')

# print(ds)

# print('\nDIMENSIONS:')
# for dim in ds.dimensions:
# 	print(f'{dim:>24}: {ds.dimensions[dim]}')
# 	try:
# 		print(ds[dim][:])
# 	except:
# 		pass
# 
# print('\nVARIABLES:')
# for var in ds.variables:
# 	print(var)

lats = ds['lat_1deg_center'][:,0]
lons = ds['lon_1deg_center'][0,:]
zc = ds['depth_center'][:]
ze = ds['depth_edge'][:]
depths = -concatenate((ze[:1], zc[:]))

d18o = ds['D18O_1deg'][:,:,:,:] # month, depth, lat, lon
d18o = concatenate((d18o[:,:1,:,:], d18o[:,:,:,:]), axis = 1)

T = ds['THETA_1deg'][:,:,:,:] # month, depth, lat, lon
T = concatenate((T[:,:1,:,:], T[:,:,:,:]), axis = 1)

d18o = d18o.filled(fill_value = nan)
d18o_mask = (isnan(d18o)).astype(int)

T = T.filled(fill_value = nan)
T_mask = (isnan(T)).astype(int)


glon, glat = meshgrid(lons, lats)
gx = cos(glon * pi / 180) * cos(glat * pi / 180)
gy = sin(glon * pi / 180) * cos(glat * pi / 180)
gz = sin(glat * pi / 180)

with open('malevich_2019_s1.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]
	
for _, r in enumerate(data):
	for k in r:
		if k not in ['corename', 'doi', 'foram', 'reference', 'species', 'subspecies']:
			r[k] = float(r[k])

for k,r in enumerate(data):
	r['Sample'] = f'{k:04.0f}'
	r['Species'] = r['species']
	r['Lat'] = r['latitude']
	r['Lon'] = r['longitude']
	r['Ref'] = r['reference']
	r['Type'] = 'planktic'

print('Extracting seawater d18O for planktics...')

with (
	open('d18Osw_estimates.csv', 'w') as fid,
	open('d18Osw_estimates_500m.csv', 'w') as fid_500m,
	):
	fid.write('Sample,Ref,d18Osw,sd18Osw,T,sT')
	fid_500m.write('Sample,Ref,d18Osw,sd18Osw,T,sT')

	d18Osw_residuals = []

	for r in tqdm(data):
		s, lon, lat, ref = r['Sample'], r['Lon'], r['Lat'], r['Ref']
# 		print(f'\tProcessing {s} from {ref}...')
		if r['Type'] == 'benthic':
			continue
		x = cos(lon * pi / 180) * cos(lat * pi / 180)
		y = sin(lon * pi / 180) * cos(lat * pi / 180)
		z = sin(lat * pi / 180)
		sqdistance = (gx-x)**2 + (gy-y)**2 + (gz-z)**2
# 		if 'Depth' in r:
# 			i = [i for i, _ in enumerate(depths) if _ >= depth][0]
# # 			print(depth, i, depths[i])
# 		else:
		sp = species_lookup[r['Species']]
		zmin = species_depth[sp]['zmin']
		zmax = species_depth[sp]['zmax']
# 		print(zmax, i, depths[i])
		i = [i for i, _ in enumerate(depths) if _ >= zmax][0]

		sqdistance += d18o_mask[0,i,:,:] * 10
		j,k = [_[0] for _ in where(sqdistance == sqdistance.min())]
	
		if GENERATE_PLOTS:
			fig = figure()
		X, Y, Tloc, M = depths[:], d18o[:,:,j,k], T[:,:,j,k], d18o_mask[0,:,j,k]
		X, Y, Tloc = X[M<1], Y[:,M<1], Tloc[:,M<1]
# 		print(X.shape, Y.shape)

		maxdepth = X[-1]

# 		if 'Depth' in r:
# 			zi = [depth]
# 		else:
		zi = linspace(zmin, zmax, int(zmax-zmin+1))
		zi_500m = linspace(0, min(500, maxdepth), int(min(500, maxdepth)+1))

		d18values, d18values_500m = [], []
		Tvalues, Tvalues_500m = [], []
		for y in Y:
			if GENERATE_PLOTS:
				plot(y, -X, 'b-', alpha = .3)
			f = interp1d(X,y)
			d18values += [f(zi)]
# 			if not 'Depth' in r:
			d18values_500m += [f(zi_500m)]

		for t in Tloc:
			f = interp1d(X,t)
			Tvalues += [f(zi)]
			Tvalues_500m += [f(zi_500m)]

		d18values = array(d18values)
		Tvalues = array(Tvalues)
		
		d18Osw_residuals += list(d18values.flatten() - d18values.mean())

# 		if 'Depth' in r:
# 			if r['Site'] not in bottoms:
# 				bottoms += [r['Site']]
# 				fid_bottom.write(f"\n{r['Site']},{d18values.mean():.2f},{(d18Osw_model_sigma**2 + d18values.std(ddof = 1)**2)**.5:.2f}")
# 			axhline(-depth, alpha = .5)			
# 			fid.write(f'\n{s},{ref},{d18values.mean():.2f},{(d18Osw_model_sigma**2 + d18values.std(ddof = 1)**2)**.5:.2f},{Tvalues.mean():.2f},{(Tsw_model_sigma**2 + Tvalues.std(ddof = 1)**2)**.5:.2f}')
# 		else:
		if GENERATE_PLOTS:
			axhspan(-zmin, -zmax, alpha = .25)

		fid.write(f'\n{s},{ref},{d18values.mean():.2f},{(d18Osw_model_sigma**2 + d18values.std(ddof = 1)**2)**.5:.2f},{Tvalues.mean():.2f},{(Tsw_model_sigma**2 + Tvalues.std(ddof = 1)**2)**.5:.2f}')
		d18values_500m = array(d18values_500m)
		Tvalues_500m = array(Tvalues_500m)
		fid_500m.write(f'\n{s},{ref},{d18values_500m.mean():.2f},{(d18Osw_model_sigma**2 + d18values_500m.std(ddof = 1)**2)**.5:.2f},{Tvalues_500m.mean():.2f},{(Tsw_model_sigma**2 + Tvalues_500m.std(ddof = 1)**2)**.5:.2f}')

		if GENERATE_PLOTS:
			title(f'{s}\n{ref}')
			savefig(f'plots/{s}_bk.pdf')
			close(fig)

if GENERATE_PLOTS:

	fig = figure(figsize = (3,2.5))
	subplots_adjust(.05, .2, .95, .95)

	R = d18Osw_residuals

	kw = dict(
		bins = linspace(-1,1,41),
		)
	hist(R, histtype = 'stepfilled', lw = 0, color = (.2,.6,1), **kw)
	hist(R, histtype = 'step', lw = 1, color = 'k', **kw)

	yticks([])
	xlabel('δ$^{18}O_{sw}$ residuals (‰)')
	gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.1f}$' if x else '$0$'))

	ymax = axis()[-1]
	axis([None, None, 0, ymax*1.7])

	rmse = (array(R)**2).mean()**.5

	res = linspace(0,1,1001)
	pres = [len([r for r in R if abs(r)<_res]) / len(R) for _res in tqdm(res)]
	f = interp1d(pres, res)
	cl95 = f(0.95)
	_, cl, _ = errorbar(0, ymax * 1.1, None, cl95, ecolor = 'k', elinewidth = 1, capthick = 1, capsize = 4, ls = 'None', marker = 'None')
	cl[0].set_marker([(1.3,0), (2,1), (0,0), (2,-1), (1.3,0)])
	cl[1].set_marker([(-1.3,0), (-2,1), (0,0), (-2,-1), (-1.3,0)])
	text(0, ymax * 1.15, f'± {cl95:.2f} (95 %)', ha = 'center', va = 'bottom', size = 10, linespacing = 1)
	fill_between([-cl95, cl95], [0]*2, [ymax * 1.1]*2, zorder = -100, color = [.9]*3, lw = 0)

	grid(False)
	text(0.02, 0.96, '(B) Variability from seawater estimates', size = 8.9, va = 'top', ha = 'left', transform = gca().transAxes)
	savefig('d18Osw_residuals.pdf')
	close(fig)

	# with open('../output/saved_values.py', 'a') as fid:
	# 	fid.write(f'dswRMSE = {rmse:.2f}\n')
	# 	fid.write(f'dswCL = {cl95:.2f}\n')

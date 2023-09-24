#! /usr/bin/env python3
'''
Assign calcification temperatures.
'''

from scipy.interpolate import interp1d
from csv import DictReader
import netCDF4 as nc
from pylab import *

d18Osw_model_sigma = 0.1
Tsw_model_sigma = 1.0

### SEAWATER d18O GRIDDED MODEL

ds = nc.Dataset('breitkreutz_2018/D18O_Breitkreuz_et_al_2018.nc')


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

with open('../1_compile_D47_data/sites.csv') as f:
	sites = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in sites:
	for k in r:
		if k not in ['Site', 'Ref']:
			r[k] = float(r[k])



print('Extracting seawater d18O for benthics...')

with (
	open('bottom_d18Osw.csv', 'w') as fid,
	):
	fid.write('Site,T,SE_T,d18Osw,SE_d18Osw')

	for s in sites:
		site, lon, lat, depth, ref = s['Site'], s['Lon'], s['Lat'], s['Depth'], s['Ref']

		print(f'\tProcessing {site} from {ref}...')

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
		
		X, Y, Tloc, M = depths[:], d18o[:,:,j,k], T[:,:,j,k], d18o_mask[0,:,j,k]
		X, Y, Tloc = X[M<1], Y[:,M<1], Tloc[:,M<1]

		maxdepth = X[-1]

		d18values, d18values_500m = [], []
		Tvalues, Tvalues_500m = [], []
		for y in Y:
			sca(ax1)
			plot(y, -X, 'b-', alpha = .1)
			f = interp1d(X,y)
			d18values += [f(depth)]

		for t in Tloc:
			sca(ax2)
			plot(t, -X, 'r-', alpha = .1)
			f = interp1d(X,t)
			Tvalues += [f(depth)]

		kw = dict(elinewidth = 1.5, alpha = 1, capsize = 5, marker = 'None', ls = 'None', capthick = 1.5)

		d18values = array(d18values)
		d18, sd18 = d18values.mean(), (d18Osw_model_sigma**2 + d18values.std(ddof = 1)**2)**.5
		sca(ax1)
		errorbar(d18, -depth, None, 1.96*sd18, ecolor = 'b', **kw)
		xlabel('d18Osw')
		ylabel('depth')
		text(.5, .97, f'{site}, {ref}', va = 'top', ha = 'center', transform = fig.transFigure)

		Tvalues = array(Tvalues)
		t, st = Tvalues.mean(), (Tsw_model_sigma**2 + Tvalues.std(ddof = 1)**2)**.5		
		sca(ax2)
		errorbar(t, -depth, None, 1.96*st, ecolor = 'r', **kw)
		xlabel('T')

		fid.write(f'\n{site},{t:.2f},{st:.2f},{d18:.2f},{sd18:.2f}')

		savefig(f'plots/benthics/{site}_bk.pdf')
		close(fig)
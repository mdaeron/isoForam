#! /usr/bin/env python3

import netCDF4 as nc
from csv import DictReader
from scipy.interpolate import interp1d
from pylab import *


with open('../1_compile_D47_data/sites.csv') as f:
	sites = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in sites:
	for k in r:
		if k not in ['Site', 'Ref']:
			r[k] = float(r[k])

seawater_chemistry = {}

try:
	for prop, netcdf in [
		('omega', 'GLODAPv2/GLODAPv2.2016b.OmegaC.nc'),
		('salinity', 'GLODAPv2/GLODAPv2.2016b.salinity.nc'),
		]:

		ds = nc.Dataset(netcdf)
	# 	print(ds)

		lats = ds['lat'][:]
		lons = ds['lon'][:]
		depths = ds['Depth'][:]

		var = netcdf.split('.')[-2]
		prp = ds[var][:,:,:]
		prp_error = ds[f'{var}_error'][:,:,:]

		prp = prp.filled(fill_value = nan)
		prp_error = prp_error.filled(fill_value = nan)

		max_depth = [[
			depths[~isnan(prp[:,i,j])]
			for j in range(len(lons))] for i in range(len(lats))]

		max_depth = array([[
			max_depth[i][j].max() if len(max_depth[i][j]) else nan
			for j in range(len(lons))] for i in range(len(lats))])

		seawater_chemistry[prop] = prp
		seawater_chemistry[f'{prop}_error'] = prp_error
		seawater_chemistry[f'{prop}_max_depth'] = max_depth

except FileNotFoundError:

	print("""
	ERROR: MISSING GLODAPv2 DATA SETS
	
	These large files (100 MB each) are not included here by default, but they can be
	obtained from <https://www.glodap.info/index.php/mapped-data-product>.
	
	The files to download are:
	
	- GLODAPv2.2016b.OmegaC.nc
	- GLODAPv2.2016b.salinity.nc

	""")
	exit()


# ds = nc.Dataset('GLODAPv2/GLODAPv2.2016b.salinity.nc')
# # print(ds)
# 
# # lats = ds['lat'][:]
# # lons = ds['lon'][:]
# # depths = ds['Depth'][:]
# 
# salinity = ds['salinity'][:,:,:]
# salinity_error = ds['salinity_error'][:,:,:]
# 
# salinity = salinity.filled(fill_value = nan)
# salinity_error = salinity_error.filled(fill_value = nan)
# 
# 
# salinity_max_depth = [[
# 	depths[~isnan(salinity[:,i,j])]
# 	for j in range(len(lons))] for i in range(len(lats))]
# 
# salinity_max_depth = array([[
# 	salinity_max_depth[i][j].max() if len(salinity_max_depth[i][j]) else nan
# 	for j in range(len(lons))] for i in range(len(lats))])
# 
# seawater_chemistry['salinity'] = salinity
# seawater_chemistry['salinity_error'] = salinity_error
# seawater_chemistry['salinity_max_depth'] = salinity_max_depth

glon, glat = meshgrid(lons, lats)
gx = cos(glon * pi / 180) * cos(glat * pi / 180)
gy = sin(glon * pi / 180) * cos(glat * pi / 180)
gz = sin(glat * pi / 180)

def get_seawater_chemistry(lat, lon, depth, prop):
	x = cos(lon * pi / 180) * cos(lat * pi / 180)
	y = sin(lon * pi / 180) * cos(lat * pi / 180)
	z = sin(lat * pi / 180)
	sqdistance = (gx-x)**2 + (gy-y)**2 + (gz-z)**2

	max_depth = seawater_chemistry[f'{prop}_max_depth']

	# penalize max_depth NaNs
	sqdistance += isnan(max_depth).astype(int) * 100
	
	if type(depth) == float:
		# penalize max_depth less than depth
		sqdistance += (max_depth < depth).astype(int) * 100

	j,k = [int(_) for _ in where(sqdistance == sqdistance.min())]
# 	distance = sqdistance[j,k]**.5 * 4e4 / 2 / pi

	if depth == 'bottom':
		i = int(where(depths == max_depth[j,k])[0])
		depth = depths[i]
		W = seawater_chemistry[prop][i,j,k]
		sW = seawater_chemistry[f'{prop}_error'][i,j,k]
	else:
		f = interp1d(depths, seawater_chemistry[prop][:,j,k])
		W = f(depth)
		sf = interp1d(depths, seawater_chemistry[f'{prop}_error'][:,j,k])
		sW = sf(depth)
	
	return W, sW, depth
	
	

with open('bottom_seawater_chemistry.csv', 'w') as fid:
	fid.write('Site,Lat,Lon,Depth,Omega_cc,SE_Omega_cc,Salinity,SE_Salinity')
	_sites = []
	for s in sites:
		if s['Site'] in _sites:
			continue
		site = s['Site']
		print(f'\tProcessing {site}')
		_sites += [site]
		W, sW, depth = get_seawater_chemistry(s['Lat'], s['Lon'], s['Depth'], 'omega')
		S, sS, depth = get_seawater_chemistry(s['Lat'], s['Lon'], s['Depth'], 'salinity')
		fid.write(f"\n{site},{s['Lat']:.2f},{s['Lon']:.2f},{s['Depth']:.0f},{W:.3f},{sW:.3f},{S:.3f},{sS:.3f}")

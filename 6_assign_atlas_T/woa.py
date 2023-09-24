#! /usr/bin/env python3

import netCDF4 as nc
from pylab import *
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d

import warnings
warnings.simplefilter('error', RuntimeWarning)

try:
	datasets = [nc.Dataset(f'../6_assign_atlas_T/woa/woa23_decav91C0_t{tindex:02.0f}_04.nc') for tindex in range(17)]
except FileNotFoundError:
	print("""
	ERROR: MISSING WOA DATA SETS
	
	These large files (2-3 GB each) are not included here by default, but they can be
	obtained from <https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/bin/woa23.pl>.
	
	The files to download are:
	
	- woa23_decav91C0_t04_04.nc
	- woa23_decav91C0_t00_04.nc
	- woa23_decav91C0_t01_04.nc
	- woa23_decav91C0_t02_04.nc
	- woa23_decav91C0_t03_04.nc
	- woa23_decav91C0_t05_04.nc
	- woa23_decav91C0_t06_04.nc
	- woa23_decav91C0_t07_04.nc
	- woa23_decav91C0_t08_04.nc
	- woa23_decav91C0_t09_04.nc
	- woa23_decav91C0_t10_04.nc
	- woa23_decav91C0_t11_04.nc
	- woa23_decav91C0_t12_04.nc
	- woa23_decav91C0_t13_04.nc
	- woa23_decav91C0_t14_04.nc
	- woa23_decav91C0_t15_04.nc
	- woa23_decav91C0_t16_04.nc
	""")
	exit()

def get_T_from_woa(lat, lon, depth, site, tindex = 'annual'):

	T, method = None, None
	
	zi = asarray(depth)
		
	if type(tindex) == str:
		tindex = {
			'annual': 0,
			'winter': 13,
			'spring': 14,
			'summer': 15,
			'fall': 16,
			}[tindex.lower()]
	
	ds = datasets[tindex]

	lats = ds['lat'][:]
	lons = ds['lon'][:]
	depths = ds['depth'][:]

	t_se = ds['t_se'][0,:,:,:]
	t_se = t_se.filled(fill_value = nan)

	t_an = ds['t_an'][0,:,:,:]
	t_an = t_an.filled(fill_value = nan)

	t_an_mask = isnan(t_an).astype(int)
		
	f_t_an = RegularGridInterpolator((depths, lats, lons), t_an)
	f_t_se = RegularGridInterpolator((depths, lats, lons), t_se)

	glon, glat = meshgrid(lons, lats)
	gx = cos(glon * pi / 180) * cos(glat * pi / 180)
	gy = sin(glon * pi / 180) * cos(glat * pi / 180)
	gz = sin(glat * pi / 180)

	Ti = array(f_t_an([[_, lat, lon] for _ in zi]))
	sTi = array(f_t_se([[_, lat, lon] for _ in zi]))
	zi, Ti, sTi = zi[~isnan(Ti)], Ti[~isnan(Ti)], sTi[~isnan(Ti)]
	
	return zi, Ti, sTi


def get_fuzzy_T_from_woa(lat, lon, depth, site, tindex = 'annual', create_plot = True, plotdir = 'plots'):

	T, method = None, None

	if type(tindex) == str:
		tindex = {
			'annual': 0,
			'winter': 13,
			'spring': 14,
			'summer': 15,
			'fall': 16,
			}[tindex.lower()]
	
	ds = datasets[tindex]

	lats = ds['lat'][:]
	lons = ds['lon'][:]
	depths = ds['depth'][:]

	t_se = ds['t_se'][0,:,:,:]
	t_se = t_se.filled(fill_value = nan)

	t_an = ds['t_an'][0,:,:,:]
	t_an = t_an.filled(fill_value = nan)

	t_an_mask = isnan(t_an).astype(int)
		
	f_t_an = RegularGridInterpolator((depths, lats, lons), t_an)
	f_t_se = RegularGridInterpolator((depths, lats, lons), t_se)

	_zi = linspace(0,5000,5001)

	glon, glat = meshgrid(lons, lats)
	gx = cos(glon * pi / 180) * cos(glat * pi / 180)
	gy = sin(glon * pi / 180) * cos(glat * pi / 180)
	gz = sin(glat * pi / 180)

	Ti = array(f_t_an([[_, lat, lon] for _ in _zi]))
	sTi = array(f_t_se([[_, lat, lon] for _ in _zi]))
	zi, Ti, sTi = _zi[~isnan(Ti)], Ti[~isnan(Ti)], sTi[~isnan(Ti)]
	
	if depth == 'bottom':
		method = 'local interpolation'
		depth = zi[-1]
		T = Ti[-1]
		sT = sTi[-1]
		
	else:
		if isnan(zi).all() or depth > zi.max():
			method = 'neareast lateral exploration'
			x = cos(lon * pi / 180) * cos(lat * pi / 180)
			y = sin(lon * pi / 180) * cos(lat * pi / 180)
			z = sin(lat * pi / 180)
			sqdistance = (gx-x)**2 + (gy-y)**2 + (gz-z)**2

			i = [i for i, _ in enumerate(depths) if _ >= depth][0]

			sqdistance += t_an_mask[i,:,:] * 10
			j,k = [int(_) for _ in where(sqdistance == sqdistance.min())]

			distance = sqdistance[j,k]**.5 * 4e4 / 2 / pi

			lzi = depths[:]
			lTi = t_an[:,j,k]
			slTi = t_se[:,j,k]
			lzi, lTi, slTi = lzi[~isnan(lTi)], lTi[~isnan(lTi)], slTi[~isnan(lTi)]
			f = interp1d(lzi, lTi)
			T = f(depth)
			sf = interp1d(lzi, slTi)
			sT = sf(depth)
		else:
			method = 'local interpolation'
			f = interp1d(zi, Ti)
			T = f(depth)
			sf = interp1d(zi, sTi)
			sT = sf(depth)

	if isnan(sT):
		sT = 0.5

	sT = max(sT, 0.1)

	if create_plot:
		fig = figure(figsize = (5,5))
		subplots_adjust(.2,.15,.95,.9)

		if not isnan(zi).all():
			plot(Ti, -zi, '-', color = [.8]*3, lw = 4, label = 'local interpolation')

		if method is not None:
			errorbar(T, -depth, None, 1.96*sT, marker = 'None', ls = 'None', ecolor = 'k', elinewidth = 1, capsize = 3, capthick = 1)
			text(T + 1.96 * sT, -depth, '    core-top depth', va = 'center', ha = 'left', size = 8)

		if method == 'neareast lateral exploration':
			plot(lTi, -lzi, '-', color = 'k', lw = .7, alpha = 0.5, label = f'nearest full-depth profile ({distance:.0f} km away)')

		legend()
		xlabel('T (°C)')
		ylabel('Depth (m)')
		text(0, 1,
			f'{site}\n',
			size = 10,
			ha = 'left',
			va = 'center',
			transform = gca().transAxes,
			)
		text(1, 1,
			f'{abs(lat):.2f} °{"S" if lat < 0 else "N"}, {abs(lon):.2f} °{"W" if lat < 0 else "E"}\n',
			size = 10,
			ha = 'right',
			va = 'center',
			transform = gca().transAxes,
			)
		savefig(f'{plotdir}/{site}.pdf')
		close(fig)
	
	return (T, sT, depth, method)
	
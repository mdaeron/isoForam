#! /usr/bin/env python3
'''
Assign calcification temperatures.
'''

from csv import DictReader
from pylab import *
from scipy.odr import Model, RealData, ODR
from woa import get_T_from_woa, get_fuzzy_T_from_woa
from tqdm import tqdm
import sys

style.use('../mydefault.mplstyle')

sys.path.append('../4_species_depth')
from species_depth import species_depth

sys.path.append('../3_species_lookup')
from species_lookup import species_lookup

DIVE_INTO_WOA = False
PROCESS_BOTTOM = True
DISCORDANT_CRITERION_FACTOR = 1.96

with open('../7_assign_Tcalcif/foram_D47_calibration_data.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in data:
	tobepopped = []
	for k in r:
		if k not in ['Sample', 'Species', 'Type', 'Site', 'Ref', 'Twoa23_vs_Tiso_species', 'Tiso_species_offset', 'Twoa23_500m_vs_Tiso_species_500m', 'Twoa23_1500m_vs_Tiso_species_1500m']:
			if r[k] == '':
				tobepopped.append(k)
			else:
				if k == 'Depth' and '-' in r[k]:
					r['Depth'] = tuple(float(_) for _ in r[k][1:-1].split('-'))
				else:
					r[k] = float(r[k])
	for k in tobepopped:
		r.pop(k)

with open('../1_compile_D47_data/sites.csv') as f:
	sites = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in sites:
	for k in r:
		if k not in ['Site', 'Ref']:
			r[k] = float(r[k])

if PROCESS_BOTTOM:
	with open('../6_assign_atlas_T/bottom_temperatures.csv', 'w') as fid:
		fid.write('Site,Lat,Lon,Depth,mean_annual_Tbottom,SE_mean_annual_Tbottom')
		for s in sites:
			depth, lat, lon, site = s['Depth'], s['Lat'], s['Lon'], s['Site']
			T, sT, depth, method = get_fuzzy_T_from_woa(lat, lon, depth, site, plotdir = 'bottom')
			fid.write(f'\n{site},{lat:.2f},{lon:.2f},{depth:.0f},{T:.2f},{sT:.2f}')
			if T is not None:
				print(f'{site:>24}: mean annual T = {T:>4.1f} ± {sT:.1f} C ({method})')
			else:
				print(f'{site} failed!')


if DIVE_INTO_WOA:
	for r in tqdm(data):
		if r['Type'] == 'benthic':
			continue
		_zmin = species_depth[species_lookup[r['Species']]]['zmin']
		_zmax = species_depth[species_lookup[r['Species']]]['zmax']

		for filename, zmin, zmax in[
			(r['Sample'], _zmin, _zmax),
			(r['Sample']+'_500m', 0, 500),
			(r['Sample']+'_1500m', 0, 1500),
			]:
			zi = linspace(zmin, zmax, int(zmax-zmin+1))
			with open(f'histograms/csv/{filename}.csv', 'w') as fid:
				Tcol = linspace(-5,35,401)[:-1] + 0.05
				fid.write('Month,' + ','.join([f'{t:.2f}' for t in Tcol]))
				for month in range(12):
					zi, Ti, sTi = get_T_from_woa(r['Lat'], r['Lon'], zi, r['Site'], tindex = month+1)
					if not (zi.min() == zmin and zi.max() == zmax):
						print(f'{filename}: WOA depth range is narrower than {zmin}-{zmax} m.')
					Ti = (Ti - 0.05).round(1) + 0.05
					counts = dict(zip(*unique(Ti, return_counts = True)))
					fid.write(f'\n{month+1},' + ','.join([f'{counts[t]}' if t in counts else '0' for t in Tcol]))



with (
	open('../6_assign_atlas_T/planktic_temperatures.csv', 'w') as fid,
	open('../6_assign_atlas_T/planktic_temperatures_500m.csv', 'w') as fid_500m,
	open('../6_assign_atlas_T/planktic_temperatures_1500m.csv', 'w') as fid_1500m,
	):
	fid.write('Sample,Species,Site,Lat,Lon,Depth,Twoa23,SD_Twoa23,Discordant')
	fid_500m.write('Sample,Species,Site,Lat,Lon,Depth,Twoa23,SD_Twoa23,Discordant')
	fid_1500m.write('Sample,Species,Site,Lat,Lon,Depth,Twoa23,SD_Twoa23,Discordant')
	_collected_ = []
	_collected_500m_ = []
	_collected_1500m_ = []
	for r in tqdm(data):
		if r['Type'] == 'planktic':
			_zmin = species_depth[r['Species']]['zmin']
			_zmax = species_depth[r['Species']]['zmax']
			for filename, zmin, zmax, _fid_, _coll_ in[
				(r['Sample'], _zmin, _zmax, fid, _collected_),
				(r['Sample']+'_500m', 0, 500, fid_500m, _collected_500m_),
				(r['Sample']+'_1500m', 0, 1500, fid_1500m, _collected_1500m_),
				]:
				with open(f'histograms/csv/{filename}.csv') as gid:
					rdata = list(DictReader(gid))
		
					fig = figure(figsize = (3.7, 2))
					subplots_adjust(.03, .25, .97, .95)

					Tcol = [k for k in rdata[0]][1:]

					## Compute tmin, tmax with fine resolution
					dT = 0.1
					bins = linspace(-5,35,int(40/dT+1))
					temps = array([float(t) for t in Tcol])
					n_cumul = bins[1:]*0

					for m in rdata:
						x = array([float(t) for t in Tcol])
						weights = array([float(m[t]) for t in Tcol])
						c = histogram(x, bins = bins, weights = weights)[0]
						n_cumul += c


					t = (bins[1:] + bins[:-1])/2
					try:
						tmin = t[n_cumul>0].min()
						tmax = t[n_cumul>0].max()
					except ValueError:
# 							print(r['Sample'])
						close(fig)
						continue

					## Compute histogram with coarse resolution
					dT = 0.5
					bins = linspace(-5,35,int(40/dT+1))
					temps = array([float(t) for t in Tcol])
					n_cumul = bins[1:]*0

					for m in rdata:
						x = array([float(t) for t in Tcol])
						weights = array([float(m[t]) for t in Tcol])
						c = histogram(x, bins = bins, weights = weights)[0]
						n_cumul += c

					t = (bins[1:] + bins[:-1])/2
					counts = hist(t, bins = bins, weights = n_cumul, histtype = 'stepfilled', lw = 0, alpha = .25, color = (0,.5,1))[0]
					
					_coll_.append(dict(
						Sample = r['Sample'],
						Species = r['Species'].split(' ')[0][0] + '. ' + ' '.join(r['Species'].split(' ')[1:]),
						Site = r['Site'],
						))
					_coll_[-1]['t'] = t*1
					_coll_[-1]['w'] = n_cumul / n_cumul.sum()
					_coll_[-1]['bins'] = bins
					
					hist(t, bins = bins, weights = n_cumul, histtype = 'step', lw = 1, alpha = 1, color = [.5]*3)
					
					score = 0 if (
						tmin > (r['Tiso_species'] + DISCORDANT_CRITERION_FACTOR * r['SE_Tiso_species'])
						or tmax < (r['Tiso_species'] - DISCORDANT_CRITERION_FACTOR * r['SE_Tiso_species'])
						) else 1

					_coll_[-1]['discordant'] = 1-score

					tavg = (t*counts).sum() / counts.sum()
					tsd = ((counts * (t-tavg)**2).sum() / (counts.sum() - 1))**.5

					_fid_.write(f"\n{r['Sample']},{r['Species']},{r['Site']},{r['Lat']:.2f},{r['Lon']:.2f},({zmin:.0f}-{zmax:.0f}),{tavg:.2f},{tsd:.2f},{'discordant' if score < 0.05 else ''}")

# 					if score < 0.05:
# 
# 						errorbar(tavg, n_cumul.max()/4, None, 1.96 * tsd,
# 							ecolor = 'b',
# 							elinewidth = 1,
# 							capsize = 3,
# 							capthick = 1,
# 							ls = 'None',
# 							marker = 'o',
# 							mfc = 'w',
# 							mec = 'b',
# 							mew = 1,
# 							)

					_coll_[-1]['Tiso'] = r['Tiso_species']
					_coll_[-1]['SE_Tiso'] = r['SE_Tiso_species']

					errorbar(r['Tiso_species'], n_cumul.max()/2, None, 1.96 * r['SE_Tiso_species'],
						ecolor = (.9,0,0) if score < 0.05 else 'k',
						elinewidth = 1,
						capsize = 3,
						capthick = 1,
						ls = 'None',
						marker = 'D',
						mfc = 'w' if score < 0.05 else 'k',
						mec = (.9,0,0) if score < 0.05 else 'k',
						mew = 1,
						ms = 4,
						)
					text(
						r['Tiso_species'],
						n_cumul.max()/2,
						'$^{18}O$ estimate\nof calcification T\n',
						va = 'bottom',
						ha = 'center',
						color = (.85,0,0) if score < 0.05 else 'k',
						size = 8.2)
			
					xmin = bins[:-1][n_cumul>0][0] - dT
					xmax = bins[1:][n_cumul>0][-1] + dT

					xtop = median([t for t,n in zip(bins[:-1], n_cumul) for _ in range(int(n))]) + dT/2

					text(
						xtop,
						0,
						'Atlas T\n',
						va = 'bottom',
						ha = 'center',
						size = 8.2)
			
					xmin = min(xmin, r['Tiso_species'] - 2.2 * r['SE_Tiso_species'])
					xmax = max(xmax, r['Tiso_species'] + 2.2 * r['SE_Tiso_species'])
					axis([xmin, xmax, None, n_cumul.max()*1.3])
					text(.5, 0.95, f"{r['Sample']} – {r['Species']}", transform = gca().transAxes, va = 'top', ha = 'center', size = 8)
					yticks([])
					xlabel('T (°C)')
					savefig(f'histograms/pdf/{filename}.pdf')
					close(fig)

for G, depth in (
	(_collected_, ''),
	(_collected_500m_, '_500m'),
	(_collected_1500m_, '_1500m'),
	):

	G = sorted(G, key = lambda _: _['Tiso'])
	
	fig = figure(figsize = (6,9))
	subplots_adjust(.03, .06, .97, .98, 0, 0)

	for k,g in enumerate(G):
		ax = subplot(len(G), 1, k+1)
		g['ax'] = ax
		
		ax.hist(
			g['t'],
			bins = g['bins'],
			weights = g['w'],
			histtype = 'stepfilled',
			lw = 0.5,
			alpha = .4,
			color = (0,.5,1),
			)
		ax.errorbar(
			g['Tiso'],
			g['w'].max()/2,
			None,
			1.96*g['SE_Tiso'],
			ecolor = (.9,0,0) if g['discordant'] else 'k',
			elinewidth = 1,
			capsize = 2,
			capthick = 1,
			ls = 'None',
			marker = 'D',
			mfc = 'w' if g['discordant'] else 'k',
			mec = (.9,0,0) if g['discordant'] else 'k',
			mew = 0.75,
			ms = 3,
			)
		ax.text(
			-2.5 if k >= 20 else 30.5,
			g['w'].max()*0.51,
			f"{g['Species']} ({g['Site']})",
			ha = 'left' if k >= 20 else 'right',
			va = 'center',
			size = 6,
			)
		yticks([])
		ax.tick_params(length = 0, labelbottom = False)		
		axis([-3,31,None, g['w'].max()*1.2])
		ax.spines['top'].set_linewidth(.5)
		ax.spines['top'].set_color([.5]*3)
		ax.spines['bottom'].set_linewidth(.5)
		ax.spines['bottom'].set_color([.5]*3)
	

	G[0]['ax'].spines['top'].set_linewidth(1)
	G[0]['ax'].spines['top'].set_color('k')
	ax.spines['bottom'].set_linewidth(1)
	ax.spines['bottom'].set_color('k')
	ax.tick_params(length = 3, labelbottom = True)
	xlabel('T (°C)')
	fig.savefig(f'planktic_T18{depth}.pdf')
	close(fig)
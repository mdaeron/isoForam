#! /usr/bin/env python3
'''
Assign calcification temperatures.
'''

from csv import DictReader
from pylab import *
from scipy.odr import Model, RealData, ODR
from tqdm import tqdm
from multiprocessing import Pool, get_context

import sys

sys.path.append('../4_species_depth')
from species_depth import species_depth

sys.path.append('../3_species_lookup')
from species_lookup import species_lookup

sys.path.append('../6_assign_atlas_T')
from woa import get_T_from_woa, get_fuzzy_T_from_woa

DIVE_INTO_WOA = False

with open('malevich_2019_s1.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]
	
with open('malevich_data_out.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in data:
	tobepopped = []
	for k in r:
		if k not in ['Sample', 'Species', 'Type', 'Site', 'Ref', 'Twoa23_vs_Tiso_species', 'Tiso_species_offset', 'Twoa23_500m_vs_Tiso_species_500m']:
			if r[k] == '':
				tobepopped.append(k)
			else:
				if k == 'Depth' and '-' in r[k]:
					r['Depth'] = tuple(float(_) for _ in r[k][1:-1].split('-'))
				else:
					r[k] = float(r[k])
	for k in tobepopped:
		r.pop(k)

if DIVE_INTO_WOA:
	
	def poolfun(r):

		print(r['Sample'])

		for filename, zmin, zmax in[
			(r['Sample'], r['zmin'], r['zmax']),
			(r['Sample']+'_500m', 0, 500),
			]:
			zi = linspace(zmin, zmax, int(zmax-zmin+1))
			with open(f'histograms/csv/{filename}.csv', 'w') as fid:
				Tcol = linspace(-5,35,401)[:-1] + 0.05
				fid.write('Month,' + ','.join([f'{t:.2f}' for t in Tcol]))
				for month in range(12):
					zi, Ti, sTi = get_T_from_woa(r['Lat'], r['Lon'], zi, r['corename'], tindex = month+1)
# 					if not (zi.min() == zmin and zi.max() == zmax):
# 						print(f'{filename}: WOA depth range is narrower than {zmin}-{zmax} m.')
					Ti = (Ti - 0.05).round(1) + 0.05
					counts = dict(zip(*unique(Ti, return_counts = True)))
					fid.write(f'\n{month+1},' + ','.join([f'{counts[t]}' if t in counts else '0' for t in Tcol]))

	if __name__ == '__main__':
		with Pool() as p:
			p.map(poolfun, data)

if __name__ == '__main__':
	with (
		open('planktic_temperatures.csv', 'w') as fid,
		open('planktic_temperatures_500m.csv', 'w') as fid_500m,
		):
		fid.write('Sample,Species,Site,Lat,Lon,Depth,Twoa23,SD_Twoa23,Discordant')
		fid_500m.write('Sample,Species,Site,Lat,Lon,Depth,Twoa23,SD_Twoa23,Discordant')
		for r in tqdm(data):
			if r['Type'] == 'planktic':
				_zmin = species_depth[r['Species']]['zmin']
				_zmax = species_depth[r['Species']]['zmax']
				for filename, zmin, zmax, _fid_ in[
					(r['Sample'], _zmin, _zmax, fid),
					(r['Sample']+'_500m', 0, 500, fid_500m),
					]:
					with open(f'histograms/csv/{filename}.csv') as gid:
						rdata = list(DictReader(gid))
		
						fig = figure()
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
						counts = hist(t, bins = bins, weights = n_cumul, histtype = 'step', lw = 2, alpha = .5, color = 'r')[0]
						
						score = 0 if (
							tmin > (r['Tiso_species'] + 1.96 * r['SE_Tiso_species'])
							or tmax < (r['Tiso_species'] - 1.96 * r['SE_Tiso_species'])
							) else 1

						tavg = (t*counts).sum() / counts.sum()
						tsd = ((counts * (t-tavg)**2).sum() / (counts.sum() - 1))**.5

						_fid_.write(f"\n{r['Sample']},{r['Species']},{r['Site']},{r['Lat']:.2f},{r['Lon']:.2f},({zmin:.0f}-{zmax:.0f}),{tavg:.2f},{tsd:.2f},{'discordant' if score < 0.05 else ''}")

# 						if score < 0.05:
# 
# 							errorbar(tavg, n_cumul.max()/4, None, 1.96 * tsd,
# 								ecolor = 'r',
# 								elinewidth = 1,
# 								capsize = 3,
# 								capthick = 1,
# 								ls = 'None',
# 								marker = 'o',
# 								mfc = 'w',
# 								mec = 'r',
# 								mew = 1,
# 								)

						errorbar(r['Tiso_species'], n_cumul.max()/4, None, 1.96 * r['SE_Tiso_species'],
							ecolor = 'b' if score < 0.05 else 'k',
							elinewidth = 1,
							capsize = 3,
							capthick = 1,
							ls = 'None',
							marker = 'o',
							mfc = 'w',
							mec = 'b' if score < 0.05 else 'k',
							mew = 1,
							)
			
						xmin = bins[:-1][n_cumul>0][0] - dT
						xmax = bins[1:][n_cumul>0][-1] + dT
			
						xmin = min(xmin, r['Tiso_species'] - 2.2 * r['SE_Tiso_species'])
						xmax = max(xmax, r['Tiso_species'] + 2.2 * r['SE_Tiso_species'])
						axis([xmin, xmax, None, None])
						title(f"{r['Sample']} – {r['Species']}")
						savefig(f'histograms/pdf/{filename}.pdf')
						close(fig)

#! /usr/bin/env python3

from pylab import *
from csv import DictReader
from scipy.interpolate import interp1d

sys.path.append('../../3_species_lookup')
from species_lookup import species_lookup

with open('meilland_2019_table_s2.csv') as f:
	data = list(DictReader(f))

for r in data:
	for k in r:
		if k == 'Station':
			continue
		r[k] = float(r[k])

depths = sorted(set([r['Depth bottom'] for r in data] + [r['Depth top'] for r in data]))

species = sorted([
	k for k in data[0]
	if k not in [
		'Station',
		'Depth top',
		'Depth bottom',
		'benthic',
		'undeterminate',
		'Total shells',
		]
	])

zi = linspace(0, max(depths), int(max(depths)+1))

with open('meilland_2019_depths.csv', 'w') as fid:
	fid.write('Species,zmin,zmax,Ref')
	for s in species:
		print(f'Processing {s}')
		fig = figure()
	
		y = zi*0
		for r in data:
			zmin = r['Depth top']
			zmax = r['Depth bottom']
			y[(zmin <= zi) & (zi <= zmax)] += r[s] / (zmax - zmin)

		subplot(121)
		plot(y, zi, 'r-')
		xticks([0])
		axvline(0, lw = 0.5, color = 'k', alpha = .25)
		xlabel('Probability distribution')
		ylabel('Depth')
		title(s)
		axis([None, None, zi[-1], zi[0]])

		x = cumsum(y)
		x -= x[0]
		x /= x[-1]

		f = interp1d(x,zi)
		intervals_95cl = []
		for x0,z0 in zip(x, zi):
			if x0 == 0:
				continue
			if x0 > 0.05:
				break
			z1 = float(f(x0 + 0.95))
			intervals_95cl.append([z0, z1, z1-z0])

		intervals_95cl = sorted(intervals_95cl, key = lambda _: _[-1])
		z_95cl = intervals_95cl[0]
	
		zmin = (z_95cl[0]//10)*10
		zmax = (z_95cl[1]//10)*10 + 10

		if 'juvenile' not in s:
			fid.write(f'\n{species_lookup[s]},{zmin:.0f},{zmax:.0f},Meilland et al. (2019)')

		axhspan(zmin, zmax, color = 'r', alpha = 0.25, lw = 0)

		subplot(122)
		plot(x, zi, 'b-')
		axhspan(*z_95cl[:2], color = 'b', alpha = 0.25, lw = 0)
		xticks([0,1])
		xlabel('Cumulative probability distribution')
		title(f'{s}: {zmin:.0f}–{zmax:.0f} m')
		axis([None, None, zi[-1], zi[0]])


		savefig(f'{s}.pdf')
	
		close(fig)


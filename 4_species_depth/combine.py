#! /usr/bin/env python3
from csv import DictReader

data = []

for file in [
# 	'Pflaumann_1999/pflaumann_1999.csv',
	'Greco_2019/greco_2019.csv',
	'Rebotim_2017/rebotim_2017.csv',
	'Meilland_2019/meilland_2019_depths.csv',
	]:

	with open(file) as f:
		data += list(DictReader(f))
	
data = sorted(data, key = lambda _: _['Species'] + _['Ref'])

species = sorted({r['Species'] for r in data})

with open('depth_data.csv', 'w') as fid:
	fid.write('Species,zmin,zmax,Ref')
	for s in species:
		G = [r for r in data if r['Species'] == s]
		if len(G) < 2:
			r = G[0]
			fid.write(f'\n{r["Species"]},{r["zmin"]},{r["zmax"]},{r["Ref"]}')
		else:
			zmin = min([int(r['zmin']) for r in G])
			zmax = max([int(r['zmax']) for r in G])
			refs = '; '.join([r['Ref'] for r in G])
			fid.write(f'\n{s},{zmin},{zmax},{refs}')

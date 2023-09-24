#! /usr/bin/env python3

from csv import DictReader
from pylab import *

style.use('../mydefault.mplstyle')

with open('../7_assign_Tcalcif/foram_D47_calibration_data.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in data:
	tobepopped = []
	for k in r:
		if k not in ['Sample', 'Species', 'Type', 'Site', 'Ref', 'Twoa23_vs_Tiso_species', 'Twoa23_500m_vs_Tiso_species_500m', 'Twoa23_1500m_vs_Tiso_species_1500m', 'Tiso_species_offset']:
			if r[k] == '':
				tobepopped.append(k)
			else:
				if k == 'Depth' and '-' in r[k]:
					r['Depth'] = tuple(float(_) for _ in r[k][1:-1].split('-'))
				else:
					r[k] = float(r[k])
	for k in tobepopped:
		r.pop(k)

data = [_ for _ in data if _['Type'] == 'planktic']
sites_with_discordants = {_['Site'] for _ in data if _['Twoa23_vs_Tiso_species'] == 'discordant'}
sites_with_concordants = {_['Site'] for _ in data if _['Twoa23_vs_Tiso_species'] == ''}
sites = sites_with_discordants.intersection(sites_with_concordants)
sites = sorted(sites)

fig = figure(figsize = (6,3))
subplots_adjust(.05, .17, .67, .95)
labels = []
for k, site in enumerate(sites):
	dX = 0
	G = [_ for _ in data if _['Site'] == site]
	for r in G:
		marker = 's' if r['Twoa23_vs_Tiso_species'] == 'discordant' else 'x'
		X = r['d13C_VPDB']
		Y = -k
		eY = 1.96*r['SE_d13C_VPDB']
		plot(X, Y, marker,
			mec = 'r' if r['Twoa23_vs_Tiso_species'] == 'discordant' else 'k',
			mfc = 'None',
			ms = 6 if r['Twoa23_vs_Tiso_species'] == 'discordant' else 5,
			zorder = 100 if r['Twoa23_vs_Tiso_species'] == 'discordant' else 200,
			mew = 1 if r['Twoa23_vs_Tiso_species'] == 'discordant' else 0.75,
			label = None if marker in labels else (
				'discordant'
				if r['Twoa23_vs_Tiso_species'] == 'discordant' else
				'concordant'
				)
			)
		labels += marker

basins = {
	'CD107 A ML 5A': 'N. Atlantic',
	'MD08-3179': 'N. Atlantic',
	'SO164-25-3': 'N. Atlantic',
	'SO213-84-2': 'Pacific',
	'SO225-53-1': 'Pacific',
	'WIND 33B': 'Indian Ocean',
	}
sitelabels = [f'{site} ({basins[site]})' for site in sites]

legend(fontsize = 9, loc = 'lower left')
gca().yaxis.tick_right()
yticks([-k for k in range(len(sites))])
xlabel('$δ^{13}C_{VPDB}$ (‰)')
gca().set_yticklabels(sitelabels)
axis([None, None, -len(sites), 1])
savefig('check_d13C.pdf')
close(fig)
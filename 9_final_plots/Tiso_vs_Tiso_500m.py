#! /usr/bin/env python3

from csv import DictReader
from pylab import *
from matplotlib import ticker

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

fig = figure(figsize = (6,3))
subplots_adjust(.15,.2,.95,.9)
ax = subplot(111)

X, sX, Y, sY = array([
	[_['Tiso_species'], _['SE_Tiso_species'], _['Tiso_species_500m'], _['SE_Tiso_species_500m']]
	for _ in data
	if _['Type'] == 'planktic'
	]).T

# plot(X, Y-X, 'wo', mec = 'k', mew = .8, ms = 4)
ax.errorbar(X, Y-X, 1.96*sX,
	ls = 'None',
	marker = 'o',
	mec = 'k',
	mew = 0.8,
	ms = 3,
	mfc = [.5]*3,
	ecolor = (0,0,0,.3),
	elinewidth = 0.8,
	capthick = 0.8,
	capsize = 1,
	zorder = 100,
	)

# axlims = ax.axis()
# xmin, xmax = min(axlims), max(axlims)
ax.axhline(0, color = 'k', lw = 0.75, alpha = .5)
# ax.plot([xmin, xmax], [xmin, xmax], 'k-', lw = 0.75, dashes = (8,2,2,2))
ax.set_xlabel('$T_{18}$ (°C) estimated over assumed living depths for species')
ax.set_ylabel('Effect of estimating $T_{18}$\nover 0–500 m depth instead\n')
# ax.axis([xmin, xmax, xmin, xmax])

fig.savefig('plots/Tiso_vs_Tiso_500m.pdf')
close(fig)
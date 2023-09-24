#! /usr/bin/env python3

from csv import DictReader
from pylab import *

style.use('../mydefault.mplstyle')

with open('../7_assign_Tcalcif/foram_D47_calibration_data.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in data:
	tobepopped = []
	for k in r:
		if k not in ['Sample', 'Species', 'Type', 'Site', 'Ref', 'Twoa23_vs_Tiso_species', 'Twoa23_500m_vs_Tiso_species_500m', 'Tiso_species_offset']:
			if r[k] == '':
				tobepopped.append(k)
			else:
				if k == 'Depth' and '-' in r[k]:
					r['Depth'] = tuple(float(_) for _ in r[k][1:-1].split('-'))
				else:
					r[k] = float(r[k])
	for k in tobepopped:
		r.pop(k)


fig = figure(figsize = (6,4))
subplots_adjust(.11, .17, .61, .92)

discordant_species = sorted({r['Species'] for r in data if r['Twoa23_vs_Tiso_species'] == 'discordant'})

kweb = dict(
	elinewidth = 1.,
	capsize = 2,
	ls = 'None',
	marker = 'None',
	)
kweb['capthick'] = kweb['elinewidth']


kwplot = dict(
	ls = 'None',
	ms = 5,
	mec = 'k',
	mew = 0.9,
	mfc = 'w',
	zorder = 1000,
	)

### DISCORDANT
G = [r for r in data if r['Type'] == 'planktic' and r['Twoa23_vs_Tiso_species'] == 'discordant']
X = [r['Tbottom_woa23'] for r in G]
Y = [r['Twoa23'] for r in G]
Z = [r['Tiso_species'] for r in G]
eX = [1.96 * r['SE_Tbottom_woa23'] for r in G]
eY = [1.96 * r['SE_Twoa23'] for r in G]
eZ = [1.96 * r['SE_Tiso_species'] for r in G]

alpha = 1
kweb['ecolor'] = [alpha*_ + (1-alpha) for _ in [1,0,0]]
ebs = errorbar(Z, Y, eY, label = 'Surface T from WOA', **kweb)
# plot([], [], '-', color = kweb['ecolor'], lw = kweb['elinewidth'], label = 'Surface T from WOA')

alpha = 1
kweb['ecolor'] = [alpha*_ + (1-alpha) for _ in [0,0,0]]
ebi = errorbar(Z, Z, eZ, label = 'Calcification T from $^{18}$O', **kweb)
# plot([], [], '-', color = kweb['ecolor'], lw = kweb['elinewidth'], label = 'Calcification T from $^{18}$O')

alpha = 1
kweb['ecolor'] = [alpha*_ + (1-alpha) for _ in [0,0,1]]
ebb = errorbar(Z, X, eX, label = 'Bottom T from WOA', **kweb)
# plot([], [], '-', color = kweb['ecolor'], lw = kweb['elinewidth'], label = 'Bottom T from WOA')


for (sp, marker) in zip(discordant_species, ['d', (4,0,0), (4,0,45), (3,0,270), (3,0,90), (3,0,180), (3,0,0)]):
	H = [r for r in G if r['Species'] == sp]
	X = [r['Tbottom_woa23'] for r in H]
	Y = [r['Twoa23'] for r in H]
	Z = [r['Tiso_species'] for r in H]
	eX = [1.96 * r['SE_Tbottom_woa23'] for r in H]
	eY = [1.96 * r['SE_Twoa23'] for r in H]
	eZ = [1.96 * r['SE_Tiso_species'] for r in H]

	kwplot['marker'] = marker
	kwplot['label'] = sp
	plot(Z, Z, **kwplot)

plot([],[],ls = 'None', marker = 'None', label = '\n')

# G = [r for r in data if r['Type'] == 'planktic' and r['Twoa23_vs_Tiso_species'] == '']
# # X = [r['Tbottom_woa23'] for r in G]
# Y = [r['Twoa23'] for r in G]
# Z = [r['Tiso_species'] for r in G]
# # eX = [1.96 * r['SE_Tbottom_woa23'] for r in G]
# eY = [1.96 * r['SE_Twoa23'] for r in G]
# eZ = [1.96 * r['SE_Tiso_species'] for r in G]
# 
# kwplot['marker'] = '+'
# kwplot['ms'] = 3
# kwplot['alpha'] = 0.3
# kwplot['label'] = 'Non-discordant'
# plot(Z, Y, **kwplot)

x1, x2 = [-7, 32]
plot([x1, x2], [x1, x2], 'k-', lw = .7, alpha = .5, dashes = (6,2,2,2))
text(
	1-0.03, 0.03,
	'only discordant planktic data shown here',
	transform = gca().transAxes,
	va = 'bottom',
	ha = 'right',
	alpha = .5,
	style = 'italic',
	size = 9,
	)

axis([x1, x2]*2)

legend(loc = 6, labelspacing = 1, fontsize = 8, bbox_to_anchor = (1.05, 0.5), frameon = False)

xlabel('Isotopic estimate of calcification temperature (°C)\nat assumed living depths')
ylabel('Other estimates of temperatures (°C)')
savefig(f'plots/Tiso_vs_Tbottom_for_discordant_planktics.pdf')
close(fig)


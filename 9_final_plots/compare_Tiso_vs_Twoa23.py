#! /usr/bin/env python3

from csv import DictReader
from york import YorkReg
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



for suffix in ['', '_500m']:
	fig = figure(figsize = (5,5))
	subplots_adjust(.15, .15, .95, .95)

	discordant_species = sorted({r['Species'] for r in data if r['Twoa23_vs_Tiso_species'] == 'discordant'})

	kweb = dict(
		ecolor = (1,.65,1) if suffix else (1,.65,.65),
		elinewidth = 1,
		capsize = 2,
		capthick = 1,
		ls = 'None',
		marker = 'None',
		)


	kwplot = dict(
		alpha = .5,
		ls = 'None',
		marker = 'o',
		ms = 5,
		mec = 'k',
		mew = 1,
		mfc = 'w',
		zorder = 1000,
		)

	# G = [r for r in data if r['Type'] == 'planktic' and r['Twoa23_vs_Tiso_species'] == '']
	# X = [r['Tiso_species'] for r in G]
	# Y = [r['Twoa23'] for r in G]
	# eX = [1.96 * r['SE_Tiso_species'] for r in G]
	# eY = [1.96 * r['SE_Twoa23'] for r in G]
	# 
	# # errorbar(X, Y, eY, eX, **kweb)
	# plot(X, Y, label = '\nconcordant\n', **kwplot)
	# 


	### NON-DISCORDANT
	if suffix:
		G = [r for r in data if r['Type'] == 'planktic' and r['Twoa23_vs_Tiso_species'] == '' or r['Twoa23_500m_vs_Tiso_species_500m'] == '']
	else:
		G = [r for r in data if r['Type'] == 'planktic' and r['Twoa23_vs_Tiso_species'] == '']

	for (sp, marker) in zip(discordant_species, 'dDs<>^v'):
		H = [r for r in G if r['Species'] == sp]
		X = [r['Tiso_species'+suffix] for r in H]
		Y = [r['Twoa23'+suffix] for r in H]
		eX = [1.96 * r['SE_Tiso_species'+suffix] for r in H]
		eY = [1.96 * r['SE_Twoa23'+suffix] for r in H]
		kwplot['marker'] = marker
		plot(X, Y, label = sp.split(' ')[0][0] + '. ' + ' '.join(sp.split(' ')[1:]), **kwplot)

	kwplot['mec'] = 'm' if suffix else 'r'
	kwplot['alpha'] = 1

	### DISCORDANT
	if suffix:
		G = [r for r in data if r['Type'] == 'planktic' and r['Twoa23_vs_Tiso_species'] == 'discordant' and r['Twoa23_500m_vs_Tiso_species_500m'] == 'discordant']
	else:
		G = [r for r in data if r['Type'] == 'planktic' and r['Twoa23_vs_Tiso_species'] == 'discordant']
	for (sp, marker) in zip(discordant_species, 'dDs<>^v'):
		H = [r for r in G if r['Species'] == sp]
		X = [r['Tiso_species'+suffix] for r in H]
		Y = [r['Twoa23'+suffix] for r in H]
		eX = [1.96 * r['SE_Tiso_species'+suffix] for r in H]
		eY = [1.96 * r['SE_Twoa23'+suffix] for r in H]

		kwplot['marker'] = marker

		errorbar(X, Y, eY, eX, **kweb)
		baz, = plot(X, Y, **kwplot)


	### NON-DISCORDANT
	G = [r for r in data if r['Type'] == 'planktic' and r['Species'] not in discordant_species]

	kwplot['mec'] = 'k'
	kwplot['alpha'] = .5
	kwplot['marker'] = 'o'

	X = [r['Tiso_species'+suffix] for r in G]
	Y = [r['Twoa23'+suffix] for r in G]
	eX = [1.96 * r['SE_Tiso_species'+suffix] for r in G]
	eY = [1.96 * r['SE_Twoa23'+suffix] for r in G]

	plot(X, Y, label = 'other species', **kwplot)

	x1, x2 = [-3, 34]
	plot([x1, x2], [x1, x2], 'k-', lw = .7, alpha = .5)

	axis([x1, x2]*2)
	legend(loc = 4, labelspacing = 0.25)
	xlabel('Isotopic estimate of calcification temperature (°C)\nat ' + ('0-500 m' if suffix else 'assumed living depths'))
	ylabel('Annual range of temperatures (°C)\nat ' + ('0-500 m' if suffix else 'assumed living depths'))
	savefig(f'plots/Tiso_vs_Twoa23_planktic{suffix}.pdf')
	close(fig)







	fig = figure(figsize = (6,3))
	subplots_adjust(.12, .15, .975, .95, .25, .25)

	for ax, T, txt in [
		(subplot(121), 'Tiso_species', 'With isotopic temperature estimates\nfor discordant data points'),
		(subplot(122), 'Twoa23', 'With atlas temperature estimates\nfor discordant data points'),
		]:

		ax.text(0.05, 0.95, txt, ha = 'left', va = 'top', transform = ax.transAxes, size = 8)

		G = [r for r in data if r['Type'] == 'planktic' and r['Twoa23_vs_Tiso_species'] == '']
		X = [(r['Tiso_species'] + 273.15)**-2 for r in G]
		Y = [r['D47_ICDES'] for r in G]
		sX = [2 * r['SE_Tiso_species'] * (r['Tiso_species'] + 273.15)**-3 for r in G]
		sY = [r['SE_D47_ICDES'] for r in G]

		A, B, CM, rchisq = YorkReg(X, Y, sX, sY)

		xmin, xmax = (273.15+32)**-2, (273.15-3)**-2
		xi = linspace(xmin, xmax)
		yi = A * xi + B
		syi = (CM[0,0]*xi**2 + 2*CM[0,1]*xi + CM[1,1])**.5

		kw_yorkreg = dict(
			lw = 0.8,
			color = 'k',
			zorder = 15,
			dashes = (6,3,2,3),
			)

		plot(xi[:-2], yi[:-2]+1.96*syi[:-2], '-', **kw_yorkreg)
		plot(xi[2:], yi[2:]-1.96*syi[2:], '-', **kw_yorkreg)

		kw = dict(
			marker = 'o',
			ls = 'None',
			mec = 'k',
			mfc = 'w',
			ms = 4,
			mew = 0.7,
			)

		ax.plot(X, Y, **kw, label = 'concordant data')

		G = [r for r in data if r['Type'] == 'planktic' and r['Twoa23_vs_Tiso_species'] == 'discordant']
		X = [(r[T] + 273.15)**-2 for r in G]
		Y = [r['D47_ICDES'] for r in G]

		kw['mfc'] = (.5,.5,.5)
		kw['zorder'] = 100

		ax.plot(X, Y, **kw, label = 'discordant data')

		sca(ax)
		Ti = [30,20,10,0]
		xticks([(273.15+_)**-2 for _ in Ti])
		ax.set_xticklabels([f"${t}\\,$°C" for t in Ti])
		if T == 'Tiso_species':
			ylabel(f'Δ$_{{47}}$ (I-CDES, ‰)')


	legend(loc = 4, markerfirst = False, labelspacing = 0.5, frameon = False, bbox_to_anchor = (0.98,0), borderpad = 0, handlelength = 0.5)

	savefig('plots/Regression_with_discordant_data.pdf')
	close(fig)
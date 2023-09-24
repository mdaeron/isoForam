#! /usr/bin/env python3
'''
Assign calcification temperatures.
'''

from csv import DictReader
from pylab import *
from tqdm import tqdm
import matplotlib.ticker as ticker

style.use('../mydefault.mplstyle')

with open('../7_assign_Tcalcif/foram_D47_calibration_data.csv') as f:
	D47_data = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in D47_data:
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
D47_data = [r for r in D47_data if r['Type'] == 'planktic']

with open('malevich_data_out.csv') as f:
	malevich_data = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in malevich_data:
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

### ALL DATA

for ref, data, titl, conc_color, warm_color, cold_color in [
	('D47',      D47_data,      'This study',             (0,0,0,.67), (.9,0,0,1), (0,.2,1,1)),
	('malevich', malevich_data, 'Malevich et al. (2019)', (0,0,0,.25), (.9,0,0,.75), (0,.2,1,.5)),
	]:
	
	G = [r for r in data if 'Twoa23' in r]

	for depth in ['', '_500m', '_1500m']:
		
		if ref == 'malevich' and depth == '_1500m':
			continue
		
		if depth:
			G = [r for r in G if r['Twoa23_vs_Tiso_species'] == 'discordant']
		
		kwarm = [_ for _, r in enumerate(G) if r['Twoa23'+depth+'_vs_Tiso_species'+depth] == 'discordant' and r['Twoa23'+depth] < r['Tiso_species'+depth]]
		kcold = [_ for _, r in enumerate(G) if r['Twoa23'+depth+'_vs_Tiso_species'+depth] == 'discordant' and r['Twoa23'+depth] > r['Tiso_species'+depth]]
		kconc = [_ for _, r in enumerate(G) if r['Twoa23'+depth+'_vs_Tiso_species'+depth] == '']

		Twoa = array([r['Twoa23'+depth] for r in G])
		eTwoa = 1.96 * array([r['SE_Twoa23'+depth] for r in G])
		Tiso = array([r['Tiso_species'+depth] for r in G])
		eTiso = 1.96 * array([r['SE_Tiso_species'+depth] for r in G])

		fig = figure(figsize = (3.7, 3.7))
		subplots_adjust(0.15, 0.15, 0.96, 0.96)

		kw = dict(ms = 4)

		if depth:
			if kwarm:
				plot(Twoa[kwarm], Tiso[kwarm], 'x', mec = warm_color, label = f'{len(kwarm)} warm samples still discordant', lw = 1.5, **kw)
			plot(Twoa[kconc], Tiso[kconc], '+', mec = conc_color, label = f'{len(kconc)} samples no longer discordant', lw = 0.7, zorder = 0, **kw)
			plot(Twoa[kcold], Tiso[kcold], 'x', mec = cold_color, label = f'{len(kcold)} {"cold " if kwarm else ""}samples still discordant', lw = 1.5, **kw)
		else:
			if kwarm:
				plot(Twoa[kwarm], Tiso[kwarm], 'x', mec = warm_color, label = f'warm discordant (N = {len(kwarm)})', lw = 1.5, **kw)
			plot(Twoa[kconc], Tiso[kconc], '+', mec = conc_color, label = f'concordant (N = {len(kconc)})', lw = 0.7, zorder = 0, **kw)
			plot(Twoa[kcold], Tiso[kcold], 'x', mec = cold_color, label = f'{"cold " if kwarm else ""}discordant (N = {len(kcold)})', lw = 1.5, **kw)

		xi = [-8, 32]
		plot(xi, xi, 'k-', dashes = (8,2,2,2), lw = 0.7, zorder = -100)
		axis([xi[0],xi[1],xi[0],xi[1]])

		gca().xaxis.set_major_locator(ticker.MultipleLocator(5))
		gca().yaxis.set_major_locator(ticker.MultipleLocator(5))
		
		legend(fontsize = 8, labelspacing = 0.3)
		if depth:
			xlabel('Atlas T over 0-500 m depth range (°C)')
			ylabel('Oxygen-18 estimate of T (°C)')
			txt = ', now assuming calcification\ndepths of 0–500 m for previously discordant data'
		else:
			xlabel('Atlas T over calcification depth range (°C)')
			ylabel('Oxygen-18 estimate of T (°C)')
			txt = ', assuming\nspecies-specific calcification depths'
		text(.99, 0.01, titl + txt, transform = gca().transAxes, va = 'bottom', ha = 'right', size = 8, weight = 'bold')
		savefig(f'../output/T18_plot_{ref}{depth}.pdf')
		close(fig)

for ref, data in [
	('D47', D47_data),
# 	('malevich', malevich_data),
	]:
	
	G = [r for r in data if 'Twoa23' in r]

	for depth in ['', '_500m', '_1500m']:
		
		if depth:			
			G = [r for r in G if r['Twoa23_vs_Tiso_species'] == 'discordant']
		
		kwarm = [_ for _, r in enumerate(G) if r['Twoa23'+depth+'_vs_Tiso_species'+depth] == 'discordant' and r['Twoa23'+depth] < r['Tiso_species'+depth]]
		kcold = [_ for _, r in enumerate(G) if r['Twoa23'+depth+'_vs_Tiso_species'+depth] == 'discordant' and r['Twoa23'+depth] > r['Tiso_species'+depth]]
		kconc = [_ for _, r in enumerate(G) if r['Twoa23'+depth+'_vs_Tiso_species'+depth] == '']

		Twoa = array([r['Twoa23'+depth] for r in G])
		eTwoa = 1.96 * array([r['SE_Twoa23'+depth] for r in G])
		Tbot = array([r['Tbottom_woa23'] for r in G])
		eTbot = 1.96 * array([r['SE_Tbottom_woa23'] for r in G])
		Tiso = array([r['Tiso_species'+depth] for r in G])
		eTiso = 1.96 * array([r['SE_Tiso_species'+depth] for r in G])

		if not depth:
			species = sorted({r['Species'] for r in G if r['Twoa23'+depth+'_vs_Tiso_species'+depth] == 'discordant'})
			sp_markers = {s: m for s,m in zip(species, ['o', 'd', (3,0,0), (3,0,180), (3,0,90), (3,0,-90), (4,0,45)])}

		fig = figure(figsize = (3.7, 3.7))
		subplots_adjust(0.15, 0.15, 0.96, 0.96)

		kweb = dict(
			marker = 'None',
			ls = 'None',
			capsize = 2,
			elinewidth = 0.75,
			capthick = 0.75,
			)

		errorbar(
			Twoa[kcold+kwarm],
			Twoa[kcold+kwarm],
			eTwoa[kcold+kwarm],
			ecolor = 'k',
			**({} if depth else dict(label = 'Calcif. T from atlas')),
			**kweb,
			)		

		errorbar(
			Twoa[kcold+kwarm],
			Tiso[kcold+kwarm],
			eTiso[kcold+kwarm],
			ecolor = (.8,0,0),
			**({} if depth else dict(label = 'Calcif. T from $^{18}O$')),
			**kweb,
			)

		errorbar(
			Twoa[kcold+kwarm],
			Tbot[kcold+kwarm],
			eTbot[kcold+kwarm],
			ecolor = (0,.2,1),
			**({} if depth else dict(label = 'Bottom T from atlas')),
			**kweb,
			)
			
		for s in species:
			gs = s.split(' ')
			gs[0] = gs[0][0] + '.'
			gs = ' '.join(gs)
			M = [r for k,r in enumerate(G) if k not in kconc and r['Species'] == s]
			if M:
				plot(
					[r['Twoa23'+depth] for r in M],
					[r['Tiso_species'+depth] for r in M],
					ls = 'None',
					marker = sp_markers[s],
					mec = (.8,0,0),
					mfc = 'w',
					mew = 0.75,
					label = f'{gs} (N = {len(M)})',
					ms = 4,
					)

		xi = [-4, 32]
		plot(xi, xi, 'k-', dashes = (8,2,2,2), lw = 0.7, zorder = -100)
		axis([xi[0],xi[1],xi[0],xi[1]])

		gca().xaxis.set_major_locator(ticker.MultipleLocator(5))
		gca().yaxis.set_major_locator(ticker.MultipleLocator(5))

		legend(fontsize = 7, labelspacing = 0.5, handlelength = 0.8)

		if depth:
			xlabel('Mean atlas T over 0-500 m depth range (°C)')
			ylabel('T (°C)')
			txt = 'This study, now assuming calcification\ndepths of 0–500 m for previously discordant data'
		else:
			xlabel('Mean atlas T over calcification depth range (°C)')
			ylabel('T (°C)')
			txt = 'This study, assuming\nspecies-specific calcification depths'
		text(.99, 0.01, txt, transform = gca().transAxes, va = 'bottom', ha = 'right', size = 8, weight = 'bold')
		savefig(f'../output/T18_plot_errorbars_{ref}{depth}.pdf')
		close(fig)

with open('../output/saved_values.py', 'a') as fid:
	N = len(malevich_data)
	Nd = len([r for r in malevich_data if r['Twoa23_vs_Tiso_species'] == 'discordant'])
	Nc = len([r for r in malevich_data if r['Twoa23_vs_Tiso_species'] == ''])
	Nd500m = len([r for r in malevich_data if r['Twoa23_500m_vs_Tiso_species_500m'] == 'discordant'])
	Nc500m = len([r for r in malevich_data if r['Twoa23_500m_vs_Tiso_species_500m'] == ''])
	fid.write(f'\nNmalevich = {N}')
	fid.write(f'\nNmalevich_discordant = {Nd}')
	fid.write(f'\nNmalevich_concordant = {Nc}')
	fid.write(f'\nNmalevich_discordant_500m = {Nd500m}')
	fid.write(f'\nNmalevich_concordant_500m = {Nc500m}')

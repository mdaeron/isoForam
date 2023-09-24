#! /usr/bin/env python3
'''
Compile species-specific offsets from Kim & O'Neil (1997) regression.
'''

from csv import DictReader
from pylab import *
from scipy.stats import t as tstudent
from matplotlib import ticker
from scipy import stats
import sys
from scipy.interpolate import interp1d

sys.path.append('../3_species_lookup')
from species_lookup import species_lookup

import warnings
warnings.simplefilter('error', RuntimeWarning)

style.use('../mydefault.mplstyle')

with open('../2_compile_d18O_data/foram_d18O_compilation.csv') as f:
	data = list(DictReader(f))

for r in data:
	for k in r:
		if k not in ['Species', 'Ref', 'Type']:
			r[k] = float(r[k])
	if ' ' not in r['Species']:
		r['Species'] = r['Species'] + ' sp.'

for r in data:
	r['klna18'] = 1000 * log((1000 + r['d18Oc']) / (1000 + r['d18Osw_VSMOW']) * 1.03092)
	r['invT'] = 1/(273.15 + r['T'])

for r in data:
	if r['Species'] == 'Globigerinoides sacculifer':
		r['Species'] = species_lookup['Globigerinoides sacculifer']

speciesref = sorted({(r['Species'], r['Ref']) for r in data})

KON97_intercepts = {}
KON97_intercepts_by_ref = {}

print('Generating d18O plots...')

for grouping, KI in [
	('by ref', KON97_intercepts_by_ref),
	('by species', KON97_intercepts),
	]:
	for step in [0,1]:
		if step:
			# in second step, combine Globorotalia truncatulinoides (s.) [fine] and [coarse]
			for r in data:
				if 'Globorotalia truncatulinoides (s.)' in r['Species']:
					r['Species'] = 'Globorotalia truncatulinoides (s.)'
				if 'Globorotalia truncatulinoides (d.)' in r['Species']:
					r['Species'] = 'Globorotalia truncatulinoides (d.)'
				if 'Globorotalia inflata' in r['Species']:
					r['Species'] = 'Globorotalia inflata'
	
		for sp, ref in [
			('Globorotalia truncatulinoides (d.)', 'Lončarić et al. (2006)'),
			('Globorotalia truncatulinoides (s.)', 'Lončarić et al. (2006)'),
			('Globorotalia inflata', 'Lončarić et al. (2006)'),
			] if step == 1 else speciesref:
		# in second step, only process Globorotalia truncatulinoides (s.)


			if grouping == 'by ref':
				G = [r for r in data if sp == r['Species'] and ref == r['Ref']]
				H = [r for r in data if sp != r['Species'] or  ref != r['Ref']]
				refs = ref
				print(f'\tPlotting {sp} from {ref}...')
			else:
				G = [r for r in data if sp == r['Species']]
				H = [r for r in data if sp != r['Species']]
				refs = sorted({r['Ref'] for r in G})
				refs = '\n'.join(refs)
				if len(G) == 0:
					print(f'\tSkipping {sp} because no records were found.')
					continue
				print(f'\tPlotting {sp} from all refs...')

			fig = figure()

			kw = dict(
				lw = 0.9,
				)
			xi = array([1/310, 1/268])

			yi = 17570 * xi - 29.13
			plot(xi, yi, 'k-', label = 'Isotopic equilibrium according to Daëron et al. (2019)', **kw)

			yi = 18030 * xi - 32.17
			plot(xi, yi, 'k-', label = 'Synthetic precipitates of Kim & O\'Neil (1997)', dashes = (8,2,2,2), **kw)

			kw = dict(
				ls = 'None',
				mfc = 'None',
				ms = 5,
				)

			plot(
				[r['invT'] for r in H],
				[r['klna18'] for r in H],
				marker = '+',
				mew = 0.75,
				label = 'other species',
				mec = (1,.8,.8),
				zorder = -100,
				**kw
				)

			plot(
				[r['invT'] for r in G],
				[r['klna18'] for r in G],
				marker = 'o',
				mew = 1.0,
				label = f'{sp}\n{refs}',
				mec = 'k',
				**kw
				)

			legend()
			axis([xi[0], xi[-1], None, None])
			Ti = [30, 20, 10, 0]
			xticks([1/(273.15+t) for t in Ti])
			gca().xaxis.set_ticklabels([f'{t:.0f} °C' for t in Ti])
			xlabel('1/T')
			ylabel('1000 ln($^{18}α$)')

			X = array([r['invT'] for r in G])
			Y = array([r['klna18'] for r in G])
			R = Y - 18030 * X
			
			KI_key = sp if grouping == 'by species' else f'{sp} - {ref}'
			KI[KI_key] = {
				'N': len(R),
				'avg': R.mean(),
				'sd': R.std(ddof = 1),
				'se': R.std(ddof = 1) / len(R)**.5,
				'Type': G[0]['Type'],
				'Grouping': 'species',
				'invT': array([r['invT'] for r in G]),
				'klna18': array([r['klna18'] for r in G]),
				'R': R[:],
				'Refs': refs.replace('\n', '; '),
				}
	
			title(f'1000 log($^{{18}}$α) = 18030 / T - {KI[KI_key]["avg"]:.02f} ± {KI[KI_key]["se"]:.02f} (1SE)')
	
			savefig(f"../2_compile_d18O_data/plots/{grouping}/{KI_key}.pdf")
			close(fig)


KON97_intercepts = { # remove all data from Lončarić et al. (2006) except G. truncatulinoides and G. inflata
	x: KON97_intercepts[x]
	for x in sorted(
		{r['Species'] for r in data if r['Ref'] != 'Lončarić et al. (2006)'}
		| {'Globorotalia truncatulinoides (s.)', 'Globorotalia truncatulinoides (d.)', 'Globorotalia inflata'})
	}

primary_observations = [k for k in KON97_intercepts]

G = [KON97_intercepts[k] for k in KON97_intercepts if 'Hoeglundina' not in k]

for t in ['benthic', 'planktic'][1:]:
	R = array([r for x in G for r in x['R'] if x['Type'] == t])
	if t == 'planktic':
		label = 'all planktics'
	else:
		label = 'all calcitic benthics'
	KON97_intercepts[label] = {
		'N': len(R),
		'avg': R.mean(),
		'sd': R.std(ddof = 1),
		'se': R.std(ddof = 1) / len(R)**.5,
		'Type': t,
		'Grouping': 'b/p',
		'R': R[:],
		'Refs': '; '.join(sorted({x['Refs'] for x in G})),
		}


for genus in sorted({k.split(' ')[0] for k in KON97_intercepts}):
	if genus != 'all':
		G = [KON97_intercepts[k] for k in KON97_intercepts if genus == k.split(' ')[0]]
		R = array([r for x in G for r in x['R']])
		KON97_intercepts[f'{genus} spp.'] = {
			'N': len(R),
			'avg': R.mean(),
			'sd': R.std(ddof = 1),
			'se': R.std(ddof = 1) / len(R)**.5,
			'Type': G[0]['Type'],
			'Grouping': 'genus',
			'R': R[:],
			'Refs': '; '.join(sorted({_ for x in G for _ in x['Refs'].split('; ')})),
			}

G = [KON97_intercepts[k] for k in KON97_intercepts if k.split(' ')[0] in ['Cibicidoides', 'Planulina'] and not k.endswith(' spp.')]
R = array([r for x in G for r in x['R']])
KON97_intercepts['Cibicidoides + Planulina'] = {
	'N': len(R),
	'avg': R.mean(),
	'sd': R.std(ddof = 1),
	'se': R.std(ddof = 1) / len(R)**.5,
	'Type': 'benthic',
	'Grouping': 'b/p',
	'R': R[:],
	'Refs': '; '.join(sorted({x['Refs'] for x in G})),
	}

if __name__ == '__main__':
	
	
	### PLOT OF ALL RESIDUALS TO THE SPECIES-SPECIFIC 18O CALIBRATIONS
	
	fig = figure(figsize = (3,2.5))
	subplots_adjust(.05, .2, .95, .95)
	
	R = []
	for k in primary_observations:
		R += [r-KON97_intercepts[k]['R'].mean() for r in KON97_intercepts[k]['R']]

	kw = dict(
		bins = linspace(-1,1,41),
		)
	hist(R, histtype = 'stepfilled', lw = 0, color = (.2,.6,1), **kw)
	hist(R, histtype = 'step', lw = 1, color = 'k', **kw)

	yticks([])
	xlabel('1000 ln($^{18}$α) residuals')
	gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.1f}$' if x else '$0$'))
	
	ymax = axis()[-1]
	axis([None, None, 0, ymax*1.7])

	rmse = (array(R)**2).mean()**.5

	res = linspace(0,1,1001)
	pres = [len([r for r in R if abs(r)<_res]) / len(R) for _res in res]
	f = interp1d(pres, res)
	cl95 = f(0.95)
	_, cl, _ = errorbar(0, ymax * 1.1, None, cl95, ecolor = 'k', elinewidth = 1, capthick = 1, capsize = 4, ls = 'None', marker = 'None')
	cl[0].set_marker([(1.3,0), (2,1), (0,0), (2,-1), (1.3,0)])
	cl[1].set_marker([(-1.3,0), (-2,1), (0,0), (-2,-1), (-1.3,0)])
	text(0, ymax * 1.15, f'± {cl95:.2f} (95 %)', ha = 'center', va = 'bottom', size = 10, linespacing = 1)
	fill_between([-cl95, cl95], [0]*2, [ymax * 1.1]*2, zorder = -100, color = [.9]*3, lw = 0)

	text(0.02, 0.96, '(B) Variability from species-specific calibrations', size = 8.9, va = 'top', ha = 'left', transform = gca().transAxes)
	grid(False)
	savefig('oxygen18_residuals.pdf')
	close(fig)
	
	with open('../output/saved_values.py', 'a') as fid:
		fid.write(f'\nklnaRMSE = {rmse:.2f}')
		fid.write(f'\nklnaCL = {cl95:.2f}')
	

	### MAKE PLOT OF ALL SPECIES OBSERVATIONS
	
	fig = figure(figsize = (7.5, 9))
	subplots_adjust(.1, .05, .97, .97, .25, .25)

	axes = {}
	planktic_species = sorted([i for i in KON97_intercepts if KON97_intercepts[i]['Grouping'] == 'species' and KON97_intercepts[i]['Type'] == 'planktic'])
	benthic_species = sorted([i for i in KON97_intercepts if KON97_intercepts[i]['Grouping'] == 'species' and KON97_intercepts[i]['Type'] == 'benthic'])

	for k, sp in enumerate(planktic_species + benthic_species):
		axes[sp] = subplot(6,4,k+1)
		ax = axes[sp]
		if k % 4 == 0:
			ylabel('1000 ln($^{18}α$)')
		Ti = [30,15,0]
		xticks([(t+273.15)**-1 for t in Ti])
		ax.set_xticklabels([f"${t}\\,$°C" for t in Ti])
		ax.yaxis.set_major_locator(ticker.MultipleLocator(3))

		X = KON97_intercepts[sp]['invT']
		Y = KON97_intercepts[sp]['klna18']
		
		kw = dict(
			ls = 'None',
			marker = '+',
			mew = .5,
			mec = 'k',
			ms = 4,
			)
		plot(X, Y, **kw)

		kw = dict(
			lw = 0.8,
			color = (.2,.6,1),
			zorder = 0,
			)

		xi = array([(t+273.15)**-1 for t in [Ti[0]+4, Ti[-1]-3]])
		yi = 17570 * xi - 29.13
		plot(xi, yi, '-', label = 'Isotopic equilibrium according to Daëron et al. (2019)', **kw)
		yi = 18030 * xi - 32.17
		plot(xi, yi, '-', label = 'Synthetic precipitates of Kim & O\'Neil (1997)', dashes = (8,2,2,2), **kw)

		axis([xi[0], xi[-1], 25.9, 35.1])
		
		sptxt = sp if sp.endswith(' sp.') else (sp.split(' ')[0][0] + '. ' + ' '.join(sp.split(' ')[1:]))
		text(0.05, 0.95, sptxt, va = 'top', ha = 'left', size = 8, weight = 'bold', transform = ax.transAxes)
		reftxt = KON97_intercepts[sp]['Refs'].split('; ')
		for k, t in enumerate(reftxt):
			reftxt[k] = t.split(' fig')[0]
		reftxt = sorted(reftxt, key = len)
		reftxt = '\n'.join(reftxt)
		text(0.95, 0.05, reftxt, va = 'bottom', ha = 'right', size = 7, transform = ax.transAxes)

	legend(loc = 'center', bbox_to_anchor = (2.5, 0.5), fontsize = 8, labelspacing = 1.5, frameon = False)

	savefig('species_specific_alpha_plots.pdf')
	close(fig)



	### WRITE CSV OF SPECIES EFFECTS

	with open('species_offsets.csv', 'w') as fid:
		fid.write('Species,Genus,Group,A,B,SE_B,SD_observations,Refs')

		for figtype in ['errorbar', 'kde']:
			fig = figure(figsize = (5,6))
			subplots_adjust(.05, .1, .95, .9)

			k = 0

			for p, (color, KIlist) in enumerate([
				('g', sorted([i for i in KON97_intercepts if KON97_intercepts[i]['Grouping'] == 'species' and KON97_intercepts[i]['Type'] == 'planktic'],
					key = lambda s: s.replace(' sp,','zzz sp,'))),
				([.6,0,.6], sorted([i for i in KON97_intercepts if KON97_intercepts[i]['Grouping'] == 'species' and KON97_intercepts[i]['Type'] == 'benthic'],
					key = lambda s: s.replace(' sp,','zzz sp,'))),
				([0,0,1], sorted([i for i in KON97_intercepts if KON97_intercepts[i]['Grouping'] == 'genus'])),
				('k', sorted([i for i in KON97_intercepts if KON97_intercepts[i]['Grouping'] == 'b/p'])),
				]):

				k -= 1
				if p > 0:
					axhline(k, color = 'k', lw = .7)

				for sp in KIlist:
					k -= 1
					if figtype == 'kde' and sp == 'all planktics':
						k -= 0.5
			
					if figtype == 'kde':
						r = KON97_intercepts[sp]['R']+32.17
						kde = stats.gaussian_kde(r)
						xi = linspace(-2, 3, 5001)
						yi = kde(xi)
						xi, yi = xi[yi>1e-2], yi[yi>1e-2]
						yi *= sqrt(len(r))/40
						fill_between(
							xi,
							k+yi,
							k-yi,
							alpha = .5,
							zorder = 5,
							color = color,
			# 				color = 'k' if KON97_intercepts[sp]['Type'] == 'benthic' else 'b',
							)
						ha = 'left' if (xi.min() + xi.max())/2 < 0.3 else 'right'
						text(
							(xi.min() - .05) if ha == 'right' else (xi.max() + .05),
							k,
							f'{species_lookup[sp] if sp in species_lookup else sp}',
							va = 'center',
							ha = ha,
							size = 7,
	# 						weight = 'bold' if sp.endswith(' spp.') or sp.startswith('all ') or ' + ' in sp else 'normal',
							style = 'normal' if sp.startswith('all ') else 'italic',
							color = 'k',
							zorder = 100,
							)
					elif figtype == 'errorbar':
						X, eX, fX = (
							KON97_intercepts[sp]['avg'] + 32.17,
							KON97_intercepts[sp]['se'] * tstudent.ppf(1 - 0.05/2, KON97_intercepts[sp]['N']-1),
							KON97_intercepts[sp]['sd'] * tstudent.ppf(1 - 0.05/2, KON97_intercepts[sp]['N']-1),
							)

						if KON97_intercepts[sp]['Grouping'] == 'species':
							genus, species, other, refs = '', species_lookup[sp] if sp in species_lookup else sp, '', KON97_intercepts[sp]['Refs']
						elif KON97_intercepts[sp]['Grouping'] == 'genus':
							genus, species, other, refs = sp.split(' ')[0] + ' spp.', '', '', KON97_intercepts[sp]['Refs']
						elif KON97_intercepts[sp]['Grouping'] == 'b/p':
							genus, species, other, refs = '', '', sp, 'see above'
						
						if not sp.endswith(' sp.'):
							fid.write(f'\n{species},{genus},{other},18030,{KON97_intercepts[sp]["avg"]:.2f},{KON97_intercepts[sp]["se"]:.3f},{KON97_intercepts[sp]["sd"]:.2f},{refs}')

						
						kw = dict(
							ls = 'None',
							marker = 'None',
				# 			ecolor = 'k' if KON97_intercepts[sp]['Type'] == 'benthic' else 'b',
							ecolor = color,
							elinewidth = 6,
							capsize = 0,
							zorder = 100,
							)
						errorbar(X, k, None, fX, alpha = .25, **kw)
						errorbar(X, k, None, eX, **kw)
	# 					text(
	# 						X, k,
	# 						f'{sp}',
	# 						va = 'center',
	# 						ha = 'center',
	# 						size = 7,
	# 						weight = 'bold' if sp.endswith(' spp.') or sp.startswith('all ') or ' + ' in sp else 'normal',
	# 						style = 'normal' if sp.startswith('all ') else 'italic',
	# 						color = color,
	# 						zorder = 100,
	# 						)

						ha = 'left' if X < 0.3 else 'right'
						text(
							(X - fX - 0.05) if ha == 'right' else (X + fX + 0.05),
							k,
							f'{species_lookup[sp] if sp in species_lookup else sp}',
							va = 'center',
							ha = ha,
							size = 8,
	# 						weight = 'bold' if sp.endswith(' spp.') or sp.startswith('all ') or ' + ' in sp else 'normal',
							style = 'normal' if sp.startswith('all ') else 'italic',
							color = 'k',
							zorder = 100,
							)
			

	
			ytop = axis()[-1]

			axvline(0, color = 'k', lw = 1, alpha = .15, zorder = -100)
			text(0.1, ytop, "Kim & O'Neil (1997)\n\n", va = 'center', ha = 'right', color = 'k', alpha = .6, size = 8)

			# Grossman & Ku (1986)
			_d18c = linspace(-5,5,2001)
			_T = 20.6 - 4.34 * _d18c
			_d18c = _d18c[abs(_T-12.5) <= 12.5]
			_T = _T[abs(_T-12.5) <= 12.5]
			_klna18 = 1000 * log((1+_d18c/1000)*1.03092)
			_B = _klna18 - 18030/(_T+273.15) + 32.17
			gku0, gku25 = _B.max(), _B.min()
			axvspan(gku25, gku0, color = 'k', alpha = .1, lw = 0, zorder = -100)
			text((gku25+gku0)/2, ytop, "Grossman & Ku (1986)\n0–25 °C\n\n\n", va = 'center', ha = 'center', color = 'k', alpha = .6, size = 8)

			# Shackleton (1974)
			_d18c = linspace(-5,5,2001)
			_T = 16.9 - 4.38 * (_d18c + 0.27) + 0.1 * (_d18c + 0.27)**2
			_d18c = _d18c[abs(_T-12.5) <= 12.5]
			_T = _T[abs(_T-12.5) <= 12.5]
			_klna18 = 1000 * log((1+_d18c/1000)*1.03092)
			_B = _klna18 - 18030/(_T+273.15) + 32.17
			shack0, shack25 = _B.max(), _B.min()
			axvspan(shack25, shack0, color = 'k', alpha = .1, lw = 0, zorder = -100)
			text((shack25+shack0)/2, ytop, "Shackleton (1974)\n0–25 °C\n\n\n\n", va = 'center', ha = 'center', color = 'k', alpha = .6, size = 8)

			# text(arag, k*0.76, "Grossman & Ku (1986) at 10 °C\n", va = 'center', ha = 'center', color = 'g', alpha = .75, rotation = 90, size = 8, linespacing = 2)
			yticks([])
			xlabel('1000 ln($^{18}$α) offset from Kim & O\'Neil (1997)')

			gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'${x:+.1f}$' if x else '$0$'))

			axis([-1.2, 2.05, None, None])

			savefig(f'species_offsets_{figtype}.pdf')
			close(fig)
	
	with open('KON97_intercepts.py', 'w') as f:
		f.write('from numpy import array\n\nKON97_intercepts = ')
		f.write(str(KON97_intercepts))

# 	for sp in sorted(KON97_intercepts, key = lambda s: s.replace(' spp.', '   spp.').replace('Cibicides + Cibicidoides + Planulina', 'zzz').replace('all ', 'zzzzzz')):
# 		if sp.endswith(' sp.'):
# 			continue
# 		if sp.startswith('all '):
# 			genus, species, binned = '', '', sp
# 			refs = 'see above'
# 		elif ' + ' in sp:
# 			genus, species, binned = '', '', sp
# 			refs = 'see above'
# 		elif sp.endswith(' spp.'):
# 			genus, species, binned = sp[:-5], '', ''
# 			refs = '; '.join(sorted(set(KON97_intercepts[sp]["Refs"].split('; '))))
# 		else:
# 			genus, species, binned = '', sp, ''
# 			refs = '; '.join(sorted(set(KON97_intercepts[sp]["Refs"].split('; '))))
# 

#! /usr/bin/env python3
'''
Assign calcification temperatures.
'''

from csv import DictReader
from pylab import *
from scipy.odr import Model, RealData, ODR

style.use('../mydefault.mplstyle')

import sys
sys.path.append('../2_compile_d18O_data')
from KON97_intercepts import KON97_intercepts

sys.path.append('../3_species_lookup')
from species_lookup import species_lookup

tobeadded = []
for k in KON97_intercepts:
	if k in species_lookup:
		if species_lookup[k] not in KON97_intercepts:
			tobeadded += [k]

for k in tobeadded:
	KON97_intercepts[species_lookup[k]] = KON97_intercepts[k]

# for _ in sorted([k for k in KON97_intercepts]):
# 	print(f'{_:45}', species_lookup[_] if _ in species_lookup else '\t\t<--- not in species_lookup')
# exit()

sys.path.append('../4_species_depth')
from species_depth import species_depth

with open('../5_assign_d18Osw/bottom_d18Osw.csv') as f:
	bottom_d18Osw = {r['Site']: {k: r[k] for k in r} for r in DictReader(f)}
for r in bottom_d18Osw.values():
	for k in r:
		if k not in ['Site']:
			r[k] = float(r[k])

with open('../5_assign_d18Osw/d18Osw_estimates.csv') as f:
	d18Osw_estimates = {r['Sample']: {k: r[k] for k in r} for r in DictReader(f)}
for r in d18Osw_estimates.values():
	for k in r:
		if k not in ['Sample', 'Ref']:
			r[k] = float(r[k])

with open('../5_assign_d18Osw/d18Osw_estimates_500m.csv') as f:
	d18Osw_estimates_500m = {r['Sample']: {k: r[k] for k in r} for r in DictReader(f)}
for r in d18Osw_estimates.values():
	for k in r:
		if k not in ['Sample', 'Ref']:
			r[k] = float(r[k])

with open('../5_assign_d18Osw/d18Osw_estimates_1500m.csv') as f:
	d18Osw_estimates_1500m = {r['Sample']: {k: r[k] for k in r} for r in DictReader(f)}
for r in d18Osw_estimates.values():
	for k in r:
		if k not in ['Sample', 'Ref']:
			r[k] = float(r[k])

try:
	with open('../6_assign_atlas_T/bottom_temperatures.csv') as f:
		bTdata = {r['Site']: r for r in DictReader(f)}
	for r in bTdata.values():
		for k in r:
			if k != 'Site':
				r[k] = float(r[k])
except FileNotFoundError:
	bTdata = {}

try:
	with open('../6_assign_atlas_T/planktic_temperatures.csv') as f:
		pTdata = {r['Sample']: r for r in DictReader(f)}
except FileNotFoundError:
	pTdata = {}

try:
	with open('../6_assign_atlas_T/planktic_temperatures_500m.csv') as f:
		pTdata_500m = {r['Sample']: r for r in DictReader(f)}
except FileNotFoundError:
	pTdata_500m = {}
	
try:
	with open('../6_assign_atlas_T/planktic_temperatures_1500m.csv') as f:
		pTdata_1500m = {r['Sample']: r for r in DictReader(f)}
except FileNotFoundError:
	pTdata_1500m = {}
	
with open('../1_compile_D47_data/foram_D47_compilation.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]
	
for r in data:
	tobepopped = []
	for k in r:
		if k not in ['Sample', 'Species', 'Type', 'Site', 'Ref']:
			if r[k] == '':
				tobepopped.append(k)
			else:
				r[k] = float(r[k])
	for k in tobepopped:
		r.pop(k)


print('Estimating calcification temperatures...')

for r in data:
	ignore_benthic = False
	s = r['Sample']
	sp = species_lookup[r['Species']]
	d18Oc, sd18Oc = r['d18O_VPDB'], r['SE_d18O_VPDB']

	if r['Type'] == 'planktic':
		d18Ow = float(d18Osw_estimates[s]['d18Osw'])
		sd18Ow = float(d18Osw_estimates[s]['sd18Osw'])
	else:
		d18Ow = float(bottom_d18Osw[r['Site']]['d18Osw'])
		sd18Ow = float(bottom_d18Osw[r['Site']]['SE_d18Osw'])

	klna18 = 1000 * log(
		(1 + d18Oc/1000)
		/ (1 + d18Ow/1000)
		* 1.03092
		)
	sklna18 = ((sd18Oc / (1 + d18Oc/1000))**2 + (sd18Ow / (1 + d18Ow/1000))**2)**.5

	if r['Type'] == 'planktic':
		d18Ow_500m = float(d18Osw_estimates_500m[s]['d18Osw'])
		sd18Ow_500m = float(d18Osw_estimates_500m[s]['sd18Osw'])
		klna18_500m = 1000 * log(
			(1 + d18Oc/1000)
			/ (1 + d18Ow_500m/1000)
			* 1.03092
			)
		sklna18_500m = ((sd18Oc / (1 + d18Oc/1000))**2 + (sd18Ow_500m / (1 + d18Ow_500m/1000))**2)**.5

		d18Ow_1500m = float(d18Osw_estimates_1500m[s]['d18Osw'])
		sd18Ow_1500m = float(d18Osw_estimates_1500m[s]['sd18Osw'])
		klna18_1500m = 1000 * log(
			(1 + d18Oc/1000)
			/ (1 + d18Ow_1500m/1000)
			* 1.03092
			)
		sklna18_1500m = ((sd18Oc / (1 + d18Oc/1000))**2 + (sd18Ow_1500m / (1 + d18Ow_1500m/1000))**2)**.5

	A = 18030
	try:
		if sp == 'Trilobatus trilobus':
			I = KON97_intercepts['Trilobatus sacculifer']
			iso_species_key = 'Trilobatus sacculifer'
		else:
			I = KON97_intercepts[sp]
			iso_species_key = sp
# 		print(f'\tFound species offset for {iso_species_key}.')
	except KeyError:
		try:
			I = KON97_intercepts[sp.split(' ')[0] + ' spp.']
			iso_species_key = sp.split(' ')[0] + ' spp.'
			print(f'\tUsing genus offset ({iso_species_key}) for {sp}.')
		except KeyError:
			
			if r['Type'] == 'benthic':
				ignore_benthic = True
# 				print(f'\tIgnored (benthic) {species_lookup[r["Species"]]}')
			else:
				I = KON97_intercepts[{
#	 				'benthic': 'all calcitic benthics',
					'planktic': 'all planktics',
					}[r['Type']]]
				print(f'\tUsing "all planktics" offset for {species_lookup[r["Species"]]}.')
				iso_species_key = 'all planktics'
	if not ignore_benthic:
		B = I['avg']
		sB = I['se']
	
		r['Tiso_species_offset'] = iso_species_key
		r['Tiso_species'] = A / (klna18 - B) - 273.15
		r['SE_Tiso_species'] = (sB**2 + sklna18**2)**.5 * A / (klna18 - B)**2
		if r['Type'] == 'planktic':
			r['Tiso_species_500m'] = A / (klna18_500m - B) - 273.15
			r['SE_Tiso_species_500m'] = (sB**2 + sklna18_500m**2)**.5 * A / (klna18_500m - B)**2
			r['Tiso_species_1500m'] = A / (klna18_1500m - B) - 273.15
			r['SE_Tiso_species_1500m'] = (sB**2 + sklna18_1500m**2)**.5 * A / (klna18_1500m - B)**2

# 		print(f'{iso_species_key:>45} B = {B+32.17:.2f}')

	if r['Type'] == 'planktic':
		a, b = 19178.76270137247, -36.54443793149699 # T-dependent Mulitza regression
		T_mulitza = a / (klna18 - b) - 273.15
		sT_mulitza = sklna18 * a / (klna18 - b)**2
		r['Tiso_Mulitza03'], r['SE_Tiso_Mulitza03'] = T_mulitza, sT_mulitza

		a, b = 18030, -32.17 # KON97
		T_kon97 = a / (klna18 - b) - 273.15
		sT_kon97 = sklna18 * a / (klna18 - b)**2
		r['Tiso_KON97'], r['SE_Tiso_KON97'] = T_kon97, sT_kon97
		
		if r['Sample'] in pTdata:
			r['Twoa23'] = float(pTdata[r['Sample']]['Twoa23'])
			r['SE_Twoa23'] = float(pTdata[r['Sample']]['SD_Twoa23'])
			r['Twoa23_vs_Tiso_species'] = pTdata[r['Sample']]['Discordant']
			r['Twoa23_500m_vs_Tiso_species_500m'] = pTdata_500m[r['Sample']]['Discordant']
			r['Twoa23_500m'] = float(pTdata_500m[r['Sample']]['Twoa23'])
			r['SE_Twoa23_500m'] = float(pTdata_500m[r['Sample']]['SD_Twoa23'])
			r['Twoa23_1500m_vs_Tiso_species_1500m'] = pTdata_1500m[r['Sample']]['Discordant']
			r['Twoa23_1500m'] = float(pTdata_1500m[r['Sample']]['Twoa23'])
			r['SE_Twoa23_1500m'] = float(pTdata_1500m[r['Sample']]['SD_Twoa23'])

	elif r['Type'] == 'benthic':
		if r['Site'] in bTdata:
			r['Twoa23'] = bTdata[r['Site']]['mean_annual_Tbottom']
			r['SE_Twoa23'] = bTdata[r['Site']]['SE_mean_annual_Tbottom']
			r['Tpub'] = r['Tatlas']
# 			r['SE_Tpub'] = r['sTatlas']
			r['SE_Tpub'] = 0.5

	if r['Site'] in bTdata:
		r['Tbottom_woa23'] = bTdata[r['Site']]['mean_annual_Tbottom']
		r['SE_Tbottom_woa23'] = bTdata[r['Site']]['SE_mean_annual_Tbottom']

	if r['Type'] == 'planktic':
		r['Tk18'] = d18Osw_estimates[s]['T']
		r['SE_Tk18'] = d18Osw_estimates[s]['sT']
		r['Tk18_500m'] = d18Osw_estimates[s]['T']
		r['SE_Tk18_500m'] = d18Osw_estimates[s]['sT']
		r['Tk18_1500m'] = d18Osw_estimates[s]['T']
		r['SE_Tk18_1500m'] = d18Osw_estimates[s]['sT']
	else:
		r['Tk18'] = float(bottom_d18Osw[r['Site']]['T'])
		r['SE_Tk18'] = float(bottom_d18Osw[r['Site']]['SE_T'])

with open('foram_D47_calibration_data.csv', 'w') as fid:
	ff = [_.split(':') for _ in [
		'Sample:s',
		'Species:s',
		'Type:s',
		'Site:s',
		'Lat:.2f',
		'Lon:.2f',
		'Depth:s',
		'd13C_VPDB:.3f',
		'SE_d13C_VPDB:.3f',
		'd18O_VPDB:.3f',
		'SE_d18O_VPDB:.3f',
		'D47_ICDES:.4f',
		'SE_D47_ICDES:.4f',
		'Tpub:.2f',
		'SE_Tpub:.2f',
		'Twoa23:.2f',
		'SE_Twoa23:.2f',
		'Twoa23_500m:.2f',
		'SE_Twoa23_500m:.2f',
		'Twoa23_1500m:.2f',
		'SE_Twoa23_1500m:.2f',
		'Tbottom_woa23:.2f',
		'SE_Tbottom_woa23:.2f',
		'Tk18:.2f',
		'SE_Tk18:.2f',
		'Tk18_500m:.2f',
		'SE_Tk18_500m:.2f',
		'Tk18_1500m:.2f',
		'SE_Tk18_1500m:.2f',
		'Tiso_species_offset:s',
		'Tiso_species:.2f',
		'SE_Tiso_species:.2f',
		'Tiso_species_500m:.2f',
		'SE_Tiso_species_500m:.2f',
		'Tiso_species_1500m:.2f',
		'SE_Tiso_species_1500m:.2f',
		'Tiso_KON97:.2f',
		'SE_Tiso_KON97:.2f',
		'Tiso_Mulitza03:.2f',
		'SE_Tiso_Mulitza03:.2f',
		'Twoa23_vs_Tiso_species:s',
		'Twoa23_500m_vs_Tiso_species_500m:s',
		'Twoa23_1500m_vs_Tiso_species_1500m:s',
		'Ref:s',
		]]
	fields = [f[0] for f in ff]

	fid.write(','.join(fields))
	for r in data:
		r['Species'] = species_lookup[r['Species']]
		r['Depth'] = str(r['Depth']) if r['Type'] == 'benthic' else f"({species_depth[r['Species']]['zmin']:.0f}-{species_depth[r['Species']]['zmax']:.0f})"
		fid.write('\n' + ','.join([f'{r[f]:{fmt}}' if f in r else '' for f,fmt in ff]))

if __name__ == '__main__':

	benthic_sites = sorted({(r['Site'], r['Tatlas']) for r in data if r['Type'] == 'benthic'}, key = lambda _: _[1])
	benthic_sites = [_[0] for _ in benthic_sites]
	
	_benthic_data_ = []


	fig = figure(figsize = (4, 6))
	subplots_adjust(.05, .09, .96, .98, .05)

	axes = [subplot(13, 1, k+1) for k in range(13)]

	k = 0
	for site in benthic_sites:
		noaxis = True
		species = sorted({
			r['Species'] for r in data
			if r['Site'] == site and r['Type'] == 'benthic'
			})

		
		j = 0
		for sp in species:
			G = [r for r in data if r['Site'] == site and r['Species'] == sp]
			Tatlas = G[0]['Tatlas']
			try:
				Twoa = G[0]['Twoa23']
				sTwoa = G[0]['SE_Twoa23']
			except:
				Twoa = 0.
				sTwoa = 0.
			N = sum([r['N'] for r in G])

			d18Oc = sum([r['d18O_VPDB']*r['N']/N for r in G])
			sd18Oc = sum([(r['SE_d18O_VPDB'] * r['N']/N)**2 for r in G])**.5
			try:
				star = ''
				A, B = 18030, KON97_intercepts[sp]['avg']
# 				print(f'Species {sp} is in KON97_intercepts')
			except KeyError:
				try:
					A, B = 18030, KON97_intercepts[sp.split(' ')[0] + ' spp.']['avg']
					star = '*'
# 					print(f'Species {sp} + spp. is in KON97_intercepts')
				except KeyError:
					print(f'Species {sp} not found')
					continue

			j += 1
			if noaxis:
				noaxis = False
				ax = axes[k]
				k += 1

			d18Ow = float(bottom_d18Osw[site]['d18Osw'])
			sd18Ow = float(bottom_d18Osw[site]['SE_d18Osw'])

			klna18 = 1000 * log(
				(1 + d18Oc/1000)
				/ (1 + d18Ow/1000)
				* 1.03092
				)
			sklna18 = (
				(sd18Oc / (1 + d18Oc/1000))**2
				+ (sd18Ow / (1 + d18Ow/1000))**2
				)**.5

			Tiso = A/(klna18 - B) - 273.15
			sTiso = sklna18 * A / (klna18 - B)**2

			ymax = j+0.8
			ax.errorbar(Tiso, j, xerr = 1.96*sTiso, ecolor = 'k', elinewidth = 1, capsize = 2, capthick = 1)
			_sp_ = sp.split(' ')[0][0] + '. ' + sp.split(' ')[1] + star
			ax.text(
				Tiso + ((-1.96*sTiso - 0.5) if Tiso > 22 else (1.96*sTiso + 0.5)),
				j,
				_sp_,
				style = 'italic',
				size = 6,
				va = 'center',
				ha = 'right' if Tiso > 22 else 'left',
				)
	
			_benthic_data_ += [{
				'Site': site,
				'Species': sp,
				'T18': Tiso,
				'SE_T18': sTiso,
				'D47': sum([r['D47_ICDES']*r['N']/N for r in G]),
				'SE_D47': sum([(r['SE_D47_ICDES'] * r['N']/N)**2 for r in G])**.5,
				'Ref': G[0]['Ref'],
				}]
	
		if not noaxis:
			ax.set_yticks([])
		
			ax.axvspan(Twoa-1.96*sTwoa, Twoa+1.96*sTwoa, color = 'k', alpha = 0.15, lw = 0)
			ax.axvline(Tatlas, color = (0,.5,1), alpha = 0.5, lw = 2)
	
			x1, x2 = ax.axis()[:2]
			x0 = (x1+x2)/2
# 			ax.set_xlim((x0-6, x0+6))
			ax.set_xlim((-4, 27))
			ax.set_ylim((0.2, ymax))
			ax.tick_params(length = 0, labelbottom = False)

	ax.set_xlabel('T (°C)')
	ax.tick_params(length = 3, labelbottom = True)

	fig.savefig('benthic_T18.pdf')
	close(fig)

	with open('benthic_T18.csv', 'w') as fid:
		fid.write(','.join(_benthic_data_[0].keys()))
		for _ in _benthic_data_:
			fid.write(f"\n{_['Site']},{_['Species']},{_['T18']:.2f},{_['SE_T18']:.2f},{_['D47']:.4f},{_['SE_D47']:.4f},{_['Ref']}")

# 	peral_sites = sorted({
# 		r['Site'] for r in data if 'Peral' in r['Ref'] and r['Type'] == 'benthic'
# 		})
# 
# 	piasecki_sites = sorted({
# 		r['Site'] for r in data if 'Piasecki' in r['Ref']
# 		})
# 
# 	isobenthics = {}
# 	Tcomparison = []
# 		
# 	for _ in sorted({r['Species'] for r in data if r['Type'] == 'benthic'}):
# 		print(_)
# 
# 	for site in peral_sites + piasecki_sites:
# 		
# 		G = [r for r in data if r['Site'] == site and r['Type'] == 'benthic']
# 		for sp, mrk in [
# 			('Uvigerina mediterranea', (3,0,90)),
# 			('Cibicides wuellerstorfi', (3,0,-90)),
# 			('Cibicidoides pachyderma', 'o'),
# 			('Planulina foveolata', (4,0,45)),
# 			('Planulina ariminensis', (4,0,0)),
# 			('Hoeglundina elegans', (3,0,180)),
# 			('Uvigerina peregrina', (3,0,0)),
# 			]:
# 			H = [r for r in G if sp == species_lookup[r['Species']]]
# 			if len(H):
# 				N = sum([r['N'] for r in H])
# 				d18Oc = sum([r['d18O_VPDB']*r['N']/N for r in H])
# 				sd18Oc = sum([(r['SE_d18O_VPDB'] * r['N']/N)**2 for r in H])**.5
# 				if site not in isobenthics:
# 					isobenthics[site] = {}
# 				try:
# 					A, B = 18030, KON97_intercepts[sp]['avg']
# 				except KeyError:
# 					A, B = 18030, KON97_intercepts[sp.split(' ')[0] + ' spp.']['avg']				
# 				d18Ow = float(bottom_d18Osw[H[0]['Site']]['d18Osw'])
# 				sd18Ow = float(bottom_d18Osw[H[0]['Site']]['SE_d18Osw'])
# 				klna18 = 1000 * log(
# 					(1 + d18Oc/1000)
# 					/ (1 + d18Ow/1000)
# 					* 1.03092
# 					)
# 					
# 				sklna18 = (
# 					(sd18Oc / (1 + d18Oc/1000))**2
# 					+ (sd18Ow / (1 + d18Ow/1000))**2
# 					)**.5
# 
# 				Tiso = A/(klna18 - B) - 273.15
# 				sTiso = sklna18 * A / (klna18 - B)**2
# 				isobenthics[site][sp] = Tiso
# # 				print(f'{site:<20} - {sp:<30}: {Tiso:.1f}  (atlas: {H[0]["Tatlas"]})')
# 				Tcomparison.append([float(H[0]["Tatlas"]), float(H[0]["sTatlas"]), Tiso, sTiso, sp, site, mrk])
# 	
# 	fig = figure(figsize = (5,5))
# 	for mrk in {_[-1] for _ in Tcomparison}:
# 		X, sX, Y, sY, sp, __, _ = zip(*[_ for _ in Tcomparison if _[-1] == mrk])
# 		errorbar(X, Y, 1.96*array(sY), ls = 'None', marker = 'None', ecolor = 'k', elinewidth = 1, capthick = 1, capsize = 2, zorder = 100)
# 		plot(X, Y, ls = 'None', marker = mrk, mfc = 'w', mec = 'k', label = sp[0], zorder = 200)
# 
# 	for site in piasecki_sites:
# 		try:
# 			X, Y, ___, __, _ = zip(*[_ for _ in Tcomparison if _[-2] == site])
# 			plot([X[0], X[0]], [min(Y), max(Y)], '-k', alpha = .5, zorder = -1000)
# 		except ValueError:
# 			pass
# 	
# 	Tmin, Tmax = -4, 27
# 	plot([Tmin, Tmax], [Tmin, Tmax], 'r-', alpha = .25)
# 	axis([Tmin, Tmax,Tmin, Tmax])
# 
# 	legend(loc = 4)
# 	xlabel('T atlas')
# 	ylabel('T iso')
# 		
# 	savefig('benthicT_comparison.pdf')
# 	close(fig)
# 
# 
# if len(bTdata):
# 	dataplot = list({
# 		(
# 			r['Site'],
# 			r['Tatlas'],
# 			r['sTatlas'],
# 			float(bTdata[r['Site']]['mean_annual_Tbottom']),
# 			float(bTdata[r['Site']]['SE_mean_annual_Tbottom']),
# 			r['Ref']
# 			)
# 		for r in data if r['Type'] == 'benthic'
# 		})
# 
# 	dataplot = [{k:v for k,v in zip(['Site', 'Told', 'sTold', 'Tnew', 'sTnew', 'Ref'], _)} for _ in dataplot]
# 
# 		
# 	fig = figure(figsize = (6,6))
# 
# 	for ref in ['Piasecki et al. (2019)', 'Peral et al. (2018)']:
# 		G = [_ for _ in dataplot if _['Ref'] == ref]
# 		X = [_['Told'] for _ in G]
# 		eX = [_['sTold']*1.96 for _ in G]
# 		Y = [_['Tnew'] for _ in G]
# 		eY = [_['sTnew']*1.96 for _ in G]
# 	
# 		errorbar(X, Y, eY, eX if 'Peral' in ref else None, ecolor = 'k', capsize = 2, elinewidth = 1, capthick = 1, ls = 'None', marker = 'None')
# 		plot(X, Y, 'o', mec = 'k', mfc = ([.75]*3) if 'Peral' in ref else 'w', label = ref)
# 
# 	xlabel('Originally published bottom temperature (°C)')
# 	ylabel('Redetermined bottom temperature (°C)')
# 
# 	plot([-3,20], [-3,20], 'k-', alpha = .5, dashes = (8,2,2,2))
# 	axis([-3,20,-3,20])
# 	legend()	
# 	savefig('Tatlas_new_vs_old.pdf')
# 	close(fig)
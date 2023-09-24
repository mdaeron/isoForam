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

# with open('../5_assign_d18Osw/bottom_d18Osw.csv') as f:
# 	bottom_d18Osw = {r['Site']: {k: r[k] for k in r} for r in DictReader(f)}
# for r in bottom_d18Osw.values():
# 	for k in r:
# 		if k not in ['Site']:
# 			r[k] = float(r[k])

with open('d18Osw_estimates.csv') as f:
	d18Osw_estimates = {r['Sample']: {k: r[k] for k in r} for r in DictReader(f)}
for r in d18Osw_estimates.values():
	for k in r:
		if k not in ['Sample', 'Ref']:
			r[k] = float(r[k])

with open('d18Osw_estimates_500m.csv') as f:
	d18Osw_estimates_500m = {r['Sample']: {k: r[k] for k in r} for r in DictReader(f)}
for r in d18Osw_estimates.values():
	for k in r:
		if k not in ['Sample', 'Ref']:
			r[k] = float(r[k])

try:
	with open('bottom_temperatures.csv') as f:
		bTdata = {r['Site']: r for r in DictReader(f)}
	for r in bTdata.values():
		for k in r:
			if k != 'Site':
				r[k] = float(r[k])
except FileNotFoundError:
	bTdata = {}

try:
	with open('planktic_temperatures.csv') as f:
		pTdata = {r['Sample']: r for r in DictReader(f)}
except FileNotFoundError:
	pTdata = {}

try:
	with open('planktic_temperatures_500m.csv') as f:
		pTdata_500m = {r['Sample']: r for r in DictReader(f)}
except FileNotFoundError:
	pTdata_500m = {}

with open('malevich_2019_s1.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]
		
for r in data:
	for k in r:
		if k not in ['corename', 'doi', 'foram', 'reference', 'species', 'subspecies']:
			r[k] = float(r[k])

for k,r in enumerate(data):
	r['Sample'] = f'{k:04.0f}'
	r['Species'] = species_lookup[r['species']]
	r['Lat'] = r['latitude']
	r['Lon'] = r['longitude']
	r['Ref'] = r['reference']
	r['Type'] = 'planktic'
	r['zmin'] = species_depth[r['Species']]['zmin']
	r['zmax'] = species_depth[r['Species']]['zmax']
	r['d18O_VPDB'] = r['d18oc']
	r['SE_d18O_VPDB'] = 0.1


print('Estimating calcification temperatures...')

for r in data:
	ignore_benthic = False
	s = r['Sample']
	sp = species_lookup[r['Species']]
	d18Oc, sd18Oc = r['d18oc'], 0.1

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

# 	elif r['Type'] == 'benthic':
# 		if r['Site'] in bTdata:
# 			r['Twoa23'] = bTdata[r['Site']]['mean_annual_Tbottom']
# 			r['SE_Twoa23'] = bTdata[r['Site']]['SE_mean_annual_Tbottom']
# 			r['Tpub'] = r['Tatlas']
# # 			r['SE_Tpub'] = r['sTatlas']
# 			r['SE_Tpub'] = 0.5

	if 'Site' in r and r['Site'] in bTdata:
		r['Tbottom_woa23'] = bTdata[r['Site']]['mean_annual_Tbottom']
		r['SE_Tbottom_woa23'] = bTdata[r['Site']]['SE_mean_annual_Tbottom']

	if r['Type'] == 'planktic':
		r['Tk18'] = d18Osw_estimates[s]['T']
		r['SE_Tk18'] = d18Osw_estimates[s]['sT']
		r['Tk18_500m'] = d18Osw_estimates[s]['T']
		r['SE_Tk18_500m'] = d18Osw_estimates[s]['sT']
	else:
		r['Tk18'] = float(bottom_d18Osw[r['Site']]['T'])
		r['SE_Tk18'] = float(bottom_d18Osw[r['Site']]['SE_T'])

with open('malevich_data_out.csv', 'w') as fid:
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
		'Tbottom_woa23:.2f',
		'SE_Tbottom_woa23:.2f',
		'Tk18:.2f',
		'SE_Tk18:.2f',
		'Tk18_500m:.2f',
		'SE_Tk18_500m:.2f',
		'Tiso_species_offset:s',
		'Tiso_species:.2f',
		'SE_Tiso_species:.2f',
		'Tiso_species_500m:.2f',
		'SE_Tiso_species_500m:.2f',
		'Tiso_KON97:.2f',
		'SE_Tiso_KON97:.2f',
		'Tiso_Mulitza03:.2f',
		'SE_Tiso_Mulitza03:.2f',
		'Twoa23_vs_Tiso_species:s',
		'Twoa23_500m_vs_Tiso_species_500m:s',
		'Ref:s',
		]]
	fields = [f[0] for f in ff]

	fid.write(','.join(fields))
	for r in data:
		r['Species'] = species_lookup[r['Species']]
		r['Depth'] = str(r['Depth']) if r['Type'] == 'benthic' else f"({species_depth[r['Species']]['zmin']:.0f}-{species_depth[r['Species']]['zmax']:.0f})"
		fid.write('\n' + ','.join([f'{r[f]:{fmt}}' if f in r else '' for f,fmt in ff]))

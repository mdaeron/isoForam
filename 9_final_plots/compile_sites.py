#! /usr/bin/env python3

from csv import DictReader

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

with open('../1_compile_D47_data/sites.csv') as f:
	sites = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in sites:
	for k in r:
		if k not in ['Site', 'Ref']:
			r[k] = float(r[k])
# sites = sorted(sites, key = lambda x: x['Site'])

with open('../5_assign_d18Osw/bottom_d18Osw.csv') as f:
	d18Osw_sites = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in d18Osw_sites:
	for k in r:
		if k not in ['Site']:
			r[k] = float(r[k])
# d18Osw_sites = sorted(d18Osw_sites, key = lambda x: x['Site'])

with open('../6_assign_atlas_T/bottom_temperatures.csv') as f:
	T_sites = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in T_sites:
	for k in r:
		if k not in ['Site']:
			r[k] = float(r[k])
# T_sites = sorted(T_sites, key = lambda x: x['Site'])

with open('../8_seawater_chemistry/bottom_seawater_chemistry.csv') as f:
	bswchem_sites = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in bswchem_sites:
	for k in r:
		if k not in ['Site']:
			r[k] = float(r[k])
# bswchem_sites = sorted(bswchem_sites, key = lambda x: x['Site'])

print('Compiling sites:')

with open('sites.csv', 'w') as f:
	f.write('Site,Lat,Lon,Depth,T,SE_T,Salinity,SE_Salinity,Omega_cc,SE_Omega_cc,d18Osw,SE_d18Osw,T_breitkreuz,SE_T_breitkreuz,Benthics/Planktics,Ref')
	for a in sites:
		b = [s for s in T_sites if s['Site'] == a['Site']][0]
		c = [s for s in bswchem_sites if s['Site'] == a['Site']][0]
		d = [s for s in d18Osw_sites if s['Site'] == a['Site']][0]
		site = a['Site']
		print(f'\tProcessing {site}')
		lat = a['Lat']
		lon = a['Lon']
		depth = a['Depth']
		ref = a['Ref']
		salinity = c['Salinity']
		se_salinity = c['SE_Salinity']
		omega = c['Omega_cc']
		se_omega = c['SE_Omega_cc']
		T = b['mean_annual_Tbottom']
		se_T = b['SE_mean_annual_Tbottom']
		Tb = d['T']
		se_Tb = d['SE_T']
		d18 = d['d18Osw']
		se_d18 = d['SE_d18Osw']

		_types = sorted({r['Type'] for r in data if r['Site'] == site})
		if len(_types) == 2:
			types = 'both'
		else:
			types = _types[0]

		f.write(f'\n{site},{lat:.2f},{lon:.2f},{depth:.0f},{T:.2f},{se_T:.2f},{salinity:.3f},{se_salinity:.3f},{omega:.3f},{se_omega:.3f},{d18:.2f},{se_d18:.2f},{Tb:.2f},{se_Tb:.2f},{types},{ref}')
		
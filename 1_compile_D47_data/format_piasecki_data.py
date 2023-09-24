#! /usr/bin/env python3
'''
Reformat new version (I-CDES reprocessed) of Piasecki et al. (2019) data
'''

from numpy import *
from csv import DictReader

with open('piasecki_samples_reprocessed.csv') as fid:
	newdata = list(DictReader(fid))

data = []
for r in newdata:
	if r['Sample']:
		try:
			d13 = array(d13)
			d18 = array(d18)
			D47 = array(D47)
			data.append(dict(
				Site = site,
				Sample = sample,
				d13C_VPDB = d13*1,
				d18O_VPDB = d18*1,
				D47_ICDES = D47*1,
				N = len(D47),
				))
		except NameError:
			pass
	if r['Site']:
		site = r['Site']
	if r['Sample']:
		sample = r['Sample']
		d13, d18, D47 = [], [], []
	d13.append(float(r['d13C VPDB (Final)']))
	d18.append(float(r['d18O VPDB (Final)']))
	D47.append(float(r['D47 I-CDES (Final)']))
else:
	d13 = array(d13)
	d18 = array(d18)
	D47 = array(D47)
	data.append(dict(
		Site = site,
		Sample = sample,
		d13C_VPDB = d13*1,
		d18O_VPDB = d18*1,
		D47_ICDES = D47*1,
		N = len(D47),
		))

newsites = {
	'13MC': '13MC-G',
	'19MC': '19MC-G',
	'50MC': '50MC-G',
	'53MC': '53MC-G',
	'89MC': '89MC-G',
	'94MC': '94MC-G',
	'GS06-144-19': 'GS06-144-19',
	'GS07-150-17-2': 'GS07-150-17-2',
	'GS07-150-22-1': 'GS07-150-22-1',
	'MP43-BC': 'MP43-BC',
	'MP46-MC': 'MP46-MC',
	'S0213-54-4': 'SO213-54-4',
	'SO213-71-2': 'SO213-71-2',
	}

Nf = 0
d13_residuals = []
d18_residuals = []
D47_residuals = []

for r in data:
	d13_residuals += [_ - r['d13C_VPDB'].mean() for _ in r['d13C_VPDB']]
	d18_residuals += [_ - r['d18O_VPDB'].mean() for _ in r['d18O_VPDB']]
	D47_residuals += [_ - r['D47_ICDES'].mean() for _ in r['D47_ICDES']]
	Nf += r['N'] - 1

d13_SD = ((array(d13_residuals)**2).sum() / Nf)**.5
d18_SD = ((array(d18_residuals)**2).sum() / Nf)**.5
D47_SD = ((array(D47_residuals)**2).sum() / Nf)**.5

for r in data:
	r['d13C_VPDB'] = r['d13C_VPDB'].mean()
	r['d18O_VPDB'] = r['d18O_VPDB'].mean()
	r['D47_ICDES'] = r['D47_ICDES'].mean()
	r['SE_d13C_VPDB'] = d13_SD / r['N']**.5
	r['SE_d18O_VPDB'] = d18_SD / r['N']**.5
	r['SE_D47_ICDES'] = D47_SD / r['N']**.5

	if 'Cpachyderma' in r['Sample']:
		r['Species'] = 'Cibicidoides pachyderma'
	elif 'Helegans' in r['Sample']:
		r['Species'] = 'Hoeglundina elegans'

	elif 'Lconvergens' in r['Sample']:
		r['Species'] = 'Lenticulina convergens'
	elif 'Lenticulinaconvergens' in r['Sample']:
		r['Species'] = 'Lenticulina convergens'
	elif 'Lenticulina_convergens' in r['Sample']:
		r['Species'] = 'Lenticulina convergens'
	elif 'lenticulinaconvergens' in r['Sample']:
		r['Species'] = 'Lenticulina convergens'

	elif 'Lenticulinaiota' in r['Sample']:
		r['Species'] = 'Lenticulina iota'
	elif 'lenticulinaiota' in r['Sample']:
		r['Species'] = 'Lenticulina iota'

	elif 'Pserata' in r['Sample']:
		r['Species'] = 'Pyrgo serrata'
	elif 'pyrgoserrata' in r['Sample']:
		r['Species'] = 'Pyrgo serrata'

	elif 'Uperegina' in r['Sample']:
		r['Species'] = 'Uvigerina peregrina'
	elif 'univ' in r['Sample']:
		r['Species'] = 'Uvigerina peregrina'
	elif 'uperigini' in r['Sample']:
		r['Species'] = 'Uvigerina peregrina'
	elif 'Upergina' in r['Sample']:
		r['Species'] = 'Uvigerina peregrina'

	elif 'Amphisterigina_radiata' in r['Sample']:
		r['Species'] = 'Amphisterigina radiata'
	elif 'pyrgo' in r['Sample']:
		r['Species'] = 'Pyrgo spp'
	elif 'Cpachy' in r['Sample']:
		r['Species'] = 'Cibicidoides pachyderma'
	elif 'Lenticulinaspp' in r['Sample']:
		r['Species'] = 'Lenticulina spp'
	elif 'Amphistegina_lessoni' in r['Sample']:
		r['Species'] = 'Amphistegina lessoni'
	elif 'umbonatus' in r['Sample']:
		r['Species'] = 'Oridorsalus umbonatus'
	elif 'Melonis_barleanum' in r['Sample']:
		r['Species'] = 'Melonis barleannum'
	elif 'Planulina_ariminensis' in r['Sample']:
		r['Species'] = 'Planulina ariminensis'
	elif 'U_med' in r['Sample']:
		r['Species'] = 'Uvigerina mediterranea'
	elif 'Cmundulus' in r['Sample']:
		r['Species'] = 'Cibicidoides mundulus'
	elif 'cibpachyderma' in r['Sample']:
		r['Species'] = 'Cibicidoides pachyderma'
	elif 'melonis' in r['Sample']:
		r['Species'] = 'Melonis spp'
	elif 'Cibicides_lobatus' in r['Sample']:
		r['Species'] = 'Cibicides lobatus'
	elif 'med' in r['Sample']:
		r['Species'] = 'Uvigerina mediterranea'

with open('piasecki_samples.csv', 'w') as fid:
	fid.write('Sample,Site,Species,N,d13C,sd13C,d18Oc,sd18Oc,D47,sD47')
	for r in data:
		if 'Species' in r:
			fid.write(f"\n{r['Sample']},{newsites[r['Site']]},{r['Species']},{r['N']},{r['d13C_VPDB']:.3f},{r['SE_d13C_VPDB']:.3f},{r['d18O_VPDB']:.3f},{r['SE_d18O_VPDB']:.3f},{r['D47_ICDES']:.4f},{r['SE_D47_ICDES']:.4f}")
		else:
			print(r['Site'], r['Sample'], r['N'])

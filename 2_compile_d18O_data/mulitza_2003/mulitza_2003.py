#! /usr/bin/env python3

import numpy as np
from csv import DictReader
from glob import glob

species = {
	'G. bull.'     : 'Globigerina bulloides',
	'G. ruber (w)' : 'Globigerinoides ruber white',
	'G. sacc.'     : 'Globigerinoides sacculifer',
	'N. pachy (l)' : 'Neogloboquadrina pachyderma (s.)',
	}

data = []
for csvfile in glob('Mulitza_2003_d18O_plank_forams/csv/*.csv'):
	with open(csvfile) as f:
		data += [{
			'Species': species[r['Species']],
			'T': r['T (∞C)'],
			'd18Osw_VSMOW': r['d18O(SMOW)'],
			'd18Oc': r['d18O'],
			} for r in DictReader(f)]

fields = [k for k in data[0]]

with open('mulitza_data.csv', 'w') as f:
	f.write(','.join(fields))
	for r in data:
		f.write('\n')
		f.write(','.join([r[k] if k in r else '' for k in fields]))


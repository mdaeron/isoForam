#! /usr/bin/env python3

import numpy as np
from csv import DictReader
from glob import glob
from pylab import *


data = []
for g in sorted(glob('Marchitto_2014_A*.csv')):
	with open(g) as f:
		data += [{k: r[k] for k in r} for r in DictReader(f)]

fields = []
for r in data:
	for k in r:
		if k not in fields:
			fields.append(k)

with open('marchitto_data.csv', 'w') as f:
	f.write(','.join(fields))
	for r in data:
		f.write('\n')
		f.write(','.join([r[k] if k in r else '' for k in fields]))

with open('marchitto_data.csv') as f:
	data = list(DictReader(f))

for r in data:
	for k in r:
		if k != 'Core':
			if r[k] == '':
				r[k] = None
			else:
				r[k] = float(r[k])

fig = figure()
foramfields = [k for k in fields if k.startswith('d18O_')]
for ff, mrk in zip(foramfields, 'osDdv^'):
	G = [r for r in data if r[ff] is not None]
	X = [r['T'] for r in G]
	Y = [r[ff] - r['d18Osw_VSMOW'] for r in G]
	plot(X, Y, mrk, alpha = .4, label = ff[5:].replace('_', '. '))

legend()
xlabel('T')
ylabel('δ$_c$ - δ$_w$')
show()
close(fig)
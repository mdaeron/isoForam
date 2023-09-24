#! /usr/bin/env python3
'''
Reprocess data from Peral et al. (2018) in I-CDES using D47crunch
'''

from D47crunch import *
from pylab import *
import os, sys, datetime
from matplotlib.patches import Ellipse
from peral_Tpubatlas import peral_Tpubatlas

### Process data with different size fractions treated as different samples:
rawdata = D47data()
	
rawdata.read('peral_2018_rawdata.csv')
rawdata.wg()
rawdata.crunch()
rawdata.standardize()
rawdata.summary(dir = 'output', verbose = True)
rawdata.table_of_sessions(dir = 'output', verbose = True)
rawdata.table_of_samples(dir = 'output', verbose = True)

### Read binned sample names from metadata:
metadata = read_csv('peral_2018_metadata.csv')
N = len(metadata)

### Define sample groups to bin together:
groups = {l['Sample']: [] for l in metadata}
for r in groups:
	groups[r] += [s for s in rawdata.samples if s.startswith(r)]

samples_new, D47_new, CM_new = rawdata.combine_samples(groups)

CM_klna = zeros((N, N))
klna = zeros((N,))
cores = []

for r in metadata:
	X = [x['d18O_VSMOW'] for s in groups[r['Sample']] for x in rawdata.samples[s]['data']]
	r['d18O_VSMOW'] = mean(X)
	r['d18c'] = (1000+r['d18O_VSMOW']) / rawdata.ALPHA_18O_ACID_REACTION / 1.03092 - 1000
	r['sd18c'] = rawdata.repeatability['r_d18O_VSMOW'] / rawdata.ALPHA_18O_ACID_REACTION / 1.03092 / len(X)**.5
# 	print(r['Sample'], r['d18Oc_VPDB'], r['d18c'])
	if r['Sample'] == '2FPA1_UviMed':
		r['d18c'] -= 0.47
	r['klna'] = 1000 * log((1000+r['d18c'])/(1000+r['d18w'])*1.03092)
	r['sklna'] = sqrt(
		1e6 * r['sd18c']**2/(1000+r['d18c'])**2
		+ 1e6 * r['sd18w']**2/(1000+r['d18w'])**2
		)
	r['T18'] = 18030 / (r['klna'] + 32.17) - 273.15
	if r['Sample'] == '2FPA1_UviMed':
		r['d18c'] += 0.47
	print(r['Sample'], r['T18'])

for i,sample1 in enumerate(samples_new):
	r = [u for u in metadata if u['Sample'] == sample1][0]
	cores += [sample1.split('_')[0]]
	klna[i] = r['klna']
	for j, sample2 in enumerate(samples_new):
		if i == j:
			CM_klna[i,i] = 1e6 * r['sd18c']**2/(1000+r['d18c'])**2 + 1e6 * r['sd18w']**2/(1000+r['d18w'])**2
			break
		else:
			core1 = sample1.split('_')[0]
			core2 = sample2.split('_')[0]
			if core1 == core2 and core1 != 'MOCOSED':
				CM_klna[i,j] = 1e6 * r['sd18w']**2/(1000+r['d18w'])**2
				CM_klna[j,i] = CM_klna[i,j]

T18 = 18030/(klna + 32.17) - 273.15

### T = 18030/(klna + 32.17) - 273.15
J = (-18030/(klna + 32.17)**2).reshape(N,1)

CM_T18 = CM_klna * (J.T * J)

def sortkey(x):
	return mean([l['T18'] for l in metadata if l['Sample'].split('_')[0] == x.split('_')[0]]) - (1 if x == 'MOCOSED_CWuel' else 0)
o = sorted(range(N), key = lambda k: sortkey(samples_new[k]))

samples_new = array(samples_new)[o]
D47_new = D47_new[o]
klna = klna[o]
T18 = T18[o]

CM_new = CM_new[o,:]
CM_new = CM_new[:,o]
CM_klna = CM_klna[o,:]
CM_klna = CM_klna[:,o]
CM_T18 = CM_T18[o,:]
CM_T18 = CM_T18[:,o]

for s,D47,sD47,k, sk, T, sT in zip(samples_new, D47_new, diag(CM_new)**.5, klna, diag(CM_klna)**.5, T18, diag(CM_T18)**.5):
	print(f'{s:<20}{D47:8.4f}{sD47:8.4f}{k:8.2f}{sk:8.2f}{T:8.2f}{sT:8.2f}')

with open('peral_samples.csv', 'w') as f:
	f.write('Sample,N,d13C_VPDB,SE_d13C_VPDB,d18O_VPDB,SE_d18O_VPDB,D47_ICDES,SE_D47_ICDES,Tatlas,sTatlas')
	for ns, D47, sD47 in zip(samples_new, D47_new, diag(CM_new)**.5):
		Tpa, sTpa = peral_Tpubatlas[ns]['Tpubatlas'], peral_Tpubatlas[ns]['SE_Tpubatlas']
		os = [s for s in rawdata.unknowns if s.startswith(ns)]
		G = [r for r in rawdata if r['Sample'] in os]
		d13C = mean([r['d13C_VPDB'] for r in G])
		SE_d13C = G[0]['SE_d13C'] / len(G)**.5
		d18Oc = mean([r['d18O_VPDB'] for r in G])
		SE_d18Oc = G[0]['SE_d18O'] / len(G)**.5
		N = len(G)
		f.write(f'\n{ns},{N},{d13C:.3f},{SE_d13C:.3f},{d18Oc:.3f},{SE_d18Oc:.3f},{D47:.4f},{sD47:.4f},{Tpa:.2f},{sTpa:.2f}')
	

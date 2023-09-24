#! /usr/bin/env python3
from csv import DictReader

with open('../4_species_depth/depth_data.csv') as f:
	data = [{k: r[k] for k in r} for r in DictReader(f)]
	
for r in data:
	for k in [k for k in r if r[k] == '']:
		r.pop(k)
	for k in r:
		if k not in ['Species', 'Ref']:
			r[k] = float(r[k])

species_depth = {
	r['Species']: r
	for r in sorted(data, key = lambda x: x['Species'])
	}

species_depth['Trilobatus trilobus'] = species_depth['Trilobatus sacculifer'].copy()
species_depth['Globigerinoides ruber'] = {
	'zmin': min(
		species_depth['Globigerinoides ruber white']['zmin'],
		species_depth['Globigerinoides ruber pink']['zmin'],
		),
	'zmax': max(
		species_depth['Globigerinoides ruber white']['zmax'],
		species_depth['Globigerinoides ruber pink']['zmax'],
		),
	}

if __name__ == '__main__':

	import sys
	sys.path.append('../3_species_lookup')
	from species_lookup import species_lookup

	with open('../1_compile_D47_data/foram_D47_compilation.csv') as f:
		species = sorted({species_lookup[r['Species']] for r in DictReader(f) if r['Type'] == 'planktic'})

	for s in species:
		if s in species_depth:
			print(f'{s:>50}: {str(species_depth[s]["zmin"]):>10} <--> {str(species_depth[s]["zmax"]):<10}')
		else:
			print(f'{s}: missing depth info.')
	

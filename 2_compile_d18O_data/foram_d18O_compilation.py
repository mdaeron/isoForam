#! /usr/bin/env python3
'''
Assign calcification temperatures.
'''

from csv import DictReader
import sys

sys.path.append('../4_species_depth')
from species_depth import species_depth
sys.path.append('../3_species_lookup')
from species_lookup import species_lookup

species = {
	'G_glutinata_fine'            : 'Globigerina glutinata',
	'G_inflata_fine'              : 'Globorotalia inflata [fine]',
	'G_inflata_coarse'            : 'Globorotalia inflata [coarse]',
	'G_trilobus_fine'             : 'Globigerinoides trilobus [fine]',
	'G_trilobus_coarse'           : 'Globigerinoides trilobus [coarse]',
	'G_ruber_fine'                : 'Globigerinoides ruber [fine]',
	'G_ruber_coarse'              : 'Globigerinoides ruber [coarse]',
	'G_truncatulinoides_d_fine'   : 'Globorotalia truncatulinoides (d.) [fine]',
	'G_truncatulinoides_d_coarse' : 'Globorotalia truncatulinoides (d.) [coarse]',
	'G_truncatulinoides_s_fine'   : 'Globorotalia truncatulinoides (s.) [fine]',
	'G_truncatulinoides_s_coarse' : 'Globorotalia truncatulinoides (s.) [coarse]',
	'C_pachyderma'                : 'Cibicidoides pachyderma',
	'P_ariminensis'               : 'Planulina ariminensis',
	'P_foveolata'                 : 'Planulina foveolata',
	'U_peregrina'                 : 'Uvigerina peregrina',
	'H_elegans'                   : 'Hoeglundina elegans',
	'C_wuellerstorfi'             : 'Cibicidoides wuellerstorfi',
	}


data = []

with open('grossman_ku_1986/grossman_ku_data.csv') as f:
	data += [r | {'Ref': 'Grossman & Ku (1986)', 'Type': 'benthic'} for r in DictReader(f) if not r['Species'].startswith('#')]

with open('loncaric_2006/loncaric_data.csv') as f:
	for r in DictReader(f):
		K = [k for k in r if k.startswith('d18O_') and r[k]]
		for k in K:
			sp = species[k[5:]] if k[5:] in species else k[5:]
			sp = species_lookup[sp]
			if sp not in species_depth:
# 				print(f"No depth info for {sp}.")
				zmin, zmax = 0, 1000
			else:
				zmin, zmax = species_depth[sp]['zmin'], species_depth[sp]['zmax']
			if float(r['Depth']) >= zmin and float(r['Depth']) <= zmax:
				data += [{
					'Species'     : species[k[5:]] if k[5:] in species else k[5:],
					'T'           : r['T'],
					'd18Osw_VSMOW': r['d18Osw_VSMOW'],
					'd18Oc'       : r[k],
					'Ref'         : 'Lončarić et al. (2006)',
					'Type'        : 'planktic',
					}]


with open('marchitto_2014/marchitto_data.csv') as f:
	for r in DictReader(f):
		K = [k for k in r if k.startswith('d18O_') and r[k]]
		for k in K:
			data += [{
				'Species'     : species[k[5:]] if k[5:] in species else k[5:],
				'T'           : r['T'],
				'd18Osw_VSMOW': r['d18Osw_VSMOW'],
				'd18Oc'       : r[k],
				'Ref'         : 'Marchitto et al. (2014)',
				'Type'        : 'benthic',
				}]

with open('keigwin_1998/keigwin_data.csv') as f:
	for r in DictReader(f):
		if r['d18Osw_VSMOW']:
			K = [k for k in r if k.startswith('d18O_') and r[k]]
			for k in K:
				data += [{
					'Species'     : species[k[5:]] if k[5:] in species else k[5:],
					'T'           : r['T'],
					'd18Osw_VSMOW': r['d18Osw_VSMOW'],
					'd18Oc'       : r[k],
					'Ref'         : 'Keigwin (1998)',
					'Type'        : 'benthic',
					}]

with open('mulitza_2003/mulitza_data.csv') as f:
	data += [r | {'Ref': 'Mulitza et al. (2003)', 'Type': 'planktic'} for r in DictReader(f)]

with open('mccorkle_2008/mccorkle_data.csv') as f:
	data += [r | {'Type': 'benthic'} for r in DictReader(f)]

with open('barras_2010/barras_2010.csv') as f:
	data += [r | {'Ref': 'Barras et al. (2010)', 'Type': 'benthic'} for r in DictReader(f)]

with open('spero_lea_1996/spero_lea_1996.csv') as f:
	data += [r | {'Ref': 'Spero & Lea (1996)', 'Type': 'planktic'} for r in DictReader(f)]

with open('spero_lea_1993/spero_lea_1993.csv') as f:
	data += [r | {'Ref': 'Spero & Lea (1993)', 'Type': 'planktic'} for r in DictReader(f) if 'sacculifer' not in r['Species']]

with open('bemis_1998/bemis_1998_orbulina_universa.csv') as f:
	bemisdata = [r for r in DictReader(f)]
	for r in bemisdata:
		for L in ['HL', 'LL']:
			data += [{
				'Species': r['Species'],
				'T': r['T'],
				'd18Osw_VSMOW': '0.0',
				'd18Oc'       : r[f'd18Oc-d18Ow {L}'],
				'Ref'         : 'Bemis et al. (1998) fig 1a',
				'Type'        : 'planktic',
				}]

with open('bemis_1998/bemis_1998_globigerina_bulloides.csv') as f:
	bemisdata = [r for r in DictReader(f)]
	for r in bemisdata:
		data += [{
			'Species': r['Species'],
			'T': r['T'],
			'd18Osw_VSMOW': '0.0',
			'd18Oc'       : r[f'd18Oc-d18Ow'],
			'Ref'         : 'Bemis et al. (1998) fig 1c',
			'Type'        : 'planktic',
			}]
			

fields = [k for k in data[0]]

with open('foram_d18O_compilation.csv', 'w') as f:
	f.write(','.join(fields))
	for r in data:
		f.write('\n')
		f.write(','.join([r[k] for k in fields]))

# for x in sorted({f'{r["Species"]} - {r["Type"]} - {r["Ref"]}' for r in data}): print(x)

# for r in data:
# 	for k in r:
# 		if k not in ['Species', 'Ref']:
# 			r[k] = float(r[k])
# 		print(r)

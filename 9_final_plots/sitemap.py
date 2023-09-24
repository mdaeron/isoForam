#! /usr/bin/env python3

from csv import DictReader
import cartopy
import cartopy.crs as ccrs
from pylab import *
import warnings


style.use('../mydefault.mplstyle')

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


with open('sites.csv') as f:
	sites = [{k: r[k] for k in r} for r in DictReader(f)]	
for r in sites:
	for k in r:
		if k not in ['Site', 'Benthics/Planktics', 'Ref']:
			r[k] = float(r[k])

for s in sites:
	G = [r for r in data if r['Site'] == s['Site'] and r['Type'] == 'planktic']
	if G:
		s['planktics'] = len(G)
		H = [r for r in G if r['Twoa23_500m_vs_Tiso_species_500m'] == 'discordant']
		s['discordants_500m'] = len(H)
		s['discordants_ratio_500m'] = len(H)/len(G)
		H = [r for r in G if r['Twoa23_vs_Tiso_species'] == 'discordant']
		s['discordants'] = len(H)
		s['discordants_ratio'] = len(H)/len(G)




fig = figure(figsize = (5.6,3.))
subplots_adjust(0.05,0,.95,1)

proj = ccrs.Robinson()

ax = plt.axes(projection = proj)
ax.set_global()
# ax.set_extent([-180+1e-6, 180-1e-6, -90, 90], ccrs.PlateCarree())
# ax.stock_img()

# ax = axes(projection = ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, color = [.9]*3)
ax.add_feature(cartopy.feature.OCEAN, color = [.7]*3)
ax.coastlines(lw = 0.75)

# y, x = zip(*list({(s['Lat'], s['Lon']) for s in sites}))
# plot(x, y, 'wD', transform = ccrs.PlateCarree(), label = 'all core tops', mec = 'k', ms = 4, mew = 0.8)


warnings.filterwarnings("ignore")

y, x = zip(*[(s['Lat'], s['Lon']) for s in sites if s['Benthics/Planktics'] == 'benthic'])
plot(x, y, 'ws', transform = ccrs.PlateCarree(), label = 'benthic', mec = 'k', ms = 5, mew = 1)

y, x = zip(*[(s['Lat'], s['Lon']) for s in sites if s['Benthics/Planktics'] == 'planktic'])
plot(x, y, 'wo', transform = ccrs.PlateCarree(), label = 'planktic', mec = 'k', ms = 5, mew = 1)

y, x = zip(*[(s['Lat'], s['Lon']) for s in sites if s['Benthics/Planktics'] == 'both'])
plot(x, y, 'wD', transform = ccrs.PlateCarree(), label = 'both', mec = 'k', ms = 4.5, mew = 1)

legend(labelspacing = 0.2, handlelength = 0.5)
savefig('plots/sitemap.pdf')
close(fig)





fig = figure(figsize = (8,5))

proj = ccrs.Robinson()

ax = plt.axes(projection = proj)
ax.set_global()
# ax.set_extent([-180+1e-6, 180-1e-6, -90, 90], ccrs.PlateCarree())
# ax.stock_img()

# ax = axes(projection = ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND, color = [.6]*3)
ax.add_feature(cartopy.feature.OCEAN, color = [.8]*3)
ax.coastlines(lw = 0.75)

# y, x = zip(*list({(s['Lat'], s['Lon']) for s in sites}))
# plot(x, y, 'wD', transform = ccrs.PlateCarree(), label = 'all core tops', mec = 'k', ms = 4, mew = 0.8)


warnings.filterwarnings("ignore")

y, x = zip(*[(s['Lat'], s['Lon']) for s in sites if s['Benthics/Planktics'] == 'benthic'])
plot(x, y, 'ws', transform = ccrs.PlateCarree(), label = 'benthic', mec = 'k', ms = 5, mew = 0.8)

y, x = zip(*[(s['Lat'], s['Lon']) for s in sites if s['Benthics/Planktics'] == 'planktic'])
plot(x, y, 'wo', transform = ccrs.PlateCarree(), label = 'planktic', mec = 'k', ms = 5, mew = 0.8)

y, x = zip(*[(s['Lat'], s['Lon']) for s in sites if s['Benthics/Planktics'] == 'both'])
plot(x, y, 'wD', transform = ccrs.PlateCarree(), label = 'both', mec = 'k', ms = 4.5, mew = 0.8)

for s in sites:
	if 'discordants' in s:
		if s['discordants']:

			text(s['Lon'], s['Lat'],
				f"{s['discordants']}     \n",
				transform = ccrs.PlateCarree(),
				color = (1,0,0),
				ha = 'center',
				va = 'bottom',
				size = 7,
				weight = 'bold',
				linespacing = 0.25,
				)

			if s['discordants_500m']:
				text(s['Lon'], s['Lat'],
					f"{s['discordants_500m']}\n",
					transform = ccrs.PlateCarree(),
					color = (0,0,1),
					ha = 'center',
					va = 'bottom',
					size = 7,
					weight = 'bold',
					linespacing = 0.25,
					)
		
			text(s['Lon'], s['Lat'],
				f"     {s['planktics']}\n",
				transform = ccrs.PlateCarree(),
				color = (1,1,1),
				ha = 'center',
				va = 'bottom',
				size = 7,
				weight = 'bold',
				linespacing = 0.25,				
				)
		

legend()
title(
	'White numbers = total planktic samples at this site;\nRed = discordant samples; Blue = discordant samples even with 0-500m depth\n',
	size = 9,
	)
savefig('plots/sitemap_with_discordant_planktics.pdf')
close(fig)

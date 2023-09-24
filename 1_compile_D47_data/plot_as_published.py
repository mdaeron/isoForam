#! /usr/bin/env python3
'''
Compile foram peral_data from original sources.
'''

from csv import DictReader
from pylab import *

style.use('../mydefault.mplstyle')

with open('peral_aspublished.csv') as f:
	peral_data = [{k: float(r[k]) if k not in ['Core', 'Species/Fraction'] else r[k] for k in r} for r in DictReader(f)]
for r in peral_data:
	r['T'] = float(r['Tpub'])
	r['sT'] = float(r['sTpub'])
	r['D47'] = float(r['D47pub'])
	r['sD47'] = float(r['sD47pub'])

with open('../1_compile_D47_data/piasecki_samples_aspublished.csv') as f:
	piasecki_data = [r for r in DictReader(f)]

for r in piasecki_data:
	r['T'] = float(r['T'])
	r['sT'] = 0
	r['D47'] = float(r['D47'])
	r['sD47'] = float(r['sD47'])
	
with open('../1_compile_D47_data/meinicke_samples.csv') as f:
	meinicke_data = [r for r in DictReader(f)]

for r in meinicke_data:
	r['T'] = float(r['pubT'])
	r['sT'] = float(r['SE_pubT'])
	r['D47'] = float(r['pubD47'])
	r['sD47'] = float(r['SE_pubD47'])

fig = figure(figsize = (4,3.5))
subplots_adjust(.16, .1, .96, .96)

kw_errorbar_foram = dict(
	ls = 'None',
	marker = 'None',
	ecolor = 'k',
	alpha = 0.15,
	capsize = 0,
	elinewidth = 1,
	zorder = 100,
	)

kw_plot_foram = dict(
	ls = 'None',
	mew = .8,
	mec = 'k',
	)

kw_lsce = dict(
	mfc = [.7]*3,
	zorder = 200,
	)

kw_bergen = dict(
	mfc = [1]*3,
	zorder = 100,
	)

kw_plot_benthic = kw_plot_foram | dict(
	marker = 's',
	ms = 5,
	)

kw_plot_planktic = kw_plot_foram | dict(
	marker = 'o',
	ms = 5,
	)


'''PERAL BENTHICS'''
G = [
	r for r in peral_data
	if ('Cibicides' in r['Species/Fraction'] or 'Uvigerina' in r['Species/Fraction'])
	]

X = [(r['T']+273.15)**-2 for r in G]
eX = [1.96 * 2 * r['sT'] * (r['T']+273.15)**-3 for r in G]
Y = [r['D47'] for r in G]
eY = [r['sD47']*1.96 for r in G]

# errorbar(X, Y, eY, eX, **kw_errorbar_foram)
plot(X, Y, label = 'Peral et al. (2018) benthics', **kw_plot_benthic, **kw_lsce)
X18 = [] + X


'''PIASECKI'''
G = [r for r in piasecki_data]

X, Y, eX, eY = [], [], [], []
for site in sorted({r['Site'] for r in G}):
	H = [r for r in G if r['Site'] == site]
	X.append((H[0]['T']+273.15)**-2)
	eX.append(1.96*2*H[0]['sT']*(H[0]['T']+273.15)**-3)
	D47 = sum([r['D47']/r['sD47']**2 for r in H]) / sum([1/r['sD47']**2 for r in H])
	sD47 = sum([1/r['sD47']**2 for r in H])**-.5
	Y.append(D47)
	eY.append(sD47*1.96)

# errorbar(X, Y, eY, eX, **kw_errorbar_foram)
plot(X, Y, label = 'Piasecki et al. (2019) benthics', **(kw_plot_benthic | kw_bergen | dict(mfc = 'k')))
X19 = [] + X


'''PERAL PLANKTICS'''
G = [
	r for r in peral_data
	if ('Cibicides' not in r['Species/Fraction'] and 'Uvigerina' not in r['Species/Fraction'])
	]

X = [(r['T']+273.15)**-2 for r in G]
eX = [1.96 * 2 * r['sT'] * (r['T']+273.15)**-3 for r in G]
Y = [r['D47'] for r in G]
eY = [r['sD47']*1.96 for r in G]

# errorbar(X, Y, eY, eX, **kw_errorbar_foram)
plot(X, Y, label = 'Peral et al. (2018) planktics', **kw_plot_planktic, **kw_lsce)
X18 += X

'''MEINICKE'''
G = [r for r in meinicke_data]

X = [(r['T']+273.15)**-2 for r in G]
eX = [1.96 * 2 * r['sT'] * (r['T']+273.15)**-3 for r in G]
Y = [r['D47'] for r in G]
eY = [r['sD47']*1.96 for r in G]	

# errorbar(X, Y, eY, eX, **kw_errorbar_foram)
plot(X, Y, label = 'Meinicke et al. (2020) planktics', **kw_plot_planktic, **kw_bergen)
X20 = [] + X

leg1 = legend(fontsize = 8, loc = 'upper left')

xi = array([min(X18)*0.995, max(X18)*1.005])
yi = 41.63e3*xi + 0.2056
a, = plot(xi, yi, 'k-', label = 'Peral et al. (2018) regression', lw = 1, dashes = (3,1))
D47 = 41.63e3*273.15**-2 + 0.2056

xi = array([min(X19)*0.995, max(X19)*1.005])
yi = 46.0e3*xi + 0.159
b, = plot(xi, yi, 'k-', label = 'Piasecki et al. (2019) regression', lw = 0.75)
T1 = ((D47 - 0.159)/46.0e3)**-.5 - 273.15
D47_1 = 46.0e3*273.15**-2 + 0.159

print(f'At 0 °C, Piasecki - Peral is {1000*(D47_1-D47):+.1f} ppm')
print(f'At 0 °C, Piasecki - Peral corresponds to {T1:+.1f} °C')

xi = array([min(X20)*0.995, max(X20)*1.005])
yi = 39.7e3*xi + 0.2259
c, = plot(xi, yi, 'k-', label = 'Meinicke et al. (2020) regression', lw = 1, dashes = (6,1.5,2,1.5))
T2 = ((D47 - 0.2259)/39.7e3)**-.5 - 273.15
D47_2 = 39.7e3*273.15**-2 + 0.2259

print(f'At 0 °C, Meinicke - Peral is {1000*(D47_2-D47):+.1f} ppm')
print(f'At 0 °C, Meinicke - Peral corresponds to {T2:+.1f} °C')

legend(
	[a,b,c],
	[_.get_label() for _ in [a,b,c]],
	loc = 'lower right',
	fontsize = 8,
# 	labelspacing = 0.4,
# 	frameon = False,
	)
	
gca().add_artist(leg1)


xmin, xmax = 304**-2, 270**-2


Ti = [30,20,10,0]
xticks([(t+273.15)**-2 for t in Ti])
gca().set_xticklabels([f"${t}\\,$°C" for t in Ti])
ylabel(f'Δ$_{{47}}$ (CDES-25, ‰)')

axis([xmin, xmax, 0.645, 0.79])
gca().xaxis.set_minor_locator(FixedLocator([(_+273.15)**-2 for _ in linspace(-2,30,17)]))
grid(visible = True, alpha = 0.2, which = 'both')
savefig('studies_as_published.pdf')
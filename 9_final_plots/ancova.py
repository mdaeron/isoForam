from pylab import *
from scipy.stats import norm
from scipy.odr import Model, RealData, ODR


def ancova( X1, Y1, sX1, sY1, X2, Y2, sX2, sY2, q = 0.05, verbose = False, plots = False) :

	X1, Y1, sX1, sY1 = list(X1), list(Y1), list(sX1), list(sY1)
	X2, Y2, sX2, sY2 = list(X2), list(Y2), list(sX2), list(sY2)
	xi = linspace(27, 35)

	def f(B, x) : return B[0]*x + B[1]
	linear = Model(f)

	data = RealData( X1, Y1, sx = sX1, sy = sY1 )
	fit = ODR(data, linear, beta0=[40e3, 0.]).run()
	A1, B1, CM1 = fit.beta[0], fit.beta[1], fit.cov_beta
	if plots:
		yi = A1*xi+B1
		syi = (CM1[0,0]*xi**2 + 2*CM1[0,1]*xi + CM1[1,1])**.5
		fill_between(xi, yi+1.96*syi, yi-1.96*syi, color = 'k', alpha = 0.2)

	data = RealData( X2, Y2, sx = sX2, sy = sY2 )
	fit = ODR(data, linear, beta0=[40e3, 0.]).run()
	A2, B2, CM2 = fit.beta[0], fit.beta[1], fit.cov_beta
	if plots:
		yi = A2*xi+B2
		syi = (CM2[0,0]*xi**2 + 2*CM2[0,1]*xi + CM2[1,1])**.5
		fill_between(xi, yi+1.96*syi, yi-1.96*syi, color = 'k', alpha = 0.2)

	slope_deviation = abs(A2-A1)/(CM1[0,0] + CM2[0,0])**.5
	p_slope = 2*norm.cdf(-slope_deviation)
	if verbose:
		print(f'Slope weighted deviation is {slope_deviation:.2f} (p_slope = {p_slope}).')

	X = [x for x in X1] + [-x for x in X2]
	Y = Y1 + Y2
	sX = sX1 + sX2
	sY = sY1 + sY2
	
	def f(B, x) : return B[0]*x*sign(x) + B[1] + B[2] * (1-sign(x))/2
	bilinear = Model(f)
	data = RealData( X, Y, sx = sX, sy = sY )
	fit = ODR(data, bilinear, beta0=[40e3, 0., 0.]).run()
	A, B, D, CM = fit.beta[0], fit.beta[1], fit.beta[2], fit.cov_beta
	if plots:
		yi = A*xi+B
		plot(xi, yi, 'k-', alpha = 0.2)
		yi += D
		plot(xi, yi, 'k-', alpha = 0.2)

	if verbose:
		print( f'Vertical offset of {1000*D:.1f} ± {1000*CM[2,2]**.5:.1f} ppm.' )

	intercept_deviation = abs(D)/CM[2,2]**.5
	p_intercept = 2*norm.cdf(-intercept_deviation)
	if verbose:
		print(f'Intercept weighted deviation is {intercept_deviation:.2f} (p_intercept = {p_intercept}).')

	return p_slope, p_intercept	

#! /usr/bin/env python3

import numpy as np
from scipy.odr import Model, RealData, ODR

def YorkReg(X, Y, sX, sY, rXY=0., max_iter = 1e3, cvg_limit=1e-15, rel_cvg=True, verbose=False, rescale = False):
	'''
	Computes straight-line regression with errors in X and Y
	according to York et al. (2004):

	Y = A + BX

	Parameters
	----------
	X : array_like
		1-D array of length N corresponding to the observed X values.
	Y : array_like
		1-D array of length N corresponding to the observed Y values.
	sX : array_like or float
		1-D array of length N corresponding to the standard errors
		assigned to the elements of X. If specified as a float, all
		elements of sX will be assigned this value.
	sY : array_like or float
		1-D array of length N corresponding to the standard errors
		assigned to the elements of Y. If specified as a float, all
		elements of sY will be assigned this value.
	rXY : array_like or float, optional
		1-D array of length N corresponding to the correlation coefficients
		between the errors in each pair of values (X[i], Y[i]) with 0 <= i < N.
		rXY[i] = cov(X[i],Y[i]) / sX[i] / sY[i], thus -1 <= rXY[i] <= +1.
		If specified as a float, all elements of rXY will be assigned this value.
	cvg_limit : float
		The threshold change in slope that triggers a new iterative step
		in computing the best-fit slope.
	rel_cvg : boolean
		If rel_cvg is True, then iteration convergence is reached when *relative*
		changes in slope are smaller than cvg_limit; otherwise, convergence is
		reached when *absolute* changes in slope are smaller than cvg_limit.
	verbose : boolean
		If True, print out information regarding the regression convergence.
	'''

	X = np.asarray(X, dtype="float")
	Y = np.asarray(Y, dtype="float")
	sX = np.asarray(sX, dtype="float")
	if sX.size == 1:
		sX = X*0+sX
	sY = np.asarray(sY, dtype="float")
	if sY.size == 1:
		sY = Y*0+sY
	rXY = np.asarray(rXY, dtype="float")
	if rXY.size == 1:
		rXY = X*0+rXY

	k = (~np.isnan(X)) & (~np.isnan(Y)) & (~np.isnan(sX)) & (~np.isnan(sY)) & (~np.isnan(rXY))

	X = X[k]
	Y = Y[k]
	sX = sX[k]
	sY = sY[k]
	rXY = rXY[k]

	if verbose:
		print(f'rXY = {rXY}')
		print(f"YorkReg -> N = {len(X)}")

	oX = sX ** -2
	oY = sY ** -2

	if verbose:
		print(f'oX = {oX}')
		print(f'oY = {oY}')

	def new_estimate(b):
		W = oX * oY / (oX + (b ** 2) * oY - 2 * b * rXY * (oX * oY) ** .5)
		avgX = sum(W * X) / sum(W)
		avgY = sum(W * Y) / sum(W)
		U = X - avgX
		V = Y - avgY
		beta = W * (U / oY + b * V / oX - (b * U + V) * rXY / (oX * oY) ** .5)
		new_b = sum(W * beta * V) / sum(W * beta * U)
		new_a = avgY - new_b * avgX
		return new_a, new_b

	# initial estimate of slope B and intercept A
	B,A = np.polyfit(X, Y, 1)
	new_A, new_B = new_estimate(B)
	if rel_cvg:
		e = abs(new_B / B - 1)
	else:
		e = abs(new_B - B)

	Niter = 0
	while e >= cvg_limit and Niter < max_iter:
		Niter += 1
		if verbose:
			print(f"YorkReg -> [{Niter}] A = {new_A}, B = {new_B}, e = {e}")
		A, B = new_A, new_B
		new_A, new_B = new_estimate(B)
		if rel_cvg:
			e = abs(new_B / B - 1)
		else:
			e = abs(new_B - B)
	
	if verbose:
		print(f"YorkReg -> [{Niter+1}] A = {new_A}, B = {new_B}, e = {e}")
	A, B = new_A, new_B

	W = oX * oY / (oX + (B ** 2) * oY - 2 * B * rXY * (oX * oY) ** .5)
	avgX = sum(W * X) / sum(W)
	avgY = sum(W * Y) / sum(W)
	U = X - avgX
	V = Y - avgY
	beta = W * (U / oY + B * V / oX - (B * U + V) * rXY / (oX * oY) ** .5)
	x = avgX + beta
	y = avgY + B * beta
	avgx = sum(W * x) / sum(W)
	avgy = sum(W * y) / sum(W)
	u = x - avgx
	v = y - avgy

	sB2 = 1 / sum(W * (u ** 2))
	sA2 = 1 / sum(W) + sB2 * avgx ** 2
	cAB = -avgx * sB2

	chi2 = sum(W * (Y - B * X - A) ** 2)
	Nf = X.size - 2
	if Nf > 0:
		rchisq = chi2/Nf
	else:
		rchisq = 0.

	if rescale:
		if rchisq > 0:
			f = max(rchisq, 1)
			if verbose:
				print(f'RMSWD = {sqrt(rchisq):.2f}')
				print(f'Scaling variance by a factor of {f:.2f}')
			sA2 *= f
			sB2 *= f
			cAB *= f

	return B, A, np.array([[sB2,cAB],[cAB,sA2]]), rchisq


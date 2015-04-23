""" peak_detect.py
	
	Module for automatically finding peaks in 1D quasi-periodic signals. 
	
	Reference: F Scholkmann, J Boss and M Wolf. An Efficient Algorithm for Automatic Peak Detection in Noisy Periodic and Quasi-Periodic Signals, Algorithms 5, 588-603 (2012)
	
	Usage:
		> peakIdx = peak_detect.getPeaks(t,x)
		> # plot the peaks
		> plot(t[peakIdx], x[peakIdx], 'or', markerfacecolor='None', markeredgewidth=2, markersize=8, markeredgecolor='r')
	
	
	v 1.0
	Christopher Smeenk
	April 2015
	

"""
import pylab as py

def calcLMS(vec):
	
	L = py.ceil(vec.size/2) - 1
	alpha = 1.
	
	M = py.rand(L,vec.size) + alpha# * py.ones((L,vec.size))
	
	for k in range(1,int(L)):
		for i in range(int(k+2), int(vec.size-k+1)):
			
			wk = 2*k  # window size
			
			if (vec[i-1] > vec[i-k-1]) and (vec[i-1]>vec[i+k-1]):
				M[k,i] = 0.
	
	
	return M
	
def getPeaks(x,y):
	
	linfit = py.polyfit(x, y, 1)
	yfit = py.polyval(linfit, x)

	ynew = y - yfit

	M = calcLMS(ynew)
	gamma_k = M.sum(1)
	lambda_m = gamma_k.argmin()
	# rescale the LMS
	Mr = M[1:lambda_m, :]

	sigma = 1./(lambda_m-1) * py.sum( py.sqrt(( Mr - 1./lambda_m* (py.ones((lambda_m-1,1))* Mr.sum(0)) )**2), 0)

	p = py.find(sigma == 0) - 1
	
	return p


def zeropad(vec, n):
	newvec = py.hstack([py.zeros((n)), vec, py.zeros((n))])
	
	return newvec


# synthetic data
"""
a = 0.5
b = 0.5
c = 0.5
d = 0.5
fs = 20. #80
f1 = 10.
f2 = 70.
f3 = 5.
t = py.linspace(0,5,1000)  # (0,5,1000)
x = a*py.sin(2.*py.pi*f1/fs*t) + b*py.sin(2.*py.pi*f2/fs*t) + c*py.sin(2.*py.pi*f3/fs*t) + d*(py.rand(t.size)-0.5)


# Code to process the data. Copied into getPeaks() above.
#linfit = py.polyfit(t, x, 1)
#yfit = py.polyval(linfit, t)

#ynew = x - yfit

#M = calcLMS(ynew)
#gamma_k = M.sum(1)
#lambda_m = gamma_k.argmin()
## rescale the LMS
#Mr = M[1:lambda_m, :]

#sigma = 1./(lambda_m-1) * py.sum( py.sqrt(( Mr - 1./lambda_m* (py.ones((lambda_m-1,1))* Mr.sum(0)) )**2), 0)

#peakIdx = py.find(sigma == 0) - 1
"""

# run the code below for debugging & visualization
"""
peakIdx = getPeaks(t,x)

py.figure()
py.plot(t,x)
py.plot(t[peakIdx], x[peakIdx], 'or', markersize=8, markerfacecolor='None', markeredgecolor='r',markeredgewidth=2)
"""

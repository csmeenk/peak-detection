import pylab as py
import peak_detect

c,t,y = py.loadtxt('monthssn.dat', unpack=1)

peakIdx = peak_detect.getPeaks(t,y)
mu = py.mean(py.diff(t[peakIdx]))
s = py.std(py.diff(t[peakIdx]))


py.figure()
py.plot(t,y)
py.plot(t[peakIdx], y[peakIdx], 'or', markersize=8, markerfacecolor='None', markeredgecolor='r',markeredgewidth=2)
py.text(1750,270, 'Period: %.4g +/- %.4g years' % ( mu, s) )

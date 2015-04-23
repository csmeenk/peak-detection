linfit = py.polyfit(t, x, 1)
yfit = py.polyval(linfit, t)
ynew = x - yfit
py.plot(t,yfit)
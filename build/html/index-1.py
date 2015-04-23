import pylab as py
a = 1.
b = 1.
c = 1.5
d = 0.5
fs = 10. #80
f1 = 10.
f2 = 70.
f3 = 5.
t = py.linspace(0,5,1000)  # (0,5,1000)
x = a*py.sin(2.*py.pi*f1/fs*t) + b*py.sin(2.*py.pi*f2/fs*t) + c*py.sin(2.*py.pi*f3/fs*t) + d*(py.rand(t.size)-0.5)
py.plot(t,x)
import numpy as np
from matplotlib import pyplot as plt
import math
from pylab import *

xmax = 10
tmax = 20
xsteps = 100
tsteps = 100
u0 = np.zeros(xsteps)
AMatrix = np.zeros((xsteps, xsteps))
uarray = np.zeros((xsteps, tsteps))
deltax = xmax / xsteps
deltat = tmax / tsteps
c = 0.5
umax = 0
tarray = np.linspace(0,tmax,num=tsteps)
aa = 0.5 - c * deltat / (2 * deltax)
bb = 0.5 + c * deltat / (2 * deltax)

# t = 0 initialisieren
x0 = np.linspace(0, xmax, num=xsteps)
for i in range(0,xsteps):
	uarray[i][0] = math.exp(-10*(x0[i]-2)**2)


AMatrix[0][2] = aa
AMatrix[0][xsteps-1] = bb
AMatrix[xsteps-1][0] = aa
AMatrix[xsteps-1][xsteps-2] = bb

for i in range(1,xsteps-1):
	AMatrix[i][i+1] = aa
	AMatrix[i][i-1] = bb



for i in range(1,tsteps-1):
	uarray[:,i] = np.dot(AMatrix, uarray[:,i-1])


umax = np.amax(AMatrix)
fig = plt.figure()
ax = plt.axes(xlim=(0, xmax), ylim=(0, umax))
	
for ind in range(0, tsteps):
	
	ax.clear()
	ax.plot(x0, uarray[:,ind])
	figtext(.002, .02, "by Frank Ehebrecht / Januar 2014")
	plt.title('t = %f'%(tarray[ind]))
	savefig('anim%03d.png'%ind, dpi=100)
# rest mit ffmpeg manuell zu film (animation.py aus matplotlib ist beim
# speichern von animationen fehlerhaft)

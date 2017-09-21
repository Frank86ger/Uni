import numpy as np
from matplotlib import pyplot as plt
import math
from pylab import *

xmax = 10.

deltax = 0.05
deltat = 0.05
xsteps = int(xmax / deltax)
tsteps = 36
tmax = deltat * tsteps

u0 = np.zeros(xsteps)
AMatrix = np.zeros((xsteps, xsteps))
uarray = np.zeros((xsteps, tsteps))

umax = 0
tarray = np.linspace(0,tmax,num=tsteps)


# t = 0 initialisieren
x0 = np.linspace(0, xmax, num=xsteps)
for i in range(0,xsteps):
	uarray[i][0] = math.exp(-(x0[i]-3)**2)



for tt in range(0, tsteps-1):
	AMatrix[0][xsteps-2] = uarray[0][tt]*(deltat/deltax)
	AMatrix[xsteps-2][xsteps-3] = uarray[xsteps-1][tt]*(deltat/deltax)
	
	for i in range(0, xsteps):
		AMatrix[i][i] = 1-uarray[i][tt]*(deltat/deltax)
	for i in range(1,xsteps-1):
		AMatrix[i+1][i] = uarray[i][tt]*(deltat/deltax)
		
	uarray[:,tt+1] = np.dot(AMatrix, uarray[:,tt])



umax = np.amax(AMatrix)
fig = plt.figure()
#ax = plt.axes(xlim=(0, xmax), ylim=(0, umax))
ax = fig.add_subplot(1,1,1)

for ind in range(0, tsteps):
	
	ax.clear()
	ax.set_xlim([0,xmax])
	ax.set_ylim([0,1.2])
	ax.plot(x0, uarray[:,ind])
	figtext(.002, .02, "by Frank Ehebrecht / Januar 2014")
	plt.title('t = %f'%(tarray[ind]))
	savefig('anim%03d.png'%ind, dpi=100)
# rest mit ffmpeg manuell zu film (animation.py aus matplotlib ist beim
# speichern von animationen fehlerhaft)

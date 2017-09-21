import numpy as np
from matplotlib import pyplot as plt
import math
from pylab import *


deltax = 0.1
deltat = 0.05
xmax = 20
#tmax = 1800 #1800!!!
tsteps = 1800

c = 0.2

xsteps = 2*int(xmax / deltax)+1
#tsteps = int(tmax / deltat)
tmax = deltat * tsteps

nthstep = 5

#bisschen unschoen
uarray = np.zeros((xsteps, tsteps))
tarray = np.zeros(tsteps)
xarray = np.zeros(xsteps)
fx = np.zeros(xsteps)
gx = np.zeros(xsteps) #gamma1
beta = np.zeros(xsteps)
AMat = np.zeros((xsteps, xsteps))
alpha = deltat/deltax
tarray = np.zeros(tsteps)

def sechh(argu):
	return 2 / (math.exp(argu)+math.exp(-argu))




for i in range(0, tsteps):
	tarray[i] = i * deltat

for i in range(0, xsteps):
	xarray[i] = i*deltax - xmax#
	fx[i] = 4*math.atan(math.exp(-xarray[i]/math.sqrt(1-c**2)))
	gx[i] = -2*c/math.sqrt(1-c**2)*1/math.cosh(xarray[i]/math.sqrt(1-c**2))
	#gx[i] = -2*c/math.sqrt(1-c**2) * sechh(xarray[i]/math.sqrt(1-c**2))
	uarray[i][0] = fx[i]
	beta[i] = math.sin(fx[i])
	AMat[i][i] = 1-alpha**2

AMat[0][1] = alpha**2
AMat[xsteps-1][xsteps-2] = alpha**2

for i in range(1,xsteps-1):
	AMat[i][i-1] = alpha**2 / 2
	AMat[i][i+1] = alpha**2 / 2

#u1
uarray[:,1] = np.multiply(deltat, gx) +  np.dot(AMat, uarray[:,0]) - np.multiply(deltat**2/2, beta)


for j in range(1,tsteps-1):
	beta = np.sin(uarray[:,j])
	uarray[:,j+1] = -uarray[:,j-1] + np.dot(np.multiply(2, AMat), uarray[:,j] - np.multiply(deltat**2, beta))


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

for ind in range(0, tsteps+1):
	
	if ind%nthstep == 0:
		ax.clear()
		ax.set_xlim([-xmax,xmax])
		ax.set_ylim([-1,7])
		ax.plot(xarray, uarray[:,ind])
		figtext(.002, .02, "by Frank Ehebrecht / Januar 2014")
		plt.title('t = %f'%(tarray[ind]))
		savefig('anim%03d.png'%int(ind/nthstep), dpi=100)
# rest mit ffmpeg manuell zu film (animation.py aus matplotlib ist beim
# speichern von animationen fehlerhaft)

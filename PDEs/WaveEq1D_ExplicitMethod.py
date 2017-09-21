import numpy as np
from matplotlib import pyplot as plt
import math
from pylab import *


deltax = 0.1
deltat = 0.05
xmax = 10
tsteps = 20

xsteps = int(xmax / deltax)
tmax = tsteps*deltat
uarray = np.zeros((xsteps+2, tsteps+1))
tarray = np.zeros(tsteps+1)

for i in range(0, tsteps+1):
	tarray[i] = i*deltat-deltat
	# i=1 => t=0


alpha = 2*deltat/deltax


x = np.zeros(xsteps+2)
fx = np.zeros(xsteps+2)

for i in range(0, xsteps+2):
	x[i] = i*deltax-deltax
	fx[i] = math.sin(x[i]*math.pi)
	# x(t=0) = x0 = x[1]

#startwerte
for i in range(1, xsteps+1):
	uarray[i][1] = fx[i]
	uarray[i][2] = (1-alpha**2)*fx[i] + alpha**2/2*(fx[i-1] + fx[i+1])
	uarray[i][0] = uarray[i][2] 


for j in range(1,tsteps):
	for i in range(0, xsteps+1):
		uarray[i][j+1] = -uarray[i][j-1] + 2*(1-alpha**2)*uarray[i][j] + alpha**2*(uarray[i+1][j]+uarray[i-1][j])


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

for ind in range(0, tsteps+1):
	
	ax.clear()
	ax.set_xlim([0,xmax])
	ax.set_ylim([-1.2,1.2])
	ax.plot(x, uarray[:,ind])
	figtext(.002, .02, "by Frank Ehebrecht / Januar 2014")
	plt.title('t = %f'%(tarray[ind]))
	savefig('anim%03d.png'%ind, dpi=100)
# rest mit ffmpeg manuell zu film (animation.py aus matplotlib ist beim
# speichern von animationen fehlerhaft)

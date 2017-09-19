import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math
import RK4_nVar as Runge

Tmax = 1
delta_DGL = 0.001
DGL_steps = math.floor(Tmax/delta_DGL)
delta_abl = 0.001
counter = 0
AA = 0
BB = 0
CC = 0
DD = 0
dett = 0
F1 = 1
F2 = 1
r_step = 1
s_step = 1
eps = 0.001


# Sachen
while((math.fabs(F1)>eps) or (math.fabs(F2)>eps)):
	counter += 1
	#DGL 3 mal loesen
	# zeit, array der dynamischen Vars = RK4.nVar(Anzahl Gleichungen, [Gleichungen], [Startwerte], Schrittweite, max Zeit)
	t1, dyns = Runge.RK4_nVar(4, ['x[1]', 'x[2]', 'x[3]', '(1+t**2)*(x[2])**2-5*(x[0])**2'], [1, 0, s_step, r_step], delta_DGL, Tmax)
	ts, dyns_s = Runge.RK4_nVar(4, ['x[1]', 'x[2]', 'x[3]', '(1+t**2)*(x[2])**2-5*(x[0])**2'], [1, 0, s_step+delta_abl, r_step], delta_DGL, Tmax)
	tr, dyns_r = Runge.RK4_nVar(4, ['x[1]', 'x[2]', 'x[3]', '(1+t**2)*(x[2])**2-5*(x[0])**2'], [1, 0, s_step, r_step+delta_abl], delta_DGL, Tmax)
	
	#Matrixelemente der Jacobi Matrix / Determinante
	AA = ((dyns_s[2][DGL_steps-1]+2) - (dyns[2][DGL_steps-1]+2)) / delta_abl
	BB = ((dyns_r[2][DGL_steps-1]+2) - (dyns[2][DGL_steps-1]+2)) / delta_abl
	CC = ((dyns_s[3][DGL_steps-1]+3) - (dyns[3][DGL_steps-1]+3)) / delta_abl
	DD = ((dyns_r[3][DGL_steps-1]+3) - (dyns[3][DGL_steps-1]+3)) / delta_abl
	dett = AA*DD-BB*CC
		
	# Naechste Schritte berechnen
	F1 = (DD*(dyns[2][DGL_steps-1]+2)-BB*(dyns[3][DGL_steps-1]+3)) / dett
	F2 = (AA*(dyns[3][DGL_steps-1]+3)-CC*(dyns[2][DGL_steps-1]+2)) / dett
	s_step -= F1
	r_step -= F2


# PLOTS
fig = plt.figure(figsize=plt.figaspect(2.))
ax1 = fig.add_subplot(1,1,1)
ax1.set_title('Abbruch nach %i Iterationen fuer eps=%f ! \n x\'\'(0)= %f und x\'\'\'(0)= %f'%(counter, eps, dyns[2][0], dyns[3][0]))
ax1.set_xlabel('Zeit t')
ax1.set_ylabel('x(t)')
ax1.plot(t1, dyns[0], 'bo')
plt.show()

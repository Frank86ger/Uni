import numpy as np
#import scipy as sc
import matplotlib.pyplot as plt
import math


# Startwerte
T = 20 * math.pi
h = 0.05
w = 2
steps = math.floor(T/h) #Anzahl Schritte
x_0 = 0
v_0 = 1
t_0 = 0


ar_x1 = np.zeros(steps) # Arrays erzeugen
ar_y1 = np.zeros(steps)
ar_t1 = np.zeros(steps)
ar_x1[0] = x_0 # Startbedingungen
ar_y1[0] = v_0
ar_t1[0] = t_0

for k in range(1,steps-1): # Heun Verfahren
	ar_t1[k] = ar_t1[0] + k * h # Zeitpunkte
	
	xtemp = ar_x1[k-1] + h * ar_y1[k-1] 	   # Hilfsschritte
	ytemp = ar_y1[k-1] - h * w**2 * ar_x1[k-1] # 
	
	ar_x1[k] = 0.5 * ar_x1[k-1] + 0.5*(xtemp + h * ytemp) # Hauptschritte
	ar_y1[k] = 0.5 * ar_y1[k-1] + 0.5*(ytemp + h * xtemp)
	

plt.plot(ar_x1, ar_y1, 'k+') # Plot
plt.show()

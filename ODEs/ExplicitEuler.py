import numpy as np
#import scipy as sc
import matplotlib.pyplot as plt
import math

# Rechte Seite der DGL
def f(x,t):
	return t - x**2

# Startwerte
T = 9
h = 0.05
x_0 = [-0.7, 0.0, 1.0, 3.0]
farbe = ['bo', 'ro', 'go', 'ko']
steps = math.floor(T/h)

# Werte berechnen
for i in range(0,4):
	ar_x1 = np.zeros(steps)
	ar_t1 = np.zeros(steps)
	ar_x1[0] = x_0[i] # Startwert
	ar_t1[0] = 0 # t_0 = 0
	for k in range(1,steps-1): # Euler Verfahren
		ar_t1[k] = ar_t1[0] + k * h
		ar_x1[k] = ar_x1[k-1] + h * f(ar_x1[k-1], ar_t1[k-1])
	plt.plot(ar_t1, ar_x1, farbe[i])

# Achsen einstellen / Plot ausgeben
plt.axis([-0.5, 9.5, -1.5, 3.5])
plt.show()

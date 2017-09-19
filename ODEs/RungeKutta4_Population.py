import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math


h = 0.1			# Schrittweite
t_0 = [0, 0, 0, 0]		# t-Startwerte
x_0 = [50, 100, 100, 120]	# x-Startwerte
y_0 = [100, 50, 100, 120]	# y-Startwerte
T = 250				# Zeitintervall [0,T]
steps = math.floor(T/h) 	# Anzahl Schritte


# Parameter
a = 0.004
b= 50
c = 0.75
d = 0.001
e = 100
f = 3.0


# Funktionen fuer RK4
def f_x(x, y):
	return a*x*(b-x-c*y)
def f_y(x,y):
	return d*y*(e-y-f*x)
def RK4_xstep(xn, yn):
	k = f_x(xn, yn)
	kges = k
	print(kges)
	k = f_x(xn+h/2*k, yn+h/2*k)
	kges = kges + 2*k
	print(kges)
	k = f_x(xn+h/2*k, yn+h/2*k)
	kges = kges + 2*k
	print(kges)
	k = f_x(xn+h*k, yn+h*k)
	kges = kges + k
	print(kges)
	return xn + h/6 * kges
def RK4_ystep(xn, yn):
	k = f_y(xn, yn)
	kges = k
	k = f_y(xn+h/2*k, yn+h/2*k)
	kges = kges + 2*k
	k = f_y(xn+h/2*k, yn+h/2*k)
	kges = kges + 2*k
	k = f_y(xn+h*k, yn+h*k)
	kges = kges + k
	return yn + h/6 * kges


# Arrays initialisieren
ar_t = np.zeros((steps, 4))
ar_x = np.zeros((steps, 4))
ar_y = np.zeros((steps, 4))


# For loop fuer die verschiedenen Startwerte
for i in range(0,4):
	# Array Startwerte setzen
	ar_t[0][i] = t_0[i]
	ar_x[0][i] = x_0[i]
	ar_y[0][i] = y_0[i]
	
	# RK4 Werte in Arrays fuellen / berechnen
	for k in range(1,steps):
		ar_t[k][i] = ar_t[0][i] + k*h
		ar_x[k][i] = RK4_xstep(ar_x[k-1][i], ar_y[k-1][i])
		ar_y[k][i] = RK4_ystep(ar_x[k-1][i], ar_y[k-1][i])

# Plot-Layout
plotit = [['ax1a', 'ax1b'],['ax2a', 'ax2b'],['ax3a', 'ax3b'],['ax4a', 'ax4b']]
fig, axes = plt.subplots(2,4)
((plotit[0][0], plotit[1][0], plotit[2][0], plotit[3][0]), (plotit[0][1], plotit[1][1], plotit[2][1], plotit[3][1])) = axes # unpack the axes

# Plotten
for i in range(0,4):
	plotit[i][0].set_title(x_0[i]) ### Startwerte. Jedoch Problem: Geht leider nicht wie bei --> print('x0=', x_0[i], ', y0=', y0[i]) <--
	plotit[i][0].set_xlabel('Zeit t')
	plotit[i][0].set_ylabel('Population x / y')
	plotit[i][0].plot(ar_t[:,i], ar_x[:,i], 'b-')
	plotit[i][0].plot(ar_t[:,i], ar_y[:,i], 'r-')
	plotit[i][0].legend(["x(t)", "y(t)"])
	
	plotit[i][1].set_xlabel('Population x')
	plotit[i][1].set_ylabel('Population y')
	plotit[i][1].plot(ar_x[:,i], ar_y[:,i], 'k-')
	plotit[i][1].legend(["y(x)"])


plt.figtext(0.2, 0.965, 'ab>de: x profitiert staerker durch seine Groesse als y im linearen Anteil.', ha='left')
plt.figtext(0.2, 0.945, 'a>d: x wird jedoch auch staerker durch seine Groesse am Wachstum gehindert als y (quadratischer Anteil)', ha='left')
plt.figtext(0.2, 0.925, 'ac=df: Beide Gruppen sind gleichermassen agressiv. ==> Wenn x=y, ist x bei kleinerer Population besser (y bei Grosser)', ha='left')
plt.show()

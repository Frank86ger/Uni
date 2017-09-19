import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math

T = 50			# Maximale Zeit
h0 = 0.01		# Start-Schrittweite
steps = 900000		# Anzahl maximaler Schritte (so gross waehlen, dass genug Platz fuer Werte vorhanden ist)
var_count = 2		# Anzahl der dynamischen Variablen
eps_tol = 1.e-10	# Epsilon-Toleranz
b_fak = 0.9		# Sicherheitsfaktor

# Startwerte der dynamischen Variablen
dyn_var_start = [5000, 120]

# Arrays fuer Werte der dynamischen Variablen
dyn_var = np.zeros((var_count, steps))
t_array = np.zeros(steps)

# Startwerte einsetzen
for i in range(0, var_count):
	dyn_var[i][0] = dyn_var_start[i]
t_array[0] = 0

# Butcher-Parameter (Fehlberg)
alpha = np.array([0, 1/4, 3/8, 12/13, 1, 1/2])
beta = np.array([[0, 0, 0, 0, 0, 0],
	[1/4, 0, 0, 0, 0, 0],
	[3/32, 9/32, 0, 0, 0, 0],
	[1932/2197, -7200/2197, 7296/2197, 0, 0, 0],
	[439/216, -8, 3680/513, -845/4104, 0, 0],
	[-8/27, 2, -3544/2565, 1859/4104, -11/40, 0]])
c = np.array([25/216, 0, 1408/2565, 2197/4104, -1/5, 0])
c_star = np.array([16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55])
# ..............................

# k-Vektoren
k_vec = np.zeros((var_count, 6))

# Epsilon-Hilfswerte
eps_vec = np.zeros(var_count)
eps_min = 0
eps_max = 0

# Vektor der h-Werte
h_vec = np.zeros(steps)
h_vec[0] = h0

# EpsSum zum Abschaetzen des Fehlers
eps_sum = np.zeros(steps)
eps_sum[0] = 0

# Funktionen f. t skalar; x vektoriell
def dyn_var_func1(t, x):
	return 2*x[0]-0.02*x[0]*x[1]
def dyn_var_func2(t, x):
	return 0.0002*x[0]*x[1]-0.8*x[1]
# Array obiger Funktionen
dyn_var_funcs = [dyn_var_func1, dyn_var_func2] #[dyn_var_func1, dyn_var_func2, ..., ]

# k-Vektor-Array berechnen (DAS IST EINE "var_count X 6"-Matrix
def k_calc(element):
	for i1 in range(0, var_count): #k-Vektor fuer jede dynamische Variable
		for i2 in range(0, 6): #6 Elemente des k-Vektors berechnen
			butch_fak = np.dot(beta[i2,:],k_vec[i1,:]) # k-Summe berechnen
			k_vec[i1][i2] = dyn_var_funcs[i1](t_array[element-1]+h_vec[element-1]*alpha[i2],
				dyn_var[:,element-1]+np.inner(h_vec[element-1] * butch_fak, np.ones(var_count)))
	# Jetzt sind alle k-Matrix eintraege berechnet -> nun: Schrittweite testen
	hnew_calc(element)

# Schrittweite testen und einstellen
def hnew_calc(element):
	for i1 in range(0, var_count):
		eps_temp = np.dot((c - c_star), k_vec[i1,:])
		eps_vec[i1] = math.fabs(eps_temp) #Array aller eps-Werte

	eps_max = np.amax(eps_vec)
	if eps_max >= eps_tol: #Falls eps zu gross
		h_vec[element-1] *= b_fak*(eps_tol/eps_max)**0.2
		k_calc(element) # Dies ist eine Rekursion!
	else: # Schrittweite zulaessig - danach neue Schrittweite
		eps_min = 0		
		for k in range(0, var_count): # groesstes eps, welches kleiner als eps_tol ist nehmen
			if (eps_vec[k] <= eps_tol) and (eps_vec[k] >= eps_min):
				eps_min = eps_vec[k]
		calc_next_step(element, eps_min) #naechsten Schritt berechnen
		h_vec[element] = h_vec[element-1] * b_fak * (eps_tol/eps_min)**0.25

# Neue Elemente berechnen
def calc_next_step(element, epsmin):	
	for i1 in range(0, var_count):
		dyn_var[i1][element] = dyn_var[i1][element-1] + h_vec[element-1] * np.dot(c_star, k_vec[i1,:])
		t_array[element] = t_array[element-1] + h_vec[element-1]
		eps_sum[element] = eps_sum[element-1] + epsmin


# !!! START !!!
counter = 0
while (t_array[counter-1] <= T):
	counter += 1
	k_calc(counter)
# !!!!!!!!!!!!!


# Bei dyn_var-arrays Nullen entfernen
dyn_var_trimmed = np.zeros((var_count, counter))
for i in range(0, var_count):
	for j in range(0, counter):
		dyn_var_trimmed[i][j] = dyn_var[i][j]

# Bei t-array Nullen entfernen
t_array_trimmed = np.zeros(counter)
for i in range(0, counter):
	t_array_trimmed[i] = t_array[i]
	
# Bei eps_sum Nullen entfernen
eps_sum_trimmed = np.zeros(counter)
for i in range(0, counter):
	eps_sum_trimmed[i] = eps_sum[i]


# PLOTS
fig, axes = plt.subplots(1,2)
(ax1, ax2) = axes

ax1.set_title('x0 = 5000 / y0 = 120')
ax1.set_xlabel('Zeit t')
ax1.set_ylabel('Population x(t) / y(t)')
ax1.plot(t_array_trimmed, dyn_var_trimmed[0], 'bo')
ax1.plot(t_array_trimmed, dyn_var_trimmed[0]+eps_sum_trimmed, 'r-') #+fehler
ax1.plot(t_array_trimmed, dyn_var_trimmed[0]-eps_sum_trimmed, 'y-') #-fehler
ax1.plot(t_array_trimmed, dyn_var_trimmed[1], 'ro')
ax1.legend(["x(t)", "x+Fehler", "x-Fehler", "y(t)"])

ax2.set_title('')
ax2.set_xlabel('Population x(t)')
ax2.set_ylabel('Population y(t)')
ax2.plot(dyn_var_trimmed[0], dyn_var_trimmed[1], 'b-')
ax2.legend(["y(x)"])
plt.show()

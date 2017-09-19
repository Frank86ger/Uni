import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


def function_builder(arg):
	def function(t, x):
		return eval(arg)
	return function
my_dynamic_function = {}

def RK4_nVar(Gleichung):
	T = 100			# Maximale Zeit
	h0 = 0.01		# Start-Schrittweite
	steps = math.floor(T/h0)# Anzahl maximaler Schritte (so gross waehlen, dass genug Platz fuer Werte vorhanden ist)
	var_count = 3		# Anzahl der dynamischen Variablen
	arg=0
	
	# Startwerte der dynamischen Variablen
	dyn_var_start = [10, 20, 30]
	
	# Arrays fuer Werte der dynamischen Variablen
	dyn_var = np.zeros((var_count, steps))
	t_array = np.zeros(steps)
	
	# Startwerte einsetzen
	for i in range(0, var_count):
		dyn_var[i][0] = dyn_var_start[i]
	t_array[0] = 0
	
	# Butcher-Parameter (RK4)
	alpha = np.array([0, 1/2, 1/2, 1])
	beta = np.array([[0, 0, 0, 0],
			[1/2, 0, 0, 0],
			[0, 1/2, 0, 0],
			[0, 0, 1, 0]])
	cvec = np.array([1/6, 1/3, 1/3, 1/6])
	# ..............................
	
	# k-Vektoren
	k_vec = np.zeros((var_count, 4))
	

	
	
	########################################
	# x-Array global definieren?
	
	
	for i in range(0, var_count):
		# argument fue die funktionen in gruen
		my_dynamic_function[i] = function_builder(Gleichung[i])
	
	#######################################
	
	
	
	
	
	# k-Vektor-Array berechnen (DAS IST EINE "var_count X 6"-Matrix
	def k_calc(element):
		for i1 in range(0, var_count): #k-Vektor fuer jede dynamische Variable
			for i2 in range(0, 4): #6 Elemente des k-Vektors berechnen
				butch_fak = np.dot(beta[i2,:],k_vec[i1,:]) # k-Summe berechnen
				k_vec[i1][i2] = my_dynamic_function[i1](t_array[element-1]+h0*alpha[i2],
					dyn_var[:,element-1]+np.inner(h0 * butch_fak, np.ones(var_count)))
		# Jetzt sind alle k-Matrix eintraege berechnet -> nun: Schrittweite testen
		calc_next_step(element)
	
	
	# Neue Elemente berechnen
	def calc_next_step(element):	
		for i1 in range(0, var_count):
			dyn_var[i1][element] = dyn_var[i1][element-1] + h0 * np.dot(cvec, k_vec[i1,:])
			t_array[element] = t_array[element-1] + h0
	
	
	# !!! START !!!
	counter = 0
	while (counter < steps-1):
		counter += 1
		k_calc(counter)
	# !!!!!!!!!!!!!
	
	return t_array, dyn_var
#### ENDE  DER  FUNKTION #########################



#func_list = ['10*(x[1]-x[0])', '25*x[0]-x[1]-x[0]*x[2]', 'x[0]*x[1]-8/3*x[2]']

ttt, dyns = RK4_nVar(['10*(x[1]-x[0])', '25*x[0]-x[1]-x[0]*x[2]', 'x[0]*x[1]-8/3*x[2]'])



# PLOTS
fig = plt.figure(figsize=plt.figaspect(2.))

ax1 = fig.add_subplot(2,2,1)
ax1.set_title('r = 25')
ax1.set_xlabel('Zeit t')
ax1.set_ylabel('Winkel phi(t)')
ax1.plot(ttt, dyns[2], 'bo')
ax1.legend('z(t)')

ax2 = fig.add_subplot(2, 2, 2, projection='3d')
ax2.plot(dyns[0], dyns[1], dyns[2], 'b-')




plt.show()

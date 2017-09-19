import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import math


tmin = -1.
tmax = 1.
h = 0.02
N_steps = math.floor((tmax-tmin)/h)+1

AMatrix = np.zeros((N_steps, N_steps))
b_vec = np.zeros(N_steps)
t_array = np.zeros(N_steps)
q_array = np.zeros(N_steps)

#t array berechnen
for i in range(0,N_steps):
	t_array[i] = tmin + i*h

# q-array berechnen
for i in range(0, N_steps):
	q_array[i]= -(1+(tmin + i*h)**2)

# ** A Matrix Elemente **
AMatrix[0][0] = -(2+h**2*q_array[0])
AMatrix[0][1] = 1
AMatrix[N_steps-1][N_steps-2] = 1 
AMatrix[N_steps-1][N_steps-1] = -(2+h**2*q_array[N_steps-1]) 

for i in range(1,N_steps-1):
	AMatrix[i][i-1] = 1
	AMatrix[i][i] = -(2+h**2*q_array[i])
	AMatrix[i][i+1] = 1
# ** **

# Loesnungvektor
b_vec[0] = -1*h**2
b_vec[N_steps-1] = -1*h**2
for i in range(1,N_steps-1):
	b_vec[i] = -1*h**2

# System loesen
x = np.linalg.solve(AMatrix, b_vec)
#print(AMatrix)


#plotten
fig = plt.figure(figsize=plt.figaspect(2.))
ax1 = fig.add_subplot(1,1,1)

ax1.set_xlabel('Zeit t')
ax1.set_ylabel('x(t)')
ax1.plot(t_array, x, 'bo')
plt.show()

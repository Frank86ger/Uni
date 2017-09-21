#/////////////////////////////////////
#  by Frank Ehebrecht im Mai 2014  ///
#/////////////////////////////////////
import os
import ctypes
import numpy as np
from ctypes import CDLL, c_double, c_int # C/++konforme Datentypen / Einbinden von C/++ Bibs
import math
import Image

# C/++ konforme Datentypen definieren
STEPS = 256
SSTEPS = STEPS*STEPS
my_C_double_type = c_double * SSTEPS
my_C_int_type = c_int * SSTEPS #neu

my_C_input = my_C_double_type()
my_C_int_input = my_C_int_type()




#---Startwerte----------------------------
#im=Image.open("128_circle.png")	# Gebiets-Geometrie
im=Image.open("256_circle.png")


EPSILON = c_double(0.3)		# epsilon-Wert
DELTA = c_double(0.0)		# delta Wert
interval = 64.0			# Intervallgroesse / quadratisch
INTERVAL = c_double(interval)	

DAMPENING = c_double(0.9)	# 3*EPSILON / Daempfung ausserhalb des Gebiets

anispeed_c = c_double(10)	# Animations-Speed (sleep zwischen frames in us)
deltaT_c = c_double(0.01)	# Zeit-Schrittweite
timesteps = c_int(16010)		# Anzahl Zeitschritte
rausch_amp = 0.3		# Rauschamplitude fuer Startbedingung

eachnth = c_int(2000)		# jedes wievielte Bild fuer die Ausgabe genutzt wird

ZMIN = c_double(-1.05)		# Z range von ZMIN
ZMAX = c_double(1.05)		# bis ZMAX

method = c_int(1)		# 1: Immer 1
plottype = c_int(2)		# Immer 2
#-----------------------------------------




temp = STEPS / interval	#Ortsaufloesung
#-----------------------------------
#Werte des Eingangsignals definieren:
#-----------------------------------
for i in range(0, SSTEPS):

	#Reines Rauschen als Startbedingung
	my_C_input[i] = (np.random.random_sample()/1 - 0.5)*2*rausch_amp


	#ZigZag Instabilitaet
	"""
	my_C_input[i] = 0.578877 * math.sin((int(i / STEPS) * 0.88357 / temp)+0)
	my_C_input[i] = my_C_input[i] + (np.random.random_sample()/1 - 0.5)*2*rausch_amp
	"""


	#Eckhaus Instabilitaet
	"""
	my_C_input[i] = .578877 * math.sin(int(i / STEPS) * 1.178 / temp)
	my_C_input[i] = my_C_input[i] + (np.random.random_sample()/1 - 0.5)*2*rausch_amp
	"""
#-----------------------------------



# Pic-Werte to C-array (Gebiet)
pix = im.load()
for i in range(0, STEPS):
		for j in range(0, STEPS):
			if(pix[i, STEPS-j-1][0] == 0):
				my_C_int_input[i*STEPS+1*j] = 0
			else:
				my_C_int_input[i*STEPS+1*j] = 1


#Bibliothek aufrufen -> Dislin wird gestartet
ctypes.cdll.LoadLibrary(os.getcwd() + "/my_SwiftHo_lib.so").my_SwiftHo_solver(my_C_input, my_C_int_input, EPSILON, deltaT_c, timesteps, anispeed_c, eachnth, ZMIN, ZMAX, method, plottype, INTERVAL, DELTA, DAMPENING)

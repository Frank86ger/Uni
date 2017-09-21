#//////////////////////////////////////
#  by Frank Ehebrecht im Juli 2014  ///
#//////////////////////////////////////
import os
import ctypes
import numpy as np
from ctypes import CDLL, c_double, c_int # C/++konforme Datentypen / Einbinden von C/++ Bibs
import math
import Image
import subprocess

# C/++ konforme Datentypen definieren
STEPS = 256
SSTEPS = STEPS*STEPS
my_C_double_type = c_double * SSTEPS
my_C_int_type = c_int * SSTEPS #neu

my_C_input = my_C_double_type()
my_C_input_cplx = my_C_double_type()
my_C_int_input = my_C_int_type()




#---Startwerte----------------------------
im=Image.open("256_tp.png")	# Gebietsgeometrie

ALPHA = c_double(0.0)		# alpha-Wert
BETA = c_double(1.5)		# beta-Wert	
INTERVAL = c_double(200.0)	# Intervallgroesse / quadratisch
I_OFFSET = c_double(-100.0)	# Offset der Skalierung von Null (fuer x und y identisch)

DAMPENING = c_double(0.9)	# Daempfung ausserhalb des Gebiets

DeltaT = c_double(0.05)		# Zeit-Schrittweite
TIMESTEPS = c_int(1000)		# Anzahl Zeitschritte
rausch_amp = 0.01		# Rauschamplitude fuer Startbedingung

EACHNTH = c_int(10)		# jedes wievielte Bild fuer die Ausgabe genutzt wird

vid_name = "vidout"		# der Name des Output-Videos (ueberschreibt vorhandene videos NICHT)
#-----------------------------------------



temp = STEPS / INTERVAL.value	#Ortsaufloesung

#-----------------------------------
#Werte des Eingangsignals definieren:
#-----------------------------------
for i in range(0, SSTEPS):
	
	# Imaginaerteil
	my_C_input_cplx[i] = 0#(np.random.random_sample()/1 - 0.5)*2*rausch_amp

	#Reines Rauschen als Startbedingung
#	my_C_input[i] = (np.random.random_sample()/1 - 0.5)*2*rausch_amp

	# 1 + Rauschen als Startbedingung
	my_C_input[i] = 1 + (np.random.random_sample()/1 - 0.5)*2*rausch_amp
	
	
	#xx = float(int(i / STEPS)) / temp
	#yy = (i % STEPS) / temp
	
#-----------------------------------

print "Calculating ... Start"
# Pic-Werte to C-array (Gebiet)
pix = im.load()
for i in range(0, STEPS):
		for j in range(0, STEPS):
			if(pix[i, STEPS-j-1][0] == 0):
				my_C_int_input[i*STEPS+1*j] = 0
			else:
				my_C_int_input[i*STEPS+1*j] = 1


#Bibliothek aufrufen -> Dislin wird gestartet
ctypes.cdll.LoadLibrary(os.getcwd() + "/my_GLG_txtout.so").my_GLG_solver(my_C_input, my_C_input_cplx, my_C_int_input, ALPHA, BETA, TIMESTEPS, DeltaT, EACHNTH, INTERVAL, DAMPENING)
print "Calculating ... Finished"

SCALE = INTERVAL.value / float(STEPS)
IMAX = INTERVAL.value + I_OFFSET.value
MAXSTEPS = (TIMESTEPS.value / EACHNTH.value) - 1

print "Creating images ... Start"
subprocess.call(['gnuplot -e \"MAXSTEPS=\'' + str(MAXSTEPS) + '\'; SCALE=\'' + str(SCALE) + '\'; IMAX=\'' + str(IMAX) + '\'; I_OFFSET=\'' + str(I_OFFSET.value) + '\'\" ' + os.getcwd() + '/plotit.gplt', '-1'], shell=True) # gnuplot Script ausfuehren

subprocess.call(['rm ' + os.getcwd() + '/GLG*', '-1'], shell=True) # Daten loeschen
print "Creating images ... Finished"
print "Creating video ... Start"
subprocess.call(['ffmpeg -r 25 -qscale 2 -i pic%04d.png ' + vid_name + '.avi', '-1'], shell=True) # Video erstellen
subprocess.call(['rm ' + os.getcwd() + '/pic*', '-1'], shell=True) # Bilder loeschen
print "Creating video ... Finished"

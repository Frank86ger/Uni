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

my_C_inputU = my_C_double_type()
my_C_inputV = my_C_double_type()
my_C_int_input = my_C_int_type()




#---Startwerte----------------------------
im=Image.open("256_all.png")	# Gebietsgeometrie

AA = c_double(3.0)		# a-Wert
BB = c_double(18.0)		# b-Wert
DU = c_double(5.0)		# DU-Wert
DV = c_double(12.0)		# DV-Wert

INTERVAL = c_double(120.0)	# Intervallgroesse / quadratisch
I_OFFSET = c_double(0.0)	# Offset der Skalierung von Null (fuer x und y identisch)

DAMPENING = c_double(0.9)	# Daempfung ausserhalb des Gebiets

DeltaT = c_double(0.05)		# Zeit-Schrittweite
TIMESTEPS = c_int(1500)		# Anzahl Zeitschritte
rausch_amp = 0.01		# Rauschamplitude fuer Startbedingung

EACHNTH = c_int(6)		# jedes wievielte Bild fuer die Ausgabe genutzt wird

vid_name = "A1-a3b18"		# der Name des Output-Videos (ueberschreibt vorhandene videos NICHT)
#-----------------------------------------



temp = STEPS / INTERVAL.value	#Ortsaufloesung

#-----------------------------------
#Werte des Eingangsignals definieren:
#-----------------------------------
for i in range(0, SSTEPS):
	
	# U = A + Rauschen als Startbedingung
	my_C_inputU[i] = AA.value + (np.random.random_sample()/1 - 0.5)*2*rausch_amp	
		
	# V = B/A + Rauschen als Startbedingung
	my_C_inputV[i] = BB.value/AA.value + (np.random.random_sample()/1 - 0.5)*2*rausch_amp	
	
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
ctypes.cdll.LoadLibrary(os.getcwd() + "/my_bruss2D_txtout.so").my_Bruss2D_solver(my_C_inputU, my_C_inputV, my_C_int_input, AA, BB, DU, DV, INTERVAL, DAMPENING, DeltaT, EACHNTH, TIMESTEPS)


print "Calculating ... Finished"
SCALE = INTERVAL.value / float(STEPS)
IMAX = INTERVAL.value + I_OFFSET.value
MAXSTEPS = (TIMESTEPS.value / EACHNTH.value) - 1

print "Creating images ... Start"
subprocess.call(['gnuplot -e \"MAXSTEPS=\'' + str(MAXSTEPS) + '\'; SCALE=\'' + str(SCALE) + '\'; IMAX=\'' + str(IMAX) + '\'; I_OFFSET=\'' + str(I_OFFSET.value) + '\'; PICNAME=\'' + str(vid_name) + '\'\" ' + os.getcwd() + '/plotit.gplt', '-1'], shell=True) # gnuplot Script ausfuehren

subprocess.call(['rm ' + os.getcwd() + '/Bruss*', '-1'], shell=True) # Daten loeschen
print "Creating images ... Finished"
print "Creating video ... Start"
subprocess.call(['ffmpeg -r 25 -qscale 2 -i pic%04d.png ' + vid_name + '.avi', '-1'], shell=True) # Video erstellen
subprocess.call(['rm ' + os.getcwd() + '/pic*', '-1'], shell=True) # Bilder loeschen
print "Creating video ... Finished"

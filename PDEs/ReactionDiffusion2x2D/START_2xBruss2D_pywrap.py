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

my_C_inputU1 = my_C_double_type()
my_C_inputU2 = my_C_double_type()
my_C_inputV1 = my_C_double_type()
my_C_inputV2 = my_C_double_type()
my_C_int_input = my_C_int_type()


#---Startwerte----------------------------
im=Image.open("256_all.png")	# Gebietsgeometrie

AA = c_double(3.0)		# a-Wert
BB = c_double(9.9)		# b-Wert
ALPHA = c_double(1.0)		# alpha-Wert
DU1 = c_double(8.33)		# DU1-Wert
DV1 = c_double(8.33)		# DV1-Wert
DU2 = c_double(46.0)		# DU2-Wert
DV2 = c_double(120.0)		# DV2-Wert

UH1 = 3.0			# uh1-Wert
UH2 = 3.0			# uh2-Wert
VH1 = 3.3			# vh1-Wert
VH2 = 3.3			# vh2-Wert

INTERVAL = c_double(200.0)	# Intervallgroesse / quadratisch
I_OFFSET = c_double(0.0)	# Offset der Skalierung von Null (fuer x und y identisch)

DAMPENING = c_double(0.9)	# Daempfung ausserhalb des Gebiets

DeltaT = c_double(0.05)		# Zeit-Schrittweite
TIMESTEPS = c_int(1000)		# Anzahl Zeitschritte
rausch_amp = 0.01		# Rauschamplitude fuer Startbedingung

EACHNTH = c_int(10)		# jedes wievielte Bild fuer die Ausgabe genutzt wird

vid_name = "A2-e"		# der Name des Output-Videos (ueberschreibt vorhandene videos NICHT)
#-----------------------------------------



temp = STEPS / INTERVAL.value	#Ortsaufloesung

#-----------------------------------
#Werte des Eingangsignals definieren:
#-----------------------------------
for i in range(0, SSTEPS):
	
	# U1 = UH1 + Rauschen als Startbedingung
	my_C_inputU1[i] = UH1 + (np.random.random_sample()/1 - 0.5)*2*rausch_amp	
	
	# U2 = UH2 + Rauschen als Startbedingung
	my_C_inputU2[i] = UH2 + (np.random.random_sample()/1 - 0.5)*2*rausch_amp	

	# U2 = UH2 + Rauschen als Startbedingung
	my_C_inputV1[i] = VH1 + (np.random.random_sample()/1 - 0.5)*2*rausch_amp
	
	# U2 = UH2 + Rauschen als Startbedingung
	my_C_inputV2[i] = VH2 + (np.random.random_sample()/1 - 0.5)*2*rausch_amp		
	
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
ctypes.cdll.LoadLibrary(os.getcwd() + "/my_2xBruss2D_txtout.so").my_2xBruss2D_solver(my_C_inputU1, my_C_inputU2, my_C_inputV1, my_C_inputV2, my_C_int_input, AA, BB, ALPHA, DU1, DV1, DU2, DV2, INTERVAL, DAMPENING, DeltaT, EACHNTH, TIMESTEPS)


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

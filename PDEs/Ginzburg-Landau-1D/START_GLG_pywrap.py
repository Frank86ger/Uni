#////////////////////////////////////////////////////////////
#  by Frank Ehebrecht im Juli 2014  			  ///
#  (bei Bedarf noch cbrange in plotit_pic.gplt einstellen ///
#////////////////////////////////////////////////////////////
import os
import ctypes
import numpy as np
from ctypes import CDLL, c_double, c_int # C/++konforme Datentypen / Einbinden von C/++ Bibs
import math
import Image
import subprocess

# C/++ konforme Datentypen definieren
STEPS = 512	# Anzahl Stuetzstellen
my_C_double_type = c_double * STEPS

my_C_input = my_C_double_type()
my_C_input_cplx = my_C_double_type()


#---Startwerte----------------------------
ALPHA = c_double(0.0)		# alpha-Wert
BETA = c_double(-2)		# beta-Wert	
INTERVAL = c_double(200.0)	# Intervallgroesse
I_OFFSET = c_double(-100.0)	# Offset der Skalierung von Null


DeltaT = c_double(0.05)		# Zeit-Schrittweite
TIMESTEPS = c_int(1000)		# Anzahl Zeitschritte
rausch_amp = 0.01		# Rauschamplitude fuer Startbedingung

EACHNTH = c_int(5)		# jedes wievielte Bild fuer die Ausgabe genutzt wird

vid_name = "A1-e2"		# der Name des Output-Videos und Bildes (ueberschreibt vorhandene videos NICHT)
#-----------------------------------------


temp = STEPS / INTERVAL.value	#Ortsaufloesung
#-----------------------------------
#Werte des Eingangsignals definieren:
#-----------------------------------
for i in range(0, STEPS):
	
	# Imaginaerteil
	my_C_input_cplx[i] = 0 #(np.random.random_sample()/1 - 0.5)*2*rausch_amp


	#Reines Rauschen als Startbedingung
#	my_C_input[i] = (np.random.random_sample()/1 - 0.5)*2*rausch_amp


	# 1 + Rauschen als Startbedingung
#	my_C_input[i] = 1.0 + (np.random.random_sample()/1 - 0.5)*2*rausch_amp
	
	
	# Benjamin Feir Instabilitaet
	"""
	x = float(i) / temp
	my_C_input[i] = math.sqrt(1.0 - (20.0* math.pi / INTERVAL.value)**2) * math.cos(20.0 * math.pi / INTERVAL.value * x) + (np.random.random_sample()/1 - 0.5)*2*rausch_amp
	my_C_input_cplx[i] = math.sqrt(1.0 - (20.0* math.pi / INTERVAL.value)**2) * math.sin(20.0 * math.pi / INTERVAL.value * x) + (np.random.random_sample()/1 - 0.5)*2*rausch_amp
	"""
	
	
	# Intermittenz
 	"""	"""
	x = float(i) / temp
	arg1 = ((x-100) + INTERVAL.value * 1.0 / 4)**2.0
	arg2 = ((x-100) - INTERVAL.value * 1.0 / 4)**2.0
	if(math.fabs(arg1) > 500 ):
		summ1 = 0
	else:
		summ1 = 1.0/math.cosh(arg1)
	if(math.fabs(arg2) > 500 ):
		summ2 = 0
	else:
		summ2 = 1.0/math.cosh(arg2)		
	my_C_input[i] = summ1 + 0.8 * summ2 + (np.random.random_sample()/1 - 0.5)*2*rausch_amp

#-----------------------------------

print "Calculating ... Start"
#Bibliothek aufrufen -> Dislin wird gestartet
ctypes.cdll.LoadLibrary(os.getcwd() + "/my_GLG_txtout.so").my_GLG_solver(my_C_input, my_C_input_cplx, ALPHA, BETA, TIMESTEPS, DeltaT, EACHNTH, INTERVAL, STEPS)
print "Calculating ... Finished"

SCALE = INTERVAL.value / float(STEPS)
IMAX = INTERVAL.value + I_OFFSET.value
MAXSTEPS = (TIMESTEPS.value / EACHNTH.value) - 1
TMAX = TIMESTEPS.value * DeltaT.value

print "Creating images ... Start"
subprocess.call(['gnuplot -e \"MAXSTEPS=\'' + str(MAXSTEPS) + '\'; SCALE=\'' + str(SCALE) + '\'; IMAX=\'' + str(IMAX) + '\'; I_OFFSET=\'' + str(I_OFFSET.value) + '\'\" ' + os.getcwd() + '/plotit_vid.gplt', '-1'], shell=True) # gnuplot Script ausfuehren

subprocess.call(['gnuplot -e \"SCALE=\'' + str(SCALE) + '\'; IMAX=\'' + str(IMAX) + '\'; I_OFFSET=\'' + str(I_OFFSET.value) + '\'; PICNAME=\'' + vid_name + '.png' + '\'; TMAX=\'' + str(TMAX) + '\'; TSCALE=\'' + str(DeltaT.value) + '\'\" ' + os.getcwd() + '/plotit_pic.gplt', '-1'], shell=True) # gnuplot Script ausfuehren

subprocess.call(['rm ' + os.getcwd() + '/GLG*', '-1'], shell=True) # Daten loeschen
print "Creating images ... Finished"
print "Creating video ... Start"
subprocess.call(['ffmpeg -r 25 -qscale 2 -i pic%04d.png ' + vid_name + '.avi', '-1'], shell=True) # Video erstellen
subprocess.call(['rm ' + os.getcwd() + '/pic*', '-1'], shell=True) # Bilder loeschen
print "Creating video ... Finished"

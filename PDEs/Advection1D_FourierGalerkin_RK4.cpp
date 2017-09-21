/**
by Frank Ehebrecht
im Mai 2014
cpplink -c advec -lfftw3
./advec
**/
#include <complex.h> //Dies unbedingt VOR fftw3.h einfuegen
#include <fftw3.h> //FFTW
#include <iostream>
#include <fstream> //Output nach txt
#include <iomanip> //Precision fuer txtOutput

#include <cmath>
#include "discpp.h" //Dislin
#include <unistd.h> //fuer sleep

using namespace std;

// PP-Parameter
#define TIMESTEPS 1024
#define SPACESTEPS 128
#define TIME 20.0
#define INTERVAL 2*M_PI
#define CC 0.3
#define ANI_SPEED 0.05

// Signal einlesen
void mFuncDef(fftw_complex* signal){
	double spacing = (INTERVAL) / SPACESTEPS;
	for(int i=0; i<SPACESTEPS; i++){
		signal[i] = exp(-2*M_PI*pow(  (spacing * i)  -M_PI/2,2));
	}
}

// Runge Kutta 4
void mRK4(fftw_complex signal[TIMESTEPS][SPACESTEPS]){
	double t_stepsize = TIME / TIMESTEPS;
	complex double k1, k2, k3, k4;
	double k = (2 * M_PI) / (INTERVAL);

	for(int t=0; t<TIMESTEPS-1; t++){
		for(int n=0; n<=SPACESTEPS/2; n++){
			k1 = -CC * I * n * k * signal[t][n];
			k2 = -CC * I * n * k * (signal[t][n] + t_stepsize/2 * k1);
			k3 = -CC * I * n * k * (signal[t][n] + t_stepsize/2 * k2);
			k4 = -CC * I * n * k * (signal[t][n] + t_stepsize * k3);
			signal[t+1][n] = signal[t][n] + t_stepsize/6 * (k1 + 2*k2 + 2*k3 + k4);	
		}
		for(int n=SPACESTEPS/2+1; n<SPACESTEPS; n++){
			k1 = -CC * I * (n-SPACESTEPS) * k * signal[t][n];
			k2 = -CC * I * (n-SPACESTEPS) * k * (signal[t][n] + t_stepsize/2 * k1);
			k3 = -CC * I * (n-SPACESTEPS) * k * (signal[t][n] + t_stepsize/2 * k2);
			k4 = -CC * I * (n-SPACESTEPS) * k * (signal[t][n] + t_stepsize * k3);
			signal[t+1][n] = signal[t][n] + t_stepsize/6 * (k1 + 2*k2 + 2*k3 + k4);
		}
	}
}

// DFT-Signal normalisieren
void mNorm(fftw_complex* signal){
	double dNorm = (1./SPACESTEPS);
	for(int i=0; i<SPACESTEPS; i++){
		signal[i] *= dNorm;
	}
}

// Mit Dislin Daten plotten
void mPlot(fftw_complex signal[TIMESTEPS][SPACESTEPS]){
	
	Dislin g;

	double aREAL[TIMESTEPS][SPACESTEPS];
	double aIMAG[TIMESTEPS][SPACESTEPS];
	double aSPACE[SPACESTEPS];
	double aEXACT[TIMESTEPS][SPACESTEPS];
	double spacing = 2 * M_PI / SPACESTEPS;
	double t_stepsize = TIME / TIMESTEPS;
	int ic;

	// sehr unschoen - weiss aber nicht, wie ich sonst die werte(typ complex) an dislin uebergeben soll :/
	for(int t=0; t<TIMESTEPS; t++){
		for(int x=0; x<SPACESTEPS; x++){
			aREAL[t][x] = creal(signal[t][x]); //Realteil
			aIMAG[t][x] = cimag(signal[t][x]); //Imaginaerteil
			aEXACT[t][x] = exp(-2*M_PI*pow(  (spacing * x - (CC*t_stepsize*t))  -M_PI/2,2)) + exp(-2*M_PI*pow(  (spacing * (x+SPACESTEPS) - (CC*t_stepsize*t))  -M_PI/2,2));
			// Exakte Loesung
		}
	}
	for(int x=0; x<SPACESTEPS; x++){
		aSPACE[x] = x * (2*M_PI)/SPACESTEPS;
	}


	// dislin parameter
	g.metafl ("cons");
	g.scrmod ("revers");
	g.disini ();
	g.pagera ();
	g.complx ();
	g.axspos (450, 1800);
	g.axslen (2200, 1200);

	g.name   ("x", "x");
	g.name   ("Amplitude", "y");

	g.labdig (-1, "x");
	g.ticks  (9, "x");
	g.ticks  (10, "y");

	g.titlin ("Advektions-Gleichung", 1);
	g.titlin ("Signal: exp(-2* pi * (x  - pi/2)^2)", 3);

	ic=g.intrgb (1.0,1.0,1.0);
	g.axsbgd (ic);

	g.graf   (0.0, 2*M_PI, 0.0, 1, 0.0, 1.0, 0, 0.5);
	g.setrgb (0.7, 0.7, 0.7);

	g.color  ("fore");
	g.height (50);
	g.title  ();
	

	g.incmrk (1);
	g.hsymbl (25);


	// Workaround fuer Animation: mit g.sendbf() Ausgabe an Bildschirm erzwingen,
	// der jeweils vorherige Schritt wird in Hintergrundsfarbe ueberschrieben
	// "black" = weiss; "white" = schwarz ....
	for(int i=0; i<TIMESTEPS-1; i++){
	
		g.incmrk (1); // Marker an
		g.polcrv("STEM"); // saeulen
		
		g.marker(16);
		g.color("black");
		g.curve(aSPACE, aREAL[i], SPACESTEPS);
		g.color("red");
		g.curve(aSPACE, aREAL[i+1], SPACESTEPS);
		
		
		g.incmrk (0); // Marker aus
		g.polcrv("LINEAR"); // linearer fit zwischen Punkten
		
		g.marker(16);
		g.color("black");
		g.curve(aSPACE, aEXACT[i], SPACESTEPS);
		g.color("white");
		g.curve(aSPACE, aEXACT[i+1], SPACESTEPS);
		
		
		g.sendbf();
		sleep(ANI_SPEED); //Pause - Ani.geschw. regulieren
	}
	g.disfin();
}


// +++ MAIN +++
int main(){


	// my_INPUT: Ortsraum; my_OUTPUT: Freqraum
	fftw_complex my_INPUT[TIMESTEPS][SPACESTEPS];
	fftw_complex my_OUTPUT[TIMESTEPS][SPACESTEPS];
	double exactsol[TIMESTEPS][SPACESTEPS];
	
	// Startsignal
	mFuncDef(my_INPUT[0]);

	// Vorwaerts-Trafo des Startsignals
	fftw_plan my_PLAN = fftw_plan_dft_1d(SPACESTEPS,		
						my_INPUT[0],	// InputDaten
						my_OUTPUT[0],	// OutputDaten
						FFTW_FORWARD,	// Vorwaerts-Trafo
						FFTW_ESTIMATE);	// in diesem Fall schnellere Variante
	fftw_execute(my_PLAN);	// oben geplantes ausfuehren
	fftw_destroy_plan(my_PLAN);

	// Normieren
	mNorm(my_OUTPUT[0]);
	
	// RK4 im Freqraum
	mRK4(my_OUTPUT);

	//Rueck-Trafo fuer alle Werte
	for(int k=1; k<TIMESTEPS; k++){
		fftw_plan my_PLAN = fftw_plan_dft_1d(SPACESTEPS,		
						my_OUTPUT[k],	// InputDaten
						my_INPUT[k],	// OutputDaten
						FFTW_BACKWARD,	// Rueck-Trafo
						FFTW_ESTIMATE);	// in diesem Fall schnellere Variante
		fftw_execute(my_PLAN);	// oben geplantes ausfuehren		
		fftw_destroy_plan(my_PLAN);
	}


	// Plot mit Dislin
	mPlot(my_INPUT);


/**+++++++++__TO TEXT FILE __++++++++++
	ofstream myfile;
	myfile.open ("out.txt");
	myfile.precision(6);
	myfile.setf( std::ios::fixed, std:: ios::floatfield );
	for(int j = 0; j<SPACESTEPS; j++){
		myfile << j << "\t" << creal(my_INPUT[200][j]) << "\t" << cimag(my_INPUT[200][j]) << "\n";
	}
	myfile.close();
_____________________________________**/

	return 0;
}

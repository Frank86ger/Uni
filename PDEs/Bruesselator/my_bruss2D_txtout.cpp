/**
by Frank Ehebrecht im Juli 2014

g++ -shared -c -fPIC my_bruss2D_txtout.cpp -o my_bruss2D_txtout.o -lfftw3
g++ -shared -Wl,-soname,my_bruss2D_txtout.so -o my_bruss2D_txtout.so my_bruss2D_txtout.o -lfftw3
**/

#include <complex.h> //Dies VOR fftw3.h einfuegen
#include <fftw3.h> //FFTW
#include <iostream>

#include <cmath> //Mathematische Funktionen
#include <unistd.h> //fuer sleep

#include <stdlib.h> //random

#include <fstream> //Output nach txt
#include <iomanip> //Precision fuer txtOutput
#include <sstream> //txt-Output richtig nummerieren koennen

# include <stdio.h>
# include <time.h>

using namespace std;

#define SSTEPS 256
#define HALFSTEPS 129
#define SSTEPS2 65536 // SSTEPS^2
#define ALLOC (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * SSTEPS * SSTEPS )
#define CPLX fftw_complex*


// Index in Reihen x Spalten Form
int mIndex(int iROW, int iCOL){
	return iROW * SSTEPS + iCOL;
}

//2D Signal Normalisieren
void mD2Norm(fftw_complex* INOUT){
	double dNorm = 1.0 / (SSTEPS2);
	for(int i=0; i<SSTEPS; i++){
		for(int j=0; j<SSTEPS; j++){
			INOUT[mIndex(i,j)] *= dNorm;
		}
	}
}

//x-Ableitung berechnen
void m2D_xDerive(fftw_complex* INOUT, double dInterval){
	double dK = 2*M_PI / dInterval;
	for(int yy=0; yy<SSTEPS; yy++){
		for(int xx=0; xx<HALFSTEPS; xx++){
			INOUT[mIndex(xx,yy)] *= I * dK * xx;
		}
		for(int xx=HALFSTEPS; xx<SSTEPS; xx++){
			INOUT[mIndex(xx,yy)] *= -1 * I * dK * (SSTEPS-xx);
		}
	}
}

//y-Ableitung berechnen
void m2D_yDerive(fftw_complex* INOUT, double dInterval){
	double dK = 2*M_PI / dInterval;
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<HALFSTEPS; yy++){
			INOUT[mIndex(xx,yy)] *= I * dK * yy;
		}
		for(int yy=HALFSTEPS; yy<SSTEPS; yy++){
			INOUT[mIndex(xx,yy)] *= -1 * I * dK * (SSTEPS-yy);
		}
	}
}

// Realteile aus Signal in Datei sFilename speichern mit Bildnummer
void m_R_data_to_file(fftw_complex* signal, string sFilename, int iBildnummer, int iShape[SSTEPS][SSTEPS]){
	
	ostringstream oss;
	oss << setfill('0') << setw(4) << iBildnummer;
	sFilename += oss.str();
	ofstream myfile;
	myfile.open (sFilename.c_str());
	myfile.precision(6);
	myfile.setf( std::ios::fixed, std:: ios::floatfield );
	for(int y=0; y<SSTEPS; y++){
		for(int x=0; x<SSTEPS; x++){
			if(iShape[x][y] == 0){
				
				myfile << creal(signal[mIndex(x,y)]) << "  ";				
				//myfile << -100 << "  ";
			}
			else{
				myfile << creal(signal[mIndex(x,y)]) << "  ";
			}
		}
		myfile << "\n";
	}
	myfile.close();
}

// OUT = IN1 + IN2
void mAdd(fftw_complex* IN1, fftw_complex* IN2, fftw_complex* OUT){
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			OUT[mIndex(xx,yy)] = IN1[mIndex(xx,yy)]+IN2[mIndex(xx,yy)];
		}
	}	
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// NichtLinearitaet
void mNonLin(fftw_plan PLAN_U_F2X, fftw_plan PLAN_V_F2X, fftw_plan PLAN_U_X2F_nl, fftw_plan PLAN_V_X2F_nl, 
	     CPLX X_U_nl, CPLX X_U, CPLX X_V, CPLX F_U_nl, int iShape[][SSTEPS], double dDamp,
	     CPLX F_U, CPLX F_V, CPLX X_V_nl, CPLX F_V_nl, double dA, double dB){
	     
	fftw_execute(PLAN_U_F2X); // FFT von U: F_U --> X_U
	fftw_execute(PLAN_V_F2X); // FFT von V: F_V --> X_V
	mD2Norm(X_U); // X_U normiert
	mD2Norm(X_V); // X_V normiert

	//Nicht-Linearitaet + Dampening
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			// NL
			X_U_nl[mIndex(xx,yy)] = dA + X_V[mIndex(xx,yy)] * X_U[mIndex(xx,yy)] * X_U[mIndex(xx,yy)];
			X_V_nl[mIndex(xx,yy)] = dB * X_U[mIndex(xx,yy)] - X_V[mIndex(xx,yy)] * X_U[mIndex(xx,yy)] * X_U[mIndex(xx,yy)];
			// Dampening fuer Shape
			if(iShape[xx][yy] == 0){
//				X_U_nl[mIndex(xx,yy)] += -1.0 * dDamp * X_U[mIndex(xx,yy)];
//				X_V_nl[mIndex(xx,yy)] += -1.0 * dDamp * X_V[mIndex(xx,yy)]; 
			}
		}
	}
	fftw_execute(PLAN_U_X2F_nl); // iFFT von U: X_U_nl --> F_U_nl (enthaelt nichtlinearitaet)
	fftw_execute(PLAN_V_X2F_nl); // iFFT von V: X_V_nl --> F_V_nl (enthaelt nichtlinearitaet)
}

// Semi-Implizites-Eulerverfahren
void mSIE(CPLX F_SIE_a, CPLX F_SIE_b, CPLX F_SIE_c, double dB, double dA, double dDU, double dDV, double deltaT, double dInterval){
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			F_SIE_a[mIndex(xx,yy)] = 1;
			F_SIE_b[mIndex(xx,yy)] = 1;
		}
	}
	m2D_xDerive(F_SIE_a, dInterval);
	m2D_xDerive(F_SIE_a, dInterval); //deldel x in a	
	m2D_yDerive(F_SIE_b, dInterval);
	m2D_yDerive(F_SIE_b, dInterval); //deldel y in b
	mAdd(F_SIE_a, F_SIE_b, F_SIE_b); //laplace in b
	
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			F_SIE_a[mIndex(xx,yy)] = dDU * F_SIE_b[mIndex(xx,yy)] - (dB + 1.0);
			F_SIE_b[mIndex(xx,yy)] = dDV * F_SIE_b[mIndex(xx,yy)];

			F_SIE_a[mIndex(xx,yy)] = 1.0 / (1.0/deltaT - F_SIE_a[mIndex(xx,yy)]);
			F_SIE_b[mIndex(xx,yy)] = 1.0 / (1.0/deltaT - F_SIE_b[mIndex(xx,yy)]);
		}
	}
}


//int main(){
extern "C" void my_Bruss2D_solver(double IN_U[SSTEPS][SSTEPS], double IN_V[SSTEPS][SSTEPS], int iShape[SSTEPS][SSTEPS], double dA, double dB, 
				double dDU, double dDV, double dInterval, double dDamp, double deltaT, int iEachNth, int iTimeSteps){
	
	// WERTE
//	double dA = 3;
//	double dB = 9;
//	double dDU = 5;
//	double dDV = 12;
	
//	double dInterval = 120;
//	double dDamp = 0.9;
//	double deltaT = 0.05;
	
//	int iEachNth = 1000;
	
	srand (time(NULL));
	
//	int iTimeSteps = 5000;
	
	// Shape
//	int iShape[SSTEPS][SSTEPS];

/**/	for(int yy=0; yy<SSTEPS; yy++){
		for(int xx=0; xx<SSTEPS; xx++){
			iShape[xx][yy] = 1;
			if (((xx-128)*(xx-128)) + ((yy-128)*(yy-128)) <= 8000.0){
				iShape[xx][yy] = 1;
			}
		}
	}


	// Speicher allocieren
	fftw_complex *X_U, *F_U, *X_U_nl, *F_U_nl, *X_V, *F_V, *X_V_nl, *F_V_nl;
	fftw_complex *F_SIE_a, *F_SIE_b, *F_SIE_c;

	X_U = ALLOC; F_U = ALLOC; X_V = ALLOC; F_V = ALLOC; X_U_nl = ALLOC; F_U_nl = ALLOC;
	X_V_nl = ALLOC; F_V_nl = ALLOC; F_SIE_a = ALLOC; F_SIE_b = ALLOC; F_SIE_c = ALLOC;
	
	
	// FFTW3 - Plaene erstellen
	fftw_plan PLAN_U_X2F = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_U, F_U, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_U_F2X = fftw_plan_dft_2d(SSTEPS, SSTEPS, F_U, X_U, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_plan PLAN_V_X2F = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_V, F_V, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_V_F2X = fftw_plan_dft_2d(SSTEPS, SSTEPS, F_V, X_V, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_plan PLAN_U_X2F_nl = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_U_nl, F_U_nl, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_V_X2F_nl = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_V_nl, F_V_nl, FFTW_FORWARD, FFTW_ESTIMATE);

////////////////////////////////////////////////////////////////////////////

	//////////// INPUTSIGNAL //////////////////
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			X_U[mIndex(xx, yy)] = IN_U[xx][yy];
			X_V[mIndex(xx, yy)] = IN_V[xx][yy];
//			X_U[mIndex(xx, yy)] = dA + (((double)(rand()%1000)/10000)-0.05) * 2.;
//			X_V[mIndex(xx, yy)] = dB/dA + (((double)(rand()%1000)/10000)-0.05) * 2.;			
	}}/////////////////////////////////////////

	// Koeffizienten fuer Semi Implizites Eulerverfahren bestimmen (a,b,c)
	mSIE(F_SIE_a, F_SIE_b, F_SIE_c, dB, dA, dDU, dDV, deltaT, dInterval);

	// Vorwaertstrafo von U und V
	fftw_execute(PLAN_U_X2F); // FFT von U: X_U --> F_U
	fftw_execute(PLAN_V_X2F); // FFT von V: X_V --> F_V

	
	// Fuer Zeitschritte tt
	for(int tt=0; tt<iTimeSteps; tt++){

		// Nichtlinearitaet (hier werden auch die Signale ruecktransformiert)
		mNonLin(PLAN_U_F2X, PLAN_V_F2X, PLAN_U_X2F_nl, PLAN_V_X2F_nl, X_U_nl, X_U, X_V, F_U_nl, iShape, dDamp, F_U, F_V, X_V_nl, F_V_nl, dA, dB);
		
		// Ausgabe alle EachNth Schritte
		if(tt%iEachNth == 0){ 
			m_R_data_to_file(X_U, "Bruss_u", tt/iEachNth, iShape);
			m_R_data_to_file(X_V, "Bruss_v", tt/iEachNth, iShape);
		}
		
		// Schritt
		for(int xx=0; xx<SSTEPS; xx++){
			for(int yy=0; yy<SSTEPS; yy++){
				F_V[mIndex(xx, yy)] = F_SIE_b[mIndex(xx, yy)] * (F_V[mIndex(xx, yy)]/deltaT + F_V_nl[mIndex(xx, yy)]);
				F_U[mIndex(xx, yy)] = F_SIE_a[mIndex(xx, yy)] * (F_U[mIndex(xx, yy)]/deltaT + F_U_nl[mIndex(xx, yy)]);
			}
		}
	}

	m_R_data_to_file(X_U, "BrussENDU", 0, iShape);
	m_R_data_to_file(X_V, "BrussENDV", 0, iShape);



////////////////////////////////////////////////////////////////////////////

	// Speicher freigeben
	fftw_free(X_U); fftw_free(F_U); fftw_free(X_U_nl); fftw_free(F_U_nl);
	fftw_free(X_V); fftw_free(F_V); fftw_free(X_V_nl); fftw_free(F_V_nl);
	fftw_free(F_SIE_a); fftw_free(F_SIE_b); fftw_free(F_SIE_c);

	fftw_destroy_plan(PLAN_U_X2F); 	fftw_destroy_plan(PLAN_U_F2X);
	fftw_destroy_plan(PLAN_V_X2F); fftw_destroy_plan(PLAN_V_F2X);
	fftw_destroy_plan(PLAN_U_X2F_nl);
	fftw_destroy_plan(PLAN_V_X2F_nl);
}

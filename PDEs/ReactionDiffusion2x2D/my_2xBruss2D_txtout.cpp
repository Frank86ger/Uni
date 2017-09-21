/**
by Frank Ehebrecht im Juli 2014

g++ -shared -c -fPIC my_2xBruss2D_txtout.cpp -o my_2xBruss2D_txtout.o -lfftw3
g++ -shared -Wl,-soname,my_2xBruss2D_txtout.so -o my_2xBruss2D_txtout.so my_2xBruss2D_txtout.o -lfftw3
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

// NL
void mNonLin(fftw_plan PLAN_U1_F2X, fftw_plan PLAN_U2_F2X, fftw_plan PLAN_V1_F2X, fftw_plan PLAN_V2_F2X,
	     CPLX X_U1, CPLX X_U2, CPLX X_V1, CPLX X_V2, CPLX X_U1_nl, CPLX X_U2_nl, CPLX X_V1_nl, CPLX X_V2_nl,
	     fftw_plan PLAN_U1_X2F_nl, fftw_plan PLAN_U2_X2F_nl, fftw_plan PLAN_V1_X2F_nl, fftw_plan PLAN_V2_X2F_nl,
	     double dAlpha, double dA, double dB, double dDamp, int iShape[SSTEPS][SSTEPS]){


	fftw_execute(PLAN_U1_F2X); // FFT von U1: F_U1 --> X_U1
	fftw_execute(PLAN_U2_F2X); // FFT von U2: F_U2 --> X_U2
	fftw_execute(PLAN_V1_F2X); // FFT von V1: F_V1 --> X_V1
	fftw_execute(PLAN_V2_F2X); // FFT von V2: F_V2 --> X_V2
	mD2Norm(X_U1); // X_U1 normiert
	mD2Norm(X_U2); // X_U2 normiert
	mD2Norm(X_V1); // X_V1 normiert
	mD2Norm(X_V2); // X_V2 normiert

	//Nicht-Linearitaet + Dampening
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			X_U1_nl[mIndex(xx,yy)] = dAlpha * X_U2[mIndex(xx,yy)] + dA + X_U1[mIndex(xx,yy)]*X_U1[mIndex(xx,yy)]*X_V1[mIndex(xx,yy)];
			X_U2_nl[mIndex(xx,yy)] = dAlpha * X_U1[mIndex(xx,yy)] + dA + X_U2[mIndex(xx,yy)]*X_U2[mIndex(xx,yy)]*X_V2[mIndex(xx,yy)];
			X_V1_nl[mIndex(xx,yy)] = dAlpha * X_V2[mIndex(xx,yy)] + dB * X_U1[mIndex(xx,yy)] - X_U1[mIndex(xx,yy)]*X_U1[mIndex(xx,yy)]*X_V1[mIndex(xx,yy)];
			X_V2_nl[mIndex(xx,yy)] = dAlpha * X_V1[mIndex(xx,yy)] + dB * X_U2[mIndex(xx,yy)] - X_U2[mIndex(xx,yy)]*X_U2[mIndex(xx,yy)]*X_V2[mIndex(xx,yy)];
			// Dampening fuer Shape
			if(iShape[xx][yy] == 0){
				X_U1_nl[mIndex(xx,yy)] += -1.0 * dDamp * X_U1[mIndex(xx,yy)];
				X_U2_nl[mIndex(xx,yy)] += -1.0 * dDamp * X_U2[mIndex(xx,yy)];
				X_V1_nl[mIndex(xx,yy)] += -1.0 * dDamp * X_V1[mIndex(xx,yy)];
				X_V2_nl[mIndex(xx,yy)] += -1.0 * dDamp * X_V2[mIndex(xx,yy)];
			}
		}
	}
	fftw_execute(PLAN_U1_X2F_nl); // iFFT von U1: X_U1_nl --> F_U1_nl
	fftw_execute(PLAN_U2_X2F_nl); // iFFT von U2: X_U2_nl --> F_U2_nl
	fftw_execute(PLAN_V1_X2F_nl); // iFFT von V1: X_V1_nl --> F_V1_nl
	fftw_execute(PLAN_V2_X2F_nl); // iFFT von V2: X_V2_nl --> F_V2_nl
}



// Semi-Implizites-Eulerverfahren (SIE)
void mSIE(CPLX F_SIE_a, CPLX F_SIE_b, CPLX F_SIE_c, CPLX F_SIE_d, double dB, double dAlpha, double dDU1, double dDU2, double dDV1, double dDV2, double deltaT, double dInterval){
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
	mAdd(F_SIE_a, F_SIE_b, F_SIE_d); //laplace in d
	
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			F_SIE_a[mIndex(xx,yy)] = dDU1 * F_SIE_d[mIndex(xx,yy)] - dAlpha - (dB + 1.0);
			F_SIE_b[mIndex(xx,yy)] = dDU2 * F_SIE_d[mIndex(xx,yy)] - dAlpha - (dB + 1.0);
			F_SIE_c[mIndex(xx,yy)] = dDV1 * F_SIE_d[mIndex(xx,yy)] - dAlpha;
			F_SIE_d[mIndex(xx,yy)] = dDV2 * F_SIE_d[mIndex(xx,yy)] - dAlpha;
			
			F_SIE_a[mIndex(xx,yy)] = 1.0 / (1.0/deltaT - F_SIE_a[mIndex(xx,yy)]);
			F_SIE_b[mIndex(xx,yy)] = 1.0 / (1.0/deltaT - F_SIE_b[mIndex(xx,yy)]);
			F_SIE_c[mIndex(xx,yy)] = 1.0 / (1.0/deltaT - F_SIE_c[mIndex(xx,yy)]);
			F_SIE_d[mIndex(xx,yy)] = 1.0 / (1.0/deltaT - F_SIE_d[mIndex(xx,yy)]);
		}
	}
}



//int main(){
extern "C" void my_2xBruss2D_solver(double IN_U1[SSTEPS][SSTEPS], double IN_U2[SSTEPS][SSTEPS], double IN_V1[SSTEPS][SSTEPS], 
				double IN_V2[SSTEPS][SSTEPS], int iShape[SSTEPS][SSTEPS], double dA, double dB,
				double dAlpha, double dDU1, double dDV1, double dDU2, double dDV2,
				double dInterval, double dDamp, double deltaT, int iEachNth, int iTimeSteps){
	
	// WERTE
//	double dA = 3;
//	double dB = 9;
//	double dAlpha = 0.1;
//	double dDU1 = 12.6;
//	double dDV1 = 27.5;
//	double dDU2 = 47.5;
//	double dDV2 = 141.5;
	
//	double dVH1 = 3;
//	double dVH2 = 3;
	
//	double dInterval = 200;
//	double dDamp = 0.9;
//	double deltaT = 0.05;
	double dInvDT = 1./deltaT;
	
//	int iEachNth = 100;
	
	srand (time(NULL));
	
//	int iTimeSteps = 1000;
	
	// Shape
//	int iShape[SSTEPS][SSTEPS];
//	for(int yy=0; yy<SSTEPS; yy++){
//		for(int xx=0; xx<SSTEPS; xx++){
//			iShape[xx][yy] = 1;
//			if (((xx-128)*(xx-128)) + ((yy-128)*(yy-128)) <= 100.0){
//				iShape[xx][yy] = 1;
//			}
//		}
//	}


	// Speicher allocieren
	fftw_complex *X_U1, *F_U1, *X_U1_nl, *F_U1_nl;
	fftw_complex *X_U2, *F_U2, *X_U2_nl, *F_U2_nl;
	fftw_complex *X_V1, *F_V1, *X_V1_nl, *F_V1_nl;
	fftw_complex *X_V2, *F_V2, *X_V2_nl, *F_V2_nl;
	fftw_complex *F_SIE_a, *F_SIE_b, *F_SIE_c, *F_SIE_d;

	X_U1 = ALLOC; F_U1 = ALLOC; X_U1_nl = ALLOC; F_U1_nl = ALLOC;
	X_U2 = ALLOC; F_U2 = ALLOC; X_U2_nl = ALLOC; F_U2_nl = ALLOC;
	X_V1 = ALLOC; F_V1 = ALLOC; X_V1_nl = ALLOC; F_V1_nl = ALLOC;
	X_V2 = ALLOC; F_V2 = ALLOC; X_V2_nl = ALLOC; F_V2_nl = ALLOC;
	
	F_SIE_a = ALLOC; F_SIE_b = ALLOC; F_SIE_c = ALLOC; F_SIE_d = ALLOC;
	
	
	// FFTW3 - Plaene erstellen
	fftw_plan PLAN_U1_X2F = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_U1, F_U1, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_U1_F2X = fftw_plan_dft_2d(SSTEPS, SSTEPS, F_U1, X_U1, FFTW_BACKWARD, FFTW_ESTIMATE);

	fftw_plan PLAN_U2_X2F = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_U2, F_U2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_U2_F2X = fftw_plan_dft_2d(SSTEPS, SSTEPS, F_U2, X_U2, FFTW_BACKWARD, FFTW_ESTIMATE);
		
	fftw_plan PLAN_V1_X2F = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_V1, F_V1, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_V1_F2X = fftw_plan_dft_2d(SSTEPS, SSTEPS, F_V1, X_V1, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_plan PLAN_V2_X2F = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_V2, F_V2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_V2_F2X = fftw_plan_dft_2d(SSTEPS, SSTEPS, F_V2, X_V2, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_plan PLAN_U1_X2F_nl = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_U1_nl, F_U1_nl, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_U2_X2F_nl = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_U2_nl, F_U2_nl, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_V1_X2F_nl = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_V1_nl, F_V1_nl, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_V2_X2F_nl = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_V2_nl, F_V2_nl, FFTW_FORWARD, FFTW_ESTIMATE);

////////////////////////////////////////////////////////////////////////////

	//////////// INPUTSIGNAL //////////////////
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			X_U1[mIndex(xx, yy)] = IN_U1[xx][yy];//dVH1 + (((double)(rand()%1000)/10000)-0.05) * 2.;
			X_U2[mIndex(xx, yy)] = IN_U2[xx][yy];//dVH1 + (((double)(rand()%1000)/10000)-0.05) * 2.;
			X_V1[mIndex(xx, yy)] = IN_V1[xx][yy];//dVH2 + (((double)(rand()%1000)/10000)-0.05) * 2.;
			X_V2[mIndex(xx, yy)] = IN_V2[xx][yy];//dVH2 + (((double)(rand()%1000)/10000)-0.05) * 2.;
	}}/////////////////////////////////////////

	// Koeffizienten fuer Semi Implizites Eulerverfahren bestimmen (a,b,c, d)
	mSIE(F_SIE_a, F_SIE_b, F_SIE_c, F_SIE_d, dB, dAlpha, dDU1, dDU2, dDV1, dDV2, deltaT, dInterval);

	// Vorwaertstrafo von U1,U2,V1,V2
	fftw_execute(PLAN_U1_X2F); // FFT von U1: X_U1 --> F_U1
	fftw_execute(PLAN_U2_X2F); // FFT von U2: X_U2 --> F_U2
	fftw_execute(PLAN_V1_X2F); // FFT von V1: X_V1 --> F_V1
	fftw_execute(PLAN_V2_X2F); // FFT von V2: X_V2 --> F_V2


	// Fuer Zeitschritte tt
	for(int tt=0; tt<iTimeSteps; tt++){

		// Nicht-Linearitaet berechnen
		mNonLin(PLAN_U1_F2X, PLAN_U2_F2X, PLAN_V1_F2X, PLAN_V2_F2X, X_U1, X_U2, X_V1, X_V2, X_U1_nl, X_U2_nl, X_V1_nl, X_V2_nl, PLAN_U1_X2F_nl, PLAN_U2_X2F_nl, PLAN_V1_X2F_nl, PLAN_V2_X2F_nl, dAlpha, dA, dB, dDamp, iShape);


		if(tt%iEachNth == 0){ 
			m_R_data_to_file(X_U1, "Bruss_u1", tt/iEachNth, iShape);
			m_R_data_to_file(X_U2, "Bruss_u2", tt/iEachNth, iShape);
			m_R_data_to_file(X_V1, "Bruss_v1", tt/iEachNth, iShape);
			m_R_data_to_file(X_V2, "Bruss_v2", tt/iEachNth, iShape);
		}


		// Schritt
		for(int xx=0; xx<SSTEPS; xx++){
			for(int yy=0; yy<SSTEPS; yy++){
				F_U1[mIndex(xx, yy)] = F_SIE_a[mIndex(xx, yy)] * (dInvDT*F_U1[mIndex(xx, yy)] + F_U1_nl[mIndex(xx, yy)]);
				F_U2[mIndex(xx, yy)] = F_SIE_b[mIndex(xx, yy)] * (dInvDT*F_U2[mIndex(xx, yy)] + F_U2_nl[mIndex(xx, yy)]);
				F_V1[mIndex(xx, yy)] = F_SIE_c[mIndex(xx, yy)] * (dInvDT*F_V1[mIndex(xx, yy)] + F_V1_nl[mIndex(xx, yy)]);
				F_V2[mIndex(xx, yy)] = F_SIE_d[mIndex(xx, yy)] * (dInvDT*F_V2[mIndex(xx, yy)] + F_V2_nl[mIndex(xx, yy)]);
			}
		}
	}


	// Rueckwaertstrafo von U und V
//	fftw_execute(PLAN_U1_F2X); // iFFT von U: F_U --> X_U
//	fftw_execute(PLAN_U2_F2X); // iFFT von U: F_U --> X_U
//	fftw_execute(PLAN_V1_F2X); // iFFT von V: F_V --> X_V
//	fftw_execute(PLAN_V2_F2X); // iFFT von V: F_V --> X_V
	
//	mD2Norm(X_U1);
//	mD2Norm(X_U2);
//	mD2Norm(X_V1);
//	mD2Norm(X_V2);

	// hier endwerte noch raus
	m_R_data_to_file(X_U1, "BrussENDU1", 0, iShape);
	m_R_data_to_file(X_U2, "BrussENDU2", 0, iShape);
	m_R_data_to_file(X_V1, "BrussENDV1", 0, iShape);
	m_R_data_to_file(X_V2, "BrussENDV2", 0, iShape);


////////////////////////////////////////////////////////////////////////////

	// Speicher freigeben
	fftw_free(X_U1); fftw_free(F_U1); fftw_free(X_U1_nl); fftw_free(F_U1_nl);
	fftw_free(X_U2); fftw_free(F_U2); fftw_free(X_U2_nl); fftw_free(F_U2_nl);
	fftw_free(X_V1); fftw_free(F_V1); fftw_free(X_V1_nl); fftw_free(F_V1_nl);
	fftw_free(X_V2); fftw_free(F_V2); fftw_free(X_V2_nl); fftw_free(F_V2_nl);
	fftw_free(F_SIE_a); fftw_free(F_SIE_b); fftw_free(F_SIE_c); fftw_free(F_SIE_d);

	fftw_destroy_plan(PLAN_U1_X2F); fftw_destroy_plan(PLAN_U1_F2X);
	fftw_destroy_plan(PLAN_U2_X2F); fftw_destroy_plan(PLAN_U2_F2X);
	fftw_destroy_plan(PLAN_V1_X2F); fftw_destroy_plan(PLAN_V1_F2X);
	fftw_destroy_plan(PLAN_V2_X2F); fftw_destroy_plan(PLAN_V2_F2X);
	fftw_destroy_plan(PLAN_U1_X2F_nl); fftw_destroy_plan(PLAN_U2_X2F_nl);
	fftw_destroy_plan(PLAN_V1_X2F_nl); fftw_destroy_plan(PLAN_V2_X2F_nl);
}


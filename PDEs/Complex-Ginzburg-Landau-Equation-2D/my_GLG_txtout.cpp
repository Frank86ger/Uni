/**
by Frank Ehebrecht im Juli 2014

g++ -shared -c -fPIC my_GLG_txtout.cpp -o my_GLG_txtout.o -lfftw3
g++ -shared -Wl,-soname,my_GLG_txtout.so -o my_GLG_txtout.so my_GLG_txtout.o -lfftw3
**/

#include <complex.h> //Dies VOR fftw3.h einfuegen
#include <fftw3.h> //FFTW
#include <iostream>

#include <cmath> //Mathematische Funktionen
//#include "discpp.h" //Dislin
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

// OUT = IN1 + IN2
void mAdd(fftw_complex* IN1, fftw_complex* IN2, fftw_complex* OUT){
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			OUT[mIndex(xx,yy)] = IN1[mIndex(xx,yy)]+IN2[mIndex(xx,yy)];
		}
	}	
}

// INPUT nach OUTPUT kopieren
void mCopy(fftw_complex* INPUT, fftw_complex* OUTPUT){
	for(int ii=0; ii<SSTEPS; ii++){
		for(int jj=0; jj<SSTEPS; jj++){
			OUTPUT[mIndex(ii,jj)] = INPUT[mIndex(ii,jj)];
		}
	}
}


// Nichtlinearitaet berechnen / F_Signal_nl wird ist Ausgabegroesse (Nichtlinearitaet in F)
void mNonLin(fftw_plan PLAN_F2X, fftw_plan PLAN_X2F_nl, fftw_complex* X_Signal_nl, fftw_complex* X_Signal, fftw_complex* F_Signal_nl, int iShape[][SSTEPS], double dDamp, double dBeta){ 
//PLAN_F2X - PLAN_X2F_nonlin - X_Signal_nonlin X_Signal dBeta, dDamp F_Signal_nonlin
	fftw_execute(PLAN_F2X); //X_Signal liegt vor
	mD2Norm(X_Signal); //X_Signal normiert

	//Nicht-Linearitaet + Dampening
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			X_Signal_nl[mIndex(xx,yy)] = -1.0 * X_Signal[mIndex(xx,yy)] * conj(X_Signal[mIndex(xx,yy)]) * X_Signal[mIndex(xx,yy)] * (1 + I*dBeta);
			if(iShape[xx][yy] == 0){ //Dampening
				X_Signal_nl[mIndex(xx,yy)] += -1.0 * dDamp * X_Signal[mIndex(xx,yy)]; 
			}
		}
	}
	fftw_execute(PLAN_X2F_nl);
	// + Orszag auf non F?
}

// ETD2-Koeffizienten. a: exp(qh) / b fuer N_n / c fuer N_(n-1)
void mETD2(fftw_complex* F_ETD2_a, fftw_complex* F_ETD2_b, fftw_complex* F_ETD2_c, double dAlpha, double deltaT, double dInterval){
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			F_ETD2_a[mIndex(xx,yy)] = 1;
			F_ETD2_b[mIndex(xx,yy)] = 1;
		}
	}
	m2D_xDerive(F_ETD2_a, dInterval);
	m2D_xDerive(F_ETD2_a, dInterval); //del x ^2 in a
	
	m2D_yDerive(F_ETD2_b, dInterval);
	m2D_yDerive(F_ETD2_b, dInterval); //del y ^2 in a
	
	mAdd(F_ETD2_a, F_ETD2_b, F_ETD2_b); //laplace in b
	
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			F_ETD2_b[mIndex(xx,yy)] = F_ETD2_b[mIndex(xx,yy)] * (1.0 + I * dAlpha) + 1.0; // ((1 + ia)lap) + 1 = q   in b
			F_ETD2_a[mIndex(xx,yy)] = cexp(F_ETD2_b[mIndex(xx,yy)]*deltaT); //exp(qh)   in a
			
			if( creal(F_ETD2_b[mIndex(xx,yy)])==0 & cimag(F_ETD2_b[mIndex(xx,yy)])==0 ){
				F_ETD2_b[mIndex(xx,yy)] = 3.0 * deltaT / 2.0;
				F_ETD2_c[mIndex(xx,yy)] = -1.0 * deltaT / 2.0;
				
			}
			else{
				F_ETD2_c[mIndex(xx,yy)] = (1.0 + deltaT * F_ETD2_b[mIndex(xx,yy)] - F_ETD2_a[mIndex(xx,yy)]) / (deltaT * F_ETD2_b[mIndex(xx,yy)] * F_ETD2_b[mIndex(xx,yy)]);
				F_ETD2_b[mIndex(xx,yy)] = ((1.0+deltaT*F_ETD2_b[mIndex(xx,yy)]) * F_ETD2_a[mIndex(xx,yy)] - 1.0 - 2.0 * deltaT * F_ETD2_b[mIndex(xx,yy)]) / (deltaT * F_ETD2_b[mIndex(xx,yy)] * F_ETD2_b[mIndex(xx,yy)]);
			}
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
				myfile << -100 << "  ";
			}
			else{
				myfile << creal(signal[mIndex(x,y)]) << "  ";
			}
		}
		myfile << "\n";
	}
	myfile.close();
}

// Imagteile aus Signal in Datei sFilename speichern mit Bildnummer
void m_C_data_to_file(fftw_complex* signal, string sFilename, int iBildnummer, int iShape[SSTEPS][SSTEPS]){
	
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
				myfile << -100 << "  ";
			}
			else{
				myfile << cimag(signal[mIndex(x,y)]) << "  ";
			}
		}
		myfile << "\n";
	}
	myfile.close();
}

// Abs von Signal in Datei sFilename speichern mit Bildnummer
void m_A_data_to_file(fftw_complex* signal, string sFilename, int iBildnummer, int iShape[SSTEPS][SSTEPS]){
	
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
				myfile << -100 << "  ";
			}
			else{
				myfile << cabs(signal[mIndex(x,y)]) << "  ";
			}
		}
		myfile << "\n";
	}
	myfile.close();
}

//int main(){
extern "C" void my_GLG_solver(double IN_real[SSTEPS][SSTEPS], double IN_imag[SSTEPS][SSTEPS], int iShape[SSTEPS][SSTEPS], double dAlpha, double dBeta, int iTimeSteps, double deltaT, int iEachNth, double dInterval, double dDamp){



	
	//double dInterval = 200;
	//double dDamp = 1;
	//double dBeta = 1.5;
	//double dAlpha = 0;
	//double deltaT = 0.05;
	
	//int iEachNth = 20;
	
	srand (time(NULL));
	
	//int iTimeSteps = 500;
	
	//int iShape[SSTEPS][SSTEPS];
	//for(int y=0; y<SSTEPS; y++){
	//	for(int x=0; x<SSTEPS; x++){
	//		iShape[x][y] = 1;
	//	}
	//}
	
	
		  
	fftw_complex *X_Signal, *F_Signal, *X_Signal_nl, *F_Signal_nl, *F_ETD2_a, *F_ETD2_b, *F_ETD2_c, *F_Signal_nl_temp;
		  
	X_Signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * SSTEPS * SSTEPS );
	F_Signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * SSTEPS * SSTEPS );
	
	X_Signal_nl = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * SSTEPS * SSTEPS );
	F_Signal_nl = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * SSTEPS * SSTEPS );
	
	F_Signal_nl_temp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * SSTEPS * SSTEPS );
	
	F_ETD2_a = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * SSTEPS * SSTEPS );
	F_ETD2_b = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * SSTEPS * SSTEPS );
	F_ETD2_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * SSTEPS * SSTEPS );
				
	
	fftw_plan PLAN_X2F = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_Signal, F_Signal, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_F2X = fftw_plan_dft_2d(SSTEPS, SSTEPS, F_Signal, X_Signal, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_plan PLAN_X2F_nl = fftw_plan_dft_2d(SSTEPS, SSTEPS, X_Signal_nl, F_Signal_nl, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_F2X_nl = fftw_plan_dft_2d(SSTEPS, SSTEPS, F_Signal_nl, X_Signal_nl, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	// INPUTSIGNAL //////////////////
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			X_Signal[mIndex(xx, yy)] = IN_real[xx][yy] + I * IN_imag[xx][yy];
			//X_Signal[mIndex(xx, yy)] = 0 + ((double)(rand()%1000)/10000)-0.05;
		}
	}
	/////////////////////////////////
	
	
	fftw_execute(PLAN_X2F); // FFT
	
	mETD2(F_ETD2_a, F_ETD2_b, F_ETD2_c, dAlpha, deltaT, dInterval); //ETD2 - Koeffizienten berechnen
	mNonLin(PLAN_F2X, PLAN_X2F_nl, X_Signal_nl, X_Signal, F_Signal_nl, iShape, dDamp, dBeta); // Nicht linearitaet berechnen
	
	// Ersten Zeitschritt berechnen
	for(int xx=0; xx<SSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			F_Signal[mIndex(xx,yy)] = F_ETD2_a[mIndex(xx,yy)] * F_Signal[mIndex(xx,yy)] + F_ETD2_b[mIndex(xx,yy)] * F_Signal_nl[mIndex(xx,yy)] + F_ETD2_c[mIndex(xx,yy)] * F_Signal_nl[mIndex(xx,yy)];
		}
	}
	mCopy(F_Signal_nl, F_Signal_nl_temp);
	m_R_data_to_file(X_Signal, "GLGreal", 0, iShape);
//	m_C_data_to_file(X_Signal, "GLGimag", 0, iShape);
	m_A_data_to_file(X_Signal, "GLGabs", 0, iShape);
	
	// Alle Zeitschritte berechnen
	for(int tt=1; tt<iTimeSteps; tt++){
	
		if(tt%iEachNth == 0){ // Ausgabe
			//system("clear");
			//cout << "Daten erzeugen: " << setprecision(3) << ((float)(tt)) / iTimeSteps * 100.0 << "%\n";
			
			m_R_data_to_file(X_Signal, "GLGreal", tt/iEachNth, iShape);
//			m_C_data_to_file(X_Signal, "GLGimag", tt/iEachNth, iShape);
			m_A_data_to_file(X_Signal, "GLGabs", tt/iEachNth, iShape);
		}

		mNonLin(PLAN_F2X, PLAN_X2F_nl, X_Signal_nl, X_Signal, F_Signal_nl, iShape, dDamp, dBeta);
		
		//Zeitschritt
		for(int xx=0; xx<SSTEPS; xx++){
			for(int yy=0; yy<SSTEPS; yy++){
				F_Signal[mIndex(xx,yy)] = F_ETD2_a[mIndex(xx,yy)] * F_Signal[mIndex(xx,yy)] + F_ETD2_b[mIndex(xx,yy)] * F_Signal_nl[mIndex(xx,yy)] + F_ETD2_c[mIndex(xx,yy)] * F_Signal_nl_temp[mIndex(xx,yy)];
			}
		}
		mCopy(F_Signal_nl, F_Signal_nl_temp);
	}

	//system("clear");
	//cout << "Daten erzeugen: 100%\n";
	
//	fftw_execute(PLAN_F2X);	
//	mD2Norm(X_Signal);
	
	// Speicher freigeben
	fftw_free(X_Signal);
	fftw_free(F_Signal);

	fftw_free(X_Signal_nl);
	fftw_free(F_Signal_nl);
	
	fftw_free(F_ETD2_a);
	fftw_free(F_ETD2_b);
	fftw_free(F_ETD2_c);

	fftw_free(F_Signal_nl_temp);
	
}

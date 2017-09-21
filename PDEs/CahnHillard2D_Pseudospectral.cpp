/**
by Frank Ehebrecht
im Mai 2014

mit main:
cpplink -c my_burgers_lib -lfftw3
./my_burgers_lib
**/
#include <complex.h> //Dies unbedingt VOR fftw3.h einfuegen
#include <fftw3.h> //FFTW
#include <iostream>

#include <cmath>
//#include "discpp.h" //Dislin
#include <unistd.h> //fuer sleep

#include <stdlib.h> // random

#include <fstream> //Output nach txt
#include <iomanip> //Precision fuer txtOutput
#include <sstream> // txt-Output richtig nummerieren koennen

using namespace std;

#define SSTEPS 128
#define CPLXSTEPS 65 // SSTEPS/2 + 1
#define SSTEPS2 16384 // SSTEPS^2


//2D Signal Normalisieren
void D2Norm(double INOUT[][SSTEPS]){
	double dNorm = 1.0 / (SSTEPS2);
	for(int i=0; i<SSTEPS; i++){
		for(int j=0; j<SSTEPS; j++){
			INOUT[i][j] *= dNorm;
		}
	}
}

// y-Ableitung
void m2D_yDerive(fftw_complex INOUT[][CPLXSTEPS], double dInterval){
	double dK = 2*M_PI / dInterval;

	for(int y=0; y<CPLXSTEPS; y++){
		for(int x=0; x<CPLXSTEPS; x++){
			INOUT[y][x] *= I * dK * x;
		}
	}
	for(int y=CPLXSTEPS; y<SSTEPS; y++){
		for(int x=0; x<CPLXSTEPS; x++){
			INOUT[y][x] *= I*dK*x;
		}
	}
}

// x-Ableitung
void m2D_xDerive(fftw_complex INOUT[][CPLXSTEPS], double dInterval){
	double dK = 2*M_PI / dInterval;
	for(int x=0; x<CPLXSTEPS; x++){
		for(int y=0; y<CPLXSTEPS; y++){
			INOUT[y][x] *= I * dK * y;
		}
		for(int y=CPLXSTEPS; y<SSTEPS; y++){
			INOUT[y][x] *= -1 * I * dK * (SSTEPS-y);
		}
	}
}

// INPUT-Array nach OUTPUT kopieren (F-Space)
double mCplxCopy(fftw_complex INPUT[][CPLXSTEPS], fftw_complex OUTPUT[][CPLXSTEPS]){
	for(int i=0; i<SSTEPS; i++){
		for(int j=0; j<CPLXSTEPS; j++){
			OUTPUT[i][j] = INPUT[i][j];
		}
	}
}

// IN1 + IN2 = OUT (X-Space)
void mRealAdd(double IN1[][SSTEPS], double IN2[][SSTEPS], double OUT[][SSTEPS]){
	for(int i=0; i<SSTEPS; i++){
		for(int j=0; j<SSTEPS; j++){
			OUT[i][j] = IN1[i][j]+IN2[i][j];
		}
	}
}

// IN1 + IN2 = OUT (F-Space)
void mCplxAdd(fftw_complex IN1[][CPLXSTEPS], fftw_complex IN2[][CPLXSTEPS], fftw_complex OUT[][CPLXSTEPS]){
	for(int i=0; i<SSTEPS; i++){
		for(int j=0; j<CPLXSTEPS; j++){
			OUT[i][j] = IN1[i][j]+IN2[i][j];
		}
	}
}

// INOUT <-> Orszag Filter
void mOrszag(fftw_complex INOUT[][CPLXSTEPS]){
	//double dFilter= 2.0*(double)CPLXSTEPS/3;
	int dFilter = trunc(2.0*(double)CPLXSTEPS/3);

	for(int x=dFilter; x<CPLXSTEPS; x++){
		for(int y=0; y<SSTEPS; y++){
			INOUT[y][x] = 0;
		}
	}
	for(int x=0; x<dFilter; x++){
		for(int y=dFilter; y<2*dFilter; y++){
			INOUT[y][x] = 0;
		}
	}
}

// Koeffizienten fuer ETD1 Verfahren berechnen: INOUT1: exp(Lh) / INOUT2: (exp(Lh)-1)/L
void mETD1(fftw_complex INOUT1[][CPLXSTEPS], fftw_complex INOUT2[][CPLXSTEPS], fftw_complex TEMP1[][CPLXSTEPS], double dInterval, double dKappa, double deltaT){
	for(int xx=0; xx<CPLXSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			INOUT1[yy][xx] = 1;
			TEMP1[yy][xx] = 1;
		}
	}
	m2D_xDerive(INOUT1, dInterval);
	m2D_xDerive(INOUT1, dInterval);
	m2D_yDerive(TEMP1, dInterval);
	m2D_yDerive(TEMP1, dInterval);
	mCplxAdd(INOUT1, TEMP1, INOUT1); //lap in INOUT1

	mCplxCopy(INOUT1, TEMP1);
	mCplxCopy(INOUT1, INOUT2);

	m2D_xDerive(INOUT2, dInterval);
	m2D_xDerive(INOUT2, dInterval);
	m2D_yDerive(TEMP1, dInterval);
	m2D_yDerive(TEMP1, dInterval);
	mCplxAdd(INOUT2, TEMP1, INOUT2); //laplap in INOUT2

	for(int xx=0; xx<CPLXSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			INOUT2[yy][xx] *= dKappa; //kappa laplap in INOUT2
			INOUT1[yy][xx] = -1 * INOUT1[yy][xx] - INOUT2[yy][xx]; // -(lap + kaplaplap)  in  INOUT1

			if( creal(INOUT1[yy][xx])==0 & cimag(INOUT1[yy][xx])==0 ){
				INOUT2[yy][xx] = 1;
			}
			else{
				INOUT2[yy][xx] = (cexp(INOUT1[yy][xx]*deltaT) - 1.0) / INOUT1[yy][xx];
			}

			INOUT1[yy][xx] = cexp(INOUT1[yy][xx]*deltaT);
		}
	}
}

// LaPlace Operator IN-Array wird erhalten / OUT wird zu LaPlace(IN)
void mLaPlace(fftw_complex IN[][CPLXSTEPS], fftw_complex OUT[][CPLXSTEPS], fftw_complex TEMP1[][CPLXSTEPS], double dInterval){
	mCplxCopy(IN, TEMP1);
	mCplxCopy(IN, OUT);
	m2D_xDerive(OUT, dInterval);
	m2D_xDerive(OUT, dInterval);
	m2D_yDerive(TEMP1, dInterval);
	m2D_yDerive(TEMP1, dInterval);
	mCplxAdd(OUT, TEMP1, OUT);
}

// Nichtlinearitaet berechnen (OHNE LaPlace) - Orszag filter wird 2x angewandt
void mNonLin(fftw_plan R2F, fftw_plan F2R, fftw_complex IN[][CPLXSTEPS], fftw_complex NonLinOut[][CPLXSTEPS], double dNonLinTemp[][SSTEPS], fftw_complex LaPlaceTemp[][CPLXSTEPS], double dInterval){
	double dNorm = 1.0 / (SSTEPS2);
	mCplxCopy(IN, NonLinOut);
	mOrszag(NonLinOut); //ORSZAG
	fftw_execute(F2R); //NonLinOut wird zerstoert
	for(int i=0; i<SSTEPS; i++){
		for(int j=0; j<SSTEPS; j++){
			dNonLinTemp[i][j] *= dNorm;
			dNonLinTemp[i][j] *= dNonLinTemp[i][j]*dNonLinTemp[i][j];
		}
	}
	fftw_execute(R2F);
	mOrszag(NonLinOut); // ORSZAG
}

// k-Matrix fuer RK4-Verfahren berechnen
void kCalc(fftw_complex kn[][CPLXSTEPS], fftw_complex nonlin[][CPLXSTEPS], fftw_complex plain[][CPLXSTEPS], fftw_complex lap[][CPLXSTEPS], double kappa){
	for(int i=0; i<SSTEPS; i++){
		for(int j=0; j<CPLXSTEPS; j++){
			kn[i][j] = nonlin[i][j] - plain[i][j] - kappa * lap[i][j];
		}
	}
}

// Fuer RK4-Verfahren (Funktionswert fuer neuen k-Schritt aendern)
void mNewOutK(fftw_complex OUTT[][CPLXSTEPS], fftw_complex NEW[][CPLXSTEPS], fftw_complex kn[][CPLXSTEPS], double dH){
	for(int xx=0; xx<CPLXSTEPS; xx++){
		for(int yy=0; yy<SSTEPS; yy++){
			NEW[yy][xx] = OUTT[yy][xx] + dH * kn[yy][xx];
		}
	}
}

// Forward FFTW INN->OUTT
void mExecutePlan(fftw_plan myPLAN){
	fftw_execute(myPLAN);
}

// iFFTW incl. Normierung. FourierWerte werden erhalten (2d_c2r unterstuetzt den FFTW_PRESERVE_INPUT-Flag nicht)
void mExeF2R(fftw_plan myPLAN, fftw_complex OUTT[][CPLXSTEPS], fftw_complex OUTT_temp[][CPLXSTEPS], double INN[][SSTEPS]){
	mCplxCopy(OUTT, OUTT_temp);
	fftw_execute(myPLAN); //OUTT --> INN, OUTT wird vernichtet
	mCplxCopy(OUTT_temp, OUTT); //copy from temp
	D2Norm(INN); // Norm input
}



// Werte aus Signal in Datei sFilename speichern mit Bildnummer
void data_to_file(double signal[][SSTEPS], string sFilename, int iBildnummer){

	ostringstream oss;
	oss << setfill('0') << setw(4) << iBildnummer;
	sFilename += oss.str();
	ofstream myfile;
	myfile.open (sFilename.c_str());
	myfile.precision(6);
	myfile.setf( std::ios::fixed, std:: ios::floatfield );
	for(int y=0; y<SSTEPS; y++){
		for(int x=0; x<SSTEPS; x++){
			myfile << signal[x][y] << "  ";
		}
		myfile << "\n";
	}
	myfile.close();
}



// +++ MAIN +++
int main(){


	// PARAMETER ++++
	double dPhase = (2*M_PI)/(SSTEPS);
	double kappa = 1.0; // Kappa
	double deltaT = 0.01; // Zeitschrittweite
	int iTimesteps = 80000; // Anzahl Zeitschritte
	double dSleep = 0.01; // Sleep zwischen DislinFrames
	int iEachNth = 400; // Jede Nte Berechnung ausgeben
	double zMin = -1.0; // zMin
	double zMax = 1.0; // zMax

	int iChooseMethod = 1; // 1: ETD1   2: RK4


	srand (time(NULL));

	double INN[SSTEPS][SSTEPS]; // X-Space
	fftw_complex OUTT[SSTEPS][CPLXSTEPS]; // F-Space

	fftw_complex OUTT_temp[SSTEPS][CPLXSTEPS]; // Tempvar -> 2d_c2r zerstoert input
	fftw_complex LaPlaceTemp[SSTEPS][CPLXSTEPS]; // LaPlaceTemp
	fftw_complex NonLinOut[SSTEPS][CPLXSTEPS]; // TempVar fuer Nichtlinearitat im cplx
	double dNonLinTemp[SSTEPS][SSTEPS]; // TempVar fuer NL in x-space


	// FFTW fuer GesamtTrafo
	fftw_plan PLAN_R2F = fftw_plan_dft_r2c_2d(SSTEPS, SSTEPS, *INN, *OUTT, FFTW_ESTIMATE);
	fftw_plan PLAN_F2R = fftw_plan_dft_c2r_2d(SSTEPS, SSTEPS, *OUTT, *INN, FFTW_ESTIMATE);

	// FFTW fuer Trafo der NichtLinearitaet
	fftw_plan PLAN_R2F_nonlin = fftw_plan_dft_r2c_2d(SSTEPS, SSTEPS, *dNonLinTemp, *NonLinOut, FFTW_ESTIMATE);
	fftw_plan PLAN_F2R_nonlin = fftw_plan_dft_c2r_2d(SSTEPS, SSTEPS, *NonLinOut, *dNonLinTemp, FFTW_ESTIMATE);


	// Initialisierung der Startwerte
	for(int x=0; x<SSTEPS; x++){
		for(int y=0; y<SSTEPS; y++){
			//INN[x][y] = 1 * sin(1*dPhase*x) * cos(1*dPhase*y);
			//INN[x][y] = ((double)(rand()%100)/100)-0.5;	// FUNZT
//			INN[x][y] = 2*((double)(rand()%100)/100)-1;
			INN[x][y] = ((double)(rand()%100)/1000)-0.05;
		}
	}

	// R2F - INN->OUTT ++++
	mExecutePlan(PLAN_R2F);
	//+++++++++++++++++++++

	switch(iChooseMethod){
		case 1:  // ETD1
		fftw_complex ETD1_a[SSTEPS][CPLXSTEPS];
		fftw_complex ETD1_b[SSTEPS][CPLXSTEPS];

		mETD1(ETD1_a, ETD1_b, LaPlaceTemp, SSTEPS, kappa, deltaT); /// ETD1 Koeffizienten berechnen
		for(int tt=0; tt<iTimesteps; tt++){
			mNonLin(PLAN_R2F_nonlin, PLAN_F2R_nonlin, OUTT, NonLinOut, dNonLinTemp, LaPlaceTemp, SSTEPS); // NonLinOut ist ^3
			mLaPlace(NonLinOut, NonLinOut, LaPlaceTemp, SSTEPS); // laplace der Nichtlinearitaet

			// Naechsten Zeitschritt mittels ETD1
			for(int xx=0; xx<CPLXSTEPS; xx++){
				for(int yy=0; yy<SSTEPS; yy++){
					OUTT[yy][xx] = ETD1_a[yy][xx] * OUTT[yy][xx] + ETD1_b[yy][xx] * NonLinOut[yy][xx];
				}
			}

			if(tt%iEachNth == 0){
				mExeF2R(PLAN_F2R, OUTT, OUTT_temp, INN);
				data_to_file(INN, "CHVals", tt/iEachNth);
			}
		}
		break;

		case 2:	// RK4
		fftw_complex laplap[SSTEPS][CPLXSTEPS];
		fftw_complex lap[SSTEPS][CPLXSTEPS];

		fftw_complex k1[SSTEPS][CPLXSTEPS];
		fftw_complex k2[SSTEPS][CPLXSTEPS];
		fftw_complex k3[SSTEPS][CPLXSTEPS];
		fftw_complex k4[SSTEPS][CPLXSTEPS];

		// ZeitSchleife.
		for(int tt=0; tt<iTimesteps; tt++){

			mLaPlace(OUTT, lap, LaPlaceTemp, SSTEPS); //laplace
			mNonLin(PLAN_R2F_nonlin, PLAN_F2R_nonlin, OUTT, NonLinOut, dNonLinTemp, LaPlaceTemp, SSTEPS); // NonLinOut ist ^3
			kCalc(k1, NonLinOut, OUTT, lap, kappa);
			mLaPlace(k1, k1, LaPlaceTemp, SSTEPS); // 1run ist korrekt
			//k1 fertig

			mNewOutK(OUTT, OUTT_temp, k1, deltaT/2);
			mLaPlace(OUTT_temp, lap, LaPlaceTemp, SSTEPS);
			mNonLin(PLAN_R2F_nonlin, PLAN_F2R_nonlin, OUTT_temp, NonLinOut, dNonLinTemp, LaPlaceTemp, SSTEPS);
			kCalc(k2, NonLinOut, OUTT_temp, lap, kappa);
			mLaPlace(k2, k2, LaPlaceTemp, SSTEPS);
			//k2 fertig

			mNewOutK(OUTT, OUTT_temp, k2, deltaT/2);
			mLaPlace(OUTT_temp, lap, LaPlaceTemp, SSTEPS);
			mNonLin(PLAN_R2F_nonlin, PLAN_F2R_nonlin, OUTT_temp, NonLinOut, dNonLinTemp, LaPlaceTemp, SSTEPS);
			kCalc(k3, NonLinOut, OUTT_temp, lap, kappa);
			mLaPlace(k3, k3, LaPlaceTemp, SSTEPS);
			//k3 fertig

			mNewOutK(OUTT, OUTT_temp, k3, deltaT);
			mLaPlace(OUTT_temp, lap, LaPlaceTemp, SSTEPS);
			mNonLin(PLAN_R2F_nonlin, PLAN_F2R_nonlin, OUTT_temp, NonLinOut, dNonLinTemp, LaPlaceTemp, SSTEPS);
			kCalc(k4, NonLinOut, OUTT_temp, lap, kappa);
			mLaPlace(k4, k4, LaPlaceTemp, SSTEPS);
			//k4 fertig

			// naechster Zeitschritt
			for(int xx=0; xx<CPLXSTEPS; xx++){
				for(int yy=0; yy<SSTEPS; yy++){
					OUTT[yy][xx] = OUTT[yy][xx] + deltaT/6 * (k1[yy][xx] + 2*k2[yy][xx] + 2*k3[yy][xx] + k4[yy][xx]);
				}
			}

			if(tt%iEachNth == 0){
				// F2R - OUTT->INN  und  Norm +++
				mExeF2R(PLAN_F2R, OUTT, OUTT_temp, INN);
				data_to_file(INN, "CHVals", tt/iEachNth);
			}
		}
		break;

		default: cout << "Bitte Verfahren auswaehlen! (1 oder 2 bei iChooseMethod!)";
	}


	// ++++ F2R - OUTT->INN und  Norm ++++++
	mExeF2R(PLAN_F2R, OUTT, OUTT_temp, INN);
	//++++++++++++++++++++++++++++++++++++++

	// Destroy Plans
	fftw_destroy_plan(PLAN_R2F);
	fftw_destroy_plan(PLAN_R2F_nonlin);
	fftw_destroy_plan(PLAN_F2R_nonlin);
	fftw_destroy_plan(PLAN_F2R);


	return 0;
}

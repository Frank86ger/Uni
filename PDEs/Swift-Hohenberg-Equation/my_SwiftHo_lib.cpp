/**
by Frank Ehebrecht
im Mai 2014

fuer Bib:
g++ -shared -c -fPIC my_SwiftHo_lib.cpp -o my_SwiftHo_lib.o -lfftw3
g++ -shared -Wl,-soname,my_SwiftHo_lib.so -o my_SwiftHo_lib.so my_SwiftHo_lib.o -lfftw3

**/
#include <complex.h> //Dies unbedingt VOR fftw3.h einfuegen
#include <fftw3.h> //FFTW
#include <iostream>

#include <cmath>
#include <unistd.h> //fuer sleep

#include <stdlib.h> // random

#include <fstream> //Output nach txt
#include <iomanip> //Precision fuer txtOutput
#include <sstream> // txt-Output richtig nummerieren koennen

using namespace std;

#define SSTEPS 256
#define CPLXSTEPS 129 // SSTEPS/2 + 1
#define SSTEPS2 65536 // SSTEPS^2

// 128 / 256 / 512
// 65  / 129 / 257
// 16384 / 65536 / 262144


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
void mETD1(fftw_complex INOUT1[][CPLXSTEPS], fftw_complex INOUT2[][CPLXSTEPS], fftw_complex TEMP1[][CPLXSTEPS], double dInterval, double dEpsi, double deltaT){
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

			INOUT2[yy][xx] = INOUT1[yy][xx] + 1.0; // lap+1
			INOUT2[yy][xx] *= -INOUT2[yy][xx]; // -(lap+1)^2
			INOUT1[yy][xx] += dEpsi; // eps - (lap+1)^2 == linearer anteil
/****/
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


void mNameUnbekannt(fftw_complex INOUT1[][CPLXSTEPS], fftw_complex INOUT2[][CPLXSTEPS], fftw_complex TEMP1[][CPLXSTEPS], double dInterval, double dEpsi, double deltaT){
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
			INOUT1[yy][xx] = INOUT2[yy][xx] + 2*INOUT1[yy][xx] + 1.0;
			INOUT1[yy][xx] = dEpsi - INOUT1[yy][xx];
			INOUT1[yy][xx] = 1.0 / (1.0 / deltaT - INOUT1[yy][xx]);			
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
void mNonLin(fftw_plan R2F, fftw_plan F2R, fftw_complex IN[][CPLXSTEPS], fftw_complex NonLinOut[][CPLXSTEPS], double dNonLinTemp[][SSTEPS], fftw_complex LaPlaceTemp[][CPLXSTEPS], double dInterval, double dDelta, int iShape[][SSTEPS], double dDampening){
	double dNorm = 1.0 / (SSTEPS2);
	double ddtemp = 0;
	mCplxCopy(IN, NonLinOut);
	mOrszag(NonLinOut); //ORSZAG
	fftw_execute(F2R); //NonLinOut wird zerstoert
	for(int i=0; i<SSTEPS; i++){
		for(int j=0; j<SSTEPS; j++){
			dNonLinTemp[i][j] *= dNorm;
			ddtemp = dNonLinTemp[i][j];
			dNonLinTemp[i][j] = dDelta*dNonLinTemp[i][j]*dNonLinTemp[i][j] - dNonLinTemp[i][j]*dNonLinTemp[i][j]*dNonLinTemp[i][j];
			
			// Im Gebiet dampening dazunehmen
			if(iShape[i][j] == 0){
				dNonLinTemp[i][j] = dNonLinTemp[i][j] - dDampening*ddtemp;
			}
		}
	}
	fftw_execute(R2F);
	mOrszag(NonLinOut); // ORSZAG
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
void data_to_file(double signal[][SSTEPS], string sFilename, int iBildnummer, int iShape[SSTEPS][SSTEPS]){
	
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
				myfile << signal[x][y] << "  ";
			}
		}
		myfile << "\n";
	}
	myfile.close();
}



// +++ MAIN +++
extern "C" void my_SwiftHo_solver(double INN[SSTEPS][SSTEPS], int iShape[SSTEPS][SSTEPS], double dEpsilon, double deltaT, int iTimesteps, double dSleep, int iEachNth, double zMin, double zMax, int iChooseMethod, int iDislinOrTxt, double myInterval, double deltaD, double dDampening){

	
	srand (time(NULL));
	
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


	// R2F - INN->OUTT ++++
	mExecutePlan(PLAN_R2F);
	//+++++++++++++++++++++
	

	switch(iChooseMethod){
		case 1:  // Name unbekannt Methode
		fftw_complex TEMP_a[SSTEPS][CPLXSTEPS];
		fftw_complex TEMP_b[SSTEPS][CPLXSTEPS];

		mNameUnbekannt(TEMP_a, TEMP_b, LaPlaceTemp, myInterval, dEpsilon, deltaT);
		for(int tt=0; tt<iTimesteps; tt++){
			// Nichtlinearitaet incl Dampening / Shape
			mNonLin(PLAN_R2F_nonlin, PLAN_F2R_nonlin, OUTT, NonLinOut, dNonLinTemp, LaPlaceTemp, myInterval, deltaD, iShape, dDampening);

			// Naechsten Zeitschritt mittels Unbekannt
			for(int xx=0; xx<CPLXSTEPS; xx++){
				for(int yy=0; yy<SSTEPS; yy++){
					OUTT[yy][xx] = TEMP_a[yy][xx] * OUTT[yy][xx] / deltaT + TEMP_a[yy][xx] * NonLinOut[yy][xx];
				}
			}

			if(tt%iEachNth == 0){ // Ausgabe
				mExeF2R(PLAN_F2R, OUTT, OUTT_temp, INN);
				switch(iDislinOrTxt){
					case 2: //als txt output
						data_to_file(INN, "SHVals", tt/iEachNth, iShape);				
					break;
					default: cout << "Keine Ausgabemethode gewaehlt!";
				}
			}		
		}
		break;
		

		case 2:	// etd1
		fftw_complex ETD1_a[SSTEPS][CPLXSTEPS];
		fftw_complex ETD1_b[SSTEPS][CPLXSTEPS];
	
		mETD1(ETD1_a, ETD1_b, LaPlaceTemp, myInterval, dEpsilon, deltaT); /// ETD1 Koeffizienten berechnen
		for(int tt=0; tt<iTimesteps; tt++){
			// Nichtlinearitaet incl Dampening / Shape
			mNonLin(PLAN_R2F_nonlin, PLAN_F2R_nonlin, OUTT, NonLinOut, dNonLinTemp, LaPlaceTemp, myInterval, deltaD, iShape, 10.0);

			// Naechsten Zeitschritt mittels ETD1
			for(int xx=0; xx<CPLXSTEPS; xx++){
				for(int yy=0; yy<SSTEPS; yy++){
					OUTT[yy][xx] = ETD1_a[yy][xx] * OUTT[yy][xx] + ETD1_b[yy][xx] * NonLinOut[yy][xx];
				}
			}
			
			if(tt%iEachNth == 0){//Ausgabe
				mExeF2R(PLAN_F2R, OUTT, OUTT_temp, INN);
				switch(iDislinOrTxt){
					case 2: //als txt output
						data_to_file(INN, "SHVals", tt/iEachNth, iShape);				
					break;
					default: cout << "Keine Ausgabemethode gewaehlt!";
				}
			}		
		}

		break;

		default: cout << "Bitte Verfahren auswaehlen! (1 oder 2 bei iChooseMethod!)";
	}
	
	
	//mExeF2R(PLAN_F2R, OUTT, OUTT_temp, INN);
	//data_to_file(INN, "SHVals", 1, iShape);				
	
	// ++++ F2R - OUTT->INN und  Norm ++++++
	mExeF2R(PLAN_F2R, OUTT, OUTT_temp, INN);
	//++++++++++++++++++++++++++++++++++++++

	// Destroy Plans
	fftw_destroy_plan(PLAN_R2F);
	fftw_destroy_plan(PLAN_R2F_nonlin);
	fftw_destroy_plan(PLAN_F2R_nonlin);
	fftw_destroy_plan(PLAN_F2R);

}

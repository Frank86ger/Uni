/**
by Frank Ehebrecht
im Mai 2014

mit main:
g++ Burgers.cpp -o Burgers.o -lfftw3
./Burgers.o
**/
#include <complex.h> //Dies unbedingt VOR fftw3.h einfuegen
#include <fftw3.h> //FFTW
#include <iostream>

#include <cmath>
//#include "discpp.h" //Dislin
//#include <unistd.h> //fuer sleep
#include <fstream> //Output nach txt
#include <iomanip> //Precision fuer txtOutput
#include <sstream> // txt-Output richtig nummerieren koennen

using namespace std;


// Signal einlesen
void mFuncDef(fftw_complex* signal, int N, double dInterval){
	double spacing = (dInterval) / N;
	for(int i=0; i<N; i++){
		signal[i] = 1*(sin(2*(spacing * i)));
	}
}

// DFT-Signal normalisieren
void mNorm(fftw_complex* signal, int N){
	double dNorm = (1./N);
	for(int i=0; i<N; i++){
		signal[i] *= dNorm;
	}
}

// Vorwaertstrafo INCL Normierung
void mFFTW(int N, fftw_complex* INPUT, fftw_complex* OUTPUT){
	fftw_plan my_PLAN = fftw_plan_dft_1d(N,		
						INPUT,	// InputDaten
						OUTPUT,	// OutputDaten
						FFTW_FORWARD,	// Vorwaerts-Trafo
						FFTW_ESTIMATE);	// in diesem Fall schnellere Variante
	fftw_execute(my_PLAN);	// oben geplantes ausfuehren
	fftw_destroy_plan(my_PLAN);
	
	mNorm(OUTPUT, N);
}

// Ruecktrafo OHNE Normierung
void mIFFTW(int N, fftw_complex* INPUT, fftw_complex* OUTPUT){
	fftw_plan my_PLAN = fftw_plan_dft_1d(N,		
						INPUT,	// InputDaten
						OUTPUT,	// OutputDaten
						FFTW_BACKWARD,	// Rueck-Trafo
						FFTW_ESTIMATE);	// in diesem Fall schnellere Variante
	fftw_execute(my_PLAN);	// oben geplantes ausfuehren		
	fftw_destroy_plan(my_PLAN);
}

// Stuetzpunkte, k-Wert, Input-Arr, Output-Arr
void mFDerive(int N, double k, fftw_complex* INPUT, fftw_complex* OUTPUT){
	for(int n=0; n<N; n++){
		OUTPUT[n] = I * n * k * INPUT[n];
	}
	for(int n=N/2+1; n<N; n++){
		OUTPUT[n] = I * (n-N) * k * INPUT[n];
	}
}

// Zwei Arrays elementweise Multiplizieren
void mMultiply(int N, fftw_complex* INPUT1, fftw_complex* INPUT2, fftw_complex* OUTPUT){
	for(int n=0; n<N; n++){
		OUTPUT[n] = INPUT1[n] * INPUT2[n];
	}
}

// RK4 - k-Wert-Schritt
void mRK_k_calc(int N, double dFac, double h, fftw_complex* INPUT1, fftw_complex* INPUT2, fftw_complex* OUTPUT){
	for(int n=0; n<N; n++){
		OUTPUT[n] = INPUT1[n] + dFac * h * INPUT2[n];
	}
}

// Die Funktion um den naechsten Zeitschritt zu berechnen
void mRK4_func(int N, double k, double dNU,fftw_complex* INPUT, fftw_complex* OUTPUT){
	for(int n=0; n<N/2; n++){
		OUTPUT[n] = -dNU * n*n * k*k * INPUT[n] + OUTPUT[n];
	}
	for(int n=N/2+1; n<N; n++){
		OUTPUT[n] = -dNU * (n-N)*(n-N) * k*k * INPUT[n] + OUTPUT[n];
	}
}

// Orszag Dealiase
void mDealiase(int N, fftw_complex* INPUT, fftw_complex* OUTPUT){
	int nmax = trunc((double)N/3); // N/2 * 2/3 -> 2/3 Regel
	for(int n=0; n<=nmax; n++){
		OUTPUT[n] = INPUT[n];
	}
	for(int n=nmax+1; n<(N-nmax); n++){
		OUTPUT[n] = 0;
	}
	for(int n=(N-nmax); n<N; n++){
		OUTPUT[n] = INPUT[n];
	}
}

// 36-Filter Dealiase
void mDealiase2(int N, fftw_complex* INPUT, fftw_complex* OUTPUT){
	for(int n=0; n<=N/2; n++){
		OUTPUT[n] = INPUT[n] * exp(-36*pow((n/(N/2)),36));
	}
	double kmax=N/2+1;
	for(int n=N/2+1; n<N; n++){
		OUTPUT[n] = INPUT[n] * exp(-36*pow(((2*kmax-n)/kmax),36));
	}
}

// Array kopieren
void mCopy(int N, fftw_complex* INPUT, fftw_complex* OUTPUT){
	for(int n=0; n<N; n++){
		OUTPUT[n] = INPUT[n];
	}
}

// Array elementweise mit Faktor multiplizieren
void mMultFak(int N, double dFak, fftw_complex* INOUT){
	for(int n=0; n<N; n++){
		INOUT[n] *= dFak;
	}
}

void data_to_file(double* signal, string sFilename, int iBildnummer, int iSteps, double dInterval){
	
	ostringstream oss;
	oss << setfill('0') << setw(4) << iBildnummer;
	sFilename += oss.str();
	ofstream myfile;
	myfile.open (sFilename.c_str());
	myfile.precision(6);
	myfile.setf( std::ios::fixed, std:: ios::floatfield );
	for(int x=0; x<iSteps; x++){
		myfile << dInterval / iSteps * x << "\t"  << signal[x] << "\n";
	}
	myfile.close();
}


// Pseudospektralverfahren mit RK4
void mPseudoSpec(fftw_complex* XSpace, fftw_complex* FSpace, int N, double dInterval, double dDeltaT, int iTimesteps, double dANiSpeed, double dNU, int* iCases, double* PltRg){

	fftw_complex FSpace_temp[N]; //temp
	fftw_complex XSpace_temp[N]; //temp
	
	double k = (2 * M_PI) / (dInterval); //k-vec
	complex double u_temp[N]; //temp
	complex double u_temp_dea[N]; //temp
	complex double RK_kvec[5][N]; //k1,2,3,4
	double RK_facs[4] = {0.0, 0.5, 0.5, 1.0}; //k1,2,3,4-Faktoren

	double aREAL[2][N]; //Realteil

	// x-skalierung
	double aSPACE[N];
	for(int x=0; x<N; x++){
		aSPACE[x] = x * (dInterval)/N;
	}	

	// Haette man sich auch sparen koennen
	complex double cplx_temp[N];
	complex double cplx_temp2[N];
	complex double cplx_temp3[N];
	

	for(int tt=0; tt<iTimesteps; tt++){ //fuer alle Zeitschritte
		for(int kn=0; kn<4; kn++){ // fuer k1,k2,k3,k4
			mRK_k_calc(N, RK_facs[kn], dDeltaT, FSpace, RK_kvec[kn], u_temp); // u_temp ist das neue signal


			switch(iCases[0]){ // Erstes Dealiasing
				case 0:
					// u_temp = u_temp_dea -> kein Dealiasing
					mCopy(N, u_temp, u_temp_dea);
					break;
				case 1:
					// Orszag Filter
					mDealiase(N, u_temp, u_temp_dea);
					break;
				case 2:
					// 36-Filter
					mDealiase2(N, u_temp, u_temp_dea);
					break;
			}


			switch(iCases[1]){ // u del u oder .5 del u^2
				case 0:
					// u del u
					mFDerive(N, k, u_temp_dea, cplx_temp); //Ableitung
					mIFFTW(N, cplx_temp, cplx_temp2); 
					mIFFTW(N, u_temp_dea, cplx_temp3);
					mMultiply(N, cplx_temp2, cplx_temp3, cplx_temp);
					mFFTW(N, cplx_temp, RK_kvec[kn+1]);
					break;
				case 1:
					// .5 del (u^2)
					mIFFTW(N, u_temp_dea, cplx_temp); // u* -> u
					mMultiply(N, cplx_temp, cplx_temp, cplx_temp3); // u^2
					mFFTW(N, cplx_temp3, cplx_temp); // (u^2)*
					mFDerive(N, k, cplx_temp, RK_kvec[kn+1]); //Ableitung
					mMultFak(N, 0.5, RK_kvec[kn+1]);
					break;
			}




			switch(iCases[2]){ // Zweites Dealiasing
				case 0:
					// kein Aliasing
					break;
				case 1:
					// Orszag-Filter
					mDealiase(N, RK_kvec[kn+1], RK_kvec[kn+1]);
					break;
				case 2:
					// 36-Filter
					mDealiase2(N, RK_kvec[kn+1], RK_kvec[kn+1]);
					break;
			}


			mRK4_func(N, k, dNU, u_temp, RK_kvec[kn+1]); //k-Vec berechnen
		}

	
		// RK4 Step ausfuehren
		for(int n=0; n<N; n++){
			FSpace_temp[n] = FSpace[n] + dDeltaT/6 * (RK_kvec[1][n] + 2*RK_kvec[2][n] + 2*RK_kvec[3][n] + RK_kvec[4][n]);
		}
	
		
		//Ruecktrafo nach XSpace
		mIFFTW(N, FSpace_temp, XSpace_temp);
		
		
		// Realteil aus komplexen Signal kopieren ... unschoen
		for(int x=0; x<N; x++){
			aREAL[0][x] = creal(XSpace[x]); //Realteil
			aREAL[1][x] = creal(XSpace_temp[x]); //Realteil
		}



		// Plot +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		if(tt%10 == 0){
			data_to_file(aREAL[0], "data", tt/10, N, dInterval);
		}
		// Plot -------------------------------------------------------------
		
		// [1] -> [0]
		for(int i=0; i<N; i++){
			FSpace[i] = FSpace_temp[i];
			XSpace[i] = XSpace_temp[i];
		}
	}
}
/**
// my_burgers_solver mittels Python aufrufen
extern "C" void my_burgers_solver(fftw_complex *XSpace, int iSpacesteps, int iTimesteps, double dInterval, double dNU, double dAniSpeed, double dDeltaT, int iCases[3], double dPlotRanges[6]) {

	fftw_complex FSpace[iSpacesteps];

	// Vorwaerts-Trafo incl Normierung fuer erstes Element
	mFFTW(iSpacesteps, XSpace, FSpace);

	// Pseudospektral-Verfahren
	mPseudoSpec(XSpace, FSpace, iSpacesteps, dInterval, dDeltaT, iTimesteps, dAniSpeed, dNU, iCases, dPlotRanges);
}
**/

/****/
// +++ MAIN +++
int main(){
	int iSpacesteps = 512;
	int iTimesteps = 1400;
	double dInterval = 2 * M_PI;
	double dNU = 0.005;
	double dAniSpeed = 1;
	double dDeltaT = 0.001;
	int iCases[3] = {1,1,1};
	double dPlotRanges[6] = {0.0, 2*M_PI, 1.0, -1.0, 1.0, 0.5};



	// XSpace: Ortsraum; FSpace: Freqraum
	fftw_complex XSpace[iSpacesteps];
	fftw_complex FSpace[iSpacesteps];
	
	// Startsignal
	mFuncDef(XSpace, iSpacesteps, dInterval);

	// Vorwaerts-Trafo incl Normierung fuer erstes Element
	mFFTW(iSpacesteps, XSpace, FSpace);



	// RK4 im Freqraum
	mPseudoSpec(XSpace, FSpace, iSpacesteps, dInterval, dDeltaT, iTimesteps, dAniSpeed, dNU, iCases, dPlotRanges);

	return 0;
}

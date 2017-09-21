/**
1D Allen-Cahn-Gleichung mit Pseudospektralverfahren und ETD1
by Frank Ehebrecht
im Mai 2014

g++ AllenCahn.cpp -o AllenCahn.o -lfftw3
./AllenCahn.o

=> Output in Form von Textdateien : "blahXXXX"
**/
#include <complex.h> //Dies unbedingt VOR fftw3.h einfuegen
#include <fftw3.h> //FFTW
#include <iostream>

#include <cmath>
//#include "discpp.h" //Dislin
#include <unistd.h> //fuer sleep

#include <fstream> //Output nach txt
#include <iomanip> //Precision fuer txtOutput
#include <sstream> // txt-Output richtig nummerieren koennen

using namespace std;


// Signal einlesen
void mFuncDef(fftw_complex* signal, int N, double dInterval, double dBeta){
	double spacing = (dInterval) / N;
	
	/**
	for(int i=0; i<N; i++){
		signal[i] = 0.05 * exp(-5.0* ((spacing * i)-0)*((spacing * i)-0));
		signal[i] = sin(spacing*i);// WEG
	}
	**/
	for(int i=0; i<=(double)N/3; i++){
		signal[i] = dBeta;
	}
	for(int i=(int)N/3; i<=N*2./3; i++){
		signal[i] = 1.0;
	}
	for(int i=(int)N*2./3; i<=N; i++){
		signal[i] = dBeta;
	}

	
}

// DFT-Signal normalisieren
void mNorm(fftw_complex* signal, int N){
	double dNorm = (1./N);
	for(int i=0; i<N; i++){
		signal[i] *= dNorm;
	}
}

// Stuetzpunkte, k-Wert, Input-Arr, Output-Arr
void mFDerive(int N, double k, fftw_complex* INPUT, fftw_complex* OUTPUT){
	for(int n=0; n<=N/2; n++){
		OUTPUT[n] = I * n * k * INPUT[n];
	}
	for(int n=N/2+1; n<N; n++){
		OUTPUT[n] = I * (n-N) * k * INPUT[n];
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
void mCplxCopy(int N, fftw_complex* INPUT, fftw_complex* OUTPUT){
	for(int n=0; n<N; n++){
		OUTPUT[n] = INPUT[n];
	}
}

// ETD1 Koeffizienten berechnen
void mETD1(fftw_complex* INOUT1, fftw_complex* INOUT2, int iSteps, double dK, double dKappa, double dT, double dBeta){
	for(int i=0; i<iSteps; i++){
		INOUT1[i] = 1.0;
	}
	mFDerive(iSteps, dK, INOUT1, INOUT1);
	mFDerive(iSteps, dK, INOUT1, INOUT1); //laplap
	for(int i=0; i<iSteps; i++){
		INOUT1[i] *= dKappa; //kap laplap
		INOUT1[i] = INOUT1[i] - dBeta; //(kap laplap - beta)

		if(creal(INOUT1[i])==0 & cimag(INOUT1[i])==0){
			INOUT2[i] = 1;
		}else{
			INOUT2[i] = (cexp(INOUT1[i]*dT) - 1) / 	INOUT1[i];
		}	
		INOUT1[i] = cexp(INOUT1[i]*dT);
	}
}

// FFTW Forward
void mExecuteForward(fftw_plan myPLAN){
	fftw_execute(myPLAN);
}

// FFTW Backward incl Normierung
void mExecuteBackward(fftw_plan myPLAN, fftw_complex* XSpace, int iSteps){
	fftw_execute(myPLAN);
	mNorm(XSpace, iSteps); // fuer reell bauen
}

// Nichtlinearen Term berechnen
void mNonLin(fftw_plan PLAN_X2F_nonlin,fftw_plan PLAN_F2X_nonlin, fftw_complex* FSpace, fftw_complex* FNonLin, fftw_complex* XNonLin, fftw_complex* FNonLin2, fftw_complex* XNonLin2, int iSpacesteps){
	double dNorm = 1.0 / iSpacesteps;
	mCplxCopy(iSpacesteps, FSpace, FNonLin);
	// Dealiasing hier? .. unnoetig
	fftw_execute(PLAN_F2X_nonlin);// FNonLin -> XNonLin
	
	for(int iii=0; iii<iSpacesteps; iii++) { XNonLin2[iii] = XNonLin[iii]; } //temp in xnonlin2
	
	
	for(int i=0; i<iSpacesteps; i++){
		XNonLin[i] *= dNorm; //Normalisieren nicht F sondern X!!
		XNonLin[i] *= XNonLin[i]; //Quadrieren
	}
	// nochma Dealiasing hier ? ... unnoetig
	fftw_execute(PLAN_X2F_nonlin);// XNonLin -> FNonLin
	mCplxCopy(iSpacesteps, FNonLin, FNonLin2); //fourier quad in FNonLin2
	
	
	for(int i=0; i<iSpacesteps; i++){
		XNonLin2[i] *= dNorm; //Normalisieren nicht F sondern X!!
		XNonLin[i] = XNonLin2[i] * XNonLin2[i] * XNonLin2[i];
	}
	fftw_execute(PLAN_X2F_nonlin);// XNonLin -> FNonLin
	
	
}

// Werte aus Signal in Datei sFilename speichern mit Bildnummer
void data_to_file(double* signal, string sFilename, int iBildnummer, int iSteps, double dInterval){
	
	ostringstream oss;
	oss << setfill('0') << setw(4) << iBildnummer;
	sFilename += oss.str();
	ofstream myfile;
	myfile.open (sFilename.c_str());
	myfile.precision(6);
	myfile.setf( std::ios::fixed, std:: ios::floatfield );
	for(int x=0; x<iSteps; x++){
		myfile << dInterval / iSteps * x -20 << "\t"  << signal[x] << "\n";
	}
	myfile.close();
}

// +++ MAIN +++
int main(){
	
	int iSpacesteps = 256;
	double dInterval = 60;
	double dDeltaT = 0.05;
	double dKappa = 1;
	int iTimeSteps = 8000;
	double dD = 1;
	double dBeta = -1.1;
	int eachNth = 10;

	double dK = (2 * M_PI) / dInterval;

	double Realout[iSpacesteps];

	// XSpace: Ortsraum; FSpace: Freqraum
	fftw_complex XSpace[iSpacesteps];
	fftw_complex FSpace[iSpacesteps];
	
	fftw_complex XNonLin[iSpacesteps];
	fftw_complex FNonLin[iSpacesteps];
	
	fftw_complex XNonLin2[iSpacesteps];
	fftw_complex FNonLin2[iSpacesteps];
	
	
	
	fftw_complex FTemp[iSpacesteps];
	
	fftw_plan PLAN_X2F = fftw_plan_dft_1d(iSpacesteps, XSpace, FSpace, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_F2X = fftw_plan_dft_1d(iSpacesteps, FSpace, XSpace, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_plan PLAN_X2F_nonlin = fftw_plan_dft_1d(iSpacesteps, XNonLin, FNonLin, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_F2X_nonlin = fftw_plan_dft_1d(iSpacesteps, FNonLin, XNonLin, FFTW_BACKWARD, FFTW_ESTIMATE);


	// Funktion definieren
	mFuncDef(XSpace, iSpacesteps, dInterval, dBeta);

	// XSpace -> FSpace
	mExecuteForward(PLAN_X2F);

	// Faktoren fuer ETD1 berechnen
	fftw_complex ETD1_a[iSpacesteps];
	fftw_complex ETD1_b[iSpacesteps];
	mETD1(ETD1_a, ETD1_b, iSpacesteps, dK, dKappa, dDeltaT, dBeta);

	//for(int kk=0; kk<256; kk++){ cout << cimag(ETD1_b[kk]) << "   "; } //WEG
	
	/** **/
	// Zeitschritte
	for(int tt=0; tt<iTimeSteps; tt++){
	
		mNonLin(PLAN_X2F_nonlin, PLAN_F2X_nonlin, FSpace, FNonLin, XNonLin, FNonLin2, XNonLin2, iSpacesteps); //Nichtlinearitaet
		

		// Naechster Zeitschritt
		for(int i=0; i<iSpacesteps; i++){
			FSpace[i] = ETD1_a[i] * FSpace[i] + ETD1_b[i] * ((1+dBeta)*FNonLin2[i] - FNonLin[i]);
		}

		
		if(tt%eachNth == 0){ // Output
			mExecuteBackward(PLAN_F2X, XSpace, iSpacesteps);
			
			for(int i=0; i<iSpacesteps; i++){
				Realout[i] = creal(XSpace[i]);
			}
			data_to_file(Realout, "blah", tt/eachNth, iSpacesteps, dInterval); //Output an File "blahXXXX"
		}
	}

	fftw_destroy_plan(PLAN_X2F);
	fftw_destroy_plan(PLAN_F2X);
	fftw_destroy_plan(PLAN_X2F_nonlin);
	fftw_destroy_plan(PLAN_F2X_nonlin);
	
	return 0;
}

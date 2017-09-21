/**
by Frank Ehebrecht im Juli 2014

compile with:
g++ -shared -c -fPIC my_GLG_txtout.cpp -o my_GLG_txtout.o -lfftw3
g++ -shared -Wl,-soname,my_GLG_txtout.so -o my_GLG_txtout.so my_GLG_txtout.o -lfftw3

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

//2D Signal Normalisieren
void mD2Norm(fftw_complex* INOUT, int iSteps){
	double dNorm = 1.0 / (iSteps);
	for(int i=0; i<iSteps; i++){
		INOUT[i] *= dNorm;
	}
}

//x-Ableitung berechnen
void mDerive(fftw_complex* INOUT, double dInterval, int iSteps){
	double dK = 2*M_PI / dInterval;
	int iHalfsteps = iSteps / 2 + 1;

	for(int xx=0; xx<iHalfsteps; xx++){
		INOUT[xx] *= I * dK * xx;
	}
	for(int xx=iHalfsteps; xx<iSteps; xx++){
		INOUT[xx] *= -1 * I * dK * (iSteps-xx);
	}
}

void mOrzag(fftw_complex* INOUT, int iSteps){
	int iCut = (int)(iSteps / 2 * 2 / 3);
	for(int i=iCut; i<2*iCut; i++){
		INOUT[i] = 0;
	}
}


// ETD2-Koeffizienten. a: exp(qh) / b fuer N_n / c fuer N_(n-1)
void mETD2(fftw_complex* F_ETD2_a, fftw_complex* F_ETD2_b, fftw_complex* F_ETD2_c, double dAlpha, double deltaT, double dInterval, int iSteps){
	for(int xx=0; xx<iSteps; xx++){
		F_ETD2_b[xx] = 1;
	}	
	mDerive(F_ETD2_b, dInterval, iSteps);
	mDerive(F_ETD2_b, dInterval, iSteps); //deldel y in b
	
	for(int xx=0; xx<iSteps; xx++){
		F_ETD2_b[xx] = F_ETD2_b[xx] * (1.0 + I * dAlpha) + 1.0; // ((1 + ia)lap) + 1 = q   in b
		F_ETD2_a[xx] = cexp(F_ETD2_b[xx]*deltaT); //exp(qh)   in a
		if( creal(F_ETD2_b[xx])==0 & cimag(F_ETD2_b[xx])==0 ){
			F_ETD2_b[xx] = 3.0 * deltaT / 2.0;
			F_ETD2_c[xx] = -1.0 * deltaT / 2.0;
		}
		else{
			F_ETD2_c[xx] = (1.0 + deltaT * F_ETD2_b[xx] - F_ETD2_a[xx]) / (deltaT * F_ETD2_b[xx] * F_ETD2_b[xx]);
			F_ETD2_b[xx] = ((1.0+deltaT*F_ETD2_b[xx]) * F_ETD2_a[xx] - 1.0 - 2.0 * deltaT * F_ETD2_b[xx]) / (deltaT * F_ETD2_b[xx] * F_ETD2_b[xx]);
		}
	}
}

// Nichtlinearitaet berechnen / F_Signal_nl wird ist Ausgabegroesse (Nichtlinearitaet in F)
void mNonLin(fftw_plan PLAN_F2X, fftw_plan PLAN_X2F_nl, fftw_complex* X_Signal_nl, fftw_complex* X_Signal, fftw_complex* F_Signal_nl, double dBeta, int iSteps){ 

	fftw_execute(PLAN_F2X); //X_Signal liegt vor
	mD2Norm(X_Signal, iSteps); //X_Signal normiert

	//Nicht-Linearitaet
	for(int xx=0; xx<iSteps; xx++){
		X_Signal_nl[xx] = -1.0 * X_Signal[xx] * conj(X_Signal[xx]) * X_Signal[xx] * (1 + I*dBeta);
	}
	fftw_execute(PLAN_X2F_nl);
	// + Orszag auf non F?
	mOrzag(F_Signal_nl, iSteps);
}

// INPUT nach OUTPUT kopieren
void mCopy(fftw_complex* INPUT, fftw_complex* OUTPUT, int iSteps){
	for(int ii=0; ii<iSteps; ii++){
		OUTPUT[ii] = INPUT[ii];
	}
}

// Realteil von Signal in Datei sFilename speichern mit Bildnummer
void m_R_data_to_file(fftw_complex* signal, string sFilename, int iBildnummer, int iSteps){
	
	string sFilename2 = sFilename + "matrix";
	
	ostringstream oss;
	oss << setfill('0') << setw(4) << iBildnummer;
	sFilename += oss.str();
	ofstream myfile;
	myfile.open (sFilename.c_str());
	myfile.precision(6);
	myfile.setf( std::ios::fixed, std:: ios::floatfield );
	for(int xx=0; xx<iSteps; xx++){
		myfile << xx << " \t" << creal(signal[xx]) << "  ";
		myfile << "\n";
	}
	myfile.close();
	
	ofstream mymatrix;

	mymatrix.open (sFilename2.c_str(), ios::app);
	mymatrix.precision(6);
	mymatrix.setf( std::ios::fixed, std::ios::floatfield );
	for(int xx=0; xx<iSteps; xx++){
		mymatrix << creal(signal[xx]) << "  ";
	}
	mymatrix << "\n";
	mymatrix.close();
}

// Betrag von Signal in Datei sFilename speichern mit Bildnummer
void m_A_data_to_file(fftw_complex* signal, string sFilename, int iBildnummer, int iSteps){

	string sFilename2 = sFilename + "matrix";
	
	ostringstream oss;
	oss << setfill('0') << setw(4) << iBildnummer;
	sFilename += oss.str();
	ofstream myfile;
	myfile.open (sFilename.c_str());
	myfile.precision(6);
	myfile.setf( std::ios::fixed, std:: ios::floatfield );
	for(int xx=0; xx<iSteps; xx++){
		myfile << xx << " \t" << cabs(signal[xx]) << "  ";
		myfile << "\n";
	}
	myfile.close();
	
	ofstream mymatrix;

	mymatrix.open (sFilename2.c_str(), ios::app);
	mymatrix.precision(6);
	mymatrix.setf( std::ios::fixed, std::ios::floatfield );
	for(int xx=0; xx<iSteps; xx++){
		mymatrix << cabs(signal[xx]) << "  ";
	}
	mymatrix << "\n";
	mymatrix.close();
}

// Phase von Signal in Datei sFilename speichern mit Bildnummer
void m_P_data_to_file(fftw_complex* signal, string sFilename, int iBildnummer, int iSteps){
	
	string sFilename2 = sFilename + "matrix";
	
	ostringstream oss;
	oss << setfill('0') << setw(4) << iBildnummer;
	sFilename += oss.str();
	ofstream myfile;
	myfile.open (sFilename.c_str());
	myfile.precision(6);
	myfile.setf( std::ios::fixed, std:: ios::floatfield );
	for(int xx=0; xx<iSteps; xx++){
		if(creal(signal[xx]) > 0){
			myfile << xx << " \t" << atan(cimag(signal[xx]) / creal(signal[xx])) << "  ";
		}
		else if(creal(signal[xx]) < 0.0 && cimag(signal[xx]) >= 0.0){
			myfile << xx << " \t" << atan(cimag(signal[xx]) / creal(signal[xx])) + M_PI << "  ";
		}
		else if(creal(signal[xx]) < 0.0 && cimag(signal[xx]) < 0.0){
			myfile << xx << " \t" << atan(cimag(signal[xx]) / creal(signal[xx])) - M_PI << "  ";
		}
		else if(creal(signal[xx]) == 0.0 && cimag(signal[xx]) > 0.0){
			myfile << xx << " \t" << M_PI / 2 << "  ";
		}
		else if(creal(signal[xx]) == 0.0 && cimag(signal[xx]) < 0.0){
			myfile << xx << " \t" << -1.0 * M_PI / 2 << "  ";
		}
		else if(creal(signal[xx]) == 0.0 && cimag(signal[xx]) == 0.0){
			myfile << xx << " \t" << 0.0 << "  ";
		}
		myfile << "\n";
	}
	myfile.close();
}


//int main(){
extern "C" void my_GLG_solver(double IN_real[1024], double IN_imag[1024], double dAlpha, double dBeta, int iTimeSteps, double deltaT, int iEachNth, double dInterval, int iSteps){


//	double dInterval = 2 * M_PI;
//	double dAlpha = 1;
//	double dBeta = 2;
//	double deltaT = 0.05;
//	int iSteps = 512;
	
//	int iEachNth = 20;
	
	srand (time(NULL));
	
//	int iTimeSteps = 1000;
		  
	fftw_complex *X_Signal, *F_Signal, *X_Signal_nl, *F_Signal_nl, *F_ETD2_a, *F_ETD2_b, *F_ETD2_c, *F_Signal_nl_temp;
		  
	X_Signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * iSteps);
	F_Signal = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * iSteps);
	
	X_Signal_nl = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * iSteps);
	F_Signal_nl = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * iSteps);
	
	F_Signal_nl_temp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * iSteps);
	
	F_ETD2_a = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * iSteps);
	F_ETD2_b = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * iSteps);
	F_ETD2_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * iSteps);
				
	
	fftw_plan PLAN_X2F = fftw_plan_dft_1d(iSteps, X_Signal, F_Signal, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_F2X = fftw_plan_dft_1d(iSteps, F_Signal, X_Signal, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	fftw_plan PLAN_X2F_nl = fftw_plan_dft_1d(iSteps, X_Signal_nl, F_Signal_nl, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_plan PLAN_F2X_nl = fftw_plan_dft_1d(iSteps, F_Signal_nl, X_Signal_nl, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	// INPUTSIGNAL //////////////////
	for(int xx=0; xx<iSteps; xx++){
		X_Signal[xx] = IN_real[xx] + I * IN_imag[xx];		
		//X_Signal[xx] = 0 + ((double)(rand()%1000)/10000)-0.05;
	}
	/////////////////////////////////

	fftw_execute(PLAN_X2F); // FFT
	mETD2(F_ETD2_a, F_ETD2_b, F_ETD2_c, dAlpha, deltaT, dInterval, iSteps); //ETD2 - Koeffizienten berechnen
	mNonLin(PLAN_F2X, PLAN_X2F_nl, X_Signal_nl, X_Signal, F_Signal_nl, dBeta, iSteps);

	m_R_data_to_file(X_Signal, "GLGreal", 0, iSteps);
	m_A_data_to_file(X_Signal, "GLGabs", 0, iSteps);
	m_P_data_to_file(X_Signal, "GLGphas", 0, iSteps);

	for(int xx=0; xx<iSteps; xx++){
		F_Signal[xx] = F_ETD2_a[xx] * F_Signal[xx] + F_ETD2_b[xx] * F_Signal_nl[xx] + F_ETD2_c[xx] * F_Signal_nl[xx];
	}
	mCopy(F_Signal_nl, F_Signal_nl_temp, iSteps);



	for(int tt=1; tt<iTimeSteps; tt++){
		
		// output
		m_R_data_to_file(X_Signal, "GLGreal", tt/iEachNth, iSteps);
		m_A_data_to_file(X_Signal, "GLGabs", tt/iEachNth, iSteps);
		m_P_data_to_file(X_Signal, "GLGphas", tt/iEachNth, iSteps);
		
		mNonLin(PLAN_F2X, PLAN_X2F_nl, X_Signal_nl, X_Signal, F_Signal_nl, dBeta, iSteps);
		for(int xx=0; xx<iSteps; xx++){
		F_Signal[xx] = F_ETD2_a[xx] * F_Signal[xx] + F_ETD2_b[xx] * F_Signal_nl[xx] + F_ETD2_c[xx] * F_Signal_nl[xx];
		}
		mCopy(F_Signal_nl, F_Signal_nl_temp, iSteps);
	}

	fftw_execute(PLAN_F2X);
	mD2Norm(X_Signal, iSteps);	
		
//	m_R_data_to_file(X_Signal, "endout", 2, iSteps);
	
	
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

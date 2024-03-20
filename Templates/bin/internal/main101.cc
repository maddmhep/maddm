#include "Pythia8/Pythia.h"
#include "math.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;
using namespace Pythia8;
using std::ifstream;

int const   r_vecsize  = 200;
int const   p_vecsize  = 200;
float       fArgonne[int(r_vecsize)*int(p_vecsize)+10];
float       fArgonne_r[int(r_vecsize)*int(p_vecsize)+10];
float       fArgonne_q[int(r_vecsize)*int(p_vecsize)+10];
double const mp=0.9382720813;
double const mn=0.93956542052;
double const mD=1.8756;
double const D_m=1.8756;

//void writeArgonne();

void writeArgonne(){
  float a1, a2, a3, a4;
  
  std::ifstream infile("../../../Cards/Argonne_WF.txt");
  std::string line; std::getline(infile, line);
  
  int count = 0;
  while (infile >> a1 >> a2 >> a3 >> a4)
    {	
      fArgonne[count] = a4+a3;
	  //fArgonne_r[count] = a1;
	  //fArgonne_q[count] = a2;
      count++;
    }
	//return fArgonne
}

void writeDbarspherical(double pcoalescence, double mDM, std::string filename, int nbins, double A, double Z){
    //writeDbarspherical(pcoalescence,mDM,outdir,nbins);
	
    float dNdE_pbar[int(nbins)];
	float T_pbar[int(nbins)];
	float log10x_pbar[int(nbins)];
	
    double eminh = -9.;
    double emaxh = 0.;
    //double nbins = 180; //Binninb modifed to account properly for the line spectra. The last bin is -0.05/-0.025/0.00+
    double DeltaBin = (emaxh-eminh)/nbins; 
    double T_pbar_min = pow(10.,eminh+DeltaBin/2.)*mDM;
    double T_pbar_max = pow(10.,emaxh-DeltaBin/2.)*mDM;
	double DeltaTBin = (log10(T_pbar_max)-log10(T_pbar_min))/nbins;
	
	std::string inputFilePath = filename+"antiprotonsPrimordial_spectrum_pythia8.dat";
	std::string outputFilePath = filename+"antideuterons_spherical_spectrum_pythia8.dat";
	
	if(A==2 && Z==1){
		outputFilePath = filename+"antideuterons_spherical_spectrum_pythia8.dat";
	}
	else if(A==3 && Z==2){
		outputFilePath = filename+"antihelions3_spherical_spectrum_pythia8.dat";
	}
	else if(A==4 && Z==2){
		outputFilePath = filename+"antihelions4_spherical_spectrum_pythia8.dat";
	}
	else{
		outputFilePath = filename+"antinucleus_spherical_spectrum_pythia8.dat";
	}
	
	double N = A-Z;
	double mA = mp*Z+mn*N;
		
	std::ifstream inputFile(inputFilePath);
	std::ofstream outputFile(outputFilePath);
	//std::ifstream infile(filename.c_str());
	
    // Read data from the input file and write to the output file
    double column1, column2;
	int count = 0;
	double log10x = 0.; 
    while (inputFile >> column1 >> column2) {
		log10x_pbar[count] = column1;
		T_pbar[count] = pow(10.,column1)*mDM;
		dNdE_pbar[count] = column2/(log(10)*T_pbar[count]);
		count++;
    }
	
	double T_A = 0.;
	double Tpbar_val = 0.;
	int binint = 0.;
	double dNdPbar_TA = 0.;
	double dNdT_A = 0.;
	double BA = 0.;
	double K_A = 0.;
	double dNdlog10x_A =0.;
	for(int t=0; t<nbins ; t++){
		T_A = T_pbar[t];
		Tpbar_val = T_A/A;
		
		if(Tpbar_val<T_pbar_min){
			dNdlog10x_A = 0.;
		    dNdPbar_TA = 0.;
		}
		else{
			K_A = sqrt( (T_A+mA)*(T_A+mA) - mA*mA );
			//BA = (mp*Z+mn*(A-Z))*pow(pcoalescence,3.)/(6.*mp*Z*mn*(A-Z)*K_Dbar);
			BA = ((mp*Z+mn*(A-Z))/(mp*Z*mn*(A-Z))) * pow(A,A) * pow(pow(pcoalescence,3.)/(24.*K_A),A-1.);
			binint = round( (log10(Tpbar_val)-log10(T_pbar_min))/DeltaTBin );
			dNdPbar_TA =  dNdE_pbar[binint-1] + (Tpbar_val-T_pbar[binint-1])/(T_pbar[binint]-T_pbar[binint-1])*(dNdE_pbar[binint]-dNdE_pbar[binint-1]);
		}
		
		dNdlog10x_A = BA*pow(dNdPbar_TA,A)*log(10.)*T_A;
		outputFile << std::setprecision(5) << std::scientific << log10(T_A/mDM) << "\t" << dNdlog10x_A << std::endl;
	}
	
    // Close the files
    inputFile.close();
    outputFile.close();
}


double WF_probability(double sigma,double d,double q) {
	double pigreek = 3.141592;
	double Prob = 8.*pow(d*d/(d*d+sigma*sigma*4),3./2.)*exp(-q*q*d*d*5.0684*5.0684/4.);
	return Prob;
}

double WF_SG_probability_pvalue(double Deltar, double Deltap, double sigma, double d) {
	//cout << "In function " << sigma << "  " << d << endl;
	double pigreek = 3.141592;
	double d_GeV = d*5.0684; //GeV^-1
    double funcr = (1.-erf(Deltar/(sqrt(2.)*sigma)));
    double funcq = (1.-erf(Deltap*d_GeV/sqrt(2.)));
	double Prob = funcr*funcq;
	return Prob;
}

double WF_SG_probability_PDF(double Deltar, double Deltap, double sigma, double d) {
	double pigreek = 3.141592;
	double d_GeV = d*5.0684; //GeV^-1
	double normalization = 4. * d_GeV/sigma / (2.*pigreek);
    double funcr = exp( - pow( Deltar/(sqrt(2.)*sigma ) ,2. ) );
    double funcq = exp( - pow( Deltap*d_GeV/(sqrt(2.)) ,2. ) );
	double Prob = funcr*funcq*normalization;
	return Prob;
}

double WF_Argonne_probability(double Deltar, double Deltap) {
	double pigreek = 3.141592;
	//double d_GeV = d*5.0684; //GeV^-1
	
	double q = Deltap/2.;
	double r_min = 0.050;   //fm
	double r_max = 19.95;   //fm
	double Dr    = 0.100;   //fm
	double p_min = 0.00125; //GeV
	double p_max = 0.49875; //GeV
	double Dp    = 0.00250; //GeV
	
	//double bin_r = round( (Deltar-r_min)/Dr );
	//double bin_p = round( (Deltap-p_min)/Dp );
	
	if(Deltar>r_max || q>p_max){
		return 0;
	}
	else if(Deltar<r_max and q<p_max){
		
		double binn_r = (Deltar+r_min)/Dr;
		double binn_p = (q+p_min)/Dp;
	
		int bin_r = int(binn_r);
		int bin_p = int(binn_p);
		
		if(bin_p==0){
			bin_p=2;
		}
		if(bin_r==0){
			bin_r=2;
		}
		if(bin_p>p_vecsize){
			bin_p=p_vecsize;
		}
		if(bin_r>r_vecsize){
			bin_r=r_vecsize;
		}
			
	  	double PDF_r1_p1 =   fArgonne[(bin_r-1)*p_vecsize+bin_p-1] ;
		double PDF_r1_p2 =   fArgonne[(bin_r-1)*p_vecsize+bin_p] ;
	  	double PDF_r2_p1 =   fArgonne[bin_r*p_vecsize+bin_p-1] ;
		double PDF_r2_p2 =   fArgonne[bin_r*p_vecsize+bin_p] ;
	
		double p1 = p_min + (bin_p-1)*Dp;
		double p2 = p_min + bin_p*Dp;
		double r1 = r_min + (bin_r-1)*Dr;
		double r2 = r_min + bin_r*Dr;
	
		double PDF_r1 =  PDF_r1_p1 + (q-p1)/(p2-p1)*(PDF_r1_p2-PDF_r1_p1);
		double PDF_r2 =  PDF_r2_p1 + (q-p1)/(p2-p1)*(PDF_r2_p2-PDF_r2_p1);
		double PDF =  PDF_r1 + (Deltar-r1)/(r2-r1)*(PDF_r2-PDF_r1);
				
		return PDF;
		
	}
}

double sourcesize_function(double p_px,double p_py,double p_pz,double n_px,double n_py,double n_pz, double p_tt, double p_xx, double p_yy, double p_zz, double n_tt, double n_xx, double n_yy, double n_zz) {
    //xx,yy,zz are in fermi!!!
	
    double mp = 0.938;
    double mn = 0.939;
    double D_m = mn+mp;
	
    double p_pT = sqrt(p_px*p_px+p_py*p_py);
    double n_pT = sqrt(n_px*n_px+n_py*n_py);
    double p_E = sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz + mp*mp);
    double n_E = sqrt(n_px*n_px+n_py*n_py+n_pz*n_pz + mn*mn);
	
    double D_px = p_px+n_px;
    double D_py = p_py+n_py;
    double D_pz = p_pz+n_pz;
    double D_p = sqrt( D_px*D_px + D_py*D_py + D_pz*D_pz );
    double D_pT = p_pT+n_pT;
    //D_pP = sqrt( D_p*D_p - D_pT*D_pT )#pz
    double D_E = sqrt( D_p*D_p + D_m*D_m );
    double D_y = 0.5*log((D_E+D_pz)/(D_E-D_pz));
    
    double beta_x = D_px/D_E;
    double beta_y = D_py/D_E;
    double beta_z = D_pz/D_E;
    double beta = sqrt( beta_x*beta_x + beta_y*beta_y + beta_z*beta_z );
    double gamma = D_E/D_m;
    
    double p_tt_prime = gamma*p_tt - gamma*beta_x*p_xx - gamma*beta_y*p_yy - gamma*beta_z*p_zz;
    double p_xx_prime = -gamma*beta_x*p_tt + (1.0+(gamma-1.0)*beta_x*beta_x/(beta*beta))*p_xx + (gamma-1.0)*(beta_x*beta_y/(beta*beta))*p_yy + (gamma-1.0)*(beta_x*beta_z/(beta*beta))*p_zz;
    double p_yy_prime = -gamma*beta_y*p_tt + (gamma-1.0)*(beta_x*beta_y/(beta*beta))*p_xx + (1.0+(gamma-1.0)*beta_y*beta_y/(beta*beta))*p_yy + (gamma-1.0)*(beta_y*beta_z/(beta*beta))*p_zz;
    double p_zz_prime = -gamma*beta_z*p_tt + (gamma-1.0)*(beta_x*beta_z/(beta*beta))*p_xx + (gamma-1.0)*(beta_y*beta_z/(beta*beta))*p_yy + (1.0+(gamma-1.0)*beta_z*beta_z/(beta*beta))*p_zz;
    
    double n_tt_prime = gamma*n_tt - gamma*beta_x*n_xx - gamma*beta_y*n_yy - gamma*beta_z*n_zz;
    double n_xx_prime = -gamma*beta_x*n_tt + (1.0+(gamma-1.0)*beta_x*beta_x/(beta*beta))*n_xx + (gamma-1.0)*(beta_x*beta_y/(beta*beta))*n_yy + (gamma-1.0)*(beta_x*beta_z/(beta*beta))*n_zz;
    double n_yy_prime = -gamma*beta_y*n_tt + (gamma-1.0)*(beta_x*beta_y/(beta*beta))*n_xx + (1.0+(gamma-1.0)*beta_y*beta_y/(beta*beta))*n_yy + (gamma-1.0)*(beta_y*beta_z/(beta*beta))*n_zz;
    double n_zz_prime = -gamma*beta_z*n_tt + (gamma-1.0)*(beta_x*beta_z/(beta*beta))*n_xx + (gamma-1.0)*(beta_y*beta_z/(beta*beta))*n_yy + (1.0+(gamma-1.0)*beta_z*beta_z/(beta*beta))*n_zz;
    
    double source_prime = sqrt( pow(p_xx_prime-n_xx_prime,2.) + pow(p_yy_prime-n_yy_prime,2.) + pow(p_zz_prime-n_zz_prime,2.));
	
	return source_prime;
}

double coalescence_function(double p_px,double p_py,double p_pz,double n_px,double n_py,double n_pz, double pcoal, double source_size, double sigma, double d, int method) {
    double mp = 0.938;
    double mn = 0.939;
    double D_m = mn+mp;
	
    double p_pT = sqrt(p_px*p_px + p_py*p_py);
    double n_pT = sqrt(n_px*n_px + n_py*n_py);
    double p_E = sqrt(p_px*p_px + p_py*p_py + p_pz*p_pz + mp*mp);
    double n_E = sqrt(n_px*n_px + n_py*n_py + n_pz*n_pz + mn*mn);
	
    double D_px = p_px+n_px;
    double D_py = p_py+n_py;
    double D_pz = p_pz+n_pz;
    double D_p = sqrt( D_px*D_px + D_py*D_py + D_pz*D_pz );
    double D_pT = p_pT+n_pT;
    //D_pP = sqrt( D_p*D_p - D_pT*D_pT )#pz
    double D_E = sqrt( D_p*D_p + D_m*D_m );
    double D_y = 0.5*log((D_E+D_pz)/(D_E-D_pz));
    double beta_x = D_px/D_E;
    double beta_y = D_py/D_E;
    double beta_z = D_pz/D_E;
    double beta = sqrt( beta_x*beta_x + beta_y*beta_y + beta_z*beta_z );
    double gamma = D_E/D_m;
    double p_E_prime = gamma*p_E - gamma*beta_x*p_px - gamma*beta_y*p_py - gamma*beta_z*p_pz;
    double p_px_prime = -gamma*beta_x*p_E + (1.0+(gamma-1.0)*beta_x*beta_x/(beta*beta))*p_px + (gamma-1.0)*(beta_x*beta_y/(beta*beta))*p_py + (gamma-1.0)*(beta_x*beta_z/(beta*beta))*p_pz;
    double p_py_prime = -gamma*beta_y*p_E + (gamma-1.0)*(beta_x*beta_y/(beta*beta))*p_px + (1.0+(gamma-1.0)*beta_y*beta_y/(beta*beta))*p_py + (gamma-1.0)*(beta_y*beta_z/(beta*beta))*p_pz;
    double p_pz_prime = -gamma*beta_z*p_E + (gamma-1.0)*(beta_x*beta_z/(beta*beta))*p_px + (gamma-1.0)*(beta_y*beta_z/(beta*beta))*p_py + (1.0+(gamma-1.0)*beta_z*beta_z/(beta*beta))*p_pz;
    double n_E_prime = gamma*n_E - gamma*beta_x*n_px - gamma*beta_y*n_py - gamma*beta_z*n_pz;
    double n_px_prime = -gamma*beta_x*n_E + (1.0+(gamma-1.0)*beta_x*beta_x/(beta*beta))*n_px + (gamma-1.0)*(beta_x*beta_y/(beta*beta))*n_py + (gamma-1.0)*(beta_x*beta_z/(beta*beta))*n_pz;
    double n_py_prime = -gamma*beta_y*n_E + (gamma-1.0)*(beta_x*beta_y/(beta*beta))*n_px + (1.0+(gamma-1.0)*beta_y*beta_y/(beta*beta))*n_py + (gamma-1.0)*(beta_y*beta_z/(beta*beta))*n_pz;
    double n_pz_prime = -gamma*beta_z*n_E + (gamma-1.0)*(beta_x*beta_z/(beta*beta))*n_px + (gamma-1.0)*(beta_y*beta_z/(beta*beta))*n_py + (1.0+(gamma-1.0)*beta_z*beta_z/(beta*beta))*n_pz;
    double D_px_prime = p_px_prime+n_px_prime;
    double D_py_prime = p_py_prime+n_py_prime;
    double D_pz_prime = p_pz_prime+n_pz_prime;
    double D_p_prime = sqrt( D_px_prime*D_px_prime + D_py_prime*D_py_prime + D_pz_prime*D_pz_prime );
    double D_E_prime = sqrt( D_p_prime*D_p_prime + D_m*D_m );
    double D_deltap = sqrt( pow(p_px_prime-n_px_prime,2.) + pow(p_py_prime-n_py_prime,2.) + pow(p_pz_prime-n_pz_prime,2.) );
    
	if(method == 1){
	    if(D_deltap<=pcoal){
			return D_E;
		}
	    else{
			return 0.0;
		}
	}
	if(method == 2){
	    if(D_deltap<=pcoal and source_size<=sigma){
			return D_E;
		}
	    else{
			return 0.0;
		}
	}
	if(method == 3){
		double Prob_WF = WF_SG_probability_PDF(source_size,D_deltap,sigma,d);
	    //double rndm = gRandom->Uniform(0.0, 1.0);
		double rndm = (float) rand()/RAND_MAX;
	    if (rndm <= Prob_WF){
    		return D_E;
		}
	    else{
			return 0.0;
		}
	}
	if(method == 31){
		double Prob_WF = WF_SG_probability_pvalue(source_size,D_deltap,sigma,d);
	    //double rndm = gRandom->Uniform(0.0, 1.0);
		double rndm = (float) rand()/RAND_MAX;
	    if (rndm <= Prob_WF){
    		return D_E;
		}
	    else{
			return 0.0;
		}
	}
	if(method == 4){
		double Prob_WF = WF_Argonne_probability(source_size,D_deltap);
	    //double rndm = gRandom->Uniform(0.0, 1.0);
		double rndm = (float) rand()/RAND_MAX;
	    if (rndm <= Prob_WF){
    		return D_E;
		}
	    else{
			return 0.0;
		}
	}
}

int main(){
    Pythia pythia;
    
    //Initialising additional parameters needed
    pythia.settings.addWord("Main:outdir", "./");
    pythia.settings.addWord("Main:mDM", "100");
    pythia.settings.addWord("Main:p_coalescence", "0.196");
    pythia.settings.addWord("Main:d_coalescence", "1.75");
    pythia.settings.addWord("Main:sigma_coalescence", "3");
    pythia.settings.addWord("Main:methodDbar", "1");
	pythia.settings.addWord("Main:nbins", "180");
    // Initialize Les Houches Event File run. List initialization information.
    pythia.readFile("../../Source/spectrum.cmnd");
  
    double pcoalescence = atof(pythia.word("Main:p_coalescence").data()); //Coalescence momentum in GeV
    double d = atof(pythia.word("Main:d_coalescence").data()); //fm
    double sigma = atof(pythia.word("Main:sigma_coalescence").data()); //fm
    double method_Dbar = atof(pythia.word("Main:methodDbar").data()); //ADD PARAMETER IN MADDM.
    double numberofbins = atof(pythia.word("Main:nbins").data()); //ADD PARAMETER IN MADDM.
    
    pythia.readString("Fragmentation:setVertices = on"); //Setting on the vertex info for baryons and mesons
    pythia.readString("PartonVertex:setVertex = on");
    pythia.readString("2112:mayDecay=on"); //Setting on decay on neutrons and antineutrons
    pythia.readString("-2112:mayDecay=on"); //Because in the Galaxy nbar decay happens
  	cout << "Setting vertex information on" << endl;
  	cout << "All Weak decays on" << endl;
	
    cout << endl;
    cout << endl;
	cout << "###########################################################################"<< endl;
	cout << " INPUT PARAMETERS FOR THE DBAR PRODUCTION, NBAR DECAY AND VERTEX INFO" << endl;
	cout << endl;
	
	/*
	Summary of the methods
	method_Dbar = 0 no Dbar, 1 pcoal, 2 pcoal+sigma, 3 Gaussian WF PDF, 31 Gaussian WF pvalue, 4 Argonne WF, 5 Spherical, 10 all methods
	*/
	
    if(method_Dbar==0){//No Dbar
    	  cout << "Not producing Dbar spectrum with MC" << endl;
    }
    if(method_Dbar==1){//Coalescence pcoal
    	  cout << "Producing Dbar spectrum with coalescence method and pcoal= " << pcoalescence << endl;
    }
	if(method_Dbar==2){//Coalescence model pcoal and sigma<3 (sharp cutoff)
  	  	cout << "Producing Dbar spectrum with coalescence method and pcoal= " << pcoalescence <<  " GeV, sigma= " << sigma << " fm" << endl;
	}
	if(method_Dbar==3){//Gaussian Wigner function
  	  	cout << "Producing Dbar spectrum with method Gaussian Wigner function and parameters d= " << d <<  " fm, sigma= " << " fm" << endl;
  	}
	if(method_Dbar==31){//Gaussian Wigner function
  	  	cout << "Producing Dbar spectrum with method Gaussian Wigner (pvalue) function and parameters d= " << d <<  " fm, sigma= " << " fm" << endl;
  	}
	if(method_Dbar==4){//Argonne Wigner function
		writeArgonne();
  	  	cout << "Producing Dbar spectrum with method Argonne Wigner function, no parameters in the model" << endl;
  	}
	if(method_Dbar==5){//Spherical approach
  	  	cout << "Producing Dbar spectrum with spherical approach and pcoal= " << pcoalescence << endl;
  	}
	if(method_Dbar==10){//Coalescence
		writeArgonne();
  	  	cout << "Producing Dbar spectrum with all the methods and parameters pcoal= " << pcoalescence <<  " GeV,  d= " << d <<  " fm, sigma= " << sigma << " fm" << endl;
  	}
	
	cout << "###########################################################################"<< endl;
    cout << endl;
    cout << endl;
  
    // Initialization
    pythia.init();
    // Allow for possibility of a few faulty events.
    int nAbort = 20;
    int iAbort = 0;
  
    // Setting for the histrogram and choice of variables
    double mDM = atof(pythia.word("Main:mDM").data()); // divide by the DM mass to get the x variable
    string outdir   = pythia.word("Main:outdir"); // Where to write the output file
  
    double eminh = -9.;
    double emaxh = 0.;
    double nbins = numberofbins; //Binninb modifed to account properly for the line spectra. The last bin is -0.05/-0.025/0.00+
    double DeltaBin = (emaxh-eminh)/nbins; 

    double nbins_Dbar = 180; //Binninb modifed to account properly for the line spectra. The last bin is -0.05/-0.025/0.00+
    double DeltaBin_Dbar = (emaxh-eminh)/nbins_Dbar; 
  
    int cont_pbar = 0;
	int cont_pbarP = 0;
    int cont_gamma = 0;
    int cont_nue = 0;
    int cont_numu = 0;
    int cont_nutau = 0;
    int cont_pos = 0;
    int cont_Dbar = 0;
	int cont_Dbar_1 = 0;
	int cont_Dbar_2 = 0;
	int cont_Dbar_3 = 0;
	int cont_Dbar_31 = 0;
	int cont_Dbar_4 = 0;
    // Histogram particle spectra
	
    Hist Dbar("antideuteron spectrum", nbins_Dbar, eminh, emaxh, false, true);
    Hist Dbar_1("antideuteron spectrum", nbins_Dbar, eminh, emaxh, false, true);
    Hist Dbar_2("antideuteron spectrum", nbins_Dbar, eminh, emaxh, false, true);
    Hist Dbar_3("antideuteron spectrum", nbins_Dbar, eminh, emaxh, false, true);
	Hist Dbar_31("antideuteron spectrum", nbins_Dbar, eminh, emaxh, false, true);
    Hist Dbar_4("antideuteron spectrum", nbins_Dbar, eminh, emaxh, false, true);

    Hist gamma("gamma spectrum", nbins, eminh, emaxh, false, true);
    Hist electron("e+- spectrum", nbins, eminh, emaxh, false, true);
    Hist antiproton("pbar spectrum", nbins, eminh, emaxh, false, true);
    Hist antiprotonP("pbar primordial spectrum", nbins, eminh, emaxh, false, true);
    Hist nue("nu_e spectrum", nbins, eminh, emaxh, false, true);
    Hist numu("nu_mu spectrum", nbins, eminh, emaxh, false, true);
    Hist nutau("nu_tau spectrum", nbins, eminh, emaxh, false, true);
    Hist rest("remaining particle spectrum", nbins, eminh, emaxh, false, true);
	
   // Begin event loop.
    int nEvent = 0;
    int maxevent = pythia.settings.mode("Main:NumberOfEvents");
  
    cout << "mDM: " << mDM << " outdir: "<<  outdir << "Nevents " << maxevent <<endl;
  
    for (int iEvent = 0; (iEvent < maxevent || maxevent <= 0) ; ++iEvent) {
      // Generate events. Quit if many failures.
      if (!pythia.next()) {
        if (++iAbort < nAbort) continue;
        cout << "Event generation aborted prematurely, owing to error!\n";
        break;
      }
      // Loop over all particles and select wished final state particles for histrograms.
  	int deuteroncheck=0; //Check ti avoid two deuterons produced with the same particles.
	int deuteroncheck_1=0;
	int deuteroncheck_2=0;
	int deuteroncheck_3=0;
	int deuteroncheck_31=0;
	int deuteroncheck_4=0;
		
  	//pythia.event.list();
	
      for (int i = 0; i < pythia.event.size(); ++i) 
        if (pythia.event[i].isFinal()) {
  		  int idAbs  = pythia.event[i].idAbs();
  		  int idap  = pythia.event[i].id();
  		  double eI  = log10((pythia.event[i].e()-pythia.event[i].m())/mDM);
		  double checkpbar_mother = 0.;
		  double checknbar_mother = 0.;
  		  //select photons
  		  if (idAbs == 22){
  			  gamma.fill(eI); 
  			  cont_gamma++;
  		  }
  		  //select positrons (electrons equivalent)
  		  else if (idap == -11){
  			  electron.fill(eI);
  			  cont_pos++;
  		  }
  		  //select antiprotons
		  else if (idap == -2212){
  			  antiproton.fill(eI);
  			  cont_pbar++;
			  
			  if (pythia.event[pythia.event[i].mother1()].tau()/1e-12<1e5){
	  			  antiprotonP.fill(eI);
	  			  cont_pbarP++;
				  /* Here we remove the antiprotons produced from the following particles: nbar0(-2112), Lambdabar0(-3122)
				  Sigmabar(-3222), Lambda_bbar0(-5122), B_sbar0(-531), B-(-521), B+(521), B0(511), B_s0(531), eta_c(441),
				  */
			  }
			  if (pythia.event[pythia.event[i].mother1()].tau()/1e-12>1e5){
				  checkpbar_mother = 1.;
				  /* Here we remove the antineutrons produced from the following particles: Lambdabar0(-3122), Sigmabar+(-3112), 
				  Sigmabar-(-3222), Lambda_b0(5122), B_sbar0(-531), B_s0(531), B-(-521), B0(511), Bbar0(-511), eta_c(441),
				  */
			  }
  		  }
		  
  		  // select various neutrinos
  		  else if (idap == 12){
  			  nue.fill(eI);
  			  cont_nue++;
  		  }
  		  else if (idap == 14){
  			  numu.fill(eI);
  			  cont_numu++;
  		  }
  		  else if (idap == 16){
  			  nutau.fill(eI);
  			  cont_nutau++;
  		  }
  		  else {
  			  rest.fill(eI);
  		  }
  		  //THE PART BELOW IS RELATED TO THE ANTIDEUTERON SPECTRUM
		  //if (idap == -2212 and deuteroncheck==0 and method_Dbar>0 and checkpbar_mother==0){
		  if (idap == -2212 and deuteroncheck==0 and method_Dbar>0){
  			  for (int j = 0; j < pythia.event.size(); ++j){ //loop over the particles of the list (Not only Final ones)
  				  int jdap  = pythia.event[j].id();
				  
				  //if the second particle is an antineutrons and it is different from the mother 
				  //of the antiprotons (i.e. the antiproton (index i) does not come from the decay 
				  //of the antineutron (index j).)
  				  if (jdap == -2112 and pythia.event[i].mother1()!=j){ 
					  
					  if (pythia.event[pythia.event[j].mother1()].tau()/1e-12>1e5){
						  checknbar_mother = 1.;
						  /* Here we remove the antineutrons produced from the weak decays.
						  Minly from following particles: nbar0(-2112), Lambdabar0(-3122)
						  Sigmabar(-3222), Lambda_bbar0(-5122), B_sbar0(-531), B-(-521), B+(521), B0(511), B_s0(531), eta_c(441),
						  */
					  }
						  
  					  double p_px = pythia.event[i].px();
  					  double p_py = pythia.event[i].py();
  					  double p_pT = sqrt( p_px*p_px + p_py*p_py);
  					  double p_pz = pythia.event[i].pz();
					  
  					  double p_tt = pythia.event[i].tProd()/1e-12; //fermi
  					  double p_xx = pythia.event[i].xProd()/1e-12; //fermi
  					  double p_yy = pythia.event[i].yProd()/1e-12; //fermi
  					  double p_zz = pythia.event[i].zProd()/1e-12; //fermi

  					  double n_px = pythia.event[j].px();
  					  double n_py = pythia.event[j].py();
  					  double n_pT = sqrt( n_px*n_px + n_py*n_py);
  					  double n_pz = pythia.event[j].pz();
				
  					  double n_tt = pythia.event[j].tProd()/1e-12; //fermi
  					  double n_xx = pythia.event[j].xProd()/1e-12; //fermi
  				      double n_yy = pythia.event[j].yProd()/1e-12; //fermi
  					  double n_zz = pythia.event[j].zProd()/1e-12; //fermi
				
  					  double source_size = 0.;
  					  if(method_Dbar==2 || method_Dbar==3 || method_Dbar==4 || method_Dbar==10){
  					    	source_size = sourcesize_function(p_px,p_py,p_pz,n_px,n_py,n_pz,p_tt,p_xx,p_yy,p_zz,n_tt,n_xx,n_yy,n_zz);
  					  }
					  
					  if(method_Dbar!=10 and method_Dbar!=5){
  					  	double D_E = coalescence_function(p_px,p_py,p_pz,n_px,n_py,n_pz,pcoalescence,source_size,d,sigma,method_Dbar);
						
          			    if(D_E>0.0){
           				 	//pythia.event.list();
							deuteroncheck = 1;
                        	double eI  = log10((D_E-D_m)/mDM);
                        	Dbar.fill(eI);
                        	cont_Dbar++;
                    	}
                	  }
					  
					  if(method_Dbar==10){
  					  	double D_E_1 = coalescence_function(p_px,p_py,p_pz,n_px,n_py,n_pz,pcoalescence,source_size,sigma,d,1);
						double D_E_2 = coalescence_function(p_px,p_py,p_pz,n_px,n_py,n_pz,pcoalescence+0.013,source_size,sigma,d,2);
						double D_E_3 = coalescence_function(p_px,p_py,p_pz,n_px,n_py,n_pz,pcoalescence,source_size,sigma,d,3);
						double D_E_31= coalescence_function(p_px,p_py,p_pz,n_px,n_py,n_pz,pcoalescence,source_size,sigma,d-0.87,31);
						double D_E_4 = coalescence_function(p_px,p_py,p_pz,n_px,n_py,n_pz,pcoalescence,source_size,sigma,d,4);
						//double D_E_4 = 0.;
						
          			    if(D_E_1>0.0 and deuteroncheck_1==0 and checkpbar_mother==0 and checknbar_mother==0){
           				 	deuteroncheck_1 = 1;
                        	double eI  = log10((D_E_1-D_m)/mDM);
                        	Dbar_1.fill(eI);
                        	cont_Dbar_1++;
                    	}
          			    if(D_E_2>0.0 and deuteroncheck_2==0){
           				 	deuteroncheck_2 = 1;
                        	double eI  = log10((D_E_2-D_m)/mDM);
                        	Dbar_2.fill(eI);
                        	cont_Dbar_2++;
                    	}
          			    if(D_E_3>0.0 and deuteroncheck_3==0){
							//pythia.event.list();
							source_size = sourcesize_function(p_px,p_py,p_pz,n_px,n_py,n_pz,p_tt,p_xx,p_yy,p_zz,n_tt,n_xx,n_yy,n_zz);
							cout << source_size << endl;
           				 	deuteroncheck_3 = 1;
                        	double eI  = log10((D_E_3-D_m)/mDM);
                        	Dbar_3.fill(eI);
                        	cont_Dbar_3++;
                    	}
          			    if(D_E_31>0.0 and deuteroncheck_31==0){
  						  	//pythia.event.list();
							deuteroncheck_31 = 1;
                        	double eI  = log10((D_E_31-D_m)/mDM);
                        	Dbar_31.fill(eI);
                        	cont_Dbar_31++;
							//cout << cont_Dbar_4 << "  " << cont_Dbar_3 << "  " << cont_Dbar_2 << "  " << cont_Dbar_1 << endl;
                    	}
          			    if(D_E_4>0.0 and deuteroncheck_4==0){
           				 	deuteroncheck_4 = 1;
                        	double eI  = log10((D_E_4-D_m)/mDM);
                        	Dbar_4.fill(eI);
                        	cont_Dbar_4++;
							//cout << cont_Dbar_4 << "  " << cont_Dbar_3 << "  " << cont_Dbar_2 << "  " << cont_Dbar_1 << endl;
                    	}
                	  }
  				}
  			} 
  		}
	
      }
    // End of event loop
      nEvent = iEvent;
    }
    //Statistic and histrograms
    pythia.stat();
    
    if(method_Dbar!=10 and method_Dbar!=0){
    	Dbar.operator*=(1./nEvent/DeltaBin_Dbar);
	}
	else if(method_Dbar==10){
		Dbar_1.operator*=(1./nEvent/DeltaBin_Dbar);
		Dbar_2.operator*=(1./nEvent/DeltaBin_Dbar);
		Dbar_3.operator*=(1./nEvent/DeltaBin_Dbar);
		Dbar_31.operator*=(1./nEvent/DeltaBin_Dbar);
		Dbar_4.operator*=(1./nEvent/DeltaBin_Dbar);
	}
    gamma.operator*=(1./nEvent/DeltaBin);
    electron.operator*=(1./nEvent/DeltaBin);
    antiproton.operator*=(1./nEvent/DeltaBin);
	antiprotonP.operator*=(1./nEvent/DeltaBin);
    nue.operator*=(1./nEvent/DeltaBin);
    numu.operator*=(1./nEvent/DeltaBin);
    nutau.operator*=(1./nEvent/DeltaBin);
  
    if(method_Dbar!=10 and method_Dbar!=0){
    	Dbar.table(outdir + "./antideuterons_spectrum_pythia8.dat", false, true, true);
	}
	else if(method_Dbar==10){
		Dbar_1.table(outdir + "./antideuterons_pcoal_spectrum_pythia8.dat", false, true, true);
		Dbar_2.table(outdir + "./antideuterons_pcoalsigma_spectrum_pythia8.dat", false, true, true);
		Dbar_3.table(outdir + "./antideuterons_GWF_spectrum_pythia8.dat", false, true, true);
		Dbar_31.table(outdir + "./antideuterons_GWF_pvalue_spectrum_pythia8.dat", false, true, true);
		Dbar_4.table(outdir + "./antideuterons_AWF_spectrum_pythia8.dat", false, true, true);
	}
    
    gamma.table(outdir + "./gammas_spectrum_pythia8.dat", false, true, true);
    electron.table(outdir + "./positrons_spectrum_pythia8.dat", false, true, true);
    antiproton.table(outdir + "./antiprotons_spectrum_pythia8.dat", false, true, true);
	antiprotonP.table(outdir + "./antiprotonsPrimordial_spectrum_pythia8.dat");
    nue.table(outdir + "./neutrinos_e_spectrum_pythia8.dat", false, true, true);
    numu.table(outdir + "./neutrinos_mu_spectrum_pythia8.dat", false, true, true);
    nutau.table(outdir + "./neutrinos_tau_spectrum_pythia8.dat", false, true, true);
    rest.table(outdir + "./restx_spectrum_pythia8.dat", false, true, true);
    
    if(method_Dbar!=10 and method_Dbar!=0){
    	cout << gamma << electron << antiproton << antiprotonP << nue << numu << nutau << Dbar << rest;
	}
	else if(method_Dbar==10){
		cout << gamma << electron << antiproton << antiprotonP << nue << numu << nutau << Dbar_1  << Dbar_2  << Dbar_3 << Dbar_31 << Dbar_4 << rest;
	}
	else if(method_Dbar==0){
		cout << gamma << electron << antiproton << nue << numu << nutau << rest;
	}
    
	/*
	//std::string outputFilePath = filename+"antideuterons_spherical_spectrum_pythia8.dat";
	std::string outputFilePath = outdir+"multiplicity_particles_pythia8.dat";
	std::ofstream outputFile(outputFilePath);
	outputFile << "#This files contains the particle multiplicities   " << std::endl;
	outputFile << "#e+  " << std::setprecision(5) << std::scientific << (double)cont_pos/maxevent << "\t" << std::endl;
	outputFile << "#gammas  " << std::setprecision(5) << std::scientific << (double)cont_pos/maxevent << "\t" << std::endl;
	outputFile << "#antiprotons  " << std::setprecision(5) << std::scientific << (double)cont_pbar/maxevent << "\t" << std::endl;
	outputFile << "#nue  " << std::setprecision(5) << std::scientific << (double)cont_nue/maxevent << "\t" << std::endl;
	outputFile << "#numu  " << std::setprecision(5) << std::scientific << (double)cont_numu/maxevent << "\t" << std::endl;
	outputFile << "#nutau  " << std::setprecision(5) << std::scientific << (double)cont_untau/maxevent << "\t" << std::endl;
	*/
  
    cout << endl;
    cout << endl;
	cout << "###########################################################################"<< endl;
    cout << "Multiplicity" << endl;
	
    if(method_Dbar!=10 and method_Dbar!=0){
    	cout << "positrons= " << (double)cont_pos/maxevent << "  gammas= " << (double)cont_gamma/maxevent << "  antiprotons= " << (double)cont_pbar/maxevent << "  antiprotons Primordial= " << (double)cont_pbarP/maxevent << "  Dbar= " << (double)cont_Dbar/maxevent << endl;
		//outputFile << "#Dbar  " << std::setprecision(5) << std::scientific << (double)cont_Dbar_1/maxevent << "\t" << std::endl;
	}
	else if(method_Dbar==10){
		cout << "positrons= " << (double)cont_pos/maxevent << "  gammas= " << (double)cont_gamma/maxevent << "  antiprotons= " << (double)cont_pbar/maxevent << "  antiprotons Primordial= " << (double)cont_pbarP/maxevent << "  Dbar pcoal= " << (double)cont_Dbar_1/maxevent << "  Dbar pcoalsigma= " << (double)cont_Dbar_2/maxevent << "  Dbar GWF= " << (double)cont_Dbar_3/maxevent << "  Dbar GWF(pvalue)= " << (double)cont_Dbar_31/maxevent << "  Dbar AWF= " << (double)cont_Dbar_4/maxevent << endl;
		//outputFile << "#Dbar pcoal  " << std::setprecision(5) << std::scientific << (double)cont_Dbar_1/maxevent << "\t" << std::endl;
		//outputFile << "#Dbar pcoal-sigma  " << std::setprecision(5) << std::scientific << (double)cont_Dbar_2/maxevent << "\t" << std::endl;
		//outputFile << "#Dbar GWF  " << std::setprecision(5) << std::scientific << (double)cont_Dbar_3/maxevent << "\t" << std::endl;
		//outputFile << "#Dbar GWF pvalue  " << std::setprecision(5) << std::scientific << (double)cont_Dbar_31/maxevent << "\t" << std::endl;
		//outputFile << "#Dbar AWF  " << std::setprecision(5) << std::scientific << (double)cont_Dbar_4/maxevent << "\t" << std::endl;
	}
	else if(method_Dbar==0){
		cout << "positrons= " << (double)cont_pos/maxevent << "  gammas= " << (double)cont_gamma/maxevent << "  antiprotons= " << (double)cont_pbar/maxevent << endl;
	}
    
    cout << "nue= " << (double)cont_nue/maxevent << "  numu= " << (double)cont_numu/maxevent << "  nutau= " << (double)cont_nutau/maxevent << endl;
	cout << "###########################################################################"<< endl;
    cout << endl;
    cout << endl;
	
	//outputFile.close();
	
	if(method_Dbar==10 || method_Dbar==5){
		writeDbarspherical(pcoalescence-0.041,mDM,outdir,nbins,2,1);
		writeDbarspherical(pcoalescence-0.041,mDM,outdir,nbins,3,2);
		writeDbarspherical(pcoalescence-0.041,mDM,outdir,nbins,4,2);
	}
    
    //Done
    return 0;
  }
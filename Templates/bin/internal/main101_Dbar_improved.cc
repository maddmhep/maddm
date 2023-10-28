#include "Pythia8/Pythia.h"
#include "math.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;
using namespace Pythia8;
using std::ifstream;


double WF_probability(double sigma,double d,double q) {
	double pigreek = 3.141592;
	double Prob = 8.*pow(d*d/(d*d+sigma*sigma*4),3./2.)*exp(-q*q*d*d*5.0684*5.0684/4.);
	return Prob;
}

double WF_SG_probability(double Deltap, double r, double sigma, double d) {
	double pigreek = 3.141592;
	double d_GeV = d*5.0684; //GeV^-1
    double funcr = (1.-erf(r/(sqrt(2.)*sigma)));
    double funcq = (1.-erf(Deltap*d_GeV/2.));
	double Prob = funcr*funcq;
	return Prob;
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

double coalescence_function(double p_px,double p_py,double p_pz,double n_px,double n_py,double n_pz, double pcoal, double source_size, double d, double sigma, int method) {
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
	    if(D_deltap<pcoal){
			return D_E;
		}
	    else{
			return 0.0;
		}
	}
	
	if(method == 2){
		double Prob_WF = WF_SG_probability(D_deltap,source_size,sigma,d);
	    //double rndm = gRandom->Uniform(0.0, 1.0);
		double rndm = (float) rand()/RAND_MAX;
		//cout << "random number. " << rndm << "  " << Prob_WF/1e-3 << "  " << D_deltap << "  " << source_size << "  " << sigma << "  " << d_D << endl;
	    if (rndm <= Prob_WF){
			//cout << "random number. " << rndm << "  " << Prob_WF/1e-3 << "  " << D_deltap << "  " << source_size << "  " << sigma << "  " << d << endl;
    		return D_E;
		}
	    else{
			return 0.0;
		}
	}
}

int main(){
  // You can always read an plain LHE file,
  // but if you ran "./configure --with-gzip" before "make"
  // then you can also read a gzipped LHE file.
  //#ifdef GZIPSUPPORT
  //bool useGzip = true;
  //#else
  //bool useGzip = false;
  //#endif
  //cout << " useGzip = " << useGzip << endl;
  // Generator. We here stick with default values, but changes
  // could be inserted with readString or readFile.
  Pythia pythia;
  // Initialize Les Houches Event File run. List initialization information.
  //pythia.readString("Beams:frameType = 4");
  //if (useGzip) pythia.readString("Beams:LHEF = ww_events.lhe.gz");
  //else         pythia.readString("Beams:LHEF = ww_events.lhe");
  //pythia.readFile("../../Source/spectrum.cmnd");
  pythia.readFile("spectrum.cmnd");
  
  //double method_Dbar = pythia.parm("Main:method_Dbar"); // 
  double method_Dbar = 2; // 
  pythia.readString("Fragmentation:setVertices = on");
  pythia.readString("PartonVertex:setVertex = on");
  
  //pythia.settings.mode("Next:numberCount",1000);
  //pythia.settings.mode("Next:numberShowInfo",1);
  //pythia.settings.mode("Next:numberShowProcess",1);
  //pythia.settings.mode("Next:numberShowEvent",1);
  
  //pythia.settings.flag("PartonLevel:ISR",false);
  //pythia.settings.flag("PartonLevel:MPI",false);
  //pythia.settings.flag("PDF:lepton",false);
  //pythia.settings.flag("TimeShower:weakShower",true);
  //pythia.settings.mode("TimeShower:weakShowerMode",0);
  //pythia.settings.mode("TimeShower:pTminWeak",0.1);
  // Initialization
  pythia.init();
  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;
  
  // Setting for the histrogram and choice of variables
  double mDM = pythia.parm("Main:spareParm1"); // divide by the DM mass to get the x variable
  string outdir   = pythia.word("Main:spareWord1"); // Where to write the output file
  
  double eminh = -9.;
  double emaxh = 0.;
  double nbins = 100;
  double nbins_Dbar = 50;
  double DeltaBin = (emaxh-eminh)/nbins; 
  double DeltaBin_Dbar = (emaxh-eminh)/nbins_Dbar; 
  
  double pcoalescence = 0.180; //Coalescence momentum in GeV
  //double d = 2.5; //fm
  double d = 2.0; //fm
  //double d = 1.8; //fm
  //double d = 1.6; //fm
  double sigma = 1.0; //fm
  double mp = 0.938;
  double mn = 0.939;
  double D_m = mn+mp;
  
  int cont_pbar = 0;
  int cont_pos = 0;
  int cont_Dbar = 0;
  // Histogram particle spectra
  Hist Dbar("antideuteron spectrum", nbins_Dbar, eminh, emaxh);
  Hist gamma("gamma spectrum", nbins, eminh, emaxh);
  Hist electron("e+- spectrum", nbins, eminh, emaxh);
  Hist proton("p spectrum", nbins, eminh, emaxh);
  Hist nue("nu_e spectrum", nbins, eminh, emaxh);
  Hist numu("nu_mu spectrum", nbins, eminh, emaxh);
  Hist nutau("nu_tau spectrum", nbins, eminh, emaxh);
  Hist rest("remaining particle spectrum", nbins, eminh, emaxh);
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
		
    for (int i = 0; i < pythia.event.size(); ++i) 
      if (pythia.event[i].isFinal()) {
	int idAbs  = pythia.event[i].idAbs();
	int idap  = pythia.event[i].id();
	double eI  = log10((pythia.event[i].e()-pythia.event[i].m())/mDM);
	
	//select photons
	if (idAbs == 22) gamma.fill(eI); 
	//select electrons (positron equivalent)
	else if (idap == -11){
		electron.fill(eI);
		cont_pos++;
	}
	//select proton
	else if (idap == -2212 or idAbs == 2112){
		proton.fill(eI);
		cont_pbar++;
	}
	// select various neutrinos
	else if (idap == 12) nue.fill(eI);
	else if (idap == 14) numu.fill(eI);
	else if (idap == 16) nutau.fill(eI);
	else {
	  rest.fill(eI);
	  //cout << "Error: stable id = " << pythia.event[i].id() << endl;
	}
	
	//THE PART BELOW IS RELATED TO THE ANTIDEUTERON SPECTRUM
	if (idap == -2212 and deuteroncheck==0){ //if the event is an antiproton
		//cout << "found pbar " << "  " << iEvent << "  "  << i << endl;
		for (int j = 0; j < pythia.event.size(); ++j){ //loop over the final particles
			int jdap  = pythia.event[j].id();
			if (jdap == -2112){ //if the second particle is an antineutrons
				//cout << "found nbar  " << "  " << iEvent << "  "  << j << endl;
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
				
				//cout << p_px << "  " << p_py << "  " << p_pz << "  " << n_px << "  " << n_py << "  " << n_pz << endl;
				
				//coalescence_function(double p_px,double p_py,double p_pz,double n_px,double n_py,double n_pz, double pcoal)
				double source_size = 0.;
				if(method_Dbar==2.){
					double source_size = sourcesize_function(p_px,p_py,p_pz,n_px,n_py,n_pz,p_tt,p_xx,p_yy,p_zz,n_tt,n_xx,n_yy,n_zz);
				}
				
				double D_E = coalescence_function(p_px,p_py,p_pz,n_px,n_py,n_pz,pcoalescence,source_size,d,sigma,method_Dbar);
                
                if(D_E>0.0){
					cout << "Dbar found. " <<  D_E << "  " << cont_Dbar << "  " << cont_pbar << endl;
                    deuteroncheck = 1;
                    double eI  = log10((D_E-D_m)/mDM);
                    Dbar.fill(eI);
                    cont_Dbar++;
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
  
  Dbar.operator*=(1./nEvent/DeltaBin_Dbar);
  gamma.operator*=(1./nEvent/DeltaBin);
  electron.operator*=(1./nEvent/DeltaBin);
  proton.operator*=(1./nEvent/DeltaBin);
  nue.operator*=(1./nEvent/DeltaBin);
  numu.operator*=(1./nEvent/DeltaBin);
  nutau.operator*=(1./nEvent/DeltaBin);
  
  Dbar.table("antideuterons_spectrum_pythia8_GWF_d2p0.dat");
  gamma.table("gammas_spectrum_pythia8.dat");
  electron.table("positrons_spectrum_pythia8.dat");
  proton.table("antiprotons_spectrum_pythia8.dat");
  nue.table("neutrinos_e_spectrum_pythia8.dat");
  numu.table("neutrinos_mu_spectrum_pythia8.dat");
  nutau.table("neutrinos_tau_spectrum_pythia8.dat");
  rest.table("restx_spectrum_pythia8.dat");
  
  cout << gamma << electron << proton << nue << numu << nutau << Dbar << rest;
  cout << cont_pos << "  " << cont_pbar << "  " << cont_Dbar << endl;
  
  //Done
  return 0;
}
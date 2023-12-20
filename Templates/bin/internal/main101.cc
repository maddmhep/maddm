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
    pythia.settings.addWord("Main:p_coalescence", "0.180");
    pythia.settings.addWord("Main:d_coalescence", "2");
    pythia.settings.addWord("Main:sigma_coalescence", "1");
    pythia.settings.addWord("Main:methodDbar", "1");
    // Initialize Les Houches Event File run. List initialization information.
    pythia.readFile("../../Source/spectrum.cmnd");
  
    double pcoalescence = atof(pythia.word("Main:p_coalescence").data()); //Coalescence momentum in GeV
    double d = atof(pythia.word("Main:d_coalescence").data()); //fm
    double sigma = atof(pythia.word("Main:sigma_coalescence").data()); //fm
    double method_Dbar = atof(pythia.word("Main:methodDbar").data()); //ADD PARAMETER IN MADDM.
  
    cout << endl;
    cout << endl;
	cout << "###########################################################################"<< endl;
	cout << " INPUT PARAMETERS FOR THE DBAR PRODUCTION, NBAR DECAY AND VERTEX INFO" << endl;
	cout << endl;
    if(method_Dbar==2){//Gaussian WF
  	  pythia.readString("Fragmentation:setVertices = on"); //Setting on the vertex info for baryons and mesons
  	  pythia.readString("PartonVertex:setVertex = off");
  	  pythia.readString("2112:mayDecay=on"); //Setting on decay on neutrons and antineutrons
  	  pythia.readString("-2112:mayDecay=on");
  	  cout << "Producing Dbar spectrum with method Gaussian Wigner function and parameters d= " << d <<  "  sigma= " << sigma << endl;
  	  cout << "Set vertex on" << endl;
  	  cout << "n/nbar decay on" << endl;
    }
    else{
  	  pythia.readString("Fragmentation:setVertices = off"); //Setting off the vertex info for baryons and mesons
  	  pythia.readString("PartonVertex:setVertex = off");
  	  if(method_Dbar==1){//Coalescence
  		  pythia.readString("2112:mayDecay=off"); //Setting off decay on neutrons and antineutrons
  		  pythia.readString("-2112:mayDecay=off"); //This avoids producing Dbar with anti-nucleons produced after nbar decay.
  	  	  cout << "Producing Dbar spectrum with coalescence method and pcoal= " << pcoalescence << endl;
  		  cout << "Set vertex off" << endl;
  		  cout << "n/nbar decay on" << endl;
  	  }
  	  if(method_Dbar==0){//No Dbar
  		  pythia.readString("2112:mayDecay=on"); //Setting off decay on neutrons and antineutrons to take the correct pbar energetics
  		  pythia.readString("-2112:mayDecay=on");
  	  	  cout << "Not producing Dbar spectrum" << endl;
  		  cout << "Set vertex off" << endl;
  		  cout << "n/nbar decay on" << endl;
  	  }
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
    double nbins = 100;
    double DeltaBin = (emaxh-eminh)/nbins; 
  
    double mp = 0.938;
    double mn = 0.939;
    double D_m = mn+mp;
  
    int cont_pbar = 0;
    int cont_gamma = 0;
    int cont_nue = 0;
    int cont_numu = 0;
    int cont_nutau = 0;
    int cont_pos = 0;
    int cont_Dbar = 0;
    // Histogram particle spectra
    Hist Dbar("antideuteron spectrum", nbins, eminh, emaxh);
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
		
  	//pythia.event.list();
	
      for (int i = 0; i < pythia.event.size(); ++i) 
        if (pythia.event[i].isFinal()) {
  		  int idAbs  = pythia.event[i].idAbs();
  		  int idap  = pythia.event[i].id();
  		  double eI  = log10((pythia.event[i].e()-pythia.event[i].m())/mDM);
	
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
  		  //select antiprotons and antineutrons (because they decay into antiprotons) 
  		  else if (idap == -2212 or idap == -2112){
  			  proton.fill(eI);
  			  cont_pbar++;
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
  		  if (idap == -2212 and deuteroncheck==0 and method_Dbar>0){ //if the event is an antiproton
  			  for (int j = 0; j < pythia.event.size(); ++j){ //loop over the final particles
				  
  				  int jdap  = pythia.event[j].id();
				  
  				  if (jdap == -2112 and pythia.event[i].mother1()!=j){ //if the second particle is an antineutrons and its different from the mother of the antiprotons (i.e. the antiproton does not come from the decay of the antineutron.)
					  
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
  					  if(method_Dbar==2.){
  					    	double source_size = sourcesize_function(p_px,p_py,p_pz,n_px,n_py,n_pz,p_tt,p_xx,p_yy,p_zz,n_tt,n_xx,n_yy,n_zz);
  					  }
  				      source_size = sourcesize_function(p_px,p_py,p_pz,n_px,n_py,n_pz,p_tt,p_xx,p_yy,p_zz,n_tt,n_xx,n_yy,n_zz);
  					  double D_E = coalescence_function(p_px,p_py,p_pz,n_px,n_py,n_pz,pcoalescence,source_size,d,sigma,method_Dbar);
                	
        			  	if(D_E>0.0){
							//cout << "Dbar found. " << iEvent << "  " << D_E << "  " << cont_Dbar << "  " << cont_pbar << "  " << deuteroncheck << "  " << source_size << endl;
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
  
    Dbar.operator*=(1./nEvent/DeltaBin);
    gamma.operator*=(1./nEvent/DeltaBin);
    electron.operator*=(1./nEvent/DeltaBin);
    proton.operator*=(1./nEvent/DeltaBin);
    nue.operator*=(1./nEvent/DeltaBin);
    numu.operator*=(1./nEvent/DeltaBin);
    nutau.operator*=(1./nEvent/DeltaBin);
  
    Dbar.table("antideuterons_spectrum_pythia8.dat", false, true);
    gamma.table("gammas_spectrum_pythia8.dat");
    electron.table("positrons_spectrum_pythia8.dat");
    proton.table("antiprotons_spectrum_pythia8.dat");
    nue.table("neutrinos_e_spectrum_pythia8.dat");
    numu.table("neutrinos_mu_spectrum_pythia8.dat");
    nutau.table("neutrinos_tau_spectrum_pythia8.dat");
    rest.table("restx_spectrum_pythia8.dat");
  
    cout << gamma << electron << proton << nue << numu << nutau << Dbar << rest;
  
    cout << endl;
    cout << endl;
	cout << "###########################################################################"<< endl;
    cout << "Multiplicity" << endl;
    cout << "positrons= " << (double)cont_pos/maxevent << "  gammas= " << (double)cont_gamma/maxevent << "  antiprotons= " << (double)cont_pbar/maxevent << "  Dbar= " << (double)cont_Dbar/maxevent << endl;
    cout << "nue= " << (double)cont_nue/maxevent << "  numu= " << (double)cont_numu/maxevent << "  nutau= " << (double)cont_nutau/maxevent << endl;
	cout << "###########################################################################"<< endl;
    cout << endl;
    cout << endl;
  
    //Done
    return 0;
  }
#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"
// #include "/T2Kflux2016/t2kflux_2016_sk_plus250kA.txt" // opens the txt data file

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
// #include "TGraph.h"
#include "TApplication.h"
// #include "TPave.h"

#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char * argv[] )
{

  //allows Canvas to open in
  TApplication *app = new TApplication("app",0,0);

  //Read file:

  TFile* fluxfile = new TFile("/home/yong/Neutrinos/Prob3++.20121225/T2Kflux2016/t2kflux_2016_sk_plus250kA.root");
  TH1D* numu = (TH1D*)fluxfile->Get("enu_sk_numu");
  TH1D* numu_osc = (TH1D*)numu->Clone("enu_sk_numu_osc");


  std::string line;
  std::ifstream dataFile("/home/yong/Neutrinos/Prob3++.20121225/T2Kflux2016/t2kflux_2016_sk_plus250kA.txt");
  int i;
  int size = 220;

  // if(dataFile.is_open()){ // this line is to get the size of data
  //   i=0;
  //   cout << "Checkpoint 1" << endl;
  //   while(std::getline(dataFile,line)){
  //
  //     if(i>2){
  //       size++; // gets the size of the data
  //     }
  //     i++;
  //   }
  //
  // }
  //

  //initialise containers for the data
  double numu_col[size],numub_col[size],nue_col[size],nueb_col[size];

  if(dataFile.is_open()){ // fill containers with flux values
    i = 0; //tracks the row getline is on
    cout << "Checkpoint 2" << endl;
    while(std::getline(dataFile,line)){

      std::stringstream ss(line);
      ss.ignore(20,'\n'); // ignores the first 20 characters

      if(i>2){ // this if loop is to avoid the first three rows of the txt file
        // cout << "Checkpoint 4" << endl;
        ss >> numu_col[i-3];
        ss >> numub_col[i-3];
        ss >> nue_col[i-3];
        ss >> nueb_col[i-3];

        // Un-comment below to view file in terminal
        // cout << "#" << i-2 << "\t" << numu_col[i-3] << "\t" << numub_col[i-3] << "\t" <<  nue_col[i-3] << "\t" <<  nueb_col[i-3] << endl;

      }
      i++;
    }
  }



   //set parameters
   bool kSquared  = true;   // using sin^2(x) variables?
   double DM2     =  2.4e-3; // delta m squared for mu <->tau
   double Theta23 =  0.5   ;
   double Theta13 =  0.025  ;
   double dm2     =  7.6e-5; // delta m squared for e <-> mu
   double Theta12 =  0.312;
   double delta   =  0 * (3.1415926/180.0); //convert to radians

   std::cout << "Using          " << std::endl
	     << "      DM2      " <<  DM2      << std::endl //delta mass squared
	     << "      Theta23  " <<  Theta23  << std::endl
	     << "      Theta13  " <<  Theta13  << std::endl
	     << "      dm2      " <<  dm2      << std::endl //delta mass squared
	     << "      Theta12  " <<  Theta12  << std::endl
	     << "      dcp      " <<  delta    << std::endl; // delta for CP violation


   //make propagator
   BargerPropagator   * bNu;

   bNu = new BargerPropagator( );
   bNu->UseMassEigenstates( false );

   //initialise histogram parameters;
   Double_t numu_oflux[size],numub_oflux[size],nue_oflux[size],nueb_oflux[size];

   int kNuBar = 1.0; // positive for neutrino, negative for antineutrino
   double BasePath = 295; //km
   double Density = 2.3; //g/cm^3

   double energy; //GeV (peak at 0.6GeV)
   double prob;

   //introduce pre-factor beta:
   double beta=1; // where 0<beta<2

   for(int bin = 1; bin <=size; bin++){

     energy = bin*0.05 - 2.5; // middle energy of each bin in GeV

     //Setting MNS matrix for neutrinos
     bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar );
     bNu->propagateLinear( 1*kNuBar, BasePath, Density ); // potentially not needed

     // numu -> nue
     prob = bNu->GetProb(2, 1);
     numu_oflux[bin-1]= (1./beta)*prob*numu_col[bin-1];

     //nue -> nue
     prob = bNu->GetProb(1, 1);
     nue_oflux[bin-1]= (1./beta)*prob*nue_col[bin-1];

     //Setting MNS matrix for anti-neutrinos
     bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, -1*kNuBar );
     bNu->propagateLinear( -1*kNuBar, BasePath, Density ); // potentially not needed

     //numub -> nueb
     prob = bNu->GetProb(-2, -1);
     numub_oflux[bin-1]= beta*prob*numub_col[bin-1];

     //nueb -> nueb
     prob = bNu->GetProb(-1, -1);
     nueb_oflux[bin-1]= beta*prob*nueb_col[bin-1];

   }


  //PLOTTING
  TCanvas *c1 = new TCanvas("c1","Canvas for Oscillated Flux ",200,10,700,500);

  TH1D* h1 = new TH1D("h", "Oscillated Flux vs Energy", 220, 0.,10.);
  h1->GetXaxis()->SetTitle("Energy/GeV");
  h1->GetYaxis()->SetTitle("Flux/cm^{2}");

  for(int j = 0; j < size; j++){
    h1->Fill(numu_oflux[j]);//fill each bin with oscillated flux values

    //maybe can Fill() with weighting using original flux and probability as weighting?
  }
  h1->Draw();

   cout << "Application running..." << endl;
   app->SetReturnFromRun(true);
   app->Run(); // need this to give options for saving and zoom etc
  //  app->Terminate();

  // TPave *pv = new TPave(0.7,0.9,0.9,0.7); // make it top right corner
  // pv->("Oscillation Parameters");
  // pv->AddEntry("gr","#delta m^{2}: ", "");
  // pv->AddEntry("gr","#delta m^{2}: ", "");
  // pv->Draw();

  // bool end;
  // cout << "End? 1/0";
  // cin>> end;
  // app->Connect("TCanvas", "Closed()", "TApplication", app, "Terminate()");

  cout << endl<<"Done Cowboy!" << endl;

  return 0;
}

// int main(int argc, char * argv[] ){
//
//   //allows Canvas to open in
//   TApplication *app = new TApplication("app",0,0);
//
//   plotEllipse(0.6,1);
//   plotEllipse(0.1,1);
//
//   cout << "Application running..." << endl;
//   app->SetReturnFromRun(true);
//   app->Run(); // need this to give options for saving and zoom etc
//   app->Terminate();
// }

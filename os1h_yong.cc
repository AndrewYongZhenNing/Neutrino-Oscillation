#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"

#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TApplication.h"


#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char * argv[] )
{

  //allows Canvas to open in
  TApplication *app = new TApplication("app",0,0);

   //set parameters
   double energy = 0.6; //GeV (peak at 0.6GeV)
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

   //Load files:
   TFile* fluxfile = new TFile("/home/yong/Neutrinos/Prob3++.20121225/T2Kflux2016/t2kflux_2016_plus250kA.root");

   TH1D* numu = (TH1D*)fluxfile->Get("enu_sk_numu"); // creates the numu histogram
   TH1D* numu_osc = (TH1D*)numu->Clone("enu_sk_numu_osc"); //creates a copy of the numu histogram

   TH1D* nue = (TH1D*)fluxfile->Get("enu_sk_nue"); // creates the numu histogram
   TH1D* nue_osc = (TH1D*)nue->Clone("enu_sk_nue_osc"); //creates a copy of the numu histogram

   TH1D* numub = (TH1D*)fluxfile->Get("enu_sk_numub"); // creates the numu histogram
   TH1D* numub_osc = (TH1D*)numu->Clone("enu_sk_numub_osc"); //creates a copy of the numu histogram

   TH1D* nueb = (TH1D*)fluxfile->Get("enu_sk_nueb"); // creates the numu histogram
   TH1D* nueb_osc = (TH1D*)nue->Clone("enu_sk_nueb_osc"); //creates a copy of the numu histogram


   int kNuBar = 1.0; // positive for neutrino, negative for antineutrino
   double BasePath = 295; //km
   double Density = 2.3; //g/cm^3

   //Setting MNS matrix for neutrinos
   bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar );
   bNu->propagateLinear( 1*kNuBar, BasePath, Density ); // potentially not needed

   double update_bin;
   double prob;

   //introduce pre-factor beta:
   double beta=1; // where 0<beta<2

   // numu -> nue
   prob = bNu->GetProb(2, 1);

   for(int bin = 1; bin <=220; bin++){
     update_bin = (1./beta)*prob*numu->GetBinContent(bin);
     numu_osc->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

   }

   //nue -> nue
   prob = bNu->GetProb(1,1);

   for(int bin = 1; bin <=220; bin++){
     update_bin = (1./beta)*prob*nue->GetBinContent(bin);
     nue_osc->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

   }

   // Setting MNS matrix for anti-neutrinos
   bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, -1*kNuBar );
   bNu->propagateLinear( -1*kNuBar, BasePath, Density ); // potentially not needed

   //numub -> nueb
   prob = bNu->GetProb(-2, -1);

   for(int bin = 1; bin <=220; bin++){
     update_bin = (1./beta)*prob*numub->GetBinContent(bin);
     numub_osc->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

   }

   //nueb -> nueb
   prob = bNu->GetProb(-1, -1);

   for(int bin = 1; bin <=220; bin++){
     update_bin = (1./beta)*prob*nueb->GetBinContent(bin);
     nueb_osc->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

   }



  //  double bin_val1,bin_val2;
  //  for(int i =1; i<=220; i++){
  //    bin_val1 = numu->GetBinContent(i);
  //    bin_val2 = numu_osc->GetBinContent(i);
  //    cout << numu_osc->GetBin(i) << "\t" << bin_val1 << "\t" << bin_val2<< endl;
  //  }


  //PLOTTING
  TCanvas *c1 = new TCanvas("c1","Canvas for Oscillated Flux ",200,10,900,600);

  //create an overlay
  THStack *hs = new THStack("hs","");
  // hs->GetXaxis()->SetTitle("Energy/GeV");
  // hs->GetYaxis()->SetTitle("Flux/cm^{2}");
  numu_osc->SetFillColor(kRed);
  numub_osc->SetFillColor(kBlue);
  nue_osc->SetFillColor(kGreen);
  nueb_osc->SetFillColor(42);
  hs->Add(numu_osc);
  hs->Add(numub_osc);
  hs->Add(nue_osc);
  hs->Add(nueb_osc);

  c1->cd(1); hs->Draw("H");
  // TH1D* h1 = new TH1D("h", "Oscillated Flux vs Energy", 220, 0.,10.);
  // h1->GetXaxis()->SetTitle("Energy/GeV");
  // h1->GetYaxis()->SetTitle("Flux/cm^{2}");
  //

  // h1->Draw();
  // numu_osc->Draw("Hist"); // "Hist" fills the histogram
  // numub_osc->Draw("Hist");

   cout << "Application running..." << endl;
   app->SetReturnFromRun(true);
   app->Run(); // need this to give options for saving and zoom etc
   app->Terminate();

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

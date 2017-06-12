#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"

#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TLegend.h"


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
   TFile* fluxfile_p = new TFile("/home/yong/Neutrinos/Prob3++.20121225/T2Kflux2016/t2kflux_2016_plus250kA.root");
   TFile* fluxfile_m = new TFile("/home/yong/Neutrinos/Prob3++.20121225/T2Kflux2016/t2kflux_2016_minus250kA.root");

   //plus250kA file:
   TH1D* numu_p = (TH1D*)fluxfile_p->Get("enu_sk_numu"); // creates the numu histogram
   TH1D* numu_osc_p = (TH1D*)numu_p->Clone("enu_sk_numu_osc"); //creates a copy of the numu histogram

   TH1D* nue_p = (TH1D*)fluxfile_p->Get("enu_sk_nue"); // creates the numu histogram
   TH1D* nue_osc_p = (TH1D*)nue_p->Clone("enu_sk_nue_osc"); //creates a copy of the numu histogram

   TH1D* numub_p = (TH1D*)fluxfile_p->Get("enu_sk_numub"); // creates the numu histogram
   TH1D* numub_osc_p = (TH1D*)numub_p->Clone("enu_sk_numub_osc"); //creates a copy of the numu histogram

   TH1D* nueb_p = (TH1D*)fluxfile_p->Get("enu_sk_nueb"); // creates the numu histogram
   TH1D* nueb_osc_p = (TH1D*)nueb_p->Clone("enu_sk_nueb_osc"); //creates a copy of the numu histogram

   //minus250kA files:
   TH1D* numu_m = (TH1D*)fluxfile_m->Get("enu_sk_numu"); // creates the numu histogram
   TH1D* numu_osc_m = (TH1D*)numu_m->Clone("enu_sk_numu_osc"); //creates a copy of the numu histogram

   TH1D* nue_m = (TH1D*)fluxfile_m->Get("enu_sk_nue"); // creates the numu histogram
   TH1D* nue_osc_m = (TH1D*)nue_m->Clone("enu_sk_nue_osc"); //creates a copy of the numu histogram

   TH1D* numub_m = (TH1D*)fluxfile_m->Get("enu_sk_numub"); // creates the numu histogram
   TH1D* numub_osc_m = (TH1D*)numub_m->Clone("enu_sk_numub_osc"); //creates a copy of the numu histogram

   TH1D* nueb_m = (TH1D*)fluxfile_m->Get("enu_sk_nueb"); // creates the numu histogram
   TH1D* nueb_osc_m = (TH1D*)nueb_m->Clone("enu_sk_nueb_osc"); //creates a copy of the numu histogram


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
     update_bin = (1./beta)*prob*numu_p->GetBinContent(bin);
     numu_osc_p->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

     update_bin = (1./beta)*prob*numu_m->GetBinContent(bin);
     numu_osc_m->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

   }

   //nue -> nue
   prob = bNu->GetProb(1,1);

   for(int bin = 1; bin <=220; bin++){
     update_bin = (1./beta)*prob*nue_p->GetBinContent(bin);
     nue_osc_p->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

     update_bin = (1./beta)*prob*nue_m->GetBinContent(bin);
     nue_osc_m->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

   }

   // Setting MNS matrix for anti-neutrinos
   bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, -1*kNuBar );
   bNu->propagateLinear( -1*kNuBar, BasePath, Density ); // potentially not needed

   //numub -> nueb
   prob = bNu->GetProb(-2, -1);

   for(int bin = 1; bin <=220; bin++){
     update_bin = (1./beta)*prob*numub_p->GetBinContent(bin);
     numub_osc_p->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

     update_bin = (1./beta)*prob*numub_m->GetBinContent(bin);
     numub_osc_m->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

   }

   //nueb -> nueb
   prob = bNu->GetProb(-1, -1);

   for(int bin = 1; bin <=220; bin++){
     update_bin = (1./beta)*prob*nueb_p->GetBinContent(bin);
     nueb_osc_p->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

     update_bin = (1./beta)*prob*nueb_m->GetBinContent(bin);
     nueb_osc_m->SetBinContent(bin,update_bin); // update the Clone copy of the histogram

   }


  //PLOTTING
  TCanvas *c1 = new TCanvas("c1","Canvas for Oscillated Flux ",200,10,900,600);


  //create an overlay
  THStack *hs = new THStack("hs","T2K Oscillated Flux minus250kA");
  // hs->GetXaxis()->SetTitle("Energy/GeV");
  // hs->GetYaxis()->SetTitle("Flux/cm^{2}");

  //plus250kA plot:
  // numu_osc_p->SetFillColor(kRed);
  // numub_osc_p->SetFillColor(kBlue);
  // nue_osc_p->SetFillColor(kGreen);
  // nueb_osc_p->SetFillColor(42);
  // hs->Add(numu_osc_p);
  // hs->Add(numub_osc_p);
  // hs->Add(nue_osc_p);
  // hs->Add(nueb_osc_p);

  //minus250kA plot:
  numu_osc_m->SetFillColor(kRed);
  numub_osc_m->SetFillColor(kBlue);
  nue_osc_m->SetFillColor(kGreen);
  nueb_osc_m->SetFillColor(42);
  hs->Add(numub_osc_m); // adding histograms to HStack depends on order; in this case: the blue one )numub_osc_m appears at the bottom
  hs->Add(numu_osc_m);
  hs->Add(nue_osc_m);
  hs->Add(nueb_osc_m);

  c1->cd(1);
  hs->Draw("H");
  hs->GetXaxis()->SetTitle("Energy/GeV");
  hs->GetXaxis()->SetRangeUser(0,1.2); //SetRangeUser sets the range of bins closest to the user defined boundaries: in this case the bin closest to 0GeV and 1.2GeV

  // hs->GetYaxis()->SetTitle("the Y axis");
  c1->Modified();


  TLegend *legend = new TLegend(0.7,0.7,0.9,0.9);
  // legend->SetHeader("Title","C");
  legend->AddEntry(numu_osc_m,"#nu_{#mu} #rightarrow #nu_{e}","f");
  legend->AddEntry(numub_osc_m,"#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","f");
  legend->AddEntry(nue_osc_m,"#nu_{e} #rightarrow #nu_{e}","f");
  legend->AddEntry(nueb_osc_m,"#bar{#nu}_{e} #rightarrow #bar{#nu}_{e}","f");
  legend->Draw();
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

  cout << endl<<"Done Cowboy!" << endl;

  return 0;
}

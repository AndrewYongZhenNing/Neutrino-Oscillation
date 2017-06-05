#include <math.h>
#include <iostream>
#include <fstream>

#include "BargerPropagator.h"

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TApplication.h"

#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>

using namespace std;
int main(int argc, char * argv[] )
// void os1_yong()
{

  //allows Canvas to open in
  TApplication *app = new TApplication("app",0,0);

   //set parameters
   bool kSquared  = true;   // using sin^2(x) variables?
   double DM2     =  2.4e-3;
   double Theta23 =  0.5   ;
   double Theta13 =  0.025  ;
   double dm2     =  7.6e-5;
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

   cout << "SIZE: " << size << endl;
   double energy = 0.6; //GeV
   int kNuBar = 1.0; // positive for neutrino, negative for antineutrino
   double BasePath = 295; //km
   double Density = 2.3; //g/cm^3

   int size = (2*M_PI/0.01) +1;
   Double_t p[size], pbar[size];

   for(delta = 0; delta < 2*M_PI; delta += 0.01){

   }
   bNu->SetMNS( Theta12,  Theta13, Theta23, dm2, DM2, delta , energy, kSquared, kNuBar ); //creates the MNS matrix given a particular set of parameters (eg values of mixing angle, CPviolating phase factor etc)
   bNu->propagateLinear( 1*kNuBar, BasePath, Density );
   //oscilate mu to electron (or the anti version)
   double prob = bNu->GetProb(2, 1);// 1 for electron, 2 for mu, 3 for tau // gets the MNS matrix element for a particular flavour oscillation given parameters set by setMNS

   std::cout << "Probability: " <<  prob << std::endl;

   // just testing with a standard TGraph example below
  //  TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
  //  Int_t n = 20;
  //  for (Int_t i=0;i<n;i++) {
  //    x[i] = i*0.1;
  //    y[i] = 10*sin(x[i]+0.2);
  //  }
  //  TGraph* gr = new TGraph(n,x,y); // TGraph only takes in pointers/array
  //  gr->Draw();
  //  c1->Update();

  cout << endl<<"Done Cowboy!" << endl;

  // app->Run();
  bool end;
  cout << "End? 1/0";
  cin>> end;
  // app->Connect("TCanvas", "Closed()", "TApplication", app, "Terminate()");

  return 0;

}

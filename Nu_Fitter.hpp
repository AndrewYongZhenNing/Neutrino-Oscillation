//
//  Nu_Fitter.cpp
//  Neutrinos
//
//  Created by Siobhan Alden & Andrew Yong on 12/06/2017.
//  Copyright © 2017 Siobhan Alden & Andrew Yong. All rights reserved.
//

#ifndef Nu_Fitter_hpp
#define Nu_Fitter_hpp

#include "BargerPropagator.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TApplication.h"


class Nu_Fitter{

public:

    Nu_Fitter(int kNuBarVar, std::string path, std::string filename1, std::string filename2, std::string filename3, std::string filename4); // overloaded constructor: takes data from two files instead
    ~Nu_Fitter();

    void make_sum(char hist_type, char phenon_mode, bool oscillate);
    // applies probability to each raw histogram and then sums them up.
    // Command: hist_type: 'f' for (fake/file) data, 'p' for prediction; phenon_mode: 'd' for disappearance; 'a' for appearanceoscillate: true for oscillation, false for no oscillation

    void print_kNu(); //tells user what value kNuBar has at the moment

    double getLLH(); // calculates the log likelihood

    void show_hist(char hist_type);
    // show the prediction histogram for analysis
    // Command: hist_type: 'f' for (fake/file) data, 'p' for prediction;

    std::vector<double> return_cparam(); // return currentsPars

    std::vector<std::string> return_sparam(); // return parsName

    void set_param(int index, double val, char vector_type);
    //redefine the value of the index'th parameter
    //Command: vector_type: 'c' to set currentPars, 'p' to set proposedPars

    double sigma_cc(int nu_mode, double bin_E);
    // returns the pre-factor when taking into account interaction cross-section
    //Command: 1: nue; -1: nueb; 2: numu; -2:numub

    void add_param(double value, std::string name); // adds new parameter into Nu_Fitter

    void get_integral(double l_bound, double u_bound); // gets the integral of _Data using lower and upper bound in energy(GeV)

    void get_integral(); // gets the integral of _Data using lower and upper bound in energy(GeV)


protected:

    bool kSquared;   // using sin^2(x) variables?
    int kNuBar; // positive for neutrino, negative for antineutrino
    double BasePath; //km
    double Density;
    int nbin;
    std::string _path;
    std::string _filename1;
    std::string _filename2;
    std::string _filename3;
    std::string _filename4;
    std::vector<double> currentPars;
    std::vector<double> proposedPars;
    std::vector<std::string> parsName;

    TH1D* _Data;
    TH1D* _Prediction;
    TH1D* _input1;// nu mu
    TH1D* _input2;// nu mu bar
    TH1D* _input3;// nu e
    TH1D* _input4;// nu e bar

    //initialises Barger propagator pointers here
    BargerPropagator   * bNu;
    BargerPropagator   * bNu2;




};


#endif /* Nu_Fitter_hpp */

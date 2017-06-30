//
//  Markov.hpp
//
//
//  Created by Siobhan Alden & Andrew Yong on 16/06/2017.
//
//

#ifndef Markov_Chain_hpp
#define Markov_Chain_hpp

#include "Nu_Fitter.hpp"
#include "Disappearance.hpp"
#include "Appearance.h"
// #include "Appearance.cpp"


#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TApplication.h"



class Markov_Chain{

public:

    Markov_Chain(std::vector<double> &obj_pars, std::vector<std::string> &parsName, int nstep, std::string filename);
    ~Markov_Chain();
    void startMH(std::vector<double>& currentPars, Nu_Fitter* oscObject, char phenon_mode);
    void startMH(std::vector<double> &currentPars, Disappearance *oscObject, bool dis);
    void startMH(std::vector<double> &currentPars, Appearance *plusObj, Appearance *minusObj);
    void set_param(int index); // if element is true, the parameter associated with the true element is varied in the Markov chain.
    void set_width(int index, double value);// fine tunes the width of a particular parameter
    void apply_constraint(); // hard-code constraints on proposedPars such that if proposedPars elements is unphysical/unmathematical, set the proposedLLH to be large

private:

    int steps;
    double currentLLH;
    double proposedLLH;
    std::vector<double> proposedPars;
    std::vector<bool> pars_bool;
    std::vector<double> width;
    std::vector<std::string> branchName;

    TFile* file;
    TTree* tree;
    TRandom3* rnd;



};


#endif /* Markov_Chain_hpp */

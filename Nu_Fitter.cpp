//
//  Nu_Fitter.cpp
//  Neutrinos
//
//  Created by Siobhan Alden & Andrew Yong on 12/06/2017.
//  Copyright © 2017 Siobhan Alden & Andrew Yong. All rights reserved.
//

#include "Nu_Fitter.hpp"


Nu_Fitter::Nu_Fitter(int kNuBarVar, std::string path, std::string filename1, std::string filename2, std::string filename3, std::string filename4){

    _path = path; // inputs the file directory
    _filename1 = filename1; // within the .root file, get each of the flux files
    _filename2 = filename2;
    _filename3 = filename3;
    _filename4 = filename4;

    // each _input file takes the raw histogram made from the flux files
    TFile* nu_file = new TFile(_path.c_str());
    _input1 = (TH1D*)nu_file->Get(_filename1.c_str());
    _input2 = (TH1D*)nu_file->Get(_filename2.c_str());
    _input3 = (TH1D*)nu_file->Get(_filename3.c_str());
    _input4 = (TH1D*)nu_file->Get(_filename4.c_str());

    _Data  = (TH1D*)_input1->Clone("Prob");// overwrite with oscillation
    _Prediction = (TH1D*)_input1->Clone("Prob");// overwrite with markov chain
    _Data->SetTitle("Some make_sum1 title");

    kSquared  = true; // hard coded for now
    BasePath = 295; //km
    Density = 2.3;
    kNuBar = kNuBarVar;// switch between neutrino and anti
    nbin  = _input1->GetNbinsX();

    //Storing oscillation parameters in vector
    currentPars.push_back(0); // energy : 0
    currentPars.push_back(2.4e-3);//DM2 : 1
    currentPars.push_back(0.5); //theta23 : 2
    currentPars.push_back(0.025); //theta13 : 3
    currentPars.push_back(7.6e-5); //dm2 : 4
    currentPars.push_back(0.312); //theta12 : 5
    currentPars.push_back(0); // deltacp : 6
    currentPars.push_back(10); // n : 7
    currentPars.push_back(1); // beta : 8

    proposedPars.push_back(0); // energy : 0
    proposedPars.push_back(2.4e-3);//DM2 : 1
    proposedPars.push_back(0.5); //theta23 : 2
    proposedPars.push_back(0.025); //theta13 : 3
    proposedPars.push_back(7.6e-5); //dm2 : 4
    proposedPars.push_back(0.312); //theta12 : 5
    proposedPars.push_back(0); // deltacp : 6
    proposedPars.push_back(10); // n : 7
    proposedPars.push_back(1); // beta : 8

    //passes the varaible name into parsName
    parsName.push_back("Energy"); // energy : 0
    parsName.push_back("DM2");//DM2 : 1
    parsName.push_back("Theta23"); //theta23 : 2
    parsName.push_back("Theta13"); //theta13 : 3
    parsName.push_back("dm2"); //dm2 : 4
    parsName.push_back("Theta12"); //theta12 : 5
    parsName.push_back("Delta"); // deltacp : 6
    parsName.push_back("n"); // n : 7
    parsName.push_back("Beta"); // beta : 8


    bNu = new BargerPropagator( );
    bNu->UseMassEigenstates( false );
    bNu2 = new BargerPropagator( );
    bNu2->UseMassEigenstates( false );

}


Nu_Fitter::~Nu_Fitter(){}


void Nu_Fitter::make_sum(char hist_type, char phenon_mode, bool oscillate){

    // if statement for probability for oscillation or probability of no oscillation
    int in_nu1, out_nu1, in_nu2, out_nu2;

    if(oscillate){ // considers neutrino oscillation

      //Appearance
        in_nu1 = 2;
        out_nu1 = 1;
      // Disappearance
        in_nu2 = 1;
        out_nu2 = 1;

    }

    else if(!oscillate){ // ensures it oscillates to itself (eg numu->numu)

        in_nu1 = 2;
        out_nu1 = 2;
        in_nu2 = 1;
        out_nu2 = 1;

    }

    // BargerPropagator   * bNu;
    // BargerPropagator   * bNu2;
    //
    // bNu = new BargerPropagator( );
    // bNu->UseMassEigenstates( false );
    // bNu2 = new BargerPropagator( );
    // bNu2->UseMassEigenstates( false );

    double E,osci_prob1,osci_prob2,osci_prob3,osci_prob4,bin_content1,bin_content2,bin_content3,bin_content4,weight;

    int count = 0;
    for ( int i = 1; i<=nbin; i++){

        count++;

        E = _input1->GetXaxis()->GetBinCenter(i);

        if(hist_type=='f'){

          // neutrino PMNS
            bNu->SetMNS( currentPars[5], currentPars[3], currentPars[2], currentPars[4], currentPars[1], currentPars[6] , E, kSquared, kNuBar );

          // antineutrino PMNS
            bNu2->SetMNS( currentPars[5], currentPars[3], currentPars[2], currentPars[4], currentPars[1], currentPars[6] , E, kSquared, -1*kNuBar );
        }

        else if(hist_type=='p'){

            bNu->SetMNS( proposedPars[5], proposedPars[3], proposedPars[2], proposedPars[4], proposedPars[1], proposedPars[6] , E, kSquared, kNuBar );


            bNu2->SetMNS( proposedPars[5], proposedPars[3], proposedPars[2], proposedPars[4], proposedPars[1], proposedPars[6] , E, kSquared, -1*kNuBar );

        }

        bNu->propagateLinear( 1*kNuBar, BasePath, Density);

        bNu2->propagateLinear( -1*kNuBar, BasePath, Density);

        // individual probabilities for each file
        osci_prob1 = bNu->GetProb(in_nu1,out_nu1);//2,1
        osci_prob2 = bNu2->GetProb(in_nu1,out_nu1);//-2,-1
        osci_prob3 = bNu->GetProb(in_nu2,out_nu2);//1,1
        osci_prob4 = bNu2->GetProb(in_nu2,out_nu2);//-1,-1
        bin_content1 = _input1->GetBinContent(i);
        bin_content2 = _input2->GetBinContent(i);
        bin_content3 = _input3->GetBinContent(i);
        bin_content4 = _input4->GetBinContent(i);

        // std::cout << "i: " << i <<"\tp1: " << osci_prob1 << "\tp2: " <<6 osci_prob2 << "\tp3: " << osci_prob3 << "\tp4: " << osci_prob4 << std::endl;



        if(hist_type == 'f'){ // applies changes to the _Data histogram

          if(phenon_mode == 'd'){ //Disappearance case
            weight = (1./currentPars[8])*sigma_cc(2,E)*osci_prob1*bin_content1 + currentPars[8]*sigma_cc(-2,E)*osci_prob2*bin_content2 + (1./currentPars[8]*sigma_cc(2,E))*osci_prob3*bin_content3 + currentPars[8]*sigma_cc(-2,E)*osci_prob4*bin_content4;
          }

          else if(phenon_mode == 'a'){ //Appearance case
            weight = (1./currentPars[8])*sigma_cc(1,E)*osci_prob1*bin_content1 + currentPars[8]*sigma_cc(-1,E)*osci_prob2*bin_content2 + (1./currentPars[8]*sigma_cc(1,E))*osci_prob3*bin_content3 + currentPars[8]*sigma_cc(-1,E)*osci_prob4*bin_content4;
          }

          else{
            std::cout << "make_sum() Error:\nInvalid phenon_mode command." << std::endl;
            weight = 0;
          }

            _Data->SetBinContent(i,weight);

            // std::cout << "i: " << i << "\tE: " << E << "\t_Data: " << _Data->GetBinContent(i) << "\tw: " << weight << std::endl;

        }

        else if(hist_type == 'p'){ // applies changes to the _Prediction histogram

          if(phenon_mode == 'd'){ //Disappearance case
            weight = (1./proposedPars[8])*sigma_cc(2,E)*osci_prob1*bin_content1 + proposedPars[8]*sigma_cc(-2,E)*osci_prob2*bin_content2 + (1./proposedPars[8]*sigma_cc(2,E))*osci_prob3*bin_content3 + proposedPars[8]*sigma_cc(-2,E)*osci_prob4*bin_content4;
          }

          else if(phenon_mode == 'a'){ //Appearance case
            weight = (1./proposedPars[8])*sigma_cc(1,E)*osci_prob1*bin_content1 + proposedPars[8]*sigma_cc(-1,E)*osci_prob2*bin_content2 + (1./proposedPars[8]*sigma_cc(1,E))*osci_prob3*bin_content3 + proposedPars[8]*sigma_cc(-1,E)*osci_prob4*bin_content4;
          }

          else{
            std::cout << "make_sum() Error:\nInvalid phenon_mode command." << std::endl;
            weight = 0;
          }

            _Prediction->SetBinContent(i,weight);

        }

        else{

            std::cout << "Invalid Hist Type" << std::endl;

        }

        // std::cout << "Data: " << _Data->GetBinContent(i) << "\tPrediciton: " << _Prediction->GetBinContent(i) <<  std::endl;

    }



}



void Nu_Fitter::print_kNu(){

    if(kNuBar ==1){
        std::cout << "Neutrino" << std::endl;
    }

    else if(kNuBar == -1){

        std::cout << "Anti-neutrino" << std::endl;
    }

}

double Nu_Fitter::getLLH(){

    double LLH = 0;
    double lambda,N;

    for(int j = 1; j<=nbin; j++){

        lambda = _Prediction->GetBinContent(j);
        N = _Data->GetBinContent(j);


        if(N>0&&lambda>0){ // to prevent nan at later bins(higher energies) with 0 entries

            LLH += lambda-N - N*log(lambda/N);
            //std::cout << " lambda " << lambda << " Data " << N << " total LLH " << LLH << " LLH " << lambda-N - N*log(lambda/N) << std::endl;

        }

    }


    return LLH;

}

void Nu_Fitter::show_hist(char hist_type){

    //allows Canvas to open in
    TApplication *app = new TApplication("app",0,0);

    //PLOTTING
    TCanvas *c1 = new TCanvas("c1","Canvas",2000,1000);

    if(hist_type == 'f'){
      _Data->Draw("HIST"); // to plot histogram without original error bar
    }

    else if(hist_type=='p'){
      _Prediction->Draw("HIST SAME C"); // to plot histogram without original error bar
    }

    else{
      std::cout << "show_hist() Error:\nInvalid command. Please close the canvas and try again." << std::endl;
    }

    c1->Update();
    app->Run(); // need this to give options for saving and zoom etc

}

std::vector<double> Nu_Fitter::return_cparam(){

    return currentPars;

}

std::vector<std::string> Nu_Fitter::return_sparam(){
  return parsName;
}

void Nu_Fitter::set_param(int index, double val, char vector_type){
    if(vector_type=='c'){
        currentPars[index] = val;
    }

    else if(vector_type == 'p'){
        proposedPars[index] = val;
    }

    else{
        std::cout << "Invalid command for vector_type. Please choose either 'c' or 'p'." << std::endl;
    }

}

double Nu_Fitter::sigma_cc(int nu_mode, double bin_E){

  double nucleons = 1.3539*10e31;
  double sig; // cm^2/nucleon

  if(nu_mode==1){ //nue mode
    sig = (0.7*10e-38)*bin_E;
  }

  else if(nu_mode==-1){ //nueb mode
    sig = (0.25*10e-38)*bin_E;
  }

  else if(nu_mode==2){ //numu mode
    sig = (0.74*10e-38)*bin_E;
  }

  else if(nu_mode==-2){
    sig = (0.26*10e-38)*bin_E;
  }

  else{
    std::cout << "sigma_cc() error:\nInvalid oscillation mode. Please choose either 1,-1,2,-2." << std::endl;
    sig=0;
  }

  return sig*nucleons; // returns cross section (cm^2)

}

void Nu_Fitter::add_param(double value, std::string name){

  currentPars.push_back(value);
  parsName.push_back(name);

  int last_element = currentPars.size()-1;
  // int last_element = currentPars.back().c;

  std::cout << "Parameter: " << name << "\tValue: " << value << "\tIndex: " << currentPars[last_element] << std::endl;
}

void Nu_Fitter::get_integral(double E_lower, double E_upper){

  double f_integral,l_bound,u_bound;

  l_bound = E_lower;
  u_bound = E_upper;

  if(l_bound < 10 && u_bound < 10){

    l_bound = l_bound/0.05;
    u_bound = u_bound/0.05;

    if(l_bound == 0){ // set it to be first bin
      l_bound = 1;
    }

    f_integral = _Data->Integral(l_bound,u_bound);

  }

  else if(l_bound >10 && u_bound >10){

    double interval = 10/0.05; // the number of bins before 10GeV
    l_bound = l_bound-10 + interval;
    u_bound = u_bound-10 + interval;
    f_integral = _Data->Integral(l_bound,u_bound);

  }

  else if(l_bound <10 && u_bound>10){

    double interval = 10/0.05;
    l_bound = l_bound/0.05;
    u_bound = u_bound-10 + interval;
    f_integral = _Data->Integral(l_bound,u_bound);

  }

  else if(l_bound < 0 && u_bound < 0){
    std::cout << "get_integral() error:\nInvalid boundaries." << std::endl;
    f_integral=0;
  }


  std::cout <<"Lower energy boundary: " << E_lower << "GeV\tUpper energy boundary: " << E_upper << "GeV\nFile integral: " << f_integral << std::endl;

}

void Nu_Fitter::get_integral(){
  double f_integral = _Data->Integral();


  std::cout <<"Lower energy boundary: 0GeV\tUpper energy boundary: 30GeV\nFile integral: " << f_integral << std::endl;
}

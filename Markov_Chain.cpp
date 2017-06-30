//
//  Markov_Chain.cpp
//
//
//  Created by Siobhan Alden & Andrew Yong on 16/06/2017.
//
//

#include "Markov_Chain.hpp"

Markov_Chain::Markov_Chain(std::vector<double> &obj_pars, std::vector<std::string> &parsName, int nstep, std::string filename){

    steps = nstep;
    file = new TFile(filename.c_str(),"RECREATE");
    tree = new TTree("t","");

    for(int i = 0; i < parsName.size();i++){ // creates the leaflist for each parameter branch
      branchName.push_back(parsName[i]+"/D");
    }

    double width_factor = 100; // default width factor

    for(int j = 0; j<obj_pars.size(); j++){ //fills the branch and initialises width vector

      TString arg1 = parsName[j].c_str();
      TString arg3 = branchName[j].c_str();

      tree->Branch(arg1, &obj_pars[j],arg3);
      width.push_back(obj_pars[j]/width_factor); // can be fine-tuned with set_width(...)

    }

    tree->Branch("LLH",&currentLLH,"LLH/D");

    tree->Fill();

    //sets all parameters to be false by default; only parameters which are set true will vary in MCMC
    for( int k = 0; k<=obj_pars.size(); k++){

        pars_bool.push_back(false);
        proposedPars.push_back(0);

    }

    rnd  = new TRandom3(0);

    //Parameters and the corresponding index:
    //E=0 : 0
    //DM2=2.4e-3 : 1
    //theta23=0.5 : 2
    //theta13=0.025 : 3
    //dm2=7.6e-3 : 4
    //theta12=0.312 : 5
    //deltacp=0 : 6
    //n=10 : 7
    //beta=1 : 8
    //coef1=1.00043 : 9
    //coef2=-0.00513335 : 10
    //coef3=-0.81076 : 11
    //coef4=-0.0217913 : 12
    //coef5=0.223944 : 13
    //coef6=0.0340821 : 14
    //coef7=-0.0569585 : 15
    //coef8=0.0133551 : 16
    //coef9=-0.000970247 : 17




}

Markov_Chain::~Markov_Chain(){};

//Disappearance MCMC:
void Markov_Chain::startMH(std::vector<double> &currentPars, Disappearance* oscObject, bool dis){

    currentLLH = oscObject->getLLH();

    for(int i=0; i<steps; i++){
        if(i%1000==0 || i==0){

            if(i==0){
              std::cout << "Initialising random walk..." << std::endl;
            }
            else{
              std::cout << "i " << i << "Running random walk..." << std::endl;
            }

        }

        for(int k = 0; k<currentPars.size();k++){
            if(pars_bool[k]){

                proposedPars[k]=rnd->Gaus(currentPars[k],width[k]);
                oscObject->set_paras_d(k, proposedPars[k],'p');

            }

        }

        oscObject->taylor('p');
        proposedLLH = oscObject->getLLH();

        apply_constraint();

        double accProb = TMath::Min(1.,TMath::Exp(currentLLH-proposedLLH));
        double fRandom = rnd->Rndm();
        std::cout << i << "\t" << currentLLH << "\t" << proposedLLH << "\t" << currentLLH-proposedLLH << std::endl;
        std::cout << "i " << i << " accprob " << accProb << " fRandom " << fRandom << " current LLH " << currentLLH << " proposed LLH " << proposedLLH <<  " parameter 1 current: " << currentPars[1] << " parameter 1 proposed:" << proposedPars[1] << std::endl;
         if ( fRandom <= accProb )
        {

            for(int m=0; m<currentPars.size(); m++)
            {

                if(pars_bool[m]){
                    //count++;
                  //  std::cout << "coef " << m << " " << currentPars[m] << " proposed " << m << " " << proposedPars[m] << std::endl;
                    currentPars[m]=proposedPars[m];

                  }

            }

            currentLLH=proposedLLH;
        }

        //std::cout << " 0 " << coef[0] << " 1 " << coef[1] << " 2 " << coef[2]  << std::endl;
        tree->Fill();


        if(i==steps-1){
            std::cout << "Metropolis-Hastings algorithm completed." << std::endl;

        }
    }

   // std::cout << count << std::endl;
    tree->Write();
    file->Close();

}

//Appearance MCMC:
void Markov_Chain::startMH(std::vector<double> &currentPars, Appearance *plusObj, Appearance *minusObj){

    std::cout << "Initialising random walk..." << std::endl;
    currentLLH = plusObj->getLLH() + minusObj->getLLH(); // add LLH from two files

    for(int i=0; i<steps; i++){

        if(i%1000==0 || i==0){
            std::cout << "Running random walk..." << std::endl;
        }

        for(int j=0; j<currentPars.size(); j++) {

            if(pars_bool[j]){

                proposedPars[j]=rnd->Gaus(currentPars[j],width[j]);
                plusObj->set_param(j, proposedPars[j],'p');
                minusObj->set_param(j, proposedPars[j],'p');

            }

            //re-make the histograms again with newly updated parameters
            plusObj->make_sum('p','a',true);
            minusObj->make_sum('p','a',true);

        }

        proposedLLH = plusObj->getLLH() + minusObj->getLLH();

        apply_constraint();

        double accProb = TMath::Min(1.,TMath::Exp(currentLLH-proposedLLH));
        double fRandom = rnd->Rndm();
        // std::cout << "i " << i << " accprob " << accProb << " fRandom " << fRandom << " current LLH " << currentLLH << " proposed LLH " << proposedLLH <<  " parameter 1 current: " << currentPars[1] << " parameter 1 proposed:" << proposedPars[1] << std::endl;

        if ( fRandom <= accProb )
        {
            for(int k=0; k<currentPars.size(); k++)
            {
                if(pars_bool[k]){
                    currentPars[k]=proposedPars[k];
                    //std::cout << "i " << i << " proposed " << proposedPars[k] <<std::endl;
                }

            }
            currentLLH=proposedLLH;

        }
        //std::cout << "i " << i << " Theta23 " << proposedPars[2] <<std::endl;

        tree->Fill();

        if(i==steps-1){
            std::cout << "Metropolis-Hastings algorithm completed." << std::endl;
        }

    }

    tree->Write();
    file->Close();


}



void Markov_Chain::set_param(int index){

    pars_bool[index] = true;

}

void Markov_Chain::set_width(int index, double value){

    width[index] = value;

}

void Markov_Chain::apply_constraint(){

  if(proposedPars[1]<0){ // if DM2<0
    proposedLLH+=1.0E6;
  }

  if(proposedPars[2]<0 || proposedPars[2]>1){ // theta13>1 or theta13<0
    proposedLLH+=1.0E6;
  }

  if(proposedPars[3]<0 || proposedPars[3]>1){ // theta13>1 or theta13<0
    proposedLLH+=1.0E6;
  }

  if(proposedPars[4]<0){ // if dm2<2
    proposedLLH+=1.0E6;
  }

  if(proposedPars[5]<0 || proposedPars[5]>1){ // theta13>1 or theta13<0
    proposedLLH+=1.0E6;
  }

  if(proposedPars[6]>TMath::Pi()){ // deltaCP>Pi
    proposedPars[6]-=2*TMath::Pi();
  }

  else if(proposedPars[6]<-1*TMath::Pi()){//deltaCP<-Pi
    proposedPars[6]+=2*TMath::Pi();
  }

  if(proposedPars[8] <0){ //beta < 0
    proposedLLH+=1.0E6;
  }



}


//Nu_Fitter MCMC:
void Markov_Chain::startMH(std::vector<double>& currentPars, Nu_Fitter* oscObject, char phenon_mode){

    currentLLH = oscObject->getLLH();

    for(int i=0; i<steps; i++){

        if(i%1000==0 || i==0){

            if(i==0){
                std::cout << "Initialising random walk..." << std::endl;
            }

            else{
                std::cout << "Running random walk..." << std::endl;
            }

        }

        for(int j=0; j<currentPars.size(); j++) {

            if(pars_bool[j]){

                proposedPars[j]=rnd->Gaus(currentPars[j],width[j]);
                oscObject->set_param(j, proposedPars[j],'p');

            }

            oscObject->make_sum('p',phenon_mode,true);

        }

        proposedLLH = oscObject->getLLH();

        double accProb = TMath::Min(1.,TMath::Exp(currentLLH-proposedLLH));
        double fRandom = rnd->Rndm();
        // std::cout << "i " << i << " accprob " << accProb << " fRandom " << fRandom << " current LLH " << currentLLH << " proposed LLH " << proposedLLH <<  " parameter 2 current: " << currentPars[2] << " parameter 2 proposed:" << proposedPars[2] << std::endl;

        if ( fRandom <= accProb )
        {
            for(int k=0; k<currentPars.size(); k++)
            {
                if(pars_bool[k]){
                    currentPars[k]=proposedPars[k];
                    //std::cout << "i " << i << " proposed " << proposedPars[k] <<std::endl;
                }

            }
            currentLLH=proposedLLH;

        }
        //std::cout << "i " << i << " Theta23 " << proposedPars[2] <<std::endl;

        tree->Fill();

        if(i==steps-1){
            std::cout << "Metropolis-Hastings algorithm completed." << std::endl;
        }

    }

    tree->Write();
    file->Close();

}

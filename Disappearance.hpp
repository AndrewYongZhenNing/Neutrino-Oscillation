//
//  Disappearance.hpp
//
//
//  Created by Siobhan Alden on 13/06/2017.
//
//

#ifndef Disappearance_hpp
#define Disappearance_hpp

#include "Nu_Fitter.hpp"
//#include "Nu_Fitter.cpp"

class Disappearance : public Nu_Fitter {

public:

    Disappearance(int kNuBarVar, std::string path, std::string filename, std::string filename2, std::string filename3, std::string filename4, int mass);
    ~Disappearance();
    void make_Prediction(char hist_type, int which);
    double series( double E);
    double fact(int n);
    double approx( double E);
    double osci_prob( double E);
    double taylorinv(double x, std::vector<double>& par);
    void taylor(char hist_type);
    //void fitinvE();
    std::vector<double> return_coef_pars();
    void set_paras_d(int index, double val, char vector_type);

private:

    double L = 295;
    double dm_sq;
    double a;
    std::vector<double> Ene, invE, prob;
    std::vector<double> coefs;
    std::vector<double> bin1;
    std::vector<double> bin2;
    std::vector<double> bin3;
    std::vector<double> bin4;



};

#endif /* Disappearance_hpp */

//
//  fubar.cpp
//
//
//  Created by Siobhan Alden & Andrew Yong on 12/06/2017.
//
//

#include <stdio.h>
#include <iostream>


#include "Nu_Fitter.hpp"
#include "Appearance.h"
#include "Markov_Chain.hpp"
#include "Disappearance.hpp"


int main(){

    // Andrew's Test

  //  Appearance plus = Appearance(1, "/home/yong/Neutrinos/Prob3++.20121225/T2Kflux2016/t2kflux_2016_plus250kA.root", "enu_sk_numu","enu_sk_numub","enu_sk_nue","enu_sk_nueb");
  //  Appearance minus = Appearance(1, "/home/yong/Neutrinos/Prob3++.20121225/T2Kflux2016/t2kflux_2016_minus250kA.root", "enu_sk_numu","enu_sk_numub","enu_sk_nue","enu_sk_nueb");
  //
  //
  //  plus.make_sum('f','a', true);// Data
  //  plus.make_sum('p','a',true);
  //  minus.make_sum('f','a', true);// Data
  //  minus.make_sum('p','a',true);
  //
  // //  plus.get_integral();
  // //  plus.get_integral(0,1.25);
  //  //
  //  std::vector<double> params = plus.return_cparam();
  //  std::vector<std::string> sparams = plus.return_sparam();
  //  //
  //  //
  //  Markov_Chain five = Markov_Chain(params, sparams, 1000000, "appearance10e6.root");

  Nu_Fitter tdb2 = Nu_Fitter(1, "/home/yong/Neutrinos/Prob3++.20121225/T2Kflux2016/t2kflux_2016_plus250kA.root", "enu_sk_numu","enu_sk_numub","enu_sk_nue","enu_sk_nueb");
  tdb2.set_param(8,0.5,'c');
  tdb2.make_sum('f','a',true);
  tdb2.make_sum('p','a',true);

   std::vector<double> params = tdb2.return_cparam();
   std::vector<std::string> sparams = tdb2.return_sparam();

  Markov_Chain t1 = Markov_Chain(params, sparams, 1000000, "tdb2.root");


   //
   //
   t1.set_param(3);//theta13
   t1.set_width(3,5E-3);
   t1.set_param(6);// delta
   t1.set_width(6,8E-2);
   t1.set_param(8); // beta
   t1.set_width(8,8E-2);
   //
   t1.startMH(params,&tdb2,'a');

   Nu_Fitter tdb3 = Nu_Fitter(1, "/home/yong/Neutrinos/Prob3++.20121225/T2Kflux2016/t2kflux_2016_plus250kA.root", "enu_sk_numu","enu_sk_numub","enu_sk_nue","enu_sk_nueb");
   tdb2.set_param(8,0.5,'c');
   tdb2.make_sum('f','a',true);
   tdb2.set_param(8,1,'p');
   tdb2.make_sum('p','a',true);

    std::vector<double> params1 = tdb3.return_cparam();
    std::vector<std::string> sparams1 = tdb3.return_sparam();

   Markov_Chain t2 = Markov_Chain(params1, sparams1, 1000000, "tdb3.root");


    //
    //
    t2.set_param(3);//theta13
    t2.set_width(3,5E-3);
    t2.set_param(6);// delta
    t2.set_width(6,8E-2);
    // t2.set_param(8); // beta
    // t2.set_width(8,8E-2);
    //
    t2.startMH(params1,&tdb3,'a');


  //  five.startMH(params,&plus,&minus);
  //  five.startMH(params, &foo, true);

}

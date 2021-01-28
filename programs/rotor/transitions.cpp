
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include "boost/math/constants/constants.hpp"
#include "eigen3/Eigen/Dense"
#include "fmt/format.h"
#include "am/wigner_gsl.h"
#include "am/am.h"
#include "lgi/lgi.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/vcs.h"
#include "sp3rlib/sp3r_operator.h"

#include <eigen3/Eigen/Eigen>

#define PI 3.14159265

namespace u3{
// HalfInt S(1), J0(2),J1(1);
  void get_transitions(const u3::SU3& x, const HalfInt& S, const HalfInt& J0, const HalfInt& J1,double hw, int A)
  {
    double hbarc=197.327; //MeVfm
    double mc2=938.92; //MeV 
    double b2=hbarc*hbarc/mc2/hw;
    double intrinsic_factor=2./A;
    // std::cout<<b2<<std::endl;

    double Qrme=b2*intrinsic_factor*sqrt(2*u3::Casimir2(x));
    // Construct the basis   
    MultiplicityTagged<int>::vector so3_branched_labels = BranchingSO3(x);
    //Bra J0, ket J1
    std::vector<int> basis0, basis1;
    for(auto L_kappa : so3_branched_labels)
      {
        int L=L_kappa.irrep;
        int kappa_max=L_kappa.tag;
        assert(kappa_max==1);
        if(am::AllowedTriangle(L,S,J0))
          basis0.push_back(L);

        if(am::AllowedTriangle(L,S,J1))
          basis1.push_back(L);
      }

    Eigen::MatrixXd Q(basis0.size(),basis1.size());
    for(int i=0; i<basis0.size(); ++i)
      for(int j=0; j<basis1.size(); ++j)
        {
          int L0 = basis0[i];
          int L1 = basis1[j];
          Q(i,j)=u3::W(x,1,L1,u3::SU3(1,1),1,2,x,1,L0,1)*am::Unitary9J(L1,S,J1,2,0,2,L0,S,J0)*Qrme;
        }
    std::cout<<J1.Str()<<"->"<<J0.Str()<<std::endl;
    std::cout<<Q<<std::endl<<std::endl;

  }
}




int main(int argc, char **argv)
{
  // if (argc<4)
  // {
  //   std::cout<<"Syntax: Z N <inputfile> <outputfile> \n <inputfile> contains list of subspace labels Nex 2Sp 2Sn 2S lambda mu"<<std::endl;
  // }
  double hw=15.0; //MeV
  
  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////
  // SU(3) caching
  u3::U3CoefInit();

  u3::SU3 x(2,1);
  int A=8;

  HalfInt S(1), J0(2),J1(2);
  u3::get_transitions(x,S,J0,J1,hw,A);
  //////////////////////////////
  J1=HalfInt(1);
  u3::get_transitions(x,S,J0,J1,hw,A);
  //////////////////////////////
  J1=HalfInt(3);
  u3::get_transitions(x,S,J0,J1,hw,A);

  ///////////////
  x=u3::SU3 (3,0);
  A=7;

  S=HalfInt(1,2), J0=HalfInt(3,2), J1=HalfInt(1,2);
  u3::get_transitions(x,S,J0,J1,hw,A);

  S=HalfInt(1,2), J0=HalfInt(3,2), J1=HalfInt(3,2);
  u3::get_transitions(x,S,J0,J1,hw,A);
}

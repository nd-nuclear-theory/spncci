/****************************************************************
  Calculate B(E2) transitions in SU(3) model analytically

  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  01/28/21 (aem) : Created
****************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "boost/math/constants/constants.hpp"
#include <Eigen/Dense>
#include "fmt/format.h"
#include "am/wigner_gsl.h"
#include "am/am.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "lgi/lgi.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/vcs.h"
#include "sp3rlib/sp3r_operator.h"

#include <Eigen/Eigen>


const double pi=boost::math::constants::pi<double>();

namespace u3{
  int wigner_sign(double x)
    {
      return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
    }

  void Q(const u3::SU3& x, const HalfInt& S, const HalfInt& J0, const HalfInt& J1,double hw)
  {
    double hbarc=197.327; //MeVfm
    double mc2=938.92; //MeV
    double b2=hbarc*hbarc/mc2/hw;

    double Qrme=sqrt(15./(16*pi))*b2*sqrt(2*u3::Casimir2(x));
    int sign=u3::wigner_sign(u3::W(x,1,x.lambda()+x.mu(),u3::SU3(1,1),1,1,x,1,x.lambda()+x.mu(),1));

    // Construct the basis
    MultiplicityTagged<unsigned int>::vector so3_branched_labels = BranchingSO3(x);
    MultiplicityTagged<unsigned int>::vector basis0, basis1;
    for(auto L_kappa : so3_branched_labels)
      {
        int L=L_kappa.irrep;
        int kappa_max=L_kappa.tag;
        assert(kappa_max==1); //Currently assume kappa_max=1 for all cases considered.
        if(am::AllowedTriangle(L,S,J0))
          for(int kappa=1; kappa<=kappa_max; ++kappa)
            basis0.emplace_back(L,kappa);

        if(am::AllowedTriangle(L,S,J1))
          for(int kappa=1; kappa<=kappa_max; ++kappa)
            basis1.emplace_back(L,kappa);
      }
    // Construct the matrix of Q
    Eigen::MatrixXd Q(basis0.size(),basis1.size());
    for(int i=0; i<basis0.size(); ++i)
      for(int j=0; j<basis1.size(); ++j)
        {
          int L0 = basis0[i].irrep;
          int kappa0=basis0[i].tag;
          int L1 = basis1[j].irrep;
          int kappa1=basis1[i].tag;
          Q(i,j)=sign*u3::W(x,kappa1,L1,u3::SU3(1,1),1,2,x,kappa0,L0,1)*am::Unitary9J(L1,S,J1,2,0,2,L0,S,J0)*Qrme;
        }

    std::cout<<J1.Str()<<"->"<<J0.Str()<<std::endl;
    std::cout<<mcutils::FormatMatrix(Q,"9.4f")<<std::endl<<std::endl;

  }

    void Qp(
      const u3::SU3& xp,const u3::SU3& xn,const u3::SU3& x,
      int rho, const HalfInt& S, const HalfInt& J0, const HalfInt& J1,double hw,
      Eigen::ArrayXXf& Q
      )
  {
    double hbarc=197.327; //MeVfm
    double mc2=938.92; //MeV
    double b2=hbarc*hbarc/mc2/hw;

    double Qrme=sqrt(15./(16*pi))*b2*sqrt(2*u3::Casimir2(xp));
    int sign=u3::wigner_sign(u3::W(xp,1,xp.lambda()+xp.mu(),u3::SU3(1,1),1,1,xp,1,xp.lambda()+xp.mu(),1));
    int rho1_max=u3::OuterMultiplicity(x,u3::SU3(1,1),x);

    // Construct the basis
    MultiplicityTagged<unsigned int>::vector so3_branched_labels = BranchingSO3(x);
    MultiplicityTagged<unsigned int>::vector basis0, basis1;
    for(auto L_kappa : so3_branched_labels)
      {
        int L=L_kappa.irrep;
        int kappa_max=L_kappa.tag;
        assert(kappa_max==1);
        if(am::AllowedTriangle(L,S,J0))
          for(int kappa=1; kappa<=kappa_max; ++kappa)
            basis0.emplace_back(L,kappa);

        if(am::AllowedTriangle(L,S,J1))
          for(int kappa=1; kappa<=kappa_max; ++kappa)
            basis1.emplace_back(L,kappa);
      }

    std::cout<<"Initial states"<<std::endl;
    for(const auto& state : basis1)
      std::cout<<state.tag<<"  "<<state.irrep<<std::endl;
    std::cout<<"Final states"<<std::endl;
    for(const auto& state : basis0)
      std::cout<<state.tag<<"  "<<state.irrep<<std::endl;

    std::cout<<std::endl;
    Q=Eigen::ArrayXXf::Zero(basis0.size(),basis1.size());
    for(int i=0; i<basis0.size(); ++i)
      for(int j=0; j<basis1.size(); ++j)
        {
          int L0 = basis0[i].irrep;
          int kappa0=basis0[i].tag;
          int L1 = basis1[j].irrep;
          int kappa1=basis1[i].tag;
          Q(i,j)=0;

          for(int rho1=1; rho1<=rho1_max; ++rho1)
          {
            Q(i,j)+=sign*u3::W(x,kappa1,L1,u3::SU3(1,1),1,2,x,kappa0,L0,rho1)
                    *u3::Z(u3::SU3(1,1),xp,x,xn,xp,1,rho,x,rho,rho1)
                    *am::Unitary9J(L1,S,J1,2,0,2,L0,S,J0)
                    *Qrme;

          }
        }
    std::cout<<J1.Str()<<"->"<<J0.Str()<<std::endl;
    std::cout<<mcutils::FormatMatrix(Q,"9.4f")<<std::endl<<std::endl;
  }
}




int main(int argc, char **argv)
{
  int max_lambda_plus_mu=39;
  u3::U3CoefInit(max_lambda_plus_mu);
  if (argc<2)
  {
    std::cout<<"Syntax: <inputfile>"<<std::endl;
    std::cout<<"  <inputfile> contains list of lambda_Z mu_Z lambda_N mu_N rho lambda mu 2S 2Ji 2Jf hw"<<std::endl;
    std::cout<<"  where (lambda_Z,mu_Z) x (lambda_N, mu_N) -> rho (lambda,mu)"<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::string inputfile=argv[1];

  std::vector<std::tuple<u3::SU3,u3::SU3,u3::SU3,int,HalfInt,HalfInt,HalfInt,double>> rme_labels;

  std::ifstream in_stream(inputfile);
  mcutils::StreamCheck(bool(in_stream),inputfile,"Failed to open input file");

  // scan input file
  std::string line;
  int line_count = 0;
  while ( std::getline(in_stream, line) )
    {
      // count line
      ++line_count;

      // set up for parsing
      std::istringstream line_stream(line);

      // parse line
      //   Nex lambda mu 2Sp 2Sn 2S count
      double hw;
      int lambda_Z, mu_Z, lambda_N, mu_N, rho, lambda, mu, SS, JJi, JJf;
      line_stream >>lambda_Z>> mu_Z>> lambda_N>> mu_N>> rho>> lambda>> mu>> SS>> JJi>> JJf >>hw;
      mcutils::ParsingCheck(line_stream, line_count, line);

      u3::SU3 x(lambda,mu),xp(lambda_Z,mu_Z), xn(lambda_N,mu_N);
      HalfInt S(SS,2), Ji(JJi,2), Jf(JJf,2);

      rme_labels.emplace_back(xp,xn,x,rho,S,Ji,Jf,hw);

    }

  // int rho;
  // u3::SU3 x,xp,xn;
  // HalfInt S,Ji,Jf;
  // double hw;
  for(auto& labels : rme_labels)
    {
      const auto& [xp,xn,x,rho,S,Ji,Jf,hw]=labels;
      std::cout<<xp.Str()<<"  "<<xn.Str()<<"  "<<x.Str()<<"  "
      <<rho<<"  "<<S<<"  "<<Ji<<"  "<<Jf<<"  "<<hw<<std::endl;
      Eigen::ArrayXXf Qp,Qn;
      std::cout<<"Qp"<<std::endl;
      u3::Qp(xp,xn,x,rho, S, Jf, Ji,hw,Qp);
      std::cout<<"Qn"<<std::endl;
      u3::Qp(xn,xp,x,rho, S, Jf, Ji,hw,Qn);
      Eigen::ArrayXXf ratio=Qp/Qn;
      std::cout<<"Ratio Qp/Qn"<<std::endl;
      std::cout<<mcutils::FormatMatrix(ratio,"9.4f")<<std::endl;
      std::cout<<"------------------------------------"<<std::endl;

    }





  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////
  // SU(3) caching


  if(false)
  {
    // // 8Li
    // u3::SU3 x(2,1), xp(3,0), xn(1,0);

    // HalfInt S(1), J0(2),J1(2);
    // u3::Q(x,S,J0,J1,hw);
    // u3::Qp(xp,xn,x,1, S, J0, J1,hw);
    // // //////////////////////////////
    // J1=HalfInt(1);
    // u3::Q(x,S,J0,J1,hw);
    // u3::Qp(xp,xn,x,1, S, J0, J1,hw);
    // // //////////////////////////////
    // J1=HalfInt(3);
    // u3::Q(x,S,J0,J1,hw);
    // u3::Qp(xp,xn,x,1, S, J0, J1,hw);
  }
  ///////////////
  if(false)
  {
    // 7Be
     double hw=15.0; //MeV
    u3::SU3 x(3,0),xp(2,0), xn(1,0);
    HalfInt S=HalfInt(1,2),J0,J1;
    // HalfInt S=HalfInt(1,2), J0=HalfInt(3,2), J1=HalfInt(1,2);
    // u3::Q(x,S,J0,J1,hw);
    // u3::Qp(xp,xn,x,1, S, J0, J1,hw);

    J0=HalfInt(3,2), J1=HalfInt(3,2);
    u3::Q(x,S,J0,J1,hw);

    Eigen::ArrayXXf Qp,Qn;
    u3::Qp(xp,xn,x,1, S, J0, J1,hw,Qp);
    u3::Qp(xn,xp,x,1, S, J0, J1,hw,Qn);
    Eigen::ArrayXXf ratio=Qp/Qn;
    std::cout<<ratio<<std::endl;
  }

  if(false)
  {
    // 9Be
     double hw=15.0; //MeV
    u3::SU3 x(3,1),xp(2,0), xn(1,1);
    HalfInt S=HalfInt(1,2),J0,J1;
    // HalfInt S=HalfInt(1,2), J0=HalfInt(3,2), J1=HalfInt(1,2);
    // u3::Q(x,S,J0,J1,hw);
    // u3::Qp(xp,xn,x,1, S, J0, J1,hw);

    J0=HalfInt(3,2), J1=HalfInt(3,2);
    u3::Q(x,S,J0,J1,hw);

    Eigen::ArrayXXf Qp,Qn;
    u3::Qp(xp,xn,x,1, S, J0, J1,hw,Qp);
    u3::Qp(xn,xp,x,1, S, J0, J1,hw,Qn);
    Eigen::ArrayXXf ratio=Qp/Qn;
    std::cout<<ratio<<std::endl;
  }




}

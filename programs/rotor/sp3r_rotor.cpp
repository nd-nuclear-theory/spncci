/****************************************************************
  Calculate rigid rotor spectrum in Sp(3,R) framework
  
  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  05/15/20 (aem) : Created
****************************************************************/

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

namespace sp3r
{
    Eigen::MatrixXd  QuadrupoleOperator(
      const sp3r::Sp3RSpace& sp3r_space, 
      const u3::U3& omegap, int kappap, int Lp,
      const u3::U3& omega, int kappa, int L,
      const vcs::MatrixCache& K_matrices
    )
    //Assumes oscillator length paramters b=1
    {
     const double pi=boost::math::constants::pi<double>();

      Eigen::MatrixXd Q;
      if(omega.N()+2 == omegap.N())
        {
          Q = u3::W(omega.SU3(),kappa,L,u3::SU3(2,0),1,2,omegap.SU3(),kappap, Lp,1)
              *sqrt(15.0/16.0/pi)*Sp3rRaisingOperator(sp3r_space, omegap, omega, K_matrices);
          std::cout<<"A "<<Sp3rRaisingOperator(sp3r_space, omegap, omega, K_matrices)<<std::endl;
        }
      if(omega.N() == omegap.N())
      {
        Q = u3::W(omega.SU3(),kappa,L,u3::SU3(1,1),1,2,omegap.SU3(),kappap, Lp,1)
            *sqrt(15.0/16.0/pi)*U3Operator(sp3r_space, omegap, omega);
            std::cout<<"C11 "<<U3Operator(sp3r_space, omegap, omega)<<std::endl;
      }
      if(omega.N()-2 == omegap.N())
        Q = u3::W(omega.SU3(),kappa,L,u3::SU3(0,2),1,2,omegap.SU3(),kappap, Lp,1)
            *sqrt(15.0/16/pi)*sp3r::Sp3rLoweringOperator(sp3r_space, omegap, omega, K_matrices);
      
      return Q; 
    }


 
}

int main(int argc, char **argv)
{
  // if (argc<4)
  // {
  //   std::cout<<"Syntax: Z N <inputfile> <outputfile> \n <inputfile> contains list of subspace labels Nex 2Sp 2Sn 2S lambda mu"<<std::endl;
  // }
  

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////
  // SU(3) caching
  u3::U3CoefInit();

  // Useful constants
  // not currently used 
  double hbarc=197.327; //MeVfm
  double mc2=938.92; //MeV 
  double hw=20.0; //MeV
  double b2=hbarc*hbarc/mc2/hw;
  std::cout<<b2<<std::endl;
  // coef=sqrt(15.0/16/pi)*b2;

  //7Be ground state band 
  lgi::NuclideType nuclei;
  nuclei[0]=4;
  nuclei[1]=3;
  HalfInt N0=lgi::Nsigma0ForNuclide(nuclei,true);
  std::cout<<"N0 "<<N0<<std::endl;

  u3::U3 sigma(N0,u3::SU3(3,0));
  HalfInt S(1,2);
  int Nmax=2;
  sp3r::Sp3RSpace sp3r_space(sigma,Nmax);
  vcs::MatrixCache K_matrices;
  vcs::GenerateKMatrices(sp3r_space,K_matrices);
  std::cout<<"Kmatrices "<<std::endl;
  for(auto it=K_matrices.begin(); it!=K_matrices.end(); ++it)
    {
        std::cout<<it->first.Str()<<std::endl;
        auto matrix=it->second;
        std::cout<<matrix<<std::endl;
        // std::cout<<matrix.inverse()<<std::endl;
    }

  std::vector<u3::U3> u3_irreps;
  u3_irreps.push_back(sigma);
  u3_irreps.emplace_back(N0+2,u3::SU3(5,0));

  std::map<std::tuple<HalfInt,int,HalfInt,int>,double> transitions;
  std::set<HalfInt> Jp_set;
  std::set<HalfInt> J_set;
  for(u3::U3 omega : u3_irreps)
    for(u3::U3 omegap : u3_irreps)
      {

        MultiplicityTagged<int>::vector L_kappa_values=u3::BranchingSO3(omega.SU3());
        MultiplicityTagged<int>::vector Lp_kappap_values=u3::BranchingSO3(omegap.SU3());
        for( auto L_kappa : L_kappa_values)
          for(auto Lp_kappap : Lp_kappap_values)
            {
              int Lp=Lp_kappap.irrep;
              int kappap_max=Lp_kappap.tag;
              int L=L_kappa.irrep;
              int kappa_max=L_kappa.tag;
              
              for(int kappa=1; kappa<=kappa_max; ++kappa)
                for(int kappap=1; kappap<=kappap_max; ++kappap)
                  {
                    HalfInt::vector J_values = am::ProductAngularMomenta(L, S);
                    HalfInt::vector Jp_values = am::ProductAngularMomenta(Lp, S);
                    Eigen::MatrixXd QRME = sp3r::QuadrupoleOperator(sp3r_space,omegap,kappap,Lp,omega,kappa,L,K_matrices);
                    
                    // Branching 
                    for(HalfInt J : J_values)
                      for(HalfInt Jp : Jp_values)
                        {
                          if(not am::AllowedTriangle(J,Jp,2))
                            continue; 

                          Eigen::MatrixXd Q = am::Unitary9J(L,S,J,2,0,2,Lp,S,Jp)*QRME;
                          int n = (omega == sigma)?1:2;
                          int np = (omegap == sigma)?1:2;
                          
                          J_set.insert(J);
                          
                          std::tuple<HalfInt,int,HalfInt,int> key(J,n,Jp,np);
                          
                          //Convert RME to Edmonds/Suhonen convention and store in map
                          // transitions[key]=Hat(Jp)*Q(0,0); //Single value for current irreps 
                          transitions[key]=Q(0,0); //Single value for current irreps 
                          
                        }
                  }
            }
      }

  // Generate quadrupole matrix elements for inband transitions in ground state band 
  for(HalfInt Ji : J_set)
    for(int n=1; n<3; n++)
      for(HalfInt Jf : J_set)
        for(int np=1; np<3; ++np)
          {
            if(am::AllowedTriangle(Ji,Jf,2) && (Ji>=Jf))
            {
              if((Ji==HalfInt(9,2) or Ji==HalfInt(11,2)) && n==1)
                continue;

              if((Jf==HalfInt(9,2) or Jf==HalfInt(11,2)) && np==1)
                continue;
              
              std::tuple<HalfInt,int,HalfInt,int> key(Ji,n,Jf,np);
              
              double E2=transitions.count(key)?transitions.at(key):std::nanf("");
              // Converting to shell model state indexing from irrep indexing 
              int ni=n;
              int nf=np;
              if(Ji==HalfInt(9,2) or Ji==HalfInt(11,2))
                ni=1;
              if(Jf==HalfInt(9,2) or Jf==HalfInt(11,2))
                nf=1;

              std::cout<<fmt::format("{:5.1f}  1  {}  nan   {:5.1f}  1  {}  nan    nan   nan   {: 6.5}   nan",float(Jf),nf,float(Ji),ni,E2)<<std::endl;
            }
          }

std::cout<<"--------------------------------"<<std::endl;

  for(HalfInt Jp : J_set)
    for(HalfInt J : J_set)
      {
        if((J==HalfInt(9,2) or J==HalfInt(11,2)))
          continue;

        if((Jp==HalfInt(9,2) or Jp==HalfInt(11,2)))
          continue;


        std::tuple<HalfInt,int,HalfInt,int> key1(Jp,2,J,2);
        std::tuple<HalfInt,int,HalfInt,int> key2(Jp,1,J,1);

        double trans1=transitions[key1];
        double trans2=transitions[key2];

        double BE2_lower=trans2*trans2/am::dim(J);
        double BE2_upper=trans1*trans1/am::dim(J);

        if(am::AllowedTriangle(J,Jp,2) && (J>=Jp))
          std::cout<<fmt::format("{:5}   {:5}  {: 6.5}    {: 6.5}    {: 6.5}",Jp.Str(),J.Str(),BE2_upper,BE2_lower,BE2_upper/BE2_lower)<<std::endl;
      }

std::cout<<"--------------------------------"<<std::endl;
 for(HalfInt Jp : J_set)
    for(HalfInt J : J_set)
      {
        
        if((Jp==HalfInt(9,2) or Jp==HalfInt(11,2)))
          continue;


        std::tuple<HalfInt,int,HalfInt,int> key1(Jp,2,J,1);
        std::tuple<HalfInt,int,HalfInt,int> key2(Jp,2,J,2);

        double trans1=transitions[key1];
        double trans2=transitions[key2];

        double BE2_inter=trans2*trans2/am::dim(J);
        double BE2_intra=trans1*trans2/am::dim(J);

        if(am::AllowedTriangle(J,Jp,2) && (J>=Jp))
          std::cout<<fmt::format("{:5}   {:5}  {: 6.5}    {: 6.5}    {: 6.5}",Jp.Str(),J.Str(),BE2_intra,BE2_inter,BE2_intra/BE2_inter)<<std::endl;
      }
        


}

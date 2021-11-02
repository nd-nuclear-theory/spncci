
/****************************************************************
  Calculate triaxial rotor spectrum in SU(3) basis

  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  05/15/20 (aem) : Created
****************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <Eigen/Dense>
#include "fmt/format.h"
#include "am/wigner_gsl.h"
#include "am/am.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"

#define PI 3.14159265

namespace rotor
{

  double GetGamma(const u3::SU3& x)
  {return atan(sqrt(3)*(x.mu()+1)/(2*x.lambda()+x.mu()+3));}

  void GetA(double alpha, double gamma, std::vector<double>&A)
  // Based on np-8-1958-237-Davydov
  {
    A.resize(3);
    for(int k=1; k<=3; ++k)
      {
        A[k-1]=alpha/mcutils::sqr(sin(gamma-2*PI/3*k));
      }
  }


  void GetLambda(const u3::SU3& x, std::vector<double>& lambda)
  // Based on eq4.5 in zpa-329-1988-33-Castanos
    {

      lambda.resize(3);
      lambda[0]=-(x.lambda()-x.mu())/3;
      lambda[1]=-(x.lambda()+2*x.mu()+3)/3;
      lambda[2]=(2*x.lambda()+x.mu()+3)/3;
    }

  void GetD(const u3::SU3& x, const std::vector<double>& lambda, std::vector<double>& D)
  // Based on eq5.5 in zpa-329-1988-33-Castanos
    {
      D.resize(3);
      for(int i=0; i<3; ++i)
        D[i]=2*lambda[i]*lambda[i]*lambda[i]+lambda[0]*lambda[1]*lambda[2];
    }

  void GetCoefs(const u3::SU3& x, const std::vector<double>& A, std::vector<double>& coefs)
  // Based on eq5.5 in zpa-329-1988-33-Castanos
    {
      assert(A.size()==3);
      coefs.resize(3);

      std::vector<double> lambda;
      rotor::GetLambda(x,lambda);

      std::vector<double> D;
      GetD(x,lambda,D);

      //First coefficient
      coefs[0]=lambda[1]*lambda[2]/(2*mcutils::sqr(lambda[0])+lambda[1]*lambda[2])*A[0];
      coefs[0]=coefs[0]+lambda[0]*lambda[2]/(2*mcutils::sqr(lambda[1])+lambda[0]*lambda[2])*A[1];
      coefs[0]=coefs[0]+lambda[0]*lambda[1]/(2*mcutils::sqr(lambda[2])+lambda[0]*lambda[1])*A[2];

      //Second coefficient
      coefs[1]=lambda[0]/(2*mcutils::sqr(lambda[0])+lambda[1]*lambda[2])*A[0];
      coefs[1]=coefs[1]+lambda[1]/(2*mcutils::sqr(lambda[1])+lambda[0]*lambda[2])*A[1];
      coefs[1]=coefs[1]+lambda[2]/(2*mcutils::sqr(lambda[2])+lambda[0]*lambda[1])*A[2];

      coefs[2]=1/(2*mcutils::sqr(lambda[0])+lambda[1]*lambda[2])*A[0];
      coefs[2]=coefs[2]+1/(2*mcutils::sqr(lambda[1])+lambda[0]*lambda[2])*A[1];
      coefs[2]=coefs[2]+1/(2*mcutils::sqr(lambda[2])+lambda[0]*lambda[1])*A[2];


    }

  // double casimirSU3=Casimir2(x);
  // std::cout<<casimirSU3<<std::endl<<std::endl;
    double CRME(const u3::SU3& x)
      {
        double rme=sqrt(2*Casimir2(x));
        return rme;
      }

    double X3RME(const u3::SU3& x, int L, int kappap, int kappa)
    // [LxQ].L
      {

        // std::cout<<"       "<<L<<"  "<<kappa<<"  "<<kappap<<"    "<<am::Unitary6J(L,2,L,1,L,1)<<"  "
        //             <<u3::W(x,kappa,L,u3::SU3(1,1),1,2,x,kappap,L,1)<<"  "
        //             <<CRME(x)<<std::endl;


        double rme=-sqrt(3)*L*(L+1)*am::Unitary6J(L,2,L,1,L,1)
                    *u3::W(x,kappa,L,u3::SU3(1,1),1,2,x,kappap,L,1)
                    *CRME(x);
        return rme;
      }


    double X4RME(const u3::SU3& x, int L, int kappap, int kappa)
    // [LxQ].[QxL]
      {
        // std::cout<<std::endl<<"Calculating X4"<<std::endl;
        // std::cout<<kappap<<"  "<<kappa<<std::endl;
        MultiplicityTagged<unsigned int>::vector Lvalues=BranchingSO3(x);
        double rme=0.0;
        for(auto L_tagged : Lvalues)
          {
            int Lbar=L_tagged.irrep;
            int kappa_bar_max=L_tagged.tag;
            for(int kappa_bar=1; kappa_bar<=kappa_bar_max; ++kappa_bar)
              {

                for(int m=1; m<=1; ++m)
                {

                // // std::cout<<"Lbar "<<Lbar<<"  "<<kappa_bar<<std::endl;
                // std::cout<<"    "
                // <<am::Unitary6J(L,m,Lbar,m,L,0)<<"  "
                // <<ParitySign(L+Lbar+m)*Hat(Lbar)/Hat(L)/Hat(m)<<std::endl;
                // // <<mcutils::sqr(am::Unitary6J(L,1,L,2,Lbar,m))<<"  "
                // std::cout<<"SU3 coefs "
                // <<L<<"  "<<Lbar<<"  "<<m<<std::endl
                // <<kappa<<" "<<kappa_bar<<"  "<<kappap<<std::endl
                // <<"  "<<u3::W(x,kappa_bar,Lbar,u3::SU3(1,1),1,2,x,kappap,L,1)<<"  "
                // <<"  "<<u3::W(x,kappa,L,u3::SU3(1,1),1,2,x,kappa_bar,Lbar,1)<<std::endl;


                rme+=ParitySign(m)*Hat(m)
                      *am::Unitary6J(L,m,Lbar,m,L,0)*L*(L+1)
                      *mcutils::sqr(am::Unitary6J(L,1,L,2,Lbar,m))
                      *u3::W(x,kappa_bar,Lbar,u3::SU3(1,1),1,2,x,kappap,L,1)
                      *u3::W(x,kappa,L,u3::SU3(1,1),1,2,x,kappa_bar,Lbar,1);
                }
              }
          }
         return rme*6*Casimir2(x);
      }


  void GetHamiltonianMatrix(
    const u3::SU3& x, const MultiplicityTagged<unsigned int>::vector& basis,
    Eigen::MatrixXd& Hrot, std::vector<double>coefs
  )
    {

      Hrot=Eigen::MatrixXd::Zero(basis.size(),basis.size());
      for(int i=0; i<basis.size(); ++i)
        {
          auto& state=basis[i];
          int L=state.irrep;
          int kappa=state.tag;
          for(int ip=0; ip<=i; ++ip)
            {
              auto& statep=basis[ip];
              int Lp=statep.irrep;
              int kappap=statep.tag;
              if(Lp!=L)
                continue;

              double ham=coefs[0]*L*(L+1)
                          +coefs[1]*rotor::X3RME(x,L,kappap,kappa)
                          +coefs[2]*rotor::X4RME(x, L,kappap,kappa);

              Hrot(i,ip)=ham;
              Hrot(ip,i)=ham;
            }
        }
    }

  //Fit parameter
  // double alpha=.108; (split=1)
// No spliting over irreps
void GetSU3RotorEnergies(const u3::SU3& x, double alpha)
{
  //Construct basis for irrep by branching to SO(3)
  MultiplicityTagged<unsigned int>::vector Lvalues=BranchingSO3(x);
  MultiplicityTagged<unsigned int>::vector basis;
  for(auto L_tagged : Lvalues)
    {
      int L=L_tagged.irrep;
      int kappa_max=L_tagged.tag;
      for(int kappa=1; kappa<=kappa_max; ++kappa)
        {
          basis.emplace_back(L,kappa);
          std::cout<<L<<"  "<<kappa<<std::endl;
        }
    }

  //Get angle based on SU(3) labels
  double gamma=rotor::GetGamma(x);
  std::cout<<"gamma is: "<<gamma/PI<<"Pi "<<std::endl;

  // Calculate coefficients
  std::vector<double> A(3);

  // A is based on Daydov model expression for moment of intertia
  rotor::GetA(alpha,gamma,A);
  // std::cout<<"A coefs"<<std::endl;
  // for(auto a : A)
  //   std::cout<<a<<std::endl;

  //Compute coefficients for terms in the Hamiltonian
  std::vector<double> coefs;
  rotor::GetCoefs(x,A,coefs);

  // std::cout<<"Rotor coefs"<<std::endl;
 //  for(auto coef : coefs)
 //    std::cout<<coef<<std::endl;
 // std::cout<<"-------"<<std::endl;

  // Construct matrix of Hrot_su3
  Eigen::MatrixXd Hrot;
  rotor::GetHamiltonianMatrix(x,basis,Hrot,coefs);


  // std::cout<<Hrot1<<std::endl<<std::endl;

  //Get eigenvalues
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Hrot);

  int eigensolver_status = eigensolver.info();
  assert(eigensolver_status==Eigen::Success);

  // save eigenresults
  auto eigenvalues = eigensolver.eigenvalues();
  auto eigenvectors = eigensolver.eigenvectors();

  std::cout<<eigenvalues<<std::endl<<std::endl;

  std::cout<<eigenvectors<<std::endl;

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
  int max_lambda_plus_mu=39;
  u3::U3CoefInit(max_lambda_plus_mu);
  double alpha=1.0;
  // If we assum a maximially asymmetric (triaxial) rotor with lambda=mu
  if(true)
  { u3::SU3 x(2,2);
    std::cout<<x.Str()<<std::endl;
    rotor::GetSU3RotorEnergies(x, alpha);
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
  }
  if(true)
  {
    u3::SU3 x(3,0);
    std::cout<<x.Str()<<std::endl;
    rotor::GetSU3RotorEnergies(x, alpha);
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
  }

  if(true)
  {
    u3::SU3 x(5,0);
    std::cout<<x.Str()<<std::endl;
    rotor::GetSU3RotorEnergies(x, alpha);
    std::cout<<"-------------------------------------------------------------------------"<<std::endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////


//   u3::SU3 x2(3,0);
//   double alpha2=alpha*(1-split);
//   double gamma2=rotor::GetGamma(x2);
//   std::cout<<"gamma is: "<<gamma2/PI<<"Pi "<<std::endl;

//   std::vector<double> A2(3);
//   rotor::GetA(alpha2,gamma2,A2);
//   std::cout<<"A coefs"<<std::endl;
//   for(auto a : A2)
//     std::cout<<a<<std::endl;

//   std::cout<<"-------"<<std::endl;
//   std::vector<double> coefs2;
//   rotor::GetCoefs(x2,A2,coefs2);
//   for(auto coef : coefs2)
//     std::cout<<coef<<std::endl;
//   std::cout<<"-------"<<std::endl;

//     // Construct matrix of Hrot_su3
//   Eigen::MatrixXd Hrot2;
//   rotor::GetHamiltonianMatrix(x,basis,Hrot2,coefs2);


//   std::cout<<Hrot2<<std::endl<<std::endl;

//   //Get eigenvalues
//   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver2(Hrot2);

//   int eigensolver_status2 = eigensolver2.info();
//   assert(eigensolver_status2==Eigen::Success);

//   // save eigenresults
//   auto eigenvalues2 = eigensolver2.eigenvalues();
//   auto eigenvectors2 = eigensolver2.eigenvectors();

//   std::cout<<eigenvalues2<<std::endl<<std::endl;

//   std::cout<<eigenvectors2<<std::endl;

//   std::cout<<"-------------------------------------------------------------------------"<<std::endl;


//   Eigen::MatrixXd Hrot=Hrot1+Hrot2;

//   std::cout<<Hrot<<std::endl<<std::endl;

//   //Get eigenvalues
//   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Hrot);

//   int eigensolver_status = eigensolver.info();
//   assert(eigensolver_status==Eigen::Success);

//   // save eigenresults
//   auto eigenvalues = eigensolver.eigenvalues();
//   auto eigenvectors = eigensolver.eigenvectors();

//   std::cout<<eigenvalues<<std::endl<<std::endl;

//   std::cout<<eigenvectors<<std::endl;

// std::cout<<"-------------------------------------------------------------------------"<<std::endl;
// std::cout<<"Testing "<<std::endl;


}

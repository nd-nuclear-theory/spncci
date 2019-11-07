
//  Calculate triaxial rotor spectrum in SU(3) basis

#include <iostream>
#include <fstream>
#include <cstdlib>


#include "eigen3/Eigen/Dense"
#include "fmt/format.h"
#include "am/wigner_gsl.h"
#include "am/am.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"


namespace rotor
{

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


      coefs[0]=0; 
      for(int i=0; i<3; i++)
        {
          if(D[i]==0)
            continue;
          // std::cout<<coefs[0]<<std::endl;
          // std::cout<<"D "<<D[i]<<"  lambda "<<lambda[i]<<"  A "<<A[i]<<std::endl;
          // std::cout<<"lambdas "<<lambda[0]*lambda[1]*lambda[2]/D[i<<std::endl;
          coefs[0]=coefs[0]+lambda[0]*lambda[1]*lambda[2]/D[i]*A[i];
          // std::cout<<coefs[0]<<std::endl;
        }

      coefs[1]=0;
      for(int i=0; i<3; i++)
        {
          if(D[i]==0)
            continue;
          coefs[1]+=lambda[i]*lambda[i]/D[i]*A[i];
        }

      coefs[2]=0;
      for(int i=0; i<3; i++)
        {
          if(D[i]==0)
            continue;
          coefs[2]+=lambda[i]/D[i]*A[i];
        }
    }

  // double casimirSU3=Casimir2(x);
  // std::cout<<casimirSU3<<std::endl<<std::endl;
    double CRME(const u3::SU3& x)
      {
        double rme=sqrt(2*Casimir2(x));
        return rme;
      }

    double X3RME(const u3::SU3& x, int L, int kappap, int kappa)
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
      {
        // std::cout<<std::endl<<"Calculating X4"<<std::endl;
        // std::cout<<kappap<<"  "<<kappa<<std::endl;
        MultiplicityTagged<int>::vector Lvalues=BranchingSO3(x);
        double rme=0.0;
        for(auto L_tagged : Lvalues)
          {
            int Lbar=L_tagged.irrep;
            int kappa_bar_max=L_tagged.tag; 
            for(int kappa_bar=1; kappa_bar<=kappa_bar_max; ++kappa_bar)
              {
                
                for(int m=1; m<=3; ++m)
                {
                // std::cout<<"Lbar "<<Lbar<<"  "<<kappa_bar<<std::endl;
                // std::cout<<"    "
                // <<sqr(am::Unitary6J(L,1,L,2,Lbar,m))<<"  "
                // <<u3::W(x,kappa_bar,Lbar,u3::SU3(1,1),1,2,x,kappap,L,1)<<"  "
                // <<u3::W(x,kappa,L,u3::SU3(1,1),1,2,x,kappa_bar,Lbar,1)<<std::endl;

                rme+=ParitySign(L+Lbar+1)*Hat(Lbar)/Hat(L)*L*(L+1)
                      *sqr(am::Unitary6J(L,1,L,2,Lbar,m))
                      // *4*(x.lambda()+x.mu()+3)*(x.lambda()+x.mu())-x.lambda()*x.mu()
                      // *6*Casimir2(x)
                      *u3::W(x,kappa_bar,Lbar,u3::SU3(1,1),1,2,x,kappap,L,1)
                      *u3::W(x,kappa,L,u3::SU3(1,1),1,2,x,kappa_bar,Lbar,1);
                }
              }
          }
         return rme*6*Casimir2(x);
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

  u3::SU3 x(2,2);
  MultiplicityTagged<int>::vector Lvalues=BranchingSO3(x);
  MultiplicityTagged<int>::vector basis;
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

  // Calculate coefficients
  std::vector<double> A(3);
  // for(int i=0; i<3; ++i)
  //   A[i]=rand()%10;

  A[0]=0;
  A[1]=-0.17;
  A[2]=-0.01;
  double a=-.065;
  double b=-.01;

  std::vector<double> coefs;
  rotor::GetCoefs(x,A,coefs);

  for(auto coef : coefs)
    std::cout<<coef<<std::endl<<std::endl;

  // Construct matrix of Hrot_su3
  Eigen::MatrixXd Hrot;
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

        // double ham=coefs[0]*L*(L+1)
        //             +coefs[1]*rotor::X3RME(x,L,kappap,kappa)
        //             +coefs[2]*rotor::X4RME(x, L,kappap,kappa);


        double ham=a*rotor::X3RME(x,L,kappap,kappa)
                    +b*rotor::X4RME(x, L,kappap,kappa);


        // std::cout<<L<<"  "<<kappap<<"  "<<kappa<<"  "<<rotor::X3RME(x,L,kappap,kappa)<<"  "
        //   <<rotor::X4RME(x, L,kappap,kappa)<<std::endl;

        Hrot(i,ip)=ham;
        Hrot(ip,i)=ham;
      }
    }
  for(auto a : A)
    std::cout<<a<<std::endl;
  std::cout<<std::endl;

  std::cout<<Hrot<<std::endl<<std::endl;

  // spncci::Vector& eigenvalues,
  // spncci::Matrix& eigenvectors,
  
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Hrot);

    int eigensolver_status = eigensolver.info();
    assert(eigensolver_status==Eigen::Success);

    // save eigenresults
    auto eigenvalues = eigensolver.eigenvalues();
    auto eigenvectors = eigensolver.eigenvectors();

    std::cout<<eigenvalues<<std::endl<<std::endl;

    std::cout<<eigenvectors<<std::endl;




// u3::SU3 x1(1,1);
// double c2=(x.lambda()+x.mu()+3)*(x.lambda()+x.mu())-x.lambda()*x.mu();
// double me=-1*u3::W(x,1,3,u3::SU3(1,1),1,1,x,1,3,1)*2*sqrt(c2/3);
// std::cout<<"me "<<me<<"  "<<sqrt(3*4)/me<<std::endl;

// // double sum=0;
// double sum1=0;
// for(auto state : basis)
//   {
//     int L=state.irrep;
//     int k=state.tag;
    
    
//     for(int L1=1; L1<=2; L1++)
//     // for(statep : basis)
//       {
//         // int Lp=statep.irrep;
//         // int kp=statep.tag;
//         sum1+=u3::W(x,k,L,x1,1,L1,x,1,3,1)*u3::W(x,k,L,x1,1,L1,x,1,3,1);
//         // sum1+=u3::W(x,k,L,x1,1,L1,x,1,3,1)*u3::W(x,k,L,x1,1,L1,x,1,3,1);

//       }
    
//   }
// std::cout<<"  sum: "<<sum1<<std::endl;
// std::cout<<u3::W(x,1,3,u3::SU3(1,1),1,2,x,1,3,1)<<std::endl;

}

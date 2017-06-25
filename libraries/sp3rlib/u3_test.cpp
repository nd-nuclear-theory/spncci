/****************************************************************
  u3_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created.
  3/8/16 (aem,mac): Add tests for U3ST.
  3/9/16 (aem,mac): Update includes and namespaces.

****************************************************************/

#include "sp3rlib/u3.h"

#include <vector>
#include <algorithm>

#include "am/am.h"

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // construction and string conversion
  ////////////////////////////////////////////////////////////////

  std::cout << u3::SU3(1,2).Str() << std::endl;
  std::cout 
    << u3::U3(0,0,0).Str()
    << u3::U3(HalfInt(1,2),HalfInt(1,2),HalfInt(1,2)).Str()
    << u3::U3(HalfInt(5,2),u3::SU3(1,0)).Str()
    << std::endl;


  ////////////////////////////////////////////////////////////////
  // comparison operators
  ////////////////////////////////////////////////////////////////

  // SU(3)
  std::vector<u3::SU3> lm_vector;
  lm_vector.push_back(u3::SU3(1,2));
  lm_vector.push_back(u3::SU3(4,3));
  lm_vector.push_back(u3::SU3(1,1));
  lm_vector.push_back(u3::SU3(4,5));
  sort(lm_vector.begin(),lm_vector.end());
  for (int i=0; i<lm_vector.size(); ++i)
    {
      std::cout << lm_vector[i].Str();
    }
  std::cout << std::endl;

  // U(3)
  std::vector<u3::U3> omega_vector;
  omega_vector.push_back(u3::U3(1,1,1));
  omega_vector.push_back(u3::U3(3,1,0));
  omega_vector.push_back(u3::U3(3,2,1));
  omega_vector.push_back(u3::U3(2,1,1));
  for (int i=0; i<omega_vector.size(); ++i)
    {
      std::cout << omega_vector[i].Str();
    }
  std::cout << std::endl;
  sort(omega_vector.begin(),omega_vector.end());
  for (int i=0; i<omega_vector.size(); ++i)
    {
      std::cout << omega_vector[i].Str();
    }
  std::cout << std::endl;

  // U(3)xSU(2)
  std::vector<u3::U3S> omegaS_vector;
  omegaS_vector.push_back(u3::U3S(u3::U3(1,1,1),3));
  omegaS_vector.push_back(u3::U3S(u3::U3(1,1,1),1));
  omegaS_vector.push_back(u3::U3S(u3::U3(3,1,0),HalfInt(3,2)));
  omegaS_vector.push_back(u3::U3S(u3::U3(3,1,0),HalfInt(1,2)));
  for (int i=0; i<omegaS_vector.size(); ++i)
    {
      std::cout << omegaS_vector[i].Str();
    }
  std::cout << std::endl;
  sort(omegaS_vector.begin(),omegaS_vector.end());
  for (int i=0; i<omegaS_vector.size(); ++i)
    {
      std::cout << omegaS_vector[i].Str();
    }
  std::cout << std::endl;

  // U(3)xSU(2)xSU(2)
  std::vector<u3::U3ST> omegaST_vector;
  omegaST_vector.push_back(u3::U3ST(u3::U3(1,1,1),3,2));
  omegaST_vector.push_back(u3::U3ST(u3::U3(1,1,1),1,1));
  omegaST_vector.push_back(u3::U3ST(u3::U3(3,1,0),HalfInt(3,2),HalfInt(1,2)));
  omegaST_vector.push_back(u3::U3ST(u3::U3(3,1,0),HalfInt(1,2),1));
  for (int i=0; i<omegaST_vector.size(); ++i)
    {
      std::cout << omegaST_vector[i].Str()<<"  ";
    }
  std::cout << std::endl;
  sort(omegaST_vector.begin(),omegaST_vector.end());
  for (int i=0; i<omegaST_vector.size(); ++i)
    {
      std::cout << omegaST_vector[i].Str()<<"  ";
    }
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // dimension, conjugation, and validation tests
  ////////////////////////////////////////////////////////////////
  std::cout<<"  "<<std::endl;
  std::cout <<"dim"<<" "<<"Conjugate"<<" "<<"ConjugationGrade"<<" "<<"Casimir2"<<" "<<"Casimir3"<<std::endl;
  for (int i=0; i<lm_vector.size(); ++i)
    {
      std::cout 
        << lm_vector[i].Str() 
        << " " << dim(lm_vector[i])
        << " " << Conjugate(lm_vector[i]).Str() 
        << " " << ConjugationGrade(lm_vector[i])
        << " " << Casimir2(lm_vector[i])
        << " " << Casimir3(lm_vector[i])
        << std::endl;
    }
  std::cout << std::endl;

  for (int i=0; i<omega_vector.size(); ++i)
    {
      std::cout 
        << omega_vector[i].Str() 
        << " " << dim(omega_vector[i])
        << " " << Conjugate(omega_vector[i]).Str() 
        << " " << ConjugationGrade(omega_vector[i])
        << " " << omega_vector[i].Valid()
        << std::endl;
    }
  std::cout << std::endl;
  std::cout 
    << u3::U3(1,3,2).Str() 
    << " " << u3::U3(1,3,2).Valid() 
    << std::endl;

  for (int i=0; i<omegaS_vector.size(); ++i)
    {
      std::cout 
        << omegaS_vector[i].Str() 
        << " " << dim(omegaS_vector[i])
        << std::endl;
    }
  std::cout << std::endl;

  ///////////////////////////////////////////////////////////////////
  //  Testing KronckerProduct, OuterMultiplicity, BranchingMultiplicity
  //  and BranchingSO3
  //  Checked against prototype Sp(3,R) libary function couple in u3.py
  //////////////////////////////////////////////////////////////////
  MultiplicityTagged<u3::SU3>::vector product;
  MultiplicityTagged<int>::vector branch;
  MultiplicityTagged<int>::vector branchR;
  HalfInt S(1,2);  // HalfInt S = HalfInt(1,2);
  HalfInt J(1,2);  // HalfInt J = HalfInt(1,2);
  HalfInt::pair r = am::ProductAngularMomentumRange(S,J);
  // Coupling (lambda1,mu1)x(lambda2,mu2) for range of lambda's and mu's 
  for(int l1=0; l1<2; l1++)
    for(int m1=1; m1<3; m1++)
      {
        u3::SU3 x1(l1,m1);
        for(int l2=3; l2<5; l2++)
          for(int m2=4; m2<5; m2++)
            {
              u3::SU3 x2(l2,m2);
              //List of product irreps from coupling
              product=KroneckerProduct(x1,x2);
              // Print (lambda1,mu1) and (lambda2,mu2)
              std::cout<<x1.Str()<<"  "<<x2.Str()<<std::endl;
              // for each product irrep, branch irrep and print product irrep 
              // followed by possible kappa L values 
              for(int i=0; i<product.size(); i+=5)  
                {
                  branch=BranchingSO3(product[i].irrep);
                  std::cout << "  " << product[i].Str() << std::endl;
                  std::cout << "Unconstrained branching" << std::endl;
                  for(int j=0; j<branch.size();j++)
                    std::cout << "    " << branch[j].irrep << "," << branch[j].tag << std::endl;
                  branchR=BranchingSO3Constrained(product[i].irrep,r);
                  std::cout << "Constrained branching" <<std::endl;
                  for(int k=0; k<branchR.size();k++)
                    std::cout << "    " << branchR[k].irrep << "," << branchR[k].tag << std::endl;
                }
            }
      }

} //main

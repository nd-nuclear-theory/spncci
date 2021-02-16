/****************************************************************
  sp3rcoef.cpp

  Anna E. McCoy
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>

#include "gsl/gsl_sf.h"  

#include "sp3rlib/sp3rcoef.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h"
 
namespace sp3r
  {
	
  // CCoef(u3::U3& n, int q, u3:: U3& np)
  // //Don't need at present.  Will come back and finish later (aem). 

  // // Coefficient of fractional parentage for expanding the raising polynomial into a
  // // polynomial of a single Jacobi coordinate and a polynomial of all other coordinates 
  // {
  //   double coef=0;
  //   int N=int(n.N());
  //   if (n.SU3()==u3::SU3(N,0))&&(np.SU3()==u3::SU3(N-q,0))
  //     coef=sqrt(
  //       gsl_sf_fact(N/2)*gsl_sf_fact(q)
  //       /(pow(2.,q/2)*gsl_sf_fact(N/2-q/2)*gsl_sf_fact(q/2))
  //       );
  //   else
  //   {
  //     if ((n.f1-2)>=n.f2)
  //       u3::U3 nb(n.f1-2,n.f2,n.f3);
  //     else if ((n.f2-2)>=f.n3)
  //       u3::U3 nb(n.f1,n.f2-2,n.f3);
  //     else if ((n.f3-2)>=0)
  //       u3::U3 nb(n.f1,n.f2,n.f3-2);
  //     // if none of the above returns true, coef is zero. 
  //     else
  //       return coef;
  //     // otherwise, calculate coef
  //  }
  // 
  
  void GenerateBCoefCache(BCoefCache& cache, int Nmax)
    {
      // TODO: Finish debugging
      std::vector<u3::U3> poly_vector=RaisingPolynomialLabels(Nmax);
      for(auto n1 : poly_vector)
        {
          MultiplicityTagged<u3::SU3>::vector n1p_set=u3::KroneckerProduct(n1.SU3(),u3::SU3(0,2));
          int N1p=int(n1.N()-2);
          u3::U3 n1p;
          bool continue_flag=true;

          for(int i=0; i<n1p_set.size(); ++i)
            {
              n1p=u3::U3(N1p,n1p_set[i].irrep);
              if((n1p.Valid())&&(n1p.SU3().lambda()%2==0)&&(n1p.SU3().mu()%2==0))
                {    
                  continue_flag=false;
                  break;
                }
            }
          if(continue_flag)
            continue;
          double coef1=vcs::BosonCreationRME(n1,n1p);
          for(auto n2 : poly_vector)
            {
              MultiplicityTagged<u3::SU3>::vector n3_set=u3::KroneckerProduct(n1.SU3(),n2.SU3());
              MultiplicityTagged<u3::SU3>::vector n3p_set=u3::KroneckerProduct(n1p.SU3(),n2.SU3());
              int N3=int(n1.N()+n2.N());
              int N3p=N3-2;
              for(auto n3_tagged :n3_set)
                {
                  if(not u3::U3::ValidLabels(N3, n3_tagged.irrep))
                    continue;
                  u3::U3 n3(N3,n3_tagged.irrep);
                  if((n3.SU3().lambda()%2!=0)||(n3.SU3().mu()%2!=0))
                    continue;
                  // If n1 or n2 are identity, the coef is 1
                  if((n1==u3::U3(0,0,0))||(n2==u3::U3(0,0,0)))
                    {
                      cache[BCoefLabels(n1,n2,n3,1)]=1;
                      continue;
                    }
                  int rho_max(n3_tagged.tag);
                  for(int rho=1; rho<=rho_max; ++rho)
                    {
                  // std::cout<<n1.Str()<<"  "<<n2.Str()<<"  "<<n3.Str()<<"  "<<rho<<std::endl;
                      BCoefLabels labels(n1,n2,n3,rho);
                      double coef2=0;
                      for(auto n3p_tagged : n3p_set)
                        {
                          u3::U3 n3p(N3p,n3p_tagged.irrep);
                          if(not n3p.Valid())
                            continue;
                          if((n3p.SU3().lambda()%2!=0)||(n3p.SU3().mu()%2!=0))
                            continue;
                          if(u3::OuterMultiplicity(n3p.SU3(),u3::SU3(2,0),n3.SU3())==0)
                            continue;
                          double coef3=vcs::BosonCreationRME(n3,n3p);
                          int rhop_max=n3p_tagged.tag;
                          for(int rhop=1; rhop<=rhop_max; ++rhop)
                            {
                              BCoefLabels labelsp(n1p,n2,n3p,rhop);
                              double coef=coef3/coef1*cache[labelsp]
                                *u3::U(u3::SU3(2,0),n1p.SU3(),n3.SU3(),n2.SU3(),n1.SU3(),1,rho,n3p.SU3(),rhop,1);
                              // std::cout<<"  "<<n1p.Str()<<"  "<<n2.Str()
                              // <<"  "<<n3p.Str()<<"  "<<rhop<<"  "<<coef1<<"  "<<coef3<<"  "<<cache[labelsp]
                              // <<"  "<<u3::U(u3::SU3(2,0),n1p.SU3(),n3.SU3(),n2.SU3(),n1.SU3(),1,rho,n3p.SU3(),rhop,1)
                              // <<"  "<<coef
                              // <<std::endl;
                              assert(cache.count(labelsp));
                              cache[labels]
                              +=coef3/coef1*cache[labelsp]
                                *u3::U(u3::SU3(2,0),n1p.SU3(),n3.SU3(),n2.SU3(),n1.SU3(),1,rho,n3p.SU3(),rhop,1);
                              // std::cout<<"    "<<cache[labels]<<std::endl;
                              // coef2+=cache[labelsp]
                              // *u3::U(u3::SU3(2,0),n1p.SU3(),n3.SU3(),n2.SU3(),n1.SU3(),1,rho,n3p.SU3(),rhop,1);
                            }
                        }
                    }
                }
            }
        }
    }
   

	}
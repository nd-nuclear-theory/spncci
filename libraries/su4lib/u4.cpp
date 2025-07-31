/****************************************************************
  u4.cpp

  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT
****************************************************************/

#include "su4lib/u4.h"

namespace u4{

  
  unsigned int OuterMultiplicity(const u4::U4& lmbda, const u4::U4& mu, const u4::U4& nu)
  {
    const auto& [lmbda1,lmbda2,lmbda3,lmbda4] = lmbda.Key();  
    const auto& [mu1,mu2,mu3,mu4] = mu.Key();  
    const auto& [nu1,nu2,nu3,nu4] = nu.Key();  

    // Converting to labels used in jmp-39-1998-5631-Pan
    // SU(4) labels for lambda
    int l1=lmbda1-lmbda4;
    int l2=lmbda2-lmbda4;
    int l3=lmbda3-lmbda4;

    // SU(4) labels for mu
    int u1=mu1-mu4;
    int u2=mu2-mu4;
    int u3=mu3-mu4;

    int v1=nu1-nu4;
    int v2=nu2-nu4;
    int v3=nu3-nu4;
    
    // in multiplicity formula nu is not just SU(4) labels of final
    // SU(4) irrep, but instead product of the two young diagrams
    // which may have columns with 4 boxes. 
    //
    // Ensuring that product has number of boxes equal to sum of boxes
    // in lmbd and mu SU(4) diagrams
    int n4 = (l1+l2+l3+u1+u2+u3-(v1+v2+v3));

    // If such a four rowed young diagram cannot be constructed
    // then the multiplicity is necessarily zero. 
    if (n4%4!=0 || n4<0)
      return 0;
    
    int v4=int(n4/4);

    v1=v1+v4;
    v2=v2+v4;
    v3=v3+v4;
    
    unsigned int mult=0;
    int eta1_min = std::max(0,v2-l1);
    int eta1_max = std::min(v2-l2,v1-l1);
    // std::cout<<"eta1 "<< eta1_min<<" "<<eta1_max<<std::endl;
    for(int eta1=eta1_min; eta1<=eta1_max; ++eta1)
    {
      int eta2_min = std::max(v3-v2+eta1,0);
      int eta2_max = std::min(eta1,u3);
      eta2_max = std::min(eta2_max,v3-v4);
      // std::cout<<"eta2 "<< eta2_min<<" "<<eta2_max<<std::endl;
      
      for(int eta2=eta2_min; eta2<=eta2_max; ++eta2)
      {
        int eta3_min = std::max(2*eta2-eta1+u2+v4-v3-u3,eta2-eta1+l3+u2-v3);
        eta3_min = std::max(eta3_min,eta2+v4-l3-u3);
        eta3_min = std::max(eta3_min,eta2+l1+l2+l3+2*u2-v1-v2-v3);

        eta3_min = std::max(eta3_min,eta1+l1+l2+u2-v1-v2);
        // std::cout<<eta3_min<<std::endl;
        eta3_min = std::max(eta3_min,0);
        eta3_min = std::max(eta3_min,int((eta2+l1+l2+l3+2*u2-v1-v2-v3)/2));
        
        int eta3_max = std::min(u2-eta1,v4-u3+eta2);
        eta3_max = std::min(eta3_max,u2-u3);
        eta3_max = std::min(eta3_max,l2-v3+u2-eta1+eta2);
        eta3_max = std::min(eta3_max,u2-eta2);
        // std::cout<<"eta3 "<< eta3_min<<" "<<eta3_max<<std::endl;
        for(int eta3=eta3_min; eta3<=eta3_max; ++eta3)
          ++mult;
      }
    }
    return mult;
  }

  MultiplicityTagged<u4::U4>::vector KroneckerProduct(const u4::U4& f, const u4::U4& g)
  {
    MultiplicityTagged<u4::U4>::vector kron_prod;
    const auto& [g1,g2,g3,g4] = g.Key();
    const auto& [f1,f2,f3,f4] = f.Key();
    auto g_tot = g1+g2+g3+g4;
    
    for(unsigned int h1=0u; h1<=g1; ++h1)
    {
      
      // Max is either dictated by number of rows in first column or maximum number of 
      // remaining elements from the first row in g plus the min of elements allowed from the second
      // of g, which is at most the number of elements from the first row of g in the first row of product tableau. 
      auto h2_max = std::min(f1+h1-f2,std::min(h1,g2)+(g1-h1));
      for(auto h2=0u; h2<=h2_max; ++h2)
      {
        auto h3_max = std::min(f2+h2-f3,g_tot-h1-h2);
        for(auto h3=0u; h3<=h3_max; ++h3)
        {
          unsigned int h4 = g_tot-h1-h2-h3;
          if(h4+f4 <=(f3+h3))
          {
            u4::U4 h(f1+h1,f2+h2,f3+h3,f4+h4);

            unsigned int mult = OuterMultiplicity(f, g, h);
            
            if(mult>0)
              kron_prod.emplace_back(u4::U4(f1+h1,f2+h2,f3+h3,f4+h4),mult);
          }
        }
      }
    }
  return kron_prod;

  }

}
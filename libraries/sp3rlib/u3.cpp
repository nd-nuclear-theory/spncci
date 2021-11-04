/****************************************************************
  u3.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/

#include "sp3rlib/u3.h"

#include <sstream>



namespace u3 
{
  std::string SU3::Str() const
  {
    return fmt::format("({},{})",lambda(),mu());
  }

  std::string U3::Str() const
  {
    // return fmt::format("{}({},{})",N(),SU3().lambda(),SU3().mu());
    std::string string = fmt::format("{}{}",N(),SU3());
    return string;
  }

  std::string U3S::Str() const
  {
    return fmt::format("{}x{}",U3(),S());
  }

  std::string U3ST::Str() const
  {
    return fmt::format("{}x{}x{}",U3(),S(),T());
  }

  unsigned int OuterMultiplicity(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x3)
  {
    unsigned int multiplicity = 0;
    // int Nx = (x1.lambda()-x1.mu()) + (x2.lambda()-x2.mu()) - (x3.lambda()-x3.mu());
    int Nx = int(x1.lambda_u()+x2.lambda_u()+x3.mu_u()) - int(x1.mu_u()+x2.mu_u()+x3.lambda_u());

    // short circuit
    if ((Nx%3)!=0)
      return 0;


    // handle conjugation for "oblate" cases
    int Mx;
    u3::SU3 y1, y2, y3;
    if (Nx>=0)
      {
      	y1 = x1; 
      	y2 = x2; 
      	y3 = x3;
      	Mx = Nx/3;
      }
    else
      {
      	y1 = Conjugate(x1);
      	y2 = Conjugate(x2);
      	y3 = Conjugate(x3);
      	Mx = -Nx/3;
      }
    
    // main calculation
    Nx = Mx+y1.mu()+y2.mu()-y3.mu();
    int Mu = std::min(y1.lambda()-Mx,y2.mu());
    if (Mu>=0)
      {
      	int Nu = std::min(y2.lambda()-Mx,y1.mu());
      	if (Nu>=0)
      	  multiplicity = std::max(std::min(Nx,Nu)-std::max(Nx-Mu,0)+1,0);
      }

    return multiplicity;
  }

  MultiplicityTagged<u3::SU3>::vector KroneckerProduct(const u3::SU3& x1, const u3::SU3& x2)
  {
    // calculate bounds on (lambda3,mu3)
    int lambda3_min = 0;  // could be further constrained
    int mu3_min = 0;  // could be further constrained
    int lambda3_max = x1.lambda()+x2.lambda()+std::min(x2.mu(),x1.lambda()+x1.mu());
    int mu3_max = x1.mu()+x2.mu()+std::min(x1.lambda(),x2.lambda());

    // allocate container for product
    MultiplicityTagged<u3::SU3>::vector product;
    int max_entries = (lambda3_max - lambda3_min + 1) * (mu3_max - mu3_min + 1);
    product.reserve(max_entries);
    // generate product
    for (unsigned int lambda3 = lambda3_min; lambda3 <= lambda3_max; ++lambda3)
      for (unsigned int mu3 = mu3_min; mu3 <= mu3_max; ++mu3)
      	{
      	  u3::SU3 x3(lambda3,mu3);
      	  unsigned int multiplicity = OuterMultiplicity(x1,x2,x3);
      	  if (multiplicity>0)
            {
       	      product.push_back(MultiplicityTagged<u3::SU3>(x3,multiplicity));
            }
      	}

    return product;
  }

  MultiplicityTagged<u3::U3>::vector KroneckerProduct(const u3::U3& omega1, const u3::U3& omega2)
  {
    // couple U(1)
    HalfInt N = omega1.N() + omega2.N();
    
    // couple SU(3)
    MultiplicityTagged<u3::SU3>::vector su3_product = KroneckerProduct(omega1.SU3(),omega2.SU3());

    // augment SU(3) entries with U(1) number
    MultiplicityTagged<u3::U3>::vector u3_product;
    u3_product.reserve(su3_product.size());
    for (const auto& x_tagged : su3_product)
      {
        if (u3::U3::ValidLabels(N,x_tagged.irrep))
          {
            u3_product.push_back({{N,x_tagged.irrep},x_tagged.tag});
          }
      }

    return u3_product;
  }

  unsigned int BranchingMultiplicitySO3(const u3::SU3& x, int L)
  {
    // std:: max cannot compare int and unsigned int
    unsigned int multiplicity 
      = std::max(0,(x.lambda()+x.mu()+2-L)/2)
      - std::max(0,(x.lambda()+1-L)/2)
      - std::max(0,(x.mu()+1-L)/2);
    return multiplicity;
  }

  MultiplicityTagged<unsigned int>::vector BranchingSO3(const u3::SU3& x)
  {

    // calculate bound on L
    int L_min = std::min(x.lambda(),x.mu())%2;
    int L_max = x.lambda()+x.mu();

    // allocate container for product
    MultiplicityTagged<unsigned int>::vector branching;
    int max_entries = L_max-L_min+1;
    branching.reserve(max_entries);

    // generate branching
    for (unsigned int L=L_min; L<=L_max; ++L)
      {
      	unsigned int multiplicity = BranchingMultiplicitySO3(x,L);
        if (multiplicity>0)
      	  branching.push_back(MultiplicityTagged<unsigned int>(L,multiplicity));
      }

    return branching;
  }

  MultiplicityTagged<unsigned int>::vector BranchingSO3Constrained(const u3::SU3& x, const HalfInt::pair& r)
  {

    // calculate bound on L
    int L_min = std::max(std::min(x.lambda(),x.mu())%2,int(r.first));
    int L_max = std::min(x.lambda()+x.mu(),int(r.second));

    // allocate container for product
    MultiplicityTagged<unsigned int>::vector branching;
    int max_entries = std::max(L_max-L_min+1,0);  // expression might be negative so impose floor of 0
    branching.reserve(max_entries);

    // generate branching
    for (unsigned int L=L_min; L<=L_max; ++L)
      {
      	unsigned int multiplicity = BranchingMultiplicitySO3(x,L);
      	if (multiplicity>0)
      	  branching.push_back(MultiplicityTagged<unsigned int>(L,multiplicity));
      }

    return branching;
  }


}  // namespace

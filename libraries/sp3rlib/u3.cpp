/****************************************************************
  u3.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "sp3rlib/u3.h"

#include <sstream>


namespace u3 
{
  
  std::string SU3::Str() const
  {
    std::ostringstream ss;

    ss << "(" << lambda << "," << mu << ")";
    return ss.str();
  }

  std::string U3::Str() const
  {
    std::ostringstream ss;

    ss << N() << SU3().Str();

    // ss << "[" << f1 << "," << f2 << "," << f3 << "]";

    return ss.str();
  }

  std::string U3S::Str() const
  {
    std::ostringstream ss;

    ss << w.Str() << "x" << S;

    return ss.str();
  }

  std::string U3ST::Str() const
  {
    std::ostringstream ss;

    ss << w.Str() << "x" << S << "x" << T;

    return ss.str();
  }

  int OuterMultiplicity(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x3)
  {
    int multiplicity = 0;
    int Nx = (x1.lambda-x1.mu) + (x2.lambda-x2.mu) - (x3.lambda-x3.mu);

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
    Nx = Mx+y1.mu+y2.mu-y3.mu;
    int Mu = std::min(y1.lambda-Mx,y2.mu);
    if (Mu>=0)
      {
	int Nu = std::min(y2.lambda-Mx,y1.mu);
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
    int lambda3_max = x1.lambda+x2.lambda+std::min(x2.mu,x1.lambda+x1.mu);  // asymmetric expression, may be further constrained
    int mu3_max = x1.mu+x2.mu+std::min(x1.lambda,x2.lambda);

    // allocate container for product

    // Is this max_entries too big?
    //
    // SU3::Couple just reserves x1.lambda+x2.lambda+x1.mu+x2.mu.
    //
    // TESTING: compare final size() with final capacity()
    //
    // Note: could trim container for product with shrink_to_fit() at end if
    // storage will be persistent enough to make this worthwhile
    // (nonbinding request in C++11)

    MultiplicityTagged<u3::SU3>::vector product;
    int max_entries = (lambda3_max - lambda3_min + 1) * (mu3_max - mu3_min + 1);
    product.reserve(max_entries);

    // generate product
    for (int lambda3 = lambda3_min; lambda3 <= lambda3_max; ++lambda3)
      for (int mu3 = mu3_min; mu3 <= mu3_max; ++mu3)
	{
	  u3::SU3 x3(lambda3,mu3);
	  int multiplicity = OuterMultiplicity(x1,x2,x3);
	  if (multiplicity>0)
 	    product.push_back(MultiplicityTagged<u3::SU3>(x3,multiplicity));
	}

    return product;
  }

  MultiplicityTagged<u3::U3>::vector KroneckerProduct(const u3::U3& w1, const u3::U3& w2)
  {
    // couple U(1)
    HalfInt N = w1.N() + w2.N();
    
    // couple SU(3)
    MultiplicityTagged<u3::SU3>::vector su3_product = KroneckerProduct(w1.SU3(),w2.SU3());

    // augment SU(3) entries with U(1) number
    MultiplicityTagged<u3::U3>::vector u3_product;
    u3_product.reserve(su3_product.size());
    for (auto x_tagged_iter = su3_product.begin(); x_tagged_iter !=su3_product.end(); ++x_tagged_iter)
      {
	MultiplicityTagged<u3::SU3> x_tagged = *x_tagged_iter;
	MultiplicityTagged<u3::U3> w_tagged(u3::U3(N,x_tagged.irrep),x_tagged.tag);
	u3_product.push_back(w_tagged);
      }

    return u3_product;
  }

  int BranchingMultiplicitySO3(const u3::SU3& x, int L)
  {
    int multiplicity 
      = std::max(0,(x.lambda+x.mu+2-L)/2)
      - std::max(0,(x.lambda+1-L)/2)
      - std::max(0,(x.mu+1-L)/2);
    return multiplicity;
  }

  MultiplicityTagged<int>::vector BranchingSO3(const u3::SU3& x)
  {

    // calculate bound on L
    int L_min = std::min(x.lambda,x.mu)%2;
    int L_max = x.lambda+x.mu;

    // allocate container for product
    MultiplicityTagged<int>::vector branching;
    int max_entries = L_max-L_min+1;
    branching.reserve(max_entries);

    // generate branching
    for (int L=L_min; L<=L_max; ++L)
      {
	int multiplicity = BranchingMultiplicitySO3(x,L);
	if (multiplicity>0)
	  branching.push_back(MultiplicityTagged<int>(L,multiplicity));
      }

    return branching;
  }

  MultiplicityTagged<int>::vector BranchingSO3Constrained(const u3::SU3& x, const HalfInt::pair& r)
  {

    // calculate bound on L
    int L_min = std::max(std::min(x.lambda,x.mu)%2,int(r.first));
    int L_max = std::min(x.lambda+x.mu,int(r.second));

    // allocate container for product
    MultiplicityTagged<int>::vector branching;
    int max_entries = L_max-L_min+1;
    branching.reserve(max_entries);

    // generate branching
    for (int L=L_min; L<=L_max; ++L)
      {
	int multiplicity = BranchingMultiplicitySO3(x,L);
	if (multiplicity>0)
	  branching.push_back(MultiplicityTagged<int>(L,multiplicity));
      }

    return branching;
  }

}  // namespace

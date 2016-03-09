/****************************************************************
  sp3r.h

  Sp(3,R) labeling and branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/9/16 (aem,mac): Created based on prototype spstates.py, sp3r.py,
    and coefficients.py.

****************************************************************/

#ifndef SP3R_H_
#define SP3R_H_

#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "am/halfint.h"
#include "sp3rlib/u3.h"

namespace sp3r 
{

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  struct Sp3R
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline SU3() 
      : lambda(0), mu(0) {}
    
    // construction from (lambda,mu)
    inline SU3(int lambda_, int mu_) 
      : lambda(lambda_), mu(mu_) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline std::pair<int,int> Key() const
    {
      return std::pair<int,int>(lambda,mu);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    // Elliott labels
    int lambda, mu;

  };

  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const SU3& x1, const SU3& x2)
  {
    return x1.Key() == x2.Key();
  }

  inline bool operator < (const SU3& x1, const SU3& x2)
  {
    return x1.Key() < x2.Key();
  }

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::SU3& x)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return (x.lambda+1)*(x.mu+1)*(x.lambda+x.mu+2)/2;
  }

  inline u3::SU3 Conjugate(const u3::SU3& x)
  // Conjugate irrep.
  {
    return u3::SU3(x.mu,x.lambda);
  }

  inline int ConjugationGrade(const u3::SU3& x)
  // Integer contribution to phase on conjugation.
  {
    return x.mu + x.lambda;
  }


}  // namespace

#endif

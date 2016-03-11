/****************************************************************
  sp_basis.h

  Sp(3,R) basis construction, indexing, and branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created.

****************************************************************/

#ifndef SP_BASIS_H_
#define SP_BASIS_H_

#include <vector>

#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // Sp(3,R) LGI enumeration
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // LGI input tabulation format:
  //
  //   Nex 2Sp 2Sn 2S lambda mu count
  //   ...
  //
  // Stored for each LGI:
  //   Nex Sp Sn S sigma
  //
  // Note that Sp and Sn are spectators for most conceivable
  // calculations but are retained for informational value.
  //
  // Each input table line results in multiple stored LGIs based on
  // the given count.
  //
  // In calculating sigma, we also need Nsigma_0 for this nucleus, since
  //
  //   Nsigma = Nsigma_0 + Nex
  
  // The Nn_max for each sigma might be determined by any number of
  // creative truncation schemes, but most simply determined based on
  // an Nmax cutoff:
  //
  //   Nn_max = Nmax - Nex
  //          = Nmax - (Nsigma - Nsigma_0)

  struct LGI
  {
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    LGI(int Nex_, const HalfInt& Sp_, const HalfInt& Sn_, const HalfInt& S_, const u3::U3& sigma_)
    : Nex(Nex_), Sp(Sp_), Sn(Sn_), S(S_), sigma(sigma_) {}


    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////
    
    int Nex;
    HalfInt Sp, Sn, S;
    u3::U3 sigma;
  };

  void GenerateLGIVector(std::vector<LGI>& basis, const std::string& lgi_filename, const HalfInt& Nsigma_0);
  // Generates vector of LGIs based on LGI input tabulation.
  //
  // Arguments:
  //   basis (vector<LGI>) : container for LGI list
  //   filename (string) : filename for LGI table file
  //   Nsigma_0 (HalfInt) : Nsigma_0 for nucleus

  void GenerateSp3RSpaces(std::map<u3::U3,sp3r::Sp3RSpace>& spaces, const std::string& lgi_filename, const HalfInt& Nsigma_0);
  // Generates vector of LGIs based on LGI input tabulation.
  //


}  // namespace

#endif

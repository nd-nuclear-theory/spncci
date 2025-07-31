/****************************************************************
  sp3rls.h

  Sp(3,R)xSU(2)>SU(2) labeling and branching restricted to definite L and S.

  SPDX-License-Identifier: MIT

  10/14/23 (cvc): Created.
****************************************************************/

#ifndef SP3RLS_H_
#define SP3RLS_H_

#include <map>
#include <string>

#include "basis/basis.h"
#include "basis/degenerate.h"
#include "basis/operator.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "sp3rlib/sp3r.h"

namespace sp3r
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Class for carrier space of Sp(3,R)>U(3)>SO(3)xSU(2)>SU(2) irrep restricted to L and S
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class SP3RLSSpace;

  class SP3RLSSPace
      : public basis::BaseSpace<SP3RLSSPace, sp3r::U3Subspace, std::tuple<u3::U3>>
  {
  public:
    SP3RLSSPace() = default;
    SP3RLSSPace(
        const u3::U3& sigma,
        unsigned int Nn_max,
        unsigned int s_const,
        unsigned int j_const,
        const bool cache_Kmatrices = true,
        // const bool subspace_labels_only = false,
        const bool branch_to_so3 = true
      );

    // accessors
    const u3::U3& sigma() const { return std::get<0>(labels()); }

    unsigned int Nn_max() const { return Nn_max_; }

    unsigned int l_const() const { return l_const_; }

    // diagnostic output
    std::string DebugStr() const;

  private:
    unsigned int l_const_;
    unsigned int Nn_max_;
  };


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sectors: U3Subspaces connected by operator w0
class SP3RLSSectors
    : public basis::BaseSectors<SP3RLSSPace>
{
 public:
  // Default constructor
  SP3RLSSectors() = default;

  // Constructor
  SP3RLSSectors(
      const SP3RLSSPace& space, const u3::U3& omega0, const bool& su3_generator = false
    );

 private:
  u3::U3 omega0_;
};



}

#endif

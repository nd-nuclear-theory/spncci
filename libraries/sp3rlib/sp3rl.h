/****************************************************************
  sp3rl.h

  Sp(3,R) labeling and branching restricted to good l.

  SPDX-License-Identifier: MIT

  10/14/23 (cvc): Created.
****************************************************************/

#ifndef SP3RL_H_
#define SP3RL_H_

#include <map>
#include <string>

#include "basis/basis.h"
#include "basis/degenerate.h"
#include "basis/hypersector.h"
#include "basis/operator.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "sp3rlib/sp3r.h"
#include "am/halfint.h"
#include "am/am.h"

#include <functional>
#include <fstream>
#include <iterator>

namespace sp3r
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Class for carrier space of Sp(3,R)>U(3)>SO(3) irrep restricted to single l
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class Sp3RLSpace;

  class Sp3RLSpace
      : public basis::BaseSpace<Sp3RLSpace, sp3r::U3Subspace, std::tuple<u3::U3>>
  {
  public:
    Sp3RLSpace() = default;
    Sp3RLSpace(
        const u3::U3& sigma,
        unsigned int Nn_max,
        unsigned int l_const, // currently passing definite l, can try more complex selection rule later
        const bool cache_Kmatrices = true,
        // const bool subspace_labels_only = false,
        const bool branch_to_so3 = true
      );

    // accessors
    const u3::U3& sigma() const { return std::get<0>(labels()); }

    unsigned int Nn_max() const { return Nn_max_; }

    unsigned int L() const { return l_const_; }

    // diagnostic output
    std::string DebugStr() const;

  private:
    unsigned int l_const_;
    unsigned int Nn_max_;
  };


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sectors: U3Subspaces connected by operator w0
class Sp3RLSectors
    : public basis::BaseSectors<Sp3RLSpace>
{
 public:
  // Default constructor
  Sp3RLSectors() = default;

  // Constructor
  Sp3RLSectors(
      const Sp3RLSpace& space, const u3::U3& omega0, const bool& su3_generator = false
    );

 private:
  u3::U3 omega0_;
};



}

#endif

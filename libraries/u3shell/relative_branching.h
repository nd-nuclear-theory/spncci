/****************************************************************
  relative_branching.h

  Branching U(3)xSU(2)xSU(2) relative rme's
                                  
  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  11/2/17 (aem): Created.
****************************************************************/
#ifndef RELATIVE_BRANCHING_H_
#define RELATIVE_BRANCHING_H_

#include "u3shell/tensor_labels.h"
#include "u3shell/upcoupling.h"

namespace u3shell
{
  typedef std::tuple<int,int,HalfInt,HalfInt> RelativeStateLabelsNLST;
  typedef std::tuple<int,HalfInt,HalfInt, RelativeStateLabelsNLST, RelativeStateLabelsNLST> RelativeStateSectorNLST;

  void BranchRelativeNLST(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    std::unordered_map<u3shell::RelativeStateSectorNLST,double,boost::hash<u3shell::RelativeStateSectorNLST>>& rmes_nlst
    );
  // Branching relative matrix elements from U3ST to NLST.  LS-branched rmes stored in unordered map

   void BranchRelativeNLSJT(
    const basis::OperatorLabelsJT& operator_labels, 
    int Nmax, int Jmax,
    const basis::RelativeSpaceLSJT& relative_space_lsjt,
    const std::unordered_map<u3shell::RelativeStateSectorNLST,double,boost::hash<u3shell::RelativeStateSectorNLST>>& rmes_nlst,
    std::array<basis::RelativeSectorsLSJT,3>& isospin_component_sectors_lsjt,
    std::array<basis::OperatorBlocks<double>,3>& isospin_component_blocks_lsjt
    );
   // Branching relative matrix elements from NLST to NLSJT.  J-branched rmes stored in block container

	void BranchRelativeRMEs(const basis::OperatorLabelsJT& operator_labels,int Nmax, int Jmax, 
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    const basis::RelativeSpaceLSJT& relative_space_lsjt,
    std::array<basis::RelativeSectorsLSJT,3>& isospin_component_sectors_lsjt,
    std::array<basis::OperatorBlocks<double>,3>& isospin_component_blocks_lsjt
    );
	// Branch relative rmes from U3ST to NLSJT


}
#endif
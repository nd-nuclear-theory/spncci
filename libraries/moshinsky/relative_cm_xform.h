/****************************************************************
  relative_cm_xform.h

  Transforming from relative to CM.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  11/22/16 (aem,mac): Created.
****************************************************************/

#ifndef RELATIVE_CM_XFORM_H_
#define RELATIVE_CM_XFORM_H_

#include <boost/functional/hash_fwd.hpp>
#include <unordered_map>

#include "am/am.h"
#include "sp3rlib/u3.h"
#include "u3shell/relative_operator.h"
#include "u3shell/tensor_labels.h"

namespace u3shell
{
  //// Defining type names for branched relative and realtive-cm states and tensors
	typedef std::tuple<int,int,HalfInt,HalfInt> RelativeStateLabelsNLST;
	typedef std::tuple<int, HalfInt,HalfInt,RelativeStateLabelsNLST,RelativeStateLabelsNLST> RelativeBraketNLST;
	typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> RelativeCMStateLabelsNLST;	
	typedef std::tuple<int,HalfInt,HalfInt,RelativeCMStateLabelsNLST,RelativeCMStateLabelsNLST> RelativeCMBraketNLST;
	// typedef std::map<u3shell::RelativeCMUnitTensorLabelsU3ST,double> RelativeCMUnitTensorExpansion;
  typedef std::unordered_map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorCache,
                  boost::hash<u3shell::RelativeUnitTensorLabelsU3ST>
                  >RelativeCMExpansion;

	void
	BranchNLST(
	  const u3shell::RelativeUnitTensorLabelsU3ST& relative_unit_tensor,
    u3::WCoefCache& w_cache,
	  std::map<RelativeBraketNLST,double>& branched_rel_unit_tensors
	  );

  inline void
  BranchNLST(
    const u3shell::RelativeUnitTensorLabelsU3ST& relative_unit_tensor,
    std::map<RelativeBraketNLST,double>& branched_rel_unit_tensors
    )
  {
    u3::WCoefCache w_cache;
    BranchNLST(relative_unit_tensor,w_cache,branched_rel_unit_tensors);
  }

  void RelativeToCMLST(
    int Nmax, 
    const std::map<RelativeBraketNLST,double>& branched_rel_unit_tensors,
    std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map);

  void UpcoupleCMU3ST(
    std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map,
    u3::WCoefCache& w_cache,
    RelativeCMUnitTensorCache& rel_cm_u3st_map
    );

  inline void UpcoupleCMU3ST(
    std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map,
    RelativeCMUnitTensorCache& rel_cm_u3st_map
    )
  {
    u3::WCoefCache w_cache;
    UpcoupleCMU3ST(rel_cm_lst_map, w_cache,rel_cm_u3st_map);
  }

  // Conversion from Relative to relative-cm at the U(3) level
  void RelativeToCMU3ST(int Nmax,  
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  std::vector<RelativeCMUnitTensorCache>& unit_tensor_rel_cm_expansions
  );

  //Conversion from relative to relative-cm by branching to LST, converting
  // to relative-cm at lst level and then upcoupling to U3ST
  void RelativeUnitTensorToRelativeCMUnitTensorU3ST(
    int Nmax,  
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
    u3::WCoefCache& w_cache,
    RelativeCMExpansion& unit_relative_cm_map
    );

  // Wrapper function in cases where cache not provided.
  inline void RelativeUnitTensorToRelativeCMUnitTensorU3ST(
    int Nmax,  
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
    RelativeCMExpansion& unit_relative_cm_map
    )
  { 
    u3::WCoefCache w_cache;
    RelativeUnitTensorToRelativeCMUnitTensorU3ST(
      Nmax,relative_unit_tensors, 
      w_cache,unit_relative_cm_map
    );
}


}
#endif
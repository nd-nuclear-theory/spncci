/****************************************************************
  import_jisp16.h
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/14/16 (aem,mac): Created to import relative jisp16 files.
  2/2/17 (mac): Fix missing include guard.
  Deprecated
****************************************************************/

#ifndef IMPORT_INTERACTION_H_
#define IMPORT_INTERACTION_H_

#include <iostream>
#include "eigen3/Eigen/Eigen"  
#include "basis/lsjt_scheme.h"

namespace u3shell
{
	std::vector<Eigen::MatrixXd>
	ImportInteraction(
		std::string filename, 
		const basis::RelativeSpaceLSJT& space,
    	const basis::RelativeSectorsLSJT& sectors,
    	std::string interaction_type
    	);
}

#endif

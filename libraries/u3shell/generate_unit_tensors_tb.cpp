/****************************************************************
  two_body_generator.cpp

  Construct files in Tomas's recoupler format for complete set of
  two-body unit tensors, in two-body truncated space.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/25/16 (mac): Created.

****************************************************************/

#include <iostream>

#include "cppformat/format.h"

#include "u3shell/indexing_u3st.h"
#include "u3shell/two_body_operator.h"

////////////////////////////////////////////////////////////////
// two-body unit tensor enumeration
////////////////////////////////////////////////////////////////

std::map<int,std::vector<TwoBodyUnitTensorLabelsU3ST>>
  GenerateTwoBodyUnitTensorLabelsU3ST(
      const u3shell::TwoBodySectorsU3ST two_body_sectors
    )
  // Generate labels for U3ST-scheme two-body unit tensors acting
  // within the given two-body sectors.
  //
  // The resulting unit tensor labels are grouped by N0, i.e., the
  // number of oscillator quanta caried by the operator.  This N0 will
  // vary from -2*Nmax to +2*Nmax.
  //
  // Arguments:
  //   Nmax (int) : maximum oscillator truncation
  //
  // Returns:
  //   (std::map<int,std::vector<TwoBodyUnitTensorLabelsU3ST>>)
  //   : map from N0 -> vector of relative unit tensor labels
{

  // iterate over sectors
  for (int sector_index = 0; sector_index < two_body_sectors.size(); ++sector_index)
    {

      // extract sector subspaces
      const u3shell::TwoBodySectorsU3ST::SectorType& sector = two_body_sectors.GetSector(sector_index);
      const u3shell::TwoBodySectorsU3ST::SubspaceType& bra_subspace = sector.bra_subspace();
      const u3shell::TwoBodySectorsU3ST::SubspaceType& ket_subspace = sector.ket_subspace();

      //

    }

}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // initialize su3lib
  u3::U3CoefInit();

  // define relevant two-body space
  const int Nmax = 4;
  u3shell::TwoBodySpaceU3ST two_body_space(Nmax);

  // generate list of unit tensors
  u3shell::TwoBodySectorsU3ST two_body_sectors(two_body_space);
  std::vector<std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST>> two_body_unit_tensor_labels;
  GenerateTwoBodyUnitTensorLabelsU3ST(
      two_body_sectors,
      two_body_unit_tensor_labels
    );
  //       int Nmax, 
  //       std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels
  //       );


  // write list of unit tensors

  // generate unit tensors

  // termination
  return 0;
}

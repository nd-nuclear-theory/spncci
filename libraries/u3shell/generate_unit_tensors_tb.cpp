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
      ////////////////////////////////
      // extract labels
      ////////////////////////////////

      // extract sector subspaces
      const u3shell::TwoBodySectorsU3ST::SectorType& sector = two_body_sectors.GetSector(sector_index);
      const u3shell::TwoBodySectorsU3ST::SubspaceType& bra_subspace = sector.bra_subspace();
      const u3shell::TwoBodySectorsU3ST::SubspaceType& ket_subspace = sector.ket_subspace();

      // extract bra subspace labels
      u3::U3 omagap;
      HalfInt Sp, Tp;
      int gp;
      std::tie(omegap,Sp,Tp,gp) =  bra_subspace.GetSubspaceLabels();

      // extract ket subspace labels
      u3::U3 omaga;
      HalfInt S, T;
      int g;
      std::tie(omega,S,T,g) =  ket_subspace.GetSubspaceLabels();

      ////////////////////////////////
      // determine SU(3)xSxT couplings
      ////////////////////////////////

      // U(3): operator "destroys" x and "creates" xp

      // SU(3) couplings
      u3::SU3 xp = omegap.SU3();
      u3::SU3 x = omega.SU3();
      MultiplicityTagged<u3::SU3>::vector x0_set
        =u3::KroneckerProduct(Conjugate(x),xp);

      // U(1) coupling
      int etap = omegap.N();
      int eta = omega.N();
      int N0 = etap - eta;

      // SxT couplings
      HalfInt::pair S0_range = am::ProductAngularMomentumRange(S,Sp);
      HalfInt::pair T0_range = am::ProductAngularMomentumRange(T,Tp);

      // parity coupling
      int g0 = (g+gp)%2;
      
      ////////////////////////////////
      // construct label set
      ////////////////////////////////

      for (MultiplicityTagged<u3::SU3> x0_rho0max : x0_set)
        // for each SU(3)
        {
          // extract SU(3) labels
          u3::SU3 x0 = x0_rho0max.irrep;
          int rho0max = x0_rho0max.multiplicity;
          
          for (HalfInt S0=S0_range.first; S0 <= S0_range.second; ++S0)
            for (HalfInt T0=T0_range.first; T0 <= T0_range.second; ++T0)
              // for each SxT
              {
                
                // collect tensorial labels
                u3shell::OperatorLabelsU3ST operator_labels(N0,x0,S0,T0,g0);

                // iterate over remaining (nontensorial) labels
                for (int bra_index = 0; bra_index < bra_space.size(); ++bra_index)
                  for (int ket_index = 0; ket_index < ket_space.size(); ++ket_index)
                    // for each <bra|ket> pair
                    {
                      // extract state labels


                    }
                    
                    

                
              }
        }

      // construct unit tensors between these subspaces

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

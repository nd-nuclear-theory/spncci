
/****************************************************************
  hyperblocks_u3s.h

  U(3)xS layer of SpNCCI basis branching.
                                  
  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  1/2/19 (aem): Created based on branching_u3s
  6/21/19 (aem): Updated ComputeManyBodyRMEs to use new seed set
    up and recurrence functions. 
****************************************************************/

#ifndef SPNCCI_SPNCCI_HYPERBLOCKS_U3S_H_
#define SPNCCI_SPNCCI_HYPERBLOCKS_U3S_H_

#include "spncci/spncci_basis.h"
#include "spncci/recurrence.h"
#include "spncci/parameters.h"
namespace spncci
{

  typedef std::vector<spncci::ObservableBabySpNCCIHypersectors> ObservableHypersectorsTable;
  typedef std::tuple<int,int,int,int> ObservableHypersectorLabels;
  typedef std::vector< std::vector<basis::OperatorHyperblocks<double>>> ObservableHyperblocksTable;


	void ContractBabySpNCCISymmetricHypersectors(
	  const spncci::LGIPair& lgi_pair,
	  int num_observables, int num_hw_values,
	  const spncci::BabySpNCCISpace& baby_spncci_space,
	  const std::vector<u3shell::ObservableSpaceU3S>& observable_spaces,
	  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
	  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors1,
	  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors2,
	  const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks1,
	  const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks2,
	  const std::vector<std::vector<u3shell::RelativeRMEsU3SSubspaces>>& observables_relative_rmes,
	  spncci::ObservableHypersectorsTable& observable_hypersectors_table,
	  spncci::ObservableHyperblocksTable& observable_hyperblocks_table
  	);

void GetBabySpNCCIHyperBlocks(
  const int observable_index,
  const int hw_index,
  const spncci::LGIPair& lgi_pair,
  std::vector<spncci::ObservableHypersectorLabels>& list_baby_spncci_hypersectors,
  basis::OperatorHyperblocks<double>& baby_spncci_observable_hyperblocks
  );

void GetOperatorTile(
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::ObservableSpaceU3S& observable_space,
  const spncci::SubspaceSpBasis& spbasis_subspace_bra,
  const spncci::SubspaceSpBasis& spbasis_subspace_ket,
  const std::vector<int>& offsets_bra_subspace,
  const std::vector<int>& offsets_ket_subspace,
  const HalfInt& J0, const HalfInt& Jp, const HalfInt& J,
  const int hw_index,
  const int observable_index,
  const spncci::LGIPair& lgi_pair,
  u3::WCoefCache& w_cache,
  const std::vector<spncci::ObservableHypersectorLabels>& list_baby_spncci_hypersectors,
  basis::OperatorHyperblocks<double>& baby_spncci_observable_hyperblocks,
  spncci::OperatorBlock& tile
  );

  void ConstructSymmetricOperatorMatrix(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    const HalfInt& J0,
    const spncci::SpaceSpBasis& spbasis_bra, //For a given J
    const spncci::SpaceSpBasis& spbasis_ket, //For a given J
    const std::vector<spncci::LGIPair>& lgi_pairs,
    int observable_index, int hw_index,
    spncci::OperatorBlock& operator_matrix
  );

void ComputeManyBodyRMEs(
  const spncci::RunParameters& run_parameters,
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::SpNCCISpace& spncci_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const std::vector<u3shell::ObservableSpaceU3S>& observable_spaces,
  const std::vector<std::vector<u3shell::RelativeRMEsU3SSubspaces>>& observables_relative_rmes,
  const spncci::KMatrixCache& k_matrix_cache,
  const spncci::KMatrixCache& kinv_matrix_cache,
  spncci::OperatorBlocks& lgi_transformations,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::LGIPair& lgi_pair
  );
  //Calculates many-body matrix elements of operators defined by observables_relative_rmes.  
  //
  // Calculates unit tensor hyperblocks for given lgi pair and it's conjugate
  // Contracts unit tensor hyperblocks with relative rmes of operators to get operator rme hyperblocks
  // Writes hyperblocks to file

//********************************************** Added by J.H. ******************************************
void ComputeOneBodyUnitTensorRMEs(
  const spncci::RunParameters& run_parameters,
  int N1vp,
  int N1vn,
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::SpNCCISpace& spncci_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::OneBodyUnitTensorSpaceU3S& one_body_unit_tensor_space,
  const spncci::KMatrixCache& k_matrix_cache,
  const spncci::KMatrixCache& kinv_matrix_cache,
  spncci::OperatorBlocks& lgi_transformations,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::LGIPair& lgi_pair,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
  spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersectors,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks2,
  spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersectors2
  );
  // Calculates many-body RMEs of one-body unit tensors.  
  // Calculates unit tensor hyperblocks for given lgi pair and it's conjugate.

void ComputeTwoBodyDensityRMEs(
  const spncci::RunParameters& run_parameters,
  int N1vp,
  int N1vn,
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::SpNCCISpace& spncci_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::TwoBodyDensitySpace& two_body_density_space,
  const spncci::KMatrixCache& k_matrix_cache,
  const spncci::KMatrixCache& kinv_matrix_cache,
  spncci::OperatorBlocks& lgi_transformations,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::LGIPair& lgi_pair
  );
  // Calculates many-body RMEs of two-body densities.  
  // Calculates tbd hyperblocks for given lgi pair and it's conjugate.

class OBUnitTensorRMELabels {
  public:
  // default constructor
  inline OBUnitTensorRMELabels() : N_ex_sigma_p_(0), lambda_sigma_p_(0), mu_sigma_p_(0), twice_Sp_p_(0), twice_Sn_p_(0), twice_S_p_(0), N_ex_omega_p_(0), lambda_omega_p_(0), mu_omega_p_(0), gamma_p_(0), upsilon_p_(0), N_ex_sigma_(0), lambda_sigma_(0), mu_sigma_(0), twice_Sp_(0), twice_Sn_(0), twice_S_(0), N_ex_omega_(0), lambda_omega_(0), mu_omega_(0), gamma_(0), upsilon_(0), lambda0_(0), mu0_(0), twice_S0_(0), eta_p_(0), eta_(0), rho0_(0), Tz_(0) {}
  // constructor
  inline OBUnitTensorRMELabels(int N_ex_sigma_p, int lambda_sigma_p, int mu_sigma_p, int twice_Sp_p, int twice_Sn_p, int twice_S_p, int N_ex_omega_p, int lambda_omega_p, int mu_omega_p, int gamma_p, int upsilon_p, int N_ex_sigma, int lambda_sigma, int mu_sigma, int twice_Sp, int twice_Sn, int twice_S, int N_ex_omega, int lambda_omega, int mu_omega, int gamma, int upsilon, int lambda0, int mu0, int twice_S0, int eta_p, int eta, int rho0, int Tz) : N_ex_sigma_p_(N_ex_sigma_p), lambda_sigma_p_(lambda_sigma_p), mu_sigma_p_(mu_sigma_p), twice_Sp_p_(twice_Sp_p), twice_Sn_p_(twice_Sn_p), twice_S_p_(twice_S_p), N_ex_omega_p_(N_ex_omega_p), lambda_omega_p_(lambda_omega_p), mu_omega_p_(mu_omega_p), gamma_p_(gamma_p), upsilon_p_(upsilon_p), N_ex_sigma_(N_ex_sigma), lambda_sigma_(lambda_sigma), mu_sigma_(mu_sigma), twice_Sp_(twice_Sp), twice_Sn_(twice_Sn), twice_S_(twice_S), N_ex_omega_(N_ex_omega), lambda_omega_(lambda_omega), mu_omega_(mu_omega), gamma_(gamma), upsilon_(upsilon), lambda0_(lambda0), mu0_(mu0), twice_S0_(twice_S0), eta_p_(eta_p), eta_(eta), rho0_(rho0), Tz_(Tz) {}
  // accessors
  inline int N_ex_sigma_p() const {
    return N_ex_sigma_p_;
  }
  inline int lambda_sigma_p() const {
    return lambda_sigma_p_;
  }
  inline int mu_sigma_p() const {
    return mu_sigma_p_;
  }
  inline int twice_Sp_p() const {
    return twice_Sp_p_;
  }
  inline int twice_Sn_p() const {
    return twice_Sn_p_;
  }
  inline int twice_S_p() const {
    return twice_S_p_;
  }
  inline int N_ex_omega_p() const {
    return N_ex_omega_p_;
  }
  inline int lambda_omega_p() const {
    return lambda_omega_p_;
  }
  inline int mu_omega_p() const {
    return mu_omega_p_;
  }
  inline int gamma_p() const {
    return gamma_p_;
  }
  inline int upsilon_p() const {
    return upsilon_p_;
  }
  inline int N_ex_sigma() const {
    return N_ex_sigma_;
  }
  inline int lambda_sigma() const {
    return lambda_sigma_;
  }
  inline int mu_sigma() const {
    return mu_sigma_;
  }
  inline int twice_Sp() const {
    return twice_Sp_;
  }
  inline int twice_Sn() const {
    return twice_Sn_;
  }
  inline int twice_S() const {
    return twice_S_;
  }
  inline int N_ex_omega() const {
    return N_ex_omega_;
  }
  inline int lambda_omega() const {
    return lambda_omega_;
  }
  inline int mu_omega() const {
    return mu_omega_;
  }
  inline int gamma() const {
    return gamma_;
  }
  inline int upsilon() const {
    return upsilon_;
  }
  inline int lambda0() const {
    return lambda0_;
  }
  inline int mu0() const {
    return mu0_;
  }
  inline int twice_S0() const {
    return twice_S0_;
  }
  inline int eta_p() const {
    return eta_p_;
  }
  inline int eta() const {
    return eta_;
  }
  inline int rho0() const {
    return rho0_;
  }
  inline int Tz() const {
    return Tz_;
  }
  // key tuple, comparisons
  typedef std::tuple<int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,int> KeyType;
  inline KeyType Key() const{
    return KeyType(N_ex_sigma_p(), lambda_sigma_p(), mu_sigma_p(), twice_Sp_p(), twice_Sn_p(), twice_S_p(), N_ex_omega_p(), lambda_omega_p(), mu_omega_p(), gamma_p(), upsilon_p(), N_ex_sigma(), lambda_sigma(), mu_sigma(), twice_Sp(), twice_Sn(), twice_S(), N_ex_omega(), lambda_omega(), mu_omega(), gamma(), upsilon(), lambda0(), mu0(), twice_S0(), eta_p(), eta(), rho0(), Tz());
  }
  inline friend bool operator == (const OBUnitTensorRMELabels& l1, const OBUnitTensorRMELabels& l2){
    return l1.Key() == l2.Key();
  }
  inline friend bool operator < (const OBUnitTensorRMELabels& l1, const OBUnitTensorRMELabels& l2){
    return l1.Key() < l2.Key();
  }
  private:
  int N_ex_sigma_p_, lambda_sigma_p_, mu_sigma_p_, twice_Sp_p_, twice_Sn_p_, twice_S_p_, N_ex_omega_p_, lambda_omega_p_, mu_omega_p_, gamma_p_, upsilon_p_, N_ex_sigma_, lambda_sigma_, mu_sigma_, twice_Sp_, twice_Sn_, twice_S_, N_ex_omega_, lambda_omega_, mu_omega_, gamma_, upsilon_, lambda0_, mu0_, twice_S0_, eta_p_, eta_, rho0_, Tz_;
};
//***************************************************************************************
}//namespace
#endif

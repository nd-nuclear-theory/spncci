/****************************************************************
  upcoupling.h
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  6/14/16 (aem,mac): Created to import relative jisp16 files.
  2/2/17 (mac):
    + Fix missing include guard.
    + Add class UpcouplingLabels.
  2/14/17 (aem): 
    Change typedef RelativeRMEsU3ST from map to unordered_map.
  2/21/17 (mac): Change interface to ReadRelativeOperatorU3ST.
  3/28/17 (aem): Change upcoupling to allow for non-isoscalar 
    interactions
  9/11/17 (aem): Update documentation and remove deprecated routines
****************************************************************/

#ifndef SPNCCI_U3SHELL_UPCOUPLING_H_
#define SPNCCI_U3SHELL_UPCOUPLING_H_

#include <iomanip>
#include <iostream>
#include <boost/functional/hash_fwd.hpp>
#include "eigen3/Eigen/Eigen" 
#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"
#include "u3shell/relative_operator.h"
#include "u3shell/unit_tensor_space_u3s.h" 
#include "sp3rlib/u3coef.h"

namespace u3shell
{
  typedef std::tuple<int,int,int> RelativeSubspaceLabelsNLST;
  typedef std::tuple<int,int,int,u3shell::RelativeSubspaceLabelsNLST,
            u3shell::RelativeSubspaceLabelsNLST> RelativeSectorNLST;


  class UpcouplingLabels
    : public std::tuple<int,int>
  // Storage of (kappa0,L0) pair.
  {
    public:

    UpcouplingLabels(int kappa0, int L0)
      : std::tuple<int,int>(kappa0,L0)
    {}

    int kappa0() const {return std::get<1>(*this);};
    int L0() const {return std::get<1>(*this);};
  };

  typedef std::unordered_map<
    std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>,
    double,
    boost::hash<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>>
    >
    RelativeRMEsU3ST;

  typedef std::unordered_map<
    std::tuple<int,int,int>,
    std::vector<double>,
    boost::hash<std::tuple<int,int,int>> 
    > 
    RelativeRMEsU3SSubspaces;


  // Mapping to store coefficients (RMEs) of relative operator in
  // terms of relative unit tensors.
  //
  // (relative_unit_tensor_labels,kappa0,L0) -> RME
  //
  // TODO: replace with map (relative_tensor_unit_labels,upcoupling_labels) -> RME

  typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> RelativeCMStateLabelsNLST;  
  typedef std::tuple<int,HalfInt,HalfInt,RelativeCMStateLabelsNLST,RelativeCMStateLabelsNLST> RelativeCMBraketNLST;
  typedef std::unordered_map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorCache,
    boost::hash<u3shell::RelativeUnitTensorLabelsU3ST>
    >RelativeCMExpansion;


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Upcoupling
  //
  // Programs calling these function need to initialize with 
  // U3CoefInit()
  //
  // Upcoupling relative operator from NLSJTg scheme to U3STg scheme
  // Relative matrix element of operators in NLSJTg are stored in 
  // std::array<basis::OperatorBlocks<double>,3> where the three componenets of the
  // array are T0=0,1,2. The vectors of Eigen::MatrixdX matrices
  // for each element of the array correspond to the basis::RelativeSectorsLSJT sectors, 
  // Each sector is labled by L'S'J'T'g' J0 g0 L S J T g for fixed T0
  // 
  // Each matrix is indexed by np,n where 2n+L=N and 2n'+L'=N'
  //
  // SU(3)XSU(2)xSU(2) reduced matrix elements stored in 
  // (Relative Tensor, k0, L0)
  // std::map<std::pair<RelativeUnitTensorLabelsU3ST,int, int>,double> rme_map;
  //
  //Generate list of U3ST operator labels which are given by 
  //	u3shell::GenerateRelativeUnitTensorLabelsU3ST(        
  //         int Nmax, 
  //         std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
  //	       )
  //
  // Iterate over vector
  // To get vector of <L0,k0max>, <L,kmax> and <Lp, Kpmax>
  //   MultiplicityTagged<int>::vector BranchingSO3Constrained(const u3::SU3& x, const HalfInt::pair& r);
  // 
  //
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Upcoupling NLSJT reduced matrix elements to NLST reduced matrix elements and stores them in rme_nlst_map
  //
  // args: rme_nlst_map is a map of the NLST reduced matrix elements
  //       space gives relative space
  //       sectors indicate sectors corresponding to defined space 
  //       sector_vector are the input matrix elements for each sector 
  // 

  void UpcouplingNLST(
      const basis::RelativeSpaceLSJT& space,
      const basis::RelativeSectorsLSJT& sectors,
      const std::vector<Eigen::MatrixXd>& sector_vector, 
      int J0, int g0, int T0, int Nmax,
      std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map
    );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Upcoupling NLST reduced matrix elements to U3ST reduced matrix elements and stores them in rme_map
  //
  // args: rme_nlst_map is a map of the NLST reduced matrix elements
  //       J0, g0,T0 and Nmax (int) are basis parameters 
  // 
  void UpcouplingU3ST(
      std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map,
      int Nmax,
      u3::WCoefCache& w_cache,
      RelativeRMEsU3ST& rme_map
    );

  // If no cache is passed as an arguemnt, creates cache for use in calculations
  void UpcouplingU3ST(
      std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map,
      int Nmax,
      RelativeRMEsU3ST& rme_map
    );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Upcoupling NLSJT reduced matrix elements to U3ST reduced matrix elements and stores them in rme_map
  //
  // 
  void Upcoupling(    
      const basis::RelativeSpaceLSJT& space,
      const std::array<basis::RelativeSectorsLSJT,3>& T0_sector_labels,
      const std::array<basis::OperatorBlocks<double>,3>& T0_sectors,
      u3::WCoefCache& w_cache,
      int J0, int g0, int T0,int Nmax,
      RelativeRMEsU3ST& rme_map
    );

  // If no cache is passed as an arguemnt, creates cache for use in calculations
  void Upcoupling(    
      const basis::RelativeSpaceLSJT& space,
      const std::array<basis::RelativeSectorsLSJT,3>& T0_sector_labels,
      const std::array<basis::OperatorBlocks<double>,3>& T0_sectors,
      int J0, int g0, int T0,int Nmax,
      RelativeRMEsU3ST& rme_map
    );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Upcoupling relative-cm reduced matrix elements to U3ST reduced matrix elements
  //
  void UpcoupleCMU3ST(
      std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map,
      u3::WCoefCache& w_cache,
      RelativeCMUnitTensorCache& rel_cm_u3st_map
    );

  // Overloading function.  Provides W coefficient cache is not provided on input.  
  void UpcoupleCMU3ST(
      std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map,
      RelativeCMUnitTensorCache& rel_cm_u3st_map
    );

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // input and output routines for U3ST decomposed operators 
  //

  void WriteRelativeOperatorU3ST(const std::string& filename, RelativeRMEsU3ST& relative_rmes, bool hermitian);  
  // Writes relative RMEs to filename
  // if hermitian=True, then hermiticity enforced by zeroing out conjugate pairs


  void ReadRelativeOperatorU3ST(const std::string& filename, RelativeRMEsU3ST& relative_rmes);

  void ReadRelativeOperatorU3ST(
    int Nmax, int N1v,
    const std::string& filename, 
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    RelativeRMEsU3SSubspaces& relative_rmes
    );

}
// 

#endif

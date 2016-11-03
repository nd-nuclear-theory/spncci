/****************************************************************
  upcoupling.h
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/14/16 (aem,mac): Created to import relative jisp16 files.
****************************************************************/
#include <iomanip>
#include <iostream>

#include "eigen3/Eigen/Eigen" 
#include "basis/lsjt_scheme.h"
#include "u3shell/relative_operator.h"
namespace u3shell
{
  typedef std::tuple<int,int,int> RelativeSubspaceLabelsNLST;
  typedef std::tuple<int,int,u3shell::RelativeSubspaceLabelsNLST,u3shell::RelativeSubspaceLabelsNLST> RelativeSectorNLST;
  typedef std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>,double>RelativeRMEsU3ST;
// Programs calling these function need to initialize with 
// U3CoefInit()

// Upcoupling relative operator from NLSJTg scheme to U3STg scheme
// Relative matrix element of operators in NLSJTg are stored in a std::vector
// of Eigen::MatrixdX matrices for each L'S'J'T'g' J0 g0 L S J T g sector
// Each matrix is indexed by np,n where 2n+L=N and 2n'+L'=N'

// SU(3)XSU(2)xSU(2) reduced matrix elements stored in 
// (Relative Tensor, k0, L0)
 // std::map<std::pair<RelativeUnitTensorLabelsU3ST,int, int>,double> rme_map;

//Generate list of U3ST operator labels which are given by 
//	u3shell::GenerateRelativeUnitTensorLabelsU3ST(        
//         int Nmax, 
//         std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
//	       )

// Iterate over vector
// To get vector of <L0,k0max>, <L,kmax> and <Lp, Kpmax>
//   MultiplicityTagged<int>::vector BranchingSO3Constrained(const u3::SU3& x, const HalfInt::pair& r);
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
    int T0, int Nmax,
    RelativeRMEsU3ST& rme_map
    );

  void Upcoupling(    
    const basis::RelativeSpaceLSJT& space,
    const basis::RelativeSectorsLSJT& sectors,
    const std::vector<Eigen::MatrixXd>& sector_vector, 
    int J0, int g0, int T0,int Nmax,
    RelativeRMEsU3ST& rme_map
    );

  void WriteRelativeOperatorU3ST(std::ostream& os, const RelativeRMEsU3ST& relative_rmes);  
  void ReadRelativeOperatorU3ST(std::istream& is, RelativeRMEsU3ST& relative_rmes);

}
// 


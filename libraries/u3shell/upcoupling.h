/****************************************************************
  upcoupling.h
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/14/16 (aem,mac): Created to import relative jisp16 files.
****************************************************************/
#include "eigen3/Eigen/Eigen" 
#include "basis/lsjt_scheme.h"
#include "u3shell/relative_operator.h"
namespace u3shell
{
  typedef std::tuple<int,int,int> RelativeSubspaceLabelsNLST;
  typedef std::tuple<int,int,u3shell::RelativeSubspaceLabelsNLST,u3shell::RelativeSubspaceLabelsNLST> RelativeSectorNLST;

// Programs calling these function need to initialize with 
// U3CoefInit()

// Upcoupling relative operator from NLSJTg scheme to U3STg scheme
// Relative matrix element of operators in NLSJTg are stored in a std::vector
// of Eigen::MatrixdX matrices for each L'S'J'T'g' J0 g0 L S J T g sector
// Each matrix is indexed by np,n where 2n+L=N and 2n'+L'=N'

// SU(3)XSU(2)xSU(2) reduced matrix elements stored in 
// (Relative Tensor, k0)
 // std::map<std::pair<RelativeUnitTensorLabelsU3ST,int>,double> rme_map;

//Generate list of U3ST operator labels which are given by 
//	u3shell::GenerateRelativeUnitTensorLabelsU3ST(        
//         int Nmax, 
//         std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
//	       )

// Iterate over vector
// To get vector of <L0,k0max>, <L,kmax> and <Lp, Kpmax>
//   MultiplicityTagged<int>::vector BranchingSO3Constrained(const u3::SU3& x, const HalfInt::pair& r);
// 
  void UpcouplingNLST(
    const basis::RelativeSpaceLSJT& space,
    const basis::RelativeSectorsLSJT& sectors,
    const std::vector<Eigen::MatrixXd>& sector_vector, 
    int J0, int g0, int T0, int Nmax,
    std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map
    );

  void UpcouplingU3ST(
    std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map,
    int J0, int g0, int T0, int Nmax,
    std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int>,double>& rme_map
    );

  void Upcoupling(    
    const basis::RelativeSpaceLSJT& space,
    const basis::RelativeSectorsLSJT& sectors,
    const std::vector<Eigen::MatrixXd>& sector_vector, 
    int J0, int g0, int T0,int Nmax,
    std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int>,double>& rme_map
    );
  
  
}
// 


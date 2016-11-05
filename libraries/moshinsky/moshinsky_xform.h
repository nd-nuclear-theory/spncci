/****************************************************************
  moshinsky.h

  U(3) Moshinsky coefficient.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created based on prototype u3states.py, u3.py, and so3.py.
  3/8/16 (aem,mac): Add U3ST structure and rename U3S structure.
  3/9/16 (aem,mac): Add KeyType typedefs.  Extract MultiplicityTagged.
  5/11/16 (aem,mac): Move to u3shell namespace.
  11/4/16 (aem): Factored out
****************************************************************/

#ifndef MOSHINSKY_H_
#define MOSHINSKY_H_

#include "eigen3/Eigen/Eigen"

// #include "sp3rlib/u3.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/u3st_scheme.h"
#include "moshinsky/relative_cm_xform.h"

namespace u3shell
{

  typedef std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> 
            IndexedRelativeUnitTensorLabelsU3ST;
  typedef std::tuple<u3shell::TwoBodyUnitTensorLabelsU3ST, int, int>
            IndexTwoBodyTensorLabelsU3ST;
  typedef std::map<IndexTwoBodyTensorLabelsU3ST,double>IndexedTwoBodyTensorRMEsU3ST;
  typedef  std::unordered_map<u3shell::RelativeUnitTensorLabelsU3ST, 
                            u3shell::TwoBodyUnitTensorCoefficientsU3ST,
                            boost::hash<u3shell::RelativeUnitTensorLabelsU3ST>
                          > TwoBodyExpansionMap;
  typedef std::map<u3shell::IndexedRelativeUnitTensorLabelsU3ST,double> 
            IndexedRelativeUnitTensorCoefficientsU3ST;

  double 
  MoshinskyCoefficient(int r, int R, int r1, int r2, const u3::SU3& x);
  //Computes and returns the SU(3) moshinsky transformation

  inline double 
  MoshinskyCoefficient(int r, int R, int r1, int r2, const u3::U3& w)
  {
    return MoshinskyCoefficient(r,R,r1,r2,w.SU3());
  }
  // Overloading Moshinsky to take integers and U3 for total symmetry


  Eigen::MatrixXd 
  MoshinskyTransform(
    const u3::SU3& x0, int etap,int eta,
    const u3shell::TwoBodySubspaceU3ST& bra_subspace, 
    const u3shell::TwoBodySubspaceU3ST& ket_subspace, 
    int rho0,
    std::string normalization="NAS"
  );
  // Generate the moshinsky transformation of unit matrix element of relative_cm unit tensor 
  // with SU(3) character x0 and relative oscillator quanta etap and eta. Returns transformed
  // matrix elements as a matrix.
  //
  // The cm oscillator quanta is fixed by conservation of oscillator quanta.  Each two-body
  // subspace has fixed total number of oscillator quanta.  
  //
  // Arguments:
  //  x0 (input) : SU(3) character of unit tensor
  //  etap (input) : number of relative oscillator quanta in bra source  subspace
  //  eta (input) : number of relative ocillator quanta in ket source subspace
  //  bra_subspace (input) : bra target subspace
  //  ket_subspace (input) : ket target subspace
  //  rho0 (input) : indexing outermultiplicity of coupling of SU(3) character of 
  //        ket_subspace and x0 to SU(3) character of bra_subspace 
  //  normalization (optional input) : if "NAS" then target RME's are normalized  
  //        anti-symmetrized, if "AS" then target RME's and anti-symmetrized but 
  //        not normalized.  Default is "NAS"   

  void 
  MoshinskyTransformUnitTensor(
    const u3shell::RelativeUnitTensorLabelsU3ST& tensor, 
    double expansion_coef, 
    u3shell::TwoBodySpaceU3ST& space,
    u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_expansion,
    std::string normalization="NAS"
  );
  // Generates expansion  of a relative unit tensor (tensor) in terms of two-body unit
  // tensors and stores them in two_body_expansion. 
  //
  // The transformation from relative to relative-cm basis is carried out internally at 
  // at U(3) level. 
  //
  // Arguments:
  //  tensor (input) : tensorial properties of relative unit tensor
  //  expansion_coef (input) : coefficient by which the relative unit tensor is multiplied. 
  //                           If considering a single unit tensor, expansion_coef=1.0;
  //  space (input) : source space
  //  two_body_expansion (output) : Container for expansion coefficients and corresponding
  //                                two-body unit tensor labels. 
  //  normalization (optional input) : if "NAS" then target RME's are normalized  
  //        anti-symmetrized, if "AS" then target RME's and anti-symmetrized but 
  //        not normalized.  Default is "NAS"

  void
  MoshinskyTransformUnitTensor(
    const u3shell::RelativeUnitTensorLabelsU3ST& tensor, 
    RelativeCMUnitTensorCache& unit_relative_cm_expansion,
    u3shell::TwoBodySpaceU3ST& space,
    u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_expansion,
    std::string normalization="NAS",
    double expansion_coef=1.0
  );
  // Generates expansion  of a relative unit tensor (tensor) in terms of two-body unit
  // tensors and stores them in two_body_expansion. 
  // 
  //
  // Arguments:
  //  tensor (input) : tensorial properties of relative unit tensor
  //  unit_relative_cm_expansion (input): expansion of relative unit tensor in terms of 
  //       relative-cm unit tensors.  Can be obained using 
  //       u3shell::RelativeUnitTensorToRelativeCMUnitTensorU3ST() in relative_cm_xform.h
  //  expansion_coef (input) : coefficient by which the relative unit tensor is multiplied. 
  //       If considering a single unit tensor, expansion_coef=1.0;
  //  space (input) : source space
  //  two_body_expansion (output) : Container for expansion coefficients and corresponding
  //       two-body unit tensor labels. 
  //  normalization (optional input) : if "NAS" then target RME's are normalized  
  //        anti-symmetrized, if "AS" then target RME's and anti-symmetrized but 
  //        not normalized.  Default is "NAS"



  void TransformRelativeTensorToTwobodyTensor(
    const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_unit_tensor_expansion, 
    u3shell::RelativeCMExpansion& unit_relative_cm_expansion,
    u3shell::TwoBodySpaceU3ST& space,
    u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_unit_tensor_expansion,
    std::string normalization="NAS"
    );


void RelativeUnitTensorToTwobodyU3ST(int Nmax,  
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  u3shell::RelativeCMExpansion& unit_relative_cm_map,
  TwoBodyExpansionMap& two_body_expansion_vector,
  std::string normalization="AS"
  );
  // Moshinsky Transform operator decomposed in terms of unit tensors to 
  // two-body nomralized anti-symmeterized space 
  // 
  // 
  // Arguments:
  //    relative_unit_tensors (input) : vector of unit tensors
  //    unit_relative_cm_map (input) : expansion of unit tensors in vector in terms of 
  //        relative-cm unit tensors
  //    two_body_expansion_vector (output)
  //    normalization (optional input) : if "NAS" then target RME's are normalized  
  //        anti-symmetrized, if "AS" then target RME's and anti-symmetrized but 
  //        not normalized.  Default is "NAS"

  
  

void
ConvertRelativeTensorToTwoBodyTensor(int Nmax,
  std::map<IndexedRelativeUnitTensorLabelsU3ST,double>& relative_rmes,
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  u3shell::RelativeCMExpansion& unit_relative_cm_map,
  IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
  std::string normalization="AS"
  );

void 
ConvertRelativeTensorToTwoBodyTensor(int Nmax,
  std::map<IndexedRelativeUnitTensorLabelsU3ST,double>& relative_rmes,
  IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
  std::string normalization="AS"
  );

} //namespace

#endif 






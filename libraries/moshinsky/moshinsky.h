/****************************************************************
  moshinsky.h

  U(3) Moshinsky coefficient.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created based on prototype u3states.py, u3.py, and so3.py.
  3/8/16 (aem,mac): Add U3ST structure and rename U3S structure.
  3/9/16 (aem,mac): Add KeyType typedefs.  Extract MultiplicityTagged.
  5/11/16 (aem,mac): Move to u3shell namespace.

****************************************************************/

#ifndef MOSHINSKY_H_
#define MOSHINSKY_H_

#include "eigen3/Eigen/Eigen"

#include "sp3rlib/u3.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/u3st_scheme.h"

namespace u3shell
{
  // double MoshinskyCoefficient(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& xr,const u3::SU3& xc,const u3::SU3& x);
  // SU(3) Moshinsky Coefficient which is equivalent to a Wigner little d function evaluated at pi/2

  // double MoshinskyCoefficient(const u3::U3& w1, const u3::U3& w2, const u3::U3& wr,const u3::U3& wc,const u3::U3& w);
  //Overleading for U3 arguements 

  double 
  MoshinskyCoefficient(int r, int R, int r1, int r2, const u3::SU3& x);
  //Overloading Moshinsky to take integer arguements for two-body and relative-center of mass arguments
  // and SU(3) total symmetry (lambda,mu)

  inline double 
  MoshinskyCoefficient(int r, int R, int r1, int r2, const u3::U3& w)
  {
    return MoshinskyCoefficient(r,R,r1,r2,w.SU3());
  }
  // Overloading Moshinsky to take integers and U3 for total symmetry


  Eigen::MatrixXd 
  MoshinskyTransform(
    const u3::SU3& x0, 
    int etap,
    int eta,
    const u3shell::TwoBodySubspaceU3ST& bra_subspace, 
    const u3shell::TwoBodySubspaceU3ST& ket_subspace, 
    int rho0,
    std::string normalization="NAS"
  );

  // u3shell::TwoBodyUnitTensorCoefficientsU3ST 
  // MoshinskyTransformUnitTensor(const u3shell::RelativeUnitTensorLabelsU3ST& tensor, u3shell::TwoBodySpaceU3ST& space);
  // // Moshinsky transform of relative unit tensor operator to twobody space and anti-symmeterizes 
 
  // The transformed coefficients are stored in the two_body_expansion container which is a 
  // std::map<u3shell::TwoBodyUnitTensorLabelsU3ST,double>

  void 
  MoshinskyTransformUnitTensor(
    const u3shell::RelativeUnitTensorLabelsU3ST& tensor, 
    double expansion_coef, 
    u3shell::TwoBodySpaceU3ST& space,
    u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_expansion,
    std::string normalization="NAS"
  );

  void
  MoshinskyTransformUnitTensor(
    const u3shell::RelativeUnitTensorLabelsU3ST& tensor, 
    std::map<u3shell::RelativeCMUnitTensorLabelsU3ST,double>& unit_relative_cm_expansion,
    u3shell::TwoBodySpaceU3ST& space,
    u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_expansion,
    std::string normalization="NAS"
  );




// Moshinsky Transform operator decomposed in terms of unit tensors to two-body nomralized anti-symmeterized space 
  
  void
  TransformRelativeTensorToTwobodyTensor(
    const RelativeUnitTensorCoefficientsU3ST& relative_unit_tensor_expansion, 
    u3shell::TwoBodySpaceU3ST& space,
    u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_unit_tensor_expansion, 
    std::string normalization="NAS"
  );
  
} //namespace

#endif 






/****************************************************************
  moshinsky_unitarity.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame
 
  SPDX-License-Identifier: MIT

  10/20/16 (aem,mac): Created to test unitarity of moshinsky 
  						        transformation.
****************************************************************/
#include "sp3rlib/u3coef.h"
#include "moshinsky/moshinsky_xform.h"
#include "u3shell/relative_operator.h"

int main(int argc, char **argv)
{
  // initialize su3lib
  u3::U3CoefInit();
  int Nmax=6;
  // Generate a space 
  u3shell::TwoBodySpaceU3ST two_body_space(Nmax);
  u3shell::RelativeCMSpaceU3ST relative_cm_space(Nmax);
  bool norm=true;
  for(int subspace_index=0; subspace_index<two_body_space.size(); ++subspace_index)
    {
      const u3shell::TwoBodySubspaceU3ST& two_body_subspace=two_body_space.GetSubspace(subspace_index);
      const u3shell::RelativeCMSubspaceU3ST& relative_cm_subspace=relative_cm_space.GetSubspace(subspace_index);
      u3::SU3 x=two_body_subspace.omega().SU3();

      Eigen::MatrixXd matrix=Eigen::MatrixXd::Zero(relative_cm_subspace.size(),two_body_subspace.size());
      Eigen::MatrixXd matrixT=Eigen::MatrixXd::Zero(relative_cm_subspace.size(),two_body_subspace.size());

      for(int relative_cm_index=0; relative_cm_index<relative_cm_subspace.size(); ++relative_cm_index)
        for(int two_body_index=0; two_body_index<two_body_subspace.size(); ++two_body_index)
          {
            const u3shell::TwoBodyStateU3ST two_body_state(two_body_subspace,two_body_index);
            const u3shell::RelativeCMStateU3ST relative_cm_state(relative_cm_subspace,relative_cm_index);
            int Nr=relative_cm_state.Nr();
            int Ncm=relative_cm_state.Ncm();
            int N1=two_body_state.N1();
            int N2=two_body_state.N2();
            double coef=norm?(1./std::sqrt(1.+KroneckerDelta(N1,N2))):1;
            matrix(relative_cm_index,two_body_index)=std::sqrt(2)*coef*u3shell::MoshinskyCoefficient( Nr, Ncm,N1, N2, x); //Match Mark
          }    

      Eigen::MatrixXd product=matrix*matrix.transpose();
      // Eigen::MatrixXd product=matrix*matrixT;

      std::cout<<"Subspace "<<subspace_index<<std::endl;
      // std::cout<<two_body_subspace.DebugStr()<<std::endl;
      std::cout<<matrix<<std::endl<<std::endl<<product<<std::endl<<std::endl;
    }
}
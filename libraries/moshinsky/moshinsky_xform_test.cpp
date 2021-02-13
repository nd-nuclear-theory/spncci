/****************************************************************
  moshinsky_xform_test.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/7/16 (aem,mac): Created to test moshinky.h, moshinsky.cpp.
  5/11/16 (aem,mac): Update namespace.

  2/13/18 (mac): WARNING: This program was being left out of the
  build due to a filename error.  It no longer builds successfully.
  It seems to depend upon a nonexistent header moshinsky.h.

****************************************************************/

#include "moshinsky/moshinsky.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/relative_operator.h"

int main(int argc, char **argv)
{
  // initialize su3lib
  u3::U3CoefInit();
  int Nmax=8;
  // Generate a space 
  u3shell::TwoBodySpaceU3ST space(Nmax);


	u3shell::RelativeStateLabelsU3ST bra(2,1,0);
	u3shell::RelativeStateLabelsU3ST ket(2,1,0);

  // Coefficient map for a single unit tensor has only one entry
  u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels(u3::SU3(2,2),1,0,bra,ket);
  u3shell::RelativeUnitTensorCoefficientsU3ST unit_tensor;
  unit_tensor[unit_tensor_labels]=1;

  u3shell::TwoBodyUnitTensorCoefficientsU3ST unit_tensor_two_body;
  u3shell::TransformRelativeTensorToTwobodyTensor(unit_tensor, space,unit_tensor_two_body);
  std::cout<<std::endl<<"Relative Unit Tensor  "<<unit_tensor_labels.Str()<<std::endl;
  for (auto key_value : unit_tensor_two_body)
  {
    // extract unit tensor labels and coefficients
    auto labels= key_value.first;
    auto coefficient = key_value.second;

    // std::cout<<labels.Str()<<std::endl
    // <<coefficient<<std::endl;
  }
  // std::cout<<"number of unit tensors"<<unit_tensor_two_body.size()<<std::endl;

  // Identity test 
  std::cout<<std::endl<<"Identity Test "<<std::endl;
  u3shell::RelativeUnitTensorCoefficientsU3ST identity;
  for (int N=0; N<=Nmax; N++)
    for (int S=0;S<=1; S++)
      for (int T=0;T<=1; T++)
        {
          if ((N+S+T)%2==1)
          {
            u3shell::RelativeStateLabelsU3ST bra(N,S,T);
            u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
            u3shell::RelativeUnitTensorLabelsU3ST tensor(u3::SU3(0,0),0,0,bra,ket);
            identity[u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,0),0,0,bra,ket)]=1;
            // std::cout<<tensor.Str()<<std::endl;
          }
        }

  // // testing moshinksy transformation 
  // u3shell::RelativeUnitTensorLabelsU3ST tensor(u3::SU3(0,0),0,0,bra,ket);
  u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion;
  u3shell::TransformRelativeTensorToTwobodyTensor(identity, space,two_body_expansion);

  for (auto key_value : two_body_expansion)
  {
    // extract unit tensor labels and coefficients
    auto labels= key_value.first;
    auto coefficient = key_value.second;

    // std::cout<<labels.Str()<<std::endl
    // <<coefficient<<std::endl;
  }

  // // Number operator test 
  // std::cout<<std::endl<<"Number operator test "<<std::endl;
  // u3shell::RelativeUnitTensorCoefficientsU3ST number_operator;
  // for (int N=0; N<=Nmax; N++)
  //   for (int S=0;S<=1; S++)
  //     for (int T=0;T<=1; T++)
  //       if ((N+S+T)%2==1)
  //         {
  //           u3shell::RelativeStateLabelsU3ST bra(N,S,T);
  //           u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
  //           number_operator[u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,0),0,0,bra,ket)]=u3shell::RelativeNumberOperator(bra,ket);
  //         }

  // testing moshinksy transformation 
  // u3shell::RelativeUnitTensorLabelsU3ST tensor(u3::SU3(0,0),0,0,bra,ket);
  // u3shell::TransformRelativeTensorToTwobodyTensor(number_operator, space, two_body_expansion);

  for (auto key_value : two_body_expansion)
  {
    // extract unit tensor labels and coefficients
    auto labels= key_value.first;
    auto rme = key_value.second;
    u3::SU3 x0(labels.x0());
    u3::SU3 xp(labels.ket().x());
    u3::SU3 x(labels.bra().x());
    int rho0=labels.rho0();

    HalfInt S0=labels.S0();
    HalfInt S=labels.ket().S();
    HalfInt Sp=labels.bra().S();
    std::tuple<HalfInt,HalfInt,HalfInt> S_tuple(S,S0,Sp);
    
    HalfInt T0=labels.T0();
    HalfInt T=labels.ket().T();
    HalfInt Tp=labels.bra().T();
    std::tuple<HalfInt,HalfInt,HalfInt> T_tuple(T,T0,Tp);
    
    int g0=labels.g0();

    MultiplicityTagged<int>::vector L_kappa=u3::BranchingSO3(x);
    MultiplicityTagged<int>::vector L0_kappa0=u3::BranchingSO3(x0);
    MultiplicityTagged<int>::vector Lp_kappap=u3::BranchingSO3(xp);
    std::map<
      std::tuple<
            std::tuple<int,int,int>,
            std::tuple<HalfInt,HalfInt,HalfInt>,
            std::tuple<HalfInt,HalfInt,HalfInt>,
            int 
            >,
      double
    > L_map;

    for(int k=0; k<L_kappa.size(); k++)
      {
        int L=L_kappa[k].irrep;
        int kappa_max=L_kappa[k].tag;
        for(int k0=0; k0<L0_kappa0.size(); k0++)
          {
            int L0=L0_kappa0[k0].irrep;
            int kappa0_max=L0_kappa0[k0].tag;
            for(int kp=0; kp<Lp_kappap.size(); kp++)
              {
                int Lp=Lp_kappap[kp].irrep;
                int kappap_max=Lp_kappap[kp].tag;
                double branched_coef=0;
                for(int kappa=1; kappa<=kappa_max; kappa++)
                  for(int kappa0=1; kappa0<=kappa0_max; kappa0++)
                    for(int kappap=1; kappap<=kappap_max; kappap++)
                      {
                        branched_coef+=u3::W(x, kappa, L, x0, kappa0, L0,xp, kappap, Lp, rho0);
                      }
                std::tuple<int,int,int> L_tuple(L,L0,Lp);

                L_map[std::make_tuple(L_tuple,S_tuple,T_tuple,g0)]+=branched_coef*rme;
              }
          }
      }

    std::cout<<labels.Str()<<std::endl
    <<rme<<std::endl;
  }


}

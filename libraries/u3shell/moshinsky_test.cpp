/****************************************************************
  moshinsky_test.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created to test moshinky.h, moshinsky.cpp.
  5/11/16 (aem,mac): Update namespace.
****************************************************************/

#include "u3shell/moshinsky.h"

int main(int argc, char **argv)
{
  // initialize su3lib
  u3::U3CoefInit();
  int Nmax=4;
  // Generate a space 
  u3shell::TwoBodySpaceU3ST space(Nmax);


  std::cout<< "Moshinsky coefficient " <<u3shell::MoshinskyCoefficient(2,2,4,0,u3::SU3(4,0))<<std::endl;

	u3shell::RelativeStateLabelsU3ST bra(2,1,0);
	u3shell::RelativeStateLabelsU3ST ket(2,1,0);

  // Coefficient map for a single unit tensor has only one entry
  u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels(u3::SU3(2,2),1,0,bra,ket);
  u3shell::RelativeUnitTensorCoefficientsU3ST unit_tensor;
  unit_tensor[unit_tensor_labels]=1;

  u3shell::TwoBodyUnitTensorCoefficientsU3ST unit_tensor_two_body=u3shell::TransformRelativeTensorToTwobodyTensor(unit_tensor, space);
  std::cout<<std::endl<<"Relative Unit Tensor  "<<unit_tensor_labels.Str()<<std::endl;
  for (auto key_value : unit_tensor_two_body)
  {

    // extract unit tensor labels and coefficients
    auto labels= key_value.first;
    auto coefficient = key_value.second;

    std::cout<<labels.Str()<<std::endl
    <<coefficient<<std::endl;
  }
  // std::cout<<"number of unit tensors"<<unit_tensor_two_body.size()<<std::endl;

  // Identity test 
  std::cout<<std::endl<<"Identity Test "<<std::endl;
  u3shell::RelativeUnitTensorCoefficientsU3ST identity;
  for (int N=0; N<=Nmax; N++)
    for (int S=0;S<=1; S++)
      for (int T=0;T<=1; T++)
        {
          u3shell::RelativeStateLabelsU3ST bra(N,S,T);
          u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
          identity[u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,0),0,0,bra,ket)]=1;
        }

  // // testing moshinksy transformation 
  // u3shell::RelativeUnitTensorLabelsU3ST tensor(u3::SU3(0,0),0,0,bra,ket);
  u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion=u3shell::TransformRelativeTensorToTwobodyTensor(identity, space);

  for (auto key_value : two_body_expansion)
  {
    // extract unit tensor labels and coefficients
    auto labels= key_value.first;
    auto coefficient = key_value.second;

    std::cout<<labels.Str()<<std::endl
    <<coefficient<<std::endl;
  }

}

/****************************************************************
  relative_operator_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include <iostream>
#include <vector>

#include "fmt/format.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/relative_operator.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

void test_relative()
{

      // declare coefficient containers
      u3shell::RelativeUnitTensorCoefficientsU3ST relative_unit_tensor_coefficients;

      // populate operator
      //
      // TO CHECK: if labels make sense
      u3shell::RelativeUnitTensorLabelsU3ST
        relative_unit_tensor_labels(
            u3::SU3(1,1),1,0,  // x0,S0,T0
            u3shell::RelativeStateLabelsU3ST(1,1,1), // bra
            u3shell::RelativeStateLabelsU3ST(1,1,1)  // ket
          );
      relative_unit_tensor_coefficients[relative_unit_tensor_labels] = 1;

      // dump operator
      std::cout << "relative_unit_tensor_coefficients" << std::endl;
      for (auto key_value : relative_unit_tensor_coefficients)
        {

          // extract unit tensor labels and coefficients
          auto labels= key_value.first;
          double coefficient = key_value.second;
        
          std::cout 
            << fmt::format(
                "  {} {:e}",
                labels.Str(),
                coefficient
              )
            << std::endl;
        }

}


int main(int argc, char **argv)
{

  u3::U3CoefInit();
  test_relative();

  if(false)
  {
    int Nmax=2;
    int N1v=1;
    u3shell::RelativeSpaceU3ST space(Nmax);
    std::map<int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>> relative_unit_tensor_labels;

    u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, N1v,relative_unit_tensor_labels);
    for(auto it=relative_unit_tensor_labels.begin(); it!=relative_unit_tensor_labels.end(); ++it)
      {
        std::cout<<it->first<<std::endl;
        std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_vec=it->second;
        for(auto tensor : unit_vec)
          std::cout<<tensor.Str()<<std::endl;
      }

    {
      int Nmax=2;
      int N1v=1;
      std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
      u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, N1v,relative_unit_tensors);
      std::cout<<std::endl<<relative_unit_tensors.size()<<std::endl;
      int index=0;
      for(auto& tensor : relative_unit_tensors)
       {
          std::cout<<index<<"   "<<tensor.Str()<<std::endl;
          index++;
        }  
    }

    {
      int Nmax=4;
      int N1v=1;
      std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
      u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, N1v,relative_unit_tensors);
      std::cout<<std::endl<<relative_unit_tensors.size()<<std::endl;
      int index=0;
      for(auto& tensor : relative_unit_tensors)
       {
          std::cout<<index<<"   "<<tensor.Str()<<std::endl;
          index++;
        }
    }
  }
  




if(true)
{
    int Nmax=20;
    int S=0;
    int T=0;
    std::cout<<"B relative operator"<<std::endl;
    for(int N=2; N<=Nmax; N++)
      {
        int Np=N-2;
        u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
        u3shell::RelativeStateLabelsU3ST bra(Np,S,T);              
        double Brme1=u3shell::Brel(bra,ket);
        double Brme2=std::sqrt((N+2)*(N+1)/2);
        std::cout<<fmt::format("{:.4f}  {:.4f} ",Brme1, Brme2)<<std::endl;
      }

    // Arel
    std::cout<<std::endl<<"A relative operator"<<std::endl;
    for(int N=0; N<=Nmax; N++)
      {
        int Np=N+2;
        u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
        u3shell::RelativeStateLabelsU3ST bra(Np,S,T);   

        double Arme1=u3shell::Arel(bra,ket);
        double Arme2=std::sqrt((N+2)*(N+1)/2);
        std::cout<<fmt::format("{:.4f}  {:.4f} ",Arme1, Arme2)<<std::endl;
      }


}



  // termination
  return 0;
}

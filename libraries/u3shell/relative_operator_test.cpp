/****************************************************************
  relative_operator_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  5/25/16 (aem): Created.
****************************************************************/
#include <iostream>
#include <vector>

#include "fmt/format.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/relative_operator.h"

#include "spncci_basis/recurrence_indexing_operator.h"

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

  u3::U3CoefInit(100);

  if(false)
    test_relative();

  if(true)
  {
    int Nmax=2;
    int N1v=1;

    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
    u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1v,relative_unit_tensor_labels);


    std::set<unsigned int> Allowed_L0_values;
    std::set<unsigned int> Allowed_S0_values;
    std::set<unsigned int> Allowed_T0_values;
    std::unordered_set<u3::U3> Allowed_w0_values;

    for(const auto& tensor : relative_unit_tensor_labels)
      {
        const auto& [N0,x0,S0,T0,g0] = tensor.operator_labels().Key();
        Allowed_w0_values.insert({N0,x0});
        Allowed_S0_values.insert(int(S0));
        Allowed_T0_values.insert(int(T0));
      }


    relative::OperatorParameters unit_tensor_parameters(
      N1v,
      Nmax,
      0,
      Allowed_w0_values,
      Allowed_L0_values,
      Allowed_S0_values,
      Allowed_T0_values
      );

    std::cout<<"Relative Sp(3,R) raising operator"<<std::endl;
    relative::spatial::OperatorSpace
      spatial_operator_space(unit_tensor_parameters);

    std::cout<<spatial_operator_space.DebugStr()<<std::endl;

    relative::spin::OperatorSpace
      spin_operator_space(unit_tensor_parameters);
    std::cout<<spin_operator_space.DebugStr()<<std::endl;

    relative::OperatorSectors operator_sectors(spatial_operator_space,spin_operator_space);
    std::cout<<operator_sectors.DebugStr()<<std::endl;


  }

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

int N1v=1;
int Nmax=2;




  // termination
  return 0;
}

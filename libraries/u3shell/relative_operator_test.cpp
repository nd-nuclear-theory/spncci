/****************************************************************
  relative_operator_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iostream>

#include "cppformat/format.h"
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

  u3shell::RelativeStateLabelsU3ST bra;
  u3shell::RelativeStateLabelsU3ST ket;

  bra=u3shell::RelativeStateLabelsU3ST(2,1,0);
  ket=u3shell::RelativeStateLabelsU3ST(2,1,0);
  std::cout<<"Number operator    "<<u3shell::RelativeNumberOperator(bra,ket)<<std::endl;
  std::cout<<"Kinetic operator   "<<u3shell::RelativeKineticEnergyOperator(bra,ket)<<std::endl;


  bra=u3shell::RelativeStateLabelsU3ST(4,1,0);
  ket=u3shell::RelativeStateLabelsU3ST(2,1,0);
  std::cout<<"Raising operator   "<<u3shell::RelativeSp3rRaisingOperator(bra,ket)<<std::endl;
  std::cout<<"Kinetic operator   "<<u3shell::RelativeKineticEnergyOperator(bra,ket)<<std::endl;


  bra=u3shell::RelativeStateLabelsU3ST(2,1,0);
  ket=u3shell::RelativeStateLabelsU3ST(4,1,0);
  std::cout<<"Lowering operator  "<<u3shell::RelativeSp3rLoweringOperator(bra,ket)<<std::endl;
  std::cout<<"Kinetic operator   "<<u3shell::RelativeKineticEnergyOperator(bra,ket)<<std::endl;



  // termination
  return 0;
}

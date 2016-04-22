/****************************************************************
  u3coef_test.cpp

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created based on prototype u3.py and 
  T. Dytrych CSU3Master.
****************************************************************/

#include "am/halfint.h"
#include "utilities/utilities.h"
#include "utilities/multiplicity_tagged.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"

void iteration_test()
{

  // label set
  std::vector<u3::UCoefLabels> label_set;


  u3::SU3 x1(2,0);
  u3::SU3 x2(3,1);
  u3::SU3 x3(4,2);
  MultiplicityTagged<u3::SU3>::vector x12_vector=KroneckerProduct(x1,x2);
  MultiplicityTagged<u3::SU3>::vector x23_vector=KroneckerProduct(x2,x3);
  for(int i=0; i<x12_vector.size(); i++)
    {
      u3::SU3 x12=x12_vector[i].irrep;
      MultiplicityTagged<u3::SU3>::vector x_vector=KroneckerProduct(x12,x3);
      for(int j=0; j<x23_vector.size(); j++)
        {
          u3::SU3 x23=x23_vector[j].irrep;
          for(int k=0; k<x_vector.size(); k++)
            {
              u3::SU3 x=x_vector[k].irrep;
              if(u3::OuterMultiplicity(x1,x23,x)>0)
                label_set.push_back(u3::UCoefLabels(x1,x2,x,x3,x12,x23));
            }
        }
    }

  for (int a=0; a<label_set.size(); a++)
    {
      std::cout << label_set[a].Str() << " " << hash_value(label_set[a]) << std::endl;
    }
}

int main(int argc, char **argv)
{

  u3::SU3 x1(2,0);
  u3::SU3 x2(2,0);
  u3::SU3 x12(4,0);
  u3::SU3 x3(2,0);
  u3::SU3 x4(1,0);
  u3::SU3 x23(4,0);
  u3::SU3 x34(3,0);
  u3::SU3 x13(4,0);
  u3::SU3 x24(3,0);
  u3::SU3 x(6,0);
  u3::SU3 xx(7,0);


  u3::U3CoefInit();
  std::cout << U(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;
  std::cout << U(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;
  std::cout << Z(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;

  std::cout << W(x1,1,2,x2,1,2,x12,1,4,1) << std::endl;
  std::cout << W(x1,1,2,x2,1,2,x12,1,2,1) << std::endl;
  std::cout << W(x1,1,2,x2,1,2,x12,1,0,1) << std::endl;

  std::cout << Unitary9LambdaMu(
    x1,  x2,  x12, 1,
    x3,  x4,  x34, 1,
    x13, x24, xx,  1,
    1,    1,   1
    )<<std::endl;

  // block access
  std::cout << "U block test" << std::endl;
  u3::UCoefLabels labels(x1, x2, x, x3, x12, x23);
  u3::UCoefBlock block(labels);
  int r12_max, r12_3_max, r23_max, r1_23_max;
  std::tie(r12_max, r12_3_max, r23_max, r1_23_max) = block.Key();
  std::cout << "multiplicities " << r12_max << " " << r12_3_max << " " << r23_max << " " << r1_23_max << std::endl;
  std::cout << block.GetCoef(1,1,1,1) << std::endl;
  
  // test iteration over many labels
  iteration_test();


}

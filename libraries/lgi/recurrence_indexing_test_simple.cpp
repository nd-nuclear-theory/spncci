/****************************************************************
    recurrence_indexing_test.cpp

  Anna E. McCoy
  INT

  SPDX-License-Identifier: MIT

 6/25/21 (aem): Created.
****************************************************************/
#include "spncci_basis/recurrence_indexing.h"

#include <cppitertools/itertools.hpp>
#include <string>

#include "am/halfint_fmt.h"
#include "fmt/format.h"
#include "lgi/lgi.h"
#include "lgi/recurrence_lgi.h"
#include "mcutils/profiling.h"
#include "u3shell/operator_indexing_spatial.h"
#include "u3ncsm/dimensions.h"
#include "sp3rlib/u3coef.h"
// namespace test
// {
// }

int main(int argc, char** argv)
{

  // su3::init();
  u3::U3CoefInit(39);
  int Nsigma_max = 2;

  nuclide::NuclideType nuclide({3,3});
  auto Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide);
  unsigned int N0 = nuclide::N0ForNuclide(nuclide);
  int N1v = nuclide::ValenceShellForNuclide(nuclide);

  auto lgi_vector = lgi::GetLGIVector(nuclide,Nsigma0,Nsigma_max);

  // uint8_t parity_bar = 0;
  // spncci::spatial::OperatorConstraints operator_constraints(N1v,Nsigma0,parity_bar);
  // spncci::spatial::RecurrenceU3Space<u3shell::spatial::OneCoordType>
  // space({{Nsigma0,{0,0}},{Nsigma0,{0,0}}},operator_constraints);

  spncci::spin::Space<lgi::LGI> spin_space(lgi_vector, Nsigma_max);
  spncci::spin::RecurrenceSpace<lgi::LGI, u3shell::spin::twobody::OperatorLabelsST> spin_recurrence_space(spin_space, spin_space);
  spncci::spatial::Space spatial_space(spin_space,Nsigma0, Nsigma_max);
  spncci::spatial::RecurrenceSpace<u3shell::spatial::OneCoordType> spatial_recurrence_space(spatial_space,spatial_space,N1v,Nsigma0);

  std::cout<<"spin space "<<spin_space.size()<<std::endl;
  std::cout<<"spatial space "<<spatial_space.size()<<std::endl;
  std::cout<<"spin recurrence "<<spin_recurrence_space.size()<<std::endl;
  std::cout<<"spatial recurrence "<<spatial_recurrence_space.size()<<std::endl;


}

/****************************************************************
  u4_branching_test.cpp

  Patrick J. Fasano
  University of Notre Dame and Lawrence Berkeley National Laboratory

  SPDX-License-Identifier: MIT

  + 03/15/21 (pjf): Created.
****************************************************************/


#include "fmt/format.h"
#include "eigen3/Eigen/Eigen"
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "su4lib/u4_branching.h"

int main()
{
  // single multiplicity from Fig. 4c' in Draayer, J. Math. Phys. 11, 3225 (1970)
  constexpr unsigned int mult = su4lib::U4toSU2SU2Multiplicity(13,9,3,0,HalfInt(5,2),HalfInt(13,2));
  std::cout << "multiplicity for (4,6,3) -> (ST)=(5/2,13/2): " << mult << std::endl << std::endl;

  // Compare to Fig. 4a in Draayer, J. Math. Phys. 11, 3225 (1970)
  // (5,0,2) -> [7,2,2,0]
  std::cout << "# (5,0,2) as [7,2,2,0]" << std::endl
            << "# S   T  mult" << std::endl;
  for (HalfInt S = HalfInt(1,2); S < 4; ++S)
    for (HalfInt T = HalfInt(1,2); T < 4; ++T)
      std::cout << fmt::format("{:g} {:g}  {:d}", S, T, su4lib::U4toSU2SU2Multiplicity(7,2,2,0,S,T)) << std::endl;
  std::cout << std::endl;

  // Compare to Fig. 4a' in Draayer, J. Math. Phys. 11, 3225 (1970)
  // (5,6,2) -> [15,10,4,2]
  std::cout << "# (5,6,2) as [15,10,4,2]" << std::endl;
  Eigen::Matrix<unsigned int, 10, 10> mult_ap;
  // n.b. the actual (ST) pairs are shifted by 1/2
  for (int S = 0; S <= 9; ++S)
    for (int T = 0; T <= 9; ++T)
      mult_ap(S,T) = su4lib::U4toSU2SU2Multiplicity(15,10,4,2,S+HalfInt(1,2),T+HalfInt(1,2));
  std::cout << mult_ap << std::endl << std::endl;

  // Compare to Fig. 4b in Draayer, J. Math. Phys. 11, 3225 (1970)
  // (5,0,3) -> [8,3,3,0]
  std::cout << "# (5,0,3) as [8,3,3,0]" << std::endl
            << "# S T  mult" << std::endl;
  auto weights = su4lib::GenerateU4toSU2SU2Weights(8,3,3,0);
  for (const auto& [S,T,m] : weights)
  {
    std::cout << fmt::format("  {:g} {:g}  {:d}", S, T, m) << std::endl;
  }
  std::cout << std::endl;

  // Compare to Fig. 4b' in Draayer, J. Math. Phys. 11, 3225 (1970)
  // (5,7,3) -> [17,12,5,2]
  std::cout << "# (5,7,3) as [17,12,5,2]" << std::endl;
  Eigen::Matrix<unsigned int, 12, 12> mult_bp;
  for (int S = 0; S <= 11; ++S)
    for (int T = 0; T <= 11; ++T)
      mult_bp(S,T) = su4lib::U4toSU2SU2Multiplicity(17,12,5,2,S,T);
  std::cout << mult_bp << std::endl << std::endl;

  // (5,7,3) -> [16,11,4,1]
  std::cout << "# (5,7,3) as [16,11,4,1]" << std::endl;
  for (int S = 0; S <= 11; ++S)
    for (int T = 0; T <= 11; ++T)
      mult_bp(S,T) = su4lib::U4toSU2SU2Multiplicity(16,11,4,1,S,T);
  std::cout << mult_bp << std::endl;
}

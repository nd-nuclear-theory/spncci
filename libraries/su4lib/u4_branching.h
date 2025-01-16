/****************************************************************
  u4_branching.h

  U(4) -> SU(2)_S x SU(2)_T branching functions

  Patrick J. Fasano
  University of Notre Dame and Lawrence Berkeley National Laboratory

  SPDX-License-Identifier: MIT

  + 03/15/21 (pjf): Created.
****************************************************************/

#ifndef SU4LIB_U4_BRANCHING_
#define SU4LIB_U4_BRANCHING_

#include <vector>
#include <tuple>
#include "am/halfint.h"

namespace su4lib
{
  inline constexpr unsigned int RacahPhi(int z)
  {
    return (z>=0 ? z*z/4 : 0);
  }

  inline constexpr unsigned int RacahOmega(int f1, int f2, HalfInt S, HalfInt T)
  // Evaluate Racah's omega function for decomposing symmetric tensor [f1,f2,0,0]
  //
  // See Racah, Rev. Mod. Phys. 21, 494 (1949)
  {
    // first check preconditions on (f1,f2) multiplicity formula
    if (f1 < f2) return RacahOmega(f2, f1, S, T);
    if ((S < 0) || (T < 0)) return 0;
    if ((2*S > (f1+f2)) || (2*T > (f1+f2))) return 0;
    if (ParitySign(2*S) != ParitySign(f1+f2)) return 0;
    if (ParitySign(2*T) != ParitySign(f1+f2)) return 0;

    // Rev. Mod. Phys. 21, 494 (1949) -- eq. 28
    // break up into positive and negative parts
    unsigned int omega_pos = 0, omega_neg = 0;
    omega_pos += RacahPhi(f2 + 2 - int(abs(T-S)));
    omega_neg += RacahPhi(f2 + 1 - int(T+S));
    omega_pos += RacahPhi(int(T+S) - f1 - 1);
    omega_neg += RacahPhi(int((T+S) - abs(T-S)) - f1 + f2 + 1)/2;

    return omega_pos - omega_neg;
  }

  inline constexpr unsigned int U4toSU2SU2Multiplicity(
      int f1, int f2, int f3, int f4, HalfInt S, HalfInt T
    )
  // Evaluate inner multiplicity for U(4) -> SU(2)_S x SU(2)_T branching.
  {
    // Racah, Rev. Mod. Phys. 21, 494 (1949) -- eq. 31
    return RacahOmega(f1-f3, f2-f4, S, T)
      - (RacahOmega(f1-f4+1, f2-f3-1, S, T) + RacahOmega(f1-f2-1, f3-f4-1, S, T));
  }

  inline std::vector<std::tuple<HalfInt,HalfInt,unsigned int>>
  GenerateU4toSU2SU2Weights(int f1, int f2, int f3, int f4)
  {
    // Wigner's P is the max S0 or T0, therefore max S and max T
    HalfInt P = HalfInt((f1-f3)+(f2-f4),2);
    std::vector<std::tuple<HalfInt,HalfInt,unsigned int>> weights;
    for (HalfInt S = P; S >= 0 ; --S)
    {
      for (HalfInt T = P; T >= 0 ; --T)
      {
        auto multiplicity = U4toSU2SU2Multiplicity(f1, f2, f3, f4, S, T);
        if (multiplicity > 0) weights.push_back({S,T,multiplicity});
      }
    }

    return weights;
  }

}  // end namespace su4lib

#endif  // SU4LIB_U4_BRANCHING_

/****************************************************************
  nuclide.h

  Nuclide defintions

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  10/15/21 (aem): Created. Definitons extracted from lgi.h
  ****************************************************************/

#ifndef NUCLIDE_H_
#define NUCLIDE_H_

#include "am/am.h"

////////////////////////////////////////////////////////////////
// Calculation of Nsigma0
////////////////////////////////////////////////////////////////
namespace nuclide
{
using NuclideType = std::array<int,2>;

inline HalfInt Nsigma0ForNuclide(const NuclideType& nuclide, bool intrinsic=false)
// Calculate Nsigma0 for nuclide.
//
// This may be thought of as the dimensionless "oscillator energy in
// the lowest Pauli-allowed configuration", including zero-point
// energy.
//
// Example:
//
//   spncci::Nsigma0ForNuclide({3,3});
//
//      => returns 11
//
// Note that this is shorthand for
//
//   spncci::Nsigma0ForNuclide(spncci::NuclideType({3,3}))
//
// made possible by automatic conversion from initializer list to
// array.
//
// Arguments:
//   nuclide (input): (N,Z) for nucleus
//
// Returns:
//   Nsigma0

{
  // each major shell eta=2*n+l (for a spin-1/2 fermion) contains (eta+1)*(eta+2) substates
  HalfInt Nsigma0 = 0;
  for (int species_index : {0,1})
    {
      int num_particles = nuclide[species_index];
      for (int eta=0; num_particles>0; ++eta)
        {
          // add contribution from particles in shell
          int shell_degeneracy = (eta+1)*(eta+2);
          int num_particles_in_shell = std::min(num_particles,shell_degeneracy);
          // want num_particles_in_shell*(eta+HalfInt(3,2)), but HalfInt does not provide multiplication
          Nsigma0 += HalfInt(num_particles_in_shell*(2*eta+3),2);

          // discard particles in shell
          num_particles -= num_particles_in_shell;
        }
    }
  // If intrinsic remove cm zero point energy 3/2
  if(intrinsic)
    Nsigma0=Nsigma0-HalfInt(3,2);
  return Nsigma0;
}


  inline int ValenceShellForNuclide(const NuclideType& nuclide)
  // Calculate valence shell for nuclide.
  //
  // This is the highest oscillator shell occupied in the lowest
  // Pauli-allowed configuration (i.e., N1v).
  //
  // Example:
  //
  //   spncci::ValenceShellForNuclide({3,3});
  //
  //      => returns 1
  //
  // Note that this is shorthand for
  //
  //   spncci::ValenceShell(spncci::NuclideType({3,3}))
  //
  // made possible by automatic conversion from initializer list to
  // array.
  //
  //
  // Arguments:
  //   nuclide (input): (N,Z) for nucleus
  //
  // Returns:
  //   N1v
  {
    // each major shell eta=2*n+l (for a spin-1/2 fermion) contains (eta+1)*(eta+2) substates
    int N1v = 0;
    for (int species_index : {0,1})
      {
        int num_particles = nuclide[species_index];
        for (int eta=0; num_particles>0; ++eta)
          {
            // add contribution from particles in shell
            int shell_degeneracy = (eta+1)*(eta+2);
            int num_particles_in_shell = std::min(num_particles,shell_degeneracy);

            // register this shell as occupied
            N1v = std::max(N1v,eta);

            // discard particles in shell
            num_particles -= num_particles_in_shell;
          }
      }

    return N1v;
  }

}// end namespace 
#endif

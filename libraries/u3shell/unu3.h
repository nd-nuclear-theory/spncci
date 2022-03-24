/****************************************************************
  un.h

  U(N) to U(3) branching and generating allowed U(3)xSU(2) state labels
                                  
  Anna E. McCoy
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/7/16 (aem): Created based on lsu3shell CUNMASTER.cpp
  9/15/16 (aem): Changed storage container of allowed U3S irreps
****************************************************************/

#ifndef UN_H_
#define UN_H_

#include <map>  
#include <vector>
#include "sp3rlib/u3.h"
#include "sp3rlib/u3s.h"
#include "sp3rlib/u3st.h"
#include <unordered_map>


namespace un
{
  typedef std::vector<int> UNLabels;
  typedef std::vector<int> BasisStateWeightVector;
  typedef std::tuple<int,int,int> SingleParticleState;
  typedef std::unordered_map< u3::U3S,int,boost::hash<u3::U3S>> SingleShellAllowedU3SIrreps;
  // typedef MultiplicityTagged<u3::U3S>::vector SingleShellAllowedU3SIrreps;
    
  // u3::U3 UNWeightToU3Labels(BasisStateWeightVector weights, std::vector<SingleParticleState> single_particle_states);


  void GenerateU3Labels(
    const UNLabels& un_labels,
    const std::vector<SingleParticleState>& single_particle_states,
    BasisStateWeightVector& weights,
    std::map<u3::U3,int>& u3_un_multiplicity_map
  );
  /*
   Evaluates all allowed U(3) patterns [N_z, N_x,N_y] of U(N) basis states and their multiplicities.
   Algorithm based on M. Hamermesh's Group Theory and it's Applications to Physical Problems chapter 11
   and Computer Physics Communications 56 (1989) 279-290. 

   INPUT:
    un_labels (UNLabels): Gelfand parent row Young shape [f] = [f_{1},f_{2},...,f_{N}] 
        that labels an U(N) irrep
    
    single_particle_states: a distribution of HO quanta in (z, x, y) directions associated 
        with each of N levels of U(N)  

    weights: A vector of weights of a basis state of U(N) irrep. As function is called recursively, 
    an empty N-dimensional vector must be provided on input.

  OUTPUT: 

  Map of U(3) pattern labels associated with a corresponding multiplicity.
  */


  // Wrapper for GenerateU3Labels providing empty vector weights that is necessary in 
  // recursive calls but not for inital call to funciton
  inline void GenerateU3Labels(
      const UNLabels& un_labels,
      const std::vector<SingleParticleState>& single_particle_states,
      std::map<u3::U3,int>& u3_un_multiplicity_map
    )
  { 
    BasisStateWeightVector weights(un_labels.size());
    GenerateU3Labels(un_labels, single_particle_states, weights, u3_un_multiplicity_map);
  }


  inline void GenerateSingleParticleStates(int n, std::vector<SingleParticleState>& single_particle_states)
  {
    for (long k = 0; k <= n; k++)
        for (long nx = k; nx >= 0; nx--)
        {
            SingleParticleState sp(n-k,nx,k-nx);
            single_particle_states.push_back(sp);
        }
  }

unsigned UNBranchingMultiplicity(const u3::U3 w, const std::map<u3::U3,int>& u3_un_multiplicity_map);


void 
GenerateAllowedSU3xSU2Irreps(
    const unsigned n, const unsigned A, 
    SingleShellAllowedU3SIrreps& allowed_irreps
  );
  /*
  Generates a list of allowed U(3)xSU(2) irreps that obey the antisymmetry constrain on A particles
  in the same major oscillator shell 

  INPUT: 
    n : major oscillator shell 
    A : number of particles in the shell

  OUTPUT: 
    allowed_irreps : a map with key-value pair of allowed U3S labels and number of occurances. 
  */


void GenerateAllowedSU3xSU2xSU2TwoBodyIrreps(
          const unsigned n,
          MultiplicityTagged<u3::U3ST>::vector& allowed_irreps
        );
// Currently not working 
} //end namespace

#endif
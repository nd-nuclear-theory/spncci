/****************************************************************
  un.h

  U(N) to U(3) branching and generating allowed U(3)xSU(2) state labels
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created based on lsu3shell CUNMASTER.cpp
****************************************************************/

#ifndef UN_H_
#define UN_H_

#include <map>  
#include <vector>
#include "sp3rlib/u3.h"



namespace un
{
  typedef std::vector<int> UNLabels;
  typedef std::vector<int> BasisStateWeightVector;
  typedef std::tuple<int,int,int> SingleParticleState;


	
	// u3::U3 UNWeightToU3Labels(BasisStateWeightVector weights, std::vector<SingleParticleState> single_particle_states);

	/*
	 * INPUT:
	 * (1) un_labels (Gelfand parent row):
	 * Young shape [f] = [f_{1},f_{2},...,f_{N}] that labels an U(N) irrep
	 * (2) single_particle_states:
	 * a distribution of HO quanta in (z, x, y) directions associated with each of
	 * N levels of U(N)  
	 * (3) weights:
	 * A weight vector of a basis state of U(N) irrep. This parameter is used in a 
	 * recursive execution. An empty N-dimensional vector must be provided on input.
	 *
	 * TASK: 
	 * evaluate all allowed U(3) patterns [N_{z}, N_{x}, N_{y}] of U(N) basis
	 * states and their multiplicities. Algorithm described in Computer Physics
	 * Communications 56 (1989) 279-290
	 *
	 * OUTPUT: 
	 * Map of [N_{z}, N_{x}, N_{y}] labels associated with a corresponding multiplicity.
	 * This map can be used to calculate allowed SU(3)xSU(2) irreps
	 */
	void GenerateU3Labels(
      const UNLabels& un_labels,
      const std::vector<SingleParticleState>& single_particle_states,
      BasisStateWeightVector& weights,
      std::map<u3::U3,int>& u3_un_multiplicity_map
    );


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

unsigned UNBranchingMultiplicity(const u3::U3 w, const std::map<u3::U3,int>& u3_un_multiplicity_map);

void GetAllowedSU3xSU2Irreps(
        const unsigned n, const unsigned A, 
        MultiplicityTagged<u3::U3S>::vector& allowed_irreps
      );

} //end namespace 

#endif
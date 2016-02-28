#ifndef CUNMASTER_H
#define CUNMASTER_H
#include <UNU3SU3/UNU3SU3Basics.h>
#include <map>
#include <set>
#include <bitset>

class CUNMaster {
	private:
	void GenerateU3Labels(	const UN::LABELS& f, 
							const U3::SPS& hosps, 
							UN::BASIS_STATE_WEIGHT_VECTOR& vWeights,
							UN::U3MULT_LIST& mU3LabelsOccurance);
	void GenerateU3LabelsUNWeights(	const UN::LABELS& f,	
									const U3::SPS& ShellSPS, 
									UN::BASIS_STATE_WEIGHT_VECTOR& vWeights,
									UN::U3MULT_LIST& mU3LabelsOccurance,
									std::map<U3::LABELS, std::set<UN::BASIS_STATE_WEIGHT_VECTOR> >&  mUNweights);
	private:
	void Weight2U3Label(	const UN::BASIS_STATE_WEIGHT_VECTOR& vWeights,
							const U3::SPS& ShellSPS, 
							U3::LABELS& vU3Labels);
	public:
/*	
 * TASK: 
 * Calculate [N_{z}, N_{x}, N_{y}] quantum numbers of U(N) irrep basis states and their 
 * occurance by generating all allowed Gel'fand patterns of a given U(N) irrep. 
 * This is needed for Draayer algorithm which calculates allowed SU(3)xSU(2) irreps in irrep of U(N).
 *
 * OUTPUT: 
 * Map of [N_{z}, N_{x}, N_{y}] labels associated with a corresponding number
 * of occurances.  This map is used to evaluate allowed SU(3)xSU(2) irreps.
 */
	void GenerateU3Labels(	const UN::LABELS& f, 
							const U3::SPS& ShellSPS, 
							UN::U3MULT_LIST& mU3LabelsOccurance) {
		UN::BASIS_STATE_WEIGHT_VECTOR Weight(f.size());
		GenerateU3Labels(f, ShellSPS, Weight, mU3LabelsOccurance);};

/*	
 * TASK: 
 * Evaluate all basis state of U(N) irrep |[f] N_{x}; N_{y}; N_{z}> and calculate how many times they occur (degeneracy)
 * If N_{z}>=N_{x}>=N_{y} store bitset representation of |[f] N_{x}; N_{y}; N_{z}> in mHWSbitsets.
 *
 * OUTPUT: 
 * 1) Map of [N_{z}, N_{x}, N_{y}] labels associated with a corresponding degeneracy. 
 * This map is used to calculate allowed SU(3)xSU(2) irreps.
 * 2) Map of [N_{z}, N_{x}, N_{y}] labels such that N_{z} >= N_{x} >= N_{y} and
 * associated vector of U(N) basis states represented as a pair of bitsets.
 * This is used to calculate HWS of SU(3)xSU(2) irreps.
 */
//	TODO: this version implements bitset as std::bitset<64>. The final version must be more generic, that is, able
//	to work with different 
	void GenerateU3LabelsUNWeights(	const UN::LABELS& f,
										const U3::SPS& ShellSPS,
										UN::U3MULT_LIST& mU3LabelsOccurance,
										std::map<U3::LABELS, std::set<UN::BASIS_STATE_WEIGHT_VECTOR> >&  mUNweights)
							{
								UN::BASIS_STATE_WEIGHT_VECTOR Weight(f.size());
								GenerateU3LabelsUNWeights(f, ShellSPS, Weight, mU3LabelsOccurance, mUNweights);
							};
	public:
// TASK:
// Return the number of times an U(3) irrep (U3Labels) occurs in U(N) irrep
// [f].  This is usually denoded as alpha, that is, [f] \alpha (N_{z} N_{x}
// N_{y}).  mU3IRs is a map of labels (NZ, NX, NY) in basis of an U(N) irrep
// [f] associated with the number of basis states of the U(N) irrep [f] that
// carry this HO quanta 
	static unsigned GetMultiplicity(const U3::LABELS U3Labels, const UN::U3MULT_LIST& mU3IRs); // returns 
/*
 * INPUT:
 * (1) n:
 * harmonic oscillator shell
 * (2) A:
 * number of fermions
 * (3) ShellStructure:
 * A distribution of HO quanta in (z, x, y) directions associated with each of
 * N=(n+1)(n+2)/2 levels of U(N) 
 *
 * TASK:
 * Use algorithm described in Computer Physics Communications 56 (1989) 279-290
 * to obtain a set of allowed SU(3)xSU(2) irreps for a system of A fermions in n-th HO shell
 *
 * OUTPUT:
 * std::vector of SU(3)xSU(2) irrep labels
 */

void GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, const U3::SPS& ShellSPS, std::vector<UN::SU3_VEC>& SpinSU3MULTLabels);
void GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, std::vector<UN::SU3_VEC>& SpinSU3MULTLabels);

void GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, const U3::SPS& ShellSPS, std::vector<std::pair<SU2::LABEL, UN::SU3_VEC> >& SpinSU3MULTLabels);
void GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, std::vector<std::pair<SU2::LABEL, UN::SU3_VEC> >& SpinSU3MULTLabels);


void GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, const U3::SPS& ShellSPS, UN::SU3xSU2_VEC& AllowedIrreps);
void GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, UN::SU3xSU2_VEC& AllowedIrreps);

void GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, const U3::SPS& ShellSPS, std::vector<SU3xSU2::LABELS>& AllowedIrreps);
void GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, std::vector<SU3xSU2::LABELS>& AllowedIrreps);
};
#endif

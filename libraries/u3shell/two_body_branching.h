/****************************************************************
  two_body_branching.h

  Branching U(3)xSU(2)xSU(2) two-body rme's
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/7/17 (aem): Created.
****************************************************************/
#ifndef TWO_BODY_BRANCHING_H_
#define TWO_BODY_BRANCHING_H_

#include "cppformat/format.h"
#include "u3shell/tensor_labels.h"
#include "moshinsky/moshinsky_xform.h"

namespace u3shell
{
	// LST
	typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> TwoBodyStateLabelsLST;
	typedef std::tuple<int,HalfInt,HalfInt,TwoBodyStateLabelsLST,TwoBodyStateLabelsLST> TwoBodyBraketLST;
	//LSJT
	typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt,HalfInt> TwoBodyStateLabelsLSJT;
	typedef std::pair<TwoBodyStateLabelsLSJT,TwoBodyStateLabelsLSJT>TwoBodyBraketLSJT;
	//JJJT
	typedef std::tuple<int, int, HalfInt,int, int, HalfInt,HalfInt,HalfInt>TwoBodyStateLabelsJJJT;
	typedef std::pair<TwoBodyStateLabelsJJJT, TwoBodyStateLabelsJJJT> TwoBodyBraketJJJT;



	void H2FormatLookUp(int Nmax,std::map<std::tuple<int,int,HalfInt>,int>& h2_lookup);
	// Creates look up table for h2 format labels for NLJ single particle states
	//
	// Arguments:
	//	Nmax (input) : Two-body truncation
	//	h2_lookup (output) : Look up table for conversion of (N,L,J) single particle
	//		labels to h2 format.

	void BranchTwoBodyNLST(
  u3shell::IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
  std::map<u3shell::TwoBodyBraketLST,double>& two_body_rmes_lst
  );
  // Branch two-body SU(3)xSU(2)xSU(2) rme's to SO(3)xSU(2)xSU(2) rme's,
  // i.e., branch from (omega,S,T) to (L,S,T) labeled states.
	// 
	// Arguments:
	//	indexed_two_body_rmes (input) : two-body rme's labeled by
	//		unit tensor labels and kappa0,L0
	//	two_body_rme_lst (output) : branched rme's 

	void BranchTwoBodyLSJT( int Jmax, int J0,
	  const std::map<TwoBodyBraketLST,double>& two_body_rme_lst,
	  std::map<TwoBodyBraketLSJT,double>&two_body_rme_lsjt
	  );
	// Branch from LST to LSJT

	void branchJJJT(
	  std::map<TwoBodyBraketLSJT,double>& two_body_rme_lsjt,
	  std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt
	  );
	// Branch from LSJT to JJJT

	void BranchTwoBodyU3STToJJJT(int Jmax, int J0,
	  u3shell::IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
	  std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt
	  );
	// Branch from omegaST to JJJT

	void PrintTwoBodyMatrixElementsJJJT(const std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt);

}
#endif
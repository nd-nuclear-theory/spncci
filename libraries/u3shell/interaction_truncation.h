/****************************************************************
  interaction_truncations.h

  Truncations for relative two-body interaction 
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  11/2/17 (aem): Created.
****************************************************************/
#ifndef INTERACTION_TRUNCATION_H_
#define INTERACTION_TRUNCATION_H_
#include <iostream>
#include <fstream>
#include <unordered_set>
#include "u3shell/relative_branching.h"
#include "cppformat/format.h"

namespace u3shell
{
  void TruncateInteractionNrel(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    int Nrel_max
    );
  // Truncate interaction by Nrel.  Interaction RME is eliminated if either of the bra or ket eta is greater than Nrel_max.


  void TruncateInteractionN0(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    int N0_max
    );
  // truncate interaction by N0.  If abs(N0) is greater than N0_max, RME is eliminated
  

    void TruncateInteractionx0(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    int x0_max
    );
  // truncate interaction by (lambda, mu).  Interaction RME is eliminated if either lambda or mu is greater than LamMu_Max.
  

  void AccumulateInteractionByU3ST(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    std::unordered_map<u3shell::OperatorLabelsU3ST,double,boost::hash<u3shell::OperatorLabelsU3ST>>& probability_by_u3st
    );
    // Get distribution of interaction of U3ST tensor componenent (N0,x0,S0,T0) normalized to 1. 


  void TruncateInteractionU3ST(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    double truncation_threshold
    );
  // truncate interaction by individual U3ST relative rmes by eliminating those which fall
  // below a certain percentage of the total interaction decompostion
  // 
  // Truncation is only carried out on the N0>=0 RMEs.  The corresponding N0<0 matrix element
  // is included based on the value of the corresponding N0>0 RME. 


  typedef std::tuple<HalfInt,u3::SU3> SpacialSymmetries;
  void AccumulateInteractionByU3(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    std::unordered_map<SpacialSymmetries,double,boost::hash<SpacialSymmetries>>& probability_by_u3
    );
    // Get distribution of interaction of U3ST tensor componenent (N0,x0) normalized to 1. 


  void TruncateInteractionByU3(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    double truncation_threshold
    );
  // truncate interaction by U3 tensor character by eliminating rmes of tensor componenents U(1)xSU(3)
  // falls below a certain percentage of the total interaction decompostion
  // 
  // Truncation is only carried out on the N0>=0 RMEs.  The corresponding N0<0 matrix element
  // is included based on the value of the corresponding N0>0 RME. 

  /*typedef vector< tuple<int, int, int, int> > u3st;
  void Contribution(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    double sum_threshhold
    );*/


  void GetRelativeRMEsU3ST(
        const std::string& interaction_filename,
        int Nmax, int Jmax,
        u3shell::RelativeRMEsU3ST& interaction_u3st
        );
  // Read relative rmes in from file and upcouple

  void PrintRelativeRMEsU3ST(const u3shell::RelativeRMEsU3ST& interaction_u3st);
  // printing out u3st rmes 
  
  void WriteRelativeRMEsU3ST(const std::string& filename, const u3shell::RelativeRMEsU3ST& interaction_u3st);
  // writing out u3st rmes 


} // end namespace

#endif


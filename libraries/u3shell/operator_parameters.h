/****************************************************************
  operator_parameters.h

  Struct defining allowed quantum numbers for an operator
                                  
  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  4/6/22 (aem): Extracted from recurrence_indexing_operator.
****************************************************************/

#ifndef OPERATOR_PARAMETERS_H_
#define OPERATOR_PARAMETERS_H_

#include <unordered_set>
#include <fstream>
#include <string>
#include "sp3rlib/u3.h"


namespace u3shell::relative
{

  static constexpr unsigned int kNone = std::numeric_limits<unsigned int>::max();

  struct OperatorParameters
  /// Provides information about the relative two-body operator
  /// Includes SU(3), angular momentum, spin, and isospin
  /// quantum numbers for selection rules.
  {
    OperatorParameters()=default;

    OperatorParameters(
      const unsigned int Nbar_max_,
      const unsigned int J0_,
      const std::unordered_set<u3::U3>& Allowed_w0_values_,
      const std::set<unsigned int>& Allowed_L0_values_,
      const std::set<uint8_t>& Allowed_S0_values_,
      const std::set<uint8_t>& Allowed_T0_values_
    )
        :
        Nbar_max{Nbar_max_},
        J0{J0_},
        Allowed_w0_values{Allowed_w0_values_},
        Allowed_L0_values{Allowed_L0_values_},
        Allowed_S0_values{Allowed_S0_values_},
        Allowed_T0_values{Allowed_T0_values_}
    {
      for(const auto& S0 : Allowed_S0_values_)
        assert(S0==0 || S0==1 || S0==2);

      for(const auto& T0 : Allowed_T0_values_)
        assert(T0==0 || T0==1 || T0==2);
    }

    // Members
    unsigned int Nbar_max;
    unsigned int J0;
    std::set<unsigned int> Allowed_L0_values;
    std::set<uint8_t> Allowed_S0_values;
    std::set<uint8_t> Allowed_T0_values;
    std::unordered_set<u3::U3> Allowed_w0_values;
  };


  inline OperatorParameters CombineParameters(std::vector<OperatorParameters> parameters)
  {
    std::unordered_set<u3::U3> Allowed_w0_values;
    std::set<unsigned int> Allowed_L0_values;
    std::set<uint8_t> Allowed_S0_values;
    std::set<uint8_t> Allowed_T0_values;
    unsigned int Nbar_max=parameters[0].Nbar_max;
    unsigned int J0 = parameters[0].J0;
    for(const auto& param : parameters)
    {
      // If J0 values do not match, then eliminate selection rules based on J0
      // by setting J0 to kNone value.
      if(J0!=param.J0)
        J0=u3shell::relative::kNone;

      Nbar_max = std::max(Nbar_max,param.Nbar_max);
      Allowed_w0_values.insert(param.Allowed_w0_values.begin(),param.Allowed_w0_values.end());
      Allowed_L0_values.insert(param.Allowed_L0_values.begin(),param.Allowed_L0_values.end());
      Allowed_S0_values.insert(param.Allowed_S0_values.begin(),param.Allowed_S0_values.end());
      Allowed_T0_values.insert(param.Allowed_T0_values.begin(),param.Allowed_T0_values.end());
    }
    OperatorParameters new_parameters(
        Nbar_max, J0, Allowed_w0_values, Allowed_L0_values, Allowed_S0_values, Allowed_T0_values
      );
    return new_parameters;
  }

  /// Write header containing information on what's in the output file
  void WriteOperatorParametersHeader(std::ofstream& output);

  /// Write operator parameters to file.
  /// If L0 or J0 is kNone, then value is given as -1.
  ///
  /// Order within file is
  ///   Nbar_max J0
  ///   Num_w0 N01 lambda01 mu01  N02 lambda02 mu02 ...
  ///   Num_L0 L01 L02 L03...
  ///   Num_S0 S01 S02 S03...
  ///   Num_T0 T01 T02 T03...
  void WriteOperatorParameters(
      const OperatorParameters& parameters, std::ofstream& output
    );


  /// Read operator parameters from file.  Option bool value indicates if file contains
  /// header or not.
  OperatorParameters ReadOperatorParametersText(std::ifstream& input,const bool header_included=false);

  //////////////////////////////////////////////////////////////////////////////////////////
  // Operator parameters for commonly used operators
  //////////////////////////////////////////////////////////////////////////////////////////
  /// Operator parameters for identity operator
  inline u3shell::relative::OperatorParameters
  IdentityParameters(const unsigned int Nbar_max)
  {return u3shell::relative::OperatorParameters(Nbar_max,0u,{{0,{0u,0u}}},{0u},{0},{0});}

  /// Operator parameters for isoscalar quadrupole operator (T0=0)
  inline OperatorParameters
  QIsoscalarParameters(const unsigned int Nbar_max)
  {return OperatorParameters(Nbar_max,2u,{{0,{1u,1u}},{-2,{0u,2u}},{2,{2u,0u}}},{2u},{0},{0});}

  /// Operator parameters for isovector quadrupole operator (T0=1)
  inline OperatorParameters
  QIsovectorParameters(const unsigned int Nbar_max)
  {return OperatorParameters(Nbar_max,2u,{{0,{1u,1u}},{-2,{0u,2u}},{2,{2u,0u}}},{2u},{0},{1});}

  /// Operator parameters for isovector quadrupole operator (T0=1)
  inline OperatorParameters
  HamiltonianParameters(const unsigned int Nbar_max, const unsigned int T0_min=0u, const unsigned int T0_max=2u)
  {
    std::set<uint8_t> Allowed_T0_values;
    for(uint8_t T0 = T0_min; T0<=T0_max; ++T0)
      Allowed_T0_values.insert(T0);
    return OperatorParameters(Nbar_max,0u,{},{},{},Allowed_T0_values);}

  /// Operator parameters for kinetic energy
  inline OperatorParameters
  KineticEnergyParameters(const unsigned int Nbar_max)
  {return OperatorParameters(Nbar_max,0u,{{0,{1u,1u}},{-2,{0u,2u}},{2,{2u,0u}}},{0u},{0},{0});}

  /// Operator parameters for kinetic energy
  inline OperatorParameters
  KSquaredParameters(const unsigned int Nbar_max)
  {return OperatorParameters(Nbar_max,0u,{{0,{0u,0u}},{-2,{0u,2u}},{2,{2u,0u}}},{0u},{0},{0});}

  /// Operator parameters for kinetic energy
  inline OperatorParameters
  RSquaredParameters(const unsigned int Nbar_max)
  {return OperatorParameters(Nbar_max,0u,{{0,{0u,0u}},{-2,{0u,2u}},{2,{2u,0u}}},{0u},{0},{0});}


}//relative namespace



#endif

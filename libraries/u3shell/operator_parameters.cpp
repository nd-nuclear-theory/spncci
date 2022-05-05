/****************************************************************
  operator_parameters.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3shell/operator_parameters.h"

namespace u3shell::relative{

  void WriteOperatorParametersHeader(std::ofstream& output)
  {
    output << fmt::format("# Nbar_max J0\n");
    output << fmt::format("# Numw0 w0values\n");
    output << fmt::format("# NumL0 L0values\n");
    output << fmt::format("# NumS0 S0values\n");
    output << fmt::format("# NumT0 T0values\n");
  }

  void WriteOperatorParameters(
    const OperatorParameters& parameters,
    std::ofstream& output
  )
  {
    // Write Nbar_max and J0 to file
    int J0out = parameters.J0==u3shell::relative::kNone?-1:parameters.J0;
    output << fmt::format("{:3d} {:3d}\n",parameters.Nbar_max,J0out);

    std::string w0_string=fmt::format("{:3d}",parameters.Allowed_w0_values.size());
    for(const auto w0 : parameters.Allowed_w0_values)
      w0_string+=fmt::format(" {:3d} {:3d} {:3d}",w0.N(),w0.SU3().lambda(),w0.SU3().mu());
    output << fmt::format("{}\n",w0_string);

    // Write L values to file
    std::string L0_string;
    if(parameters.Allowed_L0_values.count(u3shell::relative::kNone))
      L0_string = fmt::format("{:3d} {:3d}\n",1,-1);
    else
      L0_string = fmt::format("{:3d}",parameters.Allowed_L0_values.size());
      for(const auto L0 : parameters.Allowed_L0_values)
        L0_string+=fmt::format(" {:3d}",L0);

    output << fmt::format("{}\n",L0_string);

    std::string S0_string=fmt::format("{:3d}",parameters.Allowed_S0_values.size());
    for(const auto S0 : parameters.Allowed_S0_values)
      S0_string+=fmt::format(" {:3d}",S0);
    output << fmt::format("{}\n",S0_string);

    std::string T0_string=fmt::format("{:3d}",parameters.Allowed_T0_values.size());
    for(const auto T0 : parameters.Allowed_T0_values)
      T0_string+=fmt::format(" {:3d}",T0);
    output << fmt::format("{}\n",T0_string);
  }

  OperatorParameters ReadOperatorParametersText(std::ifstream& input,const bool header_included)
  {
    // Skip header lines
    std::string line;
    std::istringstream line_stream;
    int line_count = 0;
    if(header_included)
      for(int i=0; i<5; ++i) std::getline(input, line);

    // Read in parameters
    int Nbar_max, J0in;
    std::getline(input,line);
    std::istringstream(line) >> Nbar_max >> J0in;
    unsigned int J0 = (J0in==-1)?u3shell::relative::kNone: J0in;

    // Read in allowed w0 values
    int num,lambda0,mu0,N0;
    std::unordered_set<u3::U3> Allowed_w0_values;
    std::getline(input,line);
    line_stream = std::istringstream(line);
    line_stream >> num;
    for(int i=0; i<num; ++i)
    {
      line_stream >> N0 >> lambda0 >> mu0;
      Allowed_w0_values.insert({N0,{lambda0,mu0}});
    }

    int label;

    // Read in Allowed L0 values
    std::set<unsigned int> Allowed_L0_values;
    std::getline(input,line);
    line_stream = std::istringstream(line);
    line_stream >> num;
    for(int i=0; i<num; ++i)
    {
      line_stream >> label;
      if(label == -1)
        Allowed_L0_values.insert(u3shell::relative::kNone);
      else
        Allowed_L0_values.insert(label);
    }

    // Read in Allowed S0 values
    std::set<uint8_t> Allowed_S0_values;
    std::getline(input,line);
    line_stream = std::istringstream(line);
    line_stream >> num;
    for(int i=0; i<num; ++i)
    {
      line_stream >> label;
      Allowed_S0_values.insert(label);
    }

    // Read in Allowed T0 values
    std::set<uint8_t> Allowed_T0_values;
    std::getline(input,line);
    line_stream = std::istringstream(line);
    line_stream >> num;
    for(int i=0; i<num; ++i)
    {
      line_stream >> label;
      Allowed_T0_values.insert(label);
    }

    return OperatorParameters(
        Nbar_max, J0, Allowed_w0_values, Allowed_L0_values, Allowed_S0_values, Allowed_T0_values
      );
  }


}

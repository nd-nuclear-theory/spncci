/****************************************************************
  generate_lsu3shell_relative_tensors.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/8/16 (aem,mac): Created.
****************************************************************/
#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_operator.h"
#include "spncci/sp_basis.h"
#include "spncci/lgi_unit_tensor.h"
#include "u3shell/unit_tensor_expansion.h"
// Given a nucleus (protons, neutrons)
// Given Nsigma_0
// Given a Nmax truncation
// Generate a list of unit tensors
// Generate Brel
// Generate Nrel
int main(int argc, char **argv)
{
  u3::U3CoefInit();


  if(argc<4)
    {
      std::cout<<"Syntax: Protons Neutrons Nmin Nmax "<<std::endl;
      std::exit(1);
    }
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  // will be either 1 or 2; 
  int Nstep=std::stoi(argv[4]);
  assert(Nstep<=2);
  int Nmin=Nmax%Nstep;
  // Set up unit tensor model space space
  std::string model_space=fmt::format("model_space.{}_{}_Nmax{}",Z,N,Nmax);
  std::ofstream model_stream(model_space);
  model_stream<<Z<<"  "<<N<<"  "<<-1<<std::endl;
  for(int Nex=Nmin; Nex<=Nmax; Nex+=2)
    model_stream<<Nex<<"  "<<-1<<std::endl;
  model_stream.close();

  //begin control file
  // first give specifications for unit tensors, then Brel and Nrel
  std::ofstream control_stream("operators.dat");
  control_stream
    <<model_space
    <<std::endl;

  //Generate all relative unit tensors up to Nmax cutoff
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensor_labels);
  lsu3shell::GenerateLSU3ShellOperator(Nmax, relative_unit_tensor_labels);


  // Generate Brel operator up to Nmax cutoff
  std::string brel_file=fmt::format("Brel_Nmax{:02d}",Nmax);
  u3shell::RelativeUnitTensorCoefficientsU3ST Brel_operator;
  u3shell::BrelRelativeUnitTensorExpansion(Nmin,Nmax, Brel_operator);
  lsu3shell::GenerateLSU3ShellOperator(Nmax, Brel_operator, brel_file);

  //Generate Nrel operator up to Nmax cutoff
  std::string nrel_file=fmt::format("Nrel_Nmax{:02d}",Nmax);
  u3shell::RelativeUnitTensorCoefficientsU3ST Nrel_operator;
  u3shell::NintrRelativeUnitTensorExpansion(Nmin,Nmax, Nrel_operator);
  lsu3shell::GenerateLSU3ShellOperator(Nmax, Nrel_operator, nrel_file);

  int num_unit=relative_unit_tensor_labels.size();
  // number of relative operators including Brel and Nrel
  int num_ops=num_unit+2;
  control_stream<<fmt::format("  operators {} ",num_ops)<<std::endl;
  for(int i=0; i<num_unit; ++i)
    control_stream<<fmt::format("unit{:06d}",i)<<std::endl;

  control_stream<<brel_file<<std::endl;
  control_stream<<nrel_file<<std::endl;
}

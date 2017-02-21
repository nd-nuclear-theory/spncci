/****************************************************************
  generate_lsu3shell_relative_tensors.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/8/16 (aem,mac): Created.
  11/7/16 (aem): Updated documentation.
  12/2/16 (aem): Add bool for U(N)->U(3) restriction on ops.
  2/16/17 (mac): Add Arel to output.
****************************************************************/

#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_operator.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"

// For a given nucleus, Nmax and N_step, generates input files for LSU3Shell 
// RecoupleSU3Interaction which recouples operators into format
// required by SU3RME. Input files include 
// 
// Control file (operators.dat) giving
//   The model space file name
//   "operators" number of operators
//   List of front part of operator file names, e.g. if file name is 
//      "unit000001.recoupler", then "unit000001" is given in file. 
//
// Operator files include:
//    All unit tensors defined on  relative space truncated by Nmax
//    Brel
//    Arel
//    Nintr (proportional to Nrel)
//
// Commandline inputs are Z N Nmax Nstep
//  Z : number of protons
//  N : number of neutrons
//  Nmax : max number of oscillator quant in relative basis
//  Nstep : indicates if all or only same parity spaces are included
//          in model space.  Can only take values of 1 or 2. 

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  // process arguments
  if(argc<6)
    {
      std::cout<<"Syntax: Protons Neutrons Nmax Nstep N1B"<<std::endl;
      std::exit(1);
    }
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  // will be either 1 or 2; 
  int Nstep=std::stoi(argv[4]);
  int N1B=std::stoi(argv[5]);
  assert((Nstep==2)||(Nstep==1));
  int Nmin=Nmax%Nstep;
  int A=N+Z;
  int T0=0;
  int J0=-1;
  bool un_u3_restrict=false;
  if( (N==0) || (Z==0) )
    un_u3_restrict=true;

  // set up lsu3shell model space file for unit tensor calculations
  int parity=(Nstep==1)?-1:Nmin;
  lsu3shell::GenerateModelSpaceFile(Z, N, Nmax, parity);

  // Generate all relative unit tensors up to Nmax cutoff
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax+2*N1B,relative_unit_tensor_labels,J0,T0,false);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, relative_unit_tensor_labels,un_u3_restrict);

  // Generate Brel operator up to Nmax cutoff
  std::string Brel_file_name_base=fmt::format("Brel_{:02d}_Nmax{:02d}",A,Nmax);
  std::string Brel_file_name=fmt::format("{}.recoupler",Brel_file_name_base);
  u3shell::RelativeUnitTensorCoefficientsU3ST Brel_operator;
  u3shell::BrelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B,Brel_operator,A);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B,Brel_operator,Brel_file_name,un_u3_restrict);

  // Generate Arel operator up to Nmax cutoff
  std::string Arel_file_name_base=fmt::format("Arel_{:02d}_Nmax{:02d}",A,Nmax);
  std::string Arel_file_name=fmt::format("{}.recoupler",Arel_file_name_base);
  u3shell::RelativeUnitTensorCoefficientsU3ST Arel_operator;
  u3shell::ArelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B,Arel_operator,A);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B,Arel_operator,Arel_file_name,un_u3_restrict);

  // Generate Nintr operator up to Nmax cutoff
  std::string Nintr_file_name_base=fmt::format("Nrel_{:02d}_Nmax{:02d}",A,Nmax);
  u3shell::RelativeUnitTensorCoefficientsU3ST Nrel_operator;
  u3shell::NintrRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B,Nrel_operator,A);
  std::string Nintr_file_name=fmt::format("{}.recoupler",Nintr_file_name_base);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B,Nrel_operator,Nintr_file_name,un_u3_restrict);

  // generate control file for lsu3shell
  //
  // first give specifications for unit tensors, then Brel, Arel and Nintr
  std::ofstream control_stream("relative_operators.dat");
  int num_unit=relative_unit_tensor_labels.size();
  for(int i=0; i<num_unit; ++i)
    control_stream<<fmt::format("relative_unit_{:06d}",i)<<std::endl;


  control_stream<<Brel_file_name_base<<std::endl;
  control_stream<<Arel_file_name_base<<std::endl;
  control_stream<<Nintr_file_name_base<<std::endl;

  control_stream.close();
}

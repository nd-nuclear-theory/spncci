/****************************************************************
  generate_lsu3shell_relative_tensors.cpp

  For a given valence shell, Nmax and N_step, generates input files for LSU3Shell 
  RecoupleSU3Interaction which recouples operators into format
  required by SU3RME. Input files include 

  Control file (operators.dat) giving:
    basename N lambda mu 2*S


  Operators include:
     All unit tensors defined on relative space truncated by Nmax
     Brel
     Arel
     Nrel

  Command line parameters:
    Nmax : max number of oscillator quant in relative basis
    Nstep : indicates if all or only same parity spaces are included
            in model space.  Can only take values of 1 or 2. 
    ...

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  - 8/8/16 (aem,mac): Created.
  - 11/7/16 (aem): Updated documentation.
  - 12/2/16 (aem): Add bool for U(N)->U(3) restriction on ops.
  - 2/16/17 (mac): Add Arel to output.
  - 5/29/17 (mac): Extract model space file generation.
  - 5/30/17 (mac):
    + Add U3S labels to control file output for lsu3shell.
    + Switch symplectic operators from intrinsic (A-dependent)
      to relative (A-independent).

****************************************************************/

#include <fstream>
#include <ostream>  

#include "fmt/format.h"

#include "lsu3shell/lsu3shell_operator.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"

int main(int argc, char **argv)
{
  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  u3::U3CoefInit();

  // process arguments
  if(argc<5+1)
    {
      std::cout<<"Syntax: Nmax Nstep N1v J0 T0"<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  // int Z=std::stoi(argv[1]);
  // int N=std::stoi(argv[2]);
  // int A=N+Z;
  int Nmax=std::stoi(argv[1]);
  int Nstep=std::stoi(argv[2]); // will be either 1 or 2
  int N1B=std::stoi(argv[3]);
  int J0=std::stoi(argv[4]);
  int T0=std::stoi(argv[5]);
  assert((Nstep==2)||(Nstep==1));
  int Nmin=Nmax%Nstep;
  // int T0=-1;
  // // int T0=0;
  // int J0=-1;
  bool un_u3_restrict=false;
  // if( (N==0) || (Z==0) )
  //   un_u3_restrict=true;

  // set up output stream for SU3RME control file
  std::string relative_operator_filename("relative_operators.dat");
  std::ofstream control_stream(relative_operator_filename);

  ////////////////////////////////////////////////////////////////
  // unit tensors
  ////////////////////////////////////////////////////////////////

  // Generate all relative unit tensor labels (up to Nmax cutoff)
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1B,relative_unit_tensor_labels,J0,T0,false);

  int num_unit_tensors = relative_unit_tensor_labels.size();
  std::cout
    << fmt::format("number of relative tensors: {}",num_unit_tensors)
    <<std::endl;

  // generate debugging output
  std::ofstream label_stream("relative_unit_tensor_labels.dat");
  for(int unit_tensor_index=0; unit_tensor_index<num_unit_tensors; ++unit_tensor_index)
    {
      label_stream
        << fmt::format(
            "{:06d} {}",
            unit_tensor_index,
            relative_unit_tensor_labels[unit_tensor_index].Str()
          )
        << std::endl; 
    }
  label_stream.close();

  // generate unit tensor operator files
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, relative_unit_tensor_labels,un_u3_restrict);

  // generate control file entries for lsu3shell SU3RME
  for(int unit_tensor_index=0; unit_tensor_index<num_unit_tensors; ++unit_tensor_index)
    {
      const u3shell::OperatorLabelsU3ST& operator_labels = relative_unit_tensor_labels[unit_tensor_index];
      std::string file_name_base = fmt::format("relative_unit_{:06d}",unit_tensor_index);
      control_stream
        << fmt::format(
            "{:s} {:3d} {:3d} {:3d} {:3d}",
            file_name_base,
            operator_labels.N0(),operator_labels.x0().lambda(),operator_labels.x0().mu(),
            TwiceValue(operator_labels.S0())
          )
        << std::endl;
    }

  ////////////////////////////////////////////////////////////////
  // symplectic operators
  ////////////////////////////////////////////////////////////////

  // Generate Brel operator up to Nmax cutoff
  std::string Brel_file_name_base="Brel";
  std::string Brel_file_name=fmt::format("{}.recoupler",Brel_file_name_base);
  u3shell::RelativeUnitTensorCoefficientsU3ST Brel_operator;
  u3shell::BrelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B,Brel_operator);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B,Brel_operator,Brel_file_name,un_u3_restrict);

  // Generate Arel operator up to Nmax cutoff
  std::string Arel_file_name_base="Arel";
  std::string Arel_file_name=fmt::format("{}.recoupler",Arel_file_name_base);
  u3shell::RelativeUnitTensorCoefficientsU3ST Arel_operator;
  u3shell::ArelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B,Arel_operator);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B,Arel_operator,Arel_file_name,un_u3_restrict);

  // Generate Nrel operator up to Nmax cutoff
  std::string Nrel_file_name_base="Nrel";
  u3shell::RelativeUnitTensorCoefficientsU3ST Nrel_operator;
  u3shell::NrelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B,Nrel_operator);
  std::string Nrel_file_name=fmt::format("{}.recoupler",Nrel_file_name_base);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B,Nrel_operator,Nrel_file_name,un_u3_restrict);

  // generate control file entries for lsu3shell SU3RME
  control_stream
    << fmt::format(
        "{:s} {:3d} {:3d} {:3d} {:3d}",
        Brel_file_name_base,
        -2,0,2,0
      )
    << std::endl;
  control_stream
    << fmt::format(
        "{:s} {:3d} {:3d} {:3d} {:3d}",
        Arel_file_name_base,
        2,2,0,0
      )
    << std::endl;
  control_stream
    << fmt::format(
        "{:s} {:3d} {:3d} {:3d} {:3d}",
        Nrel_file_name_base,
        0,0,0,0
      )
    << std::endl;

  // ////////////////////////////////////////////////////////////////
  // // symplectic operators -- intrinsic version (A-dependent) -- DEPRECATED
  // ////////////////////////////////////////////////////////////////
  // 
  // // Generate Brel operator up to Nmax cutoff
  // std::string Brel_file_name_base=fmt::format("Brel_{:02d}_Nmax{:02d}",A,Nmax);
  // std::string Brel_file_name=fmt::format("{}.recoupler",Brel_file_name_base);
  // u3shell::RelativeUnitTensorCoefficientsU3ST Brel_operator;
  // u3shell::BrelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B,Brel_operator,A);  // TODO remove A dependence
  // lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B,Brel_operator,Brel_file_name,un_u3_restrict);
  // 
  // // Generate Arel operator up to Nmax cutoff
  // std::string Arel_file_name_base=fmt::format("Arel_{:02d}_Nmax{:02d}",A,Nmax);
  // std::string Arel_file_name=fmt::format("{}.recoupler",Arel_file_name_base);
  // u3shell::RelativeUnitTensorCoefficientsU3ST Arel_operator;
  // u3shell::ArelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B,Arel_operator,A);  // TODO remove A dependence
  // lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B,Arel_operator,Arel_file_name,un_u3_restrict);
  // 
  // // Generate Nrel operator up to Nmax cutoff
  // std::string Nrel_file_name_base=fmt::format("Nrel_{:02d}_Nmax{:02d}",A,Nmax);
  // u3shell::RelativeUnitTensorCoefficientsU3ST Nrel_operator;
  // u3shell::NrelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B,Nrel_operator,A);  // TODO remove A dependence
  // std::string Nrel_file_name=fmt::format("{}.recoupler",Nrel_file_name_base);
  // lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B,Nrel_operator,Nrel_file_name,un_u3_restrict);
  // 
  // // generate control file entries for lsu3shell SU3RME
  // control_stream
  //   << fmt::format(
  //       "{:s} {:3d} {:3d} {:3d} {:3d}",
  //       Brel_file_name_base,
  //       -2,0,2,0
  //     )
  //   << std::endl;
  // control_stream
  //   << fmt::format(
  //       "{:s} {:3d} {:3d} {:3d} {:3d}",
  //       Arel_file_name_base,
  //       2,2,0,0
  //     )
  //   << std::endl;
  // control_stream
  //   << fmt::format(
  //       "{:s} {:3d} {:3d} {:3d} {:3d}",
  //       Nrel_file_name_base,
  //       0,0,0,0
  //     )
  //   << std::endl;

  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  control_stream.close();
}

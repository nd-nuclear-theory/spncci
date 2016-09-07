/****************************************************************
  check_two_body_unit_tensors.cpp

  Read in RMEs for two body tensors and check against expected
  values.

  Example:
    check_two_body_unit_tensors 2 1

  Input files:
    lsu3shell_basis.dat -- lsu3shell tabular basis listing file
    two_body_unit_*.rme -- output of SU3RME for each operator

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/7/16 (aem,mac): Created, starting frrom
    generate_lsu3shell_two_body_tensors.cpp.
****************************************************************/

#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_rme.h"
#include "spncci/sp_basis.h"
#include "spncci/lgi_unit_tensor.h"
#include "u3shell/unit_tensor_expansion.h"
#include "u3shell/two_body_operator.h"


////////////////////////////////////////////////////////////////
// generate known matrix elements
////////////////////////////////////////////////////////////////

  void 
  GenerateTwoBodyUnitTensorMatrices(      
    const u3shell::TwoBodyUnitTensorLabelsU3ST& two_body_unit_tensor_labels,
    const lsu3shell::U3SPNBasisLSU3Labels& basis_provenance,
    const u3shell::SpaceU3SPN& space,
    u3shell::SectorsU3SPN& sectors,
    basis::MatrixVector& matrices
    )
  // Generate matrix representation of two-body unit tensor.
  //
  // The given space must be a deuteron-like two-body space.
  {
    // set up zero initialized operator
    sectors = u3shell::SectorsU3SPN(space,u3shell::OperatorLabelsU3S(two_body_unit_tensor_labels),false);
    basis::SetOperatorToZero(sectors,matrices);

    // fill in nonzero entries
    // TODO
  }

////////////////////////////////////////////////////////////////
// read matrices
////////////////////////////////////////////////////////////////

//   void 
//   ReadTwoBodyUnitTensorMatrices(      
//     const lsu3shell::LSU3BasisTable& lsu3_basis_table,
//     const u3shell::TwoBodyUnitTensorLabelsU3ST& two_body_unit_tensor_labels
//     const u3shell::SpaceU3SPN& space,
//     u3shell::SectorsU3SPN& sectors,
//     basis::MatrixVector& matrices,
//     )
//   // Read in matrix representation of two-body unit tensor.
//   {
//     sectors = u3shell::SectorsU3SPN(space,two_body_unit_tensor_labels);
//     lsu3shell::ReadLSU3ShellRMEs(
//         is,two_body_unit_tensor_labels,lsu3_basis_table,
//         space,sectors,matrices
//       );
// 
//     is_nrel << 
// 
//   }

////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  // process arguments
  if(argc<5)
    {
      std::cout<<"Syntax: Protons Neutrons Nmin Nmax 2Nsigma0"<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  // will be either 1 or 2; 
  int Nstep=std::stoi(argv[4]);
  assert(Nstep>=1);
  assert(Nstep<=2);
  int Nmin=Nmax%Nstep;
  int two_Nsigma_0=std::stoi(argv[5]);
  HalfInt Nsigma_0 = HalfInt(two_Nsigma_0,2);

  // generate list of relative unit tensors up to Nmax cutoff
  std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST> two_body_unit_tensor_labels_list;
  u3shell::GenerateTwoBodyUnitTensorLabelsU3ST(Nmax, two_body_unit_tensor_labels_list);

  // read lsu3shell basis table and construct basis mapping
  std::string lsu3shell_basis_filename("lsu3shell_basis.dat");
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3shell_basis_filename,basis_table,basis_provenance,space);
  
  // iterate over unit tensors
  int num_unit = two_body_unit_tensor_labels_list.size();
  for(int operator_index=0; operator_index<num_unit; ++operator_index)
    {
      // extract tensor information
      const u3shell::TwoBodyUnitTensorLabelsU3ST& two_body_unit_tensor_labels
        = two_body_unit_tensor_labels_list[operator_index];
      std::cout
        << fmt::format("tensor {} labels {}",operator_index,two_body_unit_tensor_labels.Str())
        << std::endl;

      // read computed rmes
      std::string rme_filename = fmt::format("two_body_unit_{:06d}",operator_index);
      std::cout << fmt::format("reading {}",rme_filename) << std::endl;
      u3shell::SectorsU3SPN sectors(space,u3shell::OperatorLabelsU3S(two_body_unit_tensor_labels),false);
      basis::MatrixVector matrices;
      std::ifstream rme_stream(rme_filename);
      lsu3shell::ReadLSU3ShellRMEs(
          rme_stream,u3shell::OperatorLabelsU3S(two_body_unit_tensor_labels),basis_table,
        space,sectors,matrices
      );

    }

}

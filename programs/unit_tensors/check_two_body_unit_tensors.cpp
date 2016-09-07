/****************************************************************
  check_two_body_unit_tensors.cpp

  Read in RMEs for two body tensors and check against expected
  values.

  Example:
    check_two_body_unit_tensors 3 3 2 1 22     // 6Li, Nmax02

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
// generate test matrix elements
////////////////////////////////////////////////////////////////

  void 
  GenerateTwo_BodyUnitTensorMatrixVector(      
    const lsu3shell::LSU3BasisTable& lsu3_basis_table,
    const u3shell::SpaceU3SPN& space, 
    basis::MatrixVector& matrix_vector 
    );
  // Generates vector of Ncm Matrix sectors in LSU3shell basis from 
  // Nrel matrix sectors. 


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

  // read and construct basis indexing
  std::string lsu3shell_basis_filename("lsu3shell_basis.dat");
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3shell_basis_filename,basis_table,basis_provenance,space);
  
  // iterate over unit tensors
  int num_unit = two_body_unit_tensor_labels_list.size();
  for(int i=0; i<num_unit; ++i)
    {
      // extract tensor information
      const u3shell::TwoBodyUnitTensorLabelsU3ST& two_body_unit_tensor_labels
        = two_body_unit_tensor_labels_list[i];
      std::cout
        << fmt::format("tensor {} labels {}",i,two_body_unit_tensor_labels.Str())
        << std::endl;


    }

}

/****************************************************************
  check_lsu3shell_rme.cpp

  WARNING: Will not build since requires old versions of
  lsu3shell::ReadLSU3ShellRMEs and/or lgi::GenerateLGIExpansion with
  stream arguments.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  12/2/16 (aem): Created.
****************************************************************/
#include <fstream>
#include <iostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_basis.h"
#include "lsu3shell/lsu3shell_operator.h"
#include "lsu3shell/lsu3shell_rme.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"

// Given two files of operator matrix elements, compare comuputed matrix elements. 
int main(int argc, char **argv)
{
  double zero_threshold=10e-8;
  double tolerance=1e-5;
  if(argc<5)
    std::cout<<" syntax: A Nsigma_0 <basis table> <operator load file>"<<std::endl;

	int A=std::stoi(argv[1]);
  int Nsigma_0=std::stoi(argv[2]);
	std::string basis_filename=argv[3];
	std::string operators=argv[4];

  // Generate basis
  lsu3shell::LSU3ShellBasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  lsu3shell::ReadLSU3ShellBasis(Nsigma_0,basis_filename,basis_table,basis_provenance,space);
  std::cout<<"Read in basis"<<std::endl;

	int N1,lambda1,mu1,S1,T1,N2,lambda2,mu2,S2,T2;
	std::string operator1_filename,operator2_filename, log_file, type1,type2;

	// Extract operator information
	std::ifstream is(operators.c_str());
	is>>type1>>N1>>lambda1>>mu1>>S1>>T1>>operator1_filename;
	is>>type2>>N2>>lambda2>>mu2>>S2>>T2>>operator2_filename;
  is>>log_file;
	is.close();

  // Construct sectors 
  u3shell::OperatorLabelsU3ST operator1(N1,u3::SU3(lambda1,mu1),S1,T1,N1%2);
  u3shell::OperatorLabelsU3ST operator2(N2,u3::SU3(lambda2,mu2),S2,T2,N1%2);
  assert(operator1==operator2);

  u3shell::SectorsU3SPN sectors(space,operator1,true);
  std::cout<<fmt::format("There are {} sectors",sectors.size())<<std::endl;

  std::ifstream is_op1(operator1_filename.c_str());
  std::ifstream is_op2(operator2_filename.c_str());

  basis::MatrixVector matrices1;
  basis::MatrixVector matrices2;
  
  std::cout<<"Get operator 1 sectors"<<std::endl;
  if(type1==type2)  
    lsu3shell::ReadLSU3ShellRMEs(is_op1,operator1,basis_table,space,sectors,matrices1);
  else if((type1=="NREL")&&(type2=="NCM")) 
    lsu3shell::GenerateNcmMatrixVector(A,is_op1,basis_table,space,matrices1);
  else
    std::cout<<"type1=type2 or type1=NREL and type2=NCM"<<std::endl;

  std::cout<< "Get operator 2 sectors"<<std::endl;
  lsu3shell::ReadLSU3ShellRMEs(is_op2,operator2,basis_table,space,sectors,matrices2);

  std::ofstream log_stream(log_file.c_str());
  
  std::cout<<"Comparing operators"<<std::endl;
  lsu3shell::CompareLSU3ShellRMEs(log_stream, basis_provenance,space,sectors,matrices2,matrices1,tolerance,true);

}

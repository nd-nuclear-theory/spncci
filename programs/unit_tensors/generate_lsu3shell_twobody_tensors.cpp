/****************************************************************
  generate_lsu3shell_twobody_tensors.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/8/16 (aem,mac): Created.
****************************************************************/
#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell_io/lsu3shell_interface.h"
#include "spncci/sp_basis.h"
#include "spncci/lgi_unit_tensor.h"
#include "u3shell/unit_tensor_expansion.h"
#include "u3shell/two_body_operator.h"
// Given a nucleus (protons, neutrons)
// Given Nsigma_0
// Given a Nmax truncation
// Generate a list of twobody unit tensors
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
	std::ofstream control_stream("twobody_operators.dat");
	control_stream
		<<model_space
		<<std::endl;

	//Generate all relative unit tensors up to Nmax cutoff
	std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST> twobody_unit_tensor_labels;
	u3shell::GenerateTwoBodyUnitTensorLabelsU3ST(Nmax, twobody_unit_tensor_labels);
	int num_unit=twobody_unit_tensor_labels.size();
	for(int i=0; i<num_unit; ++i)
		{
			u3shell::TwoBodyUnitTensorCoefficientsU3ST twobody_unit_tensor_coefficients;
			twobody_unit_tensor_coefficients[twobody_unit_tensor_labels[i]]=1;
			lsu3shell::GenerateLSU3ShellOperators(Nmax, twobody_unit_tensor_coefficients,i);
		}

	control_stream<<fmt::format("  operators {} ",num_unit)<<std::endl;
	for(int i=0; i<num_unit; ++i)
		control_stream<<fmt::format("twobody_unit{:06d}",i)<<std::endl;
}
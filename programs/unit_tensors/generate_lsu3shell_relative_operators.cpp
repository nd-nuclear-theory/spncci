/****************************************************************
  generate_lsu3shell_relative_tensors.cpp

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
	int Nmin=std::stoi(argv[3]);
	int Nmax=std::stoi(argv[4]);
	
	int m=0;
	// Set up unit tensor basis space
	std::ofstream model_stream(fmt::format("model_space.{:02d}",m));
	model_stream<<Z<<"  "<<N<<std::endl;
	for(int Nex=Nmin; Nex<=Nmax; Nex+=2)
		model_stream<<Nex<<"  "<<-1<<std::endl;
	model_stream.close();
	++m;
	//	Set up Brel lsu3shell model spaces
	// for each Nex, separate model space file
	// indexing starts at 1

	for(int Nex=Nmin; Nex<=Nmax; Nex+=2)
			{	
				std::ofstream model_stream(fmt::format("model_space.{:02d}",m));
				model_stream<<Z<<"  "<<N<<std::endl
										<<Nex<<"  "<<-1<<std::endl;
				++m;
			}

	//begin control file
	// first give specifications for unit tensors, then Brel and Nrel
	std::ofstream control_stream("CONTROL");
	control_stream
		<<fmt::format("model_space.{:02d}  model_space.{:02d}",0,0)
		<<std::endl;

	//Generate all relative unit tensors up to Nmax cutoff
	std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
	u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensor_labels);
	lsu3shell::GenerateLSU3ShellOperators(Nmax, relative_unit_tensor_labels);

 	int num_rel=relative_unit_tensor_labels.size();
	int index=num_rel;

	control_stream<<fmt::format("  operators {} {}",0,(num_rel-1))<<std::endl;

	// Generate Brel operator up to Nmax cutoff
	u3shell::RelativeUnitTensorCoefficientsU3ST Brel_operator;
	u3shell::BrelRelativeUnitTensorExpansion(Nmin,Nmax, Brel_operator);
	lsu3shell::GenerateLSU3ShellOperators(Nmax, Brel_operator, index);
	++index;

	//Generate Nrel operator up to Nmax cutoff
  u3shell::RelativeUnitTensorCoefficientsU3ST Nrel_operator;
  u3shell::NrelRelativeUnitTensorExpansion(Nmin,Nmax, Nrel_operator);
  
  lsu3shell::GenerateLSU3ShellOperators(Nmax, Nrel_operator, index);
  ++index;

	for(int Nex=Nmin+2; Nex<=Nmax; Nex+=2)
			{	
				int index_bra=Nex/2;
				int index_ket=Nex/2+1;
				// model_stream(fmt::format("model_space.{:02d}",m).c_str());
				control_stream
					<<fmt::format("model_space.{:02d}  model_space.{:02d}",index_bra,index_ket).c_str()
					<<std::endl
					<<fmt::format("operators {} {}",num_rel,num_rel)<<std::endl;
			}


  //write list of operators and corresponding operator index to file
	std::string op_filename=fmt::format("operator_list_Nmax{:02d}",Nmax);

 	std::ofstream os(op_filename.c_str());
 	os<<Nmax<<"  "<<index<<std::endl;
 	int i;
 	for(i=0; i<num_rel; ++i)
 		{
 			os<<fmt::format("operator{:06d}.recoupler",i).c_str()<<std::endl
 				<<relative_unit_tensor_labels[i].Str()<<std::endl;
 		}
 	os<<fmt::format("operator{:06d}.recoupler",i).c_str()<<std::endl
 		<<fmt::format("Brel [{} {} {} {}]", u3::U3(0,0,-2).Str(),0,0,0)<<std::endl
 		<<fmt::format("operator{:06d}.recoupler",i+1).c_str()<<std::endl
 		<<fmt::format("Nrel [{} {} {} {}]", u3::U3(0,0,0).Str(),0,0,0)<<std::endl;

 	os.close();

}
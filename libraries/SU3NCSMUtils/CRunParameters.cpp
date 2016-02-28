#include <SU3NCSMUtils/CRunParameters.h>
#include <SU3ME/CalculateMe_proton_neutron_ncsmSU3BasisJcut.h>

#include <stdexcept>

// #define NMAX 16 
// #define NMAX 14 
#define NMAX 12 

void CRunParameters::LoadInteractionTerms(int my_rank, CInteractionPPNN& interactionPPNN, CInteractionPN& interactionPN)
{
	if (my_rank == 0) { std::cout << "Loading operators ... " << std::endl; }

	interactionPPNN.LoadTwoBodyOperators(my_rank, two_body_PPNN_); // this will load all rme tables into memory
	interactionPN.AddOperators(my_rank, two_body_PN_);
	interactionPPNN.AddOneBodyOperators(my_rank, one_body_, MECalculatorData::JJ0_, MECalculatorData::MM0_);

	if (my_rank == 0) { std::cout << "Done" << std::endl; }
}

void CRunParameters::LoadRunParameters(const char* run_params_file_name)
{
	std::ifstream input_file(run_params_file_name);
	if (!input_file)
	{
		std::ostringstream error_message;
		error_message << "Could not open file with run parameters '" << run_params_file_name << "'!";
		throw std::logic_error(error_message.str());
	}

	std::string model_space_file_name;
	input_file >> model_space_file_name;

	ncsmModelSpace_.Load(model_space_file_name);

	assert(Nmax() <= 12);

	input_file >> hw_;

	if (!input_file)
	{
		std::ostringstream error_message;
		error_message << "Error while reading hw strength. The end of file '" << run_params_file_name << "' was reached unexpectedly.";
		throw std::logic_error(error_message.str());
	}

	while (true)
	{
		std::string interaction_term;
		input_file >> interaction_term;

		if (!input_file)
		{
			return;
		}

		if (interaction_term == "AB00")
		{
			double lambda;

			input_file >> lambda;
			if (!input_file)
			{
				std::ostringstream error_message;
				error_message << "Error while reading strenght of operator [A x B]^(0 0)_{0}! End of file ???";
				throw std::logic_error(error_message.str());
			}

			AddAB00(lambda);
		}

		if (interaction_term == "VCOUL")
		{
			AddVcoul();
		}
		else if (interaction_term == "TREL")
		{
			AddTrel();
		}
		else if (interaction_term == "NCM")
		{
			double lambda;

			input_file >> lambda;
			if (!input_file)
			{
				std::ostringstream error_message;
				error_message << "Error while reading Lawson term parameter! End of file ???";
				throw std::logic_error(error_message.str());
			}

			AddNcm(lambda);
		}
		else if (interaction_term == "INT")
		{
			std::string inter_file_name_prefix;
			input_file >> inter_file_name_prefix;
			
			if (!input_file)
			{
				std::ostringstream error_message;
				error_message << "Error while reading interaction file name prefix!";
				throw std::logic_error(error_message.str());
			}

			std::string ppnn_file_name(inter_file_name_prefix);
			ppnn_file_name += ".PPNN";

			std::string pn_file_name(inter_file_name_prefix);
			pn_file_name += ".PN";

			AddTwoBodyOperatorPPNN(ppnn_file_name, 1.0);
			AddTwoBodyOperatorPN(pn_file_name, 1.0);
		}
	}
}

void CRunParameters::LoadHamiltonian(int my_rank, const std::string& hamiltonian_file_name, CInteractionPPNN& interactionPPNN, CInteractionPN& interactionPN, int A)
{
	int hw;
	std::ifstream input_file(hamiltonian_file_name);
	if (!input_file)
	{
		std::ostringstream error_message;
		error_message << "Could not open file with run parameters '" << hamiltonian_file_name << "'!";
		throw std::logic_error(error_message.str());
	}

	input_file >> hw;

	if (!input_file)
	{
		std::ostringstream error_message;
		error_message << "Error while reading hw strength. The end of file '" << hamiltonian_file_name << "' was reached unexpectedly.";
		throw std::logic_error(error_message.str());
	}

	while (true)
	{
		std::string interaction_term;
		input_file >> interaction_term;

		if (!input_file)
		{
			break;
		}

		if (interaction_term == "AB00")
		{
			double lambda;

			input_file >> lambda;
			if (!input_file)
			{
				std::ostringstream error_message;
				error_message << "Error while reading strenght of operator [A x B]^(0 0)_{0}! End of file ???";
				throw std::logic_error(error_message.str());
			}

			AddAB00(lambda);
		}

		if (interaction_term == "VCOUL")
		{
			AddVcoul(hw);
		}
		else if (interaction_term == "TREL")
		{
			AddTrel(hw, A);
		}
		else if (interaction_term == "NCM")
		{
			double lambda;

			input_file >> lambda;
			if (!input_file)
			{
				std::ostringstream error_message;
				error_message << "Error while reading Lawson term parameter! End of file ???";
				throw std::logic_error(error_message.str());
			}

			AddNcm(lambda, A);
		}
		else if (interaction_term == "INT")
		{
			std::string inter_file_name_prefix;
			input_file >> inter_file_name_prefix;
			
			if (!input_file)
			{
				std::ostringstream error_message;
				error_message << "Error while reading interaction file name prefix!";
				throw std::logic_error(error_message.str());
			}

			std::string ppnn_file_name(inter_file_name_prefix);
			ppnn_file_name += ".PPNN";

			std::string pn_file_name(inter_file_name_prefix);
			pn_file_name += ".PN";

			AddTwoBodyOperatorPPNN(ppnn_file_name, 1.0);
			AddTwoBodyOperatorPN(pn_file_name, 1.0);
		}
	}

	LoadInteractionTerms(my_rank, interactionPPNN, interactionPN);
}

//TODO: 
// Generate Nmax=0 ... 12 for Vcoul, Trel, Ncm and load the appropriate
// interaction file based on Nmax of model space
void CRunParameters::AddVcoul()
{
	std::string vcoul_file_name(su3shell_data_directory);
	
#if NMAX==12
	vcoul_file_name += "/SU3_Interactions_Operators/Vcoul_10MeV/Vcoul_2b_10MeV_Nmax12_pshell.PPNN";
#elif NMAX==14
	vcoul_file_name += "/SU3_Interactions_Operators/Vcoul_10MeV/Vcoul_2b_10MeV_Nmax14_pshell.PPNN";
#elif NMAX==16
	vcoul_file_name += "/SU3_Interactions_Operators/Vcoul_10MeV/Vcoul_2b_10MeV_Nmax16_pshell.PPNN";
#endif

	AddTwoBodyOperatorPPNN(vcoul_file_name, VcoulCoeff());
}
void CRunParameters::AddVcoul(int hw)
{
	std::string vcoul_file_name(su3shell_data_directory);
	
#if NMAX==12
	vcoul_file_name += "/SU3_Interactions_Operators/Vcoul_10MeV/Vcoul_2b_10MeV_Nmax12_pshell.PPNN";
#elif NMAX==14
	vcoul_file_name += "/SU3_Interactions_Operators/Vcoul_10MeV/Vcoul_2b_10MeV_Nmax14_pshell.PPNN";
#elif NMAX==16
	vcoul_file_name += "/SU3_Interactions_Operators/Vcoul_10MeV/Vcoul_2b_10MeV_Nmax16_pshell.PPNN";
#endif

	AddTwoBodyOperatorPPNN(vcoul_file_name, VcoulCoeff(hw));
}

void CRunParameters::AddAB00(double lambda)
{
	std::string ab00_ppnn_file_name(su3shell_data_directory);
	std::string ab00_pn_file_name(su3shell_data_directory);
	std::string ab00_1b_file_name(su3shell_data_directory);
	
	ab00_ppnn_file_name += "/SU3_Interactions_Operators/AB00/AB00_2b_Nmax8_pshell.PPNN";
	ab00_pn_file_name += "/SU3_Interactions_Operators/AB00/AB00_2b_Nmax8_pshell.PN";
	ab00_1b_file_name	+= "/SU3_Interactions_Operators/AB00/AB00_1b_Nmax8";

	AddTwoBodyOperatorPPNN(ab00_ppnn_file_name, lambda);
	AddTwoBodyOperatorPN(ab00_pn_file_name, lambda);
	AddOneBodyOperator(ab00_1b_file_name, lambda);
}

void CRunParameters::AddTrel()
{
	std::string trel_ppnn_file_name(su3shell_data_directory);
	std::string trel_pn_file_name(su3shell_data_directory);
#if NMAX==12
	trel_ppnn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax12_pshell.PPNN";
	trel_pn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax12_pshell.PN";
#elif NMAX==14
	trel_ppnn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax14_pshell.PPNN";
	trel_pn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax14_pshell.PN";
#elif NMAX==16
	trel_ppnn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax16_pshell.PPNN";
	trel_pn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax16_pshell.PN";
#endif
	AddTwoBodyOperatorPPNN(trel_ppnn_file_name, TrelCoeff());
	AddTwoBodyOperatorPN(trel_pn_file_name, TrelCoeff());
}

void CRunParameters::AddTrel(int hw, int A)
{
	std::string trel_ppnn_file_name(su3shell_data_directory);
	std::string trel_pn_file_name(su3shell_data_directory);
#if NMAX==12
	trel_ppnn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax12_pshell.PPNN";
	trel_pn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax12_pshell.PN";
#elif NMAX==14
	trel_ppnn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax14_pshell.PPNN";
	trel_pn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax14_pshell.PN";
#elif NMAX==16
	trel_ppnn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax16_pshell.PPNN";
	trel_pn_file_name += "/SU3_Interactions_Operators/Trel/Trel_2b_Nmax16_pshell.PN";
#endif
	AddTwoBodyOperatorPPNN(trel_ppnn_file_name, TrelCoeff(hw, A));
	AddTwoBodyOperatorPN(trel_pn_file_name, TrelCoeff(hw, A));
}

void CRunParameters::AddNcm(double lambda)
{
	std::string bdb_ppnn_file_name(su3shell_data_directory);
	std::string bdb_pn_file_name(su3shell_data_directory);
	std::string bdb_1b_file_name(su3shell_data_directory);
#if NMAX==12
	bdb_ppnn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax12_pshell.PPNN";
	bdb_pn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax12_pshell.PN";
	bdb_1b_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b+b_1b_Nmax12";
#elif NMAX==14
	bdb_ppnn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax14_pshell.PPNN";
	bdb_pn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax14_pshell.PN";
	bdb_1b_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b+b_1b_Nmax14";
#elif NMAX==16
	bdb_ppnn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax16_pshell.PPNN";
	bdb_pn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax16_pshell.PN";
	bdb_1b_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b+b_1b_Nmax16";
#endif
	double dcoeff = lambda*NcmCoeff();

	AddTwoBodyOperatorPPNN(bdb_ppnn_file_name, dcoeff);
	AddTwoBodyOperatorPN(bdb_pn_file_name, dcoeff);
	AddOneBodyOperator(bdb_1b_file_name, dcoeff);
}

void CRunParameters::AddNcm(double lambda, int A)
{
	std::string bdb_ppnn_file_name(su3shell_data_directory);
	std::string bdb_pn_file_name(su3shell_data_directory);
	std::string bdb_1b_file_name(su3shell_data_directory);

#if NMAX==12
	bdb_ppnn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax12_pshell.PPNN";
	bdb_pn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax12_pshell.PN";
	bdb_1b_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b+b_1b_Nmax12";
#elif NMAX==14
	bdb_ppnn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax14_pshell.PPNN";
	bdb_pn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax14_pshell.PN";
	bdb_1b_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b+b_1b_Nmax14";
#elif NMAX==16
	bdb_ppnn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax16_pshell.PPNN";
	bdb_pn_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b1+b2_Plus_b1b+2_2b_Nmax16_pshell.PN";
	bdb_1b_file_name	+= "/SU3_Interactions_Operators/Ax[(B+).(B)]/b+b_1b_Nmax16";
#endif
	double dcoeff = lambda*NcmCoeff(A);

	AddTwoBodyOperatorPPNN(bdb_ppnn_file_name, dcoeff);
	AddTwoBodyOperatorPN(bdb_pn_file_name, dcoeff);
	AddOneBodyOperator(bdb_1b_file_name, dcoeff);
}

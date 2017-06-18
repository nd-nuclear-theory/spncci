/****************************************************************
  parameters.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/parameters.h"

#include <fstream>

#include "cppformat/format.h"
#include "mcutils/parsing.h"

namespace spncci
{

  RunParameters::RunParameters(int argc, char **argv)
  {
    // read from command line arguments
    //
    // TODO reorder filenames 
    if (argc<5)
      {
        std::cout << "Syntax: A twice_Nsigma0 Nsigma0_ex_max N1v Nmax num_eigenvalues <load file>"
          // <basis filename> <Nrel filename> <Brel filename> <Arel filename>" 
                  << std::endl;
        std::exit(1);
      }
    A = std::stoi(argv[1]); 
    int twice_Nsigma0= std::stoi(argv[2]);
    Nsigma0_ex_max=std::stoi(argv[3]);
    Nsigma_0=HalfInt(twice_Nsigma0,2);
    N1v=std::stoi(argv[4]);
    Nmax = std::stoi(argv[5]);
    num_eigenvalues=std::stoi(argv[6]);
    std::string load_file=argv[7];

    // std::cout<< fmt::format("{} {} {} {} {} {}",A, twice_Nsigma0, Nsigma_0, Nsigma0_ex_max, N1v, Nmax)<<std::endl;
  
    // many-body problem
    // observable_filenames = std::vector<std::string>({"hamiltonian_u3st.dat"});

    // Reading in from load life 
    int line_count=0;
    int twice_Jmin, twice_Jmax, J_step;
    double hw_min, hw_max, hw_step;
    std::string line, observable;
    std::ifstream is(fmt::format("{}.load",load_file));
  
    assert(is);
    int J0;
    while(std::getline(is,line))
      {
        std::istringstream line_stream(line);
        ++line_count;
        if(line_count==1)
          {
            line_stream >> twice_Jmin >> twice_Jmax >> J_step;
            ParsingCheck(line_stream,line_count,line);
          }
        else if(line_count==2)
          {
            line_stream >> hw_min >> hw_max >> hw_step;
            ParsingCheck(line_stream,line_count,line);
          }
        else
          {
            line_stream >> observable >> J0;
            ParsingCheck(line_stream,line_count,line);
            observable_filenames.push_back(observable);
            observable_Jvalues.push_back(J0);
          }
      }

    num_observables = observable_filenames.size();
    observable_directory="relative_observables";
    // generate list of J values 
    HalfInt Jmin(twice_Jmin,2);
    HalfInt Jmax(twice_Jmax,2);
    for(HalfInt J=Jmin; J<=Jmax; J+=J_step)
      J_values.push_back(J);

    std::cout<<"J values are: ";
    for(auto J : J_values)
      std::cout<<J<<"  ";
    std::cout<<std::endl;

    for(double hw=hw_min; hw<=hw_max; hw+=hw_step)
      hw_values.push_back(hw);

    std::cout<<"hw values are: ";
    for(auto hw : hw_values)
      std::cout<<hw<<"  ";
    std::cout<<std::endl;


    // hard-coded directory structure and filenames
    lsu3shell_rme_directory = "lsu3shell_rme";
    lsu3shell_basis_filename = lsu3shell_rme_directory + "/" + "lsu3shell_basis.dat";
    Brel_filename = lsu3shell_rme_directory + "/" + fmt::format("Brel.rme",Nsigma0_ex_max);
    Arel_filename = lsu3shell_rme_directory + "/" + fmt::format("Arel.rme",Nsigma0_ex_max);
    Nrel_filename = lsu3shell_rme_directory + "/" + fmt::format("Nrel.rme",Nsigma0_ex_max);
    relative_unit_tensor_filename_template = lsu3shell_rme_directory + "/" + "relative_unit_{:06d}.rme";

    // hard-coded eigen solver parameters   
    eigensolver_num_convergence = 2*num_eigenvalues;    // docs for SymEigsSolver say to take "ncv>=2*nev"
    eigensolver_max_iterations = 100*num_eigenvalues;
    eigensolver_tolerance = 1e-8;
  }

}  // namespace

/****************************************************************
  generate_lgi_expansions_lsu3shell.cpp

  Generates the expansion Sp(3,R) CMF lowest grade irreps in the lsu3shell
  U(3)xSU(2) proton-neutron scheme basis by solving for the simultaneous
  null space of Bintr and Ncm. 

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  2/16/22 (aem): Created.
  
****************************************************************/
#include "lgi/lgi.h"
// #include <cstdlib>

#include "LookUpContainers/CWig9lmLookUpTable.h"
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"
#include "UNU3SU3/UNU3SU3Basics.h"

#include "lgi/lgi_gen.h"
#include "lgi/dimensions.h"
#include "utilities/nuclide.h"
#include "utilities/utilities.h"
#include "mcutils/eigen.h"
#include "mcutils/io.h"
#include "fmt/format.h"
#include "u3shell/relative_operator.h"
#include "lsu3shell/lsu3shell_basis.h"


// operator_dir = ${SPNCCI_OPERATOR_DIR}/rununittensor01/
//data/lgi_set/lgi_list.dat


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int my_rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (nprocs < 2) {
    std::cerr << "Master-slave program requires at least 2 MPI processes!" << std::endl<<std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if(argc<4+1)
    if(my_rank==0)
      {
        std::cerr<<"\nSyntax: Z N Nsigma_max <operator_dir> <optional: selected_lgi_list>\n"
        <<"operator_dir: directory containing operator files\n Nrel.PN, Nrel.PPNN, Brel.PN and Brel.PPNN for given Nsigma_max\n"
        <<std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Set global variables
  //
  // If not false, interactionPN.AddOperator will hang in cases where x0 doesn't branch to S0.
  k_dependent_tensor_strenghts=false;
  // If not false, get segmentation fault when calling function Calculate_Proton_x_Identity_MeData
  // caused by WigEckSU3SO3CG::WigEckSU3SO3CG trying to calculate coupling coefficients for
  // L0 not in x0.
  precalculate_WigEckSU3SO3CG_coefficients = false;
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialization for calculating coupling and recoupling coefficients.
  su3::init();
  CWig9lmLookUpTable<RME::DOUBLE>::AllocateMemory(true);
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Read in input variables
  int Z = std::stoi(argv[1]);
  int N = std::stoi(argv[2]);
  int Nsigma_max = std::stoi(argv[3]);
 	std::string operator_dir = argv[4];

  nuclide::NuclideType nuclide({Z,N});
  bool intrinsic = true;
  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide,intrinsic);
  unsigned int N0 = nuclide::N0ForNuclide(nuclide);
  int N1v = nuclide::ValenceShellForNuclide(nuclide);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Initial basis set-up
  //
  //Generate list of CMF LGI labels and multiplicity for given nuclide and Nmax
  lgi::MultiplicityTaggedLGIVector lgi_vector = lgi::get_lgi_vector(nuclide,Nsigma0,Nsigma_max);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  static const int tag_work = 0;
  static const int tag_finished = 1;

  // Advisor part:
  if (my_rank == 0) 
    {
      char dummy;
      MPI_Status status;

      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Create list of indices of lgi which appear in selected basis
      ////////////////////////////////////////////////////////////////////////////////////////////////
      // for(int j=0; j<lgi_vector.size(); ++j)
      //   std::cout<<j<<"  "<<lgi_vector[j].Str()<<std::endl;

      std::vector<int>lgi_indices;
      if(argc == 6)
        {
          std::string selected_lgi_filename = argv[5];
          lgi::MultiplicityTaggedLGIVector selected_lgi_vector;
          lgi::ReadLGISet(selected_lgi_filename, Nsigma0,selected_lgi_vector);
          std::unordered_set<lgi::LGI> select_lgi;
          for(const auto&[lgi,dummy] : selected_lgi_vector)
            select_lgi.insert(lgi);
          for(int i=0; i<lgi_vector.size(); ++i)
            {
              if(select_lgi.count(lgi_vector[i].irrep))
                lgi_indices.push_back(i);
            }

          // for(auto k : lgi_indices)
          //   std::cout<<"  "<<k<<"  "<<lgi_vector[k].Str()<<std::endl;
        }
      else
        {
          lgi_indices.resize(lgi_vector.size());
          std::iota(lgi_indices.begin(),lgi_indices.end(),0);
        }
      ////////////////////////////////////////////////////////////////////////////////////////////////

      // Initial distribution of work over graduate students
      int i=0;
      int num_jobs=std::min(nprocs,int(lgi_indices.size()));
      for(int grad_student=1; grad_student < num_jobs; ++grad_student)
        {
          MPI_Send(&lgi_indices[i], 1, MPI_INT, grad_student, tag_work, MPI_COMM_WORLD);
          i++;
        }

      // number of working students
      int num_students=i;

      ////////////////////////////////////////////////////////////////////////////////////////////////
      // As each grad student finish their work, give them more work
      for(i; i<lgi_indices.size(); ++i)
        {  
           // send another work
           MPI_Recv(&dummy, 0, MPI_CHAR, MPI_ANY_SOURCE, tag_finished, MPI_COMM_WORLD, &status);
           MPI_Send(&lgi_indices[i], 1, MPI_INT, status.MPI_SOURCE, tag_work, MPI_COMM_WORLD);
        }
      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Wait for all grad students to finish
      while(num_students>0)
        {
          MPI_Recv(&dummy, 0, MPI_CHAR, MPI_ANY_SOURCE, tag_finished, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          num_students--;
        }
      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Tell graduate students to get some sleep
      for(int grad_student=1; grad_student < nprocs; ++grad_student)
        MPI_Send(&dummy, 0, MPI_CHAR, grad_student, tag_finished, MPI_COMM_WORLD);
    
    }

  // Grad student part 
  else
    {
      char dummy;
      MPI_Status status;   
      
      //Do lsu3shell basis initialization 
      CBaseSU3Irreps baseSU3Irreps(Z,N,Nsigma_max);
      int lgi_index;

      while(true)
        {
          MPI_Recv(&lgi_index,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

          //If all the work finished, break out
          if(status.MPI_TAG == tag_finished) break;

          // Otherwise, compute the lgi_expansion
          basis::OperatorBlock<double> lgi_expansion
            = generate_lgi_expansion(
                nuclide,Nsigma_max,N0,
                baseSU3Irreps,
                lgi_vector,lgi_index,
                operator_dir,
                MPI_COMM_SELF
                );

          std:cout<<lgi_vector[lgi_index].Str()<<std::endl;

          // Save lgi expansion to file
          std::string lgi_expansion_filename =lgi::lgi_expansion_filename(Z,N,lgi_vector[lgi_index].irrep);
          WriteOperatorBlockBinary(lgi_expansion, lgi_expansion_filename);

          // Let advisor know expansion computed
          MPI_Send(&dummy,0,MPI_CHAR,0,tag_finished,MPI_COMM_WORLD);
        }
    }

  
  /////////////////////////////////////////////////////////////////////////////////////

  // clears memory allocated for U9, U6, and Z6 coefficients
  CWig9lmLookUpTable<RME::DOUBLE>::ReleaseMemory();
  // clear memory allocated for single-shell SU(3) rmes
  CSSTensorRMELookUpTablesContainer::ReleaseMemory();  
  // clear memory allocated for SU(3)>SO(3)
  CWigEckSU3SO3CGTablesLookUpContainer::ReleaseMemory();  
  su3::finalize();
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Finalize();
  
  return EXIT_SUCCESS;

	

}

/****************************************************************
  generate_seeds_lsu3shell.cpp

  Computes SU(3)xSU(2) reduced matrix elements of relative unit tensors
  between Sp(3,R) LGI which are used as seeds in the spncci recurrence.
  RMEs are first computed in SU(3)-coupled basis using lsu3shell code and
  then transformed to Sp(3,R)-coupled basis using lgi expansions generated
  using generate_lgi_expansion_lsu3shell.  Expansions saved to files with
  filename format:
    lgi_expansion_Z{:02d}_N{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_2Sp{:02}_2Sn{:02}_2S{:02}.dat

  Writes seeds for each sigma,sigma',partiy_bar subspace to file.

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  2/18/22 (aem): Created.
****************************************************************/
#include "LookUpContainers/CWig9lmLookUpTable.h"
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"
#include "UNU3SU3/UNU3SU3Basics.h"

#include "lgi/lgi.h"
#include "lgi/recurrence_lgi.h"
#include "spncci_basis/recurrence_indexing.h"
#include "u3ncsm/dimensions.h"
#include "u3ncsm/seed_gen.h"
#include "u3shell/relative_operator.h"
#include "utilities/nuclide.h"
#include "utilities/utilities.h"

// operator_dir = ${SPNCCI_OPERATOR_DIR}/rununittensor01/

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  int my_rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (nprocs < 2) {
    std::cerr << "Master-slave program requires at least 2 MPI processes!" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if(argc<4+1)
    if(my_rank==0)
      {
        std::cerr<<"Syntax: Z N Nsigma_max <operator_dir> <optional: selected_lgi_list>"<<std::endl;
        std::cerr<<"  operator_dir: directory containing relative unit tensor operator files ending in .PN and .PPNN"<<std::endl;
        std::cerr<<"  selected_lgi_list: optional list of Sp(3,R)SpSnS irreps to include in basis.  If none given, basis is full Nsigma_max basis"<<std::endl;
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

  // Get list of LGI in basis
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  // If filename given, read in list of lgi from file
  if(argc == 6)
    {
      std::string lgi_filename = argv[5];
      lgi::ReadLGISet(lgi_filename, Nsigma0,lgi_vector);
    }
  //Otherwise, generate LGI vector by finding possible cmf LGI by counting arguments
  else
    lgi_vector = lgi::get_lgi_vector(nuclide,Nsigma0,Nsigma_max);


  if(true && my_rank==0)
    {
      for(const auto& lgi : lgi_vector)
        std::cout<<lgi.Str()<<std::endl;
    }

  // Generate list of unit tensors
  int J0=-1;
  int T0=-1;
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nsigma_max,N1v,unit_tensor_labels,J0,T0,false);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Set up recurrence space
  spncci::spin::Space<lgi::LGI> spin_space(lgi_vector, Nsigma_max);
  spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST> spin_recurrence_space(spin_space, spin_space);
  spncci::spatial::Space spatial_space(spin_space,Nsigma0, Nsigma_max);
  spncci::spatial::RecurrenceSpace spatial_recurrence_space(spatial_space,spatial_space,N1v,Nsigma0);

  // assert(spatial_recurrence_space.size()==spin_recurrence_space.size());
  //lsu3shell basis initialization
  CBaseSU3Irreps baseSU3Irreps(Z,N,Nsigma_max);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  static const int tag_work = 0;
  static const int tag_finished = 1;
  char dummy;
  MPI_Status status;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Advisor part:
  if (my_rank == 0)
    {
      // Initial distribution of work over graduate students
      // Note, we iterate over spatial recurrence subspaces because spatial basis
      // eliminates empty subspaces, while spin does not.
      int lgi_subspace_index=0;
      int num_jobs=std::min(nprocs,int(spatial_recurrence_space.size()));
      for(int grad_student=1; grad_student < num_jobs; ++grad_student)
        {
          MPI_Send(&lgi_subspace_index, 1, MPI_INT, grad_student, tag_work, MPI_COMM_WORLD);
          lgi_subspace_index++;
        }

      // number of working students
      int num_students=lgi_subspace_index;

      // As each grad student finish their work, send out additional work until
      // all lgi expansions are accounted for
      for(lgi_subspace_index; lgi_subspace_index<spatial_recurrence_space.size(); ++lgi_subspace_index)
        {
          MPI_Recv(&dummy, 0, MPI_CHAR, MPI_ANY_SOURCE, tag_finished, MPI_COMM_WORLD, &status);
          MPI_Send(&lgi_subspace_index, 1, MPI_INT, status.MPI_SOURCE, tag_work, MPI_COMM_WORLD);
        }

      // Wait for all students to finish
      while(num_students>0)
        {
          MPI_Recv(&dummy, 0, MPI_CHAR, MPI_ANY_SOURCE, tag_finished, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          num_students--;
        }

      // Tell graduate students to get some sleep
      for(int grad_student=1; grad_student < nprocs; ++grad_student)
        MPI_Send(&dummy, 0, MPI_CHAR, grad_student, tag_finished, MPI_COMM_WORLD);
    }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Grad student part
  else
    {
      int lgi_subspace_index;
      while(true)
        {
          MPI_Recv(&lgi_subspace_index,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

          //If all the work finished, break out
          if(status.MPI_TAG == tag_finished) break;

          basis::OperatorBlock<double> recurrence_seed_block
            = spncci::seeds::GenerateRecurrenceSeedBlock(
                  nuclide,Nsigma0,N1v,
                  unit_tensor_labels,
                  spin_recurrence_space,
                  spatial_recurrence_space,
                  baseSU3Irreps,
                  operator_dir,
                  lgi_subspace_index
                );

          // Write seeds to file
          const auto&[sigma_ket,sigma_bra,parity_bar] = spatial_recurrence_space.GetSubspace(lgi_subspace_index).labels();
          std::string seed_filename = spncci::seeds::seed_filename(Z,N,Nsigma0,sigma_bra,sigma_ket,parity_bar);
          utils::WriteOperatorBlockBinary(recurrence_seed_block, seed_filename);

          // Let advisor know seeds for given lgi pair computed
          MPI_Send(&dummy,0,MPI_CHAR,0,tag_finished,MPI_COMM_WORLD);
        }
    }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

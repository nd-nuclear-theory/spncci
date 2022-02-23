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
#include "lgi/lgi.h"

#include "LookUpContainers/CWig9lmLookUpTable.h"
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"
#include "UNU3SU3/UNU3SU3Basics.h"

#include "lgi/lgi_gen.h"
#include "lgi/dimensions.h"
#include "utilities/nuclide.h"
#include "mcutils/eigen.h"
#include "utilities/utilities.h"

#include "u3shell/relative_operator.h"
#include "spncci/recurrence_seeds.h"
#include "lgi/su3rme.h"

// operator_dir = ${SPNCCI_OPERATOR_DIR}/rununittensor01/
namespace spncci{
  namespace seeds{

    basis::OperatorBlock<double> GenerateRecurrenceSeedBlock(
      const nuclide::NuclideType& nuclide,
      const HalfInt& Nsigma0,
      const int N1v,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
      const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST>& spin_recurrence_space,
      const spncci::spatial::RecurrenceSpace& spatial_recurrence_space,
      const CBaseSU3Irreps& baseSU3Irreps,
      const std::string& operator_dir,
      const int recurrence_sp3r_space_index
      )
    {
      basis::OperatorBlock<double> recurrence_seed_block;

      MPI_Comm individual_comm = MPI_COMM_SELF;
      const auto&[Z,N]=nuclide;

      // Get relevant spin::RecurrenceLGISpace space and extract
      // sigma, sigma' and operator parity information
      ////////////////////////////////////////////////////////////////////////////////////
      // spatial::RecurrenceSpace() []
      // ->spatial::RecurrenceSp3RSpace() [sigma,sigma',parity_bar]
      //   ->spatial::RecurrenceNnsumSpace() [Nsum]
      //     ->spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
      //       ->spatial::RecurrenceOperatorSubspace() [x0] ->rho0_max
      //         ->spatial::RecurrenceOperatorState() [Nbar,Nbar']
      //
      // Get spatial operator subspace for given sigma_ket,sigma_bra
      // Note: for seeds we only care about the Nsum=0 subspace (with index=0) which has a single
      // spatial::RecurrenceU3Space with omega=sigma and omega' = sigma' which also has index=0
      ///////////////////////////////////////////////////////////////////////////////////
      const auto& spatial_recurrence_sp3r_space = spatial_recurrence_space.GetSubspace(recurrence_sp3r_space_index);
      const auto&[sigma_ket,sigma_bra,parity_bar] = spatial_recurrence_sp3r_space.labels();
      const spncci::spatial::RecurrenceU3Space& spatial_recurrence_u3s_space = spatial_recurrence_sp3r_space.GetSubspace(0).GetSubspace(0);
      int Nex_ket(sigma_ket.N()-Nsigma0);
      int Nex_bra(sigma_bra.N()-Nsigma0);

      int exchange_symm_bar = (parity_bar+1)%2;
      int spin_recurrence_space_index = spin_recurrence_space.LookUpSubspaceIndex({sigma_ket,sigma_bra,exchange_symm_bar});
      const auto& spin_recurrence_lgi_space = spin_recurrence_space.GetSubspace(spin_recurrence_space_index);

      // const auto&[sigma_ket,sigma_bra,exchange_symm_bar] = spin_recurrence_lgi_space.labels();

      // // If there is no corresponding
      if(spin_recurrence_space_index == -1)
        {
          std::cout<<"No corresponding spin space for "<<sigma_ket.Str()<<" "<<sigma_bra.Str()<<" "<<parity_bar<<std::endl;
          return recurrence_seed_block;
        }



      // std::cout<<"****************************************"<<std::endl;
      // std::cout<<sigma_bra.Str()<<"  "<<sigma_ket.Str()<<"  "<<parity_bar<<std::endl;
      // std::cout<<"****************************************"<<std::endl;



      // Allocate seed block
      recurrence_seed_block
        = basis::OperatorBlock<double>::Zero(
            spatial_recurrence_u3s_space.dimension(),
            spin_recurrence_lgi_space.dimension()
          );
      // std::cout<<"recurrence seed block "<<recurrence_seed_block.rows()<<" x "<<recurrence_seed_block.cols()<<std::endl;
      // Loop over all unit tensors
      //TODO: openMP parallelize over operator index.
      for(int operator_index=0; operator_index<unit_tensor_labels.size(); ++operator_index)
        {
          // Extract operator labels and apply spatial selection rules
          const auto& [operator_labels,relative_bra, relative_ket]=unit_tensor_labels[operator_index].Key();
          const auto& [dN0,x0,S0,T0,g0] = operator_labels;
          if(dN0+sigma_ket.N()!=sigma_bra.N()) continue;

          int spatial_operator_index = spatial_recurrence_u3s_space.LookUpSubspaceIndex(x0);
          if(spatial_operator_index == -1) continue;
          const auto& spatial_recurrence_operator_subspace = spatial_recurrence_u3s_space.GetSubspace(spatial_operator_index);
          // std::cout<<"spatial_operator_index "<<spatial_operator_index<<std::endl;

          const auto& [Nbar,Sbar,Tbar] = relative_ket;
          const auto& [Nbarp,Sbarp,Tbarp] = relative_bra;
          int spatial_recurrence_operator_state_index = spatial_recurrence_operator_subspace.LookUpStateIndex({Nbar,Nbarp});
          if(spatial_recurrence_operator_state_index == -1) continue;
          // std::cout<<"spatial_recurrence_operator_state_index "<<spatial_recurrence_operator_state_index<<std::endl;
          // Look up spatial and spin operator subspace indices

          int rhot_max=u3::OuterMultiplicity(sigma_ket.SU3(),x0,sigma_bra.SU3());
          // assert(rhot_max==spatial_recurrence_operator_subspace.degeneracy());
          // if(rhot_max<1) continue;

          // Definite labels for use in generate seeds
          auto spin_tensor_labels = spncci::spin::UnitTensorLabelsST(int(S0),int(T0),int(Sbar),int(Sbarp),int(Tbar),int(Tbarp));
          SU3xSU2::LABELS w0(1,x0.lambda(),x0.mu(),TwiceValue(S0));

          // Read in operator from file
          std::string operator_base_name = fmt::format("{}/relative_unit_{:06}",operator_dir,operator_index);
          std::ofstream interaction_log_file("/dev/null");
          bool log_is_on = false;
          bool generate_missing_rme = true;

          // The constructor for interactionPPNN must always be called before the constructor for
          // interactionPN otherwise a global variable doesn't get correctly allocated...F#$%*$
          CInteractionPPNN interactionPPNN(baseSU3Irreps,log_is_on,interaction_log_file);
          CInteractionPN interactionPN(baseSU3Irreps,generate_missing_rme,log_is_on,interaction_log_file);

          interactionPPNN.LoadTwoBodyOperator(operator_base_name+".PPNN");
          interactionPN.AddOperator(operator_base_name+".PN");
          interactionPPNN.TransformTensorStrengthsIntoPP_NN_structure();
          ////////////////////////////////////////////////////////////////////////////////////
          // Iterating through spin space
          ////////////////////////////////////////////////////////////////////////////////////
          // spin::RecurrenceSpace() []
          //   spin::RecurrenceLGISpace() [sigma,sigma',exchange_symm_bar]
          //   ->spin::RecurrenceSpinSpace() [S,S']
          //     -> spin::RecurrenceSpinSubspace() [Sp,Sn,Sp',Sn']/[T,T']->(gamma,gamma')
          //       ->spin::RecurrenceOperatorState() [operator_index]
          ////////////////////////////////////////////////////////////////////////////////////
          for(int spin_space_index=0; spin_space_index<spin_recurrence_lgi_space.size(); ++spin_space_index)
            {
              const auto& spin_recurrence_spin_space = spin_recurrence_lgi_space.GetSubspace(spin_space_index);
              const auto& [S_ket, S_bra] = spin_recurrence_spin_space.labels();
              if(!am::AllowedTriangle(S_ket,S0,S_bra)) continue;

              for(int spin_subspace_index=0; spin_subspace_index<spin_recurrence_spin_space.size(); ++spin_subspace_index)
                {
                  const auto& spin_recurrence_subspace = spin_recurrence_spin_space.GetSubspace(spin_subspace_index);
                  const auto& ket_upstream_labels = spin_recurrence_subspace.ket_upstream_labels();
                  const auto& bra_upstream_labels = spin_recurrence_subspace.bra_upstream_labels();
                  if(!spncci::spin::UnitTensorAllowed(ket_upstream_labels,spin_tensor_labels,bra_upstream_labels))
                    continue;

                  const auto [Sp_ket, Sn_ket] = ket_upstream_labels;
                  const auto [Sp_bra, Sn_bra] = bra_upstream_labels;

                  // Look up spin operator index
                  int spin_operator_index = spin_recurrence_subspace.LookUpStateIndex(spin_tensor_labels);
                  // std::cout<<"spin_operator_index "<<spin_operator_index<<std::endl;
                  int jdiag=0;
                  int ndiag=1;

                  // Define bra and ket basis
                  proton_neutron::ModelSpaceMapType model_space_map_ket;
                  model_space_map_ket[Nex_ket][{TwiceValue(Sp_ket),TwiceValue(Sn_ket)}][TwiceValue(S_ket)].emplace_back(1,sigma_ket.SU3().lambda(),sigma_ket.SU3().mu());
                  proton_neutron::ModelSpace ket_ncsmModelSpace(Z,N,model_space_map_ket);
                  lsu3::CncsmSU3xSU2Basis ket_ncsm_basis;
                  ket_ncsm_basis.ConstructBasis(ket_ncsmModelSpace, jdiag, ndiag, individual_comm);
                  // int dim_ket = lsu3shell::get_num_U3PNSPN_irreps(ket_ncsm_basis);

                  // Construct bra basis from model space
                  proton_neutron::ModelSpaceMapType model_space_map_bra;
                  model_space_map_bra[Nex_bra][{TwiceValue(Sp_bra),TwiceValue(Sn_bra)}][TwiceValue(S_bra)].emplace_back(1,sigma_bra.SU3().lambda(),sigma_bra.SU3().mu());
                  proton_neutron::ModelSpace bra_ncsmModelSpace(Z,N,model_space_map_bra);
                  lsu3::CncsmSU3xSU2Basis bra_ncsm_basis;
                  bra_ncsm_basis.ConstructBasis(bra_ncsmModelSpace, jdiag, ndiag, individual_comm);
                  // int dim_bra = lsu3shell::get_num_U3PNSPN_irreps(bra_ncsm_basis);

                  // std::cout<<"Calculate unit tensor RMEs in lsu3shell basis"<<std::endl;
                  basis::OperatorBlocks<double> su3_unit_tensor_rmes
                    =lsu3shell::CalculateRME(interactionPPNN,interactionPN,bra_ncsm_basis,ket_ncsm_basis,dN0,w0,rhot_max);

                  // Read in lgi expansion from file
                  std::string lgi_expansion_filename_bra =lgi::lgi_expansion_filename(Z,N,{{sigma_bra,Sp_bra,Sn_bra,S_bra},Nex_bra});
                  std::string lgi_expansion_filename_ket =lgi::lgi_expansion_filename(Z,N,{{sigma_ket,Sp_ket,Sn_ket,S_ket},Nex_ket});
                  basis::OperatorBlock<double> lgi_expansion_bra = ReadOperatorBlockBinary(lgi_expansion_filename_bra);
                  basis::OperatorBlock<double> lgi_expansion_ket = ReadOperatorBlockBinary(lgi_expansion_filename_ket);

                  // Loop over outer multiplicity and transform each block from lsu3shell basis to spncci basis.
                  // Save seeds into recurrence seed matrix
                  for(int irhot=0; irhot<rhot_max; ++irhot)
                    {
                      int rho0=irhot+1;
                      // Row index in recurrence seed matrix
                      int target_index_row
                        = spatial_recurrence_u3s_space.GetSubspaceOffset(spatial_operator_index,rho0)
                          + spatial_recurrence_operator_state_index;


                      basis::OperatorBlock<double> block
                        = lgi_expansion_bra.transpose()*su3_unit_tensor_rmes[irhot]*lgi_expansion_ket;

                      // // TEMP: write block to file for comparison with old seeds.
                      // std::string temp_filename
                      //   = fmt::format("temp/seeds_Z{:02d}_N{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_2Sp{:02d}_2Sn{:02d}_2S{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_2Sp{:02d}_2Sn{:02d}_2S{:02d}_rho{}_i{}.dat",
                      //         Z,N,Nex_bra,sigma_bra.SU3().lambda(),sigma_bra.SU3().mu(),TwiceValue(Sp_bra),TwiceValue(Sn_bra),TwiceValue(S_bra),
                      //         Nex_ket,sigma_ket.SU3().lambda(),sigma_ket.SU3().mu(),TwiceValue(Sp_ket),TwiceValue(Sn_ket),TwiceValue(S_ket),
                      //         rho0,operator_index
                      //       );

                      // WriteOperatorBlockBinary(block,temp_filename);

                      for(int i=0; i<block.rows(); i++)
                        for(int j=0; j<block.cols(); j++)
                          {
                            int gamma_bra=i+1;
                            int gamma_ket=j+1;

                            // Column index in recurrence seed matrix
                            int target_index_col
                              = spin_recurrence_lgi_space.GetSubspaceOffset(spin_space_index)
                                +spin_recurrence_spin_space.GetSubspaceOffset(spin_subspace_index,gamma_ket,gamma_bra)
                                + spin_operator_index;

                            recurrence_seed_block(target_index_row,target_index_col) = block(i,j);

                          }
                    }
                }
            }
        }
      return recurrence_seed_block;
    }

  }// namespaces
}




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
          WriteOperatorBlockBinary(recurrence_seed_block, seed_filename);

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

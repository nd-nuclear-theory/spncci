/****************************************************************
  spncci_seeds.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/

#include "spncci/spncci_seeds.h"
#include <mpi.h>
#include "fmt/format.h"
#include "basis/operator.h"
#include "lgi/su3rme.h"
#include "lgi/lgi.h"
#include "lgi/dimensions.h"
#include "spncci/eigenproblem.h"

#include "LSU3/ncsmSU3xSU2Basis.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "UNU3SU3/UNU3SU3Basics.h"
#include "LookUpContainers/CWig9lmLookUpTable.h"

namespace spncci
{
  namespace seeds
  {

    proton_neutron::ModelSpaceMapType generate_bra_model_space_Brel(
      const lgi::LGI& lgi,
      const lgi::MultiplicityTaggedLGIVector& lgi_vector
    )
      {
        const auto&[Nex,sigma,Sp,Sn,S]=lgi.Key();
        // bra model space consists subspaces connected by B[-2(0,2)]
        proton_neutron::ModelSpaceMapType model_space_map_bra;
        for(const auto&[lgi_p,dim_p] : lgi_vector)
          {
            const auto&[Nex_p,sigma_p,Sp_p,Sn_p,S_p]=lgi_p.Key();
            
            // Apply selection rule 
            if(Nex-2!=Nex_p){continue;}
            if(Sp!=Sp_p){continue;}
            if(Sn!=Sn_p){continue;}
            if(S!=S_p){continue;}          
            if(u3::OuterMultiplicity(sigma.SU3(),u3::SU3(0,2),sigma_p.SU3())==0) {continue;}
            
            //If passed all selection rules, add to bra spaces 
            SU3::LABELS w(1,sigma_p.SU3().lambda(),sigma_p.SU3().mu());
            model_space_map_bra[Nex_p][{TwiceValue(Sp_p),TwiceValue(Sn_p)}][TwiceValue(S_p)].push_back(w);
          }
        return model_space_map_bra;
      }


    void generate_lgi_expansion(
      const nuclide::NuclideType& nuclide,
      const int Nsigma_max,
      const std::string& operator_dir,
      int eigensolver_max_iterations,
      double eigensolver_tolerance,
      const MPI_Comm world_comm //set default to be MPI_COMM_WORLD
      )
    {
      // Basis set up
      bool intrinsic = true;
      HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide,intrinsic);
      unsigned int N0 = nuclide::N0ForNuclide(nuclide);
      int N1v = nuclide::ValenceShellForNuclide(nuclide);
      const auto&[Z,N]=nuclide;

      //Generate LGI vector by finding possible cmf LGI by counting arguments 
      lgi::MultiplicityTaggedLGIVector lgi_vector = lgi::get_lgi_vector(nuclide,Nsigma0,Nsigma_max);
      for(const auto& lgi : lgi_vector)
        std::cout<<lgi.Str()<<std::endl;

      // Get number of cmf U3SPN irreps for each subspace
      std::map<u3shell::U3SPN,unsigned int> cmf_basis_dimensions 
        = lsu3shell::lsu3shell_cmf_basis_dimensions(nuclide,Nsigma0,Nsigma_max);
      ////////////////////////////////////////////////////////////////////////////////////////////////
      
      su3::init();
      CWig9lmLookUpTable<RME::DOUBLE>::AllocateMemory(true);

      CBaseSU3Irreps baseSU3Irreps(Z,N,Nsigma_max);

      std::ofstream interaction_log_file("/dev/null");
      bool log_is_on = false;
      bool generate_missing_rme = true;

      std::string brel_operator_base_name = fmt::format("{}/Brel",operator_dir);
      std::string nrel_operator_base_name = fmt::format("{}/Nrel",operator_dir);

      //////////////////////////////////////////////////////////////////////////////////////////
      // Load Nrel operator used to calculate Ncm
      //////////////////////////////////////////////////////////////////////////////////////////
      CInteractionPPNN interactionPPNN_Nrel(baseSU3Irreps,log_is_on,interaction_log_file);
      CInteractionPN interactionPN_Nrel(baseSU3Irreps,generate_missing_rme,log_is_on,interaction_log_file);
      interactionPPNN_Nrel.LoadTwoBodyOperator(nrel_operator_base_name+".PPNN");
      interactionPN_Nrel.AddOperator(nrel_operator_base_name+".PN");
      interactionPPNN_Nrel.TransformTensorStrengthsIntoPP_NN_structure();
      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Calculate rmes of Brel and Ncm for each U3SpSnS subspace corresponding to an LGI
      basis::OperatorBlocks<double> cmf_basis_transformations(lgi_vector.size());
      std::cout<<"Calculate Ncm eigenvectors "<<std::endl;
      for(int l=0; l<lgi_vector.size(); ++l)
        {  
          const auto&[lgi,gamma_max] = lgi_vector[l];
          //Get labels
          const auto&[Nex,sigma,Sp,Sn,S]=lgi.Key();
          
          if(Nex==0)
            continue;
     
          // Note: lsu3shell basis construction has internal call to MPI, 
          // so we split ranks into invidual communicators, so that each 
          // MPI rank is constructing only the basis for it's only set of LGI
          //
          // Get the rank in the original communicator
          int my_rank;
          MPI_Comm_rank(world_comm, &my_rank);
          // Split the communicators such that each rank
          // belongs to its own little world
          MPI_Comm individual_comm;
          MPI_Comm_split(world_comm, my_rank, my_rank, &individual_comm);

          // With communicator split, ndiag is always 1
          int jdiag=0; 
          int ndiag=1;

          ////////////////////////////////////////////////////////////////////////////////////////////////
          // Calculate rmes for Nrel
          ////////////////////////////////////////////////////////////////////////////////////////////////
          // Temporary container for constructing the bra and ket model space 
          // model_space_map organized into {N : {SpSn : {S : <(lambda,mu)>}}} 
          proton_neutron::ModelSpaceMapType model_space_map;
          
          // Model space consists of irreps with quantum numbers Nex (lambda,mu) Sp Sn S
          // Nrel has quantum numbers 0(0,0)0 0, so bra and ket model space are the same 
          model_space_map[Nex][{TwiceValue(Sp),TwiceValue(Sn)}][TwiceValue(S)].emplace_back(1,sigma.SU3().lambda(),sigma.SU3().mu());
          proton_neutron::ModelSpace ket_ncsmModelSpace(Z,N,model_space_map);
          proton_neutron::ModelSpace bra_ncsmModelSpace(Z,N,model_space_map);
          
          // Construct ket basis from model space
          lsu3::CncsmSU3xSU2Basis ket_ncsm_basis;
          ket_ncsm_basis.ConstructBasis(ket_ncsmModelSpace, jdiag, ndiag, individual_comm);

          // Construct bra basis from model space
          lsu3::CncsmSU3xSU2Basis bra_ncsm_basis;
          bra_ncsm_basis.ConstructBasis(bra_ncsmModelSpace, jdiag, ndiag, individual_comm);


          int dim = lsu3shell::get_num_U3PNSPN_irreps(ket_ncsm_basis);
          basis::OperatorBlock<double> Ncm_matrix = double(N0+Nex)*Eigen::MatrixXd::Identity(dim,dim); 
          int dN0=0;
          // w0(rho0,lambda0,mu0,twice_S0)
          SU3xSU2::LABELS w0(1,0,0,0);
          int rhot_max=1;        
          // There is only one block in the operator_blocks returned by CalculateRME because rhot_max=1;  
          double A=Z+N;
          Ncm_matrix-=2./A*lsu3shell::CalculateRME(interactionPPNN_Nrel,interactionPN_Nrel,bra_ncsm_basis,ket_ncsm_basis,dN0,w0,rhot_max)[0];

          //Expected number of eigenvectors will eigenvalue 0
          int num_cmf_irreps = cmf_basis_dimensions[lgi.u3spn()];
          
          bool verbose=false;
          spncci::Vector eigenvalues;
          auto& eigenvectors=cmf_basis_transformations[l];
          spncci::SolveEigenproblem(
            Ncm_matrix,num_cmf_irreps,2*num_cmf_irreps,
            eigensolver_max_iterations,eigensolver_tolerance,
            eigenvalues,eigenvectors,verbose
          );
        }

      //////////////////////////////////////////////////////////////////////////////////////////
      // Load Brel used to calculated Bintr
      //////////////////////////////////////////////////////////////////////////////////////////
      CInteractionPPNN interactionPPNN_Brel(baseSU3Irreps,log_is_on,interaction_log_file);
      CInteractionPN interactionPN_Brel(baseSU3Irreps,generate_missing_rme,log_is_on,interaction_log_file);
      interactionPPNN_Brel.LoadTwoBodyOperator(brel_operator_base_name+".PPNN");
      interactionPN_Brel.AddOperator(brel_operator_base_name+".PN");
      interactionPPNN_Brel.TransformTensorStrengthsIntoPP_NN_structure();
      //////////////////////////////////////////////////////////////////////////////////////////   
      // Calculate Brel matrix and get null vectors 
      /////////////////////////////////////////////////////////////////////////////////////
      // For each LGI, calculate Brel matrix, apply Nrel eigenvector transformation
      
      for(const auto& [lgi,gamma_max] : lgi_vector)
      {
        const auto&[Nex,sigma,Sp,Sn,S]=lgi.Key();
        if (Nex==0)
          continue;
        proton_neutron::ModelSpaceMapType model_space_map_ket;
        model_space_map_ket[Nex][{TwiceValue(Sp),TwiceValue(Sn)}][TwiceValue(S)].emplace_back(1,sigma.SU3().lambda(),sigma.SU3().mu());
        proton_neutron::ModelSpaceMapType model_space_map_bra=generate_bra_model_space_Brel(lgi,lgi_vector);
        std::cout<<lgi.Str()<<std::endl;
        if (model_space_map_bra.size()==0)
          continue;
        std::cout<<"-------------------------"<<std::endl;
        for(const auto& [Nexp,temp1] : model_space_map_bra)
          for(const auto& [SpSn,temp2]: temp1)
            for(const auto& [S,temp3]: temp2)
              for(const auto& [r,lam,mu] : temp3)
                std::cout<<fmt::format("{}({},{}) {} {} {}",Nexp,lam,mu,SpSn.first,SpSn.second,S)<<std::endl;
        std::cout<<"*****************************"<<std::endl;
        proton_neutron::ModelSpace ket_ncsmModelSpace(Z,N,model_space_map_ket);
        proton_neutron::ModelSpace bra_ncsmModelSpace(Z,N,model_space_map_bra);

        // Note: lsu3shell basis construction has internal call to MPI, 
        // so we split ranks into invidual communicators, so that each 
        // MPI rank is constructing only the basis for it's only set of LGI
        //
        // Get the rank in the original communicator
        int my_rank;
        MPI_Comm_rank(world_comm, &my_rank);
        // Split the communicators such that each rank
        // belongs to its own little world
        MPI_Comm individual_comm;
        MPI_Comm_split(world_comm, my_rank, my_rank, &individual_comm);

        // With communicator split, ndiag is always 1
        int jdiag=0; 
        int ndiag=1;

        // Construct ket basis from model space
        lsu3::CncsmSU3xSU2Basis ket_ncsm_basis;
        ket_ncsm_basis.ConstructBasis(ket_ncsmModelSpace, jdiag, ndiag, individual_comm);

        // Construct bra basis from model space
        lsu3::CncsmSU3xSU2Basis bra_ncsm_basis;
        bra_ncsm_basis.ConstructBasis(bra_ncsmModelSpace, jdiag, ndiag, individual_comm);

        
          // int dim_ket = lsu3shell::get_num_U3PNSPN_irreps(ket_ncsm_basis);
          // int dim_bra = lsu3shell::get_num_U3PNSPN_irreps(ket_ncsm_basis);
          // basis::OperatorBlock<double> Ncm_matrix = double(N0+Nex)*Eigen::MatrixXd::Identity(dim,dim); 
        int dN0=-2;
        // w0(rho0,lambda0,mu0,twice_S0)
        SU3xSU2::LABELS w0(1,0,2,0);
        int rhot_max=1;  
        // since rhot_max=1, only one element of OperatorBlocks      
        auto Brel_matrix=lsu3shell::CalculateRME(interactionPPNN_Brel,interactionPN_Brel,bra_ncsm_basis,ket_ncsm_basis,dN0,w0,rhot_max)[0];
        std::cout<<Brel_matrix<<std::endl;
      }
      /////////////////////////////////////////////////////////////////////////////////////

      // clears memory allocated for U9, U6, and Z6 coefficients
      CWig9lmLookUpTable<RME::DOUBLE>::ReleaseMemory();
      // clear memory allocated for single-shell SU(3) rmes
      CSSTensorRMELookUpTablesContainer::ReleaseMemory();  
      // clear memory allocated for SU(3)>SO(3)
      CWigEckSU3SO3CGTablesLookUpContainer::ReleaseMemory();  

      su3::finalize();

    }

  }// seeds namespace 

}  // namespace

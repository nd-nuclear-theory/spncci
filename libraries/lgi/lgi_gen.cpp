/****************************************************************
  lgi_gen.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "lgi/lgi_gen.h"

#include <fstream>
#include <iostream>
#include <omp.h>

#include "LookUpContainers/CWig9lmLookUpTable.h"
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "UNU3SU3/UNU3SU3Basics.h"

#include "lsu3shell/lsu3shell_basis.h"
#include "lgi/dimensions.h"
#include "su3rme.h"
#include "lgi/null_solver.h"

namespace lgi
{
    proton_neutron::ModelSpaceMapType 
    generate_bra_model_space_Brel(
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



  basis::OperatorBlock<double> NcmMatrix(
    const int A, const double Ntot,
    const CBaseSU3Irreps& baseSU3Irreps,
    const lsu3::CncsmSU3xSU2Basis& bra_ncsm_basis,
    const lsu3::CncsmSU3xSU2Basis& ket_ncsm_basis,
    const std::string& nrel_operator_base_name
    )
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Calculate Ncm matrix 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  {
    //////////////////////////////////////////////////////////////////////////////////////////
    // Load Nrel operator used to calculate Ncm
    //////////////////////////////////////////////////////////////////////////////////////////
    // std::string nrel_operator_base_name = fmt::format("{}/Nrel",operator_dir);

    std::ofstream interaction_log_file("/dev/null");
    bool log_is_on = false;
    bool generate_missing_rme = true;

    CInteractionPPNN interactionPPNN_Nrel(baseSU3Irreps,log_is_on,interaction_log_file);
    CInteractionPN interactionPN_Nrel(baseSU3Irreps,generate_missing_rme,log_is_on,interaction_log_file);
    interactionPPNN_Nrel.LoadTwoBodyOperator(nrel_operator_base_name+".PPNN");
    interactionPN_Nrel.AddOperator(nrel_operator_base_name+".PN");
    interactionPPNN_Nrel.TransformTensorStrengthsIntoPP_NN_structure();



    // w0(rho0,lambda0,mu0,twice_S0)
    SU3xSU2::LABELS w0(1,0,0,0);
    int dN0=0;
    int rhot_max=1;        
    
    // There is only one block in the operator_blocks returned by CalculateRME because rhot_max=1;
    basis::OperatorBlock<double> Ncm 
      = -2./A*lsu3shell::CalculateRME(interactionPPNN_Nrel,interactionPN_Nrel,bra_ncsm_basis,ket_ncsm_basis,dN0,w0,rhot_max)[0];
    
    int dim=Ncm.rows();
    Ncm+=Ntot*Eigen::MatrixXd::Identity(dim,dim); 
    
    return Ncm;
  }


  basis::OperatorBlock<double> BMatrix(
    const CBaseSU3Irreps& baseSU3Irreps,
    const lsu3::CncsmSU3xSU2Basis& bra_ncsm_basis,
    const lsu3::CncsmSU3xSU2Basis& ket_ncsm_basis,
    const std::string& brel_operator_base_name
    )
     //////////////////////////////////////////////////////////////////////////////////////////   
    // Calculate Brel matrix 
    /////////////////////////////////////////////////////////////////////////////////////
    {

      std::ofstream interaction_log_file("/dev/null");
      bool log_is_on = false;
      bool generate_missing_rme = true;

      CInteractionPPNN interactionPPNN_Brel(baseSU3Irreps,log_is_on,interaction_log_file);
      CInteractionPN interactionPN_Brel(baseSU3Irreps,generate_missing_rme,log_is_on,interaction_log_file);

      interactionPPNN_Brel.LoadTwoBodyOperator(brel_operator_base_name+".PPNN");
      interactionPN_Brel.AddOperator(brel_operator_base_name+".PN");
      interactionPPNN_Brel.TransformTensorStrengthsIntoPP_NN_structure();

      int dN0=-2;
      // w0(rho0,lambda0,mu0,twice_S0)
      SU3xSU2::LABELS w0(1,0,2,0);
      int rhot_max=1;  

      return lsu3shell::CalculateRME(interactionPPNN_Brel,interactionPN_Brel,bra_ncsm_basis,ket_ncsm_basis,dN0,w0,rhot_max)[0];
     
    }


  basis::OperatorBlocks<double> generate_lgi_expansion(
    const nuclide::NuclideType& nuclide,
    const int Nsigma_max,
    const lgi::MultiplicityTaggedLGIVector& lgi_vector,
    const std::string& operator_dir,
    const MPI_Comm world_comm //set default to be MPI_COMM_WORLD
    )
  {

    unsigned int N0 = nuclide::N0ForNuclide(nuclide);
    
    su3::init();
    CWig9lmLookUpTable<RME::DOUBLE>::AllocateMemory(true);

    // Define distributions of partices over shells and do basic basis setup
    const auto&[Z,N]=nuclide;
    CBaseSU3Irreps baseSU3Irreps(Z,N,Nsigma_max);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Calculate rmes of B and Ncm for each U3SpSnS subspace corresponding to an LGI
    // Then solve for simultaneous null space of B and Ncm
    basis::OperatorBlocks<double> lgi_expansions(lgi_vector.size());
    for(int l=0; l<lgi_vector.size(); ++l)
      {  
        //Extract label information
        const auto&[lgi,gamma_max] = lgi_vector[l];
        const auto&[Nex,sigma,Sp,Sn,S]=lgi.Key();
        
        // If Nex==0, then all of the irreps in the subspace are lgi and the basis
        // for the null spaces is given by the identity matrix        
        if(Nex==0)
        {
          lgi_expansions[l]=Eigen::MatrixXd::Identity(gamma_max,gamma_max);
          continue;
        }
   
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Construct the basis 
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

        // Model space consists of irreps with quantum numbers Nex (lambda,mu) Sp Sn S
        // map used for defining model_space_map organized into {N : {SpSn : {S : <(lambda,mu)>}}} 
        proton_neutron::ModelSpaceMapType model_space_map_ket;
        model_space_map_ket[Nex][{TwiceValue(Sp),TwiceValue(Sn)}][TwiceValue(S)].emplace_back(1,sigma.SU3().lambda(),sigma.SU3().mu());
        proton_neutron::ModelSpace ket_ncsmModelSpace(Z,N,model_space_map_ket);
        
        // Construct ket basis from model space
        // num in communicator must be equal to ndiag
        lsu3::CncsmSU3xSU2Basis ket_ncsm_basis;
        ket_ncsm_basis.ConstructBasis(ket_ncsmModelSpace, jdiag, ndiag, individual_comm);
        int dim_ket = lsu3shell::get_num_U3PNSPN_irreps(ket_ncsm_basis);
        
        // Nrel has quantum numbers 0(0,0)0 0, so bra and ket model space are the same 
        // But we'll need a different bra model space for Brel
        proton_neutron::ModelSpaceMapType model_space_map_bra=generate_bra_model_space_Brel(lgi,lgi_vector);

        // If model_space_map_bra is not empty, then construct the basis for calculating B
        bool compute_B_matrix=false;
        lsu3::CncsmSU3xSU2Basis bra_ncsm_basis;
        int dim_B_bra = 0;
        if (model_space_map_bra.size()!=0)
          {
            
            proton_neutron::ModelSpace bra_ncsmModelSpace(Z,N,model_space_map_bra);
            bra_ncsm_basis.ConstructBasis(bra_ncsmModelSpace, jdiag, ndiag, individual_comm);
            dim_B_bra = lsu3shell::get_num_U3PNSPN_irreps(bra_ncsm_basis);
            compute_B_matrix=true; 
          } 
        int dim_bra = dim_ket + dim_B_bra;

        // Initialize matrix
        basis::OperatorBlock<double>& BNcm_matrix = lgi_expansions[l];
        BNcm_matrix = Eigen::MatrixXd::Zero(dim_bra,dim_ket); 
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate Ncm matrix 
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        std::string nrel_operator_base_name = fmt::format("{}/Nrel",operator_dir);
        BNcm_matrix.block(0,0,dim_ket,dim_ket)=lgi::NcmMatrix(Z+N,double(N0+Nex),baseSU3Irreps,ket_ncsm_basis,ket_ncsm_basis,nrel_operator_base_name);
        //////////////////////////////////////////////////////////////////////////////////////////   
        // Calculate Brel matrix 
        /////////////////////////////////////////////////////////////////////////////////////
        if(compute_B_matrix)
        {

          std::string brel_operator_base_name = fmt::format("{}/Brel",operator_dir);
          BNcm_matrix.block(dim_ket,0,dim_B_bra,dim_ket)
            =BMatrix(baseSU3Irreps,bra_ncsm_basis,ket_ncsm_basis,brel_operator_base_name);         
        }

        //////////////////////////////////////////////////////////////////////////////////////////   
        // Get null vectors 
        /////////////////////////////////////////////////////////////////////////////////////
        Eigen::MatrixXd null_vectors =lgi::FindNullSpaceSVD(BNcm_matrix,gamma_max,1e-6);

        std::cout<<"Null vectors"<<std::endl;
        std::cout<<null_vectors<<std::endl;
    }

    /////////////////////////////////////////////////////////////////////////////////////

    // clears memory allocated for U9, U6, and Z6 coefficients
    CWig9lmLookUpTable<RME::DOUBLE>::ReleaseMemory();
    // clear memory allocated for single-shell SU(3) rmes
    CSSTensorRMELookUpTablesContainer::ReleaseMemory();  
    // clear memory allocated for SU(3)>SO(3)
    CWigEckSU3SO3CGTablesLookUpContainer::ReleaseMemory();  

    su3::finalize();

    return lgi_expansions;
  }



  // u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1B,relative_unit_tensor_labels,J0,T0,false);


  void ComputeSeeds(
    const std::tuple<MultiplicityTagged<lgi::LGI>,MultiplicityTagged<lgi::LGI>>& lgi_pair,
    const basis::OperatorBlock<double>& lgi_expansion_bra,
    const basis::OperatorBlock<double>& lgi_expansion_ket,
    const const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels
    )
  {
    //Extract label information
    const auto&[lgi_bra,gamma_max_bra] = lgi_pair[0];
    const auto&[Nex_bra,sigma_bra,Sp_bra,Sn_bra,S_bra]=lgi_bra.Key();

    const auto&[lgi_ket,gamma_max_ket] = lgi_pair[0];
    const auto&[Nex_ket,sigma_ket,Sp_ket,Sn_ket,S_ket]=lgi_ket.Key();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Construct the basis 
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

    // Model space consists of irreps with quantum numbers Nex (lambda,mu) Sp Sn S
    // num in communicator must be equal to ndiag
    // map used for defining model_space_map organized into {N : {SpSn : {S : <(lambda,mu)>}}} 
    //
    // Construct ket basis from model space
    proton_neutron::ModelSpaceMapType model_space_map_ket;
    model_space_map_ket[Nex_ket][{TwiceValue(Sp_ket),TwiceValue(Sn_ket)}][TwiceValue(S_ket)].emplace_back(1,sigma_ket.SU3().lambda(),sigma_ket.SU3().mu());
    proton_neutron::ModelSpace ket_ncsmModelSpace(Z,N,model_space_map_ket);
    lsu3::CncsmSU3xSU2Basis ket_ncsm_basis;
    ket_ncsm_basis.ConstructBasis(ket_ncsmModelSpace, jdiag, ndiag, individual_comm);
    int dim_ket = lsu3shell::get_num_U3PNSPN_irreps(ket_ncsm_basis);
    // Construct bra basis from model space
    proton_neutron::ModelSpaceMapType model_space_map_bra;
    model_space_map_bra[Nex_bra][{TwiceValue(Sp_bra),TwiceValue(Sn_bra)}][TwiceValue(S_bra)].emplace_back(1,sigma_bra.SU3().lambda(),sigma_bra.SU3().mu());
    proton_neutron::ModelSpace bra_ncsmModelSpace(Z,N,model_space_map_bra);
    lsu3::CncsmSU3xSU2Basis bra_ncsm_basis;
    bra_ncsm_basis.ConstructBasis(bra_ncsmModelSpace, jdiag, ndiag, individual_comm);
    int dim_bra = lsu3shell::get_num_U3PNSPN_irreps(bra_ncsm_basis);

    for(int i=0; i<unit_tensor_labels.size(); ++i)
      {
        const auto& [operator_labels,relative_bra, relative_ket]=unit_tensor_labels[i].Key();
        const auto& [Nop,x0,S0,T0,g0] = operator_labels;


      }


  }


}// end namespace

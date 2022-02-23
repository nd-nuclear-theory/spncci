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
#include "UNU3SU3/UNU3SU3Basics.h"

#include "lsu3shell/lsu3shell_basis.h"
#include "lgi/dimensions.h"
#include "su3rme.h"
#include "lgi/null_solver.h"
#include "am/am.h"

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


  basis::OperatorBlock<double> generate_lgi_expansion(
    const nuclide::NuclideType& nuclide,
    const int Nsigma_max,
    const unsigned int N0,
    const CBaseSU3Irreps& baseSU3Irreps,
    const lgi::MultiplicityTaggedLGIVector& lgi_vector,
    const int lgi_index,
    const std::string& operator_dir,
    const MPI_Comm& individual_comm,
    double zero_threshold
    )
  {
    // Notes: probably want to pull out operator load into main program and pass in a argument
    // If so, need to understand how Nop selection rules applied internally.
    // Nrel and Brel would have to be combined into one operator because of funky pointer modification
    // barried deep in lsu3shell code but will have two different Nop which may be problematic.
    // Aslo want to switch to using B+Ncm rather than Bintr+Ncm.
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Calculate rmes of B and Ncm for each U3SpSnS subspace corresponding to an LGI
    // Then solve for simultaneous null space of B and Ncm
    //Extract label information
    const auto&[lgi,gamma_max] = lgi_vector[lgi_index];
    const auto&[Nex,sigma,Sp,Sn,S]=lgi.Key();

    // If Nex==0, then all of the irreps in the subspace are lgi and the basis
    // for the null spaces is given by the identity matrix
    if(Nex==0)
    {
      return Eigen::MatrixXd::Identity(gamma_max,gamma_max);
    }

    // With communicator split, ndiag is always 1
    int jdiag=0;
    int ndiag=1;
    const auto&[Z,N]=nuclide;

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
    basis::OperatorBlock<double> BNcm_matrix = Eigen::MatrixXd::Zero(dim_bra,dim_ket);
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
    basis::OperatorBlock<double> lgi_expansion = lgi::FindNullSpaceSVD(BNcm_matrix,gamma_max,zero_threshold);

    //TODO: Make optional
    if(true)
      {
        assert(mcutils::IsZero(BNcm_matrix*lgi_expansion,zero_threshold));
        assert(mcutils::IsZero(BNcm_matrix.block(dim_ket,0,dim_B_bra,dim_ket)*lgi_expansion));
        assert(mcutils::IsZero(BNcm_matrix.block(0,0,dim_ket,dim_ket)*lgi_expansion,zero_threshold));
        assert(not mcutils::IsZero(lgi_expansion,zero_threshold));
        assert(mcutils::IsZero(
          lgi_expansion.transpose()*lgi_expansion-Eigen::MatrixXd::Identity(gamma_max,gamma_max),
          zero_threshold
          ));
      }


    return lgi_expansion;
  }



  basis::OperatorBlocks<double> generate_lgi_expansion(
    const nuclide::NuclideType& nuclide,
    const int Nsigma_max,
    const lgi::MultiplicityTaggedLGIVector& lgi_vector,
    const std::string& operator_dir,
    const MPI_Comm world_comm //set default to be MPI_COMM_WORLD
    )
  {
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

    // Define distributions of partices over shells and do basic basis setup
    const auto&[Z,N]=nuclide;
    CBaseSU3Irreps baseSU3Irreps(Z,N,Nsigma_max);
    unsigned int N0 = nuclide::N0ForNuclide(nuclide);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Calculate rmes of B and Ncm for each U3SpSnS subspace corresponding to an LGI
    // Then solve for simultaneous null space of B and Ncm
    basis::OperatorBlocks<double> lgi_expansions(lgi_vector.size());
    for(int l=0; l<lgi_vector.size(); ++l)
      {  

        lgi_expansions[l] =
          generate_lgi_expansion(
            nuclide,Nsigma_max,N0,
            baseSU3Irreps,
            lgi_vector,l,
            operator_dir,
            individual_comm
            );
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


  std::vector<basis::OperatorBlocks<double>>
   ComputeSeeds(
    const nuclide::NuclideType& nuclide,
    const int Nsigma_max,
    const int N1v,
    const std::pair<MultiplicityTagged<lgi::LGI>,MultiplicityTagged<lgi::LGI>>& lgi_pair,
    const basis::OperatorBlock<double>& lgi_expansion_bra,
    const basis::OperatorBlock<double>& lgi_expansion_ket,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
    const std::string& operator_dir,
    const MPI_Comm world_comm, //set default to be MPI_COMM_WORLD
    const bool& restrict_op_J0
    )
  {
    std::vector<basis::OperatorBlocks<double>> lgi_unit_tensor_rmes(unit_tensor_labels.size());
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Set global variables
    //
    // If not false, interactionPN.AddOperator will hang in cases where x0 doesn't branch to S0.
    k_dependent_tensor_strenghts=restrict_op_J0;
    // If not false, get segmentation fault when calling function Calculate_Proton_x_Identity_MeData
    // caused by WigEckSU3SO3CG::WigEckSU3SO3CG trying to calculate coupling coefficients for
    // L0 not in x0.
    precalculate_WigEckSU3SO3CG_coefficients = false;
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialization for calculating coupling and recoupling coefficients
    su3::init();
    CWig9lmLookUpTable<RME::DOUBLE>::AllocateMemory(true);

    //Extract label information
    const auto&[Z,N]=nuclide;
    const auto&[lgi_bra,gamma_max_bra] = lgi_pair.first;
    const auto&[Nex_bra,sigma_bra,Sp_bra,Sn_bra,S_bra]=lgi_bra.Key();

    const auto&[lgi_ket,gamma_max_ket] = lgi_pair.second;
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
    // 
    // int Nsmax= std::max(Nex_ket,Nex_bra);
    CBaseSU3Irreps baseSU3Irreps(Z,N,Nsigma_max);
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
        const auto& [dN0,x0,S0,T0,g0] = operator_labels;
        const auto& [Nbar,Sbar,Tbar] = relative_ket;
        const auto& [Nbarp,Sbarp,Tbarp] = relative_bra;

        // w0(rho0,lambda0,mu0,twice_S0)
        
        SU3xSU2::LABELS w0(1,x0.lambda(),x0.mu(),TwiceValue(S0));
        int rhot_max=u3::OuterMultiplicity(sigma_ket.SU3(),x0,sigma_bra.SU3());

        bool allowed_operator = true;
        allowed_operator &= rhot_max>0;
        allowed_operator &= am::AllowedTriangle(S_ket,S0,S_bra);
        allowed_operator &= (Nex_ket+2*N1v)>=Nbar;
        allowed_operator &= (Nex_bra+2*N1v)>=Nbarp;


        ////////////////////////////////////////////////////////////////////////////////
        // Temporary for testing.  To be replaced by call to UnitTensorAllowed(...) in
        // recurrence_indexing_spin
        ////////////////////////////////////////////////////////////////////////////////
        bool spin_selection_satisfied;
        if (am::AllowedTriangle(Sp_bra,1,Sp_ket) && am::AllowedTriangle(Sn_bra,1,Sn_ket))
          spin_selection_satisfied=true;

        if (S0==0)
        {
          if (Sp_bra==Sp_ket && Sn_bra==Sn_ket)
            spin_selection_satisfied=true;
        }
        else if (S0==1)
        {
          if (Sp_bra==Sp_ket && am::AllowedTriangle(Sn_bra,1,Sn_ket))
            spin_selection_satisfied=true;

          if (am::AllowedTriangle(Sp_bra,1,Sp_ket) && Sn_bra==Sn_ket)
            spin_selection_satisfied=true;
        }
        else if (S0==2 && Tbarp==1 && Tbar==1)
        {
          if (Sp_bra==Sp_ket && am::AllowedTriangle(Sn_bra,2,Sn_ket))
            spin_selection_satisfied=true;

          if (am::AllowedTriangle(Sp_bra,2,Sp_ket) && Sn_bra==Sn_ket)
            spin_selection_satisfied=true;
        }
        else
          spin_selection_satisfied=false;


        allowed_operator&=spin_selection_satisfied;
        ////////////////////////////////////////////////////////////////////////////////



        if(not allowed_operator)
          {
            //todo: set operator to zero
            continue; 
          }

        // Read operator from file
        std::cout<<unit_tensor_labels[i].Str()<<std::endl;
        std::string operator_base_name = fmt::format("{}/relative_unit_{:06}",operator_dir,i);
        // std::cout<<"operator "<<operator_base_name<<std::endl;
        std::ofstream interaction_log_file("/dev/null");
        bool log_is_on = false;
        bool generate_missing_rme = true;

        // Apparently the constructor for interactionPPNN must always be called before the constructor for interactionPN
        // otherwise a global variable doesn't get correctly allocated....
        CInteractionPPNN interactionPPNN(baseSU3Irreps,log_is_on,interaction_log_file);
        CInteractionPN interactionPN(baseSU3Irreps,generate_missing_rme,log_is_on,interaction_log_file);
        
        interactionPPNN.LoadTwoBodyOperator(operator_base_name+".PPNN");
        interactionPN.AddOperator(operator_base_name+".PN");
        interactionPPNN.TransformTensorStrengthsIntoPP_NN_structure();

        // std::cout<<"Calculate RME"<<std::endl;
        basis::OperatorBlocks<double> su3_unit_tensor_rmes
          =lsu3shell::CalculateRME(interactionPPNN,interactionPN,bra_ncsm_basis,ket_ncsm_basis,dN0,w0,rhot_max);

        auto& lgi_pair_rmes = lgi_unit_tensor_rmes[i];
        lgi_pair_rmes.resize(rhot_max);
        for(int irhot=0; irhot<rhot_max; ++irhot)
          {
            // std::cout<<lgi_expansion_bra.transpose().rows()<<" x "<<lgi_expansion_bra.transpose().cols()<<"   "
            // <<su3_unit_tensor_rmes[irhot].rows()<<" x "<<su3_unit_tensor_rmes[irhot].cols()<<"   "
            // <<lgi_expansion_ket.rows()<<" x "<<lgi_expansion_ket.cols()<<std::endl;
            // std::cout<<mcutils::FormatMatrix(lgi_expansion_bra.transpose(),"3.2f")<<std::endl<<std::endl
            // <<mcutils::FormatMatrix(su3_unit_tensor_rmes[irhot],"3.2f")<<std::endl<<std::endl
            // <<mcutils::FormatMatrix(lgi_expansion_ket,"3.2f")<<std::endl<<std::endl;
            
            lgi_pair_rmes[irhot]=lgi_expansion_bra.transpose()*su3_unit_tensor_rmes[irhot]*lgi_expansion_ket;
          }
      }


    // clears memory allocated for U9, U6, and Z6 coefficients
    CWig9lmLookUpTable<RME::DOUBLE>::ReleaseMemory();
    // clear memory allocated for single-shell SU(3) rmes
    CSSTensorRMELookUpTablesContainer::ReleaseMemory();  
    // clear memory allocated for SU(3)>SO(3)
    CWigEckSU3SO3CGTablesLookUpContainer::ReleaseMemory();  

    su3::finalize();

    return lgi_unit_tensor_rmes;
  }


}// end namespace

/****************************************************************
  spncci_branching_u3s_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iostream>
#include <fstream>

#include "cppformat/format.h"
#include "lgi/lgi.h"
#include "spncci/spncci_basis.h"
#include "spncci/spncci_branching_u3s.h"
#include "spncci/spncci_branching_u3lsj.h"
#include "spncci/unit_tensor.h"
#include "u3shell/relative_operator.h"
#include "u3shell/upcoupling.h"

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // set up SpNCCI space
  ////////////////////////////////////////////////////////////////

  // read in LGIs
  std::string filename = "lgi_test.dat";  // test file in data/lgi_set/lgi_test.dat
  lgi::MultiplicityTaggedLGIVector multiplicity_tagged_lgi_vector;
  lgi::ReadLGISet(multiplicity_tagged_lgi_vector,filename);

  // generate SpNCCI space from LGIs
  HalfInt Nsigma_0 = HalfInt(11,1);
  int Nmax = 2;
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;  // dictionary from sigma to branching
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  spncci::GenerateSpNCCISpace(multiplicity_tagged_lgi_vector,truncator,spncci_space,sigma_irrep_map);

  ////////////////////////////////////////////////////////////////
  // construct flattened baby SpNCCI space
  ////////////////////////////////////////////////////////////////

  // put SpNCCI space into standard linearized container
  spncci::BabySpNCCISpace baby_spncci_space(spncci_space);

  // diagnostic
  std::cout << "baby_spncci_space" << std::endl;
  for (int subspace_index=0; subspace_index<baby_spncci_space.size(); ++subspace_index)
    std::cout << baby_spncci_space.GetSubspace(subspace_index).DebugStr()
              << std::endl;
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // read interaction coefficients
  ////////////////////////////////////////////////////////////////

  // Recall storage scheme for u3shell::RelativeRMEsU3ST:
  //
  // relative_rmes (u3shell::RelativeRMEsU3ST) maps
  //
  //   (unit_tensor_labels,kappa0,L0) -> RME
  //
  // where
  // 
  // unit_tensor_labels (u3shell::RelativeUnitTensorLabelsU3ST) contains
  //
  //   u3shell::OperatorLabelsU3ST operator_labels
  //   u3shell::RelativeStateLabelsU3ST bra
  //   u3shell::RelativeStateLabelsU3ST ket

    std::string interaction_filename = "unit.dat";
    std::ifstream interaction_stream(interaction_filename);
    assert(interaction_stream);
    u3shell::RelativeRMEsU3ST relative_rmes;
    u3shell::ReadRelativeOperatorU3ST(interaction_stream,relative_rmes);

    // diagnostic -- relative rme contents
    std::cout << "relative rme list" << std::endl;
    for (const auto& label_rme_pair : relative_rmes)
      {

        u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels;
        int kappa0, L0;
        std::tie(unit_tensor_labels,kappa0,L0) = label_rme_pair.first;
        double rme = label_rme_pair.second;

        // cast unit tensor labels into U3S operator labels
        u3shell::OperatorLabelsU3S operator_labels_u3s = unit_tensor_labels;
        
        std::cout << fmt::format(
            "unit tensor {} kappa0 {} L0 {} -> U3S symmetry {} rme {}",
            unit_tensor_labels.Str(),kappa0,L0,operator_labels_u3s.Str(),rme
          ) << std::endl;
      }
    std::cout << std::endl;
  
  ////////////////////////////////////////////////////////////////
  // construct fake unit tensor recurrence results
  ////////////////////////////////////////////////////////////////

  // For each allowed BabySpNCCI sector, we will fill relevant unit
  // tensor matrix with a constant matrix.

  // Recall storage scheme for spncci::UnitTensorMatricesByIrrepFamily:
  //
  // (bra_irrep_family_index,ket_irrep_family_index) (<int,int>)
  //   -> (bra_Nn,ket_Nn) (<int,int>)
  //   -> unit_tensor_sector_labels (spncci::UnitTensorU3Sector)
  //   -> sector matrix (Eigen::MatrixXd)

  spncci::UnitTensorMatricesByIrrepFamily unit_tensor_matrices;
  
  // TODO finish construction...  not worth the effort?

  for (const auto& label_rme_pair : relative_rmes)
    // for each contribution in coefficient list
    {

      // extract (unit_tensor_labels,kappa0,L0)
      u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels;
      int kappa0, L0;
      std::tie(unit_tensor_labels,kappa0,L0) = label_rme_pair.first;
      // double rme = label_rme_pair.second;

      // recombine unit tensor labels as needed in spncci::UnitTensorMatricesByIrrepFamily
      u3shell::OperatorLabelsU3S operator_labels_u3s = unit_tensor_labels;  // implicit cast to U3S labels
      int bra_Nn, ket_Nn;
      // spncci::UnitTensorU3Sector unit_tensor_sector_labels(...);
      // UnitTensorU3Sector(u3::U3 omegap, u3::U3 omega, u3shell::RelativeUnitTensorLabelsU3ST tensor, int rho0)


      // FlatKey x0_,S0_,T0_,bra_.eta(),bra_.S(),bra_.T(),ket_.eta(),ket_.S(),ket_.T()
    }


  ////////////////////////////////////////////////////////////////
  // baby SpNCCI sector construction
  ////////////////////////////////////////////////////////////////

  // for each contributing (N0,x0,S0) in coef list:
  //   generate source (sigma,Sp,Sn,S,omega) sectors (BabySpNCCISectors) for this (N0,x0,S0)
  //   generate target (omega,S) sectors (SpNCCISectorsU3) for this (N0,x0,S0)
  //   for each (kappa0,L0) within (N0,x0,S0) in coef list:
  //
  //      Alternative #1:
  //
  //      zero initialize matrices for all target sectors
  //      for each source sector:
  //        for each contributing operator ([N0,x0,S0];T0,eta_prime,...;[kappa0,L0]) in coef list:
  //          accumulate source sector matrix into corresponding block of *appropriate* target sector matrix
  //
  //      Alternative #2:
  //
  //        for each target (omega,S) sector
  //          zero initialize target sector matrix
  //          for each source (sigma,Sp,Sn,S,omega) sector contained within this target sector:
  //            for each contributing operator ([N0,x0,S0];T0,eta_prime,...;[kappa0,L0]) in coef list:
  //              accumulate source sector matrix into corresponding block of *current* target sector matrix
  //     
  //        Cons: requires us to filter source sectors by current target sector
  //
  //
  // Prerequisites:
  //
  //   - Operator symmetry (N0,x0,S0) and corresponding upcoupling (kappa0,L0) labels
  //
  //     (N0,x0,S0) -> {(kappa0,L0),...}
  //
  //   - ... ideally also with pointer to 

  // accumulate relevant symmetry and upcoupling labels for given interaction

  // Note: can make map and set unordered, but that is less pretty for debugging
  std::map<u3shell::OperatorLabelsU3S,std::set<u3shell::UpcouplingLabels>> operator_u3s_upcoupling_labels;
  for (const auto& label_rme_pair : relative_rmes)
    // for each contribution in coefficient list
    {

      // extract operator symmetry and upcoupling labels
      u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels;
      int kappa0, L0;
      std::tie(unit_tensor_labels,kappa0,L0) = label_rme_pair.first;
      u3shell::OperatorLabelsU3S operator_labels_u3s = unit_tensor_labels;  // implicit cast to U3S labels
      u3shell::UpcouplingLabels upcoupling_labels(kappa0,L0);

      // accumulate
      operator_u3s_upcoupling_labels[operator_labels_u3s].insert(upcoupling_labels);
    }

  // Alternatively, can extract from:
  // std::vector<u3shell::IndexedOperatorLabelsU3S> operator_u3s_kappa0_L0_list;
  // u3shell::GetInteractionTensorsU3S(interaction_rmes,operator_u3s_kappa0_L0_list);

  for (const auto& u3s_upcoupling_labels_set_pair : operator_u3s_upcoupling_labels)
    {
      // recover operator labels and corresponding set of upcoupling labels
      const u3shell::OperatorLabelsU3S& operator_labels_u3s = u3s_upcoupling_labels_set_pair.first;
      const std::set<u3shell::UpcouplingLabels> upcoupling_labels_set = u3s_upcoupling_labels_set_pair.second;

      // generate source sectors
      spncci::BabySpNCCISectors source_sectors(baby_spncci_space,operator_labels_u3s);
      std::cout << fmt::format("Sectors for u3s {}",operator_labels_u3s.Str()) << std::endl;
      
      for (int sector_index=0; sector_index < source_sectors.size(); ++sector_index)
        {
          const auto& source_sector = source_sectors.GetSector(sector_index);
          const auto& bra_subspace = source_sector.bra_subspace();
          const auto& ket_subspace = source_sector.ket_subspace();
          std::cout
            << fmt::format("  ({:3d},{:3d}) : {} -- {}",
                           source_sector.bra_subspace_index(),
                           source_sector.ket_subspace_index(),
                           source_sector.bra_subspace().LabelStr(),
                           source_sector.ket_subspace().LabelStr()
              )
            << std::endl;
         }
      std::cout << std::endl;
    }


  ////////////////////////////////////////////////////////////////
  // Regroup state tests -- aem old tests
  ////////////////////////////////////////////////////////////////

  if (true)
    {
      std::cout<<"Regroup test"<<std::endl;
      // build space
      spncci::SpaceU3S space(spncci_space);
      std::cout<<"irreps "<< spncci_space.size()<<std::endl;
      for(auto irrep :spncci_space)
        std::cout<<irrep.Str()<<std::endl;
      // dump subspace contents 
      std::cout<<"Regrouping"<<std::endl;     
      for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
        {
          const spncci::SubspaceU3S& subspace = space.GetSubspace(subspace_index);
          std::cout << fmt::format("subspace {}  {}",subspace_index,subspace.Str())
                    << " with dimension "<< subspace.sector_dim()
                    <<"  and subspace size "<<subspace.size()<<std::endl;
          for (int state_index=0; state_index<subspace.size(); ++state_index)
            {
              int baby_spncci_index;
              std::tie(baby_spncci_index)=subspace.GetStateLabels(state_index);
              std::cout << baby_spncci_index<< std::endl;
            } 
        }
    }
  ////////////////////////////////////////////////////////////////
  // Regroup U3S sector tests -- aem old tests
  ////////////////////////////////////////////////////////////////  
  if(false)
  {
    std::cout<<"U3S Sectors "<<std::endl;
    spncci::SpaceU3S space(spncci_space);
    // spncci::SectorLabelsU3SCache u3s_sectors;
    // To test sector construction
    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_tensor_labels;
    u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_tensor_labels);
    std::vector<u3shell::IndexedOperatorLabelsU3S> tensor_labels;
    for(auto unit_tensor : relative_tensor_labels)
      {
        // std::cout<<unit_tensor.Str()<<std::endl;
        u3shell::OperatorLabelsU3S operator_labels(unit_tensor.operator_labels());
        MultiplicityTagged<int>::vector L0_kappa0=u3::BranchingSO3(operator_labels.x0());
        int L0=L0_kappa0[0].irrep;
        int kappa0=L0_kappa0[0].tag;
        // std::cout<<"tensor "<<operator_labels.Str()<<" "<<kappa0<<"  "<<L0<<std::endl;
        tensor_labels.push_back(u3shell::IndexedOperatorLabelsU3S(operator_labels,kappa0,L0));
      }

    std::vector<spncci::SectorLabelsU3S> sector_vector;

    spncci::GetSectorsU3S(space,tensor_labels,sector_vector);
    std::cout<<"number of sectors "<<sector_vector.size()<<std::endl;
    int i=0;
    basis::MatrixVector matrix_vector(sector_vector.size());
    std::cout<<"space"<<std::endl;
    for(int i=0; i<space.size(); ++i)
      std::cout<<space.GetSubspace(i).GetSubspaceLabels().Str()<<std::endl;
    for(auto& sector : sector_vector)
      {
        // if(i<10)
        bool allowed=sector.bra_index()==sector.ket_index();
        allowed&=sector.bra_index()<3;
        if(allowed)
        {
          u3::U3S omegaS_bra=space.GetSubspace(sector.bra_index()).GetSubspaceLabels();
          u3::U3S omegaS_ket=space.GetSubspace(sector.ket_index()).GetSubspaceLabels();
          std::cout<<fmt::format(" {}",sector.Str())<<std::endl;
          std::cout<<"  "<<omegaS_bra.Str()<<"  "<< omegaS_ket.Str()<<std::endl;

          
          const spncci::SubspaceU3S& ket_subspace=space.GetSubspace(sector.ket_index());
          const spncci::SubspaceU3S& bra_subspace=space.GetSubspace(sector.bra_index());
          int sector_dim_bra=bra_subspace.sector_dim();
          int sector_dim_ket=ket_subspace.sector_dim();
          matrix_vector[i]=Eigen::MatrixXd::Zero(sector_dim_bra,sector_dim_ket);
          // std::cout<<matrix_vector[i]<<std::endl;

        }
        ++i;
      }
  }
  ////////////////////////////////////////////////////////////////
  // Regroup LS test state tests -- aem old tests
  ////////////////////////////////////////////////////////////////
  if (true)
    {
      std::cout<<"Regroup LS test"<<std::endl;
      spncci::SpaceU3S space(spncci_space);
      // build space
      spncci::SpaceLS space_ls(space,0);
      // dump subspace contents 
      std::cout<<"Regrouping"<<std::endl;     
      for (int subspace_index=0; subspace_index<space_ls.size(); ++subspace_index)
        {
          const spncci::SubspaceLS& subspace_ls = space_ls.GetSubspace(subspace_index);
          std::cout << fmt::format("subspace {}",subspace_ls.Str())
                    << "with dimension "<< subspace_ls.sector_dim()
                    <<"  "<<subspace_ls.size()<<std::endl;
          for (int state_index=0; state_index<subspace_ls.size(); ++state_index)
            {
              
              int u3_subspace_index;
              std::tie(u3_subspace_index)=subspace_ls.GetStateLabels(state_index);
              std::cout<<u3_subspace_index<<std::endl;
              // const spncci::StateLS state_ls(subspace_ls,state_index);
              // std::cout << state_ls<< std::endl;
            } 
        }
    }
} //main

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

  ////////////////////////////////////////////////////////////////
  // read interaction coefficients
  ////////////////////////////////////////////////////////////////

  std::string interaction_filename = "unit.dat";
  std::ifstream interaction_stream(interaction_filename);
  assert(interaction_stream);
  u3shell::RelativeRMEsU3ST relative_rmes;
  u3shell::ReadRelativeOperatorU3ST(interaction_stream,relative_rmes);

  // Recall:
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



  // diagnostic -- relative rme contents
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
      // FlatKey x0_,S0_,T0_,bra_.eta(),bra_.S(),bra_.T(),ket_.eta(),ket_.S(),ket_.T()
    }

  ////////////////////////////////////////////////////////////////
  // construct fake unit tensor recurrence results
  ////////////////////////////////////////////////////////////////

  // For each allowed BabySpNCCI sector, we will fill relevant unit
  // tensor matrix with constants.
  spncci::UnitTensorMatricesByIrrepFamily unit_tensor_matrices;

  // u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor_labels = 
  //       const u3shell::OperatorLabelsU3ST operator_labels,
  //       const u3shell::RelativeStateLabelsU3ST& bra,
  //       const u3shell::RelativeStateLabelsU3ST& ket
     

  // (bra_irrep_family_index,ket_irrep_family_index) (<int,int>)
  //   -> (bra_Nn,ket_Nn) (<int,int>)
  //   -> unit_tensor_sector_labels (spncci::UnitTensorU3Sector)
  //   -> sector matrix (Eigen::MatrixXd)
  
  

  ////////////////////////////////////////////////////////////////
  // baby SpNCCI sector construction
  ////////////////////////////////////////////////////////////////

  //spncci::BabySpNCCISectors scalar_sectors(baby_spncci_space,)

  ////////////////////////////////////////////////////////////////
  // Regroup state tests
  ////////////////////////////////////////////////////////////////

  if (false)
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
              const spncci::StateU3S state(subspace,state_index);
              std::cout << state.Str()<< std::endl;
            } 
        }
    }
  ////////////////////////////////////////////////////////////////
  // Regroup U3S sector tests
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
  // Regroup LS test state tests
  ////////////////////////////////////////////////////////////////
  if (false)
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
              const spncci::StateLS state_ls(subspace_ls,state_index);
              std::cout << state_ls.Str()<< std::endl;
            } 
        }
    }
} //main

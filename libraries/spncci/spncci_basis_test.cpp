/****************************************************************
  spncci_basis_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "spncci/spncci_basis.h"

#include "cppformat/format.h"
#include "lgi/lgi.h"

#include "u3shell/relative_operator.h"
#include "spncci/spncci_branching_u3s.h"
#include "spncci/spncci_branching_u3lsj.h"
int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // read in LGIs
  ////////////////////////////////////////////////////////////////

  HalfInt Nsigma_0 = HalfInt(11,1);
  // std::string filename = "libraries/spncci/sp_basis_test.dat";
  std::string filename = "lgi_test.dat";  // test file in data/lgi_set/lgi_test.dat
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  lgi::ReadLGISet(lgi_vector,filename);

  // diagnostic -- inspect LGI listing
  std::cout << "LGI vector" << std::endl;
  for (int i=0; i<lgi_vector.size(); ++i)
    std::cout << i << " " << lgi_vector[i].Str() << std::endl;
  std::cout << "********************************" << std::endl;

  ////////////////////////////////////////////////////////////////
  // build irreps
  ////////////////////////////////////////////////////////////////

  int Nmax = 2;
  spncci::SpNCCISpace sp_irrep_vector;
  spncci::SigmaIrrepMap sigma_irrep_map;  // dictionary from sigma to branching
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  spncci::GenerateSpNCCISpace(lgi_vector,truncator,sp_irrep_vector,sigma_irrep_map);

  // diagnostic -- inspect irrep listing and contents
  if(false)
  {
    std::cout << "SpNCCIIrrepFamily vector reprise" << std::endl;
    for (int i=0; i<sp_irrep_vector.size(); ++i)
      std::cout << i << " " << sp_irrep_vector[i].DebugString()
                <<"  "<<sp_irrep_vector[i].gamma_max();

    // examine irreps
    std::cout << "irreps (by sigma)" << std::endl;
    for (auto it = sigma_irrep_map.begin(); it != sigma_irrep_map.end(); ++it)
      {
        u3::U3 sigma = it->first;
        const sp3r::Sp3RSpace& irrep = it->second;

        std::cout << irrep.DebugString();
        std::cout << std::endl;
      }
    std::cout << "********************************" << std::endl;
  }

  // diagnostic -- inspect irrep listing and contents
  if(false)
  {
    // examine irreps by calling reference to irreps from vector
    std::cout << "irreps (by lgi)" << std::endl;
    for (auto sp_irrep_tag : sp_irrep_vector)
      {
        std::cout<<"Get irrep"<<std::endl;
        const spncci::SpNCCIIrrepFamily& sp_irrep = sp_irrep_tag;
        std::cout<<"Get subspace"<<std::endl;
        const sp3r::Sp3RSpace& irrep = sp_irrep.Sp3RSpace();
        std::cout<<"Debug"<<std::endl;
        std::cout << irrep.DebugString();
        std::cout << std::endl;
      }
    std::cout << "********************************" << std::endl;
  }


  ////////////////////////////////////////////////////////////////
  // indicate dimensions
  ////////////////////////////////////////////////////////////////

  if(false)
  {
    std::cout << "TotalU3Subspaces " << spncci::TotalU3Subspaces(sp_irrep_vector) << std::endl;
    std::cout << "TotalDimensionU3 " << spncci::TotalDimensionU3S(sp_irrep_vector) << std::endl;
    std::cout << "TotalDimensionU3LS " << spncci::TotalDimensionU3LS(sp_irrep_vector) << std::endl;
    std::cout << "TotalDimensionU3LSJConstrained ";
    for (HalfInt J=0; J<10; ++J)
      std::cout << J << " " << spncci::TotalDimensionU3LSJConstrained(sp_irrep_vector,J) << "    ";
    std::cout << std::endl;
    std::cout << "TotalDimensionU3LSJAll " << spncci::TotalDimensionU3LSJAll(sp_irrep_vector) << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  // Regroup state tests
  ////////////////////////////////////////////////////////////////
  if (true)
    {
      std::cout<<"Regroup test"<<std::endl;
      // build space
      spncci::SpaceU3S space(sp_irrep_vector);
      std::cout<<"irreps "<< sp_irrep_vector.size()<<std::endl;
      for(auto irrep :sp_irrep_vector)
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
  if(true)
  {
    std::cout<<"U3S Sectors "<<std::endl;
    spncci::SpaceU3S space(sp_irrep_vector);
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
      spncci::SpaceU3S space(sp_irrep_vector);
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

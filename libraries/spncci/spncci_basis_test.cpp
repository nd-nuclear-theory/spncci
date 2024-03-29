/****************************************************************
  spncci_basis_test.cpp

  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

 2/15/18 (aem): Update ReadLGISet
****************************************************************/
#include "spncci/spncci_basis.h"

#include "fmt/format.h"
#include "lgi/lgi.h"
#include "u3shell/relative_operator.h"


void GetBabySpncciBasis(spncci::BabySpNCCISpace& baby_spncci_space)
//Extracted from branching2_test.  TODO: Combine with other test examples. 
{
  // read in LGIs
  std::string filename 
  // ="/afs/crc.nd.edu/user/a/amccoy/research/codes/spncci/data/lgi_set/lgi_test.dat";
  ="data/lgi_set/lgi_test.dat";
  // = "~/research/codes/spncci/data/lgi_set/lgi_test.dat";  // test file in data/lgi_set/lgi_test.dat
  lgi::MultiplicityTaggedLGIVector multiplicity_tagged_lgi_vector;
  HalfInt Nsigma0=lgi::Nsigma0ForNuclide({3,3});
  lgi::ReadLGISet(filename, Nsigma0, multiplicity_tagged_lgi_vector);

  // generate SpNCCI space from LGIs
  // HalfInt Nsigma_0 = HalfInt(11,1);
  int Nmax = 2;
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;  // dictionary from sigma to branching
  spncci::NmaxTruncator truncator(Nsigma0,Nmax);
  spncci::GenerateSpNCCISpace(multiplicity_tagged_lgi_vector,truncator,spncci_space,sigma_irrep_map);

  ////////////////////////////////////////////////////////////////
  // construct flattened baby SpNCCI space
  ////////////////////////////////////////////////////////////////

  // put SpNCCI space into standard linearized container
  baby_spncci_space=spncci::BabySpNCCISpace(spncci_space);

  // diagnostic
  std::cout << "baby_spncci_space" << std::endl;
  for (int subspace_index=0; subspace_index<baby_spncci_space.size(); ++subspace_index)
    std::cout << baby_spncci_space.GetSubspace(subspace_index).DebugStr()
              << std::endl;
  std::cout << std::endl;


}


int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // shell counting calculations
  ////////////////////////////////////////////////////////////////
  if(false)
  {
    // std::cout << "Nsigma0ForNuclide" << std::endl;
    // std::cout
    //   << fmt::format(
    //       "3He {} 6Li {}",
    //       lgi::Nsigma0ForNuclide({2,1}),
    //       lgi::Nsigma0ForNuclide({3,3})
    //     )
    //   << std::endl;
    std::cout << "ValenceShellForNuclide" << std::endl;
    std::cout
      << fmt::format(
          "3He {} 6Li {}",
          spncci::ValenceShellForNuclide({2,1}),
          spncci::ValenceShellForNuclide({3,3})
        )
      << std::endl;
  }
  ////////////////////////////////////////////////////////////////
  // set up SpNCCI space
  ////////////////////////////////////////////////////////////////

  // read in LGIs for 6Li
  // std::string filename = "../../data/lgi_set/lgi_test.dat";  // test file in data/lgi_set/lgi_test.dat
  std::string filename = "../../data/lgi_set/lgi_test-2.dat";  // test file in data/lgi_set/lgi_test.dat
  lgi::MultiplicityTaggedLGIVector multiplicity_tagged_lgi_vector;
  HalfInt Nsigma0=lgi::Nsigma0ForNuclide({3,3});
  lgi::ReadLGISet(filename,Nsigma0,multiplicity_tagged_lgi_vector);

  if(false)
  {
    // diagnostic -- inspect LGI listing
    std::cout << "LGI set" << std::endl;
    for (int i=0; i<multiplicity_tagged_lgi_vector.size(); ++i)
      std::cout << i << " " << multiplicity_tagged_lgi_vector[i].Str() << std::endl;
    std::cout << "********************************" << std::endl;
  }

  // generate SpNCCI space from LGIs
  HalfInt Nsigma_0 = HalfInt(11,1);
  int Nmax = 2;
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;  // dictionary from sigma to branching
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  spncci::GenerateSpNCCISpace(multiplicity_tagged_lgi_vector,truncator,spncci_space,sigma_irrep_map);

  // diagnostic -- inspect irrep families
  if(false)
  {
    std::cout << "SpNCCI space" << std::endl;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {
        std::cout << "Irrep family" << std::endl;
        std::cout << "  " << spncci_irrep_family.Str() << std::endl;
        std::cout << "  " << fmt::format("Multiplicity gamma_max {}",spncci_irrep_family.gamma_max()) << std::endl;
        std::cout << "  " << "Irrep contents" << std::endl;
        std::cout << spncci_irrep_family.Sp3RSpace().DebugStr();
        std::cout << std::endl;
      }
    std::cout << "********************************" << std::endl;
  }



  // diagnostic -- inspect irrep families
  if(false)
  {
    //Alternate truncator
    // generate SpNCCI space from LGIs
    HalfInt Nsigma_0 = HalfInt(11,1);
    int Nmax = 6;
    spncci::SpNCCISpace spncci_space2;
    spncci::SigmaIrrepMap sigma_irrep_map2;  // dictionary from sigma to branching
    spncci::NlimitTruncator truncator2(Nsigma_0,Nmax,4);
    spncci::GenerateSpNCCISpace(multiplicity_tagged_lgi_vector,truncator2,spncci_space2,sigma_irrep_map2);

    std::cout << "SpNCCI space" << std::endl;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space2)
      {
        std::cout << "Irrep family" << std::endl;
        std::cout << "  " << spncci_irrep_family.Str() << std::endl;
        std::cout << "  " << fmt::format("Multiplicity gamma_max {}",spncci_irrep_family.gamma_max()) << std::endl;
        std::cout << "  " << "Irrep contents" << std::endl;
        std::cout << spncci_irrep_family.Sp3RSpace().DebugStr();
        std::cout << std::endl;
      }
    std::cout << "********************************" << std::endl;
  }


  ////////////////////////////////////////////////////////////////
  // count dimensions
  ////////////////////////////////////////////////////////////////

  if(false)
  {
    std::cout << fmt::format("  Irrep families {}",spncci_space.size()) << std::endl;
    std::cout << fmt::format("  TotalU3Subspaces {}",spncci::TotalU3Subspaces(spncci_space)) << std::endl;
    std::cout << fmt::format("  TotalDimensionU3 {}",spncci::TotalDimensionU3S(spncci_space)) << std::endl;
    std::cout << fmt::format("  TotalDimensionU3LS {}",spncci::TotalDimensionU3LS(spncci_space)) << std::endl;
    std::cout << "TotalDimensionU3LSJConstrained ";
    for (HalfInt J=0; J<10; ++J)
      std::cout << J << " " << spncci::TotalDimensionU3LSJConstrained(spncci_space,J) << "    ";
    std::cout << std::endl;
    std::cout << fmt::format("  TotalDimensionU3LSJAll {}",spncci::TotalDimensionU3LSJAll(spncci_space)) << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  // construct flattened baby SpNCCI space
  ////////////////////////////////////////////////////////////////

  // put SpNCCI space into standard linearized container
  spncci::BabySpNCCISpace baby_spncci_space(spncci_space);
  if(false)
  {
    // diagnostic
    std::cout << "baby_spncci_space" << std::endl;
    for (int subspace_index=0; subspace_index<baby_spncci_space.size(); ++subspace_index)
      std::cout << baby_spncci_space.GetSubspace(subspace_index).DebugStr()
                << std::endl;
    std::cout << std::endl;
  }
  ////////////////////////////////////////////////////////////////
  // construct baby SpNCCI Hypersectors
  ////////////////////////////////////////////////////////////////

  int N1v=1;
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1v,unit_tensor_labels,-1,-1,false);
  u3shell::RelativeUnitTensorSpaceU3S operator_space(Nmax,N1v,unit_tensor_labels);

  if(false)
  {
    spncci::BabySpNCCIHypersectors baby_spncci_hypersectors(baby_spncci_space,operator_space);
    std::cout<<"Baby SpNCCI Hypersectors"<<std::endl;
    std::cout<<baby_spncci_hypersectors.DebugStr()<<std::endl;
  }

  if(false)
  {
    std::cout<<"observable space"<<std::endl;
    // build list of observable labels 
    std::set<u3shell::IndexedOperatorLabelsU3S> observable_labels_set;
    for(auto& unit_tensor : unit_tensor_labels)
      {
        int N0=unit_tensor.N0();
        u3::SU3 x0=unit_tensor.x0();
        HalfInt S0=unit_tensor.S0();
        int kappa0_max=u3::BranchingMultiplicitySO3(x0,0);
        if(kappa0_max>0)
          {
            u3shell::OperatorLabelsU3S labels_u3s(N0,x0,S0);
            
            // kappa0=1 and L0=0
            observable_labels_set.emplace(labels_u3s,1,0);
          }
      }

    std::vector<u3shell::IndexedOperatorLabelsU3S>observable_labels;
    for(auto& labels : observable_labels_set)
      observable_labels.push_back(labels);
    // Construct observable space

    u3shell::ObservableSpaceU3S observable_space(observable_labels);
    std::cout<<observable_space.Str()<<std::endl;

    spncci::ObservableBabySpNCCIHypersectors observable_hypersectors(baby_spncci_space,observable_space);
    std::cout<<observable_hypersectors.DebugStr()<<std::endl;

  }

if(true)
  {
    std::cout<<"observable space"<<std::endl;
    // build list of observable labels 
    std::set<u3shell::IndexedOperatorLabelsU3S> observable_labels_set;
    for(auto& unit_tensor : unit_tensor_labels)
      {
        int N0=unit_tensor.N0();
        u3::SU3 x0=unit_tensor.x0();
        HalfInt S0=unit_tensor.S0();
        int kappa0_max=u3::BranchingMultiplicitySO3(x0,0);
        if(kappa0_max>0)
          {
            u3shell::OperatorLabelsU3S labels_u3s(N0,x0,S0);
            
            // kappa0=1 and L0=0
            observable_labels_set.emplace(labels_u3s,1,0);
          }
      }

    std::vector<u3shell::IndexedOperatorLabelsU3S>observable_labels;
    for(auto& labels : observable_labels_set)
      observable_labels.push_back(labels);
    // Construct observable space

    u3shell::ObservableSpaceU3S observable_space(observable_labels);
    // std::cout<<observable_space.Str()<<std::endl;

    spncci::ObservableBabySpNCCIHypersectors observable_hypersectors(baby_spncci_space,observable_space,1,0,true);
    std::cout<<observable_hypersectors.DebugStr()<<std::endl;
    // for(int observable_hypersector_index=0; 
    //     observable_hypersector_index<observable_hypersectors.size();
    //     ++observable_hypersector_index
    //   )
    //   {
    //     int baby_spncci_index_bra, baby_spncci_index_ket, operator_subspace_index, rho0;
    //     std::tie(baby_spncci_index_bra,baby_spncci_index_ket,operator_subspace_index,rho0)
    //       =list_baby_spncci_hypersectors[observable_hypersector_index];

    //     spncci::OperatorBlock& block=baby_spncci_observable_hyperblocks[observable_hypersector_index][0];
    //     int block_rows=block.rows();
    //     int block_cols=block.cols();

    //     // std::cout<<"Basis information"<<std::endl;
    //     const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra);
    //     const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket);
    //   }

  }


  if(false)
  {
    u3::U3CoefInit();
    std::cout<<"testing baby spncci reordering"<<std::endl;
    int num_irrep_families=multiplicity_tagged_lgi_vector.size();
    std::vector<int> num_subspaces;
    spncci::NumBabySpncciSubspacesInIrrepFamily(
      num_irrep_families,baby_spncci_space,num_subspaces
    );
    std::cout<<"gotten number of irrep families"<<std::endl;
    std::vector<int>ordered_subspaces;
    spncci::SortSubspacesDecending(num_subspaces,ordered_subspaces);

  }

  if(false)
  {
      HalfInt J(1,1); 
  spncci::BabySpNCCISpace baby_spncci_space;
  GetBabySpncciBasis(baby_spncci_space);
  spncci::SpaceSpBasis spbasis(baby_spncci_space,J);
  if(true)
  {
    std::cout<<"full dimension "<<spbasis.FullDimension()<<std::endl;
    std::cout<<"size "<<spbasis.size()<<std::endl;
   std::cout<<spbasis.DebugStr(true)<<std::endl;
  }

  if(true)
  {
    int full_dimension_initial=spbasis.FullDimension();
    std::cout<<full_dimension_initial<<std::endl;
    int full_dimension_final=0;
    for(int i=0; i<spbasis.size(); ++i)
    {
      const spncci::SubspaceSpBasis& subspace=spbasis.GetSubspace(i);
      int subspace_dimension=subspace.dimension();
      std::cout<<i<<":  "<<subspace_dimension<<std::endl;
      full_dimension_final+=subspace_dimension;
    } 
    std::cout<<"full dimension final "<<full_dimension_final<<std::endl;
  }
  }



} //main

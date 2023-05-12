/****************************************************************
  unit_tensor.cpp

  Anna E. McCoy
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT
****************************************************************/

#include "spncci/recurrence.h"

#include <omp.h>

#include "fmt/format.h"
#include "lgi/lgi_unit_tensors.h"
#include "mcutils/eigen.h"
#include "lgi/lgi_unit_tensors.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/sp3r_operator.h"
#include "spncci/spncci_common.h"
#include "spncci/transform_basis.h"
#include "utilities/utilities.h"

extern double zero_threshold;

namespace spncci
{

  void
  ZeroInitBlocks(int number, int rows, int cols,std::vector<basis::OperatorBlock<double>>& unit_tensor_blocks)
  {
    unit_tensor_blocks.resize(number);
    for(int i=0; i<number; ++i)
      unit_tensor_blocks[i]=Eigen::MatrixXd::Zero(rows,cols);
  }


 int FromIrrepIndexGetBabySpncciIndex(
    int irrep_family_index,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const spncci::BabySpNCCISpace& baby_spncci_space
    )
    {
      // Extract LGI labels
      const lgi::LGI& lgi=lgi_families[irrep_family_index].irrep;
      u3::U3 sigma;
      HalfInt Sp,Sn,S;
      std::tie(std::ignore,sigma,Sp,Sn,S)=lgi.Key();

      // Get baby spncci index
      spncci::BabySpNCCISubspaceLabels baby_spncci_labels(sigma,Sp,Sn,S,sigma);
      int baby_spncci_index=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labels);
      assert(baby_spncci_space.GetSubspace(baby_spncci_index).irrep_family_index()==irrep_family_index);
      return baby_spncci_index;
    }

void GetLGIUnitTensorSubspaceIndices(
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensors,
  std::set<int>& lgi_operator_subset
  )
  {
    // for each unit tensor, extract labels and identify subspace index for both the tensor and its conjugate
    for(auto unit_tensor : lgi_unit_tensors)
      {
        // Extract unit tensor labels
        u3::SU3 x0;
        HalfInt S0;
        int etap,eta;
        std::tie(x0,S0,std::ignore,etap,std::ignore,std::ignore,eta,std::ignore,std::ignore)=unit_tensor.FlatKey();

        // look up unit tensor index
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
        int operator_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
        lgi_operator_subset.insert(operator_subspace_index);

        // Add conjugate tensor for Nn=0 sectors
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
        int operator_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
        lgi_operator_subset.insert(operator_subspace_index_conj);


      }
  }


int GetUnitTensorSubspaceIndex(
  const u3::SU3& x0, HalfInt S0, int etap, int eta,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  std::set<int>& unit_tensor_subset,
  bool conjugate=false
  )
  // for a given unit tensor subspace, look up its subspace index for both the tensor and add the
  // unit tensor subspace to the accumulating set of operator subspaces appearing in the recurrence
  {
    // Construct unit tensor subspace labels
    u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels;

    if(conjugate)
      unit_tensor_subspace_labels=u3shell::UnitTensorSubspaceLabels(u3::Conjugate(x0),S0,eta,etap);
    else
      unit_tensor_subspace_labels=u3shell::UnitTensorSubspaceLabels(x0,S0,etap,eta);

    // look up subspace index
    int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);

    return unit_tensor_subspace_index;

  }

void GetCase1UnitTensors(
    const u3::SU3& x0, const HalfInt& S0, int etap,int eta,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    std::set<int>& unit_tensor_subset,
    bool conjugate=false
  )
  // Each new unit tensors must satisfy:
  //    x0 x (2,0) -> x0p
  // and
  //    (etap-2,0)x(0,eta) -> x0p
  {
    // std::cout<<"case 1"<<std::endl;

    // Get list of possible x0p values from etap and eta
    MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(u3::SU3(etap-2,0), u3::SU3(0,eta));
    // For each possible x0p
    for(auto& x0p_tagged : x0p_set)
      {
        u3::SU3 x0p(x0p_tagged.irrep);
        // std::cout<<x0p.Str()<<"  "<<etap-2<<"  "<<eta<<"  "<<x0.Str()<<std::endl;

        // If x0 x (2,0) -> x0p doesn't satisfy constraint go to next x0p
        if(u3::OuterMultiplicity(x0,u3::SU3(0,2),x0p)==0)
          continue;

        // get subspace index
        int unit_tensor_subspace_index
              =GetUnitTensorSubspaceIndex(x0p,S0,etap-2,eta,unit_tensor_space,unit_tensor_subset,conjugate);

        // If unit tensor subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(unit_tensor_subspace_index!=-1)
          unit_tensor_subset.insert(unit_tensor_subspace_index);

      }

  }

void GetCase2UnitTensors(
    const u3::SU3& x0, const HalfInt& S0, int etap,int eta,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    std::set<int>& unit_tensor_subset,
    bool conjugate=false
  )
  // Each new unit tensors must satisfy:
  //    x0 x (2,0) -> x0p
  // and
  //    (etap,0)x(0,eta+2) -> x0p
  {
    // case 2
    // std::cout<<"case 2"<<std::endl;
    // Get list of possible x0p values from etap and eta
    MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(u3::SU3(etap,0), u3::SU3(0,eta+2));
    // For each possible x0p
    for(auto& x0p_tagged : x0p_set)
      {
        u3::SU3 x0p(x0p_tagged.irrep);
        // std::cout<<x0p.Str()<<"  "<<etap<<"  "<<eta<<"  "<<x0.Str()<<std::endl;

        // If x0 x (2,0) -> x0p doesn't satisfy constraint go to next x0p
        if(u3::OuterMultiplicity(x0,u3::SU3(0,2),x0p)==0)
          continue;

        // get subspace index
        int unit_tensor_subspace_index
              =GetUnitTensorSubspaceIndex(x0p,S0,etap,eta+2,unit_tensor_space,unit_tensor_subset,conjugate);

        // If unit tensor subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(unit_tensor_subspace_index!=-1)
          unit_tensor_subset.insert(unit_tensor_subspace_index);
      }
  }

void GenerateRecurrenceUnitTensors(
    int Nmax, int N1v,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensors,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    std::map<spncci::NnPair,std::set<int>>& operator_subsets_NnpNn
  )
  // Generate vector of operator subspaces which will appear in the recurrence from
  // unit tensors with non-zero rmes between the given lgi pair
  {
    int Nrel_max=Nmax+2*N1v;

    // Get lgi unit tensor subspaces
    auto& lgi_operator_subset=operator_subsets_NnpNn[spncci::NnPair(0,0)];
    GetLGIUnitTensorSubspaceIndices(unit_tensor_space,lgi_unit_tensors,lgi_operator_subset);

    // Unit tensors subspaces are identified recursively starting from those between the lgi
    //
    // Generate unit tensors for (Nnp,0) and (0,Nn) hypersectors
    for(int Nn=0; Nn<=Nmax; Nn+=2)
      {
        std::set<int>& NnpNn_subspaces_source=operator_subsets_NnpNn[spncci::NnPair(0,Nn)];
        std::set<int>& NnpNn_subspaces_target=operator_subsets_NnpNn[spncci::NnPair(0,Nn+2)];
        std::set<int>& NnpNn_subspaces_target_conj=operator_subsets_NnpNn[spncci::NnPair(Nn+2,0)];

        for(int subspace_index : NnpNn_subspaces_source)
          {
            // Extract source unit tensor labels
            u3::SU3 x0;
            HalfInt S0;
            int etap,eta;
            std::tie(x0,S0,etap,eta)=unit_tensor_space.GetSubspace(subspace_index).labels();

            // Generate unit tensors for (0,Nn)
            bool conjugate=false;

            if(etap-2>=0)
              GetCase1UnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);

            if(eta+2<=Nrel_max)
              GetCase2UnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);


            // Generate unit tensors for (Nnp,0)
            conjugate=true;

            if(etap-2>=0)
              GetCase1UnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target_conj,conjugate);

            if(eta+2<=Nrel_max)
              GetCase2UnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target_conj,conjugate);
          }
      }

    // Generate remaining unit tensors
    for(int Nsum=0; Nsum<=2*Nmax; Nsum+=2)
      for(int Nnp=0; Nnp<=Nmax; Nnp+=2)
        {
          int Nn=Nsum-Nnp;

          // Nn must be positive integer and only need Nn<=Nnp unit tensors
          if(Nn<0 || Nn>Nnp)
            continue;

          std::set<int>& NnpNn_subspaces_source=operator_subsets_NnpNn[spncci::NnPair(Nnp,Nn)];
          // std::cout<<"Nsum group "<<Nsum<<"  "<<Nnp<<"  "<<Nn<<std::endl;

          // (Nnp,Nn)-> (Nnp,Nn+2) sectors
          if(Nn+2<=Nnp)
            {
              std::set<int>& NnpNn_subspaces_target=operator_subsets_NnpNn[spncci::NnPair(Nnp,Nn+2)];
                          // Generate unit tensors for (0,Nn)


              for(int subspace_index : NnpNn_subspaces_source)
                  {
                    // Extract source unit tensor labels
                    u3::SU3 x0;
                    HalfInt S0;
                    int etap,eta;
                    std::tie(x0,S0,etap,eta)=unit_tensor_space.GetSubspace(subspace_index).labels();

                    bool conjugate=false;

                    if(etap-2>0)
                      GetCase1UnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);

                    if(eta+2<=Nrel_max)
                      GetCase2UnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);
                  }
            }

          // Diagonal sectors, i.e., (Nnp,Nn)->(Nnp+2,Nn+2)
          if((Nnp+2<=Nmax) && (Nn+2<=Nmax))
          {
            // std::cout<<"inserting diagonal "<<NnpNn_subspaces_source.size()<<std::endl;
            std::set<int>& NnpNn_subspaces_target=operator_subsets_NnpNn[spncci::NnPair(Nnp+2,Nn+2)];
            NnpNn_subspaces_target.insert(NnpNn_subspaces_source.begin(), NnpNn_subspaces_source.end());
          }
        }// End Nnp loop
  }

  void GetLGIPairsForRecurrence(
      const lgi::MultiplicityTaggedLGIVector& lgi_families,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::SigmaIrrepMap& sigma_irrep_map,
      std::vector<spncci::LGIPair>& lgi_pairs
    )
  {

    // Organize lgi pairs by basis size -- alternative to simple loop over LGI familes above
    std::map<int, std::vector<spncci::LGIPair>, std::greater<int> > sort_map;
    for(int irrep_family_index_bra=0; irrep_family_index_bra<lgi_families.size(); ++irrep_family_index_bra)
      for(int irrep_family_index_ket=0; irrep_family_index_ket<=irrep_family_index_bra; ++irrep_family_index_ket)
        {
          int t=spncci_space[irrep_family_index_bra].Sp3RSpace().size()
                +spncci_space[irrep_family_index_ket].Sp3RSpace().size();
          sort_map[t].push_back(spncci::LGIPair(irrep_family_index_bra,irrep_family_index_ket));
        }
    // std::cout<<"sorted map"<<std::endl;
    for(auto it=sort_map.begin(); it!=sort_map.end(); ++it)
      {
        std::cout<<it->first<<std::endl;
        for(const auto& pair : it->second)
          lgi_pairs.push_back(pair);
      }

    // for(auto pair: lgi_pairs)
    //   std::cout<<pair.first<<"  "<<pair.second<<std::endl;

  }


void AddLGIPairNnSorted(
      const std::vector<int>& lgi_full_space_index_lookup,
      const spncci::SpNCCISpace& spncci_space,
      int Nmax,
      int irrep_family_index_bra,
      int irrep_family_index_ket,
      std::map<int, std::vector<spncci::LGIPair>, std::greater<int> >& sort_map
    )
  {
    int index1=lgi_full_space_index_lookup[irrep_family_index_bra];
    int index2=lgi_full_space_index_lookup[irrep_family_index_ket];

    std::string lgi_unit_tensor_filename
      =fmt::format("seeds/operators_{:06d}_{:06d}.dat",index1,index2);

    bool verbose=false;
    bool files_found=utils::FileExists(lgi_unit_tensor_filename,verbose);
    if(files_found)
      {
        int t=(Nmax-spncci_space[irrep_family_index_bra].Nex())
                +(Nmax-spncci_space[irrep_family_index_ket].Nex());
        sort_map[t].push_back(spncci::LGIPair(irrep_family_index_bra,irrep_family_index_ket));
      }

  }


void GetLGIPairsForRecurrence(
      // const lgi::MultiplicityTaggedLGIVector& lgi_families,
      const std::vector<int>& lgi_full_space_index_lookup,
      const spncci::SpNCCISpace& spncci_space,
      // const spncci::SigmaIrrepMap& sigma_irrep_map,
      int Nmax,
      std::vector<spncci::LGIPair>& lgi_pairs
    )
  {

    // Organize lgi pairs by basis size -- alternative to simple loop over LGI familes above
    std::map<int, std::vector<spncci::LGIPair>, std::greater<int> > sort_map;

    int num_irrep_families=lgi_full_space_index_lookup.size();
    for(int irrep_family_index_bra=0; irrep_family_index_bra<num_irrep_families; ++irrep_family_index_bra)
      for(int irrep_family_index_ket=0; irrep_family_index_ket<=irrep_family_index_bra; ++irrep_family_index_ket)
        {
          spncci::AddLGIPairNnSorted(
            lgi_full_space_index_lookup,spncci_space,Nmax,
            irrep_family_index_bra,irrep_family_index_ket,sort_map
          );
        }

    // std::cout<<"sorted map"<<std::endl;
    for(auto it=sort_map.begin(); it!=sort_map.end(); ++it)
      {
        // std::cout<<it->first<<std::endl;

        for(const auto& pair : it->second)
          lgi_pairs.push_back(pair);
      }
  }





void GetLGIPairsForRecurrence(
      const std::vector<int>& lgi_full_space_index_lookup,
      const spncci::SpNCCISpace& spncci_space,
      int Nmax,
      const std::vector<int>& trial_subspace,
      const std::vector<int>& test_subspace,
      std::vector<spncci::LGIPair>& lgi_pairs
    )
  {
    // number of irrep famlies in space
    int num_irrep_families=lgi_full_space_index_lookup.size();

    // Organize lgi pairs by basis size -- alternative to simple loop over LGI familes above
    std::map<int, std::vector<spncci::LGIPair>, std::greater<int> > sort_map;

    //Create pairs for calculating RMEs all irrep family pairs in trial subspace
    for(int irrep_family_index_ket : trial_subspace)
      for(int irrep_family_index_bra : trial_subspace)
        {
          if(irrep_family_index_ket<=irrep_family_index_bra)
            spncci::AddLGIPairNnSorted(
                lgi_full_space_index_lookup,spncci_space,Nmax,
                irrep_family_index_bra,irrep_family_index_bra,sort_map
              );
        }
    //Create pairs for test subspace with trial subpace
    for(int irrep_family_index_ket : trial_subspace)
      for(int irrep_family_index_bra : test_subspace)
        {
          if(irrep_family_index_ket<=irrep_family_index_bra)
            spncci::AddLGIPairNnSorted(
                lgi_full_space_index_lookup,spncci_space,Nmax,
                irrep_family_index_bra,irrep_family_index_bra,sort_map
              );
          else
            spncci::AddLGIPairNnSorted(
                lgi_full_space_index_lookup,spncci_space,Nmax,
                irrep_family_index_bra,irrep_family_index_bra,sort_map
              );
        }

    // std::cout<<"sort the lgi pairs"<<std::endl;
    for(auto it=sort_map.begin(); it!=sort_map.end(); ++it)
      {
        std::cout<<it->first<<std::endl;
        for(const auto& pair : it->second)
          lgi_pairs.push_back(pair);
      }
  }

  void
  PopulateHypersectorsWithSeeds(
    int irrep_family_index_bra, int irrep_family_index_ket,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors_Nn0,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensors,
    const std::vector<int>& rho0_values,
    basis::OperatorBlocks<double>& unit_tensor_seed_blocks,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_Nn0,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  )
  {
    int baby_spncci_index_bra
      =spncci::FromIrrepIndexGetBabySpncciIndex(irrep_family_index_bra,lgi_families,baby_spncci_space);
    int baby_spncci_index_ket
      =spncci::FromIrrepIndexGetBabySpncciIndex(irrep_family_index_ket,lgi_families,baby_spncci_space);


    const lgi::LGI& lgi_bra=lgi_families[irrep_family_index_bra].irrep;
    const lgi::LGI& lgi_ket=lgi_families[irrep_family_index_ket].irrep;
    const u3::U3& sigmap=lgi_bra.sigma();
    const u3::U3& sigma=lgi_ket.sigma();

    double conjugation_factor_base
        =ParitySign(u3::ConjugationGrade(sigmap)+lgi_bra.S()-u3::ConjugationGrade(sigma)-lgi_ket.S())
          * sqrt(1.*u3::dim(sigmap)*am::dim(lgi_bra.S())/u3::dim(sigma)/am::dim(lgi_ket.S()));

    // std::cout<<"loop over lgi unit tensors"<<std::endl;
    for(int i=0; i<lgi_unit_tensors.size();  ++i)
      {
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor=lgi_unit_tensors[i];
        int rho0=rho0_values[i];

        // Extract unit tensor labels
        u3::SU3 x0;
        HalfInt S0,T0, Sp,Tp,S,T;
        int etap, eta;
        std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=unit_tensor.FlatKey();
        // std::cout<<sigmap.Str()<<" "<<lgi_bra.S()<<"  "<<sigma.Str()<<"  "<<lgi_ket.S()<<" "<<unit_tensor.Str()<<std::endl;

        // Look up unit tensor subspace
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
        int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
        auto& subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        // Look up index of unit tensor in subspace
        int unit_tensor_state_index
              =subspace.LookUpStateIndex(std::tuple<int,int,int,int,int>(int(T0),int(Sp),int(Tp),int(S),int(T)));

        // Get Hypersector index
        int hypersector_index
            =baby_spncci_hypersectors.LookUpHypersectorIndex(
                baby_spncci_index_bra,baby_spncci_index_ket,
                unit_tensor_subspace_index,rho0
              );

        // std::cout<<hypersector_index<<"  "<<unit_tensor_state_index<<"  "<<i<<std::endl
        //   <<unit_tensor_seed_blocks[i]<<std::endl;
        unit_tensor_hyperblocks[hypersector_index][unit_tensor_state_index]=unit_tensor_seed_blocks[i];

        // Get conjugate

        // Look up conjugate unit tensor subspace
        // std::cout<<x0.Str()<<"  "<<S0<<"  "<<eta<<"  "<<etap<<std::endl;
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
        int unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
        auto& subspace_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj);

        // Look up index of unit tensor in subspace
        int unit_tensor_state_index_conj
              =subspace_conj.LookUpStateIndex(std::tuple<int,int,int,int,int>(int(T0),int(S),int(T),int(Sp),int(Tp)));

        // Get Hypersector index
        int hypersector_index_Nn0
            =baby_spncci_hypersectors_Nn0.LookUpHypersectorIndex(
                baby_spncci_index_ket,baby_spncci_index_bra,
                unit_tensor_subspace_index_conj,rho0
              );

        double conjugation_factor
        =conjugation_factor_base
          *sqrt(1.*u3::dim(u3::SU3(eta,0))*am::dim(S)*am::dim(T)/u3::dim(u3::SU3(etap,0))/am::dim(Sp)/am::dim(Tp));
        // std::cout<<"  "<<unit_tensor_state_index_conj<<"  "<<i<<"  "<<conjugation_factor<<std::endl
          // <<unit_tensor_seed_blocks[i]<<std::endl;
        unit_tensor_hyperblocks_Nn0[hypersector_index_Nn0][unit_tensor_state_index_conj]
          =conjugation_factor*unit_tensor_seed_blocks[i].transpose();

      }
    }

//************************************** Added by J.H. *******************************
int GetOneBodyUnitTensorSubspaceIndex(
  const u3::SU3& x0, HalfInt S0, int etap, int eta,
  const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
  bool conjugate=false
  )
  // for a given unit tensor subspace, look up its subspace index
  {
    // Construct unit tensor subspace labels
    u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels;

    if(conjugate)
      unit_tensor_subspace_labels=u3shell::UnitTensorSubspaceLabels(u3::Conjugate(x0),S0,eta,etap);
    else
      unit_tensor_subspace_labels=u3shell::UnitTensorSubspaceLabels(x0,S0,etap,eta);

    // look up subspace index
    int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);

    return unit_tensor_subspace_index;
  }

void GetCase1OneBodyUnitTensors(
    const u3::SU3& x0, const HalfInt& S0, int etap,int eta,
    const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
    std::set<int>& unit_tensor_subset,
    bool conjugate=false
  )
  // Each new unit tensors must satisfy:
  //    (etap-2,0)x(0,eta) -> x0p
  // and
  //    x0p x (2,0) -> x0
  {
    // Get list of possible x0p values from etap and eta
    MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(u3::SU3(etap-2,0), u3::SU3(0,eta));
    // For each possible x0p
    for(auto& x0p_tagged : x0p_set)
      {
        u3::SU3 x0p(x0p_tagged.irrep);

        // If x0p x (2,0) -> x0 doesn't satisfy constraint go to next x0p
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)
          continue;

        // get subspace index
        int unit_tensor_subspace_index
              =GetOneBodyUnitTensorSubspaceIndex(x0p,S0,etap-2,eta,unit_tensor_space,conjugate);

        // If unit tensor subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(unit_tensor_subspace_index!=-1)
          unit_tensor_subset.insert(unit_tensor_subspace_index);
      }
  }

void GetCase2OneBodyUnitTensors(
    const u3::SU3& x0, const HalfInt& S0, int etap,int eta,
    const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
    std::set<int>& unit_tensor_subset,
    bool conjugate=false
  )
  // Each new unit tensors must satisfy:
  //    (etap,0)x(0,eta+2) -> x0p
  // and
  //    x0p x (2,0) -> x0
  {
    // Get list of possible x0p values from etap and eta
    MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(u3::SU3(etap,0), u3::SU3(0,eta+2));
    // For each possible x0p
    for(auto& x0p_tagged : x0p_set)
      {
        u3::SU3 x0p(x0p_tagged.irrep);

        // If x0p x (2,0) -> x0 doesn't satisfy constraint go to next x0p
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)
          continue;

        // get subspace index
        int unit_tensor_subspace_index
              =GetOneBodyUnitTensorSubspaceIndex(x0p,S0,etap,eta+2,unit_tensor_space,conjugate);

        // If unit tensor subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(unit_tensor_subspace_index!=-1)
          unit_tensor_subset.insert(unit_tensor_subspace_index);
      }
  }

void GetCase3OneBodyUnitTensors(
    const u3::SU3& x0, const HalfInt& S0, int etap,int eta,
    const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
    std::set<int>& unit_tensor_subset,
    bool conjugate=false
  )
  // Each new unit tensors must satisfy:
  //    (etap-1,0)x(0,eta+1) -> x0p
  // and
  //    x0p x (2,0) -> x0
  {
    // Get list of possible x0p values from etap and eta
    MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(u3::SU3(etap-1,0), u3::SU3(0,eta+1));
    // For each possible x0p
    for(auto& x0p_tagged : x0p_set)
      {
        u3::SU3 x0p(x0p_tagged.irrep);

        // If x0p x (2,0) -> x0 doesn't satisfy constraint go to next x0p
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)
          continue;

        // get subspace index
        int unit_tensor_subspace_index
              =GetOneBodyUnitTensorSubspaceIndex(x0p,S0,etap-1,eta+1,unit_tensor_space,conjugate);

        // If unit tensor subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(unit_tensor_subspace_index!=-1)
          unit_tensor_subset.insert(unit_tensor_subspace_index);
      }
  }

int GetTwoBodyDensitySubspaceIndex(
  const u3::SU3& x0, HalfInt S0, int N1, int N2, int N3, int N4,
  const u3shell::TwoBodyDensitySpace& tbd_space,
  bool conjugate=false
  )
  // for a given TBD subspace, look up its subspace index
  {
    // Construct TBD subspace labels
    u3shell::TwoBodyDensitySubspaceLabels tbd_subspace_labels;

    if(conjugate)
      tbd_subspace_labels=u3shell::TwoBodyDensitySubspaceLabels(u3::Conjugate(x0),S0,N4,N3,N2,N1);
    else
      tbd_subspace_labels=u3shell::TwoBodyDensitySubspaceLabels(x0,S0,N1,N2,N3,N4);

    // look up subspace index
    int tbd_subspace_index=tbd_space.LookUpSubspaceIndex(tbd_subspace_labels);

    return tbd_subspace_index;
  }

void GetCase1TwoBodyDensities(const u3::SU3& x0,const HalfInt& S0,int N1,int N2,int N3,int N4,
    const u3shell::TwoBodyDensitySpace& tbd_space, std::set<int>& tbd_subset, bool conjugate=false){
  // Each new TBD must satisfy:
  //    (N1,0) x (N2,0) -> xf
  //    (0,N3)x(0,N4+2) -> xi
  //    xf x xi -> x0p
  //    x0p x (2,0) -> x0
  for(MultiplicityTagged<u3::SU3>& xf_tagged : KroneckerProduct(u3::SU3(N1,0), u3::SU3(N2,0))){
    for(MultiplicityTagged<u3::SU3>& xi_tagged : KroneckerProduct(u3::SU3(0,N3), u3::SU3(0,N4+2))){
      for(MultiplicityTagged<u3::SU3>& x0p_tagged : KroneckerProduct(xf_tagged.irrep, xi_tagged.irrep)){
        u3::SU3 x0p(x0p_tagged.irrep);
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)continue;
        int tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1,N2,N3,N4+2,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1,N2,N4+2,N3,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
      }
    }
  }
}

void GetCase2TwoBodyDensities(const u3::SU3& x0,const HalfInt& S0,int N1,int N2,int N3,int N4,
    const u3shell::TwoBodyDensitySpace& tbd_space, std::set<int>& tbd_subset, bool conjugate=false){
  // Each new TBD must satisfy:
  //    (N2,0) x (N1-2,0) -> xf
  //    (0,N3)x(0,N4) -> xi
  //    xf x xi -> x0p
  //    x0p x (2,0) -> x0
  for(MultiplicityTagged<u3::SU3>& xf_tagged : KroneckerProduct(u3::SU3(N2,0), u3::SU3(N1-2,0))){
    for(MultiplicityTagged<u3::SU3>& xi_tagged : KroneckerProduct(u3::SU3(0,N3), u3::SU3(0,N4))){
      for(MultiplicityTagged<u3::SU3>& x0p_tagged : KroneckerProduct(xf_tagged.irrep, xi_tagged.irrep)){
        u3::SU3 x0p(x0p_tagged.irrep);
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)continue;
        int tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N2,N1-2,N3,N4,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1-2,N2,N3,N4,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
      }
    }
  }
}

void GetCase3TwoBodyDensities(const u3::SU3& x0,const HalfInt& S0,int N1,int N2,int N3,int N4,
    const u3shell::TwoBodyDensitySpace& tbd_space, std::set<int>& tbd_subset, bool conjugate=false){
  // Each new TBD must satisfy:
  //    (N1,0) x (N2,0) -> xf
  //    (0,N3+1)x(0,N4+1) -> xi
  //    xf x xi -> x0p
  //    x0p x (2,0) -> x0
  for(MultiplicityTagged<u3::SU3>& xf_tagged : KroneckerProduct(u3::SU3(N1,0), u3::SU3(N2,0))){
    for(MultiplicityTagged<u3::SU3>& xi_tagged : KroneckerProduct(u3::SU3(0,N3+1), u3::SU3(0,N4+1))){
      for(MultiplicityTagged<u3::SU3>& x0p_tagged : KroneckerProduct(xf_tagged.irrep, xi_tagged.irrep)){
        u3::SU3 x0p(x0p_tagged.irrep);
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)continue;
        int tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1,N2,N3+1,N4+1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1,N2,N4+1,N3+1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
      }
    }
  } 
}

void GetCase4TwoBodyDensities(const u3::SU3& x0,const HalfInt& S0,int N1,int N2,int N3,int N4,
    const u3shell::TwoBodyDensitySpace& tbd_space, std::set<int>& tbd_subset, bool conjugate=false){
  // Each new TBD must satisfy:
  //    (N2,0) x (N1-1,0) -> xf
  //    (0,N4+1)x(0,N3) -> xi
  //    xf x xi -> x0p
  //    x0p x (2,0) -> x0
  for(MultiplicityTagged<u3::SU3>& xf_tagged : KroneckerProduct(u3::SU3(N2,0), u3::SU3(N1-1,0))){
    for(MultiplicityTagged<u3::SU3>& xi_tagged : KroneckerProduct(u3::SU3(0,N4+1), u3::SU3(0,N3))){
      for(MultiplicityTagged<u3::SU3>& x0p_tagged : KroneckerProduct(xf_tagged.irrep, xi_tagged.irrep)){
        u3::SU3 x0p(x0p_tagged.irrep);
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)continue;
        int tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N2,N1-1,N4+1,N3,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1-1,N2,N4+1,N3,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N2,N1-1,N3,N4+1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1-1,N2,N3,N4+1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
      }
    }
  }
}

void GetCase5TwoBodyDensities(const u3::SU3& x0,const HalfInt& S0,int N1,int N2,int N3,int N4,
    const u3shell::TwoBodyDensitySpace& tbd_space, std::set<int>& tbd_subset, bool conjugate=false){
  // Each new TBD must satisfy:
  //    (N2,0) x (N1-1,0) -> xf
  //    (0,N4+2)x(0,N3-1) -> xi
  //    xf x xi -> x0p
  //    x0p x (2,0) -> x0
  for(MultiplicityTagged<u3::SU3>& xf_tagged : KroneckerProduct(u3::SU3(N2,0), u3::SU3(N1-1,0))){
    for(MultiplicityTagged<u3::SU3>& xi_tagged : KroneckerProduct(u3::SU3(0,N4+2), u3::SU3(0,N3-1))){
      for(MultiplicityTagged<u3::SU3>& x0p_tagged : KroneckerProduct(xf_tagged.irrep, xi_tagged.irrep)){
        u3::SU3 x0p(x0p_tagged.irrep);
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)continue;
        int tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N2,N1-1,N4+2,N3-1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1-1,N2,N4+2,N3-1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);

      }
    }
  }
}

void GetCase6TwoBodyDensities(const u3::SU3& x0,const HalfInt& S0,int N1,int N2,int N3,int N4,
    const u3shell::TwoBodyDensitySpace& tbd_space, std::set<int>& tbd_subset, bool conjugate=false){
  // Each new TBD must satisfy:
  //    (N2,0) x (N1-1,0) -> xf
  //    (0,N3+1)x(0,N4) -> xi
  //    xf x xi -> x0p
  //    x0p x (2,0) -> x0
  for(MultiplicityTagged<u3::SU3>& xf_tagged : KroneckerProduct(u3::SU3(N2,0), u3::SU3(N1-1,0))){
    for(MultiplicityTagged<u3::SU3>& xi_tagged : KroneckerProduct(u3::SU3(0,N3+1), u3::SU3(0,N4))){
      for(MultiplicityTagged<u3::SU3>& x0p_tagged : KroneckerProduct(xf_tagged.irrep, xi_tagged.irrep)){
        u3::SU3 x0p(x0p_tagged.irrep);
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)continue;
        int tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N2,N1-1,N3+1,N4,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1-1,N2,N3+1,N4,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
      }
    }
  } 
}

void GetCase7TwoBodyDensities(const u3::SU3& x0,const HalfInt& S0,int N1,int N2,int N3,int N4,
    const u3shell::TwoBodyDensitySpace& tbd_space, std::set<int>& tbd_subset, bool conjugate=false){
  // Each new TBD must satisfy:
  //    (N1-1,0) x (N2-1,0) -> xf
  //    (0,N4+1)x(0,N3-1) -> xi
  //    xf x xi -> x0p
  //    x0p x (2,0) -> x0
  for(MultiplicityTagged<u3::SU3>& xf_tagged : KroneckerProduct(u3::SU3(N1-1,0), u3::SU3(N2-1,0))){
    for(MultiplicityTagged<u3::SU3>& xi_tagged : KroneckerProduct(u3::SU3(0,N4+1), u3::SU3(0,N3-1))){
      for(MultiplicityTagged<u3::SU3>& x0p_tagged : KroneckerProduct(xf_tagged.irrep, xi_tagged.irrep)){
        u3::SU3 x0p(x0p_tagged.irrep);
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)continue;
        int tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1-1,N2-1,N4+1,N3-1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N2-1,N1-1,N4+1,N3-1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
      }
    }
  } 
}

void GetCase8TwoBodyDensities(const u3::SU3& x0,const HalfInt& S0,int N1,int N2,int N3,int N4,
    const u3shell::TwoBodyDensitySpace& tbd_space, std::set<int>& tbd_subset, bool conjugate=false){
  // Each new TBD must satisfy:
  //    (N2,0) x (N1-2,0) -> xf
  //    (0,N4+1)x(0,N3-1) -> xi
  //    xf x xi -> x0p
  //    x0p x (2,0) -> x0
  for(MultiplicityTagged<u3::SU3>& xf_tagged : KroneckerProduct(u3::SU3(N2,0), u3::SU3(N1-2,0))){
    for(MultiplicityTagged<u3::SU3>& xi_tagged : KroneckerProduct(u3::SU3(0,N4+1), u3::SU3(0,N3-1))){
      for(MultiplicityTagged<u3::SU3>& x0p_tagged : KroneckerProduct(xf_tagged.irrep, xi_tagged.irrep)){
        u3::SU3 x0p(x0p_tagged.irrep);
        if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)==0)continue;
        int tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N2,N1-2,N4+1,N3-1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
	tbd_subspace_index=GetTwoBodyDensitySubspaceIndex(x0p,S0,N1-2,N2,N4+1,N3-1,tbd_space,conjugate);
        // If TBD subspace exists, add to operator_subsets for given Nnp,Nn sector
        if(tbd_subspace_index!=-1)tbd_subset.insert(tbd_subspace_index);
      }
    }
  }
}

void GetLGIOneBodyUnitTensorSubspaceIndices(
  const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
  const std::vector<u3shell::OneBodyUnitTensorLabelsU3S>& lgi_unit_tensors,
  std::set<int>& lgi_operator_subset
  )
  {
    // for each unit tensor, extract labels and identify subspace index for both the tensor and its conjugate
    for(auto unit_tensor : lgi_unit_tensors)
      {
        // Extract unit tensor labels
        u3::SU3 x0;
        HalfInt S0;
        int etap,eta;
        std::tie(x0,S0,etap,eta,std::ignore)=unit_tensor.FlatKey();

        // look up unit tensor index
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
        int operator_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
        lgi_operator_subset.insert(operator_subspace_index);

        // Add conjugate tensor for Nn=0 sectors
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
        int operator_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
        lgi_operator_subset.insert(operator_subspace_index_conj);
      }
  }

void GetLGITwoBodyDensitySubspaceIndices(
  const u3shell::TwoBodyDensitySpace& tbd_space,
  const std::vector<u3shell::TwoBodyDensityLabels>& lgi_tbds,
  std::set<int>& lgi_operator_subset
  )
  {
    // for each TBD, extract labels and identify subspace index for both the tensor and its conjugate
    for(auto tbd : lgi_tbds)
      {
        // Extract TBD labels
        u3::SU3 x0;
        HalfInt S0;
        int N1,N2,N3,N4;
        std::tie(x0,S0,N1,N2,N3,N4,std::ignore,std::ignore,std::ignore,std::ignore,std::ignore,std::ignore)=tbd.FlatKey();

        // look up TBD index
        u3shell::TwoBodyDensitySubspaceLabels tbd_subspace_labels(x0,S0,N1,N2,N3,N4);
        int operator_subspace_index=tbd_space.LookUpSubspaceIndex(tbd_subspace_labels);
        lgi_operator_subset.insert(operator_subspace_index);

        // Add conjugate tensor for Nn=0 sectors
        u3shell::TwoBodyDensitySubspaceLabels tbd_subspace_labels_conj(u3::Conjugate(x0),S0,N4,N3,N2,N1);
        int operator_subspace_index_conj=tbd_space.LookUpSubspaceIndex(tbd_subspace_labels_conj);
        lgi_operator_subset.insert(operator_subspace_index_conj);
      }
  }

void GenerateRecurrenceOneBodyUnitTensors(
    int Nmax, int N1vp, int N1vn,
    const std::vector<u3shell::OneBodyUnitTensorLabelsU3S>& lgi_unit_tensors,
    const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
    std::map<spncci::NnPair,std::set<int>>& operator_subsets_NnpNn
  )
  // Generate vector of operator subspaces which will appear in the recurrence from
  // unit tensors with non-zero rmes between the given lgi pair
  {

    int eta_max=Nmax+std::max(N1vp,N1vn);

    // Get lgi unit tensor subspaces
    auto& lgi_operator_subset=operator_subsets_NnpNn[spncci::NnPair(0,0)];
    GetLGIOneBodyUnitTensorSubspaceIndices(unit_tensor_space,lgi_unit_tensors,lgi_operator_subset);

    // Unit tensors subspaces are identified recursively starting from those between the lgi
    //
    // Generate unit tensors for (Nnp,0) and (0,Nn) hypersectors
    for(int Nn=0; Nn<=Nmax; Nn+=2)
      {
        std::set<int>& NnpNn_subspaces_source=operator_subsets_NnpNn[spncci::NnPair(0,Nn)];
        std::set<int>& NnpNn_subspaces_target=operator_subsets_NnpNn[spncci::NnPair(0,Nn+2)];
        std::set<int>& NnpNn_subspaces_target_conj=operator_subsets_NnpNn[spncci::NnPair(Nn+2,0)];

        for(int subspace_index : NnpNn_subspaces_source)
          {
            // Extract source unit tensor labels
            u3::SU3 x0;
            HalfInt S0;
            int etap,eta;
            std::tie(x0,S0,etap,eta)=unit_tensor_space.GetSubspace(subspace_index).labels();

            // Generate unit tensors for (0,Nn)
            bool conjugate=false;

            if(etap-2>=0)
              GetCase1OneBodyUnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);

            if(eta+2<=eta_max)
              GetCase2OneBodyUnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);

	    if(etap-1>=0 && eta+1<=eta_max)
              GetCase3OneBodyUnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);


            // Generate unit tensors for (Nnp,0)
            conjugate=true;

            if(etap-2>=0)
              GetCase1OneBodyUnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target_conj,conjugate);

            if(eta+2<=eta_max)
              GetCase2OneBodyUnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target_conj,conjugate);

	    if(etap-1>=0 && eta+1<=eta_max)
              GetCase3OneBodyUnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target_conj,conjugate);
          }
      }

    // Generate remaining unit tensors
    for(int Nsum=0; Nsum<=2*Nmax; Nsum+=2)
      for(int Nnp=0; Nnp<=Nmax; Nnp+=2)
        {
          int Nn=Nsum-Nnp;

          // Nn must be positive integer and only need Nn<=Nnp unit tensor RMEs
          if(Nn<0 || Nn>Nnp)
            continue;

          std::set<int>& NnpNn_subspaces_source=operator_subsets_NnpNn[spncci::NnPair(Nnp,Nn)];

          // (Nnp,Nn)-> (Nnp,Nn+2) sectors
          if(Nn+2<=Nnp)
            {
              std::set<int>& NnpNn_subspaces_target=operator_subsets_NnpNn[spncci::NnPair(Nnp,Nn+2)];

              for(int subspace_index : NnpNn_subspaces_source)
                  {
                    // Extract source unit tensor labels
                    u3::SU3 x0;
                    HalfInt S0;
                    int etap,eta;
                    std::tie(x0,S0,etap,eta)=unit_tensor_space.GetSubspace(subspace_index).labels();

                    bool conjugate=false;

                    if(etap-2>=0)
                      GetCase1OneBodyUnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);

                    if(eta+2<=eta_max)
                      GetCase2OneBodyUnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);

                    if(etap-1>=0 && eta+1<=eta_max)
                      GetCase3OneBodyUnitTensors(x0,S0,etap,eta,unit_tensor_space,NnpNn_subspaces_target,conjugate);
                  }
            }

          // Diagonal sectors, i.e., (Nnp,Nn)->(Nnp+2,Nn+2)
          if((Nnp+2<=Nmax) && (Nn+2<=Nmax))
          {
            std::set<int>& NnpNn_subspaces_target=operator_subsets_NnpNn[spncci::NnPair(Nnp+2,Nn+2)];
            NnpNn_subspaces_target.insert(NnpNn_subspaces_source.begin(), NnpNn_subspaces_source.end());
          }
        }// End Nnp loop
  }

void GenerateRecurrenceTwoBodyDensities(
    int Nmax, int N1vp, int N1vn,
    const std::vector<u3shell::TwoBodyDensityLabels>& lgi_tbds,
    const u3shell::TwoBodyDensitySpace& tbd_space,
    std::map<spncci::NnPair,std::set<int>>& operator_subsets_NnpNn
  )
  // Generate vector of operator subspaces which will appear in the recurrence from
  // TBDs with non-zero rmes between the given lgi pair
  {

    int eta_max=Nmax+std::max(N1vp,N1vn);

    // Get lgi TBD subspaces
    auto& lgi_operator_subset=operator_subsets_NnpNn[spncci::NnPair(0,0)];
    GetLGITwoBodyDensitySubspaceIndices(tbd_space,lgi_tbds,lgi_operator_subset);

    // TBDs subspaces are identified recursively starting from those between the lgi
    //
    // Generate TBDs for (Nnp,0) and (0,Nn) hypersectors
    for(int Nn=0; Nn<=Nmax; Nn+=2)
      {
        std::set<int>& NnpNn_subspaces_source=operator_subsets_NnpNn[spncci::NnPair(0,Nn)];
        std::set<int>& NnpNn_subspaces_target=operator_subsets_NnpNn[spncci::NnPair(0,Nn+2)];
        std::set<int>& NnpNn_subspaces_target_conj=operator_subsets_NnpNn[spncci::NnPair(Nn+2,0)];

        for(int subspace_index : NnpNn_subspaces_source)
          {
            // Extract source TBD labels
            u3::SU3 x0;
            HalfInt S0;
            int N1,N2,N3,N4;
            std::tie(x0,S0,N1,N2,N3,N4)=tbd_space.GetSubspace(subspace_index).labels();

            // Generate TBDs for (0,Nn)
            bool conjugate=false;

            if(N4+2<=eta_max)GetCase1TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-2>=0)GetCase2TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
	    if(N3+1<=eta_max && N4+1<=eta_max)GetCase3TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
	    if(N1-1>=0 && N4+1<=eta_max)GetCase4TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
	    if(N1-1>=0 && N3-1>=0 && N4+2<=eta_max)GetCase5TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
	    if(N1-1>=0 && N3+1<=eta_max)GetCase6TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
	    if(N1-1>=0 && N2-1>=0 && N3-1>=0 && N4+1<=eta_max)GetCase7TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
	    if(N1-2>=0 && N3-1>=0 && N4+1<=eta_max)GetCase8TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);

            // Generate TBDs for (Nnp,0)
            conjugate=true;

            if(N4+2<=eta_max)GetCase1TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-2>=0)GetCase2TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N3+1<=eta_max && N4+1<=eta_max)GetCase3TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-1>=0 && N4+1<=eta_max)GetCase4TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-1>=0 && N3-1>=0 && N4+2<=eta_max)GetCase5TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-1>=0 && N3+1<=eta_max)GetCase6TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
     if(N1-1>=0 && N2-1>=0 && N3-1>=0 && N4+1<=eta_max)GetCase7TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-2>=0 && N3-1>=0 && N4+1<=eta_max)GetCase8TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);

   	 }
    }

    // Generate remaining TBDs
    for(int Nsum=0; Nsum<=2*Nmax; Nsum+=2)
      for(int Nnp=0; Nnp<=Nmax; Nnp+=2)
        {
          int Nn=Nsum-Nnp;

          // Nn must be positive integer and only need Nn<=Nnp TBD RMEs
          if(Nn<0 || Nn>Nnp)continue;

          std::set<int>& NnpNn_subspaces_source=operator_subsets_NnpNn[spncci::NnPair(Nnp,Nn)];

          // (Nnp,Nn)-> (Nnp,Nn+2) sectors
          if(Nn+2<=Nnp)
            {
              std::set<int>& NnpNn_subspaces_target=operator_subsets_NnpNn[spncci::NnPair(Nnp,Nn+2)];

              for(int subspace_index : NnpNn_subspaces_source)
                  {
                    // Extract source TBD labels
                    u3::SU3 x0;
                    HalfInt S0;
                    int N1,N2,N3,N4;
                    std::tie(x0,S0,N1,N2,N3,N4)=tbd_space.GetSubspace(subspace_index).labels();

                    bool conjugate=false;

                    if(N4+2<=eta_max)GetCase1TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-2>=0)GetCase2TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N3+1<=eta_max && N4+1<=eta_max)GetCase3TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-1>=0 && N4+1<=eta_max)GetCase4TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-1>=0 && N3-1>=0 && N4+2<=eta_max)GetCase5TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-1>=0 && N3+1<=eta_max)GetCase6TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
     if(N1-1>=0 && N2-1>=0 && N3-1>=0 && N4+1<=eta_max)GetCase7TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
            if(N1-2>=0 && N3-1>=0 && N4+1<=eta_max)GetCase8TwoBodyDensities(x0,S0,N1,N2,N3,N4,tbd_space,NnpNn_subspaces_target,conjugate);
                  }
            }

          // Diagonal sectors, i.e., (Nnp,Nn)->(Nnp+2,Nn+2)
          if((Nnp+2<=Nmax) && (Nn+2<=Nmax))
          {
            std::set<int>& NnpNn_subspaces_target=operator_subsets_NnpNn[spncci::NnPair(Nnp+2,Nn+2)];
            NnpNn_subspaces_target.insert(NnpNn_subspaces_source.begin(), NnpNn_subspaces_source.end());
          }
        }// End Nnp loop

  }

void AddLGIPairForOneBodyRecurrence(
      const std::vector<int>& lgi_full_space_index_lookup,
      const spncci::SpNCCISpace& spncci_space,
      int Nmax,
      int irrep_family_index_bra,
      int irrep_family_index_ket,
      std::map<int, std::vector<spncci::LGIPair>, std::greater<int> >& sort_map
    )
  {
    int index1=lgi_full_space_index_lookup[irrep_family_index_bra];
    int index2=lgi_full_space_index_lookup[irrep_family_index_ket];

    std::string lgi_unit_tensor_filename
      =fmt::format("seeds/oboperators_{:06d}_{:06d}.dat",index1,index2);

    bool verbose=false;
    bool files_found=utils::FileExists(lgi_unit_tensor_filename,verbose);
    if(files_found)
      {
        int t=(Nmax-spncci_space[irrep_family_index_bra].Nex())
                +(Nmax-spncci_space[irrep_family_index_ket].Nex());
        sort_map[t].push_back(spncci::LGIPair(irrep_family_index_bra,irrep_family_index_ket));
      }

  }

void AddLGIPairForTwoBodyRecurrence(
      const std::vector<int>& lgi_full_space_index_lookup,
      const spncci::SpNCCISpace& spncci_space,
      int Nmax,
      int irrep_family_index_bra,
      int irrep_family_index_ket,
      std::map<int, std::vector<spncci::LGIPair>, std::greater<int> >& sort_map
    )
  {
    int index1=lgi_full_space_index_lookup[irrep_family_index_bra];
    int index2=lgi_full_space_index_lookup[irrep_family_index_ket];

    std::string lgi_tbd_filename
      =fmt::format("seeds/tboperators_{:06d}_{:06d}.dat",index1,index2);

    bool verbose=false;
    bool files_found=utils::FileExists(lgi_tbd_filename,verbose);
    if(files_found)
      {
        int t=(Nmax-spncci_space[irrep_family_index_bra].Nex())
                +(Nmax-spncci_space[irrep_family_index_ket].Nex());
        sort_map[t].push_back(spncci::LGIPair(irrep_family_index_bra,irrep_family_index_ket));
      }

  }

void GetLGIPairsForOneBodyRecurrence(
      const std::vector<int>& lgi_full_space_index_lookup,
      const spncci::SpNCCISpace& spncci_space,
      int Nmax,
      std::vector<spncci::LGIPair>& lgi_pairs_ob
    )
  {
    // Organize lgi pairs by basis size -- alternative to simple loop over LGI familes above
    std::map<int, std::vector<spncci::LGIPair>, std::greater<int> > sort_map;

    int num_irrep_families=lgi_full_space_index_lookup.size();
    for(int irrep_family_index_bra=0; irrep_family_index_bra<num_irrep_families; ++irrep_family_index_bra)
      for(int irrep_family_index_ket=0; irrep_family_index_ket<=irrep_family_index_bra; ++irrep_family_index_ket)
        {
          spncci::AddLGIPairForOneBodyRecurrence(
            lgi_full_space_index_lookup,spncci_space,Nmax,
            irrep_family_index_bra,irrep_family_index_ket,sort_map
          );
        }

    for(auto it=sort_map.begin(); it!=sort_map.end(); ++it)
      {
        for(const auto& pair : it->second)
          lgi_pairs_ob.push_back(pair);
      }
  }

void GetLGIPairsForTwoBodyRecurrence(
      const std::vector<int>& lgi_full_space_index_lookup,
      const spncci::SpNCCISpace& spncci_space,
      int Nmax,
      std::vector<spncci::LGIPair>& lgi_pairs_tb
    )
  {
    // Organize lgi pairs by basis size
    std::map<int, std::vector<spncci::LGIPair>, std::greater<int> > sort_map;

    int num_irrep_families=lgi_full_space_index_lookup.size();
    for(int irrep_family_index_bra=0; irrep_family_index_bra<num_irrep_families; ++irrep_family_index_bra)
      for(int irrep_family_index_ket=0; irrep_family_index_ket<=irrep_family_index_bra; ++irrep_family_index_ket)
        {
          spncci::AddLGIPairForTwoBodyRecurrence(
            lgi_full_space_index_lookup,spncci_space,Nmax,
            irrep_family_index_bra,irrep_family_index_ket,sort_map
          );
        }

    for(auto it=sort_map.begin(); it!=sort_map.end(); ++it)
      {
        for(const auto& pair : it->second)
          lgi_pairs_tb.push_back(pair);
      }
  }

  void PopulateHypersectorsWithSeedsForOneBodyRecurrence(
    int irrep_family_index_bra, int irrep_family_index_ket,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
    const spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersectors_Nn0,
    const spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersectors,
    const std::vector<u3shell::OneBodyUnitTensorLabelsU3S>& lgi_unit_tensors,
    const std::vector<int>& rho0_values,
    basis::OperatorBlocks<double>& unit_tensor_seed_blocks,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_Nn0,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  )
  {
    int baby_spncci_index_bra
      =spncci::FromIrrepIndexGetBabySpncciIndex(irrep_family_index_bra,lgi_families,baby_spncci_space);
    int baby_spncci_index_ket
      =spncci::FromIrrepIndexGetBabySpncciIndex(irrep_family_index_ket,lgi_families,baby_spncci_space);

    const lgi::LGI& lgi_bra=lgi_families[irrep_family_index_bra].irrep;
    const lgi::LGI& lgi_ket=lgi_families[irrep_family_index_ket].irrep;
    const u3::U3& sigmap=lgi_bra.sigma();
    const u3::U3& sigma=lgi_ket.sigma();

    double conjugation_factor_base
        =ParitySign(u3::ConjugationGrade(sigmap)-lgi_bra.S()-u3::ConjugationGrade(sigma)+lgi_ket.S())
          * sqrt(1.*u3::dim(sigmap)*am::dim(lgi_bra.S())/u3::dim(sigma)/am::dim(lgi_ket.S()));
    // ParitySign(n) probably is (-1)^n
    // u3::ConjugationGrade(sigma) is lambda_sigma+mu_sigma

    // std::cout<<"loop over lgi unit tensors"<<std::endl;
    for(int i=0; i<lgi_unit_tensors.size();  ++i)
      {
        const u3shell::OneBodyUnitTensorLabelsU3S& unit_tensor=lgi_unit_tensors[i];
        int rho0=rho0_values[i];

        // Extract unit tensor labels
        u3::SU3 x0;
        HalfInt S0;
        int etap, eta, Tz;
        std::tie(x0,S0,etap,eta,Tz)=unit_tensor.FlatKey();

        // Look up unit tensor subspace
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
        int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
        auto& subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        // Look up index of unit tensor in subspace
        int unit_tensor_state_index
              =subspace.LookUpStateIndex(std::tuple<int>(Tz));

        // Get Hypersector index
        int hypersector_index
            =baby_spncci_hypersectors.LookUpHypersectorIndex(
                baby_spncci_index_bra,baby_spncci_index_ket,
                unit_tensor_subspace_index,rho0
              );

        unit_tensor_hyperblocks[hypersector_index][unit_tensor_state_index]=unit_tensor_seed_blocks[i];

        // Get conjugate

        // Look up conjugate unit tensor subspace
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
        int unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
        auto& subspace_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj);

        // Look up index of unit tensor in subspace
        int unit_tensor_state_index_conj
              =subspace_conj.LookUpStateIndex(std::tuple<int>(Tz));

        // Get Hypersector index
        int hypersector_index_Nn0
            =baby_spncci_hypersectors_Nn0.LookUpHypersectorIndex(
                baby_spncci_index_ket,baby_spncci_index_bra,
                unit_tensor_subspace_index_conj,rho0
              );

        double conjugation_factor
        =ParitySign(etap+eta)*sqrt(double((eta+1)*(eta+2))/double((etap+1)*(etap+2)))*conjugation_factor_base;
        unit_tensor_hyperblocks_Nn0[hypersector_index_Nn0][unit_tensor_state_index_conj]
          =conjugation_factor*unit_tensor_seed_blocks[i].transpose();

      }
    }

  void PopulateHypersectorsWithSeedsForTwoBodyRecurrence(
    int irrep_family_index_bra, int irrep_family_index_ket,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::TwoBodyDensitySpace& tbd_space,
    const spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersectors_Nn0,
    const spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersectors,
    const std::vector<u3shell::TwoBodyDensityLabels>& lgi_tbds,
    const std::vector<int>& rho_values,
    basis::OperatorBlocks<double>& tbd_seed_blocks,
    basis::OperatorHyperblocks<double>& tbd_hyperblocks_Nn0,
    basis::OperatorHyperblocks<double>& tbd_hyperblocks,
    u3::PhiCoefCache& phi_coef_cache
  )
  {
    int baby_spncci_index_bra
      =spncci::FromIrrepIndexGetBabySpncciIndex(irrep_family_index_bra,lgi_families,baby_spncci_space);
    int baby_spncci_index_ket
      =spncci::FromIrrepIndexGetBabySpncciIndex(irrep_family_index_ket,lgi_families,baby_spncci_space);

    const lgi::LGI& lgi_bra=lgi_families[irrep_family_index_bra].irrep;
    const lgi::LGI& lgi_ket=lgi_families[irrep_family_index_ket].irrep;
    const u3::U3& sigmap=lgi_bra.sigma();
    const u3::U3& sigma=lgi_ket.sigma();

    double conjugation_factor_base
        =ParitySign(u3::ConjugationGrade(sigmap)-lgi_bra.S()-u3::ConjugationGrade(sigma)+lgi_ket.S())
          * sqrt(1.*u3::dim(sigmap)*am::dim(lgi_bra.S())/u3::dim(sigma)/am::dim(lgi_ket.S()));
    // ParitySign(n) is (-1)^n
    // u3::ConjugationGrade(sigma) is lambda_sigma+mu_sigma

    // std::cout<<"loop over lgi TBDs"<<std::endl;
    for(int i=0; i<lgi_tbds.size();  ++i)
      {
        const u3shell::TwoBodyDensityLabels& tbd=lgi_tbds[i];
        int rho=rho_values[i];

        // Extract TBD labels
        u3::SU3 x0,xf,xi;
        HalfInt S0;
        int N1,N2,N3,N4,Sf,Si,rho0,Tz;
        std::tie(x0,S0,N1,N2,N3,N4,xf,Sf,xi,Si,rho0,Tz)=tbd.FlatKey();

        // Look up TBD subspace
        u3shell::TwoBodyDensitySubspaceLabels tbd_subspace_labels(x0,S0,N1,N2,N3,N4);
        int tbd_subspace_index=tbd_space.LookUpSubspaceIndex(tbd_subspace_labels);
        auto& subspace=tbd_space.GetSubspace(tbd_subspace_index);

        // Look up index of TBD in subspace
        int tbd_state_index
              =subspace.LookUpStateIndex(std::tuple<u3::SU3,int,u3::SU3,int,int,int>(xf,Sf,xi,Si,rho0,Tz));

        // Get Hypersector index
        int hypersector_index
            =baby_spncci_hypersectors.LookUpHypersectorIndex(
                baby_spncci_index_bra,baby_spncci_index_ket,
                tbd_subspace_index,rho
              );

        tbd_hyperblocks[hypersector_index][tbd_state_index]=tbd_seed_blocks[i];

        // Get conjugate

        // Look up conjugate TBD subspace
        u3shell::TwoBodyDensitySubspaceLabels tbd_subspace_labels_conj(u3::Conjugate(x0),S0,N4,N3,N2,N1);
        int tbd_subspace_index_conj=tbd_space.LookUpSubspaceIndex(tbd_subspace_labels_conj);
        auto& subspace_conj=tbd_space.GetSubspace(tbd_subspace_index_conj);

        // Get Hypersector index
        int hypersector_index_Nn0
            =baby_spncci_hypersectors_Nn0.LookUpHypersectorIndex(
                baby_spncci_index_ket,baby_spncci_index_bra,
                tbd_subspace_index_conj,rho
              );

	int rho0max=u3::OuterMultiplicity(xf,xi,x0);
	double conjugation_factor=ParitySign(u3::ConjugationGrade(xf)+u3::ConjugationGrade(xi)-u3::ConjugationGrade(x0)+N1+N2+N3+N4
                          +rho0max-rho0)*conjugation_factor_base;
        for(int rho0p=1; rho0p<=rho0max; rho0p++){

          // Look up index of TBD in subspace
          int tbd_state_index_conj
            =subspace_conj.LookUpStateIndex(std::tuple<u3::SU3,int,u3::SU3,int,int,int>(u3::Conjugate(xi),Si,u3::Conjugate(xf),Sf,rho0p,Tz));

          tbd_hyperblocks_Nn0[hypersector_index_Nn0][tbd_state_index_conj]+=conjugation_factor
	    *u3::PhiCached(phi_coef_cache,u3::Conjugate(xi),u3::Conjugate(xf),u3::Conjugate(x0),rho0p,rho0)*tbd_seed_blocks[i].transpose();

        }

      }

    }
//************************************************************************************

void AddNn0BlocksToHyperblocks(
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors_Nn0,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_Nn0,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
)
{
  for(int hypersector_index_Nn0=0; hypersector_index_Nn0<baby_spncci_hypersectors_Nn0.size(); ++hypersector_index_Nn0)
    {
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Extracting labels from source (Nn0 sectors)
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Get hypersector indices for Nn0 sector
      // Note the labels from Nn0 sectors are conjugate labels, i.e., bra is actually ket and vice versa
      auto key=baby_spncci_hypersectors_Nn0.GetHypersector(hypersector_index_Nn0).Key();
      int unit_tensor_subspace_index_Nn0, baby_spncci_index_bra, baby_spncci_index_ket, rho0;
      std::tie(baby_spncci_index_ket,baby_spncci_index_bra,unit_tensor_subspace_index_Nn0,rho0)=key;

      // Extract unit tensor subspace labels from Nn0 tensor
      auto& unit_tensor_subspace_Nn0=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_Nn0);
      u3::SU3 x0c;
      HalfInt S0;
      int etap,eta;
      std::tie(x0c,S0,eta,etap)=unit_tensor_subspace_Nn0.labels();

      // Get bra and ket labels from Nn0 sector
      const spncci::BabySpNCCISubspace& subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra);
      const spncci::BabySpNCCISubspace& subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket);

      u3::U3 omegap=subspace_bra.omega();
      HalfInt Sp=subspace_bra.S();

      u3::U3 omega=subspace_ket.omega();
      HalfInt S=subspace_ket.S();

      // part of conjugation factor
      double conjugation_factor_base
              =ParitySign(u3::ConjugationGrade(omega)+S-u3::ConjugationGrade(omegap)-Sp)
                *sqrt(1.*u3::dim(omega)*u3::dim(u3::SU3(etap,0))*am::dim(S)
                       /u3::dim(omegap)/u3::dim(u3::SU3(eta,0))/am::dim(Sp)
                  );

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Looking up target hypersector
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Get unit tensor subspace index in hyperblocks
      u3::SU3 x0(u3::Conjugate(x0c));
      u3shell::UnitTensorSubspaceLabels unit_tensor_labels
        =u3shell::UnitTensorSubspaceLabels(x0,S0,etap,eta);

      int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
      auto& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

      // Look up hypersector
      int hypersector_index
          =baby_spncci_hypersectors.LookUpHypersectorIndex(
              baby_spncci_index_bra,baby_spncci_index_ket,unit_tensor_subspace_index,rho0
            );

      // std::cout<<"hypersector "<<hypersector_index<<"  "<<subspace_bra.LabelStr()<<" "<<subspace_ket.LabelStr()<<std::endl;
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // For each source hyperblock, identify target block and conjugate
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for(int unit_tensor_index_Nn0=0; unit_tensor_index_Nn0<unit_tensor_subspace_Nn0.size(); ++unit_tensor_index_Nn0)
        {
          int Tbp,Sbp,Sb,Tb,T0;
          std::tie(T0,Sb,Tb,Sbp,Tbp)=unit_tensor_subspace_Nn0.GetStateLabels(unit_tensor_index_Nn0);

          // Get unit tensor index
          std::tuple<int,int,int,int,int> state_labels(T0,Sbp,Tbp,Sb,Tb);
          int unit_tensor_index=unit_tensor_subspace.LookUpStateIndex(state_labels);

          // Conjugation factor

          double conjugation_factor
                  =conjugation_factor_base*sqrt(am::dim(Sbp)*am::dim(Tbp)/am::dim(Sb)/am::dim(Tb));

          // std::cout<<fmt::format("{}  {}    {} {} {}  {} {} {}  {} {} {}",
          //   hypersector_index,unit_tensor_index,x0.Str(),S0,T0,etap,Sbp,Tbp,eta,Sb,Tb)<<std::endl;
          // std::cout<<conjugation_factor*unit_tensor_hyperblocks_Nn0[hypersector_index_Nn0][unit_tensor_index_Nn0]<<std::endl;
          unit_tensor_hyperblocks[hypersector_index][unit_tensor_index]
            =conjugation_factor*unit_tensor_hyperblocks_Nn0[hypersector_index_Nn0][unit_tensor_index_Nn0].transpose();
        }
    }
}



  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Construct KBUK matrix
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
  // Was modified to account for new Sp(3,R) indexing and VCS K matrix calculation
  // Has NOT been validated
  // Computes chi[upsilon1 omega1;upsilon omega;(20)]=KBUK(upsilon1,upsilon) for given omega1,omega,sigma
  // chi is given by Eq.(E.25) in Anna's thesis
  void ConsructKBUK(
    u3::UCoefCache& u_coef_cache, int Nn, // Nn=N_omega-N_sigma
    const u3::U3& sigma, const u3::U3& omega, const u3::U3& omega1,
    const sp3r::U3Subspace& u3_subspace,
    const sp3r::U3Subspace& u3_subspace1,
    const Eigen::MatrixXd& K1, //(dim1,upsilon1_max), K1(i,upsilon1)=(K_sigma^omega1)_[n1,rho1]upsilon1, where i is index of [n1,rho1]
    const Eigen::MatrixXd& K_inv,//(upsilon_max,dim), K_inv(upsilon,i)=(K_sigma^omega1)^-1_upsilon[n,rho], where i is index of [n,rho]
    Eigen::MatrixXd& KBUK //(upislon1_max,upsilon)
    )
  {
    int upsilon_max1=u3_subspace1.upsilon_max();
    int upsilon_max=u3_subspace.upsilon_max();

    const auto& nonorthogonal_basis1 = u3_subspace1.nonorthogonal_basis();
    const auto& nonorthogonal_basis = u3_subspace.nonorthogonal_basis();
    int dim1=nonorthogonal_basis1.dimension(); // dim1 is the number of [n1,rho1]'s
    int dim=nonorthogonal_basis.dimension(); // dim is the number of [n,rho]'s

    basis::OperatorBlock<double> BU = basis::OperatorBlock<double>::Zero(dim1, dim); // BU is zero matrix of size dim1xdim
    // Now matrix BU(i,j)=2/Nn*(n||a^dagger||n1)*U[(20) n1 omega sigma;n _ rho;omega1 rho1 _], 
    // where i and j are indices for [n1,rho1] and [n,rho], is computed
    for(int u3_state_index=0; u3_state_index<nonorthogonal_basis.size(); ++u3_state_index) // loop over n
      // nonorthogonal_basis.size() is the number of n's
      {
        const auto& n = nonorthogonal_basis.GetState(u3_state_index).n(); // n is n
        const int rho_max = nonorthogonal_basis.GetState(u3_state_index).rho_max();

        // iterate over (n1,rho1)
        for (int u3_state_index1=0; u3_state_index1<nonorthogonal_basis1.size(); u3_state_index1++) // loop over n1
	  // nonorthogonal_basis1.size() is the number of n1's
          {
            const auto& n1 = nonorthogonal_basis1.GetState(u3_state_index1).n(); // n1 is n1
            const int rho1_max = nonorthogonal_basis1.GetState(u3_state_index1).rho_max();

            // u3::U3 n1(n1_rho1.irrep);
            if (u3::OuterMultiplicity(n1.SU3(), u3::SU3(2,0),n.SU3())==0) // n1 must couple with (20) to n
              continue;

            for(int rho=1; rho<=rho_max; ++rho) // loop over rho
              for(int rho1=1; rho1<=rho1_max; ++rho1) // loop over rho1
                {
                  int row=nonorthogonal_basis1.GetStateOffset(u3_state_index1,rho1);
                  int col=nonorthogonal_basis.GetStateOffset(u3_state_index,rho);
                  // BU(u3_state_index1,u3_state_index)
                  //   =2./Nn*u3boson::BosonCreationRME(n,n1)
                  //    *u3::UCached(u_coef_cache,u3::SU3(2,0),n1.SU3(),omega.SU3(),sigma.SU3(),
                  //       n.SU3(),1,n_rho.tag,omega1.SU3(),n1_rho1.tag,1
                  //     );

                  BU(row,col)
                    =2./Nn*u3boson::BosonCreationRME(n,n1)
                     *u3::UCached(u_coef_cache,u3::SU3(2,0),n1.SU3(),omega.SU3(),sigma.SU3(),
                        n.SU3(),1,rho,omega1.SU3(),rho1,1
                      );
                }
          }
      }

    // Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
//    KBUK.noalias()=K1.transpose()*BU*K_inv.transpose();
    KBUK.noalias()=K1*BU*K_inv;
  }


void 
ComputeUnitTensorHyperblocks(
  int Nmax, int N1v,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::KMatrixCache& k_matrix_map,
  const spncci::KMatrixCache& kinv_matrix_map,
  const spncci::SpNCCISpace& spncci_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
  const std::vector<std::vector<int>>& unit_tensor_hypersector_subsets,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks // Output
  // basis::OperatorHyperblocks<double> is std::vector<std::vector<OperatorBlock<double>>>
  // OperatorBlock<double> is Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
  )
// compute hyperblocks for unit tensors recursively
{
  // std::cout<<"in the recurrence"<<std::endl;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Set up for calculation
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int Nsum=2; Nsum<=2*Nmax; Nsum+=2) // Nsum=Nn+Nn', where Nn=N_omega-N_sigma and Nn'=N_omega'-N_sigma'
    {
      const std::vector<int>& unit_tensor_hypersectors=unit_tensor_hypersector_subsets[Nsum/2];

      // Parallelize here
      // unit_tensor_hyperblocks are zero initalized so there should be no race conditions in
      // writing each hyperbock to unit_tensor_hyperblocks
      // std::cout<<"omp_get_num_threads "<<omp_get_num_threads()<<std::endl;

      for(int i=0; i<unit_tensor_hypersectors.size(); ++i) // loop over hypersectors
      // for(int hypersector_index : unit_tensor_hypersectors)
      // Hypersector is given by the first index of unit_tensor_hyperblocks
      // Hypersector corresponds to sigma,Sp,Sn,S,omega,sigma',Sp',Sn',S',omega',omega0,S0,\bar{N},\bar{N}',rho0
        {
          int hypersector_index=unit_tensor_hypersectors[i];
          auto key=baby_spncci_hypersectors.GetHypersector(hypersector_index).Key();

          int unit_tensor_subspace_index, baby_spncci_subspace_indexp, baby_spncci_index, rho0;
          std::tie(baby_spncci_subspace_indexp,baby_spncci_index,unit_tensor_subspace_index,rho0)=key;

          const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
              =baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);

          const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
              =baby_spncci_space.GetSubspace(baby_spncci_index);

          const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace
              =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

          // Extract ket dimensions and mulitplicities
          int dim=baby_spncci_subspace_ket.dimension();
          int gamma_max=baby_spncci_subspace_ket.gamma_max();
          int upsilon_max=baby_spncci_subspace_ket.upsilon_max();

          // Extract bra dimensions and mulitplicities
          int dimp=baby_spncci_subspace_bra.dimension();
          int gamma_maxp=baby_spncci_subspace_bra.gamma_max();
          int upsilon_maxp=baby_spncci_subspace_bra.upsilon_max();

          // Extract Sp(3,R) space
          int irrep_family_index_bra=baby_spncci_subspace_bra.irrep_family_index();
          int irrep_family_index_ket=baby_spncci_subspace_ket.irrep_family_index();

          // std::cout<<"irrep_family_index_bra "<<irrep_family_index_bra<<" irrep_family_index_ket"<<irrep_family_index_ket<<std::endl;
          const sp3r::Sp3RSpace& irrep_bra = spncci_space[irrep_family_index_bra].Sp3RSpace();
          const sp3r::Sp3RSpace& irrep_ket = spncci_space[irrep_family_index_ket].Sp3RSpace();

          // extract subspace labels
          u3::U3 omegap,sigmap,omega,sigma; // p denotes prime. bra has primed quantum numbers
          u3::SU3 x0; // x0 is SU(3) irrep of omega0
          HalfInt S0,Sn_ket,Sp_ket,S_ket,Sn_bra,Sp_bra,S_bra;
          int etap,eta; // eta is \bar{N}, etap is \bar{N}'

          // Extracting labels
          // std::cout<<"extracting labels"<<std::endl;
          std::tie(sigmap,Sp_bra,Sn_bra,S_bra,omegap)=baby_spncci_subspace_bra.labels();
          std::tie(sigma,Sp_ket,Sn_ket,S_ket,omega)=baby_spncci_subspace_ket.labels();
          std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
          int Nn=baby_spncci_subspace_ket.Nn(); // Nn is N_omega-N_sigma

          // omega u3 subspace in irrep
          const sp3r::U3Subspace& u3_subspace=irrep_ket.LookUpSubspace(omega);
          const sp3r::U3Subspace& u3_subspacep=irrep_bra.LookUpSubspace(omegap);
          // Extracting K matrices for sp_irrep and sp_irrepp from the K_matrix_maps
          // std::cout<<"bunny1"<<std::endl;
          const vcs::MatrixCache& K_matrix_map_bra=k_matrix_map.at(sigmap); // We don't need this
          // std::cout<<"bunny2"<<std::endl;
          // Temporary fix for debug purposes
          vcs::MatrixCache null_cache; // We don't need this
          const vcs::MatrixCache& Kinv_matrix_map_bra=kinv_matrix_map.at(sigmap); // We don't need this
          // std::cout<<"bunny3"<<std::endl;
          const vcs::MatrixCache& K_matrix_map_ket=k_matrix_map.at(sigma);
          // std::cout<<"bunny4"<<std::endl;
          const vcs::MatrixCache& Kinv_matrix_map_ket=kinv_matrix_map.at(sigma);
          // std::cout<<"bunny5"<<std::endl;
          const Eigen::MatrixXd& Kp=K_matrix_map_bra.at(omegap); // We don't need this
          // std::cout<<"bunny6"<<std::endl;
          // std::cout<<sigma.Str()<<". "<<omega.Str()<<std::endl;
          const Eigen::MatrixXd& K_inv=Kinv_matrix_map_ket.at(omega);

          // Generate labels to sum over
          int rho0_max=u3::OuterMultiplicity(omega.SU3(),x0,omegap.SU3());

          // Precalculating kronecker products used in sum to calculate unit tensor matrix
          MultiplicityTagged<u3::U3>::vector omegapp_set=KroneckerProduct(omegap, u3::U3(0,0,-2));
	  // omegapp_set is vector of multiplicity tagged \bar{omega}
          MultiplicityTagged<u3::U3>::vector omega1_set=KroneckerProduct(omega, u3::U3(0,0,-2));
	  // omega1_set is vector of multiplicity tagged omega1
          MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(x0, u3::SU3(2,0));
	  // x0p_set is vector of multiplicity tagged SU(3) irreps of omega0' (N_omega0'=N_omega'-N_omega1)

          // std::cout<<"hypersector_index "<<hypersector_index<<std::endl;
           std::vector<basis::OperatorBlock<double>>& unit_tensor_blocks=unit_tensor_hyperblocks[hypersector_index];
	  // unit_tensor_blocks is vector of matrices (blocks, sectors) corresponding to given hypersector
	  // Matrix indices of block correspond to gamma,upsilon (whose number is dim) and gamma',upsilon' (whose number is dimp)
	  // Within hypersector blocks correspond to different operators (unit tensors)
	  // In my one-body case there will be at most 2 blocks in each hypersector: proton and/or neutron
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //  Calculate unit tensor matrix
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          int num_blocks=unit_tensor_blocks.size();

          for(auto& omega1_mult :omega1_set) // omega1_mult is multiplicity tagged omega1
            {
              u3::U3 omega1(omega1_mult.irrep); // omega1 is omega1

              if (not irrep_ket.ContainsSubspace(omega1)) // if omega1 is not present in Sp(3,R) irrep sigma
                continue;

              spncci::BabySpNCCISubspaceLabels baby_spncci_labels1(sigma,Sp_ket,Sn_ket,S_ket,omega1);
              int baby_spncci_subspace_index1=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labels1);
              // std::cout<<"bunny rabbit 1"<<std::endl;
              Eigen::MatrixXd K1=K_matrix_map_ket.at(omega1); // K1 is K_sigma^omega1
              const sp3r::U3Subspace& u3_subspace1=irrep_ket.LookUpSubspace(omega1);
              int upsilon_max1=u3_subspace1.upsilon_max(); // upsilon_max1 is maximal upsilon1

              int dim1=upsilon_max1*gamma_max; // gamma_max is maximal gamma (for ket, i.e., sigma)
              std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omega1;
              // std::cout<<"initializing"<<std::endl;
              // Initializing blocks for sum over omega1
              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omega1);
              // std::cout<<"here"<<std::endl;
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              // Construct KBUK matrix
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
              // fmt::print("Constructing for {}: {}  {}",sigma,omega,omega1);
              spncci::ConsructKBUK( // Computes chi[upsilon1 omega1;upsilon omega;(20)]=KBUK(upsilon1-1,upsilon-1) for given omega1,omega,sigma
                u_coef_cache, Nn,sigma, omega, omega1,
                u3_subspace,u3_subspace1,K1,K_inv,
                KBUK
              );

              // ////////////////////////////////////////////////////////////////////////////////////////////////////////
              //summing over x0'
              // std::cout<<"sum over x0"<<std::endl;
              for (auto& x0p_mult : x0p_set) // x0p_mult is multiplicity tagged SU(3) irrep of omega0'
                {
                  u3::SU3 x0p(x0p_mult.irrep); // x0p is SU(3) irrep of omega0'
                  // std::cout<<"x0p "<<x0p.Str()<<std::endl;
                  int rho0p_max=OuterMultiplicity(omega1.SU3(),x0p,omegap.SU3()); // rho0p_max is maximal rho1'

                  // summing over rho0'
                  // std::cout<<"summing over rho0p.  rho0p_max "<<rho0p_max<<std::endl;
                  for (int rho0p=1; rho0p<=rho0p_max; rho0p++) // rho0p is rho1'
                    {
                      // Zero initialize blocks accumlating sum over x0p and rho0p
                      std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_x0p;
                      ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_x0p);

                      double coef=0;
                      for(int rho0b=1; rho0b<=rho0_max; rho0b++) // rho0b is rho0'
                        {
                          //(2,0)xx0->x0p(by construction),
                          //(2,0)xomega1->omega (by construction),
                          // x0xomega->omegap, (rho0_max)
                          //omega1xx0p->omegap (rho0p_max)
                          coef+=u3::PhiCached(phi_coef_cache,omega.SU3(),x0,omegap.SU3(),rho0,rho0b)
                               *u3::UCached(u_coef_cache,x0,u3::SU3(2,0),omegap.SU3(), omega1.SU3(),x0p,1,rho0p,omega.SU3(),1,rho0b);
                        }
		      // coef is sum_{rho0'}Phi_{rho0 rho0'}[omega omega0 omega']U[omega0 (20) omega' omega1; omega0' _ rho1'; omega _ rho0']

                      // std::cout<<"computing term 3"<<std::endl;
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // third term, i.e., first term in Eq.(E.26) in Anna's thesis
                      // sum over omega'', v'' and rho0''
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // Zero initialze blocks accumulating sum over omegapp
                      std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omegapp;
                      ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omegapp);

                      // Summing over omega''
                      // std::cout<<"summing over omegapp"<<std::endl;
                      for (auto& omegapp_mult : omegapp_set) // omegapp_mult is multiplicity tagged \bar{omega}
                        {
                          // A matrix will annihilate bra
                          if(baby_spncci_subspace_bra.Nn()==0) // if omega' is sigma', i.e., Nn'=N_omega'-N_sigma'=0
                            continue;

                          u3::U3 omegapp(omegapp_mult.irrep); // omegapp is \bar{omega}
                          // std::cout<<omegapp.Str()<<std::endl;

                          if (not irrep_bra.ContainsSubspace(omegapp)) // if \bar{omega} is not present in Sp(3,R) irrep sigma'
                            continue;

                          // get hypersector index
                          spncci::BabySpNCCISubspaceLabels baby_spncci_labelspp(sigmap,Sp_bra,Sn_bra,S_bra,omegapp);
                          int baby_spncci_subspace_indexpp=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labelspp);

                          // omega'' subspace (v'')
                          sp3r::U3Subspace u3_subspacepp=irrep_bra.LookUpSubspace(omegapp);
                          int upsilon_maxpp=u3_subspacepp.upsilon_max();
                          int dimpp=upsilon_maxpp*gamma_maxp;
                          // Obtaining K matrix for omega''

////////////////////////////////////////////////////////////////////////////////////////////////////////////                          
                          Eigen::MatrixXd A 
                            = sp3r::Sp3rRaisingOperator(sigmap,u3_subspacep,u3_subspacepp,u_coef_cache);

                          // spncci::Amatrix(u_coef_cache,
                          //   u3_subspacep,u3_subspacepp,sigmap, omegap, omegapp,
                          //   K_matrix_map_bra,Kinv_matrix_map_bra,upsilon_maxp, upsilon_maxpp,
                          //   A
                          // );
////////////////////////////////////////////////////////////////////////////////////////////////////////////

                           // Zero initialze blocks accumulating sum over rho0pp and rho0bp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rhobp;

                          ZeroInitBlocks(num_blocks,dimpp,dim1,unit_tensor_blocks_rhobp);

                          // Summing over rho0bp
                          int rho0bp_max=u3::OuterMultiplicity(omega1.SU3(),x0,omegapp.SU3());
                          for(int rho0bp=1; rho0bp<=rho0bp_max; ++rho0bp) // rho0bp is \bar{rho}
                            {
                              // Get hypersector index
                              int hypersector_index3
                                =baby_spncci_hypersectors.LookUpHypersectorIndex(baby_spncci_subspace_indexpp,baby_spncci_subspace_index1,unit_tensor_subspace_index,rho0bp);
                              // std::cout<<"hypersector3 "<<hypersector_index3<<std::endl;
                              if(hypersector_index3==-1)
                                continue;

                              double coef3=0;
                              for (int rho0pp=1; rho0pp<=rho0p_max; rho0pp++) // rho0pp is \bar{rho0}
                                {
                                // omegaxx0->omegapp
                                // x0x(2,0)->x0p (construction)
                                // omegappx(2,0)->omegap (construction)
                                // omegaxx0p->omegap
                                  // std::cout<<"heres a phi"<<std::endl;
                                  coef3+=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0pp)
                                          *u3::UCached(u_coef_cache,
                                            omega1.SU3(),x0,omegap.SU3(),u3::SU3(2,0),
                                            omegapp.SU3(),rho0bp,1,x0p,1,rho0pp
                                            );
                                }// end rho0pp
				// coef3 is sum_{\bar{rho0}}Phi_{rho1' \bar{rho0}}[omega0' omega1 omega']U[omega1 omega0 omega' (20); \bar{omega} \bar{rho} _; omega0' _ \bar{rho0}]
                              // std::cout<<"sum over blocks for term 3"<<std::endl;

                              for(int b=0; b<num_blocks; ++b)
                                unit_tensor_blocks_rhobp[b]+=coef3*unit_tensor_hyperblocks[hypersector_index3][b];
                            }

                          // matrix product A*unit_tensor_block (v',v'')*(v'',v1)
                          // std::cout<<"add in A operator and sum over blocks "<<std::endl;
                          for(int b=0; b<num_blocks; ++b)
                            for(int i=0; i<gamma_maxp; ++i)
                              for(int j=0; j<gamma_max; ++j)
                                {
                                  // Get target indices
                                  int it=i*upsilon_maxp;
                                  int jt=j*upsilon_max1;
                                  // Get source indices
                                  int is=i*upsilon_maxpp;
                                  int js=j*upsilon_max1;

                                  unit_tensor_blocks_omegapp[b].block(it,jt,upsilon_maxp,upsilon_max1)
                                    +=A*unit_tensor_blocks_rhobp[b].block(is,js,upsilon_maxpp,upsilon_max1);
                                }
                        } //omegapp
                      // std::cout<<"accumulating over unit tensor blocks x0p"<<std::endl;
                      // accumulating sum over omegapp
                      for(int b=0; b<num_blocks; ++b)
                        unit_tensor_blocks_x0p[b]+=unit_tensor_blocks_omegapp[b];

                      // std::cout<<"term 1"<<std::endl;
                      // std::cout<<unit_tensor_blocks_x0p[0]<<std::endl;
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                      //first term, i.e., second term in Eq.(E.26) in Anna's thesis
                      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                      if(u3::OuterMultiplicity(u3::SU3(etap,0),u3::SU3(0,eta-2),x0p)>0)
                        {
                          assert((eta-2)>=0);
                          // look up index of subspace in unit tensor space
                          u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap,eta-2);

                          int unit_tensor_subspace_index1=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                          assert(unit_tensor_subspace_index1!=-1);

                          //(rbp,0)x(0,rb)->x0, (0,rb)x(2,0)->(0,rb-2), x0x(2,0)->x0p, (0,rb-2)x(rbp,0)->x0p
                          double
                          coef1=u3::UCached(u_coef_cache,u3::SU3(etap,0),u3::SU3(0,eta),x0p, u3::SU3(2,0),x0,1,1,u3::SU3(0,eta-2),1,1)
                                *sqrt(1.*u3::dim(x0p)*u3::dim(u3::SU3(eta,0))/(u3::dim(x0)));
			  // coef1 is U[(\bar{N}'0)(0\bar{N})omega0'(20);omega0;(0,\bar{N}-2)][dim(omega0')dim(\bar{N}0)/dim(omega0)]^1/2

                          // zero initialize blocks for accumulating first term in sum over rhobp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                            ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                          // summing over rho0bp and accumulating sectors in unit1_matrix.
                          for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp) // rho0bp is \bar{rho0}
                            {
                              int hypersector_index1=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                    baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                    unit_tensor_subspace_index1, rho0bp
                                  );

                              // Accumulate
                              if(hypersector_index1==-1)
                                continue;

                              // std::cout<<"num blocks "<<num_blocks<<std::endl;
                              for(int b=0; b<num_blocks; ++b)
                              {
                                unit_tensor_blocks_rho0bp[b]
                                  +=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                    *unit_tensor_hyperblocks[hypersector_index1][b];
                              }
                            } //end rho0bp

                          // std::cout<<"coef1 "<<coef1<<std::endl;
                          // std::cout<<unit_tensor_blocks_rho0bp[0]<<std::endl;
                          for(int b=0; b<num_blocks; ++b)
                            unit_tensor_blocks_x0p[b]+=coef1*unit_tensor_blocks_rho0bp[b];
                        }

                        // std::cout<<"term 2 "<<std::endl;
                        // std::cout<<unit_tensor_blocks_x0p[0]<<std::endl;
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // second term, i.e., third term in Eq.(E.26) in Anna's thesis
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // std::cout<<"term 2 "<<std::endl;
                        // std::cout<<"x0p "<<x0p.Str()<<std::endl;
                        if ((u3::OuterMultiplicity(u3::SU3(etap+2,0),u3::SU3(0,eta),x0p)>0) && (etap+2)<=Nmax+2*N1v)
                          {
                            // (2,0)x(rbp,0)->(rbp+2,0), (rbp,0)x(0,rb)->x0, (rbp+2,0)x(0,rb)->x0p, x0x(2,0)->x0p
                            // look up index of subspace in unit tensor space
                            u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap+2,eta);


                            int unit_tensor_subspace_index2=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                            assert(unit_tensor_subspace_index2!=-1);


                            double
                            coef2=-1*ParitySign(u3::ConjugationGrade(x0)-u3::ConjugationGrade(x0p))
                                    * u3::dim(u3::SU3(etap,0))*sqrt(u3::dim(x0p)/(6.*u3::dim(u3::SU3(eta,0))))
                                    *u3::UCached(u_coef_cache,u3::SU3(etap+2,0),u3::SU3(0,etap),x0p,x0,
                                            u3::SU3(2,0),1,1,u3::SU3(0,eta),1,1);

                            // zero initialize blocks for accumulating first term in sum over rhobp
                            std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                            for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp) // rho0bp is \bar{rho0}
                              {

                                int hypersector_index2=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                      baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                      unit_tensor_subspace_index2, rho0bp
                                    );

                                auto& bra_sub=baby_spncci_space.GetSubspace(baby_spncci_subspace_index1);
                                auto& ket_sub=baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);
                                auto& unit_sub=unit_tensor_space.GetSubspace(unit_tensor_subspace_index2);

                                // Accumulate
                                if(hypersector_index2==-1)
                                continue;

                                for(int b=0; b<num_blocks; ++b)
                                {
                                  unit_tensor_blocks_rho0bp[b]
                                    +=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                      *unit_tensor_hyperblocks[hypersector_index2][b];
                                }
                              } //end rho0bp

                            for(int b=0; b<num_blocks; ++b)
                              unit_tensor_blocks_x0p[b]+=coef2*unit_tensor_blocks_rho0bp[b];
                          }

                        for(int b=0; b<num_blocks; ++b)
                          unit_tensor_blocks_omega1[b]+=coef*unit_tensor_blocks_x0p[b];
                      }//end rho0p
                  }// end x0p sum
                // summing over n, rho, n1, rho1, v1
                for(int b=0; b<num_blocks; ++b)
                  for(int i=0; i<gamma_maxp; ++i)
                    for(int j=0; j<gamma_max; ++j)
                    {
                      int it=i*upsilon_maxp;
                      int jt=j*upsilon_max;
                      int is=i*upsilon_maxp;
                      int js=j*upsilon_max1;
                      // (v'v1) (v1 v)
                      unit_tensor_blocks[b].block(it,jt,upsilon_maxp,upsilon_max)
                        +=unit_tensor_blocks_omega1[b].block(is,js,upsilon_maxp,upsilon_max1)*KBUK;

                    }
                // std::cout<<"done with blocks"<<std::endl;
              }// end omega1_mult
        }// end hypersector index
    }// end Nsum
  // std::cout<<"end recurrence"<<std::endl;
  }



void DoRecurrenceInitialization(
  int Nmax, int N1v,
  const spncci::LGIPair& lgi_pair,
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  spncci::OperatorBlocks& lgi_transformations,
  bool transform_lgi_families,
  std::map<spncci::NnPair,std::set<int>>& unit_tensor_subspace_subsets,
  spncci::BabySpNCCIHypersectors& baby_spncci_hypersector_seeds,
  spncci::BabySpNCCIHypersectors& baby_spncci_hypersector_seeds_conj,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds_conj
  )
{
    //////////////////////////////////////////////////////////////////////
    // Extract lgi index  labels
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
    ///////////////////////////////////////////////////////////////////////////
    // Read in list of unit tensors between lgi pair and conjugates from files
    // Returned bool, files_found_test has no current use, but could be used
    //  to identify lgi pair with no non-zero rmes between them.
    // Corresponding rho0 values stored separately for later hypersector lookup
    ///////////////////////////////////////////////////////////////////////////
    // Initialize containers
    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensors;
    std::vector<int> rho0_values;

    // Get index corresponding to lgi in the full space.
    // Index may differ from lgi index in basis if space has been truncated
    int index1=lgi_full_space_index_lookup[irrep_family_index_bra];
    int index2=lgi_full_space_index_lookup[irrep_family_index_ket];

    // Read in operators
    std::string lgi_unit_tensor_filename
      =fmt::format("seeds/operators_{:06d}_{:06d}.dat",index1,index2);
    bool files_found_test=lgi::ReadUnitTensorLabels(lgi_unit_tensor_filename,lgi_unit_tensors,rho0_values);

    ///////////////////////////////////////////////////////////////////////////
    // Set up hypersectors and hyperblocks for seeds
    // Generate hypersectors from list of unit tensors and outer-mulitplicites
    //  read from files
    // Populate hyperblocks using seeds read in from file
    //////////////////////////////////////////////////////////////////////////

    // Reads in unit tensor seed blocks and stores them in a vector of blocks. Order
    // corresponds to order of (unit_tensor,rho0) pairs in corresponding operator file.
    basis::OperatorBlocks<double> unit_tensor_seed_blocks;
    std::string seed_filename
      =fmt::format("seeds/seeds_{:06d}_{:06d}.rmes",index1,index2);
    files_found_test&=lgi::ReadBlocks(seed_filename, lgi_unit_tensors.size(), unit_tensor_seed_blocks);

    // If transform_lgi_families=True, apply basis transformation to lgi
    if(transform_lgi_families)
      spncci::TransformSeeds(index1,index2,lgi_transformations,unit_tensor_seed_blocks);

    // Identify unit tensor subspaces for recurrence
    spncci::GenerateRecurrenceUnitTensors(
      Nmax,N1v,lgi_unit_tensors,
      unit_tensor_space,unit_tensor_subspace_subsets
    );

    //Generate hypersectors for unit tensor between lgi pair (lgi1,lgi2)
    baby_spncci_hypersector_seeds
      =spncci::BabySpNCCIHypersectors(
        lgi_families,baby_spncci_space,unit_tensor_space,
        unit_tensor_subspace_subsets[spncci::NnPair(0,0)],
        irrep_family_index_bra,irrep_family_index_ket
      );

    // Generate hypersectors for unit tensors between conjugate lgi pair (lgi2,lgi1)
    baby_spncci_hypersector_seeds_conj
      =spncci::BabySpNCCIHypersectors(
        lgi_families,baby_spncci_space,unit_tensor_space,
        unit_tensor_subspace_subsets[spncci::NnPair(0,0)],
        irrep_family_index_ket,irrep_family_index_bra
      );

    // Zero initialize seed hyperblocks and conjugate hyperblocks
    basis::SetHyperoperatorToZero(baby_spncci_hypersector_seeds,unit_tensor_hyperblocks_seeds);
    basis::SetHyperoperatorToZero(baby_spncci_hypersector_seeds_conj,unit_tensor_hyperblocks_seeds_conj);

    // Populate the hyperblocks and conjugate hyperblocks with the seed
    // Conjugate hyperspectors will be used in calculating Nn0 rmes in recurrence
    spncci::PopulateHypersectorsWithSeeds(
      irrep_family_index_bra, irrep_family_index_ket,lgi_families,
      baby_spncci_space,unit_tensor_space,
      baby_spncci_hypersector_seeds_conj,baby_spncci_hypersector_seeds,
      lgi_unit_tensors,rho0_values,unit_tensor_seed_blocks,
      unit_tensor_hyperblocks_seeds_conj,unit_tensor_hyperblocks_seeds
    );
    //////////////////////////////////////////////////////////////////////

}

  bool
    GenerateUnitTensorHyperblocks(
      const spncci::LGIPair& lgi_pair,
      int Nmax, int N1v,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const spncci::KMatrixCache& kinv_matrix_cache,
      std::map<spncci::NnPair,std::set<int>>& unit_tensor_subspace_subsets,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersector_seeds,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersector_seeds_conj,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds_conj,
      u3::UCoefCache& u_coef_cache,
      u3::PhiCoefCache& phi_coef_cache,
      spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
      )
  {
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;

    // std::cout<<"generate Nn0 hypersectors"<<std::endl;
    // Generate Nn=0 hypersectors to be computed by conjugation
    bool Nn0_conjugate_hypersectors=true;
    std::vector<std::vector<int>> unit_tensor_hypersector_subsets_Nn0;


    spncci::BabySpNCCIHypersectors baby_spncci_hypersectors_Nn0( // Construct Nnp=0 and Nn!=0 hypersectors
      Nmax, baby_spncci_space, unit_tensor_space,
      unit_tensor_subspace_subsets, unit_tensor_hypersector_subsets_Nn0,
      irrep_family_index_ket, irrep_family_index_bra,
      Nn0_conjugate_hypersectors
    );

    // Hypersectors for Nnp>=Nn
    Nn0_conjugate_hypersectors=false;
    std::vector<std::vector<int>> unit_tensor_hypersector_subsets;

    // (Nnp,Nn) sectors for Nnp>Nn
    // std::cout<<"generate Nnp>Nn hypersectors"<<std::endl;
    baby_spncci_hypersectors=spncci::BabySpNCCIHypersectors( // Construct Nnp>=Nn hypersectors
      Nmax,baby_spncci_space, unit_tensor_space,
      unit_tensor_subspace_subsets, unit_tensor_hypersector_subsets,
      irrep_family_index_bra,irrep_family_index_ket,
      Nn0_conjugate_hypersectors
    );

    //Zero initialize hyperblocks for both Nn=0 and Nnp>=Nn sectors
    basis::OperatorHyperblocks<double> unit_tensor_hyperblocks_Nn0;
    basis::SetHyperoperatorToZero(baby_spncci_hypersectors_Nn0,unit_tensor_hyperblocks_Nn0);
    basis::SetHyperoperatorToZero(baby_spncci_hypersectors,unit_tensor_hyperblocks);

    //Iterate over seed hypersectors and identify corresponding hypersector in hyperblocks or hyperblocks_Nn0
    //Transfer seeds into hyperblocks for recurrence
    for(int seed_hypersector_index=0; seed_hypersector_index<baby_spncci_hypersector_seeds.size(); ++seed_hypersector_index)
      {
        int bra_index,ket_index,operator_index,multiplicity_index;
        std::tie(bra_index,ket_index,operator_index,multiplicity_index)
          =baby_spncci_hypersector_seeds.GetHypersector(seed_hypersector_index).Key();
        int hypersector_index
          =baby_spncci_hypersectors.LookUpHypersectorIndex(bra_index,ket_index,operator_index,multiplicity_index);

        unit_tensor_hyperblocks[hypersector_index]=unit_tensor_hyperblocks_seeds[seed_hypersector_index];
      }

    //Transfer conjugate hypersectors into Nn0 hypersectors for recurrence
    for(int seed_hypersector_index=0; seed_hypersector_index<baby_spncci_hypersector_seeds_conj.size(); ++seed_hypersector_index)
      {
        int bra_index,ket_index,operator_index,multiplicity_index;
        std::tie(bra_index,ket_index,operator_index,multiplicity_index)
          =baby_spncci_hypersector_seeds_conj.GetHypersector(seed_hypersector_index).Key();
        int hypersector_index
          =baby_spncci_hypersectors_Nn0.LookUpHypersectorIndex(bra_index,ket_index,operator_index,multiplicity_index);

        unit_tensor_hyperblocks_Nn0[hypersector_index]=unit_tensor_hyperblocks_seeds_conj[seed_hypersector_index];
      }

    // std::cout<<"Compute Nn=0 blocks"<<std::endl;
    spncci::ComputeUnitTensorHyperblocks(
      Nmax,N1v,u_coef_cache,phi_coef_cache,
      k_matrix_cache,kinv_matrix_cache,spncci_space,baby_spncci_space,
      unit_tensor_space,baby_spncci_hypersectors_Nn0,
      unit_tensor_hypersector_subsets_Nn0,unit_tensor_hyperblocks_Nn0
    );

    // std::cout<<"Add Nn0 blocks to hyperblocks"<<std::endl;
    spncci::AddNn0BlocksToHyperblocks(
      baby_spncci_space,unit_tensor_space,
      baby_spncci_hypersectors_Nn0,baby_spncci_hypersectors,
      unit_tensor_hyperblocks_Nn0,unit_tensor_hyperblocks
    );

    // std::cout<<"Compute unit tensor hyperblocks"<<std::endl;
    spncci::ComputeUnitTensorHyperblocks(
      Nmax,N1v,u_coef_cache,phi_coef_cache,
      k_matrix_cache,kinv_matrix_cache,spncci_space,baby_spncci_space,
      unit_tensor_space,baby_spncci_hypersectors,
      unit_tensor_hypersector_subsets,unit_tensor_hyperblocks
    );
    return true;
  }

//************************************************ Added by J.H. ***********************************************
void AddNn0BlocksToOneBodyUnitTensorHyperblocks(
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersectors_Nn0,
  const spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersectors,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_Nn0,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
)
{
  for(int hypersector_index_Nn0=0; hypersector_index_Nn0<baby_spncci_hypersectors_Nn0.size(); ++hypersector_index_Nn0)
    {
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Extracting labels from source (Nn0 sectors)
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Get hypersector indices for Nn0 sector
      // Note the labels from Nn0 sectors are conjugate labels, i.e., bra is actually ket and vice versa
      auto key=baby_spncci_hypersectors_Nn0.GetHypersector(hypersector_index_Nn0).Key();
      int unit_tensor_subspace_index_Nn0, baby_spncci_index_bra, baby_spncci_index_ket, rho0;
      std::tie(baby_spncci_index_ket,baby_spncci_index_bra,unit_tensor_subspace_index_Nn0,rho0)=key;

      // Extract unit tensor subspace labels from Nn0 tensor
      auto& unit_tensor_subspace_Nn0=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_Nn0);
      u3::SU3 x0c;
      HalfInt S0;
      int etap,eta;
      std::tie(x0c,S0,eta,etap)=unit_tensor_subspace_Nn0.labels();

      // Get bra and ket labels from Nn0 sector
      const spncci::BabySpNCCISubspace& subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra);
      const spncci::BabySpNCCISubspace& subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket);

      u3::U3 omegap=subspace_bra.omega();
      HalfInt Sp=subspace_bra.S();

      u3::U3 omega=subspace_ket.omega();
      HalfInt S=subspace_ket.S();

      // conjugation factor
      double conjugation_factor
              =ParitySign(u3::ConjugationGrade(omega)+S-u3::ConjugationGrade(omegap)-Sp-etap+eta)
                *sqrt(double(u3::dim(omega)*am::dim(S)*(etap+1)*(etap+2))
                       /double(u3::dim(omegap)*am::dim(Sp)*(eta+1)*(eta+2))
                  );

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Looking up target hypersector
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Get unit tensor subspace index in hyperblocks
      u3::SU3 x0(u3::Conjugate(x0c));
      u3shell::UnitTensorSubspaceLabels unit_tensor_labels
        =u3shell::UnitTensorSubspaceLabels(x0,S0,etap,eta);

      int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
      auto& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

      // Look up hypersector
      int hypersector_index
          =baby_spncci_hypersectors.LookUpHypersectorIndex(
              baby_spncci_index_bra,baby_spncci_index_ket,unit_tensor_subspace_index,rho0
            );

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // For each source hyperblock, identify target block and conjugate
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for(int unit_tensor_index_Nn0=0; unit_tensor_index_Nn0<unit_tensor_subspace_Nn0.size(); ++unit_tensor_index_Nn0)
        {
          int Tz;
          std::tie(Tz)=unit_tensor_subspace_Nn0.GetStateLabels(unit_tensor_index_Nn0);

          // Get unit tensor index
          std::tuple<int> state_labels(Tz);
          int unit_tensor_index=unit_tensor_subspace.LookUpStateIndex(state_labels);

          unit_tensor_hyperblocks[hypersector_index][unit_tensor_index]
            =conjugation_factor*unit_tensor_hyperblocks_Nn0[hypersector_index_Nn0][unit_tensor_index_Nn0].transpose();
        }
    }
}

  void ConsructKBUK_JH(
    u3::UCoefCache& u_coef_cache, int Nn, // Nn=N_omega-N_sigma
    const u3::U3& sigma, const u3::U3& omega, const u3::U3& omega1,
    const sp3r::U3Subspace& u3_subspace,
    const sp3r::U3Subspace& u3_subspace1,
    const Eigen::MatrixXd& K1, //(dim1,upsilon1_max), K1(i,upsilon1)=(K_sigma^omega1)_[n1,rho1]upsilon1, where i is index of [n1,rho1]
    const Eigen::MatrixXd& K_inv,//(upsilon_max,dim), K_inv(upsilon,i)=(K_sigma^omega1)^-1_upsilon[n,rho], where i is index of [n,rho]
    Eigen::MatrixXd& KBUK, //(upislon1_max,upsilon)
    const bool write
    )
  {
    int upsilon_max1=u3_subspace1.upsilon_max();
    int upsilon_max=u3_subspace.upsilon_max();

    const auto& nonorthogonal_basis1 = u3_subspace1.nonorthogonal_basis();
    const auto& nonorthogonal_basis = u3_subspace.nonorthogonal_basis();
    int dim1=nonorthogonal_basis1.dimension(); // dim1 is the number of [n1,rho1]'s
    int dim=nonorthogonal_basis.dimension(); // dim is the number of [n,rho]'s

    basis::OperatorBlock<double> BU = basis::OperatorBlock<double>::Zero(dim1, dim); // BU is zero matrix of size dim1xdim
    // Now matrix BU(i,j)=2/Nn*(n||a^dagger||n1)*U[(20) n1 omega sigma;n _ rho;omega1 rho1 _], 
    // where i and j are indices for [n1,rho1] and [n,rho], is computed
    for(int u3_state_index=0; u3_state_index<nonorthogonal_basis.size(); ++u3_state_index) // loop over n
      // nonorthogonal_basis.size() is the number of n's
      {
        const auto& n = nonorthogonal_basis.GetState(u3_state_index).n(); // n is n
        const int rho_max = nonorthogonal_basis.GetState(u3_state_index).rho_max();

        // iterate over (n1,rho1)
        for (int u3_state_index1=0; u3_state_index1<nonorthogonal_basis1.size(); u3_state_index1++) // loop over n1
	  // nonorthogonal_basis1.size() is the number of n1's
          {
            const auto& n1 = nonorthogonal_basis1.GetState(u3_state_index1).n(); // n1 is n1
            const int rho1_max = nonorthogonal_basis1.GetState(u3_state_index1).rho_max();

            // u3::U3 n1(n1_rho1.irrep);
            if (u3::OuterMultiplicity(n1.SU3(), u3::SU3(2,0),n.SU3())==0) // n1 must couple with (20) to n
              continue;

            for(int rho=1; rho<=rho_max; ++rho) // loop over rho
              for(int rho1=1; rho1<=rho1_max; ++rho1) // loop over rho1
                {
                  int row=nonorthogonal_basis1.GetStateOffset(u3_state_index1,rho1);
                  int col=nonorthogonal_basis.GetStateOffset(u3_state_index,rho);
                  // BU(u3_state_index1,u3_state_index)
                  //   =2./Nn*u3boson::BosonCreationRME(n,n1)
                  //    *u3::UCached(u_coef_cache,u3::SU3(2,0),n1.SU3(),omega.SU3(),sigma.SU3(),
                  //       n.SU3(),1,n_rho.tag,omega1.SU3(),n1_rho1.tag,1
                  //     );

                  BU(row,col)
                    =2./Nn*u3boson::BosonCreationRME(n,n1)
                     *u3::UCached(u_coef_cache,u3::SU3(2,0),n1.SU3(),omega.SU3(),sigma.SU3(),
                        n.SU3(),1,rho,omega1.SU3(),rho1,1
                      );
                }
          }
      }

    // Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
    KBUK.noalias()=K1*BU*K_inv;
    if(write){
      std::cout<<"KBUK=K1*BU*K_inv, where KBUK is"<<std::endl;
      std::cout<<KBUK<<std::endl;
      std::cout<<"K1 is"<<std::endl;
      std::cout<<K1<<std::endl;
      std::cout<<"BU is"<<std::endl;
      std::cout<<BU<<std::endl;
      std::cout<<"K_inv is"<<std::endl;
      std::cout<<K_inv<<std::endl;
    }
  }

void ComputeOneBodyUnitTensorHyperblocks(
  int Nmax, int N1vp, int N1vn, int nucleon_number,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::KMatrixCache& k_matrix_map,
  const spncci::KMatrixCache& kinv_matrix_map,
  const spncci::SpNCCISpace& spncci_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersectors,
  const std::vector<std::vector<int>>& unit_tensor_hypersector_subsets,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks // Output
  // basis::OperatorHyperblocks<double> is std::vector<std::vector<OperatorBlock<double>>>
  // OperatorBlock<double> is Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
  )
// compute hyperblocks for unit tensors recursively
{
  for(int Nsum=2; Nsum<=2*Nmax; Nsum+=2) // Nsum=Nn+Nn', where Nn=N_omega-N_sigma and Nn'=N_omega'-N_sigma'
    {
      const std::vector<int>& unit_tensor_hypersectors=unit_tensor_hypersector_subsets[Nsum/2];

      for(int i=0; i<unit_tensor_hypersectors.size(); ++i) // loop over hypersectors
      // Hypersector is given by the first index of unit_tensor_hyperblocks
      // Hypersector corresponds to sigma,Sp,Sn,S,omega,sigma',Sp',Sn',S',omega',omega0,S0,eta',eta,rho0
        {
          int hypersector_index=unit_tensor_hypersectors[i];
          auto key=baby_spncci_hypersectors.GetHypersector(hypersector_index).Key();

          int unit_tensor_subspace_index, baby_spncci_subspace_indexp, baby_spncci_index, rho0;
          std::tie(baby_spncci_subspace_indexp,baby_spncci_index,unit_tensor_subspace_index,rho0)=key;

          const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
              =baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);

          const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
              =baby_spncci_space.GetSubspace(baby_spncci_index);

          const u3shell::OneBodyUnitTensorSubspaceU3S& unit_tensor_subspace
              =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

          // Extract ket dimensions and mulitplicities
          int dim=baby_spncci_subspace_ket.dimension();
          int gamma_max=baby_spncci_subspace_ket.gamma_max();
          int upsilon_max=baby_spncci_subspace_ket.upsilon_max();

          // Extract bra dimensions and mulitplicities
          int dimp=baby_spncci_subspace_bra.dimension();
          int gamma_maxp=baby_spncci_subspace_bra.gamma_max();
          int upsilon_maxp=baby_spncci_subspace_bra.upsilon_max();

          // Extract Sp(3,R) space
          int irrep_family_index_bra=baby_spncci_subspace_bra.irrep_family_index();
          int irrep_family_index_ket=baby_spncci_subspace_ket.irrep_family_index();

          const sp3r::Sp3RSpace& irrep_bra = spncci_space[irrep_family_index_bra].Sp3RSpace();
          const sp3r::Sp3RSpace& irrep_ket = spncci_space[irrep_family_index_ket].Sp3RSpace();

          // extract subspace labels
          u3::U3 omegap,sigmap,omega,sigma; // p denotes prime. bra has primed quantum numbers
          u3::SU3 x0; // x0 is Gamma0
          HalfInt S0,Sn_ket,Sp_ket,S_ket,Sn_bra,Sp_bra,S_bra;
          int etap,eta;

          // Extracting labels
          std::tie(sigmap,Sp_bra,Sn_bra,S_bra,omegap)=baby_spncci_subspace_bra.labels();
          std::tie(sigma,Sp_ket,Sn_ket,S_ket,omega)=baby_spncci_subspace_ket.labels();
          std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
          int Nn=baby_spncci_subspace_ket.Nn(); // Nn is N_omega-N_sigma

          // omega u3 subspace in irrep
          const sp3r::U3Subspace& u3_subspace=irrep_ket.LookUpSubspace(omega);
          const sp3r::U3Subspace& u3_subspacep=irrep_bra.LookUpSubspace(omegap);
          const vcs::MatrixCache& K_matrix_map_ket=k_matrix_map.at(sigma);
          const vcs::MatrixCache& Kinv_matrix_map_ket=kinv_matrix_map.at(sigma);
          const Eigen::MatrixXd& K_inv=Kinv_matrix_map_ket.at(omega);

          // Generate labels to sum over
          int rho0_max=u3::OuterMultiplicity(omega.SU3(),x0,omegap.SU3());

          // Precalculating kronecker products used in sum to calculate unit tensor matrix
          MultiplicityTagged<u3::U3>::vector omegapp_set=KroneckerProduct(omegap, u3::U3(0,0,-2));
	  // omegapp_set is vector of multiplicity tagged \bar{omega}
          MultiplicityTagged<u3::U3>::vector omega1_set=KroneckerProduct(omega, u3::U3(0,0,-2));
	  // omega1_set is vector of multiplicity tagged omega1
          MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(x0, u3::SU3(2,0));
	  // x0p_set is vector of multiplicity tagged SU(3) irreps of omega0' (N_omega0'=N_omega'-N_omega1)

          std::vector<basis::OperatorBlock<double>>& unit_tensor_blocks=unit_tensor_hyperblocks[hypersector_index];
	  // unit_tensor_blocks is vector of matrices (blocks, sectors) corresponding to given hypersector
	  // Matrix indices of block correspond to gamma,upsilon (whose number is dim) and gamma',upsilon' (whose number is dimp)
	  // Within hypersector blocks correspond to different operators (unit tensors)
	  // For one-body unit tensors there are at most 2 blocks in each hypersector: proton and/or neutron
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //  Calculate unit tensor matrix
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          int num_blocks=unit_tensor_blocks.size();

int ii=-1000;
int Nsummy;
bool write=false;
if(sigmap.N().TwiceValue()==19 && sigmap.SU3().lambda()==2 && sigmap.SU3().mu()==0 && Sp_bra.TwiceValue()==1 && Sn_bra.TwiceValue()==1 && S_bra==0 && omegap.N().TwiceValue()==19 && omegap.SU3().lambda()==2 && omegap.SU3().mu()==0 && sigma.N().TwiceValue()==19 && sigma.SU3().lambda()==2 && sigma.SU3().mu()==0 && Sp_ket.TwiceValue()==1 && Sn_ket.TwiceValue()==1 && S_ket.TwiceValue()==0 && omega.N().TwiceValue()==31 && omega.SU3().lambda()==4 && omega.SU3().mu()==2 && x0.lambda()==0 && x0.mu()==6 && S0.TwiceValue()==0 && etap==0 && eta==6){
//std::cout<<"******************************************* OB recurrence begins *************************************"<<std::endl;	
//ii=i;
Nsummy=Nsum;
//write=true;
}

          for(auto& omega1_mult :omega1_set) // omega1_mult is multiplicity tagged omega1
            {
              u3::U3 omega1(omega1_mult.irrep); // omega1 is omega1

              if (not irrep_ket.ContainsSubspace(omega1)) // if omega1 is not present in Sp(3,R) irrep sigma
                continue;

if(Nsum==Nsummy && i==ii)std::cout<<"omega1: "<<omega1.N()<<" "<<omega1.SU3().lambda()<<" "<<omega1.SU3().mu()<<std::endl;

              spncci::BabySpNCCISubspaceLabels baby_spncci_labels1(sigma,Sp_ket,Sn_ket,S_ket,omega1);
              int baby_spncci_subspace_index1=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labels1);
              Eigen::MatrixXd K1=K_matrix_map_ket.at(omega1); // K1 is K_sigma^omega1
              const sp3r::U3Subspace& u3_subspace1=irrep_ket.LookUpSubspace(omega1);
              int upsilon_max1=u3_subspace1.upsilon_max(); // upsilon_max1 is maximal upsilon1

              int dim1=upsilon_max1*gamma_max; // gamma_max is maximal gamma (for ket, i.e., sigma)
              std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omega1;
              // Initializing blocks for sum over omega1
              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omega1);

              // Construct KBUK matrix
              Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
              spncci::ConsructKBUK_JH( // Computes chi[upsilon1 omega1;upsilon omega;(20)]=KBUK(upsilon1-1,upsilon-1) for given omega1,omega,sigma
                u_coef_cache, Nn,sigma, omega, omega1,
                u3_subspace,u3_subspace1,K1,K_inv,
                KBUK,write
              );

              double coef=sqrt(1.*u3::dim(omega.SU3())/(1.*u3::dim(omega1.SU3())))/double(u3::dim(x0));

              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              // first term
              // sum over \bar{omega}, \bar{upsilon} and \bar{rho}
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              // Zero initialze blocks accumulating sum over \bar{omega}
              std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omegapp;
              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omegapp);

              // Summing over \bar{omega}
              for (auto& omegapp_mult : omegapp_set) // omegapp_mult is multiplicity tagged \bar{omega}
                {
                   // A matrix will annihilate bra
                   if(baby_spncci_subspace_bra.Nn()==0) // if omega' is sigma', i.e., Nn'=N_omega'-N_sigma'=0
                     continue;
  
                   u3::U3 omegapp(omegapp_mult.irrep); // omegapp is \bar{omega}

                   if (not irrep_bra.ContainsSubspace(omegapp)) // if \bar{omega} is not present in Sp(3,R) irrep sigma'
                     continue;

if(Nsum==Nsummy && i==ii)std::cout<<"omegapp: "<<omegapp.N()<<" "<<omegapp.SU3().lambda()<<" "<<omegapp.SU3().mu()<<std::endl;

                   // get hypersector index
                   spncci::BabySpNCCISubspaceLabels baby_spncci_labelspp(sigmap,Sp_bra,Sn_bra,S_bra,omegapp);
                   int baby_spncci_subspace_indexpp=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labelspp);

                   // \bar{omega} subspace (\bar{upsilon})
                   sp3r::U3Subspace u3_subspacepp=irrep_bra.LookUpSubspace(omegapp);
                   int upsilon_maxpp=u3_subspacepp.upsilon_max(); // upsilon_maxpp is maximal \bar{upsilon}
                   int dimpp=upsilon_maxpp*gamma_maxp;

                   Eigen::MatrixXd A 
                     = sp3r::Sp3rRaisingOperator(sigmap,u3_subspacep,u3_subspacepp,u_coef_cache);

                   // Zero initialze blocks accumulating sum over rho2, rho3 and rho0bp
                   std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rhobp;

                   ZeroInitBlocks(num_blocks,dimpp,dim1,unit_tensor_blocks_rhobp);

                   // Summing over rho0bp
                   int rho0bp_max=u3::OuterMultiplicity(omega1.SU3(),x0,omegapp.SU3());
                   for(int rho0bp=1; rho0bp<=rho0bp_max; ++rho0bp) // rho0bp is \bar{rho}
                     {
                       // Get hypersector index
                       int hypersector_index3
                         =baby_spncci_hypersectors.LookUpHypersectorIndex(baby_spncci_subspace_indexpp,baby_spncci_subspace_index1,unit_tensor_subspace_index,rho0bp);
                       if(hypersector_index3==-1)
                         continue;

if(Nsum==Nsummy && i==ii)std::cout<<"bar{rho}: "<<rho0bp<<std::endl;

                       double coef3=0;
                       for (int rho2=1; rho2<=rho0bp_max; rho2++)
                         {
               	           for (int rho3=1; rho3<=rho0bp_max; rho3++)
			     {

if(Nsum==Nsummy && i==ii)std::cout<<"rho2 rho3: "<<rho2<<" "<<rho3<<std::endl;
if(Nsum==Nsummy && i==ii)std::cout<<"coef3=coef3+Phi*Phi*U="<<coef3<<"+"<<u3::PhiCached(phi_coef_cache,omegapp.SU3(),u3::Conjugate(x0),omega1.SU3(),rho0bp,rho2)<<"*"<<PhiCached(phi_coef_cache,omega1.SU3(),u3::Conjugate(omegapp.SU3()),u3::Conjugate(x0),rho2,rho3)<<"*"<<u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(omegapp.SU3()),omega.SU3(),omega1.SU3(),u3::SU3(2,0),1,1,u3::Conjugate(x0),rho3,rho0)<<"="<<coef3+u3::PhiCached(phi_coef_cache,omegapp.SU3(),u3::Conjugate(x0),omega1.SU3(),rho0bp,rho2)*PhiCached(phi_coef_cache,omega1.SU3(),u3::Conjugate(omegapp.SU3()),u3::Conjugate(x0),rho2,rho3)*u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(omegapp.SU3()),omega.SU3(),omega1.SU3(),u3::SU3(2,0),1,1,u3::Conjugate(x0),rho3,rho0)<<std::endl;

			        coef3+=u3::PhiCached(phi_coef_cache,omegapp.SU3(),u3::Conjugate(x0),omega1.SU3(),rho0bp,rho2)
			   	  *PhiCached(phi_coef_cache,omega1.SU3(),u3::Conjugate(omegapp.SU3()),u3::Conjugate(x0),rho2,rho3)
				  *u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(omegapp.SU3()),omega.SU3(),omega1.SU3(),
				  	       u3::SU3(2,0),1,1,u3::Conjugate(x0),rho3,rho0);
			     }
		         }

if(Nsum==Nsummy && i==ii)std::cout<<"coef3=factor*coef3="<<ParitySign(u3::ConjugationGrade(x0)+u3::ConjugationGrade(omega.SU3())+u3::ConjugationGrade(omegap.SU3()))*sqrt(1.*u3::dim(omega1.SU3())*u3::dim(x0)*u3::dim(omegapp.SU3())/6)<<"*"<<coef3<<"="<<ParitySign(u3::ConjugationGrade(x0)+u3::ConjugationGrade(omega.SU3())+u3::ConjugationGrade(omegap.SU3()))*sqrt(1.*u3::dim(omega1.SU3())*u3::dim(x0)*u3::dim(omegapp.SU3())/6)*coef3<<std::endl;

		       coef3=ParitySign(u3::ConjugationGrade(x0)+u3::ConjugationGrade(omega.SU3())
	                     +u3::ConjugationGrade(omegap.SU3()))
		   	     *sqrt(1.*u3::dim(omega1.SU3())*u3::dim(x0)*u3::dim(omegapp.SU3())/6)*coef3;

if(Nsum==Nsummy && i==ii)std::cout<<"unit_tensor_blocks_rhobp[0]=unit_tensor_blocks_rhobp[0]+coef3*unit_tensor_hyperblocks[hypersector_index3][0]="<<unit_tensor_blocks_rhobp[0]<<"+"<<coef3<<"*"<<unit_tensor_hyperblocks[hypersector_index3][0]<<"="<<unit_tensor_blocks_rhobp[0]+coef3*unit_tensor_hyperblocks[hypersector_index3][0]<<std::endl;

                       // sum over blocks for term 1
                       for(int b=0; b<num_blocks; ++b)
                         unit_tensor_blocks_rhobp[b]+=coef3*unit_tensor_hyperblocks[hypersector_index3][b];
                     }

if(Nsum==Nsummy && i==ii)std::cout<<"unit_tensor_blocks_omegapp[0]=unit_tensor_blocks_omegapp[0]+A*unit_tensor_blocks_rhobp[0]="<<unit_tensor_blocks_omegapp[0]<<"+"<<A<<"*"<<unit_tensor_blocks_rhobp[0]<<"="<<unit_tensor_blocks_omegapp[0]+A*unit_tensor_blocks_rhobp[0]<<std::endl;

                   // matrix product A*unit_tensor_block (upsilon',\bar{upsilon})*(\bar{upsilon},upsilon1)
                   // add in A operator and sum over blocks
                   for(int b=0; b<num_blocks; ++b)
                     for(int i=0; i<gamma_maxp; ++i)
                       for(int j=0; j<gamma_max; ++j)
                         {
                           // Get target indices
                           int it=i*upsilon_maxp;
                           int jt=j*upsilon_max1;
                           // Get source indices
                           int is=i*upsilon_maxpp;
                           int js=j*upsilon_max1;

                           unit_tensor_blocks_omegapp[b].block(it,jt,upsilon_maxp,upsilon_max1)
                             +=A*unit_tensor_blocks_rhobp[b].block(is,js,upsilon_maxpp,upsilon_max1);
                         }

		} //omegapp

if(Nsum==Nsummy && i==ii)std::cout<<"unit_tensor_blocks_omega1[0]=unit_tensor_blocks_omegapp[0]="<<unit_tensor_blocks_omegapp[0]<<std::endl;

              // accumulating sum over omegapp
              for(int b=0; b<num_blocks; ++b)
                  unit_tensor_blocks_omega1[b]=unit_tensor_blocks_omegapp[b];

              //summing over x0'
              for (auto& x0p_mult : x0p_set) // x0p_mult is multiplicity tagged Gamma0'
                {
                  u3::SU3 x0p(x0p_mult.irrep); // x0p is Gamma0'
                                               
if(Nsum==Nsummy && i==ii)std::cout<<"Gamma0': "<<x0p.lambda()<<" "<<x0p.mu()<<std::endl;

                  int rho0p_max=OuterMultiplicity(omega1.SU3(),x0p,omegap.SU3()); // rho0p_max is maximal \bar{rho0}

                  // Zero initialize blocks accumlating sum over x0p
                  std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_x0p;
                  ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_x0p);

//                  for(int b=0; b<num_blocks; ++b)
//                    unit_tensor_blocks_x0p[b]+=unit_tensor_blocks_omegapp[b];

                  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  // second term
                  //////////////////////////////////////////////////////////////////////////////////////////////////////////
                  if(u3::OuterMultiplicity(u3::SU3(etap,0),u3::SU3(0,eta-2),x0p)>0)
                    {
                      assert((eta-2)>=0);
                      // look up index of subspace in unit tensor space
                      u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap,eta-2);

                      int unit_tensor_subspace_index1=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                      assert(unit_tensor_subspace_index1!=-1);

                      double coef1=(1.-(1./nucleon_number))*sqrt(1.*u3::dim(u3::SU3(eta,0)))*u3::dim(x0p)
                        *u3::UCached(u_coef_cache,u3::SU3(etap,0),u3::SU3(0,eta),x0p,u3::SU3(2,0),x0,1,1,u3::SU3(0,eta-2),1,1);

                      // zero initialize blocks for accumulating second term in sum over rhobp
                      std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                      ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                      // summing over rho0bp and accumulating sectors
                      for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp) // rho0bp is \bar{rho0}
                        {
                          int hypersector_index1=baby_spncci_hypersectors.LookUpHypersectorIndex(
                              baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                              unit_tensor_subspace_index1, rho0bp
                            );

                          // Accumulate
                          if(hypersector_index1==-1)
                            continue;

                          for(int b=0; b<num_blocks; ++b)
                            {
                               unit_tensor_blocks_rho0bp[b]
                                 +=u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(x0p),omega.SU3(),u3::SU3(2,0),
			                       omega1.SU3(),rho0bp,1,u3::Conjugate(x0),1,rho0)
                                   *unit_tensor_hyperblocks[hypersector_index1][b];
                            }
                          } //end rho0bp
		
if(Nsum==Nsummy && i==ii)std::cout<<"unit_tensor_blocks_x0p[0]=unit_tensor_blocks_x0p[0]+coef1*unit_tensor_blocks_rho0bp[0]="<<unit_tensor_blocks_x0p[0]<<"+"<<coef1<<"*"<<unit_tensor_blocks_rho0bp[0]<<"="<<unit_tensor_blocks_x0p[0]+coef1*unit_tensor_blocks_rho0bp[0]<<" (2nd term)"<<std::endl;

                        for(int b=0; b<num_blocks; ++b)
                          unit_tensor_blocks_x0p[b]+=coef1*unit_tensor_blocks_rho0bp[b];
                      }

                      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // third term
                      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                      if ((u3::OuterMultiplicity(u3::SU3(etap+2,0),u3::SU3(0,eta),x0p)>0) && (etap+2)<=Nmax+std::max(N1vp,N1vn))
                        {
                          // look up index of subspace in unit tensor space
                          u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap+2,eta);

                          int unit_tensor_subspace_index2=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                          assert(unit_tensor_subspace_index2!=-1);

                          double coef2=-1*(1.+(1./nucleon_number))*ParitySign(u3::ConjugationGrade(x0)+u3::ConjugationGrade(x0p))
                                  *u3::dim(u3::SU3(etap,0))*u3::dim(x0p)
                                  *u3::UCached(u_coef_cache,u3::SU3(2,0),u3::SU3(etap,0),x0p,u3::SU3(0,eta),
		         		 u3::SU3(etap+2,0),1,1,x0,1,1)/sqrt(double(u3::dim(u3::SU3(etap+2,0))));

                          // zero initialize blocks for accumulating first term in sum over rhobp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                          ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                          for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp) // rho0bp is \bar{rho0}
                            {
                              int hypersector_index2=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                    baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                    unit_tensor_subspace_index2, rho0bp
                                  );

                              // Accumulate
                              if(hypersector_index2==-1)
                                continue;

                              for(int b=0; b<num_blocks; ++b)
                                {
                                  unit_tensor_blocks_rho0bp[b]
                                    +=u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(x0p),omega.SU3(),u3::SU3(2,0),
						  omega1.SU3(),rho0bp,1,u3::Conjugate(x0),1,rho0)
                                      *unit_tensor_hyperblocks[hypersector_index2][b];
                                }
                            } //end rho0bp

if(Nsum==Nsummy && i==ii)std::cout<<"unit_tensor_blocks_x0p[0]=unit_tensor_blocks_x0p[0]+coef2*unit_tensor_blocks_rho0bp[0]="<<unit_tensor_blocks_x0p[0]<<"+"<<coef2<<"*"<<unit_tensor_blocks_rho0bp[0]<<"="<<unit_tensor_blocks_x0p[0]+coef2*unit_tensor_blocks_rho0bp[0]<<" (3rd term)"<<std::endl;

                          for(int b=0; b<num_blocks; ++b)
                            unit_tensor_blocks_x0p[b]+=coef2*unit_tensor_blocks_rho0bp[b];
                        }

                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // fourth term
                      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                      if((u3::OuterMultiplicity(u3::SU3(etap+1,0),u3::SU3(0,eta-1),x0p)>0) && (etap+1)<=Nmax+std::max(N1vp,N1vn))
                        {
                          assert((eta-1)>=0);
                          // look up index of subspace in unit tensor space
                          u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap+1,eta-1);

                          int unit_tensor_subspace_index3=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                          assert(unit_tensor_subspace_index3!=-1);

                          double coef4=0.0;
			  // x0 x (1,0) -> xpp
			  // (etap,0) x (0,eta-1) -> xpp
			  // (1,0) x xpp -> x0p
			  MultiplicityTagged<u3::SU3>::vector xpp_set=KroneckerProduct(x0,u3::SU3(1,0));
			  for (auto& xpp_mult : xpp_set) // xpp_mult is multiplicity tagged Gamma''
                            {
                              u3::SU3 xpp(xpp_mult.irrep); // xpp is Gamma''
			      if((u3::OuterMultiplicity(u3::SU3(etap,0),u3::SU3(0,eta-1),xpp)==0)
			         || (u3::OuterMultiplicity(u3::SU3(1,0),xpp,x0p)==0))
				continue;
                              coef4+=ParitySign(u3::ConjugationGrade(xpp))
				*u3::UCached(u_coef_cache,u3::SU3(etap,0),u3::SU3(0,eta),xpp,u3::SU3(1,0),x0,1,1,u3::SU3(0,eta-1),1,1)
				*u3::UCached(u_coef_cache,u3::SU3(1,0),u3::SU3(etap,0),x0p,u3::SU3(0,eta-1),u3::SU3(etap+1,0),1,1,xpp,1,1)
				*u3::UCached(u_coef_cache,x0,u3::SU3(1,0),x0p,u3::SU3(1,0),xpp,1,1,u3::SU3(2,0),1,1);
			    }
			  coef4=-1*ParitySign(u3::ConjugationGrade(x0p))*(etap+1)*u3::dim(x0p)
				  *sqrt(2.*(eta+2)/(etap+3))*coef4/nucleon_number;

                          // zero initialize blocks for accumulating second term in sum over rhobp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                            ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                          // summing over rho0bp and accumulating sectors
                          for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp) // rho0bp is \bar{rho0}
                            {
                              int hypersector_index4=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                    baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                    unit_tensor_subspace_index3, rho0bp
                                  );

                              // Accumulate
                              if(hypersector_index4==-1)
                                continue;

                              for(int b=0; b<num_blocks; ++b)
                              {
                                unit_tensor_blocks_rho0bp[b]
                                  +=u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(x0p),omega.SU3(),u3::SU3(2,0),
				     omega1.SU3(),rho0bp,1,u3::Conjugate(x0),1,rho0)
                                    *unit_tensor_hyperblocks[hypersector_index4][b];
                              }
                            } //end rho0bp

if(Nsum==Nsummy && i==ii)std::cout<<"unit_tensor_blocks_x0p[0]=unit_tensor_blocks_x0p[0]+coef4*unit_tensor_blocks_rho0bp[0]="<<unit_tensor_blocks_x0p[0]<<"+"<<coef4<<"*"<<unit_tensor_blocks_rho0bp[0]<<"="<<unit_tensor_blocks_x0p[0]+coef4*unit_tensor_blocks_rho0bp[0]<<" (4th term)"<<std::endl;

                          for(int b=0; b<num_blocks; ++b)
                            unit_tensor_blocks_x0p[b]+=coef4*unit_tensor_blocks_rho0bp[b];
                        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

if(Nsum==Nsummy && i==ii)std::cout<<"unit_tensor_blocks_omega1[0]=unit_tensor_blocks_omega1[0]+unit_tensor_blocks_x0p[0]="<<unit_tensor_blocks_omega1[0]<<"+"<<unit_tensor_blocks_x0p[0]<<"="<<unit_tensor_blocks_omega1[0]+unit_tensor_blocks_x0p[0]<<std::endl;

                        for(int b=0; b<num_blocks; ++b)
                          unit_tensor_blocks_omega1[b]+=unit_tensor_blocks_x0p[b];

                  }// end x0p sum
		
if(Nsum==Nsummy && i==ii)std::cout<<"unit_tensor_blocks_omega1[0]=unit_tensor_blocks_omega1[0]*coef="<<unit_tensor_blocks_omega1[0]<<"*"<<coef<<"="<<unit_tensor_blocks_omega1[0]*coef<<std::endl;

                for(int b=0; b<num_blocks; ++b)
                        unit_tensor_blocks_omega1[b]*=coef;

if(Nsum==Nsummy && i==ii)std::cout<<"unit_tensor_blocks[0]=unit_tensor_blocks[0]+unit_tensor_blocks_omega1[0]*KBUK="<<unit_tensor_blocks[0]<<"+"<<unit_tensor_blocks_omega1[0]<<"*"<<KBUK<<"="<<unit_tensor_blocks[0]+unit_tensor_blocks_omega1[0]*KBUK<<std::endl;

                // summing over n, rho, n1, rho1, upsilon1
                for(int b=0; b<num_blocks; ++b)
                  for(int i=0; i<gamma_maxp; ++i)
                    for(int j=0; j<gamma_max; ++j)
                    {
                      int it=i*upsilon_maxp;
                      int jt=j*upsilon_max;
                      int is=i*upsilon_maxp;
                      int js=j*upsilon_max1;
		      
                      unit_tensor_blocks[b].block(it,jt,upsilon_maxp,upsilon_max)
                        +=unit_tensor_blocks_omega1[b].block(is,js,upsilon_maxp,upsilon_max1)*KBUK;

                    }
                // done with blocks
              }// end omega1_mult
        }// end hypersector index
    }// end Nsum
  }

void DoOneBodyRecurrenceInitialization(
  int Nmax, int N1vp, int N1vn,
  const spncci::LGIPair& lgi_pair,
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::OneBodyUnitTensorSpaceU3S& one_body_unit_tensor_space,
  spncci::OperatorBlocks& lgi_transformations,
  bool transform_lgi_families,
  std::map<spncci::NnPair,std::set<int>>& unit_tensor_subspace_subsets,
  spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersector_seeds,
  spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersector_seeds_conj,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds_conj
  )
{
    // Extract lgi index  labels
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
    ///////////////////////////////////////////////////////////////////////////
    // Read in list of unit tensors between lgi pair and conjugates from files
    // Returned bool, files_found_test has no current use, but could be used
    //  to identify lgi pair with no non-zero rmes between them.
    // Corresponding rho0 values stored separately for later hypersector lookup
    ///////////////////////////////////////////////////////////////////////////
    // Initialize containers
    std::vector<u3shell::OneBodyUnitTensorLabelsU3S> lgi_unit_tensors;
    std::vector<int> rho0_values;

    // Get index corresponding to lgi in the full space.
    // Index may differ from lgi index in basis if space has been truncated
    int index1=lgi_full_space_index_lookup[irrep_family_index_bra];
    int index2=lgi_full_space_index_lookup[irrep_family_index_ket];

    // Read in operators
    std::string lgi_unit_tensor_filename
      =fmt::format("seeds/oboperators_{:06d}_{:06d}.dat",index1,index2);
    bool files_found_test=lgi::ReadOneBodyUnitTensorLabels(lgi_unit_tensor_filename,lgi_unit_tensors,rho0_values);

    ///////////////////////////////////////////////////////////////////////////
    // Set up hypersectors and hyperblocks for seeds
    // Generate hypersectors from list of unit tensors and outer-mulitplicites
    //  read from files
    // Populate hyperblocks using seeds read in from file
    //////////////////////////////////////////////////////////////////////////

    // Reads in unit tensor seed blocks and stores them in a vector of blocks. Order
    // corresponds to order of (unit_tensor,rho0) pairs in corresponding operator file.
    basis::OperatorBlocks<double> unit_tensor_seed_blocks;
    std::string seed_filename
      =fmt::format("seeds/obseeds_{:06d}_{:06d}.rmes",index1,index2);
    files_found_test&=lgi::ReadBlocks(seed_filename, lgi_unit_tensors.size(), unit_tensor_seed_blocks);

    // Identify unit tensor subspaces for recurrence
    spncci::GenerateRecurrenceOneBodyUnitTensors(
      Nmax,N1vp,N1vn,lgi_unit_tensors,
      one_body_unit_tensor_space,unit_tensor_subspace_subsets
    );

    //Generate hypersectors for unit tensors between lgi pair (lgi1,lgi2)
    baby_spncci_hypersector_seeds
      =spncci::BabySpNCCIOneBodyUnitTensorHypersectors(
        lgi_families,baby_spncci_space,one_body_unit_tensor_space,
        unit_tensor_subspace_subsets[spncci::NnPair(0,0)],
        irrep_family_index_bra,irrep_family_index_ket
      );

    // Generate hypersectors for unit tensors between conjugate lgi pair (lgi2,lgi1)
    baby_spncci_hypersector_seeds_conj
      =spncci::BabySpNCCIOneBodyUnitTensorHypersectors(
        lgi_families,baby_spncci_space,one_body_unit_tensor_space,
        unit_tensor_subspace_subsets[spncci::NnPair(0,0)],
        irrep_family_index_ket,irrep_family_index_bra
      );

    // Zero initialize seed hyperblocks and conjugate hyperblocks
    basis::SetHyperoperatorToZero(baby_spncci_hypersector_seeds,unit_tensor_hyperblocks_seeds);
    basis::SetHyperoperatorToZero(baby_spncci_hypersector_seeds_conj,unit_tensor_hyperblocks_seeds_conj);

    // Populate the hyperblocks and conjugate hyperblocks with the seed
    // Conjugate hyperspectors will be used in calculating Nn0 rmes in recurrence
    spncci::PopulateHypersectorsWithSeedsForOneBodyRecurrence(
      irrep_family_index_bra, irrep_family_index_ket,lgi_families,
      baby_spncci_space,one_body_unit_tensor_space,
      baby_spncci_hypersector_seeds_conj,baby_spncci_hypersector_seeds,
      lgi_unit_tensors,rho0_values,unit_tensor_seed_blocks,
      unit_tensor_hyperblocks_seeds_conj,unit_tensor_hyperblocks_seeds
    );
}

  bool GenerateOneBodyUnitTensorHyperblocks(
      const spncci::LGIPair& lgi_pair,
      int Nmax, int N1vp, int N1vn, int nucleon_number,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::OneBodyUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const spncci::KMatrixCache& kinv_matrix_cache,
      std::map<spncci::NnPair,std::set<int>>& unit_tensor_subspace_subsets,
      const spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersector_seeds,
      const spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersector_seeds_conj,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds_conj,
      u3::UCoefCache& u_coef_cache,
      u3::PhiCoefCache& phi_coef_cache,
      spncci::BabySpNCCIOneBodyUnitTensorHypersectors& baby_spncci_hypersectors,
      basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
      )
  {
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;

    // Generate Nn=0 hypersectors to be computed by conjugation
    bool Nn0_conjugate_hypersectors=true;
    std::vector<std::vector<int>> unit_tensor_hypersector_subsets_Nn0;

    spncci::BabySpNCCIOneBodyUnitTensorHypersectors baby_spncci_hypersectors_Nn0(
      Nmax, baby_spncci_space, unit_tensor_space,
      unit_tensor_subspace_subsets, unit_tensor_hypersector_subsets_Nn0,
      irrep_family_index_ket, irrep_family_index_bra,
      Nn0_conjugate_hypersectors
    );

    // Hypersectors for Nnp>=Nn
    Nn0_conjugate_hypersectors=false;
    std::vector<std::vector<int>> unit_tensor_hypersector_subsets;

    // (Nnp,Nn) sectors for Nnp>Nn
    // std::cout<<"generate Nnp>Nn hypersectors"<<std::endl;
    baby_spncci_hypersectors=spncci::BabySpNCCIOneBodyUnitTensorHypersectors(
      Nmax,baby_spncci_space, unit_tensor_space,
      unit_tensor_subspace_subsets, unit_tensor_hypersector_subsets,
      irrep_family_index_bra,irrep_family_index_ket,
      Nn0_conjugate_hypersectors
    );

    //Zero initialize hyperblocks for both Nn=0 and Nnp>=Nn sectors
    basis::OperatorHyperblocks<double> unit_tensor_hyperblocks_Nn0;
    basis::SetHyperoperatorToZero(baby_spncci_hypersectors_Nn0,unit_tensor_hyperblocks_Nn0);
    basis::SetHyperoperatorToZero(baby_spncci_hypersectors,unit_tensor_hyperblocks);

    //Iterate over seed hypersectors and identify corresponding hypersector in hyperblocks or hyperblocks_Nn0
    //Transfer seeds into hyperblocks for recurrence
    for(int seed_hypersector_index=0; seed_hypersector_index<baby_spncci_hypersector_seeds.size(); ++seed_hypersector_index)
      {
        int bra_index,ket_index,operator_index,multiplicity_index;
        std::tie(bra_index,ket_index,operator_index,multiplicity_index)
          =baby_spncci_hypersector_seeds.GetHypersector(seed_hypersector_index).Key();
        int hypersector_index
          =baby_spncci_hypersectors.LookUpHypersectorIndex(bra_index,ket_index,operator_index,multiplicity_index);

        unit_tensor_hyperblocks[hypersector_index]=unit_tensor_hyperblocks_seeds[seed_hypersector_index];
      }

    //Transfer conjugate hypersectors into Nn0 hypersectors for recurrence
    for(int seed_hypersector_index=0; seed_hypersector_index<baby_spncci_hypersector_seeds_conj.size(); ++seed_hypersector_index)
      {
        int bra_index,ket_index,operator_index,multiplicity_index;
        std::tie(bra_index,ket_index,operator_index,multiplicity_index)
          =baby_spncci_hypersector_seeds_conj.GetHypersector(seed_hypersector_index).Key();
        int hypersector_index
          =baby_spncci_hypersectors_Nn0.LookUpHypersectorIndex(bra_index,ket_index,operator_index,multiplicity_index);

        unit_tensor_hyperblocks_Nn0[hypersector_index]=unit_tensor_hyperblocks_seeds_conj[seed_hypersector_index];
      }
/*
//=============================================================================================
if(irrep_family_index_bra==0 && irrep_family_index_ket==0){
std::cout<<"*********************************** unit_tensor_hyperblocks_Nn0 before*************************************"<<std::endl;
      if(baby_spncci_hypersectors_Nn0.size()!=unit_tensor_hyperblocks_Nn0.size())
        std::cout<<"ERROR: baby_spncci_hypersectors_Nn0.size(),unit_tensor_hyperblocks_Nn0.size(): "
                 <<baby_spncci_hypersectors_Nn0.size()<<" "<<unit_tensor_hyperblocks_Nn0.size()<<std::endl;

      for (std::size_t hypersector_index=0; hypersector_index<baby_spncci_hypersectors_Nn0.size(); ++hypersector_index){
        auto key=baby_spncci_hypersectors_Nn0.GetHypersector(hypersector_index).Key();
        int unit_tensor_subspace_index, baby_spncci_subspace_indexp, baby_spncci_index, rho0;
        std::tie(baby_spncci_subspace_indexp,baby_spncci_index,unit_tensor_subspace_index,rho0)=key;
        const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);
        const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index);
        const u3shell::OneBodyUnitTensorSubspaceU3S& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
        int dim=baby_spncci_subspace_ket.dimension();
        int gamma_max=baby_spncci_subspace_ket.gamma_max();
        int upsilon_max=baby_spncci_subspace_ket.upsilon_max();
        if(dim!=gamma_max*upsilon_max)std::cout<<"ERROR: dim,gamma_max,upsilon_max: "<<dim<<" "<<gamma_max<<" "<<upsilon_max<<std::endl;
        int dimp=baby_spncci_subspace_bra.dimension();
        int gamma_maxp=baby_spncci_subspace_bra.gamma_max();
        int upsilon_maxp=baby_spncci_subspace_bra.upsilon_max();
        if(dimp!=gamma_maxp*upsilon_maxp)std::cout<<"ERROR: dimp,gamma_maxp,upsilon_maxp: "<<dimp<<" "<<gamma_maxp<<" "<<upsilon_maxp<<std::endl;
        u3::U3 omegap,sigmap,omega,sigma; // p denotes prime. bra has primed quantum numbers
        u3::SU3 x0;
        HalfInt S0,Sn_ket,Sp_ket,S_ket,Sn_bra,Sp_bra,S_bra;
        int etap,eta;
        std::tie(sigmap,Sp_bra,Sn_bra,S_bra,omegap)=baby_spncci_subspace_bra.labels();
        std::tie(sigma,Sp_ket,Sn_ket,S_ket,omega)=baby_spncci_subspace_ket.labels();
        std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
        std::cout<<"hypersector:"<<std::endl;
        std::cout<<"N_sigma' lambda_sigma' mu_sigma' Sp' Sn' S' N_omega' lambda_omega' mu_omega': "<<sigmap.N()<<" "<<sigmap.SU3().lambda()<<" "<<sigmap.SU3().mu()<<" "<<Sp_bra<<" "<<Sn_bra<<" "<<S_bra<<" "<<omegap.N()<<" "<<omegap.SU3().lambda()<<" "<<omegap.SU3().mu()<<std::endl;
        std::cout<<"N_sigma  lambda_sigma  mu_sigma  Sp  Sn  S  N_omega  lambda_omega  mu_omega : "<<sigma.N()<<" "<<sigma.SU3().lambda()<<" "<<sigma.SU3().mu()<<" "<<Sp_ket<<" "<<Sn_ket<<" "<<S_ket<<" "<<omega.N()<<" "<<omega.SU3().lambda()<<" "<<omega.SU3().mu()<<std::endl;
        std::cout<<"lambda0 mu0 S0 N' N rho0: "<<x0.lambda()<<" "<<x0.mu()<<" "<<S0<<" "<<etap<<" "<<eta<<" "<<rho0<<std::endl;
        if(unit_tensor_hyperblocks_Nn0[hypersector_index].size()>2)std::cout<<"ERROR: unit_tensor_hyperblocks_Nn0["<<hypersector_index<<"].size()="<<unit_tensor_hyperblocks_Nn0[hypersector_index].size()<<std::endl;

        for(int operator_index=0; operator_index<unit_tensor_hyperblocks_Nn0[hypersector_index].size(); operator_index++){
          std::cout<<"Operator index="<<operator_index<<std::endl;
          if(unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index].rows()!=dimp)std::cout<<"ERROR: unit_tensor_hyperblocks_Nn0["<<hypersector_index<<"]["<<operator_index<<"].rows(),dimp: "
                  <<unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index].rows()<<" "<<dimp<<std::endl;
          if(unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index].cols()!=dim)std::cout<<"ERROR: unit_tensor_hyperblocks_Nn0["<<hypersector_index<<"]["<<operator_index<<"].cols(),dim: "
                  <<unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index].cols()<<" "<<dim<<std::endl;

          for(int i=0; i<gamma_maxp; ++i)
            for(int j=0; j<gamma_max; ++j){
//            std::cout<<"gammap,gamma="<<i+1<<" "<<j+1<<std::endl;
              int it=i*upsilon_maxp;
              int jt=j*upsilon_max;
              // unit_tensor_hyperblocks[hypersector_index][operator_index].block(it,jt,upsilon_maxp,upsilon_max) is block corresponding to given gammap,gamma.
              // Its elements correspond to different upsilonp,upsilon.
              for(int up=0; up<upsilon_maxp; up++)
                for(int u=0; u<upsilon_max; u++){
                  std::cout<<"gamma' upsilon' gamma upsilon RME: "<<i+1<<" "<<up+1<<" "<<j+1<<" "<<u+1<<" "
                           <<unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index](it+up,jt+u)<<std::endl;
              }
          }
        }
      }
std::cout<<"***********************************************************************************************************"<<std::endl;
}
//=============================================================================================
*/
    // std::cout<<"Compute Nn=0 blocks"<<std::endl;
    spncci::ComputeOneBodyUnitTensorHyperblocks(
      Nmax,N1vp,N1vn,nucleon_number,u_coef_cache,phi_coef_cache,
      k_matrix_cache,kinv_matrix_cache,spncci_space,baby_spncci_space,
      unit_tensor_space,baby_spncci_hypersectors_Nn0,
      unit_tensor_hypersector_subsets_Nn0,unit_tensor_hyperblocks_Nn0
    );
/*
//=============================================================================================
if(irrep_family_index_bra==2 && irrep_family_index_ket==0){
std::cout<<"*********************************** unit_tensor_hyperblocks_Nn0 after**************************************"<<std::endl;
      if(baby_spncci_hypersectors_Nn0.size()!=unit_tensor_hyperblocks_Nn0.size())
        std::cout<<"ERROR: baby_spncci_hypersectors_Nn0.size(),unit_tensor_hyperblocks_Nn0.size(): "
                 <<baby_spncci_hypersectors_Nn0.size()<<" "<<unit_tensor_hyperblocks_Nn0.size()<<std::endl;

      for (std::size_t hypersector_index=0; hypersector_index<baby_spncci_hypersectors_Nn0.size(); ++hypersector_index){
        auto key=baby_spncci_hypersectors_Nn0.GetHypersector(hypersector_index).Key();
        int unit_tensor_subspace_index, baby_spncci_subspace_indexp, baby_spncci_index, rho0;
        std::tie(baby_spncci_subspace_indexp,baby_spncci_index,unit_tensor_subspace_index,rho0)=key;
        const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);
        const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index);
        const u3shell::OneBodyUnitTensorSubspaceU3S& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
        int dim=baby_spncci_subspace_ket.dimension();
        int gamma_max=baby_spncci_subspace_ket.gamma_max();
        int upsilon_max=baby_spncci_subspace_ket.upsilon_max();
        if(dim!=gamma_max*upsilon_max)std::cout<<"ERROR: dim,gamma_max,upsilon_max: "<<dim<<" "<<gamma_max<<" "<<upsilon_max<<std::endl;
        int dimp=baby_spncci_subspace_bra.dimension();
        int gamma_maxp=baby_spncci_subspace_bra.gamma_max();
        int upsilon_maxp=baby_spncci_subspace_bra.upsilon_max();
        if(dimp!=gamma_maxp*upsilon_maxp)std::cout<<"ERROR: dimp,gamma_maxp,upsilon_maxp: "<<dimp<<" "<<gamma_maxp<<" "<<upsilon_maxp<<std::endl;
        u3::U3 omegap,sigmap,omega,sigma; // p denotes prime. bra has primed quantum numbers
        u3::SU3 x0;
        HalfInt S0,Sn_ket,Sp_ket,S_ket,Sn_bra,Sp_bra,S_bra;
        int etap,eta;
        std::tie(sigmap,Sp_bra,Sn_bra,S_bra,omegap)=baby_spncci_subspace_bra.labels();
        std::tie(sigma,Sp_ket,Sn_ket,S_ket,omega)=baby_spncci_subspace_ket.labels();
        std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
        std::cout<<"hypersector:"<<std::endl;
        std::cout<<"N_sigma' lambda_sigma' mu_sigma' Sp' Sn' S' N_omega' lambda_omega' mu_omega': "<<sigmap.N()<<" "<<sigmap.SU3().lambda()<<" "<<sigmap.SU3().mu()<<" "<<Sp_bra<<" "<<Sn_bra<<" "<<S_bra<<" "<<omegap.N()<<" "<<omegap.SU3().lambda()<<" "<<omegap.SU3().mu()<<std::endl;
        std::cout<<"N_sigma  lambda_sigma  mu_sigma  Sp  Sn  S  N_omega  lambda_omega  mu_omega : "<<sigma.N()<<" "<<sigma.SU3().lambda()<<" "<<sigma.SU3().mu()<<" "<<Sp_ket<<" "<<Sn_ket<<" "<<S_ket<<" "<<omega.N()<<" "<<omega.SU3().lambda()<<" "<<omega.SU3().mu()<<std::endl;
        std::cout<<"lambda0 mu0 S0 N' N rho0: "<<x0.lambda()<<" "<<x0.mu()<<" "<<S0<<" "<<etap<<" "<<eta<<" "<<rho0<<std::endl;
        if(unit_tensor_hyperblocks_Nn0[hypersector_index].size()>2)std::cout<<"ERROR: unit_tensor_hyperblocks_Nn0["<<hypersector_index<<"].size()="<<unit_tensor_hyperblocks_Nn0[hypersector_index].size()<<std::endl;

        for(int operator_index=0; operator_index<unit_tensor_hyperblocks_Nn0[hypersector_index].size(); operator_index++){
          std::cout<<"Operator index="<<operator_index<<std::endl;
          if(unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index].rows()!=dimp)std::cout<<"ERROR: unit_tensor_hyperblocks_Nn0["<<hypersector_index<<"]["<<operator_index<<"].rows(),dimp: "
                  <<unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index].rows()<<" "<<dimp<<std::endl;
          if(unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index].cols()!=dim)std::cout<<"ERROR: unit_tensor_hyperblocks_Nn0["<<hypersector_index<<"]["<<operator_index<<"].cols(),dim: "
                  <<unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index].cols()<<" "<<dim<<std::endl;

          for(int i=0; i<gamma_maxp; ++i)
            for(int j=0; j<gamma_max; ++j){
//            std::cout<<"gammap,gamma="<<i+1<<" "<<j+1<<std::endl;
              int it=i*upsilon_maxp;
              int jt=j*upsilon_max;
              // unit_tensor_hyperblocks[hypersector_index][operator_index].block(it,jt,upsilon_maxp,upsilon_max) is block corresponding to given gammap,gamma.
              // Its elements correspond to different upsilonp,upsilon.
              for(int up=0; up<upsilon_maxp; up++)
                for(int u=0; u<upsilon_max; u++){
                  std::cout<<"gamma' upsilon' gamma upsilon RME: "<<i+1<<" "<<up+1<<" "<<j+1<<" "<<u+1<<" "
                           <<unit_tensor_hyperblocks_Nn0[hypersector_index][operator_index](it+up,jt+u)<<std::endl;
              }
          }
        }
      }
std::cout<<"***********************************************************************************************************"<<std::endl;
}
//=============================================================================================
*/
    // std::cout<<"Add Nn0 blocks to hyperblocks"<<std::endl;
    spncci::AddNn0BlocksToOneBodyUnitTensorHyperblocks(
      baby_spncci_space,unit_tensor_space,
      baby_spncci_hypersectors_Nn0,baby_spncci_hypersectors,
      unit_tensor_hyperblocks_Nn0,unit_tensor_hyperblocks
    );
/*
//=============================================================================================
if(irrep_family_index_bra==0 && irrep_family_index_ket==0){
std::cout<<"*********************************** unit_tensor_hyperblocks before**************************************"<<std::endl;
      if(baby_spncci_hypersectors.size()!=unit_tensor_hyperblocks.size())
        std::cout<<"ERROR: baby_spncci_hypersectors.size(),unit_tensor_hyperblocks.size(): "
                 <<baby_spncci_hypersectors.size()<<" "<<unit_tensor_hyperblocks.size()<<std::endl;

      for (std::size_t hypersector_index=0; hypersector_index<baby_spncci_hypersectors.size(); ++hypersector_index){
        auto key=baby_spncci_hypersectors.GetHypersector(hypersector_index).Key();
        int unit_tensor_subspace_index, baby_spncci_subspace_indexp, baby_spncci_index, rho0;
        std::tie(baby_spncci_subspace_indexp,baby_spncci_index,unit_tensor_subspace_index,rho0)=key;
        const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);
        const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index);
        const u3shell::OneBodyUnitTensorSubspaceU3S& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
        int dim=baby_spncci_subspace_ket.dimension();
        int gamma_max=baby_spncci_subspace_ket.gamma_max();
        int upsilon_max=baby_spncci_subspace_ket.upsilon_max();
        if(dim!=gamma_max*upsilon_max)std::cout<<"ERROR: dim,gamma_max,upsilon_max: "<<dim<<" "<<gamma_max<<" "<<upsilon_max<<std::endl;
        int dimp=baby_spncci_subspace_bra.dimension();
        int gamma_maxp=baby_spncci_subspace_bra.gamma_max();
        int upsilon_maxp=baby_spncci_subspace_bra.upsilon_max();
        if(dimp!=gamma_maxp*upsilon_maxp)std::cout<<"ERROR: dimp,gamma_maxp,upsilon_maxp: "<<dimp<<" "<<gamma_maxp<<" "<<upsilon_maxp<<std::endl;
        u3::U3 omegap,sigmap,omega,sigma; // p denotes prime. bra has primed quantum numbers
        u3::SU3 x0;
        HalfInt S0,Sn_ket,Sp_ket,S_ket,Sn_bra,Sp_bra,S_bra;
        int etap,eta;
        std::tie(sigmap,Sp_bra,Sn_bra,S_bra,omegap)=baby_spncci_subspace_bra.labels();
        std::tie(sigma,Sp_ket,Sn_ket,S_ket,omega)=baby_spncci_subspace_ket.labels();
        std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
        std::cout<<"hypersector:"<<std::endl;
        std::cout<<"N_sigma' lambda_sigma' mu_sigma' Sp' Sn' S' N_omega' lambda_omega' mu_omega': "<<sigmap.N()<<" "<<sigmap.SU3().lambda()<<" "<<sigmap.SU3().mu()<<" "<<Sp_bra<<" "<<Sn_bra<<" "<<S_bra<<" "<<omegap.N()<<" "<<omegap.SU3().lambda()<<" "<<omegap.SU3().mu()<<std::endl;
        std::cout<<"N_sigma  lambda_sigma  mu_sigma  Sp  Sn  S  N_omega  lambda_omega  mu_omega : "<<sigma.N()<<" "<<sigma.SU3().lambda()<<" "<<sigma.SU3().mu()<<" "<<Sp_ket<<" "<<Sn_ket<<" "<<S_ket<<" "<<omega.N()<<" "<<omega.SU3().lambda()<<" "<<omega.SU3().mu()<<std::endl;
        std::cout<<"lambda0 mu0 S0 N' N rho0: "<<x0.lambda()<<" "<<x0.mu()<<" "<<S0<<" "<<etap<<" "<<eta<<" "<<rho0<<std::endl;
        if(unit_tensor_hyperblocks[hypersector_index].size()>2)std::cout<<"ERROR: unit_tensor_hyperblocks["<<hypersector_index<<"].size()="<<unit_tensor_hyperblocks[hypersector_index].size()<<std::endl;

        for(int operator_index=0; operator_index<unit_tensor_hyperblocks[hypersector_index].size(); operator_index++){
          std::cout<<"Operator index="<<operator_index<<std::endl;
          if(unit_tensor_hyperblocks[hypersector_index][operator_index].rows()!=dimp)std::cout<<"ERROR: unit_tensor_hyperblocks["<<hypersector_index<<"]["<<operator_index<<"].rows(),dimp: "
                  <<unit_tensor_hyperblocks[hypersector_index][operator_index].rows()<<" "<<dimp<<std::endl;
          if(unit_tensor_hyperblocks[hypersector_index][operator_index].cols()!=dim)std::cout<<"ERROR: unit_tensor_hyperblocks["<<hypersector_index<<"]["<<operator_index<<"].cols(),dim: "
                  <<unit_tensor_hyperblocks[hypersector_index][operator_index].cols()<<" "<<dim<<std::endl;

          for(int i=0; i<gamma_maxp; ++i)
            for(int j=0; j<gamma_max; ++j){
//            std::cout<<"gammap,gamma="<<i+1<<" "<<j+1<<std::endl;
              int it=i*upsilon_maxp;
              int jt=j*upsilon_max;
              // unit_tensor_hyperblocks[hypersector_index][operator_index].block(it,jt,upsilon_maxp,upsilon_max) is block corresponding to given gammap,gamma.
              // Its elements correspond to different upsilonp,upsilon.
              for(int up=0; up<upsilon_maxp; up++)
                for(int u=0; u<upsilon_max; u++){
                  std::cout<<"gamma' upsilon' gamma upsilon RME: "<<i+1<<" "<<up+1<<" "<<j+1<<" "<<u+1<<" "
                           <<unit_tensor_hyperblocks[hypersector_index][operator_index](it+up,jt+u)<<std::endl;
              }
          }
        }
      }
std::cout<<"***********************************************************************************************************"<<std::endl;
}
//=============================================================================================
*/
    // std::cout<<"Compute unit tensor hyperblocks"<<std::endl;
    spncci::ComputeOneBodyUnitTensorHyperblocks(
      Nmax,N1vp,N1vn,nucleon_number,u_coef_cache,phi_coef_cache,
      k_matrix_cache,kinv_matrix_cache,spncci_space,baby_spncci_space,
      unit_tensor_space,baby_spncci_hypersectors,
      unit_tensor_hypersector_subsets,unit_tensor_hyperblocks
    );
/*
//=============================================================================================
if(irrep_family_index_bra==0 && irrep_family_index_ket==0){
std::cout<<"*********************************** unit_tensor_hyperblocks after**************************************"<<std::endl;
      if(baby_spncci_hypersectors.size()!=unit_tensor_hyperblocks.size())
        std::cout<<"ERROR: baby_spncci_hypersectors.size(),unit_tensor_hyperblocks.size(): "
                 <<baby_spncci_hypersectors.size()<<" "<<unit_tensor_hyperblocks.size()<<std::endl;

      for (std::size_t hypersector_index=0; hypersector_index<baby_spncci_hypersectors.size(); ++hypersector_index){
        auto key=baby_spncci_hypersectors.GetHypersector(hypersector_index).Key();
        int unit_tensor_subspace_index, baby_spncci_subspace_indexp, baby_spncci_index, rho0;
        std::tie(baby_spncci_subspace_indexp,baby_spncci_index,unit_tensor_subspace_index,rho0)=key;
        const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);
        const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index);
        const u3shell::OneBodyUnitTensorSubspaceU3S& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
        int dim=baby_spncci_subspace_ket.dimension();
        int gamma_max=baby_spncci_subspace_ket.gamma_max();
        int upsilon_max=baby_spncci_subspace_ket.upsilon_max();
        if(dim!=gamma_max*upsilon_max)std::cout<<"ERROR: dim,gamma_max,upsilon_max: "<<dim<<" "<<gamma_max<<" "<<upsilon_max<<std::endl;
        int dimp=baby_spncci_subspace_bra.dimension();
        int gamma_maxp=baby_spncci_subspace_bra.gamma_max();
        int upsilon_maxp=baby_spncci_subspace_bra.upsilon_max();
        if(dimp!=gamma_maxp*upsilon_maxp)std::cout<<"ERROR: dimp,gamma_maxp,upsilon_maxp: "<<dimp<<" "<<gamma_maxp<<" "<<upsilon_maxp<<std::endl;
        u3::U3 omegap,sigmap,omega,sigma; // p denotes prime. bra has primed quantum numbers
        u3::SU3 x0;
        HalfInt S0,Sn_ket,Sp_ket,S_ket,Sn_bra,Sp_bra,S_bra;
        int etap,eta;
        std::tie(sigmap,Sp_bra,Sn_bra,S_bra,omegap)=baby_spncci_subspace_bra.labels();
        std::tie(sigma,Sp_ket,Sn_ket,S_ket,omega)=baby_spncci_subspace_ket.labels();
        std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
        std::cout<<"hypersector:"<<std::endl;
        std::cout<<"N_sigma' lambda_sigma' mu_sigma' Sp' Sn' S' N_omega' lambda_omega' mu_omega': "<<sigmap.N()<<" "<<sigmap.SU3().lambda()<<" "<<sigmap.SU3().mu()<<" "<<Sp_bra<<" "<<Sn_bra<<" "<<S_bra<<" "<<omegap.N()<<" "<<omegap.SU3().lambda()<<" "<<omegap.SU3().mu()<<std::endl;
        std::cout<<"N_sigma  lambda_sigma  mu_sigma  Sp  Sn  S  N_omega  lambda_omega  mu_omega : "<<sigma.N()<<" "<<sigma.SU3().lambda()<<" "<<sigma.SU3().mu()<<" "<<Sp_ket<<" "<<Sn_ket<<" "<<S_ket<<" "<<omega.N()<<" "<<omega.SU3().lambda()<<" "<<omega.SU3().mu()<<std::endl;
        std::cout<<"lambda0 mu0 S0 N' N rho0: "<<x0.lambda()<<" "<<x0.mu()<<" "<<S0<<" "<<etap<<" "<<eta<<" "<<rho0<<std::endl;
        if(unit_tensor_hyperblocks[hypersector_index].size()>2)std::cout<<"ERROR: unit_tensor_hyperblocks["<<hypersector_index<<"].size()="<<unit_tensor_hyperblocks[hypersector_index].size()<<std::endl;

        for(int operator_index=0; operator_index<unit_tensor_hyperblocks[hypersector_index].size(); operator_index++){
          std::cout<<"Operator index="<<operator_index<<std::endl;
          if(unit_tensor_hyperblocks[hypersector_index][operator_index].rows()!=dimp)std::cout<<"ERROR: unit_tensor_hyperblocks["<<hypersector_index<<"]["<<operator_index<<"].rows(),dimp: "
                  <<unit_tensor_hyperblocks[hypersector_index][operator_index].rows()<<" "<<dimp<<std::endl;
          if(unit_tensor_hyperblocks[hypersector_index][operator_index].cols()!=dim)std::cout<<"ERROR: unit_tensor_hyperblocks["<<hypersector_index<<"]["<<operator_index<<"].cols(),dim: "
                  <<unit_tensor_hyperblocks[hypersector_index][operator_index].cols()<<" "<<dim<<std::endl;

          for(int i=0; i<gamma_maxp; ++i)
            for(int j=0; j<gamma_max; ++j){
//            std::cout<<"gammap,gamma="<<i+1<<" "<<j+1<<std::endl;
              int it=i*upsilon_maxp;
              int jt=j*upsilon_max;
              // unit_tensor_hyperblocks[hypersector_index][operator_index].block(it,jt,upsilon_maxp,upsilon_max) is block corresponding to given gammap,gamma.
              // Its elements correspond to different upsilonp,upsilon.
              for(int up=0; up<upsilon_maxp; up++)
                for(int u=0; u<upsilon_max; u++){
                  std::cout<<"gamma' upsilon' gamma upsilon RME: "<<i+1<<" "<<up+1<<" "<<j+1<<" "<<u+1<<" "
                           <<unit_tensor_hyperblocks[hypersector_index][operator_index](it+up,jt+u)<<std::endl;
              }
          }
        }
      }
std::cout<<"***********************************************************************************************************"<<std::endl;
}
//=============================================================================================
*/
    return true;
  }

void AddNn0BlocksToTwoBodyDensityHyperblocks(
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::TwoBodyDensitySpace& tbd_space,
  const spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersectors_Nn0,
  const spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersectors,
  basis::OperatorHyperblocks<double>& tbd_hyperblocks_Nn0,
  basis::OperatorHyperblocks<double>& tbd_hyperblocks,
  u3::PhiCoefCache& phi_coef_cache
)
{
  for(int hypersector_index_Nn0=0; hypersector_index_Nn0<baby_spncci_hypersectors_Nn0.size(); ++hypersector_index_Nn0)
    {
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Extracting labels from source (Nn0 sectors)
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Get hypersector indices for Nn0 sector
      // Note the labels from Nn0 sectors are conjugate labels, i.e., bra is actually ket and vice versa
      auto key=baby_spncci_hypersectors_Nn0.GetHypersector(hypersector_index_Nn0).Key();
      int tbd_subspace_index_Nn0, baby_spncci_index_bra, baby_spncci_index_ket, rho;
      std::tie(baby_spncci_index_ket,baby_spncci_index_bra,tbd_subspace_index_Nn0,rho)=key;

      // Extract TBD subspace labels from Nn0 tensor
      auto& tbd_subspace_Nn0=tbd_space.GetSubspace(tbd_subspace_index_Nn0);
      u3::SU3 x0c;
      HalfInt S0;
      int N1,N2,N3,N4;
      std::tie(x0c,S0,N4,N3,N2,N1)=tbd_subspace_Nn0.labels();

      // Get bra and ket labels from Nn0 sector
      const spncci::BabySpNCCISubspace& subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra);
      const spncci::BabySpNCCISubspace& subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket);

      u3::U3 omegap=subspace_bra.omega();
      HalfInt Sp=subspace_bra.S();

      u3::U3 omega=subspace_ket.omega();
      HalfInt S=subspace_ket.S();

      double conjugation_factor_base
              =ParitySign(u3::ConjugationGrade(omega)+S-u3::ConjugationGrade(omegap)-u3::ConjugationGrade(x0c)-Sp+N1+N2+N3+N4)
                *sqrt(double(u3::dim(omega)*am::dim(S))/double(u3::dim(omegap)*am::dim(Sp)));

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Looking up target hypersector
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Get TBD subspace index in hyperblocks
      u3::SU3 x0(u3::Conjugate(x0c));
      u3shell::TwoBodyDensitySubspaceLabels tbd_labels=u3shell::TwoBodyDensitySubspaceLabels(x0,S0,N1,N2,N3,N4);

      int tbd_subspace_index=tbd_space.LookUpSubspaceIndex(tbd_labels);
      auto& tbd_subspace=tbd_space.GetSubspace(tbd_subspace_index);

      // Look up hypersector
      int hypersector_index
          =baby_spncci_hypersectors.LookUpHypersectorIndex(
              baby_spncci_index_bra,baby_spncci_index_ket,tbd_subspace_index,rho);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // For each source hyperblock, identify target block and conjugate
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for(int tbd_index_Nn0=0; tbd_index_Nn0<tbd_subspace_Nn0.size(); ++tbd_index_Nn0)
        {
	  u3::SU3 xfc,xic;
          int Sf,Si,rho0c,Tz;
          std::tie(xic,Si,xfc,Sf,rho0c,Tz)=tbd_subspace_Nn0.GetStateLabels(tbd_index_Nn0);

	  u3::SU3 xf(u3::Conjugate(xfc));
	  u3::SU3 xi(u3::Conjugate(xic));
	  int rho0max=u3::OuterMultiplicity(xic,xfc,x0c);
	  double conjugation_factor=ParitySign(u3::ConjugationGrade(xf)+u3::ConjugationGrade(xi)+rho0max-rho0c)*conjugation_factor_base;
	  for(int rho0=1; rho0<=rho0max; rho0++){
            // Get TBD index
            std::tuple<u3::SU3,int,u3::SU3,int,int,int> state_labels(xf,Sf,xi,Si,rho0,Tz);
            int tbd_index=tbd_subspace.LookUpStateIndex(state_labels);

            tbd_hyperblocks[hypersector_index][tbd_index]+=conjugation_factor*u3::PhiCached(phi_coef_cache,xf,xi,x0,rho0,rho0c)
		    *tbd_hyperblocks_Nn0[hypersector_index_Nn0][tbd_index_Nn0].transpose();
	  }
        }
    }
}

void ComputeTwoBodyDensityHyperblocks(
  int Nmax, int N1vp, int N1vn, int nucleon_number,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::KMatrixCache& k_matrix_map,
  const spncci::KMatrixCache& kinv_matrix_map,
  const spncci::SpNCCISpace& spncci_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::TwoBodyDensitySpace& tbd_space,
  const spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersectors,
  const std::vector<std::vector<int>>& tbd_hypersector_subsets,
  basis::OperatorHyperblocks<double>& tbd_hyperblocks // Output
  // basis::OperatorHyperblocks<double> is std::vector<std::vector<OperatorBlock<double>>>
  // OperatorBlock<double> is Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>
  )
// compute hyperblocks for TBDs recursively
{
/*
  for(int Nsum=2; Nsum<=2*Nmax; Nsum+=2) // Nsum=Nn+Nn', where Nn=N_omega-N_sigma and Nn'=N_omega'-N_sigma'
    {
      const std::vector<int>& unit_tensor_hypersectors=unit_tensor_hypersector_subsets[Nsum/2];

      for(int i=0; i<unit_tensor_hypersectors.size(); ++i) // loop over hypersectors
      // Hypersector is given by the first index of unit_tensor_hyperblocks
      // Hypersector corresponds to sigma,Sp,Sn,S,omega,sigma',Sp',Sn',S',omega',omega0,S0,eta',eta,rho0
        {
          int hypersector_index=unit_tensor_hypersectors[i];
          auto key=baby_spncci_hypersectors.GetHypersector(hypersector_index).Key();

          int unit_tensor_subspace_index, baby_spncci_subspace_indexp, baby_spncci_index, rho0;
          std::tie(baby_spncci_subspace_indexp,baby_spncci_index,unit_tensor_subspace_index,rho0)=key;

          const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
              =baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);

          const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
              =baby_spncci_space.GetSubspace(baby_spncci_index);

          const u3shell::OneBodyUnitTensorSubspaceU3S& unit_tensor_subspace
              =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

          // Extract ket dimensions and mulitplicities
          int dim=baby_spncci_subspace_ket.dimension();
          int gamma_max=baby_spncci_subspace_ket.gamma_max();
          int upsilon_max=baby_spncci_subspace_ket.upsilon_max();

          // Extract bra dimensions and mulitplicities
          int dimp=baby_spncci_subspace_bra.dimension();
          int gamma_maxp=baby_spncci_subspace_bra.gamma_max();
          int upsilon_maxp=baby_spncci_subspace_bra.upsilon_max();

          // Extract Sp(3,R) space
          int irrep_family_index_bra=baby_spncci_subspace_bra.irrep_family_index();
          int irrep_family_index_ket=baby_spncci_subspace_ket.irrep_family_index();

          const sp3r::Sp3RSpace& irrep_bra = spncci_space[irrep_family_index_bra].Sp3RSpace();
          const sp3r::Sp3RSpace& irrep_ket = spncci_space[irrep_family_index_ket].Sp3RSpace();

          // extract subspace labels
          u3::U3 omegap,sigmap,omega,sigma; // p denotes prime. bra has primed quantum numbers
          u3::SU3 x0; // x0 is Gamma0
          HalfInt S0,Sn_ket,Sp_ket,S_ket,Sn_bra,Sp_bra,S_bra;
          int etap,eta;

          // Extracting labels
          std::tie(sigmap,Sp_bra,Sn_bra,S_bra,omegap)=baby_spncci_subspace_bra.labels();
          std::tie(sigma,Sp_ket,Sn_ket,S_ket,omega)=baby_spncci_subspace_ket.labels();
          std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
          int Nn=baby_spncci_subspace_ket.Nn(); // Nn is N_omega-N_sigma

          // omega u3 subspace in irrep
          const sp3r::U3Subspace& u3_subspace=irrep_ket.LookUpSubspace(omega);
          const sp3r::U3Subspace& u3_subspacep=irrep_bra.LookUpSubspace(omegap);
          const vcs::MatrixCache& K_matrix_map_ket=k_matrix_map.at(sigma);
          const vcs::MatrixCache& Kinv_matrix_map_ket=kinv_matrix_map.at(sigma);
          const Eigen::MatrixXd& K_inv=Kinv_matrix_map_ket.at(omega);

          // Generate labels to sum over
          int rho0_max=u3::OuterMultiplicity(omega.SU3(),x0,omegap.SU3());

          // Precalculating kronecker products used in sum to calculate unit tensor matrix
          MultiplicityTagged<u3::U3>::vector omegapp_set=KroneckerProduct(omegap, u3::U3(0,0,-2));
	  // omegapp_set is vector of multiplicity tagged \bar{omega}
          MultiplicityTagged<u3::U3>::vector omega1_set=KroneckerProduct(omega, u3::U3(0,0,-2));
	  // omega1_set is vector of multiplicity tagged omega1
          MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(x0, u3::SU3(2,0));
	  // x0p_set is vector of multiplicity tagged SU(3) irreps of omega0' (N_omega0'=N_omega'-N_omega1)

          std::vector<basis::OperatorBlock<double>>& unit_tensor_blocks=unit_tensor_hyperblocks[hypersector_index];
	  // unit_tensor_blocks is vector of matrices (blocks, sectors) corresponding to given hypersector
	  // Matrix indices of block correspond to gamma,upsilon (whose number is dim) and gamma',upsilon' (whose number is dimp)
	  // Within hypersector blocks correspond to different operators (unit tensors)
	  // For one-body unit tensors there are at most 2 blocks in each hypersector: proton and/or neutron
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //  Calculate unit tensor matrix
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          int num_blocks=unit_tensor_blocks.size();

          for(auto& omega1_mult :omega1_set) // omega1_mult is multiplicity tagged omega1
            {
              u3::U3 omega1(omega1_mult.irrep); // omega1 is omega1

              if (not irrep_ket.ContainsSubspace(omega1)) // if omega1 is not present in Sp(3,R) irrep sigma
                continue;

              spncci::BabySpNCCISubspaceLabels baby_spncci_labels1(sigma,Sp_ket,Sn_ket,S_ket,omega1);
              int baby_spncci_subspace_index1=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labels1);
              Eigen::MatrixXd K1=K_matrix_map_ket.at(omega1); // K1 is K_sigma^omega1
              const sp3r::U3Subspace& u3_subspace1=irrep_ket.LookUpSubspace(omega1);
              int upsilon_max1=u3_subspace1.upsilon_max(); // upsilon_max1 is maximal upsilon1

              int dim1=upsilon_max1*gamma_max; // gamma_max is maximal gamma (for ket, i.e., sigma)
              std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omega1;
              // Initializing blocks for sum over omega1
              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omega1);

              // Construct KBUK matrix
              Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
              spncci::ConsructKBUK_JH( // Computes chi[upsilon1 omega1;upsilon omega;(20)]=KBUK(upsilon1-1,upsilon-1) for given omega1,omega,sigma
                u_coef_cache, Nn,sigma, omega, omega1,
                u3_subspace,u3_subspace1,K1,K_inv,
                KBUK,write
              );

              double coef=sqrt(1.*u3::dim(omega.SU3())/(1.*u3::dim(omega1.SU3())))/double(u3::dim(x0));

              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              // first term
              // sum over \bar{omega}, \bar{upsilon} and \bar{rho}
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              // Zero initialze blocks accumulating sum over \bar{omega}
              std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omegapp;
              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omegapp);

              // Summing over \bar{omega}
              for (auto& omegapp_mult : omegapp_set) // omegapp_mult is multiplicity tagged \bar{omega}
                {
                   // A matrix will annihilate bra
                   if(baby_spncci_subspace_bra.Nn()==0) // if omega' is sigma', i.e., Nn'=N_omega'-N_sigma'=0
                     continue;
  
                   u3::U3 omegapp(omegapp_mult.irrep); // omegapp is \bar{omega}

                   if (not irrep_bra.ContainsSubspace(omegapp)) // if \bar{omega} is not present in Sp(3,R) irrep sigma'
                     continue;

                   // get hypersector index
                   spncci::BabySpNCCISubspaceLabels baby_spncci_labelspp(sigmap,Sp_bra,Sn_bra,S_bra,omegapp);
                   int baby_spncci_subspace_indexpp=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labelspp);

                   // \bar{omega} subspace (\bar{upsilon})
                   sp3r::U3Subspace u3_subspacepp=irrep_bra.LookUpSubspace(omegapp);
                   int upsilon_maxpp=u3_subspacepp.upsilon_max(); // upsilon_maxpp is maximal \bar{upsilon}
                   int dimpp=upsilon_maxpp*gamma_maxp;

                   Eigen::MatrixXd A 
                     = sp3r::Sp3rRaisingOperator(sigmap,u3_subspacep,u3_subspacepp,u_coef_cache);

                   // Zero initialze blocks accumulating sum over rho2, rho3 and rho0bp
                   std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rhobp;

                   ZeroInitBlocks(num_blocks,dimpp,dim1,unit_tensor_blocks_rhobp);

                   // Summing over rho0bp
                   int rho0bp_max=u3::OuterMultiplicity(omega1.SU3(),x0,omegapp.SU3());
                   for(int rho0bp=1; rho0bp<=rho0bp_max; ++rho0bp) // rho0bp is \bar{rho}
                     {
                       // Get hypersector index
                       int hypersector_index3
                         =baby_spncci_hypersectors.LookUpHypersectorIndex(baby_spncci_subspace_indexpp,baby_spncci_subspace_index1,unit_tensor_subspace_index,rho0bp);
                       if(hypersector_index3==-1)
                         continue;

                       double coef3=0;
                       for (int rho2=1; rho2<=rho0bp_max; rho2++)
                         {
               	           for (int rho3=1; rho3<=rho0bp_max; rho3++)
			     {
			        coef3+=u3::PhiCached(phi_coef_cache,omegapp.SU3(),u3::Conjugate(x0),omega1.SU3(),rho0bp,rho2)
			   	  *PhiCached(phi_coef_cache,omega1.SU3(),u3::Conjugate(omegapp.SU3()),u3::Conjugate(x0),rho2,rho3)
				  *u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(omegapp.SU3()),omega.SU3(),omega1.SU3(),
				  	       u3::SU3(2,0),1,1,u3::Conjugate(x0),rho3,rho0);
			     }
		         }

		       coef3=ParitySign(u3::ConjugationGrade(x0)+u3::ConjugationGrade(omega.SU3())
	                     +u3::ConjugationGrade(omegap.SU3()))
		   	     *sqrt(1.*u3::dim(omega1.SU3())*u3::dim(x0)*u3::dim(omegapp.SU3())/6)*coef3;

                       // sum over blocks for term 1
                       for(int b=0; b<num_blocks; ++b)
                         unit_tensor_blocks_rhobp[b]+=coef3*unit_tensor_hyperblocks[hypersector_index3][b];
                     }

                   // matrix product A*unit_tensor_block (upsilon',\bar{upsilon})*(\bar{upsilon},upsilon1)
                   // add in A operator and sum over blocks
                   for(int b=0; b<num_blocks; ++b)
                     for(int i=0; i<gamma_maxp; ++i)
                       for(int j=0; j<gamma_max; ++j)
                         {
                           // Get target indices
                           int it=i*upsilon_maxp;
                           int jt=j*upsilon_max1;
                           // Get source indices
                           int is=i*upsilon_maxpp;
                           int js=j*upsilon_max1;

                           unit_tensor_blocks_omegapp[b].block(it,jt,upsilon_maxp,upsilon_max1)
                             +=A*unit_tensor_blocks_rhobp[b].block(is,js,upsilon_maxpp,upsilon_max1);
                         }

		} //omegapp

              // accumulating sum over omegapp
              for(int b=0; b<num_blocks; ++b)
                  unit_tensor_blocks_omega1[b]=unit_tensor_blocks_omegapp[b];

              //summing over x0'
              for (auto& x0p_mult : x0p_set) // x0p_mult is multiplicity tagged Gamma0'
                {
                  u3::SU3 x0p(x0p_mult.irrep); // x0p is Gamma0'

                  int rho0p_max=OuterMultiplicity(omega1.SU3(),x0p,omegap.SU3()); // rho0p_max is maximal \bar{rho0}

                  // Zero initialize blocks accumlating sum over x0p
                  std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_x0p;
                  ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_x0p);

//                  for(int b=0; b<num_blocks; ++b)
//                    unit_tensor_blocks_x0p[b]+=unit_tensor_blocks_omegapp[b];

                  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  // second term
                  //////////////////////////////////////////////////////////////////////////////////////////////////////////
                  if(u3::OuterMultiplicity(u3::SU3(etap,0),u3::SU3(0,eta-2),x0p)>0)
                    {
                      assert((eta-2)>=0);
                      // look up index of subspace in unit tensor space
                      u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap,eta-2);

                      int unit_tensor_subspace_index1=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                      assert(unit_tensor_subspace_index1!=-1);

                      double coef1=(1.-(1./nucleon_number))*sqrt(1.*u3::dim(u3::SU3(eta,0)))*u3::dim(x0p)
                        *u3::UCached(u_coef_cache,u3::SU3(etap,0),u3::SU3(0,eta),x0p,u3::SU3(2,0),x0,1,1,u3::SU3(0,eta-2),1,1);

                      // zero initialize blocks for accumulating second term in sum over rhobp
                      std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                      ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                      // summing over rho0bp and accumulating sectors
                      for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp) // rho0bp is \bar{rho0}
                        {
                          int hypersector_index1=baby_spncci_hypersectors.LookUpHypersectorIndex(
                              baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                              unit_tensor_subspace_index1, rho0bp
                            );

                          // Accumulate
                          if(hypersector_index1==-1)
                            continue;

                          for(int b=0; b<num_blocks; ++b)
                            {
                               unit_tensor_blocks_rho0bp[b]
                                 +=u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(x0p),omega.SU3(),u3::SU3(2,0),
			                       omega1.SU3(),rho0bp,1,u3::Conjugate(x0),1,rho0)
                                   *unit_tensor_hyperblocks[hypersector_index1][b];
                            }
                          } //end rho0bp

                        for(int b=0; b<num_blocks; ++b)
                          unit_tensor_blocks_x0p[b]+=coef1*unit_tensor_blocks_rho0bp[b];
                      }

                      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // third term
                      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                      if ((u3::OuterMultiplicity(u3::SU3(etap+2,0),u3::SU3(0,eta),x0p)>0) && (etap+2)<=Nmax+std::max(N1vp,N1vn))
                        {
                          // look up index of subspace in unit tensor space
                          u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap+2,eta);

                          int unit_tensor_subspace_index2=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                          assert(unit_tensor_subspace_index2!=-1);

                          double coef2=-1*(1.+(1./nucleon_number))*ParitySign(u3::ConjugationGrade(x0)+u3::ConjugationGrade(x0p))
                                  *u3::dim(u3::SU3(etap,0))*u3::dim(x0p)
                                  *u3::UCached(u_coef_cache,u3::SU3(2,0),u3::SU3(etap,0),x0p,u3::SU3(0,eta),
		         		 u3::SU3(etap+2,0),1,1,x0,1,1)/sqrt(double(u3::dim(u3::SU3(etap+2,0))));

                          // zero initialize blocks for accumulating first term in sum over rhobp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                          ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                          for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp) // rho0bp is \bar{rho0}
                            {
                              int hypersector_index2=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                    baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                    unit_tensor_subspace_index2, rho0bp
                                  );

                              // Accumulate
                              if(hypersector_index2==-1)
                                continue;

                              for(int b=0; b<num_blocks; ++b)
                                {
                                  unit_tensor_blocks_rho0bp[b]
                                    +=u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(x0p),omega.SU3(),u3::SU3(2,0),
						  omega1.SU3(),rho0bp,1,u3::Conjugate(x0),1,rho0)
                                      *unit_tensor_hyperblocks[hypersector_index2][b];
                                }
                            } //end rho0bp

                          for(int b=0; b<num_blocks; ++b)
                            unit_tensor_blocks_x0p[b]+=coef2*unit_tensor_blocks_rho0bp[b];
                        }

                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // fourth term
                      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                      if((u3::OuterMultiplicity(u3::SU3(etap+1,0),u3::SU3(0,eta-1),x0p)>0) && (etap+1)<=Nmax+std::max(N1vp,N1vn))
                        {
                          assert((eta-1)>=0);
                          // look up index of subspace in unit tensor space
                          u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap+1,eta-1);

                          int unit_tensor_subspace_index3=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                          assert(unit_tensor_subspace_index3!=-1);

                          double coef4=0.0;
			  // x0 x (1,0) -> xpp
			  // (etap,0) x (0,eta-1) -> xpp
			  // (1,0) x xpp -> x0p
			  MultiplicityTagged<u3::SU3>::vector xpp_set=KroneckerProduct(x0,u3::SU3(1,0));
			  for (auto& xpp_mult : xpp_set) // xpp_mult is multiplicity tagged Gamma''
                            {
                              u3::SU3 xpp(xpp_mult.irrep); // xpp is Gamma''
			      if((u3::OuterMultiplicity(u3::SU3(etap,0),u3::SU3(0,eta-1),xpp)==0)
			         || (u3::OuterMultiplicity(u3::SU3(1,0),xpp,x0p)==0))
				continue;
                              coef4+=ParitySign(u3::ConjugationGrade(xpp))
				*u3::UCached(u_coef_cache,u3::SU3(etap,0),u3::SU3(0,eta),xpp,u3::SU3(1,0),x0,1,1,u3::SU3(0,eta-1),1,1)
				*u3::UCached(u_coef_cache,u3::SU3(1,0),u3::SU3(etap,0),x0p,u3::SU3(0,eta-1),u3::SU3(etap+1,0),1,1,xpp,1,1)
				*u3::UCached(u_coef_cache,x0,u3::SU3(1,0),x0p,u3::SU3(1,0),xpp,1,1,u3::SU3(2,0),1,1);
			    }
			  coef4=-1*ParitySign(u3::ConjugationGrade(x0p))*(etap+1)*u3::dim(x0p)
				  *sqrt(2.*(eta+2)/(etap+3))*coef4/nucleon_number;

                          // zero initialize blocks for accumulating second term in sum over rhobp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                            ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                          // summing over rho0bp and accumulating sectors
                          for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp) // rho0bp is \bar{rho0}
                            {
                              int hypersector_index4=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                    baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                    unit_tensor_subspace_index3, rho0bp
                                  );

                              // Accumulate
                              if(hypersector_index4==-1)
                                continue;

                              for(int b=0; b<num_blocks; ++b)
                              {
                                unit_tensor_blocks_rho0bp[b]
                                  +=u3::UCached(u_coef_cache,omegap.SU3(),u3::Conjugate(x0p),omega.SU3(),u3::SU3(2,0),
				     omega1.SU3(),rho0bp,1,u3::Conjugate(x0),1,rho0)
                                    *unit_tensor_hyperblocks[hypersector_index4][b];
                              }
                            } //end rho0bp

                          for(int b=0; b<num_blocks; ++b)
                            unit_tensor_blocks_x0p[b]+=coef4*unit_tensor_blocks_rho0bp[b];
                        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        for(int b=0; b<num_blocks; ++b)
                          unit_tensor_blocks_omega1[b]+=unit_tensor_blocks_x0p[b];

                  }// end x0p sum
		
                for(int b=0; b<num_blocks; ++b)
                        unit_tensor_blocks_omega1[b]*=coef;

                // summing over n, rho, n1, rho1, upsilon1
                for(int b=0; b<num_blocks; ++b)
                  for(int i=0; i<gamma_maxp; ++i)
                    for(int j=0; j<gamma_max; ++j)
                    {
                      int it=i*upsilon_maxp;
                      int jt=j*upsilon_max;
                      int is=i*upsilon_maxp;
                      int js=j*upsilon_max1;
		      
                      unit_tensor_blocks[b].block(it,jt,upsilon_maxp,upsilon_max)
                        +=unit_tensor_blocks_omega1[b].block(is,js,upsilon_maxp,upsilon_max1)*KBUK;

                    }
                // done with blocks
              }// end omega1_mult
        }// end hypersector index
    }// end Nsum
*/
  }

void DoTwoBodyRecurrenceInitialization(
  int Nmax, int N1vp, int N1vn,
  const spncci::LGIPair& lgi_pair,
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::TwoBodyDensitySpace& two_body_density_space,
  spncci::OperatorBlocks& lgi_transformations,
  bool transform_lgi_families,
  std::map<spncci::NnPair,std::set<int>>& tbd_subspace_subsets,
  spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersector_seeds,
  spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersector_seeds_conj,
  basis::OperatorHyperblocks<double>& tbd_hyperblocks_seeds,
  basis::OperatorHyperblocks<double>& tbd_hyperblocks_seeds_conj,
  u3::PhiCoefCache& phi_coef_cache
  )
{
    // Extract lgi index  labels
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
    ///////////////////////////////////////////////////////////////////////////
    // Read in list of TBDs between lgi pair and conjugates from files
    // Returned bool, files_found_test has no current use, but could be used
    //  to identify lgi pair with no non-zero rmes between them.
    // Corresponding rho values stored separately for later hypersector lookup
    ///////////////////////////////////////////////////////////////////////////
    // Initialize containers
    std::vector<u3shell::TwoBodyDensityLabels> lgi_tbds;
    std::vector<int> rho_values;

    // Get index corresponding to lgi in the full space.
    // Index may differ from lgi index in basis if space has been truncated
    int index1=lgi_full_space_index_lookup[irrep_family_index_bra];
    int index2=lgi_full_space_index_lookup[irrep_family_index_ket];

    // Read in operators
    std::string lgi_tbd_filename
      =fmt::format("seeds/tboperators_{:06d}_{:06d}.dat",index1,index2);
    bool files_found_test=lgi::ReadTwoBodyDensityLabels(lgi_tbd_filename,lgi_tbds,rho_values);

    ///////////////////////////////////////////////////////////////////////////
    // Set up hypersectors and hyperblocks for seeds
    // Generate hypersectors from list of TBDs and outer mulitplicites
    //  read from files
    // Populate hyperblocks using seeds read in from file
    //////////////////////////////////////////////////////////////////////////

    // Reads in TBD seed blocks and stores them in a vector of blocks. Order
    // corresponds to order of (TBD,rho) pairs in corresponding operator file.
    basis::OperatorBlocks<double> tbd_seed_blocks;
    std::string seed_filename
      =fmt::format("seeds/tbseeds_{:06d}_{:06d}.rmes",index1,index2);
    files_found_test&=lgi::ReadBlocks(seed_filename, lgi_tbds.size(), tbd_seed_blocks);

    // Identify TBD subspaces for recurrence
    spncci::GenerateRecurrenceTwoBodyDensities(
      Nmax,N1vp,N1vn,lgi_tbds,
      two_body_density_space,tbd_subspace_subsets
    );

    //Generate hypersectors for TBDs between lgi pair (lgi1,lgi2)
    baby_spncci_hypersector_seeds
      =spncci::BabySpNCCITwoBodyDensityHypersectors(
        lgi_families,baby_spncci_space,two_body_density_space,
        tbd_subspace_subsets[spncci::NnPair(0,0)],
        irrep_family_index_bra,irrep_family_index_ket
      );

    // Generate hypersectors for TBDs between conjugate lgi pair (lgi2,lgi1)
    baby_spncci_hypersector_seeds_conj
      =spncci::BabySpNCCITwoBodyDensityHypersectors(
        lgi_families,baby_spncci_space,two_body_density_space,
        tbd_subspace_subsets[spncci::NnPair(0,0)],
        irrep_family_index_ket,irrep_family_index_bra
      );

    // Zero initialize seed hyperblocks and conjugate hyperblocks
    basis::SetHyperoperatorToZero(baby_spncci_hypersector_seeds,tbd_hyperblocks_seeds);
    basis::SetHyperoperatorToZero(baby_spncci_hypersector_seeds_conj,tbd_hyperblocks_seeds_conj);

    // Populate the hyperblocks and conjugate hyperblocks with the seed
    // Conjugate hyperspectors will be used in calculating Nn0 rmes in recurrence
    spncci::PopulateHypersectorsWithSeedsForTwoBodyRecurrence(
      irrep_family_index_bra, irrep_family_index_ket,lgi_families,
      baby_spncci_space,two_body_density_space,
      baby_spncci_hypersector_seeds_conj,baby_spncci_hypersector_seeds,
      lgi_tbds,rho_values,tbd_seed_blocks,
      tbd_hyperblocks_seeds_conj,tbd_hyperblocks_seeds,
      phi_coef_cache
    );

}

  bool GenerateTwoBodyDensityHyperblocks(
      const spncci::LGIPair& lgi_pair,
      int Nmax, int N1vp, int N1vn, int nucleon_number,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::TwoBodyDensitySpace& tbd_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const spncci::KMatrixCache& kinv_matrix_cache,
      std::map<spncci::NnPair,std::set<int>>& tbd_subspace_subsets,
      const spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersector_seeds,
      const spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersector_seeds_conj,
      const basis::OperatorHyperblocks<double>& tbd_hyperblocks_seeds,
      const basis::OperatorHyperblocks<double>& tbd_hyperblocks_seeds_conj,
      u3::UCoefCache& u_coef_cache,
      u3::PhiCoefCache& phi_coef_cache,
      spncci::BabySpNCCITwoBodyDensityHypersectors& baby_spncci_hypersectors,
      basis::OperatorHyperblocks<double>& tbd_hyperblocks
      )
  {
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;

    // Generate Nn=0 hypersectors to be computed by conjugation
    bool Nn0_conjugate_hypersectors=true;
    std::vector<std::vector<int>> tbd_hypersector_subsets_Nn0;

    spncci::BabySpNCCITwoBodyDensityHypersectors baby_spncci_hypersectors_Nn0(
      Nmax, baby_spncci_space, tbd_space,
      tbd_subspace_subsets, tbd_hypersector_subsets_Nn0,
      irrep_family_index_ket, irrep_family_index_bra,
      Nn0_conjugate_hypersectors
    );

    // Hypersectors for Nnp>=Nn
    Nn0_conjugate_hypersectors=false;
    std::vector<std::vector<int>> tbd_hypersector_subsets;

    // (Nnp,Nn) sectors for Nnp>Nn
    // std::cout<<"generate Nnp>Nn hypersectors"<<std::endl;
    baby_spncci_hypersectors=spncci::BabySpNCCITwoBodyDensityHypersectors(
      Nmax,baby_spncci_space, tbd_space,
      tbd_subspace_subsets, tbd_hypersector_subsets,
      irrep_family_index_bra,irrep_family_index_ket,
      Nn0_conjugate_hypersectors
    );

    //Zero initialize hyperblocks for both Nn=0 and Nnp>=Nn sectors
    basis::OperatorHyperblocks<double> tbd_hyperblocks_Nn0;
    basis::SetHyperoperatorToZero(baby_spncci_hypersectors_Nn0,tbd_hyperblocks_Nn0);
    basis::SetHyperoperatorToZero(baby_spncci_hypersectors,tbd_hyperblocks);

    //Iterate over seed hypersectors and identify corresponding hypersector in hyperblocks or hyperblocks_Nn0
    //Transfer seeds into hyperblocks for recurrence
    for(int seed_hypersector_index=0; seed_hypersector_index<baby_spncci_hypersector_seeds.size(); ++seed_hypersector_index)
      {
        int bra_index,ket_index,operator_index,multiplicity_index;
        std::tie(bra_index,ket_index,operator_index,multiplicity_index)
          =baby_spncci_hypersector_seeds.GetHypersector(seed_hypersector_index).Key();
        int hypersector_index
          =baby_spncci_hypersectors.LookUpHypersectorIndex(bra_index,ket_index,operator_index,multiplicity_index);

        tbd_hyperblocks[hypersector_index]=tbd_hyperblocks_seeds[seed_hypersector_index];
      }

    //Transfer conjugate hypersectors into Nn0 hypersectors for recurrence
    for(int seed_hypersector_index=0; seed_hypersector_index<baby_spncci_hypersector_seeds_conj.size(); ++seed_hypersector_index)
      {
        int bra_index,ket_index,operator_index,multiplicity_index;
        std::tie(bra_index,ket_index,operator_index,multiplicity_index)
          =baby_spncci_hypersector_seeds_conj.GetHypersector(seed_hypersector_index).Key();
        int hypersector_index
          =baby_spncci_hypersectors_Nn0.LookUpHypersectorIndex(bra_index,ket_index,operator_index,multiplicity_index);

        tbd_hyperblocks_Nn0[hypersector_index]=tbd_hyperblocks_seeds_conj[seed_hypersector_index];
      }

    // std::cout<<"Compute Nn=0 blocks"<<std::endl;
    spncci::ComputeTwoBodyDensityHyperblocks(
      Nmax,N1vp,N1vn,nucleon_number,u_coef_cache,phi_coef_cache,
      k_matrix_cache,kinv_matrix_cache,spncci_space,baby_spncci_space,
      tbd_space,baby_spncci_hypersectors_Nn0,
      tbd_hypersector_subsets_Nn0,tbd_hyperblocks_Nn0
    );

    // std::cout<<"Add Nn0 blocks to hyperblocks"<<std::endl;
    spncci::AddNn0BlocksToTwoBodyDensityHyperblocks(
      baby_spncci_space,tbd_space,
      baby_spncci_hypersectors_Nn0,baby_spncci_hypersectors,
      tbd_hyperblocks_Nn0,tbd_hyperblocks,phi_coef_cache
    );

    // std::cout<<"Compute TBD hyperblocks"<<std::endl;
    spncci::ComputeTwoBodyDensityHyperblocks(
      Nmax,N1vp,N1vn,nucleon_number,u_coef_cache,phi_coef_cache,
      k_matrix_cache,kinv_matrix_cache,spncci_space,baby_spncci_space,
      tbd_space,baby_spncci_hypersectors,
      tbd_hypersector_subsets,tbd_hyperblocks
    );

    return true;
  }
//**************************************************************************************************************

} // End namespace




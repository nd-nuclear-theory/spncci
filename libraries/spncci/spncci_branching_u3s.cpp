/****************************************************************
  spncci_branching_u3s.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/spncci_branching_u3s.h"

#include <fstream>
#include <iostream>

#include "mcutils/parsing.h"
#include "cppformat/format.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // basis indexing in U3S scheme for spncci basis branching
  ////////////////////////////////////////////////////////////////

  SubspaceU3S::SubspaceU3S(const u3::U3S& omegaS, const SpNCCISpace& spncci_space)
  {
    labels_=omegaS;
    int state_index=0;
    int family_index=0;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {
        const int gamma_max = spncci_irrep_family.gamma_max();
        const sp3r::Sp3RSpace& irrep_space = spncci_irrep_family.Sp3RSpace();

        // if space contains omega x S
        if(irrep_space.ContainsSubspace(omegaS.U3()) && spncci_irrep_family.S()==omegaS.S())
        {
          //Construct subspace
          // dim is nu_max*gamma_max (size up subspace)
          // index is starting index in sector matrix
          int dim=gamma_max*irrep_space.LookUpSubspace(omegaS.U3()).size();
          PushStateLabels(StateLabelsType(family_index,spncci_irrep_family.sigma(),dim,state_index));

          // increment indices
          state_index+=dim;
          ++family_index;
        }
      }
    sector_size_=state_index;
  }

  std::string SubspaceU3S::Str() const
  {
    return omegaS().Str();
  }

  std::string StateU3S::Str() const
  {
    return fmt::format("[{} {}]",gamma(),sigma().Str());
  }


  SpaceU3S::SpaceU3S(SpNCCISpace& spncci_space)
  {
    for(const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      // for each SpNCCI irrep family
      {
        const int gamma_max = spncci_irrep_family.gamma_max();
        const sp3r::Sp3RSpace& irrep_space = spncci_irrep_family.Sp3RSpace();
        HalfInt S = spncci_irrep_family.S();

        // iterate through omega space
        for(int subspace_index=0; subspace_index<irrep_space.size(); ++subspace_index)
          // for each omega within SpNCCI irrep family
          {
            u3::U3 omega(irrep_space.GetSubspace(subspace_index).GetSubspaceLabels());

            // skip if omegaS already in space
            u3::U3S omegaS(omega,S);
            // std::cout<<"Subspace labels "<<omegaS.Str()<<std::endl;
            if (ContainsSubspace(omegaS))
              continue;

            // otherwise construct subspace 
            SubspaceU3S subspace(omegaS,spncci_space);
            // std::cout<<"    Subspace "<<subspace.Str()<<std::endl;
            PushSubspace(subspace);
          }
      }

    // get dimension and starting index of last subspace

    // (mac): why? for total dimension?  will restructure how store
    // subspace dimensions, not as part of labels

    const SubspaceType& subspace=GetSubspace(subspaces_.size()-1);
    int dim=std::get<2>(subspace.GetStateLabels(subspace.size()-1));
    int index=std::get<3>(subspace.GetStateLabels(subspace.size()-1));
    dimension_=dim+index;
  }

  std::string SectorLabelsU3S::Str() const
  {
    return fmt::format("( {} {}  {}{} {} : {} {}  {}", bra_index(),ket_index(), N0(), x0().Str(),S0(),kappa0(),L0(),rho0());
  }

  void GetSectorsU3S(
    const spncci::SpaceU3S& space, 
    const std::vector<u3shell::IndexedOperatorLabelsU3S>& relative_tensor_labels,
    std::vector<spncci::SectorLabelsU3S>& sector_vector
    )
  // TODO: replace with a normal sectors object
  {
    int u3s_sector_vector_index=0;
    for(auto tensor_labels:relative_tensor_labels)
      for(int j=0; j<space.size(); ++j)
        for(int i=0; i<=j; ++i)
          {
            int kappa0,L0; 
            u3shell::OperatorLabelsU3S op_labels;
            std::tie(op_labels, kappa0,L0)=tensor_labels;
            assert(kappa0!=0);
            u3::U3S omegapSp(space.GetSubspace(i).GetSubspaceLabels());
            u3::U3S omegaS(space.GetSubspace(j).GetSubspaceLabels());
            int rho0_max=u3::OuterMultiplicity(omegaS.SU3(), op_labels.x0(),omegapSp.SU3());
            // Check if allowed U(1) coupling
            if(omegaS.U3().N()+op_labels.N0()!=omegapSp.U3().N())
              continue;
            // Check if allowed SU(2) coupling 
            if(not am::AllowedTriangle(omegaS.S(), op_labels.S0(), omegapSp.S()))
              continue;
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                spncci::SectorLabelsU3S sector(i,j,op_labels, kappa0,L0,rho0);
                // should not alread by a sector in map
                // assert(not u3s_sectors.count(sector));
                // if(not u3s_sectors.count(sector))
                //   {
                //     u3s_sectors[sector]=u3s_sector_vector_index;
                //     ++u3s_sector_vector_index;
                //   }
                sector_vector.push_back(sector);
              }
          }
    // sector_vector.resize(u3s_sectors.size());
    // for(auto it=u3s_sectors.begin(); it!=u3s_sectors.end(); ++it)
    //   sector_vector[it->second]=it->first;
  }



}  // namespace

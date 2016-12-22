/****************************************************************
  sp_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include <fstream>
#include <iostream>

#include "mcutils/parsing.h"
#include "cppformat/format.h"
#include "spncci/sp_basis.h"


namespace spncci
{

  std::string SpIrrep::Str() const
  {
    std::ostringstream ss;

    ss << "[" 
       << " " << Nex()
       << " " << U3().Str()
       << " (" << Sp()
       << "," << Sn()
       << "," << S()
       << ") " << "]";
    return ss.str();
  }

  std::string SpIrrep::DebugString() const
  {
    std::ostringstream ss;

    ss << Str() << std::endl;
    ss << "  " 
      // << " Nn_max " << Nn_max  << " dimension " << dimension
       << " sp_space_ptr " << sp_space_ptr_
       << std::endl;

    return ss.str();
  }



  int NmaxTruncator::Nn_max(const u3::U3& sigma) const
  {
    HalfInt Nsigma = sigma.N();
    int value = Nmax_ - int(Nsigma - Nsigma_0_);
    return value;
  }

  void GenerateSp3RIrreps(
        const lgi::LGIVector& lgi_vector,
        const TruncatorInterface& truncator,
        SpIrrepVector& sp_irrep_vector,
        SigmaIrrepMap& sigma_irrep_map
        )
  {
    // populate SpIrrep vector with lgi's
    for(auto lgi_tag :lgi_vector)
        sp_irrep_vector.emplace_back(lgi_tag.irrep, lgi_tag.tag);

    // traverse SpIrrep vector
    for (int i=0; i<sp_irrep_vector.size(); ++i)
      {
        MultiplicityTagged<spncci::SpIrrep>& sp_irrep_tag=sp_irrep_vector[i];
        // retrieve sigma of LGI
        // u3::U3 sigma(it->sigma); // CRYPTIC
        spncci::SpIrrep& sp_irrep = sp_irrep_tag.irrep;
        u3::U3 sigma = sp_irrep.sigma();

        // find truncation
        //
        // Note: Even if a subspace is truncated into oblivion (Nn_max<0), we
        // still need to create an empty Sp3RSpace and store the
        // relevant information with the LGI, e.g., that its dimension
        // is zero.

        int Nn_max = truncator.Nn_max(sigma);

        // construct and add Sp3RSpace (if needed)

        // FAILS: since compiler anticipates the case where sigma
        // is not found as a key and map would (hypothetically)
        // need to insert an SP3RSpace using the (nonexistent)
        // default constructor sp3r::Sp3RSpace()
        //
        // sigma_irrep_map[sigma] = sp3r::Sp3RSpace(sigma,Nn_max);
        //
        // SOLUTION: add default constructor

        // FAILS in g++ 4.5: map.emplace not yet defined?
        //
        // sigma_irrep_map.emplace(
        //                sigma,  // key
        //                sigma,Nn_max  // constructor arguments for value
        //                );

        // WORKS: but cumbersome
        // sigma_irrep_map.insert(
        //                         std::make_pair(sigma,sp3r::Sp3RSpace(sigma,Nn_max))
        //                         );

        if (sigma_irrep_map.count(sigma) == 0)
          sigma_irrep_map[sigma] = sp3r::Sp3RSpace(sigma,Nn_max);
          
        // save info back to SpIrrep
        const sp3r::Sp3RSpace& irrep = sigma_irrep_map[sigma];
        sp_irrep.SaveSubspaceInfo(irrep);
      }
  }

  ////////////////////////////////////////////////////////////////
  // iteration over Sp(3,R) -> U(3) subspaces & states
  ////////////////////////////////////////////////////////////////

  int TotalU3Subspaces(const SpIrrepVector& sp_irrep_vector)
  {
    int subspaces = 0;
    for (auto sp_irrep_tag : sp_irrep_vector)
      {
        const spncci::SpIrrep& sp_irrep = sp_irrep_tag.irrep;
        const sp3r::Sp3RSpace& irrep = sp_irrep.Sp3RSpace();
        subspaces += sp_irrep_tag.tag*irrep.size();
      }

    return subspaces;
  }

  int TotalDimensionU3(const SpIrrepVector& sp_irrep_vector)
  {
    int dimension = 0;
    for (auto sp_irrep_tag : sp_irrep_vector)
      {
        const spncci::SpIrrep& sp_irrep = sp_irrep_tag.irrep;
        const sp3r::Sp3RSpace& irrep = sp_irrep.Sp3RSpace();
        for (int i_subspace = 0; i_subspace < irrep.size(); ++i_subspace)
          {
            const sp3r::U3Subspace& u3_subspace = irrep.GetSubspace(i_subspace);
            dimension += sp_irrep_tag.tag*u3_subspace.size();
          }
      }

    return dimension;
  }

  int TotalDimensionU3LS(const SpIrrepVector& sp_irrep_vector)
  {
    int dimension = 0;
    for (auto sp_irrep_tag : sp_irrep_vector)
      {
        const spncci::SpIrrep& sp_irrep = sp_irrep_tag.irrep;
        const sp3r::Sp3RSpace& irrep = sp_irrep.Sp3RSpace();
        for (int i_subspace=0; i_subspace < irrep.size(); ++i_subspace)
          {
            const sp3r::U3Subspace& u3_subspace = irrep.GetSubspace(i_subspace);

            // L branching of each irrep in subspace
            //
            // branching_vector is list of L values tagged with their
            // multiplicity.
            u3::U3 omega = u3_subspace.U3();
            MultiplicityTagged<int>::vector branching_vector = BranchingSO3(omega.SU3());
            int L_dimension = 0;
            for (int i_L=0; i_L<branching_vector.size(); ++i_L)
              L_dimension += branching_vector[i_L].tag;

            // dimension contribution of subspace
            //
            // At LS-coupled level, the dimension contribution is the
            // product of the number of U(3) reduced states and their
            // L substates.
            dimension += sp_irrep_tag.tag*u3_subspace.size()*L_dimension;
          }
      }

    return dimension;
  }

  int TotalDimensionU3LSJConstrained(const SpIrrepVector& sp_irrep_vector, const HalfInt& J)
  {
    int dimension = 0;
    for (auto sp_irrep_tag : sp_irrep_vector)
      {
        const spncci::SpIrrep& sp_irrep = sp_irrep_tag.irrep;
        HalfInt S = sp_irrep.S();
        const sp3r::Sp3RSpace& irrep = sp_irrep.Sp3RSpace();
        for (int i_subspace=0; i_subspace < irrep.size(); ++i_subspace)
          {
            const sp3r::U3Subspace& u3_subspace = irrep.GetSubspace(i_subspace);

            // L branching of each irrep in subspace
            //   subject to triangle(LSJ) constraint
            //
            // branching_vector is list of L values tagged with their
            // multiplicity.
            u3::U3 omega = u3_subspace.U3();
            MultiplicityTagged<int>::vector branching_vector = BranchingSO3Constrained(omega.SU3(),am::ProductAngularMomentumRange(S,J));
            // std::cout << " omega " << omega.Str() << " S " << S << " J " << J << std::endl;
            int L_dimension = 0;
            for (int i_L=0; i_L<branching_vector.size(); ++i_L)
              // std::cout << " L " << branching_vector[i_L].irrep << " " << branching_vector[i_L].tag << std::endl;
              L_dimension += branching_vector[i_L].tag;

            // dimension contribution of subspace
            //
            // At LS-coupled level, the dimension contribution is the
            // product of the number of U(3) reduced states and their
            // L substates.
            dimension += sp_irrep_tag.tag*u3_subspace.size()*L_dimension;

          }
      }
    return dimension;
  }

  int TotalDimensionU3LSJAll(const SpIrrepVector& sp_irrep_vector)
  {
    int dimension = 0;
    for (auto sp_irrep_tag : sp_irrep_vector)
      {
        const spncci::SpIrrep& sp_irrep = sp_irrep_tag.irrep;
        HalfInt S = sp_irrep.S();
        const sp3r::Sp3RSpace& irrep = sp_irrep.Sp3RSpace();
        for (int i_subspace=0; i_subspace < irrep.size(); ++i_subspace)
          {
            const sp3r::U3Subspace& u3_subspace = irrep.GetSubspace(i_subspace);

            // L branching of each irrep in subspace
            //
            // branching_vector is list of L values tagged with their
            // multiplicity.
            u3::U3 omega = u3_subspace.U3();
            MultiplicityTagged<int>::vector branching_vector = BranchingSO3(omega.SU3());
            for (int i_L=0; i_L<branching_vector.size(); ++i_L)
              {
                int L = branching_vector[i_L].irrep;
                int L_multiplicity = branching_vector[i_L].tag;
                int J_values = int((L+S) - abs(L-S) + 1);

                // dimension contribution of this L
                //
                // At J-coupled level, the dimension contribution of this L is the
                // product of the L multiplicity with the number of J values.
                dimension += sp_irrep_tag.tag*L_multiplicity*J_values;
              }
          }
      }

    return dimension;
  }

  std::vector< std::pair<int,int> > 
  GenerateSpIrrepPairs(spncci::SpIrrepVector sp_irrep_vector )
  {
    std::vector< std::pair<int,int> >  sp_irrep_pair_vector;
    for(int i=0; i<sp_irrep_vector.size(); ++i)
      for(int j=0; j<=i; ++j)
        {
          spncci::SpIrrep sp_irrep1=sp_irrep_vector[i].irrep;
          spncci::SpIrrep sp_irrep2=sp_irrep_vector[j].irrep;
          if (abs(sp_irrep1.Sp()-sp_irrep2.Sp())<=2)
            if (abs(sp_irrep1.Sn()-sp_irrep2.Sn())<=2)
              if (abs(sp_irrep1.S()-sp_irrep2.S())<=2)
                 sp_irrep_pair_vector.push_back(std::pair<int,int>(i,j));
        }
    return sp_irrep_pair_vector;
  }

  SubspaceU3S::SubspaceU3S(const u3::U3S& omegaS,const SpIrrepVector& sp_irrep_vector)
  {
    labels_=omegaS;
    int index=0;
    // iterate over U(3)xSU(2) irreps
    for(int i=0; i< sp_irrep_vector.size(); ++i)
      {
        auto irrep=sp_irrep_vector[i].irrep;
        // extract U(3) spaces labeled by omega
        const sp3r::Sp3RSpace& space=irrep.Sp3RSpace();
        // if space contains omega and S
        if(space.ContainsSubspace(omegaS.U3()) && irrep.S()==omegaS.S())
        {
          //Construct subspace
          //dim is nu_max*gamma_max (size up subspace)
          // index is starting index in sector matrix
          int dim=sp_irrep_vector[i].tag*space.LookUpSubspace(omegaS.U3()).size();
          PushStateLabels(StateLabelsType(i,irrep.sigma(),dim,index));
          // increment index 
          index+=dim;
        }
      }
    sector_size_=index;
  }

  std::string SubspaceU3S::Str() const
  {
    return omegaS().Str();
  }

  std::string StateU3S::Str() const
  {
    return fmt::format("[{} {}]",gamma(),sigma().Str());
  }


  SpaceU3S::SpaceU3S(SpIrrepVector& sp_irrep_vector)
  {
    // iterate over U(3)xSU(2) irreps
    for(auto irrep_tag : sp_irrep_vector)
      {
        const SpIrrep& irrep=irrep_tag.irrep;
        HalfInt S(irrep.S());
        // iterate through omega space
        const sp3r::Sp3RSpace& space=irrep.Sp3RSpace();        
        for(int i=0; i<space.size(); ++i)
          {
            u3::U3 omega(space.GetSubspace(i).GetSubspaceLabels());
            // if omegaS already in space, continue to next subspace
            u3::U3S omegaS(omega,S);
            // std::cout<<"Subspace labels "<<omegaS.Str()<<std::endl;
            if(lookup_.count(omegaS))
              continue;
            // otherwise construct subspace 
            SubspaceU3S subspace(omegaS,sp_irrep_vector);
            // std::cout<<"    Subspace "<<subspace.Str()<<std::endl;
            PushSubspace(subspace);
          }
      }
    // get dimension and starting index of last subspace
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
  {
    int u3s_sector_vector_index=0;
    SectorLabelsU3SCache u3s_sectors;
    for(auto tensor_labels:relative_tensor_labels)
      for(int j=0; j<space.size(); ++j)
        for(int i=0; i<=j; ++i)
          {
            int kappa0,L0; 
            u3shell::OperatorLabelsU3S op_labels;
            std::tie(op_labels, kappa0,L0)=tensor_labels;
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
                if(not u3s_sectors.count(sector))
                  {
                    u3s_sectors[sector]=u3s_sector_vector_index;
                    ++u3s_sector_vector_index;
                  }
              }
          }
    sector_vector.resize(u3s_sectors.size());
    for(auto it=u3s_sectors.begin(); it!=u3s_sectors.end(); ++it)
      sector_vector[it->second]=it->first;
  }

  SubspaceLS::SubspaceLS(const int& L, const HalfInt& S,const SpaceU3S& u3s_space)
  {
    labels_=std::tuple<int,HalfInt>(L,S);
    int index=0;
    // iterate over U(3)xSU(2) irreps
    for(int i=0; i<u3s_space.size(); ++i)
      {
        SubspaceU3S subspace=u3s_space.GetSubspace(i);
        u3::U3 omega(subspace.omega());
        int kappa_max=u3::BranchingMultiplicitySO3(omega.SU3(),L);
        // if space contains S and omega can branch to L
        if(kappa_max>0 && subspace.S()==S)
        {
          //Construct subspace
          int dim=subspace.sector_dim();
          PushStateLabels(StateLabelsType(omega,kappa_max,index));
          // increment index 
          index+=kappa_max*dim;
        }
      }

    sector_size_=index;
  }

  std::string SubspaceLS::Str() const
  {
    return fmt::format("[{} {}]",L(),S());
  }

  std::string StateLS::Str() const
  {
    return fmt::format("[{} {}]",omega().Str(),kappa_max());
  }


  SpaceLS::SpaceLS(const SpaceU3S& u3s_space, HalfInt J)
  {
    // iterate over U(3)xSU(2) irreps
    for(int i=0; i<u3s_space.size(); ++i)
    // for(auto u3s_subspace : u3s_space)
      {
        const SubspaceU3S& u3s_subspace=u3s_space.GetSubspace(i);
        HalfInt S(u3s_subspace.S());
        // iterate through omega space
        u3::U3 omega(u3s_subspace.omega());
        // interate over possible L values
        for(int L=int(abs(S-J)); L<=(S+J); ++L)
          {
            if(lookup_.count(std::tuple<int,HalfInt>(L,S)))
              continue;
            SubspaceLS ls_subspace(L,S,u3s_space);
            PushSubspace(ls_subspace);
          }
       }
    // get dimension and starting index of last subspace
    const SubspaceType& subspace=GetSubspace(subspaces_.size()-1);
    HalfInt S=subspace.S();
    u3::U3 omega=std::get<0>(subspace.GetStateLabels(subspace.size()-1));
    int index=std::get<2>(subspace.GetStateLabels(subspace.size()-1));
    u3::U3S omegaS(omega,S);
    int dim=u3s_space.LookUpSubspace(omegaS).sector_dim();
    dimension_=dim+index;
  }

  // SectorsU3S::SectorsU3S(SpaceU3S& space)
  // {
  //   for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
  //     for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
  //       {
  //         // retrieve subspaces
  //         const SubspaceU3S& bra_subspace = space.GetSubspace(bra_subspace_index);
  //         const SubspaceU3S& ket_subspace = space.GetSubspace(ket_subspace_index);

  //         // push sectors (taking unit multiplicity)
  //         int multiplicity_index = 1;
  //         PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace,multiplicity_index));         
  //       }
  // }

  // SectorsU3S::SectorsU3S(SpaceU3S& space, const u3shell::OperatorLabelsU3ST& operator_labels)
  // {
  //   for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
  //     for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
  //       {
  //         // retrieve subspaces
  //         const SubspaceU3S& bra_subspace = space.GetSubspace(bra_subspace_index);
  //         const SubspaceU3S& ket_subspace = space.GetSubspace(ket_subspace_index);
  //         // verify selection rules
  //         bool allowed = true;
  //         // U(1)
  //         allowed &= (ket_subspace.N() + operator_labels.N0() - bra_subspace.N() == 0);
  //         // spin & isosopin
  //         allowed &= am::AllowedTriangle(ket_subspace.S(),operator_labels.S0(),bra_subspace.S());
  //         // find SU(3) multiplicity and check SU(3) selection
  //         int multiplicity = 0;
  //         if (allowed)
  //           {
  //             multiplicity = u3::OuterMultiplicity(ket_subspace.omega().SU3(),operator_labels.x0(),bra_subspace.omega().SU3());
  //             allowed &= (multiplicity > 0);
  //           }

  //         // push sectors (tagged by multiplicity)
  //         if (allowed)
  //           for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
  //             {
  //               PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace,multiplicity_index));
  //             }
  //       }
  // }



}  // namespace

/****************************************************************
  spncci_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT
****************************************************************/

#include "spncci/spncci_basis.h"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <set>

#include "mcutils/parsing.h"
#include "fmt/format.h"
#include "am/halfint_fmt.h"

#include "sp3rlib/vcs.h"

namespace spncci
{
  int ValenceShellForNuclide(const lgi::NuclideType& nuclide)
  {
    // each major shell eta=2*n+l (for a spin-1/2 fermion) contains (eta+1)*(eta+2) substates
    int N1v = 0;
    for (int species_index : {0,1})
      {
        int num_particles = nuclide[species_index];
        for (int eta=0; num_particles>0; ++eta)
          {
            // add contribution from particles in shell
            int shell_degeneracy = (eta+1)*(eta+2);
            int num_particles_in_shell = std::min(num_particles,shell_degeneracy);

            // register this shell as occupied
            N1v = std::max(N1v,eta);

            // discard particles in shell
            num_particles -= num_particles_in_shell;
          }
      }

    return N1v;
  }


  std::string SpNCCIIrrepFamily::Str() const
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

  std::string SpNCCIIrrepFamily::DebugStr() const
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


  int NlimitTruncator::Nn_max(const u3::U3& sigma) const
  {
    int value;
    HalfInt Nsigma = sigma.N();
    if(Nsigma==Nsigma_0_)
      value=Nmax_;
    else
      value = std::max(0, Nlimit_ - int(Nsigma - Nsigma_0_));
      // std::cout<<Nsigma<<"  "<<(Nlimit_ - int(Nsigma - Nsigma_0_))<<"  "<<value<<std::endl;
    return value;
  }



  void ConstructRestrictedSp3RSpace(const u3::U3& sigma, int Nn_max, sp3r::Sp3RSpace& irrep)
    // Contruct an Sp(3,R) irrep applying necessary restrictions for A<6
  // Sp3RSpace::Sp3RSpace(const u3::U3& sigma, int Nn_max, const sp3r::RestrictedSpanakopitaType& spanakopita)
    {
      // create irrep with no restrictions
      sp3r::Sp3RSpace irrep_temp(sigma,Nn_max);

      // Generate K matrices for irrep
      vcs::MatrixCache K_matrix_cache;
      vcs::MatrixCache Kinv_matrix_cache;
      vcs::GenerateKMatrices(irrep_temp,K_matrix_cache, Kinv_matrix_cache);

      // set up container for states
      std::map<MultiplicityTagged<u3::U3>,MultiplicityTagged<u3::U3>::vector> spanakopita;
      // Iterate through K matrix identifying omega subspaces which are contained in irrep along with upsilon_max
      for(auto it=K_matrix_cache.begin(); it!=K_matrix_cache.end(); ++it)
        {
          const u3::U3& omega=it->first;
          auto& Kmatrix=it->second;
          int upsilon_max=Kmatrix.rows();
          MultiplicityTagged<u3::U3> omega_tagged(omega,upsilon_max);
          const auto& subspace=irrep_temp.GetSubspace(irrep_temp.LookUpSubspaceIndex(omega));
          MultiplicityTagged<u3::U3>::vector& states=spanakopita[omega_tagged];
          for(int i=0; i<subspace.size(); ++i)
            states.push_back(subspace.GetStateLabels(i));
        }

      // Create restricted irrep
      irrep=sp3r::Sp3RSpace(sigma,Nn_max,spanakopita);
    }

  void GenerateSpNCCISpace(
      const lgi::MultiplicityTaggedLGIVector& multiplicity_tagged_lgi_vector,
      const TruncatorInterface& truncator,
      SpNCCISpace& spncci_space,
      SigmaIrrepMap& sigma_irrep_map,
      bool restrict_sp3r_to_u3_branching
    )
  {

    // populate SpNCCI space with irrep families if gamma_max>0
    // int num_lgi=multiplicity_tagged_lgi_vector.size();
    // for(int lgi_index=0; lgi_index<multiplicity_tagged_lgi_vector.size(); ++lgi_index);
    //   {
    //     const auto& lgi_tag=multiplicity_tagged_lgi_vector[lgi_index];
    //   }
    for(auto& lgi_tag : multiplicity_tagged_lgi_vector)
    {
      // if(lgi_tag.tag>0)
        spncci_space.emplace_back(lgi_tag.irrep, lgi_tag.tag);
    }

    // build up the irrep for each irrep family
    for (spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {

        // if(spncci_irrep_family.gamma_max()==0)
        //   continue;
        // find truncation
        //
        // Note: Even if a subspace is truncated into oblivion
        // (Nn_max<0), we still need to create an empty Sp3RSpace and
        // store the relevant information with the SpNCCIIrrepFamily,
        // e.g., that its dimension is zero.

        const u3::U3& sigma = spncci_irrep_family.sigma();
        int Nn_max = truncator.Nn_max(sigma);

        // construct and add Sp3RSpace (if needed)

        // access by [] indexing
        //
        // FAILS: since compiler anticipates the case where sigma
        // is not found as a key and map would (hypothetically)
        // need to insert an SP3RSpace using the (nonexistent)
        // default constructor sp3r::Sp3RSpace()
        //
        // SOLUTION: add default constructor
        //
        // sigma_irrep_map[sigma] = sp3r::Sp3RSpace(sigma,Nn_max);

        // access by map.emplace
        //
        // FAILS in g++ 4.5: map.emplace not yet defined?
        //
        // FAILS in g++ 5.3.0: weird missing allocator errors
        //
        // sigma_irrep_map.emplace(
        //                sigma,  // key
        //                sigma,Nn_max  // constructor arguments for value
        //                );

        // access by pair insertion
        //
        // WORKS: but cumbersome
        // sigma_irrep_map.insert(
        //                         std::make_pair(sigma,sp3r::Sp3RSpace(sigma,Nn_max))
        //                         );


        if (sigma_irrep_map.count(sigma) == 0)
          {
            if(restrict_sp3r_to_u3_branching)
              {
                sp3r::Sp3RSpace irrep_restricted;

                ConstructRestrictedSp3RSpace(sigma,Nn_max,irrep_restricted);
                sigma_irrep_map[sigma]=irrep_restricted;
              }
            else
              sigma_irrep_map[sigma] = sp3r::Sp3RSpace(sigma,Nn_max,restrict_sp3r_to_u3_branching);
          }
        // save info back to SpNCCIIrrepFamily
        const sp3r::Sp3RSpace& irrep_space = sigma_irrep_map[sigma];
        spncci_irrep_family.SaveSpaceReference(irrep_space);
      }
  }

  ////////////////////////////////////////////////////////////////
  // iteration over Sp(3,R) -> U(3) subspaces & states
  ////////////////////////////////////////////////////////////////

  int TotalU3Subspaces(const SpNCCISpace& spncci_space)
  {
    int subspaces = 0;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {
        const int gamma_max = spncci_irrep_family.gamma_max();
        const sp3r::Sp3RSpace& irrep_space = spncci_irrep_family.Sp3RSpace();
        subspaces += gamma_max*irrep_space.size();
      }

    return subspaces;
  }

  int TotalDimensionU3S(const SpNCCISpace& spncci_space)
  {
    int dimension = 0;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {
        const int gamma_max = spncci_irrep_family.gamma_max();
        const sp3r::Sp3RSpace& irrep_space = spncci_irrep_family.Sp3RSpace();
        for (int subspace_index = 0; subspace_index < irrep_space.size(); ++subspace_index)
          {
            const sp3r::U3Subspace& u3_subspace = irrep_space.GetSubspace(subspace_index);
            dimension += gamma_max*u3_subspace.size();
          }
      }

    return dimension;
  }

  int TotalDimensionU3LS(const SpNCCISpace& spncci_space)
  {
    int dimension = 0;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {
        const int gamma_max = spncci_irrep_family.gamma_max();
        const sp3r::Sp3RSpace& irrep_space = spncci_irrep_family.Sp3RSpace();
        for (int subspace_index=0; subspace_index < irrep_space.size(); ++subspace_index)
          {
            const sp3r::U3Subspace& u3_subspace = irrep_space.GetSubspace(subspace_index);

            // L branching of each irrep in subspace
            //
            // branching_vector is list of L values tagged with their
            // multiplicity.
            u3::U3 omega = u3_subspace.U3();
            MultiplicityTagged<int>::vector branching_vector = BranchingSO3(omega.SU3());
            int L_dimension = 0;
            for (int L_index=0; L_index<branching_vector.size(); ++L_index)
              L_dimension += branching_vector[L_index].tag;

            // dimension contribution of subspace
            //
            // At LS-coupled level, the dimension contribution is the
            // product of the number of U(3) reduced states and their
            // L substates.
            dimension += gamma_max*u3_subspace.size()*L_dimension;
          }
      }

    return dimension;
  }

  int TotalDimensionU3LSJConstrained(const SpNCCISpace& spncci_space, const HalfInt& J)
  {
    int dimension = 0;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {
        const int gamma_max = spncci_irrep_family.gamma_max();
        const sp3r::Sp3RSpace& irrep_space = spncci_irrep_family.Sp3RSpace();
        HalfInt S = spncci_irrep_family.S();
        for (int subspace_index=0; subspace_index < irrep_space.size(); ++subspace_index)
          {
            const sp3r::U3Subspace& u3_subspace = irrep_space.GetSubspace(subspace_index);

            // L branching of each irrep in subspace
            //   subject to triangle(LSJ) constraint
            //
            // branching_vector is list of L values tagged with their
            // multiplicity.
            u3::U3 omega = u3_subspace.U3();
            MultiplicityTagged<int>::vector branching_vector = BranchingSO3Constrained(omega.SU3(),am::ProductAngularMomentumRange(S,J));
            // std::cout << " omega " << omega.Str() << " S " << S << " J " << J << std::endl;
            int L_dimension = 0;
            for (int L_index=0; L_index<branching_vector.size(); ++L_index)
              // std::cout << " L " << branching_vector[L_index].irrep << " " << branching_vector[L_index].tag << std::endl;
              L_dimension += branching_vector[L_index].tag;

            // dimension contribution of subspace
            //
            // At LS-coupled level, the dimension contribution is the
            // product of the number of U(3) reduced states and their
            // L substates.
            dimension += gamma_max*u3_subspace.size()*L_dimension;

          }
      }
    return dimension;
  }

  int TotalDimensionU3LSJAll(const SpNCCISpace& spncci_space)
  {
    int dimension = 0;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {
        const int gamma_max = spncci_irrep_family.gamma_max();
        const sp3r::Sp3RSpace& irrep_space = spncci_irrep_family.Sp3RSpace();
        HalfInt S = spncci_irrep_family.S();
        for (int subspace_index=0; subspace_index < irrep_space.size(); ++subspace_index)
          {
            const sp3r::U3Subspace& u3_subspace = irrep_space.GetSubspace(subspace_index);

            // L branching of each irrep in subspace
            //
            // branching_vector is list of L values tagged with their
            // multiplicity.
            u3::U3 omega = u3_subspace.U3();
            MultiplicityTagged<int>::vector branching_vector = BranchingSO3(omega.SU3());
            for (int L_index=0; L_index<branching_vector.size(); ++L_index)
              {
                int L = branching_vector[L_index].irrep;
                int L_multiplicity = branching_vector[L_index].tag;
                int J_values = int((L+S) - abs(L-S) + 1);

                // dimension contribution of this L
                //
                // At J-coupled level, the dimension contribution of this L is the
                // product of the L multiplicity with the number of J values.
                dimension += gamma_max*L_multiplicity*J_values;
              }
          }
      }

    return dimension;
  }

  ////////////////////////////////////////////////////////////////
  // SpNCCI irrep family pair enumeration
  ////////////////////////////////////////////////////////////////

  std::vector< std::pair<int,int> >
  GenerateSpNCCIIrrepFamilyPairs(spncci::SpNCCISpace spncci_space)
  {
    std::vector< std::pair<int,int> >  sp_irrep_family_pair_vector;
    for(int i=0; i<spncci_space.size(); ++i)
      for(int j=0; j<=i; ++j)
        {
          const spncci::SpNCCIIrrepFamily& sp_irrep_family1 = spncci_space[i];
          const spncci::SpNCCIIrrepFamily& sp_irrep_family2 = spncci_space[j];
          if (abs(sp_irrep_family1.Sp()-sp_irrep_family2.Sp())<=2)
            if (abs(sp_irrep_family1.Sn()-sp_irrep_family2.Sn())<=2)
              if (abs(sp_irrep_family1.S()-sp_irrep_family2.S())<=2)
                sp_irrep_family_pair_vector.push_back(std::pair<int,int>(i,j));
        }
    return sp_irrep_family_pair_vector;
  }

  ////////////////////////////////////////////////////////////////
  // BabySpNCCI
  ////////////////////////////////////////////////////////////////

  BabySpNCCISubspace::BabySpNCCISubspace(
      const spncci::SpNCCIIrrepFamily& spncci_irrep_family,
      int irrep_family_index,
      const sp3r::U3Subspace& u3_subspace
    )
  {

    // set labels
    labels_ = BabySpNCCISubspaceLabels(
        spncci_irrep_family.sigma(),
        spncci_irrep_family.Sp(),
        spncci_irrep_family.Sn(),
        spncci_irrep_family.S(),
        u3_subspace.U3()  // omega
      );

    // save dimension info
    gamma_max_ = spncci_irrep_family.gamma_max();
    upsilon_max_ = u3_subspace.size();
    upsilon_max_ = u3_subspace.upsilon_max();
    irrep_family_index_=irrep_family_index;
    dimension_ = gamma_max_*upsilon_max_;
  }

  std::string BabySpNCCISubspace::LabelStr() const
  {
    return fmt::format("{} ({} {} {}) {}",sigma().Str(),Sp(),Sn(),S(),omega().Str());
  }

  std::string BabySpNCCISubspace::DebugStr() const
  {
    return fmt::format("{}: gamma_max {} upsilon_max {} -> dim {}",LabelStr(),gamma_max(),upsilon_max(),size());
  }

  BabySpNCCISpace::BabySpNCCISpace(const spncci::SpNCCISpace& spncci_space)
  {
    // traverse irrep families
    for(int irrep_family_index=0; irrep_family_index<spncci_space.size(); ++irrep_family_index)
      {
        const spncci::SpNCCIIrrepFamily& spncci_irrep_family=spncci_space[irrep_family_index];
        // extract Sp(3,R) space
        const sp3r::Sp3RSpace& sp_space = spncci_irrep_family.Sp3RSpace();

        // traverse U(3) subspaces
        for (int subspace_index = 0; subspace_index < sp_space.size(); ++subspace_index)
          {
            const sp3r::U3Subspace& u3_subspace = sp_space.GetSubspace(subspace_index);
            PushSubspace(BabySpNCCISubspace(spncci_irrep_family,irrep_family_index,u3_subspace));
          }
      }
  }

  BabySpNCCISectors::BabySpNCCISectors(
      const spncci::BabySpNCCISpace& space,
      const u3shell::OperatorLabelsU3S& operator_labels
    )
    : BaseSectors(space)
  // Based loosely on u3shell::SectorsU3SPN constructor.
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const BabySpNCCISubspace& bra_subspace = space.GetSubspace(bra_subspace_index);
          const BabySpNCCISubspace& ket_subspace = space.GetSubspace(ket_subspace_index);

          // verify selection rules
          bool allowed = true;
          // U(1)
          allowed &= (ket_subspace.omega().N() + operator_labels.N0() - bra_subspace.omega().N() == 0);
          // spin
          //
          // Note: Basic two-body constaints can be placed on Sp
          // and Sn triangularity based on two-body nature of
          // operator, so delta Sp<=2 and delta Sn<=2.  However, in
          // general, the operator does not have sharp Sp0 or Sn0.
          allowed &= am::AllowedTriangle(ket_subspace.S(),operator_labels.S0(),bra_subspace.S());
          allowed &= abs(int(ket_subspace.Sp()-bra_subspace.Sp()))<=2;
          allowed &= abs(int(ket_subspace.Sn()-bra_subspace.Sn()))<=2;
          if (!allowed)
            continue;

          // find SU(3) multiplicity and check SU(3) selection
          int multiplicity = u3::OuterMultiplicity(
              ket_subspace.omega().SU3(),
              operator_labels.x0(),
              bra_subspace.omega().SU3()
            );
          allowed &= (multiplicity > 0);

          // push sectors (tagged by multiplicity)
          if (allowed)
            for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
              {
                PushSector(bra_subspace_index,ket_subspace_index,multiplicity_index);
              }
        }
  }


  BabySpNCCIHypersectors::BabySpNCCIHypersectors(
        const spncci::BabySpNCCISpace& space,
        const u3shell::RelativeUnitTensorSpaceU3S& operator_space
      )
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const BabySpNCCISubspace& bra_subspace = space.GetSubspace(bra_subspace_index);
          const BabySpNCCISubspace& ket_subspace = space.GetSubspace(ket_subspace_index);

          for(int operator_subspace_index=0; operator_subspace_index<operator_space.size(); ++operator_subspace_index)
            {
              // verify selection rules
              bool allowed = true;
              const u3shell::RelativeUnitTensorSubspaceU3S&
                operator_subspace=operator_space.GetSubspace(operator_subspace_index);

              // U(1)
              allowed &= (ket_subspace.omega().N() + operator_subspace.N0() - bra_subspace.omega().N() == 0);
              // spin
              //
              // Note: Basic two-body constaints can be placed on Sp
              // and Sn triangularity based on two-body nature of
              // operator, so (delta Sp)<=2 and (delta Sn)<=2.  However, in
              // general, the operator does not have sharp Sp0 or Sn0.
              allowed &= am::AllowedTriangle(ket_subspace.S(),operator_subspace.S0(),bra_subspace.S());
              allowed &= abs(int(ket_subspace.Sp()-bra_subspace.Sp()))<=2;
              allowed &= abs(int(ket_subspace.Sn()-bra_subspace.Sn()))<=2;
              if (!allowed)
                continue;

              // find SU(3) multiplicity and check SU(3) selection
              int multiplicity = u3::OuterMultiplicity(
                  ket_subspace.omega().SU3(),
                  operator_subspace.x0(),
                  bra_subspace.omega().SU3()
                );
              allowed &= (multiplicity > 0);

              // push sectors (tagged by multiplicity)
              if (allowed)
                for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
                  {
                    PushHypersector(HypersectorType(
                      bra_subspace_index,ket_subspace_index,operator_subspace_index,
                      bra_subspace,ket_subspace,operator_subspace,multiplicity_index));
                  }
            }
        }

  }



    BabySpNCCIHypersectors::BabySpNCCIHypersectors(
      int Nmax,
    const spncci::BabySpNCCISpace& space,
    const u3shell::RelativeUnitTensorSpaceU3S& operator_space,
    const std::map<spncci::NnPair,std::set<int>>& operator_subsets_NnpNn,
    std::vector<std::vector<int>>& unit_tensor_hypersector_subsets,
    int irrep_family_index_1, int irrep_family_index_2, bool Nn0_conjugate_hypersectors
  )
  {
    // Baby spncci hyperspector constructor for recurrence hypersectors
    //
    // If Nn0_conjugate_hypersectors=true, then only construct Nnp=0 and Nn!=0 sectors
    // else construct only the lower triangle hypersectors, excluding Nn=0.
    //
    // hypersectors are restricted based on angular momentum adddition and SU(3) coupling

    unit_tensor_hypersector_subsets.resize(Nmax+1);

    int hypersector_index=0;
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const BabySpNCCISubspace& bra_subspace = space.GetSubspace(bra_subspace_index);
          const BabySpNCCISubspace& ket_subspace = space.GetSubspace(ket_subspace_index);

          int Nnp=bra_subspace.Nn();
          int Nn=ket_subspace.Nn();

          // std::cout<<Nnp<<"  "<<Nn<<std::endl;

          // If all we want is the Nnp=0 and Nn!=0 conjugated sectors
          if(Nn0_conjugate_hypersectors)
            {
              if(Nnp!=0)
                continue;
            }
          // Otherwise, only take sectors with Nnp>=Nn
          else
          {
            if(Nn>Nnp)
              continue;
          }

          bool in_irrep_families=(
            (bra_subspace.irrep_family_index()== irrep_family_index_1)
            &&(ket_subspace.irrep_family_index()== irrep_family_index_2)
            );


          if(not in_irrep_families)
            continue;

          int Nsum=Nnp+Nn;
          // std::cout<<"Nsum "<<Nsum<<std::endl;
          const std::set<int>& operator_subset=operator_subsets_NnpNn.at(spncci::NnPair(Nnp,Nn));
          // For each operator subspace, check if its an allowed operator subspace determined
          // by SU(2) and U(3) constraints.  If allowed, push multiplicity tagged hypersectors
          for(int operator_subspace_index : operator_subset)
            {
              bool allowed_subspace = true;
              const u3shell::RelativeUnitTensorSubspaceU3S&
                operator_subspace=operator_space.GetSubspace(operator_subspace_index);

              // U(1)
              // std::cout<<"N0 "<<operator_subspace.N0()<<std::endl;
              allowed_subspace &= (ket_subspace.omega().N() + operator_subspace.N0() - bra_subspace.omega().N() == 0);
              // spin
              //
              // Note: Basic two-body constaints can be placed on Sp
              // and Sn triangularity based on two-body nature of
              // operator, so (delta Sp)<=2 and (delta Sn)<=2.  However, in
              // general, the operator does not have sharp Sp0 or Sn0.
              allowed_subspace &= am::AllowedTriangle(ket_subspace.S(),operator_subspace.S0(),bra_subspace.S());
              allowed_subspace &= abs(int(ket_subspace.Sp()-bra_subspace.Sp()))<=2;
              allowed_subspace &= abs(int(ket_subspace.Sn()-bra_subspace.Sn()))<=2;
              if (!allowed_subspace)
                continue;

              // find SU(3) multiplicity and check SU(3) selection
              int multiplicity = u3::OuterMultiplicity(
                  ket_subspace.omega().SU3(),operator_subspace.x0(),
                  bra_subspace.omega().SU3()
                );

              // push sectors (tagged by multiplicity)
              // std::cout<<"multiplicity "<<multiplicity<<std::endl;
              for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
                {
                  PushHypersector(
                    HypersectorType(
                      bra_subspace_index,ket_subspace_index,operator_subspace_index,
                      bra_subspace, ket_subspace,operator_subspace,
                      multiplicity_index
                      )
                    );
                  unit_tensor_hypersector_subsets[Nsum/2].push_back(hypersector_index);
                  ++hypersector_index;
                }
            }
        }
  }


    BabySpNCCIHypersectors::BabySpNCCIHypersectors(
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const spncci::BabySpNCCISpace& space,
    const u3shell::RelativeUnitTensorSpaceU3S& operator_space,
    const std::set<int>& operator_subset,
    // std::vector<int>& unit_tensor_hypersector_subset,
    int irrep_family_index_1, int irrep_family_index_2
  )
  {
    // Baby spncci hyperspector constructor for LGI only (Nn=Nnp=0)
    //  needed for creating seed hypersectors 
    //
    // hypersectors are restricted based on angular momentum adddition and SU(3) coupling
    // unit_tensor_hypersector_subsets.resize(Nmax+1);
    //TODO: redo so that there is a single look up of the appropriate baby spncci index
    
    // Get LGI labels 
    u3::U3 sigma1, sigma2;
    HalfInt Sp1,Sn1,S1,Sp2,Sn2,S2;
    
    const lgi::LGI& lgi_1=lgi_families[irrep_family_index_1].irrep;
    std::tie(std::ignore,sigma1,Sp1,Sn1,S1)=lgi_1.Key();
    

    const lgi::LGI& lgi_2=lgi_families[irrep_family_index_2].irrep;
    std::tie(std::ignore,sigma2,Sp2,Sn2,S2)=lgi_2.Key();
    
    // Get baby spncci subspace indices for lgi
    int bra_subspace_index=space.LookUpSubspaceIndex(BabySpNCCISubspaceLabels(sigma1,Sp1,Sn1,S1,sigma1));
    int ket_subspace_index=space.LookUpSubspaceIndex(BabySpNCCISubspaceLabels(sigma2,Sp2,Sn2,S2,sigma2));

    // retrieve subspaces for lgi
    const BabySpNCCISubspace& bra_subspace = space.GetSubspace(bra_subspace_index);
    const BabySpNCCISubspace& ket_subspace = space.GetSubspace(ket_subspace_index);

    //Initialize hypersector_index
    int hypersector_index=0;

    // For each operator subspace, check if its an allowed operator subspace determined
    // by SU(2) and U(3) constraints.  If allowed, push multiplicity tagged hypersectors
    for(int operator_subspace_index : operator_subset)
      {
        bool allowed_subspace = true;
        const u3shell::RelativeUnitTensorSubspaceU3S&
          operator_subspace=operator_space.GetSubspace(operator_subspace_index);

        // U(1) constraint
        allowed_subspace &= (ket_subspace.omega().N() + operator_subspace.N0() - bra_subspace.omega().N() == 0);
        
        // spin constraints
        //
        // Note: Basic two-body constaints can be placed on Sp
        // and Sn triangularity based on two-body nature of
        // operator, so (delta Sp)<=2 and (delta Sn)<=2.  However, in
        // general, the operator does not have sharp Sp0 or Sn0.
        allowed_subspace &= am::AllowedTriangle(ket_subspace.S(),operator_subspace.S0(),bra_subspace.S());
        allowed_subspace &= abs(int(ket_subspace.Sp()-bra_subspace.Sp()))<=2;
        allowed_subspace &= abs(int(ket_subspace.Sn()-bra_subspace.Sn()))<=2;
        if (!allowed_subspace)
          continue;

        // find SU(3) multiplicity and check SU(3) selection
        int multiplicity = u3::OuterMultiplicity(
            ket_subspace.omega().SU3(),operator_subspace.x0(),
            bra_subspace.omega().SU3()
          );

        // push sectors (tagged by multiplicity)
        for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
          {
            PushHypersector(
              HypersectorType(
                bra_subspace_index,ket_subspace_index,operator_subspace_index,
                bra_subspace, ket_subspace,operator_subspace,
                multiplicity_index
                )
              );
            // unit_tensor_hypersector_subset.push_back(hypersector_index);
            // ++hypersector_index;
          }
      }
  }



  void PrintHypersectors(
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
    )
  {
    for(int hypersector_index=0; hypersector_index<baby_spncci_hypersectors.size(); ++hypersector_index)
      {
        const auto& hypersector=baby_spncci_hypersectors.GetHypersector(hypersector_index);

        int unit_tensor_subspace_index, ket_subspace_index,bra_subspace_index, rho0;
        std::tie(bra_subspace_index, ket_subspace_index,unit_tensor_subspace_index,rho0)=hypersector.Key();

        const auto& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
        const auto& bra_subspace=baby_spncci_space.GetSubspace(bra_subspace_index);
        const auto& ket_subspace=baby_spncci_space.GetSubspace(ket_subspace_index);

        std::cout<<"hypersector "<<hypersector_index<<" "<< bra_subspace.LabelStr()<<"  "<<ket_subspace.LabelStr()
        <<"  "<<unit_tensor_subspace.LabelStr()<<rho0<<std::endl;
        for(int i=0; i<unit_tensor_subspace.size(); ++i)
        {
          int T0,Sp,Tp,S,T;
          std::tie(T0,Sp,Tp,S,T)=unit_tensor_subspace.GetStateLabels(i);
          std::cout<<fmt::format("{}  {} {}  {} {}",T0,Sp,Tp,S,T)<<std::endl;
          std::cout<<unit_tensor_hyperblocks[hypersector_index][i]<<std::endl<<std::endl;
        }
      }

  }


  ObservableBabySpNCCIHypersectors::ObservableBabySpNCCIHypersectors(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    int irrep_family_index_bra, int irrep_family_index_ket, bool restrict_lower_triangle
  )
  {
    // iterate over baby spncci subspaces
    for (int bra_subspace_index=0; bra_subspace_index<baby_spncci_space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<baby_spncci_space.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const BabySpNCCISubspace& bra_subspace = baby_spncci_space.GetSubspace(bra_subspace_index);
          const BabySpNCCISubspace& ket_subspace = baby_spncci_space.GetSubspace(ket_subspace_index);

          // Check if baby spncci subspaces contained in either the bra or the ket irrep families
          bool allowed_subspace=false;
          if(restrict_lower_triangle)
            {
              if(irrep_family_index_bra==bra_subspace.irrep_family_index()
                  &&irrep_family_index_ket==ket_subspace.irrep_family_index()
                )
                {
                  allowed_subspace=true;
                }

              allowed_subspace&=(irrep_family_index_bra>=irrep_family_index_ket);
              
            }
          else
            {
              if(irrep_family_index_bra==-1 && irrep_family_index_bra==-1 )
                allowed_subspace=true;

              else if(irrep_family_index_bra==bra_subspace.irrep_family_index()
                      &&irrep_family_index_ket==ket_subspace.irrep_family_index()
                      )
                allowed_subspace=true;

              else if(irrep_family_index_bra==ket_subspace.irrep_family_index()
                      && irrep_family_index_ket==bra_subspace.irrep_family_index()
                    )
                allowed_subspace=true;
            }
          // // If restrict_lower_triangle true, then only continue if irrep family of bra>=ket
          // if(restrict_lower_triangle && (irrep_family_index_bra<irrep_family_index_ket))
          //   allowed_subspace=false;

          // if not in bra or ket, continue;
          if(not allowed_subspace)
            continue;

          // For each observable subspace, check if its an allowed observable subspace determined
          // by SU(2) and U(3) constraints.  If allowed, push multiplicity tagged hypersectors
          for(int observable_subspace_index=0; observable_subspace_index<observable_space.size(); ++observable_subspace_index)
            {
              bool allowed_subspace = true;
              const u3shell::ObservableSubspaceU3S&
                observable_subspace=observable_space.GetSubspace(observable_subspace_index);

              // U(1)
              allowed_subspace
                &=(ket_subspace.omega().N()+observable_subspace.N0()-bra_subspace.omega().N() == 0);

              // spin
              //
              // Note: Basic two-body constaints can be placed on Sp
              // and Sn triangularity based on two-body nature of
              // observable, so (delta Sp)<=2 and (delta Sn)<=2.  However, in
              // general, the observable does not have sharp Sp0 or Sn0.
              allowed_subspace &= am::AllowedTriangle(ket_subspace.S(),observable_subspace.S0(),bra_subspace.S());
              allowed_subspace &= abs(int(ket_subspace.Sp()-bra_subspace.Sp()))<=2;
              allowed_subspace &= abs(int(ket_subspace.Sn()-bra_subspace.Sn()))<=2;
              if (!allowed_subspace)
                continue;

              // find SU(3) multiplicity and check SU(3) selection
              int multiplicity = u3::OuterMultiplicity(
                  ket_subspace.omega().SU3(),observable_subspace.x0(),
                  bra_subspace.omega().SU3()
                );

              // push sectors (tagged by multiplicity)
              for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
                {
                  PushHypersector(
                    HypersectorType(
                      bra_subspace_index,ket_subspace_index,observable_subspace_index,
                      bra_subspace, ket_subspace,observable_subspace,
                      multiplicity_index
                      )
                    );
                }
            }
        }
  }

  void NumBabySpncciSubspacesInIrrepFamily(
    int num_irrep_families,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    std::vector<int>& irrep_family_num_subspaces
  )
  {
    irrep_family_num_subspaces.resize(num_irrep_families);
    for (int subspace_index=0; subspace_index<baby_spncci_space.size(); ++subspace_index)
        {
          // retrieve subspaces
          const BabySpNCCISubspace& subspace = baby_spncci_space.GetSubspace(subspace_index);
          int irrep_family_index=subspace.irrep_family_index();
          int upsilon_max=subspace.upsilon_max();
          int gamma_max=subspace.gamma_max();
          irrep_family_num_subspaces[irrep_family_index]+=upsilon_max*gamma_max;
        }

  }

  void SortSubspacesDecending(const std::vector<int>& num_subspaces,std::vector<int>&ordered_subspaces)
  {
    ordered_subspaces.resize(num_subspaces.size());
    std::vector<int> copy_num_subspaces=num_subspaces;
    std::sort(copy_num_subspaces.begin(), copy_num_subspaces.end());
    int last_index=num_subspaces.size()-1;
    for(int j=0; j<num_subspaces.size(); ++j)
      {
        int num=num_subspaces[last_index-j];
        for(int i=0; i<num_subspaces.size(); ++i)
            if(num_subspaces[i]==num)
              ordered_subspaces[j]=i;
      }
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Branched spncci basis 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  SubspaceSpBasis::SubspaceSpBasis(
    const HalfInt& J,
    const u3shell::U3SPN& sigmaSPN,
    int irrep_family_index,
    const BabySpNCCISpace& baby_spncci_space
  )
  {

    // save labels
    // labels_ = sigmaSPN;
    labels_ = irrep_family_index;
    sigmaSPN_=sigmaSPN;
    dimension_=0;
    // scan BabySpNCCISpace for states to accumulate
    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace = baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        if(baby_spncci_subspace.irrep_family_index()!=irrep_family_index)
          continue;

        // push state
        const u3::U3& omega=baby_spncci_subspace.omega();
        const HalfInt& S=baby_spncci_subspace.S();
        int gamma_max=baby_spncci_subspace.gamma_max();
        int upsilon_max=baby_spncci_subspace.upsilon_max();

        //Should be the same for all states with same irrep_family_index
        gamma_max_=gamma_max;

        MultiplicityTagged<int>::vector so3_irreps=u3::BranchingSO3(omega.SU3());
        for(auto& tagged_L : so3_irreps)
          {
            int L=tagged_L.irrep;
            int kappa_max=tagged_L.tag;
            if(am::AllowedTriangle(L,S,J))
              {
                int degeneracy=gamma_max*upsilon_max*kappa_max;
                omegaLLabels state_labels(omega,L);
                PushStateLabels(state_labels,degeneracy);

                // record auxiliary state information
                state_kappa_max_.push_back(kappa_max);
                state_upsilon_max_.push_back(upsilon_max);
                state_baby_spncci_subspace_index_.push_back(baby_spncci_subspace_index);

                // increment dimension of subspace 
                dimension_+=degeneracy;
              }
          }
      }
  }

  std::string SubspaceSpBasis::LabelStr() const
  {
    return sigmaSPN().Str();
  }

  std::string SubspaceSpBasis::DebugStr() const
  {
    std::ostringstream os;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        const StateSpBasis state(*this,state_index);

        os << fmt::format(
            "  index {} sigmaSPN {} omega,L {}, {} multiplicity {} offset {}",
            state_index,
            state.sigmaSPN().Str(),state.omega().Str(),state.L(),
            state.degeneracy(),state.offset()
          ) << std::endl;
      }

    return os.str();
  }

  SpaceSpBasis::SpaceSpBasis(const BabySpNCCISpace& baby_spncci_space, const HalfInt& J)
  {
    J_=J;

    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        // create new subspace -- only if not already constructed for this (omega,S)
        const u3shell::U3SPN& sigmaSPN = baby_spncci_subspace.sigmaSPN();
        int irrep_family_index = baby_spncci_subspace.irrep_family_index();

        if(ContainsSubspace(irrep_family_index))
          continue;

        PushSubspace(SubspaceSpBasis(J,sigmaSPN,irrep_family_index,baby_spncci_space));
      }
  }


  SpaceSpBasis::SpaceSpBasis(const BabySpNCCISpace& baby_spncci_space, const HalfInt& J, std::set<int>irrep_family_subset)
  {
    J_=J;

    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {
  
        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_subspace_index);
  
        // create new subspace -- only if not already constructed for this (omega,S)
        const u3shell::U3SPN& sigmaSPN = baby_spncci_subspace.sigmaSPN();
        int irrep_family_index = baby_spncci_subspace.irrep_family_index();
        
        if(not irrep_family_subset.count(irrep_family_index))
          continue;

        if(ContainsSubspace(irrep_family_index))
          continue;

        PushSubspace(SubspaceSpBasis(J,sigmaSPN,irrep_family_index,baby_spncci_space));
      }
  }






  std::string SpaceSpBasis::DebugStr(bool show_subspaces) const
  {
    std::ostringstream os;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        // set up alias
        const SubspaceType& subspace = GetSubspace(subspace_index);

        os << fmt::format(
            "subspace_index {} labels {} size {} full_dimension {}",
            subspace_index,subspace.LabelStr(),subspace.size(),subspace.full_dimension()
          ) << std::endl;
        if (show_subspaces)
          os << subspace.DebugStr();

      }

    return os.str();
  }


  void GetSpBasisOffsets(
    const spncci::SpaceSpBasis& spbasis,
    std::vector<std::vector<int>>& offsets
    )
  {
    // Iterate through J-branched basis and identify starting position
    // of each irrep family in the basis listing.  Starting index stored
    // in offsets indexed by subspace index. 
    offsets.resize(spbasis.size());
    int offset=0;
    for(int subspace_index=0; subspace_index<spbasis.size(); ++subspace_index)
      {
        const spncci::SubspaceSpBasis& subspace=spbasis.GetSubspace(subspace_index);
        offsets[subspace_index].resize(subspace.size());
        const std::vector<int>& multiplicities=subspace.state_multiplicities();
        for(int state_index=0; state_index<subspace.size(); ++state_index)
          {
            offsets[subspace_index][state_index]=offset;
            int multiplicity=multiplicities[state_index];
            offset+=multiplicity;
          }
      }
    assert(spbasis.FullDimension()==offset);
  }

}  // namespace

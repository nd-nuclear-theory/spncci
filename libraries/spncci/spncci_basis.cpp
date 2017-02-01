/****************************************************************
  spncci_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/spncci_basis.h"

#include <fstream>
#include <iostream>

#include "mcutils/parsing.h"
#include "cppformat/format.h"

namespace spncci
{

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

  std::string SpNCCIIrrepFamily::DebugString() const
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

  void GenerateSpNCCISpace(
      const lgi::MultiplicityTaggedLGIVector& multiplicity_tagged_lgi_vector,
      const TruncatorInterface& truncator,
      SpNCCISpace& spncci_space,
      SigmaIrrepMap& sigma_irrep_map
    )
  {

    // populate SpNCCI space with irrep families
    for(auto& lgi_tag : multiplicity_tagged_lgi_vector)
      spncci_space.emplace_back(lgi_tag.irrep, lgi_tag.tag);

    // build up the irrep for each irrep family
    for (spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {

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
          sigma_irrep_map[sigma] = sp3r::Sp3RSpace(sigma,Nn_max);
          
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

}  // namespace

/****************************************************************
  sp_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include <fstream>
#include <iostream>

#include "mcutils/parsing.h"

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
    for(auto lgi :lgi_vector)

        sp_irrep_vector.emplace_back(lgi);

    // traverse SpIrrep vector
    for (auto it = sp_irrep_vector.begin(); it != sp_irrep_vector.end(); ++it)
      {
        // retrieve sigma of LGI
        // u3::U3 sigma(it->sigma); // CRYPTIC
        spncci::SpIrrep& sp_irrep = *it;
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
    for (auto it=sp_irrep_vector.begin(); it !=sp_irrep_vector.end(); ++it)
      {
        const spncci::SpIrrep& sp_irrep = *it;
        const sp3r::Sp3RSpace& irrep = sp_irrep.Sp3RSpace();
        subspaces += irrep.size();
      }

    return subspaces;
  }

  int TotalDimensionU3(const SpIrrepVector& sp_irrep_vector)
  {
    int dimension = 0;
    for (auto it=sp_irrep_vector.begin(); it !=sp_irrep_vector.end(); ++it)
      {
        const spncci::SpIrrep& sp_irrep = *it;
        const sp3r::Sp3RSpace& irrep = sp_irrep.Sp3RSpace();
        for (int i_subspace = 0; i_subspace < irrep.size(); ++i_subspace)
          {
            const sp3r::U3Subspace& u3_subspace = irrep.GetSubspace(i_subspace);
            dimension += u3_subspace.size();
          }
      }

    return dimension;
  }

  int TotalDimensionU3LS(const SpIrrepVector& sp_irrep_vector)
  {
    int dimension = 0;
    for (auto it=sp_irrep_vector.begin(); it !=sp_irrep_vector.end(); ++it)
      {
        const spncci::SpIrrep& sp_irrep = *it;
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
            dimension += u3_subspace.size()*L_dimension;
          }
      }

    return dimension;
  }

  int TotalDimensionU3LSJConstrained(const SpIrrepVector& sp_irrep_vector, const HalfInt& J)
  {
    int dimension = 0;
    for (auto it=sp_irrep_vector.begin(); it !=sp_irrep_vector.end(); ++it)
      {
        const spncci::SpIrrep& sp_irrep = *it;
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
            dimension += u3_subspace.size()*L_dimension;

          }
      }
    return dimension;
  }

  int TotalDimensionU3LSJAll(const SpIrrepVector& sp_irrep_vector)
  {
    int dimension = 0;
    for (auto it=sp_irrep_vector.begin(); it !=sp_irrep_vector.end(); ++it)
      {
        const spncci::SpIrrep& sp_irrep = *it;
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
                dimension += L_multiplicity*J_values;
              }
          }
      }

    return dimension;
  }

  std::vector< std::pair<int,int> > GenerateSpIrrepPairs(spncci::SpIrrepVector sp_irrep_vector )
  {
    std::vector< std::pair<int,int> >  sp_irrep_pair_vector;
    for(int i=0; i<sp_irrep_vector.size(); ++i)
      for(int j=0; j<sp_irrep_vector.size(); ++j)
      //for (int j=0; j<=i; j++)
       // need both sigma', sigma and sigma,sigma' pair so we don't have to compute conjugate unit tensors  
        {
          spncci::SpIrrep sp_irrep1=sp_irrep_vector[i];
          spncci::SpIrrep sp_irrep2=sp_irrep_vector[j];
          if (abs(sp_irrep1.Sp()-sp_irrep2.Sp())<=2)
            if (abs(sp_irrep1.Sn()-sp_irrep2.Sn())<=2)
              if (abs(sp_irrep1.S()-sp_irrep2.S())<=2)
                 sp_irrep_pair_vector.push_back(std::pair<int,int>(i,j));
        }

    return sp_irrep_pair_vector;
  }



}  // namespace

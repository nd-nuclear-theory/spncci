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
    irrep_family_index_=irrep_family_index;
    dimension_ = gamma_max_*upsilon_max_;
  }

  std::string BabySpNCCISubspace::LabelStr() const
  {
    return fmt::format("{} ({} {} {}) {}",sigma().Str(),Sp(),Sn(),S(),omega().Str());
  }

  std::string BabySpNCCISubspace::DebugStr() const
  {
    return fmt::format("{}: gamma_max {} upsilon_max {} -> dim {}",LabelStr(),gamma_max_,upsilon_max_,size());
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
                PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace,multiplicity_index));
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
    const spncci::BabySpNCCISpace& space,
    const u3shell::RelativeUnitTensorSpaceU3S& operator_space,
    std::map< spncci::NnPair, std::set<int>>& operator_subsets,
    std::vector<std::vector<int>>& unit_tensor_hypersector_subsets,
    int irrep_family_index_1, int irrep_family_index_2
  )
  {
    std::cout<<"irrep family1 "<<irrep_family_index_1<<"  irrep family2 "<<irrep_family_index_2<<std::endl;
    int hypersector_index=0;
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {

          // retrieve subspaces
          const BabySpNCCISubspace& bra_subspace = space.GetSubspace(bra_subspace_index);
          const BabySpNCCISubspace& ket_subspace = space.GetSubspace(ket_subspace_index);

          int Nnp=bra_subspace.Nn();
          int Nn=ket_subspace.Nn();

          if(Nnp>Nn)
            continue;


          bool allowed=(
            (ket_subspace.irrep_family_index()== irrep_family_index_1)
            &&(bra_subspace.irrep_family_index()== irrep_family_index_2)
            );
          allowed+=(
            (ket_subspace.irrep_family_index()== irrep_family_index_2)
            &&(bra_subspace.irrep_family_index()== irrep_family_index_1)
            );
          // Only want Nnp<=Nn, else continue;
          // Can have irrep_family_index1=bra or ket but not both
          // Can have irrep_family_index2=bra or ket but not both
          std::cout<<"hi "<<Nnp<<"  "<<Nn<<"  "<<bra_subspace.irrep_family_index()
            <<"  "<<ket_subspace.irrep_family_index()<<"  "<<allowed<<std::endl;

          // // Want only sectors with Nn>=Nnp, so if Nnp greater, we want the conjugate hypsersector
          // bool Nnp_greater=(Nnp>Nn);
          // // We only want conjugate hypersectors if irrep_family_1 and irrep_family_2 are different,
          // //  otherwise, we only want to keep the Nn>=Nnp hypersectors sectors 
          // Nnp_greater&=(irrep_family_index_1!=irrep_family_index_2);
          // if (Nnp_greater)
          //   std::cout<<"conjugate "<<bra_subspace.LabelStr()<<"  "<<ket_subspace.LabelStr()<<std::endl;
          // Key for operator subsets by NnpNn, if Nnp>Nn, want conjugate 
          spncci::NnPair NnpNn(Nnp,Nn);
          // NnpNn=Nnp_greater?spncci::NnPair(Nn,Nnp):spncci::NnPair(Nnp,Nn);
          int Nsum=Nnp+Nn;
          // if(Nnp_greater)
          //   std::cout<<"Nnp greater 1 "<<Nnp_greater<<std::endl;
          // bool allowed_hypersector=true;

          // allowed_hypersector&=Nnp_greater?
          //   (ket_subspace.irrep_family_index()== irrep_family_index_1)
          //   :(bra_subspace.irrep_family_index()== irrep_family_index_1);


          // if(Nnp_greater)
          //   std::cout<<"hi "<<ket_subspace.irrep_family_index()<<"  "<<irrep_family_index_1<<std::endl;

          // allowed_hypersector&=Nnp_greater?
          //   (bra_subspace.irrep_family_index()== irrep_family_index_2)
          //   :(ket_subspace.irrep_family_index()== irrep_family_index_2);

          // if(Nnp_greater)
          //   std::cout<<"ho  "<<allowed_hypersector<<std::endl;

          // if(not allowed_hypersector)                  
          if(not allowed)                  
            continue;

          // // if(Nnp_greater)
          // //   std::cout<<"hihi"<<std::endl;

          // int subspace_index_1=Nnp_greater?ket_subspace_index:bra_subspace_index;
          // const BabySpNCCISubspace& subspace_1=Nnp_greater?ket_subspace:bra_subspace;

          // if(Nnp_greater)
          //   std::cout<<"hoho"<<std::endl;

          // int subspace_index_2=Nnp_greater?bra_subspace_index:ket_subspace_index;
          // const BabySpNCCISubspace& subspace_2=Nnp_greater?bra_subspace:ket_subspace;          

          // if(Nnp_greater)
          //   std::cout<<"hihihi"<<std::endl;

          const std::set<int>& operator_subset=operator_subsets[NnpNn];

          // std::cout<<"Nnp greater "<<Nnp_greater<<std::endl;

          // if (Nnp_greater)
          //   std::cout<<"conjugate 2 "<<subspace_1.LabelStr()<<"  "<<subspace_2.LabelStr()<<std::endl;


          // std::cout<<"Nsum "<<Nsum<<" Nnp "<<Nnp<<" Nn "<<Nn<<" "<<operator_subset.size()<<std::endl;
          for(int operator_subspace_index : operator_subset)
            {
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
                  ket_subspace.omega().SU3(),operator_subspace.x0(),
                  bra_subspace.omega().SU3()
                );

              allowed &= (multiplicity > 0);
              // std::cout<<"multiplicity "<<multiplicity<<std::endl;
              // push sectors (tagged by multiplicity)
              if (allowed)
                for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
                  {
                    PushHypersector(
                      HypersectorType(
                        bra_subspace_index,ket_subspace_index,operator_subspace_index,
                        bra_subspace, ket_subspace,operator_subspace,
                        multiplicity_index
                        )
                      );
                    // std::cout<<subspace_1.LabelStr()<<" "<<subspace_2.LabelStr()<<" "<<operator_subspace.LabelStr()<<std::endl;
                    unit_tensor_hypersector_subsets[Nsum/2].push_back(hypersector_index);
                    ++hypersector_index;
                  }
            }
        }
  }



  ////////////////////////////////////////////////////////////////
  // precomputation of K matrices
  ////////////////////////////////////////////////////////////////

  void
  PrecomputeKMatrices(
      const spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::KMatrixCache& k_matrix_cache,
      bool intrinsic
    )
  {
    for (const auto& sigma_irrep_pair : sigma_irrep_map)
      {
        // extract sigma and irrep contents
        const u3::U3& sigma = sigma_irrep_pair.first;
        const sp3r::Sp3RSpace& sp_irrep = sigma_irrep_pair.second;

        // populate K matrix cache for this irrep
        vcs::GenerateKMatrices(sp_irrep,k_matrix_cache[sigma],intrinsic);
      }

  }

}  // namespace

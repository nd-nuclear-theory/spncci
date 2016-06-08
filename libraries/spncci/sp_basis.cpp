/****************************************************************
  sp_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include <fstream>
#include <iostream>

#include "utilities/parsing.h"

#include "spncci/sp_basis.h"


namespace spncci
{

  std::string LGI::Str() const
  {
    std::ostringstream ss;

    ss << "[" 
       << " " << Nex
       << " " << sigma.Str()
       << " (" << Sp
       << "," << Sn
       << "," << S
       << ") " << "]";
    return ss.str();
  }

  std::string LGI::DebugString() const
  {
    std::ostringstream ss;

    ss << Str() << std::endl;
    ss << "  " 
      // << " Nn_max " << Nn_max  << " dimension " << dimension
       << " irrep_ptr " << irrep_ptr
       << std::endl;

    return ss.str();
  }

  void GenerateLGIVector(LGIVectorType& lgi_vector, const std::string& lgi_filename, const HalfInt& Nsigma_0)
  {

    // open input file
    std::ifstream lgi_stream(lgi_filename.c_str());
    OpenCheck(bool(lgi_stream),lgi_filename);  // explicit cast to bool required in gcc 5

    // scan input file
    std::string line;
    int line_count = 0;
    while ( std::getline(lgi_stream, line) )
      {
        // count line
        ++line_count;

        // set up for parsing
        std::istringstream line_stream(line);

        // parse line
        //   Nex lambda mu 2Sp 2Sn 2S count
        int Nex, twice_Sp, twice_Sn, twice_S, lambda, mu, count;
        line_stream >> Nex >> twice_Sp >> twice_Sn >> twice_S >> lambda >> mu >> count;
        //line_stream >> Nex >> lambda >> mu >> twice_Sp >> twice_Sn >> twice_S >> count;
        ParsingCheck(line_stream, line_count, line);

        // conversions
        HalfInt Nsigma = Nsigma_0 + Nex;
        u3::U3 sigma(Nsigma,u3::SU3(lambda,mu));
        HalfInt Sp = HalfInt(twice_Sp,2);
        HalfInt Sn = HalfInt(twice_Sn,2);
        HalfInt S = HalfInt(twice_S,2);
      
        // replicate LGI according to count
        for (int i=1; i<=count; ++i)
          lgi_vector.emplace_back(Nex,sigma,Sp,Sn,S);
      }

    // close input file
    lgi_stream.close();
  }

  int NmaxTruncator::Nn_max(const u3::U3& sigma) const
  {
    HalfInt Nsigma = sigma.N();
    int value = Nmax_ - int(Nsigma - Nsigma_0_);
    return value;
  }

  void GenerateSp3RIrreps(
                          LGIVectorType& lgi_vector,
                          SigmaIrrepMapType& sigma_irrep_map,
                          const TruncatorInterface& truncator
                          )
  {
    // traverse LGI vector
    for (auto it = lgi_vector.begin(); it != lgi_vector.end(); ++it)
      {
        // retrieve sigma of LGI
        // u3::U3 sigma(it->sigma); // CRYPTIC
        spncci::LGI& lgi = *it;
        u3::U3 sigma = lgi.sigma;

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
          
        // save info back to LGI
        const sp3r::Sp3RSpace& irrep = sigma_irrep_map[sigma];
        lgi.SaveSubspaceInfo(irrep);

      }
  }

  ////////////////////////////////////////////////////////////////
  // iteration over Sp(3,R) -> U(3) subspaces & states
  ////////////////////////////////////////////////////////////////

  int TotalU3Subspaces(const LGIVectorType& lgi_vector)
  {
    int subspaces = 0;
    for (auto it=lgi_vector.begin(); it !=lgi_vector.end(); ++it)
      {
        const spncci::LGI& lgi = *it;
        const sp3r::Sp3RSpace& irrep = lgi.Sp3RSpace();
        subspaces += irrep.size();
      }

    return subspaces;
  }

  int TotalDimensionU3(const LGIVectorType& lgi_vector)
  {
    int dimension = 0;
    for (auto it=lgi_vector.begin(); it !=lgi_vector.end(); ++it)
      {
        const spncci::LGI& lgi = *it;
        const sp3r::Sp3RSpace& irrep = lgi.Sp3RSpace();
        for (int i_subspace = 0; i_subspace < irrep.size(); ++i_subspace)
          {
            const sp3r::U3Subspace& u3_subspace = irrep.GetSubspace(i_subspace);
            dimension += u3_subspace.size();
          }
      }

    return dimension;
  }

  int TotalDimensionU3LS(const LGIVectorType& lgi_vector)
  {
    int dimension = 0;
    for (auto it=lgi_vector.begin(); it !=lgi_vector.end(); ++it)
      {
        const spncci::LGI& lgi = *it;
        const sp3r::Sp3RSpace& irrep = lgi.Sp3RSpace();
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

  int TotalDimensionU3LSJConstrained(const LGIVectorType& lgi_vector, const HalfInt& J)
  {
    int dimension = 0;
    for (auto it=lgi_vector.begin(); it !=lgi_vector.end(); ++it)
      {
        const spncci::LGI& lgi = *it;
        HalfInt S = lgi.S;
        const sp3r::Sp3RSpace& irrep = lgi.Sp3RSpace();
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

  int TotalDimensionU3LSJAll(const LGIVectorType& lgi_vector)
  {
    int dimension = 0;
    for (auto it=lgi_vector.begin(); it !=lgi_vector.end(); ++it)
      {
        const spncci::LGI& lgi = *it;
        HalfInt S = lgi.S;
        const sp3r::Sp3RSpace& irrep = lgi.Sp3RSpace();
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

  std::vector< std::pair<int,int> > GenerateLGIPairs(spncci::LGIVectorType lgi_vector )
  {
    std::vector< std::pair<int,int> >  lgi_pair_vector;
    for(int i=0; i<lgi_vector.size(); i++)
      for (int j=0; j<=i; j++)
        {
          spncci::LGI lgi1=lgi_vector[i];
          spncci::LGI lgi2=lgi_vector[j];
          if (abs(lgi1.Sp-lgi2.Sp)<=2)
            if (abs(lgi1.Sn-lgi2.Sn)<=2)
              if (abs(lgi1.S-lgi2.S)<=2)
                {
                 
                 lgi_pair_vector.push_back(std::pair<int,int>(i,j));
                }


        }

     return lgi_pair_vector;
  }



}  // namespace

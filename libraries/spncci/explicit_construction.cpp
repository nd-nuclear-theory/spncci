/****************************************************************
  explicit_construction.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/explicit_construction.h"

#include "cppformat/format.h"
#include "mcutils/eigen.h"

namespace spncci
{

  void 
  ConstructSpNCCIBasisExplicit(
      const u3shell::SpaceU3SPN& lsu3shell_space,
      const basis::MatrixVector& lgi_expansions,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const u3shell::SectorsU3SPN& Arel_sectors,
      const basis::MatrixVector& Arel_matrices,
      basis::MatrixVector& spncci_expansions
    )
  {

    spncci_expansions.resize(baby_spncci_space.size());
    for (int subspace_index=0; subspace_index<baby_spncci_space.size(); ++subspace_index)
      {
        // extract subspace properties
        const BabySpNCCISubspace& subspace = baby_spncci_space.GetSubspace(subspace_index);
        int irrep_family_index = subspace.irrep_family_index();
        int Nex = int(subspace.omega().N()-subspace.sigma().N());

        // define aliases to the relevant lsu3shell subspaces
        int lgi_lsu3shell_subspace_index = lsu3shell_space.LookUpSubspaceIndex(subspace.sigmaSPN());
        int spncci_lsu3shell_subspace_index = lsu3shell_space.LookUpSubspaceIndex(subspace.omegaSPN());

        // diagnostics
        const u3shell::SubspaceU3SPN& lgi_lsu3shell_subspace = lsu3shell_space.GetSubspace(lgi_lsu3shell_subspace_index);
        const u3shell::SubspaceU3SPN& spncci_lsu3shell_subspace = lsu3shell_space.GetSubspace(spncci_lsu3shell_subspace_index);
        std::cout
          << fmt::format(
              "Constructing subspace: LGI sigmaSPN {} family index {} => omegaSPN {} (Nex {})",
              subspace.sigmaSPN().Str(),irrep_family_index,
              subspace.omegaSPN().Str(),spncci_lsu3shell_subspace.size(),Nex
            )
          << std::endl
          << fmt::format(
              "  lsu3shell subspace sigmaSPN: U3SPN {} subspace index {} dim {}",
              lgi_lsu3shell_subspace.U3SPN().Str(),
              lgi_lsu3shell_subspace_index,
              lgi_lsu3shell_subspace.size()
            )
          << std::endl
          << fmt::format(
              "  lsu3shell subspace omegaSPN: U3SPN {} subspace index {} dim {}",
             spncci_lsu3shell_subspace.U3SPN().Str(),
             spncci_lsu3shell_subspace_index,
             spncci_lsu3shell_subspace.size()
            )
          << std::endl;

        // define aliases to expansion matrices
        const Eigen::MatrixXd& lgi_expansion = lgi_expansions[irrep_family_index];
        Eigen::MatrixXd& spncci_expansion = spncci_expansions[subspace_index];

        // diagnostics
        // std::cout << fmt::format("  lgi_expansion ({},{})",lgi_expansion.rows(),lgi_expansion.cols()) << std::endl;

        // calculate expansion of baby SpNCCI subspace
        assert((Nex==0)||(Nex==2));
        if (Nex==0)
          // Nex=0 -- trivial expansion
          {
            spncci_expansion = lgi_expansion;
          }
        else if (Nex==2)
          // Nex=1 -- single laddering, free of outer multiplicity
          {
            // retrieve matrix for applicable Arel sector
            int Arel_sector_index = Arel_sectors.LookUpSectorIndex(
                spncci_lsu3shell_subspace_index,
                lgi_lsu3shell_subspace_index
              );
            assert(Arel_sector_index!=basis::kNone);
            const Eigen::MatrixXd& Arel_matrix = Arel_matrices[Arel_sector_index];

            // diagnostics
            // std::cout << fmt::format("  Arel sector index {} matrix ({},{})",Arel_sector_index,Arel_matrix.rows(),Arel_matrix.cols()) << std::endl;

            // retrieve applicable inverse K matrix
            //
            // Language issue: Subscripted access only works for nonconst map.
            // const Eigen::MatrixXd& k_matrix = k_matrix_cache[subspace.sigma()][subspace.omega()];
            const Eigen::MatrixXd& k_matrix = k_matrix_cache.at(subspace.sigma()).at(subspace.omega());

            assert((k_matrix.rows()==1)&&(k_matrix.cols()==1));
            double k_inverse = 1 / k_matrix(0,0);  // specialize to multiplicity-free branching

            // act with raising operator on LGI expansions
            int phase = ParitySign(
                u3::ConjugationGrade(subspace.sigma().SU3())
                +u3::ConjugationGrade(subspace.omega().SU3())
              );
            
            spncci_expansion = phase*k_inverse*Arel_matrix*lgi_expansion;
          }

      }
  }

}  // namespace

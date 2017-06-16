/****************************************************************
  computation_control.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/computation_control.h"

#include "SymEigsSolver.h"  // from spectra
#include "cppformat/format.h"
#include "lgi/lgi_solver.h"
#include "mcutils/eigen.h"

namespace spncci
{

  void ConstructBranchedObservables(
    const spncci::SpaceU3S& space_u3s,
    const std::vector<std::vector<spncci::SectorLabelsU3S>>& observable_sectors_u3s,
    const std::vector<basis::MatrixVector>& observable_matrices_u3s,
    std::map<HalfInt,spncci::SpaceLS>& spaces_lsj,
    int num_observables,
    const std::vector<HalfInt>& J_values,
    const std::vector<int>& observable_Jvalues,
    std::vector<std::map<spncci::JPair,spncci::MatrixType>>& observable_matrices
    )
  {
    // populate fully-branched many-body matrices for observables
    // map: observable -> J ->  matrix
    // std::vector<std::map<HalfInt,Eigen::MatrixXd>> observable_matrices;  
    observable_matrices.resize(num_observables);
    for (int observable_index=0; observable_index<num_observables; ++observable_index)
      {
        int J0=observable_Jvalues[observable_index];
        for (const HalfInt bra_J : J_values)
          for (const HalfInt ket_J : J_values)  
            {
              if(not am::AllowedTriangle(bra_J,J0,ket_J))
                continue;
              spncci::JPair JpJ(bra_J,ket_J);
              // set up aliases (for current observable and J space)
              const std::vector<spncci::SectorLabelsU3S>& sectors_u3s = observable_sectors_u3s[observable_index];
              const basis::MatrixVector& matrices_u3s = observable_matrices_u3s[observable_index];
              
              // determine allowed LS sectors
              const spncci::SpaceLS& bra_space_lsj = spaces_lsj[bra_J];
              const spncci::SpaceLS& ket_space_lsj = spaces_lsj[ket_J];

              Eigen::MatrixXd& observable_matrix = observable_matrices[observable_index][JpJ];

              // determine set of (L0,S0) labels for this observable (triangular with J0)
              std::vector<spncci::OperatorLabelsLS> operator_labels_ls;
              // Note: to update when J0 varies by observable
              spncci::GenerateOperatorLabelsLS(J0,operator_labels_ls);

              std::vector<spncci::SectorLabelsLS> sectors_lsj;
              spncci::GetSectorsLS(bra_space_lsj,ket_space_lsj,operator_labels_ls,sectors_lsj);

              // branch LS sectors to LSJ
              basis::MatrixVector matrices_lsj;  
              spncci::ContractAndRegroupLSJ(
                  bra_J,J0,ket_J,
                  space_u3s,sectors_u3s,matrices_u3s,
                  bra_space_lsj,ket_space_lsj,sectors_lsj,matrices_lsj
                );

              // collect LSJ sectors into J matrix
              //
              // Note: Interface needs to be generalized to handle J_bra != J_ket.
              ConstructOperatorMatrix(
                  bra_space_lsj,ket_space_lsj,sectors_lsj,matrices_lsj,
                  observable_matrix
                );
            }
      }
  }


void 
  SolveHamiltonian(
      const spncci::MatrixType& hamiltonian_matrix,
      const HalfInt& J,
      int num_eigenvalues,
      int eigensolver_num_convergence,  // whatever exactly this is...
      int eigensolver_max_iterations,
      double eigensolver_tolerance,
      std::map<HalfInt,Eigen::VectorXd>& eigenvalues,  // map: J -> eigenvalues
      std::map<HalfInt,spncci::MatrixType>& eigenvectors  // map: J -> eigenvectors
    )
  {    

    // set up aliases
    // spncci::MatrixType& hamiltonian_matrix = observable_matrices[0][J];
    std::cout << fmt::format("  Diagonalizing: J={}",J) << std::endl;

    // define eigensolver and compute
    typedef spncci::MatrixType MatrixType;  // allow for possible future switch to more compact single-precision matrix
    typedef double FloatType;
    Spectra::DenseSymMatProd<FloatType> matvec(hamiltonian_matrix);
    Spectra::SymEigsSolver<FloatType,Spectra::SMALLEST_ALGE,Spectra::DenseSymMatProd<FloatType>>
      eigensolver(
          &matvec,
          num_eigenvalues,
          eigensolver_num_convergence
        );
    eigensolver.init();
    int converged_eigenvectors = eigensolver.compute(
        eigensolver_max_iterations,
        eigensolver_tolerance,
        Spectra::SMALLEST_ALGE  // sorting rule
      );
    int eigensolver_status = eigensolver.info();
    std::cout
      << fmt::format("  Eigensolver reports: eigenvectors {} status {}",converged_eigenvectors,eigensolver_status)
      << std::endl;
    assert(converged_eigenvectors=eigensolver.eigenvalues().size());  // should this always be true?
    assert(converged_eigenvectors=eigensolver.eigenvectors().cols());  // should this always be true?

    // save eigenresults
    eigenvalues[J] = eigensolver.eigenvalues();
    eigenvectors[J] = eigensolver.eigenvectors();
    std::cout << fmt::format("  Eigenvalues (J={}):",J) << std::endl
              << mcutils::FormatMatrix(eigenvalues[J],"8.5f","    ")
              << std::endl;

    // check eigenvector norms
    Eigen::VectorXd eigenvector_norms(eigenvectors[J].cols());
    for (int eigenvector_index=0; eigenvector_index<converged_eigenvectors; ++eigenvector_index)
      {
        eigenvector_norms(eigenvector_index) = eigenvectors[J].col(eigenvector_index).norm();
        const double norm_tolerance=1e-8;
        assert(fabs(eigenvector_norms(eigenvector_index)-1)<norm_tolerance);
      }
      if (true)
        {
          std::cout << fmt::format("  Norms (J={}):",J) << std::endl
                    << mcutils::FormatMatrix(eigenvector_norms,"8.5f","    ")
                    << std::endl;
        }

      // normalize eigenvectors -- redundant with Spectra eigensolver
      for (int eigenvector_index=0; eigenvector_index<converged_eigenvectors; ++eigenvector_index)
        eigenvectors[J].col(eigenvector_index).normalize();

      // diagnostics
      if (false)
        {
          std::cout << fmt::format("  Eigenvectors -- norm (J={}):",J) << std::endl
                    << mcutils::FormatMatrix(eigenvectors[J],"8.5f","    ")
                    << std::endl;
        }
  }


}  // namespace

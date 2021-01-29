/****************************************************************
  decomposition.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/decomposition.h"
#include "mcutils/eigen.h"
#include "mcutils/profiling.h"
#include "spncci/results_output.h"

namespace spncci
{
  // Note on decomposition algorithm:
  //
  // Accumulation of probabilities is done for all the
  // eigenvectors in the subspace in parallel, by means of
  // "column-wise partial reduction operations".  See Eigen
  // documentation:
  //
  //   man/eigen/eigen-doc-3.2.8/group__TutorialReductionsVisitorsBroadcasting.html
  //
  // Note that column-wise partial reductions return a row vector.
  //
  // We wish to accumulate the probability contributions from a
  // group of rows (in the eigenvector matrix) corresponding to
  // the "substates" of a single basis "state", in the
  // terminology of basis/degenerate.h.  (These are marked by
  // the vertical side bar in the diagram below.)  We thus
  // calculate the .columnwise().squaredNorm() of these rows of
  // the eigenvector matrix.  The result is a row vector of the
  // summed squares of the amplitudes from these rows.  We must
  // accumulate these onto the appropriate row of probabilities
  // in the decomposition matrix.
  //
  // Source (eigenvectors_J):
  //   
  //                     (v1)_0     (v2)_0     (v3)_0
  //                     (v1)_1     (v2)_1     (v3)_1
  //                     ...
  //   offset ->       | (v1)_n     (v2)_n     (v3)_n
  //                /  | (v1)_n+1   (v2)_n+1   (v3)_n+1
  //   (degeneracy)-   | ...      
  //                \  | (v1)_n+d-1 (v2)_n+d-1 (v3)_n+d-1
  //                     ...
  //
  // Result of column-wise partial reduction:
  //
  //                     sum1       sum2       sum3
  //
  // Target (xxx_decompositions_J):
  //
  //                     P(v1;0)   P(v2;0)  P(v3;0)
  //                     P(v1;1)   P(v2;1)  P(v3;1)
  //                     ...
  //   accumulate to ->  P(v1;m)   P(v2;m)  P(v3;m)
  //                     ...



void CalculateNexDecompositions(
    const std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
    const std::vector<spncci::Matrix>& eigenvectors,
    std::vector<spncci::Matrix>& Nex_decompositions,
    HalfInt Nsigma0,
    int Nmax
  )
{
  for (int spj_space_index=0; spj_space_index<spaces_spbasis.size(); ++spj_space_index)
    // for each J subspace
    {
      int offset=0; 
      // set up aliases for current J subspace
      const spncci::SpaceSpBasis& spj_space = spaces_spbasis[spj_space_index];


      const spncci::Matrix& eigenvectors_J = eigenvectors[spj_space_index];
      spncci::Matrix& Nex_decompositions_J = Nex_decompositions[spj_space_index];

      // initialize decomposition matrix
      const int decomposition_size = Nmax+1;
      const int num_eigenvectors = eigenvectors_J.cols();
      Nex_decompositions_J = spncci::Matrix::Zero(decomposition_size,num_eigenvectors);

      // accumulate probability
      for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
        // for each (composite) state
        {
          // retrieve basis state information          
          const spncci::SubspaceSpBasis& spj_subspace=spj_space.GetSubspace(spj_subspace_index);
          for(int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
            {
              StateSpBasis spj_state(spj_subspace,spj_state_index);    
              int degeneracy = spj_state.degeneracy();
              int Nex = int(spj_state.omega().N()-Nsigma0);
              assert((0<=Nex)&&(Nex<=Nmax));
 
              // accumulate probability from this (composite) state
               Nex_decompositions_J.row(Nex)+=eigenvectors_J.block(offset,0,degeneracy,num_eigenvectors).colwise().squaredNorm();
              offset+=degeneracy;
            }
          

        }
          
        
    }
}


void CalculateLSDecompositions(
    const std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
    const std::vector<spncci::Matrix>& eigenvectors,
    double hw
    // std::vector<spncci::Matrix>& LS_decompositions
  )
{
  std::cout<<"Doing LS decompositions"<<std::endl;
  std::vector<spncci::Matrix> LS_decompositions(spaces_spbasis.size());
  // LS_decompositions.resize(spaces_spbasis.size());
  
  std::ofstream ls_decomposition_file;
  ls_decomposition_file.open (fmt::format("LS_decompositions_{:2.1f}.dat",hw));

  for (int spj_space_index=0; spj_space_index<spaces_spbasis.size(); ++spj_space_index)
    // for each J subspace
    {
      
      const SpaceSpBasis& spj_space=spaces_spbasis[spj_space_index];
      const spncci::Matrix& eigenvectors_J = eigenvectors[spj_space_index];
      spncci::Matrix& LS_decompositions_J = LS_decompositions[spj_space_index];
      // std::cout<<"J="<<spj_space.J()<<std::endl;
      //Get number of L subspaces and L values for given J
      std::map<std::pair<int,HalfInt>,int> LS_basis;
      int index=0;
      for(int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
        {
          const SubspaceSpBasis& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
          for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
            {
              // std::cout<<"retrieve basis state information"<<std::endl;
              StateSpBasis spj_state(spj_subspace,spj_state_index);
              int L=spj_state.L();
              HalfInt S= spj_state.S();
              std::pair<int,HalfInt>LS(L,S);
              if(LS_basis.count(LS)==0)
                {
                  LS_basis[LS]=index;
                  index++;
                }
            }
        }

      // std::cout<<"initialize decomposition matrix"<<std::endl;
      const int num_eigenvectors = eigenvectors_J.cols();
      // std::cout<<L_basis.size()<<"  "<<num_eigenvectors<<std::endl;
      LS_decompositions_J = spncci::Matrix::Zero(LS_basis.size(),num_eigenvectors);
      // std::cout<<"here"<<std::endl;
      int offset=0;
      for(int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
        {
          // std::cout<<"set up aliases for current J subspace"<<std::endl;
          const SubspaceSpBasis& spj_subspace = spj_space.GetSubspace(spj_subspace_index);

          // std::cout<<"accumulate probability"<<std::endl;
          for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
            // for each (composite) state
            {
              // retrieve basis state information
              StateSpBasis spj_state(spj_subspace,spj_state_index);
              std::pair<int,HalfInt>LS(spj_state.L(),spj_state.S());
              index=LS_basis[LS];
              int degeneracy = spj_state.degeneracy();

              // std::cout<<"accumulate probability from this (composite) state"<<std::endl;
              LS_decompositions_J.row(index) += eigenvectors_J.block(offset,0,degeneracy,num_eigenvectors).colwise().squaredNorm();
              // std::cout<<"----------------------"<<std::endl<<baby_spncci_subspace_index<<std::endl<<std::endl;
              // std::cout<<baby_spncci_decompositions_J<<std::endl<<std::endl;
              offset+=degeneracy; 
            }
        }

      mcutils::ChopMatrix(LS_decompositions_J,1e-6);
      ls_decomposition_file<<"J="<<spj_space.J().Str()<<std::endl;
      for(auto itr=LS_basis.begin(); itr!=LS_basis.end(); itr++)
        {
          int L;
          HalfInt S;
          std::tie(L,S) = itr->first;
          int index=itr->second;
          ls_decomposition_file<<L<<" "<<S<<"  "<<mcutils::FormatMatrix(LS_decompositions_J.row(index),".6f")<<std::endl;
        }
      ls_decomposition_file<<std::endl;

    }//end J loop
  ls_decomposition_file.close();

}


void CalculateBabySpNCCIDecompositions(
    const std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
    const std::vector<spncci::Matrix>& eigenvectors,
    std::vector<spncci::Matrix>& baby_spncci_decompositions,
    int baby_spncci_space_size
  )
{
  for (int spj_space_index=0; spj_space_index<spaces_spbasis.size(); ++spj_space_index)
    // for each J subspace
    {
      const SpaceSpBasis& spj_space=spaces_spbasis[spj_space_index];
      const spncci::Matrix& eigenvectors_J = eigenvectors[spj_space_index];
      spncci::Matrix& baby_spncci_decompositions_J = baby_spncci_decompositions[spj_space_index];
      // initialize decomposition matrix
      const int decomposition_size = baby_spncci_space_size;
      const int num_eigenvectors = eigenvectors_J.cols();
      baby_spncci_decompositions_J = spncci::Matrix::Zero(decomposition_size,num_eigenvectors);
      int offset=0;
      for(int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
        {
          // set up aliases for current J subspace
          const SubspaceSpBasis& spj_subspace = spj_space.GetSubspace(spj_subspace_index);

          // accumulate probability
          for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
            // for each (composite) state
            {
              // retrieve basis state information
              StateSpBasis spj_state(spj_subspace,spj_state_index);
              int degeneracy = spj_state.degeneracy();
              int baby_spncci_subspace_index = spj_state.baby_spncci_subspace_index();
              assert((0<=baby_spncci_subspace_index)&&(baby_spncci_subspace_index<baby_spncci_space_size));

              // accumulate probability from this (composite) state
              baby_spncci_decompositions_J.row(baby_spncci_subspace_index) += eigenvectors_J.block(offset,0,degeneracy,num_eigenvectors).colwise().squaredNorm();
              // std::cout<<"----------------------"<<std::endl<<baby_spncci_subspace_index<<std::endl<<std::endl;
              // std::cout<<baby_spncci_decompositions_J<<std::endl<<std::endl;
              offset+=degeneracy; 
            }
        }     
      mcutils::ChopMatrix(baby_spncci_decompositions_J,1e-10);
    }

}


  void GenerateDecompositions(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
    const spncci::RunParameters& run_parameters,
    const std::vector<spncci::Matrix>& eigenvectors,
    double hw,
    std::ofstream& results_stream
    )
  {
  
    //////////////////////////////////////////////////////////////
    // do decompositions
    ////////////////////////////////////////////////////////////////
    std::cout << "Calculate eigenstate decompositions..." << std::endl;
    mcutils::SteadyTimer timer_decompositions;
    timer_decompositions.Start();

    // decomposition matrices:
    //   - vector over J
    //   - matrix over (basis_subspace_index,eigenstate_index)
    //
    // That is, decompositions are stored as column vectors, within a
    // matrix, much like the eigenstates themselves.
    std::vector<spncci::Matrix> Nex_decompositions;
    std::vector<spncci::Matrix> baby_spncci_decompositions;
    Nex_decompositions.resize(run_parameters.J_values.size());
    baby_spncci_decompositions.resize(run_parameters.J_values.size());

    // calculate decompositions
    spncci::CalculateNexDecompositions(
      spaces_spbasis,eigenvectors,Nex_decompositions,
      run_parameters.Nsigma0,run_parameters.Nmax
    );

    spncci::CalculateBabySpNCCIDecompositions(
      spaces_spbasis,eigenvectors,baby_spncci_decompositions,
      baby_spncci_space.size()
    );
    
    //TEMP: to be replaced with full function later
    
    spncci::CalculateLSDecompositions(spaces_spbasis,eigenvectors,hw);


    // // results output: decompositions
    spncci::WriteDecompositions(
      results_stream,"Nex",".6f",spaces_spbasis,
      Nex_decompositions,run_parameters.gex
    );

    spncci::WriteDecompositions(
      results_stream,"BabySpNCCI",".4e",spaces_spbasis,
      baby_spncci_decompositions,run_parameters.gex
    );

    timer_decompositions.Stop();
    std::cout << fmt::format("  (Decompositions: {})",timer_decompositions.ElapsedTime()) << std::endl;
 
  }


  void GenerateDecompositions(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
    const spncci::RunParameters& run_parameters,
    // const std::vector<spncci::Matrix>& eigenvectors,
    double hw,
    std::ofstream& results_stream
    )
  {

    //Read in eigenvectors from files fro each J and store in matrix vector
    std::vector<spncci::Matrix> eigenvectors(run_parameters.J_values.size());
    for(int j=0; j<run_parameters.J_values.size(); ++j)
      {
        HalfInt J=run_parameters.J_values[j];
        std::string eigv_filename=fmt::format("eigenvector_{:02d}_{:02.1f}.dat",TwiceValue(J),hw);
        spncci::Matrix& eigenvectors_J=eigenvectors[j];
        spncci::ReadEigenvectors(eigv_filename,eigenvectors_J);
      }
    spncci::GenerateDecompositions(baby_spncci_space,spaces_spbasis,run_parameters,eigenvectors,hw,results_stream);

  }







}  // namespace

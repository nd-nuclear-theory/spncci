/****************************************************************
  lgi_solver.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include <fstream>
#include <omp.h>

#include "cppformat/format.h"
#include "lgi/lgi_solver.h"
#include "lgi/null_solver.h"
#include "utilities/utilities.h"

extern double zero_threshold;

namespace lgi
{

  void 
  GenerateBrelNcmMatrices(
      const u3shell::SpaceU3SPN& space, 
      const u3shell::SectorsU3SPN& Brel_sectors,
      const basis::MatrixVector& Brel_matrices,
      const u3shell::SectorsU3SPN& Ncm_sectors,
      const basis::MatrixVector& Ncm_matrices,
      basis::MatrixVector& BrelNcm_matrices 
    ) 
  // Construct the Brel+Ncm Matrices for each ket subspace and store
  // them in BrelNcm_matrices
  {    
    u3shell::U3SPN omegaSPNi, omegaSPNj;
    int i,j, start_index_i, start_index_j, group_size_i, group_size_j;
    double rme;

    // subspace_sectors is map of vectors of 
    // (dimension of Brel+Ncm matrix, indices of relevant N-2 subspace)
    //
    // Filling out subspace_sectors for Ncm. 
    // Note Ncm scalar so Ncm sectors are square.
    std::map<int,std::vector<int>> subspace_sectors;
    for(int j=0; j<space.size(); ++j)
      {
        int sub_dim=space.GetSubspace(j).size();
        subspace_sectors[j].push_back(sub_dim);
      }
    // Filling out subspace_sectors for Brel
    // Total dimension is accumulated
    for(int s=0; s<Brel_sectors.size(); ++s)
      {
        int i=Brel_sectors.GetSector(s).bra_subspace_index();
        int j=Brel_sectors.GetSector(s).ket_subspace_index();
        subspace_sectors[j].push_back(i);
        //increment size of matrix 
        int sub_dim=space.GetSubspace(i).size();
        subspace_sectors[j][0]+=sub_dim;
      }

    BrelNcm_matrices.resize(space.size());
    // for each subspace, construct the BrelNcm matrix and store in vector
    for(int j=0; j<space.size(); ++j)
      {
        // total number of rows in matrix
        int rows_total=subspace_sectors[j][0];
        auto subspace=space.GetSubspace(j);
        int columns=subspace.size();
        // std::cout<<"BrelNcm initialized"<<std::endl;
        BrelNcm_matrices[j]=Eigen::MatrixXd::Zero(rows_total,columns);
        // std::cout<<BrelNcm_matrices[j]<<std::endl<<std::endl;
        //Ncm block is square so rows=columns
        int rows=columns;
        // Because Ncm is SU(3) scalar, subspace index should equal sector index;
        // std::cout<<"Adding in Ncm"<<std::endl;
        BrelNcm_matrices[j].block(0,0,rows,columns)=Ncm_matrices[j];
        // std::cout<<Ncm_matrices[j]<<std::endl;
        // std::cout<<BrelNcm_matrices[j]<<std::endl<<std::endl;
        // increment location in full matrix
        int rows_begin=columns;
        // Fill in Brel sectors 

        for(int k=1; k<subspace_sectors[j].size(); ++k)
          {
            // subspace index for bra
            int i=subspace_sectors[j][k];
            // number of rows in sector
            rows=space.GetSubspace(i).size();
            // index of sector in vector
            int sector_index=Brel_sectors.LookUpSectorIndex(i,j,1);
            // filling in sector
            // std::cout<<"Adding Brel "<<k<<" of "<<subspace_sectors[j].size()-1<<std::endl;
            BrelNcm_matrices[j].block(rows_begin,0,rows,columns)=Brel_matrices[sector_index];
            // std::cout<<Brel_matrices[sector_index]<<std::endl<<std::endl;
            // std::cout<<BrelNcm_matrices[j]<<std::endl<<std::endl;
            rows_begin+=rows;
          }
      }
  }

  void 
    GenerateLGIExpansion(
        const u3shell::SpaceU3SPN& space, 
        const u3shell::SectorsU3SPN& Brel_sectors,
        const basis::MatrixVector& Brel_matrices,
        const u3shell::SectorsU3SPN& Ncm_sectors,
        const basis::MatrixVector& Ncm_matrices,
        HalfInt Nsigma_0,
        lgi::MultiplicityTaggedLGIVector& lgi_families,
        basis::MatrixVector& lgi_expansions
        // bool keep_empty_subspaces
      )
  {
    double threshold=10e-4;
   
    basis::MatrixVector BrelNcm_matrices;
    GenerateBrelNcmMatrices(
        space,
        Brel_sectors,Brel_matrices,Ncm_sectors,Ncm_matrices,
        BrelNcm_matrices
      );

    lgi_expansions.resize(BrelNcm_matrices.size());
    lgi_families.resize(BrelNcm_matrices.size());

    #pragma omp parallel for schedule(runtime)
    for(int i=0; i<BrelNcm_matrices.size();++i)
      {
        Eigen::MatrixXd null_vectors;
        lgi::FindNullSpaceSVD(BrelNcm_matrices[i],null_vectors,threshold);
        int nullity = null_vectors.cols();

        // save LGI labels, tagged by nullity as multiplicity
        u3shell::U3SPN labels(space.GetSubspace(i).labels());          
        int Nex=int(labels.N()-Nsigma_0);
        lgi_families[i]=MultiplicityTagged<lgi::LGI>(lgi::LGI(labels,Nex),nullity);
        // .emplace_back(lgi::LGI(labels,Nex),nullity);

        // save expansions for these LGIs
        lgi_expansions[i]=null_vectors;
      }

  }


  void 
    GenerateLGIExpansion(
        int A,
        HalfInt Nsigma_0,
        const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
        const u3shell::SpaceU3SPN& lsu3shell_space, 
        std::ifstream& is_Brel,
        std::ifstream& is_Nrel,
        lgi::MultiplicityTaggedLGIVector& lgi_families,
        basis::MatrixVector& lgi_expansions
        // bool keep_empty_subspaces
      )
  // DEPRECATED
  {

    // provide wrapper for clean, I/O-free version

    // read Brel -- DEPRECATED stream-based ReadLSU3ShellRMEs
    u3shell::OperatorLabelsU3ST Brel_labels(-2,u3::SU3(0,2),0,0,0);
    u3shell::SectorsU3SPN Brel_sectors(lsu3shell_space,Brel_labels,true);
    basis::MatrixVector Brel_matrices;
    lsu3shell::ReadLSU3ShellRMEs(
        is_Brel,Brel_labels,
        lsu3shell_basis_table,lsu3shell_space,
        Brel_sectors,Brel_matrices
      );

    // read Nrel as Ncm -- DEPRECATED stream-based GenerateNcmMatrixVector
    u3shell::OperatorLabelsU3ST Ncm_labels(0,u3::SU3(0,0),0,0,0);
    u3shell::SectorsU3SPN Ncm_sectors(lsu3shell_space,Ncm_labels,true);
    basis::MatrixVector Ncm_matrices;
    lsu3shell::GenerateNcmMatrixVector(
        A,
        is_Nrel,
        lsu3shell_basis_table,lsu3shell_space,
        Ncm_matrices
      );

    // call matrix-based version
    lgi::GenerateLGIExpansion(
        lsu3shell_space, 
        Brel_sectors,Brel_matrices,Ncm_sectors,Ncm_matrices,
        Nsigma_0,
        lgi_families,lgi_expansions
      );

  }

  void
  TransformOperatorToSpBasis(
      const u3shell::SectorsU3SPN& sectors,
      const basis::MatrixVector& basis_transformation_matrices,
      const basis::MatrixVector& lsu3shell_operator_matrices,
      basis::MatrixVector& spncci_operator_matrices
    )
  {
    // for each sector, look up bra and ket subspaces 
    spncci_operator_matrices.resize(lsu3shell_operator_matrices.size());
    
    #pragma omp parallel for schedule(runtime)
    for(int s=0; s<lsu3shell_operator_matrices.size(); ++s)
      {
        int i=sectors.GetSector(s).bra_subspace_index();
        int j=sectors.GetSector(s).ket_subspace_index();

        // get transformation matrices and transpose bra transformation matrix
        const Eigen::MatrixXd& bra=basis_transformation_matrices[i].transpose();
        const Eigen::MatrixXd& ket=basis_transformation_matrices[j];

        // transform operator to spncci basis
        spncci_operator_matrices[s]=bra*lsu3shell_operator_matrices[s]*ket;
      }
  }

}// end namespace

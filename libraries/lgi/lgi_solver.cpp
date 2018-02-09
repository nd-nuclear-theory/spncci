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
      lsu3shell::OperatorBlocks& BrelNcm_matrices 
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
        basis::MatrixVector& lgi_expansions,
        std::vector<int>& lsu3shell_index_lookup_table
      )
  {
    double threshold=10e-6;
   
    basis::MatrixVector BrelNcm_matrices;
    GenerateBrelNcmMatrices(
        space,
        Brel_sectors,Brel_matrices,Ncm_sectors,Ncm_matrices,
        BrelNcm_matrices
      );

    lsu3shell_index_lookup_table.resize(space.size());
    lgi_expansions.resize(space.size());
    // Iterates over lsu3shell subspaces.  If there exist null vectors of B+N
    // matrix, add subspace labels to list of lgi and save expansion
    for(int i=0; i<BrelNcm_matrices.size();++i)
      {
        // std::cout<<BrelNcm_matrices[i]<<std::endl;
        Eigen::MatrixXd null_vectors;
        lgi::FindNullSpaceSVD(BrelNcm_matrices[i],null_vectors,threshold);
        // std::cout<<null_vectors<<std::endl;
        int nullity = null_vectors.cols();
        // std::cout<<nullity<<std::endl;

        // if nullity>0 save lgi labels to lgi_familes 
        if(nullity>0)
          {
            // save LGI labels, tagged by nullity as multiplicity
            u3shell::U3SPN labels(space.GetSubspace(i).labels());          
            int Nex=int(labels.N()-Nsigma_0);
            // std::cout<<labels.N()<<"  "<<Nsigma_0<<"  "<<Nex<<std::endl;
            lgi_families.emplace_back(lgi::LGI(labels,Nex),nullity);

            // Record lgi index in look up table
            lsu3shell_index_lookup_table[i]=lgi_families.size()-1;
          }
        // save -1 to look up table to indicate that no lgi exists. 
        else
          lsu3shell_index_lookup_table[i]=-1;

        // Save expansion even if its a zero vector
        lgi_expansions[i]=null_vectors;
      }

  }


  void GetLGIExpansion(
      const u3shell::SpaceU3SPN& lsu3shell_space, 
      const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
      const std::string& Brel_filename,
      const std::string& Nrel_filename,
      int A, HalfInt Nsigma_0,
      lgi::MultiplicityTaggedLGIVector& lgi_families,
      lsu3shell::OperatorBlocks& lgi_expansions,
      std::vector<int>& lsu3hsell_index_lookup_table
    )
  {
    u3shell::SectorsU3SPN Bintr_sectors, Nintr_sectors;
    basis::MatrixVector Bintr_matrices, Nintr_matrices;
    bool sp3r_generators=true;

    // Read in Brel and Nrel calculated in the lsu3shell basis 
    lsu3shell::ReadLSU3ShellSymplecticOperatorRMEs(
        lsu3shell_basis_table,lsu3shell_space, 
        Brel_filename,Bintr_sectors,Bintr_matrices,
        Nrel_filename,Nintr_sectors,Nintr_matrices,
        A
      );

    // From Nintr construct Ncm
    const u3shell::SectorsU3SPN& Ncm_sectors = Nintr_sectors;
    basis::MatrixVector Ncm_matrices;
    // Note that A-1 is used here since we are generating Ncm using an
    // instric space 
    lsu3shell::GenerateLSU3ShellNcmRMEs(
        lsu3shell_space,Nintr_sectors,Nintr_matrices,
        A-1,Ncm_matrices
      );

    // Apply the lgi solver to obtain lgi_families and their expansions in the lsu3shell basis 
    lgi::GenerateLGIExpansion(
        lsu3shell_space, 
        Bintr_sectors,Bintr_matrices,Ncm_sectors,Ncm_matrices,
        Nsigma_0,lgi_families,lgi_expansions,
        lsu3hsell_index_lookup_table
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

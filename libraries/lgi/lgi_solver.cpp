/****************************************************************
  lgi_solver.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include <fstream>

#include "cppformat/format.h"
#include "lgi/lgi_solver.h"
#include "lgi/null_solver.h"
#include "utilities/utilities.h"

extern double zero_threshold;

namespace lgi
{

  void 
  GenerateBrelNcmMatrices(
      int A,
      std::ifstream& is_brel,
      std::ifstream& is_nrel,
      const lsu3shell::LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      basis::MatrixVector& BrelNcm_vector 
    ) 
  {    
    u3shell::U3SPN omegaSPNi, omegaSPNj;
    int i,j, start_index_i, start_index_j, group_size_i, group_size_j;
    double rme;

    u3shell::OperatorLabelsU3ST brel_labels(-2,u3::SU3(0,2),0,0,0);
    //generate sectors for brel.
    u3shell::SectorsU3SPN brel_sectors(space,brel_labels,true);

    // std::cout<<"Reading in Brel"<<std::endl;
    basis::MatrixVector brel_matrix_vector(space.size());
    lsu3shell::ReadLSU3ShellRMEs(is_brel,brel_labels,lsu3_basis_table,space, brel_sectors,brel_matrix_vector);
    // std::cout<<"Generating Ncm"<<std::endl;
    basis::MatrixVector ncm_matrix_vector(space.size());    
    lsu3shell::GenerateNcmMatrixVector(A,is_nrel,lsu3_basis_table,space, ncm_matrix_vector);

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
    for(int s=0; s<brel_sectors.size(); ++s)
      {
        int i=brel_sectors.GetSector(s).bra_subspace_index();
        int j=brel_sectors.GetSector(s).ket_subspace_index();
        subspace_sectors[j].push_back(i);
        //increment size of matrix 
        int sub_dim=space.GetSubspace(i).size();
        subspace_sectors[j][0]+=sub_dim;
      }

    // for each subspace, construct the BrelNcm matrix and store in vector
    for(int j=0; j<space.size(); ++j)
      {
        // total number of rows in matrix
        int rows_total=subspace_sectors[j][0];
        auto subspace=space.GetSubspace(j);
        int columns=subspace.size();
        // std::cout<<"BrelNcm initialized"<<std::endl;
        BrelNcm_vector[j]=Eigen::MatrixXd::Zero(rows_total,columns);
        // std::cout<<BrelNcm_vector[j]<<std::endl<<std::endl;
        //Ncm block is square so rows=columns
        int rows=columns;
        // Because Ncm is SU(3) scalar, subspace index should equal sector index;
        // std::cout<<"Adding in Ncm"<<std::endl;
        BrelNcm_vector[j].block(0,0,rows,columns)=ncm_matrix_vector[j];
        // std::cout<<ncm_matrix_vector[j]<<std::endl;
        // std::cout<<BrelNcm_vector[j]<<std::endl<<std::endl;
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
            int sector_index=brel_sectors.LookUpSectorIndex(i,j,1);
            // filling in sector
            // std::cout<<"Adding Brel "<<k<<" of "<<subspace_sectors[j].size()-1<<std::endl;
            BrelNcm_vector[j].block(rows_begin,0,rows,columns)=brel_matrix_vector[sector_index];
            // std::cout<<brel_matrix_vector[sector_index]<<std::endl<<std::endl;
            // std::cout<<BrelNcm_vector[j]<<std::endl<<std::endl;
            rows_begin+=rows;
          }
      }
  }

  void
  GenerateLGIExpansion(
      int A,
      HalfInt Nsigma_0,
      const lsu3shell::LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      std::ifstream& is_brel,
      std::ifstream& is_nrel,
      lgi::MultiplicityTaggedLGIVector& lgi_vector,
      basis::MatrixVector& lgi_expansion_matrix_vector,
      bool keep_empty_subspaces
    )
  // Construct Brel and Ncm matrix in lsu3shell basis and solve for null space.
  // Columns of kernel are expansion coefficients for each lgi.
  {
    double threshold=10e-4;
   
    basis::MatrixVector BrelNcm_vector(space.size());
    GenerateBrelNcmMatrices(A,is_brel,is_nrel,lsu3_basis_table, space, BrelNcm_vector);

    for(int i=0; i<BrelNcm_vector.size();++i)
      {
        Eigen::MatrixXd null_vectors;
        lgi::FindNullSpaceSVD(BrelNcm_vector[i],null_vectors,threshold);
        int nullity = null_vectors.cols();

        if ((nullity>0) || keep_empty_subspaces)
          {
            // save LGI labels, tagged by nullity as multiplicity
            u3shell::U3SPN labels(space.GetSubspace(i).GetSubspaceLabels());          
            int Nex=int(labels.N()-Nsigma_0);
            lgi_vector.emplace_back(lgi::LGI(labels,Nex),nullity);

            // save expansions for these LGIs
            lgi_expansion_matrix_vector.push_back(null_vectors);
          }
      }

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

  void 
  WriteLGILabels(const lgi::MultiplicityTaggedLGIVector& lgi_vector,   std::ofstream& os)
  {
    int Nex;
    u3::U3 sigma;
    HalfInt Sp,Sn,S;
    // std::unordered_map<lgi::LGI,int,boost::hash<lgi::LGI>> lgi_counter;
    for(auto lgi_count : lgi_vector)
      {
        std::tie(Nex,sigma,Sp,Sn,S)=lgi_count.irrep.Key();
        int count=lgi_count.tag;
        os
          <<Nex
          <<"  "<<TwiceValue(sigma.N())<<"  "<<sigma.SU3().lambda()<<"  "<<sigma.SU3().mu()
          <<"  "<<TwiceValue(Sp)<<"  "<<TwiceValue(Sn)<<"  "<<TwiceValue(S)
          <<"  "<<count<<std::endl;     
      }
    }

}// end namespace

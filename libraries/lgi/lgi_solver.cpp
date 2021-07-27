/****************************************************************
  lgi_solver.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/
#include <fstream>
#include <iostream>
#include <omp.h>

#include "fmt/format.h"
#include "lsu3shell/lsu3shell_rme.h"
#include "lgi/lgi_solver.h"
#include "lgi/null_solver.h"
#include "utilities/utilities.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"

namespace lgi
{

  void 
  GenerateBrelNcmMatrices(
      const u3shell::SpaceU3SPN& space_bra,
      const u3shell::SpaceU3SPN& space_ket, 
      const u3shell::SectorsU3SPN& Brel_sectors,
      const basis::OperatorBlocks<double>& Brel_matrices,
      const u3shell::SectorsU3SPN& Ncm_sectors,
      const basis::OperatorBlocks<double>& Ncm_matrices,
      lsu3shell::OperatorBlocks& BrelNcm_matrices 
    ) 
  // Construct the Brel+Ncm Matrices for each ket subspace and store
  // them in BrelNcm_matrices
  {    
    // Iterate through Brel and Nrel sectors to identify the relevant 
    // sectors need for constructing the Ncm+Bintr matrix for each
    // U3SPN subspace in space_ket.  Indices of bra subspaces needed
    // are accumulate in subspace_sectors for each ket subspace. 
    // Since Ncm is scalar the sector index should match ket subspace index
    // so the 0th element of the subspace_sectors vector (corresponding to Ncm subspace)
    // is instead used to accumulate the total dimension. 
    //
    // Accumulating dimensions of Ncm bra in the 0th position
    std::vector<std::vector<int>> subspace_sectors(space_ket.size());
    for(int j=0; j<space_ket.size(); ++j)
      {
        int sub_dim=space_ket.GetSubspace(j).size();
        subspace_sectors[j].push_back(sub_dim);
      }
    
    // Filling out subspace_sectors for Brel
    // Total dimension is accumulated in the 0th element of the vector
    for(int s=0; s<Brel_sectors.size(); ++s)
      {
        int i=Brel_sectors.GetSector(s).bra_subspace_index();
        int j=Brel_sectors.GetSector(s).ket_subspace_index();
        subspace_sectors[j].push_back(i);
        //increment size of matrix 
        int sub_dim=space_bra.GetSubspace(i).size();
        subspace_sectors[j][0]+=sub_dim;
      }

    BrelNcm_matrices.resize(space_ket.size());

    // for each subspace, construct the BrelNcm matrix and store in vector
    for(int j=0; j<space_ket.size(); ++j)
      {
        // total number of rows in matrix
        const auto& subspace=space_ket.GetSubspace(j);
        int total_rows=subspace_sectors[j][0];
        int cols=subspace.size();
        BrelNcm_matrices[j]=Eigen::MatrixXd::Zero(total_rows,cols);
        
        //Ncm block is square so rows=columns and subspace index should equal sector index;
        int rows=cols;
        BrelNcm_matrices[j].block(0,0,rows,cols)=Ncm_matrices[j];

        // increment location in full matrix
        int row_index=cols;

        // Fill in Brel sectors 
        for(int k=1; k<subspace_sectors[j].size(); ++k)
          {
            // subspace index for bra
            int i=subspace_sectors[j][k];
            rows=space_bra.GetSubspace(i).size();
            
            // Look up Brel sector index and add corresponding block to matrix
            int sector_index=Brel_sectors.LookUpSectorIndex(i,j,1);
            BrelNcm_matrices[j].block(row_index,0,rows,cols)=Brel_matrices[sector_index];

            // Increment row index
            row_index+=rows;
          }
      }
  }

  void 
    GenerateLGIExpansion(
        const u3shell::SpaceU3SPN& space_bra,
        const u3shell::SpaceU3SPN& space_ket, 
        const u3shell::SectorsU3SPN& Brel_sectors,
        const basis::OperatorBlocks<double>& Brel_matrices,
        const u3shell::SectorsU3SPN& Ncm_sectors,
        const basis::OperatorBlocks<double>& Ncm_matrices,
        HalfInt Nsigma_0,
        lgi::MultiplicityTaggedLGIVector& lgi_families,
        basis::OperatorBlocks<double>& lgi_expansions
      )
  {
    double threshold=1e-6;
   
    basis::OperatorBlocks<double> BrelNcm_matrices;
    GenerateBrelNcmMatrices(
        space_bra,space_ket,
        Brel_sectors,Brel_matrices,Ncm_sectors,Ncm_matrices,
        BrelNcm_matrices
      );
    
    lgi_expansions.resize(BrelNcm_matrices.size());

    for(int i=0; i<BrelNcm_matrices.size();++i)
      {
        //Solve for the null spaces of Bintr+Ncm
        Eigen::MatrixXd null_vectors;
        lgi::FindNullSpaceSVD(BrelNcm_matrices[i],null_vectors,threshold);

        int nullity = null_vectors.cols();

         // if nullity>0 save lgi labels to lgi_familes 
        if(nullity>0)
          {
            // save LGI labels, tagged by nullity as multiplicity
            u3shell::U3SPN labels(space_ket.GetSubspace(i).labels());          
            int Nex=int(labels.N()-Nsigma_0);
            lgi_families.emplace_back(lgi::LGI(labels,Nex),nullity);
          }

        // Save expansion even if its a zero vector so index matches up with su3shell basis indexing
        lgi_expansions[i]=null_vectors;
      }

  }


  void GetLGIExpansion(
      const u3shell::SpaceU3SPN& lsu3shell_space_bra,
      const u3shell::SpaceU3SPN& lsu3shell_space_ket, 
      const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table_bra,
      const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table_ket,
      const std::string& Brel_filename,
      const std::string& Nrel_filename,
      int A, HalfInt Nsigma_0,
      lgi::MultiplicityTaggedLGIVector& lgi_families,
      lsu3shell::OperatorBlocks& lgi_expansions
    )
  {
    
    bool sp3r_generators=true;
    
    // std::cout<<"Read Brel RMEs"<<std::endl;
    basis::OperatorBlocks<double> Bintr_matrices;
    u3shell::OperatorLabelsU3ST Brel_labels(-2,u3::SU3(0,2),0,0,0);
    u3shell::SectorsU3SPN Bintr_sectors(lsu3shell_space_bra,lsu3shell_space_ket,Brel_labels,true);
    lsu3shell::ReadLSU3ShellRMEs(
        Brel_filename,
        lsu3shell_basis_table_bra,lsu3shell_space_bra,
        lsu3shell_basis_table_ket,lsu3shell_space_ket,
        Brel_labels,Bintr_sectors,Bintr_matrices,
        2./A,sp3r_generators
      );

    // std::cout<<"Read Nrel RMEs" <<std::endl;
    basis::OperatorBlocks<double> Nintr_matrices;
    u3shell::OperatorLabelsU3ST Nrel_labels(0,u3::SU3(0,0),0,0,0);
    u3shell::SectorsU3SPN Nintr_sectors(lsu3shell_space_bra,lsu3shell_space_ket,Nrel_labels,true);
    lsu3shell::ReadLSU3ShellRMEs(
        Nrel_filename,
        lsu3shell_basis_table_bra,lsu3shell_space_bra,
        lsu3shell_basis_table_ket,lsu3shell_space_ket,
        Nrel_labels,Nintr_sectors,Nintr_matrices,
        2./A,sp3r_generators
      );

    // std::cout<<"From Nintr construct Ncm"<<std::endl;
    const u3shell::SectorsU3SPN& Ncm_sectors = Nintr_sectors;
    basis::OperatorBlocks<double> Ncm_matrices;
    
    // Note that A-1 is used here since we are generating Ncm using an intrinsic space 
    // Since sectors are diagonal, we only need to pass the lsu3shell_space_ket
    lsu3shell::GenerateLSU3ShellNcmRMEs(
        lsu3shell_space_ket,Nintr_sectors,Nintr_matrices,
        A-1,Ncm_matrices
      );

    // Apply the lgi solver to obtain lgi_families and their expansions in the lsu3shell basis 
    // std::cout<<"Get expansion"<<std::endl;
    lgi::GenerateLGIExpansion(
        lsu3shell_space_bra,
        lsu3shell_space_ket, 
        Bintr_sectors,Bintr_matrices,
        Ncm_sectors,Ncm_matrices,
        Nsigma_0,lgi_families,lgi_expansions
      );
  }


void WriteLGIExpansions(const std::string& filename, const lsu3shell::OperatorBlocks& lgi_expansions)
  {

    // num_rows, num_cols, rmes
    // rmes are by column then by row
    // output in binary mode 
    std::ios_base::openmode mode_argument = std::ios_base::out;
    mode_argument |= std::ios_base::binary;
    std::ofstream expansion_file;
    expansion_file.open(filename,mode_argument);

    if (!expansion_file)
     {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
     }
   
    mcutils::WriteBinary<int>(expansion_file,lgi::binary_format_code);
    // floating point precision
    mcutils::WriteBinary<int>(expansion_file,lgi::binary_float_precision);


    for(auto& block : lgi_expansions)
      {
        int num_rows=block.rows();
        int num_cols=block.cols();

        assert(num_rows==static_cast<lgi::LGIIndexType>(num_rows));
        mcutils::WriteBinary<lgi::LGIIndexType>(expansion_file,num_rows);

        assert(num_cols==static_cast<lgi::LGIIndexType>(num_cols));
        mcutils::WriteBinary<lgi::LGIIndexType>(expansion_file,num_cols);
        
        for(int j=0; j<num_cols; ++j)
          for(int i=0; i<num_rows; ++i)
            {
              auto rme=block(i,j);

              if (lgi::binary_float_precision==4)
                mcutils::WriteBinary<float>(expansion_file,rme);
              else if (lgi::binary_float_precision==8)
                mcutils::WriteBinary<double>(expansion_file,rme);
            }
      }
  }

      

void WriteLGIExpansions(const std::string& filename, const lsu3shell::OperatorBlock& lgi_expansion)
  {

    // num_rows, num_cols, rmes
    // rmes are by column then by row
    // output in binary mode 
    std::ios_base::openmode mode_argument = std::ios_base::out;
    mode_argument |= std::ios_base::binary;
    std::ofstream expansion_file;
    expansion_file.open(filename,mode_argument);

    if (!expansion_file)
     {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
     }
   
    // floating point precision
    mcutils::WriteBinary<int>(expansion_file,lgi::binary_float_precision);

    int num_rows=lgi_expansion.rows();
    int num_cols=lgi_expansion.cols();

    assert(num_rows==static_cast<lgi::LGIIndexType>(num_rows));
    mcutils::WriteBinary<lgi::LGIIndexType>(expansion_file,num_rows);

    assert(num_cols==static_cast<lgi::LGIIndexType>(num_cols));
    mcutils::WriteBinary<lgi::LGIIndexType>(expansion_file,num_cols);
    
    for(int j=0; j<num_cols; ++j)
      for(int i=0; i<num_rows; ++i)
        {
          auto rme=lgi_expansion(i,j);

          if (lgi::binary_float_precision==4)
            mcutils::WriteBinary<float>(expansion_file,rme);
          else if (lgi::binary_float_precision==8)
            mcutils::WriteBinary<double>(expansion_file,rme);
        }
  }
  
void WriteLGIExpansionsText(const std::string& filename, const lsu3shell::OperatorBlock& lgi_expansion)
  {

    // num_rows, num_cols, rmes
    // rmes are by column then by row    
    std::ofstream expansion_file;
    expansion_file.open(filename);

    if (!expansion_file)
     {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
     }
   
    int num_rows=lgi_expansion.rows();
    int num_cols=lgi_expansion.cols();
    expansion_file << num_rows << num_cols<<std::endl;
    expansion_file<<mcutils::FormatMatrix(lgi_expansion, ".8f")<<std::endl;        
  }





void ReadLGIExpansion(int num_lgi_subspaces,const std::string& filename, basis::OperatorBlocks<double>& lgi_expansions)
  {

    // num_rows, num_cols, rmes
    // rmes are by column then by row
    // output in binary mode 
    std::ios_base::openmode mode_argument = std::ios_base::in;
    mode_argument |= std::ios_base::binary;
    std::ifstream in_stream(filename, mode_argument);
    mcutils::StreamCheck(bool(in_stream),filename,"Failure opening lsu3shell rme file");
    // if (!in_stream)
    //  {
    //     std::cerr << "Could not open file '" << filename << "'!" << std::endl;
    //     return;
    //  }
   
    int binary_format_code;
    int binary_float_precision;
    mcutils::ReadBinary<int>(in_stream,binary_format_code);
    mcutils::ReadBinary<int>(in_stream,binary_float_precision);


    lgi_expansions.resize(num_lgi_subspaces);
    for(int i=0; i<num_lgi_subspaces; ++i)
      {
        basis::OperatorBlock<double>& block=lgi_expansions[i];
        lgi::LGIIndexType rows, cols;
        // Read in number of rows and cols
        mcutils::ReadBinary<lgi::LGIIndexType>(in_stream,rows);
        mcutils::ReadBinary<lgi::LGIIndexType>(in_stream,cols);

        /////////////////////////////////////////////////////////////////////////////////
        // Read in RMEs and cast to double matrix 
        if(binary_float_precision==4)
          {
            float buffer[rows*cols];
            in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
            block=Eigen::Map<Eigen::MatrixXf>(buffer,rows,cols).cast<double>();
          }
        else if (binary_float_precision==8)
          {
            double buffer[rows*cols];
            in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
            block=Eigen::Map<Eigen::MatrixXd>(buffer,rows,cols);
          }
      }
  }

}// end namespace

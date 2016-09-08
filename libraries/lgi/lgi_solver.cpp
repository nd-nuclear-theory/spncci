/****************************************************************
  lgi_solver.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "lgi/lgi_solver.h"


#include <fstream>
#include <iostream>
#include <algorithm>
//#include <functional>

#include "cppformat/format.h"

#include "u3shell/two_body_operator.h"
#include "u3shell/moshinsky.h"
// #include "spncci/sp_basis.h"
#include "basis/operator.h"
#include "lgi/lgi.h"
// #include "u3shell/u3spn_scheme.h"

namespace lgi
{
  void GenerateNcmMatrixVector(      
    const std::string& nrel_filename,
    const lsu3shell::LSU3BasisTable& lsu3_basis_table,
    const u3shell::SpaceU3SPN& space, 
    basis::MatrixVector& matrix_vector 
  )
  {

    std::ifstream is_nrel(nrel_filename.c_str());
    basis::MatrixVector nrel_matrix_vector;
    u3shell::OperatorLabelsU3S nrel_labels(0,u3::SU3(0,0),0,0);
    u3shell::SectorsU3SPN nrel_sectors(space,nrel_labels,true);
    lsu3shell::ReadLSU3ShellRMEs(is_nrel,nrel_labels,lsu3_basis_table,space, nrel_sectors,nrel_matrix_vector);

    basis::MatrixVector ncm_matrix_vector;

    for(int i=0; i<nrel_matrix_vector.size(); ++i)
    {
      auto subspace=space.GetSubspace(i);
      HalfInt N=subspace.N();
      int dim=subspace.size();
      matrix_vector[i]=Eigen::MatrixXd::Identity(dim,dim)*double(N)-nrel_matrix_vector[i];
    }
  }


  void 
  GenerateBrelNcmMatrices(
      const std::string& brel_filename,
      const std::string& nrel_filename,
      const lsu3shell::LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      basis::MatrixVector& BrelNcm_vector 
    ) 
  {    
    u3shell::U3SPN omegaSPNi, omegaSPNj;
    int i,j, start_index_i, start_index_j, group_size_i, group_size_j;
    double rme;

    u3shell::OperatorLabelsU3S brel_labels(-2,u3::SU3(0,2),0,0);
    //generate sectors for brel.
    u3shell::SectorsU3SPN brel_sectors(space,brel_labels,true);
    // subspace_sectors is map of vector of (dimension of Brel+Ncm matrix, index for Ncm, indices of relevant N-2 subspace)
    std::map<int,std::vector<int>> subspace_sectors;
    for(int j=0; j<space.size(); ++j)
      {
        int sub_dim=space.GetSubspace(j).size();
        subspace_sectors[j].push_back(sub_dim);
      }

    for(int s=0; s<brel_sectors.size(); ++s)
      {
        int i=brel_sectors.GetSector(s).bra_subspace_index();
        int j=brel_sectors.GetSector(s).ket_subspace_index();
        subspace_sectors[j].push_back(i); 
        int sub_dim=space.GetSubspace(i).size();
        subspace_sectors[j][0]+=sub_dim;
      }

    std::ifstream is_brel(brel_filename.c_str());
    basis::MatrixVector brel_matrix_vector;
    lsu3shell::ReadLSU3ShellRMEs(is_brel,brel_labels,lsu3_basis_table,space, brel_sectors,brel_matrix_vector);
    
    basis::MatrixVector ncm_matrix_vector;    
     GenerateNcmMatrixVector(nrel_filename,lsu3_basis_table,space, ncm_matrix_vector);

    int rows_begin=0;
    // for each subspace, construct the BrelNcm matrix and store in vector
    for(int j=0; j<space.size(); ++j)
      {
        // total number of rows in matrix
        int rows_total=subspace_sectors[j][0];
        auto subspace=space.GetSubspace(j);
        int columns=subspace.size();

        BrelNcm_vector[j]=Eigen::MatrixXd::Zero(rows_total,columns);
        //Ncm block is square so rows=columns
        int rows=columns;
        // Because Ncm is SU(3) scalar, subspace index should equal sector index;
        BrelNcm_vector[j].block(0,0,rows,columns)=ncm_matrix_vector[j];
        // increment location in full matrix
        rows_begin+=columns;
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
            BrelNcm_vector[j].block(rows_begin,0,rows,columns)=brel_matrix_vector[sector_index];
          }
      }
  }
  
  void GenerateLGIExpansion(
      int Nsigma_0,
      const lsu3shell::LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      const std::string& brel_filename,
      const std::string& nrel_filename,
      basis::MatrixVector& lgi_expansion_matrix_vector
      // lgi::LGIVector lgi_vector
      
    )
  // Construct Brel and Ncm matrix in lsu3shell basis and solve for null space.
  // Columns of kernel are expansion coefficients for each lgi.
  {
    basis::MatrixVector BrelNcm_vector; 
    GenerateBrelNcmMatrices(brel_filename,nrel_filename,lsu3_basis_table, space, BrelNcm_vector);

    for(int i=0; i<=BrelNcm_vector.size();++i)
      {
        lgi_expansion_matrix_vector[i]=BrelNcm_vector[i].fullPivLu().kernel();
        const u3shell::SubspaceU3SPN& subspace=space.GetSubspace(i);

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

  // void 
  // WriteLGI(const spncci::LGIVectorType& lgi_vector,   std::ofstream& os)
  // {
  //   int Nex;
  //   u3::U3 sigma;
  //   HalfInt Sp,Sn,S;
  //   std::unordered_map<spncci::LGI,int,boost::hash<spncci::LGI>> lgi_counter;
  //   for(auto a:lgi_vector)
  //     {
  //       lgi_counter[a]+=1;
  //     }
  //   for(auto b:lgi_vector)
  //     {
  //       std::tie(Nex,sigma,Sp,Sn,S)=b.Key();
  //       int count=lgi_counter[b];
  //       os
  //         <<Nex<<"  "<<TwiceValue(Sp)<<"  "<<TwiceValue(Sn)<<"  "<<TwiceValue(S)
  //         <<"  "<<sigma.SU3().lambda()<<"  "<<sigma.SU3().mu()<<"  "<<count<<std::endl;     
  //     }
  // }

  // void 
  // LGINex0Initialize(
  //     int Nsigma_0, 
  //     const LSU3BasisTable& lsu3basis_vector, 
  //     spncci::LGIVectorType& lgi_vector, 
  //     int& Nsigma_begin, std::ofstream& os
  //   )
  // // writing to file the Nex=0 LGI's for which the LGI is give by the lsu3shell irrep so the expansion coefficient is 1
  // // Currently assuming only one of each symmetry at Nex=0...my need to be adjusted later
  // {
  //   //     std::map< spncci::LGI, int> lgi_lsu3shell_map;      
  //   //     int Nex=0;
  //   //     u3::SU3 x; 
  //   //     HalfInt Sp, Sn, S;
  //   //     // for each irrep in lsu3shell basis vector, extract labels and construct 
  //   //     // LGI.  Insert into map with value index of lsu3shell state in lsu3shell basis
  //   //     for(int index=0; index<lsu3basis_vector.size(); index++)
  //   //     {
  //   //       std::tie(Nex, x, Sp, Sn, S)=lsu3basis_vector[index].Key();
  //   //       if(Nex==0)
  //   //         {
  //   //           u3::U3 sigma(Nsigma_0+Nex,x);
  //   //           lgi_lsu3shell_map[spncci::LGI(Nex,sigma,Sp,Sn,S)]=index;
  //   //         }
  //   //       else
  //   //       {
  //   //         Nsigma_begin=index;
  //   //         break;
  //   //       }
  //   //     }
  //   //     // write to file
  //   //     // lsu3_basis_size  number_of_lgis 
  //   //     //  lgi_index  num_nonzero_coefs
  //   //     //    lsu3shell_index coefficient
  //   //     //    lsu3shell_index coefficient
  //   //     //  ... 
  //   //     int lgi_index=0;
  //   //     std::string outstring=fmt::format("{:6d} {:6d}",lsu3basis_vector.size(),lgi_lsu3shell_map.size());
  //   //     os << outstring.c_str()<<std::endl;
  //   //     
  //   //     for(auto it=lgi_lsu3shell_map.begin(); it!=lgi_lsu3shell_map.end(); ++it)
  //   //       {
  //   //         os <<lgi_index<<"  "<<1<<std::endl
  //   //            <<"  "<< it->second <<"  "<<1.0<<std::endl;
  //   //         lgi_vector.push_back(it->first);
  //   //         lgi_index++;
  //   //       }
  // }

  // void WriteLSU3ShellExpansionLGI(
  //     int basis_dim, 
  //     const std::map< spncci::LGI, std::vector<std::pair<int,int>> >& lgi_lsu3shell_map,
  //     const Eigen::MatrixXd& BNcm_kernal,
  //     spncci::LGIVectorType& lgi_vector,
  //     std::string lgi_expansion_filename
  //   )
  // {
  //   int lgi_index;
  //   int num_rows=BNcm_kernal.rows();
  //   int lsu3_start_index;
  //   std::fstream s(lgi_expansion_filename.c_str());
  //   s.seekp(0,std::ios::beg);
  //   s>>lgi_index>>lsu3_start_index;
  //   s.seekp(0,std::ios::end);
  //   int column, count;
  //   for(auto it=lgi_lsu3shell_map.begin(); it!=lgi_lsu3shell_map.end(); ++it)
  //     {
  //       std::vector<std::pair<int,int>> column_count_vector=it->second;
  //       for(int i=0; i<column_count_vector.size(); ++i)
  //         {
  //           std::tie(column,count)=column_count_vector[i];
  //           s<<lgi_index<<"  "<<count<<std::endl;
  //           for(int row=0; row<num_rows; row++)
  //             {
  //               double rme=BNcm_kernal(row,column);
  //               if(fabs(rme)>10e-8)
  //                 s<<"  "<<(lsu3_start_index+row)<<"  "<<rme;
  //             } 
  //           lgi_vector.push_back(it->first);
  //           lgi_index++;
  //         }
  //     }
  //   s.seekp(0,std::ios::beg);
  //   s.write(fmt::format("{:6d} {:6d}",lgi_index,basis_dim).c_str(),13);
  //   s.close();
  // }

  // void WriteControlFile()
  // {

  // }
}// end namespace

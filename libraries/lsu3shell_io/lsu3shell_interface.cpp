/****************************************************************
  sp_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "lsu3shell_interface.h"


#include <fstream>
#include <iostream>
#include <algorithm>
//#include <functional>

#include "cppformat/format.h"

#include "u3shell/two_body_operator.h"
#include "u3shell/moshinsky.h"
#include "spncci/sp_basis.h"
#include "basis/operator.h"
// #include "u3shell/u3spn_scheme.h"

namespace lsu3shell
{
  std::string LSU3Irrep::Str() const
  {
    std::ostringstream ss;

    ss << "[" 
       << " " << Nex_
       << " " << x_.Str()
       << " " << Sp_
       << " " << Sn_
       << " " << S_
       << " " << "]";
    return ss.str();
  }


  void 
  GenerateLSU3ShellOperators(
      int Nmax, 
      const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_tensor_expansion,
      std::string filename
    )
  {
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
    // declare coefficient containers
    u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
    u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
    //moshinsky transform and accumulate coefficients
    u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
    u3shell::TransformRelativeTensorToTwobodyTensor(
      relative_tensor_expansion,
      twobody_space,
      two_body_unit_tensor_coefficients,
      "NAS"   
    );
    // for(auto it=two_body_unit_tensor_coefficients.begin(); it!=two_body_unit_tensor_coefficients.end(); ++it)
    //   std::cout<<it->first.Str()<<"  "<<it->second<<std::endl;
    // convert to biquads
    u3shell::TransformTwoBodyUnitTensorToBiquad(two_body_unit_tensor_coefficients,biquad_coefficients);
    // convert to pn scheme
    u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn);
    
    // std::string operator_stream_filename = fmt::format("operator{:06d}.recoupler",operator_index);
    // std::ofstream operator_stream(operator_stream_filename);
    std::ofstream operator_stream(filename);
    WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
    operator_stream.close();
  }

  void
  GenerateLSU3ShellOperators(
      int Nmax, 
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_tensor_labels
    )
    {
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
    for(int i=0; i<relative_tensor_labels.size(); ++i)
      {
        // declare coefficient containers
        u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
        u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
        u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
        // Moshinsky transform unit tensors to two-body operators 
        MoshinskyTransformUnitTensor(relative_tensor_labels[i], 1.0, twobody_space,two_body_unit_tensor_coefficients, "NAS");
        // convert to biquads
        u3shell::TransformTwoBodyUnitTensorToBiquad(two_body_unit_tensor_coefficients,biquad_coefficients);
        // convert biquads to pn scheme
        u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn);
        // std::cout<<fmt::format("{} operator{:06d}",relative_unit_tensor_labels.Str(),i)<<std::endl;
        std::string operator_stream_filename = fmt::format("relative_unit{:06d}.recoupler",i);
        std::ofstream operator_stream(operator_stream_filename);
        WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
        operator_stream.close();
      }
    }

void GenerateLSU3ShellOperators(
        int Nmax, 
        const u3shell::TwoBodyUnitTensorCoefficientsU3ST& twobody_tensor_expansion,
        int operator_index)
  {
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
    // declare coefficient containers
    u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
    u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
    u3shell::TransformTwoBodyUnitTensorToBiquad(twobody_tensor_expansion,biquad_coefficients);
    // convert to pn scheme
    u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn);
    std::string operator_stream_filename =fmt::format("twobody_unit{:06d}.recoupler",operator_index);
    std::ofstream operator_stream(operator_stream_filename);
    // std::cout<<"hi "<<operator_stream_filename<<std::endl;
    WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
    operator_stream.close();
  }

  void 
  ReadLSU3Vector(
    HalfInt Nsigma_0, 
    const std::string& filename, 
    LSU3BasisTable& lsu3_basis_table, 
    std::map<u3shell::U3SPN,int>& subspace_dimensions
    )
  {
    std::ifstream basis_stream(filename.c_str());
    if(not basis_stream)
      std::cout<<fmt::format("File {} did not open",filename)<<std::endl;
    std::string line;
    while(std::getline(basis_stream,line))
      {
        if(not std::isdigit(line[0]))
          continue;
        std::istringstream line_stream(line);
        // alpha_n_max and alpha_p_max correspond multiplicities arising in the coupling of protons among shells and 
        // neutrons among shells. rho0_max is outer multplicity of coupling proton to neutron. 
        int Np, Nn, Nex,twice_Sp, twice_Sn, twice_S, lambda, mu, ip,in,alpha_p_max,alpha_n_max,rho0_max,lambda_p,mu_p,lambda_n,mu_n;
        line_stream
        >>ip>>alpha_p_max>>Np>>lambda_p>>mu_p>>twice_Sp 
        >>in>>alpha_n_max>>Nn>>lambda_n>>mu_n>>twice_Sn 
        >>rho0_max>>lambda >> mu>>twice_S;

        //conversions
        Nex=Nn+Np;
        HalfInt Sp=HalfInt(twice_Sp,2);
        HalfInt Sn=HalfInt(twice_Sn,2);
        HalfInt S=HalfInt(twice_S,2);
        u3::SU3 x(lambda,mu);
        // std::cout<<fmt::format("Nex {}  Nsigma_0 {}  x {}",Nex, Nsigma_0,x.Str())<<std::endl;
        u3::U3 omega(Nsigma_0+Nex,x);
        u3::U3S omegaS(omega,S);
        u3shell::U3SPN omegaSPN(omegaS,Sp,Sn);

        int start_index=subspace_dimensions[omegaSPN];
        int dim=alpha_n_max*alpha_p_max*rho0_max;
        subspace_dimensions[omegaSPN]+=dim;
        
        LSU3BasisGroup mult_group(omegaSPN,dim,start_index);
        lsu3_basis_table.push_back(mult_group);
      } 
    basis_stream.close();
  }



  void 
  ReadLSU3ShellRMEs(
    std::ifstream& is,
    const u3shell::OperatorLabelsU3S& operator_labels,
    const LSU3BasisTable& lsu3_basis_table,
    const u3shell::SpaceU3SPN& space, 
    const u3shell::SectorsU3SPN& sectors,
    basis::MatrixVector& matrix_vector // in operator.h and initial to zero
  )
  {    
    u3shell::U3SPN omegaSPNi, omegaSPNj;
    int i,j, start_index_i, start_index_j, group_size_i, group_size_j;
    double rme;
    // u3::SU3 x0(operator_labels.x0());
    while(is)
      {
        is>>i,j;
        std::tie(omegaSPNi,group_size_i,start_index_i)=lsu3_basis_table[i];
        std::tie(omegaSPNj,group_size_j,start_index_j)=lsu3_basis_table[j];
        int i_space=space.LookUpSubspaceIndex(omegaSPNi);
        int j_space=space.LookUpSubspaceIndex(omegaSPNj);
        u3::SU3 xi(omegaSPNi.SU3());
        u3::SU3 xj(omegaSPNj.SU3());
        int rho0_max=u3::OuterMultiplicity(xj,operator_labels.x0(),xi);
        for(int gi=0; gi<group_size_i; ++gi)
          for(int gj=0; gj<group_size_j; ++gj)
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                is>>rme;
                int sector_index=sectors.LookUpSectorIndex(i_space,j_space,rho0);
                int row_index=start_index_i+gi;
                int column_index=start_index_j+gj;
                matrix_vector[sector_index](row_index,column_index)=rme;
              }
      }
  }

  void 
  GenerateBrelNrelMatrix(
    const std::string& brel_filename,
    const std::string& nrel_filename,
    const LSU3BasisTable& lsu3_basis_table,
    const u3shell::SpaceU3SPN& space, 
    basis::MatrixVector& matrix_vector 
    ) 
  {    
    u3shell::U3SPN omegaSPNi, omegaSPNj;
    int i,j, start_index_i, start_index_j, group_size_i, group_size_j;
    double rme;

    std::ifstream is_brel(brel_filename.c_str());
    u3shell::OperatorLabelsU3S brel_labels(-2,u3::SU3(0,2),0,0);
    //generate sectors for brel.
    u3shell::SectorsU3SPN sectors(space,brel_labels,true);
    while(is_brel)
      {
        is_brel>>i,j;
        std::tie(omegaSPNi,group_size_i,start_index_i)=lsu3_basis_table[i];
        std::tie(omegaSPNj,group_size_j,start_index_j)=lsu3_basis_table[j];
        int i_space=space.LookUpSubspaceIndex(omegaSPNi);
        int j_space=space.LookUpSubspaceIndex(omegaSPNj);
        int sector_index=sectors.LookUpSectorIndex(i_space,j_space,1);
        // construct matrix for both Brel and Nrel
        matrix_vector[sector_index]=
            Eigen::MatrixXd::Zero(group_size_i+group_size_j,group_size_j);
        u3::SU3 xi(omegaSPNi.SU3());
        u3::SU3 xj(omegaSPNj.SU3());

        assert(u3::OuterMultiplicity(xj,u3::SU3(0,2),xi)==1);
        for(int gi=0; gi<group_size_i; ++gi)
          for(int gj=0; gj<group_size_j; ++gj)
            { 
              is_brel>>rme;
              // write Brel to bottom part of matrix
              int row_index=start_index_j+start_index_i+gi;
              int column_index=start_index_j+gj;
              matrix_vector[sector_index](row_index,column_index)=rme;
            }
      }
    is_brel.close();
    std::ifstream is_nrel(nrel_filename.c_str());
    // u3shell::OperatorLabelsU3S nrel_labels(0,u3::SU3(0,0),0,0);
    HalfInt N=omegaSPNj.N();
    while(is_nrel)
      {
        is_nrel>>i,j;
        std::tie(omegaSPNi,group_size_i,start_index_i)=lsu3_basis_table[i];
        std::tie(omegaSPNj,group_size_j,start_index_j)=lsu3_basis_table[j];
        int i_space=space.LookUpSubspaceIndex(omegaSPNi);
        int j_space=space.LookUpSubspaceIndex(omegaSPNj);
        int sector_index=sectors.LookUpSectorIndex(i_space,j_space,1);

        assert(i_space==j_space);

        u3::SU3 xi(omegaSPNi.SU3());
        u3::SU3 xj(omegaSPNj.SU3());
        double rme_nrel,rme_n=0;
        assert(u3::OuterMultiplicity(xj,u3::SU3(0,0),xi)==1);
        for(int gi=0; gi<group_size_i; ++gi)
          for(int gj=0; gj<group_size_j; ++gj)
            { 
              {
                is_nrel>>rme_nrel;
                int row_index=start_index_i+gi;
                int column_index=start_index_j+gj;
                if(gi==gj)
                  {rme_n=double(N);}
                matrix_vector[sector_index](row_index,column_index)=rme_n-rme_nrel;
              }
            }
      }
    is_nrel.close();
    
  }

  void GenerateNcmMatrix()
  {}

  void GenerateLSU3ShellExpansionLGI(
    const LSU3BasisTable& lsu3_basis_table,
    const u3shell::SpaceU3SPN& space, 
    const std::string& brel_filename,
    const std::string& nrel_filename,
    basis::MatrixVector& lgi_expansion_matrix_vector //size of space 
  )
  // Construct Brel and Ncm matrix in lsu3shell basis and 
  // solve for null space for Nex>2 or 1. 
  
  {
    basis::MatrixVector ncm_matrix_vector;
    GenerateNcmMatrix();
    basis::MatrixVector brel_ncm_matrix_vector;
    GenerateBrelNrelMatrix(brel_filename,nrel_filename,lsu3_basis_table,space, brel_ncm_matrix_vector);
    int num_brel_sectors=brel_ncm_matrix_vector.size();
    int start_index=space.size()-num_brel_sectors;
    for(int i=0; i<=brel_ncm_matrix_vector.size();++i)
      {
        lgi_expansion_matrix_vector[i]=brel_ncm_matrix_vector[i].fullPivLu().kernel();
      }
  }


  void 
  WriteLGI(const spncci::LGIVectorType& lgi_vector,   std::ofstream& os)
  {
    int Nex;
    u3::U3 sigma;
    HalfInt Sp,Sn,S;
    std::unordered_map<spncci::LGI,int,boost::hash<spncci::LGI>> lgi_counter;
    for(auto a:lgi_vector)
      {
        lgi_counter[a]+=1;
      }
    for(auto b:lgi_vector)
      {
        std::tie(Nex,sigma,Sp,Sn,S)=b.Key();
        int count=lgi_counter[b];
        os
        <<Nex<<"  "<<TwiceValue(Sp)<<"  "<<TwiceValue(Sn)<<"  "<<TwiceValue(S)
        <<"  "<<sigma.SU3().lambda()<<"  "<<sigma.SU3().mu()<<"  "<<count<<std::endl;     
      }
  }

  void 
  LGINex0Initialize(
      int Nsigma_0, 
      const LSU3Vector& lsu3basis_vector, 
      spncci::LGIVectorType& lgi_vector, 
      int& Nsigma_begin, std::ofstream& os
    )
  // writing to file the Nex=0 LGI's for which the LGI is give by the lsu3shell irrep so the expansion coefficient is 1
  // Currently assuming only one of each symmetry at Nex=0...my need to be adjusted later
  {
//     std::map< spncci::LGI, int> lgi_lsu3shell_map;      
//     int Nex=0;
//     u3::SU3 x; 
//     HalfInt Sp, Sn, S;
//     // for each irrep in lsu3shell basis vector, extract labels and construct 
//     // LGI.  Insert into map with value index of lsu3shell state in lsu3shell basis
//     for(int index=0; index<lsu3basis_vector.size(); index++)
//     {
//       std::tie(Nex, x, Sp, Sn, S)=lsu3basis_vector[index].Key();
//       if(Nex==0)
//         {
//           u3::U3 sigma(Nsigma_0+Nex,x);
//           lgi_lsu3shell_map[spncci::LGI(Nex,sigma,Sp,Sn,S)]=index;
//         }
//       else
//       {
//         Nsigma_begin=index;
//         break;
//       }
//     }
//     // write to file
//     // lsu3_basis_size  number_of_lgis 
//     //  lgi_index  num_nonzero_coefs
//     //    lsu3shell_index coefficient
//     //    lsu3shell_index coefficient
//     //  ... 
//     int lgi_index=0;
//     std::string outstring=fmt::format("{:6d} {:6d}",lsu3basis_vector.size(),lgi_lsu3shell_map.size());
//     os << outstring.c_str()<<std::endl;
//     
//     for(auto it=lgi_lsu3shell_map.begin(); it!=lgi_lsu3shell_map.end(); ++it)
//       {
//         os <<lgi_index<<"  "<<1<<std::endl
//            <<"  "<< it->second <<"  "<<1.0<<std::endl;
//         lgi_vector.push_back(it->first);
//         lgi_index++;
//       }
   }

  void WriteLSU3ShellExpansionLGI(
        int basis_dim, 
        const std::map< spncci::LGI, std::vector<std::pair<int,int>> >& lgi_lsu3shell_map,
        const Eigen::MatrixXd& BNcm_kernal,
        spncci::LGIVectorType& lgi_vector,
        std::string lgi_expansion_filename
      )
    {
      int lgi_index;
      int num_rows=BNcm_kernal.rows();
      int lsu3_start_index;
      std::fstream s(lgi_expansion_filename.c_str());
      s.seekp(0,std::ios::beg);
      s>>lgi_index>>lsu3_start_index;
      s.seekp(0,std::ios::end);
        int column, count;
      for(auto it=lgi_lsu3shell_map.begin(); it!=lgi_lsu3shell_map.end(); ++it)
        {
          std::vector<std::pair<int,int>> column_count_vector=it->second;
          for(int i=0; i<column_count_vector.size(); ++i)
          {
            std::tie(column,count)=column_count_vector[i];
            s<<lgi_index<<"  "<<count<<std::endl;
            for(int row=0; row<num_rows; row++)
              {
                double rme=BNcm_kernal(row,column);
                if(fabs(rme)>10e-8)
                    s<<"  "<<(lsu3_start_index+row)<<"  "<<rme;
              } 
            lgi_vector.push_back(it->first);
            lgi_index++;
          }
        }
      s.seekp(0,std::ios::beg);
      s.write(fmt::format("{:6d} {:6d}",lgi_index,basis_dim).c_str(),13);
      s.close();
    }

  void WriteControlFile()
  {

  }
}// end namespace

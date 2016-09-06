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

  void ReadLSU3Vector(const std::string& filename, LSU3Vector& lsu3basis_vector)
  {
    std::ifstream basis_stream(filename.c_str());
    std::string line;
    while(std::getline(basis_stream,line))
      {
        std::istringstream line_stream(line);
        // alpha_n_max and alpha_p_max correspond multiplicities arising in the coupling of protons among shells and 
        // neutrons among shells. rho0_max is outer multplicity of coupling proton to neutron. 
        int Nex, twice_Sp, twice_Sn, twice_S, lambda, mu, ip,in,alpha_p_max,alpha_n_max,rho0_max,lambda_p,mu_p,lambda_n,mu_n;
        line_stream 
        >> Nex 
        >>ip>>alpha_p_max >>lambda_p>>mu_p>>twice_Sp 
        >>in>>alpha_n_max>>lambda_n>>mu_n>>twice_Sn 
        >>rho0_max>>lambda >> mu>>twice_S;

        //conversions
        HalfInt Sp=HalfInt(twice_Sp,2);
        HalfInt Sn=HalfInt(twice_Sn,2);
        HalfInt S=HalfInt(twice_S,2);
        u3::SU3 x(lambda,mu);
        int dim=alpha_n_max*alpha_p_max*rho0_max;
        for(int i=0; i<dim; ++i)
          lsu3basis_vector.push_back(lsu3shell::LSU3Irrep(Nex,x,Sp,Sn,S));
      } 
    basis_stream.close();
  }

  void GenerateLSU3ShellOperators(int Nmax, const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_tensor_expansion, int operator_index)
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
    
    std::string operator_stream_filename = fmt::format("operator{:06d}.recoupler",operator_index);
    std::ofstream operator_stream(operator_stream_filename);
    WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
    operator_stream.close();
  }

  void GenerateLSU3ShellOperators(int Nmax, const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_tensor_labels)
    {
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
    for(int i=0; i<relative_tensor_labels.size(); ++i)
      {
        // declare coefficient containers
        u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
        u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
        u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
        // Moshinsky transform unit tensors to two-body operators 
        MoshinskyTransformUnitTensor(relative_tensor_labels[i], 1.0, twobody_space,two_body_unit_tensor_coefficients, "AS");
        // convert to biquads
        u3shell::TransformTwoBodyUnitTensorToBiquad(two_body_unit_tensor_coefficients,biquad_coefficients);
        // convert biquads to pn scheme
        u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn);
        // std::cout<<fmt::format("{} operator{:06d}",relative_unit_tensor_labels.Str(),i)<<std::endl;
        std::string operator_stream_filename = fmt::format("operator{:06d}.recoupler",i);
        std::ofstream operator_stream(operator_stream_filename);
        WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
        operator_stream.close();
      }
    }

void GenerateLSU3ShellOperators(
        int Nmax, 
        const u3shell::TwoBodyUnitTensorCoefficientsU3ST& twobody_tensor_expansion, 
        int operator_index
      )
  {
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
    // declare coefficient containers
    u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
    u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
    u3shell::TransformTwoBodyUnitTensorToBiquad(twobody_tensor_expansion,biquad_coefficients);
    // convert to pn scheme
    u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn);
    std::string operator_stream_filename = fmt::format("operator{:06d}.recoupler",operator_index);
    std::ofstream operator_stream(operator_stream_filename);
    WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
    operator_stream.close();
  }


  void ReadLSU3ShellRMEs(
        std::ifstream& is,
        int dimp, int dim, const u3::SU3& x0,
        const lsu3shell::LSU3Vector& lsu3shell_vector_bra,
        const lsu3shell::LSU3Vector& lsu3shell_vector_ket,
        std::vector<Eigen::MatrixXd>& matrix_vector
      )
  {
    u3::SU3 x, xp;
    for(int i=0; i<dimp; ++i)
      for(int j=0; j<dim; ++j)
        {
          xp=lsu3shell_vector_bra[i].x();
          x=lsu3shell_vector_ket[j].x();
          int rho0_max=u3::OuterMultiplicity(x,x0,xp);
          for(int rho0=1; rho0<=rho0_max; ++rho0)
            {
              if (rho0<matrix_vector.size())
                matrix_vector.push_back(Eigen::MatrixXd::Zero(dimp,dim)); 
              is>>matrix_vector[rho0-1](i,j);
            }
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
  GenerateBrelNrelMatrix(
      int N_sigma,
      const std::string& brel_filename,
      const std::string& nrel_filename,
      int& dim0, int& dim2, 
      Eigen::MatrixXd& BNcmMatrix
    )
  {   
    std::ifstream is_brel(brel_filename.c_str());
    is_brel>>dim0>>dim2;
    // if((dim0+dim2)!=dim)
    //  std::cout<<"brel dimension mismatch"<<std::endl;
    BNcmMatrix=Eigen::MatrixXd::Zero(dim0+dim2,dim2);

    //Read in matrix elements of Brel
    double rme;
    for(int i=0; i<dim0; ++i)
      for(int j=0; j<dim2; ++j)
        {
          is_brel>>rme;
          BNcmMatrix(i,j)=rme;
        }
    is_brel.close();
    //Read in matrix elements of Nrel
    int dim01,dim02;
    std::ifstream is_nrel(nrel_filename.c_str());
    is_nrel>>dim01>>dim02;

    //Sanity check
    if((dim01!=dim02)||(dim01!=dim0))
      std::cout<<"dimension miss match between Brel and Nrel"<<std::endl;

    for(int i=0; i<dim2; ++i)
      for(int j=0; j<dim2; ++j)
        {
          is_nrel>>rme;
          //Ncm=Ntotal-Nre;;
          BNcmMatrix(dim0+i,j)=N_sigma-rme;
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


  void GenerateLSU3ShellExpansionLGI(
    int Nsigma_0,
    int Nsigma_min,
    int Nsigma_max, 
    std::string basis_file, 
    std::string brel_filename, 
    std::string nrel_filename,
    std::string lgi_filename,
    std::string lgi_expansion_filename
  )
  // Construct Brel and Ncm matrix in lsu3shell basis and 
  // solve for null space 
  {
    int Nex;
    u3::SU3 x; 
    HalfInt Sp, Sn, S;
    spncci::LGIVectorType lgi_vector;
    // Nsigma_begin keeps track of where in basis vector irreps with Nex=Nsigma begin
    int Nsigma_begin;
    // read in lsu3shell reduced basis states 
    LSU3Vector lsu3basis_vector;
    lsu3shell::ReadLSU3Vector(basis_file, lsu3basis_vector);
    int basis_dim=lsu3basis_vector.size();
    //if starting at Nsigma=0, use LSU3 irreps to fill out LGIVector for Nsigma=0;
    if(Nsigma_min==0)
      {
        std::ofstream os(lgi_expansion_filename.c_str());
        LGINex0Initialize(Nsigma_0, lsu3basis_vector,lgi_vector, Nsigma_begin,os);
        os.close(); 
        Nsigma_min=2;
      }
    // Otherwise read in already generated LGI's from file and locate starting position
    // in lsu3basis_vector  
    else
      {
        spncci::ReadLGISet(lgi_vector,lgi_filename, Nsigma_0);
        for(int a=0; a<lsu3basis_vector.size(); ++a)
            if(lsu3basis_vector[a].Nex()==Nsigma_min)
              {
                Nsigma_begin=a;
                continue;
              }
      }

    // For each Nsigma
    for(int Nsigma=Nsigma_min; Nsigma<=Nsigma_max; Nsigma+=2)
    {
      /////////////////////////////////////////////////////////////////////////////////
      // Construct the Brel+Nrel matrix and solve for null space 
      /////////////////////////////////////////////////////////////////////////////////
      // Assumes that Brel and Nrel rme's have already been calculated using SU3RME 
      // for each Nsigma
      // dim2 is the dimension of the Nsigma space
      // dim0 is the dimension of the Nsigma-2 space
      // dim2 and dim0 are read in from file containing Brel operator rme's
      int dim2,dim0;
      Eigen::MatrixXd BNcmMatrix;
      GenerateBrelNrelMatrix(Nsigma,brel_filename, nrel_filename,dim0, dim2, BNcmMatrix);
      //Solve for null space 
      Eigen::MatrixXd BNcm_kernal=BNcmMatrix.fullPivLu().kernel();
      ///////////////////////////////////////////////////////////////////////////////
      //Get list of LGI labels and append to lgi structure 
      // map of LGI's with value vector of pairs of corresponding column in Kernel matrix 
      // and number of non-zero rme's 
      std::map< spncci::LGI, std::vector<std::pair<int,int>> > lgi_lsu3shell_map;     
      int num_column=BNcm_kernal.cols();
      int num_rows=BNcm_kernal.rows();
      // accumulate lgi's in map
      spncci::LGI lgi;
      for(int column=0; column<num_column; column++)
        {
          int count=0;
          for(int row=0; row<num_rows; row++)
            {
              double rme=BNcm_kernal(row,column);
              if(fabs(rme)>10e-8)
                {
                  // Position in lsu3basis_vector is give by the sum of Nsigma_begin, 
                  // which keeps track of where in the vector the Nex=Nsigma subset 
                  // of the lsu3shell basis begin, and row, which indexes which of the
                  // states within subset of the basis.  
                  std::tie(Nex,x,Sp,Sn,S)=lsu3basis_vector[Nsigma_begin+row].Key();
                  u3::U3 sigma (Nex+Nsigma_0,x);
                  //lgi=spncci::LGI(Nex,sigma,Sp,Sn,S);
                  lgi=spncci::LGI(u3shell::U3SPN(u3::U3S(sigma,S),Sp,Sn),Nex);
                  count++;
                }
            }

          // LGI with same labels are collected in a vector
          lgi_lsu3shell_map[lgi].push_back(std::pair<int,int>(column,count));
        }
    //Write expansion of lgi in lsu3shell basis to file, appending to lgi_expansion_filename
    //and append LGI's to lgi_vector in canonical order   
    WriteLSU3ShellExpansionLGI(basis_dim, lgi_lsu3shell_map,BNcm_kernal, lgi_vector, lgi_expansion_filename);     
    }
    //write LGI's to file
    std::ofstream os(lgi_filename.c_str());
    WriteLGI(lgi_vector,os);
  }

  void WriteControlFile()
  {

  }
}// end namespace

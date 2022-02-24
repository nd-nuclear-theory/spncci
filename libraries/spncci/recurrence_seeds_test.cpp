/****************************************************************
  recurrence_seeds_test.cpp

  Anna E. McCoy
  INT

  SPDX-License-Identifier: MIT

 9/24/21 (aem): Created.
****************************************************************/
#include "spncci/recurrence_seeds.h"

#include <cppitertools/itertools.hpp>
#include <string>
#include <numeric>
#include <iostream>
#include <fstream>

#include "am/halfint_fmt.h"
#include "fmt/format.h"
#include "lgi/lgi.h"

#include "spncci_basis/recurrence_indexing.h"
#include "lgi/recurrence_lgi.h"

#include "lgi/lgi_unit_tensors.h"
#include "mcutils/profiling.h"
#include "mcutils/eigen.h"
#include "u3shell/relative_operator.h"
#include "utilities/utilities.h"

#include "Spectra/SymEigsSolver.h"

unsigned long int NumNonZeros(const basis::OperatorBlock<double>& block)
  {
    unsigned long int counter=0;
    for(int i=0; i<block.rows(); ++i)
      for(int j=0; j<block.cols(); ++j)
          if(std::abs(block(i,j))>1e-6)
            counter++;

    return counter;
  }


basis::OperatorBlock<double> get_temp_block(int Z, int N, lgi::LGI& bra, lgi::LGI& ket, int rho0,int operator_index)
{
  const auto&[Nex_bra,sigma_bra,Sp_bra,Sn_bra,S_bra]=bra.Key();
  const auto&[Nex_ket,sigma_ket,Sp_ket,Sn_ket,S_ket]=ket.Key();

  std::string filename
    = fmt::format("temp/seeds_Z{:02d}_N{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_2Sp{:02d}_2Sn{:02d}_2S{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_2Sp{:02d}_2Sn{:02d}_2S{:02d}_rho{}_i{}.dat",
          Z,N,Nex_bra,sigma_bra.SU3().lambda(),sigma_bra.SU3().mu(),TwiceValue(Sp_bra),TwiceValue(Sn_bra),TwiceValue(S_bra),
          Nex_ket,sigma_ket.SU3().lambda(),sigma_ket.SU3().mu(),TwiceValue(Sp_ket),TwiceValue(Sn_ket),TwiceValue(S_ket),
          rho0,operator_index
        );

  return ReadOperatorBlockBinary(filename);

}


int main(int argc, char** argv)
{
  if (argc<5)
  {
    std::cout<<"Usage: recurrence seeds_test <Z> <N> <Nmax> <N1v>"<<std::endl;
  }

  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  int N1v=std::stoi(argv[4]);

  const nuclide::NuclideType nuclide({Z,N});
  const std::string lgi_filename="lgi_families.dat";
  
  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide,true);
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  lgi::ReadLGISet(lgi_filename,Nsigma0,lgi_vector);

  // std::vector<int> lgi_full_space_index_lookup;
  // lgi::ReadLGILookUpTable(lgi_full_space_index_lookup, lgi_vector.size());

  std::vector<int> lgi_full_space_index_lookup(lgi_vector.size());
  for(int i=0; i<lgi_vector.size(); ++i)
    lgi_full_space_index_lookup[i]=i;
  
  spncci::spin::Space<lgi::LGI> spin_space(lgi_vector, Nmax);
  const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST> spin_recurrence_space(spin_space, spin_space);

  auto it =iter::imap([](MultiplicityTagged<lgi::LGI> l) { return l.irrep.U3(); }, lgi_vector)| iter::unique_everseen;
  const spncci::spatial::Space spatial_space(std::vector<u3::U3>(it.begin(), it.end()), Nsigma0, Nmax);
  const spncci::spatial::RecurrenceSpace spatial_recurrence_space(spatial_space,spatial_space,N1v,Nsigma0);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Dump spatial basis 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(false)
  {
    for(const auto& sp3r_space : spatial_recurrence_space)
      {
        std::cout<<sp3r_space.LabelStr()<<std::endl;
        for(const auto & Nsum_space : sp3r_space)
          {
            std::cout<<".."<<Nsum_space.LabelStr()<<std::endl;
            for(int u3_space_index=0; u3_space_index<Nsum_space.size(); ++u3_space_index)
            
              {
                const auto& u3_space=Nsum_space.GetSubspace(u3_space_index);
                std::cout<<fmt::format("....{}  [{}]",u3_space.LabelStr(),
                  Nsum_space.GetSubspaceDegeneracy(u3_space_index)
                  )<<std::endl;
                for(int operator_subspace_index=0; operator_subspace_index<u3_space.size(); ++operator_subspace_index)
                  { 
                    const auto& operator_subspace=u3_space.GetSubspace(operator_subspace_index);
                    std::cout<<fmt::format("......{}  [{}]",
                      operator_subspace.LabelStr(),u3_space.GetSubspaceDegeneracy(operator_subspace_index)
                      )<<std::endl;
                    for(const auto& operator_state : operator_subspace)
                      std::cout<<fmt::format("........{}",operator_state.LabelStr())<<std::endl;
                  }
              }
          }
      }
  
    std::cout<<"---------------------------------------------------"<<std::endl<<std::endl;
  }
  if(false)
  {
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Dump spin basis 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(const auto& lgi_space : spin_recurrence_space)
    {
      std::cout<<lgi_space.LabelStr()<<std::endl;
      for(const auto& recurrence_spin_space : lgi_space)
        {
          std::cout<<fmt::format("..{}",recurrence_spin_space.LabelStr())<<std::endl;
          for(int subspace_index=0; subspace_index<recurrence_spin_space.size(); ++subspace_index)
            {
              const auto& spin_subspace=recurrence_spin_space.GetSubspace(subspace_index);
              const auto& [Sp_bra,Sn_bra]=spin_subspace.bra_upstream_labels();
              const auto& [Sp_ket,Sn_ket]=spin_subspace.ket_upstream_labels();
              std::cout<<fmt::format("....({} {}) ({} {}) [{}]",
                Sp_ket,Sn_ket,Sp_bra,Sn_bra,recurrence_spin_space.GetSubspaceDegeneracy(subspace_index)
                )<<std::endl;

              for(const auto& state : spin_subspace)
                {
                  std::cout<<fmt::format("......{}",state.labels().LabelStr())<<std::endl;
                }
            }    
        }
      
    }
    std::cout<<"---------------------------------------------------"<<std::endl<<std::endl;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Dump unit tensor rme info
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(true)
  {

    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
    u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1v,relative_unit_tensor_labels);

  
  for(int irrep_family_index_bra=0; irrep_family_index_bra<lgi_vector.size(); ++irrep_family_index_bra)
    for(int irrep_family_index_ket=0; irrep_family_index_ket<lgi_vector.size(); ++irrep_family_index_ket)
      {
        std::cout<<fmt::format("{}  {}",
          lgi_vector[irrep_family_index_ket].Str(),lgi_vector[irrep_family_index_bra].Str()
          )<<std::endl;

        int index1=lgi_full_space_index_lookup[irrep_family_index_bra];
        int index2=lgi_full_space_index_lookup[irrep_family_index_ket];

        // Read in operators 

        std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensors;
        std::vector<int> rho0_values;
        std::string lgi_unit_tensor_filename=fmt::format("seeds/operators_{:06d}_{:06d}.dat",index1,index2);
        lgi::ReadUnitTensorLabels(lgi_unit_tensor_filename,lgi_unit_tensors,rho0_values);

        // Reads in unit tensor seed blocks and stores them in a vector of blocks. Order
        // corresponds to order of (unit_tensor,rho0) pairs in corresponding operator file. 
        basis::OperatorBlocks<double> unit_tensor_seed_blocks;
        std::string seed_filename=fmt::format("seeds/seeds_{:06d}_{:06d}.rmes",index1,index2);
        lgi::ReadBlocks(seed_filename, lgi_unit_tensors.size(), unit_tensor_seed_blocks);

        for(int t=0; t<lgi_unit_tensors.size(); ++t)
          { 
            //Select source block
            const basis::OperatorBlock<double>& block = unit_tensor_seed_blocks[t];
            // if(irrep_family_index_ket==irrep_family_index_bra)
            {
              int operator_index=-1;
              for(int ii=0; ii<relative_unit_tensor_labels.size(); ++ii)
                {
                  if(relative_unit_tensor_labels[ii]==lgi_unit_tensors[t])
                    {
                      operator_index=ii;
                      break;
                    }
                }

              basis::OperatorBlock<double> temp_block = get_temp_block(
                  Z,N,
                  lgi_vector[irrep_family_index_bra].irrep,
                  lgi_vector[irrep_family_index_ket].irrep,
                  rho0_values[t],operator_index
                );

              assert(block.cols()==temp_block.cols());
              assert(block.rows()==temp_block.rows());

              Eigen::JacobiSVD<basis::OperatorBlock<double>> svd1(block,Eigen::ComputeFullU);
              // svd1.setThreshold(1e-8);
              auto Svalues1=svd1.singularValues();

              Eigen::JacobiSVD<basis::OperatorBlock<double>> svd2(temp_block,Eigen::ComputeFullU);
              // svd2.setThreshold(1e-8);
              auto Svalues2=svd2.singularValues();
              if(!mcutils::IsZero(Svalues1-Svalues2,1e-6))
              {
                std::cout<<fmt::format("{}  [{}]",lgi_unit_tensors[t].Str(),rho0_values[t])<<std::endl;
                std::cout<<Svalues1<<std::endl<<std::endl<<Svalues2<<std::endl;
              }

              assert(mcutils::IsZero(Svalues1-Svalues2,1e-5));
              // std::cout<<"*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*"<<std::endl;


            }
            // std::cout<<block<<std::endl<<std::endl;
          }
        // std::cout<<"******"<<std::endl<<std::endl;
      }
  std::cout<<"---------------------------------------------------"<<std::endl<<std::endl;
  }
  if(false)
  {
    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
    u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1v,relative_unit_tensor_labels);
    for(const auto& tensor : relative_unit_tensor_labels)
      std::cout<<tensor.Str()<<std::endl;
    std::cout<<"---------------------------------------------------"<<std::endl<<std::endl;
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

  // std::cout<<"Get the seeds "<<std::endl;
  basis::OperatorBlocks<double> seed_blocks
    =spncci::recurrence::GetRecurrenceSeedsFromFile(
      lgi_vector,lgi_full_space_index_lookup,
      spatial_recurrence_space,spin_recurrence_space
      );


  std::ofstream outfile;
  outfile.open ("output.dat");
  

  unsigned long int total_rmes=0;
  unsigned long int total_non_zero=0;
  int num_zero_blocks=0;
  unsigned long int num_zeros=0;
  int num_zero_rows=0;
  int num_rows=0;
  unsigned long int num_zero_row_rmes=0;
  std::vector<float> fractions(seed_blocks.size());
  for(int b=0; b<seed_blocks.size(); ++b)
    {
      const auto& block = seed_blocks[b];
      auto num_rmes=block.rows()*block.cols();
      auto num_non_zero=NumNonZeros(block);

      for(int i=0; i<block.rows(); ++i)
        {
          int row_non_zeros=NumNonZeros(block.block(i,0,1,block.cols()));
          if(row_non_zeros==0)
          {
            num_zero_rows++;
            num_zero_row_rmes+=block.cols();
          }
        }

      num_rows+=block.rows();
      auto ratio=float(num_non_zero)/num_rmes;
      fractions[b]=ratio;
      total_rmes+=num_rmes;
      total_non_zero+=num_non_zero;
      if (num_non_zero==0)
        {
          num_zero_blocks++;
          num_zeros+=block.size();
        }

      if(true)
        outfile<<fmt::format("Num rmes: {}  Num non-zeros: {}  Percentage: {:.2f}",num_rmes,num_non_zero,ratio)<<std::endl;
    }
  outfile<<"----------------------------------------------------------------"<<std::endl;
  outfile<<fmt::format("Total rmes: {}   Total non-zeros: {}  Percentage: {:.2f}",
    total_rmes,total_non_zero,float(total_non_zero)/total_rmes)<<std::endl;

  outfile<<fmt::format("Num zero blocks: {}.  Total number of blocks: {}",num_zero_blocks,seed_blocks.size())<<std::endl;
  outfile<<fmt::format("Num zero rmes from zero blocks: {}   Total zero rmes: {}  Percentage: {:.2f}",
    num_zeros,total_rmes-total_non_zero, num_zeros/float(total_rmes-total_non_zero))<<std::endl;
  outfile<<fmt::format("Num rows of all zeros:  {}. Num rows:  {}. Num rmes:  {}  Percentage:  {:.2f}",
    num_zero_rows,num_rows,num_zero_row_rmes,float(num_zero_row_rmes)/total_rmes)<<std::endl;

  outfile.close();


  for(int index=0; index<spatial_recurrence_space.size(); ++index)
    {
      const auto& spatial_subspace = spatial_recurrence_space.GetSubspace(index);
      const auto&[sigma_ket,sigma_bra,parity_bar] = spatial_subspace.labels();
      std::string seed_filename = spncci::seeds::seed_filename(Z,N,Nsigma0,sigma_bra,sigma_ket,parity_bar);
      basis::OperatorBlock<double> seed_block1 = ReadOperatorBlockBinary(seed_filename);
      basis::OperatorBlock<double> seed_block2 = seed_blocks[index];
      std::cout<<"****************************************"<<std::endl;
      std::cout<<sigma_bra.Str()<<"  "<<sigma_ket.Str()<<"  "<<parity_bar<<std::endl;
      if(sigma_bra.N()==Nsigma0 || sigma_ket.N()==Nsigma0)
        {
          std::cout<<mcutils::FormatMatrix(seed_block1,"+.2f")<<std::endl<<std::endl;
          std::cout<<mcutils::FormatMatrix(seed_block2,"+.2f")<<std::endl<<std::endl;
        }
      if(mcutils::IsZero(seed_block1-seed_block2,1e-6))
        std::cout<<"blocks matched"<<std::endl;
      else
        std::cout<<"blocks don't match"<<std::endl;
    }

}

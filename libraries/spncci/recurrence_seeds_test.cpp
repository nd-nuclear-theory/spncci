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
#include "lgi/lgi_unit_tensors.h"
#include "mcutils/profiling.h"
#include "spncci/recurrence_indexing.h"
#include "u3shell/relative_operator.h"


unsigned long int NumNonZeros(const basis::OperatorBlock<double>& block)
  {
    unsigned long int counter=0;
    for(int i=0; i<block.rows(); ++i)
      for(int j=0; j<block.cols(); ++j)
          if(std::abs(block(i,j))>1e-6)
            counter++;

    return counter;
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

  std::vector<int> lgi_full_space_index_lookup;
  lgi::ReadLGILookUpTable(lgi_full_space_index_lookup, lgi_vector.size());
  
  spncci::spin::Space<lgi::LGI> spin_space(lgi_vector, Nmax);
  const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST> spin_recurrence_space(spin_space, spin_space);

  auto it =iter::imap([](MultiplicityTagged<lgi::LGI> l) { return l.irrep.U3(); }, lgi_vector)| iter::unique_everseen;
  const spncci::spatial::Space spatial_space(std::vector<u3::U3>(it.begin(), it.end()), Nsigma0, Nmax);
  const spncci::spatial::RecurrenceSpace spatial_recurrence_space(spatial_space,spatial_space,N1v,Nsigma0);
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Dump spatial basis 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(true)
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
  if(true)
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
            std::cout<<fmt::format("{}  [{}]",lgi_unit_tensors[t].Str(),rho0_values[t])<<std::endl;

            //Select source block 
            const auto& block = unit_tensor_seed_blocks[t];
            std::cout<<block<<std::endl<<std::endl;
          }
        std::cout<<"******"<<std::endl<<std::endl;
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
}

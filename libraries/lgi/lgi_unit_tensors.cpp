/****************************************************************
  lgi_unit_tensors.cpp

  Anna E. McCoy
  TRIUMF
  
  SPDX-License-Identifier: MIT
****************************************************************/
#include "lgi/lgi_unit_tensors.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <omp.h>  

#include "fmt/format.h"
// #include "lsu3shell/lsu3shell_rme.h"
#include "lgi/lgi_solver.h"
#include "mcutils/io.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h" //For ReadUnitTensorLabels



namespace lgi
{


  // zero tolerance 
  double zero_tolerance=1e-6;
  
  // output mode
  int binary_format_code = 1;
  int binary_float_precision=8;


  void
  TransformOperatorToSpBasis(
      const u3shell::SectorsU3SPN& sectors,
      const basis::OperatorBlocks<double>& basis_transformation_matrices,
      const basis::OperatorBlocks<double>& lsu3shell_operator_matrices,
      basis::OperatorBlocks<double>& spncci_operator_matrices
    )
  {
    // for each sector, look up bra and ket subspaces 
    spncci_operator_matrices.resize(lsu3shell_operator_matrices.size());
    
    // #pragma omp parallel for schedule(runtime)
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
  TransformOperatorToSpBasis(
    const u3shell::SpaceU3SPN& space_bra, 
    const u3shell::SpaceU3SPN& space_ket, 
    const u3shell::SectorsU3SPN& sectors,
    const basis::OperatorBlocks<double>& basis_transformation_matrices,
    const basis::OperatorBlocks<double>& lsu3shell_operator_matrices,
    std::unordered_map<u3shell::U3SPN,int,boost::hash<u3shell::U3SPN>> lgi_lookup_table,
    basis::OperatorBlocks<double>& spncci_operator_matrices
  )
  {
    // for each sector, look up bra and ket subspaces 
    spncci_operator_matrices.resize(lsu3shell_operator_matrices.size());
    
    // #pragma omp parallel for schedule(runtime)
    for(int s=0; s<lsu3shell_operator_matrices.size(); ++s)
      {
        int bra_index=sectors.GetSector(s).bra_subspace_index();
        int ket_index=sectors.GetSector(s).ket_subspace_index();
        const u3shell::SubspaceU3SPN& subspace_bra=space_bra.GetSubspace(bra_index);
        const u3shell::SubspaceU3SPN& subspace_ket=space_ket.GetSubspace(ket_index);

        int i=lgi_lookup_table[subspace_bra.U3SPN()];
        int j=lgi_lookup_table[subspace_ket.U3SPN()];

        // get transformation matrices and transpose bra transformation matrix
        const Eigen::MatrixXd& bra=basis_transformation_matrices[i].transpose();
        const Eigen::MatrixXd& ket=basis_transformation_matrices[j];

        // transform operator to spncci basis
        spncci_operator_matrices[s]=bra*lsu3shell_operator_matrices[s]*ket;
      }
  }



  void RegroupSeedBlocks(
      int unit_tensor_index,
      const u3shell::SectorsU3SPN& unit_tensor_sectors,
      const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels,
      lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels,
      basis::OperatorBlocks<double>& unit_tensor_spncci_matrices,
      std::vector<int>& lsu3shell_index_lookup_table,
      bool restrict_seeds
    )
  {    
    for(int sector_index=0; sector_index<unit_tensor_sectors.size(); ++sector_index)
      {
        // Check that the sector has non-zero rmes as defined by the zero_tolerance
        if(mcutils::IsZero(unit_tensor_spncci_matrices[sector_index],lgi::zero_tolerance))
          continue;
            
        // extract U3SPN sector information from unit tensor sectors definded in lsu3shell basis
        const typename u3shell::SectorsU3SPN::SectorType& sector = unit_tensor_sectors.GetSector(sector_index);
        int bra_subspace_index = sector.bra_subspace_index();
        int ket_subspace_index = sector.ket_subspace_index();

        int bra_lgi_index=lsu3shell_index_lookup_table[bra_subspace_index];
        int ket_lgi_index=lsu3shell_index_lookup_table[ket_subspace_index];

        // If restrict_seeds=true, keep only if bra>=ket, else keep all seeds
        bool keep=restrict_seeds?(bra_lgi_index>=ket_lgi_index):true;
        if(not keep)
          continue;

        const int rho0=sector.multiplicity_index();

        // Regroup by lgi pair, then by unit tensor
        std::pair<int,int> irrep_family_pair(bra_lgi_index,ket_lgi_index);
            
        // Labels for looking up correct block to write to file
        lgi::SeedLabels seed_labels(unit_tensor_labels,rho0,unit_tensor_index,sector_index);

        // Accumulating list 
        lgi_grouped_seed_labels[irrep_family_pair].push_back(seed_labels);
          
      }// sector index
  }



//TODO: allow for different bra and ket spaces 

  void ComputeUnitTensorSeedBlocks(
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const std::string& relative_unit_tensor_filename_template,
      const u3shell::SpaceU3SPN& lsu3shell_space, 
      const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
      const basis::OperatorBlocks<double>& lgi_expansions,
      lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels,
      std::vector<basis::OperatorBlocks<double>>& unit_tensor_spncci_matrices_array,
      std::vector<int>& lsu3shell_index_lookup_table,
      bool restrict_seeds
    )
  {
    // Compute unit tensor seed blocks from lsu3shell unit tensor rmes and lgi expansion 
    // and write blocks to file in binary format
    //
    // Container for seed labels and indices for lookup when writing to file
    // lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels

    unit_tensor_spncci_matrices_array.resize(lgi_unit_tensor_labels.size());

    #pragma omp parallel
      {
        lgi::LGIGroupedSeedLabels lgi_grouped_seed_labels_temp;
       
        // For each unit tensor 
        #pragma omp for schedule(dynamic) nowait
        for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
          {
            // Get labels of corresponding unit tensor 
            basis::OperatorBlocks<double> unit_tensor_spncci_matrices;
            const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];

            // Get file name containing lsu3shell rmes of unit tensor 
            std::string filename = fmt::format(relative_unit_tensor_filename_template,unit_tensor_index);

            // generate unit tensor sectors in lsu3shell basis
            const bool spin_scalar = false; // All spin values allowed
            u3shell::SectorsU3SPN unit_tensor_sectors(lsu3shell_space,unit_tensor_labels,spin_scalar);

            // std::cout<<"Read in lsu3shell rmes of unit tensor"<<std::endl;
            basis::OperatorBlocks<double> unit_tensor_lsu3shell_matrices;
            lsu3shell::ReadLSU3ShellRMEs(
                filename,
                lsu3shell_basis_table,lsu3shell_space,
                unit_tensor_labels,unit_tensor_sectors,
                unit_tensor_lsu3shell_matrices
              );

            // std::cout<<"Transform seed rmes to SpNCCI basis "<<std::endl;
            lgi::TransformOperatorToSpBasis(
                unit_tensor_sectors,lgi_expansions,
                unit_tensor_lsu3shell_matrices,
                unit_tensor_spncci_matrices
              );

///TODO: Check if RegroupSeedBlocks requires modification
            // std::cout<<"Re-organize unit tensor seed bocks "<<std::endl;
            // by lgi then by tensor
            lgi::RegroupSeedBlocks(
                unit_tensor_index,unit_tensor_sectors,unit_tensor_labels,
                lgi_grouped_seed_labels_temp,unit_tensor_spncci_matrices,
                lsu3shell_index_lookup_table,restrict_seeds
              );

            // Matrix copy into final container
            #pragma omp critical(matrix_copy)
              unit_tensor_spncci_matrices_array[unit_tensor_index]=unit_tensor_spncci_matrices;

          }

       #pragma omp critical(regroup)
          { 
            std::cout<<"thread "<<omp_get_thread_num()<<std::endl;
            std::cout<<"size "<<lgi_grouped_seed_labels_temp.size()<<std::endl;
            for(auto it=lgi_grouped_seed_labels_temp.begin(); it!=lgi_grouped_seed_labels_temp.end(); ++it)
              {
                auto& seed_labels=it->second;
                std::cout<<"here"<<std::endl;
                auto& accumulated_seed_labels=lgi_grouped_seed_labels[it->first];
                std::cout<<"there"<<std::endl;
                accumulated_seed_labels.insert(accumulated_seed_labels.end(),seed_labels.begin(), seed_labels.end());
                std::cout<<"thar"<<std::endl;
              }
          }
      }
  }


  void ComputeUnitTensorSeedBlocks(
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const std::string& relative_unit_tensor_filename_template,
      const u3shell::SpaceU3SPN& lsu3shell_space_bra, 
      const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table_bra,
      const u3shell::SpaceU3SPN& lsu3shell_space_ket, 
      const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table_ket,
      const MultiplicityTaggedLGIVector& lgi_vector,
      const basis::OperatorBlocks<double>& lgi_expansions,
      lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels,
      std::vector<basis::OperatorBlocks<double>>& unit_tensor_spncci_matrices_array,
      std::vector<int>& lsu3shell_index_lookup_table,
      bool restrict_seeds
    )
  {
    // Compute unit tensor seed blocks from lsu3shell unit tensor rmes and lgi expansion 
    // and write blocks to file in binary format
    //
    // Container for seed labels and indices for lookup when writing to file
    // lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels

    unit_tensor_spncci_matrices_array.resize(lgi_unit_tensor_labels.size());

    #pragma omp parallel
      {
        lgi::LGIGroupedSeedLabels lgi_grouped_seed_labels_temp;
       
        // For each unit tensor 
        #pragma omp for schedule(dynamic) nowait
        for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
          {
            // Get labels of corresponding unit tensor 
            basis::OperatorBlocks<double> unit_tensor_spncci_matrices;
            const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];

            // Get file name containing lsu3shell rmes of unit tensor 
            std::string filename = fmt::format(relative_unit_tensor_filename_template,unit_tensor_index);

            // generate unit tensor sectors in lsu3shell basis
            const bool spin_scalar = false; // All spin values allowed
            u3shell::SectorsU3SPN unit_tensor_sectors(lsu3shell_space_bra,lsu3shell_space_ket,unit_tensor_labels,spin_scalar);

            // std::cout<<"Read in lsu3shell rmes of unit tensor"<<std::endl;
            basis::OperatorBlocks<double> unit_tensor_lsu3shell_matrices;
            // bool include_sp3r_generators=false;
            double scale_factor=1.0;
            lsu3shell::ReadLSU3ShellRMEs(
                filename,
                lsu3shell_basis_table_bra,lsu3shell_space_bra,
                lsu3shell_basis_table_ket,lsu3shell_space_ket,
                unit_tensor_labels,unit_tensor_sectors,
                unit_tensor_lsu3shell_matrices,scale_factor
              );

            //Create lookup table to be able to identify which lgi the lsu3shell subspace correspond to 
            std::unordered_map<u3shell::U3SPN,int,boost::hash<u3shell::U3SPN>> lgi_lookup_table;
            int i=0;
            for(auto lgi_tagged : lgi_vector)
              {
                lgi_lookup_table[lgi_tagged.irrep]=i;
                ++i;
              }

            lgi::TransformOperatorToSpBasis(
              lsu3shell_space_bra,lsu3shell_space_ket,
                unit_tensor_sectors,lgi_expansions,
                unit_tensor_lsu3shell_matrices,
                lgi_lookup_table,
                unit_tensor_spncci_matrices
                );

////// TODO: finish updating to take different bra and ket subspaces...

            // std::cout<<"Re-organize unit tensor seed bocks "<<std::endl;
            // by lgi then by tensor
            lgi::RegroupSeedBlocks(
                unit_tensor_index,unit_tensor_sectors,unit_tensor_labels,
                lgi_grouped_seed_labels_temp,unit_tensor_spncci_matrices,
                lsu3shell_index_lookup_table,restrict_seeds
              );

            // Matrix copy into final container
            #pragma omp critical(matrix_copy)
              unit_tensor_spncci_matrices_array[unit_tensor_index]=unit_tensor_spncci_matrices;

          }

       #pragma omp critical(regroup)
          { 
            std::cout<<"thread "<<omp_get_thread_num()<<std::endl;
            for(auto it=lgi_grouped_seed_labels_temp.begin(); it!=lgi_grouped_seed_labels_temp.end(); ++it)
              {
                auto& seed_labels=it->second;
                auto& accumulated_seed_labels=lgi_grouped_seed_labels[it->first];
                accumulated_seed_labels.insert(accumulated_seed_labels.end(),seed_labels.begin(), seed_labels.end());
              }
          }
      }
  }



  void WriteHeader(std::ostream& out_file)
  {
    // write binary file header
    //
    // file format code
    mcutils::WriteBinary<int>(out_file,lgi::binary_format_code);
    // floating point precision
    mcutils::WriteBinary<int>(out_file,lgi::binary_float_precision);
  }

  void WriteUnitTensorLabels(
    const std::pair<u3shell::RelativeUnitTensorLabelsU3ST,int>& unit_tensor_labels_tagged,
    std::ostream& out_file
    )
  {
    // Extract unit tensor labels 
    u3::SU3 x0; 
    HalfInt S0,T0,Sp,Tp,S,T;
    int etap,eta;
    std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=unit_tensor_labels_tagged.first.FlatKey();
    int rho0=unit_tensor_labels_tagged.second;


    out_file<<x0.lambda()<<"  "<<x0.mu()<<"  "<<TwiceValue(S0)<<"  "<<TwiceValue(T0)<<"  "
      <<etap<<"  "<<TwiceValue(Sp)<<"  "<<TwiceValue(Tp)<<"  "
      <<eta<<"  "<<TwiceValue(S)<<"  "<<TwiceValue(T)<<"  "<<rho0
      <<std::endl;
      
  }

  void WriteMatrix(const basis::OperatorBlock<double>& block, std::ostream& out_file)
  // Writes rho0, num_rows, num_cols, rmes
  // rmes are by row then by column
  {
    int num_rows=block.rows();
    int num_cols=block.cols();

    assert(num_rows==static_cast<lgi::RMEIndexType>(num_rows));
    mcutils::WriteBinary<lgi::RMEIndexType>(out_file,num_rows);

    assert(num_cols==static_cast<lgi::RMEIndexType>(num_cols));
    mcutils::WriteBinary<lgi::RMEIndexType>(out_file,num_cols);

    
    for(int j=0; j<num_cols; ++j)
      for(int i=0; i<num_rows; ++i)
        {
          auto rme=block(i,j);

          if (lgi::binary_float_precision==4)
            mcutils::WriteBinary<float>(out_file,rme);
          else if (lgi::binary_float_precision==8)
            mcutils::WriteBinary<double>(out_file,rme);

        }
  }


  void WriteSeedsToFile(
      const lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels,
      const std::vector<basis::OperatorBlocks<double>>& unit_tensor_spncci_matrices_array
    )
  {
    std::vector<std::pair<int,int>> irrep_pairs; 
    for(auto it=lgi_grouped_seed_labels.begin(); it!=lgi_grouped_seed_labels.end(); ++it)
      irrep_pairs.push_back(it->first);

    #pragma omp parallel
    {
    #pragma omp for schedule(dynamic)
    for(int i=0; i<irrep_pairs.size(); ++i)
      {
        int lgi_bra_index, lgi_ket_index;
        std::tie(lgi_bra_index,lgi_ket_index)=irrep_pairs[i];

        // temporary container of unit tensor labels and rho0 for writing to separate file
        std::vector<std::pair<u3shell::RelativeUnitTensorLabelsU3ST,int>> unit_tensor_sector_labels;

        // output filename
        std::string seed_filename=fmt::format("seeds/seeds_{:06d}_{:06d}.rmes",lgi_bra_index,lgi_ket_index);
                
        // open output file
        // std::cout << "opening " << seed_filename << std::endl;

        // output in binary mode 
        std::ios_base::openmode mode_argument = std::ios_base::out;
        mode_argument |= std::ios_base::binary;
        std::ofstream seed_file;
        seed_file.open(seed_filename,mode_argument);

        if (!seed_file)
         {
            std::cerr << "Could not open file '" << seed_filename << "'!" << std::endl;
            assert(seed_file);
         }

        // Indicate file format code and binary precision
        WriteHeader(seed_file);

        // for each unit tensor, write non-zero sectors to file
        auto& seed_labels_list=lgi_grouped_seed_labels.at(irrep_pairs[i]);
        for(const SeedLabels& seed_labels : seed_labels_list)
          {
            u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels;
            int rho0, index1,index2; 
            std::tie(unit_tensor_labels,rho0,index1,index2)=seed_labels;

            // add unit tensor to list 
            unit_tensor_sector_labels.emplace_back(unit_tensor_labels,rho0);
                
            // get block
            const basis::OperatorBlock<double>& block=unit_tensor_spncci_matrices_array[index1][index2];

            // write block
            // #pragma omp critical (write_seeds)
              lgi::WriteMatrix(block, seed_file);                
            
          }

        seed_file.close();

        // Write unit tensor labels to file 
        std::string operator_filename=fmt::format("seeds/operators_{:06d}_{:06d}.dat",lgi_bra_index,lgi_ket_index);
                
        // open output file
        // std::cout << "opening " << operator_filename << std::endl;

        // output in binary mode 
        std::ios_base::openmode mode_argument_op = std::ios_base::out;
        std::ofstream operator_file;
        operator_file.open(operator_filename,mode_argument_op);

        if (!operator_file)
          {
            std::cerr << "Could not open file '" << operator_filename << "'!" << std::endl;
            assert(operator_file);
          }

        for(auto labels : unit_tensor_sector_labels)
            lgi::WriteUnitTensorLabels(labels,operator_file);

        operator_file.close();
      }//end for loop
    }//end parallel region
  }


  void WriteSeedsToFile(
    const basis::OperatorBlocks<double>& unit_tensor_seed_blocks,
    int lgi_bra_index, int lgi_ket_index
  )
  {
    // output filename
    std::string seed_filename=fmt::format("seeds/seeds_{:06d}_{:06d}.rmes",lgi_bra_index,lgi_ket_index);
    std::cout<<"writing to "<<seed_filename<<std::endl;
    // output in binary mode 
    std::ios_base::openmode mode_argument = std::ios_base::out;
    mode_argument |= std::ios_base::binary;
    std::ofstream seed_file;
    seed_file.open(seed_filename,mode_argument);

    if (!seed_file)
     {
        std::cerr << "Could not open file '" << seed_filename << "'!" << std::endl;
        return;
     }

    // Indicate file format code and binary precision
    lgi::WriteHeader(seed_file);

    // for each unit tensor, write non-zero sectors to file
    for(const auto& block : unit_tensor_seed_blocks)
        lgi::WriteMatrix(block, seed_file);                

    seed_file.close();
  }




  bool ReadUnitTensorLabels(
      std::string& filename,
      std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
      std::vector<int>& rho0_values
    )
    {
      // open output file
      bool file_found=false;
      // std::cout << "opening " << filename << std::endl;      
      std::ifstream in_stream(filename,std::ios_base::in);
      if(not bool(in_stream))
        {
          // std::cout<<filename+" not found."<<std::endl;
          return file_found;
        }
      // mcutils::StreamCheck(bool(in_stream),filename,"Failure opening "+filename);
      else
        file_found=true;
      // scan input file
      std::string line;
      int line_count=0;
      while ( std::getline(in_stream, line) )
        {
          // count line
          ++line_count;

          // set up for parsing
          std::istringstream line_stream(line);

          // parse line
          int lambda0,mu0,twice_S0,twice_T0,etap,twice_Sp,twice_Tp,
            eta,twice_S,twice_T,rho0;

          line_stream >>lambda0>>mu0>>twice_S0>>twice_T0
              >>etap>>twice_Sp>>twice_Tp
              >>eta>>twice_S>>twice_T>>rho0;

          mcutils::ParsingCheck(line_stream, line_count, line);
          
          // conversions
          HalfInt T0 = HalfInt(twice_T0,2);
          HalfInt S0 = HalfInt(twice_S0,2);
          HalfInt Sp = HalfInt(twice_Sp,2);
          HalfInt Tp = HalfInt(twice_Tp,2);
          HalfInt S = HalfInt(twice_S,2);
          HalfInt T = HalfInt(twice_T,2);
          // std::cout<<fmt::format("{} {} {}", Nsigma, lambda,mu)<<std::endl;
          u3::SU3 x0(lambda0,mu0);

          //Construct unit tensor labels
          u3shell::RelativeStateLabelsU3ST ket(eta,S,T);
          u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp);
          unit_tensor_labels.emplace_back(x0,S0,T0,bra,ket);
          rho0_values.push_back(rho0);
        }
      return file_found;
    }

  void ReadBlock(std::istream& in_stream, Eigen::MatrixXd& block, int float_precision)
    {
      
      //Read in number of rows and columns
      lgi::RMEIndexType rows,cols;
      mcutils::ReadBinary<lgi::RMEIndexType>(in_stream,rows);
      mcutils::ReadBinary<lgi::RMEIndexType>(in_stream,cols);

      // Read in RMEs and case to double matrix 
      //TODO Change to MatrixFloatType for spncci
      if(float_precision==4)
        {
          float buffer[rows*cols];
          in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
          block=Eigen::Map<Eigen::MatrixXf>(buffer,rows,cols).cast<double>();
        }
      else if (float_precision==8)
        {
          double buffer[rows*cols];
          in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
          block=Eigen::Map<Eigen::MatrixXd>(buffer,rows,cols);
        }
    }

  bool ReadBlocks(
      const std::string& filename, 
      int num_blocks,
      basis::OperatorBlocks<double>& blocks
    )
    {
      // std::cout << "opening " << filename << std::endl;      
      bool file_found=false;
      std::ifstream in_stream(filename,std::ios_base::in|std::ios_base::binary);
      if(not bool(in_stream))
        {
          // std::cout<<filename+" not found."<<std::endl;
          return file_found;
        }

      else
        file_found=true;

      // mcutils::StreamCheck(bool(in_stream),filename,"Failure opening "+filename);
      
      // read file header
      int format_code;
      mcutils::ReadBinary<int>(in_stream,format_code);
      assert(format_code==1);
      int float_precision;
      mcutils::ReadBinary<int>(in_stream,float_precision);
      assert((float_precision==4)||(float_precision==8));

      // Read in seed blocks 
      blocks.resize(num_blocks);
      for(int i=0; i<num_blocks; ++i)
        {
          // Read in seeds to matrix in vector
          Eigen::MatrixXd& block=blocks[i];
          ReadBlock(in_stream,block, float_precision);  
        }
      return file_found;
    }
}
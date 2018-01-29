/****************************************************************
  lgi_unit_tensors.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "lgi/lgi_unit_tensors.h"

#include <iostream>
#include <fstream>
#include "mcutils/io.h"
#include "mcutils/eigen.h"
#include "cppformat/format.h"
#include "lgi/lgi_solver.h"
#include "mcutils/parsing.h" //For ReadUnitTensorLabels
#include "mcutils/io.h"//For ReadUnitTensorLabels


namespace lgi
{


  // zero tolerance 
  double zero_tolerance=1e-6;
  
  // output mode
  int binary_format_code = 1;
  int binary_float_precision=8;

  void RegroupSeedBlocks(
      int unit_tensor_index,
      const u3shell::SectorsU3SPN& unit_tensor_sectors,
      const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels,
      lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels,
      basis::MatrixVector& unit_tensor_spncci_matrices,
      bool restrict_seeds
    )
  {
    // Create list of unit tensor sectors by lgi then by unit tensor for writing to file
    for(int sector_index=0; sector_index<unit_tensor_sectors.size(); ++sector_index)
      {
        // Check that the sector has non-zero rmes as defined by the zero_tolerance
        if(mcutils::IsZero(unit_tensor_spncci_matrices[sector_index],lgi::zero_tolerance))
          continue;
            
            // extract U3SPN sector information from unit tensor sectors definded in lsu3shell basis
            const typename u3shell::SectorsU3SPN::SectorType& sector = unit_tensor_sectors.GetSector(sector_index);
            const int bra_subspace_index = sector.bra_subspace_index();
            const int ket_subspace_index = sector.ket_subspace_index();

            // If restrict_seeds=true, keep only if bra>=ket, else keep all seeds
            bool keep=restrict_seeds?(bra_subspace_index>=ket_subspace_index):true;
            if(not keep)
              continue;

            const int rho0=sector.multiplicity_index();

            // Regroup by lgi pair, then by unit tensor
            std::pair<int,int> irrep_family_pair(bra_subspace_index,ket_subspace_index);
            
            // Labels for looking up correct block to write to file
            lgi::SeedLabels seed_labels(unit_tensor_labels,rho0,unit_tensor_index,sector_index);

            // Accumulating list 
            lgi_grouped_seed_labels[irrep_family_pair].push_back(seed_labels);
          
      }// sector index
  }

  void ComputeUnitTensorSeedBlocks(
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const std::string& relative_unit_tensor_filename_template,
      const u3shell::SpaceU3SPN& lsu3shell_space, 
      const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
      const basis::MatrixVector& lgi_expansions,
      lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels,
      std::vector<basis::MatrixVector>& unit_tensor_spncci_matrices_array,
      bool restrict_seeds
    )
  {
    // Compute unit tensor seed blocks from lsu3shell unit tensor rmes and lgi expansion 
    // and write blocks to file in binary format
    //
    // Container for seed labels and indices for lookup when writing to file
    // lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels

    unit_tensor_spncci_matrices_array.resize(lgi_unit_tensor_labels.size());

    // For each unit tensor 
    for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
      {
        // Get labels of corresponding unit tensor 
        basis::MatrixVector& unit_tensor_spncci_matrices=unit_tensor_spncci_matrices_array[unit_tensor_index];
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];

        // Get file name containing lsu3shell rmes of unit tensor 
        std::string filename = fmt::format(relative_unit_tensor_filename_template,unit_tensor_index);

        // generate unit tensor sectors in lsu3shell basis
        const bool spin_scalar = false; // All spin values allowed
        u3shell::SectorsU3SPN unit_tensor_sectors(lsu3shell_space,unit_tensor_labels,spin_scalar);

        // std::cout<<"Read in lsu3shell rmes of unit tensor"<<std::endl;
        basis::MatrixVector unit_tensor_lsu3shell_matrices;
        lsu3shell::ReadLSU3ShellRMEs(
            filename,
            lsu3shell_basis_table,lsu3shell_space,
            unit_tensor_labels,unit_tensor_sectors,unit_tensor_lsu3shell_matrices
          );

        // std::cout<<"Transform seed rmes to SpNCCI basis "<<std::endl;
        lgi::TransformOperatorToSpBasis(
            unit_tensor_sectors,lgi_expansions,
            unit_tensor_lsu3shell_matrices,unit_tensor_spncci_matrices
          );

        // std::cout<<"Re-organize unit tensor seed bocks "<<std::endl;
        // by lgi then by tensor
        lgi::RegroupSeedBlocks(
            unit_tensor_index,unit_tensor_sectors,unit_tensor_labels,
            lgi_grouped_seed_labels,unit_tensor_spncci_matrices,restrict_seeds
          );
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
      

    // if(false)
    //   {
    //     mcutils::WriteBinary<int>(out_file,x0.lambda());
    //     mcutils::WriteBinary<int>(out_file,x0.mu());
    //     mcutils::WriteBinary<int>(out_file,TwiceValue(S0));
    //     mcutils::WriteBinary<int>(out_file,TwiceValue(T0));
    //     mcutils::WriteBinary<int>(out_file,etap);
    //     mcutils::WriteBinary<int>(out_file,TwiceValue(Sp));
    //     mcutils::WriteBinary<int>(out_file,TwiceValue(Tp));
    //     mcutils::WriteBinary<int>(out_file,eta);
    //     mcutils::WriteBinary<int>(out_file,TwiceValue(S));
    //     mcutils::WriteBinary<int>(out_file,TwiceValue(T));
    //     mcutils::WriteBinary<int>(out_file,num_sectors);
    //   }
  }

  void WriteMatrix(const basis::OperatorBlock<double>& block, std::ostream& out_file)
  // Writes rho0, num_rows, num_cols, rmes
  // rmes are by row then by column
  {
    int num_rows=block.rows();
    int num_cols=block.cols();

    assert(num_rows==static_cast<RMEIndexType>(num_rows));
    mcutils::WriteBinary<RMEIndexType>(out_file,num_rows);

    assert(num_cols==static_cast<RMEIndexType>(num_cols));
    mcutils::WriteBinary<RMEIndexType>(out_file,num_cols);

    
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
      const std::vector<basis::MatrixVector>& unit_tensor_spncci_matrices_array
    )
  {

    for(auto it=lgi_grouped_seed_labels.begin(); it!=lgi_grouped_seed_labels.end(); ++it)
      {
        int lgi_bra_index, lgi_ket_index;
        std::tie(lgi_bra_index,lgi_ket_index)=it->first;

        // temporary container of unit tensor labels and rho0 for writing to separate file
        std::vector<std::pair<u3shell::RelativeUnitTensorLabelsU3ST,int>> unit_tensor_sector_labels;

        // output filename
        std::string seed_filename=fmt::format("seeds/seeds_{:06d}_{:06d}.rmes",lgi_bra_index,lgi_ket_index);
                
        // open output file
        std::cout << "opening " << seed_filename << std::endl;
        ;

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
        WriteHeader(seed_file);

        // for each unit tensor, write non-zero sectors to file
        for(const SeedLabels& seed_labels : it->second)
          {
            u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels;
            int rho0, index1,index2; 
            std::tie(unit_tensor_labels,rho0,index1,index2)=seed_labels;

            // add unit tensor to list 
            unit_tensor_sector_labels.emplace_back(unit_tensor_labels,rho0);
                
            // get block
            const basis::OperatorBlock<double>& block=unit_tensor_spncci_matrices_array[index1][index2];

            // write block
            lgi::WriteMatrix(block, seed_file);                
            
          }

        seed_file.close();

        // Write unit tensor labels to file 
        std::string operator_filename=fmt::format("seeds/operators_{:06d}_{:06d}.dat",lgi_bra_index,lgi_ket_index);
                
        // open output file
        std::cout << "opening " << operator_filename << std::endl;

        // output in binary mode 
        std::ios_base::openmode mode_argument_op = std::ios_base::out;
        std::ofstream operator_file;
        operator_file.open(operator_filename,mode_argument_op);

        if (!operator_file)
          {
            std::cerr << "Could not open file '" << operator_filename << "'!" << std::endl;
            return;
          }

        for(auto labels : unit_tensor_sector_labels)
            lgi::WriteUnitTensorLabels(labels,operator_file);

        operator_file.close();
      }
  }

  bool ReadUnitTensorLabels(
      std::string& filename,
      std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
      std::vector<int>& rho0_values
    )
    {
      // open output file
      bool file_found=false;
      std::cout << "opening " << filename << std::endl;      
      std::ifstream in_stream(filename,std::ios_base::in);
      if(not bool(in_stream))
        {
          std::cout<<filename+" not found."<<std::endl;
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
      mcutils::ReadBinary<RMEIndexType>(in_stream,rows);
      mcutils::ReadBinary<RMEIndexType>(in_stream,cols);

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
      basis::MatrixVector& blocks
    )
    {
      std::cout << "opening " << filename << std::endl;      
      bool file_found=false;
      std::ifstream in_stream(filename,std::ios_base::in|std::ios_base::binary);
      if(not bool(in_stream))
        {
          std::cout<<filename+" not found."<<std::endl;
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
/****************************************************************
  lgi_unit_tensors.h

  Code for computing unit tensor rmes between lgi families from
  unit tensor rmes computed using the LSU3Shell SU3RME code and
  writes unit tensor blocks to file

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  12/29/17 (aem): Created. Based on GetUnitTensorSeedBlocks from
    spncci/computation_control
****************************************************************/
#ifndef LGI_LGI_UNIT_TENSORS_H_
#define LGI_LGI_UNIT_TENSORS_H_

#include "lgi/lgi.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "lsu3shell/lsu3shell_rme.h"


namespace lgi
{
    // zero tolerance 
  extern double zero_tolerance;
  
  // output mode
  extern int binary_format_code;
  extern int binary_float_precision;
  typedef short unsigned int RMEIndexType;

  typedef std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int,int> SeedLabels;
  typedef std::unordered_map<std::pair<int,int>,std::vector<SeedLabels>,boost::hash<std::pair<int,int>>> 
    LGIGroupedSeedLabels;

  void ComputeUnitTensorSeedBlocks(
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const std::string& relative_unit_tensor_filename_template,
      const u3shell::SpaceU3SPN& lsu3shell_space, 
      const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
      const basis::MatrixVector& lgi_expansions,
      lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels,
      std::vector<basis::MatrixVector>& unit_tensor_spncci_matrices_array,
      bool restrict_seeds=false
    );

    // Compute seed blocks for SpNCCI recurrence from LSU3Shell rmes.  The seeds are stored in a matrix array
    // by unit tensor then by unit tensor lgi sector index.
    //
    // Corresponding to the seed blocks is a nested unordered map of 
    // {lgi_pair: {unit_tensor_labels: (rho0,unit_tensor_index,sector_index)}
    // Only those sectors for lgi_bra>=lgi_ket and which have non-zero rmes are included.
    //
    // The tuple of labels (rho0,unit_tensor_index,sector_index) allow for a simple lookup of corresponding blolck
    // in array.
    //
    // If restrict_seeds=true, then only keep seeds for which bra>=ket


  void WriteSeedsToFile(
      const lgi::LGIGroupedSeedLabels& lgi_grouped_seed_labels,
      const std::vector<basis::MatrixVector>& unit_tensor_spncci_matrices_array
    );
    // For each LGI pair (bra>=ket), write seed rmes to file "seeds_{bra_lgi_index}_ket_lgi_index}.rmes" 
    // and creates corresponding list of unit tensor labels in "operators_{bra_lgi_index}_{ket_lgi_index}".
    // In the seed file, the rmes are 
    //    by unit tensor (order given by unit_tensor_file)
    //      by rho0
    //
    // Unit tensor file:
    //    lambda0 mu0 S0 T0 etap Sp Tp eta S T num_sectors(number of non-zero sectors ~ rho0max)
    //
    // seed file:
    //    rho0, rows, columns, rmes.... 

  void ReadUnitTensorLabels(
      std::string& filename,
      std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
      std::vector<int>& rho0_values
    );
    // Reads in unit tensor labels from file and stores them in a vector of pairs of unit tensor 
    // labels and the rho0 of the corresponding unit tensor block.  
    //  Labels in file are:
    //    lambda0, mu0, 2S0, 2T0, etap, 2Sp, 2Tp, eta, 2S, 2T, rho0

  void ReadSeedBlocks(
      const std::string& filename, 
      int num_blocks,
      basis::MatrixVector& unit_tensor_seed_blocks
    );
    // Reads in unit tensor seed blocks from file and stores them in a vector of blocks
    // order of blocks corresponds to order of (unit_tensor,rho0) pairs in corresponding
    // operator file. 

}


#endif
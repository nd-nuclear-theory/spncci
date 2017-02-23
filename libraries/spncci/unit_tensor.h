/****************************************************************
  unit_tensor.h

  Unit tensor recursive evaluation
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/15/16 (aem,mac): Created.
  4/14/16 (aem): Added Np,N sector look ups and changed iteration order
  4/20/16 (aem): Defined GenerateUnitTensorU3SectorLabels function and 
                 changed interation in UnitTensorMatrixGenerator
  12/7/16 (aem): Overhall of unit tensor rme calculation
  12/21/16 (aem): Factored out conjugate sector calculation
  2/2/17 (mac): Add typedef UnitTensorMatricesByIrrepFamily.
  2/20/17 (mac): Update argument types on GenerateUnitTensorMatrix.
  2/23/17 (aem): Add conjugate sector flag to GenerateUnitTensorU3SectorLabels.
****************************************************************/

#ifndef SPNCCI_SPNCCI_UNIT_TENSOR_H_
#define SPNCCI_SPNCCI_UNIT_TENSOR_H_

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash_fwd.hpp>
#include <eigen3/Eigen/Eigen>

#include "spncci/spncci_basis.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/vcs.h"
#include "u3shell/tensor_labels.h"


namespace spncci
{

  typedef   std::map<std::pair<u3::U3,HalfInt>,int> U3SCount;

  ////////////////////////////////////////////////////////////////
  // unit tensor matrix U(3) sector labels
  ////////////////////////////////////////////////////////////////

  class UnitTensorU3Sector
  {
    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

    public:
    typedef std::tuple< u3::U3, u3::U3, u3shell::RelativeUnitTensorLabelsU3ST, int> KeyType;
    // w', w, Unit_tensor_labels,r0

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
      
    inline UnitTensorU3Sector() 
      : rho0_(0) {}
      
    // construction from labels
    inline UnitTensorU3Sector(u3::U3 omegap, u3::U3 omega, u3shell::RelativeUnitTensorLabelsU3ST tensor, int rho0)
      : omegap_(omegap), omega_(omega), tensor_(tensor), rho0_(rho0) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      return KeyType(omegap_,omega_,tensor_,rho0_);
    }

    ////////////////////////////////////////////////////////////////
    // relational operators
    ////////////////////////////////////////////////////////////////

    inline friend bool operator == (const UnitTensorU3Sector& unit1rme, const UnitTensorU3Sector& unit2rme)
    {
      return unit1rme.Key() == unit2rme.Key();
    }

    inline friend bool operator < (const UnitTensorU3Sector& unit1rme, const UnitTensorU3Sector& unit2rme)
    {
      return unit1rme.Key() < unit2rme.Key();
    }
    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    inline friend std::size_t hash_value(const UnitTensorU3Sector& tensor)
    {
      boost::hash<UnitTensorU3Sector::KeyType> hasher;
      return hasher(tensor.Key());
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
      
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    private:
    // Operator labels
    u3::U3 omegap_,omega_;
    u3shell::RelativeUnitTensorLabelsU3ST tensor_;
    int rho0_;

  };

  ////////////////////////////////////////////////////////////////
  // typedefs
  ////////////////////////////////////////////////////////////////
#ifdef HASH_UNIT_TENSOR
  typedef std::unordered_map< spncci::UnitTensorU3Sector, Eigen::MatrixXd, boost::hash<spncci::UnitTensorU3Sector> > UnitTensorSectorsCache;
#else
  typedef std::map< spncci::UnitTensorU3Sector,Eigen::MatrixXd > UnitTensorSectorsCache;
#endif

  void 
    GenerateUnitTensorU3SectorLabels(
        int N1b,
        int Nmax,
        std::pair<int,int>  irrep_family_indices,
        const spncci::SpNCCISpace& sp_irrep_vector,
        const std::map< int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>>& unit_tensor_labels_map,
        std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>>& unit_tensor_NpN_sector_map
      );
  // Generates labels for (omega' omega unit_tensor) sectors for N0>=0 unit tensors.  The N0<0
  // unit tensor can be obtained by taking the adjoint of the sector. Labels are sorted by
  // Nnp and Nn of bra and ket respectively, where Nn and Nnp are the excitation quanta within the 
  // irrep. 
  //
  // Note that N0 must equal the difference between the Nn of the states as well as the difference
  // between Nsigma, i.e., N0=Nsigma_ex_p+Nnp-Nsigma_ex-Nn, 
  //
  // Arguments: 
  //
  // N1b (input) : single particle cutoff, relative particle <= 2*N1b+Nn
  // Nmax (input) : Nmax of truncation
  // irrep_family_indices (input) : sector pair given as index pair from global list sp_irrep_vector 
  // sp_irrep_vector (input) : vector of multiplicty tagged SpIrreps
  // unit_tensor_labels_map (input) : map with list of unit tensor labels with key N0 
  // unit_tensor_NpN_sector_map (output) :  For each NpN pair key in map the corresponding value 


  typedef std::map<std::pair<int,int>,std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache>>
    UnitTensorMatricesByIrrepFamily;
  // Mapping to store matrices for unit tensors as caculated in recurrence:
  //
  // (bra_irrep_family_index,ket_irrep_family_index) (<int,int>)
  //   -> (bra_Nn,ket_Nn) (<int,int>)
  //   -> unit_tensor_sector_labels (spncci::UnitTensorU3Sector)
  //   -> sector matrix (Eigen::MatrixXd)

  void 
    GenerateUnitTensorMatrix(
        int N1b,
        int Nmax, 
        const std::pair<int,int> irrep_family_indices,
        const spncci::SpNCCISpace& sp_irrep_vector,
        u3::UCoefCache& u_coef_cache,
        u3::PhiCoefCache& phi_coef_cache,
        const spncci::KMatrixCache& k_matrix_map,
        const std::map< int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>>& unit_tensor_labels_map,
        spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices
      );
  // Uses recurrence method to compute rme's of unit tensors in unit_tensor_labels_map
  // and stores rme's in unit_tensor_rme_map. 
  //
  // Arguments:
  //   N1b (input) : [=N1v] valence shell, used to determine cutoff Nrel<=2N1b+Nex. 
  //   Nmax (input) : boson number truncation
  //      [(mac): this is actually in NCCI sense, i.e., excitation relative to Nsigma0; TODO: no need to hard
  //      code this across all irreps; rather, just use the truncation for each irrep that was built into
  //      the SpNCCI space; truncators may be more general than Nmax based]
  //   irrep_family_indices (input) : integers (irrep_family_index_bra,irrep_family_index_ket) indicate which irrep family sector to recurse.
  //   sp_irrep_vector (input) : vector of mutliplicity tagged sp irreps (TO UPDATE)
  //   u_coef_cach (input,output) : cache containing already computed U coefs.  As coefs
  //                                computed, cache is updated.
  //   phi_coef_cach (input,output) : cache containing already computed phi coefs.  As coefs
  //                                computed, cache is updated.
  //   k_matrix_map (input) : cache containing precomputed Kmatrices
  //   unit_tensor_labels_map (input) : map with list of unit tensor labels with key N0 
  //   unit_tensor_matrices (input, output) : seed rmes (input) and computed rmes (output)

  // void RegroupUnitTensorU3SSectors(
  //       bool is_new_subsector,
  //       const HalfInt& Sp, const HalfInt& S,
  //       std::pair<int,int> irrep_family_indices,
  //       const std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache>& unit_teXFnsor_rme_map,
  //       const std::pair<int,int>& sp_irrep_symmetry_sum,
  //       UnitTensorU3SSectorsCache& unit_tensor_u3S_cache
  //     );
} //namespace 

#endif

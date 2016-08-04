/****************************************************************
  unit_tensor.h

  Unit tensor recursive evaluation
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/15/16 (aem,mac): Created.
  4/14/16 (aem): Added Np,N sector look ups and changed iteration order
  4/20/16 (aem): Defined GenerateUnitTensorU3SectorLabels function and 
               changed interation in UnitTensorMatrixGenerator
****************************************************************/

#ifndef UNIT_TENSOR_H_
#define UNIT_TENSOR_H_

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash_fwd.hpp>
#include <eigen3/Eigen/Eigen>

#include "spncci/sp_basis.h"
#include "sp3rlib/u3.h"
#include "u3shell/tensor_labels.h"


namespace spncci
{

  typedef   std::map<std::pair<u3::U3,HalfInt>,int> U3SCount;

  ////////////////////////////////////////////////////////////////
  // unit tensor operator labels
  ////////////////////////////////////////////////////////////////

  class UnitTensor
  {
    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

  public:
    typedef std::tuple<u3::SU3, HalfInt, HalfInt, int, HalfInt, HalfInt, 
                        int, HalfInt, HalfInt> KeyType;
    // w0, S0, T0, rbp, Sb, Tb, r, S, T

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // default constructor
    inline UnitTensor() 
      : S0_(0), T0_(0), Sp_(0), Tp_(0), S_(0), T_(0) {}
    
    // construction from labels
    inline UnitTensor(u3::SU3 x0, HalfInt S0, HalfInt T0, 
                      int rp, HalfInt Sp, HalfInt Tp, 
                      int r, HalfInt S, HalfInt T)
      : x0_(x0), S0_(S0), T0_(T0), 
        rp_(rp), Sp_(Sp), Tp_(Tp), 
        r_(r), S_(S), T_(T) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
    inline u3::SU3 x0() const
    {
      return x0_;
    }


    inline KeyType Key() const
    {
      return KeyType(x0_,S0_,T0_,rp_,Sp_,Tp_,r_,S_,T_);
    }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    static const int quanta_label_width = 6;
    static const int spin_label_width = 2;
    inline friend std::size_t hash_value(const UnitTensor& tensor)
    {
      // TODO (Andika)
      // pack ints and Halfnts together first
      // x0, S0, T0, rp, Sp, Tp, r, S, T
      // types: U3 (3 HalfInt), HalfInt, HalfInt, int, HalfInt, HalfInt, int, HalfInt, HalfInt
      // then hash the packed values together and combine
      // U3 hashing defined in sp3rlib/u3.h
      // 
      // r and rp need 6 bits
      // TwiceValue(S0), TwiceValue(T0), etc., need only 2 bits
      // See constants defined just above...

      // TODO post-Andika: make this a lot more transparent as a sequence of "shift then or" operations
      // TODO post-Andika: fix up pseudo-Java-esque variable name conventions
      // TODO post-Andika: uniformly use hash_value instead of declaring boost hasher objects
      // TODO post-Andika: remove naive hash, or make it and #els alternative to BOOST_HASH


      int packed_Ints = (TwiceValue(tensor.S0_) << 5*spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.T0_) << 4*spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.Sp_) << 3*spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.Tp_) << 2*spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.S_) << spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.T_) << 2*quanta_label_width)
        | (tensor.rp_ << quanta_label_width) | (tensor.r_ << 0);
      
      boost::hash<int> intHasher;
      std::size_t intHash = intHasher(packed_Ints);
      std::size_t u3Hash = hash_value(tensor.x0_);
      
      #ifdef NAIVEHASH 
      // naive implementation: add hashes together
      return intHash+u3Hash;
      #endif

      #ifdef BOOSTHASH      
      // smart implementation: boost::hash_combine
      std::size_t seed = 0;
      boost::hash_combine(seed,intHash);
      boost::hash_combine(seed,u3Hash);
      return seed;
      #endif
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    // Operator labels
    u3::SU3 x0_;
    HalfInt S0_, T0_, Sp_, Tp_, S_, T_;
    int rp_, r_;

  };

  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const UnitTensor& unit1, const UnitTensor& unit2)
  {
    return unit1.Key() == unit2.Key();
  }

  inline bool operator < (const UnitTensor& unit1, const UnitTensor& unit2)
  {
    return unit1.Key() < unit2.Key();
  }

  //////////////////////////////////////////////////////////////////
  //  Additional constructors 
  //////////////////////////////////////////////////////////////////
  // Constructing UnitTensor from RelativeUnitTensorLabelsU3ST
  inline spncci::UnitTensor RelativeUnitTensor(const u3shell::RelativeUnitTensorLabelsU3ST& tensor)
  {
    u3::SU3 x0;
    HalfInt S0, T0, Sp, Tp, S, T;
    int rp,  r;
    std::tie(x0,S0,T0,rp,Sp,Tp,r,S,T)=tensor.FlatKey();
    return UnitTensor(x0,S0,T0,rp,Sp,Tp,r,S,T);
  } 


  ////////////////////////////////////////////////////////////////
  // unit tensor matrix U(3) sector labels
  ////////////////////////////////////////////////////////////////

  class UnitTensorU3Sector
  {
    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

  public:
    typedef std::tuple< u3::U3, u3::U3, spncci::UnitTensor, int> KeyType;
    // w', w, Unit tensor labels,r0

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
      
    inline UnitTensorU3Sector() 
      : rho0_(0) {}
      
    // construction from labels
    inline UnitTensorU3Sector(u3::U3 omegap, u3::U3 omega, spncci::UnitTensor tensor, int rho0)
      : omegap_(omegap), omega_(omega), tensor_(tensor), rho0_(rho0) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      return KeyType(omegap_,omega_,tensor_,rho0_);
    }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    inline friend std::size_t hash_value(const UnitTensorU3Sector& tensor)
    {
      // TODO (Andika)
      // labels omegap_, omega_, tensor_, rho0_
      // type: U3, U3, tensor, int
      // hash them all
      std::size_t omegapHash = hash_value(tensor.omegap_);
      std::size_t omegaHash = hash_value(tensor.omega_);
      std::size_t tensorHash = hash_value(tensor.tensor_);
      boost::hash<int> rhoHasher;
      std::size_t rhoHash = rhoHasher(tensor.rho0_);
     
      #ifdef NAIVEHASH 
      // naive implementation: add hashes together
      return omegapHash + omegaHash + tensorHash + rhoHash;
      #endif
      
      #ifdef BOOSTHASH
      // smart implementation: boost::hash_combine
      std::size_t seed = 0;
      boost::hash_combine(seed,omegapHash);
      boost::hash_combine(seed,omegaHash);
      boost::hash_combine(seed,tensorHash);
      boost::hash_combine(seed,rhoHash);
      return seed;
      #endif
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
    spncci::UnitTensor tensor_;
    int rho0_;

  };

  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const UnitTensorU3Sector& unit1rme, const UnitTensorU3Sector& unit2rme)
  {
    return unit1rme.Key() == unit2rme.Key();
  }

  inline bool operator < (const UnitTensorU3Sector& unit1rme, const UnitTensorU3Sector& unit2rme)
  {
    return unit1rme.Key() < unit2rme.Key();
  }

  ////////////////////////////////////////////////////////////////
  // typedefs
  ////////////////////////////////////////////////////////////////
  #ifdef HASH_UNIT_TENSOR
    typedef std::unordered_map< spncci::UnitTensorU3Sector, Eigen::MatrixXd, boost::hash<spncci::UnitTensorU3Sector> > UnitTensorSectorsCache;
  #else
    typedef std::map< spncci::UnitTensorU3Sector,Eigen::MatrixXd > UnitTensorSectorsCache;
  #endif

  ////////////////////////////////////////////////////////////////
  // unit tensor matrix U(3) sector labels
  ////////////////////////////////////////////////////////////////

  class UnitTensorU3SSector
  {
    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

  public:
    typedef std::tuple<	u3::U3, HalfInt, u3::U3, HalfInt, spncci::UnitTensor, int> KeyType;
    // w'S', wS, Unit tensor labels,r0

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
      
    inline UnitTensorU3SSector() 
      : rho0_(0) {}
      
    // construction from labels
    inline UnitTensorU3SSector(u3::U3 omegap, HalfInt Sp, u3::U3 omega, HalfInt S, spncci::UnitTensor tensor, int rho0)
      : omegap_(omegap), Sp_(Sp), omega_(omega), S_(S), tensor_(tensor), rho0_(rho0) {}

    // construction from labels
    inline UnitTensorU3SSector(HalfInt Sp, HalfInt S, UnitTensorU3Sector sector)
      : Sp_(Sp), S_(S) 
      {
        u3::U3 omegap,omega;
        spncci::UnitTensor tensor;
        int rho0;
        std::tie(omegap,omega,tensor,rho0)=sector.Key();
        omegap_=omegap;
        omega_=omega;
        tensor_=tensor;
        rho0_=rho0;
      }


    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      return KeyType(omegap_,Sp_,omega_,S_,tensor_,rho0_);
    }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    inline friend std::size_t hash_value(const UnitTensorU3SSector& sector)
    {

      boost::hash<UnitTensorU3SSector::KeyType> hasher;
        return hasher(sector.Key());
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
    HalfInt Sp_,S_;
    spncci::UnitTensor tensor_;
    int rho0_;

  };

  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const UnitTensorU3SSector& unit1rme, const UnitTensorU3SSector& unit2rme)
  {
    return unit1rme.Key() == unit2rme.Key();
  }

  inline bool operator < (const UnitTensorU3SSector& unit1rme, const UnitTensorU3SSector& unit2rme)
  {
    return unit1rme.Key() < unit2rme.Key();
  }

  typedef std::unordered_map< spncci::UnitTensorU3SSector, std::vector<Eigen::MatrixXd>, boost::hash<spncci::UnitTensorU3SSector> > UnitTensorU3SSectorsCache;






  void GenerateUnitTensors(int Nmax, std::map< int,std::vector<spncci::UnitTensor> >& unit_sym_map);
  // Generates a map containing (key, value) pair (N0, operator_labels) of the unit tensors 


  void GenerateUnitTensorMatrix(
         int N1b,
         // boson number cutoff
         int Nmax, 
         // a given spncci sector pair given as index pair  from global list lgi_vector 
         std::pair<int,int> lgi_pair,
         //coefficient cache
         u3::UCoefCache u_coef_cache,
         // Address to map with list of unit tensor labels with key N0 
         std::map< int,std::vector<spncci::UnitTensor>>& unit_sym_map,
         // Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
         std::map<
         std::pair<int,int>,
         UnitTensorSectorsCache
         >& unit_tensor_rme_map
         );

  void RegroupUnitTensorU3SSectors(
        bool is_new_subsector,
        const HalfInt& Sp, const HalfInt& S,
        std::pair<int,int> lgi_pair,
        const std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache>& unit_tensor_rme_map,
        const std::pair<int,int>& lgi_symmetry_sum,
        UnitTensorU3SSectorsCache& unit_tensor_u3S_cache
      );
} //namespace 

#endif

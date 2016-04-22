/****************************************************************
  unit_tensor.h

  Unit tensor algorithms
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/15/16 (aem,mac): Created.
  4/14/16 (aem): Added Np,N sector look ups and changed iteration order
  4/20/16 (aem): Defined GenerateUnitTensorU3SectorLabels function and 
               changed interation in UnitTensorMatrixGenerator
****************************************************************/
#ifndef UNIT_TENSOR_H_
#define UNIT_TENSOR_H_

#include "spncci/sp_basis.h"
#include "sp3rlib/vcs.h"
#include "sp3rlib/u3.h"
#include <unordered_set>
#include <boost/functional/hash_fwd.hpp>

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // unit tensor operator labels
  ////////////////////////////////////////////////////////////////

  class UnitTensor
  {
    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

  public:
    typedef std::tuple<u3::U3, HalfInt, HalfInt, int, HalfInt, HalfInt, 
                        int, HalfInt, HalfInt> KeyType;
    // w0, S0, T0, rbp, Sb, Tb, r, S, T

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // default constructor
    inline UnitTensor() 
      : S0_(0), T0_(0), Sp_(0), Tp_(0), S_(0), T_(0) {}
    
    // construction from labels
    inline UnitTensor(u3::U3 omega0, HalfInt S0, HalfInt T0, 
                      int rp, HalfInt Sp, HalfInt Tp, 
                      int r, HalfInt S, HalfInt T)
      : omega0_(omega0), S0_(S0), T0_(T0), 
        rp_(rp), Sp_(Sp), Tp_(Tp), 
        r_(r), S_(S), T_(T) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
    inline u3::U3 omega0() const
    {
      return omega0_;
    }


    inline KeyType Key() const
    {
      return KeyType(omega0_,S0_,T0_,rp_,Sp_,Tp_,r_,S_,T_);
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
      // omega0, S0, T0, rp, Sp, Tp, r, S, T
      // types: U3 (3 HalfInt), HalfInt, HalfInt, int, HalfInt, HalfInt, int, HalfInt, HalfInt
      // then hash the packed values together and combine
      // U3 hashing defined in sp3rlib/u3.h
      // 
      // r and rp need 6 bits
      // TwiceValue(S0), TwiceValue(T0), etc., need only 2 bits
      // See constants defined just above...
      int packed_Ints = (TwiceValue(tensor.S0_) << 5*spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.T0_) << 4*spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.Sp_) << 3*spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.Tp_) << 2*spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.S_) << spin_label_width+2*quanta_label_width)
        | (TwiceValue(tensor.T_) << 2*quanta_label_width)
        | (rp_ << quanta_label_width) | (r_ << 0);
      
      boost::hash<int> intHasher;
      std::size_t intHash = intHasher(packed_Ints);
      std::size_t u3Hash = u3::hash_value(tensor.omega0_);
      
      // naive implementation: add hashes together
      return intHash+u3Hash;
      
      // smart implementation: boost::hash_combine
      // std::size_t seed = 0;
      // boost::hash_combine(seed,intHash);
      // boost::hash_combine(seed,u3Hash);
      // return seed;
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    // Operator labels
    u3::U3 omega0_;
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

  ////////////////////////////////////////////////////////////////
  // unit tensor matrix U(3) sector labels
  ////////////////////////////////////////////////////////////////

  class UnitTensorU3Sector
  {
    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

  public:
    typedef std::tuple<	u3::U3, u3::U3, spncci::UnitTensor, int> KeyType;
    // w'S', wS, Unit tensor labels,r0

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
      std::size_t omegapHash = u3::hash_value(tensor.omegap_);
      std::size_t omegaHash = u3::hash_value(tensor.omega_);
      std::size_t tensorHash = hash_value(tensor.tensor_);
      boost::hash<int> rhoHasher;
      std::size_t rhoHash = rhoHasher(tensor.rho0_);
      
      // naive implementation: add hashes together
      return omegapHash + omegaHash + tensorHash + rhoHash;
      
      // smart implementation: boost::hash_combine
      // std::type_t seed = 0;
      // boost::hash_combine(seed,omegapHash);
      // boost::hash_combine(seed,omegaHash);
      // boost::hash_combine(seed,tensorHash);
      // boost::hash_combine(seed,rhoHash);
      // return seed;
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



  class UnitU3SectorPair
  {
    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

  public:
    typedef std::pair< UnitTensorU3Sector, Eigen::MatrixXd> KeyType;
    // w'S', wS, Unit tensor labels,r0

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
      
   
    // construction from labels
    inline UnitU3SectorPair(spncci::UnitTensorU3Sector labels, Eigen::MatrixXd sector)
      : labels_(labels), sector_(sector) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      return KeyType(labels_, sector_);
    }

    inline spncci::UnitTensorU3Sector labels() const
    {
      return labels_;
    }

    inline Eigen::MatrixXd sector() const
    {
      return sector_;
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
    spncci::UnitTensorU3Sector labels_;
    Eigen::MatrixXd sector_;
  };


  ////////////////////////////////////////////////////////////////
  // unit tensor operator labels
  ////////////////////////////////////////////////////////////////

  void GenerateUnitTensors(int Nmax, std::map< int,std::vector<spncci::UnitTensor> >& unit_sym_map);
  // Generates a map containing (key, value) pair (N0, operator_labels) of the unit tensors 
  //
  // 

  // inline void TestFunction(  
  //   // boson number cutoff
  //   int Nmax, 
  //   // a given spncci sector pair given as index pair  from global list lgi_vector 
  //   std::pair<int,int> lgi_pair,
  //   // Address to map with list of unit tensor labels with key N0 
  //   std::map< int,std::vector<spncci::UnitTensor>>& unit_sym_map,
  //   // Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
  //   std:: map< std::pair<int,int>, std::map< spncci::UnitTensorU3Sector,Eigen::MatrixXd> >& unit_tensor_rme_map
  //   )
  //   {
  //     U3 sp=lgi_vector[lgi_pair.first].sigma;
  //     U3 s=lgi_vector[lgi_pair.second].sigma;
  //     for (int i=0; i<unit_sym_map.size(); i++)
  //       {
  //         if (u3::OuterMultiplicity(s.SU3(),unit_sym_map[0][i].omega0.SU3(),sp.SU3())>0)
  //           std::cout<< unit_tensor_rme_map[lgi_pair][spncci::UnitTensorU3Sector(sp,s,unit_sym_map[0][i],1)]<<std::endl;
  //       }
  //   }

  void GenerateUnitTensorMatrix(
         int N1b,
         // boson number cutoff
         int Nmax, 
         // a given spncci sector pair given as index pair  from global list lgi_vector 
         std::pair<int,int> lgi_pair,
         // Address to map with list of unit tensor labels with key N0 
         std::map< int,std::vector<spncci::UnitTensor>>& unit_sym_map,
         // Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
         std::map<
         std::pair<int,int>,
         std::map< spncci::UnitTensorU3Sector,Eigen::MatrixXd >
         >& unit_tensor_rme_map
         );
} //namespace 

#endif

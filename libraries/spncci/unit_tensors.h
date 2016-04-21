/****************************************************************
  unit_tensor.h

  Unit tensor algorithms
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/15/16 (aem,mac): Created.
****************************************************************/
#ifndef UNIT_TENSOR_H_
#define UNIT_TENSOR_H_

#include "spncci/sp_basis.h"
#include "sp3rlib/vcs.h"
#include <unordered_set>

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
    typedef std::tuple<
      u3::U3, HalfInt, HalfInt, 
      int, HalfInt, HalfInt, 
      int, HalfInt, HalfInt> KeyType;
    // w0, S0, T0, rbp, Sb, Tb, r, S, T

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

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
      // hash of each parameter of KeyType
      // omega, S0, T0, rp, Sp, Tp, r, S, T
      // 
      // r and rp need 6 bits
      // TwiceValue(S0), TwiceValue(T0), etc., need only 2 bits 
      // See constants defined just above...

      return 0;
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
      // redundant function??
      // hash omegap_, omega_, tensor_, rho0_
      // hash_combine omegap_, omega_, tensor_, rho0_

      return 0;
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

  void UnitTensorMatrixGenerator(
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

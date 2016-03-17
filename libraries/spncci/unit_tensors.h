/****************************************************************
  unit_tensor.h

  Unit tensor algorithms
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/15/16 (aem,mac): Created.
****************************************************************/
#ifndef UNIT_TENSOR_H_
#define UNIT_TENSOR_H_


#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Eigenvalues>  
#include <map>
#include <tuple>
#include <vector>  

#include "am/am.h"
#include "sp3rlib/u3.h"
#include "spncci/sp_basis.h"


namespace u3
{
	struct UnitTensor
	{
		////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

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
      : S0(0), T0(0), Sp(0), Tp(0), S(0), T(0) {}
    
    // construction from (lambda,mu)
    inline UnitTensor(u3::U3 omega0_, HalfInt S0_, HalfInt T0_, 
		int rp_, HalfInt Sp_, HalfInt Tp_, 
		int r_, HalfInt S_, HalfInt T_)
		: omega0(omega0_),S0(S0_),T0(T0_), 
		rp(rp_), Sp(Sp_), Tp(Tp_), 
		r(r_), S(S_), T(T_) {}


    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      return KeyType(omega0,S0,T0,rp,Sp,Tp,r,S,T);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    // Operator labels
    u3::U3 omega0;
    HalfInt S0,T0,Sp,Tp,S,T;
    int rp,r;

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

struct UnitTensorRME
	{
		////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

		typedef std::tuple<	u3::U3, u3::U3, u3::UnitTensor, int> KeyType;
		// w'S', wS, Unit tensor labels,r0

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    
    //inline UnitTensor() 
    //  : {}
    
    // construction from (lambda,mu)
    inline UnitTensorRME(u3::U3 omegap_, u3::U3 omega_, u3::UnitTensor tensor_, int rho0_)
		: omegap(omegap_),omega(omega_), tensor(tensor_), rho0(rho0_) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      return KeyType(omegap,omega,tensor,rho0);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    // Operator labels
    u3::U3 omegap,omega;
  	u3::UnitTensor tensor;
    int rho0;

  };

  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const UnitTensorRME& unit1rme, const UnitTensorRME& unit2rme)
  {
    return unit1rme.Key() == unit2rme.Key();
  }

  inline bool operator < (const UnitTensorRME& unit1rme, const UnitTensorRME& unit2rme)
  {
    return unit1rme.Key() < unit2rme.Key();
  }

} //namespace 


void UnitSymmetryGenerator(int N1b, int Nmax, std::map< int,std::vector<u3::UnitTensor> >& unit_sym_map);
// Generates a map containing (key, value) pair (N0, operator_labels) of the unit tensors 


void UnitTensorMatrix(int Nshell, int Nsigma_0, int Nmax, std::string LGI_filename);
//
// 












#endif

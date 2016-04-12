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

      
      inline UnitTensorRME() 
      :rho0(0) {}
  

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


  void UnitSymmetryGenerator(int Nmax, std::map< int,std::vector<u3::UnitTensor> >& unit_sym_map);
  // Generates a map containing (key, value) pair (N0, operator_labels) of the unit tensors 
  //
  // 

Eigen::MatrixXd UnitTensorMatrix(
  // LGI pair sector 
  const std::pair<int,int> lgi_pair,
  // sigma' irrep
  const sp3r::Sp3RSpace& irrepp,
  // sigma irrep
  const sp3r::Sp3RSpace& irrep,
  // unit tensor labels 
  u3::UnitTensorRME unit_labels,
   // Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
   std::map<
     std::pair<int,int>, 
     std::map<u3::UnitTensorRME,Eigen::MatrixXd> 
     >& unit_tensor_rme_map
  );

  // inline void TestFunction(  
  //   // boson number cutoff
  //   int Nmax, 
  //   // a given spncci sector pair given as index pair  from global list lgi_vector 
  //   std::pair<int,int> lgi_pair,
  //   // Address to map with list of unit tensor labels with key N0 
  //   std::map< int,std::vector<u3::UnitTensor>>& unit_sym_map,
  //   // Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
  //   std:: map< std::pair<int,int>, std::map< u3::UnitTensorRME,Eigen::MatrixXd> >& unit_tensor_rme_map
  //   )
  //   {
  //     U3 sp=lgi_vector[lgi_pair.first].sigma;
  //     U3 s=lgi_vector[lgi_pair.second].sigma;
  //     for (int i=0; i<unit_sym_map.size(); i++)
  //       {
  //         if (u3::OuterMultiplicity(s.SU3(),unit_sym_map[0][i].omega0.SU3(),sp.SU3())>0)
  //           std::cout<< unit_tensor_rme_map[lgi_pair][u3::UnitTensorRME(sp,s,unit_sym_map[0][i],1)]<<std::endl;
  //       }
  //   }

  void UnitTensorMatrixGenerator(
    int N1b,
  // boson number cutoff
    int Nmax, 
  // a given spncci sector pair given as index pair  from global list lgi_vector 
    std::pair<int,int> lgi_pair,
  // Address to map with list of unit tensor labels with key N0 
    std::map< int,std::vector<u3::UnitTensor>>& unit_sym_map,
  // Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
    std::map<
      std::pair<int,int>,
      std::map< u3::UnitTensorRME,Eigen::MatrixXd >
      >& unit_tensor_rme_map
    );
} //namespace 












#endif

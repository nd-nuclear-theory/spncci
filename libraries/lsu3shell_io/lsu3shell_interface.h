/****************************************************************
  lsu3shell_basis.h

  Interface for lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
****************************************************************/

#ifndef LSU3SHELL_INTERFACE_H_
#define LSU3SHELL_INTERFACE_H_

#include <boost/functional/hash_fwd.hpp>
#include <eigen3/Eigen/Eigen>

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "u3shell/relative_operator.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/u3spn_scheme.h"

namespace lsu3shell
{

class LSU3Irrep
{
	public:

		//default constructor
		inline LSU3Irrep() :Nex_(-999) {}

		//constructor
		inline LSU3Irrep(int Nex, const u3::SU3& x, const HalfInt& Sp, const HalfInt& Sn,const HalfInt& S)
		: Nex_(Nex), x_(x), Sp_(Sp), Sn_(Sn), S_(S){}
		
		//////////////////////////////////////////////////////////////
		//accessors
		//////////////////////////////////////////////////////////////

		u3::SU3 x() const {return x_;}
		int Nex() const {return Nex_;}
		HalfInt S() const {return S_;}
		HalfInt Sp() const {return Sp_;}
		HalfInt Sn() const {return Sn_;}

		//////////////////////////////////////////////////////////////
		//key tuple, comparisons and hashing
		//////////////////////////////////////////////////////////////
		typedef std::tuple<int,u3::SU3, HalfInt, HalfInt,HalfInt> KeyType;
		inline KeyType Key() const
			{
				return KeyType(Nex_, x_, Sp_,Sn_,S_);
			}
		inline friend bool operator == (const LSU3Irrep& state1, const LSU3Irrep& state2)
			{
				return state1.Key()==state2.Key();
			}
		inline friend bool operator < (const LSU3Irrep& state1, const LSU3Irrep& state2)
			{
				return state1.Key()<state2.Key();
			}
		inline friend std::size_t hash_value(const LSU3Irrep& v)
			{
				boost::hash<LSU3Irrep::KeyType> hasher;
				return hasher(v.Key());
			}

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

	private:
		int Nex_;
		u3::SU3 x_;
		HalfInt S_,Sp_,Sn_;
};

typedef std::vector<lsu3shell::LSU3Irrep>	LSU3Vector;

// void ReadLSU3Vector(const std::string& filename, LSU3Vector& lsu3basis_vector);
// reads in lsu3 basis irreps from file and stores them in a vector

void GenerateLSU3ShellOperators(int Nmax, const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_tensor_expansion, std::string filename);
//Generate input files for LSU3shell recoupler for a relative operator
// Nmax gives the truncation for the space on which the relative unit
// tensors expansion is defined
// relative_tensor_expansion is unit tensor expansion of operator
// operator_index gives index of operator file operator.00000index.recoupler

void GenerateLSU3ShellOperators(int Nmax, const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_tensor_labels);
	// Generate input files for LSUshell recoupler for all relative unit tensors 
	// tensor's which may have non-zero matrix elements between  LGI's. 

void 
GenerateLSU3ShellOperators(
		int Nmax, 
		const u3shell::TwoBodyUnitTensorCoefficientsU3ST& twobody_tensor_expansion,
		int operator_index
	);
typedef std::tuple<u3shell::U3SPN,int,int>LSU3BasisGroup;
typedef std::vector<LSU3BasisGroup> LSU3BasisTable;

void 
ReadLSU3Vector(
    HalfInt Nsigma_0, 
    const std::string& filename, 
    LSU3BasisTable& lsu3_basis_table, 
    std::map<u3shell::U3SPN,int>& subspace_dimensions
  );

void GenerateLSU3ShellExpansionLGI(
		int Nsigma_0,
		int Nsigma_min,
		int Nsigma_max, 
		std::string basis_file, 
		std::string brel_filename, 
		std::string nrel_filename,
		std::string lgi_filename,
		std::string lgi_expansion_filename
	);

}
#endif
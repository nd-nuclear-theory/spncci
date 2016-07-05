/****************************************************************
  u3coef.h

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.

  u3coef recquires initilization of the coeffcients by calling the function 
  u3::U3CoefInit() at the start of main

                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created based on prototype u3.py and 
  T. Dytrych CSU3Master.


****************************************************************/

#ifndef U3COEF_H_
#define U3COEF_H_

#include <unordered_map>
#include <tuple>
#include <boost/functional/hash_fwd.hpp>

#include "sp3rlib/u3.h"


namespace u3
{

  ////////////////////////////////////////////////////////////////
  // direct access to su3lib FORTRAN subroutines
  ////////////////////////////////////////////////////////////////

  namespace su3lib
  {

    const size_t MAX_K = 9;

    // Subroutines of original Fortran SU(3) library
    extern "C" 
    { 
      extern void wu3r3w_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[MAX_K][MAX_K][MAX_K][MAX_K]);
      extern void wru3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
      extern void wzu3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
      extern void wu39lm_(const int&, const int& , const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, double[], const int&);
      extern void blocks_(void);
    }

  } //namespace

  ////////////////////////////////////////////////////////////////
  // wrapper functions for single-coefficient access
  ////////////////////////////////////////////////////////////////

  void U3CoefInit();

  double W(const u3::SU3& x1, int k1, int L1, const u3::SU3& x2, int k2, int L2, const u3::SU3& x3, int k3, int L3, int r0);
  // Compute SU(3) reduced coupling coefficient, referred to as Wigner coefficient 
  //
  // Provides wrapper for su3lib function wu3r3w_


  typedef std::tuple<int,int,int,int> UMultiplicityTuple;
  u3::UMultiplicityTuple UMultiplicity(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x,
                                       const u3::SU3& x3, const u3::SU3& x12, const u3::SU3& x23);
  // Calculate multiplicities for SU(3) Racah U and Z coefficients.
  //
  // Arguments:
  //   x1, x2, ... (u3::SU3): SU3 labels for recoupling coefficient
  //
  // Returns:
  //   (UMultiplicityTuple): tuple of multiplicities (r12_max,r12_3_max,r23_max,r1_23_max)

  typedef std::tuple<int,int,int,int> WMultiplicityTuple;
  u3::WMultiplicityTuple WMultiplicity(const u3::SU3& x1, int L1, const u3::SU3& x2, int L2,
                                       const u3::SU3& x3, int L3);
  // Calculate multiplicities for SU(3) Wigner coefficients.
  //
  // Arguments:
  //   x1, x2,x3 (u3::SU3): SU3 labels for coupling coefficient
  //   L1,L2, L3 (int): SO(3) labels for coupling coefficient
  //  
  // Returns:
  //   (WMultiplicityTuple): tuple of multiplicities (kappa1_max,kappa2_max,kappa3_max,rho_max)


  enum class UZMode {kU, kZ};
  double UZ(
           const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3,
           const u3::SU3& x12, int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23,
           UZMode mode
            );
  // Calculate SU(3) Racah recoupling coefficient (six-lm symbol).
  //
  // Provides wrapper for su3lib function wru3optimized_ or wzu3optimized_.
  //
  // Arguments:
  //   x1, x2, ... (u3::SU3): SU3 labels for recoupling coefficient
  //   r12, ... (int): multiplicity (rho) labels for recoupling coefficient
  //   mode (UZMode): kU for U coefficient [(1x2)x3 to 1x(2x3)], kZ for Z coefficient [(1x2)x3 to 2x(1x3)]
  //
  // Returns:
  //   (double): value of coefficient

  inline
  double U(
           const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3,
           const u3::SU3& x12, int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23
           )
  // Compute U coefficient.  See comment for UZ above.
  {
    return u3::UZ(x1,x2,x,x3,x12,r12,r12_3,x23,r23,r1_23,UZMode::kU);
  }

  inline
  double Z(
           const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3,
           const u3::SU3& x12, int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23
           )
  // Compute Z coefficient.  See comment for UZ above.
  {
    return UZ(x1,x2,x,x3,x12,r12,r12_3,x23,r23,r1_23,UZMode::kZ);
  }

  double Phi(const u3::SU3& x1,  const u3::SU3& x2,  const u3::SU3& x3, int r, int rp);
  // Compute Phi phase factor that arrises in chainging the coupling order of SU(3) irreps 

  double Unitary9LambdaMu(
                          const u3::SU3& x1,  const u3::SU3& x2,  const u3::SU3& x12, int r12,
                          const u3::SU3& x3,  const u3::SU3& x4,  const u3::SU3& x34, int r34,
                          const u3::SU3& x13, const u3::SU3& x24, const u3::SU3& x,   int r13_24,
                          int r13,     int r24,     int r12_34
                          );
  // Compute SU(3) unitary 9-(lambda,mu) symbol.
  //
  // Provides wrapper for su3lib function wu39lm_


  ////////////////////////////////////////////////////////////////
  // block storage of coefficients
  ////////////////////////////////////////////////////////////////

  class UCoefLabels
  // Class to gather and provide hashing for U coefficient labels
  {
  public:

    ////////////////////////////////////////////////////////////////
    // type definitions
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<u3::SU3,u3::SU3,u3::SU3,u3::SU3,u3::SU3,u3::SU3> KeyType;
    // tuple of SU(3) labels

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    inline UCoefLabels(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x,
                      const u3::SU3& x3, const u3::SU3& x12, const u3::SU3& x23)
      :x1_(x1), x2_(x2), x_(x), x3_(x3), x12_(x12), x23_(x23){}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
 
    inline KeyType Key() const
    {
      return KeyType(x1_, x2_, x_, x3_, x12_, x23_);
    }

    ////////////////////////////////////////////////////////////////
    // validation
    ////////////////////////////////////////////////////////////////
    
    inline bool Allowed() const
    // Checks if labels satisfy coupling constraints.
    {
      int r12_max, r12_3_max, r23_max, r1_23_max;
      std::tie(r12_max,r12_3_max,r23_max,r1_23_max) = UMultiplicity(x1_,x2_,x_,x3_,x12_,x23_);
      int r_max=r12_max*r12_3_max*r23_max*r1_23_max;
      return (r_max > 0);
    }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    inline friend std::size_t hash_value(UCoefLabels const& ucoef_labels)
    {
      // 6 labels all of SU3 type
      // SU3 type and size: ints lambda and mu in sp3rlib/u3.h
      // SU3 hash_value defined in sp3rlib/u3.h
      std::size_t ulab1 = hash_value(ucoef_labels.x1_);
      std::size_t ulab2 = hash_value(ucoef_labels.x2_);
      std::size_t ulab3 = hash_value(ucoef_labels.x_);
      std::size_t ulab4 = hash_value(ucoef_labels.x3_);
      std::size_t ulab5 = hash_value(ucoef_labels.x12_);
      std::size_t ulab6 = hash_value(ucoef_labels.x23_);

      #ifdef NAIVEHASH
      // naive implementation: sum hashes together
      std::size_t ulabsum = 0;
      return ulabsum = ulab1 + ulab2 + ulab3 + ulab4 +ulab5 + ulab6;
      #endif

      #ifdef BOOSTHASH
      // smarter implementation: boost::hashcombine
      std::size_t seed = 0;
      boost::hash_combine(seed,ulab1);
      boost::hash_combine(seed,ulab2);
      boost::hash_combine(seed,ulab3);
      boost::hash_combine(seed,ulab4);
      boost::hash_combine(seed,ulab5);
      boost::hash_combine(seed,ulab6);
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
    u3::SU3 x1_, x2_, x_, x3_, x12_, x23_;

  };

  inline bool operator == (const UCoefLabels& coef1, const UCoefLabels& coef2)
  {
    return coef1.Key() == coef2.Key();
  }

  inline bool operator < (const UCoefLabels& coef1, const UCoefLabels& coef2)
  {
    return coef1.Key() < coef2.Key();
  }

  class UCoefBlock
  // Class to store and retrieve block of U coefficients sharing same
  // SU(3) labels but with different multiplicity indices
  //
  // FUTURE: Can generalize to store blocks of Z coefficients as well,
  // like function UZ above.
  {
  public:

    ////////////////////////////////////////////////////////////////
    // type definitions
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<int,int,int,int> KeyType;
    // tuple of multiplicities

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    inline UCoefBlock()
    : r12_max_(0), r12_3_max_(0), r23_max_(0), r1_23_max_(0){}
    // Construct and store multiplicites and coefficient values

    UCoefBlock(const u3::UCoefLabels& labels);

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      return KeyType(r12_max_,r12_3_max_,r23_max_,r1_23_max_);
    }
 
    ////////////////////////////////////////////////////////////////
    // entry lookup
    ////////////////////////////////////////////////////////////////

    double GetCoef(int r12, int r12_3, int r23, int r1_23) const;

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
  
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

  private:
    // multiplicities
    int r12_max_, r12_3_max_, r23_max_, r1_23_max_;
    // coefficient values
    std::vector<double> coefs_;

  }; //end UCoefLabels

  ////////////////////////////////////////////////////////////////
  // U coefficient caching
  ////////////////////////////////////////////////////////////////

  typedef std::unordered_map<
    u3::UCoefLabels,
    u3::UCoefBlock,
    boost::hash<u3::UCoefLabels> > UCoefCache;

  extern bool g_u_cache_enabled;
  double UCached(
                 UCoefCache& cache, 
                 const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3, const u3::SU3& x12,
                 int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23
                 );
  // Cached SU(3) Racah recoupling coefficient for recoupling from (1x2)x3 to 1x(2x3). 
  //
  // Global:
  //
  //   u3::g_u_cache_enabled (bool): mode flag determining whether to
  //     use caching or calculate on the fly (for debugging and
  //     profiling)
  //
  // Arguments:
  //   cache (UCoefCache): cache to use for U coefficients
  //   x1, ...: standard U coefficient SU(3) and multiplicity labels
  //
  // Returns;
  //   (double): single coefficient value


  class WCoefLabels
  // Class to gather and provide hashing for U coefficient labels
  {
  public:

    ////////////////////////////////////////////////////////////////
    // type definitions
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<u3::SU3,int,u3::SU3,int,u3::SU3,int> KeyType;
    // tuple of SU(3) and SO(3) labels

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    inline WCoefLabels(const u3::SU3& x1, int L1, const u3::SU3& x2, int L2,
                      const u3::SU3& x3, int L3)
      :x1_(x1), L1_(L1),x2_(x2), L2_(L2), x3_(x3),L3_(L3){}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
 
    inline KeyType Key() const
    {
      return KeyType(x1_, L2_, x2_, L2_, x3_,  L3_);
    }

    ////////////////////////////////////////////////////////////////
    // validation
    ////////////////////////////////////////////////////////////////
    
    // inline bool Allowed() const
    // // Checks if labels satisfy coupling constraints.
    // {
    //   int rho_max;
    //   std::tie(r12_max,r12_3_max,r23_max,r1_23_max) = UMultiplicity(x1_,x2_,x_,x3_,x12_,x23_);
    //   int r_max=r12_max*r12_3_max*r23_max*r1_23_max;
    //   return (r_max > 0);
    // }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    inline friend std::size_t hash_value(WCoefLabels const& wcoef_labels)
    {
      // 6 labels all of SU3 type
      // SU3 type and size: ints lambda and mu in sp3rlib/u3.h
      // SU3 hash_value defined in sp3rlib/u3.h
      std::size_t wlab1 = hash_value(wcoef_labels.x1_);
      std::size_t wlab2 = hash_value(wcoef_labels.x2_);
      std::size_t wlab3 = hash_value(wcoef_labels.x3_);
      // std::size_t wlab4 = hash_value(wcoef_labels.L1_);
      // std::size_t wlab5 = hash_value(wcoef_labels.L2_);
      // std::size_t wlab6 = hash_value(wcoef_labels.L3_);

      #ifdef NAIVEHASH
      // naive implementation: sum hashes together
      std::size_t wlabsum = 0;
      return wlabsum = wlab1 + wlab2 + wlab3 + wlab4 +wlab5 + wlab6;
      #endif

      #ifdef BOOSTHASH
      // smarter implementation: boost::hashcombine
      std::size_t seed = 0;
      boost::hash_combine(seed,wlab1);
      boost::hash_combine(seed,wlab2);
      boost::hash_combine(seed,wlab3);
      boost::hash_combine(seed,wcoef_labels.L1_);
      boost::hash_combine(seed,wcoef_labels.L2_);
      boost::hash_combine(seed,wcoef_labels.L3_);
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
    u3::SU3 x1_, x2_, x3_;
    int L1_, L2_, L3_;
  };

  inline bool operator == (const WCoefLabels& coef1, const WCoefLabels& coef2)
  {
    return coef1.Key() == coef2.Key();
  }

  inline bool operator < (const WCoefLabels& coef1, const WCoefLabels& coef2)
  {
    return coef1.Key() < coef2.Key();
  }

  class WCoefBlock
  // Class to store and retrieve block of U coefficients sharing same
  // SU(3) labels but with different multiplicity indices
  //
  // FUTURE: Can generalize to store blocks of Z coefficients as well,
  // like function UZ above.
  {
  public:

    ////////////////////////////////////////////////////////////////
    // type definitions
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<int,int,int,int> KeyType;
    // tuple of multiplicities

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    inline WCoefBlock()
    : kappa1_max_(0), kappa2_max_(0), kappa3_max_(0),rho_max_(0){}
    // Construct and store multiplicites and coefficient values

    WCoefBlock(const u3::WCoefLabels& labels);

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      return KeyType(kappa1_max_,kappa2_max_,kappa3_max_,rho_max_);
    }
 
    ////////////////////////////////////////////////////////////////
    // entry lookup
    ////////////////////////////////////////////////////////////////

    double GetCoef(int kappa1, int kappa2, int kappa3, int rho) const;

    inline std::vector<double> GetCoefBlock() const
    {
      return coefs_;
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
  
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

  private:
    // multiplicities
    int kappa1_max_,kappa2_max_,kappa3_max_,rho_max_;
   

    // coefficient values
    std::vector<double> coefs_;// May need to be changed...

  };

  ////////////////////////////////////////////////////////////////
  // W coefficient caching
  ////////////////////////////////////////////////////////////////

  typedef std::unordered_map<
    u3::WCoefLabels,
    u3::WCoefBlock,
    boost::hash<u3::WCoefLabels> > WCoefCache;

  extern bool g_u_cache_enabled;
  double WCached(
                 WCoefCache& cache, 
                 const u3::SU3& x1, int kappa1, int L1, const u3::SU3& x2, int kappa2, int L2, 
                 const u3::SU3& x3, int kappa3, int L3, int rho
                 );
  // Cached SU(3) Wigner coupling coefficient for coupling from (1x2)->3. 
  //
  // Global:
  //
  //   u3::g_u_cache_enabled (bool): mode flag determining whether to
  //     use caching or calculate on the fly (for debugging and
  //     profiling)
  //
  // Arguments:
  //   cache (WCoefCache): cache to use for W coefficients
  //   x1,x2,x3 ...: standard W coefficient SU(3) and multiplicity labels
  //
  // Returns;
  //   (double): single coefficient value




} //namespace 






#endif

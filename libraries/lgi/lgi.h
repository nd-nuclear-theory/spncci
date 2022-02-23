/****************************************************************
  lgi.h

  Interface for lsu3shell basis.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  09/08/16 (aem,mac): Created.
  01/31/17 (mac): Rename LGIVector to MultiplicityTaggedLGIVector.
  02/17/17 (mac): Extract WriteLGILabels from lgi_solver and change
    to accept std::ostream for output and extract to lgi.
  10/11/17 (aem) : Extract Nsigma0ForNuclide from spncci_basis
  01/15/18 (aem) : Updated Read LGI to be consistant with write functions
  01/28/21 (aem) : Removed deprecated WriteLGIExpansion function.
  08/06/21 (pjf): Add UpstreamLabelsType for use with new recurrence indexing.
****************************************************************/
#ifndef LGI_H_
#define LGI_H_

#include <tuple>

#ifdef BASIS_HASH
#include <boost/container_hash/hash.hpp>
#endif  // BASIS_HASH

#include "am/am.h"
#include "sp3rlib/sp3r.h"
#include "u3shell/u3spn_scheme.h"
#include "lsu3shell/lsu3shell_rme.h"
namespace lgi
{
  ////////////////////////////////////////////////////////////////
  // LGI labels
  ////////////////////////////////////////////////////////////////

  // class LGI
  //
  // U3SPN labels "dressed" with easy access to the Nex quantum number
  //
  // U3SPN (parent): omega x Sp x Sn x S labels

  class LGI
    : public u3shell::U3SPN
  {
    public:

    struct UpstreamLabelsType { HalfInt Sp; HalfInt Sn; };

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying
    // null constructor

    inline LGI():Nex_(-999){};
    inline
      LGI(const u3shell::U3SPN& sigmaSPN, int Nex)
      : U3SPN(sigmaSPN), Nex_(Nex)
    {}
    //////////////////////////////////////////////////////////////
    //accessors
    //////////////////////////////////////////////////////////////

    // Note: inherits all U3SPN accessors

    int Nex() const {return Nex_;}
    u3::U3 sigma() const {return U3();}
    //U3SPN could conflict with u3shell::U3SPN
    u3shell::U3SPN u3spn() const {return {U3S(),Sp(),Sn()}; }

    UpstreamLabelsType upstream_labels() const { return {Sp(),Sn()}; }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////

    std::string Str() const;
    //////////////////////////////////////////////////////////////
    //key tuple, comparisons and hashing
    //////////////////////////////////////////////////////////////

    // basic ordering and hashing functions inherited from U3SPN

    // pseudo-key for std::tie access to LGI properties
    //
    // This Key is *not* the key actually used in sorting and hashing,
    // which is the parent U3SPN type's key.

    typedef std::tuple<int,u3::U3,HalfInt,HalfInt,HalfInt> KeyType;

    inline KeyType Key() const
    {
      return {Nex(),sigma(),Sp(),Sn(),S()};
    }


    ////////////////////////////////////////////////////////////////
    // data
    ////////////////////////////////////////////////////////////////
    private:

    // Note: U3SPN labels are inherited.

    // quick-reference information
    int Nex_;
  };

  #ifdef BASIS_HASH
  inline std::size_t hash_value(const LGI::UpstreamLabelsType& l)
  {
    return boost::hash<std::tuple<HalfInt,HalfInt>>{}({l.Sp,l.Sn});
  }

  inline std::size_t hash_value(const LGI& l)
  {
    return boost::hash<std::tuple<u3::U3S,LGI::UpstreamLabelsType>>{}({l.U3S(),l.upstream_labels()});
  }

  #endif  // BASIS_HASH

  inline bool operator<(const LGI::UpstreamLabelsType& lhs, const LGI::UpstreamLabelsType& rhs)
  {
    return std::tie(lhs.Sp,lhs.Sn) < std::tie(rhs.Sp,rhs.Sn);
  }
  inline bool operator==(const LGI::UpstreamLabelsType& lhs, const LGI::UpstreamLabelsType& rhs)
  {
    return std::tie(lhs.Sp,lhs.Sn) == std::tie(rhs.Sp,rhs.Sn);
  }
  inline bool operator>(const LGI::UpstreamLabelsType& lhs, const LGI::UpstreamLabelsType& rhs)
  {
    return std::tie(lhs.Sp,lhs.Sn) > std::tie(rhs.Sp,rhs.Sn);
  }
  inline bool UpstreamLabelsAllowed(
      const lgi::LGI::UpstreamLabelsType& ket_labels,
      const lgi::LGI::UpstreamLabelsType& bra_labels
    )
  {
    return (
        abs(ket_labels.Sp - bra_labels.Sp) <= 2
        && abs(ket_labels.Sn - bra_labels.Sn) <= 2
      );
  }


  ////////////////////////////////////////////////////////////////
  // enumeration of LGI set based on input table
  ////////////////////////////////////////////////////////////////

  // LGI input tabulation format:
  //
  //   Nex 2N lambda mu 2Sp 2Sn 2S count
  //   ...
  //
  // Each input table line results in multiple stored LGIs based on
  // the given count.
  //
  // LGI container convenience type
  //
  // STYLE: maybe LGI::vector would be more consistent
  typedef MultiplicityTagged<lgi::LGI>::vector MultiplicityTaggedLGIVector;

  std::string LGIFamilyStr(const MultiplicityTagged<lgi::LGI>& lgi_family);

  void
    WriteLGILabels(const lgi::MultiplicityTaggedLGIVector& lgi_families,std::ostream& os);

  void
    WriteLGILabels(const lgi::MultiplicityTaggedLGIVector& lgi_families, const std::string& filename);

  void
  ReadLGISet(
    const std::string& lgi_filename,
    const HalfInt& Nsigma0,
    MultiplicityTaggedLGIVector& lgi_vector
  );
  // Read in LGI from file and create vector of LGIs tagged by
  // gamma_max from tabulation in file .
  //
  // Arguments:
  //   filename (string) : filename for LGI table file
  //   Nsigma0 : minimum number of oscillator quanta for the given system of nucleons.
  //    Can be obtained from lgi::Nsigma0ForNuclide.
  //   lgi_families (MultiplicityTaggedLGIVector) : container for LGI list (OUTPUT)

 void ReadLGILookUpTable(std::vector<int>& lgi_full_space_lookup_table, int num_irrep_families);
  // Reading in and filling out table of lgi indices in basis and lgi indices in full space by
  // which the seed files are labeled.



}  // namespace lgi

namespace fmt
{
template <>
struct formatter<lgi::LGI> 
{
  char presentation = 'g';

  template <typename ParseContext>
  FMT_CONSTEXPR auto parse(ParseContext& ctx) -> decltype(ctx.begin()) 
  {
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && (*it == 'd' || *it == 'g' || *it == 'f')) presentation = *it++;

    // Check if reached the end of the range:
    if (it != end && *it != '}')
      throw format_error("invalid format");

    // Return an iterator past the end of the parsed range:
    return it;
  }

  template <typename FormatContext>
  FMT_CONSTEXPR auto format(const lgi::LGI& lgi, FormatContext& ctx) -> decltype(ctx.out()) 
  {
    if(presentation == 'f')
      return format_to(ctx.out(), "{:d}{:d}[{:f},{:f},{:f}]",lgi.Nex(),lgi.SU3(),lgi.Sp(),lgi.Sn(),lgi.S());
    else
    return format_to(ctx.out(), "{:d}{:d}[{:d},{:d},{:d}]",lgi.Nex(),lgi.SU3(),lgi.Sp(),lgi.Sn(),lgi.S());
  }
};


}// namespace fmt


#ifdef BASIS_HASH
namespace std
{
  template<> struct hash<typename lgi::LGI::UpstreamLabelsType>
  {
    inline std::size_t operator()(const lgi::LGI::UpstreamLabelsType& h) const noexcept
    {
      return boost::hash<lgi::LGI::UpstreamLabelsType>{}(h);
    }
  };

  template<> struct hash<typename lgi::LGI>
  {
    inline std::size_t operator()(const lgi::LGI& h) const noexcept
    {
      return boost::hash<lgi::LGI>{}(h);
    }
  };
}  // namespace std
# endif // BASIS_HASH

#endif

/****************************************************************
  tensor_labels.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
  
#include "u3shell/tensor_labels.h"

#include <sstream>

#include "cppformat/format.h"

namespace u3shell {

  std::string RelativeStateLabelsU3ST::Str() const
  {

    return fmt::format(
        "[{} {} {}]",
        eta(),S(),T()
      );
  }

    std::string RelativeCMStateLabelsU3ST::Str() const
  {

    return fmt::format(
        "[[{} {}]{} {} {}]",
        eta_r(),eta_cm(),x().Str(),S(),T()
      );
  }

  std::string TwoBodyStateLabelsU3ST::Str() const
  {

    return fmt::format(
        "[[{},{}]{} {} {}]",
        eta1(),eta2(),x().Str(),S(),T()
      );
  }

  std::string OperatorLabelsU3ST::Str() const
  {

    return fmt::format(
        "[{:+d}{} {} {} {}]",
        N0(),x0().Str(),S0(),T0(),g0()
      );
  }

  std::string OperatorLabelsU3S::Str() const
  {

    return fmt::format(
        "[{:+d}{} {} {}]",
        N0(),x0().Str(),S0(),g0()
      );
  }


  std::string RelativeIndexedTensorLabelsU3ST::Str() const
  {

    return fmt::format(
        "{}:{} {}",
        OperatorLabelsU3ST::Str(),
        kappa0(),L0()
      );
  }


  std::string RelativeUnitTensorLabelsU3ST::Str() const
  {

    return fmt::format(
        "Unit{}({},{})",
        OperatorLabelsU3ST::Str(),
        bra().Str(),ket().Str()
      );
  }
 
  std::string RelativeCMUnitTensorLabelsU3ST::Str() const
  {

    return fmt::format(
        "Unit{}{}({},{})",
        OperatorLabelsU3ST::Str(),rho0(),
        bra().Str(),ket().Str()
      );
  }

  std::string TwoBodyUnitTensorLabelsU3ST::Str() const
  {

    return fmt::format(
        "Unit{}{}({},{})",
        OperatorLabelsU3ST::Str(),
        rho0(),
        bra().Str(),ket().Str()
      );
  }

  std::string TwoBodyStateLabelsU3S::Str() const
  {

    return fmt::format(
        "[[{},{}]{} {}]",
        eta1(),eta2(),x().Str(),S()
      );
  }

  std::string TwoBodyUnitTensorLabelsU3S::Str() const
  {

    return fmt::format(
        "Unit{}{}({},{})",
        OperatorLabelsU3S::Str(),
        rho0(),
        bra().Str(),ket().Str()
      );
  }

 
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

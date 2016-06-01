/****************************************************************
  two_body_operator.h

  Two-body operator representation and second-quantization.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/15/16 (mac): Created.
  5/18/16 (mac): Rough in interface.

****************************************************************/

#ifndef TWO_BODY_OPERATOR_H_
#define TWO_BODY_OPERATOR_H_

#include <map>

#include "u3shell/tensor_labels.h"
#include "u3shell/relative_operator.h"

namespace u3shell
{

  ////////////////////////////////////////////////////////////////
  // coefficient storage -- two-body
  ////////////////////////////////////////////////////////////////

  typedef
    std::map<u3shell::TwoBodyUnitTensorLabelsU3ST,double>
    TwoBodyUnitTensorCoefficientsU3ST;
  // Storage of coefficients of two-body unit tensors in U(3)xSxT scheme.
  //
  // Also used to represent coefficients of other similarly-labeled
  // objects, in particular, biquads.

  struct CoefficientsPN
  // Structure to store all three pn-scheme coefficients which arise
  // at the two-body level: pppp, nnnn, and pnnp.
  {
    // constuctors
    inline CoefficientsPN()
      // Default constructor.
      : pppp(0), nnnn(0), pnnp(0)
    {}

    inline CoefficientsPN(double pppp_, double nnnn_, double pnnp_)
      // Construct from coefficients.
      : pppp(pppp_), nnnn(nnnn_), pnnp(pnnp_)
    {}

    // arithmetic assignment
    inline CoefficientsPN& operator +=(const CoefficientsPN& b)
    {
      pppp += b.pppp;
      nnnn += b.nnnn;
      pnnp += b.pnnp;
      return *this;
    }

    // data
    double pppp, nnnn, pnnp;
  };

  typedef
    std::map<u3shell::TwoBodyUnitTensorLabelsU3S,u3shell::CoefficientsPN>
    TwoBodyUnitTensorCoefficientsU3SPN;
  // Storage of all three pn-scheme coefficients for given U(3)xS
  // two-body unit tensor labels.
  //
  // Also used to represent coefficients of other similarly-labeled
  // objects, in particular, biquads.

  ////////////////////////////////////////////////////////////////
  // transformation from relative to two-body unit tensors
  ////////////////////////////////////////////////////////////////

  // TODO integrate with or replace with Anna's Moshinsky work -- this is the Moshinsky
  // xform for an arbitrary operator, if we take the Moshinsky xform
  // for a single relative unit tensor to be our "basic" calculation

  void TransformRelativeUnitTensorToTwoBodyUnitTensor
    (
     const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_unit_tensor_coefficients,
     u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_unit_tensor_coefficients
     );
  // Accumulate biquad coefficients for given two-body unit tensors.
  //
  // relative_unit_tensor_coefficients (RelativeUnitTensorCoefficientsU3ST, input)
  //   : map giving coefficients on a set of two-body U3ST unit tensors 
  // two_body_unit_tensor_coefficients (TwoBodyUnitTensorCoefficientsU3ST, input)
  //   : map giving coefficients on a resulting set of two-body U3ST unit tensors 
  //
  // Note: The map given as the second parameter, for storing the
  // output coefficients, need not initially be empty.  It is
  // permissible to accumulate onto an existing set of coefficients.
  // Coefficients for the same term will be added.

  ////////////////////////////////////////////////////////////////
  // transformation from two-body unit tensors to biquads
  ////////////////////////////////////////////////////////////////

  void TransformTwoBodyUnitTensorToBiquad(
                                          const u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_unit_tensor_coefficients,
                                          u3shell::TwoBodyUnitTensorCoefficientsU3ST& biquad_coefficients
                                          );
  // Accumulate biquad coefficients for given two-body unit tensors.
  //
  // New coefficients accumulate additively onto any existing coefficients in biquad_coefficients.
  //
  // Biquads are labeled still by the corresponding unit tensor
  // labels, and recall that all the unit tensor bra and ket labels
  // are kept in their ket-like (unconjugated) from.  Thus, e.g.,
  // a unit tensor with unit RME of the form
  //
  //   < (eta1p eta2p)^omegap || ... || (eta1 eta2)^omega >
  //
  // is labeled by
  //
  //   bra: (eta1p eta2p)^omegap  ket : (eta1 eta2)^omega
  //
  // These labels thus also label the biquad
  //
  //   (a+_eta1p a+_eta2p)^omegap x (a~_eta2 a~_eta1)^omega~
  //
  //  Note the "reversed" order of factors in the second (omega~)
  //  product.
  //
  // Arguments:
  //   two_body_unit_tensor_coefficients (TwoBodyUnitTensorCoefficientsU3ST, input)
  //     : map giving coefficients on a set of two-body U3ST unit tensors 
  //   biquad_coefficients (TwoBodyUnitTensorCoefficientsU3ST, output)
  //     : map giving coefficients on a resulting set of two-body U3ST biquads

  void TransformBiquadToPNScheme(
                                 const u3shell::TwoBodyUnitTensorCoefficientsU3ST& biquad_coefficients,
                                 u3shell::TwoBodyUnitTensorCoefficientsU3SPN& biquad_coefficients_pn
                                 );
  // Accumulate biquad coefficients for given two-body unit tensors.
  //
  // biquad_coefficients (TwoBodyUnitTensorCoefficientsU3ST, input)
  //   : map giving coefficients on a set of two-body U3ST biquads
  // biquad_coefficients (TwoBodyUnitTensorCoefficientsU3ST, output)
  //   : map giving coefficients on a resulting set of two-body U3SPN biquads
  //
  // Note: The map given as the second parameter, for storing the
  // output coefficients, need not initially be empty.  It is
  // permissible to accumulate onto an existing set of coefficients.
  // Coefficients for the same term will be added.


  ////////////////////////////////////////////////////////////////
  // operator file output for recoupler
  ////////////////////////////////////////////////////////////////

  void WriteTwoBodyOperatorRecoupler(
      std::ostream& output_stream,
      const u3shell::TwoBodyUnitTensorCoefficientsU3SPN& biquad_coefficients_pn
    );
  // Write two-body operator pn-scheme biquad coefficients in format
  // expected by TD's recoupler.
  //
  // Output format:
  //
  //   The expected input for TD's recouple consists of a sequence of
  //   label line / coefficient line pairs.  Any Pauli-forbidden
  //   like-nucleon coefficient should be set to zero.
  //
  //     label line:
  //       eta1p eta2p eta2 eta1
  //         rhop=1 lambdap mup 2Sp
  //         rho=1  lambda  mu  2S
  //         rho0   lambda0 mu0 2S0
  //
  //     coefficient line:
  //       pppp nnnn pnnp
  //
  //   Note the ordering of eta values, determined by the labeling
  //   convention discussed in the comments for
  //   TransformTwoBodyUnitTensorToBiquad.

  // Arguments:
  //   output_stream (std::ofstream) : an open text-mode output stream
  //   biquad_coefficients_pn (u3shell::TwoBodyUnitTensorCoefficientsU3SPN)
  //     : the biquad coefficients

}  // namespace

#endif

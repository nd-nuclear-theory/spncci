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


#include "sp3rlib/u3.h"
#include "u3shell/tensor.h"

namespace u3shell
{

  ////////////////////////////////////////////////////////////////
  // coefficient storage -- relative
  ////////////////////////////////////////////////////////////////

  // TODO move into separate header with relative operator definitions?
  // and perhaps upcoupling definitions?

  typedef
    std::map<u3shell::RelativeUnitTensorLabelsU3ST,double>
    RelativeUnitTensorCoefficientsU3ST;
  // Storage of coefficients of relative unit tensors in U(3)xSxT scheme.

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
    std::map<u3shell::TwoBodyUnitTensorLabelsU3SPN,u3shell::CoefficientsPN>
    TwoBodyUnitTensorCoefficientsU3SPN;
  // Storage of all three pn-scheme coefficients for given U(3)xS
  // two-body unit tensor labels.
  //
  // Also used to represent coefficients of other similarly-labeled
  // objects, in particular, biquads.

  ////////////////////////////////////////////////////////////////
  // transformation from relative to two-body unit tensors
  ////////////////////////////////////////////////////////////////

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

  // TODO integrate with or replace with Anna's Moshinsky work -- this is the Moshinsky
  // xform for an arbitrary operator, if we take the Moshinsky xform
  // for a single relative unit tensor to be our "basic" calculation



  ////////////////////////////////////////////////////////////////
  // transformation from two-body unit tensors to biquads
  ////////////////////////////////////////////////////////////////

  void TransformTwoBodyUnitTensorToBiquad(
                                          const u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_unit_tensor_coefficients,
                                          u3shell::TwoBodyUnitTensorCoefficientsU3ST& biquad_coefficients
                                          );
  // Accumulate biquad coefficients for given two-body unit tensors.
  //
  // two_body_unit_tensor_coefficients (TwoBodyUnitTensorCoefficientsU3ST, input)
  //   : map giving coefficients on a set of two-body U3ST unit tensors 
  // biquad_coefficients (TwoBodyUnitTensorCoefficientsU3ST, output)
  //   : map giving coefficients on a resulting set of two-body U3ST biquads
  //
  // Note: The map given as the second parameter, for storing the
  // output coefficients, need not initially be empty.  It is
  // permissible to accumulate onto an existing set of coefficients.
  // Coefficients for the same term will be added.

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

  // TODO replace with WriteBiquadCoefficientsPN
  class OutRecouplerStream
  // Provide output stream for writing input file to Tomas's recoupler
  // code.
  {

    ////////////////////////////////////////////////////////////////
    // constructor
    ////////////////////////////////////////////////////////////////

  public:

    inline
      OutRecouplerStream()
      : stream_(NULL) {};

    ////////////////////////////////////////////////////////////////
    // destructor
    ////////////////////////////////////////////////////////////////

    // ./libraries/u3shell/two_body_operator.h: In destructor
    // 'u3shell::OutRecouplerStream::~OutRecouplerStream()':

    //./libraries/u3shell/two_body_operator.h:157:16: warning:
    //possible problem detected in invocation of delete operator:
    //[-Wdelete-incomplete]

    //~OutRecouplerStream()
    //  {
    //    delete stream_;
    //  };

    ////////////////////////////////////////////////////////////////
    // I/O activities
    ////////////////////////////////////////////////////////////////

    void Open (const std::string& basename);
    void WriteOperator (const u3shell::TwoBodyUnitTensorCoefficientsU3SPN&);
    void Close ();

  private:
    // data
    std::ofstream* stream_;
  };


}  // namespace

#endif

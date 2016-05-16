/****************************************************************
  two_body.h

  Two-body operator second-quantization.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/15/16 (mac): Created.

****************************************************************/

#ifndef TWO_BODY_H_
#define TWO_BODY_H_

#include "sp3rlib/u3.h"
#include "u3shell/tensor.h"

namespace u3shell
{

//  ////////////////////////////////////////////////////////////////
//  // operator file output
//  ////////////////////////////////////////////////////////////////
//
//  class OutLSUOperatorStream
//  {
//
//    ////////////////////////////////////////////////////////////////
//    // constructor
//    ////////////////////////////////////////////////////////////////
//
//  public:
//
//    inline
//      OutMFDnH2Stream()
//      : stream_ppnn_(0), stream_pn_(0) {};
//
//    ////////////////////////////////////////////////////////////////
//    // destructor
//    ////////////////////////////////////////////////////////////////
//
//    ~OutMFDnH2Stream()
//      {
//        delete stream_ppnn_;
//        delete stream_pn_;
//      };
//
//    ////////////////////////////////////////////////////////////////
//    // I/O activities
//    ////////////////////////////////////////////////////////////////
//
//    void Open (const std::string& basename);
//    void WriteOperator (const & );
//    void Close ();
//
//  private:
//    // data
//    std::ofstream* stream_ppnn_, stream_pn_;
//  };


}  // namespace

#endif

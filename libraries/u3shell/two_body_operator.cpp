/****************************************************************
  two_body_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/two_body_operator.h"

#include <array>

#include "cppformat/format.h"
#include "am/am.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"

namespace u3shell {


  void TransformRelativeUnitTensorToTwoBodyUnitTensor
    (
     const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_unit_tensor_coefficients,
     u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_unit_tensor_coefficients
     )
  {

    for (const auto& key_value : relative_unit_tensor_coefficients)
      {

        // extract key and value
        const u3shell::RelativeUnitTensorLabelsU3ST& relative_unit_tensor_labels = key_value.first;
        double relative_unit_tensor_coefficient = key_value.second;

        // TODO -- or more likely to be replaced by Anna's code
      }
  }


  void TransformTwoBodyUnitTensorToBiquad(
      const u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_unit_tensor_coefficients,
      u3shell::TwoBodyUnitTensorCoefficientsU3ST& biquad_coefficients
    )
  {
    for (const auto& key_value : two_body_unit_tensor_coefficients)
      {

        // extract key and value
        const u3shell::TwoBodyUnitTensorLabelsU3ST& two_body_unit_tensor_labels = key_value.first;
        double two_body_unit_tensor_coefficient = key_value.second;

        // extract unit tensor label groups
        u3shell::OperatorLabelsU3ST::KeyType operator_labels_key;
        int rho0;
        u3shell::TwoBodyStateLabelsU3ST::KeyType bra_key, ket_key;
        std::tie(operator_labels_key,rho0,bra_key,ket_key) = two_body_unit_tensor_labels.Key();

        // extract operator label groups
        int N0;
        u3::SU3 x0;
        HalfInt S0,T0;
        int g0;
        std::tie(N0,x0,S0,T0,g0) = operator_labels_key;
        
        // extract bra labels
        int eta1p, eta2p;
        u3::SU3 xp;
        HalfInt Sp, Tp;
        std::tie(eta1p,eta2p,xp,Sp,Tp) = bra_key;

        // extract ket labels
        int eta1, eta2;
        u3::SU3 x;
        HalfInt S, T;
        std::tie(eta1,eta2,x,S,T) = ket_key;

        // calculate phase and normalization factors
        int cross_projector_grade = ConjugationGrade(x0)+ConjugationGrade(x)+ConjugationGrade(xp);
        double cross_projector_norm_factor
          = sqrt(
                 double(dim(xp)*am::dim(Sp)*am::dim(Tp))
                 / double(dim(x0)*am::dim(S0)*am::dim(T0))
                 );
        int biquad_grade = eta1+eta2+ConjugationGrade(x)+int(S)+int(T);
        double biquad_norm_factor = 1/(
                             sqrt(1+KroneckerDelta(eta1p,eta2p))
                             *sqrt(1+KroneckerDelta(eta1,eta2))
                             );

        // iterate over outer multiplicity on cross projector
        int outer_multiplicity = u3::OuterMultiplicity(x,x0,xp);
        for (int rho0bar = 1; rho0bar <= outer_multiplicity; ++rho0bar)
          {
            // calculate Phi factor
            double recoupling_phase = u3::Phi(x,x0,xp,rho0,rho0bar);

            // calculate total coefficient on biquad
            double biquad_coefficient
              = two_body_unit_tensor_coefficient
              * ParitySign(cross_projector_grade + biquad_grade)
              * recoupling_phase
              * cross_projector_norm_factor
              * biquad_norm_factor;
            
            // assemble biquad labels
            u3shell::TwoBodyUnitTensorLabelsU3ST
              biquad_labels(
                  two_body_unit_tensor_labels,  // use OperatorLabelsU3ST base class instance within two_body_unit_tensor_labels
                  rho0bar,
                  two_body_unit_tensor_labels.bra(),
                  two_body_unit_tensor_labels.ket()
                );

            // accumulate coefficient
            biquad_coefficients[biquad_labels] += biquad_coefficient;
            
          }

      }
  }

  void TransformBiquadToPNScheme(
      const u3shell::TwoBodyUnitTensorCoefficientsU3ST& biquad_coefficients,
      u3shell::TwoBodyUnitTensorCoefficientsU3SPN& biquad_coefficients_pn
    )
  {
    for (const auto& key_value : biquad_coefficients)
      {

        ////////////////////////////////
        // label extraction
        ////////////////////////////////

        // extract key and value
        const u3shell::TwoBodyUnitTensorLabelsU3ST& biquad_labels = key_value.first;
        double biquad_coefficient = key_value.second;

        // extract unit tensor label groups
        u3shell::OperatorLabelsU3ST::KeyType operator_labels_key;
        int rho0;
        u3shell::TwoBodyStateLabelsU3ST::KeyType bra_key, ket_key;
        std::tie(operator_labels_key,rho0,bra_key,ket_key) = biquad_labels.Key();

        // extract operator label groups
        int N0;
        u3::SU3 x0;
        HalfInt S0,T0;
        int g0;
        std::tie(N0,x0,S0,T0,g0) = operator_labels_key;
        
        // extract bra labels
        int eta1p, eta2p;
        u3::SU3 xp;
        HalfInt Sp, Tp;
        std::tie(eta1p,eta2p,xp,Sp,Tp) = bra_key;

        // extract ket labels
        int eta1, eta2;
        u3::SU3 x;
        HalfInt S, T;
        std::tie(eta1,eta2,x,S,T) = ket_key;

        ////////////////////////////////
        // determine pn-scheme terms
        ////////////////////////////////

        // like particle terms (pppp & nnnn):
        //   include_like (bool) : whether or not like-particle biquads are to be included
        //   norm_like (double) : overall normalization factor for like-particle biquads
        //   signs_like (std::array<int,2>) : signs (+1,-1) for like-particle biquads
        //     as {pppp,nnnn}
        //
        // like particle terms active if both bra and ket are T=1
        bool include_like = ((Tp==1)&&(T==1));
        double norm_like;
        std::array<int,2> signs_like;

        if ((Tp==1)&&(T==1)&&(T0==0))
          {
            norm_like = 1/sqrt(3);
            signs_like = {+1,+1};
          }
        else if ((Tp==1)&&(T==1)&&(T0==1))
          {
            norm_like = 1/sqrt(2);
            signs_like = {+1,-1};
          }
        else if ((Tp==1)&&(T==1)&&(T0==2))
          {
            norm_like = 1/sqrt(6);
            signs_like = {+1,+1};
          }

        // unlike particle terms (pnnp):
        //   include_unlike (bool) : whether or not unlike-particle biquads are to be included
        //   norm_unlike (double) : overall normalization factor for unlike-particle biquads
        //   signs_unlike (std::array<int,4>) : signs (+1,-1) for unlike-particle biquads
        //     as {<12|12>, <12|21>, <21|12>, <21|12>}, where recall we label biquads using 
        //     the label ordering as it appears in the corresponding RME
        //
        // unlike particle terms active except in pure isovector (Tp,T,T0)=(1,1,1) special case

        bool include_unlike = !((Tp==1)&&(T==1)&&(T0==1));
        double norm_unlike;
        std::array<int,4> signs_unlike;

        if ((Tp==0)&&(T==0)&&(T0==0))
          {
            norm_unlike = 1/2;
            signs_unlike = {-1,+1,+1,-1};
          }
        else if ((Tp==1)&&(T==0)&&(T0==1))
          {
            norm_unlike = 1/2;
            signs_unlike = {-1,+1,-1,+1};
          }
        else if ((Tp==0)&&(T==1)&&(T0==1))
          {
            norm_unlike = 1/2;
            signs_unlike = {-1,-1,+1,+1};
          }
        else if ((Tp==1)&&(T==1)&&(T0==0))
          {
            norm_unlike = 1/(2*sqrt(3));
            signs_unlike = {+1,+1,+1,+1};
          }
        else if ((Tp==1)&&(T==1)&&(T0==1))
          {
            // pass
          }
        else if ((Tp==1)&&(T==1)&&(T0==2))
          {
            norm_unlike = 1/(sqrt(2*3));
            signs_unlike = {-1,-1,-1,-1};
          }

        ////////////////////////////////
        // accumulate pn-scheme terms
        ////////////////////////////////

        // operator labels sans isospin
        u3shell::OperatorLabelsU3S operator_labels_u3s(N0,x0,S0,g0);

        // bra/ket labels sans isospin
        //
        // including both index orderings for use in unlike-particle terms
        u3shell::TwoBodyStateLabelsU3S statep12(eta1p,eta2p,xp,Sp);
        u3shell::TwoBodyStateLabelsU3S statep21(eta2p,eta1p,xp,Sp);
        u3shell::TwoBodyStateLabelsU3S state12(eta1,eta2,x,S);
        u3shell::TwoBodyStateLabelsU3S state21(eta2,eta1,x,S);
        
        // common working variables
        u3shell::TwoBodyUnitTensorLabelsU3S biquad_labels_pn;
        u3shell::CoefficientsPN coefficients_pn;

        // accumulate like-particle terms
        if (include_like)
          {
            biquad_labels_pn = u3shell::TwoBodyUnitTensorLabelsU3S(
                x0,S0,rho0,statep12,state12
              );
            coefficients_pn = u3shell::CoefficientsPN(
                norm_like*signs_like[0],
                norm_like*signs_like[1],
                0
              );
            biquad_coefficients_pn[biquad_labels_pn] += coefficients_pn;
          }

        // accumulate unlike-particle terms
        if (include_unlike)
          {

            // calculate phase factors
            int grade = eta1+eta2+ConjugationGrade(x)+int(S);
            int gradep = eta1p+eta2p+ConjugationGrade(xp)+int(Sp);

            // term <12|12>
            biquad_labels_pn = u3shell::TwoBodyUnitTensorLabelsU3S(
                operator_labels_u3s,rho0,statep12,state12
              );
            coefficients_pn = u3shell::CoefficientsPN(
                0,0,norm_unlike*signs_unlike[0]
              );
            biquad_coefficients_pn[biquad_labels_pn] += coefficients_pn;

            // term <12|21>
            biquad_labels_pn = u3shell::TwoBodyUnitTensorLabelsU3S(
                operator_labels_u3s,rho0,statep12,state21
              );
            coefficients_pn = u3shell::CoefficientsPN(
                0,0,norm_unlike*signs_unlike[1]*ParitySign(grade)
              );
            biquad_coefficients_pn[biquad_labels_pn] += coefficients_pn;


            // term <21|12>
            biquad_labels_pn = u3shell::TwoBodyUnitTensorLabelsU3S(
                operator_labels_u3s,rho0,statep21,state12
              );
            coefficients_pn = u3shell::CoefficientsPN(
                0,0,norm_unlike*signs_unlike[2]*ParitySign(gradep)
              );
            biquad_coefficients_pn[biquad_labels_pn] += coefficients_pn;

            // term <21|21>
            biquad_labels_pn = u3shell::TwoBodyUnitTensorLabelsU3S(
                operator_labels_u3s,rho0,statep21,state21
              );
            coefficients_pn = u3shell::CoefficientsPN(
                0,0,norm_unlike*signs_unlike[3]*ParitySign(grade+gradep)
              );
            biquad_coefficients_pn[biquad_labels_pn] += coefficients_pn;
          }

      } 
  }
  

  void WriteTwoBodyOperatorRecoupler(
      std::ostream& output_stream,
      const u3shell::TwoBodyUnitTensorCoefficientsU3SPN& biquad_coefficients_pn
    )
  {
      for (const auto& key_value : biquad_coefficients_pn)
        {

          // extract unit tensor labels and coefficients
          const u3shell::TwoBodyUnitTensorLabelsU3S& biquad_labels= key_value.first;
          const CoefficientsPN& coefficients_pn = key_value.second;
        
          // extract unit tensor label groups
          u3shell::OperatorLabelsU3S::KeyType operator_labels_key;
          int rho0;
          u3shell::TwoBodyStateLabelsU3S::KeyType bra_key, ket_key;
          std::tie(operator_labels_key,rho0,bra_key,ket_key) = biquad_labels.Key();

          // extract operator label groups
          //
          // Note: pn-scheme labels still contain dummy isosopin variable, to be ignored.
          int N0;
          u3::SU3 x0;
          HalfInt S0;
          int g0;
          std::tie(N0,x0,S0,g0) = operator_labels_key;
        
          // extract bra labels
          int eta1p, eta2p;
          u3::SU3 xp;
          HalfInt Sp;
          std::tie(eta1p,eta2p,xp,Sp) = bra_key;

          // extract ket labels
          int eta1, eta2;
          u3::SU3 x;
          HalfInt S;
          std::tie(eta1,eta2,x,S) = ket_key;

          // label line
          output_stream
            << fmt::format(
                "{:d} {:d} {:d} {:d}   "
                "{:d} {:d} {:d} {:d}   "
                "{:d} {:d} {:d} {:d}   "
                "{:d} {:d} {:d} {:d}   ",
                eta1p,eta2p,eta2,eta1,
                1,xp.lambda(),xp.mu(),TwiceValue(Sp),
                1,x.lambda(),x.mu(),TwiceValue(S),
                rho0,x0.lambda(),x0.mu(),TwiceValue(S0)
              )
            << std::endl;

          // coefficient line
          output_stream
            << fmt::format(
                "{:e} {:e} {:e}",
                coefficients_pn.pppp,
                coefficients_pn.nnnn,
                coefficients_pn.pnnp
              )
            << std::endl;
        }
  }


 
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

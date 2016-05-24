/****************************************************************
  two_body_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/two_body_operator.h"

#include "cppformat/format.h"
#include "am/am.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/two_body_operator.h"

namespace u3shell {


  void TransformRelativeUnitTensorToTwoBodyUnitTensor
    (
     const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_unit_tensor_coefficients,
     u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_unit_tensor_coefficients
     )
  {
    // TODO
  }


  void TransformTwoBodyUnitTensorToBiquad(
      const u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_unit_tensor_coefficients,
      u3shell::TwoBodyUnitTensorCoefficientsU3ST& biquad_coefficients
    )
  {
    for (auto key_value : two_body_unit_tensor_coefficients)
      {

        // extract unit tensor labels and coefficient from map
        u3shell::TwoBodyUnitTensorLabelsU3ST two_body_unit_tensor_labels = key_value.first;
        double two_body_unit_tensor_coefficient = key_value.second;

        // extract individual labels on unit tensor

        // Note: This more compact ostensible approach does not work...
        // "error: invalid initialization of non-const reference..."
        //
        //   std::tie(
        //       std::tie(N0,x0,S0,T0,g0),
        //       rho0,
        //       std::tie(eta1p,eta2p,xp,Sp,Tp),
        //       std::tie(eta1,eta2,x,S,T)
        //     ) = two_body_unit_tensor_labels.Key();

        u3shell::OperatorLabelsU3ST::KeyType operator_labels_key;
        int rho0;
        u3shell::TwoBodyStateLabelsU3ST::KeyType bra_key, ket_key;
        std::tie(operator_labels_key,rho0,bra_key,ket_key) = two_body_unit_tensor_labels.Key();

        int N0;
        u3::SU3 x0;
        HalfInt S0,T0;
        int g0;
        std::tie(N0,x0,S0,T0,g0) = operator_labels_key;
        
        int eta1p, eta2p;
        u3::SU3 xp;
        HalfInt Sp, Tp;
        std::tie(eta1p,eta2p,xp,Sp,Tp) = bra_key;

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

        // iterate over multiplicities from recoupling cross projector
        int outer_multiplicity = u3::OuterMultiplicity(x,x0,xp);
        //std::cout
        //  << fmt::format("{} {} {} : {}",x0.Str(),x.Str(),xp.Str(),outer_multiplicity)
        //  << std::endl;

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
                  x0,S0,T0,
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
    // TODO
  }


 
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

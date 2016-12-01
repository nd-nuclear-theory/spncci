/****************************************************************
  two_body_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/two_body_operator.h"

#include <array>
#include <cmath>

#include "cppformat/format.h"
#include "am/am.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unu3.h" 

double zero_threshold=10e-13;

namespace u3shell {

  double TwoBodyNumberOperator(
    const u3shell::TwoBodyStateLabelsU3ST& bra,
    const u3shell::TwoBodyStateLabelsU3ST& ket
  )
  {
    double rme=0;
    if (bra==ket)
      rme=(ket.eta1()+ket.eta2());
    return rme;
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

            // std::cout<<fmt::format("{} {} {} {} {}",
            //   two_body_unit_tensor_coefficient,ParitySign(cross_projector_grade + biquad_grade),
            //   recoupling_phase, cross_projector_norm_factor,
            //   biquad_norm_factor)<<std::endl;

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
        int N0, g0;
        u3::SU3 x0;
        HalfInt S0,T0;
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
        // If like particles and particles in same shell, must check if U(3)xSU(2) symmetry
        // is allowed
        bool like_allowed=true;
        if(include_like&&(eta1==eta2))
          {
            un::SingleShellAllowedU3SIrreps single_shell_allowed_irreps_ket;
            un::GenerateAllowedSU3xSU2Irreps(eta1,2,single_shell_allowed_irreps_ket);
            if(not single_shell_allowed_irreps_ket.count(u3::U3S(u3::U3(2*eta1+3,x),S)))
              like_allowed=false;
          }
        if(include_like&&(eta1p==eta2p))
          {
            un::SingleShellAllowedU3SIrreps single_shell_allowed_irreps_bra;
            un::GenerateAllowedSU3xSU2Irreps(eta1p,2,single_shell_allowed_irreps_bra);
            if(not single_shell_allowed_irreps_bra.count(u3::U3S(u3::U3(2*eta1p+3,xp),Sp)))
              like_allowed=false;
          }

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
        //   signs_unlike (std::array<int,4>) : signs (+1,-1) for unlike-particle biquads <pn|pn>
        //     as {<12|12>, <12|21>, <21|12>, <21|12>}, where recall we label biquads using 
        //     the label ordering as it appears in the corresponding RME
        //
        // unlike particle terms active except in pure isovector (Tp,T,T0)=(1,1,1) special case

        bool include_unlike = !((Tp==1)&&(T==1)&&(T0==1));
        double norm_unlike;
        std::array<int,4> signs_unlike;

        if ((Tp==0)&&(T==0)&&(T0==0))
          {
            norm_unlike = 1/2.;
            signs_unlike = {-1,+1,+1,-1};
          }
        else if ((Tp==1)&&(T==0)&&(T0==1))
          {
            norm_unlike = 1/2.;
            signs_unlike = {-1,+1,-1,+1};
          }
        else if ((Tp==0)&&(T==1)&&(T0==1))
          {
            norm_unlike = 1/2.;
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
        if (include_like&&like_allowed)
          {
            biquad_labels_pn = u3shell::TwoBodyUnitTensorLabelsU3S(
                operator_labels_u3s,rho0,statep12,state12
              );
            coefficients_pn = u3shell::CoefficientsPN(
                biquad_coefficient*norm_like*signs_like[0],
                biquad_coefficient*norm_like*signs_like[1],
                0
              );
            // std::cout
            //   << fmt::format("  like {} {:+e} {:+e}",biquad_labels_pn.Str(),coefficients_pn.pppp,coefficients_pn.nnnn)
            //   << std::endl;
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
                0,0,biquad_coefficient*norm_unlike*signs_unlike[0]
              );
            // std::cout
            //   << fmt::format("  unlike 1212 {} {:+e}",biquad_labels_pn.Str(),coefficients_pn.pnnp)
            //   << std::endl;    
            biquad_coefficients_pn[biquad_labels_pn] += coefficients_pn;

            // term <12|21>
            biquad_labels_pn = u3shell::TwoBodyUnitTensorLabelsU3S(
                operator_labels_u3s,rho0,statep12,state21
              );
            coefficients_pn = u3shell::CoefficientsPN(
                0,0,biquad_coefficient*norm_unlike*signs_unlike[1]*ParitySign(grade)
              );
            // std::cout
            //   << fmt::format("  unlike 1221 {} {:+e}",biquad_labels_pn.Str(),coefficients_pn.pnnp)
            //   << std::endl;    
            biquad_coefficients_pn[biquad_labels_pn] += coefficients_pn;


            // term <21|12>
            biquad_labels_pn = u3shell::TwoBodyUnitTensorLabelsU3S(
                operator_labels_u3s,rho0,statep21,state12
              );
            coefficients_pn = u3shell::CoefficientsPN(
                0,0,biquad_coefficient*norm_unlike*signs_unlike[2]*ParitySign(gradep)
              );
            // std::cout
            //   << fmt::format("  unlike 2112 {} {:+e}",biquad_labels_pn.Str(),coefficients_pn.pnnp)
            //   << std::endl;    
            biquad_coefficients_pn[biquad_labels_pn] += coefficients_pn;

            // term <21|21>
            biquad_labels_pn = u3shell::TwoBodyUnitTensorLabelsU3S(
                operator_labels_u3s,rho0,statep21,state21
              );
            coefficients_pn = u3shell::CoefficientsPN(
                0,0,biquad_coefficient*norm_unlike*signs_unlike[3]*ParitySign(grade+gradep)
              );
            // std::cout
            //   << fmt::format("  unlike 2121 {} {:+e}",biquad_labels_pn.Str(),coefficients_pn.pnnp)
            //   << std::endl;    
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
          if(
              (fabs(coefficients_pn.pppp)>zero_threshold)
              ||(fabs(coefficients_pn.nnnn)>zero_threshold)
              ||(fabs(coefficients_pn.pnnp)>zero_threshold)
            )
          {
            // label line
            output_stream
              << fmt::format(
                  "{:d} {:d} {:d} {:d}   "
                  "{:d} {:d} {:d} {:d}   "
                  "{:d} {:d} {:d} {:d}   "
                  "{:d} {:d} {:d} {:d}   ",
                  eta1p,eta2p,eta2,eta1,
                  1,xp.lambda(),xp.mu(),TwiceValue(Sp),
                  1,x.mu(),x.lambda(),TwiceValue(S),
                  rho0,x0.lambda(),x0.mu(),TwiceValue(S0)
                )
              << std::endl;

            // coefficient line
            output_stream
              << fmt::format(
                  "{:+e} {:+e} {:+e}",
                  coefficients_pn.pppp,
                  coefficients_pn.nnnn,
                  coefficients_pn.pnnp
                )
              << std::endl;
          }
        }
  }


 
  ////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// two-body unit tensor enumeration
////////////////////////////////////////////////////////////////
void
GenerateTwoBodyUnitTensorLabelsU3ST(
      int Nmax, 
      // const u3shell::TwoBodySectorsU3ST two_body_sectors,
      std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST>& unit_tensor_labels_set
    )
  // Generate labels for U3ST-scheme two-body unit tensors acting
  // within the given two-body sectors.
  //
  // The resulting unit tensor labels are grouped by N0, i.e., the
  // number of oscillator quanta caried by the operator, where N0 will
  // vary from -2*Nmax to +2*Nmax.  This choice of grouping is mainly
  // dictated by consistency with the treatment of relative unit
  // tensor labels in spncci, though it also does provide a "neat"
  // sorting of the labels by N0 for iteration purposes.
  //
  // Arguments:
  //   Nmax (int) : maximum oscillator truncation
  //
  // Returns:
  //   (std::map<int,std::vector<TwoBodyUnitTensorLabelsU3ST>>)
  //   : map from N0 -> vector of relative unit tensor labels
{

  // define container object
    u3shell::TwoBodySpaceU3ST two_body_space(Nmax);
    u3shell::TwoBodySectorsU3ST two_body_sectors(two_body_space);
  // iterate over sectors
  for (int sector_index = 0; sector_index < two_body_sectors.size(); ++sector_index)
    {
      ////////////////////////////////
      // extract labels
      ////////////////////////////////

      // extract sector subspaces
      const u3shell::TwoBodySectorsU3ST::SectorType& sector = two_body_sectors.GetSector(sector_index);
      const u3shell::TwoBodySectorsU3ST::SubspaceType& bra_subspace = sector.bra_subspace();
      const u3shell::TwoBodySectorsU3ST::SubspaceType& ket_subspace = sector.ket_subspace();

      // extract bra subspace labels
      u3::U3 omegap;
      HalfInt Sp, Tp;
      int gp;
      std::tie(omegap,Sp,Tp,gp) =  bra_subspace.GetSubspaceLabels();

      // extract ket subspace labels
      u3::U3 omega;
      HalfInt S, T;
      int g;
      std::tie(omega,S,T,g) =  ket_subspace.GetSubspaceLabels();

      ////////////////////////////////
      // determine SU(3)xSxT couplings
      ////////////////////////////////

      // U(3): operator "destroys" x and "creates" xp

      // SU(3) couplings
      u3::SU3 xp = omegap.SU3();
      u3::SU3 x = omega.SU3();
      MultiplicityTagged<u3::SU3>::vector x0_set=u3::KroneckerProduct(Conjugate(x),xp);

      // U(1) coupling
      int etap = int(omegap.N());
      int eta = int(omega.N());
      int N0 = etap - eta;

      // SxT couplings
      HalfInt::pair S0_range = am::ProductAngularMomentumRange(S,Sp);
      HalfInt::pair T0_range = am::ProductAngularMomentumRange(T,Tp);

      // parity coupling
      int g0 = (g+gp)%2;
      
      ////////////////////////////////
      // construct label set
      ////////////////////////////////

      for (MultiplicityTagged<u3::SU3> x0_rho0max : x0_set)
        // for each SU(3)
        {
          // extract SU(3) labels
          u3::SU3 x0 = x0_rho0max.irrep;
          int rho0max = x0_rho0max.tag;
          int k=std::min(x0.lambda(),x0.mu());
          int l=std::max(x0.lambda(),x0.mu());
          // TODO : REMOVE WHEN J0 Requirement in SU3RME fixed. 
          MultiplicityTagged::vector L0_vector=u3::BranchingSO3Constrained(x0, S0_range);
          if (L0_vector.size()==0)
            continue;
          for (HalfInt S0=S0_range.first; S0 <= S0_range.second; ++S0)
          {
          // TODO : REMOVE WHEN J0 Requirement in SU3RME fixed. 
            MultiplicityTagged<int>::vector L0_vector=u3::BranchingSO3Constrained(x0, HalfInt::pair(S0,S0));
            if (L0_vector.size()==0)
              continue;

            for (HalfInt T0=T0_range.first; T0 <= T0_range.second; ++T0)
              // for each SxT
              {
                // collect tensorial labels
                u3shell::OperatorLabelsU3ST operator_labels(N0,x0,S0,T0,g0);

                // iterate over remaining (nontensorial) labels
                for (int bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
                  {
                    int eta1p, eta2p;
                    std::tie(eta1p, eta2p)= bra_subspace.GetStateLabels(bra_index);

                    for (int ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
                      // for each <bra|ket> pair
                      {             
                        // collect state labels
                        int eta1, eta2;
                        std::tie(eta1, eta2)= ket_subspace.GetStateLabels(ket_index);

                        u3shell::TwoBodyStateLabelsU3ST bra_labels(eta1p,eta2p,xp,Sp,Tp);
                        u3shell::TwoBodyStateLabelsU3ST ket_labels(eta1,eta2,x,S,T);

                        for (int rho0=1; rho0<=rho0max; ++rho0)
                          // for each multiplicity index rho0
                          { 
                            // collect unit tensor labels
                            u3shell::TwoBodyUnitTensorLabelsU3ST unit_tensor_labels(
                                operator_labels,
                                rho0,
                                bra_labels,ket_labels
                              );
                            // push unit tensor labels onto appropriate N0 set
                            unit_tensor_labels_set.push_back(unit_tensor_labels);
                          }
                      }
                  }
            }
          }
        }

    }
}


} // namespace

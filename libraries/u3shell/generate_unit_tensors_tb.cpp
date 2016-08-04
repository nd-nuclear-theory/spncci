/****************************************************************
  generate_unit_tensors_tb.cpp

  Construct files in Tomas's recoupler format for complete set of
  two-body unit tensors, in two-body truncated space.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/25/16 (mac): Created.

****************************************************************/

#include <algorithm>
#include <iostream>
#include <fstream>

#include "cppformat/format.h"

#include "am/am.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/u3st_scheme.h"
#include "u3shell/two_body_operator.h"

////////////////////////////////////////////////////////////////
// two-body unit tensor enumeration
////////////////////////////////////////////////////////////////

std::map<int,std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST>>
  GenerateTwoBodyUnitTensorLabelsU3ST(
      const u3shell::TwoBodySectorsU3ST two_body_sectors
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
  std::map<int,std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST>> unit_tensor_labels_set_by_N0;

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
      if(T!=1)    // REMOVE
        continue; // REMOVE
      if(Tp!=1)   // REMOVE
        continue; // REMOVE
      ////////////////////////////////
      // determine SU(3)xSxT couplings
      ////////////////////////////////

      // U(3): operator "destroys" x and "creates" xp

      // SU(3) couplings
      u3::SU3 xp = omegap.SU3();
      u3::SU3 x = omega.SU3();
      MultiplicityTagged<u3::SU3>::vector x0_set
        =u3::KroneckerProduct(Conjugate(x),xp);

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
          for (HalfInt S0=S0_range.first; S0 <= S0_range.second; ++S0)
          {
           /////////////REMOVE J0=0 restriction
            if(k==0)
              {
                if((l==0)&&(S0!=0))
                  continue;
                if((l%2==0)&&(int(S0)%2!=0))
                  continue;
                if((l%2==1)&&(S0!=1))
                  continue;
              }
            if(k>2)
              continue;
            if((k%2==1)&&(S0==0))
              continue;
            /////////////////////////
            // for (HalfInt T0=T0_range.first; T0 <= T0_range.second; ++T0)
              // for each SxT
              {
                HalfInt T0=0; //REMOVE
                // collect tensorial labels
                u3shell::OperatorLabelsU3ST operator_labels(N0,x0,S0,T0,g0);

                // iterate over remaining (nontensorial) labels
                for (int bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
                  {
                    int eta1p, eta2p;
                    std::tie(eta1p, eta2p)= bra_subspace.GetStateLabels(bra_index);
                    if ((eta1p==eta2p))//&&(int(u3::ConjugationGrade(xp)+Sp+Tp)%2==0)) REMOVE
                      continue;
                    for (int ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
                      // for each <bra|ket> pair
                      {
                        // collect state labels
                        int eta1, eta2;
                        std::tie(eta1, eta2)= ket_subspace.GetStateLabels(ket_index);
                        if ((eta1==eta2))//&&(int(u3::ConjugationGrade(x)+S+T)%2==0))REMOVE
                          continue;

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
                            unit_tensor_labels_set_by_N0[N0].push_back(unit_tensor_labels);
                          }
                      }
                  }
            }
          }
        }

    }

  return unit_tensor_labels_set_by_N0;
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // configuration 
  const int max_num_output_files = 10;  // safety limit for number of outputted files
  bool verbose = true;  // debugging verbosity
  bool write_output = true;  // generate actual output files for recoupler


  // initialize su3lib
  u3::U3CoefInit();

  // generate list of unit tensor labels
  const int Nmax = 2;
  u3shell::TwoBodySpaceU3ST two_body_space(Nmax);
  u3shell::TwoBodySectorsU3ST two_body_sectors(two_body_space);
  std::map<int,std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST>> unit_tensor_labels_set_by_N0
    = GenerateTwoBodyUnitTensorLabelsU3ST(two_body_sectors);

  // flatten list of unit tensor labels
  std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST> unit_tensor_labels_set;
  for (const auto& key_value : unit_tensor_labels_set_by_N0)
    {
      int N0 = key_value.first;
      const std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST>& unit_tensor_labels_set_partial = key_value.second;
      //std::insert(
      unit_tensor_labels_set.insert(
          unit_tensor_labels_set.end(),
          unit_tensor_labels_set_partial.begin(),
          unit_tensor_labels_set_partial.end()
        );
    }

  // write list of unit tensor labels
  std::string label_stream_filename = "generate_unit_tensors_tb_labels.dat";
  std::ofstream label_stream(label_stream_filename);
  for (int unit_tensor_index=0; unit_tensor_index < unit_tensor_labels_set.size(); ++unit_tensor_index)
    {
      // extract unit tensor labels
      const u3shell::TwoBodyUnitTensorLabelsU3ST& unit_tensor_labels
        = unit_tensor_labels_set[unit_tensor_index];

      // Note: If want to ouput labels in tabular form, the individual
      // labels will need to be extracted.  However, this might better
      // be delegated to a helper function.  In the meantime, we use
      // TwoBodyUnitTensorLabelsU3ST.Str() for human-readable
      // formatting.
      //
      // // extract unit tensor label groups
      // u3shell::OperatorLabelsU3ST::KeyType operator_labels_key;
      // int rho0;
      // u3shell::TwoBodyStateLabelsU3ST::KeyType bra_key, ket_key;
      // std::tie(operator_labels_key,rho0,bra_key,ket_key) = unit_tensor_labels.Key();
      // 
      // // extract operator label groups
      // int N0;
      // u3::SU3 x0;
      // HalfInt S0,T0;
      // int g0;
      // std::tie(N0,x0,S0,T0,g0) = operator_labels_key;
      //   
      // // extract bra labels
      // int eta1p, eta2p;
      // u3::SU3 xp;
      // HalfInt Sp, Tp;
      // std::tie(eta1p,eta2p,xp,Sp,Tp) = bra_key;
      // 
      // // extract ket labels
      // int eta1, eta2;
      // u3::SU3 x;
      // HalfInt S, T;
      // std::tie(eta1,eta2,x,S,T) = ket_key;

      label_stream
        << fmt::format(
            "{:06d}  {}",
            unit_tensor_index,
            unit_tensor_labels.Str()
          )
        << std::endl; 
    }
  label_stream.close();


  // output unit tensors for recoupler

  for (int unit_tensor_index=0; unit_tensor_index < unit_tensor_labels_set.size(); ++unit_tensor_index)
    {

      // safety short circuit
      if (unit_tensor_index >= max_num_output_files)
        break;

      // head iteration (verbose output)
      if (verbose)
        {
          std::cout << std::endl;
          std::cout << fmt::format("operator {:06d}",unit_tensor_index) << std::endl;
          std::cout << std::endl;
        }

      // extract unit tensor labels
      const u3shell::TwoBodyUnitTensorLabelsU3ST& two_body_unit_tensor_labels
        = unit_tensor_labels_set[unit_tensor_index];

      // debugging
      //if (!((two_body_unit_tensor_labels.ket().T()==1)&& (two_body_unit_tensor_labels.bra().T()==1)))
      //  continue;

      // declare coefficient containers
      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
      u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
      u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;

      // populate operator
      two_body_unit_tensor_coefficients[two_body_unit_tensor_labels] = 1;

      // dump operator (verbose output)
      if (verbose)
        {
          // dump operator
          std::cout << "two_body_unit_tensor_coefficients" << std::endl;
          for (auto key_value : two_body_unit_tensor_coefficients)
            {

              // extract unit tensor labels and coefficients
              auto labels= key_value.first;
              double coefficient = key_value.second;
        
              std::cout 
                << fmt::format(
                    "  {} {:e}",
                    labels.Str(),
                    coefficient
                  )
                << std::endl;
            }
        }

      // convert to biquads
      u3shell::TransformTwoBodyUnitTensorToBiquad(
          two_body_unit_tensor_coefficients,
          biquad_coefficients
        );

      // dump operator (verbose output)
      if (verbose)
        {
          std::cout << "biquad_coefficients" << std::endl;
          for (auto key_value : biquad_coefficients)
            {

              // extract unit tensor labels and coefficients
              auto labels= key_value.first;
              double coefficient = key_value.second;
        
              std::cout 
                << fmt::format(
                    "  {} {:e}",
                    labels.Str(),
                    coefficient
                  )
                << std::endl;
            }
        }

      // convert to pn scheme
      // std::cout << "biquad_coefficients_pn conversion..." << std::endl;
      u3shell::TransformBiquadToPNScheme(
          biquad_coefficients,
          biquad_coefficients_pn
        );

      // dump operator (verbose output)
      if (verbose)
        {
          std::cout << "biquad_coefficients_pn recoupler output" << std::endl;
          WriteTwoBodyOperatorRecoupler(std::cout,biquad_coefficients_pn);
        }

      // // dump operator (verbose output)
      // if (verbose)
      //   {
      //     std::cout << "biquad_coefficients_pn" << std::endl;
      //     for (auto key_value : biquad_coefficients_pn)
      //       {

      //         // extract unit tensor labels and coefficients
      //         auto labels= key_value.first;
      //         u3shell::CoefficientsPN& coefficients_pn = key_value.second;
        
      //         std::string flag;
      //         if (abs(coefficients_pn.pnnp)>1e-5)
      //             flag = "****";
      //         std::cout 
      //           << fmt::format(
      //               "  {} {:+20e} {:+20e} {:+20e} {}",
      //               labels.Str(),
      //               coefficients_pn.pppp,
      //               coefficients_pn.nnnn,
      //               coefficients_pn.pnnp,
      //               flag
      //             )
      //           << std::endl;
      //       }
      //   }

      // dump operator in recoupler format
      if (write_output)
        {
          std::cout<<fmt::format("{} operator{:06d}",two_body_unit_tensor_labels.Str(),unit_tensor_index)<<std::endl;
          std::string operator_stream_filename = fmt::format("operator{:06d}.recoupler",unit_tensor_index);
          std::ofstream operator_stream(operator_stream_filename);
          WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
          operator_stream.close();
        }

    }

  // termination
  return 0;
}

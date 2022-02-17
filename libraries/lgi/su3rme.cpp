/****************************************************************
  su3rme.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "lgi/su3rme.h"


#include "fmt/format.h"
// #include "lgi/lgi.h"
#include "lsu3shell/lsu3shell_basis.h"

// #include "sp3rlib/sp3r.h"
// #include "spncci/spncci_basis.h"
#include "u3shell/u3spn_scheme.h"
#include "utilities/utilities.h"

#include "SU3ME/proton_neutron_ncsmSU3Basis.h"
#include "LookUpContainers/CWig9lmLookUpTable.h"
#include "LookUpContainers/lock.h"
#include "SU3ME/BaseSU3Irreps.h"
#include "SU3ME/ComputeOperatorMatrix.h"
#include "SU3ME/MeEvaluationHelpers.h"
#include "su3.h"

namespace lsu3shell
{
  unsigned int get_num_U3PNSPN_irreps(const lsu3::CncsmSU3xSU2Basis& basis)
    {
      const uint32_t number_ipin_blocks = basis.NumberOfBlocks();

      unsigned int num_irreps=0; 
      for (unsigned int ipin_block = 0; ipin_block < number_ipin_blocks; ipin_block++)
        {
          uint32_t ip = basis.getProtonIrrepId(ipin_block);
          uint32_t in = basis.getNeutronIrrepId(ipin_block);
          uint16_t alpha_p_max = basis.getMult_p(ip);
          uint16_t alpha_n_max = basis.getMult_n(in);

          int32_t ibegin = basis.blockBegin(ipin_block);
          int32_t iend = basis.blockEnd(ipin_block);
          
          // std::cout<<"alphas " <<alpha_n_max<<"  "<<alpha_p_max<<std::endl;
          // std::cout<<"begin end "<<ibegin<<" "<<iend<<std::endl;
          for (int32_t iwpn = ibegin; iwpn < iend; ++iwpn) 
            {
              SU3xSU2::LABELS omega_pn(basis.getOmega_pn(ip, in, iwpn));
              unsigned int dim=alpha_n_max*alpha_p_max*omega_pn.rho;
              num_irreps+=dim;
            }
        }
      return num_irreps;
    }

/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
template <typename tDataType>
void WriteBinary(std::ostream& os, tDataType data)
// Write binary data item to stream.
//
// Note that, if the template parameter is omitted, the data type
// of the value given for data will determine the output type, but
// explicitly giving the template parameter casts the data to the
// given data type tDataType.  Explicitly giving the template
// parameter is recommended both to document the data format in
// the output file and to avoid any ambiguity of the output type.
//
// Arguments:
//   os (input): binary stream for output
//   data (input): data value to output
//
// Ex:
//   mcutils::WriteBinary<float>(out_stream,value);
{
   os.write(reinterpret_cast<const char*>(&data), sizeof(data));
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

std::vector<double> ComputeRME(std::vector<MECalculatorData>& rmeCoeffsPNPN)
// Compute RMEs for a single pair of bra and ket irrep subspaces (i.e.,
// subspaces of equivalent irreps, distinguished by upstream quantum numbers),
// for a single operator (e.g., unit tensor).
//
// rmeCoeffsPNPN: expansion of operator in a+a+aa terms
{
   // These values must be constant for each tensor
   int32_t bra_max = rmeCoeffsPNPN[0].rmes_->m_bra_max;
   int32_t ket_max = rmeCoeffsPNPN[0].rmes_->m_ket_max;
   int32_t rhot_max = rmeCoeffsPNPN[0].rmes_->m_rhot_max;
   // tensor multiplicity is however tensor dependent
   int32_t rho0_max;

   // resulting rmes do not have rho0 dependence
   int32_t nrmes_to_compute = bra_max * ket_max * rhot_max;
   // std::cout<<"num rmes to comput "<<nrmes_to_compute<<std::endl;
   std::vector<double> result_rmes(nrmes_to_compute, 0.0);

   // iterate over PP NN and PN tensors each has the same lm0 mu0 and S0
   // quantum labels but it can have a different rho0max multiplicity
   for (size_t itensor = 0; itensor < rmeCoeffsPNPN.size(); ++itensor) {
      assert(bra_max == rmeCoeffsPNPN[itensor].rmes_->m_bra_max);
      assert(ket_max == rmeCoeffsPNPN[itensor].rmes_->m_ket_max);
      assert(rhot_max == rmeCoeffsPNPN[itensor].rmes_->m_rhot_max);

      rho0_max = rmeCoeffsPNPN[itensor].rmes_->m_tensor_max;
      int32_t rme_index = 0;
      for (int i = 0; i < bra_max; ++i) {
         for (int j = 0; j < ket_max; ++j) {
            //Summing over irho0.  For all relative unit tensors and Sp(3,R) generators rho0_max=1
            for (int irho0 = 0; irho0 < rho0_max; ++irho0) {
               // vector with rhot_max elements
               // rmes = < i ||| rho0 ||| j>[i][j][rho0] = {...}
               RME::DOUBLE* rmes = rmeCoeffsPNPN[itensor].rmes_->GetVector(i, j, irho0);
               TENSOR_STRENGTH coeff_rho0 = rmeCoeffsPNPN[itensor].coeffs_[irho0];
               for (int irhot = 0; irhot < rhot_max; ++irhot) {
                  result_rmes[rme_index + irhot] += (double)coeff_rho0 * (double)rmes[irhot];
               }
            }
            rme_index += rhot_max;
         }
      }
   }
   return result_rmes;
}


/////////////////////////////////////////////////////////////////////////////

basis::OperatorBlocks<double> CalculateRME( 
  const CInteractionPPNN& interactionPPNN,
  const CInteractionPN& interactionPN,
  const lsu3::CncsmSU3xSU2Basis& bra,
  const lsu3::CncsmSU3xSU2Basis& ket,
  int dN0,
  const SU3xSU2::LABELS& w0,
  int rhot_max
  )
  // Iterate over bra and ket irrep subspaces (i.e., subspaces of equivalent
  // irreps, distinguished by upstream quantum numbers), and for each pair of
  // irrep subspaces connected by given operator (e.g., unit tensor), calculate
  // RMEs of the given operator between that pair, and write these to file.
  //
  // Calls: ComputeRME
{
  // Should probably return matrix or vector of matrices...
  int rows = lsu3shell::get_num_U3PNSPN_irreps(bra);
  int cols = lsu3shell::get_num_U3PNSPN_irreps(ket);
  // std::cout<<"rows "<<rows<<"  cols "<<cols<<std::endl;
  basis::OperatorBlocks<double> operator_blocks(rhot_max);
  for(basis::OperatorBlock<double>& block : operator_blocks)
    block = Eigen::MatrixXd::Zero(rows,cols);

  //////////////////////////////////////////////////////////////////
  // std::cout<<"checkpoint 1"<<std::endl;
  std::vector<unsigned char> hoShells_n, hoShells_p;
  std::vector<CTensorGroup*> tensorGroupsPP, tensorGroupsNN;
  std::vector<CTensorGroup_ada*> tensorGroups_p_pn, tensorGroups_n_pn;
  std::vector<int> phasePP, phaseNN, phase_p_pn, phase_n_pn;

  unsigned char num_vacuums_J_distr_p;
  unsigned char num_vacuums_J_distr_n;
  std::vector<std::pair<CRMECalculator*, CTensorGroup::COEFF_DOUBLE*> > selected_tensorsPP,selected_tensorsNN;
  std::vector<std::pair<CRMECalculator*, unsigned int> > selected_tensors_p_pn,selected_tensors_n_pn;
  std::vector<RmeCoeffsSU3SO3CGTablePointers> rmeCoeffsPP, rmeCoeffsNN;
  std::vector<std::pair<SU3xSU2::RME*, unsigned int> > rme_index_p_pn, rme_index_n_pn;
  std::vector<MECalculatorData> rmeCoeffsPNPN;
  SU3xSU2::RME identityOperatorRMEPP, identityOperatorRMENN;

  uint16_t max_mult_p = bra.getMaximalMultiplicity_p();
  uint16_t max_mult_n = bra.getMaximalMultiplicity_n();
  InitializeIdenticalOperatorRME(identityOperatorRMEPP, max_mult_p*max_mult_p);
  InitializeIdenticalOperatorRME(identityOperatorRMENN, max_mult_n*max_mult_n);
  // std::cout<<"checkpoint 2"<<std::endl;
  SingleDistributionSmallVector distr_ip, distr_in, distr_jp, distr_jn;
  UN::SU3xSU2_VEC gamma_ip, gamma_in, gamma_jp, gamma_jn;
  SU3xSU2_SMALL_VEC vW_ip, vW_in, vW_jp, vW_jn;

  unsigned char deltaP, deltaN;
  const uint32_t number_ipin_blocks = bra.NumberOfBlocks();
  const uint32_t number_jpjn_blocks = ket.NumberOfBlocks();

  int32_t icurrentDistr_p, icurrentDistr_n;
  int32_t icurrentGamma_p, icurrentGamma_n;

  // main block-by-block calculation
   int32_t row_irrep = 0;
   // std::cout<<"for blocks"<<std::endl;
   for (unsigned int ipin_block = 0; ipin_block < number_ipin_blocks; ipin_block++) 
   {
      uint32_t ip = bra.getProtonIrrepId(ipin_block);
      uint32_t in = bra.getNeutronIrrepId(ipin_block);
      uint16_t aip_max = bra.getMult_p(ip);
      uint16_t ain_max = bra.getMult_n(in);

      int32_t Nhw_bra = bra.nhw_p(ip) + bra.nhw_n(in);
      uint16_t nrow_irreps = bra.NumberPNIrrepsInBlock(ipin_block);

      SU3xSU2::LABELS w_ip(bra.getProtonSU3xSU2(ip));
      SU3xSU2::LABELS w_in(bra.getNeutronSU3xSU2(in));

      uint16_t ilastDistr_p(std::numeric_limits<uint16_t>::max());
      uint16_t ilastDistr_n(std::numeric_limits<uint16_t>::max());

      uint32_t ilastGamma_p(std::numeric_limits<uint32_t>::max());
      uint32_t ilastGamma_n(std::numeric_limits<uint32_t>::max());

      uint32_t last_jp(std::numeric_limits<uint32_t>::max());
      uint32_t last_jn(std::numeric_limits<uint32_t>::max());

      int32_t col_irrep = 0;
      // std::cout<<"J blocks"<<std::endl;
      for (unsigned int jpjn_block = 0; jpjn_block < number_jpjn_blocks; jpjn_block++) 
      {
        uint32_t jp = ket.getProtonIrrepId(jpjn_block);
        uint32_t jn = ket.getNeutronIrrepId(jpjn_block);
        int32_t Nhw_ket = ket.nhw_p(jp) + ket.nhw_n(jn);
        uint16_t ncol_irreps = ket.NumberPNIrrepsInBlock(jpjn_block);
        if ((Nhw_ket + dN0) != Nhw_bra) 
          {
            col_irrep = col_irrep + ncol_irreps;
            continue;
          }

         SU3xSU2::LABELS w_jp(ket.getProtonSU3xSU2(jp));
         SU3xSU2::LABELS w_jn(ket.getNeutronSU3xSU2(jn));
         uint16_t ajp_max = ket.getMult_p(jp);
         uint16_t ajn_max = ket.getMult_n(jn);

         // std::cout<<"jp"<<std::endl;
         if (jp != last_jp) 
         {
            icurrentDistr_p = ket.getIndex_p<lsu3::CncsmSU3xSU2Basis::kDistr>(jp);
            icurrentGamma_p = ket.getIndex_p<lsu3::CncsmSU3xSU2Basis::kGamma>(jp);
            
            if (ilastDistr_p != icurrentDistr_p) 
              {
                distr_ip.resize(0);
                bra.getDistr_p(ip, distr_ip);
                gamma_ip.resize(0);
                bra.getGamma_p(ip, gamma_ip);
                vW_ip.resize(0);
                bra.getOmega_p(ip, vW_ip);
                distr_jp.resize(0);
                ket.getDistr_p(jp, distr_jp);
                hoShells_p.resize(0);
                deltaP = TransformDistributions_SelectByDistribution(
                     interactionPPNN, interactionPN, distr_ip, gamma_ip, vW_ip, distr_jp, hoShells_p,
                     num_vacuums_J_distr_p, phasePP, tensorGroupsPP, phase_p_pn, tensorGroups_p_pn);
              }

            if (ilastGamma_p != icurrentGamma_p || ilastDistr_p != icurrentDistr_p) 
              {
                 if (deltaP <= 4) {
                    if (!selected_tensorsPP.empty()) 
                      {
                         std::for_each(selected_tensorsPP.begin(), selected_tensorsPP.end(),
                                       CTensorGroup::DeleteCRMECalculatorPtrs());
                      }
                    selected_tensorsPP.resize(0);

                    if (!selected_tensors_p_pn.empty()) 
                      {
                        std::for_each(selected_tensors_p_pn.begin(), selected_tensors_p_pn.end(),
                                       CTensorGroup_ada::DeleteCRMECalculatorPtrs());
                      }

                    selected_tensors_p_pn.resize(0);
                    gamma_jp.resize(0);
                    ket.getGamma_p(jp, gamma_jp);
                    TransformGammaKet_SelectByGammas(
                        hoShells_p, distr_jp, num_vacuums_J_distr_p, nucleon::PROTON, phasePP,
                        tensorGroupsPP, phase_p_pn, tensorGroups_p_pn, gamma_ip, gamma_jp,
                        selected_tensorsPP, selected_tensors_p_pn);
                 }
              }

            if (deltaP <= 4) 
              {
                Reset_rmeCoeffs(rmeCoeffsPP);
                Reset_rmeIndex(rme_index_p_pn);

                vW_jp.resize(0);
                ket.getOmega_p(jp, vW_jp);
                TransformOmegaKet_CalculateRME(
                  distr_jp, gamma_ip, vW_ip, gamma_jp, num_vacuums_J_distr_p, selected_tensorsPP,
                  selected_tensors_p_pn, vW_jp, rmeCoeffsPP, rme_index_p_pn);
              }

            ilastDistr_p = icurrentDistr_p;
            ilastGamma_p = icurrentGamma_p;
         }
         // std::cout<<"jn "<<std::endl;
         if (jn != last_jn) 
           {
              icurrentDistr_n = ket.getIndex_n<lsu3::CncsmSU3xSU2Basis::kDistr>(jn);
              icurrentGamma_n = ket.getIndex_n<lsu3::CncsmSU3xSU2Basis::kGamma>(jn);

              if (ilastDistr_n != icurrentDistr_n) 
              {
                 distr_in.resize(0);
                 bra.getDistr_n(in, distr_in);
                 gamma_in.resize(0);
                 bra.getGamma_n(in, gamma_in);
                 vW_in.resize(0);
                 bra.getOmega_n(in, vW_in);

                 distr_jn.resize(0);
                 ket.getDistr_n(jn, distr_jn);
                 hoShells_n.resize(0);
                 deltaN = TransformDistributions_SelectByDistribution(
                     interactionPPNN, interactionPN, distr_in, gamma_in, vW_in, distr_jn, hoShells_n,
                     num_vacuums_J_distr_n, phaseNN, tensorGroupsNN, phase_n_pn, tensorGroups_n_pn);
              }

              if (ilastGamma_n != icurrentGamma_n || ilastDistr_n != icurrentDistr_n) 
              {
                 if (deltaN <= 4) 
                 {
                    if (!selected_tensorsNN.empty()) 
                    {
                       std::for_each(selected_tensorsNN.begin(), selected_tensorsNN.end(),
                                     CTensorGroup::DeleteCRMECalculatorPtrs());
                    }
                    selected_tensorsNN.resize(0);

                    if (!selected_tensors_n_pn.empty()) 
                    {
                       std::for_each(selected_tensors_n_pn.begin(), selected_tensors_n_pn.end(),
                                     CTensorGroup_ada::DeleteCRMECalculatorPtrs());
                    }
                    selected_tensors_n_pn.resize(0);

                    gamma_jn.resize(0);
                    ket.getGamma_n(jn, gamma_jn);
                    TransformGammaKet_SelectByGammas(
                        hoShells_n, distr_jn, num_vacuums_J_distr_n, nucleon::NEUTRON, phaseNN,
                        tensorGroupsNN, phase_n_pn, tensorGroups_n_pn, gamma_in, gamma_jn,
                        selected_tensorsNN, selected_tensors_n_pn);
                 }
              }

            if (deltaN <= 4) 
            {
               Reset_rmeCoeffs(rmeCoeffsNN);
               Reset_rmeIndex(rme_index_n_pn);

               vW_jn.resize(0);
               ket.getOmega_n(jn, vW_jn);
               TransformOmegaKet_CalculateRME(
                   distr_jn, gamma_in, vW_in, gamma_jn, num_vacuums_J_distr_n, selected_tensorsNN,
                   selected_tensors_n_pn, vW_jn, rmeCoeffsNN, rme_index_n_pn);
            }

            ilastDistr_n = icurrentDistr_n;
            ilastGamma_n = icurrentGamma_n;
         }
         //////////////////////////////////////////////////////////////////////////////////////////
         // Iterating over omega S
         //////////////////////////////////////////////////////////////////////////////////////////
         // loop over wpn that result from coupling ip x in

         // std::cout<<"loop over wpn"<<std::endl;
         int32_t ibegin = bra.blockBegin(ipin_block);
         int32_t iend = bra.blockEnd(ipin_block);
         
         for (int32_t iwpn = ibegin, i = 0; iwpn < iend; ++iwpn, ++i) 
         {
            SU3xSU2::LABELS omega_pn_I(bra.getOmega_pn(ip, in, iwpn));
            IRREPBASIS braSU3xSU2basis(bra.Get_Omega_pn_Basis(iwpn));

            int32_t jbegin = ket.blockBegin(jpjn_block);
            int32_t jend = ket.blockEnd(jpjn_block);
            for (int32_t jwpn = jbegin, j = 0; jwpn < jend; ++jwpn, ++j) 
            {
              SU3xSU2::LABELS omega_pn_J(ket.getOmega_pn(jp, jn, jwpn));
              int rhot_max = SU3::mult(omega_pn_J, w0, omega_pn_I);

              if (!SU2::mult(omega_pn_J.S2, w0.S2, omega_pn_I.S2) || !rhot_max) 
               {
                  continue;
               }

              IRREPBASIS ketSU3xSU2basis(ket.Get_Omega_pn_Basis(jwpn));

              if (deltaP + deltaN <= 4) 
               {
                  Reset_rmeCoeffs(rmeCoeffsPNPN);
                  // std::cout<<"checkpoint 3"<<std::endl;
                  if (in == jn && !rmeCoeffsPP.empty()) 
                    {
                      // create structure with < an (lmn mun)Sn ||| 1 ||| an' (lmn mun) Sn> r.m.e.
                      CreateIdentityOperatorRME(w_in, w_jn, ajn_max, identityOperatorRMENN);
                      // std::cout<<"checkpoint 3.1"<<std::endl;
                      // calculate < {(lmp mup)Sp x (lmn mun)Sn} wf Sf ||| [I_{nn} x T_{pp}]
                      //|||{(lmp' mup')Sp' x (lmn' mun')Sn'} wi Si >_{rhot}
                      Calculate_Proton_x_Identity_MeData(
                        omega_pn_I, omega_pn_J, rmeCoeffsPP,
                        identityOperatorRMENN, rmeCoeffsPNPN
                      );
                      // std::cout<<"checkpoint 3.2"<<std::endl;
                    }
                  // std::cout<<"checkpoint 4"<<std::endl;
                  if (ip == jp && !rmeCoeffsNN.empty()) 
                    {
                       // create structure with < ap (lmp mup)Sp ||| 1 ||| ap' (lmp mup) Sp> r.m.e.
                       CreateIdentityOperatorRME(w_ip, w_jp, ajp_max, identityOperatorRMEPP);
                       // calculate < {(lmp mup)Sp x (lmn mun)Sn} wf Sf ||| [T_{nn} x I_{pp}]
                       //|||{(lmp' mup')Sp' x (lmn' mun')Sn'} wi Si >_{rhot}
                       Calculate_Identity_x_Neutron_MeData(
                           omega_pn_I, omega_pn_J, identityOperatorRMEPP, rmeCoeffsNN, rmeCoeffsPNPN);
                    }
                  // std::cout<<"checkpoint 5"<<std::endl;
                  if (!rme_index_p_pn.empty() && !rme_index_n_pn.empty()) 
                    {
                       // This function assumes coefficients depend on rho0 only.
                       // It also disregard SU(3) clebsch gordan coefficients as
                       // we do not use the Wigner-Eckhart theorem.
                      CalculatePNOperatorRMECoeffs(
                        interactionPN, omega_pn_I, omega_pn_J,
                        rme_index_p_pn, rme_index_n_pn, rmeCoeffsPNPN
                      );
                    }
                  // std::cout<<"checkpoint 6"<<std::endl;
                  if (!rmeCoeffsPNPN.empty()) 
                    {
                      std::vector<double> rmes = ComputeRME(rmeCoeffsPNPN);

                      int k=0;
                      // std::cout<<"num rmes "<<rmes.size()<<std::endl;
                      for(int i_row=0; i_row<aip_max*ain_max; ++i_row)
                        for(int i_col=0; i_col<ajp_max*ajn_max; ++i_col)
                          for(int i_rhot=0; i_rhot<rhot_max; ++i_rhot)
                            {
                              int row_index = row_irrep + i;
                              int col_index = col_irrep + j;
                              // std::cout<<rmes[k]<<std::endl;
                              operator_blocks[i_rhot](row_index,col_index)=rmes[k];
                              k++;
                            }
                      // std::cout<<"kk "<<k<<std::endl;
                      // std::cout<<"k: "<<k<<"  "<<rmes.size()<<"  "<<aip_max*ain_max<<"  "<<ajp_max*ajn_max
                      // <<"  "<<rhot_max<<std::endl;
                      assert(k==rmes.size());
                      // std::cout<<"num rmes accumulated "<<k<<std::endl;
                    }
               }
            } //end ket omega S
         } //end bra omega S
         col_irrep = col_irrep + ncol_irreps;
         last_jp = jp;
         last_jn = jn;
      }
      row_irrep = row_irrep + nrow_irreps;
      
   }
   Reset_rmeCoeffs(rmeCoeffsPP);
   Reset_rmeCoeffs(rmeCoeffsNN);
   Reset_rmeCoeffs(rmeCoeffsPNPN);
   // delete arrays allocated in identityOperatorRME?? structures
   delete[] identityOperatorRMEPP.m_rme;
   delete[] identityOperatorRMENN.m_rme;
 
  return operator_blocks;
}


}// end namespace

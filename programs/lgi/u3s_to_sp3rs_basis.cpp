#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "lgi/lgi.h"
#include "lsu3shell/lsu3shell_basis.h"
// #include "lsu3shell/su3rme.h"
#include "sp3rlib/sp3r.h"
#include "spncci/spncci_basis.h"
#include "u3shell/u3spn_scheme.h"
#include "utilities/utilities.h"

#include "SU3ME/proton_neutron_ncsmSU3Basis.h"
#include "LookUpContainers/CWig9lmLookUpTable.h"
#include "LookUpContainers/lock.h"
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "SU3ME/BaseSU3Irreps.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/ComputeOperatorMatrix.h"
#include "SU3ME/InteractionPPNN.h"


#include "su3.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <regex>
#include <stack>
#include <stdexcept>
#include <vector>



typedef std::tuple<int,int,int> SpinTuple;
typedef std::tuple<SpinTuple,SpinTuple,SpinTuple,SpinTuple> CompoundSpinTuple;

  struct RunParameters
  // Stores simple parameters for run
  {
    // filenames
    std::string input_filename;
    // mode
    nuclide::NuclideType nuclide;
    int Nmax;
  };

  RunParameters ProcessArguments(int argc, char **argv)
  {
    RunParameters run_parameters;
    // usage message
    if (argc-1 < 3)
      {
        std::cout << "Syntax: get_cmf_lgi_labels Z N Nmax optional:input_filename" << std::endl;
        std::exit(EXIT_SUCCESS);
      }

    // nuclide
    int Z, N, Nmax;
    {
      std::istringstream parameter_stream(argv[1]);
      parameter_stream >> Z;
      if (!parameter_stream)
      {
        std::cerr << "Expecting numeric value for Z" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    {
      std::istringstream parameter_stream(argv[2]);
      parameter_stream >> N;
      if (!parameter_stream)
      {
        std::cerr << "Expecting numeric value for Z" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    {
      std::istringstream parameter_stream(argv[3]);
      parameter_stream >> Nmax;
      if (!parameter_stream)
      {
        std::cerr << "Expecting numeric value for Nmax" << std::endl;
        std::exit(EXIT_FAILURE);
      }

    }
    
    run_parameters.nuclide = {Z,N};
    run_parameters.Nmax=Nmax;
    
    if(argc==5)
      // input filename
      run_parameters.input_filename = argv[4];
    else
      run_parameters.input_filename = "None";
    
    return run_parameters;
  }



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


void CalculateRME(
  const std::string& rme_file_name, 
  const CInteractionPPNN& interactionPPNN,
  const CInteractionPN& interactionPN,
  const lsu3::CncsmSU3xSU2Basis& bra,
  const lsu3::CncsmSU3xSU2Basis& ket,
  int dN0,
  const SU3xSU2::LABELS& w0
  )
// Iterate over bra and ket irrep subspaces (i.e., subspaces of equivalent
// irreps, distinguished by upstream quantum numbers), and for each pair of
// irrep subspaces connected by given operator (e.g., unit tensor), calculate
// RMEs of the given operator between that pair, and write these to file.
//
// Calls: ComputeRME
{
  //////////////////////////////////////////////////////////////////
  /// Previously global
  // note: unsigned short int overflows, leading to multiplicity errors
  typedef unsigned int RMEIndexType;  

  // statistics accumulators
  bool g_zero_threshold = 1e-8;  // for output diagnostics
  //////////////////////////////////////////////////////////////////

  // set up statistics for this tensor
  int num_groups_allowed = 0;
  int num_rmes_allowed = 0;
  int num_groups_nonzero = 0;
  int num_rmes_nonzero = 0;
  double mpi_tensor_start_time = MPI_Wtime();

  // open output file
  std::ofstream out_file;
  std::ios_base::openmode mode_argument = std::ios_base::out;
  mode_argument |= std::ios_base::binary;
  out_file.open(rme_file_name, mode_argument);


   std::vector<unsigned char> hoShells_n, hoShells_p;
   std::vector<CTensorGroup*> tensorGroupsPP, tensorGroupsNN;
   std::vector<CTensorGroup_ada*> tensorGroups_p_pn, tensorGroups_n_pn;

   std::vector<int> phasePP, phaseNN, phase_p_pn, phase_n_pn;

   unsigned char num_vacuums_J_distr_p;
   unsigned char num_vacuums_J_distr_n;
   std::vector<std::pair<CRMECalculator*, CTensorGroup::COEFF_DOUBLE*> > selected_tensorsPP,
       selected_tensorsNN;
   std::vector<std::pair<CRMECalculator*, unsigned int> > selected_tensors_p_pn,
       selected_tensors_n_pn;

   std::vector<RmeCoeffsSU3SO3CGTablePointers> rmeCoeffsPP, rmeCoeffsNN;
   std::vector<std::pair<SU3xSU2::RME*, unsigned int> > rme_index_p_pn, rme_index_n_pn;

   std::vector<MECalculatorData> rmeCoeffsPNPN;

   SU3xSU2::RME identityOperatorRMEPP, identityOperatorRMENN;

   uint16_t max_mult_p = bra.getMaximalMultiplicity_p();
   uint16_t max_mult_n = bra.getMaximalMultiplicity_n();

   InitializeIdenticalOperatorRME(identityOperatorRMEPP, max_mult_p*max_mult_p);
   InitializeIdenticalOperatorRME(identityOperatorRMENN, max_mult_n*max_mult_n);

   double mpi_tensor_end_setup_time = MPI_Wtime();

   // main block-by-block calculation
   SingleDistributionSmallVector distr_ip, distr_in, distr_jp, distr_jn;
   UN::SU3xSU2_VEC gamma_ip, gamma_in, gamma_jp, gamma_jn;
   SU3xSU2_SMALL_VEC vW_ip, vW_in, vW_jp, vW_jn;

   unsigned char deltaP, deltaN;

   const uint32_t number_ipin_blocks = bra.NumberOfBlocks();
   const uint32_t number_jpjn_blocks = ket.NumberOfBlocks();

   int32_t icurrentDistr_p, icurrentDistr_n;
   int32_t icurrentGamma_p, icurrentGamma_n;

   int32_t row_irrep = 0;
   
   for (unsigned int ipin_block = 0; ipin_block < number_ipin_blocks; ipin_block++) 
   {
      //    We don't have to check whether bra.NumberOfStatesInBlock(ipin_block) > 0
      // since we are interested
      //    in reduced matrix elements. The value of J is irrelevant.
      uint32_t ip = bra.getProtonIrrepId(ipin_block);
      uint32_t in = bra.getNeutronIrrepId(ipin_block);
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

      //  Generally, ket space can be different from bra space. For this reason we loop over
      // all ket jpjn blocks.
      int32_t col_irrep = 0;
      // std::cout<<"num jpjn_blocks "<<number_jpjn_blocks<<std::endl;
      for (unsigned int jpjn_block = 0; jpjn_block < number_jpjn_blocks; jpjn_block++) 
      {
         uint32_t jp = ket.getProtonIrrepId(jpjn_block);
         uint32_t jn = ket.getNeutronIrrepId(jpjn_block);
         int32_t Nhw_ket = ket.nhw_p(jp) + ket.nhw_n(jn);
         uint16_t ncol_irreps = ket.NumberPNIrrepsInBlock(jpjn_block);
         // std::cout<<"Ns "<<Nhw_bra<<" "<<dN0<<"  "<<Nhw_ket<<std::endl;
         if ((Nhw_ket + dN0) != Nhw_bra) 
         {
            col_irrep = col_irrep + ncol_irreps;
            continue;
         }

         SU3xSU2::LABELS w_jp(ket.getProtonSU3xSU2(jp));
         SU3xSU2::LABELS w_jn(ket.getNeutronSU3xSU2(jn));

         uint16_t ajp_max = ket.getMult_p(jp);
         uint16_t ajn_max = ket.getMult_n(jn);

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

         // std::cout<<"here"<<std::endl;
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

         // loop over wpn that result from coupling ip x in
         int32_t ibegin = bra.blockBegin(ipin_block);
         int32_t iend = bra.blockEnd(ipin_block);
         // std::cout<<"loop over wpn that result from coupling ip x in"<<std::endl;
         for (int32_t iwpn = ibegin, i = 0; iwpn < iend; ++iwpn, ++i) 
         {
            SU3xSU2::LABELS omega_pn_I(bra.getOmega_pn(ip, in, iwpn));
            IRREPBASIS braSU3xSU2basis(bra.Get_Omega_pn_Basis(iwpn));

            int32_t jbegin = ket.blockBegin(jpjn_block);
            int32_t jend = ket.blockEnd(jpjn_block);
            for (int32_t jwpn = jbegin, j = 0; jwpn < jend; ++jwpn, ++j) 
            {
               SU3xSU2::LABELS omega_pn_J(ket.getOmega_pn(jp, jn, jwpn));

              // std::cout<<fmt::format("({} {})  ({} {})",
              //   omega_pn_J.lm,omega_pn_J.mu,omega_pn_I.lm,omega_pn_I.mu
              //   )<<std::endl;
              // std::cout<<fmt::format("({} {})",w0.lm,w0.mu)<<std::endl;
              // std::cout<<"mult "<<SU3::mult(omega_pn_J, w0, omega_pn_I)<<std::endl;
              // std::cout<<fmt::format("spins: {} x {} -> {}",
              //   omega_pn_J.S2,w0.S2,omega_pn_I.S2)<<std::endl;
              // std::cout<<"su2 mult "<<SU2::mult(omega_pn_J.S2, w0.S2, omega_pn_I.S2)<<std::endl;
               if (!SU2::mult(omega_pn_J.S2, w0.S2, omega_pn_I.S2) ||
                   !SU3::mult(omega_pn_J, w0, omega_pn_I)) 
               {
                  continue;
               }

               IRREPBASIS ketSU3xSU2basis(ket.Get_Omega_pn_Basis(jwpn));

               if (deltaP + deltaN <= 4) 
               {
                  Reset_rmeCoeffs(rmeCoeffsPNPN);
                  if (in == jn && !rmeCoeffsPP.empty()) 
                    {
                      // create structure with < an (lmn mun)Sn ||| 1 ||| an' (lmn mun) Sn> r.m.e.
                      CreateIdentityOperatorRME(w_in, w_jn, ajn_max, identityOperatorRMENN);
                      // calculate < {(lmp mup)Sp x (lmn mun)Sn} wf Sf ||| [I_{nn} x T_{pp}]
                      //|||{(lmp' mup')Sp' x (lmn' mun')Sn'} wi Si >_{rhot}
                      Calculate_Proton_x_Identity_MeData(
                        omega_pn_I, omega_pn_J, rmeCoeffsPP,
                        identityOperatorRMENN, rmeCoeffsPNPN
                      );
                    }

                  if (ip == jp && !rmeCoeffsNN.empty()) 
                    {
                       // create structure with < ap (lmp mup)Sp ||| 1 ||| ap' (lmp mup) Sp> r.m.e.
                       CreateIdentityOperatorRME(w_ip, w_jp, ajp_max, identityOperatorRMEPP);
                       // calculate < {(lmp mup)Sp x (lmn mun)Sn} wf Sf ||| [T_{nn} x I_{pp}]
                       //|||{(lmp' mup')Sp' x (lmn' mun')Sn'} wi Si >_{rhot}
                       Calculate_Identity_x_Neutron_MeData(
                           omega_pn_I, omega_pn_J, identityOperatorRMEPP, rmeCoeffsNN, rmeCoeffsPNPN);
                    }

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

                  if (!rmeCoeffsPNPN.empty()) 
                    {
                      std::vector<double> rmes = ComputeRME(rmeCoeffsPNPN);

                      // check if any rmes are above zero threshold
                      bool nonzero = false;
                      for (auto rme : rmes) 
                       {
                          // Caution: Important to use std::abs() rather than integer abs().
                          nonzero |= (std::abs(rme) > g_zero_threshold);
                       }

                      // count rmes
                      num_rmes_allowed += rmes.size();

                      if (nonzero) 
                       {
                          num_groups_nonzero += 1;
                          num_rmes_nonzero += rmes.size();
                       }

                      int row_index = row_irrep + i;
                      assert(row_index == static_cast<RMEIndexType>(row_index));
                      WriteBinary<RMEIndexType>(out_file, row_index);
                      int col_index = col_irrep + j;
                      assert(col_index == static_cast<RMEIndexType>(col_index));
                      WriteBinary<RMEIndexType>(out_file, col_index);
                      assert(rmes.size() == static_cast<RMEIndexType>(rmes.size()));
                      WriteBinary<RMEIndexType>(out_file, rmes.size());
                    }
               }
            }
         }
         col_irrep = col_irrep + ncol_irreps;
         last_jp = jp;
         last_jn = jn;
      }
      row_irrep = row_irrep + nrow_irreps;
   }

   double mpi_tensor_end_iteration_time = MPI_Wtime();

   Reset_rmeCoeffs(rmeCoeffsPP);
   Reset_rmeCoeffs(rmeCoeffsNN);
   Reset_rmeCoeffs(rmeCoeffsPNPN);
   // delete arrays allocated in identityOperatorRME?? structures
   delete[] identityOperatorRMEPP.m_rme;
   delete[] identityOperatorRMENN.m_rme;
 
    out_file.close();

}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{ 
  auto run_parameters = ProcessArguments(argc, argv);
  bool intrinsic = true;
  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(run_parameters.nuclide,intrinsic);
  const auto&[Z,N]=run_parameters.nuclide;
  std::string spncci_root_dir=get_spncci_project_root_dir();
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //Generate LGI vector by finding possible cmf LGI by counting arguments 
  lgi::MultiplicityTaggedLGIVector lgi_vector = lgi::get_lgi_vector(run_parameters.nuclide,Nsigma0,run_parameters.Nmax);
  std::cout<<"num lgi "<<lgi_vector.size()<<std::endl;
  std::string operator_base_name = fmt::format("{}/spncci/data/relative_operators/Brel",spncci_root_dir);
  std::cout<<operator_base_name<<std::endl;
  
  //////////////////////////////////////////////////////////////////////////////////////////
  // For some reason this cannot be part of a function
  //////////////////////////////////////////////////////////////////////////////////////////
  MPI_Init(&argc, &argv);

  su3::init();
  CWig9lmLookUpTable<RME::DOUBLE>::AllocateMemory(true);

  CBaseSU3Irreps baseSU3Irreps(Z,N,run_parameters.Nmax);

  std::ofstream interaction_log_file("/dev/null");
  bool log_is_on = false;
  bool generate_missing_rme = true;
  std::string ppnn_file_name(operator_base_name+".PPNN");
  std::string pn_file_name(operator_base_name+".PN");

  CInteractionPPNN interactionPPNN(baseSU3Irreps,log_is_on,interaction_log_file);
  CInteractionPN interactionPN(baseSU3Irreps,generate_missing_rme,log_is_on,interaction_log_file);
  interactionPPNN.LoadTwoBodyOperator(ppnn_file_name);
  interactionPN.AddOperator(pn_file_name);
  interactionPPNN.TransformTensorStrengthsIntoPP_NN_structure();
  //////////////////////////////////////////////////////////////////////////////////////////   

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Calculate rmes of Brel and Ncm for each U3SpSnS subspace corresponding to an LGI
  for(const auto&[lgi,dim] : lgi_vector)
    {  

      //Get labels
      const auto&[Nex,sigma,Sp,Sn,S]=lgi.Key();
      
      if(Nex==0)
        continue;
 
      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Setting up basis
      ////////////////////////////////////////////////////////////////////////////////////////////////
      // Temporary container for constructing the bra and ket model space 
      // model_space_map organized into {N : {SpSn : {S : <(lambda,mu)>}}} 
      // SU3_VEC is std::vector<UN::SU3> where UN::SU3(mutiplicity,lambda,mu) 
      //   where multiplicity is U(N)->U(3) braching multiplicity.
      //   Defined in lsu3shell/libraries/UNU3SU3/UNU3SU3basics.h
      std::map<int,std::map<std::pair<int,int>,std::map<int,SU3_VEC>>> model_space_map_ket;
      
      // Add irrep corresponding to lgi to map
      model_space_map_ket[Nex][{TwiceValue(Sp),TwiceValue(Sn)}][TwiceValue(S)].emplace_back(1,sigma.SU3().lambda(),sigma.SU3().mu());

      // Ket model space consists of single U(3)SpSnS subspace corresponding to lgi of interest
      proton_neutron::ModelSpace ket_ncsmModelSpace(Z,N,model_space_map_ket);
      
      // bra model space consists of U(3)SpSnS subspace corresponding to LGI plus subspaces connected by B[-2(0,2)]. 
      // Iterate over lgi_vector and add connected subspaces to map
      std::map<int,std::map<std::pair<int,int>,std::map<int,SU3_VEC>>> model_space_map_bra;
      for(const auto&[lgi_p,dim_p] : lgi_vector)
        {
          const auto&[Nex_p,sigma_p,Sp_p,Sn_p,S_p]=lgi_p.Key();
          
          // Apply selection rule 
          if(not Nex-2==Nex_p){continue;}
          if(not Sp==Sp_p){continue;}
          if(not Sn==Sn_p){continue;}
          if(not S==S_p){continue;}          
          if(u3::OuterMultiplicity(sigma.SU3(),u3::SU3(0,2),sigma_p.SU3())==0) {continue;}
          
          //If passed all selection rules, add to bra spaces 
          SU3::LABELS w(1,sigma_p.SU3().lambda(),sigma_p.SU3().mu());
          model_space_map_bra[Nex_p][{TwiceValue(Sp_p),TwiceValue(Sn_p)}][TwiceValue(S_p)].push_back(w);
        }
      std::cout<<"bra model space "<<std::endl;
      for(const auto& [N,temp1] : model_space_map_bra)
        for(const auto& [SPN,temp2] : temp1)
          for(const auto& [S,su3_vec] : temp2)
            for(const auto& [rho,lm,mu] : su3_vec)
            {
              std::cout<<fmt::format("    {} ({} {}) {}",N,lm,mu,S)<<std::endl;
            }

      // Generate bra model space with augment model_space_map
      proton_neutron::ModelSpace bra_ncsmModelSpace(Z,N,model_space_map_bra);
      
      // std::cout<<"construct ket"<<std::endl;
      // Construct ket basis from model space
      int jdiag=0; 
      int ndiag=1;
      lsu3::CncsmSU3xSU2Basis ket;
      ket.ConstructBasis(ket_ncsmModelSpace, jdiag, ndiag);

      // Construct bra basis from model space
      // int jdiag=0; ndiag=1;
      // std::cout<<"construct bra"<<std::endl;
      lsu3::CncsmSU3xSU2Basis bra;
      bra.ConstructBasis(bra_ncsmModelSpace, jdiag, ndiag);
      //////////////////////////////////////////////////////////////////////////////////////////////
      // Actual calculation
      ////////////////////////////////////////////////////////////////////////////////////////////////
      /*  Inner most loop of SU3RME calculates rmes for given pair of [omegap Sp omegan Sn S] omega S subspace
          currently by bra upstream, by ket upstream, by rho0.  But I plan to rewrite so it's 
          by rho0, by bra upstream, by ket upstream instead.
          
          Will likely pass sector blocks to rme calculator to directly accumulate in larger block

          To get dimension of big blocks, probably need to iterate over blocks and get sizes
          
          for bra proton-neutron block
            for ket proton-neutron block
                for bra omega 
                    for ket omega
                        for rho0

                          rho1  rho2  rho3  rho4
                          xxxx  xxxx  xxxx  xxxx
                          xxxx  xxxx  xxxx  xxxx
                          xxxx  xxxx  xxxx  xxxx
                          xxxx  xxxx  xxxx  xxxx

          then do null solver and or transform unit tensor rmes 

          But, if Brel+Ncm, then single omega Sp Sn S ket subspace and if unit tensor then
          single omega Sp Sn S bra subspace as well.

          For given ipin,jpjn tile, need to identify relevant section of big matrix
            -> This picks out given omegap Sp omegan Sn omega S for bra and ket
            -> Size of tile is multi x multj
            -> Different tile for each rho0 (referred to as rho_t in Tomas' code)
      */

     std::string rme_file_name(operator_base_name+".rme");
     int dN0=-2;
     // rho0, lambda0,mu0,twice_S0
     SU3xSU2::LABELS w0(1,0,2,0);
     CalculateRME(rme_file_name, interactionPPNN,interactionPN,bra,ket,dN0,w0);
    }
  // clears memory allocated for U9, U6, and Z6 coefficients
  CWig9lmLookUpTable<RME::DOUBLE>::ReleaseMemory();
  // clear memory allocated for single-shell SU(3) rmes
  CSSTensorRMELookUpTablesContainer::ReleaseMemory();  
  // clear memory allocated for SU(3)>SO(3)
  CWigEckSU3SO3CGTablesLookUpContainer::ReleaseMemory();  

  su3::finalize();

MPI_Finalize();
}//end main





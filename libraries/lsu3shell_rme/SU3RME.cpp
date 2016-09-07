#include <SU3ME/ComputeOperatorMatrix.h>
#include <LSU3/ncsmSU3xSU2Basis.h>
#include <SU3ME/proton_neutron_ncsmSU3BasisFastIteration.h>
#include <SU3NCSMUtils/CRunParameters.h>
#include <LookUpContainers/CWig9lmLookUpTable.h>
#include <LookUpContainers/lock.h>

#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <stack>

#include <boost/mpi.hpp>
#include <boost/chrono.hpp>
//	To be able to load and distribute basis from a binary archive file
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

using namespace std;

void CheckRMEs(const std::vector<MECalculatorData>& rmeCoeffsPNPN)
{
   const size_t ntotal = rmeCoeffsPNPN[0].rmes_->m_ntotal;
   for (size_t itensor = 1; itensor < rmeCoeffsPNPN.size(); ++itensor) {
      if (ntotal != rmeCoeffsPNPN[itensor].rmes_->m_ntotal) {
         std::cerr << "Error: the number of rmes for a given set of tensors is not constant!";
         exit(EXIT_FAILURE);
      }
   }
}

std::vector<double> ComputeRME(std::vector<MECalculatorData>& rmeCoeffsPNPN) {
   // These values must be constant for each tensor
   int32_t bra_max = rmeCoeffsPNPN[0].rmes_->m_bra_max;
   int32_t ket_max = rmeCoeffsPNPN[0].rmes_->m_ket_max;
   int32_t rhot_max = rmeCoeffsPNPN[0].rmes_->m_rhot_max;
   // tensor multiplicity is however tensor dependent
   int32_t rho0_max;

   // resulting rmes do not have rho0 dependence
   int32_t nrmes_to_compute = bra_max * ket_max * rhot_max;
   vector<double> result_rmes(nrmes_to_compute, 0.0);

   // iterate over PP NN and PN tensors each has the same lm0 mu0 and S0
   // quantum labels but it can have a different rho0max multiplicity
   for (size_t itensor = 0; itensor < rmeCoeffsPNPN.size(); ++itensor) {
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

void CalculateRME(const CInteractionPPNN& interactionPPNN, const CInteractionPN& interactionPN,
                  const lsu3::CncsmSU3xSU2Basis& bra, const lsu3::CncsmSU3xSU2Basis& ket) {
   std::vector<unsigned char> hoShells_n, hoShells_p;
   std::vector<CTensorGroup *> tensorGroupsPP, tensorGroupsNN;
   std::vector<CTensorGroup_ada *> tensorGroups_p_pn, tensorGroups_n_pn;

   std::vector<int> phasePP, phaseNN, phase_p_pn, phase_n_pn;

   unsigned char num_vacuums_J_distr_p;
   unsigned char num_vacuums_J_distr_n;
   std::vector<std::pair<CRMECalculator *, CTensorGroup::COEFF_DOUBLE *> > selected_tensorsPP,
       selected_tensorsNN;
   std::vector<std::pair<CRMECalculator *, unsigned int> > selected_tensors_p_pn,
       selected_tensors_n_pn;

   std::vector<RmeCoeffsSU3SO3CGTablePointers> rmeCoeffsPP, rmeCoeffsNN;
   std::vector<std::pair<SU3xSU2::RME *, unsigned int> > rme_index_p_pn, rme_index_n_pn;

   std::vector<MECalculatorData> rmeCoeffsPNPN;

   SU3xSU2::RME identityOperatorRMEPP, identityOperatorRMENN;

   InitializeIdenticalOperatorRME(identityOperatorRMEPP);
   InitializeIdenticalOperatorRME(identityOperatorRMENN);

   SingleDistribution distr_ip, distr_in, distr_jp, distr_jn;
   UN::SU3xSU2_VEC gamma_ip, gamma_in, gamma_jp, gamma_jn;
   SU3xSU2_VEC vW_ip, vW_in, vW_jp, vW_jn;

   unsigned char deltaP, deltaN;

   const uint32_t number_ipin_blocks = bra.NumberOfBlocks();
   const uint32_t number_jpjn_blocks = ket.NumberOfBlocks();

   int32_t icurrentDistr_p, icurrentDistr_n;
   int32_t icurrentGamma_p, icurrentGamma_n;

   int32_t row_irrep = 0;
   for (unsigned int ipin_block = 0; ipin_block < number_ipin_blocks; ipin_block++) {
      //		We don't have to check whether bra.NumberOfStatesInBlock(ipin_block) > 0
      //since we are interested
      //		in reduced matrix elements. The value of J is irrelevant.
      uint32_t ip = bra.getProtonIrrepId(ipin_block);
      uint32_t in = bra.getNeutronIrrepId(ipin_block);
      uint16_t nrow_irreps = bra.NumberPNIrrepsInBlock(ipin_block);

      SU3xSU2::LABELS w_ip(bra.getProtonSU3xSU2(ip));
      SU3xSU2::LABELS w_in(bra.getNeutronSU3xSU2(in));

      uint16_t ilastDistr_p(std::numeric_limits<uint16_t>::max());
      uint16_t ilastDistr_n(std::numeric_limits<uint16_t>::max());

      uint32_t ilastGamma_p(std::numeric_limits<uint32_t>::max());
      uint32_t ilastGamma_n(std::numeric_limits<uint32_t>::max());

      uint32_t last_jp(std::numeric_limits<uint32_t>::max());
      uint32_t last_jn(std::numeric_limits<uint32_t>::max());

      //	Generally, ket space can be different from bra space. For this reason we loop over
      //all ket jpjn blocks.
      int32_t col_irrep = 0;
      for (unsigned int jpjn_block = 0; jpjn_block < number_jpjn_blocks; jpjn_block++) {
         uint32_t jp = ket.getProtonIrrepId(jpjn_block);
         uint32_t jn = ket.getNeutronIrrepId(jpjn_block);
         uint16_t ncol_irreps = ket.NumberPNIrrepsInBlock(jpjn_block);

         SU3xSU2::LABELS w_jp(ket.getProtonSU3xSU2(jp));
         SU3xSU2::LABELS w_jn(ket.getNeutronSU3xSU2(jn));

         uint16_t ajp_max = ket.getMult_p(jp);
         uint16_t ajn_max = ket.getMult_n(jn);

         if (jp != last_jp) {
            icurrentDistr_p = ket.getIndex_p<lsu3::CncsmSU3xSU2Basis::kDistr>(jp);
            icurrentGamma_p = ket.getIndex_p<lsu3::CncsmSU3xSU2Basis::kGamma>(jp);

            if (ilastDistr_p != icurrentDistr_p) {
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

            if (ilastGamma_p != icurrentGamma_p || ilastDistr_p != icurrentDistr_p) {
               if (deltaP <= 4) {
                  if (!selected_tensorsPP.empty()) {
                     std::for_each(selected_tensorsPP.begin(), selected_tensorsPP.end(),
                                   CTensorGroup::DeleteCRMECalculatorPtrs());
                  }
                  selected_tensorsPP.resize(0);

                  if (!selected_tensors_p_pn.empty()) {
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

            if (deltaP <= 4) {
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

         if (jn != last_jn) {
            icurrentDistr_n = ket.getIndex_n<lsu3::CncsmSU3xSU2Basis::kDistr>(jn);
            icurrentGamma_n = ket.getIndex_n<lsu3::CncsmSU3xSU2Basis::kGamma>(jn);

            if (ilastDistr_n != icurrentDistr_n) {
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

            if (ilastGamma_n != icurrentGamma_n || ilastDistr_n != icurrentDistr_n) {
               if (deltaN <= 4) {
                  if (!selected_tensorsNN.empty()) {
                     std::for_each(selected_tensorsNN.begin(), selected_tensorsNN.end(),
                                   CTensorGroup::DeleteCRMECalculatorPtrs());
                  }
                  selected_tensorsNN.resize(0);

                  if (!selected_tensors_n_pn.empty()) {
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

            if (deltaN <= 4) {
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

         //	loop over wpn that result from coupling ip x in
         int32_t ibegin = bra.blockBegin(ipin_block);
         int32_t iend = bra.blockEnd(ipin_block);
         for (int32_t iwpn = ibegin, i = 0; iwpn < iend; ++iwpn, ++i) {
            SU3xSU2::LABELS omega_pn_I(bra.getOmega_pn(ip, in, iwpn));
            IRREPBASIS braSU3xSU2basis(bra.Get_Omega_pn_Basis(iwpn));

            int32_t jbegin = ket.blockBegin(jpjn_block);
            int32_t jend = ket.blockEnd(jpjn_block);
            for (int32_t jwpn = jbegin, j = 0; jwpn < jend; ++jwpn, ++j) {
               SU3xSU2::LABELS omega_pn_J(ket.getOmega_pn(jp, jn, jwpn));
               IRREPBASIS ketSU3xSU2basis(ket.Get_Omega_pn_Basis(jwpn));

               if (deltaP + deltaN <= 4) {
                  Reset_rmeCoeffs(rmeCoeffsPNPN);
                  if (in == jn && !rmeCoeffsPP.empty()) {
                     //	create structure with < an (lmn mun)Sn ||| 1 ||| an' (lmn mun) Sn> r.m.e.
                     CreateIdentityOperatorRME(w_in, w_jn, ajn_max, identityOperatorRMENN);
                     //	calculate <	{(lmp mup)Sp x (lmn mun)Sn} wf Sf ||| [I_{nn} x T_{pp}]
                     //|||{(lmp' mup')Sp' x (lmn' mun')Sn'} wi Si >_{rhot}
                     Calculate_Proton_x_Identity_MeData(omega_pn_I, omega_pn_J, rmeCoeffsPP,
                                                        identityOperatorRMENN, rmeCoeffsPNPN);
                  }

                  if (ip == jp && !rmeCoeffsNN.empty()) {
                     //	create structure with < ap (lmp mup)Sp ||| 1 ||| ap' (lmp mup) Sp> r.m.e.
                     CreateIdentityOperatorRME(w_ip, w_jp, ajp_max, identityOperatorRMEPP);
                     //	calculate <	{(lmp mup)Sp x (lmn mun)Sn} wf Sf ||| [T_{nn} x I_{pp}]
                     //|||{(lmp' mup')Sp' x (lmn' mun')Sn'} wi Si >_{rhot}
                     Calculate_Identity_x_Neutron_MeData(
                         omega_pn_I, omega_pn_J, identityOperatorRMEPP, rmeCoeffsNN, rmeCoeffsPNPN);
                  }

                  if (!rme_index_p_pn.empty() && !rme_index_n_pn.empty()) {
                     CalculatePNInteractionMeData(interactionPN, omega_pn_I, omega_pn_J,
                                                  rme_index_p_pn, rme_index_n_pn, rmeCoeffsPNPN);
                  }

                  if (!rmeCoeffsPNPN.empty()) {
                     vector<double> rmes = ComputeRME(rmeCoeffsPNPN);
                     std::cout << row_irrep + i << " " << col_irrep + j << " ";
                     for (auto rme : rmes)
                     {
                       std::cout << rme << " ";
                     }
                     std::cout << "size: " << rmes.size() << std::endl;
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
   Reset_rmeCoeffs(rmeCoeffsPP);
   Reset_rmeCoeffs(rmeCoeffsNN);
   Reset_rmeCoeffs(rmeCoeffsPNPN);
   //	delete arrays allocated in identityOperatorRME?? structures
   delete[] identityOperatorRMEPP.m_rme;
   delete[] identityOperatorRMENN.m_rme;
}

void CalculateRMEs(const std::string& bra_space_definition_file_name,
                   const std::string& ket_space_definition_file_name,
                   const std::string& hamiltonian_file_name) {
   int idiag = 0;
   int jdiag = 0;
   int ndiag = 1;

   InitSqrtLogFactTables();

   //	Construct BRA model space
   proton_neutron::ModelSpace bra_ncsmModelSpace;
   lsu3::CncsmSU3xSU2Basis bra;
   bra_ncsmModelSpace.Load(bra_space_definition_file_name);
   bra.ConstructBasis(bra_ncsmModelSpace, idiag, ndiag);

   //	Construct KET model space
   proton_neutron::ModelSpace ket_ncsmModelSpace;
   lsu3::CncsmSU3xSU2Basis ket;
   ket_ncsmModelSpace.Load(ket_space_definition_file_name);
   ket.ConstructBasis(ket_ncsmModelSpace, jdiag, ndiag);

   //	stringstream interaction_log_file_name;
   //	interaction_log_file_name << "interaction_loading_" << my_rank << ".log";
   //	ofstream interaction_log_file(interaction_log_file_name.str().c_str());
   ofstream interaction_log_file("/dev/null");

   int nmax = std::max(bra.Nmax(), ket.Nmax());

   CBaseSU3Irreps baseSU3Irreps(bra.NProtons(), bra.NNeutrons(), nmax);

   //	since root process will read rme files, it is save to let this
   //	process to create missing rme files (for PPNN interaction) if they do not exist
   bool log_is_on = false;
   CInteractionPPNN interactionPPNN(baseSU3Irreps, log_is_on, interaction_log_file);

   //	PN interaction is read after PPNN, and hence all rmes should be already in memory.
   //	if rme file does not exist, then false

   bool generate_missing_rme = true;
   CInteractionPN interactionPN(baseSU3Irreps, generate_missing_rme, log_is_on,
                                interaction_log_file);

   CRunParameters run_params;
   run_params.LoadHamiltonian(0, hamiltonian_file_name, interactionPPNN, interactionPN,
                              bra.NProtons() + bra.NNeutrons());

   //	interactionPN.ShowTensorGroupsAndTheirTensorStructures();

   interactionPPNN.TransformTensorStrengthsIntoPP_NN_structure();

   CalculateRME(interactionPPNN, interactionPN, bra, ket);
}

int main(int argc, char** argv) {
   boost::mpi::environment env(argc, argv);

   MECalculatorData dummy;
   int jj0 = dummy.JJ0();
   int mm0 = dummy.MM0();

   if (jj0 != 0 || mm0 != 0) {
      return EXIT_FAILURE;
   }

   boost::chrono::system_clock::time_point start;
   boost::chrono::duration<double> duration;
   int my_rank, nprocs;

   boost::mpi::communicator mpi_comm_world;
   my_rank = mpi_comm_world.rank();
   nprocs = mpi_comm_world.size();

   if (nprocs != 1) {
      if (my_rank == 0) {
         cerr << "Only one MPI process allowed at this moment!" << endl;
         mpi_comm_world.abort(EXIT_FAILURE);
      }
   }

   if (su3shell_data_directory == NULL) {
      if (my_rank == 0) {
         cerr << "System variable 'SU3SHELL_DATA' was not defined!" << endl;
      }
      mpi_comm_world.abort(EXIT_FAILURE);
   }

   if (argc != 4) {
      if (my_rank == 0) {
         cerr << "Usage: " << argv[0] << " <bra model space>  <ket model space> <operator>" << endl;
      }
      mpi_comm_world.abort(EXIT_FAILURE);
   }

   std::string bra_space_definition_file_name(argv[1]);
   std::string ket_space_definition_file_name(argv[2]);
   std::string hamiltonian_file_name(argv[3]);

   //	these two variable are needed for MFDn
   CalculateRMEs(bra_space_definition_file_name, ket_space_definition_file_name,
                 hamiltonian_file_name);

   CWig9lmLookUpTable<RME::DOUBLE>::ReleaseMemory();  // clears memory allocated for U9, U6, and Z6
                                                      // coefficients
   CSSTensorRMELookUpTablesContainer::ReleaseMemory();  // clear memory allocated for single-shell
                                                        // SU(3) rmes
   CWigEckSU3SO3CGTablesLookUpContainer::ReleaseMemory();  // clear memory allocated for SU(3)>SO(3)
                                                           // coefficients

   return EXIT_SUCCESS;
}

/****************************************************************
  unit_tensor.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/unit_tensor.h"

#include <omp.h>

#include "cppformat/format.h"
#include "mcutils/eigen.h"
#include "sp3rlib/u3coef.h"
#include "spncci/spncci_common.h"

extern double zero_threshold;

namespace spncci
{
void ZeroInitBlocks(int number, int rows, int cols,std::vector<basis::OperatorBlock<double>>& unit_tensor_blocks)
  {
    unit_tensor_blocks.resize(number);
    for(int i=0; i<number; ++i)
      unit_tensor_blocks[i]=Eigen::MatrixXd::Zero(rows,cols);
  }

// spncci::UnitTensorMatrixStatistics
//   GenerateUnitTensorMatrixStatistics(const spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices)
// {

//   // counters
//   UnitTensorMatrixStatistics statistics;

//   // level 1: traverse irrep family index pairs
//   for(const auto& irrep_family_index_pairs_entry : unit_tensor_matrices)
//     {
//       ++statistics.num_irrep_family_index_pairs;

//       // level 2: traverse Nn pairs
//       for(const auto& Nn_pairs_entry : irrep_family_index_pairs_entry.second)
//         {
//           ++statistics.num_Nn_pairs;

//           // level 3: traverse unit tensor sectors
//           for(const auto& unit_tensor_sectors_entry : Nn_pairs_entry.second)
//             {
//               ++statistics.num_unit_tensor_sectors;
              
//               // process sector matrix
//               const Eigen::MatrixXd& matrix = unit_tensor_sectors_entry.second;
//               if (not mcutils::IsZero(matrix,spncci::g_zero_tolerance))
//                 ++statistics.num_nonzero_unit_tensor_sectors;
//               statistics.num_matrix_elements += matrix.rows()*matrix.cols();
//               statistics.num_nonzero_matrix_elements += mcutils::NonzerosInMatrix(matrix,spncci::g_zero_tolerance); 
//             }
//         }
//     }
//   return statistics;
// }

void 
ComputeUnitTensorHyperblocks(
  int Nmax, int N1v,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::KMatrixCache& k_matrix_map,
  const spncci::SpNCCISpace& spncci_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
  const std::vector<std::vector<int>>& unit_tensor_hypersector_subsets,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  )
// compute hyperblocks for unit tensors recursively
{
  // std::cout<<"in the recurrence"<<std::endl;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Set up for calculation 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for(int Nsum=2; Nsum<=2*Nmax; Nsum+=2)
    {
      // std::cout<< "Nsum "<<Nsum<<std::endl;
      const std::vector<int>& unit_tensor_hypersectors=unit_tensor_hypersector_subsets[Nsum/2];
      // std::cout<<"for each of the "<<unit_tensor_hypersectors.size()<<" hypersectors "<<std::endl;
      for(int hypersector_index : unit_tensor_hypersectors)
        {
          // std::cout<<"baby_spncci_hypersectors.GetHypersector(hypersector_index)"<<std::endl;
          auto key=baby_spncci_hypersectors.GetHypersector(hypersector_index).Key();

          int unit_tensor_subspace_index, baby_spncci_subspace_indexp, baby_spncci_index, rho0;
          std::tie(baby_spncci_subspace_indexp,baby_spncci_index,unit_tensor_subspace_index,rho0)=key;
          
          const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
              =baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);

          const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
              =baby_spncci_space.GetSubspace(baby_spncci_index);
          
          const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace
              =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
          
          // If Nnp<Nn, will need to look up and use conjugate sectors 
          // std::cout<<"Nnp, Nn "<<baby_spncci_subspace_bra.Nn()<<"  "<<baby_spncci_subspace_ket.Nn()<<std::endl;
          bool conjugate=(baby_spncci_subspace_bra.Nn()>(baby_spncci_subspace_ket.Nn()-2));
          // std::cout<<"conjugate "<<conjugate<<std::endl;
          int dim=baby_spncci_subspace_ket.size();
          int gamma_max=baby_spncci_subspace_ket.gamma_max();
          int upsilon_max=baby_spncci_subspace_ket.upsilon_max();

          int dimp=baby_spncci_subspace_bra.size();
          int gamma_maxp=baby_spncci_subspace_bra.gamma_max();
          int upsilon_maxp=baby_spncci_subspace_bra.upsilon_max();
          // extract Sp(3,R) space
          int irrep_family_index_bra=baby_spncci_subspace_bra.irrep_family_index();
          int irrep_family_index_ket=baby_spncci_subspace_ket.irrep_family_index();

          const sp3r::Sp3RSpace& irrep_bra = spncci_space[irrep_family_index_bra].Sp3RSpace();
          const sp3r::Sp3RSpace& irrep_ket = spncci_space[irrep_family_index_ket].Sp3RSpace();

          // extract subspace labels 
          u3::U3 omegap,sigmap,omega,sigma;
          u3::SU3 x0;
          HalfInt S0,Sn_ket,Sp_ket,S_ket,Sn_bra,Sp_bra,S_bra;
          int etap,eta;

          // Extracting labels
          std::tie(sigmap,Sp_bra,Sn_bra,S_bra,omegap)=baby_spncci_subspace_bra.labels();
          std::tie(sigma,Sp_ket,Sn_ket,S_ket,omega)=baby_spncci_subspace_ket.labels();
          std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
          int Nn=baby_spncci_subspace_ket.Nn();

          // omega u3 subspace in irrep
          const sp3r::U3Subspace& u3_subspace=irrep_ket.LookUpSubspace(omega);
          const sp3r::U3Subspace& u3_subspacep=irrep_bra.LookUpSubspace(omegap);
          // Extracting K matrices for sp_irrep and sp_irrepp from the K_matrix_maps 
          const vcs::MatrixCache& K_matrix_map_bra=k_matrix_map.at(sigmap);
          const vcs::MatrixCache& K_matrix_map_ket=k_matrix_map.at(sigma);

          Eigen::MatrixXd Kp=K_matrix_map_bra.at(omegap);
          Eigen::MatrixXd K_inv=K_matrix_map_ket.at(omega).inverse();

          // Generate labels to sum over 
          int rho0_max=u3::OuterMultiplicity(omega.SU3(),x0,omegap.SU3());

          // Precalculating kronecker products used in sum to calculate unit tensor matrix
          MultiplicityTagged<u3::U3>::vector omegapp_set=KroneckerProduct(omegap, u3::U3(0,0,-2)); 
          MultiplicityTagged<u3::U3>::vector omega1_set=KroneckerProduct(omega, u3::U3(0,0,-2));
          MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(x0, u3::SU3(2,0));

          double conjugation_grade=1;
          if(conjugate)
            {
            conjugation_grade*=ParitySign(
                u3::ConjugationGrade(omegap)+S_bra
                +u3::ConjugationGrade(omega)+S_ket
                +u3::ConjugationGrade(x0)
              );
            // std::cout<<"conjugation_grade "<<conjugation_grade<<std::endl;
            }
          std::vector<basis::OperatorBlock<double>>& unit_tensor_blocks=unit_tensor_hyperblocks[hypersector_index];
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //  Calculate unit tensor matrix
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // std::cout<<"begin tensor calculation "<<std::endl;
          int num_blocks=unit_tensor_blocks.size();

          for(auto& omega1_mult :omega1_set)
            {
              u3::U3 omega1(omega1_mult.irrep);
              // std::cout<<"omega1 "<<omega1.Str()<<" of "<<omega1_set.size()<<std::endl;
              // std::cout<<"irreps "<<irrep_ket.sigma().Str()<<std::endl;
              if (not irrep_ket.ContainsSubspace(omega1))
                continue;
              
              spncci::BabySpNCCISubspaceLabels baby_spncci_labels1(sigma,Sp_ket,Sn_ket,S_ket,omega1);
              int baby_spncci_subspace_index1=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labels1);

              Eigen::MatrixXd K1=K_matrix_map_ket.at(omega1);
              const sp3r::U3Subspace& u3_subspace1=irrep_ket.LookUpSubspace(omega1);
              int upsilon_max1=u3_subspace1.size();

              int dim1=upsilon_max1*gamma_max;
              std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omega1;

              // std::cout<<"initializing blocks"<<std::endl;
              // Initializing blocks for sum over omega1
              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omega1);
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              // Construct KBUK matrix
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              // Input KBUK matrix
              // u3_subspace
              // ucoef cache 
              // get dims for BU from KBUK
              //
              // std::cout<<"begin constructing KBUK matrix "<<std::endl;
              Eigen::MatrixXd BU(upsilon_max1, upsilon_max);
              for(int u3_state_index=0; u3_state_index<upsilon_max; ++u3_state_index)
                {
                  MultiplicityTagged<u3::U3> n_rho=u3_subspace.GetStateLabels(u3_state_index);
                  u3::U3 n(n_rho.irrep);
                  // iterate over (n1,rho1)
                  for (int u3_state_index1=0; u3_state_index1<upsilon_max1; u3_state_index1++)
                    {
                      MultiplicityTagged<u3::U3> n1_rho1=u3_subspace1.GetStateLabels(u3_state_index1);
                      u3::U3 n1(n1_rho1.irrep);
                      if (u3::OuterMultiplicity(n1.SU3(), u3::SU3(2,0),n.SU3())>0)
                          BU(u3_state_index1,u3_state_index)
                            =2./Nn*vcs::BosonCreationRME(n,n1)
                             *u3::UCached(u_coef_cache,u3::SU3(2,0),n1.SU3(),omega.SU3(),sigma.SU3(),
                                n.SU3(),1,n_rho.tag,omega1.SU3(),n1_rho1.tag,1
                              );
                      else
                        BU(u3_state_index1,u3_state_index)=0;

                    }
                }
              Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
              KBUK.noalias()=K1*BU*K_inv;
              // std::cout<<"KBUK "<<KBUK<<std::endl;
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              //summing over x0'
              // std::cout<<"sum over x0"<<std::endl;
              for (auto& x0p_mult : x0p_set)
                {
                  u3::SU3 x0p(x0p_mult.irrep);
                  // std::cout<<"x0p"<<std::endl;
                  int rho0p_max=OuterMultiplicity(omega1.SU3(),x0p,omegap.SU3());
                  
                  // summing over rho0'
                  // std::cout<<"summing over rho0p.  rho0p_max "<<rho0p_max<<std::endl;
                  for (int rho0p=1; rho0p<=rho0p_max; rho0p++)
                    {
                      // Zero initialize blocks accumlating sum over x0p and rho0p
                      std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_x0p;
                      ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_x0p);

                      double coef=0;
                      for(int rho0b=1; rho0b<=rho0_max; rho0b++)
                        {
                          //(2,0)xx0->x0p(by construction), 
                          //(2,0)xomega1->omega (by construction),
                          // x0xomega->omegap, (rho0_max)
                          //omega1xx0p->omegap (rho0p_max)
                          coef+=u3::PhiCached(phi_coef_cache,omega.SU3(),x0,omegap.SU3(),rho0,rho0b)
                               *u3::UCached(u_coef_cache,x0,u3::SU3(2,0),omegap.SU3(), omega1.SU3(),x0p,1,rho0p,omega.SU3(),1,rho0b);
                        }
                      // std::cout<<"coef "<<coef<<std::endl;

                      // std::cout<<"computing term 3"<<std::endl;
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // third term
                      // sum over omega'', v'' and rho0''
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // Zero initialze blocks accumulating sum over omegapp
                      std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omegapp;
                      ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omegapp);

                      // Summing over omega''
                      for (auto& omegapp_mult : omegapp_set)
                        {
                          // A matrix will annihilate bra
                          if(baby_spncci_subspace_bra.Nn()==0)
                            continue;

                          u3::U3 omegapp(omegapp_mult.irrep);
                          // std::cout<<omegapp.Str()<<std::endl;

                          if (not irrep_bra.ContainsSubspace(omegapp))
                            continue;
                          
                          // get hypersector index
                          spncci::BabySpNCCISubspaceLabels baby_spncci_labelspp(sigmap,Sp_bra,Sn_bra,S_bra,omegapp);
                          int baby_spncci_subspace_indexpp=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labelspp);

                          // omega'' subspace (v'')
                          sp3r::U3Subspace u3_subspacepp=irrep_bra.LookUpSubspace(omegapp);
                          int upsilon_maxpp=u3_subspacepp.size();

                          // Obtaining K matrix for omega''
                          Eigen::MatrixXd Kpp_inv=K_matrix_map_bra.at(omegapp).inverse();
                          
                          //Constructing a^\dagger U(3) boson matrix for A matrix
                          Eigen::MatrixXd boson_matrix(upsilon_maxp,upsilon_maxpp);
                          for(int vpp=0; vpp<upsilon_maxpp; vpp++)
                            {
                              MultiplicityTagged<u3::U3> npp_rhopp=u3_subspacepp.GetStateLabels(vpp);
                              const u3::U3& npp(npp_rhopp.irrep);
                              int rhopp=npp_rhopp.tag;
                              for(int vp=0; vp<dimp; vp++)
                                {
                                  MultiplicityTagged<u3::U3> np_rhop=u3_subspacep.GetStateLabels(vp);
                                  const u3::U3& np(np_rhop.irrep);
                                  int rhop=np_rhop.tag; 
                                  // std::cout<<"A matrix"<<std::endl;
                                  // std::cout<<"getting the A matrix"<<std::endl;
                                  if (u3::OuterMultiplicity(npp.SU3(), u3::SU3(2,0),np.SU3())>0)
                                    {
                                      boson_matrix(vp,vpp)=
                                        vcs::BosonCreationRME(np,npp)
                                        *ParitySign(u3::ConjugationGrade(omegap)+u3::ConjugationGrade(omegapp))
                                        *u3::UCached(
                                            u_coef_cache,u3::SU3(2,0),npp.SU3(),omegap.SU3(),sigmap.SU3(),
                                            np.SU3(),1,rhop,omegapp.SU3(),rhopp,1);
                                    }
                                  else
                                    boson_matrix(vp,vpp)=0;
                                } //end vp
                            } //end vpp
                          // finished construction boson matrix 
                          // Matrix of symplectic raising operator A
                          Eigen::MatrixXd A=Kp*boson_matrix*Kpp_inv;
                          
                          // Zero initialze blocks accumulating sum over rho0pp and rho0bp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rhobp;
                          ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rhobp);

                          // Summing over rho0bp
                          // std::cout<<"sum over rho0bp"<<std::endl;
                          int rho0bp_max=u3::OuterMultiplicity(omega1.SU3(),x0,omegapp.SU3());
                          for(int rho0bp=1; rho0bp<=rho0bp_max; ++rho0bp)
                            {
                              // Get hypersector index 
                              int hypersector_index3
                                =baby_spncci_hypersectors.LookUpHypersectorIndex(baby_spncci_subspace_indexpp,baby_spncci_subspace_index1,unit_tensor_subspace_index,rho0bp);
                              if(hypersector_index3==-1)
                                continue;

                              double coef3=0;
                              for (int rho0pp=1; rho0pp<=rho0p_max; rho0pp++)
                                { 
                                // omegaxx0->omegapp
                                // x0x(2,0)->x0p (construction)
                                // omegappx(2,0)->omegap (construction)
                                // omegaxx0p->omegap
                                  // std::cout<<"heres a phi"<<std::endl;
                                  coef3+=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0pp)
                                          *u3::UCached(u_coef_cache,
                                            omega1.SU3(),x0,omegap.SU3(),u3::SU3(2,0),
                                            omegapp.SU3(),rho0bp,1,x0p,1,rho0pp
                                            );
                                }// end rho0pp
                              
                              for(int b=0; b<num_blocks; ++b)
                                unit_tensor_blocks_rhobp[b]+=coef3*unit_tensor_hyperblocks[hypersector_index3][b];
                            }
                        
                          // matrix product A*unit_tensor_block (v',v'')*(v'',v1)
                          for(int b=0; b<num_blocks; ++b)
                            for(int i=0; i<gamma_maxp; ++i)
                              for(int j=0; j<gamma_max; ++j)
                                {
                                  // Get target indices 
                                  int it=i*upsilon_maxp;
                                  int jt=j*upsilon_max1;
                                  // Get source indices
                                  int is=i*upsilon_maxpp;
                                  int js=j*upsilon_max1;

                                  unit_tensor_blocks_omegapp[b].block(it,jt,upsilon_maxp,upsilon_max1)
                                    +=A*unit_tensor_blocks_rhobp[b].block(is,js,upsilon_maxpp,upsilon_max1);
                                }
                        } //omegapp
                      // accumulating sum over omegapp
                      for(int b=0; b<num_blocks; ++b)
                        unit_tensor_blocks_x0p[b]+=unit_tensor_blocks_omegapp[b];

                      // std::cout<<"term 1"<<std::endl;
                      // std::cout<<unit_tensor_blocks_x0p[0]<<std::endl;
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
                      //first term 
                      //////////////////////////////////////////////////////////////////////////////////////////////////////////
                      if(u3::OuterMultiplicity(u3::SU3(etap,0),u3::SU3(0,eta-2),x0p)>0)
                        {
                          assert((eta-2)>=0);
                          // look up index of subspace in unit tensor space 
                          u3shell::UnitTensorSubspaceLabels unit_tensor_labels;
                          if(conjugate)
                            unit_tensor_labels=u3shell::UnitTensorSubspaceLabels(u3::Conjugate(x0p),S0,eta-2,etap);
                          else
                            unit_tensor_labels=u3shell::UnitTensorSubspaceLabels(x0p,S0,etap,eta-2);

                          int unit_tensor_subspace_index1=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                          assert(unit_tensor_subspace_index1!=-1);
                          
                          //(rbp,0)x(0,rb)->x0, (0,rb)x(2,0)->(0,rb-2), x0x(2,0)->x0p, (0,rb-2)x(rbp,0)->x0p                        
                          double 
                          coef1=u3::UCached(u_coef_cache,u3::SU3(etap,0),u3::SU3(0,eta),x0p, u3::SU3(2,0),x0,1,1,u3::SU3(0,eta-2),1,1)
                                *sqrt(1.*u3::dim(x0p)*u3::dim(u3::SU3(eta,0))/(u3::dim(x0)));                      
                          
                          double conjugation_factor=1;
                          if(conjugate)
                            conjugation_factor*=sqrt(1.*u3::dim(u3::SU3(etap,0))*u3::dim(omega1)*am::dim(S_ket)/(u3::dim(u3::SU3(eta-2,0))*u3::dim(omegap)*am::dim(S_bra)));

                          // zero initialize blocks for accumulating first term in sum over rhobp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                            ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                          // summing over rho0bp and accumulating sectors in unit1_matrix. 
                          // std::cout<<"rho0p_max "<<rho0p_max<<std::endl;
                          for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp)
                            {
                              int hypersector_index;

                              if(conjugate)
                                hypersector_index=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                    baby_spncci_subspace_index1,baby_spncci_subspace_indexp,
                                    unit_tensor_subspace_index1, rho0bp
                                  );
                              else
                                hypersector_index=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                    baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                    unit_tensor_subspace_index1, rho0bp
                                  );

                              // std::cout<<"hypersector "<<hypersector_index<<std::endl;                        
                              // std::cout<<unit_tensor_space.GetSubspace(unit_tensor_subspace_index1).Str()<<std::endl;
                              // std::cout<<baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp).LabelStr()<<std::endl;
                              // std::cout<<baby_spncci_space.GetSubspace(baby_spncci_subspace_index1).LabelStr()<<std::endl;
                              // Accumulate
                              // if(conjugate)
                              //   std::cout<<"hello"<<hypersector_index<<std::endl; 
                              if(hypersector_index==-1)
                                continue;

                              // std::cout<<"num blocks "<<num_blocks<<std::endl;
                              for(int b=0; b<num_blocks; ++b)
                              {
                                if(conjugate)
                                  {
                                    // Note unit_tensor_subspace_index1 corresponds to conjugate unit tenosr
                                    // So, Tbp->Tb, Sbp->Sb etc. 
                                    auto& unit_tensor_subspace1=unit_tensor_space.GetSubspace(unit_tensor_subspace_index1);
                                    int Tbp,Sbp,Sb,Tb;
                                    std::tie(std::ignore,Sbp,Tbp,Sb,Tb)=unit_tensor_subspace1.GetStateLabels(b);
                                    conjugation_factor*=sqrt(1.*am::dim(Sb)*am::dim(Tb)/(am::dim(Sbp)*am::dim(Tbp)));

                                    unit_tensor_blocks_rho0bp[b]
                                      +=conjugation_grade*conjugation_factor
                                        *u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                        *unit_tensor_hyperblocks[hypersector_index][b].transpose();
                                  }
                                else
                                  unit_tensor_blocks_rho0bp[b]
                                    +=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                      *unit_tensor_hyperblocks[hypersector_index][b];
                                
                                // std::cout<<"phi "<<u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)<<std::endl;
                                // std::cout<<unit_tensor_hyperblocks[hypersector_index][b]<<std::endl;
                              }
                            } //end rho0bp

                          // std::cout<<"coef1 "<<coef1<<std::endl;  
                          // std::cout<<unit_tensor_blocks_rho0bp[0]<<std::endl;
                          for(int b=0; b<num_blocks; ++b)
                            unit_tensor_blocks_x0p[b]+=coef1*unit_tensor_blocks_rho0bp[b];
                        } 

                        // std::cout<<"term 2 "<<std::endl;
                        // std::cout<<unit_tensor_blocks_x0p[0]<<std::endl;
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // second term 
                        //////////////////////////////////////////////////////////////////////////////////////////////////////////  
                        if ((u3::OuterMultiplicity(u3::SU3(etap+2,0),u3::SU3(0,eta),x0p)>0) && (etap+2)<=Nmax+2*N1v)
                          {
                            // (2,0)x(rbp,0)->(rbp+2,0), (rbp,0)x(0,rb)->x0, (rbp+2,0)x(0,rb)->x0p, x0x(2,0)->x0p

                            // look up index of subspace in unit tensor space 
                            u3shell::UnitTensorSubspaceLabels unit_tensor_labels;
                            if(conjugate)
                              unit_tensor_labels=u3shell::UnitTensorSubspaceLabels(u3::Conjugate(x0p),S0,eta,etap+2);
                            else  
                              unit_tensor_labels=u3shell::UnitTensorSubspaceLabels(x0p,S0,etap+2,eta);
                            

                            int unit_tensor_subspace_index2=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                            assert(unit_tensor_subspace_index2!=-1);


                            double 
                            coef2=-1*ParitySign(u3::ConjugationGrade(x0)-u3::ConjugationGrade(x0p))
                                    * u3::dim(u3::SU3(etap,0))*sqrt(u3::dim(x0p)/(6.*u3::dim(u3::SU3(eta,0))))
                                    *u3::UCached(u_coef_cache,u3::SU3(etap+2,0),u3::SU3(0,etap),x0p,x0,
                                            u3::SU3(2,0),1,1,u3::SU3(0,eta),1,1);

                          double conjugation_factor=1;
                          if(conjugate)
                            conjugation_factor*=sqrt(1.*u3::dim(u3::SU3(etap+2,0))*u3::dim(omega1)*am::dim(S_ket)/(u3::dim(u3::SU3(eta,0))*u3::dim(omegap)*am::dim(S_bra)));

                            // zero initialize blocks for accumulating first term in sum over rhobp
                            std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                            for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp)
                              {

                                int hypersector_index2;
                                
                                if(conjugate)
                                  hypersector_index2=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                      baby_spncci_subspace_index1,baby_spncci_subspace_indexp,
                                      unit_tensor_subspace_index2, rho0bp
                                    );
                                else
                                  hypersector_index2=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                      baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                      unit_tensor_subspace_index2, rho0bp
                                    );

                                if(hypersector_index==41)
                                {
                                auto& bra_sub=baby_spncci_space.GetSubspace(baby_spncci_subspace_index1);
                                auto& ket_sub=baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);
                                auto& unit_sub=unit_tensor_space.GetSubspace(unit_tensor_subspace_index2);

                                std::cout<<"hypersector "<<hypersector_index2<<" : "<<bra_sub.LabelStr()
                                <<"  "<<ket_sub.LabelStr()<<" "<<unit_sub.LabelStr()<<std::endl;
                                }

                                // Check that this was actually supposed to be baby_spncci_index1.

                                // int hypersector_index
                                //   =baby_spncci_hypersectors.LookUpHypersectorIndex(
                                //     baby_spncci_subspace_indexp,baby_spncci_index,
                                //     unit_tensor_subspace_index2, rho0bp
                                //   );
                                
                                // if(conjugate)
                                //   std::cout<<"hello there "<<hypersector_index<<std::endl;
                                // Accumulate
                                if(hypersector_index2==-1)
                                continue;

                                for(int b=0; b<num_blocks; ++b)
                                {
                                  if(hypersector_index==41)
                                    {
                                      // Note unit_tensor_subspace_index1 corresponds to conjugate unit tenosr
                                      // So, Tbp->Tb, Sbp->Sb etc. 
                                      auto& unit_tensor_subspace2=unit_tensor_space.GetSubspace(unit_tensor_subspace_index2);
                                      int Tbp,Sbp,Sb,Tb;
                                      std::tie(std::ignore,Sbp,Tbp,Sb,Tb)=unit_tensor_subspace2.GetStateLabels(b);
                                      conjugation_factor*=sqrt(1.*am::dim(Sb)*am::dim(Tb)/(am::dim(Sbp)*am::dim(Tbp)));
                                      std::cout<<"hypersector source "<<hypersector_index2<<std::endl;
                                      std::cout<<fmt::format("unit state labels {} {} {} {}",Sbp,Tbp,Sb,Tb)<<std::endl;
                                      std::cout<<"conjugation_factor "<<conjugation_factor<<"  conjugation_grade "<<conjugation_grade<<std::endl;
                                      std::cout<<unit_tensor_hyperblocks[hypersector_index2][b]<<std::endl;
                                      unit_tensor_blocks_rho0bp[b]
                                        +=conjugation_grade*conjugation_factor
                                          *u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                          *unit_tensor_hyperblocks[hypersector_index2][b].transpose();
                                    }
                                  else
                                    unit_tensor_blocks_rho0bp[b]
                                      +=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                        *unit_tensor_hyperblocks[hypersector_index2][b];
                                  
                                  // std::cout<<"phi "<<u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)<<std::endl;
                                  // std::cout<<unit_tensor_hyperblocks[hypersector_index][b]<<std::endl;
                                }

                                // for(int b=0; b<num_blocks; ++b)
                                //   unit_tensor_blocks_rho0bp[b]
                                //     +=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                //       *unit_tensor_hyperblocks[hypersector_index][b];
                              } //end rho0bp

                            for(int b=0; b<num_blocks; ++b)
                              unit_tensor_blocks_x0p[b]+=coef2*unit_tensor_blocks_rho0bp[b];
                          }

                        for(int b=0; b<num_blocks; ++b)
                          unit_tensor_blocks_omega1[b]+=coef*unit_tensor_blocks_x0p[b];
                      }//end rho0p
                  }// end x0p sum 
                // std::cout<<"summing over n, rho, n1,rho1,upsilon1"<<std::endl;
                // summing over n, rho, n1, rho1, v1
                for(int b=0; b<num_blocks; ++b)
                  for(int i=0; i<gamma_maxp; ++i)
                    for(int j=0; j<gamma_max; ++j)
                    {
                      // std::cout<<fmt::format("{} {} {}  {} {} {}  {} {} {}",
                      //   i,gamma_maxp,dimp,j, gamma_max,dim, upsilon_maxp,upsilon_max, upsilon_max1)
                      // <<std::endl;
                      int it=i*upsilon_maxp;
                      int jt=j*upsilon_max;
                      int is=i*upsilon_maxp;
                      int js=j*upsilon_max1;
                      // (v'v1) (v1 v)
                      // std::cout<<it<<"  "<<jt<<"  "<<is<<"  "<<js<<std::endl;
                      // std::cout<<unit_tensor_blocks[b]<<std::endl;
                      // std::cout<<unit_tensor_blocks[b].block(it,jt,upsilon_maxp,upsilon_max)<<std::endl<<std::endl;
                      // std::cout<<unit_tensor_blocks_omega1[b].block(is,js,upsilon_maxp,upsilon_max1)*KBUK<<std::endl<<std::endl;
                      // std::cout<<unit_tensor_blocks[b].block(it,jt,dimp,dim)+unit_tensor_blocks_omega1[b].block(is,js,dimp,dim1)*KBUK<<std::endl<<std::endl;
                      unit_tensor_blocks[b].block(it,jt,upsilon_maxp,upsilon_max)
                        +=unit_tensor_blocks_omega1[b].block(is,js,upsilon_maxp,upsilon_max1)*KBUK;
                        
                    }
                // std::cout<<"done with blocks"<<std::endl;
              }// end omega1_mult 
        }// end hypersector index 
    }// end Nsum
  // std::cout<<"end recurrence"<<std::endl;
  }
} // End namespace 
  

          

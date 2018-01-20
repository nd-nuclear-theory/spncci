/****************************************************************
  unit_tensor.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/recurrence.h"

#include <omp.h>  

#include "cppformat/format.h"
#include "mcutils/eigen.h"
#include "sp3rlib/u3coef.h"
#include "spncci/spncci_common.h"

extern double zero_threshold;

namespace spncci
{
  void 
  ZeroInitBlocks(int number, int rows, int cols,std::vector<basis::OperatorBlock<double>>& unit_tensor_blocks)
  {
    unit_tensor_blocks.resize(number);
    for(int i=0; i<number; ++i)
      unit_tensor_blocks[i]=Eigen::MatrixXd::Zero(rows,cols);
  }


 int FromIrrepIndexGetBabySpncciIndex(
    int irrep_family_index,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const spncci::BabySpNCCISpace& baby_spncci_space
    )
    {
      // Extract LGI labels 
      const lgi::LGI& lgi=lgi_families[irrep_family_index].irrep;
      u3::U3 sigma;
      HalfInt Sp,Sn,S;
      std::tie(std::ignore,sigma,Sp,Sn,S)=lgi.Key();

      // Get baby spncci index 
      spncci::BabySpNCCISubspaceLabels baby_spncci_labels(sigma,Sp,Sn,S,sigma);
      int baby_spncci_index=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labels);
      return baby_spncci_index;
    }

  void 
  PopulateHypersectorsWithSeeds(
    int irrep_family_index_bra, int irrep_family_index_ket,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors_Nn0,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensors,
    const std::vector<int>& rho0_values,
    basis::MatrixVector& unit_tensor_seed_blocks,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_Nn0,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  )
  {
    int baby_spncci_index_bra
      =spncci::FromIrrepIndexGetBabySpncciIndex(irrep_family_index_bra,lgi_families,baby_spncci_space);
    int baby_spncci_index_ket
      =spncci::FromIrrepIndexGetBabySpncciIndex(irrep_family_index_ket,lgi_families,baby_spncci_space);


    const lgi::LGI& lgi_bra=lgi_families[irrep_family_index_bra].irrep;
    const lgi::LGI& lgi_ket=lgi_families[irrep_family_index_ket].irrep;
    const u3::U3& sigmap=lgi_bra.sigma();
    const u3::U3& sigma=lgi_ket.sigma();

    double conjugation_factor_base
        =ParitySign(u3::ConjugationGrade(sigmap)+lgi_bra.S()-u3::ConjugationGrade(sigma)-lgi_ket.S())
          * sqrt(1.*u3::dim(sigma)*am::dim(lgi_ket.S())/u3::dim(sigmap)/am::dim(lgi_bra.S()));


    for(int i=0; i<lgi_unit_tensors.size();  ++i)
      {
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor=lgi_unit_tensors[i];
        int rho0=rho0_values[i];
        
        // Extract unit tensor labels 
        u3::SU3 x0;
        HalfInt S0,T0, Sp,Tp,S,T;
        int etap, eta; 
        std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=unit_tensor.FlatKey();

        // Look up unit tensor subspace
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
        int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
        auto& subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        // Look up index of unit tensor in subspace
        int unit_tensor_state_index
              =subspace.LookUpStateIndex(std::tuple<int,int,int,int,int>(int(T0),int(Sp),int(Tp),int(S),int(T)));

        // Get Hypersector index
        int hypersector_index
            =baby_spncci_hypersectors.LookUpHypersectorIndex(
                baby_spncci_index_bra,baby_spncci_index_ket,
                unit_tensor_subspace_index,rho0
              );

        unit_tensor_hyperblocks[hypersector_index][unit_tensor_state_index]=unit_tensor_seed_blocks[i];

        // Get conjugate 

        // Look up conjugate unit tensor subspace
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
        int unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
        auto& subspace_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj);

        // Look up index of unit tensor in subspace
        int unit_tensor_state_index_conj
              =subspace_conj.LookUpStateIndex(std::tuple<int,int,int,int,int>(int(T0),int(S),int(T),int(Sp),int(Tp)));

        // Get Hypersector index
        int hypersector_index_Nn0
            =baby_spncci_hypersectors_Nn0.LookUpHypersectorIndex(
                baby_spncci_index_ket,baby_spncci_index_bra,
                unit_tensor_subspace_index_conj,rho0
              );

        double conjugation_factor
        =conjugation_factor_base
          *sqrt(u3::dim(u3::SU3(eta,0))*am::dim(S)*am::dim(T)/u3::dim(u3::SU3(etap,0))/am::dim(Sp)/am::dim(Tp));
        
        unit_tensor_hyperblocks_Nn0[hypersector_index_Nn0][unit_tensor_state_index_conj]
          =conjugation_factor*unit_tensor_seed_blocks[i].transpose();

      }
    }

void AddNn0BlocksToHyperblocks(
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors_Nn0,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_Nn0,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
)
{
  for(int hypersector_index_Nn0=0; hypersector_index_Nn0<baby_spncci_hypersectors_Nn0.size(); ++hypersector_index_Nn0)
    {
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Extracting labels from source (Nn0 sectors)
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Get hypersector indices for Nn0 sector
      // Note the labels from Nn0 sectors are conjugate labels, i.e., bra is actually ket and vice versa
      auto key=baby_spncci_hypersectors_Nn0.GetHypersector(hypersector_index_Nn0).Key();
      int unit_tensor_subspace_index_Nn0, baby_spncci_index_bra, baby_spncci_index_ket, rho0;
      std::tie(baby_spncci_index_ket,baby_spncci_index_bra,unit_tensor_subspace_index_Nn0,rho0)=key;

      // Extract unit tensor subspace labels from Nn0 tensor
      auto& unit_tensor_subspace_Nn0=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_Nn0);      
      u3::SU3 x0c;
      HalfInt S0;
      int etap,eta;
      std::tie(x0c,S0,eta,etap)=unit_tensor_subspace_Nn0.labels();
      
      // Get bra and ket labels from Nn0 sector
      const spncci::BabySpNCCISubspace& subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra);
      const spncci::BabySpNCCISubspace& subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket);
      
      u3::U3 omegap=subspace_bra.omega();
      HalfInt Sp=subspace_bra.S();

      u3::U3 omega=subspace_ket.omega();
      HalfInt S=subspace_ket.S();
      
      // part of conjugation factor
      double conjugation_factor_base
              =ParitySign(u3::ConjugationGrade(omega)+S-u3::ConjugationGrade(omegap)-Sp)
                *sqrt(1.*u3::dim(omega)*u3::dim(u3::SU3(etap,0))*am::dim(S)
                       /u3::dim(omegap)/u3::dim(u3::SU3(eta,0))/am::dim(Sp)
                  );

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Looking up target hypersector
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Get unit tensor subspace index in hyperblocks 
      u3shell::UnitTensorSubspaceLabels unit_tensor_labels
        =u3shell::UnitTensorSubspaceLabels(u3::Conjugate(x0c),S0,etap,eta);

      int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
      auto& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);      

      // Look up hypersector
      int hypersector_index
          =baby_spncci_hypersectors.LookUpHypersectorIndex(
              baby_spncci_index_bra,baby_spncci_index_ket,unit_tensor_subspace_index,rho0
            );

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // For each source hyperblock, identify target block and conjugate 
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for(int unit_tensor_index_Nn0=0; unit_tensor_index_Nn0<unit_tensor_subspace_Nn0.size(); ++unit_tensor_index_Nn0)
        {
          int Tbp,Sbp,Sb,Tb,T0;
          std::tie(T0,Sb,Tb,Sbp,Tbp)=unit_tensor_subspace_Nn0.GetStateLabels(unit_tensor_index_Nn0);

          // Get unit tensor index
          std::tuple<int,int,int,int,int> state_labels(T0,Sbp,Tbp,Sb,Tb);
          int unit_tensor_index=unit_tensor_subspace.LookUpStateIndex(state_labels);

          // Conjugation factor
          double conjugation_factor
                  =conjugation_factor_base*sqrt(am::dim(Sbp)*am::dim(Tbp)/am::dim(Sb)/am::dim(Tb));

          unit_tensor_hyperblocks[hypersector_index][unit_tensor_index]
            =conjugation_factor*unit_tensor_hyperblocks_Nn0[hypersector_index_Nn0][unit_tensor_index_Nn0].transpose();
        }
    }
}



  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Construct KBUK matrix
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
  void ConsructKBUK(
    u3::UCoefCache& u_coef_cache, int Nn,
    const u3::U3& sigma, const u3::U3& omega, const u3::U3& omega1,
    const sp3r::U3Subspace& u3_subspace,
    const sp3r::U3Subspace& u3_subspace1,
    const Eigen::MatrixXd& K1, //(upsilon_max1,dim1)
    const Eigen::MatrixXd& K_inv,//(dim,upsilon_max)
    Eigen::MatrixXd& KBUK //(upislon1_max,upsilon)
    )
  {
    int upsilon_max1=u3_subspace1.upsilon_max();
    int upsilon_max=u3_subspace.upsilon_max();
    int dim1=u3_subspace1.size();
    int dim=u3_subspace.size();

    Eigen::MatrixXd BU(dim1, dim);
    for(int u3_state_index=0; u3_state_index<dim; ++u3_state_index)
      {
        MultiplicityTagged<u3::U3> n_rho=u3_subspace.GetStateLabels(u3_state_index);
        u3::U3 n(n_rho.irrep);
        // iterate over (n1,rho1)
        for (int u3_state_index1=0; u3_state_index1<dim1; u3_state_index1++)
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
    // Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
    KBUK.noalias()=K1*BU*K_inv;
    // std::cout<<"KBUK "<<KBUK<<std::endl;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////

void Amatrix(
  u3::UCoefCache& u_coef_cache,
  const sp3r::U3Subspace& u3_subspacep,
  const sp3r::U3Subspace& u3_subspacepp,
  const u3::U3& sigmap, const u3::U3& omegap, const u3::U3& omegapp,
  const vcs::MatrixCache& K_matrix_map_bra,
  const vcs::MatrixCache& Kinv_matrix_map_bra,
  int upsilon_maxp, int upsilon_maxpp,
  Eigen::MatrixXd& A
  )
{
  // underlying u3boson subspace dimensions
  int dimp=u3_subspacep.size();
  int dimpp=u3_subspacepp.size();

  Eigen::MatrixXd boson_matrix(dimp,dimpp);
  
  // Extracting K matrices 
  const Eigen::MatrixXd& Kp=K_matrix_map_bra.at(omegap);
  Eigen::MatrixXd Kpp_inv=Kinv_matrix_map_bra.at(omegapp);

  for(int vpp=0; vpp<dimpp; vpp++)
    {
      MultiplicityTagged<u3::U3> npp_rhopp=u3_subspacepp.GetStateLabels(vpp);
      const u3::U3& npp(npp_rhopp.irrep);
      int rhopp=npp_rhopp.tag;
      for(int vp=0; vp<dimp; vp++)
        {
          MultiplicityTagged<u3::U3> np_rhop=u3_subspacep.GetStateLabels(vp);
          const u3::U3& np(np_rhop.irrep);
          int rhop=np_rhop.tag; 
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
  // Matrix of symplectic raising operator A
  A=Kp*boson_matrix*Kpp_inv;
}


void 
ComputeUnitTensorHyperblocks(
  int Nmax, int N1v,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::KMatrixCache& k_matrix_map,
  const spncci::KMatrixCache& kinv_matrix_map,
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
      const std::vector<int>& unit_tensor_hypersectors=unit_tensor_hypersector_subsets[Nsum/2];

      // Parallelize here 
      // unit_tensor_hyperblocks are zero initalized so there should be no race conditions in 
      // writing each hyperbock to unit_tensor_hyperblocks
      // std::cout<<"omp_get_num_threads "<<omp_get_num_threads()<<std::endl;
      
      for(int i=0; i<unit_tensor_hypersectors.size(); ++i)  
      // for(int hypersector_index : unit_tensor_hypersectors)
        {
          int hypersector_index=unit_tensor_hypersectors[i];
          auto key=baby_spncci_hypersectors.GetHypersector(hypersector_index).Key();

          int unit_tensor_subspace_index, baby_spncci_subspace_indexp, baby_spncci_index, rho0;
          std::tie(baby_spncci_subspace_indexp,baby_spncci_index,unit_tensor_subspace_index,rho0)=key;
          
          const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
              =baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);

          const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
              =baby_spncci_space.GetSubspace(baby_spncci_index);
          
          const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace
              =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

          // Extract ket dimensions and mulitplicities 
          int dim=baby_spncci_subspace_ket.size();
          int gamma_max=baby_spncci_subspace_ket.gamma_max();
          int upsilon_max=baby_spncci_subspace_ket.upsilon_max();

          // Extract bra dimensions and mulitplicities 
          int dimp=baby_spncci_subspace_bra.size();
          int gamma_maxp=baby_spncci_subspace_bra.gamma_max();
          int upsilon_maxp=baby_spncci_subspace_bra.upsilon_max();
          
          // Extract Sp(3,R) space
          int irrep_family_index_bra=baby_spncci_subspace_bra.irrep_family_index();
          int irrep_family_index_ket=baby_spncci_subspace_ket.irrep_family_index();

          // std::cout<<"irrep_family_index_bra "<<irrep_family_index_bra<<" irrep_family_index_ket"<<irrep_family_index_ket<<std::endl;
          const sp3r::Sp3RSpace& irrep_bra = spncci_space[irrep_family_index_bra].Sp3RSpace();
          const sp3r::Sp3RSpace& irrep_ket = spncci_space[irrep_family_index_ket].Sp3RSpace();

          // extract subspace labels 
          u3::U3 omegap,sigmap,omega,sigma;
          u3::SU3 x0;
          HalfInt S0,Sn_ket,Sp_ket,S_ket,Sn_bra,Sp_bra,S_bra;
          int etap,eta;

          // Extracting labels
          // std::cout<<"extracting labels"<<std::endl;
          std::tie(sigmap,Sp_bra,Sn_bra,S_bra,omegap)=baby_spncci_subspace_bra.labels();
          std::tie(sigma,Sp_ket,Sn_ket,S_ket,omega)=baby_spncci_subspace_ket.labels();
          std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
          int Nn=baby_spncci_subspace_ket.Nn();

          // omega u3 subspace in irrep
          const sp3r::U3Subspace& u3_subspace=irrep_ket.LookUpSubspace(omega);
          const sp3r::U3Subspace& u3_subspacep=irrep_bra.LookUpSubspace(omegap);
          // Extracting K matrices for sp_irrep and sp_irrepp from the K_matrix_maps 
          // std::cout<<"bunny1"<<std::endl;
          const vcs::MatrixCache& K_matrix_map_bra=k_matrix_map.at(sigmap);
          // std::cout<<"bunny2"<<std::endl;
          // Temporary fix for debug purposes
          vcs::MatrixCache null_cache;
          const vcs::MatrixCache& Kinv_matrix_map_bra=kinv_matrix_map.at(sigmap);
          // std::cout<<"bunny3"<<std::endl;
          const vcs::MatrixCache& K_matrix_map_ket=k_matrix_map.at(sigma);
          // std::cout<<"bunny4"<<std::endl;
          const vcs::MatrixCache& Kinv_matrix_map_ket=kinv_matrix_map.at(sigma);
          // std::cout<<"bunny5"<<std::endl;
          const Eigen::MatrixXd& Kp=K_matrix_map_bra.at(omegap);
          // std::cout<<"bunny6"<<std::endl;
          // std::cout<<sigma.Str()<<". "<<omega.Str()<<std::endl;
          const Eigen::MatrixXd& K_inv=Kinv_matrix_map_ket.at(omega);

          // Generate labels to sum over 
          int rho0_max=u3::OuterMultiplicity(omega.SU3(),x0,omegap.SU3());

          // Precalculating kronecker products used in sum to calculate unit tensor matrix
          MultiplicityTagged<u3::U3>::vector omegapp_set=KroneckerProduct(omegap, u3::U3(0,0,-2)); 
          MultiplicityTagged<u3::U3>::vector omega1_set=KroneckerProduct(omega, u3::U3(0,0,-2));
          MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(x0, u3::SU3(2,0));

          // std::cout<<"hypersector_index "<<hypersector_index<<std::endl;
           std::vector<basis::OperatorBlock<double>>& unit_tensor_blocks=unit_tensor_hyperblocks[hypersector_index];
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //  Calculate unit tensor matrix
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          int num_blocks=unit_tensor_blocks.size();

          for(auto& omega1_mult :omega1_set)
            {
              u3::U3 omega1(omega1_mult.irrep);

              if (not irrep_ket.ContainsSubspace(omega1))
                continue;
                  
              spncci::BabySpNCCISubspaceLabels baby_spncci_labels1(sigma,Sp_ket,Sn_ket,S_ket,omega1);
              int baby_spncci_subspace_index1=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labels1);
              // std::cout<<"bunny rabbit 1"<<std::endl;
              Eigen::MatrixXd K1=K_matrix_map_ket.at(omega1);
              const sp3r::U3Subspace& u3_subspace1=irrep_ket.LookUpSubspace(omega1);
              int upsilon_max1=u3_subspace1.upsilon_max();

              int dim1=upsilon_max1*gamma_max;
              std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omega1;

              // Initializing blocks for sum over omega1
              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omega1);
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              // Construct KBUK matrix
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              Eigen::MatrixXd KBUK(upsilon_max1,upsilon_max);
              
              spncci::ConsructKBUK(
                u_coef_cache, Nn,sigma, omega, omega1,
                u3_subspace,u3_subspace1,K1,K_inv,
                KBUK
              );
              
              // ////////////////////////////////////////////////////////////////////////////////////////////////////////
              //summing over x0'
              // std::cout<<"sum over x0"<<std::endl;
              for (auto& x0p_mult : x0p_set)
                {
                  u3::SU3 x0p(x0p_mult.irrep);
                  // std::cout<<"x0p "<<x0p.Str()<<std::endl;
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

                      // std::cout<<"computing term 3"<<std::endl;
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // third term
                      // sum over omega'', v'' and rho0''
                      ////////////////////////////////////////////////////////////////////////////////////////////////////////
                      // Zero initialze blocks accumulating sum over omegapp
                      std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_omegapp;
                      ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_omegapp);

                      // Summing over omega''
                      // std::cout<<"summing over omegapp"<<std::endl;
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
                          int upsilon_maxpp=u3_subspacepp.upsilon_max();
                          int dimpp=upsilon_maxpp*gamma_maxp;
                          // Obtaining K matrix for omega''
                          
                          Eigen::MatrixXd A;
                          spncci::Amatrix(u_coef_cache,
                            u3_subspacep,u3_subspacepp,sigmap, omegap, omegapp,
                            K_matrix_map_bra,Kinv_matrix_map_bra,upsilon_maxp, upsilon_maxpp,
                            A
                          );

                           // Zero initialze blocks accumulating sum over rho0pp and rho0bp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rhobp;
                          
                          ZeroInitBlocks(num_blocks,dimpp,dim1,unit_tensor_blocks_rhobp);

                          // Summing over rho0bp
                          int rho0bp_max=u3::OuterMultiplicity(omega1.SU3(),x0,omegapp.SU3());
                          for(int rho0bp=1; rho0bp<=rho0bp_max; ++rho0bp)
                            {
                              // Get hypersector index 
                              int hypersector_index3
                                =baby_spncci_hypersectors.LookUpHypersectorIndex(baby_spncci_subspace_indexpp,baby_spncci_subspace_index1,unit_tensor_subspace_index,rho0bp);
                              // std::cout<<"hypersector3 "<<hypersector_index3<<std::endl;
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
                              // std::cout<<"sum over blocks for term 3"<<std::endl;

                              for(int b=0; b<num_blocks; ++b)
                                unit_tensor_blocks_rhobp[b]+=coef3*unit_tensor_hyperblocks[hypersector_index3][b];
                            }
                        
                          // matrix product A*unit_tensor_block (v',v'')*(v'',v1)
                          // std::cout<<"add in A operator and sum over blocks "<<std::endl;
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
                      // std::cout<<"accumulating over unit tensor blocks x0p"<<std::endl; 
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
                          u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap,eta-2);

                          int unit_tensor_subspace_index1=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                          assert(unit_tensor_subspace_index1!=-1);
                          
                          //(rbp,0)x(0,rb)->x0, (0,rb)x(2,0)->(0,rb-2), x0x(2,0)->x0p, (0,rb-2)x(rbp,0)->x0p                        
                          double 
                          coef1=u3::UCached(u_coef_cache,u3::SU3(etap,0),u3::SU3(0,eta),x0p, u3::SU3(2,0),x0,1,1,u3::SU3(0,eta-2),1,1)
                                *sqrt(1.*u3::dim(x0p)*u3::dim(u3::SU3(eta,0))/(u3::dim(x0)));                      
 
                          // zero initialize blocks for accumulating first term in sum over rhobp
                          std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                            ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                          // summing over rho0bp and accumulating sectors in unit1_matrix. 
                          for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp)
                            {
                              int hypersector_index1=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                    baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                    unit_tensor_subspace_index1, rho0bp
                                  );

                              // Accumulate
                              if(hypersector_index1==-1)
                                continue;

                              // std::cout<<"num blocks "<<num_blocks<<std::endl;
                              for(int b=0; b<num_blocks; ++b)
                              {
                                unit_tensor_blocks_rho0bp[b]
                                  +=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                    *unit_tensor_hyperblocks[hypersector_index1][b];
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
                        // std::cout<<"term 2 "<<std::endl;
                        // std::cout<<"x0p "<<x0p.Str()<<std::endl;
                        if ((u3::OuterMultiplicity(u3::SU3(etap+2,0),u3::SU3(0,eta),x0p)>0) && (etap+2)<=Nmax+2*N1v)
                          {
                            // (2,0)x(rbp,0)->(rbp+2,0), (rbp,0)x(0,rb)->x0, (rbp+2,0)x(0,rb)->x0p, x0x(2,0)->x0p
                            // look up index of subspace in unit tensor space 
                            u3shell::UnitTensorSubspaceLabels unit_tensor_labels(x0p,S0,etap+2,eta);
                            

                            int unit_tensor_subspace_index2=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
                            assert(unit_tensor_subspace_index2!=-1);


                            double 
                            coef2=-1*ParitySign(u3::ConjugationGrade(x0)-u3::ConjugationGrade(x0p))
                                    * u3::dim(u3::SU3(etap,0))*sqrt(u3::dim(x0p)/(6.*u3::dim(u3::SU3(eta,0))))
                                    *u3::UCached(u_coef_cache,u3::SU3(etap+2,0),u3::SU3(0,etap),x0p,x0,
                                            u3::SU3(2,0),1,1,u3::SU3(0,eta),1,1);

                            // zero initialize blocks for accumulating first term in sum over rhobp
                            std::vector<basis::OperatorBlock<double>> unit_tensor_blocks_rho0bp;
                              ZeroInitBlocks(num_blocks,dimp,dim1,unit_tensor_blocks_rho0bp);

                            for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp)
                              {

                                int hypersector_index2=baby_spncci_hypersectors.LookUpHypersectorIndex(
                                      baby_spncci_subspace_indexp,baby_spncci_subspace_index1,
                                      unit_tensor_subspace_index2, rho0bp
                                    );

                                auto& bra_sub=baby_spncci_space.GetSubspace(baby_spncci_subspace_index1);
                                auto& ket_sub=baby_spncci_space.GetSubspace(baby_spncci_subspace_indexp);
                                auto& unit_sub=unit_tensor_space.GetSubspace(unit_tensor_subspace_index2);

                                // Accumulate
                                if(hypersector_index2==-1)
                                continue;

                                for(int b=0; b<num_blocks; ++b)
                                {
                                  unit_tensor_blocks_rho0bp[b]
                                    +=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                      *unit_tensor_hyperblocks[hypersector_index2][b];
                                }
                              } //end rho0bp

                            for(int b=0; b<num_blocks; ++b)
                              unit_tensor_blocks_x0p[b]+=coef2*unit_tensor_blocks_rho0bp[b];
                          }

                        for(int b=0; b<num_blocks; ++b)
                          unit_tensor_blocks_omega1[b]+=coef*unit_tensor_blocks_x0p[b];
                      }//end rho0p
                  }// end x0p sum 
                // summing over n, rho, n1, rho1, v1
                for(int b=0; b<num_blocks; ++b)
                  for(int i=0; i<gamma_maxp; ++i)
                    for(int j=0; j<gamma_max; ++j)
                    {
                      int it=i*upsilon_maxp;
                      int jt=j*upsilon_max;
                      int is=i*upsilon_maxp;
                      int js=j*upsilon_max1;
                      // (v'v1) (v1 v)
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
  

          

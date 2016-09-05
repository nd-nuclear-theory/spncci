/****************************************************************
  relative_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "u3shell/relative_operator.h"

#include "am/am.h"
#include "cppformat/format.h"
#include "sp3rlib/vcs.h"
#include "u3shell/u3st_scheme.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/moshinsky.h"
namespace u3shell {


  void GenerateRelativeUnitTensorLabelsU3ST(
        const u3shell::RelativeSpaceU3ST& space,
        std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels
        )
  {
    for(int ket_index=0; ket_index<space.size(); ++ket_index)
      {
        int eta, S, T;
        // const u3shell::RelativeSubspaceU3ST& subspace =space.GetSubspace(ket_index);
        std::tie(eta,S,T,std::ignore)=space.GetSubspace(ket_index).Key();
        for(int bra_index=0; bra_index<space.size(); ++bra_index)
          {
            int etap,Sp,Tp;
            std::tie(etap,Sp,Tp,std::ignore)=space.GetSubspace(bra_index).Key();
            MultiplicityTagged<u3::SU3>::vector x0_list=KroneckerProduct(u3::SU3(etap,0),u3::SU3(0,eta));
            for(int S0=abs(Sp-S); S0<=(Sp+S); ++S0)
              for(int T0=abs(Tp-T); T0<=(Tp+T); ++T0)
                for(int i=0; i<x0_list.size(); ++i)
                  {
                    u3::SU3 x0=x0_list[i].irrep;
                    u3shell::RelativeStateLabelsU3ST ket(eta,S,T);
                    u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp);
                    relative_unit_tensor_labels[etap-eta].push_back(u3shell::RelativeUnitTensorLabelsU3ST(x0,S0,T0,bra,ket));
                  }
          }
      }
  }

  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax, 
        std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels
        )
  {   
    #ifdef VERBOSE
    std::cout<<"Entering GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
    #endif
    for(int N0=0; N0<=Nmax; N0+=2)
      {
        std::vector<RelativeUnitTensorLabelsU3ST> sym_vec;
        for(int Sp=0; Sp<=1; Sp++)
          for(int Tp=0; Tp<=1; Tp++)
            for(int S=0; S<=1; S++)
              for (int T=0; T<=1; T++)
                for (int S0=abs(S-Sp); S0<=(S+Sp); S0++)
                  for (int T0=abs(T-Tp); T0<=(T+Tp); T0++)
                    for(int etap=0; etap<=N0+Nmax; etap++)
                      {
                        //antisymmeterization constraint on ket 
                        if ( (etap+Sp+Tp)%2!=1 )
                          continue;
                        
                        int eta=etap-N0;
                        //antisymmeterization constraint on bra 
                        if ( (eta+S+T)%2!=1)
                          continue;

                        u3shell::RelativeStateLabelsU3ST ket(eta,S,T);
                        u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp);

                        MultiplicityTagged<u3::SU3>::vector omega0_set
                          =u3::KroneckerProduct(u3::SU3(etap,0),u3::SU3(0,eta));

                        for(int w=0; w<omega0_set.size(); w++)
                          {
                            u3::SU3 x0(omega0_set[w].irrep);
                            sym_vec.push_back(u3shell::RelativeUnitTensorLabelsU3ST(x0,S0,T0,bra,ket));
                            //std::cout<<"unit tensors  "<<spncci::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T).Str()<<std::endl;
                          }
                      }       
        relative_unit_tensor_labels[N0]=sym_vec;
      }
  #ifdef VERBOSE
  std::cout<<"Exiting GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
  #endif
  } //end function

  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax, 
        std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
        )
  {   
    #ifdef VERBOSE
    std::cout<<"Entering GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
    #endif
    int T0_max;
    for(int N0=0; N0<=Nmax; N0+=2)
      {
        for(int Sp=0; Sp<=1; Sp++)
          for(int Tp=0; Tp<=1; Tp++)
            for(int S=0; S<=1; S++)
              for (int T=0; T<=1; T++)
                for (int S0=abs(S-Sp); S0<=(S+Sp); S0++)
                  for (int T0=abs(T-Tp); T0<=(T+Tp); T0++)
                    for(int etap=0; etap<=N0+Nmax; etap++)
                      {
                        if(T0!=0)
                          continue;
                        //antisymmeterization constraint on ket 
                        if ( (etap+Sp+Tp)%2!=1 )
                          continue;
                        
                        int eta=etap-N0;
                        //antisymmeterization constraint on bra 
                        if ( (eta+S+T)%2!=1)
                          continue;

                        u3shell::RelativeStateLabelsU3ST ket(eta,S,T);
                        u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp);

                        MultiplicityTagged<u3::SU3>::vector omega0_set
                          =u3::KroneckerProduct(u3::SU3(etap,0),u3::SU3(0,eta));

                        for(int w=0; w<omega0_set.size(); w++)
                          {
                            u3::SU3 x0(omega0_set[w].irrep);
                            //Restriction to J0=0 
                            int L0_min;
                            if(std::min(x0.lambda(),x0.mu())%2==1)
                              L0_min=1;
                            else
                              L0_min=std::max(x0.lambda(),x0.mu())%2;
                            
                            int L0_max=x0.lambda()+x0.mu();
                            if((L0_max<S0)||(L0_min)>S0)
                              continue;
                            relative_unit_tensor_labels.push_back(u3shell::RelativeUnitTensorLabelsU3ST(x0,S0,T0,bra,ket));
                            //std::cout<<"unit tensors  "<<spncci::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T).Str()<<std::endl;
                          }
                      }       
      }
  #ifdef VERBOSE
  std::cout<<"Exiting GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
  #endif
  } //end function 

  double RelativeNumberOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    if (bra==ket)
      // 1.5 from the 3/2 zero point energy for a single particle
      rme=ket.eta()+1.5;
    return rme;
  }

  double RelativeSp3rRaisingOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    int eta=ket.eta();
    int etap=bra.eta();
    if(
        (bra.S()==ket.S())    // delta on spin
        && (bra.T()==bra.T()) // delta on isospin
        && ((etap-eta)==2)      //only connect states with eta+2=etap 
      )
        {
          u3::U3 sigma(HalfInt(eta%2)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omega(HalfInt(eta)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omegap(HalfInt(etap)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          MultiplicityTagged<u3::U3>n_rho(u3::U3(omega.N()-sigma.N(),0,0),1);
          MultiplicityTagged<u3::U3>np_rhop(u3::U3(omegap.N()-sigma.N(),0,0),1);
          rme=(sqrt(vcs::Omega(np_rhop.irrep, omegap)-vcs::Omega(n_rho.irrep, omega))
                    *vcs::U3BosonCreationRME(sigma,np_rhop,omegap,sigma,n_rho,omega));
        }
    return rme;
  }

  double RelativeSp3rLoweringOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {

    double rme=0;
    int eta=ket.eta();
    int etap=bra.eta();
    if((eta==0)||(eta==1))
      return rme;
    
    if(
        (bra.S()==ket.S())    // delta on spin
        && (bra.T()==bra.T()) // delta on isospin
        && ((eta-etap)==2)      //only connect states with eta+2=etap 
      )
        {
          u3::U3 sigma(HalfInt(eta%2)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omega(HalfInt(eta)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omegap(HalfInt(etap)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          MultiplicityTagged<u3::U3>n_rho(u3::U3(omega.N()-sigma.N(),0,0),1);
          MultiplicityTagged<u3::U3>np_rhop(u3::U3(omegap.N()-sigma.N(),0,0),1);
          rme=parity(etap-eta)
                    *sqrt(1.*u3::dim(omega)/u3::dim(omegap))
                    *sqrt(vcs::Omega(n_rho.irrep, omega)-vcs::Omega(np_rhop.irrep, omegap))
                    *vcs::U3BosonCreationRME(sigma,n_rho,omega,sigma,np_rhop,omegap);
        }
    return rme;
  }

  double RelativeKineticEnergyOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    if (bra.eta()==ket.eta())
      rme=u3shell::RelativeNumberOperator(bra,ket);
    if (bra.eta()==(ket.eta()+2))
      rme=-sqrt(1.5)*u3shell::RelativeSp3rRaisingOperator(bra,ket);
    if (bra.eta()==(ket.eta()-2))
      rme=-sqrt(1.5)*u3shell::RelativeSp3rLoweringOperator(bra,ket);

    return rme;
  }

} // namespace

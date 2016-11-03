/****************************************************************
  relative_to_twobody.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
****************************************************************/
#include <fstream>

#include "cppformat/format.h"
#include "basis/lsjt_operator.h"

#include "am/am.h"
#include "am/wigner_gsl.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/import_interaction.h"
#include "u3shell/relative_operator.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/moshinsky.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/upcoupling.h"


void H2FormatLookUp(int Nmax,std::map<std::tuple<int,int,HalfInt>,int>& h2_lookup)
  {
    int index=1;
    for(int N=0; N<=Nmax; ++N)
      for(int L=(N%2); L<=N; L+=2)
        for(HalfInt J=abs(L-HalfInt(1,2)); J<=(L+HalfInt(1,2)); ++J)
          {
            std::tuple<int,int,HalfInt> key(N,L,J);
            h2_lookup[key]=index;
            std::cout<<fmt::format("{} {} {}  {}",N, L, J, index)<<std::endl;
            ++index;
          }
  }

typedef std::tuple<int,int,HalfInt,HalfInt> RelativeStateLabelsNLST;
typedef std::tuple<int, HalfInt,HalfInt,RelativeStateLabelsNLST,RelativeStateLabelsNLST>RelativeBraketNLST;

void
BranchNLST(
  const u3shell::RelativeUnitTensorLabelsU3ST& relative_unit_tensor,
  std::map<RelativeBraketNLST,double>& branched_rel_unit_tensors
  )
{
  int N0, etap,eta,eta_cm;
  HalfInt S0,T0,Sp,Tp,S,T;
  u3::SU3 x0;
  std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=relative_unit_tensor.FlatKey();
  MultiplicityTagged<int>::vector L0_set=u3::BranchingSO3(x0);
  for(auto Lk0 : L0_set)
    {
      int L0=Lk0.irrep;
      int kappa0_max=Lk0.tag;
      for(int kappa0=1; kappa0<=kappa0_max; ++kappa0)
        for(int Lp=etap%2; Lp<=etap; Lp+=2)
          for(int L=eta%2; L<=eta; L+=2)
            {
              if(not am::AllowedTriangle(L,L0,Lp))
                continue;
              RelativeStateLabelsNLST bra_nlst(etap,Lp,Sp,Tp);
              RelativeStateLabelsNLST ket_nlst(eta,L,S,T);
              RelativeBraketNLST braket_nlst(L0,S0,T0,bra_nlst,ket_nlst);
              branched_rel_unit_tensors[braket_nlst]
                +=u3::W(u3::SU3(eta,0),1,L,x0,kappa0,L0,u3::SU3(etap,0),1,Lp,1);
            }
    }
}

// Nr,Lr,Ncm,Lcm,L,S,T
typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> RelativeCMStateLabelsNLST;
typedef std::tuple<int,HalfInt,HalfInt,RelativeCMStateLabelsNLST,RelativeCMStateLabelsNLST> RelativeCMBraketNLST;

void RelativeToCMLST(
  int Nmax, 
  const std::map<RelativeBraketNLST,double>& branched_rel_unit_tensors,
  std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map)
// Coupling on center of mass as SO(3) level
{
  for(auto it=branched_rel_unit_tensors.begin(); it!=branched_rel_unit_tensors.end(); ++it)
    {
      int L0,Lr,Lrp,N,Np;
      HalfInt S,T,Sp,Tp,S0,T0;
      RelativeStateLabelsNLST ket_rel,bra_rel;
      RelativeBraketNLST rel_tensor=it->first;
      double rme_rel=it->second;
      std::tie(L0,S0,T0,bra_rel,ket_rel)=rel_tensor;
      std::tie(Np,Lrp,Sp,Tp)=bra_rel;
      std::tie(N,Lr,S,T)=ket_rel;
      int Ncm_max=Nmax;//-std::max(Np,N);
      for(int Ncm=0; Ncm<=Ncm_max; ++Ncm)
        for(int Lcm=Ncm%2; Lcm<=Ncm; ++Lcm)
          for(int Lp=abs(Lrp-Lcm); Lp<=(Lrp+Lcm); ++Lp)
            for(int L=abs(Lr-Lcm); L<=(Lr+Lcm); ++L)
                {
                  RelativeCMStateLabelsNLST bra_rel_cm(Np,Lrp,Ncm,Lcm,Lp,Sp,Tp);
                  RelativeCMStateLabelsNLST ket_rel_cm(N,Lr,Ncm,Lcm,L,S,T);
                  RelativeCMBraketNLST braket_rel_cm(L0,S0,T0,bra_rel_cm,ket_rel_cm);
                  rel_cm_lst_map[braket_rel_cm]
                    // =am::Wigner6J(Lrp,Lp,Lcm,L,Lr,L0)
                    // *parity(Lrp+Lcm+L+L0)*Hat(L)*Hat(Lrp);
                    =am::Unitary9J(L0,Lr,Lrp,0,Lcm,Lcm,L0,L,Lp)*rme_rel;
                    std::cout<<
                    fmt::format("{} {}   {} {}   {} {}   {:3f}   {:3f}",
                    L0,Lcm, Lrp,Lp,Lr,L,
                    am::Unitary9J(L0,0,L0,Lr,Lcm,L,Lrp,Lcm,Lp),rme_rel)<<std::endl;
                }
    }
}

typedef std::map<u3shell::RelativeCMUnitTensorLabelsU3ST,double> RelativeCMUnitTensorExpansion;

void UpcoupleCMU3ST(
  std::map<RelativeCMBraketNLST,double> rel_cm_lst_map,
  RelativeCMUnitTensorExpansion& rel_cm_u3st_map
  )
{
  int Nrp,Lrp,Ncm,Lcm,Nr,Lr,Lp,L,L0;
  HalfInt Sp,Tp,S,T,S0,T0;
  RelativeCMStateLabelsNLST bra_cm,ket_cm;
  // Sum over Lr, Lcm
  // Sum over Lrp, Lcm
  // Sum over L, Lp with dim factors taken into account 
  for(auto it=rel_cm_lst_map.begin(); it!=rel_cm_lst_map.end(); ++it)
    {
      std::tie(L0,S0,T0,bra_cm,ket_cm)=it->first;
      double rme=it->second;
      std::tie(Nrp,Lrp,Ncm,Lcm,Lp,Sp,Tp)=bra_cm;
      std::tie(Nr,Lr,Ncm,Lcm,L,S,T)=ket_cm;

      MultiplicityTagged<u3::SU3>::vector  x_set=u3::KroneckerProduct(u3::SU3(Nr,0),u3::SU3(Ncm,0));
      MultiplicityTagged<u3::SU3>::vector  xp_set=u3::KroneckerProduct(u3::SU3(Nrp,0),u3::SU3(Ncm,0));
      for(int i=0; i<x_set.size(); ++i)
        {
          u3::SU3 x=x_set[i].irrep;
          int kappa_max=u3::BranchingMultiplicitySO3(x,L);
          for(int ip=0; ip<xp_set.size(); ++ip)
            {
              u3::SU3 xp=xp_set[ip].irrep;
              MultiplicityTagged<u3::SU3>::vector x0_set=u3::KroneckerProduct(xp,u3::Conjugate(x));
              int kappap_max=u3::BranchingMultiplicitySO3(xp,Lp); 
              for(int kappa=1; kappa<=kappa_max; ++kappa)
                for(int kappap=1; kappap<=kappap_max; ++kappap)
                  {
                    for(int i0=0; i0<x0_set.size();  i0++)
                      {
                        u3::SU3 x0=x0_set[i0].irrep;
                        int rho0_max=u3::OuterMultiplicity(x,x0,xp);
                        int kappa0_max=u3::BranchingMultiplicitySO3(x0,L0);
                        u3shell::RelativeCMStateLabelsU3ST bra(Nrp,Ncm,xp,Sp,Tp);
                        u3shell::RelativeCMStateLabelsU3ST ket(Nr,Ncm,x,S,T);
                        // if(rel_cm_u3st_map.count(u3shell::RelativeCMUnitTensorLabelsU3ST(x0,S0,T0,1,bra,ket))!=0)
                        //   continue;
                        // if (L0!=std::max(x0.lambda(),x0.mu()))
                        if(L0!=0)
                          continue;
                        for(int kappa0=1; kappa0<=kappa0_max; ++kappa0)
                          for(int rho0=1; rho0<=rho0_max; ++rho0)
                              {
                                // u3shell::RelativeCMStateLabelsU3ST bra(Nrp,Ncm,xp,Sp,Tp);
                                // u3shell::RelativeCMStateLabelsU3ST ket(Nr,Ncm,x,S,T);
                                u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st(x0,S0,T0,rho0,bra,ket);
                                
                                rel_cm_u3st_map[braket_u3st]
                                +=u3::dim(x0)*am::dim(Lp)/u3::dim(xp)/am::dim(L0)
                                  *u3::W(u3::SU3(Nrp,0),1,Lrp,u3::SU3(Ncm,0),1,Lcm,xp,kappap,Lp,1)
                                  *u3::W(u3::SU3(Nr,0),1,Lr,u3::SU3(Ncm,0),1,Lcm,x,kappa,L,1)
                                  *u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0)
                                  *rme;
                                std::cout<<fmt::format("{} {} {} {} {} {} {} || {} {} {} || {} {} {} {} {} {} {}",
                                  Nrp,Lrp,Ncm,Lcm,xp.Str(), kappap,Lp, x0.Str(), kappa0,L0, Nr,Lr,Ncm,Lcm,x.Str(), kappa, L)<<std::endl;
                                std::cout<<fmt::format("   {:12f} {:12f} {:12f} {:12f} {:12f}",
                                  u3::W(u3::SU3(Nrp,0),1,Lrp,u3::SU3(Ncm,0),1,Lcm,xp,kappap,Lp,1),
                                  u3::W(u3::SU3(Nr,0),1,Lr,u3::SU3(Ncm,0),1,Lcm,x,kappa,L,1),
                                  u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0),
                                  rme,
                                  u3::dim(x0)*am::dim(Lp)/u3::dim(xp)/am::dim(L0)
                                  *u3::W(u3::SU3(Nrp,0),1,Lrp,u3::SU3(Ncm,0),1,Lcm,xp,kappap,Lp,1)
                                  *u3::W(u3::SU3(Nr,0),1,Lr,u3::SU3(Ncm,0),1,Lcm,x,kappa,L,1)
                                  *u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0)
                                  *rme
                                  )<<std::endl;
                              }
                      }
                  }
            }
        }
    }
}

void RelativeUnitTensorToRelativeCMUnitTensorU3ST(int Nmax,  
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  std::map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorExpansion>& unit_relative_cm_map)
// Coupling on center of mass by branching and re-upcoupling
{
  for(int i=0; i<relative_unit_tensors.size(); ++i)
    {
      const u3shell::RelativeUnitTensorLabelsU3ST& tensor=relative_unit_tensors[i];

      std::map<RelativeBraketNLST,double> branched_rel_unit_tensors;
      BranchNLST(tensor,branched_rel_unit_tensors);

      std::map<RelativeCMBraketNLST,double> rel_cm_lst_map;
      RelativeToCMLST(Nmax, branched_rel_unit_tensors,rel_cm_lst_map);

      RelativeCMUnitTensorExpansion rel_cm_u3st_map;
      UpcoupleCMU3ST(rel_cm_lst_map,rel_cm_u3st_map);

      unit_relative_cm_map[tensor]=rel_cm_u3st_map;
    }
}

void RelativeToCMU3ST(int Nmax,  
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  std::vector<RelativeCMUnitTensorExpansion>& unit_tensor_rel_cm_expansions
  )
// Coupling on center of mass at U(3) level
{
  int N0, Np,N, kappa0,L0;
  HalfInt S0,T0,Sp,Tp,S,T;
  u3::SU3 x0;
  u3shell::RelativeUnitTensorLabelsU3ST tensor;
  for(auto tensor:relative_unit_tensors)
    {
      RelativeCMUnitTensorExpansion relative_cm_u3st_map;
      std::tie(x0,S0,T0,Np,Sp,Tp,N,S,T)=tensor.FlatKey();
      u3::SU3 xr(N,0);
      u3::SU3 xrp(Np,0);
      for(int Ncm=0; Ncm<=Nmax; Ncm++)
        {
          u3::SU3 x_cm(Ncm,0);
          MultiplicityTagged<u3::SU3>::vector x_set=u3::KroneckerProduct(xr,x_cm);
          MultiplicityTagged<u3::SU3>::vector xp_set=u3::KroneckerProduct(xrp,x_cm);
          for(auto ip : xp_set)
            for(auto i: x_set)
              {
                u3::SU3 x(i.irrep);
                u3::SU3 xp(ip.irrep);
                int rho0_max=u3::OuterMultiplicity(x,x0,xp);
                u3shell::RelativeCMStateLabelsU3ST bra(Np,Ncm,xp,Sp,Tp);
                u3shell::RelativeCMStateLabelsU3ST ket(N,Ncm,x,S,T);

                for(int rho0=1; rho0<=rho0_max; ++rho0)
                {
                  u3shell::RelativeCMUnitTensorLabelsU3ST tensor_cm(x0,S0,T0,rho0,bra,ket);
                  relative_cm_u3st_map[tensor_cm]=u3::U(x0,xr,xp,x_cm,xrp,1,1,x,1,rho0);
                }
              }
        }
      unit_tensor_rel_cm_expansions.push_back(relative_cm_u3st_map);
    }
}

typedef  std::map<u3shell::RelativeUnitTensorLabelsU3ST, u3shell::TwoBodyUnitTensorCoefficientsU3ST> TwoBodyExpansionMap;

void RelativeUnitTensorToTwobodyU3ST(int Nmax,  
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  std::map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorExpansion>& unit_relative_cm_map,
  TwoBodyExpansionMap& two_body_expansion_vector
  )
{
  u3shell::TwoBodySpaceU3ST space(Nmax);
  for(auto tensor : relative_unit_tensors)
    {
      RelativeCMUnitTensorExpansion& unit_relative_cm_expansion=unit_relative_cm_map[tensor];
      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion;
      MoshinskyTransformUnitTensor(tensor, unit_relative_cm_expansion,space, two_body_expansion,"AS");
      two_body_expansion_vector[tensor]=two_body_expansion;
    }
}

void
RelativeUnitTensorToTwobodyU3ST(int Nmax,  
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  TwoBodyExpansionMap& two_body_expansion_vector
  )
// Via U(3) coupling to center of mass 
{
  u3shell::TwoBodySpaceU3ST space(Nmax);
  for(auto tensor : relative_unit_tensors)
    {
      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion;
      double expansion_coef=1.0;
      MoshinskyTransformUnitTensor( tensor, expansion_coef,space,two_body_expansion, "AS");
      two_body_expansion_vector[tensor]=two_body_expansion;
    }
}

typedef std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> IndexedRelativeUnitTensorLabelsU3ST;
typedef  std::map<u3shell::RelativeUnitTensorLabelsU3ST, u3shell::TwoBodyUnitTensorCoefficientsU3ST> TwoBodyExpansionMap;
typedef std::tuple<u3shell::TwoBodyUnitTensorLabelsU3ST, int, int>IndexTwoBodyTensorLabelsU3ST;
typedef std::map<IndexTwoBodyTensorLabelsU3ST,double>IndexedTwoBodyTensorRMEsU3ST;

void
ContractRelativeAndTwoBodyUnitTensorRME(
  std::map<IndexedRelativeUnitTensorLabelsU3ST,double>& relative_rmes,
  TwoBodyExpansionMap& two_body_expansion_map,
  IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes
  )
{
  int kappa0,L0;
  u3shell::RelativeUnitTensorLabelsU3ST tensor;
  for(auto it=relative_rmes.begin(); it!=relative_rmes.end(); ++it)
  {
    IndexedRelativeUnitTensorLabelsU3ST indexed_tensor=it->first;
    double rel_rme=it->second;
    std::tie(tensor,kappa0,L0)=it->first;//indexed_tensor;
    u3shell::TwoBodyUnitTensorCoefficientsU3ST& 
      two_body_expansion=two_body_expansion_map[tensor];
    // For each two-body operator, multiply by rel_rme and add to accumlating tb_rme
    // Summing over Nr and Ncm with delta condition on S and T
    for(auto it2=two_body_expansion.begin(); it2!=two_body_expansion.end();  ++it2)
      {
        u3shell::TwoBodyUnitTensorLabelsU3ST two_body_tensor=it2->first;
        double coef=it2->second;
        IndexTwoBodyTensorLabelsU3ST indexed_two_body_tensor(two_body_tensor,kappa0,L0);
        indexed_two_body_rmes[indexed_two_body_tensor]+=rel_rme*coef;
      }
  }
}


typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> TwoBodyStateLabelsLST;
typedef std::tuple<int,HalfInt,HalfInt,TwoBodyStateLabelsLST,TwoBodyStateLabelsLST> TwoBodyBraketLST;

void BranchTwoBodyNLST(
  IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
  std::map<TwoBodyBraketLST,double>& two_body_rmes_lst
  )
{
  for(auto it=indexed_two_body_rmes.begin(); it!=indexed_two_body_rmes.end(); ++it)
    {
      u3shell::TwoBodyUnitTensorLabelsU3ST tensor_u3st;
      int kappa0,L0;

      IndexTwoBodyTensorLabelsU3ST indexed_two_body_tensor=it->first;
      double rme_u3st=it->second;
      std::tie(tensor_u3st,kappa0,L0)=indexed_two_body_tensor;
      u3::SU3 x0=tensor_u3st.x0();
      HalfInt S0=tensor_u3st.S0();
      HalfInt T0=tensor_u3st.T0();
      int rho0=tensor_u3st.rho0();
      u3shell::TwoBodyStateLabelsU3ST bra=tensor_u3st.bra();
      u3shell::TwoBodyStateLabelsU3ST ket=tensor_u3st.ket();
      int N1,N2,N1p,N2p;
      HalfInt Sp,S,Tp,T;
      u3::SU3 xp,x;
      std::tie(N1p,N2p,xp,Sp,Tp)=bra.Key();
      std::tie(N1,N2,x,S,T)=ket.Key();
      MultiplicityTagged<int>::vector L_branch=BranchingSO3(x);
      MultiplicityTagged<int>::vector Lp_branch=BranchingSO3(xp);
      for(auto lp: Lp_branch)
        {
          int Lp=lp.irrep;
          int kappap_max=lp.tag;
          for(auto l: L_branch)
            {
              int L=l.irrep;
              int kappa_max=l.tag;
              // std::cout<<etap%2<<"  "<<etap<<"  "<<eta%2<<eta<<"  "<<std::endl;
              for(int kappap=1; kappap<=kappap_max; ++kappap)
                for(int kappa=1; kappa<=kappa_max; ++kappa)
                  for(int L1p=N1p%2; L1p<=N1p; L1p+=2)
                    for(int L2p=N2p%2; L2p<=N2p; L2p+=2)
                      for(int L1=N1%2; L1<=N1; L1+=2)  
                        for(int L2=N2%2; L2<=N2; L2+=2)
                          {     
                            if((abs(L1-L2)>L)||((L1+L2)<L)) //am::triangular
                              continue;
                            if((abs(L1p-L2p)>Lp)||((L1p+L2p)<Lp))
                              continue;
                            TwoBodyStateLabelsLST bra(N1p,L1p,N2p,L2p,Lp,Sp,Tp);
                            TwoBodyStateLabelsLST ket(N1, L1,N2,L2,L,S,T);                            
                            TwoBodyBraketLST braket(L0,S0,T0,bra,ket);
                            int n1=(N1-L1)/2, n2=(N2-L2)/2,n1p=(N1p-L1p)/2, n2p=(N2p-L2p)/2;
                            double rme_lst=rme_u3st*u3::W(u3::SU3(N1p,0),1,L1p,u3::SU3(N2p,0),1,L2p,xp,kappap,Lp,1)
                                            *u3::W(u3::SU3(N1,0),1,L1,u3::SU3(N2,0),1,L2,x,kappa,L,1)
                                            *u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0)
                                            *parity(n1+n2+n1p+n2p);
                 
                            two_body_rmes_lst[braket]+=rme_lst;
                          }
            }
        }

    }

}

typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt,HalfInt> TwoBodyStateLabelsLSJT;
typedef std::pair<TwoBodyStateLabelsLSJT,TwoBodyStateLabelsLSJT>TwoBodyBraketLSJT;

void BranchTwoBodyLSJT( int Jmax, int J0,
  const std::map<TwoBodyBraketLST,double>& two_body_rme_lst,
  std::map<TwoBodyBraketLSJT,double>&two_body_rme_lsjt
  )
{
  for(auto it=two_body_rme_lst.begin(); it!=two_body_rme_lst.end(); it++)
    {
      double rme_lst=it->second;
      if (fabs(rme_lst)<10e-10)
        continue;
      TwoBodyStateLabelsLST bra_lst,ket_lst;
      HalfInt S0,T0,S,T,Sp,Tp;
      int eta1,eta2,eta1p,eta2p,L1,L2,L1p,L2p,L,Lp,L0;
      std::tie(L0,S0,T0,bra_lst,ket_lst)=it->first;
      std::tie(eta1,L1,eta2,L2,L,S,T)=ket_lst;
      std::tie(eta1p,L1p,eta2p,L2p,Lp,Sp,Tp)=bra_lst;

      int J_max=std::min(Jmax,int(L+S));
      int Jp_max=std::min(Jmax,int(Lp+Sp));
      for(HalfInt J=abs(L-S); J<=J_max; ++J)
          for(HalfInt Jp=abs(Lp-Sp); Jp<=Jp_max; ++Jp)
            {
              if((J0<abs(Jp-J)) || (J0>(Jp+J)))
                continue;

              double so3coef=am::Unitary9J(L, S, J, L0,S0,J0,Lp,Sp,Jp);
              TwoBodyStateLabelsLSJT state(eta1,L1,eta2,L2,L,S,J,T);
              TwoBodyStateLabelsLSJT statep(eta1p,L1p,eta2p,L2p,Lp,Sp,Jp,Tp);
              TwoBodyBraketLSJT braket(statep,state);
              two_body_rme_lsjt[braket]+=so3coef*rme_lst;
              }
    }
}

typedef std::tuple<int, int, HalfInt,int, int, HalfInt,HalfInt,HalfInt>TwoBodyStateLabelsJJJT;
typedef std::pair<TwoBodyStateLabelsJJJT, TwoBodyStateLabelsJJJT> TwoBodyBraketJJJT;

void branchJJJT(
  std::map<TwoBodyBraketLSJT,double>&two_body_rme_lsjt,
  std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt
  )
{
 for(auto it=two_body_rme_lsjt.begin(); it!=two_body_rme_lsjt.end(); ++it)
  {
    HalfInt S0,T0,S,T,Sp,Tp,Jp,J,J1p,J2p,J1,J2;
    int N1,N2,N1p,N2p,L1,L2,L1p,L2p,L,Lp,L0;
    TwoBodyStateLabelsLSJT bra,ket;
    std::tie(bra,ket)=it->first;
    std::tie(N1,L1,N2,L2,L,S,J,T)=ket;
    std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp)=bra;
    double rme=it->second;
    // std::cout<<fmt::format("[{} {} {} {}] {} {} {} {} | |[{} {} {} {}] {} {} {} {}     {}",
    //   N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,N1,L1,N2,L2,L,S,J,T,rme)<<std::endl;

    for(J1p=abs(L1p-HalfInt(1,2)); J1p<=(L1p+HalfInt(1,2)); ++J1p)
      for(J2p=abs(L2p-HalfInt(1,2)); J2p<=(L2p+HalfInt(1,2)); ++J2p)
        for(J1=abs(L1-HalfInt(1,2)); J1<=(L1+HalfInt(1,2)); ++J1)
          for(J2=abs(L2-HalfInt(1,2)); J2<=(L2+HalfInt(1,2)); ++J2)
            {
              double norm_factor=1.0;
              if((N1==N2)&&(L1==L2)&&(J1==J2))
                norm_factor*=1/sqrt(2);
              if((N1p==N2p)&&(L1p==L2p)&&(J1p==J2p))
                norm_factor*=1/sqrt(2);

              if(not am::AllowedTriangle(J1,J2,J))
                continue;
              if(not am::AllowedTriangle(J1p,J2p,Jp))
                continue;
              double coefJp=am::Unitary9J(L1p,HalfInt(1,2),J1p,L2p,HalfInt(1,2),J2p,Lp,Sp,Jp);
              double coefJ=am::Unitary9J(L1,HalfInt(1,2),J1,L2,HalfInt(1,2),J2,L,S,J);

              // std::cout<<J1p<<" "<<J2p<<J1<<" "<<J2<<"    "<<rme<<"  "<<norm_factor<<"  "<<coefJp<<"  "<<coefJ<<std::endl;

              TwoBodyStateLabelsJJJT braJ(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp);
              TwoBodyStateLabelsJJJT ketJ(N1,L1,J1,N2,L2,J2,J,T);
              TwoBodyBraketJJJT braket(braJ, ketJ);
              two_body_rme_jjjt[braket]+=rme*coefJp*coefJ*norm_factor;
            }
  }
 
}

void VerboseTest(int Nmax, int Jmax, int J0,std::map<IndexedRelativeUnitTensorLabelsU3ST,double>& number_relative_rmes)
{
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensors);

  std::map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorExpansion> unit_relative_cm_map;
  RelativeUnitTensorToRelativeCMUnitTensorU3ST(Nmax,relative_unit_tensors,unit_relative_cm_map);
  // RelativeUnitTensorToRelativeCMUnitTensorU3ST(2*Nmax,relative_unit_tensors,unit_relative_cm_map);


  std::cout<<std::endl<<"Relative-CM unit tensor expansion coefficients"<<std::endl;
  for(auto it=unit_relative_cm_map.begin(); it!=unit_relative_cm_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST tensor=it->first;
      RelativeCMUnitTensorExpansion expansion=it->second;
      std::cout<<tensor.Str()<<std::endl;
      for(auto it2=expansion.begin(); it2!=expansion.end(); ++it2)
        {
          u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st=it2->first;
          double coef=it2->second;
          if(fabs(coef)>10e-8)
            std::cout<<"  "<<braket_u3st.Str()<<"  "<<coef<<std::endl;
        }
    }

 
  TwoBodyExpansionMap two_body_expansion_map;
  RelativeUnitTensorToTwobodyU3ST(
      Nmax, relative_unit_tensors,
      unit_relative_cm_map,
      two_body_expansion_map
    );

  std::cout<<std::endl<<"Two-Body unit tensor expansion coefficients"<<std::endl;
  for(auto tensor : relative_unit_tensors)
    {
      u3shell::TwoBodyUnitTensorCoefficientsU3ST& expansion=two_body_expansion_map[tensor];
      std::cout<<tensor.Str()<<std::endl;
      for(auto it=expansion.begin(); it!=expansion.end(); ++it)
        {
          u3shell::TwoBodyUnitTensorLabelsU3ST tb_tensor=it->first;
          double coef=it->second;
          std::cout<<"  "<<tb_tensor.Str()<<"  "<<coef<<std::endl;
        }
    }

  std::cout<<std::endl<<"Contracting reltive rme's with two-body expansion of relative unit tensors "<<std::endl;
  IndexedTwoBodyTensorRMEsU3ST indexed_two_body_rmes;
  ContractRelativeAndTwoBodyUnitTensorRME(number_relative_rmes,two_body_expansion_map,indexed_two_body_rmes);
  for(auto it=indexed_two_body_rmes.begin(); it!=indexed_two_body_rmes.end(); ++it)
    {
      IndexTwoBodyTensorLabelsU3ST indexed_tensor=it->first;
      double rme=it->second;

      int kappa0, L0; 
      u3shell::TwoBodyUnitTensorLabelsU3ST two_body_tensor;
      std::tie(two_body_tensor,kappa0,L0)=indexed_tensor;
      std::cout<<two_body_tensor.Str()<<"  "<<kappa0<<"  "<<L0<<"    "<<rme<<std::endl;

    }

  std::map<TwoBodyBraketLST,double> two_body_rmes_lst;

  std::cout<<std::endl<<"Branching"<<std::endl;
  BranchTwoBodyNLST(indexed_two_body_rmes,two_body_rmes_lst);

  for(auto it=two_body_rmes_lst.begin(); it!=two_body_rmes_lst.end(); ++it)
    {
      TwoBodyBraketLST braket=it->first;
      double rme=it->second;
      TwoBodyStateLabelsLST bra,ket;
      int L0;
      HalfInt S0,T0;
      std::tie(L0,S0,T0,bra,ket)=braket;
      int N1,L1,N2,L2,L,N1p,L1p,N2p,L2p,Lp;
      HalfInt Sp,Tp,S,T;
      std::tie(N1,L1,N2,L2,L,S,T)=ket;
      std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Tp)=bra;
      std::cout<<fmt::format("{} {} {} {} {} {} {} || {} {} {} || {} {} {} {} {} {} {}   {}",
        N1p,L1p,N2p,L2p,Lp,Sp,Tp,L0,S0,T0,N1,L1,N2,L2,L,S,T,rme)<<std::endl;
    }

  std::cout<<"Branch to NLSJT"<<std::endl;

  std::map<TwoBodyBraketLSJT,double> two_body_rme_lsjt;

  BranchTwoBodyLSJT(4, 0, two_body_rmes_lst,two_body_rme_lsjt);
  for(auto it=two_body_rme_lsjt.begin(); it!=two_body_rme_lsjt.end(); ++it)
  {
    TwoBodyBraketLSJT braket=it->first;
    double rme=it->second;
    TwoBodyStateLabelsLSJT bra,ket;
    std::tie(bra,ket)=braket;
    int N1,N2,N1p,N2p,L1,L2,L1p,L2p,L,Lp;
    HalfInt Sp,Tp,Jp,S,T,J;
    std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp)=bra;
    std::tie(N1,L1,N2,L2,L,S,J,T)=ket;
    std::cout<<fmt::format("[{} {} {} {}] {} {} {} {} | |[{} {} {} {}] {} {} {} {}     {}",
      N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,N1,L1,N2,L2,L,S,J,T,rme)<<std::endl;
  }
  std::map<std::tuple<int,int,HalfInt>,int> h2_lookup;
  H2FormatLookUp(Nmax,h2_lookup);

  std::cout<<std::endl<<"Branch to JJJT"<<std::endl;
  std::map<TwoBodyBraketJJJT,double> two_body_rme_jjjt;
  branchJJJT(two_body_rme_lsjt,two_body_rme_jjjt);
  for(auto it=two_body_rme_jjjt.begin(); it!=two_body_rme_jjjt.end(); ++it)
    {
      TwoBodyStateLabelsJJJT bra,ket;
      std::tie(bra,ket)=it->first;
      double rme=it->second;
      int N1p,N2p,N1,N2,L1,L2,L1p,L2p;
      HalfInt J1,J2,J1p,J2p, Jp,J,Tp,T;
      std::tie(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp)=bra;
      std::tie(N1,L1,J1,N2,L2,J2,J,T)=ket;

      // std::cout<<fmt::format("[{} {} {} {} {} {}]{} {} || [{} {} {} {} {} {}]{} {}     {} ",
      //   N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp, N1,L1,J1,N2,L2,J2,J,T,rme
      //   )<<std::endl;
      if(Tp==0)
        continue;
      std::tuple<int,int,HalfInt> lookup1p(N1p,L1p,J1p);
      std::tuple<int,int,HalfInt> lookup2p(N2p,L2p,J2p);
      std::tuple<int,int,HalfInt> lookup1(N1,L1,J1);
      std::tuple<int,int,HalfInt> lookup2(N2,L2,J2);

      int a1=h2_lookup[lookup1];
      int a2=h2_lookup[lookup2];
      int a1p=h2_lookup[lookup1p];
      int a2p=h2_lookup[lookup2p];
      std::cout<<fmt::format("{} {} {} {}   {} 11   {}", a1p,a2p,a1,a2,TwiceValue(J),rme)<<std::endl;
    }

}

void BranchTwoBodyU3STToJJJT(int Jmax, int J0,
  IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
  std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt
  )
{
  std::map<TwoBodyBraketLST,double> two_body_rmes_lst;
  BranchTwoBodyNLST(indexed_two_body_rmes,two_body_rmes_lst);

  std::map<TwoBodyBraketLSJT,double> two_body_rme_lsjt;
  BranchTwoBodyLSJT(Jmax, J0,two_body_rmes_lst,two_body_rme_lsjt);

  branchJJJT(two_body_rme_lsjt, two_body_rme_jjjt);
}


void 
ConvertRelativeTensorToTwoBodyTensor(int Nmax,
  std::map<IndexedRelativeUnitTensorLabelsU3ST,double>& relative_rmes,
  IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes
  )
{
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensors);

  std::map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorExpansion> 
    unit_relative_cm_map;
  std::cout<<"Relative to CM"<<std::endl;  
  RelativeUnitTensorToRelativeCMUnitTensorU3ST(Nmax, relative_unit_tensors, unit_relative_cm_map);
  std::cout<<"CM to two-body"<<std::endl;
  TwoBodyExpansionMap two_body_expansion_map;
  RelativeUnitTensorToTwobodyU3ST(Nmax, relative_unit_tensors,
   // unit_relative_cm_map,
   two_body_expansion_map);

  ContractRelativeAndTwoBodyUnitTensorRME(relative_rmes,two_body_expansion_map,indexed_two_body_rmes);
}




//////////

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  int Nmax=4; 
  int Jmax=4; 
  int J0=0;
   std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  int Nop=Nmax;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nop, relative_unit_tensors);
  //////// Identity////////////////////////////////////////////////////////////////////////
  std::map<IndexedRelativeUnitTensorLabelsU3ST,double> identity_relative_rmes;
  for(auto tensor : relative_unit_tensors)
    {
      if(
        (tensor.x0()==u3::SU3(0,0))
        &&(tensor.S0()==0)
        &&(tensor.T0()==0)
        &&(tensor.bra()==tensor.ket())
        )
        {
          IndexedRelativeUnitTensorLabelsU3ST index_tensor(tensor,1,0);
          identity_relative_rmes[index_tensor]=1;
        }
    }
  // IndexedTwoBodyTensorRMEsU3ST indexed_two_body_rmes;
  // ContractRelativeAndTwoBodyUnitTensorRME(identity_relative_rmes,two_body_expansion_map,indexed_two_body_rmes);
  ////Relative Number Operator/////////////////////////////////////////////////////////////////////


  std::map<IndexedRelativeUnitTensorLabelsU3ST,double> number_relative_rmes;
  for(auto tensor : relative_unit_tensors)
    {
      if(
        (tensor.x0()==u3::SU3(0,0))
        &&(tensor.S0()==0)
        &&(tensor.T0()==0)
        &&(tensor.bra()==tensor.ket())
        )
        {
          IndexedRelativeUnitTensorLabelsU3ST index_tensor(tensor,1,0);
          number_relative_rmes[index_tensor]=tensor.bra().eta();
        }
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::string interaction_file="/Users/annamccoy/projects/spncci/libraries/u3shell/test/ksqr_Nmax04_rel.dat";
  // std::string interaction_file="/Users/annamccoy/projects/spncci/data/jisp16_Nmax20_hw20.0_rel.dat";

  // Defining containers for reading in interaction
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  basis::OperatorLabelsJT operator_labels;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  //Read interaction and store sector information in relative_component_sectors
  // and matrix elements in relative_component_matrices
  basis::ReadRelativeOperatorLSJT(
    interaction_file,relative_lsjt_space,operator_labels,
    relative_component_sectors,relative_component_matrices, true
    );

  // Extract out T0=0 sectors and matrix elements
  const basis::MatrixVector& sector_vector=relative_component_matrices[0];
  const basis::RelativeSectorsLSJT& relative_lsjt_sectors=relative_component_sectors[0];

  //upcouple to LST
  int g0=0, T0=0;
  std::cout<<"Upcoupling to NLST"<<std::endl;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  u3shell::UpcouplingNLST(relative_lsjt_space,relative_lsjt_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);

  // Upcouple to U(3) level
  u3shell::RelativeRMEsU3ST rme_map;
  std::cout<<"Upcoupling to U3ST"<<std::endl;
  u3shell::UpcouplingU3ST(rme_nlst_map, T0, Nmax, rme_map);
  for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST op_labels;
      int kappa0,L0;
      std::tie(op_labels, kappa0,L0)=it->first;
      double rme=it->second;
      double check=u3shell::RelativeKineticEnergyOperator(op_labels.bra(), op_labels.ket());
      if(fabs(rme)>10e-10)
        std::cout<<fmt::format("{} {} {}   {}   {}",op_labels.Str(), kappa0,L0,rme,check)<<std::endl;
    }



  u3shell::RelativeStateLabelsU3ST bra(0,0,1);
  u3shell::RelativeStateLabelsU3ST ket(2,0,1);
  u3::SU3 x0(0,2);
  u3shell::RelativeUnitTensorLabelsU3ST relative_tensor(x0,0,0,bra,ket);
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors_2;
  relative_unit_tensors_2.push_back(relative_tensor);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  // u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensors);

  std::cout<<std::endl;
  std::map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorExpansion> unit_relative_cm_map;
  std::cout<<"Relative to CM"<<std::endl;
  RelativeUnitTensorToRelativeCMUnitTensorU3ST(Nmax, relative_unit_tensors_2,unit_relative_cm_map);
  // for(auto it=unit_relative_cm_map.begin(); it!=unit_relative_cm_map.end(); ++it)
  //   {
  //     u3shell::RelativeUnitTensorLabelsU3ST tensor=it->first;
  //     RelativeCMUnitTensorExpansion expansion=it->second;
  //     std::cout<<tensor.Str()<<std::endl;
  //     for(auto it2=expansion.begin(); it2!=expansion.end(); ++it2)
  //     {
          
  //         u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st=it2->first;
  //         double coef=it2->second;
  //         if(fabs(coef)>10e-8)
  //           std::cout<<"  "<<braket_u3st.Str()<<"  "<<coef<<std::endl;
  //     }
  //   }

  std::vector<RelativeCMUnitTensorExpansion> unit_tensor_rel_cm_expansions;
  RelativeToCMU3ST(Nmax, relative_unit_tensors_2,unit_tensor_rel_cm_expansions);
  for(int i=0; i<unit_tensor_rel_cm_expansions.size(); ++i)
    {

      u3shell::RelativeUnitTensorLabelsU3ST tensor=relative_unit_tensors_2[i];
      RelativeCMUnitTensorExpansion expansion1=unit_tensor_rel_cm_expansions[i];
      RelativeCMUnitTensorExpansion expansion2=unit_relative_cm_map[tensor];
      std::cout<<tensor.Str()<<std::endl;
      std::cout<<"size "<<expansion1.size()<<"  "<<expansion2.size()<<std::endl;
      
      for(auto it=expansion2.begin(); it!=expansion2.end(); ++it)
        {
          u3shell::RelativeCMUnitTensorLabelsU3ST tensor_cm=it->first;
          double rme2=it->second;
          double rme1=expansion1[tensor_cm];
          if((fabs(rme2)<10e-10)&&(fabs(rme1)<10e-10))
            continue;
          std::cout<<"  "<<tensor_cm.Str()<<"  "<<rme1<<"  "<<rme2<<std::endl;
        }
    }



  // TwoBodyExpansionMap two_body_expansion_map;
  // RelativeUnitTensorToTwobodyU3ST(Nmax,relative_unit_tensors,two_body_expansion_map);
  // RelativeUnitTensorToTwobodyU3ST(Nmax,relative_unit_tensors, unit_relative_cm_map,two_body_expansion_map);
  // for(auto it=two_body_expansion_map.begin(); it!=two_body_expansion_map.end(); ++it)
  // {
  //   u3shell::RelativeUnitTensorLabelsU3ST rel_tensor=it->first;
  //   u3shell::TwoBodyUnitTensorCoefficientsU3ST expansion=it->second;
  //   std::cout<<rel_tensor.Str()<<std::endl;
  //   for(auto i=expansion.begin(); i!=expansion.end(); ++i)
  //     {
  //       u3shell::TwoBodyUnitTensorLabelsU3ST tb_tensor(i->first);
  //       double coef=i->second;
  //       std::cout<<"  "<<tb_tensor.Str()<<"  "<<coef<<std::endl;
  //     }
  // }

  // IndexedTwoBodyTensorRMEsU3ST indexed_two_body_rmes;
  // ContractRelativeAndTwoBodyUnitTensorRME(rme_map, two_body_expansion_map,indexed_two_body_rmes);
  // for(auto it=indexed_two_body_rmes.begin(); it!=indexed_two_body_rmes.end(); ++it)
  //   {
  //     int kappa0,L0;
  //     u3shell::TwoBodyUnitTensorLabelsU3ST tb_tensor;
  //     std::tie(tb_tensor,kappa0,L0)=it->first;
  //     double rme=it->second;
  //     if(fabs(rme)>10e-10)
  //       std::cout<<fmt::format("{} {} {}   {}",tb_tensor.Str(), kappa0,L0,rme)<<std::endl;
  //   }

  // std::cout<<"Convering to Two-body"<<std::endl;
  // IndexedTwoBodyTensorRMEsU3ST indexed_two_body_rmes2;
  // ConvertRelativeTensorToTwoBodyTensor(Nmax,rme_map,indexed_two_body_rmes2);
  // // for(auto it=indexed_two_body_rmes2.begin(); it!=indexed_two_body_rmes2.end(); ++it)
  // //     {
  // //       int kappa0,L0;
  // //       u3shell::TwoBodyUnitTensorLabelsU3ST tb_tensor;
  // //       std::tie(tb_tensor,kappa0,L0)=it->first;
  // //       double rme=it->second;
  // //       if(fabs(rme)>10e-10)
  // //         std::cout<<fmt::format("{} {} {}   {}",tb_tensor.Str(), kappa0,L0,rme)<<std::endl;
  // //     }
  // std::cout<<"Branching "<<std::endl;
  // std::map<TwoBodyBraketJJJT,double> two_body_rme_jjjt;
  // BranchTwoBodyU3STToJJJT(Jmax, J0,indexed_two_body_rmes2,two_body_rme_jjjt);


  // std::map<std::tuple<int,int,HalfInt>,int> h2_lookup;
  // H2FormatLookUp(Nmax,h2_lookup);



  // int a, b, c, d, JJ, TT;
  // double trme;
  // std::string line;
  // std::map<std::tuple<int,int,int,int,int,int>,double> test_map;
  // std::string file="/Users/annamccoy/projects/spncci/libraries/u3shell/test/tbme-Trel.dat";
  // std::ifstream stream(file.c_str());
  // while(std::getline(stream,line))
  // {
  //   std::istringstream(line)>>a>>b>>c>>d>>JJ>>TT>>trme;
  //   std::tuple<int,int,int,int,int,int> key(a,b,c,d,JJ,TT);
  //   test_map[key]=trme/10.;
  // }



  // for(auto it=two_body_rme_jjjt.begin(); it!=two_body_rme_jjjt.end(); ++it)
  //   {
  //     TwoBodyStateLabelsJJJT bra,ket;
  //     std::tie(bra,ket)=it->first;
  //     if(bra>ket)
  //       continue;
  //     double rme=it->second;
  //     int N1p,N2p,N1,N2,L1,L2,L1p,L2p;
  //     HalfInt J1,J2,J1p,J2p, Jp,J,Tp,T;
  //     std::tie(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp)=bra;
  //     std::tie(N1,L1,J1,N2,L2,J2,J,T)=ket;

  //     if(Tp==0)
  //       continue;
  //     ///////////////////
  //     if(J!=0)
  //       continue;
  //     if((N1+N2)%2==1)
  //       continue;

  //     std::tuple<int,int,HalfInt> lookup1p(N1p,L1p,J1p);
  //     std::tuple<int,int,HalfInt> lookup2p(N2p,L2p,J2p);
  //     std::tuple<int,int,HalfInt> lookup1(N1,L1,J1);
  //     std::tuple<int,int,HalfInt> lookup2(N2,L2,J2);

  //     int a1=h2_lookup[lookup1];
  //     int a2=h2_lookup[lookup2];
  //     if(a1>a2)
  //       continue;
  //     int a1p=h2_lookup[lookup1p];
  //     int a2p=h2_lookup[lookup2p];
  //     if(a1p>a2p)
  //       continue;
  //     if (fabs(rme)>10e-10) 
  //       {
  //         std::tuple<int,int,int,int,int,int> key(a1p,a2p,a1,a2,TwiceValue(J),11);
  //         double trme=test_map[key];
  //         std::cout<<fmt::format("{:3} {:3} {:3} {:3}   {:3}   11   {:12f}  {:12f}   {:12f}", a1p,a2p,a1,a2,TwiceValue(J),rme,trme, fabs(rme/trme))<<std::endl;

  //       }
  //   }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // Direct route 
  // std::vector< u3shell::TwoBodyUnitTensorCoefficientsU3ST > two_body_expansion_vector2;
  // RelativeUnitTensorToTwobodyU3ST(Nmax,  relative_unit_tensors, two_body_expansion_vector2);

  // std::cout<<"  "<<std::endl;
  // for(int i=0; i<relative_unit_tensors.size(); ++i)
  //   {
  //     u3shell::RelativeUnitTensorLabelsU3ST tensor(relative_unit_tensors[i]);
  //     u3shell::TwoBodyUnitTensorCoefficientsU3ST& expansion=two_body_expansion_vector2[i];
  //     std::cout<<tensor.Str()<<std::endl;
  //     for(auto it=expansion.begin(); it!=expansion.end(); ++it)
  //       {
  //         u3shell::TwoBodyUnitTensorLabelsU3ST tb_tensor=it->first;
  //         double coef=it->second;
  //         std::cout<<"  "<<tb_tensor.Str()<<"  "<<coef<<std::endl;
  //       }
  //   }




}



  // for(auto relative_tensor:relative_unit_tensors)
  // {
  // std::map<RelativeBraketNLST,double> branched_rel_unit_tensors;
  // BranchNLST(relative_tensor, branched_rel_unit_tensors);
  // for(auto it=branched_rel_unit_tensors.begin(); it!=branched_rel_unit_tensors.end(); ++it)
  //   {
  //     int L0;
  //     HalfInt S0,T0;
  //     RelativeStateLabelsNLST bra,ket;
  //     std::tie(L0,S0,T0,bra,ket)=it->first;
  //     double rme=it->second;
  //     int Np,N,Lp,L;
  //     HalfInt Sp,S,Tp,T;
  //     std::tie(Np,Lp,Sp,Tp)=bra;
  //     std::tie(N,L,S,T)=ket;
  //     // std::cout<<fmt::format("{} {} {} {} || {} {} {} || {} {} {} {}    {}",Np,Lp,Sp,Tp,L0,S0,T0,N,L,S,T,rme)<<std::endl;
  //   }

  // std::map<RelativeCMBraketNLST,double> rel_cm_lst_map;
  // RelativeToCMLST(Nmax, branched_rel_unit_tensors,rel_cm_lst_map);
  // for(auto it=rel_cm_lst_map.begin(); it!=rel_cm_lst_map.end(); ++it)
  //   {
  //     int Nr,Ncm,Lr,Lcm,Nrp,Lrp,L0,L,Lp;
  //     HalfInt Sp,Tp,S,T,S0,T0;
  //     RelativeCMStateLabelsNLST bra,ket;
  //     std::tie(L0,S0,T0,bra,ket)=it->first;
  //     std::tie(Nr,Lr,Ncm,Lcm,L,S,T)=ket;
  //     std::tie(Nrp,Lrp,Ncm,Lcm,Lp,Sp,Tp)=bra;
  //     double rme=it->second;
  //     // std::cout<<fmt::format("[{} {} {} {}]{} {} {}|| {} {} {} ||[{} {} {} {}]{} {} {}    {}",
  //     //   Nrp,Lrp,Ncm,Lcm,Lp,Sp,Tp,L0,S0,T0,Nr,Lr,Ncm,Lcm,L,S,T,rme)<<std::endl;

  //   }

  // RelativeCMUnitTensorExpansion rel_cm_u3st_map;
  // UpcoupleCMU3ST(rel_cm_lst_map,rel_cm_u3st_map);
  // for(auto it=rel_cm_u3st_map.begin(); it!=rel_cm_u3st_map.end(); ++it)
  //   {
  //     u3shell::RelativeCMUnitTensorLabelsU3ST tensor=it->first;
  //     double rme=it->second;
  //     // std::cout<<tensor.Str()<<"  "<<rme<<std::endl;

  //   }

  // // std::map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorExpansion> unit_relative_cm_map;
  // // unit_relative_cm_map[relative_tensor]=rel_cm_u3st_map;

  // std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_tensor_vec;
  // relative_tensor_vec.push_back(relative_tensor);

  // std::map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorExpansion> unit_relative_cm_map;
  // RelativeUnitTensorToRelativeCMUnitTensorU3ST(Nmax,  relative_tensor_vec,unit_relative_cm_map);

  // TwoBodyExpansionMap two_body_expansion_map;
  // RelativeUnitTensorToTwobodyU3ST(Nmax,relative_tensor_vec, unit_relative_cm_map,two_body_expansion_map);
  // for(auto it=two_body_expansion_map.begin(); it!=two_body_expansion_map.end(); ++it)
  // {
  //   u3shell::RelativeUnitTensorLabelsU3ST rel_tensor=it->first;
  //   u3shell::TwoBodyUnitTensorCoefficientsU3ST expansion=it->second;
  //   std::cout<<rel_tensor.Str()<<std::endl;
  //   for(auto i=expansion.begin(); i!=expansion.end(); ++i)
  //     {
  //       u3shell::TwoBodyUnitTensorLabelsU3ST tb_tensor(i->first);
  //       double coef=i->second;
  //       std::cout<<"  "<<tb_tensor.Str()<<"  "<<coef<<std::endl;
  //     }
  // }
  // }
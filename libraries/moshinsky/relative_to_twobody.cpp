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
// #include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
// #include "u3shell/import_interaction.h"
#include "u3shell/relative_operator.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/moshinsky.h"
#include "moshinsky/relative_cm_xform.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/upcoupling.h"


typedef std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> IndexedRelativeUnitTensorLabelsU3ST;
// typedef  std::map<u3shell::RelativeUnitTensorLabelsU3ST, u3shell::TwoBodyUnitTensorCoefficientsU3ST> TwoBodyExpansionMap;
typedef  std::unordered_map<u3shell::RelativeUnitTensorLabelsU3ST, 
                            u3shell::TwoBodyUnitTensorCoefficientsU3ST,
                            boost::hash<u3shell::RelativeUnitTensorLabelsU3ST>
                          > TwoBodyExpansionMap;

typedef std::tuple<u3shell::TwoBodyUnitTensorLabelsU3ST, int, int>IndexTwoBodyTensorLabelsU3ST;
typedef std::map<IndexTwoBodyTensorLabelsU3ST,double>IndexedTwoBodyTensorRMEsU3ST;
typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt,HalfInt> TwoBodyStateLabelsLSJT;
typedef std::pair<TwoBodyStateLabelsLSJT,TwoBodyStateLabelsLSJT>TwoBodyBraketLSJT;
typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> TwoBodyStateLabelsLST;
typedef std::tuple<int,HalfInt,HalfInt,TwoBodyStateLabelsLST,TwoBodyStateLabelsLST> TwoBodyBraketLST;


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

void RelativeUnitTensorToTwobodyU3ST(int Nmax,  
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  u3shell::RelativeCMExpansion& unit_relative_cm_map,
  TwoBodyExpansionMap& two_body_expansion_vector
  )
{
  u3shell::TwoBodySpaceU3ST space(Nmax);
  for(auto tensor : relative_unit_tensors)
    {
      u3shell::RelativeCMUnitTensorCache& unit_relative_cm_expansion=unit_relative_cm_map[tensor];
      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion;
      u3shell::MoshinskyTransformUnitTensor(tensor, unit_relative_cm_expansion,space, two_body_expansion,"AS");
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
      u3shell::MoshinskyTransformUnitTensor( tensor, expansion_coef,space,two_body_expansion, "AS");
      two_body_expansion_vector[tensor]=two_body_expansion;
    }
}

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
      if(fabs(rme_u3st)<10e-10)
        continue;
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
                            // if((N1==2)&&(N2==2))
                            //   std::cout<<fmt::format("{} {} {} {} {} {} {}   {} {} {} {} {} {} {}   {} {} {} {} {}", 
                            //     N1p,L1p,N2p,L2p,Lp,Sp,Tp,N1, L1,N2,L2,L,S,T,
                            //     rme_u3st,
                            //     u3::W(u3::SU3(N1p,0),1,L1p,u3::SU3(N2p,0),1,L2p,xp,kappap,Lp,1),
                            //     u3::W(u3::SU3(N1,0),1,L1,u3::SU3(N2,0),1,L2,x,kappa,L,1),
                            //     u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0),
                            //     parity(n1+n2+n1p+n2p)
                            //     )<<std::endl;

                            two_body_rmes_lst[braket]+=rme_lst;
                          }
            }
        }
    }
}

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

  u3shell::RelativeCMExpansion unit_relative_cm_map;
  u3shell::RelativeUnitTensorToRelativeCMUnitTensorU3ST(Nmax,relative_unit_tensors,unit_relative_cm_map);

  // std::cout<<std::endl<<"Relative-CM unit tensor expansion coefficients"<<std::endl;
  // for(auto it=unit_relative_cm_map.begin(); it!=unit_relative_cm_map.end(); ++it)
  //   {
  //     u3shell::RelativeUnitTensorLabelsU3ST tensor=it->first;
  //     u3shell::RelativeCMUnitTensorCache expansion=it->second;
  //     std::cout<<tensor.Str()<<std::endl;
  //     for(auto it2=expansion.begin(); it2!=expansion.end(); ++it2)
  //       {
  //         u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st=it2->first;
  //         double coef=it2->second;
  //         if(fabs(coef)>10e-8)
  //           std::cout<<"  "<<braket_u3st.Str()<<"  "<<coef<<std::endl;
  //       }
  //   }

 
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

  u3shell::RelativeCMExpansion 
    unit_relative_cm_map;
  RelativeUnitTensorToRelativeCMUnitTensorU3ST(Nmax, relative_unit_tensors, unit_relative_cm_map);

  TwoBodyExpansionMap two_body_expansion_map;
  RelativeUnitTensorToTwobodyU3ST(Nmax, relative_unit_tensors,unit_relative_cm_map,two_body_expansion_map);

  ContractRelativeAndTwoBodyUnitTensorRME(relative_rmes,two_body_expansion_map,indexed_two_body_rmes);
}

void 
GetInteractionMatrix(
    std::string interaction_file, 
    basis::RelativeSpaceLSJT& relative_lsjt_space,
    basis::RelativeSectorsLSJT& relative_lsjt_sectors,
    basis::MatrixVector& sector_vector
  )
{
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
  sector_vector=relative_component_matrices[0];
  relative_lsjt_sectors=relative_component_sectors[0];
  }

void PrintTwoBodyMatrixElementsJJJT(const std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt)
{
  for(auto it=two_body_rme_jjjt.begin(); it!=two_body_rme_jjjt.end(); ++it)
    {
      TwoBodyStateLabelsJJJT bra,ket;
      std::tie(bra,ket)=it->first;
      if(bra>ket)
        continue;
      double rme=it->second;
      int N1p,N2p,N1,N2,L1,L2,L1p,L2p;
      HalfInt J1,J2,J1p,J2p, Jp,J,Tp,T;
      std::tie(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp)=bra;
      std::tie(N1,L1,J1,N2,L2,J2,J,T)=ket;
      if(fabs(rme)>10e-10)
        std::cout<<fmt::format("{} {} {} {} {} {} {} {}   {} {} {} {} {} {} {} {}   {}",
          N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp,N1,L1,J1,N2,L2,J2,J,T,rme)<<std::endl;
    }
}
void PrintTwoBodyIndexedRMEU3ST(const IndexedTwoBodyTensorRMEsU3ST& two_body_rmes)
{
  for(auto it=two_body_rmes.begin(); it!=two_body_rmes.end(); ++it)
      {
        int kappa0,L0;
        u3shell::TwoBodyUnitTensorLabelsU3ST tb_tensor;
        std::tie(tb_tensor,kappa0,L0)=it->first;
        double rme=it->second;
        if(fabs(rme)>10e-10)
          std::cout<<fmt::format("{} {} {}   {}",tb_tensor.Str(), kappa0,L0,rme)<<std::endl;
      }
}
//////////

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  int Nmax=6; 
  int Jmax=4; 
  int J0=0;
  
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensors);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //// Identity Operator 
  ///////////////////////////////////////////////////////////////////////////////////////////////////  
//   std::map<IndexedRelativeUnitTensorLabelsU3ST,double> identity_relative_rmes;

//   for(auto tensor : relative_unit_tensors)
//     {
//       if(
//         (tensor.x0()==u3::SU3(0,0))
//         &&(tensor.S0()==0)
//         &&(tensor.T0()==0)
//         &&(tensor.bra()==tensor.ket())
//         )
//         {
//           IndexedRelativeUnitTensorLabelsU3ST index_tensor(tensor,1,0);
//           identity_relative_rmes[index_tensor]=1;
//         }
//     }
//   IndexedTwoBodyTensorRMEsU3ST id_two_body_rmes;
//   ConvertRelativeTensorToTwoBodyTensor(Nmax,identity_relative_rmes, id_two_body_rmes);
//   // PrintTwoBodyIndexedRMEU3ST(id_two_body_rmes);

//   std::map<TwoBodyBraketLST,double> two_body_rmes_lst;
//   BranchTwoBodyNLST(id_two_body_rmes,two_body_rmes_lst);
//   // std::cout<<"Branch NLST"<<std::endl;
//   // for(auto it=two_body_rmes_lst.begin(); it!=two_body_rmes_lst.end(); ++it)
//   //   {
//   //     TwoBodyBraketLST braket=it->first;
//   //     double rme=it->second;
//   //     TwoBodyStateLabelsLST bra,ket;
//   //     int L0;
//   //     HalfInt S0,T0;
//   //     std::tie(L0,S0,T0,bra,ket)=braket;
//   //     int N1,L1,N2,L2,L,N1p,L1p,N2p,L2p,Lp;
//   //     HalfInt Sp,Tp,S,T;
//   //     std::tie(N1,L1,N2,L2,L,S,T)=ket;
//   //     std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Tp)=bra;
//   //     if(fabs(rme)>10e-10)
//   //     std::cout<<fmt::format("{} {} {} {} {} {} {} || {} {} {} || {} {} {} {} {} {} {}   {}",
//   //       N1p,L1p,N2p,L2p,Lp,Sp,Tp,L0,S0,T0,N1,L1,N2,L2,L,S,T,rme)<<std::endl;
//   //   }

//   std::map<TwoBodyBraketJJJT,double> id_two_body_rme_jjjt;
//   BranchTwoBodyU3STToJJJT(Jmax, J0, id_two_body_rmes, id_two_body_rme_jjjt);

//   // std::cout<<"Identity test "<<std::endl;
//   // PrintTwoBodyMatrixElementsJJJT(id_two_body_rme_jjjt);
//   ///////////////////////////////////////////////////////////////////////////////////////////////////
//   ////Relative Number Operator/////////////////////////////////////////////////////////////////////
//   ///////////////////////////////////////////////////////////////////////////////////////////////////
//   std::map<IndexedRelativeUnitTensorLabelsU3ST,double> number_relative_rmes;
//   for(auto tensor : relative_unit_tensors)
//     {
//       if(
//         (tensor.x0()==u3::SU3(0,0))
//         &&(tensor.S0()==0)
//         &&(tensor.T0()==0)
//         &&(tensor.bra()==tensor.ket())
//         )
//         {
//           IndexedRelativeUnitTensorLabelsU3ST index_tensor(tensor,1,0);
//           number_relative_rmes[index_tensor]=tensor.bra().eta();
//         }
//     }
  
//   IndexedTwoBodyTensorRMEsU3ST num_two_body_rmes;
//   ConvertRelativeTensorToTwoBodyTensor(Nmax,number_relative_rmes, num_two_body_rmes);
//   // std::cout<<"Relative Number operator two-body expansion"<<std::endl;
//   // PrintTwoBodyIndexedRMEU3ST(num_two_body_rmes);

// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // std::string interaction_file="/Users/annamccoy/projects/spncci/libraries/moshinsky/test/ksqr_Nmax06_rel.dat";
  std::string interaction_file="/Users/annamccoy/projects/spncci/data/jisp16_Nmax20_hw20.0_rel.dat";
  basis::RelativeSectorsLSJT relative_lsjt_sectors;
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  basis::MatrixVector sector_vector;
  GetInteractionMatrix(interaction_file, relative_lsjt_space,relative_lsjt_sectors,sector_vector);

  //upcouple to LST
  int g0=0, T0=0;
  std::cout<<"Upcoupling to NLST"<<std::endl;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  u3shell::UpcouplingNLST(relative_lsjt_space,relative_lsjt_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);

  // Upcouple to U(3) level
  u3shell::RelativeRMEsU3ST rme_map;
  std::cout<<"Upcoupling to U3ST"<<std::endl;
  u3shell::UpcouplingU3ST(rme_nlst_map, T0, Nmax, rme_map);
  //Check for Kinetic energy only
  // for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
  //   {
  //     u3shell::RelativeUnitTensorLabelsU3ST op_labels;
  //     int kappa0,L0;
  //     std::tie(op_labels, kappa0,L0)=it->first;
  //     double rme=it->second;
  //     double check=u3shell::RelativeKineticEnergyOperator(op_labels.bra(), op_labels.ket());
  //     if(fabs(rme)>10e-10)
  //       std::cout<<fmt::format("{} {} {}   {}   {}",op_labels.Str(), kappa0,L0,rme,check)<<std::endl;
  //   }

//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   u3shell::RelativeStateLabelsU3ST bra(0,0,1);
//   u3shell::RelativeStateLabelsU3ST ket(2,0,1);
//   u3::SU3 x0(0,2);
//   u3shell::RelativeUnitTensorLabelsU3ST relative_tensor(x0,0,0,bra,ket);
//   std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors_2;
//   relative_unit_tensors_2.push_back(relative_tensor);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout<<"Convering to Two-body"<<std::endl;
  IndexedTwoBodyTensorRMEsU3ST k2_indexed_two_body_rmes;
  ConvertRelativeTensorToTwoBodyTensor(Nmax,rme_map,k2_indexed_two_body_rmes);
  // PrintTwoBodyIndexedRMEU3ST(k2_indexed_two_body_rmes);  

  std::map<TwoBodyBraketLST,double> k2_two_body_rmes_lst;
  BranchTwoBodyNLST(k2_indexed_two_body_rmes,k2_two_body_rmes_lst);

  std::map<TwoBodyBraketLSJT,double> k2_two_body_rme_lsjt;
  BranchTwoBodyLSJT(Jmax, J0,k2_two_body_rmes_lst,k2_two_body_rme_lsjt);

  // Read in matrix elements from file
  // int N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,gp,N1,L1,N2,L2,L,S,J,T,g;
  double trme;
  std::string line;
  // std::map<TwoBodyBraketLSJT,double> test_map;
  // // std::string file="/Users/annamccoy/projects/spncci/libraries/moshinsky/test/ksqr_Nmax06_lsjt_AS.dat";
  // std::string file="/Users/annamccoy/projects/spncci/libraries/moshinsky/test/jisp16_Nmax20_hw20.0_lsjt_AS.dat";

  // std::ifstream stream(file.c_str());
  // while(std::getline(stream,line))
  // {
  //   std::istringstream(line)>>T0>> N1p>>L1p>>N2p>>L2p>>Lp>>Sp>>Jp>>Tp>>gp>>N1>>L1>>N2>>L2>>L>>S>>J>>T>>g>>trme;
  //   TwoBodyStateLabelsLSJT bra( N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp);
  //   TwoBodyStateLabelsLSJT ket(N1,L1,N2,L2,L,S,J,T);
  //   TwoBodyBraketLSJT key(bra,ket);
  //   test_map[key]=trme;
  //   double rme=k2_two_body_rme_lsjt[key];
  //   if((N1p+N2p)>Nmax)
  //     continue;
  //   if((N1+N2)>Nmax)
  //     continue;
  //   if(fabs(trme)<10e-10)
  //     continue;
  //   if(fabs(trme-rme)>10e-10)
  //     std::cout<<fmt::format("[{} {} {} {}] {} {} {} {} | |[{} {} {} {}] {} {} {} {} {:12f} {:12f}",
  //       N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,N1,L1,N2,L2,L,S,J,T,test_map[key],rme
  //       )<<std::endl;

  // }
  // stream.close();
 




  std::cout<<"Branching "<<std::endl;
  std::map<TwoBodyBraketJJJT,double> two_body_rme_jjjt;
  BranchTwoBodyU3STToJJJT(Jmax, J0,k2_indexed_two_body_rmes,two_body_rme_jjjt);

  std::map<std::tuple<int,int,HalfInt>,int> h2_lookup;
  H2FormatLookUp(Nmax,h2_lookup);

  // Read in matrix elements from file
  int a, b, c, d, JJ, TT;
  // double trme;
  // std::string line;
  std::map<std::tuple<int,int,int,int,int,int>,double> test_map_jj;
  // std::string file_jj="/Users/annamccoy/projects/spncci/libraries/moshinsky/test/tbme-Trel.dat";
  std::string file_jj="/Users/annamccoy/projects/spncci/libraries/moshinsky/test/JISP16-tb-6-20.dat";

  std::ifstream stream_jj(file_jj.c_str());
  while(std::getline(stream_jj,line))
  {
    std::istringstream(line)>>a>>b>>c>>d>>JJ>>TT>>trme;
    std::tuple<int,int,int,int,int,int> key(a,b,c,d,JJ,TT);
    // test_map_jj[key]=trme/10.;//Kinetic energy
    test_map_jj[key]=trme;
  }

  for(auto it=two_body_rme_jjjt.begin(); it!=two_body_rme_jjjt.end(); ++it)
    {
      TwoBodyStateLabelsJJJT bra,ket;
      std::tie(bra,ket)=it->first;
      if(bra>ket)
        continue;
      double rme=it->second;
      int N1p,N2p,N1,N2,L1,L2,L1p,L2p;
      HalfInt J1,J2,J1p,J2p, Jp,J,Tp,T;
      std::tie(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp)=bra;
      std::tie(N1,L1,J1,N2,L2,J2,J,T)=ket;

      if(Tp==0)
        continue;
      ///////////////////
      // if(J!=0)
      //   continue;
      if((N1+N2)%2==1)
        continue;

      std::tuple<int,int,HalfInt> lookup1p(N1p,L1p,J1p);
      std::tuple<int,int,HalfInt> lookup2p(N2p,L2p,J2p);
      std::tuple<int,int,HalfInt> lookup1(N1,L1,J1);
      std::tuple<int,int,HalfInt> lookup2(N2,L2,J2);

      int a1=h2_lookup[lookup1];
      int a2=h2_lookup[lookup2];
      if(a1>a2)
        continue;
      int a1p=h2_lookup[lookup1p];
      int a2p=h2_lookup[lookup2p];
      if(a1p>a2p)
        continue;
      if (fabs(rme)>10e-10) 
        {
          std::tuple<int,int,int,int,int,int> key(a1p,a2p,a1,a2,TwiceValue(J),22);
          double trme=test_map_jj[key];
          if(fabs(rme-trme)>10e-6)
            std::cout<<fmt::format("{:3} {:3} {:3} {:3}   {:3}   22   {:12f}  {:12f}   {:12f}", 
              a1p,a2p,a1,a2,TwiceValue(J),rme,trme, fabs(rme/trme))<<std::endl;
        }
    }
}
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

  // u3shell::RelativeCMUnitTensorCache rel_cm_u3st_map;
  // UpcoupleCMU3ST(rel_cm_lst_map,rel_cm_u3st_map);
  // for(auto it=rel_cm_u3st_map.begin(); it!=rel_cm_u3st_map.end(); ++it)
  //   {
  //     u3shell::RelativeCMUnitTensorLabelsU3ST tensor=it->first;
  //     double rme=it->second;
  //     // std::cout<<tensor.Str()<<"  "<<rme<<std::endl;

  //   }

  // // RelativeCMExpansion unit_relative_cm_map;
  // // unit_relative_cm_map[relative_tensor]=rel_cm_u3st_map;

  // std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_tensor_vec;
  // relative_tensor_vec.push_back(relative_tensor);

  // RelativeCMExpansion unit_relative_cm_map;
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
/****************************************************************
  interaction_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
****************************************************************/
#include <fstream>


#include "cppformat/format.h"
#include "basis/lsjt_operator.h"

#include "am/wigner_gsl.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/import_interaction.h"
#include "u3shell/moshinsky.h"
#include "u3shell/relative_operator.h"
#include "u3shell/upcoupling.h"


// eta1,L1, eta2,L2,L,S,J,T
typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt,HalfInt> TwoBodyStateLabelsLSJT;
typedef std::pair<TwoBodyStateLabelsLSJT,TwoBodyStateLabelsLSJT>TwoBodyRMELablesLSJT;

typedef std::tuple<u3shell::TwoBodyStateLabelsU3ST,int,int>TwoBodyStateLabelsU3STBranched;
typedef std::tuple<TwoBodyStateLabelsU3STBranched,TwoBodyStateLabelsU3STBranched,int,HalfInt,HalfInt>
  TwoBodyBraketU3STBranched;

typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> TwoBodyStateLabelsLSTBranched;
typedef std::tuple<TwoBodyStateLabelsLSTBranched,TwoBodyStateLabelsLSTBranched,int, HalfInt,HalfInt>TwoBodyBraketLSTBranched;


typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt,HalfInt> TwoBodyStateLabelsLSJTBranched;
typedef std::pair<TwoBodyStateLabelsLSJTBranched,TwoBodyStateLabelsLSJTBranched>TwoBodyBraketBranched;


typedef std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,u3shell::TwoBodyStateLabelsU3ST,u3shell::TwoBodyStateLabelsU3ST>TwoBodyRMEU3ST;
typedef std::tuple<u3shell::TwoBodyStateLabelsU3ST,u3shell::TwoBodyStateLabelsU3ST,int> TwoBodyBraket;
typedef std::tuple<u3::SU3,HalfInt,HalfInt,int,u3shell::TwoBodyStateLabelsU3ST,u3shell::TwoBodyStateLabelsU3ST, int> RelativeOperatorRMEU3ST;

typedef std::tuple<int, int, HalfInt,int, int, HalfInt,HalfInt,HalfInt>TwoBodyStateLabelsJJJTBranched;
typedef std::pair<TwoBodyStateLabelsJJJTBranched, TwoBodyStateLabelsJJJTBranched> TwoBodyBraketJJJT;


void KineticEnergyRMEMAP(
  const u3shell::RelativeSpaceU3ST& space,
  std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int>,double>& rme_map)
{
 
  double rme=0;
  for(int ket_index=0; ket_index<space.size(); ++ket_index)
    {
      int eta, S,T;
      std::tie(eta,S,T,std::ignore)=space.GetSubspace(ket_index).Key();
      for(int bra_index=0; bra_index<space.size(); ++bra_index)
        {
          int etap,Sp,Tp;
          std::tie(etap,Sp,Tp,std::ignore)=space.GetSubspace(bra_index).Key();
          u3shell::RelativeStateLabelsU3ST ket(eta,S,T);
          u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp);
          if(S!=Sp)
            continue;
          if(T!=Tp)
            continue;
          if((etap-eta)==2)
            {
              if(u3::OuterMultiplicity(u3::SU3(eta,0),u3::SU3(2,0),u3::SU3(etap,0))!=0)
                {
                  rme=u3shell::RelativeKineticEnergyOperator(bra,ket);
                  u3shell::RelativeUnitTensorLabelsU3ST labels(u3::SU3(2,0),0,0,bra,ket);
                  std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int> key(labels,1);
                  rme_map[key]=rme;
                }
            }

          if(etap==eta)
            {
            if(u3::OuterMultiplicity(u3::SU3(eta,0),u3::SU3(0,0),u3::SU3(etap,0))!=0)
                {
                  rme=u3shell::RelativeKineticEnergyOperator(bra,ket);
                  u3shell::RelativeUnitTensorLabelsU3ST labels(u3::SU3(0,0),0,0,bra,ket);
                  std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int> key(labels,1);
                  rme_map[key]=rme;
                }
            }
          if(etap-eta==-2)
          {
            if(u3::OuterMultiplicity(u3::SU3(eta,0),u3::SU3(0,2),u3::SU3(etap,0))!=0)
                {
                  rme=u3shell::RelativeKineticEnergyOperator(bra,ket);
                  u3shell::RelativeUnitTensorLabelsU3ST labels(u3::SU3(0,2),0,0,bra,ket);
                  std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int> key(labels,1);
                  rme_map[key]=rme;
                }
          }
        }
    }
}

void ContractRMEUnitTensors(int g0, int T0,
  const std::map< TwoBodyBraket,std::map<u3shell::RelativeUnitTensorLabelsU3ST,double> >& twobody_rme_u3st_map,
  std::map< std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> ,double >& relative_rme_map, 
  u3::WCoefCache& cache,
  std::map<RelativeOperatorRMEU3ST,double>& operator_rme_u3st_map
  )
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // iterating over two-body matrix elements 
  u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor;
  u3shell::TwoBodyStateLabelsU3ST bra_u3st,ket_u3st;
  int rho0;
  // std::map<RelativeOperatorRMEU3ST,double>operator_rme_u3st_map;
  // sum over etap,eta, Sp,S, Tp, T, T0
  // resuts in <omega' S'||V^(omega0 S0 kappa0)||omega S>_rho0
  for(auto it1=twobody_rme_u3st_map.begin(); it1!= twobody_rme_u3st_map.end(); ++it1)
    {
      //extract state labels 
      std::tie(bra_u3st,ket_u3st,rho0)=it1->first;
      const std::map<u3shell::RelativeUnitTensorLabelsU3ST,double>& unit_rme_map=it1->second;
      for(auto it2=unit_rme_map.begin(); it2!=unit_rme_map.end(); ++it2)
        {
          relative_unit_tensor=it2->first;
          u3::SU3 x0=relative_unit_tensor.x0();
          HalfInt S0=relative_unit_tensor.S0();
          HalfInt T0=relative_unit_tensor.T0();
          auto N0=relative_unit_tensor.N0();

          MultiplicityTagged<int>::vector L0_branch=BranchingSO3(x0);

          // std::cout<<N0<<std::endl;
          
          double unit_rme=it2->second;
          for(auto it:L0_branch)
            {
              int L0=it.irrep;
              int kappa0_max=it.tag;
              for(int kappa0=1; kappa0<=kappa0_max; ++kappa0)
                {
                  RelativeOperatorRMEU3ST relative_operator_labels(x0,S0,T0,kappa0,L0,bra_u3st,ket_u3st,rho0);
                  std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> relative_operator(relative_unit_tensor,kappa0,L0);
                  if(relative_rme_map.count(relative_operator)==0)
                    continue;
                  double relative_rme=relative_rme_map[relative_operator];
                  operator_rme_u3st_map[relative_operator_labels]+=unit_rme*relative_rme;
                  // if(N0==-2)
                  //   std::cout<<fmt::format("({}|| {} {} {} {}||{}){}  {}  {}",bra_u3st.Str(),x0.Str(),S0,T0,kappa0,ket_u3st.Str(),rho0,unit_rme,relative_rme)<<std::endl;

                }
            }
        }
    }

  std::cout<<"summing complete"<<std::endl;
}
void branchU3LST(int J0,
  const std::map<RelativeOperatorRMEU3ST,double>& operator_rme_u3st_map,
  u3::WCoefCache& cache,
  std::map<TwoBodyBraketU3STBranched,double>&rme_u3st_branched
  )                      
{
  // sum over omega0 kappa0 rho0
  //std::map<TwoBodyBraketU3STBranched,double> rme_u3st_branched;
  for(auto it=operator_rme_u3st_map.begin(); it!=operator_rme_u3st_map.end(); ++it)
    {
      double rme=it->second;
      if (fabs(rme)<10e-13)
        continue;

      u3::SU3 x0;
      HalfInt S0,T0;
      int rho0,kappa0;
      u3shell::TwoBodyStateLabelsU3ST bra_u3st,ket_u3st;
      std::tie(x0,S0,T0,kappa0,bra_u3st,ket_u3st,rho0)=(it->first);
      u3::SU3 x=ket_u3st.x();
      u3::SU3 xp=bra_u3st.x();
      int eta1p=bra_u3st.eta1(), eta2p=bra_u3st.eta2(), eta1=ket_u3st.eta1(), eta2=ket_u3st.eta2();
      // if((eta1p==1)&&(eta2p==3)&&(eta1==1)&&(eta2==1))
      //   std::cout<<fmt::format("({}|| {} {} {} {}||{}){}  {}",
      //       bra_u3st.Str(),x0.Str(),S0,T0,kappa0,ket_u3st.Str(),rho0,rme)<<std::endl;

      MultiplicityTagged<int>::vector L_branch=u3::BranchingSO3(x);
      MultiplicityTagged<int>::vector Lp_branch=u3::BranchingSO3(xp);
      MultiplicityTagged<int>::vector L0_branch=u3::BranchingSO3(x0);

      for(int l0=0; l0<L0_branch.size(); ++l0)
        {
          int L0=L0_branch[l0].irrep;
          if((L0>(J0+S0))||(L0<abs(S0-J0)))
            continue;
          int kappa0_max=L0_branch[l0].tag;
          if (kappa0>kappa0_max)
            continue;
          for(int l=0; l<L_branch.size(); ++l)
            for(int lp=0; lp<Lp_branch.size(); ++lp)
              {
                int L=L_branch[l].irrep;
                int kappa_max=L_branch[l].tag;
                int Lp=Lp_branch[lp].irrep;
                int kappap_max=Lp_branch[lp].tag;
                for(int kappa=1; kappa<=kappa_max; ++kappa)
                  for(int kappap=1; kappap<=kappap_max; ++kappap)
                    {
                      
                      double wcoef0=WCached(cache,x, kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0);
                      TwoBodyStateLabelsU3STBranched ket(ket_u3st,kappa,L);
                      TwoBodyStateLabelsU3STBranched bra(bra_u3st,kappap,Lp);
                      TwoBodyBraketU3STBranched braket(bra,ket,L0,S0,T0);
                      rme_u3st_branched[braket]+=wcoef0*rme;
                      // if((eta1p==1)&&(eta2p==3)&&(eta1==1)&&(eta2==1))
                      //   std::cout<<fmt::format("    ({} {} {}||{} {} {}|| {} {} {})   {}  {}",bra_u3st.Str(),kappap,Lp,L0,S0, T0,ket_u3st.Str(),kappa,L,wcoef0*rme,rme_u3st_branched[braket])<<std::endl;
                    }
              }
        }
    }
  // u3shell::TwoBodyStateLabelsU3ST bra(0,2,u3::SU3(2,0),1,1), ket(1,3,u3::SU3(2,1),1,1);
  // TwoBodyStateLabelsU3STBranched bra_branched(bra,1,2), ket_branched(ket,1,2);
  // TwoBodyBraketU3STBranched braket(bra_branched,ket_branched,0,0,0);
  // std::cout<<rme_u3st_branched[braket]<<std::endl;

}

void branchLST( 
  const std::map<TwoBodyBraketU3STBranched,double>&rme_u3st_branched,
  u3::WCoefCache& cache,
  std::map<TwoBodyBraketLSTBranched,double>& rme_lst_branched
  )

{
  // std::map<TwoBodyBraketLSTBranched,double> rme_lst_branched;
  for(auto it=rme_u3st_branched.begin(); it!=rme_u3st_branched.end(); ++it)
    {

      double rme=it->second;
      if (fabs(rme)<10e-10)
        continue;

      u3shell::TwoBodyStateLabelsU3ST bra_u3st,ket_u3st;
      TwoBodyStateLabelsU3STBranched bra_u3lst,ket_u3lst;
      HalfInt S0,S,T,Sp,Tp, T0;
      u3::SU3 x,xp; 
      int eta1,eta2,eta1p,eta2p, L,Lp,L0, kappap,kappa;

      std::tie(bra_u3lst,ket_u3lst,L0,S0,T0)=(it->first);

      std::tie(bra_u3st,kappap,Lp)=bra_u3lst;
      std::tie(ket_u3st,kappa,L)=ket_u3lst;
      // std::cout<<fmt::format("({} {} {}||V||{} {} {})  {}",
      //     bra_u3st.Str(),kappap, Lp,ket_u3st.Str(),kappa,L,rme)
      //   <<std::endl;
      std::tie(eta1p,eta2p,xp,Sp,Tp)=bra_u3st.Key();
      std::tie(eta1,eta2,x,S,T)=ket_u3st.Key();
      // std::cout<<fmt::format("({} {} {}||{} {} {}||{} {} {})",
      //   bra_u3st.Str(),kappap,Lp,L0,S0,T0,ket_u3st.Str(),kappa,L)
      // <<std::endl;

      MultiplicityTagged<int>::vector L1_branch=u3::BranchingSO3(u3::SU3(eta1,0));
      MultiplicityTagged<int>::vector L2_branch=u3::BranchingSO3(u3::SU3(eta2,0));
      MultiplicityTagged<int>::vector L1p_branch=u3::BranchingSO3(u3::SU3(eta1p,0));
      MultiplicityTagged<int>::vector L2p_branch=u3::BranchingSO3(u3::SU3(eta2p,0));
       
      //int phase=1;  
        for(int l1=0; l1<L1_branch.size(); ++l1)
          for(int l2=0; l2<L2_branch.size(); ++l2)
            {
              int L1=L1_branch[l1].irrep;
              int L2=L2_branch[l2].irrep;
              double norm_factor=1;
              if((eta1==eta2)&&(L1!=L2))
                    norm_factor=sqrt(2);

              double wcoef=WCached(cache,u3::SU3(eta1,0),1,L1,u3::SU3(eta2,0),1,L2,x,kappa,L,1);
              if(fabs(wcoef)<10e-10)
                continue;
              for(int l1p=0; l1p<L1p_branch.size(); ++l1p)
                for(int l2p=0; l2p<L2p_branch.size();++l2p)
                  {
                    int L1p=L1p_branch[l1p].irrep;
                    int L2p=L2p_branch[l2p].irrep;
                    double norm_factor_p=1;
                    if((eta1p==eta2p)&&(L2p!=L1p))
                          norm_factor_p=sqrt(2);

                    TwoBodyStateLabelsLSTBranched ket_lst(eta1,L1,eta2,L2,L,S,T);
                    TwoBodyStateLabelsLSTBranched bra_lst(eta1p,L1p,eta2p,L2p,Lp,Sp,Tp);
                    TwoBodyBraketLSTBranched braket_lst(bra_lst,ket_lst,L0,S0,T0);
                    double wpcoef=WCached(cache,u3::SU3(eta1p,0),1,L1p, u3::SU3(eta2p,0),1,L2p,xp,kappap,Lp,1);
                    int n1=(eta1-L1)/2, n2=(eta2-L2)/2, n1p=(eta1p-L1p)/2, n2p=(eta2p-L2p)/2;
                    
                    double rme_branched=norm_factor*norm_factor_p*wcoef*wpcoef*rme;
                    // std::cout<<fmt::format("([{} {} {} {}] {} {} {}||{} {} {}||[{} {} {} {}] {} {} {})  {} {} {} {} {} {}",
                    //                 eta1p,L1p,eta2p,L2p,Lp,Sp,Tp,L0,S0,T0,eta1,L1,eta2,L2,L,S,T,norm_factor,parity(eta1+eta2+eta1p+eta2p),wcoef,wpcoef,rme,norm_factor*parity(eta1+eta2+eta1p+eta2p)*wcoef*wpcoef*rme)
                    //     <<std::endl; 
                    
                    rme_lst_branched[braket_lst]+=rme_branched*parity(n1+n2+n1p+n2p);
                  }
            }
    }
}

void branchLSJT( int Jmax, int J0,
  const std::map<TwoBodyBraketLSTBranched,double>& rme_lst_branched,
  std::map<TwoBodyBraketBranched,double>&rme_lsjt_branched
  )
{
  for(auto it=rme_lst_branched.begin(); it!=rme_lst_branched.end(); it++)
    {
      double rme_lst=it->second;
      if (fabs(rme_lst)<10e-10)
        continue;
      TwoBodyStateLabelsLSTBranched bra_lst,ket_lst;
      HalfInt S0,T0,S,T,Sp,Tp;
      int eta1,eta2,eta1p,eta2p,L1,L2,L1p,L2p,L,Lp,L0;
      std::tie(bra_lst,ket_lst,L0,S0,T0)=it->first;
      std::tie(eta1,L1,eta2,L2,L,S,T)=ket_lst;
      std::tie(eta1p,L1p,eta2p,L2p,Lp,Sp,Tp)=bra_lst;
      // if(fabs(rme_lst)>10e-10)
      // std::cout<<fmt::format("hi  ([{} {} {} {}] {} {} {}||{} {} {}||[{} {} {} {}] {} {} {})  {}",
      //               eta1p,L1p,eta2p,L2p,Lp,Sp,Tp,L0,S0,T0,eta1,L1,eta2,L2,L,S,T, rme_lst)
      //   <<std::endl;

      int J_max=std::min(Jmax,int(L+S));
      int Jp_max=std::min(Jmax,int(Lp+Sp));
      for(HalfInt J=abs(L-S); J<=J_max; ++J)
          for(HalfInt Jp=abs(Lp-Sp); Jp<=Jp_max; ++Jp)
            {
              if((J0<abs(Jp-J)) || (J0>(Jp+J)))
                continue;

              double so3coef=am::Unitary9J(L, S, J, L0,S0,J0,Lp,Sp,Jp);
              TwoBodyStateLabelsLSJTBranched state(eta1,L1,eta2,L2,L,S,J,T);
              TwoBodyStateLabelsLSJTBranched statep(eta1p,L1p,eta2p,L2p,Lp,Sp,Jp,Tp);
              TwoBodyBraketBranched braket(statep,state);
              rme_lsjt_branched[braket]+=so3coef*rme_lst;
              }
    }
}

void branchJJJT(
  const std::map<TwoBodyBraketBranched,double>&rme_lsjt_branched,
  std::map<TwoBodyBraketJJJT,double>& branched_rme_jjjt
  )
{
  // std::map<TwoBodyBraketJJJT,double> branched_rme_jjjt;
  for(auto it=rme_lsjt_branched.begin(); it!=rme_lsjt_branched.end(); ++it)
    {
      HalfInt S0,T0,S,T,Sp,Tp,Jp,J,J1p,J2p,J1,J2;
      int eta1,eta2,eta1p,eta2p,L1,L2,L1p,L2p,L,Lp,L0;
      TwoBodyStateLabelsLSJTBranched bra,ket;
      std::tie(bra,ket)=it->first;
      std::tie(eta1,L1,eta2,L2,L,S,J,T)=ket;
      std::tie(eta1p,L1p,eta2p,L2p,Lp,Sp,Jp,Tp)=bra;
      double rme=it->second;
      // std::cout<<fmt::format("([{} {}  {} {}] {} {} {} {}||  ||[{} {}  {} {}] {} {} {} {})  {}",
      //   eta1p,L1p,eta2p,L2p,Lp,Sp,Jp,Tp,eta1,L1,eta2,L2,L,S,J,T,rme)<<std::endl;
      for(J1p=abs(L1p-HalfInt(1,2)); J1p<=(L1p+HalfInt(1,2)); ++J1p)
        for(J2p=abs(L2p-HalfInt(1,2)); J2p<=(L2p+HalfInt(1,2)); ++J2p)
          for(J1=abs(L1-HalfInt(1,2)); J1<=(L1+HalfInt(1,2)); ++J1)
            for(J2=abs(L2-HalfInt(1,2)); J2<=(L2+HalfInt(1,2)); ++J2)
              {
                double norm_factor=1;
                if((eta1==L1)&&(eta2==L2)&&(J1!=J2))
                  norm_factor=sqrt(2);
                double norm_factor_p=1;
                if((eta1p==L1p)&&(eta2p==L2p)&&(J1p!=J2p))
                  norm_factor_p=sqrt(2);
                double coefJp=am::Unitary9J(L1p,HalfInt(1,2),J1p,L2p,HalfInt(1,2),J2p,Lp,Sp,Jp);
                double coefJ=am::Unitary9J(L1,HalfInt(1,2),J1,L2,HalfInt(1,2),J2,L,S,J);

                TwoBodyStateLabelsJJJTBranched braJ(eta1p,L1p,J1p,eta2p,L2p,J2p,Jp,Tp);
                TwoBodyStateLabelsJJJTBranched ketJ(eta1,L1,J1,eta2,L2,J2,J,T);
                TwoBodyBraketJJJT braket(braJ, ketJ);
                branched_rme_jjjt[braket]+=rme*coefJp*coefJ*norm_factor*norm_factor_p;
              }
    }
}



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

void CheckRME(std::istream& is, std::map<TwoBodyBraketBranched,double>& rme_lsjt_branched)
{
    int N1p,N2p,N1,N2,L1p,L2p,L1,L2,L,Lp,gp,g;
    int Sp,Tp,Jp,S,T,J,T0;
    double rme_mark;
    std::string line;
    std::cout<<"Checking"<<std::endl;
    while(std::getline(is,line))
    {
      std::istringstream(line)>>T0>>N1p>>L1p>>N2p>>L2p>>Lp>>Sp>>Jp>>Tp>>gp
                              >>N1>>L1>>N2>>L2>>L>>S>>J>>T>>g>>rme_mark;

      TwoBodyStateLabelsLSJTBranched bra(N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp);
      TwoBodyStateLabelsLSJTBranched ket(N1,L1,N2,L2,L,S,J,T);  
      TwoBodyBraketBranched braket(bra,ket);
      if(rme_lsjt_branched.count(braket)!=0)
        {
          double rme_anna=rme_lsjt_branched[braket];
          if(fabs(rme_anna-rme_mark)>10e-8)
          {
            double factor=rme_mark/rme_anna;
            std::cout<<fmt::format("([{} {} {} {}] {} {} {} {} ||  ||[{} {} {} {}] {} {} {} {})  {:10f}  {:10f}  {:10f}",
              N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,N1,L1,N2,L2,L,S,J,T,rme_anna,rme_mark,factor)<<std::endl;
          }
        } 
      else
        if(fabs(rme_mark)>10e-10)
          std::cout<<fmt::format("missing ([{} {} {} {}] {} {} {} {} ||  ||[{} {} {} {}] {} {} {} {})  {:10f}"
            ,N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,N1,L1,N2,L2,L,S,J,T,rme_mark) 
          <<std::endl;      
    }
}

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  u3::WCoefCache cache;
  
  int Nmax=4;
  int Jmax=Nmax+1;
  int J0=0;
  int g0=0;
  int T0=0;
std::cout<< "Nmax "<<Nmax<<std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative unit tensors 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
  u3shell::RelativeSpaceU3ST relative_u3st_space(Nmax);
  u3shell::TwoBodySpaceU3ST  twobody_u3st_space(Nmax);
  std::map<TwoBodyBraketBranched,double>branched_rme_map;

  std::map<int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>> relative_unit_tensor_labels;
  std::map< TwoBodyBraket,std::map<u3shell::RelativeUnitTensorLabelsU3ST,double> >twobody_rme_u3st_map;
  double expansion_coef=1.0;

  u3shell::GenerateRelativeUnitTensorLabelsU3ST(relative_u3st_space,relative_unit_tensor_labels);
  std::cout<<"Generating complete  "<<relative_unit_tensor_labels.size()<<std::endl;
  for(auto it0=relative_unit_tensor_labels.begin(); it0!=relative_unit_tensor_labels.end(); ++it0)
    {
      // std::cout<<(it0->first)<<std::endl;
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& tensor_list=(it0->second);//relative_unit_tensor_labels[N0];
      for(int t=0; t<tensor_list.size(); ++t)
        {
          u3shell::RelativeUnitTensorLabelsU3ST tensor=tensor_list[t];
          u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion;
          u3shell::MoshinskyTransformUnitTensor(tensor,1.0,twobody_u3st_space,two_body_expansion);
          // iterate through two-body expansion
          for(auto it=two_body_expansion.begin(); it!=two_body_expansion.end(); it++)
            {
              u3shell::TwoBodyUnitTensorLabelsU3ST twobody_tensor=it->first;
              TwoBodyBraket braket(twobody_tensor.bra(), twobody_tensor.ket(),twobody_tensor.rho0());
              twobody_rme_u3st_map[braket][tensor]=it->second;

            }
        }
    }
  std::cout<<"Moshinky transform complete"<<std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // u3shell::TwoBodyStateLabelsU3ST bra,ket;
  // int rho0;
  // int eta1,eta2,eta1p,eta2p;
  // u3::SU3 x,xp;
  // HalfInt S,T,Sp,Tp;

  // for(auto it=twobody_rme_u3st_map.begin(); it!=twobody_rme_u3st_map.end(); ++it)
  //   {
  //     std::tie(bra,ket,rho0)=it->first;
  //     std::tie(eta1,eta2,x,S,T)=ket.Key();
  //     std::tie(eta1p,eta2p,xp,Sp,Tp)=bra.Key();
  //     const std::map<u3shell::RelativeUnitTensorLabelsU3ST,double> & tensor_map=it->second;
  //     std::cout<<fmt::format("([{} {}]{} {} {}||U||[{} {}]{} {} {}){})",eta1p,eta2p,xp.Str(),Sp,Tp,eta1,eta2,x.Str(),S,T,rho0)<<std::endl;
  //     for(auto it2=tensor_map.begin(); it2!=tensor_map.end(); ++it2)
  //       {
  //         std::cout<<fmt::format("   {}  {}",it2->first.Str(),it2->second)<<std::endl;
  //       }
  //   }
////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int Nmax_relative=8;
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax_relative, Jmax);
  std::string interaction_file;
  std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>,double> relative_rme_map; 
  
  //// Identity test
  // basis::RelativeSectorsLSJT relative_lsjt_sectors(relative_lsjt_space, J0,T0, g0);
  // interaction_file="NONE";
  // basis::MatrixVector sector_vector=u3shell::ImportInteraction(interaction_file, relative_lsjt_space, relative_lsjt_sectors, "Identity");

  // Kinetic Energy Test
  
  // interaction_file="NONE";
  // basis::MatrixVector sector_vector=u3shell::ImportInteraction(interaction_file, relative_lsjt_space, relative_lsjt_sectors, "Kinetic");
  // u3shell::Upcoupling(relative_lsjt_space,relative_lsjt_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);
  // std::cout<<"upcoupling complete"<<std::endl;
  // basis::MatrixVector sector_vector;
  // KineticEnergyRMEMAP(relative_u3st_space,relative_rme_map);


  // From File
  // set up operator
  //  interaction_file="/Users/annamccoy/projects/shell/libraries/moshinsky/test/symmunit_Nmax04_rel.dat";
  // interaction_file="/Users/annamccoy/projects/shell/libraries/moshinsky/test/ksqr_Nmax04_rel.dat";

  interaction_file="/Users/annamccoy/projects/spncci/data/jisp16_Nmax20_hw20.0_rel.dat";
  
  basis::OperatorLabelsJT operator_labels;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  
  basis::ReadRelativeOperatorLSJT(
      interaction_file,
      relative_lsjt_space,      
      operator_labels,
      relative_component_sectors,
      relative_component_matrices, 
      true);


  for(int i=0; i<relative_component_matrices[0].size(); ++i)
    {
      basis::RelativeStateLSJT bra(relative_component_sectors[0].GetSector(i).bra_subspace(),0);
      basis::RelativeStateLSJT ket(relative_component_sectors[0].GetSector(i).ket_subspace(),0);
      if (bra.L()==ket.L())
      {
        for(int n=0; n<relative_component_matrices[0][i].cols(); ++n)
          for(int np=0; np<n; ++np)
            {
              double rme=relative_component_matrices[0][i](np,n);
              relative_component_matrices[0][i](n,np)=rme;
              // std::cout<<rme<<std::endl;
            }
      }
    }

  for(int i=0; i<relative_component_matrices[0].size(); ++i)
    {
      basis::RelativeStateLSJT bra(relative_component_sectors[0].GetSector(i).bra_subspace(),0);
      basis::RelativeStateLSJT ket(relative_component_sectors[0].GetSector(i).ket_subspace(),0);

      std::cout<<fmt::format("{} {} {} {} {}", bra.L(), ket.L(), ket.S(), ket.J(), ket.T())
      <<std::endl;
      std::cout<<relative_component_matrices[0][i]<<std::endl;
    }

  const basis::MatrixVector& sector_vector=relative_component_matrices[0];
  const basis::RelativeSectorsLSJT& relative_lsjt_sectors=relative_component_sectors[0];
  u3shell::Upcoupling(relative_lsjt_space,relative_lsjt_sectors,sector_vector,J0,g0,T0,20,relative_rme_map);
  std::cout<<"upcoupling complete"<<std::endl;


  for(auto it=relative_rme_map.begin(); it!=relative_rme_map.end(); ++it)
  {
    int kappa0,L0; 
    u3shell::RelativeUnitTensorLabelsU3ST tensor; 
    std::tie(tensor,kappa0,L0)=it->first;
    std::cout<<tensor.Str()<<"  "<<it->second<<std::endl;
  }

  std::map<RelativeOperatorRMEU3ST,double> operator_rme_u3st_map;
  ContractRMEUnitTensors(g0,T0,twobody_rme_u3st_map,relative_rme_map,cache,operator_rme_u3st_map);
  std::cout<<"Contracting complete"<<std::endl;

  std::map<TwoBodyBraketU3STBranched,double>rme_u3st_branched;
  branchU3LST(J0,operator_rme_u3st_map,cache,rme_u3st_branched);                      
  std::cout<<"U3ST branching complete"<<std::endl;

  std::map<TwoBodyBraketLSTBranched,double> rme_lst_branched;
  branchLST(rme_u3st_branched,cache,rme_lst_branched);
  std::cout<<"LST branching complete"<<std::endl;

  std::map<TwoBodyBraketBranched,double> rme_lsjt_branched;
  branchLSJT(Jmax,J0, rme_lst_branched,rme_lsjt_branched);
  std::cout<<"LSJT branching complete"<<std::endl; 

  // std::cout<<"printing"<<std::endl;
  // for(auto it=rme_lsjt_branched.begin(); it!=rme_lsjt_branched.end(); ++it)
  //   {
  //     int N1,N2,N1p,N2p,L1,L2,L1p,L2p,Lp,L;
  //     HalfInt S,Sp,T,J,Jp,Tp;
  //     TwoBodyStateLabelsLSJTBranched bra,ket;
  //     std::tie(bra,ket)=it->first;
  //     std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp)=bra;
  //     std::tie(N1,L1,N2,L2,L,S,J,T)=ket;
  //     double rme=it->second;
  //     std::cout<<fmt::format("([{} {} {} {}] {} {} {} {} ||  ||[{} {} {} {}] {} {} {} {})  {:6}",
  //       N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,N1,L1,N2,L2,L,S,J,T,rme)<<std::endl;
  //   }
  // std::ifstream is("/Users/annamccoy/projects/shell/libraries/moshinsky/test/symmunit_Nmax04_lsjt_NAS.dat");
  // std::ifstream is("/Users/annamccoy/projects/shell/libraries/moshinsky/test/ksqr_Nmax04_lsjt_NAS.dat");
  // std::ifstream is("/Users/annamccoy/projects/shell/libraries/moshinsky/test/jisp16_Nmax20_hw20.0_lsjt_NAS.dat");

  // if(!is)
  //   std::cout<<"Didn't open"<<std::endl;
  // CheckRME(is, rme_lsjt_branched);
  // std::cout<<"Finished checking"<<std::endl;


  std::map<TwoBodyBraketJJJT,double> branched_rme_jjjt;
  branchJJJT(rme_lsjt_branched,branched_rme_jjjt);

  std::map<std::tuple<int,int,HalfInt>,int> h2_lookup;
  H2FormatLookUp(Nmax, h2_lookup);

  
  for(auto it=branched_rme_jjjt.begin(); it!=branched_rme_jjjt.end(); ++it)
    {
      int N1,N2,N1p,N2p,L1,L2,L1p,L2p;
      HalfInt J1,J2,J1p,J2p,T,J,Jp,Tp;
      TwoBodyStateLabelsJJJTBranched bra,ket;
      //typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt,HalfInt>
      std::tie(bra,ket)=it->first;
      std::tie(N1,L1,J1,N2,L2,J2,J,T)=ket;
      std::tie(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp)=bra;

      // parity even
      if((N1p+N2p)%2==1)
        continue;

      double me=it->second;
      if(fabs(me)>10e-10) 
        if((T==1)&&(J==0))
          {
            int b1=h2_lookup[std::tuple<int,int,HalfInt>(N1p,L1p,J1p)];
            int b2=h2_lookup[std::tuple<int,int,HalfInt>(N2p,L2p,J2p)];
            int k1=h2_lookup[std::tuple<int,int,HalfInt>(N1,L1,J1)];
            int k2=h2_lookup[std::tuple<int,int,HalfInt>(N2,L2,J2)];
            if(b1>k1)
              continue;
            if((b1==k1)&&(b2>k2))
              continue;
            std::cout<<fmt::format("{} {} {} {} {} {}  {}",
              b1,b2,k1,k2,J,T,me)<<std::endl;
            // std::cout<<fmt::format("([{} {} {}; {} {} {}] {} {} ||  ||[{} {} {}; {} {} {}] {} {})  {}",
            //   N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp,N1,L1,J1,N2,L2,J2,J,T,me)<<std::endl;
          }
  }
}

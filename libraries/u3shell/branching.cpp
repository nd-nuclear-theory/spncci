/****************************************************************
  interaction_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
****************************************************************/
#include "cppformat/format.h"

#include "am/wigner_gsl.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/import_interaction.h"
#include "u3shell/moshinsky.h"
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

void branchLST(  int Jmax, int J0, int g0, int T0,
  const std::map< TwoBodyBraket,std::map<u3shell::RelativeUnitTensorLabelsU3ST,double> >& twobody_rme_u3st_map,
  std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int>,double>& relative_rme_map, 
  u3::WCoefCache& cache,
  std::map<TwoBodyBraketBranched,double>&branched_rme_map
  )
{

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // iterating over two-body matrix elements 
  u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor;
  u3shell::TwoBodyStateLabelsU3ST bra_u3st,ket_u3st;
  int rho0;
  std::map<RelativeOperatorRMEU3ST,double>operator_rme_u3st_map;
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
          
          double unit_rme=it2->second;
          
          for(int kappa0=1; kappa0<=9; ++kappa0)
            {
              RelativeOperatorRMEU3ST relative_operator_labels(x0,S0,T0,kappa0,bra_u3st,ket_u3st,rho0);
              std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int> relative_operator(relative_unit_tensor,kappa0);
              if(relative_rme_map.count(relative_operator)==0)
                continue;
              operator_rme_u3st_map[relative_operator_labels]+=unit_rme*relative_rme_map[relative_operator];
            }
        }
    }
  std::cout<<"summing complete"<<std::endl;
  // sum over omega0 kappa0 rho0
  std::map<TwoBodyBraketU3STBranched,double> rme_u3st_branched;
  for(auto it=operator_rme_u3st_map.begin(); it!=operator_rme_u3st_map.end(); ++it)
    {
      double rme=it->second;
      if ((fabs(rme)<10e-10)||(fabs(rme)>10e10))
        continue;
      u3::SU3 x0;
      HalfInt S0,T0;
      int rho0,kappa0;
      std::tie(x0,S0,T0,kappa0,bra_u3st,ket_u3st,rho0)=(it->first);
      u3::SU3 x=ket_u3st.x();
      u3::SU3 xp=bra_u3st.x();

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
                    }
              }
        }
    }


  std::map<TwoBodyBraketLSTBranched,double> rme_lst_branched;
  for(auto it=rme_u3st_branched.begin(); it!=rme_u3st_branched.end(); ++it)
    {

      double rme=it->second;
      if ((fabs(rme)<10e-10)||(fabs(rme)>10e10))
        continue;
      // std::cout<<rme<<std::endl;
      TwoBodyStateLabelsU3STBranched bra_u3lst,ket_u3lst;
      HalfInt S0,T0,S,T,Sp,Tp;
      u3::SU3 x,xp; 
      int eta1,eta2,eta1p,eta2p, L,Lp,L0, kappap,kappa;

      std::tie(bra_u3lst,ket_u3lst,L0,S0,T0)=(it->first);

      std::tie(bra_u3st,kappap,Lp)=ket_u3lst;
      std::tie(ket_u3st,kappa,L)=ket_u3lst;
      // if(fabs(rme)>10e-10)
      //   std::cout<<fmt::format("({} {} {}||V||{} {} {})  {}",
      //     bra_u3st.Str(),kappap, Lp,ket_u3st.Str(),kappa,L,rme)
      //   <<std::endl;


      std::tie(eta1,eta2,x,S,T)=ket_u3st.Key();
      std::tie(eta1p,eta2p,xp,Sp,Tp)=bra_u3st.Key();

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
              if((eta1==eta2))
                {
                  if(L2<L1)
                    continue;
                  if(L2>L1)
                    norm_factor=sqrt(2);
                }
              // if((eta1==eta2)&&(L2<L1))
              //   {
                  // phase=parity(eta1+eta2+x.lambda()+x.mu()+L1+L2+L+int(S+T));
                  // int temp=L1;
                  // L1=L2;
                  // L2=temp;
                // }
              double wcoef=WCached(cache,u3::SU3(eta1,0),1,L1,u3::SU3(eta2,0),1,L2,x,kappa,L,1);
              if(fabs(wcoef)<10e-10)
                continue;
              for(int l1p=0; l1p<L1p_branch.size(); ++l1p)
                for(int l2p=0; l2p<L2p_branch.size();++l2p)
                  {
                    double norm_factor_p=1;
                    int L1p=L1p_branch[l1p].irrep;
                    int L2p=L2p_branch[l2p].irrep;
                    if((eta1p==eta2p))
                      {
                        if(L2p<L1p)
                         continue;
                        if(L2p>L1p)
                          norm_factor_p=sqrt(2);
                      }
                      //  {
                      // phase*=parity(eta1p+eta2p+xp.lambda()+xp.mu()+L1p+L2p+Lp+int(Sp+Tp));
                      // int temp=L1p;
                      // L1p=L2p;
                      // L2p=temp;
                        
                      //  }
                    TwoBodyStateLabelsLSTBranched ket_lst(eta1,L1,eta2,L2,L,S,T);
                    TwoBodyStateLabelsLSTBranched bra_lst(eta1p,L1p,eta2p,L2p,Lp,Sp,Tp);
                    TwoBodyBraketLSTBranched braket_lst(bra_lst,ket_lst,L0,S0,T0);
                    double wpcoef=WCached(cache,u3::SU3(eta1p,0),1,L1p, u3::SU3(eta2p,0),1,L2p,xp,kappap,Lp,1);
                    double rme_branched=norm_factor*norm_factor_p*parity(eta1+eta2+eta1p+eta2p)*wcoef*wpcoef*rme;
                    //if((eta1p==3)&&(eta2==3)&&(L1==1)&&(L2==3)&&(L==3))
                    // std::cout<<fmt::format("([{} {} {} {}] {} {} {}||{} {} {}||[{} {} {} {}] {} {} {})  {} {} {} {} {} {}",
                    //                 eta1p,L1p,eta2p,L2p,Lp,Sp,Tp,L0,S0,T0,eta1,L1,eta2,L2,L,S,T,norm_factor,parity(eta1+eta2+eta1p+eta2p),wcoef,wpcoef,rme,norm_factor*parity(eta1+eta2+eta1p+eta2p)*wcoef*wpcoef*rme)
                    //     <<std::endl; 
                    
                    rme_lst_branched[braket_lst]+=rme_branched;
                  }
            }
    }
  for(auto it=rme_lst_branched.begin(); it!=rme_lst_branched.end(); it++)
    {
      double rme_lst=it->second;
      if ((fabs(rme_lst)<10e-10)||(fabs(rme_lst)>10e10))
        continue;
      TwoBodyStateLabelsLSTBranched bra_lst,ket_lst;
      HalfInt S0,T0,S,T,Sp,Tp;
      int eta1,eta2,eta1p,eta2p,L1,L2,L1p,L2p,L,Lp,L0;
      std::tie(bra_lst,ket_lst,L0,S0,T0)=it->first;
      std::tie(eta1,L1,eta2,L2,L,S,T)=ket_lst;
      std::tie(eta1p,L1p,eta2p,L2p,Lp,Sp,Tp)=bra_lst;
      // if(fabs(rme_lst>10e-10)||fabs(rme_lst)>10e10)
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
              TwoBodyStateLabelsLSJTBranched state(eta1,L1,eta2,L2,L,S,T,J);
              TwoBodyStateLabelsLSJTBranched statep(eta1p,L1p,eta2p,L2p,Lp,Sp,Tp,J);
              TwoBodyBraketBranched braket(statep,state);
              branched_rme_map[braket]+=so3coef*rme_lst;
              }
    }
}


int main(int argc, char **argv)
{
  u3::U3CoefInit();
  u3::WCoefCache cache;
  
  int Nmax=8;
  int Jmax=Nmax+2;
  int J0=0;
  int g0=0;
  int T0=0;
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
  int N0max=relative_unit_tensor_labels.size()-1;
  for(int N0=0; N0<=N0max; ++N0)
    {
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& tensor_list=relative_unit_tensor_labels[N0];
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
  //// Identity test
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  basis::RelativeSectorsLSJT relative_lsjt_sectors(relative_lsjt_space, J0,T0, g0);
  std::vector<Eigen::MatrixXd> sector_vector;
  std::string interaction_file;
  std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int>,double> relative_rme_map; 
  interaction_file="NONE";
  sector_vector=u3shell::ImportInteraction(interaction_file, relative_lsjt_space, relative_lsjt_sectors, "Identity");
  u3shell::Upcoupling(relative_lsjt_space,relative_lsjt_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);
  std::cout<<"upcoupling complete"<<std::endl;
  branchLST(Jmax,J0, g0, T0, twobody_rme_u3st_map, relative_rme_map,cache,  branched_rme_map);
  for(auto it=branched_rme_map.begin(); it!=branched_rme_map.end(); ++it)
    {
      int N1,N2,N1p,N2p, L1,L2,L1p,L2p,L,Lp;
      HalfInt S,T,J,Sp,Jp,Tp;
      TwoBodyStateLabelsLSJT bra,ket;
      //typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt,HalfInt>
      std::tie(bra,ket)=it->first;
      std::tie(N1,L1,N2,L2,L,S,T,J)=ket;
      std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Tp,Jp)=bra;

      double me=it->second;
      if((fabs(me)>10e-10) &&(fabs(me)<10e10))
        std::cout<<fmt::format("([{} {} {} {}] {} {} {} {}||  ||[{} {} {} {}] {} {} {} {})  {}",
          N1,L1,N2,L2,L,S,T,J,N1p,L1p,N2p,L2p,Lp,Sp,Tp,Jp,me)<<std::endl;
  }
}

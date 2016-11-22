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
#include "sp3rlib/u3coef.h"
#include "u3shell/relative_operator.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/moshinsky_xform.h"
#include "moshinsky/relative_cm_xform.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/upcoupling.h"


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

void BranchTwoBodyNLST(
  u3shell::IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
  std::map<TwoBodyBraketLST,double>& two_body_rmes_lst
  )
{
  for(auto it=indexed_two_body_rmes.begin(); it!=indexed_two_body_rmes.end(); ++it)
    {
      u3shell::TwoBodyUnitTensorLabelsU3ST tensor_u3st;
      int kappa0,L0;

      u3shell::IndexTwoBodyTensorLabelsU3ST indexed_two_body_tensor=it->first;
      double rme_u3st=it->second;
      if(fabs(rme_u3st)<10e-13)
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
      if (fabs(rme_lst)<10e-13)
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

void BranchTwoBodyU3STToJJJT(int Jmax, int J0,
  u3shell::IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
  std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt
  )
{
  std::map<TwoBodyBraketLST,double> two_body_rmes_lst;
  BranchTwoBodyNLST(indexed_two_body_rmes,two_body_rmes_lst);

  std::map<TwoBodyBraketLSJT,double> two_body_rme_lsjt;
  BranchTwoBodyLSJT(Jmax, J0,two_body_rmes_lst,two_body_rme_lsjt);

  branchJJJT(two_body_rme_lsjt, two_body_rme_jjjt);
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
      if(fabs(rme)>10e-13)
        std::cout<<fmt::format("{} {} {} {} {} {} {} {}   {} {} {} {} {} {} {} {}   {}",
          N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp,N1,L1,J1,N2,L2,J2,J,T,rme)<<std::endl;
    }
}
void PrintTwoBodyIndexedRMEU3ST(const u3shell::IndexedTwoBodyTensorRMEsU3ST& two_body_rmes)
{
  for(auto it=two_body_rmes.begin(); it!=two_body_rmes.end(); ++it)
      {
        int kappa0,L0;
        u3shell::TwoBodyUnitTensorLabelsU3ST tb_tensor;
        std::tie(tb_tensor,kappa0,L0)=it->first;
        double rme=it->second;
        if(fabs(rme)>10e-13)
          std::cout<<fmt::format("{} {} {}   {}",tb_tensor.Str(), kappa0,L0,rme)<<std::endl;
      }
}
//////////

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  u3::WCoefCache w_cache;

  int Nmax=4; 
  int Jmax=4; 
  int J0=0;
      
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensors);

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //// Identity Operator 
  ///////////////////////////////////////////////////////////////////////////////////////////////////  
  std::map<u3shell::IndexedRelativeUnitTensorLabelsU3ST,double> identity_relative_rmes;

  for(auto tensor : relative_unit_tensors)
    {
      if(
        (tensor.x0()==u3::SU3(0,0))
        &&(tensor.S0()==0)
        &&(tensor.T0()==0)
        &&(tensor.bra()==tensor.ket())
        )
        {
          u3shell::IndexedRelativeUnitTensorLabelsU3ST index_tensor(tensor,1,0);
          identity_relative_rmes[index_tensor]=1;
        }
    }
  u3shell::IndexedTwoBodyTensorRMEsU3ST id_two_body_rmes;
  u3shell::ConvertRelativeTensorToTwoBodyTensor(Nmax,identity_relative_rmes, id_two_body_rmes);
  PrintTwoBodyIndexedRMEU3ST(id_two_body_rmes);

  std::map<TwoBodyBraketLST,double> two_body_rmes_lst;
  BranchTwoBodyNLST(id_two_body_rmes,two_body_rmes_lst);
  // std::cout<<"Branch NLST"<<std::endl;
  // for(auto it=two_body_rmes_lst.begin(); it!=two_body_rmes_lst.end(); ++it)
  //   {
  //     TwoBodyBraketLST braket=it->first;
  //     double rme=it->second;
  //     TwoBodyStateLabelsLST bra,ket;
  //     int L0;
  //     HalfInt S0,T0;
  //     std::tie(L0,S0,T0,bra,ket)=braket;
  //     int N1,L1,N2,L2,L,N1p,L1p,N2p,L2p,Lp;
  //     HalfInt Sp,Tp,S,T;
  //     std::tie(N1,L1,N2,L2,L,S,T)=ket;
  //     std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Tp)=bra;
  //     if(fabs(rme)>10e-10)
  //     std::cout<<fmt::format("{} {} {} {} {} {} {} || {} {} {} || {} {} {} {} {} {} {}   {}",
  //       N1p,L1p,N2p,L2p,Lp,Sp,Tp,L0,S0,T0,N1,L1,N2,L2,L,S,T,rme)<<std::endl;
  //   }

  std::map<TwoBodyBraketJJJT,double> id_two_body_rme_jjjt;
  BranchTwoBodyU3STToJJJT(Jmax, J0, id_two_body_rmes, id_two_body_rme_jjjt);

  std::cout<<"Identity JJJT "<<std::endl;
  PrintTwoBodyMatrixElementsJJJT(id_two_body_rme_jjjt);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // std::string interaction_file="/afs/crc.nd.edu/user/a/amccoy/projects/spncci/data/jisp16_Nmax20_hw20.0_rel.dat";
  std::string interaction_file="../../data/jisp16_Nmax20_hw20.0_rel.dat";
  basis::RelativeSectorsLSJT relative_lsjt_sectors;
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  basis::MatrixVector sector_vector;
  u3shell::GetInteractionMatrix(interaction_file, relative_lsjt_space,relative_lsjt_sectors,sector_vector);

  //upcouple to LST
  int g0=0, T0=0;
  std::cout<<"Upcoupling"<<std::endl;
  u3shell::RelativeRMEsU3ST rme_map;
  u3shell::Upcoupling(relative_lsjt_space,relative_lsjt_sectors, sector_vector, w_cache, J0, g0,T0, Nmax,rme_map);

  std::cout<<"Convering to Two-body"<<std::endl;
  u3shell::IndexedTwoBodyTensorRMEsU3ST j16_indexed_two_body_rmes;
  u3shell::ConvertRelativeTensorToTwoBodyTensor(Nmax,rme_map,j16_indexed_two_body_rmes);
  // PrintTwoBodyIndexedRMEU3ST(k2_indexed_two_body_rmes);  

  std::map<TwoBodyBraketLST,double> j16_two_body_rmes_lst;
  BranchTwoBodyNLST(j16_indexed_two_body_rmes,j16_two_body_rmes_lst);

  std::map<TwoBodyBraketLSJT,double> j16_two_body_rme_lsjt;
  BranchTwoBodyLSJT(Jmax, J0, j16_two_body_rmes_lst,j16_two_body_rme_lsjt);

  // Read in matrix elements from file
  // int N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,gp,N1,L1,N2,L2,L,S,J,T,g;
  double trme;
  std::string line;
  // std::map<TwoBodyBraketLSJT,double> test_map;
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
  BranchTwoBodyU3STToJJJT(Jmax, J0,j16_indexed_two_body_rmes,two_body_rme_jjjt);

  std::map<std::tuple<int,int,HalfInt>,int> h2_lookup;
  H2FormatLookUp(Nmax,h2_lookup);

  // Read in matrix elements from file
  int a, b, c, d, JJ, TT;
  // double trme;
  // std::string line;
  std::map<std::tuple<int,int,int,int,int,int>,double> test_map_jj;
  std::string file_jj="test/JISP16-tb-6-20.dat";

  std::ifstream stream_jj(file_jj.c_str());
  while(std::getline(stream_jj,line))
  {
    std::istringstream(line)>>a>>b>>c>>d>>JJ>>TT>>trme;
    std::tuple<int,int,int,int,int,int> key(a,b,c,d,JJ,TT);
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
      if (fabs(rme)>10e-13) 
        {
          std::tuple<int,int,int,int,int,int> key(a1p,a2p,a1,a2,TwiceValue(J),22);
          double trme=test_map_jj[key];
          if(fabs(rme-trme)>10e-8)
            std::cout<<fmt::format("{:3} {:3} {:3} {:3}   {:3}   22   {:15.13f}  {:15.13f}   {:15.13f}", 
              a1p,a2p,a1,a2,TwiceValue(J),rme,trme, fabs(rme-trme))<<std::endl;
        }
    }
}
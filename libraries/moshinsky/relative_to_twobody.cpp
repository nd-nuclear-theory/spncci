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
#include "moshinsky/relative_cm_xform.h"
#include "u3shell/upcoupling.h"
#include "u3shell/two_body_branching.h"

extern double zero_threshold;

void PrintTwoBodyIndexedRMEU3ST(const u3shell::IndexedTwoBodyTensorRMEsU3ST& two_body_rmes)
{
  for(auto it=two_body_rmes.begin(); it!=two_body_rmes.end(); ++it)
      {
        int kappa0,L0;
        u3shell::TwoBodyUnitTensorLabelsU3ST tb_tensor;
        std::tie(tb_tensor,kappa0,L0)=it->first;
        double rme=it->second;
        if(fabs(rme)>zero_threshold)
          std::cout<<fmt::format("{} {} {}   {}",tb_tensor.Str(), kappa0,L0,rme)<<std::endl;
      }
}
//////////

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  u3::WCoefCache w_cache;
  zero_threshold=1e-6;

  int Nmax=4; 
  int Jmax=4; 
  int J0=0;
      
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensors);
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //// Identity Operator 
  ///////////////////////////////////////////////////////////////////////////////////////////////////  
  u3shell::RelativeRMEsU3ST identity_relative_rmes;

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

  std::map<u3shell::TwoBodyBraketLST,double> two_body_rmes_lst;
  u3shell::BranchTwoBodyNLST(id_two_body_rmes,two_body_rmes_lst);
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
  //     if(fabs(rme)>zero_threshold)
  //     std::cout<<fmt::format("{} {} {} {} {} {} {} || {} {} {} || {} {} {} {} {} {} {}   {}",
  //       N1p,L1p,N2p,L2p,Lp,Sp,Tp,L0,S0,T0,N1,L1,N2,L2,L,S,T,rme)<<std::endl;
  //   }

  std::map<u3shell::TwoBodyBraketJJJT,double> id_two_body_rme_jjjt;
  u3shell::BranchTwoBodyU3STToJJJT(Jmax, J0, id_two_body_rmes, id_two_body_rme_jjjt);

  std::cout<<"Identity JJJT "<<std::endl;
  u3shell::PrintTwoBodyMatrixElementsJJJT(id_two_body_rme_jjjt);
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

  std::map<u3shell::TwoBodyBraketLST,double> j16_two_body_rmes_lst;
  u3shell::BranchTwoBodyNLST(j16_indexed_two_body_rmes,j16_two_body_rmes_lst);

  std::map<u3shell::TwoBodyBraketLSJT,double> j16_two_body_rme_lsjt;
  u3shell::BranchTwoBodyLSJT(Jmax, J0, j16_two_body_rmes_lst,j16_two_body_rme_lsjt);

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
  //   if(fabs(trme)<zero_threshold)
  //     continue;
  //   if(fabs(trme-rme)>zero_threshold)
  //     std::cout<<fmt::format("[{} {} {} {}] {} {} {} {} | |[{} {} {} {}] {} {} {} {} {:12f} {:12f}",
  //       N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,N1,L1,N2,L2,L,S,J,T,test_map[key],rme
  //       )<<std::endl;

  // }
  // stream.close();
 

  std::cout<<"Branching "<<std::endl;
  std::map<u3shell::TwoBodyBraketJJJT,double> two_body_rme_jjjt;
  u3shell::BranchTwoBodyU3STToJJJT(Jmax, J0,j16_indexed_two_body_rmes,two_body_rme_jjjt);

  std::map<std::tuple<int,int,HalfInt>,int> h2_lookup;
  u3shell::H2FormatLookUp(Nmax,h2_lookup);

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
      u3shell::TwoBodyStateLabelsJJJT bra,ket;
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
      if (fabs(rme)>zero_threshold) 
        {
          std::tuple<int,int,int,int,int,int> key(a1p,a2p,a1,a2,TwiceValue(J),22);
          double trme=test_map_jj[key];
          if(fabs(rme-trme)>zero_threshold)
            std::cout<<fmt::format("{:3} {:3} {:3} {:3}   {:3}   22   {:15.13f}  {:15.13f}   {:15.13f}", 
              a1p,a2p,a1,a2,TwiceValue(J),rme,trme, fabs(rme-trme))<<std::endl;
        }
    }
}
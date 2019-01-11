/****************************************************************
  relative_to_twobody.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
****************************************************************/
#include <fstream>

#include "fmt/format.h"
#include "basis/lsjt_operator.h"

#include "am/am.h"
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "am/wigner_gsl.h"
#include "basis/jjjt_scheme.h"
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

void RelativeToTwoBody(
  const std::string interaction_filename, int Nmax, int N1v,
  int Jmax, int J0, int T0, int g0,
  std::map<u3shell::TwoBodyBraketJJJT,double>& two_body_rme_jjjt)
{

  u3::WCoefCache w_cache;
  // std::string interaction_file="/afs/crc.nd.edu/user/a/amccoy/projects/spncci/data/jisp16_Nmax20_hw20.0_rel.dat";
  // std::string interaction_file="../../data/jisp16_Nmax20_hw20.0_rel.dat";
  basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Jmax);

  basis::OperatorLabelsJT operator_labels;
  std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
  std::array<basis::MatrixVector,3> isospin_component_matrices_lsjt;
  basis::ReadRelativeOperatorLSJT(
      interaction_filename,relative_space_lsjt,operator_labels,
      isospin_component_sectors_lsjt,isospin_component_matrices_lsjt, true
    );

  std::cout<<"Upcoupling"<<std::endl;
  u3shell::RelativeRMEsU3ST rme_map;
  u3shell::Upcoupling(
    relative_space_lsjt,
    isospin_component_sectors_lsjt,
    isospin_component_matrices_lsjt,
    w_cache, J0, g0,T0, Nmax,rme_map
    );
  std::cout<<rme_map.size()<<std::endl;
  // for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
  //   {
  //     u3shell::RelativeUnitTensorLabelsU3ST tensor;
  //     int kappa0, L0;
  //     std::tie(tensor,kappa0,L0)=it->first;
  //     std::cout<<tensor.Str()<<" "<<kappa0<<" "<<L0<<"  "<<it->second<<std::endl;
  //   }
  std::cout<<"Convering to Two-body"<<std::endl;
  u3shell::IndexedTwoBodyTensorRMEsU3ST j16_indexed_two_body_rmes;
  u3shell::ConvertRelativeTensorToTwoBodyTensor(Nmax,N1v,rme_map,j16_indexed_two_body_rmes);
  // PrintTwoBodyIndexedRMEU3ST(k2_indexed_two_body_rmes);

  std::cout<<"branch nlst"<<std::endl;
  std::map<u3shell::TwoBodyBraketLST,double> j16_two_body_rmes_lst;
  u3shell::BranchTwoBodyNLST(j16_indexed_two_body_rmes,j16_two_body_rmes_lst);

  std::cout<<"branch lsjt"<<std::endl;
  std::map<u3shell::TwoBodyBraketLSJT,double> j16_two_body_rme_lsjt;
  u3shell::BranchTwoBodyLSJT(Jmax, J0, j16_two_body_rmes_lst,j16_two_body_rme_lsjt);

  std::cout<<"Branching "<<std::endl;
  // std::map<u3shell::TwoBodyBraketJJJT,double> two_body_rme_jjjt;
  u3shell::BranchTwoBodyU3STToJJJT(Jmax, J0,j16_indexed_two_body_rmes,two_body_rme_jjjt);
}


void CompareToMFDN(
  int Nmax,
  const std::map<u3shell::TwoBodyBraketJJJT,double>& two_body_rme_jjjt,
  std::string file_jj
  )
{
  std::map<std::tuple<int,int,HalfInt>,int> h2_lookup;
  u3shell::H2FormatLookUp(Nmax,h2_lookup);

  // Read in matrix elements from file
  int a, b, c, d, JJ, TT;
  double trme;
  std::string line;
  std::map<std::tuple<int,int,int,int,int,int>,double> test_map_jj;
  // std::string file_jj="test/JISP16-tb-6-20.dat";

  std::ifstream stream_jj(file_jj.c_str());

  while(std::getline(stream_jj,line))
  {
    std::istringstream(line)>>a>>b>>c>>d>>JJ>>TT>>trme;
    std::tuple<int,int,int,int,int,int> key(a,b,c,d,JJ,TT);
    test_map_jj[key]=trme;
  }


  std::vector<u3shell::pn_rmes>  two_body_rmes_jjjpn;
  u3shell::ConvertJJJTToPN(Nmax,two_body_rme_jjjt,two_body_rmes_jjjpn);

  for(int i=0; i<3; ++i)
    for(auto it=two_body_rmes_jjjpn[i].begin(); it!=two_body_rmes_jjjpn[i].end(); ++it)
      {
        int J0;
        u3shell::TwoBodyStateLabelsJJJPN bra,ket;
        std::tie(J0,bra,ket)=it->first;
        double rme=it->second;

        int a1,a2,a1p,a2p,Jp,J;
        std::tie(a1,a2,J)=ket;
        std::tie(a1p,a2p,Jp)=bra;

        int type;
        if(i==0)
          type=11;
        if(i==1)
          type=22;
        if(i==2)
          type=12;

        // file format for comparison currently assumes J0=0
        assert(J==Jp);
        std::tuple<int,int,int,int,int,int> key(a1p,a2p,a1,a2,TwiceValue(J),type);
        double trme=test_map_jj[key];
        // if(fabs(rme-trme)>zero_threshold)
          std::cout<<fmt::format("{:3} {:3} {:3} {:3}   {:3}   {:3}   {:15.8f}  {:15.8f}   {:5.6f}",
            a1p,a2p,a1,a2,TwiceValue(J),type,rme,trme, fabs(rme/trme))<<std::endl;

      }

  // for(auto it=two_body_rme_jjjt.begin(); it!=two_body_rme_jjjt.end(); ++it)
  //   {
  //     int T0, J0;
  //     u3shell::TwoBodyStateLabelsJJJT bra,ket;
  //     std::tie(J0,T0,bra,ket)=it->first;
  //     if(bra>ket)
  //       continue;
  //     double rme=it->second;
  //     int N1p,N2p,N1,N2,L1,L2,L1p,L2p;
  //     HalfInt J1,J2,J1p,J2p, Jp,J,Tp,T;
  //     std::tie(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp)=bra;
  //     std::tie(N1,L1,J1,N2,L2,J2,J,T)=ket;

  //     if(Tp==0)
  //       continue;

  //     // threshold set in mfdn two-body rme file
  //     if((N1+N2)>6)
  //       continue;

  //     // threshold set in mfdn two-body rme file
  //     if((N1p+N2p)>6)
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
  //     if (fabs(rme)>zero_threshold)
  //       {
  //         std::tuple<int,int,int,int,int,int> key(a1p,a2p,a1,a2,TwiceValue(J),11);
  //         double trme=test_map_jj[key];
  //         // if(fabs(rme-trme)>zero_threshold)
  //           std::cout<<fmt::format("{:3} {:3} {:3} {:3}   {:3}   11   {:15.8f}  {:15.8f}   {:5.6f}",
  //             a1p,a2p,a1,a2,TwiceValue(J),rme,trme, fabs(rme/trme))<<std::endl;
  //       }
  // }
}


  // ///////////////////////////////////////////////////////////////////////////////////////////////////
  // //// Identity Operator  TODO make into function
  // ///////////////////////////////////////////////////////////////////////////////////////////////////
  // u3shell::RelativeRMEsU3ST identity_relative_rmes;

  // for(auto tensor : relative_unit_tensors)
  //   {
  //     if(
  //       (tensor.x0()==u3::SU3(0,0))
  //       &&(tensor.S0()==0)
  //       &&(tensor.T0()==0)
  //       &&(tensor.bra()==tensor.ket())
  //       )
  //       {
  //         u3shell::IndexedRelativeUnitTensorLabelsU3ST index_tensor(tensor,1,0);
  //         identity_relative_rmes[index_tensor]=1;
  //       }
  //   }
  // u3shell::IndexedTwoBodyTensorRMEsU3ST id_two_body_rmes;
  // u3shell::ConvertRelativeTensorToTwoBodyTensor(Nmax,N1v,identity_relative_rmes, id_two_body_rmes);
  // PrintTwoBodyIndexedRMEU3ST(id_two_body_rmes);

  // std::map<u3shell::TwoBodyBraketLST,double> two_body_rmes_lst;
  // u3shell::BranchTwoBodyNLST(id_two_body_rmes,two_body_rmes_lst);
  // // std::cout<<"Branch NLST"<<std::endl;
  // // for(auto it=two_body_rmes_lst.begin(); it!=two_body_rmes_lst.end(); ++it)
  // //   {
  // //     TwoBodyBraketLST braket=it->first;
  // //     double rme=it->second;
  // //     TwoBodyStateLabelsLST bra,ket;
  // //     int L0;
  // //     HalfInt S0,T0;
  // //     std::tie(L0,S0,T0,bra,ket)=braket;
  // //     int N1,L1,N2,L2,L,N1p,L1p,N2p,L2p,Lp;
  // //     HalfInt Sp,Tp,S,T;
  // //     std::tie(N1,L1,N2,L2,L,S,T)=ket;
  // //     std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Tp)=bra;
  // //     if(fabs(rme)>zero_threshold)
  // //     std::cout<<fmt::format("{} {} {} {} {} {} {} || {} {} {} || {} {} {} {} {} {} {}   {}",
  // //       N1p,L1p,N2p,L2p,Lp,Sp,Tp,L0,S0,T0,N1,L1,N2,L2,L,S,T,rme)<<std::endl;
  // //   }

  // std::map<u3shell::TwoBodyBraketJJJT,double> id_two_body_rme_jjjt;
  // u3shell::BranchTwoBodyU3STToJJJT(Jmax, J0, id_two_body_rmes, id_two_body_rme_jjjt);

  // std::cout<<"Identity JJJT "<<std::endl;
  // u3shell::PrintTwoBodyMatrixElementsJJJT(id_two_body_rme_jjjt);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  // u3::WCoefCache w_cache;
  zero_threshold=1e-6;

  int Nmax=4;
  int N1v=1;
  int Jmax=Nmax+2;
  int J0=2;

  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, N1v,relative_unit_tensors);


  // std::string interaction_filename="../../data/relative_interactions/jisp16_Nmax20_hw20.0_rel.dat";
  // std::string interaction_filename="../../data/relative_interactions/coulomb_test_Nmax20_steps500_rel.dat";
  // std::string interaction_filename="../../data/relative_interactions/nnloopt_Nmax20_hw40.0_caveat-nmax30_rel.dat";
  std::string interaction_filename="../../data/relative_interactions/quadrupole_test_Nmax6_total_rel.dat";

  // need to add sym link of data to nuclty directory rel
  // std::string interaction_filename="../../data/JISP16_Nmax20_hw20.0_rel.dat";
  // std::string interaction_filename="test/coulomb_Nmax20_rel.dat";

  int T0=-1;
  int g0=0;

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //// LSJT
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  u3::WCoefCache w_cache;
  // std::string interaction_file="/afs/crc.nd.edu/user/a/amccoy/projects/spncci/data/jisp16_Nmax20_hw20.0_rel.dat";
  // std::string interaction_file="../../data/jisp16_Nmax20_hw20.0_rel.dat";
  basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Jmax);

  basis::OperatorLabelsJT operator_labels;
  std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
  std::array<basis::MatrixVector,3> isospin_component_matrices_lsjt;
  basis::ReadRelativeOperatorLSJT(
      interaction_filename,relative_space_lsjt,operator_labels,
      isospin_component_sectors_lsjt,isospin_component_matrices_lsjt, true
    );

  std::cout<<"Upcoupling"<<std::endl;
  u3shell::RelativeRMEsU3ST rme_map;
  u3shell::Upcoupling(
    relative_space_lsjt,
    isospin_component_sectors_lsjt,
    isospin_component_matrices_lsjt,
    w_cache, J0, g0,T0, Nmax,rme_map
    );

  std::cout<<rme_map.size()<<std::endl;
  if(false)
    {
      for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
        {
          u3shell::RelativeUnitTensorLabelsU3ST tensor;
          int kappa0, L0;
          std::tie(tensor,kappa0,L0)=it->first;
          std::cout<<tensor.Str()<<" "<<kappa0<<" "<<L0<<"  "<<it->second<<std::endl;
        }
    }

  std::cout<<"Convering to Two-body"<<std::endl;
  u3shell::IndexedTwoBodyTensorRMEsU3ST j16_indexed_two_body_rmes;
  u3shell::ConvertRelativeTensorToTwoBodyTensor(Nmax,N1v,rme_map,j16_indexed_two_body_rmes);

  if(false)
    {
      PrintTwoBodyIndexedRMEU3ST(j16_indexed_two_body_rmes);
    }

  std::cout<<"branch nlst"<<std::endl;
  std::map<u3shell::TwoBodyBraketLST,double> j16_two_body_rmes_lst;
  u3shell::BranchTwoBodyNLST(j16_indexed_two_body_rmes,j16_two_body_rmes_lst);


  std::cout<<"branch lsjt"<<std::endl;
  std::map<u3shell::TwoBodyBraketLSJT,double> j16_two_body_rme_lsjt;
  u3shell::BranchTwoBodyLSJT(Jmax, J0, j16_two_body_rmes_lst,j16_two_body_rme_lsjt);

  if(true)
    {
      for(auto it=j16_two_body_rme_lsjt.begin(); it!=j16_two_body_rme_lsjt.end(); ++it)
        {
          int T0, J0;
          u3shell::TwoBodyStateLabelsLSJT bra,ket;
          std::tie(J0,T0,bra,ket)=it->first;
          // if(bra>ket)
          //   continue;
          double rme=it->second;
          int N1p,N2p,N1,N2,L1,L2,L1p,L2p,Lp,L;
          HalfInt Sp,S,Jp,J,Tp,T;
          std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp)=bra;
          std::tie(N1,L1,N2,L2,L,S,J,T)=ket;

          // T0  N1' l1' N2' l2' L' S' J' T' g'  N1 l1 N2 l2 L S J T g  JT-RME
          std::cout<<fmt::format("{}   {} {}  {} {}  {} {} {} {} {}   {} {}  {} {}  {} {} {} {} {}   {:13.6f} ",
            T0,
            N1p,L1p, N2p,L2p,Lp,Sp,Jp,Tp,(N1p+N2p)%2,
            N1,L1,N2,L2,L,S,J,T,(N1+N2)%2,
            rme
            )<<std::endl;
        }
    }

  // std::cout<<"Branching "<<std::endl;
  std::map<u3shell::TwoBodyBraketJJJT,double> two_body_rme_jjjt;
  u3shell::BranchTwoBodyU3STToJJJT(Jmax, J0,j16_indexed_two_body_rmes,two_body_rme_jjjt);

  if(false)
    {

      for(auto it=two_body_rme_jjjt.begin(); it!=two_body_rme_jjjt.end(); ++it)
        {
          int T0, J0;
          u3shell::TwoBodyStateLabelsJJJT bra,ket;
          std::tie(J0,T0,bra,ket)=it->first;
          // if(bra>ket)
          //   continue;
          double rme=it->second;
          int N1p,N2p,N1,N2,L1,L2,L1p,L2p;
          HalfInt J1,J2,J1p,J2p, Jp,J,Tp,T;
          std::tie(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp)=bra;
          std::tie(N1,L1,J1,N2,L2,J2,J,T)=ket;

          // T0  N1' l1' j1' N2' l2' J' T' g'  N1 l1 j1 N2 l2 j2 J T g  JT-RME
          std::cout<<fmt::format("{}   {}  {}  {:3.1f}   {}  {}  {:3.1f}   {}  {}  {}    {}  {}  {:3.1f}   {}  {}  {:3.1f}   {}  {}  {}    {:13.6f} ",
            T0,
            N1p,L1p,float(J1p), N2p,L2p,float(J2p), float(Jp),float(Tp),(N1p+N2p)%2,
            N1,L1,float(J1), N2,L2,float(J2),J,T,(N1+N2)%2,
            rme
            )<<std::endl;
        }
    }









//   // transforming from relative to two-body jjjt branched matrix elements
//   std::map<u3shell::TwoBodyBraketJJJT,double> two_body_rme_jjjt;
//   RelativeToTwoBody(
//     interaction_filename, Nmax, N1v,
//     Jmax, J0, T0, g0,two_body_rme_jjjt);

//   std::map<std::tuple<int,int,HalfInt>,int> h2_lookup;
//   u3shell::H2FormatLookUp(Nmax,h2_lookup);



//   for(auto it=two_body_rme_jjjt.begin(); it!=two_body_rme_jjjt.end(); ++it)
//     {
//       int T0, J0;
//       u3shell::TwoBodyStateLabelsJJJT bra,ket;
//       std::tie(J0,T0,bra,ket)=it->first;
//       if(bra>ket)
//         continue;
//       double rme=it->second;
//       int N1p,N2p,N1,N2,L1,L2,L1p,L2p;
//       HalfInt J1,J2,J1p,J2p, Jp,J,Tp,T;
//       std::tie(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp)=bra;
//       std::tie(N1,L1,J1,N2,L2,J2,J,T)=ket;

// // T0  N1' l1' j1' N2' l2' J' T' g'  N1 l1 j1 N2 l2 j2 J T g  JT-RME
//       std::cout<<fmt::format("{}   {} {} {}  {} {} {}  {} {} {}   {} {} {}  {} {} {}  {} {} {}   {:13.6f} ",
//         T0,
//         N1p,L1p,J1p, N2p,L2p,J2p, Jp,Tp,(N1p+N2p)%2,
//         N1,L1,J1, N2,L2,J2,J,T,(N1+N2)%2,
//         rme
//         )<<std::endl;
//     }



  // // std::string file_jj="test/JISP16-tb-6-20.dat";
  // std::string file_jj="test/VC-tb-6-20.dat";
  // CompareToMFDN(Nmax,two_body_rme_jjjt,file_jj);


}

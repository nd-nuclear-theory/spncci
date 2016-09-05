/****************************************************************
  interaction_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
****************************************************************/
#include "cppformat/format.h"

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
#include "u3shell/u3st_scheme.h"

typedef std::tuple<int,int,u3::SU3,int,int,HalfInt,HalfInt> TwoBodyStateLabelsU3STBranched;
typedef std::tuple<u3::SU3,int,int,HalfInt,int,TwoBodyStateLabelsU3STBranched,TwoBodyStateLabelsU3STBranched>TwoBodyRMELablesU3STBranched;

typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt,HalfInt> TwoBodyStateLabelsLSJT;
// eta1,L1, eta2,L2,L,S,J,T
typedef std::pair<TwoBodyStateLabelsLSJT,TwoBodyStateLabelsLSJT>TwoBodyRMELablesLSJT;


int main(int argc, char **argv)
{
  u3::U3CoefInit();
  int Nmax=2;
  int J0=0;
  int g0=0;
	int T0=0;
  int Jmax=4;
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);
  basis::RelativeSectorsLSJT relative_sectors(relative_space, J0, T0, g0);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // JISP16
	// std::string interaction_file="data/Vrel_JISP16_bare_Jmax4.hw20";
  // std::vector<Eigen::MatrixXd> sector_vector=u3shell::ImportInteraction(interaction_file, space, sectors, "JISP");
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Kinetic eneryg test 
  // std::string interaction_file="NONE";
  // std::vector<Eigen::MatrixXd> sector_vector=u3shell::ImportInteraction(interaction_file, relative_space, relative_sectors, "Kinetic");
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Identity test 
  std::string interaction_file="NONE";
  std::vector<Eigen::MatrixXd> sector_vector=u3shell::ImportInteraction(interaction_file, relative_space, relative_sectors, "Identity");
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  u3shell::UpcouplingNLST(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);
  
  // std::cout<<"UpcouplingNLST"<<std::endl;
  for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
    {
      int L0, S0, L,S,T, Lp, Sp, Tp;
      u3shell::RelativeSubspaceLabelsNLST bra, ket;
      std::tie(L0,S0,bra,ket)=it->first;
      std::tie(L,S,T)=ket;
      std::tie(Lp,Sp,Tp)=bra;
      // std::cout<<fmt::format("{} {} ({},{},{}) ({},{},{})", L0,S0,Lp,Sp,Tp,L,S,T)<<std::endl<<it->second<<std::endl;
    }
  std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>,double> relative_rme_map;

  u3shell::Upcoupling(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);
  // std::cout<<"UpcouplingU3ST"<<std::endl;
  for (auto it=relative_rme_map.begin(); it!= relative_rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST labels;
      int kappa0, L0;
      std::tie(labels,kappa0, L0)=it->first;
      double coef=it->second;

      // if (fabs(coef)>10e-13)
      //   std::cout<<labels.Str()<<"  "<<kappa0<<"  "<<L0<<std::endl<<it->second<<std::endl<<std::endl;
    }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Generate a RME expansion for all relative unit tensors and store for easy look-up.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  u3shell::TwoBodySpaceU3ST two_body_space(Nmax);
  std::map<u3shell::RelativeUnitTensorLabelsU3ST,u3shell::TwoBodyUnitTensorCoefficientsU3ST> relative_two_body_expansion_map;
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensor_labels);
  double expansion_coef=1;
  for (int i=0; i<relative_unit_tensor_labels.size(); ++i)
    {
      u3shell::RelativeUnitTensorLabelsU3ST tensor(relative_unit_tensor_labels[i]);
      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion;
      u3shell::MoshinskyTransformUnitTensor(tensor,expansion_coef,two_body_space,two_body_expansion);
      relative_two_body_expansion_map[tensor]=two_body_expansion;
    }
  // for(auto it=relative_two_body_expansion_map.begin(); it!=relative_two_body_expansion_map.end(); ++it )
  // {
  //   std::cout<<(it->first).Str()<<std::endl;;
  //   for(auto it2=(it->second).begin(); it2!=(it->second).end(); it2++)
  //     std::cout<<(it2->first).Str()<<"  "<<it2->second<<std::endl;
  // }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Combining unit tensors and relative rme's of operator
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::map<std::pair<TwoBodyStateLabelsLSJT,TwoBodyStateLabelsLSJT>,double> branched_rme;
  // iterate through rme's
  std::map<TwoBodyRMELablesU3STBranched,double> two_body_rme_u3st_branched;

  /// For testing 
  std::map<u3shell::TwoBodyUnitTensorLabelsU3ST,double> test_map;

  for(auto rit=relative_rme_map.begin(); rit!=relative_rme_map.end(); rit++)
    {
      if(fabs(rit->second)<=10e-14)
        continue;
      u3shell::RelativeUnitTensorLabelsU3ST relative_tensor;
///
///      
      int kappa0,L0,g;
      std::tie(relative_tensor,kappa0, L0)=rit->first;
      // std::cout<<relative_tensor.Str()<<std::endl;

      double relative_rme=rit->second;
      const u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_tensors=relative_two_body_expansion_map[relative_tensor];
      u3::SU3 x0(relative_tensor.x0());
      HalfInt S0=relative_tensor.S0();
      // MultiplicityTagged<int>::vector L0_kappa0_max=BranchingSO3(x0);

      for(auto tbit=two_body_tensors.begin(); tbit!=two_body_tensors.end(); ++tbit)
        {
          u3shell::TwoBodyUnitTensorLabelsU3ST two_body_tensor=tbit->first;
          
          // std::cout<<two_body_tensor.Str()<<std::endl;
          
          double unit_tensor_rme=tbit->second;
          // std::cout<<"unit_tensor "<<two_body_tensor.Str()<<"  "<<unit_tensor_rme<<std::endl;
          int eta1,eta2, eta1p,eta2p;
          HalfInt S,T, Sp, Tp;
          u3::SU3 x, xp;
          std::tie(eta1,eta2,x,S,T)=two_body_tensor.ket().Key();
          std::tie(eta1p,eta2p,xp,Sp,Tp)=two_body_tensor.bra().Key();
          int rho0=two_body_tensor.rho0();
          MultiplicityTagged<int>::vector L_kappa_max=BranchingSO3(x);
          MultiplicityTagged<int>::vector Lp_kappap_max=BranchingSO3(xp);
          test_map[two_body_tensor]+=relative_rme*unit_tensor_rme;
          if((L0>(S0+J0))||L0<abs(S0-L0))
            continue;
          for(int lp=0; lp<Lp_kappap_max.size(); ++lp)
            {
              int Lp=Lp_kappap_max[lp].irrep;
              int kappap_max=Lp_kappap_max[lp].tag;
              for(int l=0; l<L_kappa_max.size(); ++l)
                {
                  int L=L_kappa_max[l].irrep;
                  int kappa_max=L_kappa_max[l].tag;
                  for(int kappap=1; kappap<=kappap_max; ++kappap)
                    for(int kappa=1; kappa<=kappa_max; ++kappa)
                      {
                        double u3coef1=u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0);
                        // std::cout<<u3coef1<<"  "<<relative_rme<<"  "<<unit_tensor_rme<<std::endl;
                        TwoBodyStateLabelsU3STBranched ket(eta1,eta2,x,kappa,L,S,T);
                        TwoBodyStateLabelsU3STBranched bra(eta1p,eta2p,xp,kappap,Lp,Sp,Tp);
                        TwoBodyRMELablesU3STBranched key(x0,kappa0,L0,S0,T0,bra,ket);
                        // testing 
                        TwoBodyStateLabelsU3STBranched ket1(1,1,u3::SU3(2,0),1,0,0,1);
                        TwoBodyRMELablesU3STBranched key1(u3::SU3(0,0),1,0,0,0,ket1,ket1);
                        // if(key==key1)
                          // std::cout<<relative_rme<<"  "<<unit_tensor_rme<<"  "<<u3coef1<<std::endl;
                        // 
                        two_body_rme_u3st_branched[key]+=relative_rme*unit_tensor_rme*u3coef1;
                      }
                }
            }
        }
      }
      //for(auto it=test_map.begin(); it!= test_map.end(); ++it)
      //  std::cout<<it->first.Str()<<"  "<<it->second<<std::endl;
      for(auto it=two_body_rme_u3st_branched.begin(); it!= two_body_rme_u3st_branched.end(); ++it)
        {
          int kappa0,L0,T0,kappa,L,kappap,Lp,eta1,eta2,eta1p,eta2p;
          HalfInt Sp,S,Tp,T,S0;
          TwoBodyStateLabelsU3STBranched bra,ket;
          u3::SU3 x0,x,xp;
          std::tie(x0,kappa0,L0,S0,T0,bra,ket)=it->first;
          std::tie(eta1,eta2,x,kappa,L,S,T)=ket;
          std::tie(eta1p,eta2p,xp,kappap,Lp,Sp,Tp)=bra;
          HalfInt::pair Jp_range=am::ProductAngularMomentumRange(Lp, Sp);
          HalfInt::pair J_range=am::ProductAngularMomentumRange(L, S);
          double rme_u3st_branched=it->second;
          // std::cout<<fmt::format("[{} {} {} {} {}]({} {} {} {} {} {} {}, {} {} {} {} {} {} {}):  {}", x0.Str(),
          //   kappa0,L0,S0,T0,eta1p,eta2p,xp.Str(),kappap,Lp,Sp,Tp,eta1,eta2,x.Str(),kappa,L,S,T,it->second)
          // <<std::endl;

          MultiplicityTagged<int>::vector L1_kappa1_max=BranchingSO3(u3::SU3(eta1,0));
          MultiplicityTagged<int>::vector L2_kappa2_max=BranchingSO3(u3::SU3(eta2,0));
          MultiplicityTagged<int>::vector L1p_kappa1p_max=BranchingSO3(u3::SU3(eta1p,0));
          MultiplicityTagged<int>::vector L2p_kappa2p_max=BranchingSO3(u3::SU3(eta2p,0));
                            
          for(HalfInt J=J_range.first; J<=J_range.second; ++J)
            for(HalfInt Jp=Jp_range.first; Jp<=Jp_range.second; ++Jp)
              {
                double so3_coef=am::Unitary9J(L,S,J, L0,S0,J0, Lp,Sp,Jp);

                for(int l1p=0; l1p<L1p_kappa1p_max.size(); ++l1p)
                  {
                    int L1p(L1p_kappa1p_max[l1p].irrep);
                    for(int l2p=0; l2p<L2p_kappa2p_max.size(); ++l2p)
                      {
                        int L2p(L2p_kappa2p_max[l2p].irrep);
                        if((Lp<abs(L1p-L2p))||(Lp>(L1p+L2p)))
                          continue;
                        for(int l1=0;l1<L1_kappa1_max.size(); ++l1)
                          {
                            int L1(L1_kappa1_max[l1].irrep);
                            for(int l2=0; l2<L2_kappa2_max.size(); ++l2)
                              {
                                int L2(L2_kappa2_max[l2].irrep);
                                if((L<abs(L1-L2))||(L>(L1+L2)))
                                  continue;
                                double u3coef2=u3::W(u3::SU3(eta1,0),1,L1,u3::SU3(eta2,0),1,L2,x,kappa,L,1)
                                          *u3::W(u3::SU3(eta1p,0),1,L1p,u3::SU3(eta2p,0),1,L2p,xp,kappap,Lp,1)
                                          *parity((eta1-L1)/2+(eta2-L2)/2+(eta1p-L1p)/2+(eta2p-L2p)/2);
                                double rme=rme_u3st_branched*u3coef2*so3_coef;
                                TwoBodyStateLabelsLSJT state=std::make_tuple(eta1,L1,eta2,L2,L,S,J,T);
                                TwoBodyStateLabelsLSJT statep=std::make_tuple(eta1p,L1p,eta2p,L2p,Lp,Sp,Jp,Tp);
                                TwoBodyRMELablesLSJT key=std::make_pair(statep,state);
                                branched_rme[key]+=rme;
                              }
                          }

                      } 
                  }
              }
          
      }
  for(auto bit=branched_rme.begin(); bit!=branched_rme.end(); ++bit)
    {
       
      if(fabs(bit->second)>10e-14)
        std::cout<<bit->second<<std::endl;
    }
}




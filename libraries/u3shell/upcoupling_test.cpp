/****************************************************************
  upcoupling_test.cpp

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
#include "u3shell/u3st_scheme.h"

void IdentityTest(
  int Nmax, int Jmax, int J0, int T0, int g0,
  u3shell::RelativeRMEsU3ST& relative_rme_map
 )
{
  std::cout<<"Identity test"<<std::endl;

  basis::RelativeSpaceLSJT relative_space(Nmax, Jmax);
  basis::RelativeSectorsLSJT relative_sectors(relative_space, J0,T0, g0);
  std::vector<Eigen::MatrixXd> sector_vector;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;


  std::string interaction_file="NONE";
  sector_vector
    =u3shell::ImportInteraction(interaction_file, relative_space, relative_sectors, "Identity");
  //////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"UpcouplingNLST"<<std::endl;
  u3shell::UpcouplingNLST(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);
  for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
    {
      int L0, S0, L,S,T, Lp, Sp, Tp;
      u3shell::RelativeSubspaceLabelsNLST bra, ket;
      std::tie(L0,S0,bra,ket)=it->first;
      std::tie(L,S,T)=ket;
      std::tie(Lp,Sp,Tp)=bra;
      Eigen::MatrixXd sectorNLST=it->second;
      if(fabs(sectorNLST.sum())>10e-10)
        std::cout<<fmt::format("{} {} ({},{},{}) ({},{},{})", L0,S0,Lp,Sp,Tp,L,S,T)
        <<std::endl<<sectorNLST<<std::endl;
    }

  u3shell::Upcoupling(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);
  std::cout<<"UpcouplingU3ST"<<std::endl;
  for (auto it=relative_rme_map.begin(); it!= relative_rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST labels;
      int kappa0, L0;
      std::tie(labels,kappa0,L0)=it->first;
      double coef=it->second;
      if (fabs(coef)>10e-13)
        std::cout<<labels.Str()<<"  "<<kappa0<<"  "<<L0<<std::endl<<it->second<<std::endl<<std::endl;
    }
}

void
KineticCheck(
  int Nmax, int Jmax, int J0, int T0, int g0,
  u3shell::RelativeRMEsU3ST& relative_rme_map
  )
{
  std::cout<<"Kinetic Energy"<<std::endl;
  basis::RelativeSpaceLSJT relative_space(Nmax, Jmax);
  basis::RelativeSectorsLSJT relative_sectors(relative_space, J0,T0, g0);
  std::vector<Eigen::MatrixXd> sector_vector;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;

  std::string interaction_file="NONE";
  sector_vector=u3shell::ImportInteraction(interaction_file, relative_space, relative_sectors, "Kinetic");
  u3shell::UpcouplingNLST(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);

  std::cout<<"UpcouplingNLST"<<std::endl;
  for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
    {
      int L0, S0, L,S,T, Lp, Sp, Tp;
      u3shell::RelativeSubspaceLabelsNLST bra, ket;
      std::tie(L0,S0,bra,ket)=it->first;
      std::tie(L,S,T)=ket;
      std::tie(Lp,Sp,Tp)=bra;
      Eigen::MatrixXd sectorNLST=it->second;
      if(fabs(sectorNLST.sum())>10e-8)
        std::cout<<fmt::format("{} {} ({},{},{}) ({},{},{})", L0,S0,Lp,Sp,Tp,L,S,T)<<std::endl<<it->second<<std::endl;
    }

  std::cout<<"UpcouplingU3ST"<<std::endl;
  u3shell::Upcoupling(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);

  for (auto it=relative_rme_map.begin(); it!=relative_rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST labels;
      int kappa0, L0;
      std::tie(labels,kappa0,L0)=it->first;
      u3::SU3(labels.x0());
      u3shell::RelativeStateLabelsU3ST kett(labels.ket());
      u3shell::RelativeStateLabelsU3ST brat(labels.bra());
      double coefout=it->second;
      // double RelativeKineticEnergyOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)

      if (fabs(coefout)>10e-8)
      {
        double Trme=RelativeKineticEnergyOperator(brat,kett);
        std::cout<<labels.Str()
        <<"  "<<kappa0<<"  "<<L0
        <<std::endl
        <<it->second
        <<"  "
        <<Trme
        <<std::endl;     
      }

    }
}
void
JISPCheck(
  int Nmax, int Jmax, int J0, int T0, int g0,
  std::string interaction_file,
  std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>,double>& relative_rme_map

)
{  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  basis::OperatorLabelsJT operator_labels;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  
  basis::ReadRelativeOperatorLSJT(
    interaction_file,relative_lsjt_space,operator_labels,
    relative_component_sectors,relative_component_matrices, true
    );

  // for(int i=0; i<relative_component_matrices[0].size(); ++i)
  //   {
  //     basis::RelativeStateLSJT bra(relative_component_sectors[0].GetSector(i).bra_subspace(),0);
  //     basis::RelativeStateLSJT ket(relative_component_sectors[0].GetSector(i).ket_subspace(),0);
  //     std::cout<<fmt::format("{} {} {} {} {}", bra.L(), ket.L(), ket.S(), ket.J(), ket.T())<<std::endl;
  //     std::cout<<relative_component_matrices[0][i]<<std::endl;
  //  }

  const basis::MatrixVector& sector_vector=relative_component_matrices[0];
  const basis::RelativeSectorsLSJT& relative_lsjt_sectors=relative_component_sectors[0];
  u3shell::Upcoupling(relative_lsjt_space,relative_lsjt_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);
  std::cout<<"upcoupling complete"<<std::endl;
}

typedef std::tuple<int,int,u3::SU3,HalfInt, HalfInt> RelativeCMU3STLabels;

typedef std::tuple<u3::SU3,HalfInt, HalfInt,int, int,RelativeCMU3STLabels, RelativeCMU3STLabels,int> RelativeCMU3STBraket;
typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> RelativeCMLSTLabels;
typedef std::tuple<int,HalfInt,HalfInt,RelativeCMLSTLabels,RelativeCMLSTLabels> RelativeCMLSTBraket;

void
RelativeToCMU3ST(int Nmax,  
  const std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>,double>& relative_rme_map,
  std::map<RelativeCMU3STBraket,double>& relative_cm_u3st_map)
{
  int N0, etap,eta,eta_cm, kappa0,L0;
  HalfInt S0,T0,Sp,Tp,S,T;
  u3::SU3 x0;
  u3shell::RelativeUnitTensorLabelsU3ST tensor;
  for(auto it=relative_rme_map.begin(); it!=relative_rme_map.end(); ++it)
    {
      std::tie(tensor,kappa0,L0)=it->first;
      std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=tensor.FlatKey();
      if (etap>Nmax)
        continue;
      if(eta>Nmax)
        continue;
      u3::SU3 xr(eta,0);
      u3::SU3 xrp(etap,0);
      double rme=it->second;
      std::cout<<fmt::format("{} {} {}  {}",tensor.Str(),kappa0,L0,rme)<<std::endl;

      // for(eta_cm=0; eta_cm<=Nmax;eta_cm++)
      for(eta_cm=0; eta_cm<=2;eta_cm++)
        {
          u3::SU3 x_cm(eta_cm,0);
          MultiplicityTagged<u3::SU3>::vector x_set=u3::KroneckerProduct(xr,x_cm);
          MultiplicityTagged<u3::SU3>::vector xp_set=u3::KroneckerProduct(xrp,x_cm);
          for(auto ip : xp_set)
            for(auto i: x_set)
              {
                u3::SU3 x(i.irrep);
                u3::SU3 xp(ip.irrep);
                int rho0_max=u3::OuterMultiplicity(x,x0,xp);
                RelativeCMU3STLabels bra(etap,eta_cm,xp,Sp,Tp);
                RelativeCMU3STLabels ket(eta, eta_cm,x,S,T);
                for(int rho0=1; rho0<=rho0_max; ++rho0)
                {
                  RelativeCMU3STBraket braket(x0,S0,T0,kappa0,L0,bra,ket,rho0);
                  relative_cm_u3st_map[braket]=u3::U(x0,xr,xp,x_cm,xrp,1,1,x,1,rho0)*rme;
                                    // std::cout
                  // <<fmt::format("{} {} {}",x0.Str(),u3::SU3(0,0).Str(),x0.Str())<<std::endl
                  // <<fmt::format("{} {} {}",xr.Str(),x_cm.Str(),x.Str())<<std::endl
                  // <<fmt::format("{} {} {}",xrp.Str(),x_cm.Str(),xp.Str())<<std::endl
                  // <<u3::U(x0,xr,xp,x_cm,xrp,1,1,x,1,rho0)<<std::endl;
                  // std::cout<<fmt::format("U({} {} {} {} {} {} {}",x0.Str(),xr.Str(),xp.Str(),x_cm.Str(),xrp.Str(),x.Str(),rho0)
                  //   <<std::endl;
                }
              }
        }
    }
    std::cout<<"Finished upgrading to relative-CM"<<std::endl;
}
void
CMBranchLST(int Nmax,  
  const std::map<RelativeCMU3STBraket,double>& relative_cm_u3st_map,
  std::map<RelativeCMLSTBraket,double>& relative_cm_lst_map
  )
{
  int N0, etap,eta,eta_cm, kappa0, L0;
  HalfInt S0,T0,Sp,Tp,S,T;
  u3::SU3 x0;
  u3shell::RelativeUnitTensorLabelsU3ST tensor;
  for(auto it=relative_cm_u3st_map.begin(); it!=relative_cm_u3st_map.end(); ++it)
    {
      RelativeCMU3STLabels bra,ket;
      int rho0;
      u3::SU3 x,xp;
      std::tie(x0,S0,T0,kappa0,L0,bra,ket,rho0)=it->first;
      double rme=it->second;
      std::tie(etap,eta_cm,xp,Sp,Tp)=bra;
      std::tie(eta, eta_cm,x,S,T)=ket;
      
      // std::cout<<fmt::format("{} {} {} {} {}  {}  {} {} {} {} {}  {} {} {} {} {}  {}",
      //   x0.Str(),kappa0,L0,S0,T0,rho0,etap,eta_cm,xp.Str(),Sp,Tp,eta, eta_cm,x.Str(),S,T,rme)<<std::endl;

      u3::SU3 xr(eta,0);
      u3::SU3 xrp(etap,0);
      u3::SU3 x_cm(eta_cm,0);  
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
              for(int kappap=1; kappap<=kappap_max; ++kappap)
                for(int kappa=1; kappa<=kappa_max; ++kappa)
                  for(int Lrp=etap%2; Lrp<=etap; Lrp+=2)
                    for(int Lr=eta%2; Lr<=eta; Lr+=2)  
                      for(int L_cm=eta_cm%2; L_cm<=eta_cm; L_cm+=2)
                        {     
                          if((abs(Lr-L_cm)>L)||((Lr+L_cm)<L))
                            continue;
                          if((abs(Lrp-L_cm)>Lp)||((Lrp+L_cm)<Lp))
                            continue;
                          RelativeCMLSTLabels bra(etap,Lrp,eta_cm,L_cm,Lp,Sp,Tp);
                          RelativeCMLSTLabels ket(eta,Lr,eta_cm,L_cm,L,S,T);
                          if (bra>ket)
                            continue;
                          RelativeCMLSTBraket braket(L0,S0,T0,bra,ket);
                          int n=(eta-Lr)/2, np=(etap-Lrp)/2;
                          double rme_lst=rme*u3::W(xrp,1,Lrp,x_cm,1,L_cm,xp,kappap,Lp,1)
                                          *u3::W(xr,1,Lr,x_cm,1,L_cm,x,kappa,L,1)
                                          *u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0)
                                          *parity(n+np);
                          double c1=u3::W(xrp,1,Lrp,x_cm,1,L_cm,xp,kappap,Lp,1);
                          double c2=u3::W(xr,1,Lr,x_cm,1,L_cm,x,kappa,L,1);
                          double c3=u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0);
                          int c4=parity(n+np);
                          // std::cout<<fmt::format("{} {} {} {}",c1,c2,c3,c4)<<std::endl;
                          relative_cm_lst_map[braket]+=rme_lst;
                          // std::cout<<fmt::format("{} {} {}  {} {} {}  {}  {}  {}", Lrp,L_cm,Lp,Lr,L_cm,L,rme, rme_lst,relative_cm_lst_map[braket])
                          // <<std::endl;
                        }
            }
        }
    }
}

typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt, HalfInt> RelativeCMLSJTLabels;
typedef std::tuple<int,HalfInt,RelativeCMLSJTLabels,RelativeCMLSJTLabels> RelativeCMLSJTBraket;

void
CMBranchLSJT(
    int Nmax,int J0, const std::map<RelativeCMLSTBraket,double>& relative_cm_lst_map,
    std::map<RelativeCMLSJTBraket,double>& relative_cm_lsjt_map
 )
{
  int eta,etap,eta_cm, Lr,L_cm,Lrp,L,Lp,L0;
  HalfInt Sp,Tp,S0,T0,S,T;
  RelativeCMLSTLabels bra,ket;
  for(auto it=relative_cm_lst_map.begin(); it!=relative_cm_lst_map.end(); ++it)
    {
      if(fabs(it->second)>10e-10)
        {
          std::tie(L0,S0,T0,bra,ket)=it->first;
          std::tie(etap,Lrp,eta_cm,L_cm,Lp,Sp,Tp)=bra;
          std::tie(eta,Lr,eta_cm,L_cm,L,S,T)=ket;
          for(HalfInt Jp=abs(Lp-Sp); Jp<=(Lp+Sp); ++Jp)
            for(HalfInt J=abs(L-S); J<=(L+S); ++J)
              {
                double coef=am::Unitary9J(L,S,J,L0,S0,J0,Lp,Sp,Jp);
                RelativeCMLSJTLabels braj(etap,Lrp,eta_cm,L_cm,Lp,Sp,Jp,Tp);
                RelativeCMLSJTLabels ketj(eta,Lr,eta_cm,L_cm,L,S,J,T);
                RelativeCMLSJTBraket braketj(J0,T0,braj,ketj);
                relative_cm_lsjt_map[braketj]+=coef*(it->second);
              }
        }
    }

}

void
RelativeToCMLST(
  int Nmax,int T0,
  const std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map,
  std::map<RelativeCMLSTBraket,double>& rel_cm_lst_map)
{
    // Convert to LST CM space
  u3shell::RelativeSubspaceLabelsNLST bra_nlst,ket_nlst;
  int Lr,Lrp,L0;
  HalfInt S,T,Sp,Tp,S0;
  for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
    {
      std::tie(L0,S0,bra_nlst,ket_nlst)=it->first;
      std::tie(Lr,S,T)=ket_nlst;
      std::tie(Lrp,Sp,Tp)=bra_nlst;
      const Eigen::MatrixXd& sector(it->second);
      int nmax=sector.cols()-1;

      for (int np=0; np<=nmax; ++np)
        {
          int Nrp=2*np+Lrp;
          if(Nrp>Nmax)
            continue;
          // u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
          for (int n=0; n<=nmax; ++n)
            {
              int Nr=2*n+Lr;
              if(Nr>Nmax)
                continue;
              // u3shell::RelativeStateLabelsU3ST ket(N,S,T);
              //Extract rme
              double rme_nlst=sector(np,n);
              if(fabs(rme_nlst)<10e-10)
                continue;
              std::cout<<fmt::format("{} {}  {} {} {} {}  {} {} {} {}  {}",L0,S0,Nrp,Lrp,Sp,Tp,Nr,Lr,S,T,rme_nlst)<<std::endl;
              for(int Ncm=0; Ncm<=Nmax; ++Ncm)
                for(int Lcm=Ncm%2; Lcm<=Ncm; Lcm+=2)
                  {
                    if((Ncm+Nr)>Nmax)
                      continue;
                    if((Ncm+Nrp)>Nmax)
                      continue;
                    for(int L=abs(Lr-Lcm); L<=(Lr+Lcm); ++L)
                      for(int Lp=abs(Lrp-Lcm); Lp<=(Lrp+Lcm); ++Lp)
                        {
                          RelativeCMLSTLabels bra_cm(Nrp,Lrp,Ncm,Lcm,Lp,Sp,Tp);
                          RelativeCMLSTLabels ket_cm(Nr,Lr,Ncm,Lcm,L,S,T);
                          RelativeCMLSTBraket braket_cm(L0,S0,T0,bra_cm,ket_cm);
                          rel_cm_lst_map[braket_cm]=am::Unitary9J(L0,0,L0,Lr,Lcm,L,Lrp,Lcm,Lp)*rme_nlst;
                          std::cout<<fmt::format("   {} {}  {} {}  {}", Ncm, Lcm, Lp,L,rel_cm_lst_map[braket_cm])<<std::endl;
                        }
                  }
            }
        }
    }
  // Branch and compare.
  std::map<RelativeCMLSJTBraket,double> rel_cm_lsjt_map;
  RelativeCMLSTLabels bra_cm_lst, ket_cm_lst;
  int Lcm, Ncm, L, Lp, Nr,Nrp;
  HalfInt T00; 

}



void
ReadWriteCheck(
    u3shell::RelativeRMEsU3ST& relative_rme_map,
    std::string filename
)
{
   std::ostringstream os;
  WriteRelativeOperatorU3ST(os,relative_rme_map);
  std::ofstream stream;
  stream.open(filename.c_str());
  stream<<os.str();
  // std::cout<<os.str();
  stream.close();

  std::ifstream is(filename.c_str());
  if(!is)
    std::cout<<"Didn't open"<<std::endl;
  u3shell::RelativeRMEsU3ST k2_relative_rmes;
  u3shell::ReadRelativeOperatorU3ST(is, k2_relative_rmes);
}

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  int Nmax=4;
  int Jmax=Nmax+2;
  int J0=0;
  int g0=0;
	int T0=0;

  u3shell::RelativeRMEsU3ST id_relative_rme_map;
  IdentityTest(Nmax,Jmax,J0,T0, g0, id_relative_rme_map);

  // u3shell::RelativeRMEsU3ST ke_relative_rme_map;
  // KineticCheck(Nmax,Jmax,J0,T0,g0,ke_relative_rme_map);
  // std::string filename="Trel_upcouled";
  // ReadWriteCheck(ke_relative_rme_map,filename);

  // std::string interaction_file="/Users/annamccoy/projects/spncci/data/jisp16_Nmax20_hw20.0_rel.dat";
  std::string interaction_file="/Users/annamccoy/projects/spncci/libraries/u3shell/test/symmunit_Nmax04_rel.dat";

  std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>,double> relative_rme_map; 
  JISPCheck(Nmax, 3, J0, T0, g0,interaction_file, relative_rme_map);
  std::map<RelativeCMLSTBraket,double> relative_cm_lst_map;

  // CMBranchLST(2, relative_rme_map,relative_cm_lst_map);
  // int eta,etap,eta_cm, Lr,L_cm,Lrp,L,Lp,L0;
  // HalfInt Sp,Tp,S0,T00,S,T;
  // RelativeCMLSTLabels bra,ket;
  // for(auto it=relative_cm_lst_map.begin(); it!=relative_cm_lst_map.end(); ++it)
  //   {
  //     if(fabs(it->second)>10e-10)
  //       {
  //         std::tie(L0,S0,T00,bra,ket)=it->first;
  //         std::tie(etap,Lrp,eta_cm,L_cm,Lp,Sp,Tp)=bra;
  //         std::tie(eta,Lr,eta_cm,L_cm,L,S,T)=ket;
  //         std::cout<<fmt::format("({} {} {} {}; {} {} {}|| {} {} {} || {} {} {} {}; {} {} {})  {}",
  //                 etap,Lrp,eta_cm,L_cm,Lp,Sp,Tp,L0,S0,T00, eta, Lr,eta_cm,L_cm,L,S,T,it->second)
  //         <<std::endl;
  //       }
  //   }

  // HalfInt Jp,J;
  // int etar,etarp;
  // RelativeCMLSJTLabels braj,ketj;
  // std::map<RelativeCMLSJTBraket,double> relative_cm_lsjt_map;
  // CMBranchLSJT(Nmax,J0,relative_cm_lst_map,relative_cm_lsjt_map);
  // for(auto it=relative_cm_lsjt_map.begin(); it!=relative_cm_lsjt_map.end(); ++it)
  //   {
  //     if(fabs(it->second)>10e-10)
  //     {
  //       std::tie(J0,T00,braj,ketj)=it->first;
  //       std::tie(etarp,Lrp,eta_cm,L_cm,Lp,Sp,Jp,Tp)=braj;
  //       std::tie(etar,Lr,eta_cm,L_cm,L,S,J,T)=ketj;
  //       int g=(etar+eta_cm)%2;
  //       int gp=(etarp+eta_cm)%2;

  //       std::cout<<fmt::format("{} {} {} {} {} {} {} {} {}  {} {} {} {} {} {} {} {} {}   {}",
  //         etarp,Lrp,eta_cm,L_cm,Lp,Sp,Jp,Tp,gp,
  //         etar,Lr,eta_cm,L_cm,L,S,J,T,g,it->second)
  //       <<std::endl;
  //     }
  //   }

  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  basis::OperatorLabelsJT operator_labels;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  
  basis::ReadRelativeOperatorLSJT(
    interaction_file,relative_lsjt_space,operator_labels,
    relative_component_sectors,relative_component_matrices, true
    );

  const basis::MatrixVector& sector_vector=relative_component_matrices[0];
  const basis::RelativeSectorsLSJT& relative_lsjt_sectors=relative_component_sectors[0];
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;

  u3shell::UpcouplingNLST(relative_lsjt_space,relative_lsjt_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);
  std::cout<<"UpcouplingNLST complete"<<std::endl;

  // Convert to LST CM space
  std::map<RelativeCMLSTBraket,double> rel_cm_lst_map;
  RelativeToCMLST(Nmax,T0, rme_nlst_map, rel_cm_lst_map);

  std::map<RelativeCMLSJTBraket,double> rel_cm_lsjt_map;
  CMBranchLSJT(Nmax, J0, rel_cm_lst_map,rel_cm_lsjt_map);
  std::cout<<"finished branching"<<std::endl;

  // Compare.
  RelativeCMLSTLabels bra_cm_lst, ket_cm_lst;
  int Lcm, Ncm, L, Lp, Nr,Nrp,Lr,Lrp,L0;
  HalfInt T00,S,T,Sp,Tp,S0,J,Jp; 
  RelativeCMLSJTLabels bra_j, ket_j;
  for(auto it=rel_cm_lsjt_map.begin(); it!=rel_cm_lsjt_map.end(); ++it)
    {
      std::tie(J0,T00,bra_j,ket_j)=it->first;
      std::tie(Nr,Lr,Ncm,Lcm,L,S,J,T)=ket_j;
      std::tie(Nrp,Lrp,Ncm,Lcm,Lp,Sp,Jp,Tp)=bra_j;
      double rme_j=it->second;
      if(fabs(rme_j)>10e-10)  
        std::cout<<fmt::format("{} {} {} {} {} {} {} {} {}  {} {} {} {} {} {} {} {} {}   {}",
        Nrp,Lrp,Ncm,Lcm,Lp,Sp,Jp,Tp, (Nrp+Ncm)%2,Nr,Lr,Ncm,Lcm,L,S,J,T,(Nr+Ncm)%2,rme_j)
        <<std::endl;
    }
  std::map<RelativeCMLSTBraket,double> relative_cm_lst_map_2;

  u3shell::RelativeRMEsU3ST rme_map;
  u3shell::UpcouplingU3ST(rme_nlst_map,J0,g0, T0, Nmax, rme_map);

  std::map<RelativeCMU3STBraket,double> relative_cm_u3st_map;
  RelativeToCMU3ST(Nmax,rme_map,relative_cm_u3st_map);
  RelativeCMU3STLabels bra_u3,ket_u3;
  u3::SU3 x0,x,xp;
  int rho0, kappa0;
  for(auto it=relative_cm_u3st_map.begin(); it!=relative_cm_u3st_map.end(); ++it)
    {
      std::tie(x0,S0,T00,kappa0, L0,bra_u3,ket_u3,rho0)=it->first;
      double rme_u3=it->second;
      std::tie(Nr,Ncm,x,S,T)=ket_u3;
      std::tie(Nrp,Ncm,xp,Sp,Tp)=bra_u3;
      std::cout<<fmt::format("{} {} {} {} {} {}  {} {} {} {} {}  {} {} {} {} {}   {}",
        x0.Str(),S0, T0, kappa0, L0, rho0, Nrp,Ncm,xp.Str(),Sp,Tp,Nr,Ncm,x.Str(),S,T,rme_u3)
      <<std::endl;
    }

  CMBranchLST(Nmax, relative_cm_u3st_map,relative_cm_lst_map_2);
  int eta,etap,eta_cm,L_cm;
  // HalfInt Sp,Tp,S0,T00,S,T;
  RelativeCMLSTLabels bra,ket;
  std::cout<<"branching"<<std::endl;
  for(auto it=relative_cm_lst_map_2.begin(); it!=relative_cm_lst_map_2.end(); ++it)
    {
      if(fabs(it->second)>10e-10)
        {
          std::tie(L0,S0,T00,bra,ket)=it->first;
          std::tie(etap,Lrp,eta_cm,L_cm,Lp,Sp,Tp)=bra;
          std::tie(eta,Lr,eta_cm,L_cm,L,S,T)=ket;
          std::cout<<fmt::format("({} {} {} {}; {} {} {}|| {} {} {} || {} {} {} {}; {} {} {})  {}",
                  etap,Lrp,eta_cm,L_cm,Lp,Sp,Tp,L0,S0,T00, eta, Lr,eta_cm,L_cm,L,S,T,it->second)
          <<std::endl;
        }
    }


  std::map<RelativeCMLSJTBraket,double> rel_cm_lsjt_map_2;
  CMBranchLSJT(Nmax, J0, relative_cm_lst_map_2,rel_cm_lsjt_map_2);
  std::cout<<"finished branching"<<std::endl;

  // Compare.
  // RelativeCMLSTLabels bra_cm_lst, ket_cm_lst;
  // int Lcm, Ncm, L, Lp, Nr,Nrp,Lr,Lrp,L0;
  // HalfInt T00,S,T,Sp,Tp,S0,J,Jp; 
  // RelativeCMLSJTLabels bra_j, ket_j;
  for(auto it=rel_cm_lsjt_map_2.begin(); it!=rel_cm_lsjt_map_2.end(); ++it)
    {
      std::tie(J0,T00,bra_j,ket_j)=it->first;
      std::tie(Nr,Lr,Ncm,Lcm,L,S,J,T)=ket_j;
      std::tie(Nrp,Lrp,Ncm,Lcm,Lp,Sp,Jp,Tp)=bra_j;
      double rme_j=it->second;
      if(fabs(rme_j)>10e-10)  
        std::cout<<fmt::format("{} {} {} {} {} {} {} {} {}  {} {} {} {} {} {} {} {} {}   {}",
        Nrp,Lrp,Ncm,Lcm,Lp,Sp,Jp,Tp, (Nrp+Ncm)%2,Nr,Lr,Ncm,Lcm,L,S,J,T,(Nr+Ncm)%2,rme_j)
        <<std::endl;
    }


  std::cout<<u3::W(u3::SU3(3,0),1,1,u3::SU3(0,3),1,1,u3::SU3(0,0),1,0,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(3,0),1,1,u3::SU3(0,3),1,1,u3::SU3(2,2),1,0,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(2,0),1,0,u3::SU3(1,0),1,1,u3::SU3(3,0),1,1,1)<<std::endl;

  std::cout<<u3::W(u3::SU3(3,0),1,1,u3::SU3(1,1),1,1,u3::SU3(2,2),1,0,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(1,1),1,1,u3::SU3(0,3),1,1,u3::SU3(2,2),1,0,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(1,1),1,1,u3::SU3(1,1),1,1,u3::SU3(2,2),1,0,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(2,0),1,0,u3::SU3(1,0),1,1,u3::SU3(1,1),1,1,1)<<std::endl;
  std::cout<<std::endl;
  std::cout<<u3::W(u3::SU3(3,0),1,1,u3::SU3(2,2),1,0,u3::SU3(3,0),1,1,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(3,0),1,1,u3::SU3(0,0),1,0,u3::SU3(3,0),1,1,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(1,1),1,1,u3::SU3(0,0),1,0,u3::SU3(1,1),1,1,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(1,1),1,1,u3::SU3(2,2),1,0,u3::SU3(1,1),1,1,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(1,1),1,1,u3::SU3(2,2),1,0,u3::SU3(3,0),1,1,1)<<std::endl;
  std::cout<<u3::W(u3::SU3(3,0),1,1,u3::SU3(2,2),1,0,u3::SU3(1,1),1,1,1)<<std::endl;



}
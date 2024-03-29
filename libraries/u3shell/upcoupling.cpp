/****************************************************************
  upcoupling.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3shell/upcoupling.h"

#include <fstream>
#include <unordered_set>
#include "am/am.h"
#include "am/wigner_gsl.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "u3shell/tensor_labels.h"

extern double zero_threshold;

namespace u3shell
{
  // Upcoupling to SO(3) x SU(2) x SU(2)
  void UpcouplingNLST(
      const basis::RelativeSpaceLSJT& space,
      const basis::RelativeSectorsLSJT& sectors,
      const std::vector<Eigen::MatrixXd>& sector_vector, 
      int J0, int g0, int T0, int Nmax,
      std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map
    )
  {
  	std::cout<<zero_threshold<<std::endl;
    for (int index=0; index<sectors.size(); ++index)
      {
        if(sector_vector[index].size()==0)
          continue; 
        Eigen::MatrixXd sector=sector_vector[index];
        auto sector_labels(sectors.GetSector(index));
        auto bra_lsjt(sector_labels.bra_subspace().labels());
        auto ket_lsjt(sector_labels.ket_subspace().labels());

        bool diagonal_sector=(bra_lsjt==ket_lsjt);
        // If the sector is diagonal, then fill in lower triangle of the sector
        if(diagonal_sector)
          {
            int nmax=sector.cols()-1;
            for(int n=0; n<=nmax; ++n)
              for(int np=0; np<n; ++np)
                {
                  double rme=sector(np,n);
                  sector(n,np)=rme;
                }
          }
  
        int Lp,L,Sp,S,Tp,T,Jp,J,g,gp;
        std::tie(L,S,J,T,g)=ket_lsjt;
        std::tie(Lp,Sp,Jp,Tp,gp)=bra_lsjt;
        u3shell::RelativeSubspaceLabelsNLST bra_nlst(Lp,Sp,Tp);              
        u3shell::RelativeSubspaceLabelsNLST ket_nlst(L,S,T);
  
        // Iterate over possible S0 and L0 values and accumulate tensor components 
        for (int S0=0; S0<=2; ++S0)            
          for (int L0=abs(S0-J0); L0<=(S0+J0); ++L0)
            {
              RelativeSectorNLST key(L0,S0,T0,bra_nlst,ket_nlst);
              // upcoupling factor
              double so3_coef=am::Unitary9J(L,S,J, L0,S0,J0, Lp,Sp,Jp)
                              *am::dim(Jp)*am::dim(S0)*am::dim(L0)/am::dim(J0)/am::dim(Sp)/am::dim(Lp);

              double so3_coef_conj=0;
              RelativeSectorNLST key_conj;
              // See notes on conjugation in basis:lsjt_operator.h 
              if(not diagonal_sector)
                {
                  so3_coef_conj=am::Unitary9J(Lp,Sp,Jp,L0,S0,J0,L,S,J)
                    *am::dim(J)*am::dim(S0)*am::dim(L0)/am::dim(J0)/am::dim(S)/am::dim(L)
                    *ParitySign(Tp-T)*Hat(Tp)/Hat(T)
                    *ParitySign(Jp-J)*Hat(Jp)/Hat(J);
                  key_conj=RelativeSectorNLST(L0,S0,T0,ket_nlst,bra_nlst);
                }
              
              // if((L==20)&&(S==1)&& T==0 && L==Lp && S==Sp && T==Tp && L0==2 && S0==1 && T0==0)
              //   {
              //     std::cout<<Jp<<"  "<<J<<"  "<<so3_coef<<"  "<<so3_coef_conj<<std::endl;
              //     std::cout<<so3_coef*sector<<std::endl<<so3_coef_conj*sector.transpose()<<std::endl<<std::endl;
              //   }
              // If lst sector is not in map, initialize to value
              if(not rme_nlst_map.count(key))
                rme_nlst_map[key]=so3_coef*sector;
              else
                rme_nlst_map[key]+=so3_coef*sector;

              if(not diagonal_sector)
                {
                  if(not rme_nlst_map.count(key_conj))
                    rme_nlst_map[key_conj]=so3_coef_conj*sector.transpose();
                  else
                    rme_nlst_map[key_conj]+=so3_coef_conj*sector.transpose();
                }
              // if((L==20)&&(S==1)&& T==0 && L==Lp && S==Sp && T==Tp && L0==2 && S0==1 && T0==0)
              //   std::cout<<rme_nlst_map[key_conj]<<"  "<<rme_nlst_map[key]<<std::endl;
            }
      }

      // Enforcing Hermiticity on RMEs
      u3shell::RelativeSubspaceLabelsNLST bra_nlst,ket_nlst;
      int L,S,T,Lp,Sp,Tp,L0,S0;
      T0 = 0;
      for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
      {
        std::tie(L0,S0,T0,bra_nlst,ket_nlst)=it->first;
        Eigen::MatrixXd& sector(it->second);
        int nmax=sector.cols()-1;
        int npmax=sector.rows()-1;

        for (int np=0; np<=npmax; ++np)
          {
            int Np=2*np+Lp;
            if(Np>Nmax)
              continue;
            u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
            for (int n=0; n<=nmax; ++n)
              {
                int N=2*n+L;
                if(N>Nmax)
                  continue;
                u3shell::RelativeStateLabelsU3ST ket(N,S,T);
                //Extract rme
                double rme_nlst=sector(np,n);
                if (fabs(rme_nlst)<=zero_threshold) {
                  RelativeSectorNLST key(L0,S0,T0,bra_nlst,ket_nlst);
                  RelativeSectorNLST key_conj=RelativeSectorNLST(L0,S0,T0,ket_nlst,bra_nlst);
                  // Set both RME and conjugate to 0 is one is 0
                  rme_nlst_map[key](np,n) = 0;
                  rme_nlst_map[key_conj](n,np) = 0;
                  continue; 
                }
              }
          }
      }
  }

  // Upcoupling from SO(3) x SU(2) x SU(2) to SU(3) x SU(2) x SU(2)
  void UpcouplingU3ST(
      std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map,
      int Nmax,
      u3::WCoefCache& w_cache,
      RelativeRMEsU3ST& rme_map
    )
  {
    u3shell::RelativeSubspaceLabelsNLST bra_nlst,ket_nlst;
    int L,S,T,Lp,Sp,Tp,L0,S0,T0;
    for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
      {
        std::tie(L0,S0,T0,bra_nlst,ket_nlst)=it->first;
        std::tie(L,S,T)=ket_nlst;
        std::tie(Lp,Sp,Tp)=bra_nlst;
        const Eigen::MatrixXd& sector(it->second);
        int nmax=sector.cols()-1;
        int npmax=sector.rows()-1;

        for (int np=0; np<=npmax; ++np)
          {
            int Np=2*np+Lp;
            if(Np>Nmax)
              continue;
            u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
            for (int n=0; n<=nmax; ++n)
              {
                int N=2*n+L;
                if(N>Nmax)
                  continue;
                u3shell::RelativeStateLabelsU3ST ket(N,S,T);
                //Extract rme
                double rme_nlst=sector(np,n);
                //if (fabs(rme_nlst)<=zero_threshold)//REMOVE 
                  //continue;
                // generate list of allowed x0's from coupling bra and ket
                MultiplicityTagged<u3::SU3>::vector x0_set=u3::KroneckerProduct(bra.x(),u3::Conjugate(ket.x()));
                for(int i=0; i<x0_set.size(); i++)
                  {
                    u3::SU3 x0(x0_set[i].irrep);
                    if ((x0.lambda()+x0.mu())>=82)
                      {
                        std::cout<<fmt::format("x0={} is out of bounds for W coefficeint.  RME={}",x0.Str(),rme_nlst)<<std::endl;
                        continue;
                      }
                    u3shell::RelativeUnitTensorLabelsU3ST operator_labels(x0,S0,T0,bra,ket);
                    int kappa0_max=u3::BranchingMultiplicitySO3(x0, L0);
                    for(int kappa0=1; kappa0<=kappa0_max; ++kappa0)
                      {
                        double u3_coef=u3::WCached(w_cache,ket.x(), 1,L, x0, kappa0, L0, bra.x(),1,Lp,1)
                          *u3::dim(x0)*am::dim(Lp)/(u3::dim(bra.x())*am::dim(L0));
                        std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(operator_labels,kappa0,L0);
                        rme_map[key]+=u3_coef*rme_nlst*ParitySign(n+np);
                        /*if (fabs(rme_map[key]) <= zero_threshold) {
                          rme_map[key] = 0;
                        }*/
                      }
                  }
              }
          }
      }
    // Zero out rmes 
    u3shell::RelativeUnitTensorLabelsU3ST operator_labels;
    int kappa0,Np,N;
    u3::SU3 x0;
    double rme;
    for(auto it=rme_map.begin(); it!=rme_map.end(); it++)
      {
        std::tie(operator_labels,kappa0,L0)=it->first;
        rme=it->second;
        if (fabs(rme)<=zero_threshold)
          {
            rme_map[it->first]=0.0;

            // enforcing hermiticity
            u3shell::RelativeStateLabelsU3ST bra(N,S,T), ket(Np,Sp,Tp);
            RelativeUnitTensorLabelsU3ST unit_tensor_conj(u3::Conjugate(x0),S0,T0,bra,ket);
            std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key_conj(unit_tensor_conj,kappa0,L0);
            rme_map[key_conj]=0.0;
          }
      }
  }

  // If no cache is passed as an arguemnt, creates cache for use in calculations
  void UpcouplingU3ST(
      std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map,
      int Nmax,
      RelativeRMEsU3ST& rme_map
    )
  {
    u3::WCoefCache w_cache;
    UpcouplingU3ST(rme_nlst_map,Nmax,w_cache,rme_map);
  }

  void Upcoupling(    
      const basis::RelativeSpaceLSJT& space,
      const std::array<basis::RelativeSectorsLSJT,3>& T0_sector_labels,
      const std::array<basis::MatrixVector,3>& T0_sectors,
      u3::WCoefCache& w_cache, 
      int J0, int g0, int T0,int Nmax,
      RelativeRMEsU3ST& rme_map
    )
  {

    // if T0=-1, then upcouple all possible values of T0, i.e., T0=0,1,2;
    // otherwise upcouple just the desired T0 sector. 
    int T0_min=(T0==-1)?0:T0;
    int T0_max=(T0==-1)?2:T0;

    std::map<RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
    for(int T00=T0_min; T00<=T0_max; ++T00)
      {
        const basis::RelativeSectorsLSJT& sector_labels=T0_sector_labels[T00];
        const basis::OperatorBlocks<double>& sectors=T0_sectors[T00];
        u3shell::UpcouplingNLST(space,sector_labels,sectors,J0,g0,T00,Nmax,rme_nlst_map);
      }
    
    u3shell::UpcouplingU3ST(rme_nlst_map,Nmax,w_cache,rme_map);
  }

  void Upcoupling(    
      const basis::RelativeSpaceLSJT& space,
      const std::array<basis::RelativeSectorsLSJT,3>& T0_sector_labels,
      const std::array<basis::OperatorBlocks<double>,3>& T0_sectors,
      int J0, int g0, int T0,int Nmax,
      RelativeRMEsU3ST& rme_map
    )
  {
    u3::WCoefCache w_cache;
    Upcoupling(space, T0_sector_labels, T0_sectors, w_cache,J0, g0, T0, Nmax,rme_map);
  }

 void UpcoupleCMU3ST(
      std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map,
      u3::WCoefCache& w_cache,
      RelativeCMUnitTensorCache& rel_cm_u3st_map
    )
  {
    int Nrp,Lrp,Ncm,Lcm,Nr,Lr,Lp,L,L0;
    HalfInt Sp,Tp,S,T,S0,T0;
    RelativeCMStateLabelsNLST bra_cm,ket_cm;
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

                          for(int kappa0=1; kappa0<=kappa0_max; ++kappa0)
                            for(int rho0=1; rho0<=rho0_max; ++rho0)
                              {
                                u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st(x0,S0,T0,rho0,bra,ket); 
                                rel_cm_u3st_map[braket_u3st]
                                  +=am::dim(Lp)/u3::dim(xp)
                                  //u3::dim(x0)*am::dim(Lp)/u3::dim(xp)/am::dim(L0)
                                  *u3::WCached(w_cache,u3::SU3(Nrp,0),1,Lrp,u3::SU3(Ncm,0),1,Lcm,xp,kappap,Lp,1)
                                  *u3::WCached(w_cache,u3::SU3(Nr,0),1,Lr,u3::SU3(Ncm,0),1,Lcm,x,kappa,L,1)
                                  *u3::WCached(w_cache,x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0)
                                  *rme;
                              }
                        }
                    }
              }
          }
      }
  }

  void UpcoupleCMU3ST(
      std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map,
      RelativeCMUnitTensorCache& rel_cm_u3st_map
    )
  {
    u3::WCoefCache w_cache;
    UpcoupleCMU3ST(rel_cm_lst_map, w_cache,rel_cm_u3st_map);
  }

  void WriteRelativeOperatorU3ST(const std::string& filename, RelativeRMEsU3ST& relative_rmes, bool hermitian)
  {
    u3shell::RelativeUnitTensorLabelsU3ST labels;
    int kappa0,L0,etap,eta;
    HalfInt Sp,Tp,S,T,S0,T0;
    u3::SU3 x0;
    double rme;

    // Zero out rmes 
    for(auto it=relative_rmes.begin(); it!=relative_rmes.end(); it++)
      {
        std::tie(labels,kappa0,L0)=it->first;
        std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=labels.FlatKey();
        rme=it->second;
        if (fabs(rme)<zero_threshold)
          {
            relative_rmes[it->first]=0.0;

            // enforcing hermiticity
            if(hermitian)
              {
                u3shell::RelativeStateLabelsU3ST bra(eta,S,T), ket(etap,Sp,Tp);
                RelativeUnitTensorLabelsU3ST unit_tensor_conj(u3::Conjugate(x0),S0,T0,bra,ket);
                std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key_conj(unit_tensor_conj,kappa0,L0);
                 relative_rmes[key_conj]=0.0;
              }
          }
      }

    std::ofstream os(filename);
    // Label Header
    os 
      << "# RELATIVE U3ST" << std::endl
      << "#   N' S' T'   N S T   lamda mu S0 T0 kappa0 L0   RME" << std::endl;
    for(auto it=relative_rmes.begin(); it!=relative_rmes.end(); it++)
      {
        std::tie(labels,kappa0,L0)=it->first;
        std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=labels.FlatKey();

        rme=it->second;
        // if (fabs(rme)<zero_threshold)
        //   rme=0.0;
        const int width=3;
        const int precision=16;

        os << std::setprecision(precision);
        os
          << " " << std::setw(width) << etap
          << " " << std::setw(width) << Sp
          << " " << std::setw(width) << Tp
          << " " << "  "
          << " " << std::setw(width) << eta
          << " " << std::setw(width) << S
          << " " << std::setw(width) << T
          << " " << "  "
          << " " << std::setw(width) << x0.lambda()
          << " " << std::setw(width) << x0.mu()
          << " " << std::setw(width) << S0 
          << " " << std::setw(width) << T0
          << " " << std::setw(width) << kappa0
          << " " << std::setw(width) << L0
          << " " << "  "
          << " " << std::showpoint << std::scientific << rme
          << std::endl;
      }
  }

  void ReadRelativeOperatorU3ST(const std::string& filename, RelativeRMEsU3ST& relative_rmes)
  {
    //OBSOLETE, only used in test programs
    // open file
    std::ifstream in_stream(filename);
    mcutils::StreamCheck(bool(in_stream),filename,"Failure opening relative rme file");
    std::cout << fmt::format("Reading relative rmes from {}...",filename) << std::endl;

    // process stream
    //
    int etap,eta,lambda0,mu0,kappa0, L0;
    int Sp,Tp,S,T,S0,T0;
    double rme;
    std::string line;

    int line_count = 0;
    while(std::getline(in_stream,line))
      {
        ++line_count;
        // if(line_count<3)
        // {
        //   std::cout<<line<<std::endl;
        //   continue;
        // }
        std::istringstream line_stream(line);
        line_stream>>etap>>Sp>>Tp>>eta>>S>>T>>lambda0>>mu0>>S0>>T0>>kappa0>>L0>>rme;
        mcutils::ParsingCheck(line_stream,line_count,line);
        if (fabs(rme)>zero_threshold)
          {
            u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp), ket(eta,S,T);
            RelativeUnitTensorLabelsU3ST unit_tensor(u3::SU3(lambda0,mu0),S0,T0,bra,ket);
            std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(unit_tensor,kappa0,L0);
            relative_rmes[key]=rme;
          }
      }

    // close file
    in_stream.close();
  }


  void ReadRelativeOperatorU3ST(
    int Nmax, int N1v,
    const std::string& filename, 
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    RelativeRMEsU3SSubspaces& relative_rmes
    )
  {
    int etap,eta,lambda0,mu0,kappa0, L0;
    int Sp,Tp,S,T,S0,T0;
    double rme;
    std::string line;
    int eta_max=Nmax+2*N1v;

    std::ifstream in_stream(filename);
    mcutils::StreamCheck(bool(in_stream),filename,"Failure opening relative rme file");
    std::cout << fmt::format("Reading relative rmes from {}...",filename) << std::endl;

    int line_count = 0;
    while(std::getline(in_stream,line))
      {
        ++line_count;

        if(line_count<3)
        {
          // std::cout<<line<<std::endl;
          continue;
        }

        // std::cout<<line<<std::endl;
        std::istringstream line_stream(line);
        line_stream>>etap>>Sp>>Tp>>eta>>S>>T>>lambda0>>mu0>>S0>>T0>>kappa0>>L0>>rme;
        mcutils::ParsingCheck(line_stream,line_count,line);
        if((etap>eta_max)||(eta>eta_max)||(abs(etap-eta)>Nmax))
          continue;
        if (fabs(rme)>zero_threshold)
          {
            // Look up unit tensor subspace and state index for relative unit tensor
            std::tuple<u3::SU3, HalfInt, int, int> subspace_labels(u3::SU3(lambda0,mu0),S0,etap,eta);
            int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(subspace_labels);
            const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
            u3shell::RelativeUnitTensorStateLabelsU3S state_labels(T0,Sp,Tp,S,T);
            int unit_tensor_state_index=unit_tensor_subspace.LookUpStateIndex(state_labels);

            // Storing RME
            std::tuple<int,int,int> key(unit_tensor_subspace_index,kappa0,L0);
            
            // Check if subspace+op labels already in map, if not, add and resize corresponding vector rmes
            if(not relative_rmes.count(key))
              relative_rmes[key].resize(unit_tensor_subspace.size());

            // Add rme
            relative_rmes[key][unit_tensor_state_index]=rme;
          }
      }
  }


  typedef std::unordered_map<RelativeUnitTensorLabelsU3ST,std::tuple<int,int,double>,boost::hash<RelativeUnitTensorLabelsU3ST>>
                    InteractionCache;


  double SU3Casimir(const u3::SU3& x0)
    {
      int lm=x0.lambda();
      int mu=x0.mu();
      return 2./3*(lm*lm+mu*mu+lm*mu+3*lm+3*mu);
    }

  void U3DecomposeRelativeOperatorU3ST(
    int Nmax,
    double hbar_omega,
    const std::string& filename
    )
  {
    int etap,eta,lambda0,mu0,kappa0, L0;
    int Sp,Tp,S,T,S0,T0;
    double rme;
    std::string line;
    std::map<std::pair<int,u3::SU3>,double> u3_decomposition;
    std::ifstream in_stream(filename);
    mcutils::StreamCheck(bool(in_stream),filename,"Failure opening relative rme file");
    std::cout << fmt::format("Reading relative rmes from {}...",filename) << std::endl;

    int line_count = 0;
    double total; 
    while(std::getline(in_stream,line))
      {
        std::istringstream(line)>>etap>>Sp>>Tp>>eta>>S>>T>>lambda0>>mu0>>S0>>T0>>kappa0>>L0>>rme;
        if (fabs(rme)>zero_threshold)
          {
            int N0=etap-eta;
            u3::SU3 x0(lambda0,mu0);
            std::pair<int,u3::SU3> key(N0,x0);
            u3_decomposition[key]+=rme*rme;
            total+=rme*rme;
          }
      }
    std::ofstream os(fmt::format("u3_decomposition_hw{:.1f}_Nmax{:02d}_u3st.dat",hbar_omega,Nmax));

    for(auto it=u3_decomposition.begin(); it!=u3_decomposition.end(); it++)
      {
        u3::SU3 x0; 
        int N0;
        std::tie(N0,x0)=it->first;
        double casimir=SU3Casimir(x0);
        double amplitude=it->second/total;
        const int width=3;
        const int precision=16;
        os << std::setprecision(precision);
        os
          << " " << std::setw(width) << N0
          << " " << std::setw(width) << x0.lambda()
          << " " << std::setw(width) << x0.mu()
          << " " << "  "
          << " " << std::showpoint << casimir
          << " " << std::showpoint << std::scientific << amplitude
          << std::endl;
      }


  }




}//end namespace

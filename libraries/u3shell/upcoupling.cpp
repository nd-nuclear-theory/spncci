/****************************************************************
  upcoupling.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/14/16 (aem,mac): Created to import relative jisp16 files.
****************************************************************/
#include "u3shell/upcoupling.h"

#include "cppformat/format.h"

#include "am/wigner_gsl.h"
// #include "sp3rlib/u3.h"
#include "u3shell/tensor_labels.h"

extern double zero_threshold;

namespace u3shell
{

  void UpcouplingNLST(
    const basis::RelativeSpaceLSJT& space,
    const basis::RelativeSectorsLSJT& sectors,
    const std::vector<Eigen::MatrixXd>& sector_vector, 
    int J0, int g0, int T0, int Nmax,
    std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map
    )
  {
    for (int index=0; index<sectors.size(); ++index)
      {
        if(sector_vector[index].size()==0)
          continue;
        Eigen::MatrixXd sector=sector_vector[index];
        auto sector_labels(sectors.GetSector(index));
        auto bra_lsjt(sector_labels.bra_subspace().GetSubspaceLabels());
        auto ket_lsjt(sector_labels.ket_subspace().GetSubspaceLabels());
        int Lp,L,Sp,S,Tp,T,Jp,J,g,gp;
        std::tie(L,S,J,T,g)=ket_lsjt;
        std::tie(Lp,Sp,Jp,Tp,gp)=bra_lsjt;
        u3shell::RelativeSubspaceLabelsNLST bra_nlst(Lp,Sp,Tp);              
        u3shell::RelativeSubspaceLabelsNLST ket_nlst(L,S,T);
        // For now we assume T==Tp for all operators.  If this is not the case, then need to change
        // symmetry relations 
        assert(T==Tp);
        // if((T0<abs(T-Tp))||(T0>(T+Tp)))
        //   continue;
        bool diagonal_sector=(bra_lsjt==ket_lsjt);
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
        for (int S0=0; S0<=2; ++S0)  
          for (int L0=abs(S0-J0); L0<=(S0+J0); ++L0)
            {
              RelativeSectorNLST key=make_tuple(L0,S0,bra_nlst,ket_nlst);
              double so3_coef=am::Unitary9J(L,S,J, L0,S0,J0, Lp,Sp,Jp)
                                *(2.*Jp+1)*(2.*S0+1)*(2.*L0+1)/(2.*J0+1)/(2.*Sp+1)/(2.*Lp+1);
              
              if (fabs(so3_coef)<zero_threshold)
                continue;

              double so3_coef_conj=0;
              RelativeSectorNLST key_conj;
              if(not diagonal_sector)
              {
                so3_coef_conj=am::Unitary9J(Lp,Sp,Jp,L0,S0,J0,L,S,J)
                                *(2.*J+1)*(2.*S0+1)*(2.*L0+1)/(2.*J0+1)/(2.*S+1)/(2.*L+1);
                key_conj=make_tuple(L0,S0,ket_nlst,bra_nlst);
              }
              
              if(not rme_nlst_map.count(key))
                {
                  rme_nlst_map[key]=so3_coef*sector;
                  if(not diagonal_sector)
                    rme_nlst_map[key_conj]=so3_coef_conj*sector.transpose();
                }
              else
                {
                  rme_nlst_map[key]+=so3_coef*sector;
                  if(not diagonal_sector)
                    rme_nlst_map[key_conj]+=so3_coef_conj*sector.transpose();
                }
            }
      }
    // std::cout<<"zero out"<<std::endl;
    // std::vector<RelativeSectorNLST> remove_list;
    // for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end();++it)
    //   {
    //     int del=1;
    //     Eigen::MatrixXd& sector_out(it->second);
    //     int num_cols=sector_out.cols();
    //     int num_rows=sector_out.rows();
    //     for (int i=0; i<num_rows; ++i)
    //       for (int j=0; j<num_cols; ++j)
    //         if(fabs(sector_out(i,j))<zero_threshold)
    //           rme_nlst_map[it->first](i,j)=0.0;
    //         else
    //           del=0;
        
    //     if(del==1)
    //       remove_list.push_back(it->first);
    //   }
    // for(int r=0; r<remove_list.size(); ++r)
    //   rme_nlst_map.erase(remove_list[r]);
  }

  void UpcouplingU3ST(
    std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map,
    int T0, int Nmax,
    u3::WCoefCache& w_cache,
    RelativeRMEsU3ST& rme_map
    )
  {
    u3shell::RelativeSubspaceLabelsNLST bra_nlst,ket_nlst;
    int L,S,T,Lp,Sp,Tp,L0,S0;
    for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
      {
        std::tie(L0,S0,bra_nlst,ket_nlst)=it->first;
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
                if (fabs(rme_nlst)<=zero_threshold)//REMOVE
                  continue;
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
                                      *u3::dim(x0)*(2.*Lp+1)/(u3::dim(bra.x())*(2.*L0+1));
                        std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(operator_labels,kappa0,L0);
                        rme_map[key]+=u3_coef*rme_nlst*parity(n+np);
                      }
                  }
              }
          }
      }
    // // Removing zeros 
    // std::vector<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>> remove_list;
    // for(auto it=rme_map.begin(); it!=rme_map.end(); it++)
    //   {
    //     if(fabs(it->second)<zero_threshold)
    //       remove_list.push_back(it->first);
    //   }
    // for(auto tensor : remove_list)
    //   rme_map.erase(tensor);
  }

  // If no cache is passed as an arguemnt, creates cache for use in calculations
  void UpcouplingU3ST(
    std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map,
    int T0, int Nmax,
    RelativeRMEsU3ST& rme_map
    )
    {
      u3::WCoefCache w_cache;
      UpcouplingU3ST(rme_nlst_map,T0, Nmax,w_cache,rme_map);
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

  void WriteRelativeOperatorU3ST(std::ostream& os, const RelativeRMEsU3ST& relative_rmes)
    {
      u3shell::RelativeUnitTensorLabelsU3ST labels;
      int kappa0,L0,etap,eta;
      HalfInt Sp,Tp,S,T,S0,T0;
      u3::SU3 x0;
      double rme;
      for(auto it=relative_rmes.begin(); it!=relative_rmes.end(); it++)
        {
          std::tie(labels,kappa0,L0)=it->first;
          std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=labels.FlatKey();

          rme=it->second;
          if (fabs(rme)<zero_threshold)
            rme=0.0;
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
  void ReadRelativeOperatorU3ST(std::istream& is, RelativeRMEsU3ST& relative_rmes)
  {
    int etap,eta,lambda0,mu0,kappa0, L0;
    int Sp,Tp,S,T,S0,T0;
    double rme;
    std::string line;

    while(std::getline(is,line))
    {
      std::istringstream(line)>>etap>>Sp>>Tp>>eta>>S>>T>>lambda0>>mu0>>S0>>T0>>kappa0>>L0>>rme;
      if (fabs(rme)>zero_threshold)
        {
          u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp), ket(eta,S,T);
          RelativeUnitTensorLabelsU3ST unit_tensor(u3::SU3(lambda0,mu0),S0,T0,bra,ket);
          std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(unit_tensor,kappa0,L0);
          relative_rmes[key]=rme;
        }
    }
  }

  typedef std::unordered_map<RelativeUnitTensorLabelsU3ST,std::tuple<int,int,double>,boost::hash<RelativeUnitTensorLabelsU3ST>>
    InteractionCache;

 void ReadRelativeOperatorU3ST(std::istream& is, InteractionCache& relative_rmes)
  {
    int etap,eta,lambda0,mu0,kappa0, L0;
    int Sp,Tp,S,T,S0,T0;
    double rme;
    std::string line;

    while(std::getline(is,line))
    {
      std::istringstream(line)>>etap>>Sp>>Tp>>eta>>S>>T>>lambda0>>mu0>>S0>>T0>>kappa0>>L0>>rme;
      if (fabs(rme)>zero_threshold)
        {
          u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp), ket(eta,S,T);
          RelativeUnitTensorLabelsU3ST key(u3::SU3(lambda0,mu0),S0,T0,bra,ket);
          std::tuple<int,int,double> value(kappa0,L0,rme);
          relative_rmes[key]=value;
        }
    }
  }


  void Upcoupling(    
    const basis::RelativeSpaceLSJT& space,
    const basis::RelativeSectorsLSJT& sectors,
    const std::vector<Eigen::MatrixXd>& sector_vector,
    u3::WCoefCache& w_cache, 
    int J0, int g0, int T0,int Nmax,
    RelativeRMEsU3ST& rme_map
    )
  {
    std::map<RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
    u3shell::UpcouplingNLST(space,sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);
    u3shell::UpcouplingU3ST(rme_nlst_map,T0,Nmax,w_cache,rme_map);
  }

  void Upcoupling(    
    const basis::RelativeSpaceLSJT& space,
    const basis::RelativeSectorsLSJT& sectors,
    const std::vector<Eigen::MatrixXd>& sector_vector, 
    int J0, int g0, int T0,int Nmax,
    RelativeRMEsU3ST& rme_map
    )
  {
    u3::WCoefCache w_cache;
    Upcoupling(space,sectors, sector_vector, w_cache,J0, g0, T0, Nmax,rme_map);
  }

  void GetInteractionMatrix(
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

}//end namespace
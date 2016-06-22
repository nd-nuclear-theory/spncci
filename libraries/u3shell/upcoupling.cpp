/****************************************************************
  upcoupling.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/14/16 (aem,mac): Created to import relative jisp16 files.
****************************************************************/
#include "u3shell/upcoupling.h"

#include "am/wigner_gsl.h"
#include "cppformat/format.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/import_interaction.h" 

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
        const Eigen::MatrixXd& sector=sector_vector[index];
        auto sector_labels(sectors.GetSector(index));
        auto bra_lsjt(sector_labels.bra_subspace().GetSubspaceLabels());
        auto ket_lsjt(sector_labels.ket_subspace().GetSubspaceLabels());
        int Lp,L,Sp,S,Tp,T,Jp,J,g,gp;
        std::tie(L,S,J,T,g)=ket_lsjt;
        std::tie(Lp,Sp,Jp,Tp,gp)=bra_lsjt;
        u3shell::RelativeSubspaceLabelsNLST bra_nlst=std::make_tuple(Lp,Sp,Tp);              
        u3shell::RelativeSubspaceLabelsNLST ket_nlst=std::make_tuple(L,S,T);
        if((T0<abs(T-Tp))||(T0>(T+Tp)))
          continue;
        for (int S0=0; S0<=2; ++S0)  
          for (int L0=abs(S0-J0); L0<=(S0+J0); ++L0)
            {
              double so3_coef=am::Unitary9J(L,S,J, L0,S0,J0, Lp,Sp,Jp)
                                *(2.*Jp+1)*(2.*S0+1)*(2.*L0+1)/(2.*J0+1)/(2.*Sp+1)/(2.*Lp+1);
              if (so3_coef==0)
                continue;
              RelativeSectorNLST key=make_tuple(L0,S0,bra_nlst,ket_nlst);
              if(not rme_nlst_map.count(key))
                {
                rme_nlst_map[key]=so3_coef*sector;
                }
              else
                {
                  rme_nlst_map[key]+=so3_coef*sector;
                }
            }
      }
    // std::cout<<"zero out"<<std::endl;
    
    for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end();++it)
      {
        int del=1;
        Eigen::MatrixXd& sector_out(it->second);
        int num_cols=sector_out.cols();
        int num_rows=sector_out.rows();
        for (int i=0; i<num_rows; ++i)
          for (int j=0; j<num_cols; ++j)
            if((fabs(sector_out(i,j))<10e-10)||(fabs(sector_out(i,j)>10e8)))
              rme_nlst_map[it->first](i,j)=0.0;
            else
              del=0;
        if(del==1)
          rme_nlst_map.erase(it->first);
        // else
        //   std::cout<<rme_nlst_map[it->first]<<std::endl;

     }
    // std::cout<<"Exiting"<<std::endl;
  }

  void UpcouplingU3ST(
    std::map<RelativeSectorNLST,Eigen::MatrixXd>& rme_nlst_map,
    int J0, int g0, int T0, int Nmax,
    std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int>,double>& rme_map
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

        for (int np=0; np<=nmax; ++np)
          {
            int Np=2*np+Lp;
            if(Np>Nmax)
              continue;
            u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
            for (int n=0; n<=nmax; ++n)
              {
                int N=2*n+L;
                if(Np>Nmax)
                  continue;
                u3shell::RelativeStateLabelsU3ST ket(N,S,T);
                //Extract rme
                double rme_nlst=sector(np,n);
                if (fabs(rme_nlst)<=10e-10)
                  continue;
                if (fabs(rme_nlst>=10e10))
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
                        double u3_coef=u3::W(ket.x(), 1,L, x0, kappa0, L0, bra.x(),1,Lp,1)*u3::dim(x0)*(2.*Lp+1)/u3::dim(bra.x())/(2.*L0+1);
                        std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int> key(operator_labels,kappa0);
                        rme_map[key]+=u3_coef*rme_nlst;
                      }
                  }
              }
          }
      }
  }

 
  void Upcoupling(    
    const basis::RelativeSpaceLSJT& space,
    const basis::RelativeSectorsLSJT& sectors,
    const std::vector<Eigen::MatrixXd>& sector_vector, 
    int J0, int g0, int T0,int Nmax,
    std::map<std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int>,double>& rme_map
    )
  {
    std::map<RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
    u3shell::UpcouplingNLST(space,sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);
    u3shell::UpcouplingU3ST(rme_nlst_map,J0,g0,T0,Nmax,rme_map);
  }
}//end namespace
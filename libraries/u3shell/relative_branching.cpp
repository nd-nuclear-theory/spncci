/****************************************************************
  relative_branching.cpp
****************************************************************/
#include "u3shell/relative_branching.h"

#include "am/wigner_gsl.h"
#include "cppformat/format.h"

extern double zero_threshold;

namespace u3shell
{
  void BranchRelativeNLST(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    std::unordered_map<u3shell::RelativeStateSectorNLST,double,boost::hash<u3shell::RelativeStateSectorNLST>>& rmes_nlst
    )
  {

    std::cout<<"branching "<<std::endl;
    u3::WCoefCache w_cache;
 
    for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
      {
      // extract lables
        u3shell::RelativeUnitTensorLabelsU3ST tensor;
        int kappa0, L0, N, Np;
        HalfInt S0, T0, Sp,Tp,S,T;
        u3::SU3 x0;
        std::tie(tensor,kappa0,L0)=it->first;
        std::tie(x0,S0,T0,Np,Sp,Tp,N,S,T)=tensor.FlatKey();
        double rme=it->second;

        // branch to L and sum over SU(3) labels
        MultiplicityTagged<int>::vector L0_branch=BranchingSO3(x0);
        for(int Lp=Np%2; Lp<=Np; Lp+=2)
          for(int L=N%2; L<=N; L+=2)
            {
              if(not am::AllowedTriangle(Lp,L0,L))
                continue;

              int n((N-L)/2);
              int np((Np-Lp)/2);
              u3shell::RelativeStateLabelsNLST bra_nlst(Np,Lp,Sp,Tp);              
              u3shell::RelativeStateLabelsNLST ket_nlst(N,L,S,T);
              u3shell::RelativeStateSectorNLST sector(L0,S0,T0,bra_nlst,ket_nlst);
              // double coef=u3::WCached(w_cache,u3::SU3(N,0), 1,L, x0, kappa0, L0, u3::SU3(Np,0),1,Lp,1);
              rmes_nlst[sector]+=parity(n+np)*u3::WCached(w_cache,u3::SU3(N,0), 1,L, x0, kappa0, L0, u3::SU3(Np,0),1,Lp,1)*rme;
              // std::cout<<fmt::format("{} {} {} {}  {} {} {} {}  {} {} {} {} {}  {}  {} ",Np,Lp,Sp,Tp,N,L,S,T, x0, kappa0,L0,S0,T0, )
            }
      }
  }


  void BranchRelativeNLSJT(
    const basis::OperatorLabelsJT& operator_labels, 
    int Nmax, int Jmax,
    const basis::RelativeSpaceLSJT& relative_space_lsjt,
    const std::unordered_map<u3shell::RelativeStateSectorNLST,double,boost::hash<u3shell::RelativeStateSectorNLST>>& rmes_nlst,
    std::array<basis::RelativeSectorsLSJT,3>& isospin_component_sectors_lsjt,
    std::array<basis::MatrixVector,3>& isospin_component_blocks_lsjt
    )

  {  // setting up J sectors 
    basis::ConstructZeroOperatorRelativeLSJT(
        operator_labels,relative_space_lsjt,
        isospin_component_sectors_lsjt,isospin_component_blocks_lsjt
      );


    // branch to J and populate lsjt blocks
    for(auto it=rmes_nlst.begin(); it!=rmes_nlst.end(); ++it)
      {
        // extract labels
        int L0,L,Lp,N,Np;
        HalfInt S0,T0,Sp,Tp,S,T;
        u3shell::RelativeStateLabelsNLST bra_nlst, ket_nlst;
        std::tie(L0,S0,T0,bra_nlst,ket_nlst)=it->first;
        std::tie(N,L,S,T)=ket_nlst;
        std::tie(Np,Lp,Sp,Tp)=bra_nlst;
        HalfInt J0=operator_labels.J0;

        double rme=it->second;

        // Check if N valid for given Nmax truncation
        if( N>Nmax || Np>Nmax )
          continue;

        // Check angular momentum coupling
        if(not am::AllowedTriangle(L0,S0,J0))
          continue;

        // Get index for block vector
        int T00(T0);

        // Get indices for matrix element in block
        int np((Np-Lp)/2);
        int n((N-L)/2);

        // For Jp and J, populate corresponding block
        for(HalfInt Jp=abs(Lp-Sp); Jp<=Lp+Sp; ++Jp)
          for(HalfInt J=abs(L-S); J<=L+S; ++J)
            {
              if(J>Jmax || Jp >Jmax)
                continue;

              if(not am::AllowedTriangle(J,J0,Jp))
                continue;

              int subspace_index_bra = relative_space_lsjt.LookUpSubspaceIndex(
                basis::RelativeSubspaceLSJTLabels(Lp,Sp,Jp,Tp,Lp%2)
                );

              int subspace_index_ket = relative_space_lsjt.LookUpSubspaceIndex(
                basis::RelativeSubspaceLSJTLabels(L,S,J,T,L%2)
                );

              // only store canonical ordered subspaces because others obtained by symmetry
              if(subspace_index_bra>subspace_index_ket)
                continue;

              // only store upper triangle of diagonal sectors
              if((subspace_index_bra==subspace_index_ket)&&(np>n))
                continue;

              // Look up block index for block vector
              int sector_index = isospin_component_sectors_lsjt[T00].LookUpSectorIndex(subspace_index_bra,subspace_index_ket);
              auto& block=isospin_component_blocks_lsjt[T00][sector_index];
              block(np,n)+=am::Unitary9J(L,S,J, L0,S0,J0, Lp,Sp,Jp)*rme;
            }
      }
  }

  void BranchRelativeRMEs(const basis::OperatorLabelsJT& operator_labels,int Nmax, int Jmax, 
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    const basis::RelativeSpaceLSJT& relative_space_lsjt,
    std::array<basis::RelativeSectorsLSJT,3>& isospin_component_sectors_lsjt,
    std::array<basis::MatrixVector,3>& isospin_component_blocks_lsjt
    )
  {  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // branching interaction
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout<<"branching LST"<<std::endl;
    
    std::unordered_map<u3shell::RelativeStateSectorNLST,double,boost::hash<u3shell::RelativeStateSectorNLST>> rmes_nlst;
    u3shell::BranchRelativeNLST( interaction_u3st,rmes_nlst);

    // Generatesectors for each T0    
    for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
      isospin_component_sectors_lsjt[T0]= basis::RelativeSectorsLSJT(relative_space_lsjt,operator_labels.J0,T0,operator_labels.g0);

    std::cout<<"branching LSJT  "<<std::endl;
    BranchRelativeNLSJT(operator_labels, Nmax, Jmax,relative_space_lsjt,rmes_nlst,
      isospin_component_sectors_lsjt,isospin_component_blocks_lsjt
    );

  }


} // end namespace



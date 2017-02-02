/****************************************************************
  spncci_branching_u3s.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/spncci_branching_u3s.h"

#include <fstream>
#include <iostream>

#include "cppformat/format.h"
#include "mcutils/parsing.h"
#include "spncci/unit_tensor.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // baby SpNCCI
  ////////////////////////////////////////////////////////////////

  BabySpNCCISubspace::BabySpNCCISubspace(const spncci::SpNCCIIrrepFamily& spncci_irrep_family, const sp3r::U3Subspace& u3_subspace)
  {

    // set labels
    labels_ = BabySpNCCISubspaceLabels(
        spncci_irrep_family.sigma(),
        spncci_irrep_family.Sp(),
        spncci_irrep_family.Sn(),
        spncci_irrep_family.S(),
        u3_subspace.U3()  // omega
      );

    // save dimension info
    gamma_max_ = spncci_irrep_family.gamma_max();
    upsilon_max_ = u3_subspace.size();
    dimension_ = gamma_max_*upsilon_max_;
  }

  std::string BabySpNCCISubspace::LabelStr() const
  {
    return fmt::format("{} ({} {} {}) {}",sigma().Str(),Sp(),Sn(),S(),omega().Str());
  }

  std::string BabySpNCCISubspace::DebugStr() const
  {
    return fmt::format("{}: gamma_max {} upsilon_max {} -> dim {}",LabelStr(),gamma_max_,upsilon_max_,size());
  }

  BabySpNCCISpace::BabySpNCCISpace(const spncci::SpNCCISpace& spncci_space)
  {

    // traverse irrep families
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {

        // extract Sp(3,R) space
        const sp3r::Sp3RSpace& sp_space = spncci_irrep_family.Sp3RSpace();

        // traverse U(3) subspaces
        for (int subspace_index = 0; subspace_index < sp_space.size(); ++subspace_index)
          {
            const sp3r::U3Subspace& u3_subspace = sp_space.GetSubspace(subspace_index);
            assert(u3_subspace.size()!=0);
            PushSubspace(BabySpNCCISubspace(spncci_irrep_family,u3_subspace));
          }
      }
  }

  BabySpNCCISectors::BabySpNCCISectors(
      const spncci::BabySpNCCISpace& space,
      const u3shell::OperatorLabelsU3S& operator_labels
    )
  // Based loosely on u3shell::SectorsU3SPN constructor.
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const BabySpNCCISubspace& bra_subspace = space.GetSubspace(bra_subspace_index);
          const BabySpNCCISubspace& ket_subspace = space.GetSubspace(ket_subspace_index);

          // verify selection rules
          bool allowed = true;
          // U(1)
          allowed &= (ket_subspace.omega().N() + operator_labels.N0() - bra_subspace.omega().N() == 0);
          // spin
          //
          // Note: Basic constaints could/should also be placed on Sp
          // and Sn triangularity based on two-body nature of
          // operator, so delta Sp<=2 and delta Sn<=2.  However, in
          // general, the operator does not have sharp Sp0 or Sn0.
          allowed &= am::AllowedTriangle(ket_subspace.S(),operator_labels.S0(),bra_subspace.S());
          // allowed &= am::AllowedTriangle(ket_subspace.Sp(),operator_labels.S0(),bra_subspace.Sp());
          // allowed &= am::AllowedTriangle(ket_subspace.Sn(),operator_labels.S0(),bra_subspace.Sn());
          if (!allowed)
            continue;

          // find SU(3) multiplicity and check SU(3) selection
          int multiplicity = u3::OuterMultiplicity(
              ket_subspace.omega().SU3(),
              operator_labels.x0(),
              bra_subspace.omega().SU3()
            );
          allowed &= (multiplicity > 0);

          // push sectors (tagged by multiplicity)
          if (allowed)
            for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
              {
                PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace,multiplicity_index));
              }
        }
  }

  ////////////////////////////////////////////////////////////////
  // basis indexing in U3S scheme for spncci basis branching
  ////////////////////////////////////////////////////////////////

  SubspaceU3S::SubspaceU3S(const u3::U3S& omegaS, const SpNCCISpace& spncci_space)
  {
    labels_=omegaS;
    int state_index=0;
    int family_index=0;
    for (const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      {
        const int gamma_max = spncci_irrep_family.gamma_max();
        const sp3r::Sp3RSpace& irrep_space = spncci_irrep_family.Sp3RSpace();

        // if space contains omega x S
        if(irrep_space.ContainsSubspace(omegaS.U3()) && spncci_irrep_family.S()==omegaS.S())
        {
          //Construct subspace
          // dim is nu_max*gamma_max (size up subspace)
          // index is starting index in sector matrix
          int dim=gamma_max*irrep_space.LookUpSubspace(omegaS.U3()).size();
          PushStateLabels(StateLabelsType(family_index,spncci_irrep_family.sigma(),dim,state_index));

          // increment indices
          state_index+=dim;
          ++family_index;
        }
      }
    sector_size_=state_index;
  }

  std::string SubspaceU3S::Str() const
  {
    return omegaS().Str();
  }

  std::string StateU3S::Str() const
  {
    return fmt::format("[{} {}]",gamma(),sigma().Str());
  }


  SpaceU3S::SpaceU3S(SpNCCISpace& spncci_space)
  {
    for(const spncci::SpNCCIIrrepFamily& spncci_irrep_family : spncci_space)
      // for each SpNCCI irrep family
      {
        const int gamma_max = spncci_irrep_family.gamma_max();
        const sp3r::Sp3RSpace& irrep_space = spncci_irrep_family.Sp3RSpace();
        HalfInt S = spncci_irrep_family.S();

        // iterate through omega space
        for(int subspace_index=0; subspace_index<irrep_space.size(); ++subspace_index)
          // for each omega within SpNCCI irrep family
          {
            u3::U3 omega(irrep_space.GetSubspace(subspace_index).GetSubspaceLabels());

            // skip if omegaS already in space
            u3::U3S omegaS(omega,S);
            // std::cout<<"Subspace labels "<<omegaS.Str()<<std::endl;
            if (ContainsSubspace(omegaS))
              continue;

            // otherwise construct subspace 
            SubspaceU3S subspace(omegaS,spncci_space);
            // std::cout<<"    Subspace "<<subspace.Str()<<std::endl;
            PushSubspace(subspace);
          }
      }

    // get dimension and starting index of last subspace

    // (mac): why? for total dimension?  will restructure how store
    // subspace dimensions, not as part of labels

    const SubspaceType& subspace=GetSubspace(subspaces_.size()-1);
    int dim=std::get<2>(subspace.GetStateLabels(subspace.size()-1));
    int index=std::get<3>(subspace.GetStateLabels(subspace.size()-1));
    dimension_=dim+index;
  }

  std::string SectorLabelsU3S::Str() const
  {
    return fmt::format("( {} {}  {}{} {} : {} {}  {}", bra_index(),ket_index(), N0(), x0().Str(),S0(),kappa0(),L0(),rho0());
  }

  void GetSectorsU3S(
    const spncci::SpaceU3S& space, 
    const std::vector<u3shell::IndexedOperatorLabelsU3S>& relative_tensor_labels,
    std::vector<spncci::SectorLabelsU3S>& sector_vector
    )
  // TODO: replace with a normal sectors object
  {
    int u3s_sector_vector_index=0;
    for(auto tensor_labels:relative_tensor_labels)
      for(int j=0; j<space.size(); ++j)
        for(int i=0; i<=j; ++i)
          {
            int kappa0,L0; 
            u3shell::OperatorLabelsU3S op_labels;
            std::tie(op_labels, kappa0,L0)=tensor_labels;
            assert(kappa0!=0);
            u3::U3S omegapSp(space.GetSubspace(i).GetSubspaceLabels());
            u3::U3S omegaS(space.GetSubspace(j).GetSubspaceLabels());
            int rho0_max=u3::OuterMultiplicity(omegaS.SU3(), op_labels.x0(),omegapSp.SU3());
            // Check if allowed U(1) coupling
            if(omegaS.U3().N()+op_labels.N0()!=omegapSp.U3().N())
              continue;
            // Check if allowed SU(2) coupling 
            if(not am::AllowedTriangle(omegaS.S(), op_labels.S0(), omegapSp.S()))
              continue;
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                spncci::SectorLabelsU3S sector(i,j,op_labels, kappa0,L0,rho0);
                // should not alread by a sector in map
                // assert(not u3s_sectors.count(sector));
                // if(not u3s_sectors.count(sector))
                //   {
                //     u3s_sectors[sector]=u3s_sector_vector_index;
                //     ++u3s_sector_vector_index;
                //   }
                sector_vector.push_back(sector);
              }
          }
    // sector_vector.resize(u3s_sectors.size());
    // for(auto it=u3s_sectors.begin(); it!=u3s_sectors.end(); ++it)
    //   sector_vector[it->second]=it->first;
  }

  ////////////////////////////////////////////////////////////////
  // Anna's contraction
  ////////////////////////////////////////////////////////////////

  typedef std::map< std::pair<int,int>, std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache >
                      > UnitTensorSectorMatrices;


  void 
  ContractAndRegroupAEM(
      int Nmax, int N1b,
      const std::vector<spncci::SectorLabelsU3S>& sector_labels_vector,
      const u3shell::RelativeRMEsU3ST& interaction_rme_cache,
      const spncci::SpaceU3S& space,
      UnitTensorSectorMatrices& unit_tensor_sector_cache,
      basis::MatrixVector& matrix_vector
    )
  // Args:
  //  Nmax (input) : Basis truncation parameter
  //  N1b (input) : Basis single particle cutoff for Nmax=0
  //  sector_labels_vector (input) : vector of sector labels 
  //  interaction_rme_cache (input) : Container holding interaction rme's keyed
  //     by RelativeUnitTensorU3ST labels
  //  space (input) : space of omegaS subspaces
  //  unit_tensor_sector_cache (input) : nested container holding unit tensor rmes
  //  matrix_vector (output) : vector of U3S sectors indexed by U3S labels and kappa0,L0
  {
    // Contracting interaction matrix elements with unit tensor matrix elements
    std::cout<<"contracting"<<std::endl;
    // Initial sectors to zero matrices
    matrix_vector.resize(sector_labels_vector.size());
    for(int s=0; s<matrix_vector.size(); ++s)
      {
        const spncci::SectorLabelsU3S& sector=sector_labels_vector[s];
        const spncci::SubspaceU3S& ket_subspace=space.GetSubspace(sector.ket_index());
        const spncci::SubspaceU3S& bra_subspace=space.GetSubspace(sector.bra_index());
        int sector_dim_bra=bra_subspace.sector_dim();
        int sector_dim_ket=ket_subspace.sector_dim();
        matrix_vector[s]=Eigen::MatrixXd::Zero(sector_dim_bra,sector_dim_ket);
      }

    // iterate over interaction get unit tensor,kappa0,L0
    // iterate over sectors, get i,j, omega'S', rho0, omegaS->Nn etc. 
    // get corresponding unit tensor sector   
    for(auto it=interaction_rme_cache.begin(); it!=interaction_rme_cache.end(); ++it)
      {
        // Extract labels 
        int kappa0,L0;
        u3shell::RelativeUnitTensorLabelsU3ST tensor_u3st;
        std::tie(tensor_u3st,kappa0,L0)=it->first;
        double interaction_rme=it->second;
        //Check that unit tensor has rme between states in Nmax truncated basis
        int rp=tensor_u3st.bra().eta();
        int r=tensor_u3st.ket().eta();
        if((r>Nmax+2*N1b)||(rp>Nmax+2*N1b))
          continue;
        // Iterate over U3 sectors and find target sectors
        #pragma omp parallel for schedule(runtime)
        for(int s=0; s<sector_labels_vector.size(); ++s)
          {
        const spncci::SectorLabelsU3S& sector=sector_labels_vector[s];
          
        // Checking if sector is target sector
        bool allowed=sector.operator_labels()==u3shell::OperatorLabelsU3S(tensor_u3st.operator_labels());
        allowed&=sector.kappa0()==kappa0;
        allowed&=(sector.L0()==L0);
        if(not allowed)
          continue;
          
        //get subspace labels
        const spncci::SubspaceU3S& ket_subspace=space.GetSubspace(sector.ket_index());
        const spncci::SubspaceU3S& bra_subspace=space.GetSubspace(sector.bra_index());
        const u3::U3& omegap=bra_subspace.GetSubspaceLabels().U3();
        const u3::U3& omega=ket_subspace.GetSubspaceLabels().U3();
        int rho0=sector.rho0();
        // iterate through U3S subspaces, states correspond to lgi sets
        for(int i=0; i<bra_subspace.size(); ++i)
          for(int j=0; j<ket_subspace.size(); ++j)
            {
        // extracting subsector information
        // (dimp,dim)->size of subsector
        // (indexp,index)-> position of upper left corner of subsector in sector
        int dimp, dim, indexp,index,gammap,gamma;
        u3::U3 sigma,sigmap;
        std::tie(gammap,sigmap,dimp,indexp)=bra_subspace.GetStateLabels(i);
        std::tie(gamma,sigma,dim,index)=ket_subspace.GetStateLabels(j);

        //Keys for looking up subsector in unit tensor cache
        std::pair<int,int> lgi_pair(gammap,gamma);
        std::pair<int,int> NnpNn(int(omegap.N()-sigmap.N()),int(omega.N()-sigma.N()));              
        spncci::UnitTensorU3Sector unit_sector(omegap,omega,tensor_u3st,rho0);
        // Get cache containing unit tensor sector 
        if(not unit_tensor_sector_cache.count(lgi_pair))
          continue;
        if(not unit_tensor_sector_cache[lgi_pair].count(NnpNn)) 
          continue;
        spncci::UnitTensorSectorsCache& cache=unit_tensor_sector_cache[lgi_pair][NnpNn];
                
        //REMOVE : turn back on when finished debugging
        // Diagonistic: check if expected unit tensor sectors are found
        // if(not cache.count(unit_sector))
        // for(auto t=cache.begin(); t!=cache.end(); ++t)
        // std::cout<<t->first.Str()<<std::endl;
        #pragma omp critical
        {
          if(cache.count(unit_sector))
            matrix_vector[s].block(indexp,index,dimp,dim)+=interaction_rme*cache[unit_sector];
        }
      }
  }
}
  }


}  // namespace

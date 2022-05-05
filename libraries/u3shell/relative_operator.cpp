/****************************************************************
  relative_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3shell/relative_operator.h"

#include<algorithm>

#include "am/am.h"
#include "fmt/format.h"
#include "sp3rlib/vcs.h"
#include "u3shell/u3st_scheme.h"
#include "u3shell/two_body_operator.h"
#include "boost/math/constants/constants.hpp"

namespace u3shell {

  bool J0Allowed(const u3::SU3& x0, int S0, int J0)
    {
      for(int L0=abs(S0-J0); L0<=(S0+J0); ++L0)
        {
          if(u3::BranchingMultiplicitySO3(x0,L0)>0)
            return true;
        }
      return false;
    }

  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax,
        int N1v,
        std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels,
        int J0,
        int T00,
        bool restrict_positive_N0
        )
  {
    #ifdef VERBOSE
    std::cout<<"Entering GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
    #endif

    bool restrict_J0 = (J0!=-1);
    int eta_max=Nmax+2*N1v;
    int N0_min=restrict_positive_N0?0:-1*eta_max;


    for(int N0=N0_min; N0<=Nmax; N0+=2)
      {
        for( int etap=0; etap<=eta_max; ++etap)
          {
            int eta=etap-N0;
            if((eta<0)||(eta>eta_max))
              continue;
            // Get allowed x0 values
            MultiplicityTagged<u3::SU3>::vector x0_set
              =u3::KroneckerProduct(u3::SU3(etap,0),u3::SU3(0,eta));

            for(int Sp=0; Sp<=1; Sp++)
              for(int Tp=0; Tp<=1; Tp++)
                for(int S=0; S<=1; S++)
                  for (int T=0; T<=1; T++)
                    for (int S0=abs(S-Sp); S0<=(S+Sp); S0++)
                      {
                        //antisymmeterization constraint on ket
                        if ( (etap+Sp+Tp)%2!=1 )
                          continue;
                        //antisymmeterization constraint on bra
                        if ( (eta+S+T)%2!=1)
                          continue;

                      // std::cout<<"hi"<<std::endl;
                        int T0_min=(T00==-1)?abs(Tp-T):T00;
                        int T0_max=(T00==-1)?(Tp+T):T00;
                        // std::cout<<T0_min<<"  "<<T0_max<<std::endl;
                        for(int T0=T0_min; T0<=T0_max; ++T0)
                        {
                          if(not am::AllowedTriangle(T,Tp,T0))
                            continue;

                          u3shell::RelativeStateLabelsU3ST ket(eta,S,T);
                          u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp);
                          // std::cout<<fmt::format("{} {} {}   {} {} {}",etap,Sp,Tp,eta,S,T)<<std::endl;
                          for(int w=0; w<x0_set.size(); w++)
                            {
                              u3::SU3 x0(x0_set[w].irrep);
                              // If restrict on J0 and J0 allowed or not restricted
                              if((restrict_J0 && J0Allowed(x0,S0,J0)) || (not restrict_J0))
                                  relative_unit_tensor_labels.emplace_back(x0,S0,T0,bra,ket);

                              //std::cout<<"unit tensors  "<<spncci::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T).Str()<<std::endl;
                            }
                        }
                      }
          }
      }
  #ifdef VERBOSE
  std::cout<<"Exiting GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
  #endif
  } //end function

  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax,
        int N1v,
        std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels,
        int J0,
        int T00,
        bool restrict_positive_N0
      )
  {
    std::vector<RelativeUnitTensorLabelsU3ST> temp_vector;
    GenerateRelativeUnitTensorLabelsU3ST(Nmax, N1v,temp_vector,J0,T00,restrict_positive_N0);
    for (auto& tensor : temp_vector)
    {
      // std::cout<<"tensor "<<tensor.Str()<<std::endl;
      relative_unit_tensor_labels[tensor.N0()].push_back(tensor);
    }
  }


  double Nrel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    if (bra==ket)
      rme=ket.eta();
    return rme;
  }

  double Hrel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    if (bra==ket)
      rme=ket.eta()+3/2.;
    return rme;
  }

  //RME expression based on McCoy thesis 2018
  double Arel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0.;
    int eta=ket.eta();
    int etap=bra.eta();
    if(
        (bra.S()==ket.S())    // delta on spin
        && (bra.T()==bra.T()) // delta on isospin
        && ((etap-eta)==2)      //only connect states with eta+2=etap
      )
        rme=std::sqrt((eta+2)*(eta+1)/2);

    return rme;
  }

//RME expression based on McCoy thesis 2018
  double Brel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {

    double rme=0.0;
    int eta=ket.eta();
    int etap=bra.eta();
    if((eta==0)||(eta==1))
      return rme;

    if(
        (bra.S()==ket.S())    // delta on spin
        && (bra.T()==bra.T()) // delta on isospin
        && ((eta-etap)==2)      //only connect states with eta+2=etap
      )
        rme=std::sqrt((eta+2)*(eta+1)/2);

    return rme;
  }

  //RME expression based on McCoy thesis 2018
  double Crel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
    {
      double rme=0.0;
      int eta=ket.eta();
      int etap=bra.eta();
      if((eta==0))
        return rme;

      if(
          (bra.S()==ket.S())    // delta on spin
          && (bra.T()==bra.T()) // delta on isospin
          && (eta==etap)      //only connect states with eta+2=etap
        )
          rme=std::sqrt(4.*(eta*eta+3*eta)/3.);

      return rme;

    }

  double K2rel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    if (bra.eta()==ket.eta())
      // 1.5 from the 3/2 zero point energy for a single particle
      rme=u3shell::Nrel(bra,ket)+1.5;
    if (bra.eta()==(ket.eta()+2))
      rme=-sqrt(1.5)*u3shell::Arel(bra,ket);
    if (bra.eta()==(ket.eta()-2))
      rme=-sqrt(1.5)*u3shell::Brel(bra,ket);

    return rme;
  }

  double Qrel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    rme+=std::sqrt(3)*u3shell::Crel(bra,ket);
    rme+=std::sqrt(3)*u3shell::Arel(bra,ket);
    rme+=std::sqrt(3)*u3shell::Brel(bra,ket);
    return rme;
  }


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// New code for updated recurrence
////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace relative
{

RelativeOperator::RelativeOperator(
    const std::vector<OperatorParameters>& parameter_set,
    const std::vector<RMEFunction>& rme_functions,
    const std::vector<double>& coefficients
  )
{
  // Generate sectors
  parameters_ = u3shell::relative::CombineParameters(parameter_set);

  sectors_ = OperatorSectors{
    std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(parameters_),
    std::make_shared<const u3shell::spin::twobody::OperatorSpace>(parameters_),
    parameters_.J0
  };

  // Compute RMEs

  rmes_=std::vector<double>(sectors_.num_elements(), 0.0);
  for (int i = 0; i < rme_functions.size(); ++i)
  {
    auto temp = RelativeOperatorRMEs(sectors_, rme_functions[i], coefficients[i]);
    rmes_ = rmes_
            + RelativeOperatorRMEs(sectors_, rme_functions[i], coefficients[i]);
  }
}

RelativeOperator::RelativeOperator(const std::string& filename, const std::string& filetype)
{
  if(filetype=="text")
  {
    std::string line;
    std::ifstream input(filename);
    mcutils::StreamCheck(bool(input),filename,fmt::format("{} failed to open",filename));

    // skip header
    for(int i=0; i<6; ++i)
      std::getline(input,line);

    parameters_ = u3shell::relative::ReadOperatorParametersText(input);
    sectors_ = ConstructOperatorSectors(parameters_);

    // Read in rmes
    std::size_t index = 0;
    rmes_.resize(sectors_.num_elements());
    while(std::getline(input,line))
    {
      std::istringstream line_stream(line);
      int lambda0,mu0,S0,T0,kappa0,L0,Nbarp,Sbarp,Tbarp,Nbar,Sbar,Tbar;
      double rme;
      line_stream>>lambda0>>mu0>>S0>>T0>>kappa0>>L0>>Nbarp>>Sbarp>>Tbarp>>Nbar>>Sbar>>Tbar>>rme;
      rmes_[index]=rme;
      ++index;
    }
    input.close();
  }
}


  double IdentityRME(
    const TensorLabelsU3ST& tensor_labels,
    const StateLabelsNST& bra,
    const StateLabelsNST& ket,
    const unsigned int kappa0,
    const unsigned int L0
  )
  {
    const auto&[x0,S0,T0] = tensor_labels;
    if(S0==0 && T0==0 && x0==u3::SU3(0u,0u) && kappa0==1 && L0 == 0)
      if(bra[0]==ket[0] && bra[1]==ket[1] && bra[2]==ket[2])
        return 1.0;

    return 0.0;
  }

  double QuadrupoleRME(
    const TensorLabelsU3ST& tensor_labels,
    const StateLabelsNST& bra,
    const StateLabelsNST& ket,
    const unsigned int kappa0,
    const unsigned int L0
  )
  {
    double rme = 0.0;
    const auto&[x0,S0,T0] = tensor_labels;

    if(bra[1]==ket[1] && bra[2]==ket[2] && S0==0 && T0!=2 && kappa0==1 && L0 == 2)
    {
      int T = ket[2];
      double isospin_coefficient=(T0==1)?2*sqrt(T*(T+1)):1;
      if(x0==u3::SU3(2u,0u))
        rme = std::sqrt(3)*spatial::onecoord::Arme(bra[0],ket[0]);
      if(x0==u3::SU3(0u,2u))
        rme = std::sqrt(3)*spatial::onecoord::Brme(bra[0],ket[0]);
      if(x0==u3::SU3(1u,1u))
        rme = std::sqrt(3)*spatial::onecoord::Crme(bra[0],ket[0]);

      double coef = sqrt(5./(16*boost::math::constants::pi<double>()));
      rme*=coef*isospin_coefficient;
    }

    return rme;
  }


  double KSquaredRME(
    const TensorLabelsU3ST& tensor_labels,
    const StateLabelsNST& bra,
    const StateLabelsNST& ket,
    const unsigned int kappa0,
    const unsigned int L0
  )
  {
    double rme = 0.0;
    const auto&[x0,S0,T0] = tensor_labels;

    if(bra[1]==ket[1] && bra[2]==ket[2] && S0==0 && T0==0 && kappa0==1 && L0 == 0)
    {
      if(x0==u3::SU3(2u,0u))
        rme = -std::sqrt(1.5)*spatial::onecoord::Arme(bra[0],ket[0]);
      if(x0==u3::SU3(0u,2u))
        rme = -std::sqrt(1.5)*spatial::onecoord::Brme(bra[0],ket[0]);
      if(x0==u3::SU3(0u,0u))
        rme = spatial::onecoord::Hrme(bra[0],ket[0]);
    }

    return rme;
  }

  double RSquaredRME(
    const TensorLabelsU3ST& tensor_labels,
    const StateLabelsNST& bra,
    const StateLabelsNST& ket,
    const unsigned int kappa0,
    const unsigned int L0
  )
  {
    double rme = 0.0;
    const auto&[x0,S0,T0] = tensor_labels;

    if(bra[1]==ket[1] && bra[2]==ket[2] && S0==0 && T0==0 && kappa0==1 && L0 == 0)
    {
      if(x0==u3::SU3(2u,0u))
        rme = std::sqrt(1.5)*spatial::onecoord::Arme(bra[0],ket[0]);
      if(x0==u3::SU3(0u,2u))
        rme = std::sqrt(1.5)*spatial::onecoord::Brme(bra[0],ket[0]);
      if(x0==u3::SU3(0u,0u))
        rme = spatial::onecoord::Hrme(bra[0],ket[0]);
    }

    return rme;
  }



  std::vector<double> RelativeOperatorRMEs(
    const OperatorSectors& sectors,
    const RMEFunction& rme_function,
    const double coef
  )
  {
    std::vector<double> rme_array(sectors.num_elements(),0.0);
    for (size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
    {
      const auto& sector = sectors.GetSector(sector_index);
      const auto sector_offset = sectors.GetSectorOffset(sector_index);

      const auto& spatial_subspace = sector.bra_subspace();
      const auto& L0 = spatial_subspace.L0();

      const auto& spin_subspace = sector.ket_subspace();
      const auto& S0 = spin_subspace.S0();
      const auto& exchange_symm_bar = spin_subspace.exchange_symm_bar();

      const int parity_bar = (exchange_symm_bar + 1) % 2;

      for (std::size_t spin_state_index = 0; spin_state_index < spin_subspace.size();
           ++spin_state_index)
      {
        const auto& [S00, T0, Sbar, Sbarp, Tbar, Tbarp] =
            spin_subspace.GetState(spin_state_index).labels();

        assert(S00 == S0);

        for(size_t x0subspace_index=0; x0subspace_index<spatial_subspace.size(); ++x0subspace_index)
        {
          const auto& x0_subspace = spatial_subspace.GetSubspace(x0subspace_index);
          const auto& x0 = x0_subspace.x0();
          const auto& N0 = x0_subspace.N0();
          for (unsigned int kappa0 = 1;kappa0 <= spatial_subspace.GetSubspaceDegeneracy(x0subspace_index); ++kappa0)
          {
            const auto x0kappa0_offset =
                spatial_subspace.GetSubspaceOffset(x0subspace_index, kappa0);

            for (std::size_t spatial_state_index=0; spatial_state_index<x0_subspace.size(); ++spatial_state_index)
            {
              const auto& spatial_state = x0_subspace.GetState(spatial_state_index);
              const auto& Nbar = spatial_state.Nbar();
              const auto& Nbarp = spatial_state.Nbarp();

              size_t array_index =
                  sector_offset
                  + sector.element_offset(spin_state_index, x0kappa0_offset, spatial_state_index);

              rme_array[array_index]
                = coef*rme_function({x0, S0, T0}, {Nbarp, Sbarp, Tbarp}, {Nbar, Sbar, Tbar}, kappa0, L0);
            }
          }
        }
      }
    }
    return rme_array;
  }


  void WriteRelativeOperatorText(const RelativeOperator& op, const std::string& filename)
  {
    std::ofstream output(filename);
    // Write header info
    u3shell::relative::WriteOperatorParametersHeader(output);
    output << fmt::format("# lambda0 mu0 S0 T0   kappa0 L0   Nbarp Sbarp Tbarp   Nbar Sbar Tbar  rme\n");

    u3shell::relative::WriteOperatorParameters(op.parameters(), output);
    const auto& sectors = op.sectors();

    for (size_t sector_index = 0; sector_index < sectors.size(); ++sector_index)
    {
      const auto& sector = sectors.GetSector(sector_index);
      const auto sector_offset = sectors.GetSectorOffset(sector_index);

      const auto& spatial_subspace = sector.bra_subspace();
      const auto& L0 = spatial_subspace.L0();

      const auto& spin_subspace = sector.ket_subspace();
      const auto& S0 = spin_subspace.S0();
      const auto& exchange_symm_bar = spin_subspace.exchange_symm_bar();

      const int parity_bar = (exchange_symm_bar + 1) % 2;

      for (std::size_t spin_index = 0; spin_index < spin_subspace.size();++spin_index)
      {
        const auto& [S00, T0, Sbar, Sbarp, Tbar, Tbarp] =
            spin_subspace.GetState(spin_index).labels();

        assert(S00 == S0);

        for(size_t x0subspace_index=0; x0subspace_index<spatial_subspace.size(); ++x0subspace_index)
        {
          const auto& x0_subspace = spatial_subspace.GetSubspace(x0subspace_index);
          const auto& x0 = x0_subspace.x0();
          const auto& N0 = x0_subspace.N0();
          for (unsigned int kappa0 = 1;kappa0 <= spatial_subspace.GetSubspaceDegeneracy(x0subspace_index); ++kappa0)
          {
            const auto x0kappa0_offset = spatial_subspace.GetSubspaceOffset(x0subspace_index, kappa0);

            for (std::size_t spatial_index=0; spatial_index<x0_subspace.size(); ++spatial_index)
            {
              const auto& spatial_state = x0_subspace.GetState(spatial_index);
              const auto& Nbar = spatial_state.Nbar();
              const auto& Nbarp = spatial_state.Nbarp();

              size_t array_index =
                  sector_offset + sector.element_offset(spin_index, x0kappa0_offset, spatial_index);

              double rme = op.rmes()[array_index];

              output << fmt::format("{:2} {:2} {:2} {:2}   {:2} {:2}   {:2} {:2} {:2}   {:2} {:2} {:2}  {:.8e}\n",
                x0.lambda(),x0.mu(),S0,T0,kappa0,L0,Nbarp,Sbarp,Tbarp,Nbar,Sbar,Tbar,rme);
            }
          }
        }
      }
    }

    output.close();

  }

  //NOTE:  Not currently used except in testing
  double RelativeOperator::ReducedMatrixElement(
      const TensorLabelsU3ST& tensor_labels,
      const StateLabelsNST& bra,
      const StateLabelsNST& ket,
      const unsigned int kappa0,
      const unsigned int L0
    )
    {
      const auto& [x0,S0,T0] = tensor_labels;
      const auto& [Nbarp,Sbarp,Tbarp] = bra;
      const auto& [Nbar,Sbar,Tbar] = ket;

      const auto& spatial_space = sectors().bra_space();
      const std::size_t parity_space_index = spatial_space.LookUpSubspaceIndex(Nbar % 2);
      if(parity_space_index!=basis::kNone)
      {
        const auto& parity_space = spatial_space.GetSubspace(parity_space_index);
        const std::size_t N0_space_index = parity_space.LookUpSubspaceIndex(Nbarp-Nbar);
        if(N0_space_index!=basis::kNone)
        {
          const auto& N0_space = parity_space.GetSubspace(N0_space_index);
          const std::size_t L0_subspace_index = N0_space.LookUpSubspaceIndex(L0);
          if(L0_subspace_index!=basis::kNone)
          {
            const auto& L0_space = N0_space.GetSubspace(L0_subspace_index);
            const std::size_t x0_index = L0_space.LookUpSubspaceIndex(x0);
            if(x0_index!=basis::kNone)
            {
              const auto& x0_subspace = L0_space.GetSubspace(x0_index);
              const std::size_t x0_kappa0_index = L0_space.GetSubspaceOffset(x0_index, kappa0);
              const std::size_t Nbar_index = x0_subspace.LookUpStateIndex(Nbar);
              if (Nbar_index != basis::kNone)
              {
                const auto& spin_space = sectors().ket_space();
                const std::size_t S0_subspace_index = spin_space.LookUpSubspaceIndex({(Sbar+Tbar)%2,S0});
                if(S0_subspace_index!=basis::kNone)
                {
                  const auto spin_label_index =
                      u3shell::spin::twobody::LookUpIndex(
                          {S0, T0, Sbar, Sbarp, Tbar, Tbarp}
                        );
                  const auto spin_index = spin_space.GetSubspace(S0_subspace_index).LookUpStateIndex(spin_label_index);
                  const auto sector_index = sectors().LookUpSectorIndex({parity_space_index,N0_space_index,L0_subspace_index,S0_subspace_index});


                  if(sector_index != basis::kNone)
                  {
                    const auto sector_offset =
                        sectors().GetSectorOffset(sector_index);

                    //// Querry Patrick: Get segfault if I do it this way??? Why?
                    // const auto& sector = sectors().GetSector(sector_index);
                    // const auto element_offset = sector.element_offset(
                    //     spin_index, x0_kappa0_index, Nbar_index
                    //   );

                    const auto element_offset = sectors().GetSector(sector_index).element_offset(
                        spin_index, x0_kappa0_index, Nbar_index
                      );

                    return rmes_[sector_offset + element_offset];
                  }
                }
              }
            }
          }
        }
      }
      // If rme not found in sectors, return zero
      std::string Lkappaprint = L0==u3shell::relative::kNone?"":fmt::format("{} {}",kappa0,L0);
      fmt::print("rme for {} {} {} ({} {} {}, {} {} {}) {} not found.  Returning 0.0.\n",
        x0,S0,T0,Nbarp,Sbarp,Tbarp,Nbar,Sbar,Tbar,Lkappaprint
        );
      return 0.0;
    }



  }  // namespace relative


} // namespace

/****************************************************************
  interaction_truncations.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  11/6/17 (aem): Created.
****************************************************************/
#include "u3shell/interaction_truncation.h"
#include "fmt/format.h"

namespace u3shell
{
  void TruncateInteractionNrel(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    int Nrel_max
    )
  // Truncate interaction by Nrel.  Interaction RME is eliminated if either of the bra or ket eta is greater than Nrel_max.
  {
    for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
        {
          u3shell::OperatorLabelsU3ST key; 
          u3shell::RelativeUnitTensorLabelsU3ST tensor;
          int kappa0, L0;
          std::tie(tensor,kappa0,L0)=it->first;
          int Nrel=std::max(tensor.bra().eta(), tensor.ket().eta());

          if(Nrel<=Nrel_max)
            truncated_interaction_u3st[it->first]=it->second;
        }
  }

  void TruncateInteractionN0(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    int N0_max
    )
  // truncate interaction by N0.  If abs(N0) is greater than N0_max, RME is eliminated
  {
    for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
        {
          u3shell::OperatorLabelsU3ST key; 
          u3shell::RelativeUnitTensorLabelsU3ST tensor;
          int kappa0, L0;
          std::tie(tensor,kappa0,L0)=it->first;
          
          if(abs(tensor.N0())<=N0_max)
            truncated_interaction_u3st[it->first]=it->second;
        }
  }

  void TruncateInteractionx0(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    int x0_max
    )
  // truncate interaction by (lambda, mu).  Interaction RME is eliminated if either lambda or mu is greater than LamMu_Max.
  {
    for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
        {
          u3shell::OperatorLabelsU3ST key; 
          u3shell::RelativeUnitTensorLabelsU3ST tensor;
          int kappa0, L0;
          std::tie(tensor,kappa0,L0)=it->first;
          int U3max=std::max(tensor.x0().lambda(), tensor.x0().mu());
          
          if(U3max<=x0_max)
            truncated_interaction_u3st[it->first]=it->second;
        }
  }


  void AccumulateInteractionByU3ST(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    std::unordered_map<u3shell::OperatorLabelsU3ST,double,boost::hash<u3shell::OperatorLabelsU3ST>>& probability_by_u3st
    )
    // Get distribution of interaction of U3ST tensor componenent (N0,x0,S0,T0) normalized to 1. 
    {
      double total;
      for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
        {
          u3shell::RelativeUnitTensorLabelsU3ST tensor;
          int kappa0, L0;
          std::tie(tensor,kappa0,L0)=it->first;

          const u3shell::OperatorLabelsU3ST& operator_labels=tensor.operator_labels();
          int N0=tensor.N0();
          if(N0>=0)
          {
            double squared_amplitude=(it->second)*(it->second);
            probability_by_u3st[operator_labels]+=squared_amplitude;
            total+=squared_amplitude;
          }
        }
      for(auto it=probability_by_u3st.begin(); it!=probability_by_u3st.end(); ++it)
        probability_by_u3st[it->first]*=1./total;

    }


  void TruncateInteractionU3ST(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    double truncation_threshold
    )
  // truncate interaction by individual U3ST relative rmes by eliminating those which fall
  // below a certain percentage of the total interaction decompostion
  // 
  // Truncation is only carried out on the N0>=0 RMEs.  The corresponding N0<0 matrix element
  // is included based on the value of the corresponding N0>0 RME. 
  {
    // accumulate probability by U3ST (N0,x0,S0,T0)
    std::unordered_map<u3shell::OperatorLabelsU3ST,double,boost::hash<u3shell::OperatorLabelsU3ST>> probability_by_u3st;
    u3shell::AccumulateInteractionByU3ST(interaction_u3st, probability_by_u3st);

    std::unordered_set<u3shell::OperatorLabelsU3ST,boost::hash<u3shell::OperatorLabelsU3ST>> allowed_u3st;
    for(auto it=probability_by_u3st.begin(); it!=probability_by_u3st.end(); ++it)
      {
        if(fabs(it->second)>truncation_threshold)
          {
          allowed_u3st.insert(it->first);
          std::cout<<it->first.Str()<<"  "<<it->second<<std::endl;
          }
      }

    // identify allowed tensors based on u3st symmetry
    for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
        {
          u3shell::OperatorLabelsU3ST key; 
          u3shell::RelativeUnitTensorLabelsU3ST tensor;
          int kappa0, L0;
          std::tie(tensor,kappa0,L0)=it->first;
          if(tensor.N0()>=0)
            key=tensor.operator_labels();
          else
            {
              int N0=-1*tensor.N0();
              u3::SU3 x0(u3::Conjugate(tensor.x0()));
              key=u3shell::OperatorLabelsU3ST(N0,x0,tensor.S0(), tensor.T0(), tensor.g0());
            }
          if(allowed_u3st.count(key))
            truncated_interaction_u3st[it->first]=it->second;
        }
  }

  typedef std::tuple<HalfInt,u3::SU3> SpacialSymmetries;
  void AccumulateInteractionByU3(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    std::unordered_map<SpacialSymmetries,double,boost::hash<SpacialSymmetries>>& probability_by_u3
    )
    {
      double total;
      for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
        {
          u3shell::RelativeUnitTensorLabelsU3ST tensor;
          int kappa0, L0;
          std::tie(tensor,kappa0,L0)=it->first;
          u3::SU3 x0=tensor.x0(); 
          int N0=tensor.N0();
          if(N0>=0)
          {
            double squared_amplitude=(it->second)*(it->second);
            probability_by_u3[u3shell::SpacialSymmetries(N0,x0)]+=squared_amplitude;
            total+=squared_amplitude;
          }
        }
      for(auto it=probability_by_u3.begin(); it!=probability_by_u3.end(); ++it)
        probability_by_u3[it->first]*=1./total;

    }

  void TruncateInteractionByU3(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    double truncation_threshold
    )
  {
    std::unordered_map<SpacialSymmetries,double,boost::hash<SpacialSymmetries>> distribution_over_u1su3;
    u3shell::AccumulateInteractionByU3(interaction_u3st, distribution_over_u1su3);
    
    std::unordered_set<u3shell::SpacialSymmetries,boost::hash<u3shell::SpacialSymmetries>> allowed_u1su3;
    for(auto it=distribution_over_u1su3.begin(); it!=distribution_over_u1su3.end(); ++it)
      {
        // std::cout<<it->second<<"  "<<truncation_threshold<<std::endl;
        if(fabs(it->second)>truncation_threshold)
          allowed_u1su3.insert(it->first);
      }

    for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
      {
        u3shell::RelativeUnitTensorLabelsU3ST tensor;
        int kappa0, L0;
        std::tie(tensor,kappa0,L0)=it->first;
        u3::SU3 x0=tensor.x0(); 
        int N0=tensor.N0();
        u3shell::SpacialSymmetries key=N0<0?u3shell::SpacialSymmetries(-N0,u3::Conjugate(x0)):u3shell::SpacialSymmetries(N0,x0);
        if(allowed_u1su3.count(key))
          truncated_interaction_u3st[it->first]=it->second;
      }
  }


  /*void Contribution(
    const u3shell::RelativeRMEsU3ST& interaction_u3st,
    u3shell::RelativeRMEsU3ST& truncated_interaction_u3st,
    double sum_threshhold
    )
  // truncate interaction by N0.  If abs(N0) is greater than N0_max, RME is eliminated
  {
    u3st keys;
    double sum = 0;
    map <u3st, int> u3st_sums;
    for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
        {
          tuple <int, int, int, int> u3group;
          u3shell::OperatorLabelsU3ST key; 
          u3shell::RelativeUnitTensorLabelsU3ST tensor;
          int kappa0, L0;
          std::tie(tensor,kappa0,L0)=it->first;
          u3group = make_tuple(tensor.x0().lambda(), tensor.x0().mu(), tensor.S0(), tensor.T0())
          if(std::find(keys.begin(), keys.end(), u3group) != keys.end()) {*/
            /* keys contains u3group */
            /*keys.push_back(u3group);
            u3st_sums[u3group] = abs(it->second);
          } else {*/
            /* keys does not contain u3group */
            /*u3st_sums[u3group] += abs(it->second);
          } 
          sum += abs(it->second);
        }
  }*/
  

  void GetRelativeRMEsU3ST(
        const std::string& interaction_filename,
        int Nmax, int Jmax,
        u3shell::RelativeRMEsU3ST& interaction_u3st
        )
    {
      // Read in the interaction from file
      basis::RelativeSpaceLSJT relative_space_lsjt;
      std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
      std::array<basis::MatrixVector,3> isospin_component_blocks_lsjt;
      basis::RelativeOperatorParametersLSJT operator_labels;
      
      // Reads in relative operator and fills out isospin_component_sectors_lsjt, 
      // isospin_component_blocks_lsjt, operator_labels and relative_space_lsjt.
      basis::ReadRelativeOperatorLSJT(
        interaction_filename,relative_space_lsjt,operator_labels,
        isospin_component_sectors_lsjt, isospin_component_blocks_lsjt, true
        );

      // upcouple interaction
      std::cout<<"upcoupling "<<std::endl; 
      for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
        {
          u3shell::Upcoupling(
            relative_space_lsjt,
            isospin_component_sectors_lsjt,
            isospin_component_blocks_lsjt,
            operator_labels.J0, operator_labels.g0, T0,Nmax, interaction_u3st);
        }
    }

  void PrintRelativeRMEsU3ST(const u3shell::RelativeRMEsU3ST& interaction_u3st)
  // printing out u3st rmes 
  {
    for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
      {
        int kappa0,L0;
        u3shell::OperatorLabelsU3ST op;
        std::tie(op,kappa0,L0)=it->first;
        std::cout<<fmt::format("{:20} {} {}  {:20.8f} ",op.Str(),kappa0,L0,it->second)<<std::endl;
      }
  }

  void WriteRelativeRMEsU3ST(const std::string& filename, const u3shell::RelativeRMEsU3ST& interaction_u3st)
  // writing out u3st rmes 
  {
    std::ofstream os(filename);
    os << "#   N0 lambda mu   S0 T0 g0   kappa0 L0   RME" << std::endl;
    for(auto it=interaction_u3st.begin(); it!=interaction_u3st.end(); ++it)
      {
        int kappa0,L0;
        u3shell::OperatorLabelsU3ST op;
        std::tie(op,kappa0,L0)=it->first;
        //std::cout<<fmt::format("{:20} {} {}  {:20.8f} ",op.Str(),kappa0,L0,it->second)<<std::endl;
        const int width=3;
        const int precision=16;
        os << std::setprecision(precision);
        os
          << " " << std::setw(width) << op.N0()
          << " " << std::setw(width) << op.x0().lambda()
          << " " << std::setw(width) << op.x0().mu()
          << " " << "  "
          << " " << std::setw(width) << op.S0()
          << " " << std::setw(width) << op.T0()
          << " " << std::setw(width) << op.g0()
          << " " << "  "
          << " " << std::setw(width) << kappa0
          << " " << std::setw(width) << L0
          << " " << "  "
          << " " << std::showpoint << std::scientific << it->second
          << std::endl;
      }
  }
} // end namespace



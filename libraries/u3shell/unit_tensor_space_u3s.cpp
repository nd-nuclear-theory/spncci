/****************************************************************
  unit_tensor_space_u3s.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/unit_tensor_space_u3s.h"

#include <sstream>

#include "am/am.h"
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "fmt/format.h"



namespace u3shell {

  RelativeUnitTensorSubspaceU3S::RelativeUnitTensorSubspaceU3S(
    u3::SU3 x0, HalfInt S0, int etap, int eta,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels
    )
  {
    // set values
    labels_ = SubspaceLabelsType(x0,S0,etap,eta);
    // std::cout<<fmt::format("{} {} {} {}",x0.Str(),S0,etap,eta)<<std::endl;
    for(auto& tensor : unit_tensor_labels)
      {
        if(
            (x0==tensor.x0())
            &&(S0==int(tensor.S0()))
            &&(etap==tensor.bra().eta())
            &&(eta==tensor.ket().eta())
          )
          {

            int T0=int(tensor.T0());
            int Sp=int(tensor.bra().S());
            int Tp=int(tensor.bra().T());
            int S=int(tensor.ket().S());
            int T=int(tensor.ket().T());

            // std::cout<<tensor.Str()<<std::endl;
            // std::cout<<fmt::format("{} {} {} {} {}",T0,Sp,Tp,S,T)<<std::endl;

            PushStateLabels(StateLabelsType(T0,Sp,Tp,S,T));
          }
      }
  }

  std::string RelativeUnitTensorSubspaceU3S::Str() const
  {

    return fmt::format(
        "[{} {} {} {}]",
        x0().Str(),S0(),etap(),eta()
      );
  }

  std::string RelativeUnitTensorSubspaceU3S::LabelStr() const
  {

    return Str();
  }

  RelativeUnitTensorSpaceU3S::RelativeUnitTensorSpaceU3S(
    int Nmax, int N1v,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels
    )
  {
    // save Nmax
    Nmax_ = Nmax;
    N1v_=N1v;
    int eta_max=Nmax+2*N1v;
    // for each N in 0..Nmax
    for (int N0=-Nmax; N0<=Nmax; N0+=2)
      // for each etap
      for(int etap=0; etap<=eta_max; ++etap)
        {
          int eta=etap-N0;
          if((eta<0)||(eta>eta_max))
            continue;
          MultiplicityTagged<u3::SU3>::vector
            x0_set=u3::KroneckerProduct(u3::SU3(etap,0), u3::SU3(0,eta));

          for(auto& x0_tagged : x0_set)
            {
              u3::SU3 x0(x0_tagged.irrep);
             // for each S0 in 0..2
              for (int S0=0; S0<=2; ++S0)
                {
                  // construct subspace
                  RelativeUnitTensorSubspaceU3S
                    subspace(x0,S0,etap,eta,unit_tensor_labels);
                  // push subspace if nonempty
                  if (subspace.size()!=0)
                    PushSubspace(subspace);
                }
            }
        }
  }


  RelativeUnitTensorSpaceU3S::RelativeUnitTensorSpaceU3S(
    int Nmax, int N1v,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
    const std::map< std::pair<int,int>, UnitTensorSubspaceLabelsSet>&
        NnpNn_organized_unit_tensor_subspaces
    )
  {
    // save Nmax
    Nmax_ = Nmax;
    N1v_=N1v;

    for(auto it=NnpNn_organized_unit_tensor_subspaces.begin(); it!=NnpNn_organized_unit_tensor_subspaces.end(); ++it)
      {
        const std::pair<int,int>& NnpNn=it->first;
        const UnitTensorSubspaceLabelsSet& unit_tensor_subspace_labels_set=it->second;
        for(auto& unit_tensor_subspace_labels : unit_tensor_subspace_labels_set)
          {
            u3::SU3 x0;
            HalfInt S0;
            int etap,eta;
            std::tie(x0,S0,etap,eta)=unit_tensor_subspace_labels;
            // construct subspace
            RelativeUnitTensorSubspaceU3S subspace(x0,S0,etap,eta,unit_tensor_labels);
            // push subspace if nonempty
            if (subspace.size()!=0)
              PushSubspace(subspace);
          }
      }
  }

  std::string RelativeUnitTensorSpaceU3S::Str() const
  {

    std::ostringstream os;
    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);

        std::string subspace_string
          = fmt::format(
              "index {:3} {} size {:3}",
              subspace_index,
              subspace.Str(),
              subspace.size()
            );

        os << subspace_string << std::endl;
      }
    return os.str();
  }

   ObservableSubspaceU3S::ObservableSubspaceU3S(
    int N0, u3::SU3 x0, HalfInt S0, int kappa0, int L0
    )
  {
    // set values
    labels_ = SubspaceLabelsType(N0,x0,S0,kappa0,L0);
    PushStateLabels(1);

  }


  std::string ObservableSubspaceU3S::Str() const
  {

    return fmt::format(
        "[{} {} {} {} {}]",
        N0(),x0().Str(),S0(),kappa0(),L0()
      );
  }

  std::string ObservableSubspaceU3S::LabelStr() const
  {

    return Str();
  }

  ObservableSpaceU3S::ObservableSpaceU3S(
    const std::vector<u3shell::IndexedOperatorLabelsU3S>& observable_labels
  )
  {

    for(auto& tensor : observable_labels)
      {
        u3shell::OperatorLabelsU3S u3s_labels;
        int kappa0,L0;
        std::tie(u3s_labels,kappa0,L0)=tensor;

        // Extract labels
        int N0=u3s_labels.N0();
        u3::SU3 x0=u3s_labels.x0();
        HalfInt S0=u3s_labels.S0();

        ObservableSubspaceU3S
          subspace(N0,x0,S0,kappa0,L0);

        // push subspace if nonempty
        if (subspace.size()!=0)
          PushSubspace(subspace);
      }

  }


  std::string ObservableSpaceU3S::Str() const
  {

    std::ostringstream os;
    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);

        std::string subspace_string
          = fmt::format(
              "index {:3} {} size {:3}",
              subspace_index,
              subspace.Str(),
              subspace.size()
            );

        os << subspace_string << std::endl;
      }
    return os.str();
  }
}

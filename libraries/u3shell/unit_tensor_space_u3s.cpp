/****************************************************************
  unit_tensor_space_u3s.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/

#include "u3shell/unit_tensor_space_u3s.h"

#include <sstream>

#include "am/am.h"
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "fmt/format.h"



namespace u3shell {

  RelativeUnitTensorSubspaceU3S::RelativeUnitTensorSubspaceU3S(
    u3::SU3 x0, HalfInt S0, unsigned int etap, unsigned int eta,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels
    )
      : BaseSubspace{SubspaceLabelsType(x0,S0,etap,eta)}
  {
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

//********************************************** Added by J.H. *********************************************
  OneBodyUnitTensorSubspaceU3S::OneBodyUnitTensorSubspaceU3S(
    u3::SU3 x0, HalfInt S0, unsigned int etap, unsigned int eta,
    const std::vector<u3shell::OneBodyUnitTensorLabelsU3S>& one_body_unit_tensor_labels
    )
      : BaseSubspace{SubspaceLabelsType(x0,S0,etap,eta)}
  {
    for(auto& tensor : one_body_unit_tensor_labels)
      {
        if(
            (x0==tensor.x0())
            &&(S0==int(tensor.S0()))
            &&(etap==tensor.Nbra())
            &&(eta==tensor.Nket())
          )
          {
            int Tz=int(tensor.Tz());
            PushStateLabels(StateLabelsType(Tz));
          }
      }
  }

  OneBodyUnitTensorSpaceU3S::OneBodyUnitTensorSpaceU3S(
    int Nmax, int N1vp, int N1vn,
    const std::vector<u3shell::OneBodyUnitTensorLabelsU3S>& one_body_unit_tensor_labels
    )
  {
    // save Nmax
    Nmax_ = Nmax;
    N1vp_=N1vp;
    N1vn_=N1vn;
    int eta_max=Nmax+std::max(N1vp,N1vn);
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
             // for each S0 in 0..1
              for (int S0=0; S0<=1; ++S0)
                {
                  // construct subspace
                  OneBodyUnitTensorSubspaceU3S
                    subspace(x0,S0,etap,eta,one_body_unit_tensor_labels);
                  // push subspace if nonempty
                  if (subspace.size()!=0)
                    PushSubspace(subspace);
                }
            }
        }
  }

  TwoBodyDensitySubspace::TwoBodyDensitySubspace(
    u3::SU3 x0, HalfInt S0, unsigned int N1, unsigned int N2, unsigned int N3, unsigned int N4,
    const std::vector<u3shell::TwoBodyDensityLabels>& two_body_density_labels
    )
      : BaseSubspace{SubspaceLabelsType(x0,S0,N1,N2,N3,N4)}
  {
    for(auto& tensor : two_body_density_labels)
      {
        if(
            (x0==tensor.x0())
            &&(S0==int(tensor.S0()))
            &&(N1==tensor.N1())
            &&(N2==tensor.N2())
	    &&(N3==tensor.N3())
            &&(N4==tensor.N4())
          )
          {
	    u3::SU3 xf(tensor.xf());
	    int Sf=tensor.Sf();
	    u3::SU3 xi(tensor.xi());
            int Si=tensor.Si();
	    int rho0=tensor.rho0();
            int Tz=tensor.Tz();
            PushStateLabels(StateLabelsType(xf,Sf,xi,Si,rho0,Tz));
          }
      }
  }

  TwoBodyDensitySpace::TwoBodyDensitySpace(
    int Nmax, int N1vp, int N1vn,
    const std::vector<u3shell::TwoBodyDensityLabels>& two_body_density_labels
    )
  {
    // save Nmax
    Nmax_ = Nmax;
    N1vp_=N1vp;
    N1vn_=N1vn;
    int eta_max=Nmax+std::max(N1vp,N1vn);
    for(int N0=-Nmax; N0<=Nmax; N0+=2){
      // for each N1,N2,N3
      for(int N1=0; N1<=eta_max; ++N1){
        for(int N2=0; N2<=eta_max; ++N2){
	  for(int N3=0; N3<=eta_max; ++N3){
            int N4=N1+N2-N3-N0; // N0=N1+N2-N3-N4
            if((N4<0)||(N4>eta_max))continue;
            for(MultiplicityTagged<u3::SU3>& xf_tagged : u3::KroneckerProduct(u3::SU3(N1,0),u3::SU3(N2,0))){
              for(MultiplicityTagged<u3::SU3>& xi_tagged : u3::KroneckerProduct(u3::SU3(0,N3),u3::SU3(0,N4))){
                for(MultiplicityTagged<u3::SU3>& x0_tagged : u3::KroneckerProduct(xf_tagged.irrep,xi_tagged.irrep)){
                  u3::SU3 x0(x0_tagged.irrep);
                  // for each S0 in 0..2
                  for(int S0=0; S0<=2; ++S0){
                    // construct subspace
                    TwoBodyDensitySubspace subspace(x0,S0,N1,N2,N3,N4,two_body_density_labels);
                    // push subspace if nonempty
                    if(subspace.size()!=0)PushSubspace(subspace);
                  }
	        }
	      }
            }
	  }
	}
      }
    }
  }
//**********************************************************************************************************

   ObservableSubspaceU3S::ObservableSubspaceU3S(
    int N0, u3::SU3 x0, HalfInt S0, int kappa0, int L0
    )
      : BaseSubspace{SubspaceLabelsType(N0,x0,S0,kappa0,L0)}
  {
    // set values
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

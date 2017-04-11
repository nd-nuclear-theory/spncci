/****************************************************************
  unit_tensor_space_u3s.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/unit_tensor_space_u3s.h"

#include <sstream>

#include "am/am.h"
#include "cppformat/format.h"



namespace u3shell {

  RelativeUnitTensorSubspaceU3S::RelativeUnitTensorSubspaceU3S(
    u3::SU3 x0, HalfInt S0, int etap, int eta,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels
    )
  {
    // set values
    labels_ = SubspaceLabelsType(x0,S0,etap,eta);
    std::cout<<fmt::format("{} {} {} {}",x0.Str(),S0,etap,eta)<<std::endl;
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
}
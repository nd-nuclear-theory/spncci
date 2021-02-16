/****************************************************************
  relative_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3shell/relative_operator.h"

#include "am/am.h"
#include "fmt/format.h"
#include "sp3rlib/vcs.h"
#include "u3shell/u3st_scheme.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/moshinsky_xform.h"
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
    int N0_min=restrict_positive_N0?0:-1*Nmax;
    int eta_max=Nmax+2*N1v;

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
      if((eta==0)||(eta==1))
        return rme;
      
      if(
          (bra.S()==ket.S())    // delta on spin
          && (bra.T()==bra.T()) // delta on isospin
          && (eta==etap)      //only connect states with eta+2=etap 
        )
          rme=std::sqrt(4.*(eta*eta+3*eta)/3);
      
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



} // namespace

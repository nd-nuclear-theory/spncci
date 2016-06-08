/****************************************************************
  relative_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "am/am.h"
#include "cppformat/format.h"
#include "sp3rlib/vcs.h"
#include "u3shell/relative_operator.h"

namespace u3shell {


  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax, 
        std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels
        )
  {   
    #ifdef VERBOSE
    std::cout<<"Entering GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
    #endif
    for(int N0=0; N0<=Nmax; N0+=2)
      {
        std::vector<RelativeUnitTensorLabelsU3ST> sym_vec;
        for(int Sp=0; Sp<=1; Sp++)
          for(int Tp=0; Tp<=1; Tp++)
            for(int S=0; S<=1; S++)
              for (int T=0; T<=1; T++)
                for (int S0=abs(S-Sp); S0<=(S+Sp); S0++)
                  for (int T0=abs(T-Tp); T0<=(T+Tp); T0++)
                    for(int etap=0; etap<=N0+Nmax; etap++)
                      {
                        //antisymmeterization constraint on ket 
                        if ( (etap+Sp+Tp)%2!=1 )
                          continue;
                        
                        int eta=etap-N0;
                        //antisymmeterization constraint on bra 
                        if ( (eta+S+T)%2!=1)
                          continue;

                        u3shell::RelativeStateLabelsU3ST ket(eta,S,T);
                        u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp);

                        MultiplicityTagged<u3::SU3>::vector omega0_set
                          =u3::KroneckerProduct(u3::SU3(etap,0),u3::SU3(0,eta));

                        for(int w=0; w<omega0_set.size(); w++)
                          {
                            u3::SU3 x0(omega0_set[w].irrep);
                            sym_vec.push_back(u3shell::RelativeUnitTensorLabelsU3ST(x0,S0,T0,bra,ket));
                            //std::cout<<"unit tensors  "<<spncci::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T).Str()<<std::endl;
                          }
                      }       
        relative_unit_tensor_labels[N0]=sym_vec;
      }
  #ifdef VERBOSE
  std::cout<<"Exiting GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
  #endif
  } //end function


 

  double RelativeNumberOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    if (bra==ket)
      rme=ket.eta();
    return rme;
  }

  double RelativeSp3rRaisingOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    int eta=ket.eta();
    int etap=bra.eta();
    if(
        (bra.S()==ket.S())    // delta on spin
        && (bra.T()==bra.T()) // delta on isospin
        && (etap-eta==2)      //only connect states with eta+2=etap 
      )
        {
          u3::U3 sigma(HalfInt(eta%2)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omega(HalfInt(eta)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omegap(HalfInt(etap)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          MultiplicityTagged<u3::U3>n_rho(u3::U3(eta,0,0),1);
          MultiplicityTagged<u3::U3>np_rhop(u3::U3(etap,0,0),1);

          rme=(sqrt(vcs::Omega(np_rhop.irrep, omegap)-vcs::Omega(n_rho.irrep, omega))
                    *vcs::U3BosonCreationRME(sigma,np_rhop,omegap,sigma,n_rho,omega));
          // std::cout
          // <<sqrt(vcs::Omega(np_rhop.irrep, omegap)-vcs::Omega(n_rho.irrep, omega))
          // <<"  "<<vcs::U3BosonCreationRME(sigma,np_rhop,omegap,sigma,n_rho,omega)
          // <<"  "<<rme
          // <<std::endl;
        }
    // std::cout<<rme<<std::endl;    
    return rme;
  }

  double RelativeSp3rLoweringOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    int eta=ket.eta();
    int etap=bra.eta();
    if(
        (bra.S()==ket.S())    // delta on spin
        && (bra.T()==bra.T()) // delta on isospin
        && (eta-etap==2)      //only connect states with eta+2=etap 
      )
        {
          u3::U3 sigma(HalfInt(eta%2)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omega(HalfInt(eta)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omegap(HalfInt(etap)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          MultiplicityTagged<u3::U3>n_rho(u3::U3(eta,0,0),1);
          MultiplicityTagged<u3::U3>np_rhop(u3::U3(etap,0,0),1);

          rme=parity(u3::ConjugationGrade(omegap))//.lambda()+omegap.mu+omega.lambda+omega.mu)
                    *sqrt(1.*u3::dim(omega)/u3::dim(omegap))
                    *sqrt(vcs::Omega(n_rho.irrep, omega)-vcs::Omega(np_rhop.irrep, omegap))
                    *vcs::U3BosonCreationRME(sigma,n_rho,omega,sigma,np_rhop,omegap);
        }
    return rme;
  }

  double RelativeKineticEnergyOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {
    double rme=0;
    if (bra.eta()==ket.eta())
      rme=u3shell::RelativeNumberOperator(bra,ket);
    if (bra.eta()==(ket.eta()+2))
      rme=sqrt(1.5)*u3shell::RelativeSp3rRaisingOperator(bra,ket);
    if (bra.eta()==(ket.eta()-2))
      rme=sqrt(1.5)*u3shell::RelativeSp3rLoweringOperator(bra,ket);

    return rme;
  }

} // namespace

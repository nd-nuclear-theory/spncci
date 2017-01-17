/****************************************************************
  relative_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "u3shell/relative_operator.h"

#include "am/am.h"
#include "cppformat/format.h"
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
        std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels,
        int J0,
        int T00,
        bool restrict_positive_N0
        )
  {   
    #ifdef VERBOSE
    std::cout<<"Entering GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
    #endif
    
    bool restrict_J0=true;
    if(J0==-1)
      restrict_J0=false;

    int N0_min=restrict_positive_N0?0:-1*Nmax;
    for(int N0=N0_min; N0<=Nmax; N0+=2)
      {
        // std::cout<<" N0="<<N0<<std::endl;
        for(int Sp=0; Sp<=1; Sp++)
          for(int Tp=0; Tp<=1; Tp++)
            for(int S=0; S<=1; S++)
              for (int T=0; T<=1; T++)
                for (int S0=abs(S-Sp); S0<=(S+Sp); S0++)
                  for(int etap=0; etap<=Nmax; etap++)
                    {
                      // std::cout<<"hi"<<std::endl;
                      int T0_min=(T00==-1)?abs(Tp-T):T00;
                      int T0_max=(T00==-1)?(Tp+T):T00;
                      // std::cout<<T0_min<<"  "<<T0_max<<std::endl;
                      for(int T0=T0_min; T0<=T0_max; ++T0)
                      {
                        if(not am::AllowedTriangle(T,Tp,T0))
                          continue;
                        //antisymmeterization constraint on ket 
                        if ( (etap+Sp+Tp)%2!=1 )
                          continue;  
                        int eta=etap-N0;
                        if((eta<0)||(eta>Nmax))
                          continue;
                        //antisymmeterization constraint on bra 
                        if ( (eta+S+T)%2!=1)
                          continue;

                        u3shell::RelativeStateLabelsU3ST ket(eta,S,T);
                        u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp);
                        MultiplicityTagged<u3::SU3>::vector omega0_set
                          =u3::KroneckerProduct(u3::SU3(etap,0),u3::SU3(0,eta));
                        // std::cout<<fmt::format("{} {} {}   {} {} {}",etap,Sp,Tp,eta,S,T)<<std::endl;
                        for(int w=0; w<omega0_set.size(); w++)
                          {
                            u3::SU3 x0(omega0_set[w].irrep);
                            // If restrict on J0 and J0 allowed or not restricted
                            if((restrict_J0 && J0Allowed(x0,S0,J0)) || (not restrict_J0))
                                relative_unit_tensor_labels.emplace_back(x0,S0,T0,bra,ket);
                            
                            //std::cout<<"unit tensors  "<<spncci::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T).Str()<<std::endl;
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
        std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels,
        int J0,
        int T0,
        bool restrict_positive_N0
      )
  {
    std::vector<RelativeUnitTensorLabelsU3ST> temp_vector;
    GenerateRelativeUnitTensorLabelsU3ST(Nmax, temp_vector,J0,T0,restrict_positive_N0);
    for (auto tensor : temp_vector)
    {
      // std::cout<<"tensor "<<tensor.Str()<<std::endl;
      relative_unit_tensor_labels[tensor.N0()].push_back(tensor);
    }
  }



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
        && ((etap-eta)==2)      //only connect states with eta+2=etap 
      )
        {
          u3::U3 sigma(HalfInt(eta%2)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omega(HalfInt(eta)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omegap(HalfInt(etap)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          MultiplicityTagged<u3::U3>n_rho(u3::U3(omega.N()-sigma.N(),0,0),1);
          MultiplicityTagged<u3::U3>np_rhop(u3::U3(omegap.N()-sigma.N(),0,0),1);
          rme=(sqrt(vcs::Omega(np_rhop.irrep, omegap)-vcs::Omega(n_rho.irrep, omega))
                    *vcs::U3BosonCreationRME(sigma,np_rhop,omegap,sigma,n_rho,omega));
        }
    return rme;
  }

  double RelativeSp3rRaisingPolynomial(
    const u3::U3& n0, 
    const u3shell::RelativeStateLabelsU3ST& bra, 
    const u3shell::RelativeStateLabelsU3ST& ket, 
    sp3r::BCoefCache& bcoef_cache
    )
  {
    if(n0==u3::U3(0,0,0))
      return 1.0;

    double rme=0;
    int eta=ket.eta();
    int etap=bra.eta();
    if(
        (bra.S()==ket.S())    // delta on spin
        && (bra.T()==bra.T()) // delta on isospin
        && (eta+n0.N()==etap) // U(1) product rule
      )
        {
          u3::U3 sigma(HalfInt(eta%2)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omega(HalfInt(eta)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omegap(HalfInt(etap)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 n(omega.N()-sigma.N(),0,0);
          u3::U3 np(omegap.N()-sigma.N(),0,0);
          MultiplicityTagged<u3::U3>n_rho(n,1);
          MultiplicityTagged<u3::U3>np_rhop(np,1);
          rme=sqrt(vcs::Omega(np, omegap)-vcs::Omega(n, omega))
                *u3::U(n0.SU3(), n.SU3(),omegap.SU3(),sigma.SU3(),np.SU3(),1,1,omega.SU3(),1,1)
                *bcoef_cache[sp3r::BCoefLabels(n0,n,np,1)];
        }
    return rme;
  }




  double RelativeSp3rLoweringOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)
  {

    double rme=0;
    int eta=ket.eta();
    int etap=bra.eta();
    if((eta==0)||(eta==1))
      return rme;
    
    if(
        (bra.S()==ket.S())    // delta on spin
        && (bra.T()==bra.T()) // delta on isospin
        && ((eta-etap)==2)      //only connect states with eta+2=etap 
      )
        {
          u3::U3 sigma(HalfInt(eta%2)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omega(HalfInt(eta)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          u3::U3 omegap(HalfInt(etap)+HalfInt(1,2),HalfInt(1,2),HalfInt(1,2));
          MultiplicityTagged<u3::U3>n_rho(u3::U3(omega.N()-sigma.N(),0,0),1);
          MultiplicityTagged<u3::U3>np_rhop(u3::U3(omegap.N()-sigma.N(),0,0),1);
          rme=parity(etap-eta)
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
      // 1.5 from the 3/2 zero point energy for a single particle
      rme=u3shell::RelativeNumberOperator(bra,ket)+1.5;
    if (bra.eta()==(ket.eta()+2))
      rme=-sqrt(1.5)*u3shell::RelativeSp3rRaisingOperator(bra,ket);
    if (bra.eta()==(ket.eta()-2))
      rme=-sqrt(1.5)*u3shell::RelativeSp3rLoweringOperator(bra,ket);

    return rme;
  }

} // namespace

/****************************************************************
  moshinsky.cpp
                                
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/moshinsky.h"

#include <cmath>

namespace u3shell
{
  double MoshinskyCoefficient(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& xr,const u3::SU3& xc,const u3::SU3& x)
  // SU(3) Moshinsky Coefficient which is equivalent to a Wigner little d function evaluated at pi/2
  {
		
    HalfInt J(x.lambda(),2);
    HalfInt Mp(x1.lambda()-x2.lambda(),2);
    HalfInt M(xr.lambda()-xc.lambda(),2); 
	

    double moshinsky_coef=0;
    int Kmax=std::max(
                      std::max(int(J+M),int(J-M)),int(J-Mp));
		
    for(int K=0; K<=Kmax; K++)
      moshinsky_coef=moshinsky_coef+ParitySign(K)*Choose(int(J+M),int(J-Mp-K))*Choose(int(J-M),K);
    moshinsky_coef=moshinsky_coef*(ParitySign(J-Mp)
                                   *sqrt(Factorial(int(J+Mp))*Factorial(int(J-Mp))
                                         /(pow(2.,TwiceValue(J))*Factorial(int(J+M))*Factorial(int(J-M)))
                                        ));
    return moshinsky_coef;
  }

  double MoshinskyCoefficient(const u3::U3& w1, const u3::U3& w2, const u3::U3& wr,const u3::U3& wc,const u3::U3& w)
  //Overleading for U3 arguements 
  {
    return MoshinskyCoefficient(w1.SU3(), w2.SU3(), wr.SU3(), wc.SU3(), w.SU3());
  }

  double MoshinskyCoefficient(int r1, int r2, int r, int R, const u3::SU3& x)
  //Overloading Moshinsky to take integer arguements for two-body and relative-center of mass arguments
  // and SU(3) total symmetry (lambda,mu)
  {
    HalfInt J(x.lambda(),2);
    HalfInt Mp(r1-r2,2);
    HalfInt M(r-R,2) ;

    double moshinsky_coef=0;
    int Kmax=std::max(std::max(int(J+M),int(J-M)),int(J-Mp));

    for(int K=0; K<=Kmax; K++)
      moshinsky_coef=moshinsky_coef+ParitySign(K)*Choose(int(J+M),int(J-Mp-K))*Choose(int(J-M),K);		
			
    moshinsky_coef=moshinsky_coef*(
                                   ParitySign(J-Mp)
                                   *sqrt(Factorial(int(J+Mp))*Factorial(int(J-Mp))
                                         /(pow(2.,TwiceValue(J))*Factorial(int(J+M))*Factorial(int(J-M)))
                                         ));
    return moshinsky_coef;	
  }

  double MoshinskyCoefficient(int r1, int r2, int r, int R, const u3::U3& w)
  // Overloading Moshinsky to take integers and U3 for total symmetry
  {
    return MoshinskyCoefficient(r1, r2, r, R, w.SU3());
  }

  void MoshinskyTransformation(const u3shell::RelativeUnitTensorLabelsU3ST& tensor, int Nmax)
  {
    //extract operator labels 
    int etap=tensor.bra().eta();
    int eta=tensor.ket().eta();
    u3::SU3 x0(tensor.x0());

    // Generate a space 
    u3shell::TwoBodySpaceU3ST space(Nmax);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sectors are generated subject to the constraints that 
    // omega.N()+N0=omegap.N()
    // omega x x0 -> omegap
    // S x S0 ->Sp
    // T x T0 ->Tp
    u3shell::TwoBodySectorsU3ST sector_labels_list(space,tensor);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Iterate over sectors
    for (int sector_index=0; sector_index<sector_labels_list.size(); sector_index++)
    {
      // for a given sector, extract bra and ket subspace information
      auto sector_labels=sector_labels_list.GetSector(sector_index);
      
      auto bra(sector_labels.bra_subspace());
      u3::U3 omegap(bra.omega());
      int Np=int(omegap.N());
      
      auto ket(sector_labels.ket_subspace());
      u3::U3 omega(ket.omega());
      int N=int(omega.N()); //for Two-body state, N must be int. A/2=2/2

      int rho0=sector_labels.multiplicity_index();
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // eta_cm is fixed by unit tensor plus constraints
      //    (eta,0)x(eta_cm)->omega and eta+eta_cm=N
      //    (etap,0)x(eta_cm)->omegap and etap+eta_cm=Np
      //    and eta_cm>0
      int eta_cm=Np-etap;
      /// Checking valid value of eta_cm
      if (
        eta_cm!=(N-eta)
        ||(eta_cm<0)
        ||(u3::OuterMultiplicity(u3::SU3(eta,0),u3::SU3(eta_cm,0),omega.SU3())==0)
        ||(u3::OuterMultiplicity(u3::SU3(etap,0),u3::SU3(eta_cm,0),omegap.SU3())==0)
        )
        continue;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //bra subspace
      const u3shell::TwoBodySubspaceU3ST& bra_subspace = space.GetSubspace(sector_labels.bra_subspace_index());
      //bra_subspace
      const u3shell::TwoBodySubspaceU3ST& ket_subspace = space.GetSubspace(sector_labels.ket_subspace_index());
      // dimension of the subspaces 
      int dimb=bra_subspace.size();
      int dimk=ket_subspace.size();
      // set up two-body sector for given operator,(omegap,Sp,Tp),(omega,S,T), rho0
      Eigen::MatrixXd sector(dimk,dimb); 

      Eigen::MatrixXd bra_moshinky(dimb,1);
      Eigen::MatrixXd ket_moshinky(1,dimk);
      
      // iterate over bra subspace
      for (int bra_state_index=0; bra_state_index<bra_subspace.size(); ++bra_state_index)
        {
          const u3shell::TwoBodyStateU3ST bra_state(bra_subspace,bra_state_index);
          int eta1p=bra_state.N1();
          int eta2p=bra_state.N2();
          bra_moshinky(bra_state_index,0)=MoshinskyCoefficient(eta1p, eta2p, etap, eta_cm, omegap.SU3());
        } 

      // iterate over ket subspace
      for (int ket_state_index=0; ket_state_index<ket_subspace.size(); ++ket_state_index)
        {
          const u3shell::TwoBodyStateU3ST ket_state(ket_subspace,ket_state_index);
          int eta1=ket_state.N1();
          int eta2=ket_state.N2();

          ket_moshinky(0,ket_state_index)=MoshinskyCoefficient(eta1, eta2, eta, eta_cm, omega);
        }  
      // summing sector terms over Ncm
      sector=
                    u3::U(x0,u3::SU3(eta,0),omegap.SU3(),u3::SU3(eta_cm,0),u3::SU3(etap,0),1,1,omega.SU3(),1,rho0)
                    *(bra_moshinky*ket_moshinky);

      // Storing sector in matrix  
      std::cout<<bra.Str()<<"||"<<tensor.Str()<<"||"<<ket.Str()<<std::endl<<sector<<std::endl;
    }
  }
} //namespace

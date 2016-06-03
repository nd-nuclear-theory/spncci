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


  void MoshinskyTransformUnitTensor(
        const u3shell::RelativeUnitTensorLabelsU3ST& tensor, 
        double expansion_coef, 
        u3shell::TwoBodySpaceU3ST& space,
        u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_expansion
      )
  {
          //extract operator labels 
    int etap=tensor.bra().eta();
    int eta=tensor.ket().eta();
    HalfInt S=tensor.ket().S();
    HalfInt T=tensor.ket().T();
    HalfInt Sp=tensor.bra().S();
    HalfInt Tp=tensor.bra().T();
    u3::SU3 x0(tensor.x0());
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

      const u3shell::TwoBodySubspaceU3ST& bra_subspace(sector_labels.bra_subspace());
      const u3shell::TwoBodySubspaceU3ST& ket_subspace(sector_labels.ket_subspace());

      if ((bra_subspace.S()!=Sp) || (bra_subspace.T()!=Tp) || (ket_subspace.S()!=S) || (ket_subspace.T()!=T))
        continue;

      u3::U3 omegap(bra_subspace.omega());
      //u3::U3 omegap(sector_labels.bra_subspace().omega());
      int Np=int(omegap.N());
      u3::SU3 xp(omegap.SU3());

      u3::U3 omega(ket_subspace.omega());
      int N=int(omega.N()); //for Two-body state, N must be int. A/2=2/2
      u3::SU3 x(omega.SU3());

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
      // dimension of the subspaces 
      int dimb=bra_subspace.size();
      int dimk=ket_subspace.size();
      // set up two-body sector for given operator,(omegap,Sp,Tp),(omega,S,T), rho0
      Eigen::MatrixXd sector(dimk,dimb); 

      Eigen::MatrixXd bra_moshinky_12(dimb,1);
      Eigen::MatrixXd ket_moshinky_12(1,dimk);
      // Changed particles 
      Eigen::MatrixXd bra_moshinky_21(dimb,1);
      Eigen::MatrixXd ket_moshinky_21(1,dimk);

      std::vector<u3shell::TwoBodyStateLabelsU3ST>bra_labels;
      std::vector<u3shell::TwoBodyStateLabelsU3ST>ket_labels;

      //ucoefficient from factoring out center of mass 
      //the expansion_coef corresponds to the reduce matrix element of the operator in the relative space
      double relative_coef=(u3::U(x0,u3::SU3(eta,0),xp,u3::SU3(eta_cm,0),u3::SU3(etap,0),1,1,x,1,rho0)
                    *expansion_coef);

      // iterate over bra subspace
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int bra_state_index=0; bra_state_index<dimb; ++bra_state_index)
        {
          const u3shell::TwoBodyStateU3ST bra_state(bra_subspace,bra_state_index);
          int eta1p=bra_state.N1();
          int eta2p=bra_state.N2();

          bra_labels.push_back(
            u3shell::TwoBodyStateLabelsU3ST(eta1p, eta2p, xp, bra_subspace.S(), bra_subspace.T())
            );
          
          // Phase arising from exchanging particle 1 and particle 2
          bool phase_even=(eta1p+eta2p+u3::ConjugationGrade(xp)+int(Sp+Tp))%2==0;
          int exchange_phase=phase_even?1:-1;
          
          // overall factor for the bra
          // includes the relative coeffcients of convenience 
          // the factor of 1/sqrt(2) comes from the anti-symmeterization of the states
          // the factor of 1/sqrt(1+delta) comes from the normalization for particles in the same shell
          double coef=relative_coef/std::sqrt(2.*(1+KroneckerDelta(eta1p,eta2p)));

          bra_moshinky_12(bra_state_index,0)=coef*MoshinskyCoefficient(eta1p, eta2p, etap, eta_cm, xp);
          bra_moshinky_21(bra_state_index,0)=exchange_phase*coef*MoshinskyCoefficient(eta2p, eta1p, etap, eta_cm, xp);
        } 

      // iterate over ket subspace
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int ket_state_index=0; ket_state_index<dimk; ++ket_state_index)
        {
          const u3shell::TwoBodyStateU3ST ket_state(ket_subspace,ket_state_index);
          int eta1=ket_state.N1();
          int eta2=ket_state.N2();

          // Phase arising from exchanging particle 1 and particle 2
          bool phase_even=(eta1+eta2+u3::ConjugationGrade(x)+int(S+T))%2==0;
          int exchange_phase=phase_even?1:-1;

          // overall factor for the ket
          // the factor of 1/sqrt(2) comes from the anti-symmeterization of the states
          // the factor of 1/sqrt(1+delta) comes from the normalization for particles in the same shell
          double coef=1./std::sqrt(2.*(1+KroneckerDelta(eta1,eta2)));

          ket_labels.push_back(
            u3shell::TwoBodyStateLabelsU3ST(eta1, eta2, omega.SU3(), ket_subspace.S(), ket_subspace.T())
            );

          ket_moshinky_12(0,ket_state_index)=coef*MoshinskyCoefficient(eta1, eta2, eta, eta_cm, x);
          ket_moshinky_21(0,ket_state_index)=exchange_phase*coef*MoshinskyCoefficient(eta2, eta1, eta, eta_cm, x);
        }  
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // summing sector terms over Ncm and antisymmeterizing
      sector=(
        bra_moshinky_12*ket_moshinky_12
        -bra_moshinky_21*ket_moshinky_12
        -bra_moshinky_12*ket_moshinky_21
        +bra_moshinky_21*ket_moshinky_21
        );
      // Accumluating the moshinsky transformed coeffcients 
      for (int i=0; i<dimb; i++)
        {
        for (int j=0; j<dimk; j++)
          {
            double coefficient=sector(i,j);
            //std::cout<<coefficient<<"  "<<(coefficient>.00005)<<std::endl;
            if (fabs(coefficient)>1.0e-10)
            {
              TwoBodyUnitTensorLabelsU3ST tboperator(tensor,rho0,bra_labels[i],ket_labels[j]);
              // std::cout<<"twobody operator  "<<tboperator.Str()<<std::endl;
              two_body_expansion[tboperator]+=coefficient;                
            }
          }
        }
    }
  }

  u3shell::TwoBodyUnitTensorCoefficientsU3ST 
  TransformRelativeTensorToTwobodyTensor(const RelativeUnitTensorCoefficientsU3ST& relative_unit_tensor_exansion, u3shell::TwoBodySpaceU3ST& space)
  {    
    //container for acculumating two-body coefficients
    u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion;

    for (auto rel_key_value : relative_unit_tensor_exansion)
    {
      u3shell::RelativeUnitTensorLabelsU3ST tensor(rel_key_value.first);
      double expansion_coef=rel_key_value.second;

      MoshinskyTransformUnitTensor(tensor, expansion_coef, space, two_body_expansion);
    }
    return two_body_expansion;
  }
} //namespace

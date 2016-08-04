/****************************************************************
  moshinsky.cpp
                                
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include <cmath>
#include "cppformat/format.h"


#include "sp3rlib/u3coef.h"
#include "u3shell/moshinsky.h"
//#include <eigen3/Eigen/Eigenvalues>  

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



  Eigen::MatrixXd MoshinskyTransform(
        const u3::SU3& x0, 
        int etap,
        int eta,
        const u3shell::TwoBodySubspaceU3ST& bra_subspace, 
        const u3shell::TwoBodySubspaceU3ST& ket_subspace, 
        int rho0
      )
  {
    HalfInt Sp=bra_subspace.S();
    HalfInt Tp=bra_subspace.T();
    u3::U3 omegap(bra_subspace.omega());
    //u3::U3 omegap(sector_labels.bra_subspace().omega());
    int Np=int(omegap.N());
    u3::SU3 xp(omegap.SU3());

    HalfInt S=ket_subspace.S();
    HalfInt T=ket_subspace.T();
    u3::U3 omega(ket_subspace.omega());
    int N=int(omega.N()); //for Two-body state, N must be int. A/2=2/2
    u3::SU3 x(omega.SU3());


    // dimension of the subspaces 
    int dimb=bra_subspace.size();
    int dimk=ket_subspace.size();
    Eigen::MatrixXd sector(dimb,dimk);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // eta_cm is fixed by unit tensor plus constraints
    //    (eta,0)x(eta_cm)->omega and eta+eta_cm=N
    //    (etap,0)x(eta_cm)->omegap and etap+eta_cm=Np
    //    and eta_cm>0
    int eta_cm=Np-etap;
    /// Checking valid value of eta_cm
    if (
      eta_cm==(N-eta)
      &&(eta_cm>=0)
      &&(u3::OuterMultiplicity(u3::SU3(eta,0),u3::SU3(eta_cm,0),x)!=0)
      &&(u3::OuterMultiplicity(u3::SU3(etap,0),u3::SU3(eta_cm,0),xp)!=0)
      )
    {
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // // set up two-body sector for given operator,(omegap,Sp,Tp),(omega,S,T), rho0
      // Eigen::MatrixXd sector(dimk,dimb); 
      Eigen::MatrixXd bra_moshinky_12(dimb,1);
      Eigen::MatrixXd ket_moshinky_12(1,dimk);

      //ucoefficient from factoring out center of mass 
      double relative_coef=u3::U(x0,u3::SU3(eta,0),xp,u3::SU3(eta_cm,0),u3::SU3(etap,0),1,1,x,1,rho0);
                            // *parity(eta+x.lambda()+x.mu()+etap+xp.lambda()+xp.mu());
    
      // iterate over bra subspace
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int bra_state_index=0; bra_state_index<dimb; ++bra_state_index)
        {
          const u3shell::TwoBodyStateU3ST bra_state(bra_subspace,bra_state_index);
          int eta1p=bra_state.N1();
          int eta2p=bra_state.N2();
          
          // Phase arising from exchanging particle 1 and particle 2
          bool phase_even=(eta1p+eta2p+u3::ConjugationGrade(xp)+int(Sp+Tp))%2==0;
          int exchange_phase=phase_even?1:-1;
          
          // overall factor for the bra
          // the factor of 1/sqrt(1+delta) comes from the normalization for particles in the same shell
          double coef=1./std::sqrt(1.+KroneckerDelta(eta1p,eta2p));
          // std::cout<<fmt::format("{} {} {} {} {}  {}", eta1p, eta2p, etap, eta_cm, xp.Str(), MoshinskyCoefficient(eta1p, eta2p, etap, eta_cm, xp))<<std::endl;
          bra_moshinky_12(bra_state_index,0)=coef*MoshinskyCoefficient(eta_cm, etap, eta1p, eta2p, xp);
          // bra_moshinky_12(bra_state_index,0)=coef*MoshinskyCoefficient( eta1p, eta2p, etap,eta_cm, xp);

        } 

      // iterate over ket subspace
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int ket_state_index=0; ket_state_index<dimk; ++ket_state_index)
        {
          const u3shell::TwoBodyStateU3ST ket_state(ket_subspace,ket_state_index);
          int eta1=ket_state.N1();
          int eta2=ket_state.N2();

          // overall factor for the ket
          // the factor of 1/sqrt(1+delta) comes from the normalization for particles in the same shell
          double coef=1./std::sqrt(1.+KroneckerDelta(eta1,eta2));
          // ket_moshinky_12(0,ket_state_index)=coef*MoshinskyCoefficient( eta1, eta2, eta, eta_cm,x);
          ket_moshinky_12(0,ket_state_index)=coef*MoshinskyCoefficient(eta_cm,eta,eta1, eta2,  x);

        }  
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // summing sector terms over Ncm and antisymmeterizing
      sector=(
        2*relative_coef*(bra_moshinky_12*ket_moshinky_12)
        );
    }
    return sector;
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

    OperatorLabelsU3ST operator_labels(etap-eta,x0,tensor.S0(), tensor.T0(), tensor.g0());

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sectors are generated subject to the constraints that 
    // omega.N()+N0=omegap.N()
    // omega x x0 -> omegap
    // S x S0 ->Sp
    // T x T0 ->Tp
    u3shell::TwoBodySectorsU3ST sector_labels_list(space,tensor);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Iterate over two-body sectors
    for (int sector_index=0; sector_index<sector_labels_list.size(); sector_index++)
    {
      // for a given sector, extract bra and ket subspace information
      auto sector_labels=sector_labels_list.GetSector(sector_index);

      const u3shell::TwoBodySubspaceU3ST& bra_subspace(sector_labels.bra_subspace());
      const u3shell::TwoBodySubspaceU3ST& ket_subspace(sector_labels.ket_subspace());
      int rho0=sector_labels.multiplicity_index();

      if ((bra_subspace.S()!=Sp) || (bra_subspace.T()!=Tp) || (ket_subspace.S()!=S) || (ket_subspace.T()!=T))
        continue;

      // dimension of the subspaces 
      int dimb=bra_subspace.size();
      int dimk=ket_subspace.size();
      // set up two-body sector for given operator,(omegap,Sp,Tp),(omega,S,T), rho0
      Eigen::MatrixXd sector(dimk,dimb);

      u3::SU3 xp(bra_subspace.omega().SU3());
      u3::SU3 x(ket_subspace.omega().SU3());
      HalfInt Np(bra_subspace.omega().N());
      HalfInt N(ket_subspace.omega().N());

      int eta_cm=int(Np-etap);
      /// Checking valid value of eta_cm
      if (
        eta_cm!=(N-eta)
        ||(eta_cm<0)
        ||(u3::OuterMultiplicity(u3::SU3(eta,0),u3::SU3(eta_cm,0),x)==0)
        ||(u3::OuterMultiplicity(u3::SU3(etap,0),u3::SU3(eta_cm,0),xp)==0)
        )
        continue;

      sector=MoshinskyTransform(x0, etap, eta, bra_subspace, ket_subspace, rho0);

      // Accumluating the moshinsky transformed coeffcients 
      for (int i=0; i<dimb; ++i)
      {
        const u3shell::TwoBodyStateU3ST bra_state(bra_subspace,i);
        int eta1p=bra_state.N1();
        int eta2p=bra_state.N2();
        u3shell::TwoBodyStateLabelsU3ST bra(eta1p, eta2p, xp, Sp, Tp);
        for (int j=0; j<dimk; ++j)
          {
            const u3shell::TwoBodyStateU3ST ket_state(ket_subspace,j);
            int eta1=ket_state.N1();
            int eta2=ket_state.N2();
            u3shell::TwoBodyStateLabelsU3ST ket(eta1, eta2, x, S, T);

            double coefficient=expansion_coef*sector(i,j);      
            TwoBodyUnitTensorLabelsU3ST tboperator(operator_labels,rho0,bra,ket);

            two_body_expansion[tboperator]+=coefficient;                
            // if((eta1p==0)&&(eta2p==4)&&(eta1==1)&&(eta2==1))
            //   std::cout<<operator_labels.Str()<<"  "<<bra.Str()<<"  "<<ket.Str()<<"   "<<two_body_expansion[tboperator]<<std::endl;
          }

        }
      // remove unit tensors with coefficient zero
      std::vector<TwoBodyUnitTensorLabelsU3ST> delete_list;
      for(auto key_value : two_body_expansion)
      // for(auto key_value=two_body_expansion.begin(); key_value!=end; ++key_value)
      {
        if(fabs(key_value.second)<10e-10)
          delete_list.push_back(key_value.first);
      }
      // for(int i=0; i<delete_list.size(); ++i)
      //   {
      //     auto key=delete_list[i];
      //     two_body_expansion.erase(key);
      //   }
    }
  }

  u3shell::TwoBodyUnitTensorCoefficientsU3ST 
  TransformRelativeTensorToTwobodyTensor(const RelativeUnitTensorCoefficientsU3ST& relative_unit_tensor_expansion, u3shell::TwoBodySpaceU3ST& space)
  {    
    //container for acculumating two-body coefficients
    u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion;

    for (auto rel_key_value : relative_unit_tensor_expansion)
    {
      u3shell::RelativeUnitTensorLabelsU3ST tensor(rel_key_value.first);
      double expansion_coef=rel_key_value.second;

      MoshinskyTransformUnitTensor(tensor, expansion_coef, space, two_body_expansion);
    }
    return two_body_expansion;
  }
} //namespace

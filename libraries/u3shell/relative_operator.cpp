/****************************************************************
  relative_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/relative_operator.h"

#include "cppformat/format.h"
#include "am/am.h"

namespace u3shell {

  // void GenerateRelativeUnitTensorLabelsU3ST(
  //       int Nmax, 
  //       std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels
  //       )
  // {   
  //   #ifdef VERBOSE
  //   std::cout<<"Entering GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
  //   #endif
  //   for(int N0=0; N0<=Nmax; N0+=2)
  //     {
  //       std::vector<RelativeUnitTensorLabelsU3ST> sym_vec;
  //       for(int Sp=0; Sp<=1; Sp++)
  //         for(int Tp=0; Tp<=1; Tp++)
  //           for(int S=0; S<=1; S++)
  //             for (int T=0; T<=1; T++)
  //               for (int S0=abs(S-Sp); S0<=(S+Sp); S0++)
  //                 for (int T0=abs(T-Tp); T0<=(T+Tp); T0++)
  //                   for(int etap=0; etap<=N0+Nmax; etap++)
  //                     {
  //                       //antisymmeterization constraint on ket 
  //                       if ( (etap+Sp+Tp)%2!=1 )
  //                         continue;
                        
  //                       int eta=etap-N0;
  //                       //antisymmeterization constraint on bra 
  //                       if ( (eta+S+T)%2!=1)
  //                         continue;

  //                       u3shell::RelativeStateLabelsU3ST ket(eta,S,T);
  //                       u3shell::RelativeStateLabelsU3ST bra(etap,Sp,Tp);

  //                       MultiplicityTagged<u3::SU3>::vector omega0_set
  //                         =u3::KroneckerProduct(u3::SU3(etap,0),u3::SU3(0,eta));

  //                       for(int w=0; w<omega0_set.size(); w++)
  //                         {
  //                           u3::SU3 x0(omega0_set[w].irrep);
  //                           sym_vec.push_back(u3shell::RelativeUnitTensorLabelsU3ST(x0,S0,T0,bra,ket));
  //                           //std::cout<<"unit tensors  "<<spncci::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T).Str()<<std::endl;
  //                         }
  //                     }       
  //       relative_unit_tensor_labels[N0]=sym_vec;
  //     }
  // #ifdef VERBOSE
  // std::cout<<"Exiting GenerateRelativeUnitTensorLabelsU3ST"<<std::endl;
  // #endif
  // } //end function


 
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

/****************************************************************
  generate_unit_tensor_statistics.cpp

  For checking number of unit tensor rmes for spacial, spin and spacial+spin

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  2/19/21 (aem): Created.
****************************************************************/
#include <iostream>
#include <fstream>
#include "fmt/format.h"
#include "u3shell/relative_operator.h"
#include "u3shell/tensor_labels.h"
#include "lgi/lgi.h"
namespace u3shell{

void GenerateRelativeUnitTensorLabels(
    int Nmax, int N1v,
    std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
    )
  {   
    
    // Assuming parity conserving, N0 must be even
    int Nrel_max=Nmax+2*N1v;
    for(int N=0; N<=Nrel_max; ++N)
      for(int Np=0; Np<=Nrel_max; ++Np)
        {
          int N0=Np-N;
          if (N0%2!=0)
            continue;
            
          // Get allowed x0 values
          MultiplicityTagged<u3::SU3>::vector x0_set
            =u3::KroneckerProduct(u3::SU3(Np,0),u3::SU3(0,N));

           for(int Sp=0; Sp<=1; Sp++)
              for(int Tp=0; Tp<=1; Tp++)
                for(int S=0; S<=1; S++)
                  for (int T=0; T<=1; T++)
                    for (int S0=abs(S-Sp); S0<=(S+Sp); S0++)
                      {
                        //antisymmeterization constraint on ket 
                        if ( (Np+Sp+Tp)%2!=1 )
                          continue;  
                        //antisymmeterization constraint on bra 
                        if ( (N+S+T)%2!=1)
                          continue;

                        int T0_min=abs(Tp-T);
                        int T0_max=Tp+T;
                        
                        for(int T0=T0_min; T0<=T0_max; ++T0)
                          {
                            if(not am::AllowedTriangle(T,Tp,T0))
                              continue;

                            u3shell::RelativeStateLabelsU3ST ket(N,S,T);
                            u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
                            
                            for(int w=0; w<x0_set.size(); w++)
                              {
                                u3::SU3 x0(x0_set[w].irrep);
                                relative_unit_tensor_labels.emplace_back(x0,S0,T0,bra,ket);
                              }
                          }
                      }
        }
  }

void GenerateSpacialRelativeUnitTensorLabels(
    int Nmax, int N1v,
    std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
    )
  {   
    
    // Assuming parity conserving, N0 must be even
    int Nrel_max=Nmax+2*N1v;
    for(int N=0; N<=Nrel_max; ++N)
      for(int Np=0; Np<=Nrel_max; ++Np)
        {
          int N0=Np-N;
          if (N0%2!=0)
            continue;
            
          // Get allowed x0 values
          MultiplicityTagged<u3::SU3>::vector x0_set
            =u3::KroneckerProduct(u3::SU3(Np,0),u3::SU3(0,N));
          
          int Sp=0,S=0,T=0,Tp=0,T0=0,S0=0;


          u3shell::RelativeStateLabelsU3ST ket(N,S,T);
          u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
                            
          for(int w=0; w<x0_set.size(); w++)
            {
              u3::SU3 x0(x0_set[w].irrep);
              relative_unit_tensor_labels.emplace_back(x0,S0,T0,bra,ket);
            }
        }
  }


void GenerateSpinRelativeUnitTensorLabels(
    int Nmax, int N1v,
    std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
    )
  {   
    
    // Assuming parity conserving, N0 must be even
    int N=0, Np=0, N0=0;
    u3::SU3 x0(0,0);

   for(int Sp=0; Sp<=1; Sp++)
      for(int Tp=0; Tp<=1; Tp++)
        for(int S=0; S<=1; S++)
          for (int T=0; T<=1; T++)
            for (int S0=abs(S-Sp); S0<=(S+Sp); S0++)
              {
                //parity conserving constraint
                int etap=(Sp+Tp)%2;
                int eta=(S+T)%2;
                
                if ((etap-eta)%2!=0)
                  continue;
                
                int T0_min=abs(Tp-T);
                int T0_max=Tp+T;
                
                for(int T0=T0_min; T0<=T0_max; ++T0)
                  {
                    if(not am::AllowedTriangle(T,Tp,T0))
                      continue;

                    u3shell::RelativeStateLabelsU3ST ket(N,S,T);
                    u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
                    relative_unit_tensor_labels.emplace_back(x0,S0,T0,bra,ket);
                  }
              }
 
  }


typedef std::tuple<HalfInt,HalfInt,HalfInt>SpinTuple;

}//namespace 


int main(int argc, char **argv)
  {
    if(argc<7)
      {
        std::cout<<"Syntax: Z N Nmax N1v lgi_filename"<<std::endl;
        std::exit(EXIT_FAILURE);
      }
  
    lgi::NuclideType nuclide;
    nuclide[0] = std::stoi(argv[1]);
    nuclide[1] = std::stoi(argv[2]);
    int Nsmax = std::stoi(argv[3]);
    int Nmax = std::stoi(argv[4]);
    int N1v = std::stoi(argv[5]); 
    std::string lgi_filename = argv[6]; 
  
    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
    u3shell::GenerateRelativeUnitTensorLabels(Nmax, N1v,relative_unit_tensor_labels);
    std::cout<<"Num unit tensors: "<<relative_unit_tensor_labels.size()<<std::endl; 
    // for(auto& unit_tensor : relative_unit_tensor_labels)
    //   std::cout<<unit_tensor.Str()<<std::endl;


    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> spacial_relative_unit_tensor_labels;
    u3shell::GenerateSpacialRelativeUnitTensorLabels(Nmax, N1v,spacial_relative_unit_tensor_labels);
    std::cout<<"Num spacial unit tensors: "<<spacial_relative_unit_tensor_labels.size()<<std::endl; 
    // for(auto& unit_tensor : spacial_relative_unit_tensor_labels)
    //   std::cout<<unit_tensor.Str()<<std::endl;



    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> spin_relative_unit_tensor_labels;
    u3shell::GenerateSpinRelativeUnitTensorLabels(Nmax, N1v,spin_relative_unit_tensor_labels);
    std::cout<<"Num spin unit tensors: "<<spin_relative_unit_tensor_labels.size()<<std::endl; 
    // for(auto& unit_tensor : spin_relative_unit_tensor_labels)
    //   std::cout<<unit_tensor.Str()<<std::endl;


    //Read in list of RMEs and figure out how many 
    
    bool intrinsic=true;
    HalfInt N0=lgi::Nsigma0ForNuclide(nuclide, intrinsic);
    
    lgi::MultiplicityTaggedLGIVector lgi_vector;
    lgi::ReadLGISet(lgi_filename, N0,lgi_vector);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Iterate over lgi pairs and get # of spin, spacial and combined tiles and # RME (numbers)
    long int num_spin_tiles=0, num_spin_rmes=0;
    long int num_spacial_tiles=0, num_spacial_rmes=0;
    long int num_tiles=0, num_rmes=0;
    for (auto& tagged_lgi1 : lgi_vector)
      {

        int Nex1;
        u3::U3 sigma1;
        HalfInt Sp1,Sn1,S1;
        std::tie(Nex1,sigma1,Sp1,Sn1,S1)=tagged_lgi1.irrep.Key();

        if (Nex1>Nsmax)
          continue;

        int gamma1=tagged_lgi1.tag;
        int Nn1_max=int(Nmax+N0-sigma1.N());
        sp3r::Sp3RSpace irrep1(sigma1,  Nn1_max);
        for (auto& tagged_lgi2 : lgi_vector)
          {

            if(tagged_lgi1.irrep<tagged_lgi2.irrep)
              continue;

            int Nex2;
            u3::U3 sigma2;
            HalfInt Sp2,Sn2,S2;
            std::tie(Nex2,sigma2,Sp2,Sn2,S2)=tagged_lgi2.irrep.Key();
            int gamma2=tagged_lgi2.tag;
            int Nn2_max=int(Nmax+N0-sigma2.N());
            sp3r::Sp3RSpace irrep2(sigma2,  Nn2_max);

            if (Nex2>Nsmax)
              continue;

            //////////////////////////////////////////////////////////////////////////////////
            for (auto const& tensor : spin_relative_unit_tensor_labels)
              {
                HalfInt S0=tensor.S0();
                if(am::AllowedTriangle(S2,S1,S0) && (abs(Sp1-Sp2)<=2) && (abs(Sn1-Sn2)<=2) && (sigma1==sigma2))
                  {
                    num_spin_rmes += gamma1*gamma2;
                    num_spin_tiles += 1;
                  }     
              }
            // std::cout<<"spin done"<<std::endl;
            //////////////////////////////////////////////////////////////////////////////////


            for(int subspace1_index=0; subspace1_index<irrep1.size(); ++subspace1_index)
              for(int subspace2_index=0; subspace2_index<irrep2.size(); ++subspace2_index)
                {
                  const auto& subspace1=irrep1.GetSubspace(subspace1_index);
                  const auto& subspace2=irrep2.GetSubspace(subspace2_index);
                  
                  u3::U3 omega1=subspace1.U3();
                  int upsilon1=subspace1.upsilon_max();
                  
                  u3::U3 omega2=subspace2.U3();
                  int upsilon2=subspace2.upsilon_max();

                  //////////////////////////////////////////////////////////////////////////////////
                  for (auto const& tensor : spacial_relative_unit_tensor_labels)
                    {
                      u3::SU3 x0=tensor.x0();
                      bool allowed=(u3::OuterMultiplicity(omega1.SU3(),x0,omega2.SU3())>0);
                      allowed&=(omega1.N()+tensor.N0())==omega2.N();
                      allowed&=(S1==S2);
                      allowed&=(Sp1==Sp2);
                      allowed&=(Sn1==Sn2);
                      if(allowed)
                        {
                          num_spacial_rmes += gamma1*gamma2*upsilon1*upsilon2;
                          num_spacial_tiles += 1;
                        }

                    }

                  //////////////////////////////////////////////////////////////////////////////////
                  
                  for (auto const& tensor : relative_unit_tensor_labels)
                    {
                      u3::SU3 x0;
                      HalfInt S0;
                      std::tie(x0,S0,std::ignore,std::ignore,std::ignore,std::ignore,std::ignore,std::ignore,std::ignore)=tensor.FlatKey();
                      if((u3::OuterMultiplicity(omega1.SU3(),x0,omega2.SU3())>0)&&((omega1.N()+tensor.N0())==omega2.N()))
                        if(am::AllowedTriangle(S2,S1,S0) && (abs(Sp1-Sp2)<=2) && (abs(Sn1-Sn2)<=2))
                          {
                            num_rmes += gamma1*gamma2*upsilon1*upsilon2;
                            num_tiles += 1;
                          }

                    }
                  //////////////////////////////////////////////////////////////////////////////////

                }
          // std::cout<<"spacial done"<<std::endl;


          }
      }

    std::cout<<fmt::format("Spin     tiles: {:15}    rmes: {:15d}",num_spin_tiles,num_spin_rmes)<<std::endl;
    std::cout<<fmt::format("Spacial  tiles: {:15}    rmes: {:15d}",num_spacial_tiles,num_spacial_rmes)<<std::endl;
    std::cout<<fmt::format("Total    tiles: {:15}    rmes: {:15d}",num_tiles,num_rmes)<<std::endl;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




  }




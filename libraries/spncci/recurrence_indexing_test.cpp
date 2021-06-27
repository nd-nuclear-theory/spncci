/****************************************************************
    recurrence_indexing_test.cpp

  Anna E. McCoy
  INT

  SPDX-License-Identifier: MIT

 6/25/21 (aem): Created.
****************************************************************/
#include "spncci/recurrence_indexing.h"
#include "am/halfint_fmt.h"
#include "fmt/format.h"
#include "lgi/lgi.h"


int main(int argc, char **argv)
{

        //   spncci::ValenceShellForNuclide({2,1}),

  ////////////////////////////////////////////////////////////////
  // set up SpinSpace
  ////////////////////////////////////////////////////////////////

  // read in LGIs for 6Li
  std::string filename = "../../data/lgi_set/lgi_test.dat";  // test file in data/lgi_set/lgi_test.dat
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  HalfInt Nsigma0=lgi::Nsigma0ForNuclide({3,3});
  lgi::ReadLGISet(filename,Nsigma0,lgi_vector);
  auto spin_space=spncci::spin::Space(lgi_vector);

  if(true)
  {
    // diagnostic -- inspect LGI listing
    std::cout << "LGI set" << std::endl;
    for (int i=0; i<lgi_vector.size(); ++i)
      std::cout << i << " " << lgi_vector[i].Str() << std::endl;
    std::cout << "********************************" << std::endl;

    int num_lgi_spaces=spin_space.size();
    std::cout<<fmt::format("spin::Space dimension: {}",spin_space.dimension())<<std::endl;
    for(int i=0; i<num_lgi_spaces; ++i)
      {
        const spncci::spin::LGISpace& lgi_space=spin_space.GetSubspace(i);
        int num_subspaces=lgi_space.size();
        const u3::U3& sigma=lgi_space.sigma(); 
        std::cout<<fmt::format("Sigma Space : {}    Size: {:4d}   Dimension: {:4d}",sigma.Str(),num_subspaces,lgi_space.dimension())<<std::endl;
        
        for(int j=0; j<num_subspaces; ++j)
          {
            const spncci::spin::SpNCCISpinSubspace& spin_subspace=lgi_space.GetSubspace(j);
            const HalfInt& S=spin_subspace.S();
            int num_spin_states=spin_subspace.size();
            std::cout<<fmt::format("    Spin Subspace : {}       Size: {:4d}   Dimension: {:4d}",S,num_spin_states,spin_subspace.dimension())<<std::endl;
            for (int k=0; k<num_spin_states; ++k)
              {
                spncci::spin::SpinState spin_state(spin_subspace, k);
                const HalfInt& Sp=spin_state.Sp();
                const HalfInt& Sn=spin_state.Sn();
                const int gamma=spin_state.degeneracy();

                int state_index=spin_space.GetSubspaceOffset(i)  // offset to given sigma subspace
                      + spin_space.GetSubspace(i).GetSubspaceOffset(j) 
                      + spin_space.GetSubspace(i).GetSubspace(j).GetStateOffset(k, gamma);

                
                std::cout<<fmt::format("       Spin State : {} {}      gamma: {:4d}",Sp,Sn,gamma)<<std::endl;
                std::cout<<fmt::format("         index: {}",state_index)<<std::endl;

              }            

          }
      }

  }



}
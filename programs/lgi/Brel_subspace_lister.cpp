
/****************************************************************
  Program that generates bra model space for Brel given an input model space
  
  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  03/30/19 (aem) : Created
****************************************************************/


#include <iostream>
#include <fstream>

#include "fmt/format.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "utilities/nuclide.h"


typedef std::tuple<int,int,int,int,u3::SU3> U3SubspaceLabel;
int main(int argc, char **argv)
{
  if (argc<4)
  {
    std::cout<<"Syntax: Z N <inputfile> <outputfile> \n <inputfile> contains list of subspace labels Nex 2Sp 2Sn 2S lambda mu"<<std::endl;
  }
  

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////
  // SU(3) caching
  u3::U3CoefInit();

  std::array<int,2> nuclide;
  nuclide[0]=std::stoi(argv[1]);
  nuclide[1]=std::stoi(argv[2]);
  int A=nuclide[0]+nuclide[1];
  std::string subspace_filename_ket=argv[3];
  std::string subspace_filename_bra=argv[4];
  
  HalfInt N0 = nuclide::Nsigma0ForNuclide(nuclide);
  // Read in labels 
  std::string line;
  std::ifstream stream_in;
  stream_in.open(subspace_filename_ket);
  std::vector<U3SubspaceLabel> subspace_labels_ket;
  while(std::getline(stream_in,line))
    {
      // parse line
      std::istringstream line_stream(line);
      int Nex,twice_Sp,twice_Sn,twice_S, lambda, mu;
      line_stream>>Nex>>twice_Sp>>twice_Sn>>twice_S>>lambda>>mu;
      subspace_labels_ket.push_back(U3SubspaceLabel(Nex,twice_Sp,twice_Sn,twice_S,u3::SU3(lambda,mu)));
    }

  //Iterate through ket labels and generate set of bra labels for Brel and Nrel
    std::set<U3SubspaceLabel> subspace_labels_bra;
  for(auto label: subspace_labels_ket)
    {
      // Add label to bra set 
      subspace_labels_bra.insert(label);

      //Extract labels and generate list of subspaces connected by Brel
      int Nex,twice_Sp,twice_Sn,twice_S;
      u3::SU3 x;
      std::tie(Nex,twice_Sp,twice_Sn,twice_S,x)=label;
      
      if(Nex!=0)
        {
          MultiplicityTagged<u3::SU3>::vector bra_su3_labels=KroneckerProduct(x, u3::SU3(0,2));
          for(auto xp_tagged : bra_su3_labels)
            {
              const u3::SU3& xp=xp_tagged.irrep;
              u3::U3 wp=u3::U3(N0+Nex-2, xp);
              if(wp.Valid() && wp.f3()>A/2)
              {
                U3SubspaceLabel labelp(Nex-2,twice_Sp,twice_Sn,twice_S,xp);
                subspace_labels_bra.insert(labelp);
              }

            }
        }
    }
  stream_in.close();

  //Write bra labels to file
  std::ofstream stream_out;
  stream_out.open(subspace_filename_bra);
  for(auto label : subspace_labels_bra)
    {
      //Extract labels and generate list of subspaces connected by Brel
      int Nex,twice_Sp,twice_Sn,twice_S;
      u3::SU3 x;
      std::tie(Nex,twice_Sp,twice_Sn,twice_S,x)=label;
      
      //Write labels to file
      stream_out<<fmt::format("{:2d}   {:2d} {:2d} {:2d}   {:2d} {:2d}",Nex,twice_Sp,twice_Sn,twice_S,x.lambda(),x.mu())<<std::endl;

    }
  stream_out.close();

}

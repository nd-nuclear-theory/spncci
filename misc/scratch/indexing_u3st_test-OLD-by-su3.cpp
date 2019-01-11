/****************************************************************
  indexing_u3st_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "fmt/format.h"

#include "u3shell/indexing_u3st.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // relativ-cm state tests
  ////////////////////////////////////////////////////////////////


  if (true)
    {
      std::cout << "Relative-cm subspace" << std::endl;

      // subspace construction

      // RelativeSubspaceLSJT(0,0,0,0,0,0);  // should violate assertion due to T
      // RelativeSubspaceLSJT(0,0,0,1,0,7);  // should violate assertion due to Nmax
      u3shell::RelativeCMSubspaceU3ST subspace(0,u3::SU3(2,0),1,0,0,6);  // Ncm x S T g Nmax
      std::cout << fmt::format("size {}",subspace.size()) << std::endl;
      

      // index-based looping
      for (int k=0; k<subspace.Dimension(); ++k)
	{
	  u3shell::RelativeCMStateU3ST state(subspace,k);
          std::cout << fmt::format("index {} N {}",state.index(),state.N()) << std::endl;
	};

      std::cout << "State lookup" << std::endl;
      for (int Nr=2; Nr<=2; Nr+=2)
	{
	  u3shell::RelativeCMStateU3ST state(subspace,u3shell::RelativeCMSubspaceU3ST::StateLabelsType(Nr));
	  std::cout << fmt::format("Nr {} -> index {} N {}",Nr,state.index(),state.N()) << std::endl;
	}


      for (int lm=0; lm<=10; ++lm)
        for (int mu=0; mu<=10; ++mu)
          for (int S=0; S<=1; ++S)
            for (int Ncm=0; Ncm<=8; ++Ncm)
            {
              int T = 0;
              u3::SU3 x(lm,mu);
              int g = (Ncm + S + T + 1)%2;
              int Nmax = 10+g;
              u3shell::RelativeCMSubspaceU3ST subspace(Ncm,x,S,T,g,Nmax);  // Ncm x S T g Nmax
              if (subspace.size()>0)
                {
                  std::cout 
                    << fmt::format(
                                   "Ncm {} x {} S {} T {} g {} size {}",
                                   subspace.Ncm(),subspace.SU3().Str(),subspace.S(),subspace.T(),subspace.g(),subspace.size()
                                   )
                    << std::endl;
                }
            }

    }

  ////////////////////////////////////////////////////////////////
  // two-body state tests
  ////////////////////////////////////////////////////////////////
#if 0
  if (false)
    {
      std::cout << "Two-body subspace" << std::endl;

      // subspace construction

      TwoBodySubspaceLSJT subspace(0,0,0,1,0,6);

      // index-based looping
      for (int k = 0; k < subspace.Dimension(); ++k)
	{
	  TwoBodyStateLSJT state(subspace,k);
	  std::cout << state.index() 
		    << " " << state.N() 
		    << " " << state.N1() << " " << state.l1() << " " << state.N2() << " " << state.l2() 
		    << std::endl; // << " " << orbital.GetIndex1()
	};


      // spaces and sectors

      // first set up relative spaces

      std::cout << "Two-body space" << std::endl;
      int Nmax = 3;
      TwoBodySpaceLSJT space(Nmax);
      space.Print(std::cout);

      //      for (int s=0; s<space.size(); ++s)
      //	{
      //	  const TwoBodySubspaceLSJT& subspace = space.GetSubspace(s);
      //	  std::cout 
      //	    << std::setw(3) << s 
      //	    << std::setw(3) << subspace.L() 
      //	    << std::setw(3) << subspace.S() 
      //	    << std::setw(3) << subspace.J() 
      //	    << std::setw(3) << subspace.T() 
      //	    << std::setw(3) << subspace.g()
      //	    << std::setw(3) << subspace.Nmax()
      //	    << std::setw(3) << subspace.Dimension()
      //	    << std::endl;
      //	}

      // then set up allowed sectors
      std::cout << "Two-body interaction sectors" << std::endl;
      int J0 = 0;  // also try J0=2 for quadrupole operator
      int g0 = 0;
      TwoBodySectorsLSJT sectors(space,J0,g0);

      for (int sector_index=0; sector_index < sectors.size(); ++sector_index)
	{
	  int s2 = sectors.GetSector(sector_index).index2();
	  const TwoBodySubspaceLSJT& subspace2 = space.GetSubspace(s2);
	  int s1 = sectors.GetSector(sector_index).index1();
	  const TwoBodySubspaceLSJT& subspace1 = space.GetSubspace(s1);

	  std::cout 
	    << " sector " 
	    << std::setw(3) << sector_index 
	    << "     "
	    << " index "
	    << std::setw(3) << s2
	    << " LSJTg "
	    << std::setw(3) << subspace2.L() 
	    << std::setw(3) << subspace2.S() 
	    << std::setw(3) << subspace2.J() 
	    << std::setw(3) << subspace2.T() 
	    << std::setw(3) << subspace2.g()
	    // << std::setw(3) << subspace2.Nmax()
	    << " dim "
	    << std::setw(3) << subspace2.Dimension()
	    << "     "
	    << " index "
	    << std::setw(3) << s1
	    << " LSJTg "
	    << std::setw(3) << subspace1.L() 
	    << std::setw(3) << subspace1.S() 
	    << std::setw(3) << subspace1.J() 
	    << std::setw(3) << subspace1.T() 
	    << std::setw(3) << subspace1.g()
	    // << std::setw(3) << subspace1.Nmax()
	    << " dim "
	    << std::setw(3) << subspace1.Dimension()
	    << std::endl;
	}
    }
#endif
  


  // termination
  return 0;
}

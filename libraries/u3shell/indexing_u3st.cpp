/****************************************************************
  indexing_u3st.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include "u3shell/indexing_u3st.h"

namespace u3shell {


  ////////////////////////////////////////////////////////////////
  // relative-cm relative LSJT scheme implementation 
  ////////////////////////////////////////////////////////////////

  RelativeCMSubspaceU3ST::RelativeCMSubspaceU3ST(int Ncm, const u3::SU3& x, int S, int T, int g, int Nmax)
  {
    // set values
    labels_ = SubspaceLabelsType(Ncm,x,S,T,g);
    Nmax_ = Nmax;

    // validate
    assert(ValidLabels()); 

    // set up indexing

    // iterate over total oscillator quanta
    for (int N = g; N <= Nmax; N +=2)
      {
	int Nr = N - Ncm;

        // impose coupling constraint
        if (u3::OuterMultiplicity(u3::SU3(Ncm,0),u3::SU3(Nr,0),x)==0)
          continue;

        // keep surviving states
        PushStateLabels(StateLabelsType(Nr)); 
      }
      

  }

  bool RelativeCMSubspaceU3ST::ValidLabels() const
  {
    bool valid = true;

    // impose antisymmetry
    valid &= ((Ncm()+S()+T()+g())%2==1);

    // truncation
    valid &= ((Nmax()%2)==g());

    return valid;
  }

  ////////////////////////////////////////////////////////////////
  // two-body LSJT scheme implementation 
  ////////////////////////////////////////////////////////////////

#if 0
  TwoBodySubspaceLSJT::TwoBodySubspaceLSJT(int L, int S, int J, int T, int g, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(L,S,J,T,g);
    Nmax_ = Nmax;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over total oscillator quanta
    for (int N = g; N <= Nmax; N +=2)
      // iterate over oscillator (Nl) orbitals for particle 1
      for (int N1 = 0; N1 <= N; ++N1)
	for (int l1 = N1%2; l1 <= N1; l1 +=2) 
	  {
	    // iterate over oscillator (Nl) orbitals for particle 2
	    // subject to given total N
	    int N2 = N - N1;
	    for (int l2 = N2%2; l2 <= N2; l2 +=2) 
	      {
		// impose triangularity
		if (!(AllowedTriangle(l1,l2,L)))
		  continue;

		// impose antisymmetry
		if ((N1==N2)&&(l1==l2)&&(!((L+S+T)%2==1)))
		  continue;

		// keep surviving states
		PushStateLabels(StateLabelsType(N1,l1,N2,l2)); 
	      }
	  }
  }

  bool TwoBodySubspaceLSJT::ValidLabels() const
  {
    bool valid = true;
      
    // triangularity
    valid &= AllowedTriangle(L(),S(),J());

    // truncation
    valid &= ((Nmax()%2)==g());

    return valid;
  }


  TwoBodySpaceLSJT::TwoBodySpaceLSJT(int Nmax)
  {
    // save Nmax
    Nmax_ = Nmax;

    // iterate over L
    for (int L=0; L<=Nmax; ++L)
      {
	// iterate over S
	for (int S=0; S<=1; ++S)
	  {
	    // iterate over J
	    // imposing triangularity (LSJ)
	    for (int J=abs(L-S); J<=L+S; ++J)
	      {

		// iterate over T
		for (int T=0; T<=1; ++T)
		  {

		    // iterate over g
		    for (int g=0; g<=1; ++g)

		      {
			
			// downshift Nmax to match parity of subspace
			// required to pass label validity tests
			int Nmax_subspace = Nmax - (Nmax-g)%2;
		    
			TwoBodySubspaceLSJT subspace(L,S,J,T,g,Nmax_subspace);
			// std::cout 
			//    << std::setw(3) << L 
			// 	  << std::setw(3) << S 
			// 	  << std::setw(3) << T 
			// 	  << std::setw(3) << J 
			// 	  << std::setw(3) << g 
			// 	  << std::setw(3) << Nmax_for_subspace 
			// 	  << std::setw(3) << subspace.Dimension()
			// 	  << std::endl;
			if (subspace.Dimension()!=0)
			  PushSubspace(subspace);
		      }
		  }
	      }
	  }
      }
  }

  void TwoBodySpaceLSJT::Print(std::ostream& os) const
  {

    const int lw = 3;

    for (int s=0; s<size(); ++s)
      {
	const SubspaceType& subspace = GetSubspace(s);
	os
	  << " " << "index"
	  << " " << std::setw(lw) << s 
	  << " " << "LSJTg"
	  << " " << std::setw(lw) << subspace.L() 
	  << " " << std::setw(lw) << subspace.S() 
	  << " " << std::setw(lw) << subspace.J() 
	  << " " << std::setw(lw) << subspace.T() 
	  << " " << std::setw(lw) << subspace.g()
	  << " " << "Nmax"
	  << " " << std::setw(lw) << subspace.Nmax()
	  << " " << "dim"
	  << " " << std::setw(lw) << subspace.Dimension()
	  << " " << std::endl;
      }
  }


  TwoBodySectorsLSJT::TwoBodySectorsLSJT(TwoBodySpaceLSJT& space, int J0, int g0)
  {
    for (int s2=0; s2<space.size(); ++s2)
      for (int s1=0; s1<space.size(); ++s1)
	{
	  // verify triangularity and p 
	  int J2 = space.GetSubspace(s2).J();
	  int J1 = space.GetSubspace(s1).J();
	  int g2 = space.GetSubspace(s2).g();
	  int g1 = space.GetSubspace(s1).g();

	  if ( AllowedTriangle(J2,J0,J1) && ((g2+g0+g1)%2==0))
	    PushSector(Sector(s2,s1));
	}
  }
#endif
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

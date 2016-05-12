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

  TwoBodySubspaceU3ST::TwoBodySubspaceU3ST(u3::SU3 x, int S, int T, int g, int Nmax)
  {

    // set values
    labels_ = SubspaceLabelsType(x,S,T,g);
    Nmax_ = Nmax;

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over total oscillator quanta
    for (int N = g; N <= Nmax; N +=2)
      // iterate over oscillator quanta for particle 1 
      for (int N1 = 0; N1 <= N; ++N1)
    	  {
	       // oscillator quanta for particle 2 subject to given total N
    	    int N2 = N - N1;
    	    
      		// impose coupling constraint
      		if (u3::OuterMultiplicity(u3::SU3(N1,0),u3::SU3(N2,0),x)==0)
      		  continue;

      		// impose antisymmetry constraint
      		if ((N1==N2)&&(x.lambda()+x.mu()+S+T)%2==1))
      		  continue;

      		// keep surviving states
      		PushStateLabels(StateLabelsType(N1,N2)); 
      	}
	}

  bool TwoBodySubspaceU3ST::ValidLabels() const
  {
    bool valid = true; 
    // truncation
    valid &= ((Nmax()%2)==g());

    return valid;
  }


  TwoBodySpaceU3ST::TwoBodySpaceU3ST(int Nmax)
  {
    // save Nmax
    Nmax_ = Nmax;

    // iterate over lambda, max value N1max+N2max=2Nmax
    for (int lambda=0; lambda<=2*Nmax; ++lambda)
      //iterate over mu
      for (int mu=0; mu<=Nmax; ++mu)      	
        // iterate over S
      	for (int S=0; S<=1; ++S)
      		// iterate over T
      		for (int T=0; T<=1; ++T)
    		    // iterate over g
    		    for (int g=0; g<=1; ++g)  
		          {
          			// downshift Nmax to match parity of subspace
          			// required to pass label validity tests
          			int Nmax_subspace = Nmax - (Nmax-g)%2;
		    
			          TwoBodySubspaceU3ST subspace(u3::SU3(lambda,mu),S,T,g,Nmax_subspace);
			
          			if (subspace.Dimension()!=0)
          			  PushSubspace(subspace);

		          }
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

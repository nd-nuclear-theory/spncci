/****************************************************************
  indexing_u3st.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <sstream>

#include "cppformat/format.h"
#include "u3shell/indexing_u3st.h"

namespace u3shell {

  TwoBodySubspaceU3ST::TwoBodySubspaceU3ST(u3::U3 omega, int S, int T, int g)
  {

    // set values
    labels_ = SubspaceLabelsType(omega,S,T,g);

    // validate subspace labels
    assert(ValidLabels()); 

    // set up indexing
    // iterate over oscillator quanta for particle 1 
    u3::SU3 x = omega.SU3();
    for (int N1 = 0; N1 <= N(); ++N1)
      {
        // oscillator quanta for particle 2 subject to given total N
        int N2 = N() - N1;

        // impose coupling constraint
        if (u3::OuterMultiplicity(u3::SU3(N1,0),u3::SU3(N2,0),x)==0)
          continue;

        // impose antisymmetry constraint
        if ((N1==N2)&&!((ConjugationGrade(x)+S+T+1)%2==0))
          continue;

        // keep surviving states
        PushStateLabels(StateLabelsType(N1)); 
      }

  }

  bool TwoBodySubspaceU3ST::ValidLabels() const
  {
    bool valid = true; 
    // truncation
    valid &= (N()+g())%2==0;

    return valid;
  }

  std::string TwoBodySubspaceU3ST::Str() const
  {

    return fmt::format(
                       "[{} {} {} {}]",
                       omega().Str(),S(),T(),g()
                       );
  }

  TwoBodySpaceU3ST::TwoBodySpaceU3ST(int Nmax)
  {
    // save Nmax
    Nmax_ = Nmax;

    // for each N in 0..Nmax
    for (int N=0; N<=Nmax; ++N)
      // for lambda in 0..N
      for (int lambda=0; lambda<=N; ++lambda)
        //for mu in 0..floor(N/2)
        for (int mu=0; mu<=N/2; ++mu)      	
        // for each S in 0..1
          for (int S=0; S<=1; ++S)
            // for each T in 0..1
            for (int T=0; T<=1; ++T)
              
              {
                // g is fixed by g~N
                int g = N%2;
                
                // create U(3) label
                u3::SU3 x(lambda,mu);
                
                // check validity of U(3) before attempting construction
                if (!u3::U3::ValidLabels(N,x))
                  continue;

                // construct subspace
                u3::U3 omega(N,x);
                TwoBodySubspaceU3ST subspace(omega,S,T,g);
			
                // push subspace if nonempty
                if (subspace.Dimension()!=0)
                  PushSubspace(subspace);
              }
  }


  std::string TwoBodySpaceU3ST::Str() const
  {

    std::ostringstream os;
    for (int s=0; s<size(); ++s)
      {
	const SubspaceType& subspace = GetSubspace(s);
        
        std::string subspace_string 
          = fmt::format(
                        "index {:3} {} size {:3}",
                        s,
                        subspace.Str(),
                        subspace.size()
                        );
                        
        os << subspace_string << std::endl;
      }
    return os.str();
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

/****************************************************************
  u3st_scheme.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/

#include "u3shell/u3st_scheme.h"

#include <sstream>

#include "am/am.h"
#include "fmt/format.h"

namespace u3shell {

  RelativeSubspaceU3ST::RelativeSubspaceU3ST(int N, int S, int T, int g)
    : BaseSubspace{SubspaceLabelsType(N,S,T,g)}
  {
    // validate subspace labels
    assert(ValidLabels());
    PushStateLabels(1);
  }

  bool RelativeSubspaceU3ST::ValidLabels() const
  {
    bool valid = true;
    // truncation
    valid &= (N()+g())%2==0;

    return valid;
  }

  std::string RelativeSubspaceU3ST::Str() const
  {

    return fmt::format(
        "[{} {} {} {}]",
        N(),S(),T(),g()
      );
  }

  RelativeSpaceU3ST::RelativeSpaceU3ST(int Nmax)
  {
    // save Nmax
    Nmax_ = Nmax;

    // for each N in 0..Nmax
    for (int N=0; N<=Nmax; ++N)
      // for each S in 0..1
      for (int S=0; S<=1; ++S)
        // for each T in 0..1
        for (int T=0; T<=1; ++T)
          {
            // g is fixed by g~N
            int g = N%2;
            if((N+S+T)%2==1)
              {
                // construct subspace
                RelativeSubspaceU3ST subspace(N,S,T,g);

                // push subspace if nonempty
                if (subspace.size()!=0)
                  PushSubspace(subspace);
              }
          }
  }


  std::string RelativeSpaceU3ST::Str() const
  {

    std::ostringstream os;
    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);

        std::string subspace_string
          = fmt::format(
              "index {:3} {} size {:3}",
              subspace_index,
              subspace.Str(),
              subspace.size()
            );

        os << subspace_string << std::endl;
      }
    return os.str();
  }

  RelativeSectorsU3ST::RelativeSectorsU3ST(RelativeSpaceU3ST& space)
    : BaseSectors(space)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const RelativeSubspaceU3ST& bra_subspace = space.GetSubspace(bra_subspace_index);
          const RelativeSubspaceU3ST& ket_subspace = space.GetSubspace(ket_subspace_index);

          // push sectors (taking unit multiplicity)
          int multiplicity_index = 1;
          PushSector(bra_subspace_index,ket_subspace_index,multiplicity_index);

        }
  }

  RelativeSectorsU3ST::RelativeSectorsU3ST(RelativeSpaceU3ST& space, const OperatorLabelsU3ST& operator_labels)
    : BaseSectors(space)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const RelativeSubspaceU3ST& bra_subspace = space.GetSubspace(bra_subspace_index);
          const RelativeSubspaceU3ST& ket_subspace = space.GetSubspace(ket_subspace_index);


          // verify selection rules
          bool allowed = true;
          // U(1)
          allowed &= (ket_subspace.N() + operator_labels.N0() - bra_subspace.N() == 0);
          // spin & isosopin
          allowed &= am::AllowedTriangle(ket_subspace.S(),operator_labels.S0(),bra_subspace.S());
          allowed &= am::AllowedTriangle(ket_subspace.T(),operator_labels.T0(),bra_subspace.T());
          // parity
          //MARK?
          //allowed &= ((ket_subspace.g() + operator_labels.g0() - bra_subspace.g())%2 == 0);
          // find SU(3) multiplicity and check SU(3) selection
          int multiplicity = 0;
          if (allowed)
            {
              multiplicity = u3::OuterMultiplicity(ket_subspace.SU3(),operator_labels.x0(),bra_subspace.SU3());
              allowed &= (multiplicity > 0);
              // std::cout << fmt::format("{}x{}->{}: {} {}",ket_subspace.omega().SU3().Str(),operator_labels.x0.Str(),bra_subspace.omega().SU3().Str(),multiplicity,allowed)
              //           << std::endl;
            }

          // push sectors (tagged by multiplicity)
          if (allowed)
            for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
              {
                PushSector(bra_subspace_index,ket_subspace_index,multiplicity_index);
              }
        }
  }

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
  RelativeCMSubspaceU3ST::RelativeCMSubspaceU3ST(u3::U3 omega, int S, int T, int g)
    : BaseSubspace{SubspaceLabelsType(omega,S,T,g)}
  {

    // validate subspace labels
    assert(ValidLabels());

    // set up indexing
    // iterate over oscillator quanta for particle 1
    u3::SU3 x = omega.SU3();
    for (int Nr = 0; Nr <= N(); ++Nr)
      {
        // oscillator quanta for particle 2 subject to given total N
        int Ncm = N() - Nr;

        // impose coupling constraint
        if (u3::OuterMultiplicity(u3::SU3(Nr,0),u3::SU3(Ncm,0),x)==0)
          continue;

        // impose antisymmetry constraint
        if ((Nr+S+T)%2==0)
          continue;

        // keep surviving states
        PushStateLabels(StateLabelsType(Nr,Ncm));
      }
  }

  bool RelativeCMSubspaceU3ST::ValidLabels() const
  {
    bool valid = true;
    // truncation
    valid &= (N()+g())%2==0;

    return valid;
  }

  std::string RelativeCMSubspaceU3ST::Str() const
  {

    return fmt::format(
        "[{} {} ({},{})]",
        omega().Str(),S(),T(),g()
      );
  }

  RelativeCMSpaceU3ST::RelativeCMSpaceU3ST(int Nmax)
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
                RelativeCMSubspaceU3ST subspace(omega,S,T,g);

                // push subspace if nonempty
                if (subspace.size()!=0)
                  PushSubspace(subspace);
              }
  }


  std::string RelativeCMSpaceU3ST::Str() const
  {

    std::ostringstream os;
    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
  const SubspaceType& subspace = GetSubspace(subspace_index);

        std::string subspace_string
          = fmt::format(
              "index {:3} {} size {:3}",
              subspace_index,
              subspace.Str(),
              subspace.size()
            );

        os << subspace_string << std::endl;
      }
    return os.str();
  }

  RelativeCMSectorsU3ST::RelativeCMSectorsU3ST(RelativeCMSpaceU3ST& space)
    : BaseSectors(space)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
  {
          // retrieve subspaces
          const RelativeCMSubspaceU3ST& bra_subspace = space.GetSubspace(bra_subspace_index);
          const RelativeCMSubspaceU3ST& ket_subspace = space.GetSubspace(ket_subspace_index);

          // push sectors (taking unit multiplicity)
          int multiplicity_index = 1;
          PushSector(bra_subspace_index,ket_subspace_index,multiplicity_index);

  }
  }

  RelativeCMSectorsU3ST::RelativeCMSectorsU3ST(RelativeCMSpaceU3ST& space, const OperatorLabelsU3ST& operator_labels)
    : BaseSectors(space)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
  {
          // retrieve subspaces
          const RelativeCMSubspaceU3ST& bra_subspace = space.GetSubspace(bra_subspace_index);
          const RelativeCMSubspaceU3ST& ket_subspace = space.GetSubspace(ket_subspace_index);


    // verify selection rules
          bool allowed = true;
          // U(1)
          allowed &= (ket_subspace.N() + operator_labels.N0() - bra_subspace.N() == 0);
          // spin & isosopin
          allowed &= am::AllowedTriangle(ket_subspace.S(),operator_labels.S0(),bra_subspace.S());
          allowed &= am::AllowedTriangle(ket_subspace.T(),operator_labels.T0(),bra_subspace.T());
          // parity
          allowed &= ((ket_subspace.g() + operator_labels.g0() - bra_subspace.g())%2 == 0);
          // find SU(3) multiplicity and check SU(3) selection
          int multiplicity = 0;
          if (allowed)
            {
              multiplicity = u3::OuterMultiplicity(ket_subspace.omega().SU3(),operator_labels.x0(),bra_subspace.omega().SU3());
              allowed &= (multiplicity > 0);
              // std::cout << fmt::format("{}x{}->{}: {} {}",ket_subspace.omega().SU3().Str(),operator_labels.x0.Str(),bra_subspace.omega().SU3().Str(),multiplicity,allowed)
              //           << std::endl;
            }

          // push sectors (tagged by multiplicity)
    if (allowed)
            for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
              {
                PushSector(bra_subspace_index,ket_subspace_index,multiplicity_index);
              }
  }
  }



  TwoBodySubspaceU3ST::TwoBodySubspaceU3ST(u3::U3 omega, int S, int T, int g)
    : BaseSubspace{SubspaceLabelsType(omega,S,T,g)}
  {

    // validate subspace labels
    assert(ValidLabels());

    // set up indexing
    // iterate over oscillator quanta for particle 1
    u3::SU3 x = omega.SU3();
    for (int N1 = 0; N1 <= N(); ++N1)
      {
        // oscillator quanta for particle 2 subject to given total N
        int N2 = N() - N1;

        // impose canonical ordering one single-particle states
        if (!(N1<=N2))
          continue;

        // impose coupling constraint
        if (u3::OuterMultiplicity(u3::SU3(N1,0),u3::SU3(N2,0),x)==0)
          continue;

        // impose antisymmetry constraint
        if ((N1==N2)&&!((ConjugationGrade(x)+S+T+1)%2==0))
          continue;

        // keep surviving states
        PushStateLabels(StateLabelsType(N1,N2));
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
        omega().Str(),S(),T(),g());
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
                if (subspace.size()!=0)
                  PushSubspace(subspace);
              }
  }


  std::string TwoBodySpaceU3ST::Str() const
  {

    std::ostringstream os;
    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
	const SubspaceType& subspace = GetSubspace(subspace_index);

        std::string subspace_string
          = fmt::format(
              "index {:3} {} size {:3}",
              subspace_index,
              subspace.Str(),
              subspace.size()
            );

        os << subspace_string << std::endl;
      }
    return os.str();
  }

  TwoBodySectorsU3ST::TwoBodySectorsU3ST(TwoBodySpaceU3ST& space)
    : BaseSectors(space)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
	{
          // retrieve subspaces
          const TwoBodySubspaceU3ST& bra_subspace = space.GetSubspace(bra_subspace_index);
          const TwoBodySubspaceU3ST& ket_subspace = space.GetSubspace(ket_subspace_index);

          // push sectors (taking unit multiplicity)
          int multiplicity_index = 1;
          PushSector(bra_subspace_index,ket_subspace_index,multiplicity_index);

	}
  }

  TwoBodySectorsU3ST::TwoBodySectorsU3ST(TwoBodySpaceU3ST& space, const OperatorLabelsU3ST& operator_labels)
    : BaseSectors(space)
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
	{
          // retrieve subspaces
          const TwoBodySubspaceU3ST& bra_subspace = space.GetSubspace(bra_subspace_index);
          const TwoBodySubspaceU3ST& ket_subspace = space.GetSubspace(ket_subspace_index);


	  // verify selection rules
          bool allowed = true;
          // U(1)
          allowed &= (ket_subspace.N() + operator_labels.N0() - bra_subspace.N() == 0);
          // spin & isosopin
          allowed &= am::AllowedTriangle(ket_subspace.S(),operator_labels.S0(),bra_subspace.S());
          allowed &= am::AllowedTriangle(ket_subspace.T(),operator_labels.T0(),bra_subspace.T());
          // parity
          allowed &= ((ket_subspace.g() + operator_labels.g0() - bra_subspace.g())%2 == 0);
          // find SU(3) multiplicity and check SU(3) selection
          int multiplicity = 0;
          if (allowed)
            {
              multiplicity = u3::OuterMultiplicity(ket_subspace.omega().SU3(),operator_labels.x0(),bra_subspace.omega().SU3());
              allowed &= (multiplicity > 0);
              // std::cout << fmt::format("{}x{}->{}: {} {}",ket_subspace.omega().SU3().Str(),operator_labels.x0.Str(),bra_subspace.omega().SU3().Str(),multiplicity,allowed)
              //           << std::endl;
            }

          // push sectors (tagged by multiplicity)
	  if (allowed)
            for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
              {
                PushSector(bra_subspace_index,ket_subspace_index,multiplicity_index);
              }
	}
  }


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

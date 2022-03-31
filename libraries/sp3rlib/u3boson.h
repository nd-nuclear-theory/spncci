  /****************************************************************
  u3boson.h

  Define matrix elements for operators in U(3)-boson basis

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/6/22 (aem): Created. Extracted from vcs.h
  3/7/22 (aem): Add U3BosonSpace class.

****************************************************************/
#ifndef U3BOSON_H_
#define U3BOSON_H_

#include <unordered_map>
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "basis/operator.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"


namespace u3boson
{
  std::vector<u3::U3> RaisingPolynomialLabels(int Nn_max);
  // Generate full set of raising polynomial U3 labels up to given Nn_max.
  //
  // Labels are generated in "canonical" order, defined as
  // lexicographic by N(lambda,mu).
  //
  // Arguments:
  //   Nn_max : maximum excitation quanta of raising polynomial
  //     labels
  //
  // Returns:
  //   Raising polynomial labels


  class U3Subspace;
  class U3State;
  class U3BosonSpace;
  // Class for U(3) boson space, which is a direct product of the space for a unitary U(3) irrep
  // and the space for a unitary irrep of a Weyl (boson) algebra. See, e.g, jpa-18-1985-939-Rowe
  // or ps-91-2016-033003-Rowe.
  //
  // u3boson::U3BosonSpace [sigma]
  // -> U3Subspace [omega]
  //    -> U3State [n,rho]
  //
  // Note: U3Subspace is BaseDegenerateSubspace.  To get number of "states", use .size() if want
  // Number of actual states (n,rho), use .dimension()

  // Subspace constructor
  class U3Subspace
      : public basis::
            BaseDegenerateSubspace<U3Subspace, std::tuple<u3::U3>, U3State, std::tuple<u3::U3>>
  {
   public:
    U3Subspace() = default;

    U3Subspace(
        const u3::U3& omega,
        const MultiplicityTagged<u3::U3>::vector& nrho_vector
      )
        : BaseDegenerateSubspace{omega}
    {
      for (const auto& [n, rho_max] : nrho_vector) PushStateLabels(n, rho_max);
    }

    u3::U3 omega() const { return std::get<0>(labels()); }
    std::string DebugStr() const;
    std::string LabelStr() const {return omega().Str();}

   private:
  };

  // State constructor
  class U3State
      : public basis::BaseDegenerateState<U3Subspace>
  {
   public:
    // pass-through constructors

    U3State(const SubspaceType& subspace, std::size_t index)
    // Construct state by index.
        : basis::BaseDegenerateState<U3Subspace>(subspace, index)
    {}

    U3State(
        const SubspaceType& subspace,
        const typename SubspaceType::StateLabelsType& state_labels
      )
    // Construct state by reverse lookup on labels.
        : basis::BaseDegenerateState<U3Subspace>(subspace, state_labels)
    {}

    // pass-through accessors for subspace labels
    u3::U3 n() const { return std::get<0>(labels()); }
    int rho_max() const { return subspace().GetStateDegeneracy(index()); }

    // private:
  };

  // Space constructor
  class U3BosonSpace
      : public basis::BaseSpace<U3BosonSpace, U3Subspace>
  {
    public:
      U3BosonSpace() = default;
      U3BosonSpace(const u3::U3& sigma, const int Nn_max);

      u3::U3 sigma() const { return std::get<0>(sigma_); }
      int Nn_max() const {return Nn_max_;}
      std::string DebugStr() const;

    private:
      unsigned int Nn_max_;
      std::tuple<u3::U3> sigma_;

  };


  class U3BosonSectors
    : public basis::BaseSectors<U3BosonSpace>
  {
    public:

    // Default constructor
    U3BosonSectors() = default;

    // Constructor
    U3BosonSectors(
        const U3BosonSpace& space,
        const u3::U3& omega0
      );

  private:
    u3::U3 omega0_;

  };


  double BosonCreationRME(const u3::U3& np, const u3::U3& n);
  // SU(3) Reduced matrix element of a^\dagger boson creation operator
  // in a Weyl-boson basis
  //
  // Based on protoype u3boson.py  Formula is given by:
  //   G. Rosensteel and D. J. Rowe, J. Math. Phys. 24 (1983) 2461.
  //
  // Returns:
  //    rme: (double) reduced matrix element of boson creation operator.

  double U3BosonCreationRME(
    const u3::U3& sigmap, const u3::U3& np, unsigned int rhop, const u3::U3& omegap,
    const u3::U3& sigma,  const u3::U3& n,  unsigned int rho,  const u3::U3& omega
  );
  // SU(3) Reduced matrix element of a^\dagger boson creation operator in the
  // U(3)-boson basis (defined in U3BosonSpace).
  //
  // Based on G. Rosensteel and D. J. Rowe. J. Math Phys. 24 (1983) 2461.


  double U3BosonCreationRME(
  const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop, const u3::U3& omegap,
  const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
  );
  // Convenient overload.  Will be deprecated with new recurrence algorithm.

  double U3BosonAnnihilationRME(
    const u3::U3& sigmap, const u3::U3& np, unsigned int rhop, const u3::U3& omegap,
    const u3::U3& sigma,  const u3::U3& n,  unsigned int rho,  const u3::U3& omega
  );
  // SU(3) Reduced matrix element of boson annihlation operator in the
  // U(3)-boson basis (defined in U3BosonSpace).
  //
  // Based on G. Rosensteel and D. J. Rowe. J. Math Phys. 24 (1983) 2461.


  basis::OperatorBlock<double>
  U3BosonCreationOperator(
    const u3::U3& sigma,
    const u3boson::U3Subspace& bra_subspace,
    const u3boson::U3Subspace& ket_subspace,
    u3::UCoefCache& u_coef_cache
    );

  // Returns SU(3) reduced matrix elements of boson creation operator
  //  for a given pair of U(3) subspaces.


  basis::OperatorBlock<double>
  U3BosonAnnihilationOperator(
    const u3::U3& sigma,
    const u3boson::U3Subspace& bra_subspace,
    const u3boson::U3Subspace& ket_subspace,
    u3::UCoefCache& u_coef_cache
    );

  // Returns SU(3) reduced matrix elements of boson annihilation operator
  //  for a given pair of U(3) subspaces.

}

#endif

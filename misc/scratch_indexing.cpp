  // This basis is for *distinguishable* particle states, subject to
  // the symmetry constraints which apply to *identical*-particle
  // states.  That is, although the antisymmetry constraint for identical
  // orbitals (omega+S+T~1) is applied to the quantum numbers when
  // enumerating the basis, the states
  //    
  //   |(N1,N2)...>  and  |(N2,N1)...>
  //
  // are still counted separately in the basis.  No lexicographical
  // ordering constraint N1<=N2 on the two single-particle states is
  // imposed.  The basis is therefore redundant.  This simplifies
  // implementation of double contractions over "particle 1" and
  // "particle 2" indices as matrix multiplication, without the
  // calling code having to worry about "swapping" single particle
  // states within the two-body state and applying the relevant phase
  // (~omega+S+T+g+1).
  //

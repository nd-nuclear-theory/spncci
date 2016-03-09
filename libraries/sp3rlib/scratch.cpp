  // overload for irreps containing SU(3) information, e.g., U3, U3S, and U3ST
  //
  // Only reports SU(3) info...
  template <typename IRREP>
    inline int OuterMultiplicity(const IRREP& w1, const IRREP& w2, const IRREP& w3)
    {
      return OuterMultiplicity(w1.SU3(),w2.SU3(),w3.SU3());
    };


    HalfInt::pair L_range(std::min(x.lambda,x.mu)%2,x.lambda+x.mu);
    L_range = am::AngularMomentumRangeIntersection(L_range,r);
    int L_min = int(L_range.first);
    int L_max = int(L_range.second);


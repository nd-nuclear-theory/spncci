        // Note: This more compact ostensible approach does not work...
        // "error: invalid initialization of non-const reference..."
        //
        //   std::tie(
        //       std::tie(N0,x0,S0,T0,g0),
        //       rho0,
        //       std::tie(eta1p,eta2p,xp,Sp,Tp),
        //       std::tie(eta1,eta2,x,S,T)
        //     ) = two_body_unit_tensor_labels.Key();


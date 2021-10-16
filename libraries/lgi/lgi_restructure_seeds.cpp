  void 
  restructure_seeds(
    int irrep_family_index_bra, int irrep_family_index_ket, // Index in lgi_vector
    const lgi::MultiplicityTaggedLGIVector& lgi_families,   // LGI vector
    const std::vector<int>& lgi_full_space_index_lookup,
    spncci::OperatorBlocks& lgi_transformations,  // May not be initialized if transfrom_lgi_families=False
    bool transform_lgi_families
  )
  {
    //////////////////////////////////////////////////////////////////////
    // Extract lgi index  labels 
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
    ///////////////////////////////////////////////////////////////////////////
    // Read in list of unit tensors between lgi pair and conjugates from files
    // Returned bool, files_found_test has no current use, but could be used 
    //  to identify lgi pair with no non-zero rmes between them.
    // Corresponding rho0 values stored separately for later hypersector lookup
    ///////////////////////////////////////////////////////////////////////////
    // Initialize containers 
    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensors;
    std::vector<int> rho0_values;

    // Get index corresponding to lgi in the full space.
    // Index may differ from lgi index in basis if space has been truncated
    int index1=lgi_full_space_index_lookup[irrep_family_index_bra];
    int index2=lgi_full_space_index_lookup[irrep_family_index_ket];

    // Read in operators 
    std::string lgi_unit_tensor_filename
      =fmt::format("seeds/operators_{:06d}_{:06d}.dat",index1,index2);
    
    std::vector<int> rho0_values;
    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensors;
    bool files_found_test=lgi::ReadUnitTensorLabels(lgi_unit_tensor_filename,lgi_unit_tensors,rho0_values);

    ///////////////////////////////////////////////////////////////////////////
    // Reads in unit tensor seed blocks and stores them in a vector of blocks. Order
    // corresponds to order of (unit_tensor,rho0) pairs in corresponding operator file. 
    //////////////////////////////////////////////////////////////////////////

    // Reads in unit tensor seed blocks and stores them in a vector of blocks. Order
    // corresponds to order of (unit_tensor,rho0) pairs in corresponding operator file. 
    basis::OperatorBlocks<double> unit_tensor_seed_blocks;
    std::string seed_filename
      =fmt::format("seeds/seeds_{:06d}_{:06d}.rmes",index1,index2);
    basis::OperatorBlocks<double> unit_tensor_seed_blocks;
    files_found_test&=lgi::ReadBlocks(seed_filename, lgi_unit_tensors.size(), unit_tensor_seed_blocks);

    // If transform_lgi_families=True, apply basis transformation to lgi
    // TODO: figure out how this works...
    if(transform_lgi_families)
      spncci::TransformSeeds(index1,index2,lgi_transformations,unit_tensor_seed_blocks);
  
    const auto& [lgi_bra,gamma_max_bra]=lgi_families[irrep_family_index_bra];
    const auto& [lgi_ket,gamma_max_ket]=lgi_families[irrep_family_index_ket];

    const auto& [sigma_bra,Sp_bra,Sn_bra,S_bra] = lgi_bra.Key();
    const auto& [sigma_ket,Sp_ket,Sn_ket,S_ket] = lgi_ket.Key();

    // TODO: Will need to be set once indexing set up. 
    int spatial_dim=1;
    int spin_dim=1;
    Eigen::MatrixX::Zero seed_matrix(spatial_dim,spin_dim);

    // std::cout<<"loop over lgi unit tensors"<<std::endl;
    for(int i=0; i<lgi_unit_tensors.size();  ++i)
      {
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor=lgi_unit_tensors[i];
        int rho0=rho0_values[i];
        
        // Extract unit tensor labels 
        const auto & [x0,S0,T0,etap,Sp,Tp,eta,S,T]=unit_tensor.FlatKey();
        const auto& seed_block=unit_tensor_seed_blocks[i];

        // Regroup labels to look up in seed matrix 
        // int spatial_index(sigma,sigma',x0,etap,eta)
        // int spin_index(S_ket,S_bra,Sp_ket,Sp_bra,Sn_ket,Sn_bra,gamma_start,gammap_start,S0,T0,Sp,S,Tp,T)
        // seed_matrix.block(spatial_index,spin_index,1,gamma_max_ket*gamma_max_bra)
        //    = Map<Matrix1d>(seed_block.data(), seed_block.size());

      }
    }

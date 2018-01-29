/****************************************************************
  explicit_construction.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/explicit_construction.h"

#include <omp.h>
#include "cppformat/format.h"
#include "mcutils/eigen.h"
#include "lgi/lgi_solver.h"
#include "lgi/lgi_unit_tensors.h"

namespace spncci
{
  typedef std::vector<basis::MatrixVector> PolynomialMatrices; 

  void
  GenerateSpRaisingPolynomials(
      const spncci::SpNCCIIrrepFamily& sp_irrep_family,
      const u3shell::SpaceU3SPN& lsu3shell_space,
      const u3shell::SectorsU3SPN& Arel_sectors,
      const basis::MatrixVector& Arel_matrices,
      PolynomialMatrices& polynomial_matrices
    )
  {
    
    const sp3r::Sp3RSpace& u3_subspaces=sp_irrep_family.Sp3RSpace();
    u3shell::U3SPN sigmaSPN(sp_irrep_family);
    int lgi_lsu3shell_subspace_index=lsu3shell_space.LookUpSubspaceIndex(sigmaSPN);
    // index starts with 1 because we don't care about lgi nex=0 subspace (i=0 subspace)
    polynomial_matrices.resize(u3_subspaces.size());
    for(int w=1; w<u3_subspaces.size(); ++w)
      {
        const sp3r::U3Subspace&  omega_subspace=u3_subspaces.GetSubspace(w);
        u3::U3 omega(omega_subspace.labels());
        MultiplicityTagged<u3::SU3>::vector omegapp_list=KroneckerProduct(omega.SU3(), u3::SU3(0,2));
        // get index for omega in lsu3shell basis
        int omega_lsu3shell_subspace_index
          =lsu3shell_space.LookUpSubspaceIndex(u3shell::U3SPN(omega,sigmaSPN.Sp(),sigmaSPN.Sn(),sigmaSPN.S()));

        // resize matrix vector for omega polynomials 
        polynomial_matrices[w].resize(omega_subspace.size());

        // iterate over n_rho "states" in u3 subspace 
        for(int i=0; i<omega_subspace.size(); ++i)
          {
            MultiplicityTagged<u3::U3> n_rho(omega_subspace.GetStateLabels(i));
            const u3::U3& n=n_rho.irrep;
            bool init=true;
            for(int w_bar=0; w_bar<omegapp_list.size(); ++w_bar)
              {
                
                u3::U3 omega_bar(omega.N()-2,omegapp_list[w_bar].irrep);
                // check if omega_bar is a valid U3 state and in irrep
                if(not omega_bar.Valid())
                  continue;

                int omega_bar_index=u3_subspaces.LookUpSubspaceIndex(omega_bar);

                if(omega_bar_index==basis::kNone)
                  continue;

                // get index for omega in lsu3shell basis
                int omega_bar_lsu3shell_subspace_index
                  =lsu3shell_space.LookUpSubspaceIndex(u3shell::U3SPN(omega_bar,sigmaSPN.Sp(),sigmaSPN.Sn(),sigmaSPN.S()));

                // get Aintr matrix 
                int Arel_sector_index
                  =Arel_sectors.LookUpSectorIndex(omega_lsu3shell_subspace_index,omega_bar_lsu3shell_subspace_index);

                const Eigen::MatrixXd& A=Arel_matrices[Arel_sector_index];
                
                // special case: if Nn=2, then polynomial is A
                if((omega.N()-sigmaSPN.N())==2)
                  { 
                    polynomial_matrices[w][i]=A;
                    continue;
                  }

                const sp3r::U3Subspace& omega_bar_subspace=u3_subspaces.GetSubspace(omega_bar_index);

                // initial polynomial matrix
                if(init)
                  {
                    int rows=A.rows();
                    int cols=polynomial_matrices[omega_bar_index][0].cols();
                    polynomial_matrices[w][i]=Eigen::MatrixXd::Zero(rows,cols);
                    init=false;
                  }

                // sum over n_bar,rho_bar
                for(int j=0; j<omega_bar_subspace.size(); ++j)
                  {
                    u3::U3 n_bar(omega_bar_subspace.GetStateLabels(j).irrep);
                    int rho_bar=omega_bar_subspace.GetStateLabels(j).tag;

                    double 
                    coef=2./int(n.N())*vcs::BosonCreationRME(n,n_bar)
                          *u3::U(sigmaSPN.U3().SU3(),n_bar.SU3(), omega.SU3(),u3::SU3(2,0),
                                  omega_bar.SU3(),rho_bar,1,n.SU3(),1,1);
                  
                    // std::cout<<polynomial_matrices[omega_bar_index][j].rows()<<std::endl<<std::endl;
                    // std::cout<<"A "<<A.rows()<<" "<<A.cols()<<std::endl<<A<<std::endl;
                    // std::cout<<polynomial_matrices[w][i].rows()<<std::endl<<std::endl;
                    polynomial_matrices[w][i]+=coef*A*polynomial_matrices[omega_bar_index][j];
                  }
              }
          }
      }
  }


  void 
  ConstructSpNCCIBasisExplicit(
      const u3shell::SpaceU3SPN& lsu3shell_space,
      const spncci::SpNCCISpace& sp_irrep_families,
      const lgi::MultiplicityTaggedLGIVector& lgi_families,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const spncci::KMatrixCache& kinv_matrix_cache,
      const u3shell::SectorsU3SPN& Arel_sectors,
      const basis::MatrixVector& Arel_matrices,
      basis::MatrixVector& spncci_expansions,
      bool restrict_sp3r_u3_branching
    )
  {
    // Get expansion from file
    basis::MatrixVector lgi_expansions; 
    std::string lgi_expansion_filename="lgi_expansions.dat";
    lgi::ReadBlocks(lgi_expansion_filename,lgi_families.size(),lgi_expansions);

    u3::UCoefCache u_coef_cache;
    // Generate raising polynomials
    // std::cout<<"generate polynomial matrices"<<std::endl;
    std::vector<PolynomialMatrices> polynomial_matrices_by_lgi_family(sp_irrep_families.size()); 
    for(int p=0; p<sp_irrep_families.size(); ++p)
      {
        const spncci::SpNCCIIrrepFamily& sp_irrep_family=sp_irrep_families[p];
        spncci::GenerateSpRaisingPolynomials(
          sp_irrep_family,lsu3shell_space,
          Arel_sectors,Arel_matrices,
          polynomial_matrices_by_lgi_family[p]
          );
      }
    // std::cout<<"get spncci expansions "<<std::endl;
    spncci_expansions.resize(baby_spncci_space.size());
    for (int subspace_index=0; subspace_index<baby_spncci_space.size(); ++subspace_index)
      {
        // extract subspace properties
        const BabySpNCCISubspace& subspace = baby_spncci_space.GetSubspace(subspace_index);
        int irrep_family_index = subspace.irrep_family_index();
        int Nex = int(subspace.omega().N()-subspace.sigma().N());

        // extract sp3r properties
        const sp3r::Sp3RSpace& u3_subspaces=sp_irrep_families[irrep_family_index].Sp3RSpace();
        int omega_index=u3_subspaces.LookUpSubspaceIndex(subspace.omega());
        const auto& u3_subspace=u3_subspaces.GetSubspace(omega_index); 
        basis::MatrixVector& raising_polynomials=polynomial_matrices_by_lgi_family[irrep_family_index][omega_index];

        // define aliases to the relevant lsu3shell subspaces
        int lgi_lsu3shell_subspace_index = lsu3shell_space.LookUpSubspaceIndex(subspace.sigmaSPN());
        int spncci_lsu3shell_subspace_index = lsu3shell_space.LookUpSubspaceIndex(subspace.omegaSPN());

        // diagnostics
        const u3shell::SubspaceU3SPN& lgi_lsu3shell_subspace = lsu3shell_space.GetSubspace(lgi_lsu3shell_subspace_index);
        const u3shell::SubspaceU3SPN& spncci_lsu3shell_subspace = lsu3shell_space.GetSubspace(spncci_lsu3shell_subspace_index);
        // std::cout
        //   << fmt::format(
        //       "Constructing subspace: LGI sigmaSPN {} family index {} => omegaSPN {} (Nex {})",
        //       subspace.sigmaSPN().Str(),irrep_family_index,
        //       subspace.omegaSPN().Str(),spncci_lsu3shell_subspace.size(),Nex
        //     )
        //   << std::endl
        //   << fmt::format(
        //       "  lsu3shell subspace sigmaSPN: U3SPN {} subspace index {} dim {}",
        //       lgi_lsu3shell_subspace.U3SPN().Str(),
        //       lgi_lsu3shell_subspace_index,
        //       lgi_lsu3shell_subspace.size()
        //     )
        //   << std::endl
        //   << fmt::format(
        //       "  lsu3shell subspace omegaSPN: U3SPN {} subspace index {} dim {}",
        //      spncci_lsu3shell_subspace.U3SPN().Str(),
        //      spncci_lsu3shell_subspace_index,
        //      spncci_lsu3shell_subspace.size()
        //     )
        //   << std::endl;

        // define aliases to expansion matrices
        // std::cout<<"get lgi expansion "<<std::endl;
        // const Eigen::MatrixXd& lgi_expansion = lgi_expansions[irrep_family_index];
        const Eigen::MatrixXd& lgi_expansion = lgi_expansions[lgi_lsu3shell_subspace_index];
        Eigen::MatrixXd& spncci_expansion = spncci_expansions[subspace_index];

        // diagnostics
        // std::cout << fmt::format("  lgi_expansion ({},{})",lgi_expansion.rows(),lgi_expansion.cols()) << std::endl;

        // calculate expansion of baby SpNCCI subspace
        // assert((Nex==0)||(Nex==2));
        if (Nex==0)
          // Nex=0 -- trivial expansion
          {
            spncci_expansion = lgi_expansion;
          }
        else
          {
            int upsilon_max=subspace.upsilon_max();
            int u3_subspace_dim=u3_subspace.size();
            int gamma_max=subspace.gamma_max();
            int rows=raising_polynomials[0].rows();
            int cols=lgi_expansion.cols();

            spncci_expansion=Eigen::MatrixXd::Zero(rows,upsilon_max*gamma_max);
            Eigen::MatrixXd k_matrix_inverse = kinv_matrix_cache.at(subspace.sigma()).at(subspace.omega());
            
            int phase = ParitySign(
                u3::ConjugationGrade(subspace.sigma().SU3())
                +u3::ConjugationGrade(subspace.omega().SU3())
              );
            assert(raising_polynomials.size()==u3_subspace_dim);
            for(int u=0; u<upsilon_max; ++u)
            {
              Eigen::MatrixXd omega_expansion_temp=Eigen::MatrixXd::Zero(rows,cols);
              for(int u_bar=0; u_bar<u3_subspace_dim; ++u_bar)
              {
                // std::cout<<u<<" "<<u_bar<<" "<<k_matrix_inverse.rows()<<" "<<k_matrix_inverse.cols()
                // <<" "<<raising_polynomials[u_bar].rows()<<"  "<<raising_polynomials[u_bar].rows()<<std::endl;
                omega_expansion_temp+=phase*k_matrix_inverse(u_bar,u)*raising_polynomials[u_bar]*lgi_expansion;
              }

              // Populating spncci_expansion with omega expansions. 
              // omega_expansion_temp
              // u=upsilon-1 
              // g=gamma-1
              // index in expansion is upsilon_max(gamma-1)+(upsilon-1)
              for(int g=0; g<gamma_max; ++g)
                {
                  int column_index=g*upsilon_max+u;
                  spncci_expansion.block(0,column_index,rows,1)=omega_expansion_temp.block(0,g,rows,1);
                }

            }
          }
      }
  }


  void 
  ComputeUnitTensorSectorsExplicit(
    const u3::U3& sigmap, const u3::U3& sigma,
    const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const u3shell::SpaceU3SPN& lsu3shell_space,
    const u3shell::SectorsU3SPN& lsu3shell_operator_sectors,
    basis::MatrixVector& lsu3shell_operator_matrices,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const basis::MatrixVector& spncci_expansions,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
    )
    {
      // #pragma omp parallel
        
      // #pragma omp for schedule(runtime)
     
      // for each of the unit tensors in lsu3shell basis

      for(int s=0; s<lsu3shell_operator_sectors.size(); ++s)
        {
          // extract U3SPN labels
          auto& lsu3shell_sector=lsu3shell_operator_sectors.GetSector(s);
          int lsu3shell_bra_index=lsu3shell_sector.bra_subspace_index();
          int lsu3shell_ket_index=lsu3shell_sector.ket_subspace_index();
          int rho0=lsu3shell_sector.multiplicity_index();

          u3shell::U3SPN lsu3shell_ket_subspace_labels=lsu3shell_space.GetSubspace(lsu3shell_ket_index).labels();
          u3shell::U3SPN lsu3shell_bra_subspace_labels=lsu3shell_space.GetSubspace(lsu3shell_bra_index).labels();

          // look up unit tensor subspace index 
          u3shell::UnitTensorSubspaceLabels unit_tensor_labels(unit_tensor.x0(),unit_tensor.S0(),unit_tensor.bra().eta(),unit_tensor.ket().eta());
          int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels);
          
          // look up unit tensor index
          std::tuple<int,int,int,int,int> unit_tensor_state_labels(
              int(unit_tensor.T0()), int(unit_tensor.bra().S()),
              int(unit_tensor.bra().T()),int(unit_tensor.ket().S()),int(unit_tensor.ket().T())
            );
          int unit_tensor_index=unit_tensor_space.GetSubspace(unit_tensor_subspace_index).LookUpStateIndex(unit_tensor_state_labels);
          
          // Look up baby spncci index
          spncci::BabySpNCCISubspaceLabels baby_spncci_labels_bra(
              sigmap,
              lsu3shell_bra_subspace_labels.Sp(),
              lsu3shell_bra_subspace_labels.Sn(),
              lsu3shell_bra_subspace_labels.S(),
              lsu3shell_bra_subspace_labels.U3()
            );
          int baby_spncci_subspace_index_bra=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labels_bra);
            
          // If baby spncci index=-1, then omega is not in sigma irrep
          // so go to next omega.   
          if(baby_spncci_subspace_index_bra==-1)
            continue;

          spncci::BabySpNCCISubspaceLabels baby_spncci_labels_ket(
              sigma,
              lsu3shell_ket_subspace_labels.Sp(),
              lsu3shell_ket_subspace_labels.Sn(),
              lsu3shell_ket_subspace_labels.S(),
              lsu3shell_ket_subspace_labels.U3()
            );

            int baby_spncci_subspace_index_ket=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_labels_ket);
            if(baby_spncci_subspace_index_ket==-1)
              continue;

            // look up hyper sector index
            int hypersector_index
                  =baby_spncci_hypersectors.LookUpHypersectorIndex(
                    baby_spncci_subspace_index_bra,baby_spncci_subspace_index_ket,
                    unit_tensor_subspace_index, rho0
                  );
            
            // If not in hypersector set, continue.
            if(hypersector_index==-1)
              continue;

            // std::cout<<fmt::format("unit tensor subspace {}, unit tensor {}, baby spncci bra {} ket{} hypersector{}",
            //   unit_tensor_subspace_index,unit_tensor_index,baby_spncci_subspace_index_bra, baby_spncci_subspace_index_ket,hypersector_index)<<std::endl;

            // Get bra and ket lsu3shell expansion and compute unit tensor block
            const Eigen::MatrixXd& bra_expansion=spncci_expansions[baby_spncci_subspace_index_bra];
            const Eigen::MatrixXd& ket_expansion=spncci_expansions[baby_spncci_subspace_index_ket];
            Eigen::MatrixXd temp=bra_expansion.transpose()*lsu3shell_operator_matrices[s]*ket_expansion;
            // std::cout<<temp<<std::endl;

            // std::cout<<fmt::format(" bra : {} x {}  operator : {} x {}  ket : {} x {}",
            //   bra_expansion.rows(), bra_expansion.cols(), 
            //   lsu3shell_operator_matrices[s].rows(), lsu3shell_operator_matrices[s].cols(),
            //   ket_expansion.rows(), ket_expansion.cols()
            //   )<<std::endl;

            // Store unit tensor block in hypersector structure
            unit_tensor_hyperblocks[hypersector_index][unit_tensor_index]+=temp;
            // std::cout<<"hypersector "<<hypersector_index<<"  tensor "<<unit_tensor_index<<std::endl;
            // std::cout<<"temp "<<temp<<" in hypersectors "<<unit_tensor_hyperblocks[hypersector_index][unit_tensor_index]<<std::endl<<std::endl;
        }          
    }

  void 
  CheckUnitTensorRecurrence(
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_explicit
    )
    {
      // Comparing recurrence to explicit hyperblocks and printing error message if difference between
      // hyperblocks exceeds tolerance of 1e-4   
      bool errors=false;
      for(int i=0; i<unit_tensor_hyperblocks.size(); ++i)
        for(int j=0; j<unit_tensor_hyperblocks[i].size(); ++j)
          {
            auto& hypersector=baby_spncci_hypersectors.GetHypersector(i);
            int bra, ket, tensor, rho0;
            std::tie(bra,ket,tensor,rho0)=hypersector.Key();
            auto& bra_subspace=baby_spncci_space.GetSubspace(bra);
            auto& ket_subspace=baby_spncci_space.GetSubspace(ket);
            auto& tensor_subspace=unit_tensor_space.GetSubspace(tensor);
            int Sp,Tp,S,T,T0;
            std::tie(T0,Sp,Tp,S,T)=tensor_subspace.GetStateLabels(j);
            
            const Eigen::MatrixXd matrix1=unit_tensor_hyperblocks[i][j];

            // std::cout<<bra_subspace.LabelStr()<<"  "<<ket_subspace.LabelStr()<<"  "<<tensor_subspace.LabelStr()<<"  "
            //              << rho0<<std::endl;
            // std::cout<<"   "<<T0<<"  "<<Sp<<"  "<<Tp<<"  "<<S<<"  "<<T<<std::endl;

            const Eigen::MatrixXd matrix2=unit_tensor_hyperblocks_explicit[i][j];

            if(not mcutils::IsZero(matrix1-matrix2, 1e-4))
              {
                errors=true;
                std::cout<<"hyperblock "<<i<<" sub-block "<<j<<" is not correct"<<std::endl;
                std::cout<<bra_subspace.LabelStr()<<"  "<<ket_subspace.LabelStr()<<"  "<<tensor_subspace.LabelStr()<<"  "
                         << rho0<<std::endl;
                std::cout<<"the matrix should be "<<bra_subspace.size()<<" x "<<ket_subspace.size()<<std::endl;
                std::cout<<"gammma_max: "<<ket_subspace.gamma_max()<<" upsilon_max "<<ket_subspace.upsilon_max()<<std::endl;
                std::cout<<"matrix1"<<std::endl<<matrix1<<std::endl<<"matrix2"
                         <<std::endl<<matrix2<<std::endl;
              }
          }
      // If no error found, print no errors.
      assert(not errors);
      if(not errors)
        std::cout<<"no errors"<<std::endl;

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Checking the other direction 

      // Comparing recurrence to explicit hyperblocks and printing error message if difference between
      // hyperblocks exceeds tolerance of 1e-4   
      errors=false;
      for(int i=0; i<unit_tensor_hyperblocks_explicit.size(); ++i)
        for(int j=0; j<unit_tensor_hyperblocks_explicit[i].size(); ++j)
          {
            auto& hypersector=baby_spncci_hypersectors.GetHypersector(i);
            int bra, ket, tensor, rho0;
            std::tie(bra,ket,tensor,rho0)=hypersector.Key();
            auto& bra_subspace=baby_spncci_space.GetSubspace(bra);
            auto& ket_subspace=baby_spncci_space.GetSubspace(ket);
            auto& tensor_subspace=unit_tensor_space.GetSubspace(tensor);
            int Sp,Tp,S,T,T0;
            std::tie(T0,Sp,Tp,S,T)=tensor_subspace.GetStateLabels(j);
            
            const Eigen::MatrixXd matrix1=unit_tensor_hyperblocks_explicit[i][j];

            // std::cout<<bra_subspace.LabelStr()<<"  "<<ket_subspace.LabelStr()<<"  "<<tensor_subspace.LabelStr()<<"  "
            //              << rho0<<std::endl;
            // std::cout<<"   "<<T0<<"  "<<Sp<<"  "<<Tp<<"  "<<S<<"  "<<T<<std::endl;

            const Eigen::MatrixXd matrix2=unit_tensor_hyperblocks[i][j];

            if(not mcutils::IsZero(matrix1-matrix2, 1e-4))
              {
                errors=true;
                std::cout<<"hyperblock "<<i<<" sub-block "<<j<<" is not correct"<<std::endl;
                std::cout<<bra_subspace.LabelStr()<<"  "<<ket_subspace.LabelStr()<<"  "<<tensor_subspace.LabelStr()<<"  "
                         << rho0<<std::endl;
                std::cout<<"the matrix should be "<<bra_subspace.size()<<" x "<<ket_subspace.size()<<std::endl;
                std::cout<<"gammma_max: "<<ket_subspace.gamma_max()<<" upsilon_max "<<ket_subspace.upsilon_max()<<std::endl;
                std::cout<<"matrix1"<<std::endl<<matrix1<<std::endl<<"matrix2"
                         <<std::endl<<matrix2<<std::endl;
              }
          }
      // If no error found, print no errors.
      assert(not errors);
      if(not errors)
        std::cout<<"no errors"<<std::endl;

    }


  void ExplicitlyComputeUnitTensorHyperBlocks(
      const u3::U3& sigmap, const u3::U3& sigma,
      const spncci::RunParameters& run_parameters,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const basis::MatrixVector& spncci_expansions,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_explicit
    )
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // computing unit tensors
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     // read lsu3shell basis (regroup into U3SPN subspaces)
    lsu3shell::LSU3ShellBasisTable lsu3shell_basis_table;
    lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
    u3shell::SpaceU3SPN lsu3shell_space;
    lsu3shell::ReadLSU3ShellBasis(
        run_parameters.Nsigma0,run_parameters.lsu3shell_basis_filename,
        lsu3shell_basis_table,lsu3shell_basis_provenance,lsu3shell_space
      );


    // Zero initialize hypersectors 
    // basis::OperatorHyperblocks<double> unit_tensor_hyperblocks_explicit;
    basis::SetHyperoperatorToZero(baby_spncci_hypersectors,unit_tensor_hyperblocks_explicit);
    
    const std::string relative_unit_tensor_filename_template=run_parameters.relative_unit_tensor_filename_template;

    // for each unit tensor
    for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
      {
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor = lgi_unit_tensor_labels[unit_tensor_index];

        // get unit tensor labels 
        u3::SU3 x0; 
        HalfInt S0,T0,Sp,Tp,S,T;
        int etap,eta;
        std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=unit_tensor.FlatKey();

        // Look up unit tensor subspace
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
        int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
        auto& subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        // Get unit tensor index in subspace 
        int unit_tensor_state_index
          =subspace.LookUpStateIndex(std::tuple<int,int,int,int,int>(int(T0), int(Sp),int(Tp),int(S),int(T)));

        // Construct lsu3shell sectors for unit tensor
        const bool spin_scalar = false;
        u3shell::SectorsU3SPN unit_tensor_sectors;
        unit_tensor_sectors = u3shell::SectorsU3SPN(lsu3shell_space,unit_tensor,spin_scalar);
        
        // read in lsu3shell rms for unit tensor 
        basis::MatrixVector unit_tensor_lsu3shell_blocks;
        std::string filename = fmt::format(relative_unit_tensor_filename_template,unit_tensor_index);
        lsu3shell::ReadLSU3ShellRMEs(
            filename,
            lsu3shell_basis_table,lsu3shell_space,
            unit_tensor,unit_tensor_sectors,unit_tensor_lsu3shell_blocks
          );

        // Compute unit tensor hyperblocks from lsu3shell rmes using explicit basis construction
        spncci::ComputeUnitTensorSectorsExplicit(
            sigmap, sigma, unit_tensor,unit_tensor_space,
            lsu3shell_space,unit_tensor_sectors,unit_tensor_lsu3shell_blocks,
            baby_spncci_space, spncci_expansions,baby_spncci_hypersectors,
            unit_tensor_hyperblocks_explicit
          );

        // Compute conjugate unit tensor hyperblocks if not diagonal sectors 
        if(not (sigmap==sigma))
          {
            spncci::ComputeUnitTensorSectorsExplicit(
                sigma, sigmap, unit_tensor,unit_tensor_space,
                lsu3shell_space,unit_tensor_sectors,unit_tensor_lsu3shell_blocks,
                baby_spncci_space, spncci_expansions,baby_spncci_hypersectors,
                unit_tensor_hyperblocks_explicit
              );
          }
      } // end unit tensor index 
  }

void ExplicitBasisConstruction(
      const spncci::RunParameters& run_parameters,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      spncci::KMatrixCache& k_matrix_cache, 
      spncci::KMatrixCache& kinv_matrix_cache,
      bool restrict_sp3r_to_u3_branching,
      basis::MatrixVector& spncci_expansions
    )
  {

    std::cout << "Read lsu3shell basis..." << std::endl;
    // read lsu3shell basis (regroup into U3SPN subspaces)
    lsu3shell::LSU3ShellBasisTable lsu3shell_basis_table;
    lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
    u3shell::SpaceU3SPN lsu3shell_space;
    lsu3shell::ReadLSU3ShellBasis(
        run_parameters.Nsigma0,run_parameters.lsu3shell_basis_filename,
        lsu3shell_basis_table,lsu3shell_basis_provenance,lsu3shell_space
      );

    std::cout << "Solve for LGIs..." << std::endl;

    // lgi::MultiplicityTaggedLGIVector lgi_families;
    // basis::MatrixVector lgi_expansions;
    
    lgi::MultiplicityTaggedLGIVector lgi_families;
    std::string lgi_filename="lgi_families.dat";
    lgi::ReadLGISet(lgi_filename, run_parameters.Nsigma0,lgi_families);

    // lgi::GetLGIExpansion(
    //     lsu3shell_space,lsu3shell_basis_table,
    //     run_parameters.Brel_filename,run_parameters.Nrel_filename,
    //     run_parameters.A, run_parameters.Nsigma0,
    //     lgi_families, lgi_expansions
    //   );

    u3shell::SectorsU3SPN Aintr_sectors;
    basis::MatrixVector Aintr_matrices;
    lsu3shell::ReadLSU3ShellSymplecticRaisingOperatorRMEs(
        lsu3shell_basis_table,lsu3shell_space, 
        run_parameters.Arel_filename,Aintr_sectors,Aintr_matrices,
        run_parameters.A
      );



    spncci::ConstructSpNCCIBasisExplicit(
        lsu3shell_space,spncci_space,lgi_families,baby_spncci_space,
        k_matrix_cache,kinv_matrix_cache,Aintr_sectors,Aintr_matrices,spncci_expansions,
        restrict_sp3r_to_u3_branching
      );

  }

  void CheckHyperBlocks(
      int irrep_family_index_bra, int irrep_family_index_ket,
      const spncci::RunParameters& run_parameters,
      const spncci::SpNCCISpace& spncci_space,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const basis::MatrixVector& spncci_expansions,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
   )
  {
      basis::OperatorHyperblocks<double> unit_tensor_hyperblocks_explicit;

      const u3::U3& sigmap=spncci_space[irrep_family_index_bra].sigma();
      const u3::U3& sigma =spncci_space[irrep_family_index_ket].sigma();


      spncci::ExplicitlyComputeUnitTensorHyperBlocks(
          sigmap, sigma,run_parameters, 
          unit_tensor_space,lgi_unit_tensor_labels,
          baby_spncci_space, spncci_expansions,
          baby_spncci_hypersectors,
          unit_tensor_hyperblocks_explicit
        );


        std::cout<<"print hypersectors explicit"<<std::endl;
        spncci::PrintHypersectors(
          baby_spncci_space,unit_tensor_space, 
          baby_spncci_hypersectors,unit_tensor_hyperblocks_explicit
          );



      std::cout<<"checking hypersectors"<<std::endl;
      spncci::CheckUnitTensorRecurrence(
       baby_spncci_space,unit_tensor_space,baby_spncci_hypersectors,
       unit_tensor_hyperblocks,unit_tensor_hyperblocks_explicit
      );
  }



}  // namespace

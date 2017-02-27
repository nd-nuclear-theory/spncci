/****************************************************************
  explicit_construction.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/explicit_construction.h"

#include <omp.h>
#include "cppformat/format.h"
#include "mcutils/eigen.h"

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
        u3::U3 omega(omega_subspace.GetSubspaceLabels());
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
      const basis::MatrixVector& lgi_expansions,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const u3shell::SectorsU3SPN& Arel_sectors,
      const basis::MatrixVector& Arel_matrices,
      basis::MatrixVector& spncci_expansions
    )
  {
   
    u3::UCoefCache u_coef_cache;
    // Generate raising polynomials
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
        const Eigen::MatrixXd& lgi_expansion = lgi_expansions[irrep_family_index];
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
            int gamma_max=subspace.gamma_max();
            int rows=raising_polynomials[0].rows();
            int cols=lgi_expansion.cols();
            
            spncci_expansion=Eigen::MatrixXd::Zero(rows,upsilon_max*gamma_max);
            Eigen::MatrixXd k_matrix_inverse = k_matrix_cache.at(subspace.sigma()).at(subspace.omega()).inverse();
            
            int phase = ParitySign(
                u3::ConjugationGrade(subspace.sigma().SU3())
                +u3::ConjugationGrade(subspace.omega().SU3())
              );

            assert(raising_polynomials.size()==upsilon_max);
            for(int u=0; u<upsilon_max; ++u)
            {
              Eigen::MatrixXd omega_expansion_temp=Eigen::MatrixXd::Zero(rows,cols);
              for(int u_bar=0; u_bar<raising_polynomials.size(); ++u_bar)
              {
                omega_expansion_temp+=phase*k_matrix_inverse(u,u_bar)*raising_polynomials[u_bar]*lgi_expansion;
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
    const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor,
    const u3shell::SpaceU3SPN& lsu3shell_space,
    const u3shell::SectorsU3SPN& lsu3shell_operator_sectors,
    basis::MatrixVector& lsu3shell_operator_matrices,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const basis::MatrixVector& spncci_expansions,
    spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_sectors_explicit
    )
    {
      // #pragma omp parallel
      {
        // #pragma omp for schedule(runtime)
        for(int s=0; s<lsu3shell_operator_sectors.size(); ++s)
          {
            // std::cout<<"sector "<<s<<std::endl;
            auto& lsu3shell_sector=lsu3shell_operator_sectors.GetSector(s);
            int lsu3shell_bra_index=lsu3shell_sector.bra_subspace_index();
            int lsu3shell_ket_index=lsu3shell_sector.ket_subspace_index();
            int rho0=lsu3shell_sector.multiplicity_index();

            u3shell::U3SPN lsu3shell_ket_subspace_labels=lsu3shell_space.GetSubspace(lsu3shell_ket_index).GetSubspaceLabels();
            u3shell::U3SPN lsu3shell_bra_subspace_labels=lsu3shell_space.GetSubspace(lsu3shell_bra_index).GetSubspaceLabels();

            u3::U3 omegap(lsu3shell_bra_subspace_labels.U3());
            u3::U3 omega(lsu3shell_ket_subspace_labels.U3());

            spncci::UnitTensorU3Sector unit_U3Sector(omegap,omega,unit_tensor,rho0);

            // std::cout<<"for each bra "<<std::endl;
            for(int baby_spncci_bra_index=0; baby_spncci_bra_index<baby_spncci_space.size(); ++baby_spncci_bra_index)
              {
                const spncci::BabySpNCCISubspace& baby_spncci_bra_subspace=baby_spncci_space.GetSubspace(baby_spncci_bra_index);

                // check if baby spncci subspace matches the lsu3shell subspace 
                bool allowed_bra;      
                // std::cout<<omegap.Str()<<"  "<<baby_spncci_bra_subspace.omega().Str()<<std::endl;
                allowed_bra=(baby_spncci_bra_subspace.omega()==omegap);
                // std::cout<<"allowed 1 "<<allowed_bra<<std::endl;
                allowed_bra&=(baby_spncci_bra_subspace.Sp()==lsu3shell_bra_subspace_labels.Sp());
                // std::cout<<"allowed 2 "<<allowed_bra<<std::endl;
                allowed_bra&=(baby_spncci_bra_subspace.Sn()==lsu3shell_bra_subspace_labels.Sn());
                // std::cout<<"allowed 3 "<<allowed_bra<<std::endl;
                allowed_bra&=(baby_spncci_bra_subspace.S()==lsu3shell_bra_subspace_labels.S());
                // std::cout<<"allowed 4 "<<allowed_bra<<std::endl;
                if(not allowed_bra)
                  continue;
                // std::cout<<"for each ket "<<std::endl;
                for(int baby_spncci_ket_index=0; baby_spncci_ket_index<baby_spncci_space.size(); ++baby_spncci_ket_index)
                  {
                    const spncci::BabySpNCCISubspace& baby_spncci_ket_subspace=baby_spncci_space.GetSubspace(baby_spncci_ket_index);

                    // check if baby spncci subspace matches the lsu3shell subspace 
                    bool allowed_ket;      
                    allowed_ket=(baby_spncci_ket_subspace.omega()==omega);
                    allowed_ket&=(baby_spncci_ket_subspace.Sp()==lsu3shell_ket_subspace_labels.Sp());
                    allowed_ket&=(baby_spncci_ket_subspace.Sn()==lsu3shell_ket_subspace_labels.Sn());
                    allowed_ket&=(baby_spncci_ket_subspace.S()==lsu3shell_ket_subspace_labels.S());
                  
                    if(not allowed_ket)
                      continue;

                    // std::cout<<"get lgi information"<<std::endl;
                    // lgi information
                    int irrep_family_bra_index=baby_spncci_bra_subspace.irrep_family_index();
                    int irrep_family_ket_index=baby_spncci_ket_subspace.irrep_family_index();
                    std::pair<int,int> lgi_family_pair(irrep_family_bra_index,irrep_family_ket_index);

                    // NnN information
                    int Nnp=int(omegap.N()-baby_spncci_bra_subspace.sigma().N());
                    int Nn=int(omega.N()-baby_spncci_ket_subspace.sigma().N());
                    std::pair<int,int>NnpNn(Nnp,Nn);
                    
                    // Get omega and omega' expansions 
                    // std::cout<<"computing temp"<<std::endl;
                    const Eigen::MatrixXd& bra_expansion=spncci_expansions[baby_spncci_bra_index];
                    const Eigen::MatrixXd& ket_expansion=spncci_expansions[baby_spncci_ket_index];
                    Eigen::MatrixXd temp=bra_expansion.transpose()*lsu3shell_operator_matrices[s]*ket_expansion;
                    // std::cout<<temp<<std::endl;

                    if(not CheckIfZeroMatrix(temp))
                      // #pragma omp critical
                      {
                        unit_tensor_sectors_explicit[lgi_family_pair][NnpNn][unit_U3Sector]=temp;
                      }
                  }
              }
          }
        }          
    }

}  // namespace

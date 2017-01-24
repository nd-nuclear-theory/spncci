/****************************************************************
  explicit_basis_construction.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/17/17 (aem): Created.
****************************************************************/
// Generating the symplectic basis by explicit laddering with
// raising polynomials 

#include <fstream>
#include <unordered_set>
#include "cppformat/format.h"

#include "lgi/lgi_solver.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h"
#include "spncci/sp_basis.h"
#include "u3shell/unit_tensor_expansion.h"
#include "utilities/utilities.h"



int main(int argc, char **argv)
{
  u3::U3CoefInit();
  double zero_threshold=1e-6;

  if(argc<8)
    std::cout<<"Syntax : A  twice_Nsigma_0  Nmax  N1B  <basis>  <nrel>  <brel>"<<std::endl;

	// Extract input values
  int A=std::stoi(argv[1]);
  int twice_Nsigma_0=std::stoi(argv[2]);
  HalfInt Nsigma_0=HalfInt(twice_Nsigma_0,2);
  int Nmax=std::stoi(argv[3]);
  int N1B=std::stoi(argv[4]);
  std::string lsu3_filename = argv[5];
  std::string nrel_filename = argv[6];
  std::string brel_filename = argv[7];
  std::string arel_filename = argv[8];
  
  // Set up lsu3shell basis
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);

  // Get lgi expansions
  std::cout<<"Generating LGI expansion"<<std::endl;
  basis::MatrixVector lgi_expansion_matrix_vector;
  std::ifstream is_nrel(nrel_filename.c_str());
  std::ifstream is_brel(brel_filename.c_str());
  lgi::LGIVector lgi_vector;
  bool keep_empty_subspaces=true;
  lgi::GenerateLGIExpansion(A,Nsigma_0,basis_table,space, is_brel,is_nrel,
    lgi_vector,lgi_expansion_matrix_vector,keep_empty_subspaces);
  is_nrel.close();
  is_brel.close();

  // for(auto& matrix : lgi_expansion_matrix_vector)
  //   std::cout<<matrix<<std::endl<<std::endl;

  // Read in Arel
  basis::MatrixVector Arel_matrices;
  u3shell::OperatorLabelsU3ST arel_labels(2,u3::SU3(2,0),0,0,0);
  u3shell::SectorsU3SPN arel_sectors(space,arel_labels,true);

  std::ifstream is_arel(arel_filename.c_str());
  if(not is_arel)
    std::cout<<fmt::format("{} not found",arel_filename)<<std::endl;

  lsu3shell::ReadLSU3ShellRMEs(
    is_arel,arel_labels,basis_table,space,
    arel_sectors,Arel_matrices);
  is_arel.close();
  //Read in Brel
  basis::MatrixVector Brel_matrices;
  u3shell::OperatorLabelsU3ST brel_labels(-2,u3::SU3(0,2),0,0,0);
  u3shell::SectorsU3SPN brel_sectors(space,brel_labels,true);

  std::ifstream is_brel2(brel_filename.c_str());
  if(not is_brel2)
    std::cout<<fmt::format("{} not found",brel_filename)<<std::endl;  

  lsu3shell::ReadLSU3ShellRMEs(
    is_brel2,brel_labels,basis_table,space,
    brel_sectors,Brel_matrices);
  is_brel2.close();



  for(auto& matrix :Arel_matrices)
    std::cout<<matrix<<std::endl<<std::endl;

  for(int i=0; i<arel_sectors.size(); ++i)
  {
    auto& arel_sector=arel_sectors.GetSector(i);
    int lgi_index=arel_sector.ket_subspace_index();
    int w_index=arel_sector.bra_subspace_index();
    Eigen::MatrixXd& Arel=Arel_matrices[i];
    Eigen::MatrixXd& LGIs=lgi_expansion_matrix_vector[lgi_index];
    // std::cout<<"Arel "<<std::endl;
    // std::cout<<Arel<<std::endl;
    // std::cout<<"LGI "<<std::endl;
    // std::cout<<LGIs<<std::endl;
    int j=brel_sectors.LookUpSectorIndex(lgi_index,w_index,1);
    Eigen::MatrixXd& Brel=Brel_matrices[j];
    // std::cout<<"Brel "<<std::endl;
    // std::cout<<Brel<<std::endl;
    Eigen::MatrixXd omega_matrix=Arel*LGIs;
    // std::cout<<"omegas "<<std::endl;
    // std::cout<<omega_matrix<<std::endl;
    // std::cout<<"Brel on omega_matrix "<<std::endl;
    // std::cout<<Brel*omega_matrix<<std::endl;
    std::cout<<space.GetSubspace(lgi_index).GetSubspaceLabels().Str()<<std::endl;
    std::cout<<LGIs.transpose()*Brel*omega_matrix<<std::endl<<std::endl;

  }

}//End Program


//   // for(int i=0; i<lgi_vector.size(); ++i)
//   //   {
//   //     std::cout<<"LGI "<<lgi_vector[i].Str()<<std::endl;
//   //     std::cout<<lgi_expansion_matrix_vector[i]<<std::endl;
//   //     assert(lgi_expansion_matrix_vector[i].rows()!=0);
//   //   }
//   // Get Spncci vector
//   std::cout<<"Generating irreps"<<std::endl;
//   spncci::SigmaIrrepMap sigma_irrep_map;
//   spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
//   spncci::SpIrrepVector sp_irrep_vector;
//   spncci::GenerateSp3RIrreps(lgi_vector,truncator,sp_irrep_vector,sigma_irrep_map);
//   std::cout<<"Sp irrep vector size "<<sp_irrep_vector.size()<<std::endl;

//   //Generate K matrices
//   std::cout<<"Generating K matrices"<<std::endl;
//   std::unordered_map<u3::U3,vcs::MatrixCache,boost::hash<u3::U3>> k_matrix_map;
//   std::unordered_set<u3::U3,boost::hash<u3::U3> >sigma_set;
//   for(int l=0; l<sp_irrep_vector.size(); l++)
//     sigma_set.insert(sp_irrep_vector[l].irrep.sigma());

//   for( const auto& s : sigma_set)
//     {
//       vcs::MatrixCache K_map;
//       vcs::GenerateKMatrices(sigma_irrep_map[s], K_map);
//       // for(auto it=K_map.begin(); it!=K_map.end(); ++it)
//       //   {
//       //     Eigen::MatrixXd& Kmatrix=it->second;
//       //     NormalizeMatrix(Kmatrix, "rows");
//       //     std::cout<<Kmatrix<<std::endl;
//       //     std::cout<<it->second<<std::endl;
//       //   }
//       k_matrix_map[s]=K_map;
//     }












//   std::vector<basis::MatrixVector> sp_irrep_expansion(sp_irrep_vector.size());
//   //Set up irrep matrix indexing and subspace
//   std::cout<<"zero initialize matrices"<<std::endl;
//   for(int s=0; s<sp_irrep_vector.size(); ++s)
//   	{
//   		auto& lgi_family=sp_irrep_vector[s]; 
//   		const spncci::SpIrrep& sp_irrep=lgi_family.irrep;
//   		const u3::U3& sigma=sp_irrep.sigma();
//   		//TODO may need to extract these labels from LGI rather than sp vector
//   		HalfInt Sp=sp_irrep.Sp();
//   		HalfInt Sn=sp_irrep.Sn();
//   		HalfInt S=sp_irrep.S();
//   		const sp3r::Sp3RSpace& irrep=sp_irrep.Sp3RSpace();
//   		int family_size=lgi_family.tag;
//   		for(int i=0; i<irrep.size(); ++i)
//   		{
//   			u3::U3 omega(irrep.GetSubspace(i).GetSubspaceLabels());
//   			int upsilon_max=irrep.GetSubspace(i).size();
//         int cols=upsilon_max*family_size;
//         u3shell::U3SPN u3spn_labels(u3::U3S(omega,S),Sp,Sn);
//   			int subspace_index=space.LookUpSubspaceIndex(u3spn_labels);
//         int rows=space.GetSubspace(subspace_index).size();
//         sp_irrep_expansion[s].push_back(Eigen::MatrixXd::Zero(rows,cols));
//   		}
//   	}

//   //For each raising polynomial, read in matrix
//   std::cout<<"Iterating over polynomials"<<std::endl;
//   std::vector<u3::U3>raising_polynomials=sp3r::RaisingPolynomialLabels(2);
//   for(int i=0; i<raising_polynomials.size(); ++i)
//     {
//       u3::U3 n(raising_polynomials[i]);
//       u3shell::OperatorLabelsU3ST prel_operator_labels(int(n.N()),n.SU3(),0,0,0);
//       u3shell::SectorsU3SPN prel_sectors(space,prel_operator_labels,true);
//       std::ifstream is_operator(fmt::format("Prel_{:02d}_Nmax{:02d}_{:06d}.rme",A,Nmax,i));

//       if(not is_operator)
//         {
//           std::cout<<fmt::format("Prel_{:02d}_Nmax{:02d}_{:06d}.rme not found",A, Nmax, i)<<std::endl;
//           continue;
//         }
//       basis::MatrixVector lsu3shell_operator_matrices;
//       std::cout<<"Reading in sectors for "<<n.Str()<<std::endl;
//       lsu3shell::ReadLSU3ShellRMEs(is_operator,prel_operator_labels, basis_table,space, prel_sectors,lsu3shell_operator_matrices);
      
//       // iterate over sectors 
//       // Get lgi expansion
//       // std::cout<<"Looping over p sectors"<<std::endl;
//       for(int j=0; j< prel_sectors.size(); ++j)
//         {
//           // std::cout<<"entering p sector "<<j<<std::endl;
//           int bra_index=prel_sectors.GetSector(j).bra_subspace_index();
//           int ket_index=prel_sectors.GetSector(j).ket_subspace_index();
//           int rho=prel_sectors.GetSector(j).multiplicity_index();
//           // Identify the irrep
//           auto& lgi_family=sp_irrep_vector[ket_index]; 
//           const spncci::SpIrrep& sp_irrep=lgi_family.irrep;
//           const u3::U3& sigma=sp_irrep.sigma();
//           const sp3r::Sp3RSpace& irrep=sp_irrep.Sp3RSpace();
//           int family_size=lgi_family.tag;
//           if(family_size==0)
//             continue;
//           // Get omega labels for lookup in irrep
//           // Get sigma to identify irrep
//           u3::U3 omega(space.GetSubspace(bra_index).GetSubspaceLabels().U3());
//           std::cout<<"bra "<<space.GetSubspace(bra_index).GetSubspaceLabels().Str()<<std::endl;
//           // Check that space is still aligned with lgi vector
//           u3::U3 sigma2(space.GetSubspace(ket_index).GetSubspaceLabels().U3());
//           std::cout<<"ket "<<space.GetSubspace(ket_index).GetSubspaceLabels().Str()<<std::endl;
//           // sigma from lgi list and sigma2 from lookup. 
//           assert(sigma==sigma2);
//           // Index in sp_irrep_expansion
//           int omega_index=irrep.LookUpSubspaceIndex(omega);
//           // std::cout<<"somewhere in the middle"<<std::endl;
//           // multiplicity of excited states in each irrep
//           const sp3r::U3Subspace& u3_subspace=irrep.GetSubspace(omega_index);
//           // Get Kmatrix look up index for n and rho 
//           int u=u3_subspace.LookUpStateIndex(MultiplicityTagged<u3::U3>(n,rho));
//           // Get K matrices for lgi family
//           vcs::MatrixCache& k_matrices=k_matrix_map[sigma];
//           // Get lgi expansion 
//           Eigen::MatrixXd& lgi_expansion=lgi_expansion_matrix_vector[ket_index];
//           // std::cout<<lgi_expansion<<std::endl;
//           int num_lgi=lgi_expansion.cols();
//           if (num_lgi==0)
//             continue;
//           Eigen::MatrixXd& prel_operator=lsu3shell_operator_matrices[j];
//           Eigen::MatrixXd w_expansion=prel_operator*lgi_expansion;
//           basis::MatrixVector temp;
//           temp.push_back(prel_operator.transpose()*prel_operator);
//           //temp.push_back(w_expansion.transpose()*w_expansion);
//           ZeroOutMatrix(temp,1e-8);
//           int rows=w_expansion.rows();
//           Eigen::MatrixXd& Kmatrix=k_matrices[omega];
//           // iterating over each lgi in lgi family
//           int u_max=u3_subspace.size();
//           for(int l=0; l<num_lgi; ++l)
//           {
//             for(int up=0; up<u_max; ++up)
//               {
//                 // target index 
//                 int column=u_max*l+up;
//                  sp_irrep_expansion[ket_index][omega_index].block(0,column,rows,1)
//                     +=Kmatrix(up,u)*w_expansion.block(0,l,rows,1);
//                 // std::cout<<"ho"<<std::endl;
//               }
//           }
//           // std::cout<<"end lgi loop"<<std::endl;
//         }
//     }
// }

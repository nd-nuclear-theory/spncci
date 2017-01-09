/****************************************************************
  lsu3shell_rme.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "lsu3shell/lsu3shell_rme.h"


#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
//#include <functional>

#include "cppformat/format.h"

extern double zero_threshold;

namespace lsu3shell
{
  //TODO Need to account for rho multiplicity in bra and ket. 
  void 
  ReadLSU3ShellRMEs(
      std::ifstream& is,
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::SectorsU3SPN& sectors,
      basis::MatrixVector& matrix_vector 
    )
  {    
    int i,j;
    double rme;
    basis::SetOperatorToZero(sectors,matrix_vector);
    std::string line;
    while(std::getline(is,line))
      {
        // skip initial header line
        if(not std::isdigit(line[0]))
          continue;
        // extract bra/ket lsu3shell basis multiplicity group indices
        std::istringstream line_stream(line);
        line_stream >> i >> j;
        // std::cout<<i<<" "<<j<<std::endl;
        // retrieve lsu3shell basis multiplicity group information
        u3shell::U3SPN omegaSPNi, omegaSPNj;
        // std::tie(omegaSPNi,group_size_i,start_index_i)=lsu3_basis_table[i];
        // std::tie(omegaSPNj,group_size_j,start_index_j)=lsu3_basis_table[j];
        const LSU3BasisGroupData& group_i = lsu3_basis_table[i];
        const LSU3BasisGroupData& group_j = lsu3_basis_table[j];

        u3::SU3 xi(group_i.omegaSPN.SU3());
        u3::SU3 xj(group_j.omegaSPN.SU3());
        // std::cout<<fmt::format("{}  {}  {}", group_i.omegaSPN.Str(), operator_labels.Str(),group_j.omegaSPN.Str())<<std::endl;
        int rho0_max=u3::OuterMultiplicity(xj,operator_labels.x0(),xi);
        // std::cout<<group_i.dim<<"  "<<group_j.dim<<"  "<<rho0_max<<std::endl;
        // extract and store matrix elements
        int i_space=space.LookUpSubspaceIndex(group_i.omegaSPN);
        int j_space=space.LookUpSubspaceIndex(group_j.omegaSPN);
        for(int gi=0; gi<group_i.dim; ++gi)
          for(int gj=0; gj<group_j.dim; ++gj)
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                // std::cout<<"getting rme"<<std::endl;
                line_stream >> rme;
                if(fabs(rme)<zero_threshold)
                  continue;
                // std::cout<<fmt::format("{} {}  {} {} {}  {}",i,j,i_space,j_space,rho0,rme)<<std::endl;
                // Note: Since rho0 is most rapidly varying index in sector enumeration, we could just 
                // calculate the sector_index by offsetting from the sector with rho0=1.
                int sector_index=sectors.LookUpSectorIndex(i_space,j_space,rho0);
                assert(sector_index!=-1);
                int row_index=group_i.start_index+gi;
                int column_index=group_j.start_index+gj;
                // std::cout<<fmt::format("sector {} row {} column {} matrix ({},{})  {}",
                //   sector_index, row_index,column_index, matrix_vector[sector_index].rows(),
                //   matrix_vector[sector_index].cols(),rme)<<std::endl;
                // std::cout<<"sector index "<<sector_index<<std::endl;
                matrix_vector[sector_index](row_index,column_index)=rme;
              }
      }
  }

  bool 
  CompareLSU3ShellRMEs(
      std::ostream& log_stream,
      const U3SPNBasisLSU3Labels& basis_provenance,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::SectorsU3SPN& sectors,
      const basis::MatrixVector& matrices1,
      const basis::MatrixVector& matrices2,
      double tolerance,
      bool verbose
    )
  {  

    // initialize statistics
    int entries_compared = 0;
    double max_residual = 0.;
    double total_sqr_residual = 0.;
    bool success = true;

    // iterate over sectors
    for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
      {

        // retrieve sector
        const u3shell::SectorsU3SPN::SectorType& sector = sectors.GetSector(sector_index);
        if (verbose)
          log_stream
            << fmt::format(
                "sector {} bra {} ket {} dim {}x{}",
                sector_index,
                sector.bra_subspace().Str(),
                sector.ket_subspace().Str(),
                sector.bra_subspace().size(),
                sector.ket_subspace().size()
              )
            << std::endl;

        // iterate over matrix elements
        for (int bra_index=0; bra_index<sector.bra_subspace().size(); ++bra_index)
          for (int ket_index=0; ket_index<sector.ket_subspace().size(); ++ket_index)
            {
              // retrieve matrix elements
              double rme1 = matrices1[sector_index](bra_index,ket_index);
              double rme2 = matrices2[sector_index](bra_index,ket_index);

              // compare matrix elements
              double residual = std::fabs(rme2-rme1);
              ++entries_compared;
              max_residual = std::max(max_residual,residual);
              total_sqr_residual += sqr(residual);
              bool entries_agree = (residual <= tolerance);
              // std::cout<<std::endl<<fmt::format("residual {}  tolerance {}  bool {}", residual,tolerance,entries_agree)<<std::endl<<std::endl;
              success &= entries_agree;

              // write entry diagnostics
              if (!entries_agree)
                {
                log_stream
                  << fmt::format(
                      "  FAIL: {}:({},{}) rme1 {:e} rme2 {:e} residual {:e}",
                      sector_index,bra_index,ket_index,
                      rme1,rme2,residual
                    )
                  << std::endl;
                const lsu3shell::LSU3BasisGroupLabels& bra_labels = basis_provenance[sector.bra_subspace_index()][bra_index];
                const lsu3shell::LSU3BasisGroupLabels& ket_labels = basis_provenance[sector.ket_subspace_index()][ket_index];
                log_stream
                  << fmt::format(
                      "  bra ip {} in {} Np {} Nn {} ; ket ip {} in {} Np {} Nn {}",
                      bra_labels.ip,bra_labels.in,bra_labels.Np,bra_labels.Nn,
                      ket_labels.ip,ket_labels.in,ket_labels.Np,ket_labels.Nn
                    )
                  << std::endl;
                }
            }
      }

    // generate global diagnostics
    double rms_residual = std::sqrt(total_sqr_residual/entries_compared);
    log_stream
      << fmt::format(
          "entries {} rms_residual {:e} max_residual {:e} total_sqr_residual {:e}",
          entries_compared,rms_residual,max_residual,total_sqr_residual
        )
      << std::endl<<std::endl;
        

    return success;
  }

  void GenerateNcmMatrixVector(
    int A,      
    std::ifstream& is_nrel,
    const lsu3shell::LSU3BasisTable& lsu3_basis_table,
    const u3shell::SpaceU3SPN& space, 
    basis::MatrixVector& matrix_vector 
  )
  {
    assert(is_nrel.is_open());

    // Read in Nrel matrix elements and populate sectors
    u3shell::OperatorLabelsU3ST nrel_labels(0,u3::SU3(0,0),0,0,0);
    basis::MatrixVector nrel_matrix_vector;
    u3shell::SectorsU3SPN nrel_sectors(space,nrel_labels,true);
    lsu3shell::ReadLSU3ShellRMEs(
        is_nrel,nrel_labels,lsu3_basis_table,space, 
        nrel_sectors,nrel_matrix_vector
      );

    // Resize vector
    matrix_vector.resize(nrel_matrix_vector.size());

    // Iterate over Nrel subspaces and populate Ncm sectors in matrix_vector
    for(int i=0; i<nrel_matrix_vector.size(); ++i)
      {
        auto subspace=space.GetSubspace(i);
        // eigenvalue of N is given by total number of oscilator quanta minus
        // the zero point energy boson 3A/2. 
        HalfInt N=subspace.N()-3.*A/2;
        // std::cout<<N<<std::endl<<nrel_matrix_vector[i]<<std::endl;
        int dim=subspace.size();
        // std::cout<< Eigen::MatrixXd::Identity(dim,dim)*double(N) <<"     "<<nrel_matrix_vector[i]<<std::endl;
        // std::cout<<"nrel"<<std::endl<<nrel_matrix_vector[i]<<std::endl;
    
        // Ncm=N-Nrel
        matrix_vector[i]=Eigen::MatrixXd::Identity(dim,dim)*double(N)-nrel_matrix_vector[i];
      }
  }



////////////////////////////////////////////////////////////////
}// end namespace

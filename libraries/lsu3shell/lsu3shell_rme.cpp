/****************************************************************
  lsu3shell_rme.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "lsu3shell/lsu3shell_rme.h"


#include <fstream>
#include <iostream>
#include <algorithm>
//#include <functional>

#include "cppformat/format.h"

namespace lsu3shell
{
  void 
  ReadLSU3ShellRMEs(
      std::ifstream& is,
      const u3shell::OperatorLabelsU3S& operator_labels,
      const LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::SectorsU3SPN& sectors,
      basis::MatrixVector& matrix_vector 
    )
  {    
    int i,j, start_index_i, start_index_j, group_size_i, group_size_j;
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
        line_stream >> i, j;

        // retrieve lsu3shell basis multiplicity group information
        u3shell::U3SPN omegaSPNi, omegaSPNj;
        // std::tie(omegaSPNi,group_size_i,start_index_i)=lsu3_basis_table[i];
        // std::tie(omegaSPNj,group_size_j,start_index_j)=lsu3_basis_table[j];
        const LSU3BasisGroupData& group_i = lsu3_basis_table[i];
        const LSU3BasisGroupData& group_j = lsu3_basis_table[j];

        u3::SU3 xi(group_i.omegaSPN.SU3());
        u3::SU3 xj(group_j.omegaSPN.SU3());
        int rho0_max=u3::OuterMultiplicity(xj,operator_labels.x0(),xi);

        // extract and store matrix elements
        int i_space=space.LookUpSubspaceIndex(group_i.omegaSPN);
        int j_space=space.LookUpSubspaceIndex(group_j.omegaSPN);
        for(int gi=0; gi<group_size_i; ++gi)
          for(int gj=0; gj<group_size_j; ++gj)
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                line_stream >> rme;

                // Note: Since rho0 is most rapidly varying index in sector enumeration, we could just 
                // calculate the sector_index by offsetting from the sector with rho0=1.
                int sector_index=sectors.LookUpSectorIndex(i_space,j_space,rho0);
                int row_index=group_i.start_index+gi;
                int column_index=group_j.start_index+gj;
                matrix_vector[sector_index](row_index,column_index)=rme;
              }
      }
  }

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
}// end namespace

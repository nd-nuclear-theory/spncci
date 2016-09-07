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
    u3shell::U3SPN omegaSPNi, omegaSPNj;
    int i,j, start_index_i, start_index_j, group_size_i, group_size_j;
    double rme;
    basis::SetOperatorToZero(sectors,matrix_vector);

    while(is)
      {
        is>>i,j;
        std::tie(omegaSPNi,group_size_i,start_index_i)=lsu3_basis_table[i];
        std::tie(omegaSPNj,group_size_j,start_index_j)=lsu3_basis_table[j];
        int i_space=space.LookUpSubspaceIndex(omegaSPNi);
        int j_space=space.LookUpSubspaceIndex(omegaSPNj);
        u3::SU3 xi(omegaSPNi.SU3());
        u3::SU3 xj(omegaSPNj.SU3());
        int rho0_max=u3::OuterMultiplicity(xj,operator_labels.x0(),xi);
        for(int gi=0; gi<group_size_i; ++gi)
          for(int gj=0; gj<group_size_j; ++gj)
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                is>>rme;
                int sector_index=sectors.LookUpSectorIndex(i_space,j_space,rho0);
                int row_index=start_index_i+gi;
                int column_index=start_index_j+gj;
                matrix_vector[sector_index](row_index,column_index)=rme;
              }
      }
  }

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
}// end namespace

/****************************************************************
  su3rme.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "lsu3shell/su3rme.h"


namespace lsu3shell
{
  unsigned int get_num_ncsmsU3xSU3Basis_irreps(const lsu3::CncsmSU3xSU2Basis& basis)
    {
      const uint32_t number_ipin_blocks = basis.NumberOfBlocks();

      unsigned int num_irreps=0; 
      for (unsigned int ipin_block = 0; ipin_block < number_ipin_blocks; ipin_block++)
        {
          int32_t ibegin = basis.blockBegin(ipin_block);
          int32_t iend = basis.blockEnd(ipin_block);
          for (int32_t iwpn = ibegin; iwpn < iend; ++iwpn) 
            {
              num_irreps++;
            }
        }
      
      return num_irreps;
    }
}// end namespace

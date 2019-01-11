/****************************************************************
  make_6j_table_q3.cpp

  Generate table of "222" 6-j coefficients needed for calculation of
  Q-invariant q3 in analysis postprocessor.

  The coefficient is:

    {2 2     2        }
    {J J_bar J_bar_bar} 

  where each pair must be triangular with 2.  This triangularity
  condition also enforces that the J's must either all be integer or
  all be half-integer, as expected for J values drawn from the same
  many-body space.

  Output is only generated for J values in canonical order
  J<=J_bar<=J_bar_bar, since all other orderings give identical values
  for the coefficient by 6-j symmetry.

  Syntax:
    make_6j_table_q3 Jmax

  Example:
    make_6j_table_q3 20 > make_6j_table_q3_Jmax20.dat

  Output format:
    J J_bar J_bar_bar coeff
    ...

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  7/16/17 (mac): Created.
****************************************************************/

#include <iostream>

#include "am/am.h"
#include "am/halfint.h"
#include "am/wigner_gsl.h"

#include "fmt/format.h"


int main(int argc, char **argv)
{

  // process arguments
  if(argc<1+1)
    {
      std::cout << "Syntax: J_max" << std::endl;
      std::exit(EXIT_SUCCESS);
    }
  int Jmax = std::atoi(argv[1]);

  // tabulation loop
  for (HalfInt J=0; J<=Jmax; J+=HalfInt(1,2))
    for (HalfInt J_bar=J; J_bar<=Jmax; J_bar+=1)
      for (HalfInt J_bar_bar=J_bar; J_bar_bar<=Jmax; J_bar_bar+=1)
        {
          // short circuit on triangularity
          if (!(
                  am::AllowedTriangle(2,J,J_bar)
                  &&am::AllowedTriangle(2,J,J_bar_bar)
                  &&am::AllowedTriangle(2,J_bar,J_bar_bar)
                ))
            continue;

          // tabulate
          double coef = am::Wigner6J(2,2,2,J,J_bar,J_bar_bar);
          std::cout
            << fmt::format("{:4.1f} {:4.1f} {:4.1f} {:+11.8f}",float(J),float(J_bar),float(J_bar_bar),coef)
            << std::endl;
        }

}

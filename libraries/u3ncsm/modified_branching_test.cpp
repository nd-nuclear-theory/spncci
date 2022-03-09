/****************************************************************
  dimensions_test.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  11/2/21 (aem): Created.
****************************************************************/
#include <fstream>
#include <ostream>  

#include "fmt/format.h"
#include "lgi/lgi.h"
#include "u3ncsm/dimensions.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "utilities/utilities.h"
#include "utilities/nuclide.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3boson.h"
#include "sp3rlib/vcs.h"


void test_modified_branching_rule(const nuclide::NuclideType& nuclide, int Nmax=20)
  {
    bool intrinsic = true;
    HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide, intrinsic);

    // Generate list of lsu3shell center of mass free U(3)SpSnS irreps
    std::map<u3shell::U3SPN, unsigned int>
    lgi_dimensions = lsu3shell::lsu3shell_cmf_basis_dimensions(nuclide,Nsigma0,Nmax);

    // Loop the the basis.  For each lgi, ladder and remove U(3)SpSnS irreps
    // which are obtained by laddering
    for(const auto& [lgi,multiplicity] : lgi_dimensions)
      {
        if(multiplicity==0)
          continue;

        const auto& sigma = lgi.U3();
        auto Nn_max = int(Nmax-(sigma.N()-Nsigma0));

        // Sp(3,R) irrep must be unitary for branching to work.
        assert(sp3r::IsUnitary(sigma));

        // Construct U3 boson irrep and generate K matrices
        vcs::U3BosonSpace u3boson_space(sigma,Nn_max);
        auto K_matrices = vcs::GetKMatrices(sigma,u3boson_space);

        for(const auto& [omega,KK] : K_matrices)
          {
            if(sigma==omega)
              continue;

            u3shell::U3SPN labels(omega,lgi.Sp(),lgi.Sn(),lgi.S());
            int upsilon_max = KK[0].cols();

            // Debug check
            if(upsilon_max*multiplicity > lgi_dimensions[labels])
              {
                fmt::print("{}  {} {} {}\n",
                  labels.Str(),upsilon_max,multiplicity,lgi_dimensions[labels]
                );
              }

            assert(upsilon_max*multiplicity<=lgi_dimensions[labels]);

            // decrement
            lgi_dimensions[labels]-=upsilon_max*multiplicity;
          }
      }

    std::cout<<"----------------------------------------"<<std::endl;
    for(const auto& [lgi,multiplicity] : lgi_dimensions)
      if(multiplicity>0)
        {
          fmt::print("{}  {}\n",lgi.Str(),multiplicity);
        }

    std::cout<<"----------------------------------------"<<std::endl;
  }

int main(int argc, char **argv)
{
  std::vector<nuclide::NuclideType> nuclide_list = {{1,1},{2,1},{2,2},{2,3}};
  int Nmax=4;
  u3::U3CoefInit(39);

  for(const auto& nuclide : nuclide_list)
    test_modified_branching_rule(nuclide, Nmax);

}//end main

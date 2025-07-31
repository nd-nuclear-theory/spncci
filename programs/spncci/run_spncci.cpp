/****************************************************************
  run_spncci.cpp

  // Working template for final run...putting the pieces together.

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
  
  3/10/22 (aem): Created.
****************************************************************/

#include "fmt/format.h"

#include "lgi/lgi.h"
#include "sp3rlib/u3coef.h"
#include "spncci_basis/recurrence_indexing.h"

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  int Nmax, Z, N;

  nuclide::NuclideType nuclide({Z,N});
  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide);
  int N1v = nuclide::ValenceShellForNuclide(nuclide);

  // Either generate or read from file
  MultiplicityTaggedLGIVector lgi_vector = get_lgi_vector(nuclide, Nsigma0,Nmax);

  // std::cout<<"setting up recurrence indexing "<<std::endl;
  spncci::spin::Space<lgi::LGI> spin_space(lgi_vector, Nmax);
  const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST>
    spin_recurrence_space(spin_space, spin_space);


    // spatial::RecurrenceSpace() []
    // ->spatial::RecurrenceSp3RSpace() [sigma,sigma',parity_bar]
    //   ->spatial::RecurrenceNnsumSpace() [Nsum]
    //     ->spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
    //       ->spatial::RecurrenceOperatorSubspace() [x0] ->rho0_max
    //         ->spatial::RecurrenceOperatorState() [Nbar,Nbar']

  auto it =iter::imap([](MultiplicityTagged<lgi::LGI> l) { return l.irrep.U3(); }, lgi_vector)| iter::unique_everseen;
  const spncci::spatial::Space spatial_space(std::vector<u3::U3>(it.begin(), it.end()), Nsigma0, Nmax);
  const spncci::spatial::RecurrenceSpace spatial_recurrence_space(spatial_space,spatial_space,N1v,Nsigma0);


  // Loop over RecurrenceSp3Spaces, i.e., over sigma,sigma',parity_bar
  for(int i=0; i<spatial_recurrence_space.size(); ++spatial_recurrence_space)
    {
      const auto& recurrence_sp3r_space = spatial_recurrence_space.GetSubspace(i);
      const auto& [sigma_ket,sigma_bra,parity_bar] = recurrence_sp3r_space.labels();

      // Get seeds
      std::string seed_filename = spncci::seeds::seed_filename(Z,N,Nsigma0,sigma_bra,sigma_ket,parity_bar);
      basis::OperatorBlock<double> seeds = utils::ReadOperatorBlockBinary(seed_filename)

      // Construct recurrence tiles
      // Compute rmes
      // Contract
      // Branch

    }

  // Diagonalize Matrix

}

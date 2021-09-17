/****************************************************************
    recurrence_indexing_test.cpp

  Anna E. McCoy
  INT

  SPDX-License-Identifier: MIT

 6/25/21 (aem): Created.
****************************************************************/
#include "spncci/recurrence_indexing.h"

#include <string>

#include "am/halfint_fmt.h"
#include "cppitertools/imap.hpp"
#include "cppitertools/unique_everseen.hpp"
#include "mcutils/profiling.h"
#include "fmt/format.h"
#include "lgi/lgi.h"


int main(int argc, char** argv)
{
  //   spncci::ValenceShellForNuclide({2,1}),

  ////////////////////////////////////////////////////////////////
  // set up SpinSpace
  ////////////////////////////////////////////////////////////////

  // read in LGIs for 6Li
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <lgi file> <Nmax>" << std::endl;
    return 1;
  }
  std::string filename = argv[1];
  // "../../data/lgi_set/lgi_test.dat";  // test file in data/lgi_set/lgi_test.dat
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  HalfInt Nsigma0 = lgi::Nsigma0ForNuclide({3, 3});
  int N1v = 1;
  int Nmax = std::stoi(argv[2]);
  lgi::ReadLGISet(filename, Nsigma0, lgi_vector);

  constexpr bool check_lgi = true;
  constexpr bool check_spin = false;
  constexpr bool check_spin_recurrence = false;
  constexpr bool check_spatial = true;
  constexpr bool check_spatial_recurrence = true;

  if constexpr (check_lgi)
  {
    std::cout << "LGI set" << std::endl;
    for (std::size_t i = 0; i < lgi_vector.size(); ++i)
      std::cout << i << " " << lgi_vector[i].Str() << std::endl;

    std::cout << "********************************" << std::endl;
  }

  // Checking spncci::spin::Space indexing
  if constexpr (check_spin)
  {
    // diagnostic -- inspect LGI listing
    auto spin_space = spncci::spin::Space<lgi::LGI>(lgi_vector, Nmax);

    int num_lgi_spaces = spin_space.size();
    std::cout << fmt::format("spin::Space dimension: {}", spin_space.dimension())
              << std::endl;
    for (std::size_t i = 0; i < num_lgi_spaces; ++i)
    {
      const auto& spin_lgi_space = spin_space.GetSubspace(i);
      int num_subspaces = spin_lgi_space.size();
      const u3::U3& sigma = spin_lgi_space.sigma();
      std::cout << fmt::format(
          "spin::LGISpace : {}    Size: {:4d}   Dimension: {:4d}",
          sigma.Str(),
          num_subspaces,
          spin_lgi_space.dimension()
        ) << std::endl;

      for (std::size_t j = 0; j < num_subspaces; ++j)
      {
        const auto& spin_subspace = spin_lgi_space.GetSubspace(j);
        const HalfInt& S = spin_subspace.S();
        std::size_t num_spin_states = spin_subspace.size();
        std::cout << fmt::format(
            "    Spin Subspace : {}       Size: {:4d}   Dimension: {:4d}",
            S,
            num_spin_states,
            spin_subspace.dimension()
          ) << std::endl;
        for (std::size_t k = 0; k < num_spin_states; ++k)
        {
          auto spin_state = spin_subspace.GetState(k);
          const HalfInt& Sp = spin_state.labels().Sp;
          const HalfInt& Sn = spin_state.labels().Sn;
          const int gamma_max = spin_state.degeneracy();

          // TODO: fix
          int state_index =
              spin_space.GetSubspaceOffset(i)  // offset to given sigma subspace
              + spin_space.GetSubspace(i).GetSubspaceOffset(j)
              + spin_space.GetSubspace(i).GetSubspace(j).GetStateOffset(
                  k, gamma_max
                );

          std::cout << fmt::format(
              "       Spin State : {} {}      gamma: {:4d}", Sp, Sn, gamma_max
            ) << std::endl;
          std::cout << fmt::format("         index: {}", state_index)
                    << std::endl;
        }
      }
    }
  }

  // Checking spncci::spin::RecurrenceSpace indexing
  if constexpr (check_spin_recurrence)
  {
    std::cout << "___________________________________________________"
              << std::endl;

    spncci::spin::Space<lgi::LGI> spin_space(lgi_vector, Nmax);
    std::vector<std::tuple<u3::SU3, int, int>> spatial_unit_tensors;
    spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST>
        recurrence_space(spin_space, spin_space);
    std::cout << fmt::format("spin::Space dimension: {}", spin_space.dimension())
              << std::endl;

    std::size_t total_dimension = 0;
    std::cout << fmt::format(
        "spin::RecurrenceSpace dimension: {}", recurrence_space.dimension()
      ) << std::endl;

    for (std::size_t i = 0; i < recurrence_space.size(); ++i)
    {
      const auto& recurrence_lgi_space = recurrence_space.GetSubspace(i);

      const u3::U3& sigma_ket = recurrence_lgi_space.sigma_ket();
      const u3::U3& sigma_bra = recurrence_lgi_space.sigma_bra();
      std::cout << fmt::format(
          "  spin::RecurrenceLGISpace {} {} {}  size: {:4d}  dimension: {:4d}",
          sigma_ket.Str(),
          sigma_bra.Str(),
          std::get<2>(recurrence_lgi_space.labels()),
          recurrence_lgi_space.size(),
          recurrence_lgi_space.dimension()
        ) << std::endl;

      for (std::size_t j = 0; j < recurrence_lgi_space.size(); ++j)
      {
        const auto& recurrence_spin_space = recurrence_lgi_space.GetSubspace(j);
        const auto& [S_ket, S_bra] = recurrence_spin_space.labels();

        std::cout << fmt::format(
            "    spin::RecurrenceSpinSpace {} {}  size:  {:4d}  dimension:  "
            "{:4d}",
            S_ket,
            S_bra,
            recurrence_spin_space.size(),
            recurrence_spin_space.dimension()
          ) << std::endl;

        for (std::size_t k = 0; k < recurrence_spin_space.size(); ++k)
        {
          const auto& recurrence_spin_subspace =
              recurrence_spin_space.GetSubspace(k);
          const int degeneracy = recurrence_spin_space.GetSubspaceDegeneracy(k);
          const int ket_degeneracy =
              recurrence_spin_space.GetKetSubspaceDegeneracy(k);
          const int bra_degeneracy =
              recurrence_spin_space.GetBraSubspaceDegeneracy(k);

          const auto [Sp_ket, Sn_ket] =
              recurrence_spin_subspace.ket_upstream_labels();
          const auto [Sp_bra, Sn_bra] =
              recurrence_spin_subspace.bra_upstream_labels();
          std::cout << fmt::format(
              "      spin::RecurrenceSpinSubspace ({},{})  ({},{})   "
              "bra_subspace_degeneracy: "
              "{:3d} ket_subspace_degeneracy: {:3d}   degeneracy: {:4d}",
              Sp_ket,
              Sn_ket,
              Sp_bra,
              Sn_bra,
              ket_degeneracy,
              bra_degeneracy,
              degeneracy
            ) << std::endl;
          std::cout << fmt::format(
              "         size: {:3d} dimension: {:3d}   offset: {:4d}",
              recurrence_spin_subspace.size(),
              recurrence_spin_subspace.dimension(),
              recurrence_spin_space.GetSubspaceOffset(k, degeneracy)
            ) << std::endl;
          for (std::size_t l = 0; l < recurrence_spin_subspace.size(); ++l)
          {
            const auto& recurrence_operator_state =
                recurrence_spin_subspace.GetState(l);
            std::cout << fmt::format(
                "        spin::RecurrenceOperatorState {:02x}   "
                "S0: {}  T0:{}  Sbar: {}  Sbarp: {}  Tbar: {}  Tbarp: {}\n",
                recurrence_operator_state.labels().id(),
                recurrence_operator_state.labels().S0(),
                recurrence_operator_state.labels().T0(),
                recurrence_operator_state.labels().Sbar(),
                recurrence_operator_state.labels().Sbarp(),
                recurrence_operator_state.labels().Tbar(),
                recurrence_operator_state.labels().Tbarp()
              );
          }
          total_dimension += degeneracy * recurrence_spin_subspace.dimension();
        }
      }
    }
    std::cout << fmt::format("total_dimension: {}", total_dimension) << std::endl;
  }

  // Checking spncci::spatial::Space indexing
  if constexpr (check_spatial)
  {
    std::cout << "Nmax = " << Nmax << std::endl;
    auto it =
        iter::imap(
            [](MultiplicityTagged<lgi::LGI> l) { return l.irrep.U3(); }, lgi_vector
          )
        | iter::unique_everseen;
    spncci::spatial::Space spatial_space(
        std::vector<u3::U3>(it.begin(), it.end()), Nsigma0, Nmax
      );
    std::cout << fmt::format(
        "spatial::Space dimension: {}", spatial_space.dimension()
      ) << std::endl;
    for (std::size_t i = 0; i < spatial_space.size(); ++i)
    {
      const auto& spatial_lgi_space = spatial_space.GetSubspace(i);
      const u3::U3& sigma = spatial_lgi_space.sigma();
      int num_spatial_lgi_subspace = spatial_lgi_space.size();
      std::cout << fmt::format(
          "spatial::LGISpace : {}    Size: {:4d}   Dimension: {:4d}",
          sigma.Str(),
          num_spatial_lgi_subspace,
          spatial_lgi_space.dimension()
        ) << std::endl;
      for (std::size_t j = 0; j < num_spatial_lgi_subspace; ++j)
      {
        const auto& u3subspace = spatial_lgi_space.GetSubspace(j);
        const u3::U3& omega = u3subspace.omega();
        const int& upsilon_max = u3subspace.upsilon_max();
        std::cout << fmt::format(
            "    U3Subspace : {}   size: {:4d}   upsilon_max : {:4d}   "
            "dimension: {:4d}",
            omega.Str(),
            u3subspace.size(),
            upsilon_max,
            u3subspace.dimension()
          ) << std::endl;
        for (std::size_t k = 0; k < u3subspace.size(); ++k)
        {
          const auto u3state = u3subspace.GetState(k);
          const u3::U3& n = u3state.n();
          const int& rho_max = u3state.rho_max();
          std::cout << fmt::format("       U3State : {} {} ", n.Str(), rho_max)
                    << std::endl;

          int state_index =
              spatial_space.GetSubspaceOffset(i)
              + spatial_space.GetSubspace(i).GetSubspaceOffset(
                  j
                )  // omega offset upsilon=0
              + k;
        }
      }
    }
  }

  // Checking spncci::spatial::RecurrenceSpace indexing
  if constexpr (check_spatial_recurrence)
  {
    std::cout << "___________________________________________________"
              << std::endl;

    auto it =
        iter::imap(
            [](MultiplicityTagged<lgi::LGI> l) { return l.irrep.U3(); }, lgi_vector
          )
        | iter::unique_everseen;
    const spncci::spatial::Space spatial_space(
        std::vector<u3::U3>(it.begin(), it.end()), Nsigma0, Nmax
      );
    std::vector<std::tuple<u3::SU3, int, int>> spatial_unit_tensors;
    // spncci::spatial::RecurrenceSpace recurrence_space(
    //     spatial_space, spatial_space, N1v, Nsigma0
    //   );
    // std::cout << fmt::format(
    //     "spatial::Space dimension: {}", spatial_space.dimension()
    //   ) << std::endl;

    // std::cout << fmt::format(
    //     "spatial::RecurrenceSpace dimension: {}", recurrence_space.dimension()
    //   ) << std::endl;

    double lgi_size_acc = 0., lgi_size2_acc = 0.;
    std::size_t lgi_size_max = 0;
    double lgi_dim_acc = 0., lgi_dim2_acc = 0.;
    std::size_t lgi_dim_max = 0;

    std::size_t u3_count = 0;
    double u3_size_acc = 0., u3_size2_acc = 0.;
    std::size_t u3_size_max = 0;
    double u3_dim_acc = 0., u3_dim2_acc = 0.;
    std::size_t u3_dim_max = 0;

    // for (std::size_t i = 0; i < recurrence_space.size(); ++i)
    // {
    //   const auto& recurrence_lgi_space = recurrence_space.GetSubspace(i);
    std::size_t recurrence_space_size = 0;

    #pragma omp parallel for collapse(3) default(none) schedule(dynamic) \
            shared(spatial_space, Nsigma0, N1v, std::cout) \
            reduction(+: recurrence_space_size, u3_count, \
                lgi_size_acc, lgi_size2_acc, lgi_dim_acc, lgi_dim2_acc, \
                u3_size_acc, u3_size2_acc, u3_dim_acc, u3_dim2_acc \
              ) \
            reduction(max: lgi_size_max, lgi_dim_max, u3_size_max, u3_dim_max)
    for (const auto parity_bar : {0, 1})
    for (std::size_t lgi_space_ket_index = 0; lgi_space_ket_index < spatial_space.size(); ++ lgi_space_ket_index)
    for (std::size_t lgi_space_bra_index = 0; lgi_space_bra_index < spatial_space.size(); ++ lgi_space_bra_index)
    {
      std::ostringstream buffer{};
      const auto& lgi_space_ket = spatial_space.GetSubspace(lgi_space_ket_index);
      const auto& lgi_space_bra = spatial_space.GetSubspace(lgi_space_bra_index);
      mcutils::SteadyTimer timer{};
      timer.Start();
      const spncci::spatial::RecurrenceLGISpace recurrence_lgi_space(lgi_space_ket, lgi_space_bra, spncci::spatial::UnitTensorConstraintParameters(N1v, Nsigma0, parity_bar));
      timer.Stop();
      if (recurrence_lgi_space.dimension() == 0) continue;
      ++recurrence_space_size;

      // collect statistics
      const auto lgi_size = recurrence_lgi_space.size();
      lgi_size_acc += static_cast<double>(lgi_size);
      lgi_size2_acc += static_cast<double>(lgi_size) * lgi_size;
      lgi_size_max = std::max(lgi_size_max, lgi_size);
      const auto lgi_dim = recurrence_lgi_space.dimension();
      lgi_dim_acc += static_cast<double>(lgi_dim);
      lgi_dim2_acc += static_cast<double>(lgi_dim) * lgi_dim;
      lgi_dim_max = std::max(lgi_dim_max, lgi_dim);

      const u3::U3& sigma_ket = recurrence_lgi_space.sigma_ket();
      const u3::U3& sigma_bra = recurrence_lgi_space.sigma_bra();
      const uint8_t& parity_bar = recurrence_lgi_space.parity_bar();
      buffer << fmt::format(
          "  spatial::RecurrenceLGISpace {:>10s} {:>10s} {}  size: {:4d}  "
          "dimension: {:8d}  recurrence size:  {:10d}  time: {}",
          sigma_ket.Str(),
          sigma_bra.Str(),
          (parity_bar == 0 ? '+' : '-'),
          recurrence_lgi_space.size(),
          recurrence_lgi_space.dimension(),
          (recurrence_lgi_space.dimension()
           - recurrence_lgi_space.GetSubspace(0).dimension())
              * recurrence_lgi_space.GetSubspace(0).dimension(),
          timer.ElapsedTime()
        ) << std::endl;

      for (std::size_t j = 0; j < recurrence_lgi_space.size(); ++j)
      {
        const auto& recurrence_Nnsum_space = recurrence_lgi_space.GetSubspace(j);
        const int& Nnsum = recurrence_Nnsum_space.Nnsum();

        buffer << fmt::format(
            "    spatial::RecurrenceNnsumSpace {:2d}   size:  {:4d}  "
            "dimension:  "
            "{:4d}  final block size:  {:6d}",
            Nnsum,
            recurrence_Nnsum_space.size(),
            recurrence_Nnsum_space.dimension(),
            recurrence_Nnsum_space.dimension()
                * recurrence_lgi_space.LookUpSubspace({0}).dimension()
          ) << std::endl;

        for (std::size_t k = 0; k < recurrence_Nnsum_space.size(); ++k)
        {
          const auto& recurrence_u3space = recurrence_Nnsum_space.GetSubspace(k);

          // collect statistics
          u3_count += 1;
          const auto u3_size = recurrence_u3space.size();
          u3_size_acc += static_cast<double>(u3_size);
          u3_size2_acc += static_cast<double>(u3_size) * u3_size;
          u3_size_max = std::max(u3_size_max, u3_size);
          const auto u3_dim = recurrence_u3space.dimension();
          u3_dim_acc += static_cast<double>(u3_dim);
          u3_dim2_acc += static_cast<double>(u3_dim) * u3_dim;
          u3_dim_max = std::max(u3_dim_max, u3_dim);
          continue;

          const u3::U3& omega_ket = recurrence_u3space.omega_ket();
          const u3::U3& omega_bra = recurrence_u3space.omega_bra();
          const int& upsilon_max_ket = recurrence_Nnsum_space.upsilon_max_ket(k);
          const int& upsilon_max_bra = recurrence_Nnsum_space.upsilon_max_bra(k);
          const int& degeneracy = recurrence_Nnsum_space.GetSubspaceDegeneracy(k);
          buffer << fmt::format(
              "      spatial::RecurrenceU3Space {:>10s} {:>10s}   size:  {:2d} "
              " "
              "dimension:  {:3d}  offset: {:5d}   upsilon_max_ket: {:2d} "
              "upsilon_max_bra: {:2d}   degeneracy: {:4d}",
              omega_ket.Str(),
              omega_bra.Str(),
              recurrence_u3space.size(),
              recurrence_u3space.dimension(),
              recurrence_Nnsum_space.GetSubspaceOffset(k, 1),
              upsilon_max_ket,
              upsilon_max_bra,
              degeneracy
            ) << std::endl;
          // for (std::size_t l = 0; l < recurrence_u3space.size(); ++l)
          // {
          //   const auto& operator_subspace = recurrence_u3space.GetSubspace(l);
          //   buffer << fmt::format(
          //       "        spatial::RecurrenceOperatorSubspace {}   degeneracy: "
          //       "{:4d}",
          //       operator_subspace.x0().Str(),
          //       recurrence_u3space.GetSubspaceDegeneracy(l)
          //     ) << std::endl;
          //   for (std::size_t m = 0; m < operator_subspace.size(); ++m)
          //   {
          //     const auto operator_state = operator_subspace.GetState(m);
          //     buffer << fmt::format(
          //         "          spatial::RecurrenceOperatorState Nbar: {} Nbarp: "
          //         "{}",
          //         operator_state.Nbar(),
          //         operator_state.Nbar_p()
          //       ) << std::endl;
          //   }
          // }
        }
      }
      #pragma omp critical
      std::cout << buffer.str() << std::flush;
    }
    const auto size = recurrence_space_size;
    std::cout << fmt::format(
        "spatial::RecurrenceLGISpace size: {} (avg) {} (std) {} (max)",
        lgi_size_acc / size,
        std::sqrt(
            lgi_size2_acc / size - lgi_size_acc * lgi_size_acc / (size * size)
          ),
        lgi_size_max
      ) << std::endl;
    std::cout << fmt::format(
        "spatial::RecurrenceLGISpace dimension: {} (avg) {} (std) {} (max)",
        lgi_dim_acc / size,
        std::sqrt(lgi_dim2_acc / size - lgi_dim_acc * lgi_dim_acc / (size * size)),
        lgi_dim_max
      ) << std::endl;
    std::cout << fmt::format(
        "spatial::RecurrenceU3Space size: {} (avg) {} (std) {} (max)",
        u3_size_acc / u3_count,
        std::sqrt(
            u3_size2_acc / u3_count - u3_size_acc * u3_size_acc / (u3_count * u3_count)
          ),
        u3_size_max
      ) << std::endl;
    std::cout << fmt::format(
        "spatial::RecurrenceU3Space dimension: {} (avg) {} (std) {} (max)",
        u3_dim_acc / u3_count,
        std::sqrt(u3_dim2_acc / u3_count - u3_dim_acc * u3_dim_acc / (u3_count * u3_count)),
        u3_dim_max
      ) << std::endl;
  }
}

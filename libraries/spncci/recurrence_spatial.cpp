/****************************************************************
  recurrence_spatial.cpp

  Patrick J. Fasano
  University of Notre Dame and LBNL

  SPDX-License-Identifier: MIT
****************************************************************/

#include "spncci/recurrence_spatial.h"

#include <Eigen/Dense>
#include <cppitertools/enumerate.hpp>
#include <utility>
#include <vector>

#include "sp3rlib/u3coef.h"
#include "spncci_basis/recurrence_indexing.h"
#include "spncci/vcs_cache.h"
#include "u3shell/tensor_labels.h"

namespace spncci::recurrence
{

SpatialRecurrenceMatrix::SpatialRecurrenceMatrix(
    std::shared_ptr<const RecurrenceSp3RSpace> space_ptr
  )
    : recurrence_space_ptr_{space_ptr},
      recurrence_blocks_{recurrence_space().size()},
      lgi_recurrence_dimension_{recurrence_space().GetSubspace(/*Nnsum=*/0).dimension()}
{
#ifndef NDEBUG
  recurrence_done_ = std::vector<bool>(recurrence_space().size());
#endif

  for (const auto& Nnsum_space : recurrence_space())
  {
    auto& recurrence_block = GetRecurrenceBlock(Nnsum_space.Nnsum());
    recurrence_block.resize(Nnsum_space.size());
    for (const auto& [i, u3_space] : iter::enumerate(Nnsum_space))
    {
      recurrence_block.at(i) = basis::OperatorBlock<double>::Zero(
          u3_space.dimension(), lgi_recurrence_dimension_
        );
    }
#ifndef NDEBUG
    recurrence_done_.at(Nnsum_space.Nnsum() / 2) = false;
#endif
  }
}

// Functions only defined for this file.  Avoids possible overlaps in names
// in other libraries.
inline namespace
{

basis::OperatorBlock<double> aDaggerBosonMatrix(
    const u3::U3& sigma,
    const u3boson::U3Subspace& u3_subspace_bra,
    const u3boson::U3Subspace& u3_subspace_ket
  )
{
  const auto& bra_dimension = u3_subspace_bra.dimension();
  const auto& omega_bra = u3_subspace_bra.omega();
  const auto& ket_dimension = u3_subspace_ket.dimension();
  const auto& omega_ket = u3_subspace_ket.omega();
  basis::OperatorBlock<double> boson_matrix{bra_dimension, ket_dimension};

  for (std::size_t index_bra = 0; index_bra < bra_dimension; ++index_bra)
  {
    for (std::size_t index_ket = 0; index_ket < ket_dimension; ++index_ket)
    {
      const auto& u3_state_bra = u3_subspace_bra.GetState(index_bra);
      const auto& u3_state_ket = u3_subspace_ket.GetState(index_ket);
      const auto& n_bra = u3_state_bra.n();
      const auto& n_ket = u3_state_ket.n();
      for (int rho_bra = 1; rho_bra <= u3_state_bra.rho_max();
           ++rho_bra, ++index_bra)
      {
        for (int rho_ket = 1; rho_ket <= u3_state_ket.rho_max();
             ++rho_ket, ++index_ket)
        {
          boson_matrix(index_bra, index_ket) =
              ParitySign(
                  u3::ConjugationGrade(omega_ket) + u3::ConjugationGrade(omega_bra)
                )
              // clang-format off
              * u3::U(
                  u3::SU3{2, 0}, n_ket.SU3(), omega_bra.SU3(), sigma.SU3(),
                  n_bra.SU3(), 1, rho_bra,
                  omega_ket.SU3(), rho_ket, 1
                )
              // clang-format on
              * u3boson::BosonCreationRME(n_bra, n_ket);
        }
      }
    }
  }
  return boson_matrix;
}

inline basis::OperatorBlock<double> AMatrix(
    const u3::U3& sigma,
    const spncci::spatial::U3Subspace& u3_subspace_bra,
    const spncci::spatial::U3Subspace& u3_subspace_ket
  )
{
  return u3_subspace_bra.K_matrix()
         * aDaggerBosonMatrix(
              sigma,
              u3_subspace_bra.nonorthogonal_basis(),
              u3_subspace_ket.nonorthogonal_basis()
            )
         * u3_subspace_ket.Kinv_matrix();
}

inline basis::OperatorBlock<double> BMatrix(
    const u3::U3& sigma,
    const spncci::spatial::U3Subspace& u3_subspace_bra,
    const spncci::spatial::U3Subspace& u3_subspace_ket,
    const vcs::MatrixCache& K_matrix_cache,
    const vcs::MatrixCache& Kinv_matrix_cache
  )
{
  return ParitySign(
             u3::ConjugationGrade(u3_subspace_ket.omega())
             - u3::ConjugationGrade(u3_subspace_bra.omega())
           )
         * std::sqrt(
             u3::dim(u3_subspace_ket.omega()) / u3::dim(u3_subspace_bra.omega())
           )
         * (u3_subspace_ket.K_matrix()
            * aDaggerBosonMatrix(
                sigma,
                u3_subspace_ket.nonorthogonal_basis(),
                u3_subspace_bra.nonorthogonal_basis()
              )
            * u3_subspace_bra.Kinv_matrix())
               .transpose();
}

inline basis::OperatorBlock<double> ChiMatrix(
    const u3::U3& sigma,
    const spncci::spatial::U3Subspace& u3_subspace_bra,
    const spncci::spatial::U3Subspace& u3_subspace_ket
  )
{
  return ParitySign(
             u3::ConjugationGrade(u3_subspace_ket.omega())
             + u3::ConjugationGrade(u3_subspace_bra.omega())
           )
         * (2. / float(u3_subspace_bra.omega().N() - sigma.N()))
         * u3_subspace_bra.Kinv_matrix().transpose()
         * aDaggerBosonMatrix(
             sigma,
             u3_subspace_bra.nonorthogonal_basis(),
             u3_subspace_ket.nonorthogonal_basis()
           )
         * u3_subspace_ket.K_matrix().transpose();
}

};  // namespace


basis::OperatorBlock<double> UMatrix1(
    const u3::U3& target_omega_bra,
    const u3::U3& target_omega_ket,
    const u3::U3& source_omega_bra,
    const u3::U3& source_omega_ket,
    const SpatialRecurrenceMatrix::RecurrenceU3Space& target_U3_space,
    const SpatialRecurrenceMatrix::RecurrenceU3Space& source_U3_space
  )
{
  basis::OperatorBlock<double> u_matrix = basis::OperatorBlock<double>::Zero(
      target_U3_space.dimension(), source_U3_space.dimension()
    );

  for (auto&& [target_x0_index, target_x0_subspace] :
       iter::enumerate(target_U3_space))
  {
    auto target_r0_max = target_U3_space.GetSubspaceDegeneracy(target_x0_index);
    auto source_x0_index =
        source_U3_space.LookUpSubspaceIndex(target_x0_subspace.labels());
    const auto& source_x0_subspace = source_U3_space.GetSubspace(source_x0_index);
    auto source_r0_max = source_U3_space.GetSubspaceDegeneracy(source_x0_index);
    for (std::size_t target_r0 = 1; target_r0 <= target_r0_max; ++target_r0)
      for (std::size_t source_r0 = 1; source_r0 <= source_r0_max; ++source_r0)
      {
        auto target_offset =
            target_U3_space.GetSubspaceOffset(target_x0_index, target_r0);
        auto source_offset =
            source_U3_space.GetSubspaceOffset(source_x0_index, source_r0);
        for (auto&& [target_op_index, target_op] :
             iter::enumerate(target_x0_subspace))
        {
          auto source_op_index =
              source_x0_subspace.LookUpStateIndex(target_op.labels());
          u_matrix(
              target_offset + target_op_index, source_offset + source_op_index
            ) =
              u3::U(
                  u3::SU3{2u, 0u},
                  source_omega_ket,
                  target_omega_bra,
                  target_x0_subspace.x0(),
                  target_omega_ket,
                  1,
                  target_r0,
                  source_x0_subspace.x0(),// Querry PJF: should this be source_ket?
                  source_r0,
                  1
                );
        }
      }
  }

  return u_matrix;
}


basis::OperatorBlock<double> UMatrix2(
    const u3::U3& target_omega_bra,
    const u3::U3& target_omega_ket,
    const u3::U3& source_omega_bra,
    const u3::U3& source_omega_ket,
    const SpatialRecurrenceMatrix::RecurrenceU3Space& target_U3_space,
    const SpatialRecurrenceMatrix::RecurrenceU3Space& source_U3_space
  )
{
  basis::OperatorBlock<double> u_matrix = basis::OperatorBlock<double>::Zero(
      target_U3_space.dimension(), source_U3_space.dimension()
    );

  for (auto&& [target_x0_index, target_x0_subspace] :
       iter::enumerate(target_U3_space))
    for (auto&& [source_x0_index, source_x0_subspace] :
         iter::enumerate(source_U3_space))
    {
      if (!u3::OuterMultiplicity(
              {2u, 0u}, target_x0_subspace.x0(), source_x0_subspace.x0()
            ))
        continue;
      auto target_r0_max = target_U3_space.GetSubspaceDegeneracy(target_x0_index);
      auto source_r0_max = source_U3_space.GetSubspaceDegeneracy(source_x0_index);
      for (std::size_t target_r0 = 1; target_r0 <= target_r0_max; ++target_r0)
        for (std::size_t source_r0 = 1; source_r0 <= source_r0_max; ++source_r0)
        {
          auto target_offset =
              target_U3_space.GetSubspaceOffset(target_x0_index, target_r0);
          auto source_offset =
              source_U3_space.GetSubspaceOffset(source_x0_index, source_r0);
          for (auto&& [target_op_index, target_op] :
               iter::enumerate(target_x0_subspace))
          {
            const auto& l = target_op.labels();
            const auto& Nbar = target_op.Nbar();
            const auto& Nbarp = target_op.Nbarp();
            auto source_op_index =
                source_x0_subspace.LookUpStateIndex({Nbar - 2});
            u_matrix(
                target_offset + target_op_index, source_offset + source_op_index
              ) =
                u3::U(
                    // clang-format off
                    source_omega_ket, u3::SU3{2u, 0u}, target_omega_bra, target_x0_subspace.x0(),
                    target_omega_ket, 1, target_r0,
                    source_x0_subspace.x0(), source_r0, 1
                    // clang-format on
                  )
                * u3::U(
                    // clang-format off
                    u3::SU3{Nbarp, 0u}, u3::SU3{0u, Nbar}, source_x0_subspace.x0(), u3::SU3{2u, 0u},
                    target_x0_subspace.x0(), 1, 1,
                    u3::SU3{0u, Nbar}, 1, 1
                    // clang-format on
                  );
          }
        }
    }

  return u_matrix;
}

basis::OperatorBlocks<double> UNbarMatrix(
    const u3::U3& target_omega_bra,
    const u3::U3& target_omega_ket,
    const u3::U3& source_omega_bra,
    const u3::U3& source_omega_ket,
    const spncci::spatial::RecurrenceOperatorSectors<SpatialRecurrenceMatrix::OperatorStateLabelType>&
      recurrence_operator_sectors
  )
{
  basis::OperatorBlocks<double>UNbar_matrices;
  basis::SetOperatorToZero(recurrence_operator_sectors,UNbar_matrices);

  for(const auto&& [sector_index,sector] : iter::enumerate(recurrence_operator_sectors))
  {
    basis::OperatorBlock<double> block = UNbar_matrices[sector_index];
    auto target_r0_max = sector.bra_subspace_degeneracy();
    auto source_r0_max = sector.ket_subspace_degeneracy();
    const auto& source_subspace = sector.source_subspace();
    const auto& target_subspace = sector.target_subspace();

    for(int source_r0=1; source_r0<=source_r0_max; ++source_r0)
    {
      for(int target_r0=1; target_r0<=target_r0_max; ++target_r0)
      {
        auto ucoef = u3::U(
              u3::SU3{2u, 0u},
              source_omega_ket,
              target_omega_bra,
              target_subspace.x0(),
              target_omega_ket,
              1,
              target_r0,
              source_omega_ket,
              source_r0,
              1
            );

        for(std::size_t source_state_index=0; source_state_index<source_subspace.size(); source_state_index++)
          for(std::size_t target_state_index=0; target_state_index<target_subspace.size(); target_state_index++)
          {
            auto target_offset = (target_r0-1)*target_subspace.size()+target_state_index;
            auto source_offset = (source_r0-1)*source_subspace.size()+source_state_index;
            const auto target_Nbar = target_subspace.GetState(target_state_index).Nbar();
            const auto source_Nbar = source_subspace.GetState(source_state_index).Nbar();
            if(target_Nbar==source_Nbar)
            {
              block(target_offset,source_offset)=ucoef;
            }
          }
        }

      }
    }
}

basis::OperatorBlock<double> GetKronProd(const basis::OperatorBlock<double>& Matrix1, const basis::OperatorBlock<double>& Matrix2)
{
    // Kronecker product of ChiMatrix x AMatrix
  auto rows1=Matrix1.rows();
  auto cols1=Matrix1.cols();
  auto rows2=Matrix2.rows();
  auto cols2=Matrix2.cols();
  basis::OperatorBlock<double> ProductMatrix{rows1+rows2,cols1+cols2};
  for(auto i1=0; i1<rows1; i1++)
    for(auto j1=0; j1<cols1; j1++)
    {
      ProductMatrix.block(i1*rows2,j1*cols2,rows2,cols2)=Matrix1(i1,j1)*Matrix2;
    }
  return ProductMatrix;
}



void ComputeNnSum4RecurrenceTerm(
      const spncci::spatial::Sp3RSpace& sp3r_space_bra,
      const spncci::spatial::Sp3RSpace& sp3r_space_ket,
      const spncci::spatial::RecurrenceU3Space<SpatialRecurrenceMatrix::OperatorStateLabelType>&
        target_recurrence_u3_space,
      const spncci::spatial::RecurrenceU3Space<SpatialRecurrenceMatrix::OperatorStateLabelType>&
        source_recurrence_u3_space,
      const spncci::spatial::RecurrenceOperatorSectors<SpatialRecurrenceMatrix::OperatorStateLabelType>&
        recurrence_operator_sectors,
      const basis::OperatorBlock<double>& source_recurrence_u3_tile,
      basis::OperatorBlock<double>& target_recurrence_u3_tile
  )
{
      const auto& sigmap = sp3r_space_bra.sigma();
      const auto& sigma = sp3r_space_ket.sigma();
      const auto& omega1 = source_recurrence_u3_space.omega_ket();
      const auto& omega2 = source_recurrence_u3_space.omega_bra();
      const auto& omega =  target_recurrence_u3_space.omega_ket();
      const auto& omegap = target_recurrence_u3_space.omega_bra();


  // Case 1: Nn>Nn'
  if((omega.N()-sigma.N())>=(omegap.N()-sigmap.N()))
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // If Nn >= Nn'
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// This is where recurrence terms change depending on Nnsum-4 or Nnsum-2
    basis::OperatorBlock<double>
      chi_matrix =
        ChiMatrix(
            sigma,
            sp3r_space_bra.LookUpSubspace(omega),
            sp3r_space_ket.LookUpSubspace(omega1)
          );

    basis::OperatorBlock<double>
      A_matrix =
        AMatrix(
            sigmap,
            sp3r_space_bra.LookUpSubspace(omegap),
            sp3r_space_ket.LookUpSubspace(omega2)

          );

    // Kronecker product of ChiMatrix x AMatrix
    basis::OperatorBlock<double> chiA_matrix = GetKronProd(chi_matrix,A_matrix);


    // Vector of matrice indexed by operator sector index (correpsonds to omega0).  Each operator block is indexed by
    // (columns) r0_source, Nbar_source,  and (rows) r0_target, Nbar_target.
    const auto& UNbar_matrices
      = UNbarMatrix(
        omegap,
        omega,
        omega2,
        omega1,
        recurrence_operator_sectors
      );

    // Appy delta_Nsum=4 reccurence terms to Nsum-4 recurrence matrix
    for(auto target_upsilon_index=0; target_upsilon_index<chiA_matrix.rows(); target_upsilon_index++)
      for(auto source_upsilon_index=0; source_upsilon_index<chiA_matrix.rows(); source_upsilon_index++)
        for(std::size_t i=0; i<recurrence_operator_sectors.size(); ++i)
        {
          const auto& sector = recurrence_operator_sectors.GetSector(i);

          std::size_t target_operator_offset=target_recurrence_u3_space.GetSubspaceOffset(sector.target_subspace_index(),target_upsilon_index);
          std::size_t source_operator_offset=source_recurrence_u3_space.GetSubspaceOffset(sector.source_subspace_index(),source_upsilon_index);

          const basis::OperatorBlock<double>& UNbar_matrix = UNbar_matrices[i];
          target_recurrence_u3_tile.block(target_operator_offset,0,UNbar_matrix.rows(),target_recurrence_u3_tile.cols())
            += chiA_matrix(target_upsilon_index,source_upsilon_index)
                *UNbar_matrix
                *source_recurrence_u3_tile.block(source_upsilon_index,0,UNbar_matrix.cols(),source_recurrence_u3_tile.cols());

        }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  }
  // Term 2 Nn'>Nn
  else
  {

  }
}


void SpatialRecurrenceMatrix::GenerateRecurrenceBlock(unsigned int Nnsum)
{
  auto& recurrence_block = GetRecurrenceBlock(Nnsum);
  const auto& lgi_recurrence_subspace =
      recurrence_space().GetSubspace(/*Nnsum=*/0).GetSubspace(0);
  const auto& bra_lgi = lgi_recurrence_subspace.omega_bra();
  const auto& ket_lgi = lgi_recurrence_subspace.omega_ket();
  const auto& recurrence_Nnsum_space = recurrence_space().GetSubspace(Nnsum / 2);
  // initialize Nnsum=0 block to identity
  if (Nnsum == 0)
  {
    recurrence_block[0] = basis::OperatorBlock<double>::Identity(
        lgi_recurrence_dimension_, lgi_recurrence_dimension_
      );
    return;
  }

  if (Nnsum >= 4)
  {
    auto recurrence_u3_sectors =
        spncci::spatial::RecurrenceU3Sectors(recurrence_space(), Nnsum, Nnsum - 4);
    for (const auto& recurrence_u3_sector : recurrence_u3_sectors)
    {
      const auto& source_u3_space = recurrence_u3_sector.source_subspace();
      const auto source_u3_offset =
          recurrence_space()
              .GetSubspace(Nnsum / 2 - 2)
              .GetSubspaceOffset(recurrence_u3_sector.source_subspace_index());
      const auto& target_u3_space = recurrence_u3_sector.target_subspace();
      const auto target_u3_offset =
          recurrence_space().GetSubspace(Nnsum).GetSubspaceOffset(
              recurrence_u3_sector.target_subspace_index()
            );

      const auto& omega1 = source_u3_space.omega_ket();
      const auto& omega2 = source_u3_space.omega_bra();
      const auto& omega = target_u3_space.omega_ket();
      const auto& omegap = target_u3_space.omega_bra();

      const auto& chi_matrix =
          ChiMatrix(
              ket_lgi,
              recurrence_space().ket_space().LookUpSubspace(omega),
              recurrence_space().ket_space().LookUpSubspace(omega1)
            )
              .eval();
      const auto& A_matrix =
          AMatrix(
              ket_lgi,
              recurrence_space().ket_space().LookUpSubspace(omega),
              recurrence_space().ket_space().LookUpSubspace(omega1)

            )
              .eval();
      const auto& u1_matrix = UMatrix1(
          omegap, omega, omega2, omega1, target_u3_space, source_u3_space
        );

      assert(
          chi_matrix.rows() * A_matrix.rows() * u1_matrix.rows()
          == target_u3_space.dimension()
        );
      assert(
          chi_matrix.cols() * A_matrix.cols() * u1_matrix.cols()
          == source_u3_space.dimension()
        );

      basis::OperatorBlock<double> recurrence_tile{
          target_u3_space.dimension(), source_u3_space.dimension()
        };
    }
  }

  if (Nnsum >= 2)
  {
    auto recurrence_u3_sectors =
        spncci::spatial::RecurrenceU3Sectors(recurrence_space(), Nnsum, Nnsum - 2);
    for (const auto& [i, recurrence_u3_sector] :
         iter::enumerate(recurrence_u3_sectors))
    {}
  }

#ifndef NDEBUG
  recurrence_done_[Nnsum / 2] = true;
#endif
}

};  // namespace spncci::recurrence

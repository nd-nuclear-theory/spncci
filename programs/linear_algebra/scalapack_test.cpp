/****************************************************************
  scalapack_test.cpp

  Test matrix distribution and diagonalization with SCALAPACK
  and MPI.

  Input file:
  Pr Pc block_size dimension
  i j rows cols

  Anna E. McCoy and Patrick J. Fasano
  University of Notre Dame

  01/11/19 (pjf): Create.
****************************************************************/

#include <mpi.h>
#include <mkl.h>
#include <mkl_scalapack.h>
#include <mkl_blacs.h>
#include <mkl_pblas.h>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "eigen3/Eigen/Dense"
#include "fmt/format.h"


// constexpr int k_done = 0xbaff1e;
constexpr int k_done=999;

std::pair<int,int> calculateProcessIndices(
    int global_i, int global_j,
    int process_rows, int process_cols,
    int row_block_size, int col_block_size
  )
{
  return {(global_i/row_block_size)%process_rows, (global_j/col_block_size)%process_cols};
}

int getRankForIndices(
    int global_i, int global_j,
    int process_rows, int process_cols,
    int row_block_size, int col_block_size
  )
{
  int pr, pc;
  std::tie(pr, pc) = calculateProcessIndices(
      global_i, global_j,
      process_rows, process_cols,
      row_block_size, col_block_size
    );
  return (pr*process_cols + pc);
}

std::pair<int,int> getProcessIndicesForRank(
    int rank,
    int process_rows, int process_cols
  )
{
  return {(rank-rank%process_cols)/process_rows, rank%process_cols};
}

std::pair<int,int> calculateProcessLocalBlockIndices(
    int global_i, int global_j,
    int process_rows, int process_cols,
    int row_block_size, int col_block_size
  )
{
  return {global_i/(row_block_size*process_rows), global_j/(col_block_size*process_cols)};
}

std::pair<int,int> calculateBlockLocalElementIndices(
    int global_i, int global_j,
    int process_rows, int process_cols,
    int row_block_size, int col_block_size
  )
{
  return {global_i%row_block_size, global_j%col_block_size};
}

std::pair<int,int> calculateProcessLocalElementIndices(
    int global_i, int global_j,
    int process_rows, int process_cols,
    int row_block_size, int col_block_size
  )
{
  int l, m, x, y;
  std::tie(l,m) = calculateProcessLocalBlockIndices(
      global_i, global_j,
      process_rows, process_cols,
      row_block_size, col_block_size
    );
  std::tie(x,y) = calculateBlockLocalElementIndices(
      global_i, global_j,
      process_rows, process_cols,
      row_block_size, col_block_size
    );
  return {l*row_block_size + x, m*col_block_size + y};
}

void constructMatrix(int dimension, Eigen::MatrixXf& matrix)
{
  matrix = Eigen::MatrixXf::Zero(dimension, dimension);
  for (int i=0; i<dimension; ++i)
    for(int j=0; j<dimension; ++j)
    {
      // matrix elements of (b^\dagger + b)^2
      float val=0;
      if (i==j+2)
        val = std::sqrt((j+1)*(j+2));
      else if (i==j)
        val = (2*j+1);
      else if (i==j-2)
        val = std::sqrt(j*(j-1));
      matrix(i,j) = val;
    }

}

int main(int argc, char *argv[])
{
  // initialize MPI and OpenMP
  int num_mpi_ranks, mpi_rank;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_mpi_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Get_processor_name(processor_name, &name_len);

  std::ifstream input_stream("scalapack_test.in");
  int Pr, Pc, block_size, dimension;
  input_stream >> Pr >> Pc >> block_size >> dimension;
  int start_row, start_col, num_row, num_col;
  input_stream >> start_row >> start_col >> num_row >> num_col;
  assert(num_mpi_ranks == Pr*Pc+1);
  assert(start_row<dimension);
  assert(start_col<dimension);
  assert(start_row+num_row<=dimension);
  assert(start_col+num_col<=dimension);

  // set up BLACS -- process grid, etc.
  int my_pc, my_pr, dummy;
  MKL_INT zero = 0, one = 1;
  MKL_INT ictxt=MPI_COMM_WORLD, info, negone=-1, ten=10;
  char order = 'R';  // row-major process assignment
  MKL_INT mypnum, nprocs;
  blacs_pinfo(&mypnum, &nprocs);
  blacs_get(&ictxt, &zero, &ictxt);
  std::cout << "y" << mpi_rank << std::endl << std::flush;
  blacs_gridinit(&ictxt, &order, &Pr, &Pc);
  blacs_gridinfo(&ictxt, &dummy, &dummy, &my_pr, &my_pc);
  // don't use Pr and Pc from gridinfo, but from earlier input

  std::cout << fmt::format(
                  "{:d} (pr,pc)=({:d},{:d})",
                  mpi_rank, my_pr, my_pc
                )
            << std::endl << std::flush;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Comm role_comm;
  MPI_Comm_split(MPI_COMM_WORLD, mpi_rank==(num_mpi_ranks-1), mpi_rank, &role_comm);

  if (mpi_rank == (num_mpi_ranks-1))
  {
    // construct and distribute matrix from master (last MPI rank)
    Eigen::MatrixXf matrix;
    constructMatrix(dimension, matrix);
    std::cout << matrix << std::endl;

    int end_row = start_row + num_row - 1;  // last row to be sent
    int end_col = start_col + num_col - 1;  // last col to be sent
    int size_rows=0, size_cols=0;

    // loop over blocks of the matrix
    for (int start_i = start_row; start_i <= end_row; start_i += size_rows)
      for (int start_j = start_col; start_j <= end_col; start_j += size_cols)
      {
        int target_rank = getRankForIndices(start_i, start_j, Pr, Pc, block_size, block_size);
        std::cout << "send (" << start_i << "," << start_j << ")->" << target_rank << std::endl;
        // we send from the start index to the edge of the block, or to the
        // end of the submatrix
        size_rows = std::min(block_size-(start_i%block_size), end_row-start_i+1);
        size_cols = std::min(block_size-(start_j%block_size), end_col-start_j+1);

        // determine size of packed data
        int int_pack_size, row_pack_size, total_pack_size;
        MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, &int_pack_size);
        MPI_Pack_size(size_rows, MPI_FLOAT, MPI_COMM_WORLD, &row_pack_size);
        total_pack_size = 4*int_pack_size + size_cols*row_pack_size;

        // allocate buffer for pack
        std::unique_ptr<char> buffer(new char[total_pack_size]);

        // pack index information
        int position = 0;
        MPI_Pack(&start_i, 1, MPI_INTEGER, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);
        MPI_Pack(&start_j, 1, MPI_INTEGER, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);
        MPI_Pack(&size_rows, 1, MPI_INTEGER, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);
        MPI_Pack(&size_cols, 1, MPI_INTEGER, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);

        // pack submatrix data
        for (int j = start_j; j < start_j+size_cols; ++j)
        {
          MPI_Pack(&(matrix.data()[dimension*j+start_i]), size_rows, MPI_FLOAT, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);
        }

        // send packed data to destination node
        MPI_Send(buffer.get(), position, MPI_PACKED, target_rank, 0, MPI_COMM_WORLD);
      }

    // inform all workers that matrix distribution is finished
    for (int target_rank=0; target_rank<(num_mpi_ranks-1); ++target_rank)
      MPI_Send(&k_done, 1, MPI_INTEGER, target_rank, k_done, MPI_COMM_WORLD);
  }
  else
  {
    // allocate buffer
    int int_pack_size, row_pack_size, max_pack_size;
    MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, &int_pack_size);
    MPI_Pack_size(block_size, MPI_FLOAT, MPI_COMM_WORLD, &row_pack_size);
    max_pack_size = 4*int_pack_size + block_size*row_pack_size;
    std::unique_ptr<char> buffer(new char[max_pack_size]);

    // determine needed storage size and allocate local matrix
    int local_rows = numroc(&dimension, &block_size, &my_pr, &zero, &Pr);
    int local_cols = numroc(&dimension, &block_size, &my_pc, &zero, &Pc);
    int local_storage_size = local_rows*local_cols;
    Eigen::MatrixXf local_storage = Eigen::MatrixXf::Zero(local_rows, local_cols);

    // receive matrix
    while (true)
    {
      // figure out what type of message this is
      int source, tag;
      std::unique_ptr<MPI_Status> status(new MPI_Status);
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status.get());

      if (status->MPI_TAG == k_done)
      {
        // done receiving, break out
        MPI_Recv(buffer.get(), 1, MPI_INTEGER, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status.get());
        break;
      }
      else
      {
        // receive next chunk
        MPI_Recv(buffer.get(), max_pack_size, MPI_PACKED, status->MPI_SOURCE, status->MPI_TAG, MPI_COMM_WORLD, status.get());
        int start_i, start_j, size_rows, size_cols, position=0;
        MPI_Unpack(buffer.get(), max_pack_size, &position, &start_i, 1, MPI_INTEGER, MPI_COMM_WORLD);
        MPI_Unpack(buffer.get(), max_pack_size, &position, &start_j, 1, MPI_INTEGER, MPI_COMM_WORLD);
        MPI_Unpack(buffer.get(), max_pack_size, &position, &size_rows, 1, MPI_INTEGER, MPI_COMM_WORLD);
        MPI_Unpack(buffer.get(), max_pack_size, &position, &size_cols, 1, MPI_INTEGER, MPI_COMM_WORLD);

        int s, t;
        for (int j = start_j; j < start_j+size_cols; ++j)
        {
          // extract pack in place
          std::tie(s,t) = calculateProcessLocalElementIndices(start_i, j, Pr, Pc, block_size, block_size);
          MPI_Unpack(buffer.get(), max_pack_size, &position, local_storage.data()+(local_rows*t+s), size_rows, MPI_FLOAT, MPI_COMM_WORLD);
        }
      }
    }

    // check matrix distribution
    std::stringstream line_stream;
    line_stream << "local " << mpi_rank << std::endl;
    line_stream << local_storage << std::endl;

    // set up ScaLAPACK environment
    std::unique_ptr<MKL_INT> descA(new MKL_INT[9]);  // descriptor for matrix
    std::unique_ptr<MKL_INT> descZ(new MKL_INT[9]);  // descriptor for eigenvectors
    descinit(descA.get(),
        &dimension, &dimension, &block_size, &block_size,
        &zero, &zero, &ictxt, &local_rows, &info
      );
    descinit(descZ.get(),
        &dimension, &dimension, &block_size, &block_size,
        &zero, &zero, &ictxt, &local_rows, &info
      );

    // construct storage for eigensystem and work space
    char jobz='V', uplo='L';
    MKL_INT lwork;
    Eigen::VectorXf eigval_storage = Eigen::VectorXf::Zero(dimension);
    Eigen::MatrixXf eigvec_storage = Eigen::MatrixXf::Zero(local_rows, local_cols);
    std::vector<float> work(2);

    // dummy call to determine needed size for work
    lwork = -1;
    pssyev(&jobz, &uplo, &dimension,
        local_storage.data(), &one, &one, descA.get(),
        eigval_storage.data(), eigvec_storage.data(), &one, &one, descZ.get(),
        work.data(), &lwork, &info);

    // resize to actual needed work size
    lwork = int(work[0]);
    work.resize(lwork);

    // do actual diagonalization
    pssyev(&jobz, &uplo, &dimension,
        local_storage.data(), &one, &one, descA.get(),
        eigval_storage.data(), eigvec_storage.data(), &one, &one, descZ.get(),
        work.data(), &lwork, &info);

    // output matrix
    line_stream << "eigval " << mpi_rank << std::endl;
    line_stream << eigval_storage.transpose() << std::endl;
    line_stream << "eigvec " << mpi_rank << std::endl;
    line_stream << eigvec_storage << std::endl;

    // dump accumulated diagonistic info
    for (int i=0; i<num_mpi_ranks-1; ++i)
    {
      MPI_Barrier(role_comm);
      if (mpi_rank==i)
        std::cout << line_stream.str() << std::endl;
    }
  }

  // Finalize the MPI environment.
  MPI_Comm_free(&role_comm);
  MPI_Finalize();

  return EXIT_SUCCESS;
}

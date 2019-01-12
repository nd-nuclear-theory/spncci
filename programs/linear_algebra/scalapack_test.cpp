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

#include <cmath>
#include <iostream>
#include <fstream>

#include "eigen3/Eigen/Dense"
#include "fmt/format.h"

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
  return {global_i/(row_block_size*process_rows), global_j/(col_block_size*process_cols)};
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
  assert(start_row+num_row<dimension);
  assert(start_col+num_col<dimension);

  int my_pc, my_pr;
  std::tie(my_pr, my_pc) = getProcessIndicesForRank(mpi_rank, Pr, Pc);

  std::cout << fmt::format(
                  "{:d} (pr,pc)=({:d},{:d})",
                  mpi_rank, my_pr, my_pc
                )
            << std::endl << std::flush;

  if (mpi_rank == num_mpi_ranks)
  {
    Eigen::MatrixXf matrix = Eigen::MatrixXf::Zero(dimension, dimension);
    for (int i=0; i<dimension; ++i)
      for(int j=0; j<dimension; ++j)
      {
        matrix(i,j) = i*dimension + j;
      }

    int end_row = start_row + num_row - 1;
    int end_col = start_col + num_col - 1;
    for (int start_i = start_row; start_i <= end_row; start_i += block_size)
      for (int start_j = start_col; start_j <= end_col; start_j += block_size)
      {
        int target_rank = getRankForIndices(start_i, start_j, Pr, Pc, block_size, block_size);
        int size_rows = std::min(block_size, end_row-start_i+1);
        int size_cols = std::min(block_size, end_col-start_j+1);

        int num_matrix_elements = size_rows*size_cols;
        int int_pack_size, row_pack_size, total_pack_size;
        MPI_Pack_size(1, MPI_INTEGER, MPI_COMM_WORLD, &int_pack_size);
        MPI_Pack_size(size_rows, MPI_FLOAT, MPI_COMM_WORLD, &row_pack_size);
        total_pack_size = 4*int_pack_size + size_cols*row_pack_size;

        std::unique_ptr<char> buffer(new char[total_pack_size]);
        int position;
        MPI_Pack(&start_i, 1, MPI_INTEGER, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);
        MPI_Pack(&start_j, 1, MPI_INTEGER, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);
        MPI_Pack(&size_rows, 1, MPI_INTEGER, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);
        MPI_Pack(&size_cols, 1, MPI_INTEGER, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);
        for (int j = start_j; j < start_j+size_cols; ++j)
          MPI_Pack(&(matrix.data()[j]), size_rows, MPI_FLOAT, buffer.get(), total_pack_size, &position, MPI_COMM_WORLD);

        MPI_Send(buffer.get(), position, MPI_PACKED, target_rank, 0, MPI_COMM_WORLD);
      }
    int k_done = -1;
    MPI_Send(&k_done, 1, MPI_INTEGER, )
  }
  else
  {
    // determine needed storage size
    int zero = 0;
    int local_rows = numroc(&dimension, &block_size, &my_pr, &zero, &Pr);
    int local_cols = numroc(&dimension, &block_size, &my_pc, &zero, &Pc);
    int local_storage_size = local_rows*local_cols;
    std::unique_ptr<float> local_storage(new float[local_storage_size]);


  }

  // Finalize the MPI environment.
  MPI_Finalize();

  return EXIT_SUCCESS;
}

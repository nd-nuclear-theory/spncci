/****************************************************************
  scalapack_test.cpp

  Test matrix distribution and diagonalization with SCALAPACK
  and MPI.

  Input file:
  Pr Pc block_size dimension
  i j rows cols

  Anna E. McCoy and Patrick J. Fasano
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  01/11/19 (pjf): Create.
****************************************************************/

#include <omp.h>
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

struct RunParameters
{
  int Pr;
  int Pc;
  int block_size;
  int dimension;
};

struct LocalParameters
{
  int mpi_rank, num_mpi_ranks;
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  int pr, pc, rows, cols;
  MKL_INT ctxt;
};

struct SubmatrixDescriptor
{
  int i;
  int j;
  int rows;
  int cols;
};

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

void initializeRunParameters(RunParameters& params)
{
  std::ifstream input_stream("scalapack_test.in");
  input_stream >> params.Pr >> params.Pc >> params.block_size >> params.dimension;
}

void initializeParallel(int argc, char *argv[], RunParameters& params, LocalParameters& local)
{
  // initialize MPI and OpenMP
  int name_len;
  int required=MPI_THREAD_MULTIPLE;
  int provided;
  MPI_Init_thread(&argc, &argv, required, &provided);
  assert(required<=provided);
  MPI_Comm_size(MPI_COMM_WORLD, &local.num_mpi_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &local.mpi_rank);
  MPI_Get_processor_name(local.processor_name, &name_len);
  assert(local.num_mpi_ranks == params.Pr*params.Pc);

  MKL_INT zero = 0, one = 1, dummy;
  MKL_INT info, negone=-1, ten=10;
  char order = 'R';  // row-major process assignment
  MKL_INT mypnum, nprocs;
  local.ctxt=MPI_COMM_WORLD;
  blacs_pinfo(&mypnum, &nprocs);
  blacs_get(&local.ctxt, &zero, &local.ctxt);
  blacs_gridinit(&local.ctxt, &order, &params.Pr, &params.Pc);
  blacs_gridinfo(&local.ctxt, &dummy, &dummy, &local.pr, &local.pc);
  // don't use Pr and Pc from gridinfo, but from earlier input

  local.rows = numroc(&params.dimension, &params.block_size, &local.pr, &zero, &params.Pr);
  local.cols = numroc(&params.dimension, &params.block_size, &local.pc, &zero, &params.Pc);
}

void constructMatrix(SubmatrixDescriptor& descr, Eigen::MatrixXf& matrix)
{
  matrix = Eigen::MatrixXf::Zero(descr.rows, descr.cols);
  for (int i=descr.i; i<descr.i+descr.rows; ++i)
    for(int j=descr.j; j<descr.j+descr.cols; ++j)
    {
      // matrix elements of (b^\dagger + b)^2
      float val=0;
      if (i==j+2)
        val = std::sqrt((j+1)*(j+2));
      else if (i==j)
        val = (2*j+1);
      else if (i==j-2)
        val = std::sqrt(j*(j-1));
      matrix(i-descr.i,j-descr.j) = val;
      // matrix(i-start_i, j-start_j) = i+100*j;
    }

}

void sendMatrixChunk(
    const RunParameters& params,
    const LocalParameters& local,
    const SubmatrixDescriptor& chunk,
    Eigen::MatrixXf& chunk_matrix
  )
{
  // descriptor for block to be sent to other rank
  SubmatrixDescriptor block;

  // loop over blocks of the matrix
  for (block.i = chunk.i; block.i < chunk.i+chunk.rows; block.i += block.rows)
    for (block.j = chunk.j; block.j < chunk.j+chunk.cols; block.j += block.cols)
    {
      int target_rank = getRankForIndices(
          block.i, block.j,
          params.Pr, params.Pc,
          params.block_size, params.block_size
        );
      // we send from the start index to the edge of the block, or to the
      // end of the submatrix
      block.rows = std::min(
          params.block_size-(block.i%params.block_size),
          (chunk.i+chunk.rows)-block.i
        );
      block.cols = std::min(
          params.block_size-(block.j%params.block_size),
          (chunk.j+chunk.cols)-block.j
        );
      std::cout << fmt::format(
          "{:d} send ({:2d},{:2d})+({:2d},{:2d})->{:d}\n",
          local.mpi_rank, block.i, block.j, block.rows, block.cols, target_rank
        );

      // determine size of packed data
      int int_pack_size, row_pack_size, total_pack_size;
      MPI_Pack_size(4, MPI_INTEGER, MPI_COMM_WORLD, &int_pack_size);
      MPI_Pack_size(block.rows, MPI_FLOAT, MPI_COMM_WORLD, &row_pack_size);
      total_pack_size = int_pack_size + block.cols*row_pack_size;

      // allocate buffer for pack
      std::vector<unsigned char> buffer(total_pack_size);

      // pack index information
      int position = 0;
      MPI_Pack(
          &block, 4, MPI_INTEGER,
          buffer.data(), total_pack_size, &position,
          MPI_COMM_WORLD
        );

      // pack submatrix data
      //
      // (x,y) denote the indices in chunk_matrix
      // therefore, x is the offset between the current row
      // and the first row of the chunk. Likewise, y is
      // the offset between the current column and the first column
      // of the chunk
      int x = block.i - chunk.i;
      for (int y=block.j-chunk.j; y < block.cols+(block.j-chunk.j); ++y)
      {
        MPI_Pack(
            chunk_matrix.data()+(x+chunk.rows*y), block.rows, MPI_FLOAT,
            buffer.data(), total_pack_size, &position,
            MPI_COMM_WORLD
          );
      }

      // send packed data to destination node
      MPI_Send(buffer.data(), position, MPI_PACKED, target_rank, 0, MPI_COMM_WORLD);
    }
}

void receiveMatrixChunk(
    RunParameters& params,
    LocalParameters& local,
    MPI_Status& status,
    std::vector<unsigned char>& buffer,
    Eigen::MatrixXf& local_storage
  )
{
  // receive next chunk
  int max_pack_size = buffer.size();
  MPI_Recv(
      buffer.data(), max_pack_size, MPI_PACKED,
      status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status
    );
  int position=0;
  SubmatrixDescriptor block;
  MPI_Unpack(buffer.data(), max_pack_size, &position, &block, 4, MPI_INTEGER, MPI_COMM_WORLD);
  int s,t;
  std::tie(s,t) = calculateProcessLocalElementIndices(
      block.i, block.j,
      params.Pr, params.Pc, params.block_size, params.block_size
    );
  std::cout << fmt::format(
      "{:d} recv ({:2d},{:2d})+({:2d},{:2d})<-{:d} @ ({:2d},{:2d})\n",
      local.mpi_rank, block.i, block.j, block.rows, block.cols,
      status.MPI_SOURCE, s, t
    );

  // int s, t;
  assert(t+block.cols<=local_storage.cols());
  const int& i = block.i;
  for (int j = block.j; j < block.j+block.cols; ++j)
  {
    // extract pack in place
    std::tie(s,t) = calculateProcessLocalElementIndices(
        i, j,
        params.Pr, params.Pc, params.block_size, params.block_size
      );
    assert(s+block.rows<=local_storage.rows());
    MPI_Unpack(buffer.data(), max_pack_size, &position, local_storage.data()+(s+local.rows*t), block.rows, MPI_FLOAT, MPI_COMM_WORLD);
  }
}

int main(int argc, char *argv[])
{
  RunParameters params;
  LocalParameters local;

  initializeRunParameters(params);
  initializeParallel(argc, argv, params, local);

  // determine needed storage size and allocate local matrix
  Eigen::MatrixXf local_storage = Eigen::MatrixXf::Zero(local.rows, local.cols);
  // omp_set_num_threads(2);

  std::cout << fmt::format(
                  "{:d} (pr,pc)=({:d},{:d}) (rows,cols)=({:d},{:d})",
                  local.mpi_rank, local.pr, local.pc, local.rows, local.cols
                )
            << std::endl << std::flush;

  MPI_Barrier(MPI_COMM_WORLD);

  #pragma omp parallel
  {
    #pragma omp single nowait
    // one thread is dedicated to receiving matrix blocks and placing them
    // into local storage
    {
      std::cout << fmt::format(
          " rank {:d}, thread {:d} -- listening\n",
          local.mpi_rank, omp_get_thread_num()
        ) << std::flush;
      // allocate buffer
      int int_pack_size, row_pack_size, max_pack_size;
      MPI_Pack_size(4, MPI_INTEGER, MPI_COMM_WORLD, &int_pack_size);
      MPI_Pack_size(params.block_size, MPI_FLOAT, MPI_COMM_WORLD, &row_pack_size);
      max_pack_size = int_pack_size + params.block_size*row_pack_size;
      std::vector<unsigned char> buffer(max_pack_size);
      int done_count=0;

      // receive matrix
      while (true)
      {
        // figure out what type of message this is
        int source, tag;
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == k_done)
        {
          // done receiving, break out
          MPI_Recv(buffer.data(), 1, MPI_INTEGER, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
          if (local.mpi_rank == 0)
          {
            ++done_count;
            std::cout << fmt::format("done from {:d}, {:d}/{:d}\n", status.MPI_SOURCE, done_count, local.num_mpi_ranks);
            if (done_count == local.num_mpi_ranks)
            {
              std::cout << fmt::format("***all done***\n");

              // inform all workers that matrix distribution is finished
              for (int target_rank=1; target_rank<local.num_mpi_ranks; ++target_rank)
              {
                MPI_Send(&k_done, 1, MPI_INTEGER, target_rank, k_done, MPI_COMM_WORLD);
              }
              break;
            }
          }
          else  // (local.mpi_rank != 0)
          {
            break;
          }
        }
        else
        {
          receiveMatrixChunk(params, local, status, buffer, local_storage);
        }
      }
    }  // end omp single

    #pragma omp single
    // For the present example, matrix construction is single-threaded.
    // In general, one should make matrix generation multi-threaded.
    {
      std::cout << fmt::format(" rank {:d}, thread {:d} -- distributing\n", local.mpi_rank, omp_get_thread_num());
      // construct and distribute matrix
      Eigen::MatrixXf matrix;
      int cols_per_proc = params.dimension/local.num_mpi_ranks;
      int start_j = cols_per_proc*local.mpi_rank;
      int num_cols = (local.mpi_rank==local.num_mpi_ranks-1 ? params.dimension-start_j : cols_per_proc);

      SubmatrixDescriptor matrix_descr = {0, start_j, params.dimension, num_cols};
      constructMatrix(matrix_descr, matrix);
      // dump accumulated diagonistic info
      // for (int i=0; i<local.num_mpi_ranks; ++i)
      // {
      //   MPI_Barrier(MPI_COMM_WORLD);
      //   if (local.mpi_rank==i)
      //     std::cout << matrix << std::endl << std::endl;
      //   MPI_Barrier(MPI_COMM_WORLD);
      // }

      sendMatrixChunk(params, local, matrix_descr, matrix);

      // inform master that we've finished distributing matrix
      MPI_Send(&k_done, 1, MPI_INTEGER, 0, k_done, MPI_COMM_WORLD);
    }  // end omp single
  }  // end omp parallel

  MPI_Barrier(MPI_COMM_WORLD);

  // check matrix distribution
  std::stringstream line_stream;
  line_stream << "local " << local.mpi_rank << std::endl;
  line_stream << local_storage << std::endl;

  // set up ScaLAPACK environment
  MKL_INT zero=0, one=1, info;
  std::unique_ptr<MKL_INT> descA(new MKL_INT[9]);  // descriptor for matrix
  std::unique_ptr<MKL_INT> descZ(new MKL_INT[9]);  // descriptor for eigenvectors
  descinit(descA.get(),
      &params.dimension, &params.dimension, &params.block_size, &params.block_size,
      &zero, &zero, &local.ctxt, &local.rows, &info
    );
  descinit(descZ.get(),
      &params.dimension, &params.dimension, &params.block_size, &params.block_size,
      &zero, &zero, &local.ctxt, &local.rows, &info
    );

  // construct storage for eigensystem and work space
  char jobz='V', uplo='L';
  MKL_INT lwork;
  Eigen::VectorXf eigval_storage = Eigen::VectorXf(params.dimension);
  Eigen::MatrixXf eigvec_storage = Eigen::MatrixXf(local.rows, local.cols);
  std::vector<float> work(2);

  // dummy call to determine needed size for work
  lwork = -1;
  pssyev(&jobz, &uplo, &params.dimension,
      local_storage.data(), &one, &one, descA.get(),
      eigval_storage.data(), eigvec_storage.data(), &one, &one, descZ.get(),
      work.data(), &lwork, &info);

  // resize to actual needed work size
  lwork = int(work[0]);
  work.resize(lwork);

  // do actual diagonalization
  pssyev(&jobz, &uplo, &params.dimension,
      local_storage.data(), &one, &one, descA.get(),
      eigval_storage.data(), eigvec_storage.data(), &one, &one, descZ.get(),
      work.data(), &lwork, &info);

  // output matrix
  line_stream << "eigval " << local.mpi_rank << std::endl;
  line_stream << eigval_storage.transpose() << std::endl;
  line_stream << "eigvec " << local.mpi_rank << std::endl;
  line_stream << eigvec_storage << std::endl;

  // dump accumulated diagonistic info
  for (int i=0; i<local.num_mpi_ranks; ++i)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    if (local.mpi_rank==i)
      std::cout << line_stream.str() << std::endl;
  }


  // Finalize the MPI environment.
  MPI_Finalize();

  return EXIT_SUCCESS;
}

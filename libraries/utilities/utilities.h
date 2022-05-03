/****************************************************************
  utilities.h

  Define arithmetic shorthands.

  Mark A. Caprio and Anna E. McCoy
  University of Notre Dame

  SPDX-License-Identifier: MIT

  Created by Mark A. Caprio on 2/17/11.
  Drawing upon libmc/mcutils C code.
  2/23/11 (mac): Renamed from mc_arithmetic to arithmetic.
  3/9/16 (mac): Imported into spncci project as utilities.h.
  3/26/17 (mac): Add deprecation notes.
  6/11/19 (aem): Add FileExists function to check if file exits
  
****************************************************************/

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <fstream>
#include <unistd.h>
#include "fmt/format.h"
// #include "gsl/gsl_sf_gamma.h"
#include <Eigen/Eigen>
#include "basis/operator.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"

// extern double zero_threshold;

// Based on https://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}


namespace utils
{
inline bool FileExists(std::string filename, bool verbose)
  {
    int res = access(filename.c_str(), R_OK);
    if (res < 0)
      {
        if (errno == ENOENT)
          {
            if(verbose)
              std::cout<<filename<<" does not exit"<<std::endl;
          }
        else if (errno == EACCES)
          {
            if(verbose)
              std::cout<<filename<<" not readible"<<std::endl;
            exit(EXIT_FAILURE);
          }
        else
          {
            if(verbose)
              std::cout<<"unknown error associated with "<<filename<<std::endl;
          }
        return false;
      }

    else
      return true;
  }


 inline std::string get_spncci_project_root_dir()
  {
    const char *tmp = std::getenv("SPNCCI_PROJECT_ROOT_DIR");
      std::string spncci_root_dir(tmp ? tmp : "");
      if (spncci_root_dir.empty()) {
          std::cerr << "[ERROR] SPNCCI_PROJECT_ROOT_DIR not defined." << std::endl;
          exit(EXIT_FAILURE);
      }
    return spncci_root_dir;
  }


  // template <typename tDataType>
  //   void WriteBinary(std::ostream& os, const tDataType &data)

  inline void WriteOperatorBlockBinary(const basis::OperatorBlock<double>& block, std::ofstream& output)
    {
          // write number of rows and columns
          int num_rows=block.rows();
          int num_cols=block.cols();
          mcutils::WriteBinary<int>(output,num_rows);
          mcutils::WriteBinary<int>(output,num_cols);
          // dump matrix to file.  Order is column major (Eigen default for .data())
          mcutils::WriteBinary<double>(output,block.data(),num_rows*num_cols);

    }


  inline void WriteOperatorBlockBinary(const basis::OperatorBlock<double>& block, const std::string& filename)
    {
          std::ios_base::openmode mode_argument = std::ios_base::out | std::ios::app | std::ios_base::binary;
          std::ofstream output;
          output.open(filename,mode_argument);
          WriteOperatorBlockBinary(block, output);
          output.close();

    }

  inline basis::OperatorBlock<double> ReadOperatorBlockBinary(const std::string& filename)
    {
      int rows,cols;
      std::ios_base::openmode mode_argument = std::ios_base::in;
      mode_argument |= std::ios_base::binary;
      std::ifstream is(filename, mode_argument);
      mcutils::StreamCheck(bool(is),filename,fmt::format("Failure opening {}",filename));
      mcutils::ReadBinary<int>(is,rows);
      mcutils::ReadBinary<int>(is,cols);
      double buffer[rows*cols];

      mcutils::ReadBinary(is,buffer,rows*cols);
      basis::OperatorBlock<double> block=Eigen::Map<Eigen::MatrixXd>(buffer,rows,cols);
      return block;

    }

  //   // num_rows, num_cols, rmes
  //   // rmes are by column then by row
  //   // output in binary mode
  //   std::ios_base::openmode mode_argument = std::ios_base::in;
  //   mode_argument |= std::ios_base::binary;
  //   std::ifstream in_stream(filename, mode_argument);
  //   mcutils::StreamCheck(bool(in_stream),filename,"Failure opening lsu3shell rme file");
  //   // if (!in_stream)
  //   //  {
  //   //     std::cerr << "Could not open file '" << filename << "'!" << std::endl;
  //   //     return;
  //   //  }

  //   int binary_format_code;
  //   int binary_float_precision;
  //   mcutils::ReadBinary<int>(in_stream,binary_format_code);
  //   mcutils::ReadBinary<int>(in_stream,binary_float_precision);


  //   lgi_expansions.resize(num_lgi_subspaces);
  //   for(int i=0; i<num_lgi_subspaces; ++i)
  //     {
  //       basis::OperatorBlock<double>& block=lgi_expansions[i];
  //       lgi::LGIIndexType rows, cols;
  //       // Read in number of rows and cols

  //       /////////////////////////////////////////////////////////////////////////////////
  //       // Read in RMEs and cast to double matrix
  //       if(binary_float_precision==4)
  //         {
  //           float buffer[rows*cols];
  //           in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
  //           block=Eigen::Map<Eigen::MatrixXf>(buffer,rows,cols).cast<double>();
  //         }
  //       else if (binary_float_precision==8)
  //         {
  //           double buffer[rows*cols];
  //           in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
  //           block=Eigen::Map<Eigen::MatrixXd>(buffer,rows,cols);
  //         }
  //     }
  // }










// // ONLYIF(cond,x) evaluates and returns x only if cond is true
// #define ONLYIF(cond,x) ( (cond) ? (x) : 0)

// // sqr(x) returns the arithmetic square of x by self-multiplication
// //   Note: Use of inline template avoids double evaluation of x which
// //   would occur in a macro implementation.

// template <typename T>
// inline
// T sqr(const T& x)
// {
//   return x*x;
// }

// template <typename T>
// inline
// int KroneckerDelta(const T& x, const T& y)
// // Return Kronecker delta of variables x and y.
// //
// // That is, returns 1 if x==y, 0 otherwise.
// {
//   return int(x==y);  // Is int(true) guaranteed to be 1?
// }



// inline int Choose(int x, int y)
// {
//   int choose=0;
//   if ((x>=y)&&(y>=0))
//     {
//       gsl_sf_result result;
//       gsl_sf_choose_e(x,y,&result);
//       choose=result.val;
//     }
//   return choose;
// }

// inline double Factorial(int x)
// {
//   return gsl_sf_fact(x);
// }

// inline int parity(const int i)
// {
//   if( (i%2)==0)
//     return 1;
//   else
//     return -1;
// }


// inline bool CheckIfZeroMatrix(const Eigen::MatrixXd& matrix, double zero_threshold=0)
//     {
//       int rows=matrix.rows();
//       int cols=matrix.cols();
//       for(int j=0; j<cols; ++j)
//         for(int i=0; i<rows; ++i)
//           {
//             if(fabs(matrix(i,j))>zero_threshold)
//                 return false;
//           }
//       return true;
//     }

// inline void ZeroOutMatrix(basis::OperatorBlocks<double>& matrix_vector,double threshold)
//   {
//     for(auto& matrix : matrix_vector)
//       for(int i=0; i<matrix.rows(); ++i)
//         for(int j=0; j<matrix.cols(); ++j)
//           {
//             if(fabs(matrix(i,j))<threshold)
//               matrix(i,j)=0;
//           }
//   }

// inline void NormalizeMatrix(Eigen::MatrixXd& matrix, std::string type)
//   {
//     int rows=matrix.rows();
//     int cols=matrix.cols();
//     if(type=="column")
//       {
//         for(int j=0; j<cols; ++j)
//           {
//             int sum=0;
//             for(int i=0; i<rows; ++i)
//                 sum+=matrix(i,j)*matrix(i,j);
//             matrix.block(0,j,rows,1)*=std::sqrt(sum);
//           }
//       }
//     if (type=="row")
//       {
//         for(int i=0; i<rows; ++i)
//           {
//             int sum=0;
//             for(int j=0; j<cols; ++j)
//               sum+=matrix(i,j)*matrix(i,j);
//             matrix.block(i,0,1,cols)*=std::sqrt(sum);
//           }
//       }
//   }
}
#endif

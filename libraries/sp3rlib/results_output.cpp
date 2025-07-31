/****************************************************************
  results_output.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/08/24 (cvc): Created.

****************************************************************/

#include "sp3rlib/results_output.h"

namespace sp3r
{

  void StartNewSection(std::ostream& out_stream, const std::string& title)
  {
    out_stream << std::endl;
    out_stream << fmt::format("[{}]",title) << std::endl;
  }

void WriteEigenvalues(
    std::ostream& out_stream,
    const std::vector<std::pair<HalfInt,Halfint>>& SJ_values,
    const std::vector<sp3r::Vector>& eigenvalues
  )
{
  StartNewSection(out_stream,"Energies");
  out_stream << "# J S i E" << std::endl;
  for (int subspace_index=0; subspace_index<SJ_values.size(); ++subspace_index)
    {
      // Get operator and coefficient
      const auto& [S,J] = SJ_values[subspace_index];
      const Eigen::VectorXd& eigenvalues_SJ = eigenvalues[subspace_index];

      // iterate over eigenvalues in J subspace
      for (int eigenstate_index=0; eigenstate_index<eigenvalues_SJ.size(); ++eigenstate_index)
        out_stream
          << fmt::format("{:4.1f} {:4.1f} {:3d} {:+9.4f}",double(J),double(S),eigenstate_index,eigenvalues_SJ[eigenstate_index])
          << std::endl;
    }
}


void WriteEigenvectors(
  sp3r::Matrix& eigenvectors,
  const HalfInt& S,
  const HalfInt& J,
  std::ofstream& out_file,
  const int& binary_float_precision
  )
  {
    int rows=eigenvectors.rows();
    int cols=eigenvectors.cols();
    mcutils::WriteBinary<int>(out_file,TwiceValue(J));
    mcutils::WriteBinary<int>(out_file,TwiceValue(S));
    mcutils::WriteBinary<int>(out_file,binary_float_precision);
    mcutils::WriteBinary<int>(out_file,rows);
    mcutils::WriteBinary<int>(out_file,cols);

    int size=rows*cols;

    // write matrix.  Order is column major (Eigen default)
    if(binary_float_precision==4)
      {
        Eigen::MatrixXf buffer_matrix=eigenvectors.cast<float>();
        out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*binary_float_precision);
      }

    else if (binary_float_precision==8)
      {
        Eigen::MatrixXd buffer_matrix=eigenvectors;
        out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*binary_float_precision);

      }
  }


void WriteEigenvectors(
  sp3r::Matrix& eigenvectors,
  const HalfInt& S,
  const HalfInt& J,
  const std::string& filename,
  const int& binary_float_precision
  )
  {
    std::ios_base::openmode mode_argument = std::ios_base::out | std::ios_base::binary;
    std::ofstream out_file;
    out_file.open(filename,mode_argument);

    if (!out_file)
      {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
      }
    sp3r::WriteEigenvectors(eigenvectors,S,J,out_file,binary_float_precision);
  }



void ReadEigenvectors(
    const std::string& filename,
    sp3r::Matrix& eigenvectors
  )
  {
    std::ios_base::openmode mode_argument = std::ios_base::in | std::ios_base::binary;
    std::ifstream in_stream;
    in_stream.open(filename,mode_argument);

    if (!in_stream)
      {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
      }


    int twice_J,twice_S,binary_float_precision, rows, cols;
    mcutils::ReadBinary<int>(in_stream,twice_J);
    mcutils::ReadBinary<int>(in_stream,twice_S);
    mcutils::ReadBinary<int>(in_stream,binary_float_precision);
    mcutils::ReadBinary<int>(in_stream,rows);
    mcutils::ReadBinary<int>(in_stream,cols);

    HalfInt J=HalfInt(twice_J,2);
    HalfInt S=HalfInt(twice_S,2);

    // Read matrix.  Order is column major (Eigen default)
    if(lgi::binary_float_precision==4)
      {
        float buffer[rows*cols];
        in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
        eigenvectors
            =Eigen::Map<Eigen::MatrixXf>(buffer,rows,cols).cast<double>();
      }
    else if (lgi::binary_float_precision==8)
      {
        double buffer[rows*cols];
        in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
        eigenvectors
          =Eigen::Map<Eigen::MatrixXd>(buffer,rows,cols);
      }
  }



}
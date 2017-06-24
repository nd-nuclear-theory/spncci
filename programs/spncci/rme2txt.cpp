/****************************************************************
  rme2txt.cpp

  Dump binary RME file to text format.

  Output format:
    i j num_rmes
      rme0
      rme1
      ...
    ...

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/12/17 (mac): Created.
  6/23/17 (mac): Update to read rme binary file format v1,
    using code from lsu3shell_rme.
****************************************************************/

#include <fstream>
#include <iostream>
#include <string>

#include "cppformat/format.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"

typedef short unsigned int RMEIndexType;

int main(int argc, char **argv)
{

  // process arguments
  if(argc<1+1)
    {
      std::cout<<"Syntax : <filename>"<<std::endl;
      std::exit(EXIT_SUCCESS);
    }
  std::string in_stream_filename(argv[1]);

  // open file
  std::ifstream in_stream(in_stream_filename,std::ios_base::in|std::ios_base::binary);
  mcutils::StreamCheck(bool(in_stream),in_stream_filename,"Failure opening lsu3shell rme file");

  // read file header
  int format_code;
  mcutils::ReadBinary<int>(in_stream,format_code);
  assert(format_code==1);
  int float_precision;
  mcutils::ReadBinary<int>(in_stream,float_precision);
  assert((float_precision==4)||(float_precision==8));
  std::cout
    << fmt::format("RME input: filename {}, format_code {}, float_precision {}",in_stream_filename,format_code,float_precision)
    << std::endl;

  // process file data
  while (in_stream)
    {

      // read bra/ket lsu3shell basis multiplicity group indices
      RMEIndexType i, j;
      mcutils::ReadBinary<RMEIndexType>(in_stream,i);
      mcutils::ReadBinary<RMEIndexType>(in_stream,j);

      // quit if this has brought us past end of file
      if (!in_stream)
        break;

      // read multiplicity given in file
      RMEIndexType num_rmes;
      mcutils::ReadBinary<RMEIndexType>(in_stream,num_rmes);

      // write entry header
      std::cout << fmt::format("{:d} {:d} {:d}",i,j,num_rmes) << std::endl;

      // dump rmes
      for (int rme_index=0; rme_index<num_rmes; ++rme_index)
        {
          // read rme
          double rme;
          if (float_precision==4)
            {
              float rme_float;
              mcutils::ReadBinary<float>(in_stream,rme_float);
              rme = rme_float;
            }
          else if (float_precision==8)
            {
              mcutils::ReadBinary<double>(in_stream,rme);
            }

          // write rme
          std::cout
            << fmt::format("  {:+19.12e}",rme)
            << std::endl;
        }
    }

}

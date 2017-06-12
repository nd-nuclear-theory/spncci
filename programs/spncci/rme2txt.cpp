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
****************************************************************/

#include <fstream>
#include <iostream>
#include <string>

#include "cppformat/format.h"
#include "mcutils/io.h"

int main(int argc, char **argv)
{

  // process arguments
  if(argc<1+1)
    {
      std::cout<<"Syntax : <filename>"<<std::endl;
      std::exit(EXIT_SUCCESS);
    }
  std::string in_stream_filename(argv[1]);

  // process file
  std::ifstream in_stream(in_stream_filename,std::ios_base::in|std::ios_base::binary);
  while (in_stream)
    {

      // read indexing
      int i, j, num_rmes;
      mcutils::ReadBinary<int>(in_stream,i);
      mcutils::ReadBinary<int>(in_stream,j);
      mcutils::ReadBinary<int>(in_stream,num_rmes);

      // quit if this has brought us past end of file
      if (!in_stream)
        break;

      // dump entry
      std::cout << fmt::format("{:d} {:d} {:d}",i,j,num_rmes) << std::endl;
      for (int rme_index=0; rme_index<num_rmes; ++rme_index)
        {
          float rme;
          mcutils::ReadBinary<float>(in_stream,rme);
          std::cout << fmt::format("  {:12.6f}",rme) << std::endl;
        }
    }
  

}

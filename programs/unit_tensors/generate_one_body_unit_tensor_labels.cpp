/**************************************************************
Generates one-body unit tensor labels and writes them into file
obd.dat to be used by SU3RME_Alexis and SU3RME_Alexis_digest.
obd.dat contains lines with 5 integers:

  N' N lambda0 mu0 2*S0
***************************************************************/
#include <iostream>
#include <fstream>
#include "u3shell/relative_operator.h"

int main(int argc, char **argv)
{
  if (argc != 4) {
    std::cout << "Usage: " << argv[0] << " Nmax N1vp N1vn" << std::endl;
    return EXIT_FAILURE;
  }

  int Nmax = std::stoi(argv[1]);
  int N1vp = std::stoi(argv[2]);
  int N1vn = std::stoi(argv[3]);

  // Generate one-body unit tensor labels
  std::vector<u3shell::OneBodyUnitTensorLabelsU3S> unit_tensor_labels;
  u3shell::GenerateOneBodyUnitTensorLabelsU3S(Nmax, N1vp, N1vn, unit_tensor_labels);

  // Write one-body unit tensor labels into file obd.dat
  std::ofstream file("obd.dat");
  if (!file) {
    std::cout << "generate_one_body_unit_tensor_labels: Could not open file obd.dat" << std::endl;
    return EXIT_FAILURE;
  }
  int skip;
  if (N1vn>N1vp) {
    skip = 1;
  } else {
    skip = -1;
  }
  for (u3shell::OneBodyUnitTensorLabelsU3S labels : unit_tensor_labels) {
    if (labels.Tz() == skip) continue;
    file << labels.Nbra() << " " << labels.Nket() << " " << labels.x0().lambda() << " "
	 << labels.x0().mu() << " " << labels.S0().TwiceValue() << std::endl;
  }

  return EXIT_SUCCESS;
}

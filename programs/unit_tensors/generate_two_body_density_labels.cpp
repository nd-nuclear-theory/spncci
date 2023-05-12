/**************************************************************
Generates two-body density labels and writes them into file
tbd.dat to be used by SU3RME_Alexis_2B and SU3RME_Alexis_2B_digest.
tbd.dat contains lines with 13 integers:

  N1 N2 N3 N4 lmf muf SSf lmi mui SSi lm0 mu0 SS0
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

  // Generate two-body density labels
  std::vector<u3shell::TwoBodyDensityLabels> tbd_labels;
  u3shell::GenerateTwoBodyDensityLabels(Nmax, N1vp, N1vn, tbd_labels);

  int skip;
  if (N1vn>N1vp) {
    skip = 1;
  } else {
    skip = -1;
  }
  std::vector<std::array<int,13>> tensors;
  for (u3shell::TwoBodyDensityLabels labels : tbd_labels) {
    if (labels.Tz() == skip) continue;
    std::array<int,13> tensor={labels.N1(),labels.N2(),labels.N3(),labels.N4(),
                               labels.xf().lambda(),labels.xf().mu(),2*labels.Sf(),
                               labels.xi().lambda(),labels.xi().mu(),2*labels.Si(),
                               labels.x0().lambda(),labels.x0().mu(),labels.S0().TwiceValue()};
    auto it=std::find(tensors.begin(), tensors.end(), tensor);
    if(it==tensors.end())tensors.push_back(tensor);
  }

  // Write two-body density labels into file tbd.dat
  std::ofstream file("tbd.dat");
  if (!file) {
    std::cout << "ERROR: generate_two_body_density_labels: Could not open file tbd.dat" << std::endl;
    return EXIT_FAILURE;
  }
  for (std::array<int,13> tensor : tensors) {
    file << tensor[0] << " " << tensor[1] << " " << tensor[2] << " " << tensor[3] << " "
	 << tensor[4] << " " << tensor[5] << " " << tensor[6] << " "
	 << tensor[7] << " " << tensor[8] << " " << tensor[9] << " "
	 << tensor[10] << " " << tensor[11] << " " << tensor[12] << std::endl;
  }

  return EXIT_SUCCESS;
}

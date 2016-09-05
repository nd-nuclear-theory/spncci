#include <SU3ME/ModelSpaceExclusionRules.h>
#include <LSU3/ncsmSU3xSU2Basis.h>
#include <SU3ME/global_definitions.h>
#include <vector>
#include <stack>
#include <ctime>
using namespace std;

void IterateOverBasisTabularOutput(const lsu3::CncsmSU3xSU2Basis &basis) {
  uint16_t ap_max, an_max;
  uint16_t irrep_dim(0);
  uint32_t idim(0);
  uint64_t firstStateId = basis.getFirstStateId();

  uint32_t ip;
  uint32_t in;

  SingleDistribution distr_p;
  SingleDistribution distr_n;
  UN::SU3xSU2_VEC gamma_p;
  UN::SU3xSU2_VEC gamma_n;
  SU3xSU2_VEC omega_p, omega_n;
  SU3xSU2::LABELS omega_pn;

  //	loop over (ip, in) pairs
  for (int ipin_block = 0; ipin_block < basis.NumberOfBlocks(); ipin_block++) {
    if (!basis.NumberOfStatesInBlock(ipin_block)) {
      continue;
    }
    ip = basis.getProtonIrrepId(ipin_block);
    in = basis.getNeutronIrrepId(ipin_block);

    int N = basis.nhw_p(ip) + basis.nhw_n(in);

    ap_max = basis.getMult_p(ip);
    an_max = basis.getMult_n(in);
    SU3xSU2::LABELS irrep_p(basis.getProtonSU3xSU2(ip));
    SU3xSU2::LABELS irrep_n(basis.getProtonSU3xSU2(in));

    int lmp(irrep_p.lm);
    int mup(irrep_p.mu);
    int ssp(irrep_p.S2);

    int lmn(irrep_n.lm);
    int mun(irrep_n.mu);
    int ssn(irrep_n.S2);

    for (int iwpn = basis.blockBegin(ipin_block);
         iwpn < basis.blockEnd(ipin_block); ++iwpn) {
      omega_pn = basis.getOmega_pn(ip, in, iwpn);
      int SS = omega_pn.S2;
      int lm = omega_pn.lm;
      int mu = omega_pn.mu;
      int rho0_max = omega_pn.rho;
      std::cout << N << " " << ip << " " << ap_max << " " << lmp << " " << mup << " " << ssp << " ";
      std::cout << in << " " << an_max << " " << lmn << " " << mun << " " << ssn << " ";
      std::cout << rho0_max << " " << lm << " " << mu << " " << SS << std::endl;
    }
  }
}

int main(int argc, char **argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0]
         << " <model space definition> " << endl;
    cout << "Implicitly ndiag=1 and idiag=0." << endl;
    return EXIT_FAILURE;
  }

  int ndiag = 1;
  int idiag = 0;

  proton_neutron::ModelSpace ncsmModelSpace(argv[1]);
  if (ncsmModelSpace.JJ() != 255)
  {
     std::cerr << "The value of 2J has to be equal to 255! The current value: "
               << (int)ncsmModelSpace.JJ() << std::endl;
    return EXIT_FAILURE;
  }
  lsu3::CncsmSU3xSU2Basis basis(ncsmModelSpace, idiag, ndiag);
  IterateOverBasisTabularOutput(basis);
}

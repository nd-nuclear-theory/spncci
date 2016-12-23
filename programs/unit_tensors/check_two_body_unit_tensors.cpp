/****************************************************************
  check_two_body_unit_tensors.cpp

  Read in RMEs for two body tensors and check against expected
  values.

  Syntax:
    check_two_body_unit_tensors Z=1 N=1 Nmax Nstep 2*Nsigma_0=6

  Example:
    check_two_body_unit_tensors 1 1 0 1 6

  Input files:
    lsu3shell_basis.dat -- lsu3shell tabular basis listing file
    two_body_unit_*.rme -- output of SU3RME for each operator

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/7/16 (aem,mac): Created, starting frrom
    generate_lsu3shell_two_body_tensors.cpp.
****************************************************************/

#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_rme.h"
#include "spncci/sp_basis.h"
#include "sp3rlib/u3coef.h"
// #include "spncci/lgi_unit_tensor.h"
#include "u3shell/unit_tensor_expansion.h"
#include "u3shell/two_body_operator.h"


////////////////////////////////////////////////////////////////
// generate known matrix elements
////////////////////////////////////////////////////////////////

  void 
  GenerateTwoBodyUnitTensorMatrices(      
    const u3shell::TwoBodyUnitTensorLabelsU3ST& two_body_unit_tensor_labels,
    const lsu3shell::U3SPNBasisLSU3Labels& basis_provenance,
    const u3shell::SpaceU3SPN& space,
    u3shell::SectorsU3SPN& sectors,
    basis::MatrixVector& matrices
    )
  // Generate matrix representation of two-body unit tensor.
  //
  // The given space must be a deuteron-like two-body space.
  {
    // set up zero initialized operator
    sectors = u3shell::SectorsU3SPN(space,u3shell::OperatorLabelsU3S(two_body_unit_tensor_labels),false);
    basis::SetOperatorToZero(sectors,matrices);

    // unpack unit tensor labels
    //
    // Notation: b = "bar", p = "prime"

    u3shell::OperatorLabelsU3ST::KeyType tensor_labels;
    int rho0;
    u3shell::TwoBodyStateLabelsU3ST::KeyType bra_labels, ket_labels;

    std::tie(tensor_labels,rho0,bra_labels,ket_labels)
      = two_body_unit_tensor_labels.Key();

    int N0; u3::SU3 x0; HalfInt S0, T0; int g0;
    std::tie(N0,x0,S0,T0,g0)
      = tensor_labels;

    int eta1bp, eta2bp; u3::SU3 xbp; HalfInt Sbp, Tbp;
    std::tie(eta1bp,eta2bp,xbp,Sbp,Tbp) = bra_labels;
    HalfInt Nbp = eta1bp+eta2bp+3;  // U(1) label includes zero-point
    u3::U3S omegaSbp = u3::U3S(u3::U3(Nbp,xbp),Sbp);
    int gbp = (eta1bp+eta2bp) % 2;

    int eta1b, eta2b; u3::SU3 xb; HalfInt Sb, Tb;
    std::tie(eta1b,eta2b,xb,Sb,Tb) = ket_labels;
    HalfInt Nb = eta1b+eta2b+3;  // U(1) label includes zero-point
    u3::U3S omegaSb = u3::U3S(u3::U3(Nb,xb),Sb);
    int gb = (eta1b+eta2b) % 2;


    // iterate over sectors
    for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
      {
        // retrieve sector
        const u3shell::SectorsU3SPN::SectorType& sector = sectors.GetSector(sector_index);
        const int rho0p = sector.multiplicity_index();

        // short-circuit -- check sector symmetry labels
        //
        // Do they match those for the nonzero matrix element of the
        // unit tensor?
        //
        // Note: Can ignore bra/ket Sp & Sn labels (trivial on
        // deuteron) and U(3) N labels
        bool matches = true;
        matches &= (rho0p == rho0);
        matches &= (sector.bra_subspace().U3S() == omegaSbp);
        matches &= (sector.ket_subspace().U3S() == omegaSb);
        if (!matches)
          continue;

        // iterate over matrix elements
        for (int bra_index=0; bra_index<sector.bra_subspace().size(); ++bra_index)
          for (int ket_index=0; ket_index<sector.ket_subspace().size(); ++ket_index)
            {

              // retrieve state shell labels
              const lsu3shell::LSU3BasisGroupLabels& bra_labels = basis_provenance[sector.bra_subspace_index()][bra_index];
              const lsu3shell::LSU3BasisGroupLabels& ket_labels = basis_provenance[sector.ket_subspace_index()][ket_index];
              int eta1p = bra_labels.Np;
              int eta2p = bra_labels.Nn;
              int eta1 = ket_labels.Np;
              int eta2 = ket_labels.Nn;

              // canonicalize labels
              //
              //   for comparison with U3ST scheme unit "source matrix
              //   element" labels
              //
              // U3ST two-body state interchange phase factor is
              //
              //   ~ 1+g+omega+S+T
              //
              // Though our target matrix element is not isospin
              // scheme, we need to anticipate the swap on the
              // isospin-scheme source matrix element from which we
              // are calculating it.

              double canonicalization_phase = +1;
              if (eta1p > eta2p)
                {
                  canonicalization_phase *= ParitySign(1+gbp+ConjugationGrade(xbp)+Sbp+Tbp);
                  std::swap(eta1p,eta2p);
                }
              if (eta1 > eta2)
                {
                  canonicalization_phase *= ParitySign(1+gb+ConjugationGrade(xb)+Sb+Tb);
                  std::swap(eta1,eta2);
                }

              // short circuit -- check target rme labels for match with unit tensor source rme
              bool eta_matches = (
                  (std::pair<int,int>(eta1p,eta2p) == std::pair<int,int>(eta1bp,eta2bp))
                  && (std::pair<int,int>(eta1,eta2) == std::pair<int,int>(eta1b,eta2b))
                );
              if (!eta_matches)
                continue;

              // calculate rme

              // calculate isospin Clebsch-Gordan coefficient
              //
              // adapted from shell project moshinsky_xform.cpp

              // static const double kPPCoefficients[] = {+1,+sqrt(1/2.),+sqrt(1/10.)};
              // static const double kNNCoefficients[] = {+1,-sqrt(1/2.),+sqrt(1/10.)};
              static const double kPNCoefficients11[] = {+1,0,-sqrt(2./5.)};
              static const double kPNCoefficient10 = 1.;
              static const double kPNCoefficient01 = -sqrt(1/3.);
              static const double kPNCoefficient00 = 1.;

              double isospin_coefficient;
              int Tp = int(Tbp);
              int T = int(Tb);
              if ((Tp==1)&&(T==1))
                isospin_coefficient = kPNCoefficients11[int(T0)];
              else if ((Tp==0)&&(T==1))
                isospin_coefficient = kPNCoefficient01;
              else if ((Tp==1)&&(T==0))
                isospin_coefficient = kPNCoefficient10;
              else if ((Tp==0)&&(T==0))
                isospin_coefficient = kPNCoefficient00;

              // calculate antisymmetry normalization factors
              //
              // Includes overall factor of 1/2, as in shell project
              // moshinsky_xform.cpp.

              double pn_normalization_factor = 0.5;
              if (eta1p==eta2p)
                pn_normalization_factor *= sqrt(2.);
              if (eta1==eta2)
                pn_normalization_factor *= sqrt(2.);

              // combine factors
              //
              // RME relations are given on "SU(3) omegaST <->
              // omegaSpn scheme" notes page 4
              double rme = pn_normalization_factor * canonicalization_phase*isospin_coefficient;
              
              // save matrix element
              matrices[sector_index](bra_index,ket_index) = rme;
            }
      }

  }

////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  // process arguments
  if(argc<5)
    {
      std::cout<<"Syntax: Protons Neutrons Nmin Nmax 2Nsigma0"<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  int Z=std::stoi(argv[1]);
  assert(Z==1);
  int N=std::stoi(argv[2]);
  assert(N==1);
  int Nmax=std::stoi(argv[3]);
  // will be either 1 or 2; 
  int Nstep=std::stoi(argv[4]);
  assert((1<=Nstep)&&(Nstep<=2));
  int Nmin=Nmax%Nstep;
  int two_Nsigma_0=std::stoi(argv[5]);
  assert(two_Nsigma_0==2*3);  // deuteron Nsigma_0 = 2*(3/2) = 3
  HalfInt Nsigma_0 = HalfInt(two_Nsigma_0,2);

  // generate list of relative unit tensors up to Nmax cutoff
  std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST> two_body_unit_tensor_labels_list;
  u3shell::GenerateTwoBodyUnitTensorLabelsU3ST(Nmax, two_body_unit_tensor_labels_list);

  // read lsu3shell basis table and construct basis mapping
  std::string lsu3shell_basis_filename("lsu3shell_basis.dat");
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3shell_basis_filename,basis_table,basis_provenance,space);
  
  // iterate over unit tensors
  int num_unit = two_body_unit_tensor_labels_list.size();
  for(int operator_index=0; operator_index<num_unit; ++operator_index)
    {
      // extract tensor information
      const u3shell::TwoBodyUnitTensorLabelsU3ST& two_body_unit_tensor_labels
        = two_body_unit_tensor_labels_list[operator_index];
      std::cout
        << fmt::format("tensor {} labels {}",operator_index,two_body_unit_tensor_labels.Str())
        << std::endl;

      // read computed rmes
      std::string rme_filename = fmt::format("two_body_unit_{:06d}",operator_index);
      std::cout << fmt::format("reading {}",rme_filename) << std::endl;
      u3shell::SectorsU3SPN sectors(space,u3shell::OperatorLabelsU3S(two_body_unit_tensor_labels),false);
      basis::MatrixVector matrices_input;
      std::ifstream rme_stream(rme_filename);
      lsu3shell::ReadLSU3ShellRMEs(
          rme_stream,
          u3shell::OperatorLabelsU3ST(two_body_unit_tensor_labels),basis_table,
          space,sectors,matrices_input
      );

      // calculate reference RMEs
      u3shell::SectorsU3SPN sectors_dummy;  // duplicates sectors established above
      basis::MatrixVector matrices_reference;
      GenerateTwoBodyUnitTensorMatrices(      
          two_body_unit_tensor_labels,basis_provenance,
          space,sectors,matrices_reference
        );

      // compare RMEs
      double epsilon = 1e-6;
      lsu3shell::CompareLSU3ShellRMEs(
          std::cout,
          basis_provenance,
          space,sectors,matrices_input,matrices_reference,
          epsilon,true
        );
    }

}

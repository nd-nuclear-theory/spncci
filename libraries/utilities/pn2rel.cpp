/****************************************************************
  pn2rel.cpp

  Populate Hamiltonian-like operator from Petr Navratil ("PN") format
  relative operator file.

  See lsjt_operator.h for documentation of operator storage and the
  relative operator file format.

  commandline arguments: Nmax Jmax source_filename output_filename

  Language: C++11

  Anna E. McCoy
  TRIUMF

  5/2/18 (aem): Created, based upon jpv2rel.cpp.
****************************************************************/
#include <fstream>

#include "basis/lsjt_operator.h"
#include "basis/jjjpn_scheme.h"  // for TwoBodySpecies enum typedef
#include "cppformat/format.h"
#include "mcutils/eigen.h"  // for debugging output
#include "mcutils/parsing.h"


namespace utils
{

  void ReadPNOperator(
      const std::string& source_filename,
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::OperatorLabelsJT& operator_labels,
      const basis::RelativeSectorsLSJT& sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      bool verbose
    )
  {
    std::cout<<"begin"<<std::endl;
    // validate operator labels
    assert(
        (operator_labels.J0==0) && (operator_labels.g0==0)
        && (operator_labels.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian)
      );

    // open stream for reading
    std::cout
      << "Reading relative operator file (Petr Navratil format)..." << std::endl
      << "  Filename: " << source_filename << std::endl;
    std::ifstream is(source_filename);
    // StreamCheck(bool(is),source_filename,"Failed to open relative operator file");

    // read source file
    std::string line;
    int line_count = 0;
    int np,Lp,n,L,S,Tz,J;
    double me;
    // int identifier;


    while (std::getline(is,line))
      {
        ++line_count;
        std::istringstream line_stream(line);
        
        // read rest of header line
        line_stream >> np >> Lp >> n >> L >> S >> J >> Tz >> me;
        // ParsingCheck(line_stream,line_count,line);

        if(Lp>L)
          {
            std::swap(np,n);
            std::swap(Lp,L);
          }


        // deduce implied sectors labels
        int T = (L+S+1)%2; // isospin forced by L+S+T~1
        int Tp = (Lp+S+1)%2;
        int g = L%2;  // parity forced by L~g
        int gp = Lp%2;

        std::cout<<fmt::format("{} {}  {} {}   {} {} {}    {}",np,Lp,n,L,S,J,Tz,me)<<std::endl;
        // look up corresponding sector in our internal representation
        int subspace_index_bra = relative_space.LookUpSubspaceIndex(
            basis::RelativeSubspaceLSJTLabels(Lp,S,J,Tp,gp)
          );

        int subspace_index_ket = relative_space.LookUpSubspaceIndex(
            basis::RelativeSubspaceLSJTLabels(L,S,J,T,g)
          );
        
        // std::cout<<"subspce indices "<< subspace_index_bra<<"  "<<subspace_index_ket<<std::endl;

        // short circuit if subspace falls outside our target truncation
        if ((subspace_index_bra==basis::kNone)||(subspace_index_ket==basis::kNone))
          {
            std::cout << "ERROR: Input sector contains LSJT subspace not present in target truncation" << std::endl;
            std::exit(EXIT_FAILURE);
          }

        // set up references for convenience
        // const basis::RelativeSectorsLSJT& sectors = relative_component_sectors[Tz+1];
        basis::OperatorBlocks<double>& matrices = relative_component_matrices[Tz+1];

        
        int sector_index = sectors.LookUpSectorIndex(subspace_index_bra,subspace_index_ket);
        // std::cout<<"get sector "<<sector_index<<std::endl;
        const basis::RelativeSectorsLSJT::SectorType& sector = sectors.GetSector(sector_index);

        // look up target matrix dimensions
        int dimension_bra = sector.bra_subspace().size();
        int dimension_ket = sector.ket_subspace().size();
        // std::cout<<"dimensions "<<dimension_bra<<"  "<<dimension_ket<<std::endl;
        // std::cout<<matrices[sector_index].rows()<<"  "<<matrices[sector_index].cols()<<std::endl;
        // std::cout<<np<<"  "<<n<<std::endl;
        // save matrix element
        // std::cout<<"saving for sector "<<sector_index<<std::endl;
        if ((np<dimension_bra)&&(n<dimension_ket))
          matrices[sector_index](np,n) = me; 
        
      }

  }
  ////////////////////////////////////////////////////////////////

  void ReadPNOperatorPN(
      const std::string& source_filename,
      const basis::RelativeSpaceLSJT& relative_space,
      const basis::RelativeOperatorParametersLSJT& operator_parameters,
      const std::array<basis::RelativeSectorsLSJT,3>& relative_component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& relative_component_matrices,
      bool verbose
    )
  {

    // validate operator labels
    assert(
        (operator_parameters.J0==0) && (operator_parameters.g0==0)
        && (operator_parameters.T0_min==0) && (operator_parameters.T0_max==2)
        && (operator_parameters.symmetry_phase_mode==basis::SymmetryPhaseMode::kHermitian)
      );

    // Relation of JT-reduced matrix elements to pp/nn/pn matrix elements
    //
    //   < T || A^{T0} || T >  vs.  < TTz | A | TTz>
    //
    // For T=0 sectors:
    //
    //   <0||A0||0> = <00|A|00>
    //
    // For T=1 sectors:
    // 
    // {<1||A0||1>, <1||A1||1>, <1||A2||1>}
    // = 1/3. * {
    //           {1,1,1},
    //           {sqrt(9./2.),-sqrt(9./2.),0},
    //           {sqrt(5./2.),sqrt(5./2.),-sqrt(10.)}
    //         }
    //   * {<1+1|A|1+1>, <1-1|A|1-1>, <10|A|10>}
    //
    // We have listed Tz sectors in the order pp/pn/nn to match
    // Tz=-1,0,1 ordering.
    
    // transformation matrix
    //
    // indexed by (T,Tz)
    static Eigen::Matrix3d kIsospinCoefficientMatrixTzToTForT1;
    kIsospinCoefficientMatrixTzToTForT1
      << 1, 1, 1,
      std::sqrt(9./2.),  0, -std::sqrt(9./2.),
      std::sqrt(5./2.), -std::sqrt(10.), std::sqrt(5./2.) ;
    kIsospinCoefficientMatrixTzToTForT1 *= 1/3.;

    // set up storage for input matrix elements
    std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors_input;
    std::array<basis::OperatorBlocks<double>,3> relative_component_matrices_input;
    basis::ConstructZeroOperatorRelativeLSJT(
        basis::RelativeOperatorParametersLSJT(operator_parameters,operator_parameters.Nmax,operator_parameters.Jmax),
        relative_space,relative_component_sectors_input,relative_component_matrices_input
      );

    //Hack to get correct zero initialized matrices and sectors 
    // PN matrix elements all treated as T0=0 sectors 
    const basis::RelativeSectorsLSJT& input_sectors = relative_component_sectors[0];
    relative_component_matrices_input[1]=relative_component_matrices_input[0];
    relative_component_matrices_input[2]=relative_component_matrices_input[0];
  
    //For proton-neutron RMEs 'pretend' operator is isoscalar
    //The three components of relative_component_matrices_input will be 0: Tz=-1, 1: Tz=0 and 2: Tz=1
    const basis::OperatorLabelsJT operator_labels_isoscalar(0,0,0,0,basis::SymmetryPhaseMode::kHermitian);
    utils::ReadPNOperator(
        source_filename,
        relative_space,
        operator_labels_isoscalar,
        input_sectors,
        relative_component_matrices_input,
        verbose
      );

    std::cout<<"finished reading in matrix elements"<<std::endl;

    // accumulate matrix elements
    const int num_sectors = relative_component_sectors_input[0].size();
    for (int sector_index=0; sector_index < num_sectors; ++sector_index)
      // for each source "isoscalar operator" sector
      {
        // set up aliases for convenience
        const typename basis::RelativeSectorsLSJT::SectorType& 
          input_sector = input_sectors.GetSector(sector_index);

        // extract sector isospin labels
          // std::cout<<"isospin labels"<<std::endl;
        int bra_T = input_sector.bra_subspace().T();
        int ket_T = input_sector.ket_subspace().T();
        if (bra_T!=ket_T)
          continue;  // short circuit known vanishing T-changing sectors
        int T = ket_T;
        // std::cout<<"for Tz"<<std::endl;
        for(int Tz=-1; Tz<=1; ++Tz)
          {
            const Eigen::MatrixXd& 
              input_matrix = relative_component_matrices_input[Tz+1][sector_index];
          
            
            if (T==0)
              // sector with (T'T)=(0,0)
              {
                // std::cout<<"T=0"<<std::endl;
                // T=0 sectors only relevant in Tz=0 file; they are zeroed out in Tz!=0 files
                if(Tz==0)
                  {
                    // Simple copy to corresponding target sector.  Though we
                    // can write this as an accumulation for consistency.
                    int T0=0;
                    relative_component_matrices[T0][sector_index] += input_matrix;
                  }
              }
            else if (T==1)
              // sector with (T'T)=(1,1)
              {
                // std::cout<<"T=1"<<std::endl;
                for (int T0=0; T0<=2; ++T0)
                  {
                    // look up target sector
                    int target_sector_index 
                      = relative_component_sectors[T0].LookUpSectorIndex(input_sector.bra_subspace_index(),input_sector.ket_subspace_index());
                    assert(sector_index!=basis::kNone);
                    const typename basis::RelativeSectorsLSJT::SectorType& target_sector
                      = relative_component_sectors[T0].GetSector(target_sector_index);

                    // accumulate matrix for sector
                    Eigen::MatrixXd& matrix
                      = relative_component_matrices[T0][target_sector_index];
                    // std::cout<<"helo"<<std::endl;
                    double isospin_coefficient = kIsospinCoefficientMatrixTzToTForT1(T0,Tz+1);
                    // std::cout<<matrix<<std::endl;
                    // std::cout<<isospin_coefficient<<std::endl;
                    // std::cout<<input_matrix<<std::endl;
                    matrix += isospin_coefficient * input_matrix;

                    // debugging output
                    // std::cout
                    //   << fmt::format(
                    //       "input sector Tz={:+d} {:s}x{:s} -> target sector T0={:d} {:s}x{:s}   coefficient {:.4f}",
                    //       Tz,
                    //       input_sector.bra_subspace().LabelStr(),input_sector.ket_subspace().LabelStr(),
                    //       T0,
                    //       target_sector.bra_subspace().LabelStr(),target_sector.ket_subspace().LabelStr(),
                    //       isospin_coefficient
                    //     )
                    //   << std::endl;
                    // std::cout
                    //   << mcutils::FormatMatrix(input_matrix,".7f","   ")
                    //   << std::endl;
                    // std::cout << "  ... accumulating target ... " << std::endl;
                    // std::cout
                    //   << mcutils::FormatMatrix(matrix,".7f","   ")
                    //   << std::endl;
                  }
              }
          }
      }
  }

  ////////////////////////////////////////////////////////////////
} // namespace



////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // header
  if(argc<3)
    {
      std::cout<<"Syntax: Nmax Jmax source_filename output_filename"<<std::endl;
      std::exit(EXIT_FAILURE);
    }

  int Nmax=std::stoi(argv[1]);
  int Jmax=std::stoi(argv[2]);
  std::string source_filename=argv[3];
  std::string output_filename=argv[4];

  std::cout << std::endl;
  std::cout << "np2rel -- PN to relative file conversion" << std::endl;
  std::cout << std::endl;

  // set up zero operator
  std::cout << "Operator setup..." << std::endl;
  basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);
  basis::OperatorLabelsJT operator_labels(0,0,0,2,basis::SymmetryPhaseMode::kHermitian);
  basis::RelativeOperatorParametersLSJT operator_parameters(operator_labels,Nmax,Jmax);
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> relative_component_matrices;
  basis::ConstructZeroOperatorRelativeLSJT(
      operator_parameters,relative_space,relative_component_sectors,relative_component_matrices
    );

  // // operator diagnostics
  // std::cout << "  Truncation:"
  //           << " Nmax " << parameters.Nmax
  //           << " Jmax " << parameters.Jmax
  //           << std::endl;
  // std::cout << "  Matrix elements:";
  // for (int T0=0; T0<=2; ++T0)
  //   std::cout << " " << basis::UpperTriangularEntries(relative_component_sectors[T0]);
  // std::cout << std::endl;
  // std::cout << "  Allocated:";
  // for (int T0=operator_labels.T0_min; T0<=operator_labels.T0_max; ++T0)
  //   std::cout << " " << basis::AllocatedEntries(relative_component_matrices[T0]);
  // std::cout << std::endl;


  // Read in Petr Navratil formatted matrix elements 
  utils::ReadPNOperatorPN(
    source_filename,relative_space,operator_parameters,
    relative_component_sectors,relative_component_matrices,
    false  // verbose
    );
  
  // write operator
  basis::WriteRelativeOperatorLSJT(
      output_filename,
      relative_space,
      operator_labels,
      relative_component_sectors,
      relative_component_matrices,
      true  // verbose
    );

  // termination
  return 0;
}

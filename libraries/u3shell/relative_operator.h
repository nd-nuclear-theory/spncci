/****************************************************************
  relative_operator.h

  Relative operator representation and enumeration.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  5/25/16 (mac): Created with code from two_body_operator and
                 tensor_labels.
  12/6/16 (aem): Added optional parameters to GenerateRelativeUnitTensors
  9/13/19 (aem): Rename relative operators to be consistent with other libraries/programs
   5/3/22 (aem): Add RelativeOperator class and supporting functions
****************************************************************/
#ifndef RELATIVE_OPERATOR_H_
#define RELATIVE_OPERATOR_H_

#include "sp3rlib/sp3r.h"
#include "utilities/utilities.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/u3st_scheme.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/operator_indexing_spin.h"
#include "u3shell/operator_indexing_spatial.h"
#include "u3shell/operator_indexing_sectors.h"
#include "u3shell/relative_upcoupling.h"

namespace u3shell
{

  ////////////////////////////////////////////////////////////////////////
  // coefficient container for expansion in terms of relative unit tensors
  ////////////////////////////////////////////////////////////////////////

  // typedef
  //   std::map<u3shell::RelativeUnitTensorLabelsU3ST,double>
  //   RelativeUnitTensorCoefficientsU3ST;

  typedef
    std::unordered_map<
      u3shell::RelativeUnitTensorLabelsU3ST,
      double,
      boost::hash<u3shell::RelativeUnitTensorLabelsU3ST>
    > RelativeUnitTensorCoefficientsU3ST;


  ////////////////////////////////////////////////////////////////
  // generation of unit tensor label lists
  ///////////////////////////////////////////////////////////////

  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax,
        int N1v,
        std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels,
        int J0=-1,
        int T00=-1,
        bool restrict_positive_N0=false
        );

  /// Generate labels for U3ST-scheme relative tensors acting within
  /// the relative space of a given Nmax truncation, J0 and T0 and stores
  /// label in vector
  ///
  ///  Tensors in vector are ordered by:
  ///    N0, Sp,Tp,S,T,S0,T0,etap (eta=etap-N0)
  ///
  ///  Tensors are subject to trianglarity constraints on
  ///    (Sp,S0,S) and (Tp,T0,T)
  ///
  ///  and parity constraint on bra and ket
  ///    (eta+S+T)~1 and (etap+Sp+Tp)~1
  ///
  ///  Arguments:
  ///    space (input) : space on which the unit tensors are represented
  ///    relative_unit_tensor_labels (output) : container for unit tensor labels grouped by N0
  ///    J0 (optional) : Angular momentum of operator. If default value of -1, operators are not
  ///                    restricted to those which can branch to J0.
  ///    T0 (optional) : Isospin componenet of operator.  If default value of -1, then T0 takes
  ///                    all values that satisfy triangularity constraint.
  ///    bool (optional) : if true restricts, N0 or operator to positive values only


  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax,
        int N1v,
        std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels,
        int J0=-1,
        int T00=-1,
        bool restrict_positive_N0=false
      );
  // Overload of above function where container is map with unit tensors sorted by N0.

//******************************** Added by J.H. *************************************
  void GenerateOneBodyUnitTensorLabelsU3S(
        int Nmax,
        int N1vp,
	int N1vn,
        std::vector<OneBodyUnitTensorLabelsU3S>& one_body_unit_tensor_labels
        );

  /// Generates labels for one-body unit tensors acting within
  /// the space of a given Nmax truncation and stores
  /// labels in vector
  ///
  ///  Arguments:
  ///    Nmax (input) : Nmax
  ///    N1vp,N1vn (input) : N of valence shell for protons and neutrons
  ///    one_body_unit_tensor_labels (output) : container for one-body unit tensor labels

  void GenerateTwoBodyDensityLabels(
        int Nmax,
        int N1vp,
        int N1vn,
        std::vector<TwoBodyDensityLabels>& two_body_density_labels
        );

  /// Generates labels for two-body densities acting within
  /// the space of a given Nmax truncation and stores
  /// labels in vector
  ///
  ///  Arguments:
  ///    Nmax (input) : Nmax
  ///    N1vp,N1vn (input) : N of valence shell for protons and neutrons
  ///    two_body_density_labels (output) : container for two-body density labels
//************************************************************************************

  double Nrel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);
  // U(3) RME of relative number operator between relative harmonic oscillator states (bra and ket)

  double Arel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);
  // U(3) RME of Sp(3,R) raising operator between relative harmonic oscillator states (bra and ket)

  double Brel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);
  // U(3) RME of Sp(3,R) lowering operator between relative harmonic oscillator states (bra and ket)

  double Crel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);
  // U(3) RME of tensor of relative U(3) generators between relative harmonic oscillator states (bra and ket)

  double K2rel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);
  // U(3) RME of k^2 between relative harmonic oscillator states (bra and ket)

  double Qrel(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);
  // U(3) RME of mass quadrupole operator between relative harmonic oscillator states (bra and ket)

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// New code for updated recurrence
////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace spatial::onecoord
{
  /// SU(3) RME of the Sp(3,R) raising operator in a single coordiante
  /// harmonic oscillator basis.
  ///
  /// Args:
  ///   Nbarp : Harmonic oscillator quanta of bra
  ///   Nbar  : Harmonic oscillator quanta of ket
  ///
  /// Returns:
  ///   SU(3) reduced rme.
  inline double Arme(const unsigned int Nbarp, const unsigned int Nbar)
  {
    if(Nbarp==Nbar+2)
      return std::sqrt((Nbar+2)*(Nbar+1)/2);
    else
      return 0.0;
  }

  /// SU(3) RME of the Sp(3,R) lowering operator in a single coordiante
  /// harmonic oscillator basis.
  ///
  /// Args:
  ///   Nbarp : Harmonic oscillator quanta of bra
  ///   Nbar  : Harmonic oscillator quanta of ket
  ///
  /// Returns:
  ///   SU(3) reduced rme.
  inline double Brme(const unsigned int Nbarp, const unsigned int Nbar)
    {
      if(Nbarp+2==Nbar)
        return std::sqrt((Nbar+2)*(Nbar+1)/2);
      else
        return 0.0;
    }

  /// SU(3) RME of the SU(3) generators in a single coordiante
  /// harmonic oscillator basis.
  ///
  /// Args:
  ///   Nbarp : Harmonic oscillator quanta of bra
  ///   Nbar  : Harmonic oscillator quanta of ket
  ///
  /// Returns:
  ///   SU(3) reduced rme.
  inline double Crme(const unsigned int Nbarp, const unsigned int Nbar)
    {
      if(Nbarp==Nbar)
        return std::sqrt(4.*(Nbar*Nbar+3*Nbar)/3.);
      else
        return 0.0;
    }


  /// SU(3) RME of the harmonic oscillator Hamiltonian in a single
  /// coordiante harmonic oscillator basis.
  ///
  /// Args:
  ///   Nbarp : Harmonic oscillator quanta of bra
  ///   Nbar  : Harmonic oscillator quanta of ket
  ///
  /// Returns:
  ///   SU(3) reduced rme.
  inline double Hrme(const unsigned int Nbarp, const unsigned int Nbar)
    {
      if(Nbarp==Nbar)
        return Nbar+1.5;
      else
        return 0.0;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// New code part of spncci rewrite
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace relative
{
  using StateLabelsNST = std::array<unsigned int,3>;
  using TensorLabelsU3ST = std::tuple<u3::SU3,uint8_t,uint8_t>;
  using RMEFunction = std::function<double(TensorLabelsU3ST,StateLabelsNST,StateLabelsNST,unsigned int, unsigned int)>;

  std::vector<double> RelativeOperatorRMEs(
      const OperatorSectors& sectors,
      const RMEFunction& rme_function,
      const double coef=1.0
    );
  /// Function what computes all of the SU(3)xSU(2)xSU(2) reduced RMEs for sectors using
  /// the function rme_function.
  ///
  /// Args:
  ///   sectors: Sectors for operator in a relative harmonic oscillator basis
  ///   rme_function: Function for computing RMEs of operator in a relative harmonic oscillator
  ///     basis.

  class RelativeOperator
  {

  public:
    RelativeOperator()=default;
    RelativeOperator(const OperatorParameters& parameters,
      const OperatorSectors& sectors,
    const std::vector<double>& rmes
    )
    : parameters_{parameters}, sectors_{sectors},rmes_{rmes}
    {}

    // Simple constructor with one parameter and one rme function
    RelativeOperator(
      const OperatorParameters& parameters,
      const RMEFunction& rme_function,
      const double coef=1.0
    )
    : parameters_{parameters}
    {
       // Generate sectors
      sectors_ = ConstructOperatorSectors(parameters_);

      // Compute RMEs
      rmes_ =  RelativeOperatorRMEs(sectors_,rme_function);
    }

    /// Constructor that can take combines multiple operators and parameter sets
    ///   and combine them into a single operator
    /// Args:
    ///   parameter_sets : vector of parameters consisting the terms of the operators
    ///        e.g., the parameter set for H = T + V would consist of
    ///           {KineticEnergyParameters(Nbar_max), HamiltonianParameters(Nbar_max)}.
    ///
    ///       The combined set of parameters will be used to define the sectors whose rme's
    ///       will be calculated.
    ///
    ///   rme_functions : functions which can be used to compute the rmes of the operator.
    ///       The rme's for all sectors will be computed for each operator.
    ///
    ///   coefficients : coefficients used for combining the different functions. Length of coefficient
    ///       vector must be the same as the length if the rme_functions vector.
    ///
    ///   The operator defined by,
    ///     rme_functions={A,B}
    ///     coefficients = {alpha,bet}
    ///
    ///   is given by
    ///       C =alpha*A + beta*B
    RelativeOperator(
      const std::vector<OperatorParameters>& parameter_set,
      const std::vector<RMEFunction>& rme_functions,
      const std::vector<double>& coefficients
    );

    /// Construct operator from file
    RelativeOperator(const std::string& filename, const std::string& filetype="text");

    /// Constructor that reads in lsjt branched matrix elements and upcouples
    /// the matrix elements to obtain U(3)xSU(2)xSU(2) reduced matrix elements.
    RelativeOperator(
      const unsigned int Nbar_max,
      const OperatorParameters& parameters,
      const std::string& input_filename,
      const unsigned int input_Nbar_max,
      const unsigned int input_Jmax
    ): parameters_{parameters}
    {

      sectors_ = ConstructOperatorSectors(parameters_);

      rmes_ = u3shell::relative::UpcoupleU3ST(
          Nbar_max, sectors_,input_filename, input_Nbar_max, input_Jmax
        );

    }

    // Accessors
    inline const OperatorParameters& parameters() const { return parameters_; }
    inline const OperatorSectors& sectors() const {return sectors_;}
    inline const std::vector<double>& rmes() const {return rmes_;}

    double ReducedMatrixElement(
      const TensorLabelsU3ST& tensor_labels,
      const StateLabelsNST& bra,
      const StateLabelsNST& ket,
      const unsigned int kappa0=1,
      const unsigned int L0 = u3shell::relative::kNone
    );//NOTE:  Not currently used except in testing

  private:
    OperatorParameters parameters_;
    OperatorSectors sectors_;
    std::vector<double> rmes_;
  };

  /// Write Relative operator to file
  /// Writes operator parameters and rmes
  void WriteRelativeOperatorText(const RelativeOperator& op, const std::string& filename);


  /// Functions of type RMEFunction for calculating RMEs for differnt operators
  /// in the relative harmonic oscillator basis.
  ///
  /// Args:
  ///   tensor_labels: array containing (x0,S0,T0) where x0=$(\lambda_0,\mu_0)$
  ///   bra: array containing N,S,T for the bra of the RME.
  ///   ket: array containing N,S,T for the ket of the RME.

  double IdentityRME(
    const TensorLabelsU3ST& tensor_labels,
    const StateLabelsNST& bra,
    const StateLabelsNST& ket,
    const unsigned int kappa0=1,
    const unsigned int L0 = 0
  );
  /// Computes SU(3) reduced matrix element of the Identity operator


  double QuadrupoleRME(
    const TensorLabelsU3ST& tensor_labels,
    const StateLabelsNST& bra,
    const StateLabelsNST& ket,
    const unsigned int kappa0=1,
    const unsigned int L0 = 2
  );
  /// Computes SU(3) reduced matrix element of the Quadrupole operator


  double KSquaredRME(
    const TensorLabelsU3ST& tensor_labels,
    const StateLabelsNST& bra,
    const StateLabelsNST& ket,
    const unsigned int kappa0=1,
    const unsigned int L0 = 0
  );
  /// Computes SU(3) reduced matrix element of the $k^2$ operator

  double RSquaredRME(
    const TensorLabelsU3ST& tensor_labels,
    const StateLabelsNST& bra,
    const StateLabelsNST& ket,
    const unsigned int kappa0=1,
    const unsigned int L0 = 0
  );
  /// Computes SU(3) reduced matrix element of the $r^2$ operator



}// relative
  
}  // namespace

#endif

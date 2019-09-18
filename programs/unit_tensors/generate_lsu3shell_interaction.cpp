/****************************************************************
  generate_lsu3shell_interaction.cpp

  For a given nucleus, Nmax and N_step, generates input files for LSU3Shell 
  RecoupleSU3Interaction for interaction decomposed in terms of 
  SU(3)xSU(2) PN biquads which recouples interaction into format
  required by LSU3Shell. Input files include 

  Control file (operators.dat) giving:
    basename N lambda mu 2*S

  Command line parameters:
    Nmax : max number of oscillator quant in relative basis
    Nstep : indicates if all or only same parity spaces are included
            in model space.  Can only take values of 1 or 2. 
    ...

  Anna E. McCoy and Mark A. Caprio
  TRIUMF and University of Notre Dame

  - 9/12/19 (aem): Created.
 
****************************************************************/

#include <fstream>
#include <ostream>  

#include <ctime> //for testing
#include <boost/algorithm/string.hpp>

#include "fmt/format.h"

#include "lsu3shell/lsu3shell_operator.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"

namespace u3shell
{
    //Copied form generate_relative_u3st_operators.cpp and removed A and coef dependence
    // TODO move functions from both programs into a library
    void Interaction(
      int Nmax, int Jmax, int J0, int T0, int g0, 
      std::string& interaction_filename,
      u3shell::RelativeRMEsU3ST& interaction_u3st
    )
    {     
      // Read in the interaction from file
      basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Jmax);
      std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
      std::array<basis::MatrixVector,3> isospin_component_matrices_lsjt;

      basis::RelativeOperatorParametersLSJT operator_labels;
      basis::ReadRelativeOperatorLSJT(
        interaction_filename,relative_space_lsjt,operator_labels,
        isospin_component_sectors_lsjt, isospin_component_matrices_lsjt, true
        );

      // upcouple interaction
      u3shell::Upcoupling(
        relative_space_lsjt,
        isospin_component_sectors_lsjt,
        isospin_component_matrices_lsjt,
        J0, g0, T0,Nmax, interaction_u3st);
    }
}

int main(int argc, char **argv)
{
  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  u3::U3CoefInit();

  // process arguments
  if(argc<4+1)
    {
      std::cout<<"Syntax: Nmax Jmax T0 <interaction filename>"<<std::endl;
      std::exit(EXIT_FAILURE);
    }

  int Nmax=std::stoi(argv[1]);  //Nmax of interaction
  int Jmax=std::stoi(argv[2]);  //Jmax of interaction
  int T0=std::stoi(argv[3]);    //If T0=-1, then T0=0,1,2
  std::string filename=argv[4]; //interaction filename

  //Interations are angular momentum scalars and parity conserving
  int J0=0;
  int g0=0;

  // Unnecessary for interaction but relevant for other opertators. 
  int Nstep=1;

  bool un_u3_restrict=false;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Get biquad expansion of relative unit tensors 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Generate all relative unit tensor labels (up to Nmax cutoff)
  // Artificually set N1v to zero since Nmax in this case already accounts for N1v

  int N1v=0;
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1v,relative_unit_tensor_labels,J0,T0,false);

  int num_unit_tensors = relative_unit_tensor_labels.size();
  std::cout
    << fmt::format("number of relative tensors: {}",num_unit_tensors)
    <<std::endl;


  // generate unit tensor operator files
  
  // lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, relative_unit_tensor_labels,un_u3_restrict);

  u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
  int i=0;

  // container for computed biquad expansion of unit tensor 
  std::vector<u3shell::TwoBodyUnitTensorCoefficientsU3SPN> 
    biquad_coefficients_pn_list(relative_unit_tensor_labels.size());

  
  double duration;
  std::cout<<"generating biquads for relative unit tensors"<<std::endl; 
  std::clock_t start_biquads;
  start_biquads=std::clock();
  #pragma omp parallel for schedule(runtime) 
  for(int i=0; i<relative_unit_tensor_labels.size(); ++i)
    {
      const u3shell::RelativeUnitTensorLabelsU3ST& tensor=relative_unit_tensor_labels[i];
      // declare coefficient containers
      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
      u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
      u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
      
      // Moshinsky transform unit tensors to two-body operators 
      MoshinskyTransformTensor(tensor, twobody_space, 
        two_body_unit_tensor_coefficients, "NAS");
      
      // convert to biquads
      u3shell::TransformTwoBodyUnitTensorToBiquad(two_body_unit_tensor_coefficients,biquad_coefficients);
      // convert biquads to pn scheme
      u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn,un_u3_restrict);
      
      biquad_coefficients_pn_list[i]=biquad_coefficients_pn;
    }
  duration = (std::clock() - start_biquads)/(double) CLOCKS_PER_SEC;
  std::cout<<fmt::format("finished generating biquads.  Time: {:.4f} ",duration)<<std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Contract with interaction    
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"Reading in interaction and upcoupling"<<std::endl;  
  u3shell::RelativeRMEsU3ST interaction_u3st;
  u3shell::Interaction(Nmax, Jmax, J0, T0, g0, filename,interaction_u3st);

  //Container for biquad coefficients of interaction.  Stored by kappa0, by rho0.
  u3shell::TwoBodyOperatorCoefficientsU3SPN interaction_biquad_coefficieints_pn;

  std::cout<<"Starting interaction unit tensor contraction"<<std::endl;
  std::clock_t start_contraction;
  start_contraction=std::clock();  

  //Iterate over relative unit tensors converted to biquads
  for(int i=0; i<relative_unit_tensor_labels.size(); ++i)
    {
      // Extract relative unit tensor labels and coefficients 
      u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels=relative_unit_tensor_labels[i];
      u3shell::TwoBodyUnitTensorCoefficientsU3SPN& biquad_coefficients_pn=biquad_coefficients_pn_list[i];

      // Reorganize biquads to correct ordering for recoupler
      u3shell::TwoBodyUnitTensorCoefficientsRecouplerU3SPN biquad_coefficients_pn_recoupler;
      u3shell::RegroupBiquadsForRecoupler(biquad_coefficients_pn,biquad_coefficients_pn_recoupler);

      // Get tensor labels from relative unit tensor labels to look up relative interaction rmes 
      const u3shell::OperatorLabelsU3ST& operator_labels=unit_tensor_labels.operator_labels();
      const u3::SU3& x0=operator_labels.x0();
      const HalfInt& S0=operator_labels.S0();
      
      // Since J0=0, L0=S0
      int L0(S0);
      int kappa0_max=u3::BranchingMultiplicitySO3(x0, L0);

      //For each kappa0 value, check if there is a non-zero interaction rme.  If so, contract with expanded unit tensor 
      for(int kappa0=1; kappa0<=kappa0_max; ++kappa0)
        {
          std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(unit_tensor_labels,kappa0,L0);
          if (interaction_u3st.count(key))
            {
              double rme=interaction_u3st[key];
              
              for(auto& key_value : biquad_coefficients_pn_recoupler)
                {
                  const u3shell::TwoBodyUnitTensorLabelsU3S& biquad_labels= key_value.first;
                  std::vector<u3shell::CoefficientsPN>& coefficients_vector=key_value.second;
                  
                  std::vector<std::vector<u3shell::CoefficientsPN>>& interaction_coefficients_array=interaction_biquad_coefficieints_pn[biquad_labels];
                  int rho0_max=coefficients_vector.size();
                  
                  //If coefficient array not allocated, allocate
                  if(interaction_coefficients_array.size()==0)
                    {
                      
                      interaction_coefficients_array.resize(kappa0_max);
                      for(int k=0; k<kappa0_max; ++k)
                        interaction_coefficients_array[k].resize(rho0_max);
                    } 
                  for(int r=0; r<rho0_max; ++r)
                    {
                      u3shell::CoefficientsPN& coefficients=coefficients_vector[r];
                      u3shell::CoefficientsPN new_coefficients(coefficients.pppp*rme, coefficients.nnnn*rme, coefficients.pnnp*rme);
                      interaction_coefficients_array[kappa0-1][r]+=new_coefficients;
                    }
                }
            }
        }
    }// end contraction 

  duration = (std::clock() - start_contraction)/(double) CLOCKS_PER_SEC;
  std::cout<<fmt::format("Finished contraction.  Time: {:.4f}",duration)<<std::endl;
  
  // Write biquads to file
  //Strips path to file from filename so only the actual filename remains
  std::vector<std::string> split_filename; 
  boost::split(split_filename, filename, boost::is_any_of("/")); 
  std::string output_filename=split_filename[split_filename.size()-1]+"_recoupler";

  std::cout<<"writing to file "<<output_filename<<std::endl;
  u3shell::WriteTwoBodyInteractionRecoupler(output_filename,interaction_biquad_coefficieints_pn);

}

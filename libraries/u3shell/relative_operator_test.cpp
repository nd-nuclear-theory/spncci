/****************************************************************
  relative_operator_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  5/25/16 (aem): Created.
****************************************************************/
#include <iostream>
#include <vector>
#include <memory>
#include "fmt/format.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/relative_operator.h"
#include "u3shell/upcoupling.h"
////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  u3::U3CoefInit(100);

  // Validating OperatorParameters
  if(true)
  {
    unsigned int Nmax=2;
    unsigned int N1v=1;
    unsigned int Nbar_max = Nmax+2*N1v;

    ////////////////////////////////////////////////////////////////
    // Checking predefined operators
    ////////////////////////////////////////////////////////////////
    // Identity
    auto identity_parameters = u3shell::relative::IdentityParameters(Nbar_max);
    assert(identity_parameters.Nbar_max==4);
    assert(identity_parameters.J0==0);
    assert(identity_parameters.Allowed_w0_values.size()==1);
    assert(*identity_parameters.Allowed_w0_values.begin()==u3::U3(0,{0u,0u}));
    assert(*identity_parameters.Allowed_L0_values.begin()==0);
    assert(*identity_parameters.Allowed_T0_values.begin()==0);
    assert(*identity_parameters.Allowed_S0_values.begin()==0);

    auto rsquared_parameters = u3shell::relative::RSquaredParameters(Nbar_max);
    assert(rsquared_parameters.J0==0);
    assert(rsquared_parameters.Allowed_w0_values.size()==3);
    assert(rsquared_parameters.Allowed_w0_values.count(u3::U3(-2,{0u,2u})));
    assert(rsquared_parameters.Allowed_w0_values.count(u3::U3(2,{2u,0u})));
    assert(rsquared_parameters.Allowed_w0_values.count(u3::U3(0,{0u,0u})));
    assert(*rsquared_parameters.Allowed_L0_values.begin()==0);
    assert(*rsquared_parameters.Allowed_S0_values.begin()==0);
    assert(*rsquared_parameters.Allowed_T0_values.begin()==0);


    auto quadrupoleT1_parameters = u3shell::relative::QIsovectorParameters(Nbar_max);
    assert(quadrupoleT1_parameters.Allowed_w0_values.size()==3);
    assert(quadrupoleT1_parameters.Allowed_w0_values.count(u3::U3(-2,{0u,2u})));
    assert(quadrupoleT1_parameters.Allowed_w0_values.count(u3::U3(2,{2u,0u})));
    assert(quadrupoleT1_parameters.Allowed_w0_values.count(u3::U3(0,{1u,1u})));
    assert(*quadrupoleT1_parameters.Allowed_L0_values.begin()==2);
    assert(*quadrupoleT1_parameters.Allowed_S0_values.begin()==0);
    assert(*quadrupoleT1_parameters.Allowed_T0_values.begin()==1);
    assert(quadrupoleT1_parameters.J0==2);


    //Combining parameters
    u3shell::relative::OperatorParameters combined_parameters =
        u3shell::relative::CombineParameters(
            {identity_parameters, rsquared_parameters, quadrupoleT1_parameters}
          );

    assert(combined_parameters.J0==u3shell::relative::kNone);
    assert(combined_parameters.Allowed_w0_values.size()==4);
    assert(combined_parameters.Allowed_w0_values.count(u3::U3(-2,{0u,2u})));
    assert(combined_parameters.Allowed_w0_values.count(u3::U3(2,{2u,0u})));
    assert(combined_parameters.Allowed_w0_values.count(u3::U3(0,{0u,0u})));
    assert(combined_parameters.Allowed_w0_values.count(u3::U3(0,{1u,1u})));
    assert(*combined_parameters.Allowed_S0_values.begin()==0);
    assert(combined_parameters.Allowed_L0_values.size()==2);
    assert(combined_parameters.Allowed_L0_values.count(0));
    assert(combined_parameters.Allowed_L0_values.count(2));
    assert(combined_parameters.Allowed_T0_values.size()==2);
    assert(combined_parameters.Allowed_T0_values.count(0));
    assert(combined_parameters.Allowed_T0_values.count(1));


    // Operator parameters for all possible unit tensors
    // Generate list of relative unit tensors for a given Nmax and N1v.
    std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
    u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1v,relative_unit_tensor_labels);
    // for(const auto& tensor : relative_unit_tensor_labels)
    //   fmt::print("{}\n",tensor.Str());

    std::set<unsigned int> Allowed_L0_values={u3shell::relative::kNone};
    std::set<uint8_t> Allowed_S0_values;
    std::set<uint8_t> Allowed_T0_values;
    std::unordered_set<u3::U3> Allowed_w0_values;

    for(const auto& tensor : relative_unit_tensor_labels)
      {
        const auto& [N0,x0,S0,T0,g0] = tensor.operator_labels().Key();
        Allowed_w0_values.insert({N0,x0});
        Allowed_S0_values.insert(int(S0)); //No direct conversion from HaltInt to uint8_t
        Allowed_T0_values.insert(int(T0));
      }

    u3shell::relative::OperatorParameters unit_tensor_parameters(
      Nmax+2*N1v,
      u3shell::relative::kNone,
      Allowed_w0_values,
      Allowed_L0_values,
      Allowed_S0_values,
      Allowed_T0_values
      );

    fmt::print("Nbar max {}\n",unit_tensor_parameters.Nbar_max);

    ////////////////////////////////////////////////////////////////
    // Checking construction of corresponding sectors
    ///////////////////////////////////////////////////////////////
    // Generate the corresponding spatial and spin operator spaces.
    u3shell::spatial::onecoord::OperatorSpace
      spatial_operator_space(unit_tensor_parameters);
    // std::cout<<spatial_operator_space.DebugStr()<<std::endl;

    u3shell::spin::twobody::OperatorSpace
      spin_operator_space(unit_tensor_parameters);
    // std::cout<<spin_operator_space.DebugStr()<<std::endl;

    // Generate corresponding sectors
    u3shell::relative::OperatorSectors operator_sectors(
      std::make_shared<u3shell::spatial::onecoord::OperatorSpace>(spatial_operator_space),
      std::make_shared<u3shell::spin::twobody::OperatorSpace>(spin_operator_space)
    );
    unsigned int num_tensors = 0;
    assert(operator_sectors.num_elements()==relative_unit_tensor_labels.size());
    // fmt::print("Num elements {}\n",operator_sectors.num_elements());
    // fmt::print("Num unit tensors {}\n",relative_unit_tensor_labels.size());
    // std::cout<<operator_sectors.DebugStr()<<std::endl;

  }



////////////////////////////////////////////////////////////////
// Checking relative operator construction
////////////////////////////////////////////////////////////////
if(true)
{
  fmt::print("Testing relative operator construction\n");
  unsigned int Nbar_max=4;

  ////////////////////////////////////////////////////////////////
  // Quadrupole operator
  ////////////////////////////////////////////////////////////////
  auto quadrupoleT0_parameters = u3shell::relative::QIsoscalarParameters(Nbar_max);
  auto quadrupoleT1_parameters = u3shell::relative::QIsovectorParameters(Nbar_max);

  auto sectors_T0 = u3shell::relative::OperatorSectors{
    std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(quadrupoleT0_parameters),
    std::make_shared<const u3shell::spin::twobody::OperatorSpace>(quadrupoleT0_parameters)
  };

  auto sectors_T1 = u3shell::relative::OperatorSectors{
    std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(quadrupoleT1_parameters),
    std::make_shared<const u3shell::spin::twobody::OperatorSpace>(quadrupoleT1_parameters)
  };

  std::vector<double> quadrupole_test0 = u3shell::relative::RelativeOperatorRMEs(sectors_T0, u3shell::relative::QuadrupoleRME);
  std::vector<double> quadrupole_test1 = u3shell::relative::RelativeOperatorRMEs(sectors_T1, u3shell::relative::QuadrupoleRME);

  u3shell::relative::RelativeOperator quadruple_operator_T0(
    quadrupoleT0_parameters,
    u3shell::relative::QuadrupoleRME
  );

  u3shell::relative::RelativeOperator quadruple_operator_T1(
    quadrupoleT1_parameters,
    u3shell::relative::QuadrupoleRME
  );

  // Check that the stored rmes are the same.
  for(int i=0; i<quadrupole_test0.size(); ++i)
  {
    assert(fabs(quadrupole_test0[i]-quadruple_operator_T0.rmes()[i])<1e-10);
  }

  for(int i=0; i<quadrupole_test1.size(); ++i)
  {
    assert(fabs(quadrupole_test1[i]-quadruple_operator_T1.rmes()[i])<1e-10);
  }
  fmt::print("isoscalar and isovector quadrupole operators validated\n");

  // Combining the isoscalar and isovector
  auto combined_parameters =
    u3shell::relative::CombineParameters({quadrupoleT0_parameters,quadrupoleT1_parameters});
  auto spatial_ptr
    = std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(combined_parameters);
  auto spin_ptr
    = std::make_shared<const u3shell::spin::twobody::OperatorSpace>(combined_parameters);
  auto sectors = u3shell::relative::OperatorSectors{
    std::make_shared<const u3shell::spatial::onecoord::OperatorSpace>(combined_parameters),
    std::make_shared<const u3shell::spin::twobody::OperatorSpace>(combined_parameters)
  };

  std::vector<double> quadrupole_test = u3shell::relative::RelativeOperatorRMEs(sectors, u3shell::relative::QuadrupoleRME);

  u3shell::relative::RelativeOperator quadruple_operator(
    combined_parameters,
    u3shell::relative::QuadrupoleRME
  );

  for(int i=0; i<quadrupole_test.size(); ++i)
  {
    assert(fabs(quadrupole_test[i]-quadruple_operator.rmes()[i])<1e-10);
  }

  u3shell::relative::RelativeOperator quadruple_operator2(
    {quadrupoleT0_parameters,quadrupoleT1_parameters},
    {u3shell::relative::QuadrupoleRME},
    {1.0}
  );

  // Compare rmes vs test function
  for(int i=0; i<quadrupole_test.size(); ++i)
  {
    assert(fabs(quadrupole_test[i]-quadruple_operator2.rmes()[i])<1e-10);
  }
  fmt::print("isoscalar+isovector quadrupole operator validated\n");


  ////////////////////////////////////////////////////////////////
  // Checking construction via upcoupling
  ////////////////////////////////////////////////////////////////

  // Parameters for input from file.  Nmax must be >= Nbar_max.
  std::string filename = "Daejeon16_Nmax40_hw15.0_rel.dat";
      std::string inputfile = fmt::format(
          "{}/spncci/data/relative_interactions/{}",
          utils::get_spncci_project_root_dir(),
          filename
        );

  int Nmax = 40;
  int Jmax=Nmax+1;
  u3shell::relative::RelativeOperator daejeon_operator(Nbar_max,u3shell::relative::HamiltonianParameters(Nbar_max),inputfile,Nmax,Jmax);


  bool verbose = true;
  assert(utils::FileExists(inputfile, verbose));
  basis::RelativeSpaceLSJT relative_space_lsjt(Nbar_max, Jmax);
  std::array<basis::RelativeSectorsLSJT, 3> sectors_lsjt;
  std::array<basis::OperatorBlocks<double>, 3> blocks_lsjt;
  basis::RelativeOperatorParametersLSJT op_labels_lsjt;

  basis::ReadRelativeOperatorLSJT(
      inputfile,
      relative_space_lsjt,
      op_labels_lsjt,
      sectors_lsjt,
      blocks_lsjt,
      true
    );

  // For testing
  u3shell::RelativeRMEsU3ST rme_map;
  u3shell::Upcoupling(relative_space_lsjt,sectors_lsjt, blocks_lsjt,0,0,-1,Nbar_max,rme_map);

  for(const auto& [key, rme1] : rme_map)
    {
      const auto&[labels,kappa0,L0] = key;
      const auto&[x0,S0,T0,Nbarp,Sbarp,Tbarp,Nbar,Sbar,Tbar]=labels.FlatKey();
      u3shell::relative::StateLabelsNST ket{int(Nbar),int(Sbar),int(Tbar)};
      u3shell::relative::StateLabelsNST bra{int(Nbarp),int(Sbarp),int(Tbarp)};
      double rme2 = daejeon_operator.ReducedMatrixElement({x0,int(S0),int(T0)},bra,ket,kappa0,L0);
      // fmt::print("{} {} {}   {} {} {}   {} {} {}  {} {}\n",x0,S0,T0,Nbarp,Sbarp,Tbarp,Nbar,Sbar,Tbar,kappa0,L0);
      // std::cout<<rme1<<" "<<rme2<<std::endl;
      assert(fabs(rme1-rme2)<1e-8);
    }

  std::string output_filename = "relative_operator_daejeon_test_Nmax04.dat";
  std::cout<<"write to file"<<std::endl;
  u3shell::relative::WriteRelativeOperatorText(daejeon_operator,output_filename);
  std::cout<<"read from file"<<std::endl;
  u3shell::relative::RelativeOperator daejeon_operator2(output_filename);
  const auto& rmes1 = daejeon_operator.rmes();
  const auto& rmes2 = daejeon_operator2.rmes();
  assert(rmes1.size()==rmes2.size());
  for(int i=0; i<rmes1.size(); ++i)
  {
    assert(fabs(rmes1[i]-rmes2[i])<1e-7);
  }

  ////////////////////////////////////////////////////////////////
  // Constructing Hamiltonian operator
  ////////////////////////////////////////////////////////////////
  auto hamiltonian_parameters
  = u3shell::relative::CombineParameters({
    u3shell::relative::KSquaredParameters(Nbar_max),
    u3shell::relative::HamiltonianParameters(Nbar_max,0u,0u)
    });

  auto hamiltonian_sectors = u3shell::relative::ConstructOperatorSectors(hamiltonian_parameters);

  std::vector<double> hamiltonian_rmes = u3shell::relative::UpcoupleU3ST(
      Nbar_max, hamiltonian_sectors, inputfile, Nmax, Jmax
    );

  hamiltonian_rmes =
      hamiltonian_rmes
      + u3shell::relative::RelativeOperatorRMEs(
          hamiltonian_sectors, u3shell::relative::KSquaredRME, 15.0 / 4
        );

  u3shell::relative::RelativeOperator hamiltonian_operator(hamiltonian_parameters,hamiltonian_sectors,hamiltonian_rmes);

}



  // termination
  return 0;
}

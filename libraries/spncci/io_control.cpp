/****************************************************************
  io_control.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/io_control.h"

#include "cppformat/format.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // reading lsu3shell RMEs
  ////////////////////////////////////////////////////////////////

  void
    ReadLSU3ShellSymplecticOperatorRMEs(
        const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
        const u3shell::SpaceU3SPN& lsu3shell_space, 
        const std::string& Brel_filename, u3shell::SectorsU3SPN& Bintr_sectors, basis::MatrixVector& Bintr_matrices,
        const std::string& Nrel_filename, u3shell::SectorsU3SPN& Nintr_sectors, basis::MatrixVector& Nintr_matrices,
        int A
      )
  {
    // read Brel => Bintr
    u3shell::OperatorLabelsU3ST Brel_labels(-2,u3::SU3(0,2),0,0,0);
    Bintr_sectors = u3shell::SectorsU3SPN(lsu3shell_space,Brel_labels,true);
    bool sp3r_generators=true;
    lsu3shell::ReadLSU3ShellRMEs(
      sp3r_generators,
        Brel_filename,
        lsu3shell_basis_table,lsu3shell_space,
        Brel_labels,Bintr_sectors,Bintr_matrices,
        2./A
      );
    
    // read Nrel => Nintr
    u3shell::OperatorLabelsU3ST Nrel_labels(0,u3::SU3(0,0),0,0,0);
    Nintr_sectors = u3shell::SectorsU3SPN(lsu3shell_space,Nrel_labels,true);
    lsu3shell::ReadLSU3ShellRMEs(
        sp3r_generators,
        Nrel_filename,
        lsu3shell_basis_table,lsu3shell_space,
        Nrel_labels,Nintr_sectors,Nintr_matrices,
        2./A
      );

  }


  void
    ReadLSU3ShellSymplecticRaisingOperatorRMEs(
        const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
        const u3shell::SpaceU3SPN& lsu3shell_space, 
        const std::string& Arel_filename, u3shell::SectorsU3SPN& Aintr_sectors, basis::MatrixVector& Aintr_matrices,
        int A
      )
  {    
    // read Arel => Aintr
    u3shell::OperatorLabelsU3ST Arel_labels(2,u3::SU3(2,0),0,0,0);
    Aintr_sectors = u3shell::SectorsU3SPN(lsu3shell_space,Arel_labels,true);
    bool sp3r_generators=true;
    lsu3shell::ReadLSU3ShellRMEs(
        sp3r_generators,
        Arel_filename,
        lsu3shell_basis_table,lsu3shell_space,
        Arel_labels,Aintr_sectors,Aintr_matrices,
        2./A
      );

  }


  void
  ReadRelativeObservables(
    int Nmax, int N1v, const std::vector<double>& hw_values,
    const std::string& observable_directory,const std::vector<std::string>& observable_filenames, 
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    std::vector<std::vector<u3shell::RelativeRMEsU3SSubspaces>>& observables_relative_rmes,
    std::vector<std::vector<u3shell::IndexedOperatorLabelsU3S>>& relative_observable_labels
    )
  {
  std::cout << "Read observable relative rmes..." << std::endl;
  int num_observables=observable_filenames.size();

  // Resizing container holding observable labels
  relative_observable_labels.resize(num_observables); 

  // Setting up array storing relative rmes, array indexed first by hbar_omega index, then by observerable index
  observables_relative_rmes.resize(hw_values.size());
  for(int h=0; h<hw_values.size(); ++h)
    observables_relative_rmes[h].resize(num_observables);

  // for each observable, loop over all hw values, read in rmes from file and generate set of symmetry labels
  for (int observable_index=0; observable_index<num_observables; ++observable_index)
    {      
      // temporary container for accumulating set of symmetry labels over hw values for each operator
      std::unordered_set<u3shell::IndexedOperatorLabelsU3S, boost::hash<u3shell::IndexedOperatorLabelsU3S>> symmetries_u3s;
      
      // for each value of hbar omega
      for(int h=0; h<hw_values.size(); ++h)
        {
          double hw=hw_values[h];

          std::string observable_filename
            =fmt::format("{}/{}_hw{:2.1f}_Nmax{:02d}_u3st.dat", 
              observable_directory,observable_filenames[observable_index],hw,Nmax
              );

          std::cout << fmt::format("  Reading {}...",observable_filename)<< std::endl;

          // Read in and store relative obserevable symmetries and rmes in array
          u3shell::RelativeRMEsU3SSubspaces& relative_rmes=observables_relative_rmes[h][observable_index];
          u3shell::ReadRelativeOperatorU3ST(Nmax, N1v,observable_filename,unit_tensor_space,relative_rmes);
          
          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // turn into function
          // // Print out observable rmes 
          // for(auto it=relative_rmes.begin(); it!=relative_rmes.end(); ++it)
          //   {
          //     int unit_tensor_subspace_index,kappa0,L0;
          //     std::tie(unit_tensor_subspace_index,kappa0,L0)=it->first;

          //     const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace
          //       =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
              
          //     const std::vector<double>& rmes=it->second;
              
          //     std::cout<<unit_tensor_subspace.LabelStr()<<"  "<<kappa0<<"  "<<L0<<std::endl;
          //     for(int unit_tensor_index=0; unit_tensor_index<unit_tensor_subspace.size(); ++unit_tensor_index)
          //         {
          //           int T0, S,T,Sp,Tp;
          //           std::tie(T0,Sp,Tp,S,T)=unit_tensor_subspace.GetStateLabels(unit_tensor_index);
          //           std::cout<<fmt::format("{}   {} {}  {} {}  {}",T0,Sp,Tp,S,T,rmes[unit_tensor_index])<<std::endl;
          //         }     
          //   }
          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
          
          // loop over relative observable labels and accumulate a set of labels for each observable over hbar_omega values
          for(auto it=relative_rmes.begin(); it!=relative_rmes.end(); ++it)
            {
              int unit_tensor_subspace_index, kappa0,L0,etap,eta;
              HalfInt S0;
              u3::SU3 x0;
              std::tie(unit_tensor_subspace_index,kappa0,L0)=it->first;
              std::tie(x0,S0,etap,eta)=unit_tensor_space.GetSubspace(unit_tensor_subspace_index).labels();
              symmetries_u3s.insert(u3shell::IndexedOperatorLabelsU3S(u3shell::OperatorLabelsU3S(etap-eta,x0,S0),kappa0,L0));
            }
        }

      // Transfer accumulated set of observable labels to array for external use
      for(auto tensor : symmetries_u3s)
        relative_observable_labels[observable_index].push_back(tensor);
    }
  }



}  // namespace

/****************************************************************
  io_control.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/io_control.h"
#include <fstream>
#include <omp.h>  
#include "fmt/format.h"

namespace spncci
{


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

  void ReadBlock(std::istream& in_stream, Eigen::MatrixXd& block)
    {
      
      //Read in number of rows and columns
      int rows,cols;
      mcutils::ReadBinary<int>(in_stream,rows);
      mcutils::ReadBinary<int>(in_stream,cols);

      // Read in RMEs and case to double matrix 
      //TODO Change to MatrixFloatType for spncci
      // std::cout<<rows<<"  "<<cols<<"  "<<lgi::binary_float_precision<<std::endl;
      if(lgi::binary_float_precision==4)
        {
          float buffer[rows*cols];
          in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
          block=Eigen::Map<Eigen::MatrixXf>(buffer,rows,cols).cast<double>();
        }
      else if (lgi::binary_float_precision==8)
        {
          double buffer[rows*cols];
          in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
          block=Eigen::Map<Eigen::MatrixXd>(buffer,rows,cols);
        }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void WriteHyperBlock(
    const basis::OperatorHyperblocks<double>& baby_spncci_observable_hyperblocks,
    const std::string& filename,
    int irrep_family_index_bra,int irrep_family_index_ket
    )
  {
    std::ios_base::openmode mode_argument = std::ios_base::out | std::ios::app | std::ios_base::binary;
    std::ofstream out_file;
    out_file.open(filename,mode_argument);

    if (!out_file)
      {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
      }

    int num_hyperblocks=baby_spncci_observable_hyperblocks.size(); 
    
    // Write irrep family indices and number of hyperblocks 
    mcutils::WriteBinary<int>(out_file,irrep_family_index_bra);
    mcutils::WriteBinary<int>(out_file,irrep_family_index_ket);
    mcutils::WriteBinary<int>(out_file,num_hyperblocks);

    // for each block, 
    for(int hypersector_index=0; hypersector_index<num_hyperblocks; ++hypersector_index)
      {
        const basis::OperatorBlock<double>& block
            =baby_spncci_observable_hyperblocks[hypersector_index][0];

        int num_rows=block.rows();
        int num_cols=block.cols();

        // write number of rows and columns
        mcutils::WriteBinary<int>(out_file,num_rows);
        mcutils::WriteBinary<int>(out_file,num_cols);
      
        int size=num_rows*num_cols;

        // write matrix.  Order is column major (Eigen default)
        if(lgi::binary_float_precision==4)
            {
              Eigen::MatrixXf buffer_matrix=block.cast<float>();
              out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);
              
            }  
          
        else if (lgi::binary_float_precision==8)
          {
            Eigen::MatrixXd buffer_matrix=block;
            out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);

            // if(irrep_family_index_bra==1 && irrep_family_index_ket==0)
            //   {
            //     std::ios_base::openmode mode_argument = std::ios_base::out | std::ios_base::binary;
            //     std::ofstream test_file;
            //     test_file.open("test",mode_argument);

            //     double buff[size];
            //     buff=*block.data();

            //     std::cout<<block<<std::endl<<"buffer block"<<std::endl<<buffer_matrix<<std::endl;
            //     std::cout<<"size "<<sizeof(buff)<<std::endl;
                
            //     test_file.write(reinterpret_cast<char*>(&buff),sizeof(buf));
            //     test_file.close();

            //     std::ifstream test_in;
            //     test_in.open("test",std::ios_base::in | std::ios_base::binary);

            //     double buffer[size];
            //     std::cout<<"size2 "<<sizeof(buffer)<<std::endl;
                
            //     Eigen::MatrixXd block2;
            //     test_in.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
            //     block2=Eigen::Map<Eigen::MatrixXd>(buffer,num_rows,num_cols);
            //     std::cout<<"block2"<<std::endl<<block2<<std::endl;
            //     test_in.close();
            //   }

          }


      }
    out_file.close();
  
    // std::ifstream file(filename,std::ios::binary | std::ios::ate);
    // std::cout<<"file size "<<file.tellg()<<std::endl;
    // file.close();
  }


void ReadObservableHyperblocks(
  // int observable_index, int hw_index,
  std::istream& in_stream,
  spncci::LGIPair& lgi_pair,
  basis::OperatorHyperblocks<double>& baby_spncci_observable_hyperblocks  
  )

{
  
  // Read lgi family 
  int irrep_family_index_bra, irrep_family_index_ket, num_hyperblocks;
  mcutils::ReadBinary<int>(in_stream,irrep_family_index_bra);
  mcutils::ReadBinary<int>(in_stream,irrep_family_index_ket);

  assert(lgi_pair.first==irrep_family_index_bra);
  assert(lgi_pair.second==irrep_family_index_ket);
  lgi_pair=spncci::LGIPair(irrep_family_index_bra,irrep_family_index_ket);

  // std::cout<<irrep_family_index_bra<<"  "<<irrep_family_index_ket<<"  "<<omp_get_thread_num()<<std::endl;
  // Read number of hyperblocks 
  mcutils::ReadBinary<int>(in_stream,num_hyperblocks);
  baby_spncci_observable_hyperblocks.resize(num_hyperblocks);
  // std::cout<<" Read in each hyperblock "<< num_hyperblocks<<std::endl;
  for(int hypersector_index=0; hypersector_index<num_hyperblocks; ++hypersector_index)
    {
      baby_spncci_observable_hyperblocks[hypersector_index].resize(1);
      spncci::ReadBlock(in_stream, baby_spncci_observable_hyperblocks[hypersector_index][0]);
    }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void WriteHyperSectors(
    const spncci::ObservableBabySpNCCIHypersectors& baby_spncci_observable_hypersectors,
    const std::string& filename,
    int irrep_family_index_bra,int irrep_family_index_ket
    )
  {
    std::ios_base::openmode mode_argument = std::ios_base::out | std::ios::app | std::ios_base::binary;
    std::ofstream out_file;
    out_file.open(filename,mode_argument);

    if (!out_file)
      {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
      }

    int num_hypersectors=baby_spncci_observable_hypersectors.size(); 
    
    // Write irrep family indices and number of hyperblocks 
    mcutils::WriteBinary<int>(out_file,irrep_family_index_bra);
    mcutils::WriteBinary<int>(out_file,irrep_family_index_ket);
    mcutils::WriteBinary<int>(out_file,num_hypersectors);

    // for each block, 
    for(int hypersector_index=0; hypersector_index<num_hypersectors; ++hypersector_index)
      {
        const auto& baby_spncci_hypersector
            =baby_spncci_observable_hypersectors.GetHypersector(hypersector_index);

        int baby_spncci_index_bra, baby_spncci_index_ket, operator_subspace_index, rho0;
        std::tie(baby_spncci_index_bra,baby_spncci_index_ket,operator_subspace_index,rho0)
            =baby_spncci_hypersector.Key();

        // write number of rows and columns
        mcutils::WriteBinary<int>(out_file,baby_spncci_index_bra);
        mcutils::WriteBinary<int>(out_file,baby_spncci_index_ket);
        mcutils::WriteBinary<int>(out_file,operator_subspace_index);
        mcutils::WriteBinary<int>(out_file,rho0);
      }

    out_file.close();
  }



void ReadObservableHypersectors(
  std::istream& in_stream,
  spncci::LGIPair& lgi_pair,
  std::vector<spncci::ObservableHypersectorLabels>& list_baby_spncci_hypersectors,
  int& num_hypersectors
)
{
  
  // Read lgi family 
  int irrep_family_index_bra, irrep_family_index_ket;
  mcutils::ReadBinary<int>(in_stream,irrep_family_index_bra);
  mcutils::ReadBinary<int>(in_stream,irrep_family_index_ket);
  lgi_pair=spncci::LGIPair(irrep_family_index_bra,irrep_family_index_ket);

  // std::cout<<irrep_family_index_bra<<"  "<<irrep_family_index_ket<<"  "<<omp_get_thread_num()<<std::endl;
  // Read number of hyperblocks 
  mcutils::ReadBinary<int>(in_stream,num_hypersectors);
  // std::cout<<"number of hypersectors "<<num_hypersectors<<std::endl;
  list_baby_spncci_hypersectors.resize(num_hypersectors);
  // std::cout<<" Read in each hyperblock "<< num_hyperblocks<<std::endl;
  for(int hypersector_index=0; hypersector_index<num_hypersectors; ++hypersector_index)
    {
      int bra_index,ket_index,operator_index,rho0;
      mcutils::ReadBinary<int>(in_stream,bra_index);
      mcutils::ReadBinary<int>(in_stream,ket_index);
      mcutils::ReadBinary<int>(in_stream,operator_index);
      mcutils::ReadBinary<int>(in_stream,rho0);
      list_baby_spncci_hypersectors[hypersector_index]
        =spncci::ObservableHypersectorLabels(bra_index,ket_index,operator_index,rho0);
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void WriteBabySpncciObservableRMEs(
    const spncci::LGIPair& lgi_pair,
    const spncci::ObservableHypersectorsTable& baby_spncci_observable_hypersectors_table,
    spncci::ObservableHyperblocksTable& observable_hyperblocks_table
    )
  {
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
            
    
    int thread_num=omp_get_thread_num();
    
    // For each observable 
    for(int observable_index=0; observable_index<observable_hyperblocks_table.size(); ++observable_index)
      {
        const auto& observable_hyperblocks_by_lgi_by_hw=observable_hyperblocks_table[observable_index];
        const spncci::ObservableBabySpNCCIHypersectors& baby_spncci_observable_hypersectors
                =baby_spncci_observable_hypersectors_table[observable_index];
        
        std::string hypersectors_filename
          =fmt::format("hyperblocks/observable_hypersectors_{:02d}_{:02d}.rmes",observable_index,thread_num);
        spncci::WriteHyperSectors(
          baby_spncci_observable_hypersectors,hypersectors_filename,
          irrep_family_index_bra,irrep_family_index_ket
        );

        // for each hw value
        for(int hw_index=0; hw_index<observable_hyperblocks_by_lgi_by_hw.size(); ++hw_index)
          {
            const basis::OperatorHyperblocks<double>& baby_spncci_observable_hyperblocks
              =observable_hyperblocks_by_lgi_by_hw[hw_index];


            
            // One file per observable, per hw value
            std::string filename=fmt::format("hyperblocks/observable_hyperblocks_{:02d}_{:02d}_{:02d}.rmes",observable_index,hw_index,thread_num);
            
            // #pragma omp critical (write_observables)
            // Now each thread writes to a separate file
            spncci::WriteHyperBlock(
              baby_spncci_observable_hyperblocks,filename,
              irrep_family_index_bra,irrep_family_index_ket
            );
          }
      }        
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Old version
    void WriteBabySpncciObservableRMEs(
    const spncci::LGIPair& lgi_pair,
    spncci::ObservableHyperblocksTable& observable_hyperblocks_table
    )
  {
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;

    // For each observable 
    for(int observable_index=0; observable_index<observable_hyperblocks_table.size(); ++observable_index)
      {
        const auto& observable_hyperblocks_by_lgi_by_hw=observable_hyperblocks_table[observable_index];
        
        // for each hw value
        for(int hw_index=0; hw_index<observable_hyperblocks_by_lgi_by_hw.size(); ++hw_index)
          {
            const basis::OperatorHyperblocks<double>& baby_spncci_observable_hyperblocks
              =observable_hyperblocks_by_lgi_by_hw[hw_index];

            // One file per observable, per hw value
            int thread_num=omp_get_thread_num();
            std::string filename=fmt::format("hyperblocks/observable_hyperblocks_{:02d}_{:02d}_{:02d}.rmes",observable_index,hw_index,thread_num);
            
            // #pragma omp critical (write_observables)
            // Now each thread writes to a separate file
            spncci::WriteHyperBlock(
              baby_spncci_observable_hyperblocks,filename,
              irrep_family_index_bra,irrep_family_index_ket
            );
          }
      }        
  }



}  // namespace

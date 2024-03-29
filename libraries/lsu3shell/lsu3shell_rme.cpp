/****************************************************************
  lsu3shell_rme.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT
****************************************************************/
#include "lsu3shell/lsu3shell_rme.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "fmt/format.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"

namespace lsu3shell
{

  // global mode setting for rme I/O
  bool g_rme_binary_format = true;
  typedef unsigned int RMEIndexType;

  // debugging flags
  bool g_verbose_rme_listing = false;
  bool g_verbose_rme_block_statistics = false;

  void 
  ReadLSU3ShellRMEsText(
      const std::string& filename,
      const LSU3ShellBasisTable& lsu3_basis_table_bra,
      const u3shell::SpaceU3SPN& space_bra, 
      const LSU3ShellBasisTable& lsu3_basis_table_ket,
      const u3shell::SpaceU3SPN& space_ket, 
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      lsu3shell::OperatorBlocks& blocks,
      double scale_factor
    )
  {
    // open file
    std::ifstream in_stream(filename);
    mcutils::StreamCheck(bool(in_stream),filename,"Failure opening lsu3shell rme file");

    // process stream
    basis::SetOperatorToZero(sectors,blocks);
    std::string line;
    while(std::getline(in_stream,line))
      {
        // skip initial header line
        if(not std::isdigit(line[0]))
          continue;

        // read bra/ket lsu3shell basis multiplicity group indices
        int i,j;
        std::istringstream line_stream(line);
        line_stream >> i >> j;
        // std::cout<<i<<" "<<j<<std::endl;

        // retrieve lsu3shell basis multiplicity group information
        u3shell::U3SPN omegaSPNi, omegaSPNj;
        const LSU3ShellBasisGroupData& group_i = lsu3_basis_table_bra[i];
        const LSU3ShellBasisGroupData& group_j = lsu3_basis_table_ket[j];

        u3::SU3 xi(group_i.omegaSPN.SU3());
        u3::SU3 xj(group_j.omegaSPN.SU3());
        // std::cout<<fmt::format("{}  {}  {}", group_i.omegaSPN.Str(), operator_labels.Str(),group_j.omegaSPN.Str())<<std::endl;
        int rho0_max=u3::OuterMultiplicity(xj,operator_labels.x0(),xi);
        // std::cout<<group_i.dim<<"  "<<group_j.dim<<"  "<<rho0_max<<std::endl;
        // extract and store matrix elements
        int i_space=space_bra.LookUpSubspaceIndex(group_i.omegaSPN);
        int j_space=space_ket.LookUpSubspaceIndex(group_j.omegaSPN);

        // warn if file provides an rme of multiplicity zero
        //
        // This still falls short of verifying the the multiplicity
        // given in the file, but at least it's something...
        if (rho0_max==0)
          {
            std::cout
              << "WARN: input rme entry present for rme with expected zero rho0_max"
              << std::endl
              << fmt::format("i {} j {}:  {}  {}  {}",i,j,group_i.omegaSPN.Str(),operator_labels.Str(),group_j.omegaSPN.Str())
              << std::endl
              << fmt::format("{}: {}",filename,line)
              <<std::endl;
          }

        // verbose rme diagnostics
        if (g_verbose_rme_listing)
          {
            std::cout
              << fmt::format("{}: i {:3d} j {:3d} rho0_max {:3d}:  {}  {}  {}",filename,i,j,rho0_max,group_i.omegaSPN.Str(),operator_labels.Str(),group_j.omegaSPN.Str())
              << std::endl;
          }

        // extract and store matrix elements
        for(int gi=0; gi<group_i.dim; ++gi)
          for(int gj=0; gj<group_j.dim; ++gj)
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                double rme;
                // std::cout<<"getting rme"<<std::endl;
                line_stream >> rme;

                // verbose rme diagnostics
                if (g_verbose_rme_listing)
                  {
                    std::cout
                      << fmt::format("  {:e}",rme)
                      << std::endl;
                  }

                // std::cout<<fmt::format("{} {}  {} {} {}  {}",i,j,i_space,j_space,rho0,rme)<<std::endl;
                // Note: Since rho0 is most rapidly varying index in sector enumeration, we could just 
                // calculate the sector_index by offsetting from the sector with rho0=1.
                int sector_index=sectors.LookUpSectorIndex(i_space,j_space,rho0);
                assert(sector_index!=basis::kNone);
                int row_index=group_i.start_index+gi;
                int column_index=group_j.start_index+gj;
                blocks[sector_index](row_index,column_index)=scale_factor*rme;
              }        
      }

    // close file
    in_stream.close();
  };

  void 
  ReadLSU3ShellRMEsBinary(
      bool sp3r_generators,
      const std::string& filename,
      const LSU3ShellBasisTable& lsu3_basis_table_bra,
      const u3shell::SpaceU3SPN& space_bra, 
      const LSU3ShellBasisTable& lsu3_basis_table_ket,
      const u3shell::SpaceU3SPN& space_ket, 
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      lsu3shell::OperatorBlocks& blocks,
      double scale_factor
    )
  {
    // open file
    std::ifstream in_stream(filename,std::ios_base::in|std::ios_base::binary);
    mcutils::StreamCheck(bool(in_stream),filename,"Failure opening lsu3shell rme file");

    // read file header
    int format_code;
    mcutils::ReadBinary<int>(in_stream,format_code);
    assert(format_code==1);
    int float_precision;
    mcutils::ReadBinary<int>(in_stream,float_precision);
    assert((float_precision==4)||(float_precision==8));
    // std::cout
    //   << fmt::format("RME input: filename {}, format_code {}, float_precision {}",filename,format_code,float_precision)
    //   << std::endl;

    // allocate matrices for operator
    basis::SetOperatorToZero(sectors,blocks);

    // read rmes
    while(in_stream)
      {
        // read bra/ket lsu3shell basis multiplicity group indices
        RMEIndexType i, j;
        mcutils::ReadBinary<RMEIndexType>(in_stream,i);
        mcutils::ReadBinary<RMEIndexType>(in_stream,j);

        // quit if this has brought us past end of file
        if (!in_stream)
          break;

        // retrieve lsu3shell basis multiplicity group information
        u3shell::U3SPN omegaSPNi, omegaSPNj;
        assert((i<lsu3_basis_table_bra.size())&&(j<lsu3_basis_table_ket.size()));
        const LSU3ShellBasisGroupData& group_i = lsu3_basis_table_bra[i];
        const LSU3ShellBasisGroupData& group_j = lsu3_basis_table_ket[j];
        u3::SU3 xi(group_i.omegaSPN.SU3());
        u3::SU3 xj(group_j.omegaSPN.SU3());
        // std::cout<<fmt::format("{}  {}  {}", group_i.omegaSPN.Str(), operator_labels.Str(),group_j.omegaSPN.Str())<<std::endl;
        int rho0_max=u3::OuterMultiplicity(xj,operator_labels.x0(),xi);
        // std::cout<<group_i.dim<<"  "<<group_j.dim<<"  "<<rho0_max<<std::endl;
        int i_subspace_index=space_bra.LookUpSubspaceIndex(group_i.omegaSPN);
        int j_subspace_index=space_ket.LookUpSubspaceIndex(group_j.omegaSPN);
        assert((i_subspace_index!=basis::kNone)&&(j_subspace_index!=basis::kNone));

        // verify multiplicity given in file
        RMEIndexType num_rmes = group_i.dim*group_j.dim*rho0_max;
        if(rho0_max==0)
          std::cout<<fmt::format("{} {}  {}  {}  {}", i,j,group_i.omegaSPN.Str(), operator_labels.Str(),group_j.omegaSPN.Str())<<std::endl;
        mcutils::VerifyBinary<RMEIndexType>(
            in_stream,num_rmes,
            fmt::format("Unexpected value encountered reading binary rme file {}",filename),"rho0_max"
          );
        
        // verbose rme diagnostics
        if (g_verbose_rme_listing)
          {
            std::cout
              << fmt::format("{}: i {:3d} j {:3d} rho0_max {:3d}:  {}  {}  {}",filename,i,j,rho0_max,group_i.omegaSPN.Str(),operator_labels.Str(),group_j.omegaSPN.Str())
              << std::endl;
          }

        // extract and store matrix elements
        for(int gi=0; gi<group_i.dim; ++gi)
          for(int gj=0; gj<group_j.dim; ++gj)
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                // read rme
                double rme;
                if (float_precision==4)
                  {
                    float rme_float;
                    mcutils::ReadBinary<float>(in_stream,rme_float);
                    rme = rme_float;
                  }
                else if (float_precision==8)
                  {
                    mcutils::ReadBinary<double>(in_stream,rme);
                  }

                // verbose rme diagnostics
                if (g_verbose_rme_listing)
                  {
                    std::cout
                      << fmt::format("  {:e}",rme)
                      << std::endl;
                  }

                if(sp3r_generators)
                  {
                    std::stringstream ss;
                    ss << rme;
                    ss >> rme;
                  }

                // Note: Since rho0 is most rapidly varying index in sector enumeration, we could just 
                // calculate the sector_index by offsetting from the sector with rho0=1.
                int sector_index=sectors.LookUpSectorIndex(i_subspace_index,j_subspace_index,rho0);
                assert(sector_index!=basis::kNone);
                int row_index=group_i.start_index+gi;
                int column_index=group_j.start_index+gj;
                blocks[sector_index](row_index,column_index)=scale_factor*rme;
              }        
      }
    // close file
    in_stream.close();
  };

  void 
  ReadLSU3ShellRMEs(
      bool sp3r_generators,
      const std::string& filename,
      const LSU3ShellBasisTable& lsu3_basis_table_bra,
      const u3shell::SpaceU3SPN& space_bra, 
      const LSU3ShellBasisTable& lsu3_basis_table_ket,
      const u3shell::SpaceU3SPN& space_ket, 
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      lsu3shell::OperatorBlocks& blocks,
      double scale_factor
    )
  {

      // LSU3ShellBasisTable lsu3_basis_table_bra=lsu3_basis_table;
      // LSU3ShellBasisTable lsu3_basis_table_ket=lsu3_basis_table;
      // u3shell::SpaceU3SPN space_bra=space;
      // u3shell::SpaceU3SPN space_ket=space;

    // std::ios_base::openmode mode_argument = std::ios_base::in;
    // if (g_rme_binary_format)
    //   mode_argument |= std::ios_base::binary;

    if (g_rme_binary_format)
      // binary format
      ReadLSU3ShellRMEsBinary(sp3r_generators,filename,lsu3_basis_table_bra,space_bra,lsu3_basis_table_ket,space_ket,operator_labels,sectors,blocks,scale_factor);
    else
      // text format
      ReadLSU3ShellRMEsText(filename,lsu3_basis_table_bra,space_bra,lsu3_basis_table_ket,space_ket,operator_labels,sectors,blocks,scale_factor);

    if (g_verbose_rme_block_statistics)
      {
        std::cout << fmt::format("OPERATOR {}: sectors {} ",filename,blocks.size()) << std::endl;
        for(int sector_index=0; sector_index<blocks.size(); ++sector_index)
          {
            auto& block = blocks[sector_index];
            std::cout
              << fmt::format(
                  "  sector_index {:d} {:d}x{:d}: norm {:e}",
                  sector_index,block.rows(),block.cols(),block.norm()
                )
              << std::endl;
          }
      }                         
                                 
  }


  void 
  ReadLSU3ShellRMEs(
      bool sp3r_generators,
      const std::string& filename,
      const LSU3ShellBasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      lsu3shell::OperatorBlocks& blocks,
      double scale_factor
    )
  //Wrapper function for when bra space and ket space are the same.
  {
      const LSU3ShellBasisTable& lsu3_basis_table_bra=lsu3_basis_table;
      const LSU3ShellBasisTable& lsu3_basis_table_ket=lsu3_basis_table;
      const u3shell::SpaceU3SPN& space_bra=space;
      const u3shell::SpaceU3SPN& space_ket=space;
      ReadLSU3ShellRMEs(
        sp3r_generators,filename,
        lsu3_basis_table_bra,space_bra, 
        lsu3_basis_table_ket,space_ket,
        operator_labels,sectors,
        blocks,scale_factor
      );
  }

  void 
  ReadLSU3ShellRMEs(
      const std::string& filename,
      const LSU3ShellBasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      lsu3shell::OperatorBlocks& blocks,
      double scale_factor
    )
  {
    LSU3ShellBasisTable lsu3_basis_table_bra=lsu3_basis_table;
    LSU3ShellBasisTable lsu3_basis_table_ket=lsu3_basis_table;
    u3shell::SpaceU3SPN space_bra=space;
    u3shell::SpaceU3SPN space_ket=space;

    bool sp3r_generators=false;
    lsu3shell::ReadLSU3ShellRMEs(
      sp3r_generators,filename,lsu3_basis_table,space, 
      operator_labels, sectors,blocks,scale_factor
    );

  }

  void
    ReadLSU3ShellSymplecticOperatorRMEs(
        const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
        const u3shell::SpaceU3SPN& lsu3shell_space, 
        const std::string& Brel_filename, u3shell::SectorsU3SPN& Bintr_sectors, lsu3shell::OperatorBlocks& Bintr_matrices,
        const std::string& Nrel_filename, u3shell::SectorsU3SPN& Nintr_sectors, lsu3shell::OperatorBlocks& Nintr_matrices,
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
        const std::string& Arel_filename, u3shell::SectorsU3SPN& Aintr_sectors, lsu3shell::OperatorBlocks& Aintr_matrices,
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


  bool 
  CompareLSU3ShellRMEs(
      std::ostream& log_stream,
      const U3SPNBasisLSU3Labels& basis_provenance,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::SectorsU3SPN& sectors,
      const lsu3shell::OperatorBlocks& matrices1,
      const lsu3shell::OperatorBlocks& matrices2,
      double tolerance,
      bool verbose
    )
  {  

    // initialize statistics
    int entries_compared = 0;
    double max_residual = 0.;
    double total_sqr_residual = 0.;
    bool success = true;

    // iterate over sectors
    for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
      {

        // retrieve sector
        const u3shell::SectorsU3SPN::SectorType& sector = sectors.GetSector(sector_index);
        if (verbose)
          log_stream
            << fmt::format(
                "sector {} bra {} ket {} dim {}x{}",
                sector_index,
                sector.bra_subspace().LabelStr(),
                sector.ket_subspace().LabelStr(),
                sector.bra_subspace().size(),
                sector.ket_subspace().size()
              )
            << std::endl;

        // iterate over matrix elements
        for (int bra_index=0; bra_index<sector.bra_subspace().size(); ++bra_index)
          for (int ket_index=0; ket_index<sector.ket_subspace().size(); ++ket_index)
            {
              // retrieve matrix elements
              double rme1 = matrices1[sector_index](bra_index,ket_index);
              double rme2 = matrices2[sector_index](bra_index,ket_index);

              // compare matrix elements
              double residual = std::fabs(rme2-rme1);
              ++entries_compared;
              max_residual = std::max(max_residual,residual);
              total_sqr_residual += mcutils::sqr(residual);
              bool entries_agree = (residual <= tolerance);
              // std::cout<<std::endl<<fmt::format("residual {}  tolerance {}  bool {}", residual,tolerance,entries_agree)<<std::endl<<std::endl;
              success &= entries_agree;

              // write entry diagnostics
              if (!entries_agree)
                {
                  log_stream
                    << fmt::format(
                        "  FAIL: {}:({},{}) rme1 {:e} rme2 {:e} residual {:e}",
                        sector_index,bra_index,ket_index,
                        rme1,rme2,residual
                      )
                    << std::endl;
                  const lsu3shell::LSU3ShellBasisGroupLabels& bra_labels = basis_provenance[sector.bra_subspace_index()][bra_index];
                  const lsu3shell::LSU3ShellBasisGroupLabels& ket_labels = basis_provenance[sector.ket_subspace_index()][ket_index];
                  log_stream
                    << fmt::format(
                        "  bra ip {} in {} Np {} Nn {} ; ket ip {} in {} Np {} Nn {}",
                        bra_labels.ip,bra_labels.in,bra_labels.Np,bra_labels.Nn,
                        ket_labels.ip,ket_labels.in,ket_labels.Np,ket_labels.Nn
                      )
                    << std::endl;
                }
            }
      }

    // generate global diagnostics
    double rms_residual = std::sqrt(total_sqr_residual/entries_compared);
    log_stream
      << fmt::format(
          "entries {} rms_residual {:e} max_residual {:e} total_sqr_residual {:e}",
          entries_compared,rms_residual,max_residual,total_sqr_residual
        )
      << std::endl<<std::endl;
        

    return success;
  }

  void GenerateLSU3ShellNcmRMEs(
      const u3shell::SpaceU3SPN& space,
      const u3shell::SectorsU3SPN& Nrel_sectors,
      const lsu3shell::OperatorBlocks& Nrel_matrices,
      int A,
      lsu3shell::OperatorBlocks& Ncm_matrices
    )
  {

    // populate matrices for Ncm
    Ncm_matrices.resize(Nrel_sectors.size());
    for(int sector_index=0; sector_index<Nrel_sectors.size(); ++sector_index)
      {
        // define alias for subspace
        //
        // Note: Sectors here are diagonal and in 1-1 correspondence
        // with subspaces.  We do not actually need the space as an
        // argument, but keep it for uniformity.
        assert(Nrel_sectors.GetSector(sector_index).IsDiagonal());
        const auto& subspace = Nrel_sectors.GetSector(sector_index).ket_subspace();

        // calculate Ncm
        //
        // Obtain as Ncm=N-Nrel, where eigenvalue of N is given by
        // total number of oscilator quanta minus the zero point
        // energy 3A/2.
        
        HalfInt N = subspace.N()-HalfInt(3*A,2);
        int dim = subspace.size();
        Ncm_matrices[sector_index]
          = Eigen::MatrixXd::Identity(dim,dim)*double(N)-Nrel_matrices[sector_index];
      }

  }

  ////////////////////////////////////////////////////////////////
}// end namespace

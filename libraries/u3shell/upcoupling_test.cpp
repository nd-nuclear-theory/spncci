/****************************************************************
  upcoupling_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
****************************************************************/
#include <fstream>

#include "cppformat/format.h"
#include "basis/lsjt_operator.h"

#include "am/am.h"
// #include "am/wigner_gsl.h"
// #include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/import_interaction.h"
#include "u3shell/relative_operator.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/relative_cm_xform.h" 
// #include "moshinsky/moshinsky_xform.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/upcoupling.h"
#include "u3shell/u3st_scheme.h"

// typedefs for going from Relative to Relative-Center-of-Mass
typedef std::tuple<int,int,u3::SU3,HalfInt, HalfInt> RelativeCMU3STLabels;
typedef std::tuple<u3::SU3,HalfInt, HalfInt,int, int,RelativeCMU3STLabels, RelativeCMU3STLabels,int> RelativeCMU3STBraket;
typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> RelativeCMLSTLabels;
typedef std::tuple<int,HalfInt,HalfInt,RelativeCMLSTLabels,RelativeCMLSTLabels> RelativeCMLSTBraket;

//TODO remove import interaction and replace with output file from shell 
void IdentityTest(
  int Nmax, int Jmax, int J0, int T0, int g0,
  u3shell::RelativeRMEsU3ST& relative_rme_map
 )
// Checking up coupling using identity function found in ImportInteraction
// Needs to be updated to use shell identity file
{
  std::cout<<"Identity test"<<std::endl;

  basis::RelativeSpaceLSJT relative_space(Nmax, Jmax);
  basis::RelativeSectorsLSJT relative_sectors(relative_space, J0,T0, g0);
  std::vector<Eigen::MatrixXd> sector_vector;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;


  std::string interaction_file="NONE";
  sector_vector
    =u3shell::ImportInteraction(interaction_file, relative_space, relative_sectors, "Identity");
  //////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"UpcouplingNLST"<<std::endl;
  u3shell::UpcouplingNLST(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);
  for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
    {
      int L0, S0, L,S,T, Lp, Sp, Tp;
      u3shell::RelativeSubspaceLabelsNLST bra, ket;
      std::tie(L0,S0,bra,ket)=it->first;
      std::tie(L,S,T)=ket;
      std::tie(Lp,Sp,Tp)=bra;
      Eigen::MatrixXd sectorNLST=it->second;
      if(fabs(sectorNLST.sum())>10e-10)
        std::cout<<fmt::format("{} {} ({},{},{}) ({},{},{})", L0,S0,Lp,Sp,Tp,L,S,T)
        <<std::endl<<sectorNLST<<std::endl;
    }

  u3shell::Upcoupling(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);
  std::cout<<"UpcouplingU3ST"<<std::endl;
  for (auto it=relative_rme_map.begin(); it!= relative_rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST labels;
      int kappa0, L0;
      std::tie(labels,kappa0,L0)=it->first;
      double coef=it->second;
      if (fabs(coef)>10e-13)
        std::cout<<labels.Str()<<"  "<<kappa0<<"  "<<L0<<std::endl<<it->second<<std::endl<<std::endl;
    }
}

void
KineticCheck()
// Checking upcoupling using kinetic energy (k^2) using function in import_interaction
// Function given analytically 
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int Nmax=6;
  int Jmax=4;
  int J0=0;
  int T0=0;
  int g0=0;
  std::cout<<"Ksqr check"<<std::endl;
  std::string interaction_file="/Users/annamccoy/projects/spncci/libraries/moshinsky/test/ksqr_Nmax06_rel.dat";
  basis::RelativeSectorsLSJT relative_lsjt_sectors;
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  basis::MatrixVector sector_vector;
  u3shell::GetInteractionMatrix(interaction_file, relative_lsjt_space,relative_lsjt_sectors,sector_vector);

  //upcouple to LST
  std::cout<<"Upcoupling to NLST"<<std::endl;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  u3shell::UpcouplingNLST(relative_lsjt_space,relative_lsjt_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);

  // Upcouple to U(3) level
  u3shell::RelativeRMEsU3ST rme_map;
  std::cout<<"Upcoupling to U3ST"<<std::endl;
  u3shell::UpcouplingU3ST(rme_nlst_map, T0, Nmax, rme_map);
  for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST op_labels;
      int kappa0,L0;
      std::tie(op_labels, kappa0,L0)=it->first;
      double rme=it->second;
      double check=u3shell::RelativeKineticEnergyOperator(op_labels.bra(), op_labels.ket());
      if(fabs(rme)>10e-13)
        std::cout<<fmt::format("{} {} {}   {}   {}",op_labels.Str(), kappa0,L0,rme,check)<<std::endl;
    }
}

void
ReadWriteCheck(
    u3shell::RelativeRMEsU3ST& relative_rme_map,
    std::string filename
  )
  {
     std::ostringstream os;
    WriteRelativeOperatorU3ST(os,relative_rme_map);
    std::ofstream stream;
    stream.open(filename.c_str());
    stream<<os.str();
    // std::cout<<os.str();
    stream.close();

    std::ifstream is(filename.c_str());
    if(!is)
      std::cout<<"Didn't open"<<std::endl;
    u3shell::RelativeRMEsU3ST k2_relative_rmes;
    u3shell::ReadRelativeOperatorU3ST(is, k2_relative_rmes);
  }

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  int Nmax=10;
  int Jmax=Nmax+2;
  int J0=0;
  int g0=0;
	int T0=0;

  u3shell::RelativeRMEsU3ST id_relative_rme_map;
  IdentityTest(Nmax,Jmax,J0,T0, g0, id_relative_rme_map);

  u3shell::RelativeRMEsU3ST ke_relative_rme_map;
  KineticCheck();

  std::string filename="Trel_upcouled";
  ReadWriteCheck(ke_relative_rme_map,filename);
}
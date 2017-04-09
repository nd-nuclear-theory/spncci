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

double zero_threshold=1e-6;

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
  // basis::RelativeSectorsLSJT relative_sectors(relative_space, J0,T0, g0);
  std::vector<Eigen::MatrixXd> sector_vector; 

  std::string interaction_file="../moshinsky/test/identity_Nmax04_rel.dat";
  // basis::RelativeSectorsLSJT relative_sectors_lsjt;
  basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Jmax);
  
  std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
  std::array<basis::MatrixVector,3> isospin_component_matrices_lsjt;
  basis::OperatorLabelsJT operator_labels;

  basis::ReadRelativeOperatorLSJT(
    interaction_file,relative_space_lsjt,operator_labels,
    isospin_component_sectors_lsjt, isospin_component_matrices_lsjt, true
    );

  //////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"UpcouplingNLST"<<std::endl;

  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  const basis::RelativeSectorsLSJT& sectors_lsjt=isospin_component_sectors_lsjt[T0];
  const basis::MatrixVector& matrices_lsjt=isospin_component_matrices_lsjt[T0];
  u3shell::UpcouplingNLST(relative_space_lsjt,sectors_lsjt,matrices_lsjt,J0,g0,T0,Nmax,rme_nlst_map);

  for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
    {
      int L0, S0, L,S,T, Lp, Sp, Tp;
      u3shell::RelativeSubspaceLabelsNLST bra, ket;
      std::tie(L0,S0,T0,bra,ket)=it->first;
      std::tie(L,S,T)=ket;
      std::tie(Lp,Sp,Tp)=bra;
      Eigen::MatrixXd sectorNLST=it->second;
      if(fabs(sectorNLST.sum())>zero_threshold)
        std::cout<<fmt::format("{} {} ({},{},{}) ({},{},{})", L0,S0,Lp,Sp,Tp,L,S,T)
        <<std::endl<<sectorNLST<<std::endl;
    }

  u3shell::Upcoupling(
    relative_space,isospin_component_sectors_lsjt,isospin_component_matrices_lsjt,
    J0,g0,T0,Nmax,relative_rme_map
    );

  std::cout<<"UpcouplingU3ST"<<std::endl;
  for (auto it=relative_rme_map.begin(); it!= relative_rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST labels;
      int kappa0, L0;
      std::tie(labels,kappa0,L0)=it->first;
      double coef=it->second;
      if (fabs(coef)>zero_threshold)
        std::cout<<labels.Str()<<"  "<<kappa0<<"  "<<L0<<std::endl<<it->second<<std::endl<<std::endl;
    }
}

void
KineticCheck(u3shell::RelativeRMEsU3ST& rme_map)
// Checking upcoupling using kinetic energy (k^2) using function in import_interaction
// Function given analytically 
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int Nmax=20;
  int Jmax=4;
  int J0=0;
  int T0=0;
  int g0=0;
  std::cout<<"Ksqr check"<<std::endl;
  std::string interaction_file="../moshinsky/test/ksqr_Nmax06_rel.dat";
  basis::RelativeSectorsLSJT relative_lsjt_sectors;
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  
  std::array<basis::RelativeSectorsLSJT,3> T0_sector_labels_lsjt;
  std::array<basis::MatrixVector,3> T0_sectors_lsjt;
  basis::OperatorLabelsJT operator_labels;


  basis::ReadRelativeOperatorLSJT(
    interaction_file,relative_lsjt_space,operator_labels,
    T0_sector_labels_lsjt, T0_sectors_lsjt, true
    // relative_component_sectors,relative_component_matrices, true
    );


  
  // u3shell::GetInteractionMatrix(interaction_file, relative_lsjt_space,relative_lsjt_sectors,sector_vector);

  //upcouple to LST
  std::cout<<"Upcoupling to NLST"<<std::endl;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  const basis::RelativeSectorsLSJT& sector_labels_lsjt=T0_sector_labels_lsjt[T0];
  const basis::MatrixVector& sectors_lsjt=T0_sectors_lsjt[T0];
  u3shell::UpcouplingNLST(relative_lsjt_space,sector_labels_lsjt,sectors_lsjt,J0,g0,T0,Nmax,rme_nlst_map);

  // Upcouple to U(3) level
  // u3shell::RelativeRMEsU3ST rme_map;
  std::cout<<"Upcoupling to U3ST"<<std::endl;
  u3shell::UpcouplingU3ST(rme_nlst_map, Nmax, rme_map);
  for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST op_labels;
      int kappa0,L0;
      std::tie(op_labels, kappa0,L0)=it->first;
      double rme=it->second;
      double check=u3shell::RelativeKineticEnergyOperator(op_labels.bra(), op_labels.ket());
      if(fabs(rme)>zero_threshold)
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
    std::cout<<os.str();
    stream.close();

    std::ifstream is(filename.c_str());
    if(!is)
      std::cout<<"Didn't open"<<std::endl;
    u3shell::RelativeRMEsU3ST relative_rmes2;
    u3shell::ReadRelativeOperatorU3ST(is, relative_rmes2);
    for(auto it=relative_rmes2.begin(); it!=relative_rmes2.end(); ++it)
      {
        int kappa0,L0;
        u3shell::RelativeUnitTensorLabelsU3ST tensor;
        std::tie(tensor,kappa0,L0)=it->first;
        if(not relative_rme_map.count(it->first))
        {
          std::cout<<fmt::format("[{} {} {}] {} not in write map", tensor.Str(),kappa0,L0,it->second)<<std::endl;
          continue;
        }
        if(fabs(relative_rme_map[it->first]-relative_rmes2[it->first])<zero_threshold)
          std::cout<<fmt::format("[{} {} {}]", tensor.Str(),kappa0,L0)
          <<"  "<<relative_rme_map[it->first]
          <<"  "<<relative_rmes2[it->first]
          <<"  "<<fabs(relative_rme_map[it->first]-relative_rmes2[it->first])
          <<std::endl;
      }
    for(auto it=relative_rme_map.begin(); it!=relative_rme_map.end(); ++it)
      {
        if(fabs(it->second)>zero_threshold)
          if(not relative_rmes2.count(it->first))
          {
            int kappa0,L0;
            u3shell::RelativeUnitTensorLabelsU3ST tensor;
            std::tie(tensor,kappa0,L0)=it->first;
            std::cout<<fmt::format("[{} {} {}] {} not in read map", tensor.Str(),kappa0,L0,it->second)<<std::endl;
          }
        
      }
  }

int main(int argc, char **argv)
{
  // double zero_threshold=10e-6;
  u3::U3CoefInit();
  int Nmax=10;
  int Jmax=Nmax+2;
  int J0=0;
  int g0=0;
	int T0=0;
  u3shell::RelativeRMEsU3ST id_relative_rme_map;
  IdentityTest(Nmax,Jmax,J0,T0, g0, id_relative_rme_map);

  u3shell::RelativeRMEsU3ST ke_relative_rme_map;
  KineticCheck(ke_relative_rme_map);

  std::string filename="Trel_upcouled";
  ReadWriteCheck(ke_relative_rme_map,filename);

  std::string id_filename="Id_upcouled";
  std::ofstream os(id_filename.c_str());
  u3shell::WriteRelativeOperatorU3ST(os,id_relative_rme_map);
  os.close();
}
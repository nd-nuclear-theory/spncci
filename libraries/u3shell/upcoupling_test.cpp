/****************************************************************
  upcoupling_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
  5/14/19 (aem): Updated basis::OperatorLabelsJT to 
    basis::RelativeOperatorParametersLSJT for reading in operators
****************************************************************/
#include <fstream>

#include "fmt/format.h"
#include "am/am.h"
#include "mcutils/eigen.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/relative_cm_xform.h" 
#include "u3shell/tensor_labels.h"
#include "u3shell/upcoupling.h"
#include "u3shell/u3st_scheme.h"

double zero_threshold=1e-6;

// typedefs for going from Relative to Relative-Center-of-Mass
typedef std::tuple<int,int,u3::SU3,HalfInt, HalfInt> RelativeCMU3STLabels;
typedef std::tuple<u3::SU3,HalfInt, HalfInt,int, int,RelativeCMU3STLabels, RelativeCMU3STLabels,int> RelativeCMU3STBraket;
typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> RelativeCMLSTLabels;
typedef std::tuple<int,HalfInt,HalfInt,RelativeCMLSTLabels,RelativeCMLSTLabels> RelativeCMLSTBraket;


void IdentityTest(
  int Nmax, int Jmax, int J0, int T0, int g0,
  u3shell::RelativeRMEsU3ST& relative_rme_map
 )
// Checks that the identity operator upcouples correctly
{
  std::cout<<"Identity test"<<std::endl;

  basis::RelativeSpaceLSJT relative_space(Nmax, Jmax);
  // basis::RelativeSectorsLSJT relative_sectors(relative_space, J0,T0, g0);
  std::vector<Eigen::MatrixXd> sector_vector; 

  std::string interaction_file="../moshinsky/test/identity_Nmax04_rel.dat";
  // basis::RelativeSectorsLSJT relative_sectors_lsjt;
  basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Jmax);
  
  std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
  std::array<basis::OperatorBlocks<double>,3> isospin_component_matrices_lsjt;
  basis::RelativeOperatorParametersLSJT operator_labels;

  basis::ReadRelativeOperatorLSJT(
    interaction_file,relative_space_lsjt,operator_labels,
    isospin_component_sectors_lsjt, isospin_component_matrices_lsjt, true
    );

  //////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"UpcouplingNLST"<<std::endl;

  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  const basis::RelativeSectorsLSJT& sectors_lsjt=isospin_component_sectors_lsjt[T0];
  const basis::OperatorBlocks<double>& matrices_lsjt=isospin_component_matrices_lsjt[T0];
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
  // int Nmax=20;
  int Nmax=20;
  int Jmax=Nmax+1;
  int J0=0;
  int T0=0;
  int g0=0;
  std::cout<<"Ksqr check"<<std::endl;
  std::string interaction_file="../moshinsky/test/ksqr_Nmax16_rel.dat";
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  
  std::array<basis::RelativeSectorsLSJT,3> T0_sector_labels_lsjt;
  std::array<basis::OperatorBlocks<double>,3> T0_sectors_lsjt;
  basis::RelativeOperatorParametersLSJT operator_labels;
  basis::ReadRelativeOperatorLSJT(
    interaction_file,relative_lsjt_space,operator_labels,
    T0_sector_labels_lsjt, T0_sectors_lsjt, true
    // relative_component_sectors,relative_component_matrices, true
    );

  //upcouple to LST
  std::cout<<"Upcoupling to NLST"<<std::endl;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  const basis::RelativeSectorsLSJT& sector_labels_lsjt=T0_sector_labels_lsjt[T0];
  const basis::OperatorBlocks<double>& sectors_lsjt=T0_sectors_lsjt[T0];
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
      double check=u3shell::K2rel(op_labels.bra(), op_labels.ket());
      if(fabs(rme)>zero_threshold)
        std::cout<<fmt::format("{:40} {} {}   {:10.4f}   {:10.4f}  {:5.2f}",op_labels.Str(), kappa0,L0,rme,check,rme/check)<<std::endl;
    }
}

void
ReadWriteCheck(
    u3shell::RelativeRMEsU3ST& relative_rme_map,
    std::string filename
  )
  {
    std::string filename2="upcoupling_test.dat";
    bool hermitian=true;
    WriteRelativeOperatorU3ST(filename2,relative_rme_map,hermitian);
    std::ostringstream os(filename2);
    std::ofstream stream;
    stream.open(filename.c_str());
    stream<<os.str();
    std::cout<<os.str();
    stream.close();

    u3shell::RelativeRMEsU3ST relative_rmes2;
    u3shell::ReadRelativeOperatorU3ST(filename, relative_rmes2);
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

void UpcoupleQmass(int Nmax, int Jmax)
  {
    
    u3shell::RelativeRMEsU3ST rme_map;
    basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
    // std::vector<std::string> file_end={"pp","nn"};
    std::vector<std::string> file_end={"total"};
    for(int i=0; i<file_end.size(); ++i)
    {
      std::string filename=fmt::format("../../data/relative_interactions/quadrupole_test_Nmax6_{}_rel.dat",file_end[i]);
      std::cout<<filename<<std::endl;
      std::array<basis::RelativeSectorsLSJT,3> T0_sector_labels_lsjt;
      std::array<basis::OperatorBlocks<double>,3> T0_sectors_lsjt;
      basis::RelativeOperatorParametersLSJT operator_labels;
      basis::ReadRelativeOperatorLSJT(
        filename,relative_lsjt_space,operator_labels,
        T0_sector_labels_lsjt, T0_sectors_lsjt, true
        );
 
      Upcoupling(    
        relative_lsjt_space,
        T0_sector_labels_lsjt,
        T0_sectors_lsjt,
        2, 0, -1, Nmax,
        rme_map
      );
    }
  
    for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
      {
        u3shell::RelativeUnitTensorLabelsU3ST op_labels;
        int kappa0,L0;
        std::tie(op_labels, kappa0,L0)=it->first;
        double rme=it->second;
        double check=u3shell::Qrel(op_labels.bra(), op_labels.ket());        
        if(fabs(rme)>zero_threshold)
          std::cout<<fmt::format("{} {} {}   {}  {}  {}",op_labels.Str(), kappa0,L0,rme,check, check/rme)<<std::endl;
      }

  }


void UpcoupleQisovector(int Nmax, int Jmax)
  {
    
    u3shell::RelativeRMEsU3ST rme_map;
    basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
    // std::vector<std::string> file_end={"pp","nn"};
    std::vector<std::string> file_end={"total"};
    for(int i=0; i<file_end.size(); ++i)
    {
      std::string filename="../../data/relative_interactions/quadrupole-isovector.dat";
      // std::string filename="../../data/relative_interactions/quadrupole-isoscalar.dat";
      std::cout<<filename<<std::endl;
      std::array<basis::RelativeSectorsLSJT,3> T0_sector_labels_lsjt;
      std::array<basis::OperatorBlocks<double>,3> T0_sectors_lsjt;
      basis::RelativeOperatorParametersLSJT operator_labels;
      basis::ReadRelativeOperatorLSJT(
        filename,relative_lsjt_space,operator_labels,
        T0_sector_labels_lsjt, T0_sectors_lsjt, true
        );
 
      Upcoupling(    
        relative_lsjt_space,
        T0_sector_labels_lsjt,
        T0_sectors_lsjt,
        2, 0, -1, Nmax,
        rme_map
      );
    }
  
    for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
      {
        u3shell::RelativeUnitTensorLabelsU3ST op_labels;
        int kappa0,L0;
        std::tie(op_labels, kappa0,L0)=it->first;
        double rme=it->second;
        double check=u3shell::Qrel(op_labels.bra(), op_labels.ket());        
        if(fabs(rme)>zero_threshold)
          std::cout<<fmt::format("{} {} {}   {}  {}  {}",op_labels.Str(), kappa0,L0,rme,check, check/rme)<<std::endl;
      }

  }




void
QCheck()
// Checking upcoupling using kinetic energy (k^2) using function in import_interaction
// Function given analytically 
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // int Nmax=20;
  int Nmax=2;
  int Jmax=2;
  int J0=2;
  int g0=0;
  
  u3shell::RelativeRMEsU3ST rme_map;
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  // std::vector<std::string> file_end={"pp","nn","pn"};
  std::vector<std::string> file_end={"total"};
  std::string filename_base="../../data/relative_interactions/quadrupole_test_Nmax6_{}_rel.dat";

  //upcouple to LST
  std::cout<<"Upcoupling to NLST"<<std::endl;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_ppnn;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_pn;

  for(int i=0; i<1; ++i)
    {
      std::string filename=fmt::format(filename_base,file_end[i]);

      std::array<basis::RelativeSectorsLSJT,3> T0_sector_labels_lsjt;
      std::array<basis::OperatorBlocks<double>,3> T0_sectors_lsjt;
      basis::RelativeOperatorParametersLSJT operator_labels;
      basis::ReadRelativeOperatorLSJT(
        filename,relative_lsjt_space,operator_labels,
        T0_sector_labels_lsjt, T0_sectors_lsjt, true
        );

      for(int T0=0; T0<3; ++T0)
        {
          std::cout<<"T0="<<T0<<std::endl;
          const basis::RelativeSectorsLSJT& sector_labels_lsjt=T0_sector_labels_lsjt[T0];
          const basis::OperatorBlocks<double>& sectors_lsjt=T0_sectors_lsjt[T0];
          u3shell::UpcouplingNLST(relative_lsjt_space,sector_labels_lsjt,sectors_lsjt,J0,g0,T0,Nmax,rme_nlst_ppnn);
          for(auto it=rme_nlst_ppnn.begin(); it!=rme_nlst_ppnn.end(); ++it)
            if(not mcutils::IsZero(it->second,1e-6))
              std::cout<<it->second<<std::endl<<std::endl;
        }

    }

    // std::string filename=fmt::format(filename_base,file_end[2]);

    // std::array<basis::RelativeSectorsLSJT,3> T0_sector_labels_lsjt;
    // std::array<basis::OperatorBlocks<double>,3> T0_sectors_lsjt;
    // basis::OperatorLabelsJT operator_labels;

    // basis::ReadRelativeOperatorLSJT(
    //   filename,relative_lsjt_space,operator_labels,
    //   T0_sector_labels_lsjt, T0_sectors_lsjt, true
    //   );

    // for(int T0=0; T0<3; ++T0)
    //   {
    //     const basis::RelativeSectorsLSJT& sector_labels_lsjt=T0_sector_labels_lsjt[T0];
    //     const basis::OperatorBlocks<double>& sectors_lsjt=T0_sectors_lsjt[T0];
    //     u3shell::UpcouplingNLST(relative_lsjt_space,sector_labels_lsjt,sectors_lsjt,J0,g0,T0,Nmax,rme_nlst_pn);
    //   }
    // for(auto it=rme_nlst_pn.begin(); it!=rme_nlst_pn.end(); ++it)
    //   {
    //     if(not mcutils::IsZero((it->second-rme_nlst_ppnn[it->first]),1e-7))
    //       std::cout<<"sectors"<<std::endl<<it->second<<std::endl<<std::endl<<rme_nlst_ppnn[it->first]<<std::endl<<std::endl;
    //   }



}


int main(int argc, char **argv)
{
  // double zero_threshold=10e-6;
  u3::U3CoefInit();
  // int Nmax=10;
 //  int Nmax=6;
 //  int N1v=1;
 //  int Jmax=Nmax+2;
 //  int J0=0;
 //  int g0=0;
	// int T0=-1;



  // // UpcoupleQmass(Nmax,Jmax);

  // // QCheck();

  // u3shell::RelativeRMEsU3ST id_relative_rme_map;
  // IdentityTest(Nmax,Jmax,J0,T0, g0, id_relative_rme_map);

  // u3shell::RelativeRMEsU3ST ke_relative_rme_map;
  // KineticCheck(ke_relative_rme_map);

  int Nmax=2;
  int Jmax=4;
  UpcoupleQisovector(Nmax, Jmax);

  // std::string filename="Trel_upcouled";
  // ReadWriteCheck(ke_relative_rme_map,filename);

  // std::string id_filename="Id_upcouled";
  // std::ofstream os(id_filename.c_str());
  // u3shell::WriteRelativeOperatorU3ST(os,id_relative_rme_map);
  // os.close();

  // std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_tensor_labels;

  // u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, 1,unit_tensor_labels,J0,T0,false);
  // u3shell::RelativeUnitTensorSpaceU3S unit_tensor_space(Nmax,1,unit_tensor_labels);

  
  // for(auto& tensor : unit_tensor_labels)
  //   std::cout<<tensor.Str()<<std::endl;


  // for(int i=0; i<unit_tensor_space.size(); ++i)
  //   {
  //     std::cout<<unit_tensor_space.GetSubspace(i).Str()<<std::endl;
  //   }  


  // // will throw errors unless you change Nmax in Kinetic check to match Nmax above. 
  // u3shell::RelativeRMEsU3SSubspaces relative_rmes;
  // ReadRelativeOperatorU3ST(Nmax, N1v,filename, unit_tensor_space,relative_rmes);
  // for(auto it=relative_rmes.begin(); it!=relative_rmes.end(); ++it)
  //   {
  //     int index,kappa0,L0;
  //     std::tie(index,kappa0,L0)=it->first;
  //     std::cout<<index<<"  "<<kappa0<<"  "<<L0<<std::endl;
  //     std::cout<<unit_tensor_space.GetSubspace(index).Str()<<std::endl;
  //     std::vector<double>rmes(it->second);
  //     for(int i=0; i<rmes.size(); ++i)
  //       {
  //         std::cout<<rmes[i]<<std::endl;
  //       }
  //   }



}
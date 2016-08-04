/****************************************************************
  upcoupling_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
****************************************************************/
#include <fstream>

#include "cppformat/format.h"

#include "am/am.h"
#include "am/wigner_gsl.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/import_interaction.h"
#include "u3shell/relative_operator.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/moshinsky.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/upcoupling.h"
#include "u3shell/u3st_scheme.h"

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  int Nmax=14;
  int Jmax=Nmax+2;
  int J0=0;
  int g0=0;
	int T0=0;
  basis::RelativeSpaceLSJT relative_space(Nmax, Jmax);
  basis::RelativeSectorsLSJT relative_sectors(relative_space, J0,T0, g0);
  std::vector<Eigen::MatrixXd> sector_vector;
  u3shell::RelativeRMEsU3ST relative_rme_map;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  
  std::string interaction_file;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // JISP16
	// std::string interaction_file="data/Vrel_JISP16_bare_Jmax4.hw20";
  // std::vector<Eigen::MatrixXd> sector_vector=u3shell::ImportInteraction(interaction_file, space, sectors, "JISP");
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// Identity test 
  std::cout<<"Identity test"<<std::endl;
  Jmax=Nmax+2;
  interaction_file="NONE";
  sector_vector=u3shell::ImportInteraction(interaction_file, relative_space, relative_sectors, "Identity");
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  relative_rme_map.clear();
  rme_nlst_map.clear();
  u3shell::UpcouplingNLST(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);
  
  std::cout<<"UpcouplingNLST"<<std::endl;
  for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
    {
      int L0, S0, L,S,T, Lp, Sp, Tp;
      u3shell::RelativeSubspaceLabelsNLST bra, ket;
      std::tie(L0,S0,bra,ket)=it->first;
      std::tie(L,S,T)=ket;
      std::tie(Lp,Sp,Tp)=bra;
      Eigen::MatrixXd sectorNLST=it->second;
      if(fabs(sectorNLST.sum())>10e-10)
	      std::cout<<fmt::format("{} {} ({},{},{}) ({},{},{})", L0,S0,Lp,Sp,Tp,L,S,T)<<std::endl<<sectorNLST<<std::endl;
    }

  u3shell::Upcoupling(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);
  std::cout<<"UpcouplingU3ST"<<std::endl;
  for (auto it=relative_rme_map.begin(); it!= relative_rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST labels;
      int kappa0;
      std::tie(labels,kappa0)=it->first;
      double coef=it->second;
      if (fabs(coef)>10e-13)
        std::cout<<labels.Str()<<"  "<<kappa0<<std::endl<<it->second<<std::endl<<std::endl;
    }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"Kinetic Energy"<<std::endl;
  //// Kinetic eneryg test 
  interaction_file="NONE";
  sector_vector=u3shell::ImportInteraction(interaction_file, relative_space, relative_sectors, "Kinetic");
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // for(int i=0; i<sector_vector.size(); ++i)
  // 	std::cout<<sector_vector[i]<<std::endl;

  relative_rme_map.clear();
  rme_nlst_map.clear();

  u3shell::UpcouplingNLST(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);
  // std::cout<<"UpcouplingNLST"<<std::endl;
  for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
    {
      int L0, S0, L,S,T, Lp, Sp, Tp;
      u3shell::RelativeSubspaceLabelsNLST bra, ket;
      std::tie(L0,S0,bra,ket)=it->first;
      std::tie(L,S,T)=ket;
      std::tie(Lp,Sp,Tp)=bra;
      Eigen::MatrixXd sectorNLST=it->second;
      // if(fabs(sectorNLST.sum())>10e-8)
	     //  std::cout<<fmt::format("{} {} ({},{},{}) ({},{},{})", L0,S0,Lp,Sp,Tp,L,S,T)<<std::endl<<it->second<<std::endl;
    }

  u3shell::Upcoupling(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);
  // std::cout<<"UpcouplingU3ST"<<std::endl;
  std::ostringstream os;
  WriteRelativeOperatorU3ST(os,relative_rme_map);
  std::ofstream stream;
  stream.open("Trel_upcouled");
  stream<<os.str();
  // std::cout<<os.str();
  stream.close();

  std::ifstream is("Trel_upcouled");
  if(!is)
    std::cout<<"Didn't open"<<std::endl;
  u3shell::RelativeRMEsU3ST k2_relative_rmes;
  u3shell::ReadRelativeOperatorU3ST(is, k2_relative_rmes);


  for (auto it=k2_relative_rmes.begin(); it!= k2_relative_rmes.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST labels;
      int kappa0;
      std::tie(labels,kappa0)=it->first;
      u3::SU3(labels.x0());
      u3shell::RelativeStateLabelsU3ST kett(labels.ket());
      u3shell::RelativeStateLabelsU3ST brat(labels.bra());
      double coefout=it->second;
      // double RelativeKineticEnergyOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)

      if (fabs(coefout)>10e-8)
      {
        double Trme=RelativeKineticEnergyOperator(brat,kett);
        std::cout<<labels.Str()
        <<"  "<<kappa0
        <<std::endl
        <<it->second
  			<<"  "
  			<<Trme
  			<<std::endl;     
      }

    }

//   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   // JISP16
//   std::cout<<"JISP16"<<std::endl;
//   interaction_file="data/Vrel_JISP16_bare_Jmax4.hw20";
//   sector_vector=u3shell::ImportInteraction(interaction_file, relative_space, relative_sectors, "JISP");
//   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   // for(int i=0; i<sector_vector.size(); ++i)
//   //  std::cout<<sector_vector[i]<<std::endl;

//   relative_rme_map.clear();
//   rme_nlst_map.clear();

//   u3shell::UpcouplingNLST(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,rme_nlst_map);
// //  std::cout<<"UpcouplingNLST"<<std::endl;
//   for(auto it=rme_nlst_map.begin(); it!=rme_nlst_map.end(); ++it)
//     {
//       int L0, S0, L,S,T, Lp, Sp, Tp;
//       u3shell::RelativeSubspaceLabelsNLST bra, ket;
//       std::tie(L0,S0,bra,ket)=it->first;
//       std::tie(L,S,T)=ket;
//       std::tie(Lp,Sp,Tp)=bra;
//       Eigen::MatrixXd sectorNLST=it->second;
//       //if(fabs(sectorNLST.sum())>10e-8)
//       //  std::cout<<fmt::format("{} {} ({},{},{}) ({},{},{})", L0,S0,Lp,Sp,Tp,L,S,T)<<std::endl<<it->second<<std::endl;
//     }

//   u3shell::Upcoupling(relative_space,relative_sectors,sector_vector,J0,g0,T0,Nmax,relative_rme_map);
// //  std::cout<<"UpcouplingU3ST"<<std::endl;
//   std::map<std::tuple<u3::SU3,HalfInt,HalfInt>,double> distribution;
//   for (auto it=relative_rme_map.begin(); it!= relative_rme_map.end(); ++it)
//     {
//       u3shell::RelativeUnitTensorLabelsU3ST labels;
//       int kappa0;
//       std::tie(labels,kappa0)=it->first;
//       u3::SU3 x0(labels.x0());
//       HalfInt S0=labels.S0();
//       HalfInt T0=labels.T0();
//       u3shell::RelativeStateLabelsU3ST kett(labels.ket());
//       u3shell::RelativeStateLabelsU3ST brat(labels.bra());
//       double coefout=it->second;
//       // double RelativeKineticEnergyOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket)

//       // if ((fabs(coefout)>1)&&(fabs(coefout)<10e10))
//       // {
//       //   std::cout<<labels.Str()<<"  "<<kappa0<<std::endl<<it->second<<std::endl;     
//       // }
//       distribution[std::tuple<u3::SU3,HalfInt,HalfInt>(x0,S0,T0)]+=sqr(coefout);
//     }
//     for(auto it=distribution.begin(); it!=distribution.end(); ++it)
//       { 
//         u3::SU3 x0;
//         HalfInt S0,T0;
//         std::tie(x0,S0,T0)=it->first;
//         std::cout<<fmt::format("[{} {} {}] {}",x0.Str(),S0,T0,it->second)<<std::endl;
//       }


}
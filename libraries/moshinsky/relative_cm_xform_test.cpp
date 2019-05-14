/****************************************************************
  relative_to_twobody.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
****************************************************************/
#include <fstream>
#include <iostream>
#include "fmt/format.h"
#include "basis/lsjt_operator.h"

#include "am/am.h"
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "am/wigner_gsl.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/relative_operator.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/upcoupling.h"
#include "moshinsky/relative_cm_xform.h"

extern double zero_threshold;

typedef std::tuple<int,int,HalfInt,HalfInt, HalfInt> RelativeStateLabelsNLSJT;
typedef std::tuple<RelativeStateLabelsNLSJT,RelativeStateLabelsNLSJT>RelativeBraketNLSJT;

typedef std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> IndexedRelativeUnitTensorLabelsU3ST;
// typedef  std::map<u3shell::RelativeUnitTensorLabelsU3ST, u3shell::TwoBodyUnitTensorCoefficientsU3ST> TwoBodyExpansionMap;

// typedef std::map<u3shell::RelativeCMUnitTensorLabelsU3ST,double> RelativeCMUnitTensorCache;
typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt, HalfInt> RelativeCMLSJTLabels;
typedef std::tuple<int,HalfInt,RelativeCMLSJTLabels,RelativeCMLSJTLabels> RelativeCMLSJTBraket;


void BranchNLSJT(
	std::map<u3shell::RelativeBraketNLST,double>& rel_unit_tensors_lst,
	std::map<RelativeBraketNLSJT,double>& rel_unit_tensors_lsjt,
	HalfInt J0
	)
{
	int N,L,Np,Lp,L0;
	HalfInt Sp,Tp,S,T,T0,S0,J;
	u3shell::RelativeStateLabelsNLST bra_lst,ket_lst;
	for(auto it=rel_unit_tensors_lst.begin(); it!=rel_unit_tensors_lst.end(); ++it)
		{
			std::tie(L0,S0,T0,bra_lst,ket_lst)=it->first;
			double rme=it->second;
			std::tie(N,L,S,T)=ket_lst;
			std::tie(Np,Lp,Sp,Tp)=bra_lst;
			for(HalfInt Jp=abs(Lp-Sp); Jp<=(Lp+Sp); ++Jp)
				for(HalfInt J=abs(L-S); J<=(L+S); ++J)
					{
						if(not am::AllowedTriangle(J,J0,Jp))
							continue;
						RelativeStateLabelsNLSJT bra_lsjt(Np,Lp,Sp,Jp,Tp);
						RelativeStateLabelsNLSJT ket_lsjt(N,L,S,J,T);
						RelativeBraketNLSJT braket(bra_lsjt,ket_lsjt);
						rel_unit_tensors_lsjt[braket]+=am::Unitary9J(L,S,J,L0,S0,J0,Lp,Sp,Jp)
																						*parity((N-L)/2+(Np-Lp)/2)*rme;
					}
		}
}

void
CMBranchLSJT(
    int Nmax,int J0, const std::map<u3shell::RelativeCMBraketNLST,double>& rel_cm_lst_map,
    std::map<RelativeCMLSJTBraket,double>& relative_cm_lsjt_map
 )
// Branching from SO(3)xSU(2) RelxCM to SU(2) RelxCM
{
  int eta,etap,eta_cm, Lr,L_cm,Lrp,L,Lp,L0;
  HalfInt Sp,Tp,S0,T0,S,T;
  u3shell::RelativeCMStateLabelsNLST bra,ket;
  for(auto it=rel_cm_lst_map.begin(); it!=rel_cm_lst_map.end(); ++it)
    {
      if(fabs(it->second)>zero_threshold)
        {
          std::tie(L0,S0,T0,bra,ket)=it->first;
          std::tie(etap,Lrp,eta_cm,L_cm,Lp,Sp,Tp)=bra;
          std::tie(eta,Lr,eta_cm,L_cm,L,S,T)=ket;
          for(HalfInt Jp=abs(Lp-Sp); Jp<=(Lp+Sp); ++Jp)
            for(HalfInt J=abs(L-S); J<=(L+S); ++J)
              {
                double coef=am::Unitary9J(L,S,J,L0,S0,J0,Lp,Sp,Jp);
                RelativeCMLSJTLabels braj(etap,Lrp,eta_cm,L_cm,Lp,Sp,Jp,Tp);
                RelativeCMLSJTLabels ketj(eta,Lr,eta_cm,L_cm,L,S,J,T);
                RelativeCMLSJTBraket braketj(J0,T0,braj,ketj);
                relative_cm_lsjt_map[braketj]+=parity((eta-Lr)/2+(etap-Lrp)/2)*coef*(it->second);
              }
        }
    }

}

  void WriteRelativeOperatorParametersLSJT(
      std::ostream& os,
      HalfInt J0, int T0_min, int T0_max, int g0, int Nmax, int Jmax
      )
  {
    int version = 1;
    os
      << "# RELATIVE LSJT" << std::endl
      << "#   version" << std::endl
      << "#   J0 g0 T0_min T0_max symmetry_phase_mode  [P0=(-)^g0]" << std::endl
      << "#   Nmax Jmax" << std::endl
      << "#   T0   N' L' S' J' T'   N L S J T   JT-RME" << std::endl
      << " " << version << std::endl
      << " " << J0 << " " << g0
      << " " << T0_min << " " << T0_max
      << " " << 0 << std::endl
      << " " << Nmax << " " << Jmax << std::endl;
  }

  void WriteRelativeOperatorComponentLSJT(
      std::ostream& os,
      int T0, int Nmax, int Jmax,
			std::map<RelativeBraketNLSJT,double>& rel_tensors_lsjt
    )
  {
	  basis::OperatorLabelsJT operator_labels;
	  operator_labels.J0=0;
	  operator_labels.g0=0;
	  operator_labels.symmetry_phase_mode=basis::SymmetryPhaseMode::kHermitian;
	  operator_labels.T0_min=0;
	  operator_labels.T0_max=0;

  	basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);

  // populate operator containers
	  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
	  std::array<basis::OperatorBlocks<double>,3> relative_component_matrices; //DUMMY TO USE MARK CODE
	  basis::ConstructIdentityOperatorRelativeLSJT(
	      operator_labels,
	      relative_space,
	      relative_component_sectors,
	      relative_component_matrices
    	);
     const basis::RelativeSectorsLSJT& sectors=relative_component_sectors[0];

    // iterate over sectors
    for (int sector_index = 0; sector_index < sectors.size(); ++sector_index)
      {
        // extract sector
				basis::RelativeSectorsLSJT::SectorType sector = sectors.GetSector(sector_index);
				basis::RelativeSectorsLSJT::SubspaceType bra_subspace = sector.bra_subspace();
				basis::RelativeSectorsLSJT::SubspaceType ket_subspace = sector.ket_subspace();

        // verify that sector is canonical
        //
        // This is a check that the caller's sector construction
        // followed the specification that only "upper triangle"
        // sectors are stored.
        assert(sector.bra_subspace_index()<=sector.ket_subspace_index());

				// iterate over matrix elements
				for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
				  for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
				    {
              // diagonal sector: restrict to upper triangle
              if (sector.IsDiagonal())
                if (!(bra_index<=ket_index))
                  continue;

              // define states
				      const basis::RelativeStateLSJT bra(bra_subspace,bra_index);
				      const basis::RelativeStateLSJT ket(ket_subspace,ket_index);
				      RelativeStateLabelsNLSJT bra2(bra.N(),bra.L(),bra.S(),bra.J(),bra.T());
				      RelativeStateLabelsNLSJT ket2(ket.N(),ket.L(),ket.S(),ket.J(),ket.T());
				      RelativeBraketNLSJT braket1(bra2,ket2);
				      RelativeBraketNLSJT braket2(ket2,bra2);
					    const int width = 3;
					    const int precision = 14;  // less than 16 to provide some roundoff and avoid ugliness on doubles
					    double matrix_element=0.0;
					    if(rel_tensors_lsjt.count(braket1)!=0)
					    	matrix_element=rel_tensors_lsjt[braket1];
					    if(rel_tensors_lsjt.count(braket2)!=0)
					    	matrix_element=rel_tensors_lsjt[braket2];

			        os << std::setprecision(precision);
						  os
								<< " " << std::setw(width) << T0
								<< " " << "  "
								<< " " << std::setw(width) << bra.N()
								<< " " << std::setw(width) << bra.L()
								<< " " << std::setw(width) << bra.S()
								<< " " << std::setw(width) << bra.J()
								<< " " << std::setw(width) << bra.T()
								<< " " << "  "
								<< " " << std::setw(width) << ket.N()
								<< " " << std::setw(width) << ket.L()
								<< " " << std::setw(width) << ket.S()
								<< " " << std::setw(width) << ket.J()
								<< " " << std::setw(width) << ket.T()
								<< " " << "  "
								<< " " << std::showpoint << std::scientific << matrix_element
								<< std::endl;
					  }
			}
	}

void
ReadRelativeCMOperatorNLSJT(std::istream& is,
		std::map<RelativeCMLSJTBraket,double>& relative_cm_lsjt_map,
		int J0
	)
{
	int Nrp,Lrp,Ncmp,Lcmp,Lp,Nr,Lr,Ncm,Lcm,L,g,gp;
	int Sp,Jp,Tp,S,J,T,T0;
	double rme;
	std::string line;

	while(std::getline(is,line))
		{
			std::istringstream(line)>>T0>>Nrp>>Lrp>>Ncmp>>Lcmp>>Lp>>Sp>>Jp>>Tp>>gp
			>>Nr>>Lr>>Ncm>>Lcm>>L>>S>>J>>T>>g>>rme;


			if(fabs(rme)>zero_threshold)
				{
					// std::cout<<fmt::format(" {} {} {} {} {} {} {} {}  {} {} {} {} {} {} {} {}  {} ",
					// 	Nrp,Lrp,Ncm,Lcm,Lp,Sp,Jp,Tp,Nr,Lr,Ncm,Lcm,L,S,J,T,rme)<<std::endl;
					assert(Ncm==Ncmp);
					assert(Lcm==Lcmp);
					RelativeCMLSJTLabels bra(Nrp,Lrp,Ncm,Lcm,Lp,Sp,Jp,Tp);
					RelativeCMLSJTLabels ket(Nr,Lr,Ncm,Lcm,L,S,J,T);
					RelativeCMLSJTBraket braket;
					if(bra>ket)
						braket=RelativeCMLSJTBraket(J0,T0,ket,bra);
					else
						braket=RelativeCMLSJTBraket(J0,T0,bra,ket);
					relative_cm_lsjt_map[braket]=rme;
				}
		}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// In order to compare with shell relative-cm transformation, the code must be run once, then the
// shell code using the input file generated here, then run this code again with same input tensor.
////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  int Nmax=4;
  int Jmax=4;
  int J0=0,Sp,Tp,S,T;
 	int Np,N,lambda,mu,S0,T00;
 	std::cin >>Np>>Sp>>Tp>>N>>S>>T>>lambda>>mu>>S0>>T00>>Nmax;

  // Setting up tensor
  u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
  u3shell::RelativeStateLabelsU3ST ket(N,S,T);
  u3::SU3 x0(lambda,mu);
  u3shell::RelativeUnitTensorLabelsU3ST relative_tensor(x0,S0,T00,bra,ket);

  // // Can be used to generate list of tensor labels to be fed into input stream for this program
  // // std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors_out;
  // //
  // u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensors_out);
  // for(auto tensor : relative_unit_tensors_out)
  //  std::cout<<tensor.bra().eta()<<" "<<tensor.bra().S()<<" "<<tensor.bra().T()
  //  <<" "<<tensor.ket().eta()<<" "<<tensor.ket().S()<<" "<<tensor.ket().T()
  //  <<" "<<tensor.x0().lambda()<<" "<<tensor.x0().mu()<<" "<<tensor.S0()<<" "<<tensor.T0()
  //  <<" "<<Nmax<<std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Branching unit tensor and writing to input file for shell moshinsky transformation
  // for comparison of J-branched rne's of tensor on relative-cm space
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::string interaction_file="/Users/annamccoy/projects/spncci/libraries/moshinsky/test/symmunit_Nmax04_rel.dat";

  std::map<u3shell::RelativeBraketNLST,double> rel_tensors_lst;
 	u3shell::BranchNLST(relative_tensor, rel_tensors_lst);
	std::cout<<"Branch NLST"<<std::endl;
 	for(auto it=rel_tensors_lst.begin(); it!=rel_tensors_lst.end(); ++it)
 		{
 			int Np,Lp,N,L,L0;
 			HalfInt Sp,Tp,S,T,S0,T0;
 			u3shell::RelativeStateLabelsNLST bra,ket;
 			std::tie(L0,S0,T0,bra,ket)=it->first;
 			double rme=it->second;
 			std::tie(Np,Lp,Sp,Tp)=bra;
 			std::tie(N,L,S,T)=ket;
 			std::cout<<fmt::format("{} {} {}  {} {} {} {}   {} {} {} {}   {}",L0,S0,T0,Np,Lp,Sp,Tp,N,L,S,T,rme)<<std::endl;
 		}

 	std::cout<<std::endl<<"Branch NLSJT"<<std::endl;
	std::map<RelativeBraketNLSJT,double> rel_tensors_lsjt;
 	BranchNLSJT(rel_tensors_lst,rel_tensors_lsjt,J0);
 	for(auto it=rel_tensors_lsjt.begin(); it!=rel_tensors_lsjt.end(); ++it)
 		{
 			int Np,Lp,N,L;
 			HalfInt Sp,Tp,Jp,S,T,J;
 			RelativeStateLabelsNLSJT bra,ket;
 			std::tie(bra,ket)=it->first;
 			double rme=it->second;
 			std::tie(Np,Lp,Sp,Jp,Tp)=bra;
 			std::tie(N,L,S,J,T)=ket;
 			std::cout<<fmt::format("{} {} {} {} {}  {} {} {} {} {}  {}", Np,Lp,Sp,Jp,Tp,N,L,S,J,T,rme)<<std::endl;
 		}

	// Writing to shell input file
 	int T0_min=0, T0_max=0,  g0=0, T0=0;
 	std::ofstream os(interaction_file.c_str());
 	WriteRelativeOperatorParametersLSJT(os,J0, T0_min, T0_max, g0, Nmax, Jmax);
 	WriteRelativeOperatorComponentLSJT(os,T0,Nmax,Jmax,rel_tensors_lsjt);
 	os.close();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Checking upcoupling from interaction file generated from unit tensor
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Defining containers for reading in interaction
  basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Jmax);
    // basis::OperatorLabelsJT operator_labels;
  basis::RelativeOperatorParametersLSJT operator_labels;
  std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors;
  std::array<basis::OperatorBlocks<double>,3> isospin_component_matrices;
  //Read interaction and store sector information in relative_component_sectors
  // and matrix elements in relative_component_matrices
  basis::ReadRelativeOperatorLSJT(
    interaction_file,relative_space_lsjt,operator_labels,
    isospin_component_sectors,isospin_component_matrices, true
    );

  // Extract out T0=0 sectors and matrix elements
  const basis::OperatorBlocks<double>& relative_matrices_lsjt=isospin_component_matrices[0];
  const basis::RelativeSectorsLSJT& relative_sectors_lsjt=isospin_component_sectors[0];

  // upcouple to LST
  std::cout<<"Upcoupling to NLST"<<std::endl;
  std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  u3shell::UpcouplingNLST(relative_space_lsjt,relative_sectors_lsjt,relative_matrices_lsjt,int(J0),g0,T0,Nmax,rme_nlst_map);

  // Upcouple to U(3) level
  u3shell::RelativeRMEsU3ST rme_map;
  std::cout<<"Upcoupling to U3ST"<<std::endl;
  u3shell::UpcouplingU3ST(rme_nlst_map, Nmax, rme_map);
  // Should be a single unit tensor with rme=1
  for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST op_labels;
      int kappa0,L0;
      std::tie(op_labels, kappa0,L0)=it->first;
      double rme=it->second;
      if(fabs(rme)>zero_threshold)
        std::cout<<fmt::format("{} {} {}   {}",op_labels.Str(), kappa0,L0,rme)<<std::endl;
    }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Coupling to center of mass at LST level and branching to check against
  // output from shell code
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::map<u3shell::RelativeCMBraketNLST,double> rel_cm_lst_map;
	u3shell::RelativeToCMLST(Nmax, rel_tensors_lst,rel_cm_lst_map);

  std::map<RelativeCMLSJTBraket,double> relative_cm_lsjt_map;
	CMBranchLSJT(Nmax,int(J0), rel_cm_lst_map,relative_cm_lsjt_map);

	// Comparison
  std::string rel_cm_file="/Users/annamccoy/projects/spncci/libraries/moshinsky/test/symmunit_Nmax04_rcmlsjt.dat";
  std::ifstream is(rel_cm_file.c_str());

  std::map<RelativeCMLSJTBraket,double> relative_cm_lsjt_map_2;
	ReadRelativeCMOperatorNLSJT(is, relative_cm_lsjt_map_2,int(J0));

	for(auto it=relative_cm_lsjt_map.begin(); it!=relative_cm_lsjt_map.end(); it++)
	{
		HalfInt T0,J0;
		int Nrp,Lrp,Nr,Lr,Ncm,Lcm,Lp,L;
		HalfInt Sp,Tp,Jp,S,T,J;
		RelativeCMLSJTLabels bra, ket;
		std::tie(J0,T0,bra,ket)=it->first;
		double rme=it->second;
		std::tie(Nrp,Lrp,Ncm,Lcm,Lp,Sp,Jp,Tp)=bra;
		std::tie(Nr,Lr,Ncm,Lcm,L,S,J,T)=ket;
		if((Nr+Ncm)>Nmax)
			continue;
		if((Nrp+Ncm)>Nmax)
			continue;
		double rme2=relative_cm_lsjt_map_2[it->first];
		if(fabs(rme)>zero_threshold)
			std::cout<<fmt::format("{} {} {} {} {} {} {} {}   {} {} {} {} {} {} {} {}   {:12f} {:12f} {:12f}",
			Nrp,Lrp,Ncm,Lcm,Lp,Sp,Jp,Tp,Nr,Lr,Ncm,Lcm,L,S,J,T,rme,rme2,fabs(rme-rme2))<<std::endl;

	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////  U3ST
////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  relative_unit_tensors.push_back(relative_tensor);

	u3shell::RelativeCMExpansion rel_cm_u3st_map;
	u3shell::RelativeUnitTensorToRelativeCMUnitTensorU3ST(Nmax, relative_unit_tensors, rel_cm_u3st_map);

  for(auto it=rel_cm_u3st_map.begin(); it!=rel_cm_u3st_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST tensor=it->first;
      u3shell::RelativeCMUnitTensorCache expansion=it->second;
      std::cout<<tensor.Str()<<std::endl;

      for(auto it2=expansion.begin(); it2!=expansion.end(); ++it2)
      {

          u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st=it2->first;
          double coef=it2->second;
          if(fabs(coef)>zero_threshold)
            std::cout<<"  "<<braket_u3st.Str()<<"  "<<coef<<std::endl;
      }
    }

  // //Transformation to relative-cm at U(3) level, currently not working
  // //
  // RelativeCMUnitTensorCache expansion;
  // UpcoupleCMU3ST(rel_cm_lst_map,expansion);
  // for(auto it2=expansion.begin(); it2!=expansion.end(); ++it2)
  // {

  //     u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st=it2->first;
  //     double coef=it2->second;
  //     // if(fabs(coef)>zero_threshold)
  //       std::cout<<"  "<<braket_u3st.Str()<<"  "<<coef<<std::endl;
  // }





}

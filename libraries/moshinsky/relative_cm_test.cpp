/****************************************************************
  relative_to_twobody.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/15/16 (aem,mac): Created.
****************************************************************/
#include <fstream>
#include <iostream>
#include "cppformat/format.h"
#include "basis/lsjt_operator.h"

#include "am/am.h"
#include "am/wigner_gsl.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/import_interaction.h"
#include "u3shell/relative_operator.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/moshinsky.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/upcoupling.h"

typedef std::tuple<int,int,HalfInt,HalfInt> RelativeStateLabelsNLST;
typedef std::tuple<int, HalfInt,HalfInt,RelativeStateLabelsNLST,RelativeStateLabelsNLST>RelativeBraketNLST;

typedef std::tuple<int,int,HalfInt,HalfInt, HalfInt> RelativeStateLabelsNLSJT;
typedef std::tuple<RelativeStateLabelsNLSJT,RelativeStateLabelsNLSJT>RelativeBraketNLSJT;

typedef std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> IndexedRelativeUnitTensorLabelsU3ST;
typedef  std::map<u3shell::RelativeUnitTensorLabelsU3ST, u3shell::TwoBodyUnitTensorCoefficientsU3ST> TwoBodyExpansionMap;
typedef std::tuple<u3shell::TwoBodyUnitTensorLabelsU3ST, int, int>IndexTwoBodyTensorLabelsU3ST;
typedef std::map<IndexTwoBodyTensorLabelsU3ST,double>IndexedTwoBodyTensorRMEsU3ST;

void
BranchNLST(
  const u3shell::RelativeUnitTensorLabelsU3ST& relative_unit_tensor,
  std::map<RelativeBraketNLST,double>& branched_rel_unit_tensors
  )
{
  int N0, etap,eta,eta_cm;
  HalfInt S0,T0,Sp,Tp,S,T;
  u3::SU3 x0;
  std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=relative_unit_tensor.FlatKey();
  MultiplicityTagged<int>::vector L0_set=u3::BranchingSO3(x0);
  for(auto Lk0 : L0_set)
    {
      int L0=Lk0.irrep;
      int kappa0_max=Lk0.tag;
      for(int kappa0=1; kappa0<=kappa0_max; ++kappa0)
        for(int Lp=etap%2; Lp<=etap; Lp+=2)
          for(int L=eta%2; L<=eta; L+=2)
            {
              if(not am::AllowedTriangle(L,L0,Lp))
                continue;
              RelativeStateLabelsNLST bra_nlst(etap,Lp,Sp,Tp);
              RelativeStateLabelsNLST ket_nlst(eta,L,S,T);
              RelativeBraketNLST braket_nlst(L0,S0,T0,bra_nlst,ket_nlst);
              branched_rel_unit_tensors[braket_nlst]
                +=u3::W(u3::SU3(eta,0),1,L,x0,kappa0,L0,u3::SU3(etap,0),1,Lp,1);
            }
    }
}

void BranchNLSJT(
	std::map<RelativeBraketNLST,double>& rel_unit_tensors_lst,
	std::map<RelativeBraketNLSJT,double>& rel_unit_tensors_lsjt,
	HalfInt J0
	)
{
	int N,L,Np,Lp,L0;
	HalfInt Sp,Tp,S,T,T0,S0,J;
	RelativeStateLabelsNLST bra_lst,ket_lst;
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

typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt> RelativeCMStateLabelsNLST;
typedef std::tuple<int,HalfInt,HalfInt,RelativeCMStateLabelsNLST,RelativeCMStateLabelsNLST> RelativeCMBraketNLST;

void RelativeToCMLST(
  int Nmax, 
  const std::map<RelativeBraketNLST,double>& branched_rel_unit_tensors,
  std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map)
// Coupling on center of mass as SO(3) level
{
  for(auto it=branched_rel_unit_tensors.begin(); it!=branched_rel_unit_tensors.end(); ++it)
    {
      int L0,Lr,Lrp,N,Np;
      HalfInt S,T,Sp,Tp,S0,T0;
      RelativeStateLabelsNLST ket_rel,bra_rel;
      RelativeBraketNLST rel_tensor=it->first;
      double rme_rel=it->second;
      std::tie(L0,S0,T0,bra_rel,ket_rel)=rel_tensor;
      std::tie(Np,Lrp,Sp,Tp)=bra_rel;
      std::tie(N,Lr,S,T)=ket_rel;
      int Ncm_max=Nmax-std::max(Np,N);
      for(int Ncm=0; Ncm<=Ncm_max; ++Ncm)
        for(int Lcm=Ncm%2; Lcm<=Ncm; Lcm+=2)
          for(int Lp=abs(Lrp-Lcm); Lp<=(Lrp+Lcm); ++Lp)
            for(int L=abs(Lr-Lcm); L<=(Lr+Lcm); ++L)
                {
                  RelativeCMStateLabelsNLST bra_rel_cm(Np,Lrp,Ncm,Lcm,Lp,Sp,Tp);
                  RelativeCMStateLabelsNLST ket_rel_cm(N,Lr,Ncm,Lcm,L,S,T);
                  RelativeCMBraketNLST braket_rel_cm(L0,S0,T0,bra_rel_cm,ket_rel_cm);
                  rel_cm_lst_map[braket_rel_cm]
                    =am::Unitary9J(L0,0,L0,Lr,Lcm,L,Lrp,Lcm,Lp)
                    *parity(Lr+Lrp+L+Lp)*rme_rel;
                }
    }
}
typedef std::tuple<int,int,int,int,int,HalfInt,HalfInt, HalfInt> RelativeCMLSJTLabels;
typedef std::tuple<int,HalfInt,RelativeCMLSJTLabels,RelativeCMLSJTLabels> RelativeCMLSJTBraket;

void
CMBranchLSJT(
    int Nmax,int J0, const std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map,
    std::map<RelativeCMLSJTBraket,double>& relative_cm_lsjt_map
 )
// Branching from SO(3)xSU(2) RelxCM to SU(2) RelxCM
{
  int eta,etap,eta_cm, Lr,L_cm,Lrp,L,Lp,L0;
  HalfInt Sp,Tp,S0,T0,S,T;
  RelativeCMStateLabelsNLST bra,ket;
  for(auto it=rel_cm_lst_map.begin(); it!=rel_cm_lst_map.end(); ++it)
    {
      if(fabs(it->second)>10e-10)
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

typedef std::map<u3shell::RelativeCMUnitTensorLabelsU3ST,double> RelativeCMUnitTensorExpansion;

void UpcoupleCMU3ST(
  std::map<RelativeCMBraketNLST,double>& rel_cm_lst_map,
  RelativeCMUnitTensorExpansion& rel_cm_u3st_map
  )
{
  int Nrp,Lrp,Ncm,Lcm,Nr,Lr,Lp,L,L0;
  HalfInt Sp,Tp,S,T,S0,T0;
  RelativeCMStateLabelsNLST bra_cm,ket_cm;
  for(auto it=rel_cm_lst_map.begin(); it!=rel_cm_lst_map.end(); ++it)
    {
      std::tie(L0,S0,T0,bra_cm,ket_cm)=it->first;
      double rme=it->second;
      std::tie(Nrp,Lrp,Ncm,Lcm,Lp,Sp,Tp)=bra_cm;
      std::tie(Nr,Lr,Ncm,Lcm,L,S,T)=ket_cm;

      MultiplicityTagged<u3::SU3>::vector  x_set=u3::KroneckerProduct(u3::SU3(Nr,0),u3::SU3(Ncm,0));
      MultiplicityTagged<u3::SU3>::vector  xp_set=u3::KroneckerProduct(u3::SU3(Nrp,0),u3::SU3(Ncm,0));
      for(int i=0; i<x_set.size(); ++i)
        {
          u3::SU3 x=x_set[i].irrep;
          int kappa_max=u3::BranchingMultiplicitySO3(x,L);
          for(int ip=0; ip<xp_set.size(); ++ip)
            {
              u3::SU3 xp=xp_set[ip].irrep;
              MultiplicityTagged<u3::SU3>::vector x0_set=u3::KroneckerProduct(xp,u3::Conjugate(x));
              int kappap_max=u3::BranchingMultiplicitySO3(xp,Lp); 
              for(int kappa=1; kappa<=kappa_max; ++kappa)
                for(int kappap=1; kappap<=kappap_max; ++kappap)
                  {
                    for(int i0=0; i0<x0_set.size();  i0++)
                      {
                        u3::SU3 x0=x0_set[i0].irrep;
                        int rho0_max=u3::OuterMultiplicity(x,x0,xp);
                        int kappa0_max=u3::BranchingMultiplicitySO3(x0,L0);
                        u3shell::RelativeCMStateLabelsU3ST bra(Nrp,Ncm,xp,Sp,Tp);
                        u3shell::RelativeCMStateLabelsU3ST ket(Nr,Ncm,x,S,T);

                        for(int kappa0=1; kappa0<=kappa0_max; ++kappa0)
                          for(int rho0=1; rho0<=rho0_max; ++rho0)
                              {
                                u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st(x0,S0,T0,rho0,bra,ket); 
                                rel_cm_u3st_map[braket_u3st]
                                	+=am::dim(Lp)/u3::dim(xp)
                                //u3::dim(x0)*am::dim(Lp)/u3::dim(xp)/am::dim(L0)
                                  *u3::W(u3::SU3(Nrp,0),1,Lrp,u3::SU3(Ncm,0),1,Lcm,xp,kappap,Lp,1)
                                  *u3::W(u3::SU3(Nr,0),1,Lr,u3::SU3(Ncm,0),1,Lcm,x,kappa,L,1)
                                  *u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0)
                                  *rme;
                              }
                      }
                  }
            }
        }
    }
}

void RelativeToCMU3ST(int Nmax,  
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  std::vector<RelativeCMUnitTensorExpansion>& unit_tensor_rel_cm_expansions
  )
// Coupling on center of mass at U(3) level
{
  int N0, Np,N, kappa0,L0;
  HalfInt S0,T0,Sp,Tp,S,T;
  u3::SU3 x0;
  u3shell::RelativeUnitTensorLabelsU3ST tensor;
  for(auto tensor:relative_unit_tensors)
    {
      RelativeCMUnitTensorExpansion relative_cm_u3st_map;
      std::tie(x0,S0,T0,Np,Sp,Tp,N,S,T)=tensor.FlatKey();
      u3::SU3 xr(N,0);
      u3::SU3 xrp(Np,0);
      for(int Ncm=0; Ncm<=Nmax; Ncm++)
        {
          u3::SU3 x_cm(Ncm,0);
          MultiplicityTagged<u3::SU3>::vector x_set=u3::KroneckerProduct(xr,x_cm);
          MultiplicityTagged<u3::SU3>::vector xp_set=u3::KroneckerProduct(xrp,x_cm);
          for(auto ip : xp_set)
            for(auto i: x_set)
              {
                u3::SU3 x(i.irrep);
                u3::SU3 xp(ip.irrep);
                int rho0_max=u3::OuterMultiplicity(x,x0,xp);
                u3shell::RelativeCMStateLabelsU3ST bra(Np,Ncm,xp,Sp,Tp);
                u3shell::RelativeCMStateLabelsU3ST ket(N,Ncm,x,S,T);

                for(int rho0=1; rho0<=rho0_max; ++rho0)
                {
                  u3shell::RelativeCMUnitTensorLabelsU3ST tensor_cm(x0,S0,T0,rho0,bra,ket);
                  relative_cm_u3st_map[tensor_cm]
                  	=//parity(x.lambda()+x.mu()+xp.lambda()+xp.mu()+Np+N)
                  		u3::U(x0,xr,xp,x_cm,xrp,1,1,x,1,rho0);
                }
              }
        }
      unit_tensor_rel_cm_expansions.push_back(relative_cm_u3st_map);
    }
}

void RelativeUnitTensorToRelativeCMUnitTensorU3ST(int Nmax,  
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
  std::map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorExpansion>& unit_relative_cm_map)
// Coupling on center of mass by branching and re-upcoupling
{
  for(int i=0; i<relative_unit_tensors.size(); ++i)
    {
      const u3shell::RelativeUnitTensorLabelsU3ST& tensor=relative_unit_tensors[i];

      std::map<RelativeBraketNLST,double> branched_rel_unit_tensors;
      BranchNLST(tensor,branched_rel_unit_tensors);

      std::map<RelativeCMBraketNLST,double> rel_cm_lst_map;
      RelativeToCMLST(Nmax, branched_rel_unit_tensors,rel_cm_lst_map);
      
      RelativeCMUnitTensorExpansion rel_cm_u3st_map;
      UpcoupleCMU3ST(rel_cm_lst_map,rel_cm_u3st_map);

      unit_relative_cm_map[tensor]=rel_cm_u3st_map;
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
	  std::array<basis::MatrixVector,3> relative_component_matrices; //DUMMY TO USE MARK CODE
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


			if(fabs(rme)>10e-8)
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
////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  int Nmax=4; 
  int Jmax=4; 
  int J0=0,Sp,Tp,S,T;
 	int Np,N,lambda,mu,S0,T00;
 	std::cin >>Np>>Sp>>Tp>>N>>S>>T>>lambda>>mu>>S0>>T00>>Nmax;

  u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
  u3shell::RelativeStateLabelsU3ST ket(N,S,T);
  u3::SU3 x0(lambda,mu);
  u3shell::RelativeUnitTensorLabelsU3ST relative_tensor(x0,S0,T00,bra,ket);

  std::map<RelativeBraketNLST,double> rel_tensors_lst;
 	BranchNLST(relative_tensor, rel_tensors_lst);
	// std::cout<<"Branch NLST"<<std::endl;
 	for(auto it=rel_tensors_lst.begin(); it!=rel_tensors_lst.end(); ++it)
 		{
 			int Np,Lp,N,L,L0;
 			HalfInt Sp,Tp,S,T,S0,T0;
 			RelativeStateLabelsNLST bra,ket;
 			std::tie(L0,S0,T0,bra,ket)=it->first;
 			double rme=it->second;
 			std::tie(Np,Lp,Sp,Tp)=bra;
 			std::tie(N,L,S,T)=ket;
 			std::cout<<fmt::format("{} {} {}  {} {} {} {}   {} {} {} {}   {}",L0,S0,T0,Np,Lp,Sp,Tp,N,L,S,T,rme)<<std::endl;
 		}

 	// std::cout<<std::endl<<"Branch NLSJT"<<std::endl;
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

	// For comparison with Mark
  std::string interaction_file="/Users/annamccoy/projects/spncci/libraries/moshinsky/test/symmunit_Nmax04_rel.dat";

 	int T0_min=0, T0_max=0,  g0=0, T0=0;
 	std::ofstream os(interaction_file.c_str());
 	WriteRelativeOperatorParametersLSJT(os,J0, T0_min, T0_max, g0, Nmax, Jmax);
 	WriteRelativeOperatorComponentLSJT(os,T0,Nmax,Jmax,rel_tensors_lsjt);
 	os.close();

  // Defining containers for reading in interaction
  basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
  basis::OperatorLabelsJT operator_labels;
  std::array<basis::RelativeSectorsLSJT,3> relative_component_sectors;
  std::array<basis::MatrixVector,3> relative_component_matrices;
  //Read interaction and store sector information in relative_component_sectors
  // and matrix elements in relative_component_matrices
  basis::ReadRelativeOperatorLSJT(
    interaction_file,relative_lsjt_space,operator_labels,
    relative_component_sectors,relative_component_matrices, true
    );

  // Extract out T0=0 sectors and matrix elements
  const basis::MatrixVector& sector_vector=relative_component_matrices[0];
  const basis::RelativeSectorsLSJT& relative_lsjt_sectors=relative_component_sectors[0];

  //upcouple to LST
  // std::cout<<"Upcoupling to NLST"<<std::endl;
  // std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> rme_nlst_map;
  // u3shell::UpcouplingNLST(relative_lsjt_space,relative_lsjt_sectors,sector_vector,int(J0),g0,T0,Nmax,rme_nlst_map);

  // // Upcouple to U(3) level
  // u3shell::RelativeRMEsU3ST rme_map;
  // std::cout<<"Upcoupling to U3ST"<<std::endl;
  // u3shell::UpcouplingU3ST(rme_nlst_map, T0, Nmax, rme_map);
  // for(auto it=rme_map.begin(); it!=rme_map.end(); ++it)
  //   {
  //     u3shell::RelativeUnitTensorLabelsU3ST op_labels;
  //     int kappa0,L0;
  //     std::tie(op_labels, kappa0,L0)=it->first;
  //     double rme=it->second;
  //     // double check=u3shell::RelativeKineticEnergyOperator(op_labels.bra(), op_labels.ket());
  //     if(fabs(rme)>10e-10)
  //       std::cout<<fmt::format("{} {} {}   {}",op_labels.Str(), kappa0,L0,rme)<<std::endl;
  //   }

  std::map<RelativeCMBraketNLST,double> rel_cm_lst_map;
	RelativeToCMLST(Nmax, rel_tensors_lst,rel_cm_lst_map);

  std::map<RelativeCMLSJTBraket,double> relative_cm_lsjt_map;
	CMBranchLSJT(Nmax,int(J0), rel_cm_lst_map,relative_cm_lsjt_map);

	
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
			if(fabs(rme)>10e-8)
				std::cout<<fmt::format("{} {} {} {} {} {} {} {}   {} {} {} {} {} {} {} {}   {:12f} {:12f} {:12f}",
				Nrp,Lrp,Ncm,Lcm,Lp,Sp,Jp,Tp,Nr,Lr,Ncm,Lcm,L,S,J,T,rme,rme2,fabs(rme-rme2))<<std::endl;

		}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////  U3ST
////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors;
  relative_unit_tensors.push_back(relative_tensor);

	std::map<u3shell::RelativeUnitTensorLabelsU3ST,RelativeCMUnitTensorExpansion> rel_cm_u3st_map;
	RelativeUnitTensorToRelativeCMUnitTensorU3ST(Nmax, relative_unit_tensors, rel_cm_u3st_map);


  // RelativeCMUnitTensorExpansion expansion;
  // UpcoupleCMU3ST(rel_cm_lst_map,expansion);
  // for(auto it2=expansion.begin(); it2!=expansion.end(); ++it2)
  // {
      
  //     u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st=it2->first;
  //     double coef=it2->second;
  //     // if(fabs(coef)>10e-8)
  //       std::cout<<"  "<<braket_u3st.Str()<<"  "<<coef<<std::endl;
  // }


  for(auto it=rel_cm_u3st_map.begin(); it!=rel_cm_u3st_map.end(); ++it)
    {
      u3shell::RelativeUnitTensorLabelsU3ST tensor=it->first;
      RelativeCMUnitTensorExpansion expansion=it->second;
      std::cout<<tensor.Str()<<std::endl;

      for(auto it2=expansion.begin(); it2!=expansion.end(); ++it2)
      {
          
          u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st=it2->first;
          double coef=it2->second;
          if(fabs(coef)>10e-8)
            std::cout<<"  "<<braket_u3st.Str()<<"  "<<coef<<std::endl;
      }
    }

  // std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensors_out;

  // u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax, relative_unit_tensors_out);
  // for(auto tensor : relative_unit_tensors_out)
  // 	std::cout<<tensor.bra().eta()<<" "<<tensor.bra().S()<<" "<<tensor.bra().T()
  // 	<<" "<<tensor.ket().eta()<<" "<<tensor.ket().S()<<" "<<tensor.ket().T()
  // 	<<" "<<tensor.x0().lambda()<<" "<<tensor.x0().mu()<<" "<<tensor.S0()<<" "<<tensor.T0()
  // 	<<" "<<Nmax<<std::endl;



 //  std::vector<RelativeCMUnitTensorExpansion> unit_tensor_rel_cm_expansions;
	// RelativeToCMU3ST(Nmax, relative_unit_tensors, unit_tensor_rel_cm_expansions);
	// for(int i=0; i<relative_unit_tensors.size(); ++i)
	// 	{
	// 		u3shell::RelativeUnitTensorLabelsU3ST tensor=relative_unit_tensors[i];
	// 		std::cout<<tensor.Str()<<std::endl;
	// 		RelativeCMUnitTensorExpansion& expansion1=unit_tensor_rel_cm_expansions[i];
	// 		RelativeCMUnitTensorExpansion& expansion2=rel_cm_u3st_map[tensor];
	// 		for(auto it=expansion1.begin(); it!=expansion1.end(); ++it)
	// 			{
	// 				double rme1=it->second;
	// 				double rme2=expansion2[it->first];
	// 				if(fabs(rme1-rme2)>10e-8)	
	// 					std::cout<<it->first.Str()<<"  "<<rme1<<"  "<<rme2<<std::endl;
	// 			}
	// 	}


}
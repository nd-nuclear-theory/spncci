/****************************************************************
  moshinsky_coefficient_test.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created to test moshinky.h, moshinsky.cpp.
  5/11/16 (aem,mac): Update namespace.
****************************************************************/
#include <iomanip>
#include <iostream>
#include <fstream>

#include "cppformat/format.h"

#include "sp3rlib/u3coef.h"
#include "moshinsky/moshinsky_xform.h"
#include "u3shell/relative_operator.h"
#include "am/wigner_gsl.h"

void CompairisonTests(int Nmax)
	{
		for (int N=0; N<=Nmax; ++N)
			for(int N1=0; N1<=N; ++N1)
				for(int Nr=0; Nr<=N; ++Nr)
					{
						int N2=N-N1;
						int Ncm=N-Nr;
						auto tb_coupled=u3::KroneckerProduct(u3::SU3(N1,0),u3::SU3(N2,0));
						for(auto a: tb_coupled)
							{
								u3::SU3 x(a.irrep);
								if(u3::OuterMultiplicity(u3::SU3(Nr,0), u3::SU3(Ncm,0),x)>0)
									{
										double mb=u3shell::MoshinskyCoefficient(Nr,Ncm,N1,N2,x);
										std::cout<<fmt::format("({},0) ({},0) ({},0) ({},0) {}  {}",Nr,Ncm,N1,N2,x.Str(),mb)<<std::endl;
									}
							}
					}
	}

void SymmetryTests(int Nmax)
	{
		for (int N1=0; N1<Nmax; ++N1)
			for(int N2=0; N2<Nmax; ++N2)
				for(int Nr=0; Nr<Nmax; ++Nr)
					{
						int Ncm=N1+N2-Nr;
						auto tb_coupled=u3::KroneckerProduct(u3::SU3(N1,0),u3::SU3(N2,0));
						for(auto a: tb_coupled)
							{
								u3::SU3 x(a.irrep);
								if(u3::OuterMultiplicity(u3::SU3(Nr,0), u3::SU3(Ncm,0),x)>0)
									{
										double mb1=u3shell::MoshinskyCoefficient(Nr,Ncm,N1,N2,x);
										double mb2=parity(Nr+N1+N2+x.lambda()+x.mu())*u3shell::MoshinskyCoefficient(Nr,Ncm,N2,N1,x);
										if(fabs(mb1-mb2)>10e-10)
											std::cout<<fmt::format("({} {}; {} | {} {}; {})  {}  {}",Nr,Ncm,x.Str(),N1,N2,x.Str(),mb1,mb2)<<std::endl;
									}
							}
					}
	}

typedef std::tuple<int,int,int,int,int,int,int,int,int> BranchedMoshinsky;
void BranchTest(int Nmax, const std::string& filename)
{
	int N1,L1,N2,L2,Nr,Lr,Ncm,Lcm,L;
	double rme;
	std::map<BranchedMoshinsky,double>moshisky_check;
	std::ifstream is(filename.c_str());
	assert(is);
	std::string line;
	while(std::getline(is,line))
		{
			std::istringstream(line)>>Nr>>Lr>>Ncm>>Lcm>>N1>>L1>>N2>>L2>>L>>rme;
			BranchedMoshinsky key(Nr,Lr,Ncm,Lcm,N1,L1,N2,L2,L);
			moshisky_check[key]=rme;
		}

	std::map<BranchedMoshinsky,double>moshisky_map;
	for (int N=0; N<=Nmax; ++N)
		for(int N1=0; N1<=N; ++N1)
			for(int Nr=0; Nr<=N; ++Nr)
				{
					int N2=N-N1;
					int Ncm=N-Nr;
					MultiplicityTagged<u3::SU3>::vector tb_coupled=u3::KroneckerProduct(u3::SU3(N1,0),u3::SU3(N2,0));
					for(int a=0; a<tb_coupled.size(); ++a)
						{
							u3::SU3 x(tb_coupled[a].irrep);
							if(u3::OuterMultiplicity(u3::SU3(Nr,0), u3::SU3(Ncm,0),x)>0)
								{
									double mb=u3shell::MoshinskyCoefficient(Nr,Ncm,N1,N2,x);
									for(int L1=N1%2; L1<=N1; L1+=2)
										for(int L2=N2%2; L2<=N2; L2+=2)
											for(int Lr=Nr%2; Lr<=Nr; Lr+=2)
												for(int Lcm=Ncm%2; Lcm<=Ncm; Lcm+=2)
													for(int L=abs(L1-L2); L<=L1+L2; ++L)
														{
															if ((L<abs(Lr-Lcm))||(L>(Lr+Lcm)))
																continue;
															BranchedMoshinsky labels(Nr,Lr,Ncm,Lcm,N1,L1,N2,L2,L);
															double mb_branched=0;
															int kappa_max=u3::BranchingMultiplicitySO3(x,L);
															int n1=(N1-L1)/2;
															int n2=(N2-L2)/2;
															int nr=(Nr-Lr)/2;
															int ncm=(Ncm-Lcm)/2;
															for(int kappa=1; kappa<=kappa_max; ++kappa)
																mb_branched+=mb*u3::W(u3::SU3(N1,0),1,L1,u3::SU3(N2,0),1,L2,x,kappa,L,1)
																							 *u3::W(u3::SU3(Nr,0),1,Lr,u3::SU3(Ncm,0),1,Lcm,x,kappa,L,1);

															moshisky_map[labels]+=parity(n1+n2+nr+ncm)*mb_branched;
														}
								}
						}
				}
	for(auto it=moshisky_map.begin(); it!=moshisky_map.end(); ++it)
		{
			std::tie(Nr,Lr,Ncm,Lcm,N1,L1,N2,L2,L)=it->first;
			int n1=(N1-L1)/2;
			int n2=(N2-L1)/2;
			double rme_check=moshisky_check[it->first];
			if(fabs(rme_check-it->second)>10e-10)
				std::cout<<fmt::format("{} {} {} {} {} {} {} {} {}  {} {}", Nr,Lr,Ncm,Lcm,N1,L1,N2,L2,L, it->second, rme_check)<<std::endl;
		}
}

int main(int argc, char **argv)
{
	u3::U3CoefInit();
	// CompairisonTests(6);
	// SymmetryTests(20);
	std::string filename="/Users/annamccoy/projects/spncci/data/moshinsky_bracket_table_Nmax6.out";
	// BranchTest(6,filename);

	// checking so3 coupling rule
	// int Smax=10;
	// for(int S=0; S<=Smax; ++S)
	// 	for(int Sp=0; Sp<=Smax; ++Sp)
	// 		for(int S0=abs(S-Sp); S0<=(S+Sp); ++S0)
	// 			{
	// 				double coef1=am::Unitary9J(S0,0,S0,S,0,S,Sp,0,Sp);
	// 				// double coef2=Hat(Sp)/(1.*Hat(S0)*Hat(S));
	// 				if (fabs(coef1-1)>10e-10)
	// 					std::cout<<fmt::format("{} {} {}  {}",S,S0,Sp,coef1)<<std::endl;
	// 			}




}
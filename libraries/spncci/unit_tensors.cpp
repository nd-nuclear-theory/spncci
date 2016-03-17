/****************************************************************
  unit_tensor.cpp

  Unit tensor algorithms
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/15/16 (aem,mac): Created.
****************************************************************/
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include <map>
#include <sstream>
#include <vector>

#include "am/am.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/vcs.h"
#include "sp3rlib/sp3r.h"
#include "spncci/sp_basis.h"
#include "spncci/unit_tensors.h"

  void UnitSymmetryGenerator(int N1b, int Nmax, std::map< int,std::vector<u3::UnitTensor> >& unit_sym_map)
// Generates a map containing (key, value) pair (N0, operator_labels) of the unit tensors 
// Generated for rbp>=rb.  To get the other half, use conjugation 
  {		
  	for(int N0=0; N0<=Nmax; N0+=2)
  	{
  		std::vector<u3::UnitTensor> sym_vec;
  		for(int rp=N0; rp<=N0+Nmax; rp+=2)
  		{
  			int r=rp-N0;
  			MultiplicityTagged<u3::U3>::vector omega0_set=u3::KroneckerProduct(u3::U3(rp,0,0),u3::U3(0,0,-r));
  			for(int w=0; w<omega0_set.size(); w++)
  			{
  				u3::U3 omega0(omega0_set[w].irrep);
  				for(HalfInt S=0; S<=1; S++)
  					for(HalfInt Sp=0; Sp<=1; Sp++)
  						for(HalfInt T=0; T<=1; T++)
  							for(HalfInt Tp=0; Tp<=1; Tp++)
  								for (HalfInt S0=abs(S-Sp); S0<=(S+Sp); S0++)
  									for (HalfInt T0=abs(T-Tp); T0<=(T+Tp); T0++)
  										sym_vec.push_back(u3::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T));
  								}
  							}
  							unit_sym_map[N0]=sym_vec;
  						}
	} //end function


//	//Generate list of LGI's
// 		spncci::LGIVectorType lgi_vector;	
// 		spncci::GenerateLGIVector(lgi_vector,LGI_filename,Nsigma_0);
//		spncci::GenerateSp3RIrreps(lgi_vector,sigma_irrep_map,truncator);
//  Generate list of LGI pairs in canonical ordering ket then bra. 
//		std::pair<spncci::LGI,spncci::LGI> lgi_pair
//		
//  For given nuclei and Nmax cutoff, generate possible unit tensor symmetries. 
// 		std::map< int,std::vector<u3::UnitTensor>> unit_sym_map;
// 		UnitSymmetryGenerator(N1b, Nmax, unit_sym_map);




	void UnitTensorMatrixGenerator(
	// boson number cutoff
		int Nmax, 
	// a given spncci sector pair  
		std::pair<spncci::LGI,spncci::LGI> lgi_pair,
	// Address to map with list of unit tensor labels with key N0 
		std::map< int,std::vector<u3::UnitTensor>>& unit_sym_map,
	// Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
		std:: map<std::pair<spncci::LGI,spncci::LGI>,std::map< u3::UnitTensorRME,Eigen::MatrixXd>> unit_tensor_rme_map
		)
//
// 
	{
		// initial declarations 		
		u3::U3 omega0;
		HalfInt S0, T0, Sbp, Tbp, Sb, Tb ;
		int rbp,rb;

	  // Extracting LGI labels from pair
		spncci::LGI lgip=lgi_pair.first;
		spncci::LGI  lgi=lgi_pair.second;
		//std::tie (lgip,lgi)=lgi_pair;

	  //Generate Sp(3,R) irreps 
		sp3r::Sp3RSpace irrepp(lgip.sigma,Nmax-lgip.Nex);
		sp3r::Sp3RSpace irrep(lgi.sigma, Nmax-lgi.Nex);

		std::map< u3::U3,std::map<u3::U3,Eigen::MatrixXd> > K_matrix_map;

	  // Calculating K matrices for each sigma in LGI set and storing in map K_matrix_map with key sigma
		for (int k = 0; k<irrep.size();++k)
			{
				std::map<u3::U3,Eigen::MatrixXd> K_map;
				vcs::GenerateKMatrices(irrep, K_map);
				K_matrix_map[lgi.sigma]=K_map;
			}
		for (int kp = 0; kp<irrepp.size();++kp)
			{
				std::map<u3::U3,Eigen::MatrixXd> K_map;
				vcs::GenerateKMatrices(irrep, K_map);
				K_matrix_map[lgip.sigma]=K_map;
			}      	      		
		////////////////////////////////////////////////////////////////////////////////////
		// Looping over omega' and omega subspaces 
		////////////////////////////////////////////////////////////////////////////////////
		// omega subspace
		for(int i=0; i<irrep.size(); i++ )
			{
				sp3r::U3Subspace u3_subspace=irrep.GetSubspace(i);
				u3::U3 omega=u3_subspace.GetSubspaceLabels();
				int dim=u3_subspace.size();
	    	// Looking up and inverting K 
				Eigen::MatrixXd K_inv=K_matrix_map[lgi.sigma][omega].inverse();
				MultiplicityTagged<u3::U3>::vector omega1_set=KroneckerProduct(omega, u3::U3(0,0,-2));
				//omega' subspace 
				for(int ip=0; ip<irrepp.size(); ip++ )
					{
						sp3r::U3Subspace u3_subspacep=irrep.GetSubspace(ip);
						u3::U3 omegap=u3_subspacep.GetSubspaceLabels();
				      // may not need to define
						int dimp=u3_subspacep.size();
						std::vector<u3::UnitTensor> operator_set=unit_sym_map[int(omegap.N()-omega.N())];
						Eigen::MatrixXd Kp=K_matrix_map[lgip.sigma][omegap];
				    	// For A RME
						MultiplicityTagged<u3::U3>::vector omegapp_set=KroneckerProduct(omegap, u3::U3(0,0,-2));

				      // Iterating over the operator labels 
						for (int w=0; w<operator_set.size(); w++)
							{						      		
								std::tie (omega0, S0, T0, rbp, Sbp, Tbp, rb, Sb, Tb) = operator_set[w].Key();
								int rho0_max=OuterMultiplicity(omega.SU3(),omega0.SU3(),omegap.SU3());
								MultiplicityTagged<u3::U3>::vector omega0p_set=KroneckerProduct(omega, u3::U3(2,0,0));
					      		// Iterating over outer multiplicity
								for (int rho0=1; rho0<=rho0_max; rho0++)
									{
						      				// Initializing unit tensor matrix with dim. v' v
										Eigen::MatrixXd unit_tensor_matrix(dimp,dim);

						      				// summing over omega1
										for (int w1=0; w1<omega1_set.size(); w1++)
										{	
											u3::U3 omega1=omega1_set[w1].irrep;
											sp3r::U3Subspace u3_subspace1=irrep.LookUpSubspace(omega1);
						      						// Look up K1 matrix (dim v1, v1)
											Eigen::MatrixXd K1=K_matrix_map[lgi.sigma][omega1];
											int dim1=u3_subspace1.size();
						      						/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								      				// Matrix of B*U coefs with dim v1 and v
								      				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
											Eigen::MatrixXd BU(dim1,dim);
								      				// iterating over (n,rho)
											for (int m=0; m<dim; m++)
											{

													MultiplicityTagged<u3::U3> n_rho=u3_subspace.GetStateLabels(m);
													u3::U3 n(n_rho.irrep);
													for (int m1=0; m1<dim1; m1++)
														{
															MultiplicityTagged<u3::U3> n1_rho1=u3_subspace1.GetStateLabels(m);
															u3::U3 n1(n1_rho1.irrep);
															BU(m1,m)=(vcs::BosonCreationRME(n,n1)
																*u3::U(u3::SU3(2,0),n1.SU3(),omega.SU3(),lgi.sigma.SU3(),n.SU3(),1,n_rho.tag,omega1.SU3(),n1_rho1.tag,1));
														}	
												}								        
								              /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
											Eigen::MatrixXd KBUK=K1*BU*K_inv;

						      						// summing over omega0'
											for (int w0p=0; w0p<omega0p_set.size(); w0p++)
												{
													u3::U3 omega0p=omega0p_set[w0p].irrep;



													int rho0p_max=OuterMultiplicity(omega1.SU3(),omega0p.SU3(),omegap.SU3());
							      							// summing over rho0'
													for (int rho0p=1; rho0p<=rho0p_max; rho0p++)
														{
								      									// Coefficients 
															double coef1=
															(
																u3::U(
																	omega0.SU3(),u3::SU3(2,0),omegap.SU3(), omega1.SU3(),
																	omega0p.SU3(),1,rho0p,omega.SU3(),1,rho0
																	)
																*u3::U(
																	u3::SU3(rbp,0),u3::SU3(0,rb),omega0p.SU3(), u3::SU3(2,0), 
																	omega0.SU3(),1,1,u3::SU3(0,rb-2),1,1
																	)
																*sqrt(
																	1.*u3::dim(omega0p)*Factorial(rb)
																	/(Factorial(2)*Factorial(rb-2)*u3::dim(omega0))
																	)
																);

															double coef2=
															(-1
																*u3::U(
																	omega0.SU3(),u3::SU3(2,0),omegap.SU3(), omega1.SU3(),
																	omega0p.SU3(),1,rho0p,omega.SU3(),1,rho0
																	)
																*u3::U(
																	u3::SU3(2,0),u3::SU3(rbp,0),omega0p.SU3(), u3::SU3(0,rb), 
																	u3::SU3(0,rb-2),1,1,omega0.SU3(),1,1
																	)
																*sqrt(
																	1.*Factorial(rbp+2)*u3::dim(omega0p)*u3::dim(u3::SU3(rbp,0))
																	/(
																		Factorial(2)*Factorial(rbp)*u3::dim(omega0)*u3::dim(u3::SU3(rbp+2,0))
																		)
																	)
																);

															double coef3=u3::U(
																omega0.SU3(),u3::SU3(2,0),omegap.SU3(), omega1.SU3(),
																omega0p.SU3(),1,rho0p,omega.SU3(),1,rho0
																);

								      									// Getting unit tensors for terms 1 and 2
															u3::UnitTensorRME unit1_labels(omegap,omega1,u3::UnitTensor(omega0p,S0,T0,rbp,Sbp,Tbp,rb-2,Sb,Tb),rho0p);
															u3::UnitTensorRME unit2_labels(omegap,omega1,u3::UnitTensor(omega0p,S0,T0,rbp+2,Sbp,Tbp,rb,Sb,Tb),rho0p);

								      									////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								      									// Term 3, sum over omega'', v'' and rho0''
								      									////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
															Eigen::MatrixXd UM(dimp,dim1);
								      									// Summing over omega''
															for (int wpp=0; wpp<omegapp_set.size(); wpp++)
																{
																	u3::U3 omegapp(omegapp_set[wpp].irrep);
									      											// omega'' subspace (v'')
																	sp3r::U3Subspace u3_subspacepp=irrep.LookUpSubspace(omegapp);
																	int dimpp=u3_subspacepp.size();
									      											// initial matrix of U()*a^\dagger
																	Eigen::MatrixXd UA(dimp,dimpp);
									      											// Obtaining K matrix for omega''
																	Eigen::MatrixXd Kpp_inv=K_matrix_map[lgip.sigma][omegapp].inverse();
																	int rho0pp_max=u3::OuterMultiplicity(omegapp.SU3(),omega0.SU3(),omegap.SU3());

																	for (int rho0pp=1; rho0pp<=rho0pp_max; rho0pp++)
																		{
																			u3::UnitTensorRME unit3_labels(
																				omegapp, 
																				omega1, 
																				u3::UnitTensor(omega0,S0,T0,rbp,Sbp,Tbp,rb,Sb,Tb),rho0pp
																				);

										      						Eigen::MatrixXd unit3_matrix=unit_tensor_rme_map[lgi_pair][unit3_labels];
																			for(int vpp=0; vpp<dimpp; vpp++)
																				{
																					u3::U3 npp;
																					int rhopp;
											      															//MultiplicityTagged<u3::U3> npp_rhopp=u3_subspacepp.GetStateLabels(vpp)
																					std::tie (npp,rhopp)=u3_subspacepp.GetStateLabels(vpp).Key();

																				}
																		}

																}

								      									//unit_tensor_matrix+=coef1*unit_tensor_rme_map[unit1_labels]

														}
												}
										}
									}
							}
					}
			}
		
	} // End function


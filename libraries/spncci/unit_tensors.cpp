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
//#include <sstream>
#include <vector>

#include "am/am.h"
#include "am/halfint.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/vcs.h"
#include "sp3rlib/sp3r.h"
#include "spncci/sp_basis.h"
#include "spncci/unit_tensors.h"

extern spncci::LGIVectorType lgi_vector;

extern std::map< u3::U3,std::map<u3::U3,Eigen::MatrixXd> > K_matrix_map;
	
  //spncci::GenerateLGIVector(lgi_vector,filename,Nsigma_0);
namespace u3
{

	std::string UnitTensor::Str() const
	{
		std::ostringstream ss;
		ss<<omega0.Str()<<"x"<<S0<<"x"<<T0<<"("<<rp<<" "<<Sp<<" "<<Tp<<", "<<r<<" "<<S<<" "<<T<<")";
		return ss.str();
	}

	std::string UnitTensorRME::Str() const
	{
		std::ostringstream ss;
		ss<<omegap.Str()<<" "<<omega.Str()<<" "<<tensor.Str()<<" "<<rho0;
		return ss.str();
	}



  void UnitSymmetryGenerator(int Nmax, std::map< int,std::vector<u3::UnitTensor> >& unit_sym_map)
// Generates a map containing (key, value) pair (N0, operator_labels) of the unit tensors 
// Generated for rp>=r.  To get the other half, use conjugation 
  {		
  	
  	for(int N0=0; N0<=Nmax; N0+=2)
	  	{
	  		std::vector<u3::UnitTensor> sym_vec;
	  		for(int Sp=0; Sp<=1; Sp++)
	  			for(int Tp=0; Tp<=1; Tp++)
	  			  for(int S=0; S<=1; S++)
	  					for (int T=0; T<=1; T++)
		  					for (int S0=abs(S-Sp); S0<=(S+Sp); S0++)
		  						for (int T0=abs(T-Tp); T0<=(T+Tp); T0++)
		  							for(int rp=0; rp<=N0+Nmax; rp++)
			  							{
			  								if ( (rp+Sp+Tp)%2!=1 )
													continue;
			  								
			  								int r=rp-N0;
			  								if ( (r+S+T)%2!=1)
			  									continue;

			  								MultiplicityTagged<u3::U3>::vector omega0_set
			  									=u3::KroneckerProduct(u3::U3(rp,0,0),u3::U3(0,0,-r));
			  								for(int w=0; w<omega0_set.size(); w++)
				  								{
				  									u3::U3 omega0(omega0_set[w].irrep);
				  									sym_vec.push_back(u3::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T));
				  									//std::cout<<"unit tensors  "<<u3::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T).Str()<<std::endl;
				  								}
			  							}	  		
  			unit_sym_map[N0]=sym_vec;
  		}
	} //end function




Eigen::MatrixXd UnitTensorMatrix(
	// LGI pair sector 
	const std::pair<int,int> lgi_pair,
	// sigma' irrep
	const sp3r::Sp3RSpace& irrepp,
	// sigma irrep
	const sp3r::Sp3RSpace& irrep,
	// unit tensor labels 
	u3::UnitTensorRME unit_labels,
   // Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
   std::map< u3::UnitTensorRME,Eigen::MatrixXd >& unit_tensor_rme_map
	)
{

	// initial declarations 		
	u3::U3 omega0;
	HalfInt S0, T0, Sbp, Tbp, Sb, Tb ;
	int rbp,rb;
	// v',v
	int N1b=2;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//  Calculate unit tensor matrix
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	spncci::LGI  lgip=lgi_vector[lgi_pair.first];
	spncci::LGI  lgi=lgi_vector[lgi_pair.second];
	u3::U3 omegap=unit_labels.omegap;
	u3::U3 omega=unit_labels.omega;
	sp3r::U3Subspace u3_subspacep=irrepp.LookUpSubspace(omegap);
	sp3r::U3Subspace u3_subspace=irrep.LookUpSubspace(omega);
	int rho0=unit_labels.rho0;
	int Nn=int(omega.N()-lgi.sigma.N());
	int Nnp=int(omegap.N()-lgip.sigma.N());


	int dimp=u3_subspacep.size();
	int dim=u3_subspace.size();
	assert(dimp!=0 && dim!=0);
	Eigen::MatrixXd unit_tensor_matrix=Eigen::MatrixXd::Zero(dimp,dim);
	u3::UnitTensor unit_tensor=unit_labels.tensor;
	std::tie (omega0, S0, T0, rbp, Sbp, Tbp, rb, Sb, Tb) = unit_tensor.Key();


	Eigen::MatrixXd Kp=K_matrix_map[lgip.sigma][omegap];
	Eigen::MatrixXd K_inv=K_matrix_map[lgi.sigma][omega].inverse();

	MultiplicityTagged<u3::U3>::vector omegapp_set=KroneckerProduct(omegap, u3::U3(0,0,-2)); 
	MultiplicityTagged<u3::U3>::vector omega0p_set=KroneckerProduct(omega0, u3::U3(2,0,0));
	MultiplicityTagged<u3::U3>::vector omega1_set=KroneckerProduct(omega, u3::U3(0,0,-2));

	// summing over omega1
	for (int w1=0; w1<omega1_set.size(); w1++)
		{	
			u3::U3 omega1=omega1_set[w1].irrep;
			
			//check that omega1 in irrep  
			if (not irrep.ContainsSubspace(omega1))
      		continue;

      // omega1 sector
			sp3r::U3Subspace u3_subspace1=irrep.LookUpSubspace(omega1);
			int dim1=u3_subspace1.size();

		  // Look up K1 matrix (dim v1, v1)
			Eigen::MatrixXd K1=K_matrix_map[lgi.sigma][omega1];
			
			// Initializing unit tensor matrix with dim. v' v1
			Eigen::MatrixXd unit_matrix= Eigen::MatrixXd::Zero(dimp,dim1);

			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Matrix of B*U coefs with dim v1 and v
			///////////////////////////////////////////////////////////////////////////////////////////////////////////
			//std::cout<<"hi"<<std::endl;
			Eigen::MatrixXd BU(dim1,dim);
			//iterating over (n,rho)
			for (int m=0; m<dim; m++)
				{
					MultiplicityTagged<u3::U3> n_rho=u3_subspace.GetStateLabels(m);
					u3::U3 n(n_rho.irrep);

					for (int m1=0; m1<dim1; m1++)
						{
							MultiplicityTagged<u3::U3> n1_rho1=u3_subspace1.GetStateLabels(m1);
							u3::U3 n1(n1_rho1.irrep);

							if (
								u3::OuterMultiplicity(n1.SU3(), u3::SU3(2,0),n.SU3())>0
								)
								BU(m1,m)=(
									vcs::BosonCreationRME(n,n1)
			 						*u3::U(u3::SU3(2,0),n1.SU3(),omega.SU3(),lgi.sigma.SU3(),n.SU3(),1,n_rho.tag,omega1.SU3(),n1_rho1.tag,1)
									);
							
							else
								{
									BU(m1,m)=0;
									//continue???
								}
						}	
				}								        
			Eigen::MatrixXd KBUK(dim1,dim);
			KBUK=K1*BU*K_inv;
			////////////////////////////////////////////////////////////////////////////////////////////////////////

		  //summing over omega0'
			for (int w0p=0; w0p<omega0p_set.size(); w0p++)
				{
					u3::U3 omega0p=omega0p_set[w0p].irrep;
					int rho0p_max=OuterMultiplicity(omega1.SU3(),omega0p.SU3(),omegap.SU3());
				  
				  // summing over rho0'
					for (int rho0p=1; rho0p<=rho0p_max; rho0p++)
						{
							
							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							// third term
							// sum over omega'', v'' and rho0''
							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							
							double coef3=u3::U(
										omega0.SU3(),u3::SU3(2,0),omegap.SU3(), omega1.SU3(),
										omega0p.SU3(),1,rho0p,omega.SU3(),1,rho0
										);
							//Initilize 3rd-term-unit-tensor matrix
							Eigen::MatrixXd unit3_matrix=Eigen::MatrixXd::Zero(dimp,dim1);
							if ( rbp<=(Nnp-2+N1b) && rb<=(Nn-2+N1b) )
								{	
								  // Summing over omega''
									for (int wpp=0; wpp<omegapp_set.size(); wpp++)
										{
											
											u3::U3 omegapp(omegapp_set[wpp].irrep);
										  if (not irrepp.ContainsSubspace(omegapp))
			      						continue;
			      					if (
			      						unit_tensor_rme_map.count(
			      						u3::UnitTensorRME(omegapp,omega1,u3::UnitTensor(omega0,S0,T0,rbp,Sbp,Tbp,rb,Sb,Tb),1))==0
			      						)
			      						continue;						

										  // omega'' subspace (v'')
											sp3r::U3Subspace u3_subspacepp=irrepp.LookUpSubspace(omegapp);
											int dimpp=u3_subspacepp.size();
										  // Obtaining K matrix for omega''
											Eigen::MatrixXd Kpp_inv=K_matrix_map[lgip.sigma][omegapp].inverse();
											// Initialize matrix of a^\dagger for A matrix
											Eigen::MatrixXd boson_matrix(dimp,dimpp);
									    //Constructing a^\dagger matrix

											for(int vpp=0; vpp<dimpp; vpp++)
												{
													MultiplicityTagged<u3::U3> npp_rhopp=u3_subspacepp.GetStateLabels(vpp);

													for(int vp=0; vp<dimp; vp++)
														{
															MultiplicityTagged<u3::U3> np_rhop=u3_subspacep.GetStateLabels(vp);

															if (u3::OuterMultiplicity(npp_rhopp.irrep.SU3(), u3::SU3(2,0),np_rhop.irrep.SU3())>0)
																boson_matrix(vp,vpp)=
																	vcs::U3BosonCreationRME(lgip.sigma, np_rhop, omegap, lgip.sigma, npp_rhopp,omegapp);
															else
																boson_matrix(vp,vpp)=0;

														}
												}
											Eigen::MatrixXd unit3pp_matrix=Eigen::MatrixXd::Zero(dimpp,dim1);
											int rho0pp_max=u3::OuterMultiplicity(omega1.SU3(),omega0.SU3(),omegapp.SU3());
											// Summing over rho0''

											for (int rho0pp=1; rho0pp<=rho0pp_max; rho0pp++)
												{
													// Retriving unit tensor matrix 
													u3::UnitTensorRME unit3_labels(
														omegapp, 
														omega1, 
														u3::UnitTensor(omega0,S0,T0,rbp,Sbp,Tbp,rb,Sb,Tb),rho0pp
														);

													assert(unit_tensor_rme_map.count(unit3_labels)>0);
													
													unit3pp_matrix+=
													u3::U(u3::SU3(2,0),omega0.SU3(),omegap.SU3(), omega1.SU3(),
														omega0p.SU3(),1,rho0p,omegapp.SU3(),rho0pp, 1)
														*unit_tensor_rme_map[unit3_labels];

												} //end rho0pp
											// matrix product (v',v')*(v',v'')*(v'',v1)
											unit3_matrix+=Kp*boson_matrix*Kpp_inv*unit3pp_matrix;
										} // end wpp
									}
								unit_matrix+=coef3*unit3_matrix;

							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							//first term 
							//////////////////////////////////////////////////////////////////////////////////////////////////////////	
						  if (u3::OuterMultiplicity(u3::SU3(rbp,0),u3::SU3(0,rb-2),omega0p.SU3())>0 && (rb-2)>=0)
								{

									
									u3::UnitTensorRME unit1_labels(omegap,omega1,u3::UnitTensor(omega0p,S0,T0,rbp,Sbp,Tbp,rb-2,Sb,Tb),rho0p);
									if (unit_tensor_rme_map.count(unit1_labels)!=0)
										{

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

											unit_matrix+=coef1*unit_tensor_rme_map[unit1_labels];
										}
								}
							
							//////////////////////////////////////////////////////////////////////////////////////////////////////////
							// second term 
							//////////////////////////////////////////////////////////////////////////////////////////////////////////	
						  if (
						  			(u3::OuterMultiplicity(u3::SU3(rbp+2,0),u3::SU3(0,rb),omega0p.SU3())>0)
								  	&&
								  	rb<=(omega1.N()-lgi.sigma.N()+N1b)
								  	&&
								  	(rbp+2)<=(omegap.N()-lgip.sigma.N()+N1b)
								  )
								{
									u3::UnitTensorRME unit2_labels(omegap,omega1,u3::UnitTensor(omega0p,S0,T0,rbp+2,Sbp,Tbp,rb,Sb,Tb),rho0p);
									double coef2;
									if(unit_tensor_rme_map.count(unit2_labels)>0)
										{
											

											coef2=
												(-1
													*u3::U(
														omega0.SU3(),u3::SU3(2,0),omegap.SU3(), omega1.SU3(),
														omega0p.SU3(),1,rho0p,omega.SU3(),1,rho0
														)
													*u3::U(
														u3::SU3(2,0),u3::SU3(rbp,0),omega0p.SU3(), u3::SU3(0,rb), 
														u3::SU3(rbp+2,0),1,1,omega0.SU3(),1,1
														)
													*sqrt(
														1.*Factorial(rbp+2)*u3::dim(omega0p)*u3::dim(u3::SU3(rbp,0))
														/(
															Factorial(2)*Factorial(rbp)*u3::dim(omega0)*u3::dim(u3::SU3(rbp+2,0))
															)
														)
													);
											unit_matrix+=coef2*unit_tensor_rme_map[unit2_labels];
										}

								}
							//////////////////////////////////////////////////////////////////////////////////////////////////////////
						} //end rho0p
				} //end sum over w0p
				// summing over n, rho, n1, rho1, v1
				unit_tensor_matrix+=unit_matrix*KBUK;				
		}// end sum over omega1
		assert(unit_tensor_matrix.cols()!=0 && unit_tensor_matrix.rows()!=0);
		return unit_tensor_matrix;
	} // End function


	void UnitTensorMatrixGenerator(
	// single particle cutoff, relative particle <= 2*N1b+Nn
		int N1b,
	// boson number cutoff
		int Nmax, 
	// a given spncci sector pair given as index pair  from global list lgi_vector 
		std::pair<int,int> lgi_pair,
	// Address to map with list of unit tensor labels with key N0 
		std::map< int,std::vector<u3::UnitTensor>>& unit_sym_map,
	// Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
     std::map< u3::UnitTensorRME,Eigen::MatrixXd >& unit_tensor_rme_map

		)
	// Generates all unit tensor matrix matrices between states in the irreps of lgi_pair
	// The unit tensors are stored in the map unit_tensor_rme_map which has key lgi_pair and 
	// value map(matrix labels for w'w sector, matrix) 
	{
		// initial declarations 		
		u3::U3 omega0;
		HalfInt S0, T0, Sbp, Tbp, Sb, Tb ;
		int rbp,rb;

	  // Extracting LGI labels from pair
		spncci::LGI  lgip=lgi_vector[lgi_pair.first];
		spncci::LGI  lgi=lgi_vector[lgi_pair.second];

	  //Generate Sp(3,R) irreps 
		sp3r::Sp3RSpace irrepp(lgip.sigma,Nmax-lgip.Nex);
		sp3r::Sp3RSpace irrep(lgi.sigma, Nmax-lgi.Nex);

	  // Calculating K matrices for each sigma in LGI set and storing in map K_matrix_map with key sigma
		for (int k = 0; k<irrep.size(); k++)
			{
				std::map<u3::U3,Eigen::MatrixXd> K_map;
				vcs::GenerateKMatrices(irrep, K_map);
				K_matrix_map[lgi.sigma]=K_map;
			}
		for (int kp = 0; kp<irrepp.size(); kp++)
			{
				std::map<u3::U3,Eigen::MatrixXd> K_map;
				vcs::GenerateKMatrices(irrepp, K_map);
				K_matrix_map[lgip.sigma]=K_map;
			}      	      		
		////////////////////////////////////////////////////////////////////////////////////
		// Looping over omega' and omega subspaces 
		////////////////////////////////////////////////////////////////////////////////////
		//  omega' subspace
		for(int ip=0; ip<irrepp.size(); ip++ )
			{
				//sp3r::U3Subspace& u3_subspace=irrep.GetSubspace(i);
				u3::U3 omegap=irrepp.GetSubspace(ip).GetSubspaceLabels();		 
				int Nnp=int(omegap.N()-lgip.sigma.N());
				
				//omega subspace
				for(int i=0; i<irrep.size(); i++ )
					{
						//sp3r::U3Subspace& u3_subspacep=irrepp.GetSubspace(ip);
						u3::U3 omega=irrep.GetSubspace(i).GetSubspaceLabels();
						if (TwiceValue(omegap.N())==TwiceValue(omega.N()) && (not (omegap==omega)))
							continue;

						int Nn=int(omega.N()-lgi.sigma.N());
						// Get set of operator labels for given omega'omega sector
						std::vector<u3::UnitTensor>& operator_set=unit_sym_map[abs(int(omegap.N()-omega.N()))];
						// Iterating over the operator labels 						

						for (int w=0; w<operator_set.size(); w++)
							{						      		
								u3::UnitTensor unit_tensor=operator_set[w];

								std::tie (omega0, S0, T0, rbp, Sbp, Tbp, rb, Sb, Tb) = unit_tensor.Key();
								int rho0_max;

								if ((omegap.N()-omega.N())<0)
									{
										if( (rbp > (Nn+N1b)) || (rb > (Nnp+N1b)) )
											continue;

										rho0_max=OuterMultiplicity(omegap.SU3(),omega0.SU3(),omega.SU3());
									}
								else 
									{
										if( (rb > (Nn+N1b)) || (rbp > (Nnp+N1b)) )
											continue;

										rho0_max=OuterMultiplicity(omega.SU3(),omega0.SU3(),omegap.SU3());
									}	
								// Iterating over outer multiplicity
								for (int rho0=1; rho0<=rho0_max; rho0++)
									{
										
										// if LGI, unit tensor matrix is already calculated 
										if (omegap.N()==lgip.sigma.N() && omega.N()==lgi.sigma.N())
												continue;
										u3::UnitTensorRME unit_labels;
										Eigen::MatrixXd temp_matrix;
										
										// In the special case that omegap.N()!=sigmap.N() but omega.N()==sigma.N(), then to calculate we
										// need to calculate the conjugate transpose of the unit tensor matrix and then invert and multiply 
										// by factor to obtain desired matrix
										if (omegap.N()!=lgip.sigma.N() && omega.N()==lgi.sigma.N())
											{
												u3::UnitTensorRME unit_map_key(omegap,omega,unit_tensor,rho0);
												unit_labels
													=u3::UnitTensorRME(
																							omega,omegap,
																							UnitTensor(u3::Conjugate(omega0),S0,T0,rb,Sb,Tb,rbp,Sbp,Tbp),
																							rho0
																						);

												temp_matrix=u3::UnitTensorMatrix(
														std::pair<int,int>(lgi_pair.second,lgi_pair.first), irrep, irrepp,
														unit_labels, 
														unit_tensor_rme_map 
														);
												
												unit_tensor_rme_map[unit_map_key]=
													ParitySign(rbp+rb+ConjugationGrade(omega)+ConjugationGrade(omegap))
													*sqrt(1.*dim(u3::SU3(rbp,0))*dim(omega)/(dim(u3::SU3(rb,0))*dim(omegap)))
													*temp_matrix.transpose();
											}


										else 
											{
												// In the case that omegap.N()< omega.N(), then the correct unit tensor symmetry will be the 
												//  conjugate symmetry	
												if ( (omegap.N()-omega.N())<0 )
													// Operator labels are conjugates 
													unit_labels=
														u3::UnitTensorRME(
																							omegap,omega,
																							u3::UnitTensor(Conjugate(omega0),S0, T0,rb,Sb,Tb,rbp,Sbp,Tbp),
																							rho0
																						);
												else
													unit_labels=u3::UnitTensorRME(omegap,omega,unit_tensor,rho0);
												/////////////////////
												temp_matrix=u3::UnitTensorMatrix(lgi_pair, irrepp, irrep, unit_labels,unit_tensor_rme_map);
												if (temp_matrix.any())
													unit_tensor_rme_map[unit_labels]
														=u3::UnitTensorMatrix(lgi_pair, irrepp, irrep, unit_labels,unit_tensor_rme_map);
												// else
													//std::cout<<unit_labels.Str()<<unit_tensor_rme_map[unit_labels]
													// <<unit_tensor_rme_map.count(unit_labels)<<std::endl;

											}
									}
							}
					}
			}
	}

} // End namespace 

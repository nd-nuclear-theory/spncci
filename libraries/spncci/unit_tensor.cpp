/****************************************************************
  unit_tensor.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/


#include <omp.h>
#include "spncci/unit_tensor.h"

extern spncci::LGIVectorType lgi_vector;
extern std::map< u3::U3, vcs::MatrixCache > K_matrix_map;

namespace spncci
{
  typedef std::pair<UnitTensorU3Sector, Eigen::MatrixXd> UnitTensorU3SectorPair;

  std::string UnitTensor::Str() const
  {
    std::ostringstream ss;
    ss << omega0_.Str() << "x" << S0_ << "x" << T0_
       << "(" << rp_ << " " << Sp_ << " " << Tp_
       << ", " << r_ << " " << S_ << " " << T_ << ")";
    return ss.str();
  }

  std::string UnitTensorU3Sector::Str() const
  {
    std::ostringstream ss;
    ss << omegap_.Str() << " " << omega_.Str() << " " << tensor_.Str() << " " << rho0_;
    return ss.str();
  }

  std::pair<int,int> GetNSectorIndices(const int Nmax, const int irrep_size, const int Nn, std::vector<int>& NPartition)
  {
    int i_min, i_max;
    // Starting index
    if (Nn==0)
      i_min=0;
    else
      i_min=NPartition[Nn/2];
    // Ending index+1.   
    if (NPartition.size()==(Nn+2)/2)
      i_max=irrep_size-1;
    else 
      i_max=NPartition[(Nn+2)/2]-1;
    return std::pair<int,int>(i_min,i_max);
  }


  void GenerateUnitTensors(int Nmax, std::map< int,std::vector<spncci::UnitTensor> >& unit_sym_map)
  // Generates a map containing (key, value) pair (N0, operator_labels) of the unit tensors 
  // Generated for rp>=r.  To get the other half, use conjugation 
  {		
  	#ifdef VERBOSE
    std::cout<<"Entering GenerateUnitTensors"<<std::endl;
    #endif
    for(int N0=0; N0<=Nmax; N0+=2)
      {
        std::vector<spncci::UnitTensor> sym_vec;
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

                        MultiplicityTagged<u3::SU3>::vector omega0_set
                        //  =u3::KroneckerProduct(u3::U3(rp,0,0),u3::U3(0,0,-r));
                          =u3::KroneckerProduct(u3::SU3(rp,0),u3::SU3(0,r));
                        for(int w=0; w<omega0_set.size(); w++)
                          {
                            //
                            u3::U3 omega0(N0,omega0_set[w].irrep);
                            //
                            // u3::U3 omega0(omega0_set[w].irrep);
                            sym_vec.push_back(spncci::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T));
                            //std::cout<<"unit tensors  "<<spncci::UnitTensor(omega0,S0,T0,rp,Sp,Tp,r,S,T).Str()<<std::endl;
                          }
                      }	  		
        unit_sym_map[N0]=sym_vec;
      }
  #ifdef VERBOSE
  std::cout<<"Exiting GenerateUnitTensors"<<std::endl;
  #endif
  } //end function

  ////////////////////////////////////////////////////////////////////////////////////
  void GenerateUnitTensorU3SectorLabels(
                                        // single particle cutoff, relative particle <= 2*N1b+Nn
                                        int N1b,
                                        // boson number cutoff
                                        int Nmax,
                                        // a given spncci sector pair given as index pair from global list lgi_vector 
                                        std::pair<int,int>  lgi_pair,
                                        // Address to map with list of unit tensor labels with key N0 
                                        std::map< int,std::vector<spncci::UnitTensor>>& unit_sym_map,
                                        // For each NpN pair key in map the corresponding value is a vector of UnitTensorU3Sectors. 
                                        std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>>& unit_tensor_NpN_sector_map
                                        )
  // Generates labels of all sectors of unit tensor matrix matrices between states in the irreps of lgi_pair
  //
  // TODO replace argument lgi_pair with (lgip,lgi)
  {   
    #ifdef VERBOSE
    std::cout<<"Entering GenerateU3SectorLabels"<<std::endl;
    #endif

    // initial declarations     
    u3::U3 omega0;
    HalfInt S0, T0, Sbp, Tbp, Sb, Tb ;
    int rbp,rb;
    spncci::UnitTensorU3Sector unit_U3Sectors;

    // Extracting LGI labels from pair
    const spncci::LGI& lgip=lgi_vector[lgi_pair.first];
    const spncci::LGI& lgi=lgi_vector[lgi_pair.second];
    u3::U3 sigmap=lgip.sigma;
    u3::U3 sigma=lgi.sigma; 

    const sp3r::Sp3RSpace& irrepp=lgip.Sp3RSpace();
    const sp3r::Sp3RSpace& irrep=lgi.Sp3RSpace();

    int irrep_size=irrep.size();
    int irrepp_size=irrepp.size();

    // partition irreps by Nn and Nnp.  Each int in vector corresponds to the start of the next N space 
    std::vector<int> NpPartition=sp3r::PartitionIrrepByNn(irrepp, Nmax);
    std::vector<int> NPartition=sp3r::PartitionIrrepByNn(irrep, Nmax);
    ////////////////////////////////////////////////////////////////////////////////////
    // Looping over omega' and omega subspaces 
    ////////////////////////////////////////////////////////////////////////////////////

    for (int Nsum=2; Nsum<=2*Nmax; Nsum+=2)
      for (int Nnp=0; Nnp<=std::min(Nsum,Nmax); Nnp+=2)
        {
          if((Nnp+lgip.Nex)>Nmax)
            continue;
          int Nn=Nsum-Nnp;

          if ((Nn+lgi.Nex)>Nmax)
            continue; 

          std::pair<int,int> NpN_pair(Nnp,Nn);
          // Selecting section of spaces to iterate over
          int ip_min, ip_max, i_min, i_max;
          std::tie (i_min,i_max)=GetNSectorIndices(Nmax, irrep_size, Nn, NPartition);
          std::tie (ip_min,ip_max)=GetNSectorIndices(Nmax, irrepp_size, Nnp, NpPartition);
          // Get set of operator labels for given omega'omega sector
          std::vector<spncci::UnitTensor>& operator_set=unit_sym_map[abs(lgip.Nex-lgi.Nex+Nnp-Nn)];
          //  omega' subspace
          for(int ip=ip_min; ip<=ip_max; ip++ )
            {
              u3::U3 omegap=irrepp.GetSubspace(ip).GetSubspaceLabels();    
              //omega subspace
              for(int i=i_min; i<=i_max; i++ )
                {
                  u3::U3 omega=irrep.GetSubspace(i).GetSubspaceLabels();
                  // Iterating over the operator labels             
                  for (int w=0; w<operator_set.size(); w++)
                    {                     
                      spncci::UnitTensor unit_tensor=operator_set[w];
                      //unpack unit_tensor labels 
                      std::tie (omega0, S0, T0, rbp, Sbp, Tbp, rb, Sb, Tb) = unit_tensor.Key();
                      
                      //Checking angular momentum constraint 
                      if (abs(lgi.S+lgip.S)<S0)
                        continue;
                      ///////////////////////////////////////////////////////////////////////////////////////
                      // defining rho0_max for the two cases 
                      ///////////////////////////////////////////////////////////////////////////////////////
                      // if Np<N, then the operator labels need to be conjugated 
                      bool N_greater=(omegap.N()-omega.N())<0;

                      int r1=N_greater?rbp:rb;
                      int r2=N_greater?rb:rbp;
                      if (r1>(Nn+N1b+lgi.Nex)|| r2>(Nnp+N1b+lgip.Nex))
                        continue;

                      u3::SU3 lm=N_greater?omegap.SU3():omega.SU3();
                      u3::SU3 lmp=N_greater?omega.SU3():omegap.SU3();
                      int rho0_max=OuterMultiplicity(lm,omega0.SU3(),lmp);
                      ///////////////////////////////////////////////////////////////////////////////////////
                      // Iterating over outer multiplicity
                      for (int rho0=1; rho0<=rho0_max; rho0++)
                        {                   
                          unit_U3Sectors=N_greater?
                            spncci::UnitTensorU3Sector(omegap,omega,spncci::UnitTensor(Conjugate(omega0),S0, T0,rb,Sb,Tb,rbp,Sbp,Tbp),rho0)
                            :spncci::UnitTensorU3Sector(omegap,omega,unit_tensor,rho0);

                          unit_tensor_NpN_sector_map[NpN_pair].push_back(unit_U3Sectors);
                        }
                    }
                }
            }
        }    

    #ifdef VERBOSE
    std::cout<<"Exiting GenerateU3SectorLabels"<<std::endl;
    #endif
  }
  ////////////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd UnitTensorMatrix(
                                   u3::UCoefCache& u_coef_cache,
                                   // LGI pair sector 
                                  const spncci::LGI& lgip,
                                  const spncci::LGI& lgi,
                                   // vector of addresses to relevant Np,N sectors of unit tensor matrix
                                   spncci::UnitTensorSectorsCache& sector_NpN2,
                                   spncci::UnitTensorSectorsCache& sector_NpN4,
                                   // sigma' irrep
                                   const sp3r::Sp3RSpace& irrepp,
                                   // sigma irrep
                                   const sp3r::Sp3RSpace& irrep,
                                   // unit tensor labels 
                                   spncci::UnitTensorU3Sector unit_labels
                                   )
  {
    #ifdef VERBOSE
    std::cout<<"Entering UnitTensorMatrix"<<std::endl;
    #endif
    // initial declarations 		
    u3::U3 omegap, omega, omega0;
    HalfInt S0, T0, Sbp, Tbp, Sb, Tb ;
    UnitTensor tensor;
    int rbp, rb, rho0;
    // v',v
    int N1b=2;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Set up for calculation 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Extracting labels 

    std::tie(omegap,omega,tensor,rho0) = unit_labels.Key();
    sp3r::U3Subspace u3_subspacep = irrepp.LookUpSubspace(omegap);
    sp3r::U3Subspace u3_subspace  = irrep.LookUpSubspace(omega);

    int Nn=int(omega.N()-lgi.sigma.N());
    int Nnp=int(omegap.N()-lgip.sigma.N());

    int dimp=u3_subspacep.size();
    int dim=u3_subspace.size();
    assert(dimp!=0 && dim!=0);
	
    // unpacking the unit tensor labels 
    std::tie(omega0, S0, T0, rbp, Sbp, Tbp, rb, Sb, Tb) = tensor.Key();

    // Extracting K matrices for lgi and lgip from the K_matrix_maps 
    vcs::MatrixCache& K_matrix_map_lgi=K_matrix_map[lgi.sigma];
    vcs::MatrixCache& K_matrix_map_lgip=K_matrix_map[lgip.sigma];

    Eigen::MatrixXd Kp=K_matrix_map_lgip[omegap];
    Eigen::MatrixXd K_inv=K_matrix_map_lgi[omega].inverse();

    // Precalculating kronecker products used in sum to calculate unit tensor matrix
    MultiplicityTagged<u3::U3>::vector omegapp_set=KroneckerProduct(omegap, u3::U3(0,0,-2)); 
    MultiplicityTagged<u3::U3>::vector omega0p_set=KroneckerProduct(omega0, u3::U3(2,0,0));
    MultiplicityTagged<u3::U3>::vector omega1_set=KroneckerProduct(omega, u3::U3(0,0,-2));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Calculate unit tensor matrix
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Eigen::MatrixXd unit_tensor_matrix=Eigen::MatrixXd::Zero(dimp,dim);

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
        Eigen::MatrixXd K1=K_matrix_map_lgi[omega1];
			
        // Initializing unit tensor matrix with dim. v' v1
        Eigen::MatrixXd unit_matrix= Eigen::MatrixXd::Zero(dimp,dim1);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Matrix of B*U coefs with dim v1 and v
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
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

                if (u3::OuterMultiplicity(n1.SU3(), u3::SU3(2,0),n.SU3())>0)
                {
                  int r12_max, r12_3_max, r23_max, r1_23_max;
                  std::tie(r12_max,r12_3_max,r23_max,r1_23_max) = UMultiplicity(u3::SU3(2,0),n1.SU3(),omega.SU3(),lgi.sigma.SU3(),n.SU3(),omega1.SU3());
                  BU(m1,m)=(
                            vcs::BosonCreationRME(n,n1)
                            *u3::UCached(u_coef_cache,u3::SU3(2,0),n1.SU3(),omega.SU3(),lgi.sigma.SU3(),n.SU3(),1,n_rho.tag,omega1.SU3(),n1_rho1.tag,1)
                            );
                }
                else
                  {
                    BU(m1,m)=0;
                  }
              }	
          }								        
        Eigen::MatrixXd KBUK(dim1,dim);
        KBUK.noalias()=K1*BU*K_inv;
        // std::cout<<"KBUK "<<KBUK<<std::endl;
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        //summing over omega0'
        for (int w0p=0; w0p<omega0p_set.size(); w0p++)
          {
            u3::U3 omega0p=omega0p_set[w0p].irrep;
            int rho0p_max=OuterMultiplicity(omega1.SU3(),omega0p.SU3(),omegap.SU3());
				  
            // summing over rho0'
            for (int rho0p=1; rho0p<=rho0p_max; rho0p++)
              {
                // std::pair<int,int>NpN4_pair(Nnp-2,Nn-2);
                //////////////////////////////////////////////////////////////////////////////////////////////////////////
                // third term
                // sum over omega'', v'' and rho0''
                //////////////////////////////////////////////////////////////////////////////////////////////////////////
                double coef3=u3::UCached(u_coef_cache,
                                   omega0.SU3(),u3::SU3(2,0),omegap.SU3(), omega1.SU3(),
                                   omega0p.SU3(),1,rho0p,omega.SU3(),1,rho0
                                   );
                //Initilize 3rd-term-unit-tensor matrix
                Eigen::MatrixXd unit3_matrix=Eigen::MatrixXd::Zero(dimp,dim1);

                if ( rbp<=(Nnp-2+N1b+lgip.Nex) && rb<=(Nn-2+N1b+lgi.Nex) )
                  {	
                   
                    // Summing over omega''
                    for (int wpp=0; wpp<omegapp_set.size(); wpp++)
                      {
											 

                        u3::U3 omegapp(omegapp_set[wpp].irrep);
                        if (not irrepp.ContainsSubspace(omegapp))
                          continue;
                        spncci::UnitTensorU3Sector unit_tensor_u3_sector_1(omegapp,omega1,spncci::UnitTensor(omega0,S0,T0,rbp,Sbp,Tbp,rb,Sb,Tb),1);
                        if ((sector_NpN4).count(unit_tensor_u3_sector_1)==0)
                          continue;	
                        // omega'' subspace (v'')
                        sp3r::U3Subspace u3_subspacepp=irrepp.LookUpSubspace(omegapp);
                        int dimpp=u3_subspacepp.size();
                        // Obtaining K matrix for omega''
                        Eigen::MatrixXd Kpp_inv=K_matrix_map_lgip[omegapp].inverse();
                        // Initialize matrix of a^\dagger for A matrix
                        Eigen::MatrixXd boson_matrix(dimp,dimpp);
                        //Constructing a^\dagger matrix
                        for(int vpp=0; vpp<dimpp; vpp++)
                          {
                            MultiplicityTagged<u3::U3> npp_rhopp=u3_subspacepp.GetStateLabels(vpp);
                            const u3::U3& npp=npp_rhopp.irrep;
                            const int& rhopp=npp_rhopp.tag;

                            for(int vp=0; vp<dimp; vp++)
                              {
                                MultiplicityTagged<u3::U3> np_rhop=u3_subspacep.GetStateLabels(vp);
                                const u3::U3& np=np_rhop.irrep;
                                const int& rhop=np_rhop.tag; 
                                if (u3::OuterMultiplicity(npp.SU3(), u3::SU3(2,0),np.SU3())>0)
                                  boson_matrix(vp,vpp)=
                                    //vcs::U3BosonCreationRME(lgip.sigma, np_rhop, omegap, lgip.sigma, npp_rhopp,omegapp);
                                    vcs::BosonCreationRME(np,npp)
                                    *u3::UCached(u_coef_cache,lgip.sigma.SU3(), npp.SU3(), omegap.SU3(), u3::SU3(2,0), 
                                           omegapp.SU3(), rhopp, 1, np.SU3(), 1, rhop);

                                else
                                  boson_matrix(vp,vpp)=0;

                              }
                          }

                        //std::cout<<"boson matrix "<<boson_matrix<<std::endl;
                        Eigen::MatrixXd unit3pp_matrix=Eigen::MatrixXd::Zero(dimpp,dim1);
                        int rho0pp_max=u3::OuterMultiplicity(omega1.SU3(),omega0.SU3(),omegapp.SU3());
                        // Summing over rho0''

                        for (int rho0pp=1; rho0pp<=rho0pp_max; rho0pp++)
                          {
                            // Retriving unit tensor matrix 
                            spncci::UnitTensorU3Sector unit3_labels(omegapp,omega1,spncci::UnitTensor(omega0,S0,T0,rbp,Sbp,Tbp,rb,Sb,Tb),rho0pp);

                            assert(sector_NpN4.count(unit3_labels)>0);

                            unit3pp_matrix+=
                              u3::UCached(u_coef_cache,u3::SU3(2,0),omega0.SU3(),omegap.SU3(), omega1.SU3(),
                                    omega0p.SU3(),1,rho0p,omegapp.SU3(),rho0pp, 1)
                              *sector_NpN4[unit3_labels];

                          } //end rho0pp
                        // matrix product (v',v')*(v',v'')*(v'',v1)
                        unit3_matrix+=Kp*boson_matrix*Kpp_inv*unit3pp_matrix;
                      } // end wpp
                  }
                // std::cout<<"unit3_matrix "<<unit3_matrix<<std::endl;
                unit_matrix+=coef3*unit3_matrix;

                // std::cout<<"term 3  "<<unit_matrix<<std::endl;							
                //std::pair<int,int>NpN2_pair(Nnp,Nn-2);
                //////////////////////////////////////////////////////////////////////////////////////////////////////////
                //first term 
                //////////////////////////////////////////////////////////////////////////////////////////////////////////	
                if (u3::OuterMultiplicity(u3::SU3(rbp,0),u3::SU3(0,rb-2),omega0p.SU3())>0 && (rb-2)>=0)
                  {
                    spncci::UnitTensorU3Sector unit1_labels(omegap,omega1,spncci::UnitTensor(omega0p,S0,T0,rbp,Sbp,Tbp,rb-2,Sb,Tb),rho0p);
                    if ( sector_NpN2.count(unit1_labels)!=0)
                      {
                        double coef1=(
                                      u3::UCached(u_coef_cache,omega0.SU3(),u3::SU3(2,0),omegap.SU3(), omega1.SU3(),
                                            omega0p.SU3(),1,rho0p,omega.SU3(),1,rho0)
                                      *u3::UCached(u_coef_cache,u3::SU3(rbp,0),u3::SU3(0,rb),omega0p.SU3(), u3::SU3(2,0), 
                                             omega0.SU3(),1,1,u3::SU3(0,rb-2),1,1)
                                      *sqrt(1.*u3::dim(omega0p)*Choose(rb,2)/u3::dim(omega0))
                                      );											
                        unit_matrix+=coef1*sector_NpN2[unit1_labels];
                      }
                  }
                // std::cout<< "term1  "<<unit_matrix<<std::endl;
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
                    spncci::UnitTensorU3Sector unit2_labels(omegap,omega1,spncci::UnitTensor(omega0p,S0,T0,rbp+2,Sbp,Tbp,rb,Sb,Tb),rho0p);
                    double coef2;
                    if(sector_NpN2.count(unit2_labels)>0)
                      {
                        coef2=
                          (-1
                           *u3::UCached(u_coef_cache,
                                  omega0.SU3(),u3::SU3(2,0),omegap.SU3(), omega1.SU3(),
                                  omega0p.SU3(),1,rho0p,omega.SU3(),1,rho0
                                  )
                           *u3::UCached(u_coef_cache,
                                  u3::SU3(2,0),u3::SU3(rbp,0),omega0p.SU3(), u3::SU3(0,rb), 
                                  u3::SU3(rbp+2,0),1,1,omega0.SU3(),1,1
                                  )
                           *sqrt(
                                 1.*Choose(rbp+2,2)*u3::dim(omega0p)*u3::dim(u3::SU3(rbp,0))
                                 /(u3::dim(omega0)*u3::dim(u3::SU3(rbp+2,0)))
                                 )
                           );
                        // std::cout<<"unit  "<<unit_matrix<<std::endl<<"sector  "<<sector_NpN2[unit2_labels]<<std::endl;
                        unit_matrix+=coef2*sector_NpN2[unit2_labels];
                      }
                  }
                  // std::cout<<"term1"<<std::endl;
                //////////////////////////////////////////////////////////////////////////////////////////////////////////
              } //end rho0p

          } //end sum over w0p
        // summing over n, rho, n1, rho1, v1
        unit_tensor_matrix+=unit_matrix*KBUK;


      }// end sum over omega1
    assert(unit_tensor_matrix.cols()!=0 && unit_tensor_matrix.rows()!=0);
    // std::cout<<unit_tensor_matrix <<std::endl;
    #ifdef VERBOSE
    std::cout<<"Exiting UnitTensorMatrix"<<std::endl;
    #endif
    return unit_tensor_matrix;
  } // End function


  void GenerateUnitTensorU3Sector(
                                  u3::UCoefCache& u_coef_cache,
                                  const spncci::UnitTensorU3Sector& unit_tensor_u3_sector, 
                                  // LGI pair sector 
                                  const spncci::LGI& lgip,
                                  const spncci::LGI& lgi,                                 
                                  // vector of addresses to relevant Np,N sectors of unit tensor matrix
                                  // Eigen doesn't like const 
                                  spncci::UnitTensorSectorsCache& sector_NpN2,
                                  spncci::UnitTensorSectorsCache& sector_NpN4,
                                  // sigma' irrep
                                  const sp3r::Sp3RSpace& irrepp,
                                  // sigma irrep
                                  const sp3r::Sp3RSpace& irrep,
                                  bool Nn_zero,
                                  std::vector< UnitTensorU3SectorPair >& unit_tensor_u3_sector_pairs)
  {
    //calculate unit tensor matrix.   
    // int zerocout=0;
    /////////////////////////////////////////////////////////////////////////////////////
    #ifdef VERBOSE
    std::cout<<"Entering GenerateUnitTensorU3Sector"<<std::endl;
    #endif
    Eigen::MatrixXd temp_matrix;
    // In the special case that omegap.N()!=sigmap.N() but omega.N()==sigma.N(), then to calculate we
    // need to calculate the conjugate transpose of the unit tensor matrix and then invert and multiply 
    // by factor to obtain desired matrix
    if (Nn_zero)
      {
        u3::U3 omegap,omega,omega0;
        int rp, r,rho0;
        HalfInt S0, T0, Sp, Tp, S, T;
        spncci::UnitTensor unit_tensor;

        std::tie (omegap,omega,unit_tensor,rho0)=unit_tensor_u3_sector.Key();
        std::tie (omega0,S0,T0,rp,Sp,Tp,r,S,T)=unit_tensor.Key();
        spncci::UnitTensorU3Sector unit_tensor_calc_u3_sector
          =spncci::UnitTensorU3Sector(omega,omegap,UnitTensor(u3::Conjugate(omega0),S0,T0,r,S,T,rp,Sp,Tp),rho0);

        //  Call UnitTensorMatrix function to calculate the Unit Tensor sub matrix for the v'v 
        //  corresponding to omega' and omega
        temp_matrix=spncci::UnitTensorMatrix(
                                             u_coef_cache,
                                             lgi,lgip,
                                             sector_NpN2,sector_NpN4,
                                             irrep,irrepp,
                                             unit_tensor_calc_u3_sector
                                             );
        
        double coef=ParitySign(rp+r+ConjugationGrade(omega)+ConjugationGrade(omegap))
              *sqrt(1.*dim(u3::SU3(rp,0))*dim(omega)/(dim(u3::SU3(r,0))*dim(omegap)));

        // if the matrix has non-zero entries,
        if (temp_matrix.any())
          {
            // apply symmtry factors, transpose the matrix and 
            
            unit_tensor_u3_sector_pairs.push_back(UnitTensorU3SectorPair(unit_tensor_u3_sector,coef*temp_matrix.transpose()));
          }
        // else
          // zerocout++;
      }
    // otherwise, directly apply the algorithm
    else 
      {
        //  Call UnitTensorMatrix function to calculate the Unit Tensor sub matrix for the v'v 
        //  corresponding to omega' and omega
        temp_matrix=spncci::UnitTensorMatrix(
                                             u_coef_cache,
                                             lgip,lgi,
                                             sector_NpN2,sector_NpN4,
                                             irrepp,irrep,
                                             unit_tensor_u3_sector
                                             );
      
    // If temp_matrix is non-zero, add unit tensor sub matrix into the unit_tensor_rme_map
        if (temp_matrix.any())
          {
            unit_tensor_u3_sector_pairs.push_back(UnitTensorU3SectorPair(unit_tensor_u3_sector,temp_matrix));
          }
        }
        // std::cout<<"zero count  "<<zerocout<<std::endl;
  #ifdef VERBOSE
  std::cout<<"Number of pairs  "<<unit_tensor_u3_sector_pairs.size()<<std::endl;
  std::cout<<"Exiting GenerateUnitTensorU3Sector"<<std::endl;
  #endif
  }


void GenerateNpNSector(const std::pair<int,int> NpN_pair, 
                        const spncci::LGI& lgip,
                        const spncci::LGI& lgi,
                        u3::UCoefCache& u_coef_cache,
                        std::map<
                                std::pair<int,int>,
                                spncci::UnitTensorSectorsCache
                                >& unit_tensor_rme_map,
                        std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector> >& unit_tensor_NpN_sector_map)
{
          ////////////////////// NSectors///////////////////////////////////////////////////////
          const sp3r::Sp3RSpace& irrepp=lgip.Sp3RSpace();
          const sp3r::Sp3RSpace& irrep=lgi.Sp3RSpace();

          //std::pair<int,int> NpN_pair(Nnp,Nn);
          const std::vector<spncci::UnitTensorU3Sector>& unit_U3Sector_vector=unit_tensor_NpN_sector_map[NpN_pair];
          // taking care of two possible cases
          int Nnp=NpN_pair.first;
          int Nn=NpN_pair.second;
          bool Nn_zero=(Nnp!=0 && Nn==0);

          std::pair<int,int> NpN2=Nn_zero?
            std::pair<int,int>(Nn,Nnp-2)
            :std::pair<int,int>(Nnp,Nn-2);
          
          std::pair<int,int> NpN4=Nn_zero?
            NpN4=std::pair<int,int>(Nn-2,Nnp-2)
            :std::pair<int,int>(Nnp-2,Nn-2);

          UnitTensorSectorsCache& sector_NpN2=unit_tensor_rme_map[NpN2];
          UnitTensorSectorsCache& sector_NpN4=unit_tensor_rme_map[NpN4];


          int sector_count = 0;  // debugging variable
          #ifdef VERBOSE
          std::cout<<"Begin generating sectors "<< unit_U3Sector_vector.size()<<std::endl;
          #endif

          #pragma omp parallel reduction(+:sector_count)
          {
            
            #ifdef VERBOSE_OMP
            #pragma omp single
            std::cout << "omp_get_num_threads " << omp_get_num_threads() << std::endl;
            #endif

            // private storage of generated sectors
            std::vector< spncci::UnitTensorU3SectorPair > u3sector_pairs;
            // generate sectors
            int dist=0;

            #pragma omp for schedule(runtime)
            for (int i=0; i<unit_U3Sector_vector.size(); i++)
              {
                const spncci::UnitTensorU3Sector& unit_tensor_u3_sector=unit_U3Sector_vector[i];
                GenerateUnitTensorU3Sector(u_coef_cache, unit_tensor_u3_sector, lgip, lgi, sector_NpN2, sector_NpN4, irrepp, irrep, Nn_zero, u3sector_pairs);
              }

            // save out sectors
            #pragma omp critical
            {
              #ifdef VERBOSE_OMP
              std::cout << "  Saving sectors from thread " << omp_get_thread_num() << std::endl;
              #endif

              unit_tensor_rme_map[NpN_pair].insert(u3sector_pairs.begin(),u3sector_pairs.end());
              // sector_count += u3sector_pairs.size();
              // for (int j=0; j<u3sector_pairs.size(); j++)
              //   {
              //     std::cout << " " << u3sector_pairs[j].second<<std::endl;

              //     unit_tensor_rme_map[NpN_pair].insert(u3sector_pairs[j]);
              //   }
            }
          }  // omp parallel
          #ifdef VERBOSE
          std::cout<<"Finish generating sectors "<<std::endl;
          #endif
          // remove invalid NpN4 key that was inserted into map.          
          if ((Nnp+Nn)==2)
            unit_tensor_rme_map.erase(NpN4);
}

  void GenerateUnitTensorMatrix(
                                // single particle cutoff, relative particle <= 2*N1b+Nn
                                int N1b,
                                // boson number cutoff
                                int Nmax, 
                                // a given spncci sector pair given as index pair  from global list lgi_vector 
                                std::pair<int,int> lgi_pair,
                                // Address to map with list of unit tensor labels with key N0 
                                std::map< int,std::vector<spncci::UnitTensor>>& unit_sym_map,
                                // Address to map of map unit tensor matrix elements keyed by unit tensor labels for key LGI pair
                                std::map<
                                std::pair<int,int>,
                                spncci::UnitTensorSectorsCache
                                >& unit_tensor_rme_map
                                )
  // Generates all unit tensor matrix matrices between states in the irreps of lgi_pair
  // The unit tensors are stored in the map of a map unit_tensor_rme_map which has key lgi_pair to 
  // a map with key std::pair<Nnp,Nn> and value map(matrix labels for w'w sector, matrix) 
  { 
    #ifdef VERBOSE
    std::cout<<"Entering GenerateUnitTensorMatrix"<<std::endl;
    #endif
    std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector> > unit_tensor_NpN_sector_map;


    // extract LGI labels from pair
    const spncci::LGI& lgip=lgi_vector[lgi_pair.first];
    const spncci::LGI& lgi=lgi_vector[lgi_pair.second];
    const sp3r::Sp3RSpace& irrepp=lgip.Sp3RSpace();
    const sp3r::Sp3RSpace& irrep=lgi.Sp3RSpace();

    // Generate list of labels for the unit tensor u(3) sectors of the unit tensor matrices
    // function also generates list of UCoefLabels for all the U coeffients that will be precalculate 
    // and stored in a hash table 
    GenerateUnitTensorU3SectorLabels(N1b,Nmax,lgi_pair,unit_sym_map,unit_tensor_NpN_sector_map);


    // collect full list of unit tensor U(3) sectors for use in U coefficient caching
    std::vector<spncci::UnitTensorU3Sector> unit_tensor_u3_sector_vector;
    std::vector<u3::UCoefLabels> u_coef_labels_vector;

    for (auto it = unit_tensor_NpN_sector_map.begin(); it != unit_tensor_NpN_sector_map.end(); ++it)
      {
        const std::vector<spncci::UnitTensorU3Sector>& unit_tensor_u3_sectors_for_NpN = it->second;
        unit_tensor_u3_sector_vector.insert(
                                      unit_tensor_u3_sector_vector.end(),
                                      unit_tensor_u3_sectors_for_NpN.begin(),
                                      unit_tensor_u3_sectors_for_NpN.end()
                                      );
      }    
    int num_unit_tensor_sectors=0;
    u3::UCoefCache u_coef_cache;
    ////////////////////////////////////////////////////////////////////////////////////
    // Looping over NpN subspaces 
    ////////////////////////////////////////////////////////////////////////////////////
    for (int Nsum=2; Nsum<=2*Nmax; Nsum+=2)
      for (int Nnp=0; Nnp<=std::min(Nsum,Nmax); Nnp+=2)
        {
          if((Nnp+lgip.Nex)>Nmax)
            continue;

          int Nn=Nsum-Nnp;
          if ((Nn+lgi.Nex)>Nmax)
            continue;
              
          GenerateNpNSector( NpN_pair,lgip,lgi,u_coef_cache,unit_tensor_rme_map,unit_tensor_NpN_sector_map); 
          int num=unit_tensor_rme_map[NpN_pair].size();
          num_unit_tensor_sectors+=num;              
        } 

    std::cout<<"number of unit tensor sectors "<< num_unit_tensor_sectors<<std::endl;

    #ifdef VERBOSE
    std::cout<<"Exiting GenerateUnitTensorMatrix"<<std::endl; 
    #endif     
  }  // end function
        
} // End namespace 
  
          

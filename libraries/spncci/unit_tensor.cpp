/****************************************************************
  unit_tensor.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/unit_tensor.h"

#include <omp.h>

#include "cppformat/format.h"
#include "sp3rlib/u3coef.h"

extern double zero_threshold;

namespace spncci
{
  typedef std::pair<UnitTensorU3Sector, Eigen::MatrixXd> UnitTensorU3SectorPair;

  std::string UnitTensorU3Sector::Str() const
  {
    std::ostringstream ss;
    ss << omegap_.Str() << " " << omega_.Str() << " " << tensor_.Str() << " " << rho0_;
    return ss.str();
  }

  std::pair<int,int> GetNSectorIndices(const int Nmax, const int irrep_size, const int Nn, std::vector<int>& NPartition)
  // 
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

  ////////////////////////////////////////////////////////////////////////////////////
  void 
  GenerateUnitTensorU3SectorLabels(
    int N1b,
    int Nmax,
    std::pair<int,int>  sp_irrep_pair,
    const spncci::SpNCCISpace& sp_irrep_vector,
    std::map< int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>>& unit_tensor_labels_map,
    std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>>& unit_tensor_NpN_sector_map
    )
  {   
    #ifdef VERBOSE
    std::cout<<"Entering GenerateU3SectorLabels"<<std::endl;
    #endif

    // Extracting SpNCCIIrrepFamily labels from pair
    const spncci::SpNCCIIrrepFamily& sp_irrepp=sp_irrep_vector[sp_irrep_pair.first];
    const spncci::SpNCCIIrrepFamily& sp_irrep=sp_irrep_vector[sp_irrep_pair.second];
    u3::U3 sigmap=sp_irrepp.sigma();
    u3::U3 sigma=sp_irrep.sigma(); 
    const sp3r::Sp3RSpace& irrepp=sp_irrepp.Sp3RSpace();
    const sp3r::Sp3RSpace& irrep=sp_irrep.Sp3RSpace();

    int irrep_size=irrep.size();
    int irrepp_size=irrepp.size();

    // partition irreps by Nn and Nnp.  Each int in vector corresponds to the start of the next N space 
    std::vector<int> NpPartition=sp3r::PartitionIrrepByNn(irrepp, Nmax);
    std::vector<int> NPartition=sp3r::PartitionIrrepByNn(irrep, Nmax);
    ////////////////////////////////////////////////////////////////////////////////////
    // Looping over omega' and omega subspaces 
    ////////////////////////////////////////////////////////////////////////////////////
    // std::cout<<"iterating"<<std::endl;
    bool conj_sector=sp_irrep_pair.first>sp_irrep_pair.second;
    int Nnp_max=conj_sector?0:Nmax;

    // Loop over Nnp+Nn starting from 2, (Nnp+Nn=0 accounted for elsewhere)
    // std::cout<<"Loop over omega subspaces"<<std::endl;
    for (int Nsum=2; Nsum<=2*Nmax; Nsum+=2)
      for (int Nnp=0; Nnp<=std::min(Nsum,Nnp_max); Nnp+=2)
        {
          int Nn=Nsum-Nnp;

          if((Nnp+sp_irrepp.Nex())>Nmax)
            continue;

          if ((Nn+sp_irrep.Nex())>Nmax)
            continue; 

          int N0=sp_irrepp.Nex()+Nnp-sp_irrep.Nex()-Nn;

          std::pair<int,int> NpN_pair(Nnp,Nn);
          // Selecting section of spaces to iterate over
          int ip_min, ip_max, i_min, i_max;
          // std::cout<<"iterator"<<std::endl;
          std::tie (i_min,i_max)=GetNSectorIndices(Nmax, irrep_size, Nn, NPartition);
          std::tie (ip_min,ip_max)=GetNSectorIndices(Nmax, irrepp_size, Nnp, NpPartition);
          // Get set of operator labels for given omega'omega sector
          // taking into account different Nsigmas 
          std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& N0_operator_set
            =unit_tensor_labels_map[N0];

          // iterate over omega' subspace

          for(int ip=ip_min; ip<=ip_max; ip++ )
            {
              u3::U3 omegap=irrepp.GetSubspace(ip).GetSubspaceLabels();    
              // iterate over omega subspace
              for(int i=i_min; i<=i_max; i++ )
                {
                  u3::U3 omega=irrep.GetSubspace(i).GetSubspaceLabels();
                  // Iterating over the operator labels             
                  for (auto& unit_tensor : N0_operator_set)
                    {                     
                      //unpack unit_tensor labels 
                      u3::SU3 x0(unit_tensor.x0());
                      // HalfInt S0(unit_tensor.S0());
                      // HalfInt T0(unit_tensor.T0());
                      int rbp=unit_tensor.bra().eta();
                      int rb=unit_tensor.ket().eta();

                      if(rb>(2*N1b+Nn+sp_irrep.Nex()))
                        continue;
                      if(rbp>(2*N1b+Nnp+sp_irrepp.Nex()))
                        continue;
                      int rho0_max=OuterMultiplicity(omega.SU3(),x0,omegap.SU3());

                      ///////////////////////////////////////////////////////////////////////////////////////
                      // Iterating over outer multiplicity
                      for (int rho0=1; rho0<=rho0_max; rho0++)
                        {                   
                          spncci::UnitTensorU3Sector unit_U3Sector(omegap,omega,unit_tensor,rho0);
                          unit_tensor_NpN_sector_map[NpN_pair].push_back(unit_U3Sector);
                        }
                    }
                }
            }
        }    
      // std::cout<<"end subspaces"<<std::endl;
    #ifdef VERBOSE
    std::cout<<"Exiting GenerateU3SectorLabels"<<std::endl;
    #endif
  }

////////////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd 
    UnitTensorMatrix(
    u3::UCoefCache& u_coef_cache,
    u3::PhiCoefCache& phi_coef_cache,
    std::unordered_map<u3::U3,vcs::MatrixCache, boost::hash<u3::U3>> k_matrix_map,
     // SpNCCIIrrepFamily pair sector 
    const spncci::SpNCCIIrrepFamily& sp_irrepp,
    const spncci::SpNCCIIrrepFamily& sp_irrep,
    const std::pair<int,int>& lgi_mult,
     // vector of addresses to relevant Np,N sectors of unit tensor matrix
     spncci::UnitTensorSectorsCache& sector_NpN2,
     spncci::UnitTensorSectorsCache& sector_NpN4,
     spncci::UnitTensorU3Sector unit_labels
     )
  {
    #ifdef VERBOSE
    std::cout<<"Entering UnitTensorMatrix"<<std::endl;
    #endif
    // v',v
    int N1b=2;


    // std::cout<<"unit labels "<<unit_labels.Str()<<std::endl;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Set up for calculation 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Extracting labels 
    u3::U3 omegap, omega;
    u3shell::RelativeUnitTensorLabelsU3ST tensor;
    int rho0;
    std::tie(omegap,omega,tensor,rho0) = unit_labels.Key();
    int multp=lgi_mult.first;
    int mult=lgi_mult.second;

    const sp3r::Sp3RSpace& irrepp=sp_irrepp.Sp3RSpace();
    const sp3r::Sp3RSpace& irrep=sp_irrep.Sp3RSpace();
    sp3r::U3Subspace u3_subspacep = irrepp.LookUpSubspace(omegap);
    sp3r::U3Subspace u3_subspace  = irrep.LookUpSubspace(omega);

    int Nn=int(omega.N()-sp_irrep.sigma().N());
    int Nnp=int(omegap.N()-sp_irrepp.sigma().N());

    int dimp=u3_subspacep.size();
    int dim=u3_subspace.size();
    assert(dimp!=0 && dim!=0);

    // unpacking the unit tensor labels 
    u3::SU3 x0=tensor.x0();
    HalfInt S0=tensor.S0();
    HalfInt T0=tensor.T0();
    int rbp=tensor.bra().eta();
    int rb=tensor.ket().eta();

    int rho0_max=u3::OuterMultiplicity(omega.SU3(),x0,omegap.SU3());
 
    // Extracting K matrices for sp_irrep and sp_irrepp from the K_matrix_maps 
    vcs::MatrixCache& K_matrix_map_sp_irrep=k_matrix_map[sp_irrep.sigma()];
    vcs::MatrixCache& K_matrix_map_sp_irrepp=k_matrix_map[sp_irrepp.sigma()];

    Eigen::MatrixXd Kp=K_matrix_map_sp_irrepp[omegap];
    Eigen::MatrixXd K_inv=K_matrix_map_sp_irrep[omega].inverse();

    // Precalculating kronecker products used in sum to calculate unit tensor matrix
    MultiplicityTagged<u3::U3>::vector omegapp_set=KroneckerProduct(omegap, u3::U3(0,0,-2)); 
    MultiplicityTagged<u3::U3>::vector omega1_set=KroneckerProduct(omega, u3::U3(0,0,-2));
    MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(x0, u3::SU3(2,0));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Calculate unit tensor matrix
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Eigen::MatrixXd unit_tensor_matrix=Eigen::MatrixXd::Zero(dimp*multp,dim*mult);

    // summing over omega1
    for (auto& omega1_tagged : omega1_set)
      {	
        u3::U3 omega1(omega1_tagged.irrep);			
        //check that omega1 in irrep  
        if (not irrep.ContainsSubspace(omega1))
          continue;
        // omega1 sector
        sp3r::U3Subspace u3_subspace1=irrep.LookUpSubspace(omega1);
        int dim1=u3_subspace1.size();
        // Look up K1 matrix (dim v1, v1)
        Eigen::MatrixXd K1=K_matrix_map_sp_irrep[omega1];
        // Initializing unit tensor matrix with dim. v' v1
        Eigen::MatrixXd unit_matrix
          =Eigen::MatrixXd::Zero(dimp*multp,dim1*mult);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Matrix of B*U coefs with dim v1 and v
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
        Eigen::MatrixXd BU(dim1,dim);
        //iterating over (n,rho)
        for (int m=0; m<dim; m++)
          {
            MultiplicityTagged<u3::U3> n_rho=u3_subspace.GetStateLabels(m);
            u3::U3 n(n_rho.irrep);
            // iterate over (n1,rho1)
            for (int m1=0; m1<dim1; m1++)
              {
                MultiplicityTagged<u3::U3> n1_rho1=u3_subspace1.GetStateLabels(m1);
                u3::U3 n1(n1_rho1.irrep);
                // check allowed couping
                // (2,0)xn1->n, n1xsigma->omega1 (rho1), omega1x(2,0)->omega, nxsigma->omega(rho) 
                if (u3::OuterMultiplicity(n1.SU3(), u3::SU3(2,0),n.SU3())>0)
                    BU(m1,m)=vcs::BosonCreationRME(n,n1)
                             *u3::UCached(u_coef_cache,u3::SU3(2,0),n1.SU3(),omega.SU3(),sp_irrep.sigma().SU3(),
                                          n.SU3(),1,n_rho.tag,omega1.SU3(),n1_rho1.tag,1);
                else
                  BU(m1,m)=0;
              }	
          }								        
        Eigen::MatrixXd KBUK(dim1,dim);
        KBUK.noalias()=K1*BU*K_inv;

        // std::cout<<"KBUK "<<KBUK<<std::endl;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        //summing over x0'
        for (auto& x0p_tagged : x0p_set)
          {
            u3::SU3 x0p(x0p_tagged.irrep);
            Eigen::MatrixXd unit_matrix_x0=Eigen::MatrixXd::Zero(dimp*multp,dim1*mult);

            // std::cout<<"x0 "<<x0p.Str()<<std::endl;
            int rho0p_max=OuterMultiplicity(omega1.SU3(),x0p,omegap.SU3());
  				  // std::cout<<omega1.Str()<<"  "<<x0p.Str()<<"  "<<rho0p_max<<"  "<<rho0_max<<std::endl;
            // summing over rho0'
            for (int rho0p=1; rho0p<=rho0p_max; rho0p++)
              {
                double coef=0;
                for(int rho0b=1; rho0b<=rho0_max; rho0b++)
                  {
                    //(2,0)xx0->x0p(by construction), 
                    //(2,0)xomega1->omega (by construction),
                    // x0xomega->omegap, (rho0_max)
                    //omega1xx0p->omegap (rho0p_max)
                    // std::cout<<"rhos "<<rho0p<<"  "<<rho0b<<std::endl;
                  
                    // std::cout<<"coef sum"<< omega1.Str()<<"  "
                    // <<u3::UCached(u_coef_cache,x0,u3::SU3(2,0),omegap.SU3(), omega1.SU3(),x0p,1,rho0p,omega.SU3(),1,rho0b)<<std::endl;
                    coef+=u3::PhiCached(phi_coef_cache,omega.SU3(),x0,omegap.SU3(),rho0,rho0b)
                         *u3::UCached(u_coef_cache,x0,u3::SU3(2,0),omegap.SU3(), omega1.SU3(),x0p,1,rho0p,omega.SU3(),1,rho0b);
                  }
                ////////////////////////////////////////////////////////////////////////////////////////////////////////
                // third term
                // sum over omega'', v'' and rho0''
                ////////////////////////////////////////////////////////////////////////////////////////////////////////
                //Initilize 3rd-term-unit-tensor matrix
                Eigen::MatrixXd unit3_matrix=Eigen::MatrixXd::Zero(dimp*multp,dim1*mult);

                // Summing over omega''
                for (auto& omegapp_tagged : omegapp_set)
                  {
                    // A matrix will annihilate bra
                    if(Nnp==0)
                      continue;
                    u3::U3 omegapp(omegapp_tagged.irrep);
                    // std::cout<<"omegapp "<<omegapp.Str()<<std::endl;
                    if (not irrepp.ContainsSubspace(omegapp))
                      continue;
                    
                   // quick trap to check if sectors are found in NpN4, if not found for rho0b'=1, then continue      
                   spncci::UnitTensorU3Sector 
                      unit_tensor_u3_sector_1(omegapp,omega1,tensor,1);
                    if ((sector_NpN4).count(unit_tensor_u3_sector_1)==0)
                      continue;	

                    // omega'' subspace (v'')
                    sp3r::U3Subspace u3_subspacepp=irrepp.LookUpSubspace(omegapp);
                    int dimpp=u3_subspacepp.size();
                    // Obtaining K matrix for omega''
                    Eigen::MatrixXd Kpp_inv=K_matrix_map_sp_irrepp[omegapp].inverse();
                    //Constructing a^\dagger U(3) boson matrix for A matrix
                    Eigen::MatrixXd boson_matrix(dimp,dimpp);
                    for(int vpp=0; vpp<dimpp; vpp++)
                      {
                        MultiplicityTagged<u3::U3> npp_rhopp=u3_subspacepp.GetStateLabels(vpp);
                        const u3::U3& npp=npp_rhopp.irrep;
                        const int& rhopp=npp_rhopp.tag;
                        // std::cout<<"  "<<npp.Str()<<"  "<<rhopp<<std::endl;
                        for(int vp=0; vp<dimp; vp++)
                          {
                            MultiplicityTagged<u3::U3> np_rhop=u3_subspacep.GetStateLabels(vp);
                            const u3::U3& np=np_rhop.irrep;
                            const int& rhop=np_rhop.tag; 
                            // std::cout<<"A matrix"<<std::endl;
                            if (u3::OuterMultiplicity(npp.SU3(), u3::SU3(2,0),np.SU3())>0)
                              boson_matrix(vp,vpp)=
                                vcs::BosonCreationRME(np,npp)
                                *ParitySign(u3::ConjugationGrade(omegap)+u3::ConjugationGrade(omegapp))
                                *u3::UCached(u_coef_cache,u3::SU3(2,0),npp.SU3(),omegap.SU3(),sp_irrepp.sigma().SU3(),
                                                          np.SU3(),1,rhop,omegapp.SU3(),rhopp,1);
                                //*u3::UCached(u_coef_cache,sp_irrepp.sigma().SU3(),npp.SU3(),omegap.SU3(),u3::SU3(2,0), 
                                              // omegapp.SU3(),rhopp,1,np.SU3(),1,rhop);
                            else
                              boson_matrix(vp,vpp)=0;
                          } //end vp
                      } //end vpp
                    // std::cout<<"boson matrix pieces"<<std::endl;
                    // std::cout<<Kp<<"  "<<boson_matrix<<"  "<<Kpp_inv<<std::endl;
                    Eigen::MatrixXd A=Kp*boson_matrix*Kpp_inv;
                    // Unit tensor matrix 
                    Eigen::MatrixXd unit3pp_matrix=Eigen::MatrixXd::Zero(dimpp*multp,dim1*mult);
                    int rho0pp_max=u3::OuterMultiplicity(omega1.SU3(),x0,omegapp.SU3());
                    // Summing over rho0''
                    // std::cout<<"summing over rho0"<<std::endl;
                    for(int rho0pp=1; rho0pp<=rho0pp_max; ++rho0pp)
                      {
                        double coef3=0;
                        for (int rho0bp=1; rho0bp<=rho0p_max; rho0bp++)
                          // omegaxx0->omegapp
                          // x0x(2,0)->x0p (construction)
                          // omegappx(2,0)->omegap (construction)
                          // omegaxx0p->omegap
                            coef3+=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)
                                        *u3::UCached(u_coef_cache,
                                          omega1.SU3(),x0,omegap.SU3(),u3::SU3(2,0),
                                          omegapp.SU3(),rho0pp,1,x0p,1,rho0bp
                                          );    

                        // Retriving unit tensor matrix 
                        spncci::UnitTensorU3Sector unit3_labels(omegapp,omega1,tensor,rho0pp);
                        if(sector_NpN4.count(unit3_labels)>0)
                        // assert(sector_NpN4.count(unit3_labels)>0);
                          unit3pp_matrix+=coef3*sector_NpN4[unit3_labels];
                      } //end rho0pp
                    // matrix product (v',v'')*(v'',v1)

                    for(int i=0; i<multp; ++i)
                      for(int j=0; j<mult; ++j)
                        {
                          // Get target indices 
                          int it=i*dimp;
                          int jt=j*dim1;
                          // Get source indices
                          int is=i*dimpp;
                          int js=j*dim1;

                          unit3_matrix.block(it,jt,dimp,dim1)+=A*unit3pp_matrix.block(is,js,dimpp,dim1);
                          // std::cout<<A<<"  "<<unit3pp_matrix<<std::endl;
                        }
                        // unit3_matrix+=Kp*boson_matrix*Kpp_inv*unit3pp_matrix; //Now A
                  } // end omegapp
                // std::cout<<"unit3_matrix "<<unit3_matrix<<std::endl;
                unit_matrix_x0+=unit3_matrix;
                // std::cout<<"term 3  "<<unit_matrix_x0<<std::endl;							
                // //////////////////////////////////////////////////////////////////////////////////////////////////////////
                //first term 
                //////////////////////////////////////////////////////////////////////////////////////////////////////////
                if(u3::OuterMultiplicity(u3::SU3(rbp,0),u3::SU3(0,rb-2),x0p)>0)
                {
                  // std::cout<<rbp<<" "<<rb-2<<"  "<<std::endl;
                  assert((rb-2)>=0);
                  //(rbp,0)x(0,rb)->x0, (0,rb)x(2,0)->(0,rb-2), x0x(2,0)->x0p, (0,rb-2)x(rbp,0)->x0p
                  double 
                  coef1=u3::UCached(u_coef_cache,u3::SU3(rbp,0),u3::SU3(0,rb),x0p, u3::SU3(2,0),x0,1,1,u3::SU3(0,rb-2),1,1)
                        *sqrt((rb+2)*(rb+1.)*u3::dim(x0p)/(2.*u3::dim(x0)));                      
                  
                  u3shell::RelativeStateLabelsU3ST ket(tensor.ket().eta()-2,tensor.ket().S(),tensor.ket().T());
                  
                  // zero initialize unit1_matrix depending on N0 sign
                  Eigen::MatrixXd unit1_matrix=Eigen::MatrixXd::Zero(dimp*multp,dim1*mult);

                  // summing over rho0bp and accumulating sectors in unit1_matrix. 
                  for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp)
                    {
                      spncci::UnitTensorU3Sector unit1_labels;
                      unit1_labels=spncci::UnitTensorU3Sector(omegap,omega1,u3shell::RelativeUnitTensorLabelsU3ST(x0p,S0,T0,tensor.bra(),ket),rho0bp);
                      // Accumulate
                      // std::cout<<"unit 1 labels "<<unit1_labels.Str()<<"  "<<sector_NpN2.count(unit1_labels)<<std::endl;
                      if (sector_NpN2.count(unit1_labels)!=0)  
                          unit1_matrix+=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)*sector_NpN2[unit1_labels];
                    } //end rho0bp

                  // accumulate term 1 sectors in unit matrix sector
                  // std::cout<<"coef1 "<<coef1
                  // <<"  "<<u3::UCached(u_coef_cache,u3::SU3(rbp,0),u3::SU3(0,rb),x0p, u3::SU3(2,0),x0,1,1,u3::SU3(0,rb-2),1,1)
                  // <<"  "<<sqrt((rb+2)*(rb+1.)*u3::dim(x0p)/(2.*u3::dim(x0)))<<std::endl;
                  unit_matrix_x0+=coef1*unit1_matrix;
                  // std::cout<< "unit 1  "<<unit1_matrix<<std::endl;

                } 
                // if(omp_get_thread_num()==7)
                  // std::cout<< "term 1  "<<unit_matrix_x0<<std::endl;
                //////////////////////////////////////////////////////////////////////////////////////////////////////////
                // second term 
                //////////////////////////////////////////////////////////////////////////////////////////////////////////	
                if (u3::OuterMultiplicity(u3::SU3(rbp+2,0),u3::SU3(0,rb),x0p)>0)
                  {
                    // (2,0)x(rbp,0)->(rbp+2,0), (rbp,0)x(0,rb)->x0, (rbp+2,0)x(0,rb)->x0p, x0x(2,0)->x0p
                    double coef2=-1*(rbp+2)*(rbp+1)*sqrt(u3::dim(x0p)/(2.*(rbp+4)*(rbp+3)*u3::dim(x0)))
                            *u3::UCached(u_coef_cache,u3::SU3(2,0),u3::SU3(rbp,0),x0p,u3::SU3(0,rb),
                                          u3::SU3(rbp+2,0),1,1,x0,1,1);

                    spncci::UnitTensorU3Sector unit2_labels;
                    u3shell::RelativeStateLabelsU3ST bra(tensor.bra().eta()+2,tensor.bra().S(),tensor.bra().T());

                    // zero initialize unit1_matrix depending on N0 sign
                    Eigen::MatrixXd unit2_matrix=Eigen::MatrixXd::Zero(dimp*multp,dim1*mult);
                    // std::cout<<fmt::format("{} {} {} {}",dimp,multp,dim1,mult)<<std::endl;
                    // std::cout<<"tensor "<<tensor.Str()<<std::endl;
                    for(int rho0bp=1; rho0bp<=rho0p_max; ++rho0bp)
                      {
                        unit2_labels=spncci::UnitTensorU3Sector(omegap,omega1,u3shell::RelativeUnitTensorLabelsU3ST(x0p,S0,T0,bra,tensor.ket()),rho0bp);
                        // std::cout<<"unit2_labels "<<unit2_labels.Str()<<std::endl;

                        if(sector_NpN2.count(unit2_labels)>0)
                        {
                          // std::cout<<"unit2"<<std::endl<<unit2_matrix<<std::endl;
                          // std::cout<<"sector2"<<std::endl<<sector_NpN2[unit2_labels]<<std::endl;
                          unit2_matrix+=u3::PhiCached(phi_coef_cache,x0p,omega1.SU3(),omegap.SU3(),rho0p,rho0bp)*sector_NpN2[unit2_labels];

                        }
                      }
                    // std::cout<<"accumulate"<<std::endl;
                    // accumulate term 2 sectors in unit matrix sector
                    unit_matrix_x0+=coef2*unit2_matrix;
                  }
                    // std::cout<<"term 2  "<<unit_matrix_x0<<std::endl;
                ////////////////////////////////////////////////////////////////////////////////////////////////////////
                unit_matrix+=coef*unit_matrix_x0;
              } //end rho0p
          } //end sum over x0p
          // std::cout<<"unit matrix "<<unit_matrix<<std::endl;
        // summing over n, rho, n1, rho1, v1
        for(int i=0; i<multp; ++i)
          for(int j=0; j<mult; ++j)
          {
            int it=i*dimp;
            int jt=j*dim;
            int is=i*dimp;
            int js=j*dim1;
            // (v'v1) (v1 v)
            // std::cout<<unit_labels.Str()<<std::endl;
            // std::cout<<"unit_matrix"<<std::endl;
            // std::cout<<unit_matrix<<std::endl<<std::endl;
            // std::cout<<"KBUK"<<std::endl;
            // std::cout<<KBUK<<std::endl<<std::endl;
            // std::cout<<"unit_tensor_matrix"<<std::endl;
            // std::cout<<unit_tensor_matrix<<std::endl<<std::endl;
            // std::cout<<it<<" "<<jt<<" "<<dimp<<" "<<dim<<" "<<is<<" "<<js<<" "<<dimp<<"  "<<dim1<<std::endl;
            unit_tensor_matrix.block(it,jt,dimp,dim)+=unit_matrix.block(is,js,dimp,dim1)*KBUK;
          }
      }// end sum over omega1
      // std::cout<<multp<<"  "<<dimp<<"   "<<mult<<" "<<dim<<std::endl;
      // std::cout<<unit_tensor_matrix<<std::endl<<std::endl;;
    assert(unit_tensor_matrix.cols()!=0 && unit_tensor_matrix.rows()!=0);
    // std::cout<<unit_tensor_matrix <<std::endl;
    #ifdef VERBOSE
    std::cout<<"Exiting UnitTensorMatrix"<<std::endl;
    #endif
    return unit_tensor_matrix;
  } // End function

void 
GenerateNpNSector(
  const std::pair<int,int> NpN_pair, 
  const spncci::SpNCCIIrrepFamily& sp_irrepp,
  const spncci::SpNCCIIrrepFamily& sp_irrep,
  const std::pair<int,int>& lgi_multiplicities,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  std::unordered_map<u3::U3,vcs::MatrixCache, boost::hash<u3::U3>> k_matrix_map,
  std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache>& unit_tensor_rme_map,
  std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector> >& unit_tensor_NpN_sector_map
  )
{
  ////////////////////// NSectors///////////////////////////////////////////////////////
 const std::vector<spncci::UnitTensorU3Sector>& unit_U3Sector_vector
          =unit_tensor_NpN_sector_map[NpN_pair];

  int Nnp=NpN_pair.first;
  int Nn=NpN_pair.second;

  std::pair<int,int> NpN2=std::pair<int,int>(Nnp,Nn-2);
  std::pair<int,int> NpN4=std::pair<int,int>(Nnp-2,Nn-2);

  UnitTensorSectorsCache& sector_NpN2=unit_tensor_rme_map[NpN2];
  UnitTensorSectorsCache& sector_NpN4=unit_tensor_rme_map[NpN4];

  // debugging variable
  int sector_count = 0;  
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
    std::vector<spncci::UnitTensorU3SectorPair> u3_sector_pairs;

    #pragma omp for schedule(runtime)
    for (int i=0; i<unit_U3Sector_vector.size(); i++)
      {
        const spncci::UnitTensorU3Sector& unit_tensor_u3_sector=unit_U3Sector_vector[i];
        // std::cout<<"unit u3 sector "<<unit_tensor_u3_sector.Str()<<std::endl;

        Eigen::MatrixXd temp_matrix
          =spncci::UnitTensorMatrix(
            u_coef_cache,phi_coef_cache,k_matrix_map,sp_irrepp,sp_irrep, lgi_multiplicities,
            sector_NpN2,sector_NpN4,unit_tensor_u3_sector
            );
        // std::cout<<"temp_matrix"<<std::endl;
        // std::cout<<temp_matrix<<std::endl;
        // If temp_matrix is non-zero, add unit tensor sub matrix into the unit_tensor_rme_map
        // if (temp_matrix.any())
        
        if (not CheckIfZeroMatrix(temp_matrix, zero_threshold))
            u3_sector_pairs.push_back(UnitTensorU3SectorPair(unit_tensor_u3_sector,temp_matrix));
      }
    // save out sectors
    #pragma omp critical
    {
      #ifdef VERBOSE_OMP
      std::cout << "  Saving sectors from thread " << omp_get_thread_num() << std::endl;
      #endif
      unit_tensor_rme_map[NpN_pair].insert(u3_sector_pairs.begin(),u3_sector_pairs.end());
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
  if ((Nnp==0)||(Nn==0))
    unit_tensor_rme_map.erase(NpN4);
}

  void 
  GenerateUnitTensorMatrix(
    int N1b,
    int Nmax, 
    std::pair<int,int> sp_irrep_pair,
    const spncci::SpNCCISpace& sp_irrep_vector,
    u3::UCoefCache& u_coef_cache,
    u3::PhiCoefCache& phi_coef_cache,
    std::unordered_map<u3::U3,vcs::MatrixCache, boost::hash<u3::U3>> k_matrix_map,
    // std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>>& unit_tensor_NpN_sector_map,
    std::map< int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>>& unit_tensor_labels_map,
    spncci::UnitTensorMatricesByIrrepFamily& sp_irrep_unit_tensor_rme_map
    // std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache>& unit_tensor_rme_map
    )
  // Generates all unit tensor matrix matrices between states in the irreps of sp_irrep_pair
  // The unit tensors are stored in the map of a map unit_tensor_rme_map which has key 
  // sp_irrep_pair to a map with key std::pair<Nnp,Nn> and value map(matrix labels for 
  // w'w sector, matrix) 
  { 
    #ifdef VERBOSE
    std::cout<<"Entering GenerateUnitTensorMatrix"<<std::endl;
    #endif

    std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache>& unit_tensor_rme_map
        =sp_irrep_unit_tensor_rme_map[sp_irrep_pair];

    // extract SpNCCIIrrepFamily labels from pair
    const spncci::SpNCCIIrrepFamily& sp_irrepp=sp_irrep_vector[sp_irrep_pair.first];
    const spncci::SpNCCIIrrepFamily& sp_irrep=sp_irrep_vector[sp_irrep_pair.second];
    const sp3r::Sp3RSpace& irrepp=sp_irrepp.Sp3RSpace();
    const sp3r::Sp3RSpace& irrep=sp_irrep.Sp3RSpace();
    const HalfInt& Sp=sp_irrepp.S();
    const HalfInt& S=sp_irrep.S();
    std::pair<int,int> lgi_multiplicities(
      sp_irrep_vector[sp_irrep_pair.first].gamma_max(),
      sp_irrep_vector[sp_irrep_pair.second].gamma_max()
      );
    // std::cout<<sp_irrep_vector[sp_irrep_pair.first].Str()<<"  "<< sp_irrep_vector[sp_irrep_pair.second].Str()<<std::endl;
    // std::cout<<lgi_multiplicities.first<<"  "<<lgi_multiplicities.second<<std::endl;
    ////////////////////////////////////////////////////////////////////////////////////
    // Looping over NpN subspaces 
    ////////////////////////////////////////////////////////////////////////////////////    
    int num_unit_tensor_sectors=0;
    int Np_truncate=Nmax-sp_irrepp.Nex();
    int N_truncate=Nmax-sp_irrep.Nex();
    if(Np_truncate<0 || N_truncate<0)
      return;
    ////////////////////////////////////////////////////////////////////////////////////
    // Compute sector for conjugate subspaces (Nn=0, Nnp!=0)
    ////////////////////////////////////////////////////////////////////////////////////    
    std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>> unit_tensor_NpN_sector_map_conj;
    // Temporary container
    // std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache> unit_tensor_rme_map_conj;
    // Get reverse pair
    std::pair<int,int> sp_irrep_pair_conj(sp_irrep_pair.second,sp_irrep_pair.first);
    // Get labels, should only be for Nnp=0 since iconj>jconj
    std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache>& unit_tensor_rme_map_conj
      =sp_irrep_unit_tensor_rme_map[sp_irrep_pair_conj];

    GenerateUnitTensorU3SectorLabels(
      N1b,Nmax,sp_irrep_pair_conj,sp_irrep_vector,
      unit_tensor_labels_map,unit_tensor_NpN_sector_map_conj);
    // Swap multiplicity labels 
    std::pair<int,int> lgi_mult_conj(
      sp_irrep_vector[sp_irrep_pair.second].gamma_max(),
      sp_irrep_vector[sp_irrep_pair.first].gamma_max()
      );
    // For each NnpNn sector
    for(int Nn=2; Nn<=Np_truncate; ++Nn)
      {
        std::pair<int,int> NpN_pair(0,Nn);
        spncci::GenerateNpNSector(
            NpN_pair,sp_irrep,sp_irrepp,lgi_mult_conj,
            u_coef_cache, phi_coef_cache, k_matrix_map,
            unit_tensor_rme_map_conj,unit_tensor_NpN_sector_map_conj
            ); 
      }
    ////////////////////////////////////////////////////////////////////////////////////
    // Populate map with conjugated sectors
    ////////////////////////////////////////////////////////////////////////////////////          
    for(auto it=unit_tensor_rme_map_conj.begin(); it!=unit_tensor_rme_map_conj.end(); ++it)
      {
        std::pair<int,int>NnpNn(it->first.second,it->first.first);
        if(NnpNn==std::pair<int,int>(0,0))
          continue;
        const spncci::UnitTensorSectorsCache& cache=it->second;
        for(auto it2=cache.begin(); it2!=cache.end(); ++it2)
          {
            u3::U3 omegap,omega;
            u3shell::RelativeUnitTensorLabelsU3ST tensor;
            int rho0;
            std::tie(omega,omegap,tensor,rho0)=it2->first.Key();
            //Swapping unit tensor labels;
            int rp=tensor.ket().eta();
            int r=tensor.bra().eta();
            HalfInt Sb=tensor.bra().S();
            HalfInt Tb=tensor.bra().T();
            HalfInt Sbp=tensor.ket().S();
            HalfInt Tbp=tensor.ket().T();
            // Conjugation phase
            //Conjugating matrix element
            double coef=ParitySign(rp+r+ConjugationGrade(omega)+ConjugationGrade(omegap)+ConjugationGrade(tensor.x0()))
                        *sqrt(1.*dim(u3::SU3(rp,0))*dim(omega)/(dim(u3::SU3(r,0))*dim(omegap)));
            coef*=ParitySign(S+tensor.S0()+Sp+Sb+Sbp+Tb+Tbp)
                  *sqrt(1.*am::dim(S)*am::dim(Sbp)*am::dim(Tbp)/am::dim(Sp)/am::dim(Sb)/am::dim(Tb));
            // un-conjugated sector labels
            spncci::UnitTensorU3Sector sector(omegap,omega,u3shell::Conjugate(tensor),rho0);
            // std::cout<<it2->second<<std::endl;
            Eigen::MatrixXd temp=it2->second;
            unit_tensor_rme_map[NnpNn][sector]=coef*temp.transpose();
            // std::cout<<sector.Str()<<std::endl;
            // std::cout<<temp<<"  "<<coef<<"  "<<unit_tensor_rme_map[NnpNn][sector]<<std::endl;
          }
      }
    ////////////////////////////////////////////////////////////////////////////////////
    // Nn>0
    ////////////////////////////////////////////////////////////////////////////////////    
    std::map<std::pair<int,int>,std::vector<spncci::UnitTensorU3Sector>> unit_tensor_NpN_sector_map;
    GenerateUnitTensorU3SectorLabels(
      N1b,Nmax,sp_irrep_pair,sp_irrep_vector,
      unit_tensor_labels_map,unit_tensor_NpN_sector_map);

    for (int Nsum=2; Nsum<=2*Nmax; Nsum+=2)
      for (int Nn=2; Nn<=std::min(Nsum,N_truncate); Nn+=2)
        {
          // Check Nmax constrain
          int Nnp=Nsum-Nn;
          if(Nnp==0)
            continue;
          if (Nnp>Np_truncate)
            continue;
         //////////////////////////////////////////////////////////////////////////////     
          //  Compute rmes for NSectors
          //////////////////////////////////////////////////////////////////////////////    
          std::pair<int,int> NpN_pair(Nnp,Nn);              
          GenerateNpNSector(
            NpN_pair,sp_irrepp,sp_irrep,lgi_multiplicities,
            u_coef_cache, phi_coef_cache,k_matrix_map,
            unit_tensor_rme_map,unit_tensor_NpN_sector_map
            ); 
          // diagnostic output
          int num=unit_tensor_rme_map[NpN_pair].size();
          num_unit_tensor_sectors+=num;
          // std::cout<<Nnp<<"  "<<Nn<<"  "<<"num  "<< num_unit_tensor_sectors<<std::endl;
        } 
    ////////////////////////////////////////////////////////////////////////////////////    

  #ifdef VERBOSE
  std::cout<<"number of unit tensor sectors "<< num_unit_tensor_sectors<<std::endl;
  std::cout<<"Exiting GenerateUnitTensorMatrix"<<std::endl; 
  #endif     
}  // end function


} // End namespace 
  
          

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <map>
#include <mpi.h>

#include <su3.h>
#include "wigxjpf.h"

std::vector<std::array<int,3>> KroneckerProduct(const int lm1, const int mu1, const int lm2, const int mu2){
  std::vector<std::array<int,3>> product;
  for(int lm=0; lm<=lm1+lm2+std::min(mu2,lm1+mu1); lm++){
    for(int mu=0; mu<=mu1+mu2+std::min(lm1,lm2); mu++){
      int rhomax=su3::mult(lm1,mu1,lm2,mu2,lm,mu);
      if(rhomax!=0)product.push_back({lm,mu,rhomax});
    }
  }
  return product;
}

int su3dim(const int lm, const int mu){
  return ((lm+1)*(mu+1)*(lm+mu+2))/2;
}

int main(int argc, char **argv){

int JJ_bra=2;
int JJ_ket=2;
int N1vp=1;
int N1vn=1;
int N0=2;

if(argc!=2){
  std::cout<<"Usage: "<<argv[0]<<" Nmax"<<std::endl;
  return EXIT_FAILURE;
}

const int Nmax=std::stoi(argv[1]);

//**********************************************
// Read wave function
//**********************************************
std::map<std::array<int,4>,std::map<int,std::map<std::array<int,8>,std::map<int,double>>>> amplitude_by_labels;
// amplitude_by_labels[lm_omega,mu_omega,L,kappamax][SS][gamma,Nex_sigma,lm_sigma,mu_sigma,SSp,SSn,upsilon,Nex_omega][kappa]
std::ifstream wffile("wavefunction.dat");
if(!wffile){
  std::cout<<"Could not open wavefunction.dat file!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int gamma,Nex_sigma,lm_sigma,mu_sigma,SSp,SSn,SS,upsilon,Nex_omega,lm_omega,mu_omega,kappa,L;
  double amplitude;
  wffile>>gamma>>Nex_sigma>>lm_sigma>>mu_sigma>>SSp>>SSn>>SS>>upsilon>>Nex_omega>>lm_omega>>mu_omega>>kappa>>L>>amplitude;
  if(wffile){
    int kappamax=su3::kmax(lm_omega,mu_omega,L);
    std::array<int,4> key1={lm_omega,mu_omega,L,kappamax};
    std::array<int,8> key2={gamma,Nex_sigma,lm_sigma,mu_sigma,SSp,SSn,upsilon,Nex_omega};
    std::map<std::array<int,4>,std::map<int,std::map<std::array<int,8>,std::map<int,double>>>>::iterator it1=amplitude_by_labels.find(key1);
    if(it1==amplitude_by_labels.end()){
      std::map<int,std::map<std::array<int,8>,std::map<int,double>>> value1;
      std::map<std::array<int,8>,std::map<int,double>> value2;
      std::map<int,double> value3;
      value3[kappa]=amplitude;
      value2[key2]=value3;
      value1[SS]=value2;
      amplitude_by_labels[key1]=value1;
    }else{
      std::map<int,std::map<std::array<int,8>,std::map<int,double>>>::iterator it2=amplitude_by_labels[key1].find(SS);
      if(it2==amplitude_by_labels[key1].end()){
        std::map<std::array<int,8>,std::map<int,double>> value2;
        std::map<int,double> value3;
        value3[kappa]=amplitude;
        value2[key2]=value3;
        amplitude_by_labels[key1][SS]=value2;
      }else{
        std::map<std::array<int,8>,std::map<int,double>>::iterator it3=amplitude_by_labels[key1][SS].find(key2);
	if(it3==amplitude_by_labels[key1][SS].end()){
          std::map<int,double> value3;
          value3[kappa]=amplitude;
          amplitude_by_labels[key1][SS][key2]=value3;
        }else{
          amplitude_by_labels[key1][SS][key2][kappa]=amplitude;
	}
      }
    }
  }else{
    break;
  }
}
wffile.close();

//**********************************************
// Read TBDRMEs between SpNCCI basis states
//**********************************************
std::map<std::array<int,37>,double> proton_TBDRME_by_labels,neutron_TBDRME_by_labels,pn_TBDRME_by_labels;
std::ifstream btbdfile("basistbdrmes.dat");
if(!btbdfile){
  std::cout<<"Could not open basistbdrmes.dat file!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int gammap,Nex_sigmap,lm_sigmap,mu_sigmap,SSp_bra,SSn_bra,SS_bra,upsilonp,Nex_omegap,lm_omegap,mu_omegap,
      gamma,Nex_sigma,lm_sigma,mu_sigma,SSp_ket,SSn_ket,SS_ket,upsilon,Nex_omega,lm_omega,mu_omega,
      N1,N2,N3,N4,lmf,muf,Sf,lmi,mui,Si,lm0,mu0,SS0,rho0,rho,Tz;
  double rme;
  btbdfile>>gammap>>Nex_sigmap>>lm_sigmap>>mu_sigmap>>SSp_bra>>SSn_bra>>SS_bra>>upsilonp>>Nex_omegap>>lm_omegap>>mu_omegap
          >>gamma>>Nex_sigma>>lm_sigma>>mu_sigma>>SSp_ket>>SSn_ket>>SS_ket>>upsilon>>Nex_omega>>lm_omega>>mu_omega
          >>N1>>N2>>N3>>N4>>lmf>>muf>>Sf>>lmi>>mui>>Si>>lm0>>mu0>>SS0>>rho0>>rho>>Tz>>rme;
  if(btbdfile){
    std::array<int,37> key={gammap,Nex_sigmap,lm_sigmap,mu_sigmap,SSp_bra,SSn_bra,SS_bra,upsilonp,Nex_omegap,lm_omegap,mu_omegap,
                            gamma,Nex_sigma,lm_sigma,mu_sigma,SSp_ket,SSn_ket,SS_ket,upsilon,Nex_omega,lm_omega,mu_omega,
                            N1,N2,N3,N4,lmf,muf,2*Sf,lmi,mui,2*Si,lm0,mu0,SS0,rho0,rho};
    if(Tz==1){
      proton_TBDRME_by_labels[key]=rme;
    }else if(Tz==-1){
      neutron_TBDRME_by_labels[key]=rme;
    }else{
      pn_TBDRME_by_labels[key]=rme;
    }
  }else{
    break;
  }
}
btbdfile.close();

//**********************************************
// Calculate TBDRMEs
//**********************************************
int etamax=std::max(N1vp,N1vn)+Nmax;
int Ntotmax=N0+Nmax;

std::vector<std::array<int,4>> N1N2l1l2;
for(int N1=0; N1<=etamax; N1++){
for(int N2=0; N2<=std::min(etamax,Ntotmax-N1); N2++){
for(int l1=N1%2; l1<=N1; l1+=2){
for(int l2=N2%2; l2<=N2; l2+=2){
  N1N2l1l2.push_back({N1,N2,l1,l2});
}
}
}
}
double factor0=sqrt(double((JJ_bra+1)*(JJ_ket+1)));
/*
int rank,size;
MPI_Status status;
MPI_Comm comm;
comm=MPI_COMM_WORLD;
MPI_Init(NULL, NULL);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
*/
int rank,nproc;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&nproc);

//std::cout<<"Hello from rank "<<rank<<std::endl;

su3::init(100);

//for(int rank=0; rank<N1N2l1l2.size(); rank++){
 int N1=N1N2l1l2[rank][0];
 int N2=N1N2l1l2[rank][1];
 int l1=N1N2l1l2[rank][2];
 int l2=N1N2l1l2[rank][3];
// std::cout<<"Rank "<<rank<<" computes N1 N2 l1 l2 = "<<N1<<" "<<N2<<" "<<l1<<" "<<l2<<std::endl;

 std::map<std::array<int,11>,std::array<double,3>> tbdrmes;
 // Key is {N3,N4,l3,l4,jj1,jj2,jj3,jj4,JJf,JJi,JJ0}
 // Value is {proton, neutron, proton-neutron}
 for(int N3=0; N3<=etamax; N3++){
  int N4min=std::max(N1+N2-N3-Nmax,0);
  if((N1+N2-N3-N4min)%2!=0)N4min++; // so that TBD operator doesn't change parity
  for(int N4=N4min; N4<=std::min(std::min(etamax,Ntotmax-N3),N1+N2-N3+Nmax); N4+=2){
   for(int l3=N3%2; l3<=N3; l3+=2){
    for(int l4=N4%2; l4<=N4; l4+=2){
     for(int jj1=std::abs(2*l1-1); jj1<=2*l1+1; jj1+=2){
      for(int jj2=std::abs(2*l2-1); jj2<=2*l2+1; jj2+=2){
       for(int jj3=std::abs(2*l3-1); jj3<=2*l3+1; jj3+=2){
        for(int jj4=std::abs(2*l4-1); jj4<=2*l4+1; jj4+=2){
         for(int JJf=std::abs(jj1-jj2); JJf<=jj1+jj2; JJf+=2){
          for(int JJi=std::abs(jj3-jj4); JJi<=jj3+jj4; JJi+=2){
           for(int JJ0=std::max(std::abs(JJf-JJi),std::abs(JJ_ket-JJ_bra)); JJ0<=std::min(JJf+JJi,JJ_ket+JJ_bra); JJ0+=2){
            tbdrmes[{N3,N4,l3,l4,jj1,jj2,jj3,jj4,JJf,JJi,JJ0}]={0.0,0.0,0.0};
	   }
	  }
	 }
        }
       }
      }
     }
    }
   }
  }
 }
 for(std::array<int,3> xf : KroneckerProduct(N1,0,N2,0)){
  int phase1=N1+N2+xf[0]+xf[1];
  std::vector<double> su3cgs1;
  su3cgs1.reserve(N1*N2*std::max(xf[0],xf[1]));
  for(int Lf=0; Lf<=xf[0]+xf[1]; Lf++){
   int kappafmax=su3::kmax(xf[0],xf[1],Lf);
   if(kappafmax==0 || std::abs(l1-l2)>Lf || Lf>l1+l2)continue;
   su3::wu3r3w(N1,0,N2,0,xf[0],xf[1],l1,l2,Lf,kappafmax,1,1,1,su3cgs1);
   double factor1=factor0*sqrt(double(2*Lf+1));
   for(int N3=0; N3<=etamax; N3++){
    int N4min=std::max(N1+N2-N3-Nmax,0);
    if((N1+N2-N3-N4min)%2!=0)N4min++; // so that TBD operator doesn't change parity
    for(int N4=N4min; N4<=std::min(std::min(etamax,Ntotmax-N3),N1+N2-N3+Nmax); N4+=2){
     for(std::array<int,3> xi : KroneckerProduct(0,N3,0,N4)){
      int phase2=phase1+N3+N4+xi[0]+xi[1];
      std::vector<double> su3cgs2;
      su3cgs2.reserve(N3*N4*std::max(xi[0],xi[1]));
      for(int l3=N3%2; l3<=N3; l3+=2){
       for(int l4=N4%2; l4<=N4; l4+=2){
        for(int Li=0; Li<=xi[0]+xi[1]; Li++){
         int kappaimax=su3::kmax(xi[0],xi[1],Li);
         if(kappaimax==0 || std::abs(l3-l4)>Li || Li>l3+l4)continue;
         su3::wu3r3w(0,N3,0,N4,xi[0],xi[1],l3,l4,Li,kappaimax,1,1,1,su3cgs2);
         double factor2=factor1*sqrt(double(2*Li+1));
         for(std::array<int,3> x0 : KroneckerProduct(xf[0],xf[1],xi[0],xi[1])){
          int phase3=phase2+x0[0]+x0[1]+x0[2];
          std::vector<double> su3cgs3;
          su3cgs3.reserve(std::max(xf[0],xf[1])*std::max(xi[0],xi[1])*std::max(x0[0],x0[1]));
          std::vector<double> dzu3;
          su3::wzu3(xf[0],xf[1],0,0,x0[0],x0[1],xi[0],xi[1],xf[0],xf[1],xi[0],xi[1],1,x0[2],1,x0[2],dzu3);
          for(int L0=0; L0<=x0[0]+x0[1]; L0++){
           int kappa0max=su3::kmax(x0[0],x0[1],L0);
           if(kappa0max==0 || std::abs(Lf-Li)>L0 || L0>Lf+Li)continue;
           su3::wu3r3w(xf[0],xf[1],xi[0],xi[1],x0[0],x0[1],Lf,Li,L0,kappa0max,kappaimax,kappafmax,x0[2],su3cgs3);
	   double factor3=factor2*sqrt(double(2*L0+1));
           for(std::map<std::array<int,4>,std::map<int,std::map<std::array<int,8>,std::map<int,double>>>>::iterator
	       itbra1=amplitude_by_labels.begin(); itbra1!=amplitude_by_labels.end(); itbra1++){
            int lm_omega_bra=itbra1->first[0];
	    int mu_omega_bra=itbra1->first[1];
	    int L_bra=itbra1->first[2];
	    int kappamax_bra=itbra1->first[3];
	    double factor4=factor3*sqrt(double(2*L_bra+1));
	    int phase4=phase3+lm_omega_bra+mu_omega_bra;
	    for(std::map<std::array<int,4>,std::map<int,std::map<std::array<int,8>,std::map<int,double>>>>::iterator
                itket1=amplitude_by_labels.begin(); itket1!=amplitude_by_labels.end(); itket1++){
             int L_ket=itket1->first[2];
	     if(std::abs(L_ket-L0)>L_bra || L_bra>L_ket+L0)continue;
             int lm_omega_ket=itket1->first[0];
             int mu_omega_ket=itket1->first[1];
	     int rhomax=su3::mult(lm_omega_ket,mu_omega_ket,x0[0],x0[1],lm_omega_bra,mu_omega_bra);
	     if(rhomax==0)continue;
             int kappamax_ket=itket1->first[3];
             std::vector<double> su3cgs4;
             su3cgs4.reserve(std::max(lm_omega_ket,mu_omega_ket)*std::max(x0[0],x0[1])*std::max(lm_omega_bra,mu_omega_bra));
             su3::wu3r3w(lm_omega_ket,mu_omega_ket,x0[0],x0[1],lm_omega_bra,mu_omega_bra,L_ket,L0,L_bra,
	                 kappamax_bra,kappa0max,kappamax_ket,rhomax,su3cgs4);
	     int phase5=phase4-lm_omega_ket-mu_omega_ket;
             for(int jj1=std::abs(2*l1-1); jj1<=2*l1+1; jj1+=2){
              for(int jj2=std::abs(2*l2-1); jj2<=2*l2+1; jj2+=2){
               for(int JJf=std::abs(jj1-jj2); JJf<=jj1+jj2; JJf+=2){
	        for(int SSf=0; SSf<=2; SSf+=2){
                 if(std::abs(2*Lf-SSf)>JJf || JJf>2*Lf+SSf)continue;
                 double factor5=factor4*sqrt(double((jj1+1)*(jj2+1)*(JJf+1)*(SSf+1)))*wig9jj(2*l1,2*l2,2*Lf,1,1,SSf,jj1,jj2,JJf);
                 for(int jj3=std::abs(2*l3-1); jj3<=2*l3+1; jj3+=2){
                  for(int jj4=std::abs(2*l4-1); jj4<=2*l4+1; jj4+=2){
                   for(int JJi=std::abs(jj3-jj4); JJi<=jj3+jj4; JJi+=2){
                    for(int SSi=0; SSi<=2; SSi+=2){
                     if(std::abs(2*Li-SSi)>JJi || JJi>2*Li+SSi)continue;
                     double factor6=factor5*sqrt(double((jj3+1)*(jj4+1)*(JJi+1)*(SSi+1)))*wig9jj(2*l3,2*l4,2*Li,1,1,SSi,jj3,jj4,JJi);
                     for(int JJ0=std::max(std::abs(JJf-JJi),std::abs(JJ_ket-JJ_bra)); JJ0<=std::min(JJf+JJi,JJ_ket+JJ_bra); JJ0+=2){
                      std::array<int,11> keyt={N3,N4,l3,l4,jj1,jj2,jj3,jj4,JJf,JJi,JJ0};
                      for(int SS0=std::abs(SSf-SSi); SS0<=SSf+SSi; SS0+=2){
                       if(std::abs(2*L0-SS0)>JJ0 || JJ0>2*L0+SS0)continue;
                       double factor7=factor6*sqrt(double(SS0+1))*wig9jj(2*Lf,2*Li,2*L0,SSf,SSi,SS0,JJf,JJi,JJ0);
                       for(std::map<int,std::map<std::array<int,8>,std::map<int,double>>>::iterator
		           itbra2=itbra1->second.begin(); itbra2!=itbra1->second.end(); itbra2++){
	                int SS_bra=itbra2->first;
		        double factor8=factor7*sqrt(double(SS_bra+1));
		        for(std::map<int,std::map<std::array<int,8>,std::map<int,double>>>::iterator
                            itket2=itket1->second.begin(); itket2!=itket1->second.end(); itket2++){
                         int SS_ket=itket2->first;
		         if(std::abs(SS_ket-SS0)>SS_bra || SS_bra>SS_ket+SS0)continue;
                         double factor9=factor8*wig9jj(2*L_ket,SS_ket,JJ_ket,2*L0,SS0,JJ0,2*L_bra,SS_bra,JJ_bra);
		         double factord=sqrt(double((SS_ket+1)*su3dim(lm_omega_ket,mu_omega_ket))
		 	    	            /double((SS_bra+1)*su3dim(lm_omega_bra,mu_omega_bra)));
		         int phase=phase5+(SS_ket-SS_bra)/2;
		         for(std::map<std::array<int,8>,std::map<int,double>>::iterator
                             itbra3=itbra2->second.begin(); itbra3!=itbra2->second.end(); itbra3++){
                          int gamma_bra=itbra3->first[0];
		          int Nex_sigma_bra=itbra3->first[1];
		          int lm_sigma_bra=itbra3->first[2];
		          int mu_sigma_bra=itbra3->first[3];
		          int SSp_bra=itbra3->first[4];
		          int SSn_bra=itbra3->first[5];
		          int upsilon_bra=itbra3->first[6];
		          int Nex_omega_bra=itbra3->first[7];
		          int Nn_bra=Nex_omega_bra-Nex_sigma_bra;
		          for(std::map<std::array<int,8>,std::map<int,double>>::iterator
                              itket3=itket2->second.begin(); itket3!=itket2->second.end(); itket3++){
                           int gamma_ket=itket3->first[0];
                           int Nex_sigma_ket=itket3->first[1];
                           int lm_sigma_ket=itket3->first[2];
                           int mu_sigma_ket=itket3->first[3];
                           int SSp_ket=itket3->first[4];
                           int SSn_ket=itket3->first[5];
                           int upsilon_ket=itket3->first[6];
                           int Nex_omega_ket=itket3->first[7];
		           int Nn_ket=Nex_omega_ket-Nex_sigma_ket;
                           for(int rho=1; rho<=su3::mult(lm_omega_ket,mu_omega_ket,x0[0],x0[1],lm_omega_bra,mu_omega_bra); rho++){
			    for(int rho0=1; rho0<=x0[2]; rho0++){
                             double TBDRME_p;
 		             double TBDRME_n;
		             double TBDRME_pn;
                             if(Nn_bra>=Nn_ket){
                              std::array<int,37> key={gamma_bra,Nex_sigma_bra,lm_sigma_bra,mu_sigma_bra,SSp_bra,SSn_bra,
                                                      SS_bra,upsilon_bra,Nex_omega_bra,lm_omega_bra,mu_omega_bra,
                                                      gamma_ket,Nex_sigma_ket,lm_sigma_ket,mu_sigma_ket,SSp_ket,SSn_ket,
                                                      SS_ket,upsilon_ket,Nex_omega_ket,lm_omega_ket,mu_omega_ket,
                                                      N1,N2,N3,N4,xf[0],xf[1],SSf,xi[0],xi[1],SSi,x0[0],x0[1],SS0,rho0,rho};
                              TBDRME_p=proton_TBDRME_by_labels[key];
                              TBDRME_n=neutron_TBDRME_by_labels[key];
                              TBDRME_pn=pn_TBDRME_by_labels[key];
                             }else{
		              TBDRME_p=0.0;
                              TBDRME_n=0.0;
                              TBDRME_pn=0.0;
		              int sign=1;
			      int indz=rho0-1;
                              for(int rho0p=1; rho0p<=x0[2]; rho0p++){
		               sign=-sign;
                               std::array<int,37> key={gamma_ket,Nex_sigma_ket,lm_sigma_ket,mu_sigma_ket,SSp_ket,SSn_ket,
                                                       SS_ket,upsilon_ket,Nex_omega_ket,lm_omega_ket,mu_omega_ket,
		                                       gamma_bra,Nex_sigma_bra,lm_sigma_bra,mu_sigma_bra,SSp_bra,SSn_bra,
                                                       SS_bra,upsilon_bra,Nex_omega_bra,lm_omega_bra,mu_omega_bra,
                                                       N4,N3,N2,N1,xi[1],xi[0],SSi,xf[1],xf[0],SSf,x0[1],x0[0],SS0,rho0p,rho};
    	                       double factorz=sign*dzu3[indz];
                               TBDRME_p+=factorz*proton_TBDRME_by_labels[key];
                               TBDRME_n+=factorz*neutron_TBDRME_by_labels[key];
	                       TBDRME_pn+=factorz*pn_TBDRME_by_labels[key];
			       indz+=x0[2];
                              }
                              TBDRME_p*=factord;
		              TBDRME_n*=factord;
		              TBDRME_pn*=factord;
		              if(phase%2!=0){
		               TBDRME_p=-TBDRME_p;
		               TBDRME_n=-TBDRME_n;
		               TBDRME_pn=-TBDRME_pn;
		              }
		             }
                             for(int ikappa_bra=0; ikappa_bra<kappamax_bra; ikappa_bra++){
		              double factor10=factor9*itbra3->second[ikappa_bra+1];
                              for(int ikappa_ket=0; ikappa_ket<kappamax_ket; ikappa_ket++){
                               double factor11=factor10*itket3->second[ikappa_ket+1];
			       for(int ikappaf=0; ikappaf<kappafmax; ikappaf++){
		                double factor12=factor11*su3cgs1[ikappaf];
			        for(int ikappai=0; ikappai<kappaimax; ikappai++){
                                 double factor13=factor12*su3cgs2[ikappai];
			         for(int ikappa0=0; ikappa0<kappa0max; ikappa0++){
			          double factor14=factor13
			 	      *su3cgs3[rho0-1+x0[2]*ikappaf+x0[2]*kappafmax*ikappai+x0[2]*kappafmax*kappaimax*ikappa0]
				      *su3cgs4[rho-1+rhomax*ikappa_ket+rhomax*kappamax_ket*ikappa0+rhomax*kappamax_ket*kappa0max*ikappa_bra];
		                  tbdrmes[keyt][0]+=factor14*TBDRME_p;
                                  tbdrmes[keyt][1]+=factor14*TBDRME_n;
                                  tbdrmes[keyt][2]+=factor14*TBDRME_pn;
                                 }
                                }
			       }
			      }
			     }
	                    }
			   }
			  }
			 }
		        }
	               }
	              }
                     }
		    }
	           }
	          }
                 }
	        }
	       }
              }
             }
	    }
	   }
	  }
         }
        }
       }
      }
     }
    }
   }
  }
 }
 for(std::map<std::array<int,11>,std::array<double,3>>::iterator it=tbdrmes.begin(); it!=tbdrmes.end(); it++){
   std::cout<<N1<<" "<<l1<<" "<<it->first[4]<<" "<<N2<<" "<<l2<<" "<<it->first[5]<<" "<<it->first[0]<<" "<<it->first[2]<<" "
	    <<it->first[6]<<" "<<it->first[1]<<" "<<it->first[3]<<" "<<it->first[7]<<" "<<it->first[8]<<" "<<it->first[9]<<" "
	    <<it->first[10]<<" "<<it->second[0]<<" "<<it->second[1]<<" "<<it->second[2]<<std::endl;
 }
//}
su3::finalize();

MPI_Finalize();

return EXIT_SUCCESS;
}

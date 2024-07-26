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

double like_int(const int N1,const int l1,const int jj1,const int N2,const int l2,const int jj2,const int N3,const int l3,const int jj3,
	    const int N4,const int l4,const int jj4,const int JJ,std::map<std::array<int,13>,double>& v){
  std::array<int,13> key={N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJ};
  std::map<std::array<int,13>,double>::iterator it=v.find(key);
  if(it!=v.end())return v[key];

  key={N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,JJ};
  it=v.find(key);
  if(it!=v.end()){
    int phase=1+(JJ-jj1-jj2)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }

  key={N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,JJ};
  it=v.find(key);
  if(it!=v.end()){
    int phase=1+(JJ-jj3-jj4)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }

  key={N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,JJ};
  it=v.find(key);
  if(it!=v.end()){
    int phase=(jj1+jj2+jj3+jj4)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }
  
  key={N3,l3,jj3,N4,l4,jj4,N1,l1,jj1,N2,l2,jj2,JJ};
  it=v.find(key);
  if(it!=v.end())return v[key];
  
  key={N4,l4,jj4,N3,l3,jj3,N1,l1,jj1,N2,l2,jj2,JJ};
  it=v.find(key);
  if(it!=v.end()){
    int phase=1+(JJ-jj3-jj4)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }
  
  key={N3,l3,jj3,N4,l4,jj4,N2,l2,jj2,N1,l1,jj1,JJ};
  it=v.find(key);
  if(it!=v.end()){
    int phase=1+(JJ-jj1-jj2)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }
  
  key={N4,l4,jj4,N3,l3,jj3,N2,l2,jj2,N1,l1,jj1,JJ};
  it=v.find(key);
  if(it!=v.end()){
    int phase=(jj1+jj2+jj3+jj4)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }

  return 0.0;
}

double pn_int(const int N1,const int l1,const int jj1,const int N2,const int l2,const int jj2,const int N3,const int l3,const int jj3,
	      const int N4,const int l4,const int jj4,const int JJ,std::map<std::array<int,13>,double>& vpn){
  std::array<int,13> key={N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJ};
  std::map<std::array<int,13>,double>::iterator it=vpn.find(key);
  if(it!=vpn.end())return vpn[key];
 
  key={N3,l3,jj3,N4,l4,jj4,N1,l1,jj1,N2,l2,jj2,JJ};
  it=vpn.find(key);
  if(it!=vpn.end())return vpn[key];
  
  return 0.0;
}

int main(int argc, char **argv){

if(argc!=4){
  std::cout<<"Usage: "<<argv[0]<<" Nmax N1vp N1vn"<<std::endl;
  std::cout<<"where Nmax is Nmax for SpNCCI calculation of target"<<std::endl;
  return EXIT_FAILURE;
}

int Nmax = std::stoi(argv[1]);

//********************************************************
// Read OB rhos
//********************************************************
std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>> obrho_by_labels;
// obrho_by_labels[lm1',mu1',lm1,mu1][SS1',SS1][kappa1',L1',kappa1,L1][N,N',lm0,mu0,S0,rho][proton,neutron]
std::ifstream obrhos("obrho.dat");
if(!obrhos){
  std::cout<<"Could not open file obrho.dat!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int N,Np,lm0,mu0,S0,rho,lmp,mup,kappap,Lp,SSp,lm,mu,kappa,L,SS;
  obrhos>>N>>Np>>lm0>>mu0>>S0>>rho>>lmp>>mup>>kappap>>Lp>>SSp>>lm>>mu>>kappa>>L>>SS;
  if(obrhos){
    std::array<double,2> obrho;
    obrhos>>obrho[0]>>obrho[1];
    std::array<int,4> key1={lmp,mup,lm,mu};
    std::array<int,2> key2={SSp,SS};
    std::array<int,4> key3={kappap,Lp,kappa,L};
    std::array<int,6> key4={N,Np,lm0,mu0,S0,rho};
    std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>>::iterator it1=obrho_by_labels.find(key1);
    if(it1==obrho_by_labels.end()){
      std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>> value1;
      std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>> value2;
      std::map<std::array<int,6>,std::array<double,2>> value3;
      value3[key4]=obrho;
      value2[key3]=value3;
      value1[key2]=value2;
      obrho_by_labels[key1]=value1;
    }else{
      std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>::iterator
	it2=obrho_by_labels[key1].find(key2);
      if(it2==obrho_by_labels[key1].end()){
        std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>> value2;
        std::map<std::array<int,6>,std::array<double,2>> value3;
        value3[key4]=obrho;
        value2[key3]=value3;
	obrho_by_labels[key1][key2]=value2;
      }else{
        std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>::iterator it3=obrho_by_labels[key1][key2].find(key3);
	if(it3==obrho_by_labels[key1][key2].end()){
          std::map<std::array<int,6>,std::array<double,2>> value3;
          value3[key4]=obrho;
	  obrho_by_labels[key1][key2][key3]=value3;
        }else{
          obrho_by_labels[key1][key2][key3][key4]=obrho;
        }
      }
    }
  }else{
    break;
  }
}
obrhos.close();

//*******************************************************************
// Read interaction
//*******************************************************************
std::map<std::array<int,13>,double> vp,vn,vpn;
// Labels are N1 l1 2j1 N2 l2 2j2 N3 l3 2j3 N4 l4 2j4 JJ
std::ifstream h2("h2v0.dat");
if(!h2){
  std::cout<<"Could not open h2v0.dat file!"<<std::endl;
  return EXIT_FAILURE;
}
int i,ii,Nobmax;
h2>>i;
h2>>i;
h2>>Nobmax;
h2>>i;
h2>>i>>ii;
std::vector<std::array<int,3>> oblabels; // ob labels N,l,2j
for(int N=0; N<=Nobmax; N++){
  for(int l=N%2; l<=N; l=l+2){
    for(int jj=std::abs(2*l-1); jj<=2*l+1; jj=jj+2){
      std::array<int,3> ob={N,l,jj};
      oblabels.push_back(ob);
    }
  }
}
while(true){
  int k1,k2,k3,k4,JJ,species;
  h2>>k1>>k2>>k3>>k4>>JJ>>species;
  if(h2){
    double me;
    h2>>me;
    int N1=oblabels[k1-1][0];
    int l1=oblabels[k1-1][1];
    int jj1=oblabels[k1-1][2];
    int N2=oblabels[k2-1][0];
    int l2=oblabels[k2-1][1];
    int jj2=oblabels[k2-1][2];
    int N3=oblabels[k3-1][0];
    int l3=oblabels[k3-1][1];
    int jj3=oblabels[k3-1][2];
    int N4=oblabels[k4-1][0];
    int l4=oblabels[k4-1][1];
    int jj4=oblabels[k4-1][2];
    double d=1.0;
    if(k1==k2)d*=sqrt(2.0);
    if(k3==k4)d*=sqrt(2.0);
    if(species==11){
      vp[{N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJ}]=d*me;
    }else if(species==22){
      vn[{N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJ}]=d*me;
    }else{
      vpn[{N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJ}]=me;
    }
  }else{
    break;
  }
}
h2.close();

//*****************************************************************
// Calculate direct part of potential kernel
//*****************************************************************
int N1vp = std::stoi(argv[2]);
int N1vn = std::stoi(argv[3]);
int etamaxp=Nmax+N1vp;
int etamaxn=Nmax+N1vn;
int etamax=std::max(etamaxp,etamaxn);
int Nclustermax=Nmax+1-Nmax%2;

int lm1pmy=0;
int mu1pmy=1;
int kappa1pmy=1;
int L1pmy=1;
int SS1pmy=0;
int Npmy=0;
int lmpmy=0;
int mupmy=1;
int kappapmy=1;
int Lpmy=1;
int SSpmy=1;
int lm1my=0;
int mu1my=1;
int kappa1my=1;
int L1my=1;
int SS1my=0;
int Nmy=0;
int lmmy=0;
int mumy=1;
int kappamy=1;
int Lmy=1;
int SSmy=1;
int JJmy=1;

bool write=false;

std::cout<<"Number of lm1p,mu1p,lm1,mu1: "<<obrho_by_labels.size()<<std::endl;

int rank,nproc;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&nproc);

std::cout<<"Hello from rank "<<rank<<std::endl;
   
su3::init(100);

//for(int ind=0; ind<obrho_by_labels.size(); ind++){
 std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>>::iterator it1=obrho_by_labels.begin();
 std::advance(it1,rank);
 int lm1p=it1->first[0];
 int mu1p=it1->first[1];
 int lm1=it1->first[2];
 int mu1=it1->first[3];
 std::cout<<"Rank "<<rank<<" computes lm1p mu1p lm1 mu1 = "<<lm1p<<" "<<mu1p<<" "<<lm1<<" "<<mu1<<std::endl;
 std::map<std::array<int,15>,std::map<std::array<int,4>,std::array<double,2>>> pot_dir;
 // Key is {SS1',N',lm',mu',kappa',L',SS',SS1,N,lm,mu,kappa,L,SS,JJ}
 // Second key is {kappa1',L1',kappa1,L1}
 // Value is {proton projectile, neutron projectile}
 for(int Np=0; Np<=Nclustermax; Np++){
  for(std::array<int,3> xp : KroneckerProduct(lm1p,mu1p,Np,0)){
   for(int Lp=0; Lp<=xp[0]+xp[1]; Lp++){
    int kappamaxp=su3::kmax(xp[0],xp[1],Lp);
    for(int kappap=1; kappap<=kappamaxp; kappap++){
     for(int N=0; N<=Nclustermax; N++){
      for(std::array<int,3> x : KroneckerProduct(lm1,mu1,N,0)){
       for(int L=0; L<=x[0]+x[1]; L++){
        int kappamax=su3::kmax(x[0],x[1],L);
        for(int kappa=1; kappa<=kappamax; kappa++){
         for(std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>::iterator
             it2=it1->second.begin(); it2!=it1->second.end(); it2++){
          int SS1p=it2->first[0];
          int SS1=it2->first[1];
          for(int SSp=std::abs(SS1p-1); SSp<=SS1p+1; SSp+=2){
           for(int JJ=std::abs(2*Lp-SSp); JJ<=2*Lp+SSp; JJ+=2){
            for(int SS=std::abs(SS1-1); SS<=SS1+1; SS+=2){
             if(std::abs(2*L-SS)>JJ || JJ>2*L+SS)continue;
	     std::map<std::array<int,4>,std::array<double,2>> mes;
	     for(std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>::iterator
                 it3=it2->second.begin(); it3!=it2->second.end(); it3++){
              mes[it3->first]={0.0,0.0};
             }
	     pot_dir[{SS1p,Np,xp[0],xp[1],kappap,Lp,SSp,SS1,N,x[0],x[1],kappa,L,SS,JJ}]=mes;
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
 for(int Na=0; Na<=etamax; Na++){
  int Nbmin=std::max(0,Na-Nmax);
  if((Na-Nbmin)%2!=0)Nbmin++; // so that the OBD operator doesn't change parity
  for(int Nb=Nbmin; Nb<=std::min(etamax,Na+Nmax); Nb+=2){
   for(std::array<int,3> x0 : KroneckerProduct(Na,0,0,Nb)){
    int rho0max=su3::mult(lm1,mu1,x0[0],x0[1],lm1p,mu1p);
    if(rho0max==0)continue;
    std::vector<double> su3cgs1,su3cgs2;
    su3cgs1.reserve(std::max(lm1,mu1)*std::max(x0[0],x0[1])*std::max(lm1p,mu1p));
    su3cgs2.reserve(Na*Nb*std::max(x0[0],x0[1]));
    for(int L0=0; L0<=x0[0]+x0[1]; L0++){
     int kappa0max=su3::kmax(x0[0],x0[1],L0);
     if(kappa0max==0)continue;
     for(int L1ppp=0; L1ppp<=lm1p+mu1p; L1ppp++){
      int kappa1pppmax=su3::kmax(lm1p,mu1p,L1ppp);
      if(kappa1pppmax==0)continue;
      for(int L1pp=0; L1pp<=lm1+mu1; L1pp++){
       int kappa1ppmax=su3::kmax(lm1,mu1,L1pp);
       if(kappa1ppmax==0 || std::abs(L1pp-L0)>L1ppp || L1ppp>L1pp+L0)continue;
       su3::wu3r3w(lm1,mu1,x0[0],x0[1],lm1p,mu1p,L1pp,L0,L1ppp,kappa1pppmax,kappa0max,kappa1ppmax,rho0max,su3cgs1);
       for(int Np=0; Np<=Nclustermax; Np++){
        for(std::array<int,3> xp : KroneckerProduct(lm1p,mu1p,Np,0)){
         std::vector<double> su3cgs3;
	 su3cgs3.reserve(std::max(lm1p,mu1p)*Np*std::max(xp[0],xp[1]));
	 for(int Lp=0; Lp<=xp[0]+xp[1]; Lp++){
          int kappapmax=su3::kmax(xp[0],xp[1],Lp);
	  if(kappapmax==0)continue;
          for(int lp=0; lp<=Np; lp++){
	   if(su3::kmax(Np,0,lp)==0 || std::abs(L1ppp-lp)>Lp || Lp>L1ppp+lp)continue;
	   su3::wu3r3w(lm1p,mu1p,Np,0,xp[0],xp[1],L1ppp,lp,Lp,kappapmax,1,kappa1pppmax,1,su3cgs3);
           for(int N=0; N<=Nclustermax; N++){
            for(std::array<int,3> x : KroneckerProduct(lm1,mu1,N,0)){
	     std::vector<double> su3cgs4;
             su3cgs4.reserve(std::max(lm1,mu1)*N*std::max(x[0],x[1]));
	     for(int L=0; L<=x[0]+x[1]; L++){
	      int kappamax=su3::kmax(x[0],x[1],L);
              if(kappamax==0)continue;
	      double factor1=-sqrt(double((2*L1ppp+1)*(2*L0+1)*(2*Lp+1)*(2*L+1)));
              if(Nb%2!=0)factor1=-factor1;
              for(int l=0; l<=N; l++){
               if(su3::kmax(N,0,l)==0 || std::abs(L1pp-l)>L || L>L1pp+l)continue;
	       su3::wu3r3w(lm1,mu1,N,0,x[0],x[1],L1pp,l,L,kappamax,1,kappa1ppmax,1,su3cgs4);
               for(int la=0; la<=Na; la++){
                if(su3::kmax(Na,0,la)==0)continue;
		for(int lb=0; lb<=Nb; lb++){
                 if(su3::kmax(0,Nb,lb)==0 || std::abs(la-lb)>L0 || L0>la+lb)continue;
		 su3::wu3r3w(Na,0,0,Nb,x0[0],x0[1],la,lb,L0,kappa0max,1,1,1,su3cgs2);
                 for(std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>::iterator
                     it2=it1->second.begin(); it2!=it1->second.end(); it2++){
                  int SS1p=it2->first[0];
                  int SS1=it2->first[1];
                  for(int SS0=0; SS0<=2; SS0+=2){
                   if(std::abs(SS1-SS0)>SS1p || SS1p>SS1+SS0)continue;
                   for(int JJ0=std::abs(2*L0-SS0); JJ0<=2*L0+SS0; JJ0+=2){
                    for(int II1ppp=std::abs(2*L1ppp-SS1p); II1ppp<=2*L1ppp+SS1p; II1ppp+=2){
                     for(int II1pp=std::abs(2*L1pp-SS1); II1pp<=2*L1pp+SS1; II1pp+=2){
                      if(std::abs(II1pp-JJ0)>II1ppp || II1ppp>II1pp+JJ0)continue;
                      double factor2=factor1*(II1ppp+1)*(II1pp+1)*(JJ0+1)*sqrt(double((SS0+1)*(SS1p+1)))
                                            *wig9jj(2*L1pp,SS1,II1pp,2*L0,SS0,JJ0,2*L1ppp,SS1p,II1ppp);
                      for(int SSp=std::abs(SS1p-1); SSp<=SS1p+1; SSp+=2){
                       for(int JJ=std::abs(2*Lp-SSp); JJ<=2*Lp+SSp; JJ+=2){
                        for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp+=2){
                         if(std::abs(II1ppp-jjp)>JJ || JJ>II1ppp+jjp)continue;
                         double factor3=factor2*sqrt(double((jjp+1)*(SSp+1)))
                                               *wig9jj(2*L1ppp,2*lp,2*Lp,SS1p,1,SSp,II1ppp,jjp,JJ);
                         for(int SS=std::abs(SS1-1); SS<=SS1+1; SS+=2){
                          if(std::abs(2*L-SS)>JJ || JJ>2*L+SS)continue;
                          for(int jj=std::abs(2*l-1); jj<=2*l+1; jj+=2){
                           if(std::abs(II1pp-jj)>JJ || JJ>II1pp+jj || std::abs(JJ0-jjp)>jj || jj>JJ0+jjp)continue;
                           double factor4=factor3*sqrt(double((jj+1)*(SS+1)))*wig6jj(II1pp,JJ0,II1ppp,jjp,JJ,jj)
                                                 *wig9jj(2*L1pp,2*l,2*L,SS1,1,SS,II1pp,jj,JJ);
                           for(int jja=std::abs(2*la-1); jja<=2*la+1; jja+=2){
                            for(int jjb=std::abs(2*lb-1); jjb<=2*lb+1; jjb+=2){
                             if(std::abs(jja-jjb)>JJ0 || JJ0>jja+jjb)continue;
                             double factor5=factor4*sqrt(double((jja+1)*(jjb+1)))*wig9jj(2*la,2*lb,2*L0,1,1,SS0,jja,jjb,JJ0);
		             for(int JJ0p=std::abs(jja-jjp); JJ0p<=jja+jjp; JJ0p+=2){
			      if(std::abs(jjb-jj)>JJ0p || JJ0p>jjb+jj)continue;
			      double factor6=factor5*(JJ0p+1)*wig6jj(jja,jjb,JJ0,jj,jjp,JJ0p);
			      if(((jjb+II1pp+JJ+JJ0p)/2)%2!=0)factor6=-factor6;
			      double vmep=like_int(Na,la,jja,Np,lp,jjp,Nb,lb,jjb,N,l,jj,JJ0p,vp);
			      double vmen=like_int(Na,la,jja,Np,lp,jjp,Nb,lb,jjb,N,l,jj,JJ0p,vn);
			      double vmepn=pn_int(Na,la,jja,Np,lp,jjp,Nb,lb,jjb,N,l,jj,JJ0p,vpn);
			      int ind=-1;
			      for(int ikappa1ppp=0; ikappa1ppp<kappa1pppmax; ikappa1ppp++){
			       for(int ikappa0=0; ikappa0<kappa0max; ikappa0++){
			        double factor7=factor6*su3cgs2[ikappa0];
				for(int ikappa1pp=0; ikappa1pp<kappa1ppmax; ikappa1pp++){
                                 for(int rho0=1; rho0<=rho0max; rho0++){
				  ind++;
                                  double factor8=factor7*su3cgs1[ind];
				  std::array<int,6> keyrho={Na,Nb,x0[0],x0[1],SS0/2,rho0};
			          for(int kappap=1; kappap<=kappapmax; kappap++){
			           double factor9=factor8*su3cgs3[ikappa1ppp+kappa1pppmax*(kappap-1)];
				   for(int kappa=1; kappa<=kappamax; kappa++){
				    double factor10=factor9*su3cgs4[ikappa1pp+kappa1ppmax*(kappa-1)];
                                    std::array<int,15> key1={SS1p,Np,xp[0],xp[1],kappap,Lp,SSp,SS1,N,x[0],x[1],kappa,L,SS,JJ};
                                    for(std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>::iterator
                                        it3=it2->second.begin(); it3!=it2->second.end(); it3++){
                                     std::array<int,4> key2=it3->first;
		     	             if(Np<=etamaxp && N<=etamaxp && Na<=etamaxp && Nb<=etamaxp){
                                      pot_dir[key1][key2][0]+=factor10*vmep*it3->second[keyrho][0];
                                      pot_dir[key1][key2][1]+=factor10*vmepn*it3->second[keyrho][0];
				     }
                                     if(Np<=etamaxn && N<=etamaxn && Na<=etamaxn && Nb<=etamaxn){
				      //pot_dir[key1][key2][0]+= in contruction
                                      pot_dir[key1][key2][1]+=factor10*vmen*it3->second[keyrho][1];
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
    }
   }
  }
 }
 for(std::map<std::array<int,15>,std::map<std::array<int,4>,std::array<double,2>>>::iterator it1=pot_dir.begin(); it1!=pot_dir.end(); it1++){
  int SS1p=it1->first[0];
  int Np=it1->first[1];
  int lmp=it1->first[2];
  int mup=it1->first[3];
  int kappap=it1->first[4];
  int Lp=it1->first[5];
  int SSp=it1->first[6];
  int SS1=it1->first[7];
  int N=it1->first[8];
  int lm=it1->first[9];
  int mu=it1->first[10];
  int kappa=it1->first[11];
  int L=it1->first[12];
  int SS=it1->first[13];
  int JJ=it1->first[14];
  for(std::map<std::array<int,4>,std::array<double,2>>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
   std::cout<<lm1p<<" "<<mu1p<<" "<<it2->first[0]<<" "<<it2->first[1]<<" "<<SS1p<<" "<<Np<<" "
	    <<lmp<<" "<<mup<<" "<<kappap<<" "<<Lp<<" "<<SSp<<" "
	    <<lm1<<" "<<mu1<<" "<<it2->first[2]<<" "<<it2->first[3]<<" "<<SS1<<" "<<N<<" "
	    <<lm<<" "<<mu<<" "<<kappa<<" "<<L<<" "<<SS<<" "<<JJ<<" "<<it2->second[0]<<" "<<it2->second[1]<<std::endl;
  }
 }
//}
su3::finalize();

MPI_Finalize();

return EXIT_SUCCESS;
}

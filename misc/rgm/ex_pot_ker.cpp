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

if(argc!=5){
  std::cout<<"Usage: "<<argv[0]<<" Nmax N1vp N1vn N0"<<std::endl;
  std::cout<<"where Nmax is Nmax for SpNCCI calculation of target"<<std::endl;
  return EXIT_FAILURE;
}

int Nmax = std::stoi(argv[1]);

//****************************************************************
// Read TB rhos
//****************************************************************
std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>> tbrho_by_labels;
// tbrho_by_labels[lm1',mu1',lm1,mu1][SS1',SS1][kappa1',L1',kappa1,L1][N1,N2,N3,N4,lmf,muf,Sf,lmi,mui,Si,lm0,mu0,S0,rho0,rho][proton,neutron,proton-neutron]

std::ifstream tbrhos("tbrho.dat");
if(!tbrhos){
  std::cout<<"Could not open file with TB rhos!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int N1,N2,N3,N4,lmf,muf,Sf,lmi,mui,Si,lm0,mu0,S0,rho0,rho,lmp,mup,kappap,Lp,SSp,lm,mu,kappa,L,SS;
  tbrhos>>N1>>N2>>N3>>N4>>lmf>>muf>>Sf>>lmi>>mui>>Si>>lm0>>mu0>>S0>>rho0>>rho>>lmp>>mup>>kappap>>Lp>>SSp>>lm>>mu>>kappa>>L>>SS;
  if(tbrhos){
    std::array<double,3> tbrho;
    tbrhos>>tbrho[0]>>tbrho[1]>>tbrho[2];
    std::array<int,4> key1={lmp,mup,lm,mu};
    std::array<int,2> key2={SSp,SS};
    std::array<int,4> key3={kappap,Lp,kappa,L};
    std::array<int,15> key4={N1,N2,N3,N4,lmf,muf,Sf,lmi,mui,Si,lm0,mu0,S0,rho0,rho};
    std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>>::iterator it1=tbrho_by_labels.find(key1);
    if(it1==tbrho_by_labels.end()){
      std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>> value1;
      std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>> value2;
      std::map<std::array<int,15>,std::array<double,3>> value3;
      value3[key4]=tbrho;
      value2[key3]=value3;
      value1[key2]=value2;
      tbrho_by_labels[key1]=value1;
    }else{
      std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>::iterator
	it2=tbrho_by_labels[key1].find(key2);
      if(it2==tbrho_by_labels[key1].end()){
        std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>> value2;
        std::map<std::array<int,15>,std::array<double,3>> value3;
        value3[key4]=tbrho;
        value2[key3]=value3;
	tbrho_by_labels[key1][key2]=value2;
      }else{
        std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>::iterator it3=tbrho_by_labels[key1][key2].find(key3);
	if(it3==tbrho_by_labels[key1][key2].end()){
          std::map<std::array<int,15>,std::array<double,3>> value3;
          value3[key4]=tbrho;
	  tbrho_by_labels[key1][key2][key3]=value3;
        }else{
          tbrho_by_labels[key1][key2][key3][key4]=tbrho;
        }
      }
    }
  }else{
    break;
  }
}
tbrhos.close();

//***************************************************************************
// Read interaction
//***************************************************************************
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

//*****************************************************************************
// Calculate exchange part of potential kernel
//*****************************************************************************
int N1vp = std::stoi(argv[2]);
int N1vn = std::stoi(argv[3]);
int N0 = std::stoi(argv[4]);
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
bool write=true;

std::vector<std::array<int,5>> x1px1Np;
for(std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>>::iterator it1=tbrho_by_labels.begin(); it1!=tbrho_by_labels.end(); it1++){
  int lm1p=it1->first[0];
  int mu1p=it1->first[1];
  int lm1=it1->first[2];
  int mu1=it1->first[3];
  for(int Np=0; Np<=std::min(Nclustermax,etamax); Np++){
    x1px1Np.push_back({lm1p,mu1p,lm1,mu1,Np});
  }
}
//std::cout<<"Number of outer loop iterations: "<<x1px1Np.size()<<std::endl;

int rank,nproc;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&nproc);

//std::cout<<"Hello from rank "<<rank<<std::endl;

su3::init(50);

//for(std::array<int,5> outer : x1px1Np){
 int lm1p=x1px1Np[rank][0];
 int mu1p=x1px1Np[rank][1];
 int lm1=x1px1Np[rank][2];
 int mu1=x1px1Np[rank][3];
 std::array<int,4> x1px1={lm1p,mu1p,lm1,mu1};
 int Np=x1px1Np[rank][4];
// std::cout<<"Rank "<<rank<<" computes lm1p mu1p lm1 mu1 Np = "<<lm1p<<" "<<mu1p<<" "<<lm1<<" "<<mu1<<" "<<Np<<std::endl;
 std::map<std::array<int,14>,std::map<std::array<int,4>,std::array<double,2>>> pot_ex;
 // First key is {SS1',lm',mu',kappa',L',SS',SS1,N,lm,mu,kappa,L,SS,JJ}
 // Second key is {kappa1',L1',kappa1,L1}
 // Value is {proton projectile, neutron projectile}
 for(int N=0; N<=Nclustermax; N++){
  for(std::array<int,3> x : KroneckerProduct(lm1,mu1,N,0)){
   for(int L=0; L<=x[0]+x[1]; L++){
    int kappamax=su3::kmax(x[0],x[1],L);
    if(kappamax==0)continue;
    for(int kappa=1; kappa<=kappamax; kappa++){
     for(std::array<int,3> xp : KroneckerProduct(lm1p,mu1p,Np,0)){
      for(int Lp=0; Lp<=xp[0]+xp[1]; Lp++){
       int kappapmax=su3::kmax(xp[0],xp[1],Lp);
       for(int kappap=1; kappap<=kappapmax; kappap++){
        for(std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>::iterator
	    it2=tbrho_by_labels[x1px1].begin(); it2!=tbrho_by_labels[x1px1].end(); it2++){
         int SS1p=it2->first[0];
         int SS1=it2->first[1];
         for(int SSp=std::abs(SS1p-1); SSp<=SS1p+1; SSp+=2){
          for(int JJ=std::abs(2*Lp-SSp); JJ<=2*Lp+SSp; JJ+=2){
           for(int SS=std::abs(SS1-1); SS<=SS1+1; SS+=2){
            if(std::abs(2*L-SS)>JJ || JJ>2*L+SS)continue;
	    std::map<std::array<int,4>,std::array<double,2>> me;
            for(std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>::iterator
                it3=it2->second.begin(); it3!=it2->second.end(); it3++){
	     me[it3->first]={0.0,0.0};
            }
	    pot_ex[{SS1p,xp[0],xp[1],kappap,Lp,SSp,SS1,N,x[0],x[1],kappa,L,SS,JJ}]=me;
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
 for(int Nd=0; Nd<=std::min(etamax,N0+Nmax-Np); Nd++){
  for(std::array<int,3> xd : KroneckerProduct(0,Np,0,Nd)){
   std::vector<double> su3cgs1;
   su3cgs1.reserve(Np*Nd*std::max(xd[0],xd[1]));
   for(int Na=0; Na<=etamax; Na++){
    int Nad=Na-Np-Nd;
    int Nbmin=std::max(0,-Nad-Nmax);
    if((Nad+Nbmin)%2!=0)Nbmin++; // so that the TBD operator doesn't change parity
    for(int Nb=Nbmin; Nb<=std::min(std::min(etamax,N0+Nmax-Na),Nmax-Nad); Nb+=2){
     for(std::array<int,3> xab : KroneckerProduct(Na,0,Nb,0)){
      std::vector<double> su3cgs2;
      su3cgs2.reserve(Na*Nb*std::max(xab[0],xab[1]));
      for(std::array<int,3> x0 : KroneckerProduct(xab[0],xab[1],xd[0],xd[1])){
       int rhomax=su3::mult(lm1,mu1,x0[0],x0[1],lm1p,mu1p);
       if(rhomax==0)continue;
       std::vector<double> su3cgs3,su3cgs4;
       su3cgs3.reserve(std::max(xab[0],xab[1])*std::max(xd[0],xd[1])*std::max(x0[0],x0[1]));
       su3cgs4.reserve(std::max(lm1,mu1)*std::max(x0[0],x0[1])*std::max(lm1p,mu1p));
       for(int Lab=0; Lab<=xab[0]+xab[1]; Lab++){
	int kappaabmax=su3::kmax(xab[0],xab[1],Lab);
	if(kappaabmax==0)continue;
        for(int Ld=0; Ld<=xd[0]+xd[1]; Ld++){
         int kappadmax=su3::kmax(xd[0],xd[1],Ld);
	 if(kappadmax==0)continue;
         for(int L0=0; L0<=x0[0]+x0[1]; L0++){
          int kappa0max=su3::kmax(x0[0],x0[1],L0);
          if(kappa0max==0 || std::abs(Lab-Ld)>L0 || L0>Lab+Ld)continue;
          su3::wu3r3w(xab[0],xab[1],xd[0],xd[1],x0[0],x0[1],Lab,Ld,L0,kappa0max,kappadmax,kappaabmax,x0[2],su3cgs3);
          for(int L1ppp=0; L1ppp<=lm1p+mu1p; L1ppp++){
           int kappa1pppmax=su3::kmax(lm1p,mu1p,L1ppp);
	   if(kappa1pppmax==0)continue;
           for(int L1pp=0; L1pp<=lm1+mu1; L1pp++){
            int kappa1ppmax=su3::kmax(lm1,mu1,L1pp);
            if(kappa1ppmax==0 || std::abs(L1pp-L0)>L1ppp || L1ppp>L1pp+L0)continue;
	    su3::wu3r3w(lm1,mu1,x0[0],x0[1],lm1p,mu1p,L1pp,L0,L1ppp,kappa1pppmax,kappa0max,kappa1ppmax,rhomax,su3cgs4);
            for(std::array<int,3> xp : KroneckerProduct(lm1p,mu1p,Np,0)){
	     std::vector<double> su3cgs5;
             su3cgs5.reserve(std::max(lm1p,mu1p)*Np*std::max(xp[0],xp[1]));
             for(int Lp=0; Lp<=xp[0]+xp[1]; Lp++){
              int kappapmax=su3::kmax(xp[0],xp[1],Lp);
	      if(kappapmax==0)continue;
              for(int lp=0; lp<=Np; lp++){
               if(su3::kmax(Np,0,lp)==0 || std::abs(L1ppp-lp)>Lp || Lp>L1ppp+lp)continue;
	       su3::wu3r3w(lm1p,mu1p,Np,0,xp[0],xp[1],L1ppp,lp,Lp,kappapmax,1,kappa1pppmax,1,su3cgs5);
               for(int N=0; N<=Nclustermax; N++){
                for(std::array<int,3> x : KroneckerProduct(lm1,mu1,N,0)){
	         std::vector<double> su3cgs6;
                 su3cgs6.reserve(std::max(lm1,mu1)*N*std::max(x[0],x[1]));
                 for(int L=0; L<=x[0]+x[1]; L++){
                  int kappamax=su3::kmax(x[0],x[1],L);
		  if(kappamax==0)continue;
		  double factor1=sqrt(double((2*Lab+1)*(2*Ld+1)*(2*L0+1)*(2*L1ppp+1)*(2*Lp+1)*(2*L+1)));
		  if((Np+Nd)%2!=0)factor1=-factor1;
                  for(int l=0; l<=N; l++){
                   if(su3::kmax(N,0,l)==0 || std::abs(L1pp-l)>L || L>L1pp+l)continue;
		   su3::wu3r3w(lm1,mu1,N,0,x[0],x[1],L1pp,l,L,kappamax,1,kappa1ppmax,1,su3cgs6);
                   for(int ld=0; ld<=Nd; ld++){
                    if(su3::kmax(0,Nd,ld)==0 || std::abs(lp-ld)>Ld || Ld>lp+ld)continue;
		    su3::wu3r3w(0,Np,0,Nd,xd[0],xd[1],lp,ld,Ld,kappadmax,1,1,1,su3cgs1);
                    for(int la=0; la<=Na; la++){
                     if(su3::kmax(Na,0,la)==0)continue;
                     for(int lb=0; lb<=Nb; lb++){
                      if(su3::kmax(Nb,0,lb)==0 || std::abs(la-lb)>Lab || Lab>la+lb)continue;
		      su3::wu3r3w(Na,0,Nb,0,xab[0],xab[1],la,lb,Lab,kappaabmax,1,1,1,su3cgs2);
                      for(int SSab=0; SSab<=2; SSab+=2){
                       for(int SSd=0; SSd<=2; SSd+=2){
                        for(int SS0=std::abs(SSab-SSd); SS0<=SSab+SSd; SS0+=2){
                         for(int JJab=std::abs(2*Lab-SSab); JJab<=2*Lab+SSab; JJab+=2){
                          for(int JJd=std::abs(2*Ld-SSd); JJd<=2*Ld+SSd; JJd+=2){
                           for(int JJ0=std::abs(2*L0-SS0); JJ0<=2*L0+SS0; JJ0+=2){
                            if(std::abs(JJab-JJd)>JJ0 || JJ0>JJab+JJd)continue;
                            double factor2=factor1*(JJab+1)*(JJd+1)*(JJ0+1)*sqrt(double((SSab+1)*(SSd+1)*(SS0+1)))
                                           *wig9jj(2*Lab,2*Ld,2*L0,SSab,SSd,SS0,JJab,JJd,JJ0);
                            if(((JJab+JJd+JJ0)/2)%2!=0)factor2=-factor2;
                            for(std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>::iterator it2=tbrho_by_labels[x1px1].begin(); it2!=tbrho_by_labels[x1px1].end(); it2++){
                             int SS1p=it2->first[0];
                             int SS1=it2->first[1];
                             for(int II1ppp=std::abs(2*L1ppp-SS1p); II1ppp<=2*L1ppp+SS1p; II1ppp+=2){
                              for(int II1pp=std::abs(2*L1pp-SS1); II1pp<=2*L1pp+SS1; II1pp+=2){
                               if(std::abs(II1pp-JJ0)>II1ppp || II1ppp>II1pp+JJ0)continue;
                               double factor3=factor2*(II1ppp+1)*(II1pp+1)*sqrt(double(SS1p+1))
                                              *wig9jj(2*L1pp,SS1,II1pp,2*L0,SS0,JJ0,2*L1ppp,SS1p,II1ppp);
                               for(int SSp=std::abs(SS1p-1); SSp<=SS1p+1; SSp+=2){
                                for(int JJ=std::abs(2*Lp-SSp); JJ<=2*Lp+SSp; JJ+=2){
                                 for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp+=2){
                                  if(std::abs(II1ppp-jjp)>JJ || JJ>II1ppp+jjp)continue;
                                  double factor4=factor3*(jjp+1)*sqrt(double(SSp+1))*wig9jj(2*L1ppp,2*lp,2*Lp,SS1p,1,SSp,II1ppp,jjp,JJ);
                                  for(int SS=std::abs(SS1-1); SS<=SS1+1; SS+=2){
                                   if(std::abs(2*L-SS)>JJ || JJ>2*L+SS)continue;
                                   for(int jj=std::abs(2*l-1); jj<=2*l+1; jj+=2){
                                    if(std::abs(II1pp-jj)>JJ || JJ>II1pp+jj || std::abs(JJ0-jjp)>jj || jj>JJ0+jjp)continue;
                                    double factor5=factor4*sqrt(double((SS+1)*(jj+1)))*wig6jj(II1pp,JJ0,II1ppp,jjp,JJ,jj)
                                                   *wig9jj(2*L1pp,2*l,2*L,SS1,1,SS,II1pp,jj,JJ);
                                    if(((jj+II1pp+JJ)/2)%2!=0)factor5=-factor5;
                                    for(int jjd=std::abs(2*ld-1); jjd<=2*ld+1; jjd+=2){
                                     if(std::abs(jjp-jjd)>JJd || JJd>jjp+jjd || std::abs(jjd-JJab)>jj || jj>jjd+JJab)continue;
                                     double factor6=factor5*sqrt(double(jjd+1))*wig6jj(jjd,JJd,jjp,JJ0,jj,JJab)
                                                    *wig9jj(2*lp,2*ld,2*Ld,1,1,SSd,jjp,jjd,JJd);
                                     for(int jja=std::abs(2*la-1); jja<=2*la+1; jja+=2){
                                      for(int jjb=std::abs(2*lb-1); jjb<=2*lb+1; jjb+=2){
                                       if(std::abs(jja-jjb)>JJab || JJab>jja+jjb)continue;
                                       double factor7=factor6*sqrt(double((jja+1)*(jjb+1)))*wig9jj(2*la,2*lb,2*Lab,1,1,SSab,jja,jjb,JJab);
                                       double vmep=like_int(Na,la,jja,Nb,lb,jjb,Nd,ld,jjd,N,l,jj,JJab,vp);
                                       double vmen=like_int(Na,la,jja,Nb,lb,jjb,Nd,ld,jjd,N,l,jj,JJab,vn);
                                       double vmepn=pn_int(Na,la,jja,Nb,lb,jjb,Nd,ld,jjd,N,l,jj,JJab,vpn);
				       int ind3=-1;
				       for(int ikappa0=0; ikappa0<kappa0max; ikappa0++){
				        for(int ikappad=0; ikappad<kappadmax; ikappad++){
				         double factor8=factor7*su3cgs1[ikappad];
				         for(int ikappaab=0; ikappaab<kappaabmax; ikappaab++){
				          double factor9=factor8*su3cgs2[ikappaab];
	                                  for(int rho0=1; rho0<=x0[2]; rho0++){
				           ind3++;
				           double factor10=factor9*su3cgs3[ind3];
					   int ind5=-1;
					   for(int kappap=1; kappap<=kappapmax; kappap++){
				            for(int ikappa1ppp=0; ikappa1ppp<kappa1pppmax; ikappa1ppp++){
					     ind5++;
					     double factor11=factor10*su3cgs5[ind5];
					     int ind6=-1;
                                             for(int kappa=1; kappa<=kappamax; kappa++){
					      std::array<int,14> key1={SS1p,xp[0],xp[1],kappap,Lp,SSp,SS1,N,x[0],x[1],kappa,L,SS,JJ};
					      for(int ikappa1pp=0; ikappa1pp<kappa1ppmax; ikappa1pp++){
					       ind6++;
					       double factor12=factor11*su3cgs6[ind6];
					       for(int rho=1; rho<=rhomax; rho++){
					        double factor13=factor12
			                *su3cgs4[rho-1+rhomax*ikappa1pp+rhomax*kappa1ppmax*ikappa0+rhomax*kappa1ppmax*kappa0max*ikappa1ppp];
		                  std::array<int,15> keyrho={Na,Nb,Np,Nd,xab[0],xab[1],SSab/2,xd[0],xd[1],SSd/2,x0[0],x0[1],SS0/2,rho0,rho};
                                                for(std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>::iterator
                                                    it3=it2->second.begin(); it3!=it2->second.end(); it3++){
                                                 std::array<int,4> key2=it3->first;
                                                 if(Na<=etamaxp && Nb<=etamaxp && Np<=etamaxp && Nd<=etamaxp)
						   pot_ex[key1][key2][0]+=0.5*factor13*vmep*it3->second[keyrho][0];
						 if(Na<=etamaxn && Nb<=etamaxn && Np<=etamaxn && Nd<=etamaxn)
					           pot_ex[key1][key2][1]+=0.5*factor13*vmen*it3->second[keyrho][1];
						 //if(Na<=etamaxn && Nb<=etamaxp && Np<=etamaxp && Nd<=etamaxn)
						 // pot_ex[key1][key2][0]+= in construction
						 if(Na<=etamaxp && Nb<=etamaxn && Np<=etamaxn && Nd<=etamaxp)
					           pot_ex[key1][key2][1]+=factor13*vmepn*it3->second[keyrho][2]; 
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
 for(std::map<std::array<int,14>,std::map<std::array<int,4>,std::array<double,2>>>::iterator it1=pot_ex.begin(); it1!=pot_ex.end(); it1++){
  int SS1p=it1->first[0];
  int lmp=it1->first[1];
  int mup=it1->first[2];
  int kappap=it1->first[3];
  int Lp=it1->first[4];
  int SSp=it1->first[5];
  int SS1=it1->first[6];
  int N=it1->first[7];
  int lm=it1->first[8];
  int mu=it1->first[9];
  int kappa=it1->first[10];
  int L=it1->first[11];
  int SS=it1->first[12];
  int JJ=it1->first[13];
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

#include <cstdio>
#include <fstream>
#include <istream>
#include <iostream>
#include <omp.h>

#include "spncci/computation_control.h"
#include "wigxjpf.h"

double like_int(const int N1,const int l1,const int jj1,const int N2,const int l2,const int jj2,const int N3,const int l3,const int jj3,
	    const int N4,const int l4,const int jj4,const int J,std::map<std::array<int,13>,double>& v){
  std::array<int,13> key={N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,J};
  std::map<std::array<int,13>,double>::iterator it=v.find(key);
  if(it!=v.end())return v[key];

  key={N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,J};
  it=v.find(key);
  if(it!=v.end()){
    int phase=1+J-(jj1+jj2)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }

  key={N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,J};
  it=v.find(key);
  if(it!=v.end()){
    int phase=1+J-(jj3+jj4)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }

  key={N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,J};
  it=v.find(key);
  if(it!=v.end()){
    int phase=(jj1+jj2+jj3+jj4)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }
  
  key={N3,l3,jj3,N4,l4,jj4,N1,l1,jj1,N2,l2,jj2,J};
  it=v.find(key);
  if(it!=v.end())return v[key];
  
  key={N4,l4,jj4,N3,l3,jj3,N1,l1,jj1,N2,l2,jj2,J};
  it=v.find(key);
  if(it!=v.end()){
    int phase=1+J-(jj3+jj4)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }
  
  key={N3,l3,jj3,N4,l4,jj4,N2,l2,jj2,N1,l1,jj1,J};
  it=v.find(key);
  if(it!=v.end()){
    int phase=1+J-(jj1+jj2)/2;
    if(phase%2==0){
      return v[key];
    }else{
      return -v[key];
    }
  }
  
  key={N4,l4,jj4,N3,l3,jj3,N2,l2,jj2,N1,l1,jj1,J};
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
	      const int N4,const int l4,const int jj4,const int J,std::map<std::array<int,13>,double>& vpn){
  std::array<int,13> key={N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,J};
  std::map<std::array<int,13>,double>::iterator it=vpn.find(key);
  if(it!=vpn.end())return vpn[key];
 
  key={N3,l3,jj3,N4,l4,jj4,N1,l1,jj1,N2,l2,jj2,J};
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
std::cout<<"Nmax = "<<Nmax<<" (this is Nmax for SpNCCI calculation of target)"<<std::endl;

/*
double wcoef=u3::W(u3::SU3(4,3),1,2,u3::SU3(2,0),1,2,u3::SU3(5,2),1,2,1);
std::cout<<"Wigner coef: "<<wcoef<<std::endl;
double ucoef=u3::UCached(u_coef_cache,u3::SU3(6,4),u3::SU3(2,0),u3::SU3(5,5),u3::SU3(2,0),u3::SU3(6,5),1,1,u3::SU3(4,0),1,1);
std::cout<<"U coef: "<<ucoef<<std::endl;
double ninej=wig9jj(1,2,3,
                    4,6,8,
                    3,6,9);
std::cout<<"9j coef: "<<ninej<<std::endl;
*/

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

std::map<std::array<int,13>,double> vp,vn,vpn;
// Labels are N1 l1 2j1 N2 l2 2j2 N3 l3 2j3 N4 l4 2j4 J

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
      vp[{N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJ/2}]=d*me;
    }else if(species==22){
      vn[{N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJ/2}]=d*me;
    }else{
      vpn[{N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJ/2}]=me;
    }
  }else{
    break;
  }
}
h2.close();

int N1vp = std::stoi(argv[2]);
int N1vn = std::stoi(argv[3]);
int N0 = std::stoi(argv[4]);
int etamaxp=Nmax+N1vp;
int etamaxn=Nmax+N1vn;
int etamax=std::max(etamaxp,etamaxn);
int Nclustermax=Nmax+1-Nmax%2;
if(false){//
u3::UCoefCache u_coef_cache;
spncci::InitializeSpNCCI();

//*********************************************************
// Norm kernel
//*********************************************************
/*
std::map<std::array<int,12>,std::map<std::array<int,3>,std::array<double,2>>> PkA;
// First key is {lm1',mu1',kappa1',L1',2S1',N',lm1,mu1,kappa1,L1,2S1,N}
// Second key is {lm,mu,2S}
// Values are proton and neutron

for(std::map<std::array<int,16>,std::array<double,2>>::iterator it=obrho_by_labels.begin(); it!=obrho_by_labels.end(); it++){
  int lm1p=it->first[6];
  int mu1p=it->first[7];
  int kappa1p=it->first[8];
  int L1p=it->first[9];
  int SS1p=it->first[10];
  int Np=it->first[1];
  int lm1=it->first[11];
  int mu1=it->first[12];
  int kappa1=it->first[13];
  int L1=it->first[14];
  int SS1=it->first[15];
  int N=it->first[0];
  int lm0=it->first[2];
  int mu0=it->first[3];
  int S0=it->first[4];
  int rho=it->first[5];
  double factor=sqrt(double((SS1p+1)*(2*S0+1)*u3::dim(u3::SU3(lm0,mu0)))/double(u3::dim(u3::SU3(N,0))));
  std::array<double,2> factorpn;
  factorpn[0]=factor*it->second[0];
  factorpn[1]=factor*it->second[1];
  std::array<int,12> key={lm1p,mu1p,kappa1p,L1p,SS1p,Np,lm1,mu1,kappa1,L1,SS1,N};
  std::map<std::array<int,12>,std::map<std::array<int,3>,std::array<double,2>>>::iterator iter=PkA.find(key);
  if(iter==PkA.end()){
    std::map<std::array<int,3>,std::array<double,2>> map;
    for(auto x : u3::KroneckerProduct(u3::SU3(lm1,mu1),u3::SU3(N,0))){
      if(u3::OuterMultiplicity(u3::SU3(lm1p,mu1p),u3::SU3(Np,0),x.irrep)==0)continue;
      int lm=x.irrep.lambda();
      int mu=x.irrep.mu();
      double u=u3::UCached(u_coef_cache,u3::SU3(lm1,mu1),u3::SU3(lm0,mu0),x.irrep,u3::SU3(Np,0),u3::SU3(lm1p,mu1p),rho,1,u3::SU3(N,0),1,1);
      for(int SS=std::abs(SS1-1); SS<=SS1+1; SS=SS+2){
        if(SS<std::abs(SS1p-1) || SS>SS1p+1)continue;
	std::array<int,3> key2={lm,mu,SS};
	std::array<double,2> array;
	int phase=Np-N+(SS1+1+SS)/2+lm0+mu0;
	int sign;
	if(phase%2==0){
          sign=1;
	}else{
          sign=-1;
	}
	for(int i=0; i<=1; i++){
          array[i]=sign*factorpn[i]*wig6jj(SS1,2*S0,SS1p,1,SS,1)*u;
	}
	map[key2]=array;
      }
    }
    PkA[key]=map;
  }else{
    for(auto x : u3::KroneckerProduct(u3::SU3(lm1,mu1),u3::SU3(N,0))){
      if(u3::OuterMultiplicity(u3::SU3(lm1p,mu1p),u3::SU3(Np,0),x.irrep)==0)continue;
      int lm=x.irrep.lambda();
      int mu=x.irrep.mu();
      double u=u3::UCached(u_coef_cache,u3::SU3(lm1,mu1),u3::SU3(lm0,mu0),x.irrep,u3::SU3(Np,0),u3::SU3(lm1p,mu1p),rho,1,u3::SU3(N,0),1,1);
      for(int SS=std::abs(SS1-1); SS<=SS1+1; SS=SS+2){
        if(SS<std::abs(SS1p-1) || SS>SS1p+1)continue;
	std::array<int,3> key2={lm,mu,SS};
	int phase=Np-N+(SS1+1+SS)/2+lm0+mu0;
	int sign;
	if(phase%2==0){
          sign=1;
	}else{
          sign=-1;
	}
	for(int i=0; i<=1; i++){
          PkA[key][key2][i]+=sign*factorpn[i]*wig6jj(SS1,2*S0,SS1p,1,SS,1)*u;
	}
      }
    }
  }
}

std::ofstream PkAfile("PiA.dat");
if(!PkAfile){
  std::cout<<"Could not open PkA file"<<std::endl;
  return EXIT_FAILURE;
}
PkAfile<<"lm1' mu1' kappa1' L1' 2S1' N' lm1 mu1 kappa1 L1 2S1 N lm mu 2S proton neutron"<<std::endl;
for(std::map<std::array<int,12>,std::map<std::array<int,3>,std::array<double,2>>>::iterator it1=PkA.begin(); it1!=PkA.end(); it1++){
  int lm1p=it1->first[0];
  int mu1p=it1->first[1];
  int kappa1p=it1->first[2];
  int L1p=it1->first[3];
  int SS1p=it1->first[4];
  int Np=it1->first[5];
  int lm1=it1->first[6];
  int mu1=it1->first[7];
  int kappa1=it1->first[8];
  int L1=it1->first[9];
  int SS1=it1->first[10];
  int N=it1->first[11];
  for(std::map<std::array<int,3>,std::array<double,2>>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
    PkAfile<<lm1p<<" "<<mu1p<<" "<<kappa1p<<" "<<L1p<<" "<<SS1p<<" "<<Np<<" "<<lm1<<" "<<mu1<<" "<<kappa1<<" "<<L1<<" "<<SS1<<" "<<N<<" "
           <<it2->first[0]<<" "<<it2->first[1]<<" "<<it2->first[2]<<" "<<it2->second[0]<<" "<<it2->second[1]<<std::endl;
  }
}
PkAfile.close();
*/
std::map<std::array<int,15>,std::array<double,2>> norm_ex;
// Key is {lm1',mu1',kappa1',L1',SS1',N',lm1,mu1,kappa1,L1,SS1,N,lm,mu,SS}
// Value is {proton, neutron}

for(std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>>::iterator it1=obrho_by_labels.begin(); it1!=obrho_by_labels.end(); it1++){
 int lm1p=it1->first[0];
 int mu1p=it1->first[1];
 int lm1=it1->first[2];
 int mu1=it1->first[3];
 u3::SU3 x1p(lm1p,mu1p);
 u3::SU3 x1(lm1,mu1);
 for(int Np=0; Np<=Nclustermax; Np++){
  for(int N=0; N<=Nclustermax; N++){
   for(auto x : u3::KroneckerProduct(x1,u3::SU3(N,0))){
    if(u3::OuterMultiplicity(x1p,u3::SU3(Np,0),x.irrep)==0)continue;
    for(auto x0 : u3::KroneckerProduct(u3::SU3(N,0),u3::SU3(0,Np))){
     for(int rho=1; rho<=u3::OuterMultiplicity(x1,x0.irrep,x1p); rho++){
      double factor1=sqrt(double(2*u3::dim(x0.irrep))/double((N+1)*(N+2)))
	             *u3::UCached(u_coef_cache,x1,x0.irrep,x.irrep,u3::SU3(Np,0),x1p,rho,1,u3::SU3(N,0),1,1);
      if((Np-N+x0.irrep.lambda()+x0.irrep.mu())%2!=0)factor1=-factor1;
      for(std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>::iterator
          it2=it1->second.begin(); it2!=it1->second.end(); it2++){
       int SS1p=it2->first[0];
       int SS1=it2->first[1];
       for(int SS=std::abs(SS1-1); SS<=SS1+1; SS=SS+2){
        if(std::abs(SS1p-1)>SS || SS>SS1p+1)continue;
	for(int S0=0; S0<=1; S0++){
         if(std::abs(SS1-2*S0)>SS1p || SS1p>SS1+2*S0)continue;
	 double factor2=factor1*sqrt(double((SS1p+1)*(2*S0+1)))*wig6jj(SS1,2*S0,SS1p,1,SS,1);
	 if(((SS1+1+SS)/2)%2!=0)factor2=-factor2;
	 std::array<int,6> key4={N,Np,x0.irrep.lambda(),x0.irrep.mu(),S0,rho};
         for(std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>::iterator
             it3=it2->second.begin(); it3!=it2->second.end(); it3++){
          int kappa1p=it3->first[0];
	  int L1p=it3->first[1];
	  int kappa1=it3->first[2];
          int L1=it3->first[3];
          std::array<int,15> key={lm1p,mu1p,kappa1p,L1p,SS1p,Np,lm1,mu1,kappa1,L1,SS1,N,x.irrep.lambda(),x.irrep.mu(),SS};
          std::map<std::array<int,15>,std::array<double,2>>::iterator iter=norm_ex.find(key);
          if(iter==norm_ex.end()){
           std::array<double,2> mes;
           mes[0]=factor2*it3->second[key4][0];
	   mes[1]=factor2*it3->second[key4][1];
	   norm_ex[key]=mes;
          }else{
           norm_ex[key][0]+=factor2*it3->second[key4][0];
           norm_ex[key][1]+=factor2*it3->second[key4][1];
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

u_coef_cache.clear();

std::ofstream PkAfile("PiA.dat");
if(!PkAfile){
  std::cout<<"Could not open PkA file"<<std::endl;
  return EXIT_FAILURE;
}
//PkAfile<<"lm1' mu1' kappa1' L1' 2S1' N' lm1 mu1 kappa1 L1 2S1 N lm mu 2S proton neutron"<<std::endl;
for(std::map<std::array<int,15>,std::array<double,2>>::iterator it=norm_ex.begin(); it!=norm_ex.end(); it++){
 PkAfile<<it->first[0]<<" "<<it->first[1]<<" "<<it->first[2]<<" "<<it->first[3]<<" "<<it->first[4]<<" "<<it->first[5]<<" "
	<<it->first[6]<<" "<<it->first[7]<<" "<<it->first[8]<<" "<<it->first[9]<<" "<<it->first[10]<<" "<<it->first[11]<<" "
	<<it->first[12]<<" "<<it->first[13]<<" "<<it->first[14]<<" "<<it->second[0]<<" "<<it->second[1]<<std::endl;
}
PkAfile.close();
}//
//*********************************************************
// Direct part of potential kernel
//*********************************************************

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

std::map<std::array<int,19>,std::map<std::array<int,4>,std::array<double,2>>> pot_dir;
// First key is {lm1',mu1',SS1',N',lm',mu',kappa',L',SS',lm1,mu1,SS1,N,lm,mu,kappa,L,SS,JJ}
// Second key is {kappa1',L1',kappa1,L1}
// Value is {proton projectile, neutron projectile}

int number=0;
for(std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>>::iterator it1=obrho_by_labels.begin(); it1!=obrho_by_labels.end(); it1++){
 int lm1p=it1->first[0];
 int mu1p=it1->first[1];
 int lm1=it1->first[2];
 int mu1=it1->first[3];
 u3::SU3 x1p(lm1p,mu1p);
 u3::SU3 x1(lm1,mu1);
 for(int Np=0; Np<=Nclustermax; Np++){
  for(auto xp : u3::KroneckerProduct(x1p,u3::SU3(Np,0))){
   for(auto Lptagged : u3::BranchingSO3(xp.irrep)){
    int Lp=int(Lptagged.irrep);
    for(int kappap=1; kappap<=Lptagged.tag; kappap++){
     for(int N=0; N<=Nclustermax; N++){
      for(auto x : u3::KroneckerProduct(x1,u3::SU3(N,0))){
       for(auto Ltagged : u3::BranchingSO3(x.irrep)){
        int L=int(Ltagged.irrep);
        for(int kappa=1; kappa<=Ltagged.tag; kappa++){
         for(std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>::iterator
             it2=it1->second.begin(); it2!=it1->second.end(); it2++){
          int SS1p=it2->first[0];
          int SS1=it2->first[1];
          for(int SSp=std::abs(SS1p-1); SSp<=SS1p+1; SSp=SSp+2){
           for(int JJ=std::abs(2*Lp-SSp); JJ<=2*Lp+SSp; JJ=JJ+2){
            for(int SS=std::abs(SS1-1); SS<=SS1+1; SS=SS+2){
             if(std::abs(2*L-SS)>JJ || JJ>2*L+SS)continue;
             std::array<int,19> key1={lm1p,mu1p,SS1p,Np,xp.irrep.lambda(),xp.irrep.mu(),kappap,Lp,SSp,
                                      lm1,mu1,SS1,N,x.irrep.lambda(),x.irrep.mu(),kappa,L,SS,JJ};
	     std::map<std::array<int,4>,std::array<double,2>> mes;
	     for(std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>::iterator
                 it3=it2->second.begin(); it3!=it2->second.end(); it3++){
              std::array<int,4> key2=it3->first;
              number++;
              mes[key2]={0.0,0.0};
             }
	     pot_dir[key1]=mes;
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
std::cout<<"Number of MEs between SA-RGM basis states: "<<number<<std::endl;
std::cout<<"Number of outer loop iterations: "<<obrho_by_labels.size()<<std::endl;

// By default, variables are shared except loop iteration counters, which are private
// Private variables are lm1p,mu1p,lm1,mu1,x1p,x1,Nbmin,L0,L1pp,L1ppp,factor1,Lp,lp,factor2,L,l,factor3,la,lb,factor4,SS1p,SS1,keyrho,factor5,factor6,factor7,factor8,factor9,vmep,vmen,vmepn,key1,key2
#pragma omp parallel
{
//int num_threads=omp_get_num_threads();
spncci::InitializeSpNCCI();
//std::cout<<"Number of threads: "<<num_threads<<std::endl;
#pragma omp for
for(int ind=0; ind<obrho_by_labels.size(); ind++){
 std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>>::iterator it1=obrho_by_labels.begin();
 std::advance(it1,ind);
 int lm1p=it1->first[0];
 int mu1p=it1->first[1];
 int lm1=it1->first[2];
 int mu1=it1->first[3];
 u3::SU3 x1p(lm1p,mu1p);
 u3::SU3 x1(lm1,mu1);
/*bool write1;
if(lm1p==lm1pmy && mu1p==mu1pmy && lm1==lm1my && mu1==mu1my){
write1=true;
}else{
write1=false;
}*/
 for(int Na=0; Na<=etamax; Na++){
  int Nbmin=std::max(0,Na-Nmax);
  if((Na-Nbmin)%2!=0)Nbmin++; // so that the OBD operator doesn't change parity
  for(int Nb=Nbmin; Nb<=std::min(etamax,Na+Nmax); Nb=Nb+2){
   for(auto x0 : u3::KroneckerProduct(u3::SU3(Na,0),u3::SU3(0,Nb))){
    for(int rho0=1; rho0<=u3::OuterMultiplicity(x1,x0.irrep,x1p); rho0++){
     for(auto L0tagged : u3::BranchingSO3(x0.irrep)){
      int L0=int(L0tagged.irrep);
      for(auto L1pptagged : u3::BranchingSO3(x1)){
       int L1pp=int(L1pptagged.irrep);
       for(auto L1ppptagged : u3::BranchingSO3(x1p)){
	int L1ppp=int(L1ppptagged.irrep);
        if(std::abs(L1pp-L0)>L1ppp || L1ppp>L1pp+L0)continue;
        for(int kappa0=1; kappa0<=L0tagged.tag; kappa0++){
         for(int kappa1ppp=1; kappa1ppp<=L1ppptagged.tag; kappa1ppp++){
          for(int kappa1pp=1; kappa1pp<=L1pptagged.tag; kappa1pp++){
           double factor1=-sqrt(double((2*L1ppp+1)*(2*L0+1)))*u3::W(x1,kappa1pp,L1pp,x0.irrep,kappa0,L0,x1p,kappa1ppp,L1ppp,rho0);
	   if(Nb%2!=0)factor1=-factor1;
/*if(write && write1){
std::cout<<"Na Nb lm0 mu0 rho0 L0 L1pp L1ppp kappa0 kappa1ppp kappa1pp: "<<Na<<" "<<Nb<<" "<<x0.irrep.lambda()<<" "<<x0.irrep.mu()<<" "<<rho0<<" "<<L0<<" "<<L1pp<<" "<<L1ppp<<" "<<kappa0<<" "<<kappa1ppp<<" "<<kappa1pp<<std::endl;
std::cout<<"factor1=-(-1)^Nb*sqrt((2*L1ppp+1)*(2*L0+1))*<x1 kappa1pp L1pp; x0 kappa0 L0||x1p kappa1ppp L1ppp>_rho0="<<factor1<<std::endl;
}*/
           for(int Np=0; Np<=Nclustermax; Np++){
	    for(auto xp : u3::KroneckerProduct(x1p,u3::SU3(Np,0))){
	     for(auto Lptagged : u3::BranchingSO3(xp.irrep)){
	      int Lp=int(Lptagged.irrep);
	      for(int kappap=1; kappap<=Lptagged.tag; kappap++){
               for(auto lptagged : u3::BranchingSO3(u3::SU3(Np,0))){
		int lp=int(lptagged.irrep);
                if(std::abs(L1ppp-lp)>Lp || Lp>L1ppp+lp)continue;
                double factor2=factor1*sqrt(double(2*Lp+1))*u3::W(x1p,kappa1ppp,L1ppp,u3::SU3(Np,0),1,lp,xp.irrep,kappap,Lp,1);
/*bool write2;
if(Np==Npmy && xp.irrep.lambda()==lmpmy && xp.irrep.mu()==mupmy && Lp==Lpmy && kappap==kappapmy){
write2=true;
}else{
write2=false;
}
if(write && write1 && write2){
std::cout<<"lp: "<<lp<<std::endl;
std::cout<<"factor2=factor1*sqrt(2*Lp+1)*<x1p kappa1ppp L1ppp; (Np,0)lp||xp kappap Lp>="<<factor2<<std::endl;
}*/
		for(int N=0; N<=Nclustermax; N++){
                 for(auto x : u3::KroneckerProduct(x1,u3::SU3(N,0))){
                  for(auto Ltagged : u3::BranchingSO3(x.irrep)){
	           int L=int(Ltagged.irrep);
                   for(int kappa=1; kappa<=Ltagged.tag; kappa++){
                    for(auto ltagged : u3::BranchingSO3(u3::SU3(N,0))){
		     int l=int(ltagged.irrep);
                     if(std::abs(L1pp-l)>L || L>L1pp+l)continue;
                     double factor3=factor2*sqrt(double(2*L+1))*u3::W(x1,kappa1pp,L1pp,u3::SU3(N,0),1,l,x.irrep,kappa,L,1);
/*bool write3;
if(N==Nmy && x.irrep.lambda()==lmmy && x.irrep.mu()==mumy && L==Lmy && kappa==kappamy){
write3=true;
}else{
write3=false;
}
if(write && write1 && write2 && write3){
std::cout<<"l: "<<l<<std::endl;
std::cout<<"factor3=factor2*sqrt(2*L+1)*<x1 kappa1pp L1pp; (N,0)l||x kappa L>="<<factor3<<std::endl;
}*/
		     for(auto latagged : u3::BranchingSO3(u3::SU3(Na,0))){
		      int la=int(latagged.irrep);
                      for(auto lbtagged : u3::BranchingSO3(u3::SU3(0,Nb))){
		       int lb=int(lbtagged.irrep);
                       if(std::abs(la-lb)>L0 || L0>la+lb)continue;
		       double factor4=factor3*u3::W(u3::SU3(Na,0),1,la,u3::SU3(0,Nb),1,lb,x0.irrep,kappa0,L0,1);
/*if(write && write1 && write2 && write3){
std::cout<<"la lb: "<<la<<" "<<lb<<std::endl;
std::cout<<"factor4=factor3*<(Na,0)la; (0,Nb)lb||x0 kappa0 L0>="<<factor4<<std::endl;
}*/
                       for(std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>>::iterator
                           it2=it1->second.begin(); it2!=it1->second.end(); it2++){
                        int SS1p=it2->first[0];
                        int SS1=it2->first[1];
		        for(int S0=0; S0<=1; S0++){
                         if(std::abs(SS1-2*S0)>SS1p || SS1p>SS1+2*S0)continue;
			 std::array<int,6> keyrho={Na,Nb,x0.irrep.lambda(),x0.irrep.mu(),S0,rho0};
			 for(int J0=std::abs(L0-S0); J0<=L0+S0; J0++){
                          for(int II1ppp=std::abs(2*L1ppp-SS1p); II1ppp<=2*L1ppp+SS1p; II1ppp=II1ppp+2){
                           for(int II1pp=std::abs(2*L1pp-SS1); II1pp<=2*L1pp+SS1; II1pp=II1pp+2){
                            if(std::abs(II1pp-2*J0)>II1ppp || II1ppp>II1pp+2*J0)continue;
                            double factor5=factor4*(II1ppp+1)*(II1pp+1)*(2*J0+1)*sqrt(double((2*S0+1)*(SS1p+1)))
				           *wig9jj(2*L1pp,SS1,II1pp,2*L0,2*S0,2*J0,2*L1ppp,SS1p,II1ppp);
/*bool write4;
if(SS1p==SS1pmy && SS1==SS1my){
write4=true;
}else{
write4=false;
}
if(write && write1 && write2 && write3 && write4){
std::cout<<"S0 J0 II1ppp II1pp: "<<S0<<" "<<J0<<" "<<II1ppp<<" "<<II1pp<<std::endl;
std::cout<<"factor5=factor4*(II1ppp+1)*(II1pp+1)*(2*J0+1)*sqrt((2*S0+1)*(SS1p+1))*wig9j(L1pp,S1,I1pp,L0,S0,J0,L1ppp,S1p,I1ppp)="<<factor5<<std::endl;
}*/
			    for(int SSp=std::abs(SS1p-1); SSp<=SS1p+1; SSp=SSp+2){
			     for(int JJ=std::abs(2*Lp-SSp); JJ<=2*Lp+SSp; JJ=JJ+2){
			      for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp=jjp+2){
                               if(std::abs(II1ppp-jjp)>JJ || JJ>II1ppp+jjp)continue;
			       double factor6=factor5*sqrt(double((jjp+1)*(SSp+1)))
				              *wig9jj(2*L1ppp,2*lp,2*Lp,SS1p,1,SSp,II1ppp,jjp,JJ);
/*bool write5;
if(SSp==SSpmy && JJ==JJmy){
write5=true;
}else{
write5=false;
}
if(write && write1 && write2 && write3 && write4 && write5){
std::cout<<"jjp: "<<jjp<<std::endl;
std::cout<<"factor6=factor5*sqrt((jjp+1)*(SSp+1))*wig9j(L1ppp,lp,Lp,S1p,1/2,Sp,I1ppp,jp,J)="<<factor6<<std::endl;
}*/
			       for(int SS=std::abs(SS1-1); SS<=SS1+1; SS=SS+2){
                                if(std::abs(2*L-SS)>JJ || JJ>2*L+SS)continue;
				for(int jj=std::abs(2*l-1); jj<=2*l+1; jj=jj+2){
                                 if(std::abs(II1pp-jj)>JJ || JJ>II1pp+jj || std::abs(2*J0-jjp)>jj || jj>2*J0+jjp)continue;
				 double factor7=factor6*sqrt(double((jj+1)*(SS+1)))*wig6jj(II1pp,2*J0,II1ppp,jjp,JJ,jj)
					        *wig9jj(2*L1pp,2*l,2*L,SS1,1,SS,II1pp,jj,JJ);
/*bool write6;
if(SS==SSmy){
write6=true;
}else{
write6=false;
}
if(write && write1 && write2 && write3 && write4 && write5 && write6){
std::cout<<"jj: "<<jj<<std::endl;
std::cout<<"factor7=factor6*sqrt((jj+1)*(SS+1))*wig6j(I1pp,J0,I1ppp,jp,J,j)*wig9j(L1pp,l,L,S1,1/2,S,I1pp,j,J)="<<factor7<<std::endl;
}*/
                                 for(int jja=std::abs(2*la-1); jja<=2*la+1; jja=jja+2){
				  for(int jjb=std::abs(2*lb-1); jjb<=2*lb+1; jjb=jjb+2){
                                   if(std::abs(jja-jjb)>2*J0 || 2*J0>jja+jjb)continue;
				   double factor8=factor7*sqrt(double((jja+1)*(jjb+1)))*wig9jj(2*la,2*lb,2*L0,1,1,2*S0,jja,jjb,2*J0);
/*if(write && write1 && write2 && write3 && write4 && write5 && write6){
std::cout<<"jja jjb: "<<jja<<" "<<jjb<<std::endl;
std::cout<<"factor8=factor7*sqrt((jja+1)*(jjb+1))*wig9j(la,lb,L0,1/2,1/2,S0,ja,jb,J0)="<<factor8<<std::endl;
}*/
				   for(int J0p=std::abs(jja-jjp)/2; J0p<=(jja+jjp)/2; J0p++){
			            if(std::abs(jjb-jj)>2*J0p || 2*J0p>jjb+jj)continue;
				    double factor9=factor8*(2*J0p+1)*wig6jj(jja,jjb,2*J0,jj,jjp,2*J0p);
				    if(((jjb+II1pp+JJ)/2+J0p)%2!=0)factor9=-factor9;
/*if(write && write1 && write2 && write3 && write4 && write5 && write6){
std::cout<<"J0p: "<<J0p<<std::endl;
std::cout<<"factor9=factor8*(-1)^(jb+J0p+I1pp+J)*(2*J0p+1)*wig6j(ja,jb,J0,j,jp,J0p)="<<factor9<<std::endl;
}*/
				    double vmep=like_int(Na,la,jja,Np,lp,jjp,Nb,lb,jjb,N,l,jj,J0p,vp);
				    double vmen=like_int(Na,la,jja,Np,lp,jjp,Nb,lb,jjb,N,l,jj,J0p,vn);
				    double vmepn=pn_int(Na,la,jja,Np,lp,jjp,Nb,lb,jjb,N,l,jj,J0p,vpn);
                                    std::array<int,19> key1={lm1p,mu1p,SS1p,Np,xp.irrep.lambda(),xp.irrep.mu(),kappap,Lp,SSp,
					                     lm1,mu1,SS1,N,x.irrep.lambda(),x.irrep.mu(),kappa,L,SS,JJ};
                                    for(std::map<std::array<int,4>,std::map<std::array<int,6>,std::array<double,2>>>::iterator
                                        it3=it2->second.begin(); it3!=it2->second.end(); it3++){
                                     std::array<int,4> key2=it3->first;
		     	             if(Np<=etamaxp && N<=etamaxp && Na<=etamaxp && Nb<=etamaxp){
                                      pot_dir[key1][key2][0]+=factor9*vmep*it3->second[keyrho][0];
                                      pot_dir[key1][key2][1]+=factor9*vmepn*it3->second[keyrho][0];
				     }
                                     if(Np<=etamaxn && N<=etamaxn && Na<=etamaxn && Nb<=etamaxn){
				      //pot_dir[key1][key2][0]+= in contruction
                                      pot_dir[key1][key2][1]+=factor9*vmen*it3->second[keyrho][1];
				     }
/*bool write7;
if(it3->first[0]==kappa1pmy && it3->first[1]==L1pmy && it3->first[2]==kappa1my && it3->first[3]==L1my){
write7=true;
}else{
write7=false;
}
if(write && write1 && write2 && write3 && write4 && write5 && write6 && write7)std::cout<<"pot_dir[key1][key2][1]=pot_dir[key1][key2][1]+factor9*vmepn*rho+factor9*vmen*rho=pot_dir[key1][key2][1]+"<<factor9<<"*"<<vmepn<<"*"<<it3->second[keyrho][0]<<"+"<<factor9<<"*"<<vmen<<"*"<<it3->second[keyrho][1]<<"="<<pot_dir[key1][key2][1]<<std::endl;*/
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
} // End of OpenMP parallel region
std::ofstream ViAPiAfile("ViAPiA.dat");
if(!ViAPiAfile){
  std::cout<<"Could not open ViAPiA file"<<std::endl;
  return EXIT_FAILURE;
} 
//ViAPiAfile<<"lm1' mu1' kappa1' L1' SS1' N' lm' mu' kappa' L' SS' lm1 mu1 kappa1 L1 SS1 N lm mu kappa L SS JJ prot_proj neut_proj"<<std::endl;
for(std::map<std::array<int,19>,std::map<std::array<int,4>,std::array<double,2>>>::iterator it1=pot_dir.begin(); it1!=pot_dir.end(); it1++){
  int lm1p=it1->first[0];
  int mu1p=it1->first[1];
  int SS1p=it1->first[2];
  int Np=it1->first[3];
  int lmp=it1->first[4];
  int mup=it1->first[5];
  int kappap=it1->first[6];
  int Lp=it1->first[7];
  int SSp=it1->first[8];
  int lm1=it1->first[9];
  int mu1=it1->first[10];
  int SS1=it1->first[11];
  int N=it1->first[12];
  int lm=it1->first[13];
  int mu=it1->first[14];
  int kappa=it1->first[15];
  int L=it1->first[16];
  int SS=it1->first[17];
  int JJ=it1->first[18];
  for(std::map<std::array<int,4>,std::array<double,2>>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
    ViAPiAfile<<lm1p<<" "<<mu1p<<" "<<it2->first[0]<<" "<<it2->first[1]<<" "<<SS1p<<" "<<Np<<" "
	      <<lmp<<" "<<mup<<" "<<kappap<<" "<<Lp<<" "<<SSp<<" "
	      <<lm1<<" "<<mu1<<" "<<it2->first[2]<<" "<<it2->first[3]<<" "<<SS1<<" "<<N<<" "
	      <<lm<<" "<<mu<<" "<<kappa<<" "<<L<<" "<<SS<<" "<<JJ<<" "<<it2->second[0]<<" "<<it2->second[1]<<std::endl;
  }
}
ViAPiAfile.close();
return EXIT_SUCCESS;//
//**************************************************
// Exchange part of potential kernel
//**************************************************
/*
std::map<std::array<int,19>,std::map<std::array<int,4>,std::array<double,2>>> pot_ex;
// First key is {lm1',mu1',SS1',N',lm',mu',kappa',L',SS',lm1,mu1,SS1,N,lm,mu,kappa,L,SS,JJ}
// Second key is {kappa1',L1',kappa1,L1}
// Value is {proton projectile, neutron projectile}
for(int N=0; N<=std::min(Nclustermax,etamax); N++){
 for(int Na=0; Na<=std::min(etamax,N0+Nmax-N); Na++){
  for(int Nc=0; Nc<=etamax; Nc++){
   int Nac=N+Na-Nc;
   int Ndmin=std::max(0,Nac-Nmax);
   if((Nac-Ndmin)%2!=0)Ndmin++; // so that the TBD operator doesn't change parity
   for(int Nd=Ndmin; Nd<=std::min(std::min(etamax,N0+Nmax-Nc),Nac+Nmax); Nd=Nd+2){
    for(auto xa : u3::KroneckerProduct(u3::SU3(Na,0),u3::SU3(N,0))){
     for(auto xcd : u3::KroneckerProduct(u3::SU3(0,Nc),u3::SU3(0,Nd))){
      for(auto x0 : u3::KroneckerProduct(xa.irrep,xcd.irrep)){
       for(auto Latagged : u3::BranchingSO3(xa.irrep)){
	int La=int(Latagged.irrep);
        for(auto Lcdtagged : u3::BranchingSO3(xcd.irrep)){
         int Lcd=int(Lcdtagged.irrep);
         for(auto L0tagged : u3::BranchingSO3(x0.irrep)){
          int L0=int(L0tagged.irrep);
          if(std::abs(La-Lcd)>L0 || L0>La+Lcd)continue;
	  for(int rho0=1; rho0<=x0.tag; rho0++){
           for(int kappaa=1; kappaa<=Latagged.tag; kappaa++){
            for(int kappacd=1; kappacd<=Lcdtagged.tag; kappacd++){
             for(int kappa0=1; kappa0<=L0tagged.tag; kappa0++){
              double factor1=sqrt(double((2*La+1)*(2*Lcd+1)*(2*L0+1)))
		             *u3::W(xa.irrep,kappaa,La,xcd.irrep,kappacd,Lcd,x0.irrep,kappa0,L0,rho0);
	      if((Nc+Nd)%2!=0)factor1=-factor1;
              for(std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>>::iterator it1=tbrho_by_labels.begin(); it1!=tbrho_by_labels.end(); it1++){
               int lm1p=it1->first[0]; 
               int mu1p=it1->first[1];
               int lm1=it1->first[2];
               int mu1=it1->first[3];
               u3::SU3 x1p(lm1p,mu1p);
               u3::SU3 x1(lm1,mu1);
	       for(int rho=1; rho<=u3::OuterMultiplicity(x1,x0.irrep,x1p); rho++){
                for(auto L1ppptagged : u3::BranchingSO3(x1p)){
                 int L1ppp=int(L1ppptagged.irrep);
		 for(auto L1pptagged : u3::BranchingSO3(x1)){
                  int L1pp=int(L1pptagged.irrep);
		  if(std::abs(L1pp-L0)>L1ppp || L1ppp>L1pp+L0)continue;
		  for(int kappa1ppp=1; kappa1ppp<=L1ppptagged.tag; kappa1ppp++){
                   for(int kappa1pp=1; kappa1pp<=L1pptagged.tag; kappa1pp++){
                    double factor2=factor1*sqrt(double(2*L1ppp+1))*u3::W(x1,kappa1pp,L1pp,x0.irrep,kappa0,L0,x1p,kappa1ppp,L1ppp,rho);
		    for(int Np=0; Np<=Nclustermax; Np++){
		     for(auto xp : u3::KroneckerProduct(x1p,u3::SU3(Np,0))){
		      for(auto Lptagged : u3::BranchingSO3(xp.irrep)){
		       int Lp=int(Lptagged.irrep);
		       for(int kappap=1; kappap<=Lptagged.tag; kappap++){
		        for(auto lptagged : u3::BranchingSO3(u3::SU3(Np,0))){
		         int lp=int(lptagged.irrep);
			 if(std::abs(L1ppp-lp)>Lp || Lp>L1ppp+lp)continue;
			 double factor3=factor2*sqrt(double(2*Lp+1))*u3::W(x1p,kappa1ppp,L1ppp,u3::SU3(Np,0),1,lp,xp.irrep,kappap,Lp,1);
			 for(auto x : u3::KroneckerProduct(x1,u3::SU3(N,0))){
                          for(auto Ltagged : u3::BranchingSO3(x.irrep)){
                           int L=int(Ltagged.irrep);
                           for(int kappa=1; kappa<=Ltagged.tag; kappa++){
                            for(auto ltagged : u3::BranchingSO3(u3::SU3(N,0))){
                             int l=int(ltagged.irrep);
                             if(std::abs(L1pp-l)>L || L>L1pp+l)continue;
                             double factor4=factor3*sqrt(double(2*L+1))*u3::W(x1,kappa1pp,L1pp,u3::SU3(N,0),1,l,x.irrep,kappa,L,1);
			     for(auto latagged : u3::BranchingSO3(u3::SU3(Na,0))){
			      int la=int(latagged.irrep);
			      if(std::abs(la-l)>La || La>la+l)continue;
			      double factor5=factor4*u3::W(u3::SU3(Na,0),1,la,u3::SU3(N,0),1,l,xa.irrep,kappaa,La,1);
                              for(auto lctagged : u3::BranchingSO3(u3::SU3(0,Nc))){
                               int lc=int(lctagged.irrep);
			       for(auto ldtagged : u3::BranchingSO3(u3::SU3(0,Nd))){
                                int ld=int(ldtagged.irrep);
				if(std::abs(lc-ld)>Lcd || Lcd>lc+ld)continue;
				double factor6=factor5*u3::W(u3::SU3(0,Nc),1,lc,u3::SU3(0,Nd),1,ld,xcd.irrep,kappacd,Lcd,1);
				for(int Sa=0; Sa<=1; Sa++){
                                 for(int Scd=0; Scd<=1; Scd++){
                                  for(int S0=std::abs(Sa-Scd); S0<=Sa+Scd; S0++){
			           std::array<int,15> keyrho={Na,N,Nc,Nd,xa.irrep.lambda(),xa.irrep.mu(),Sa,
					              xcd.irrep.lambda(),xcd.irrep.mu(),Scd,x0.irrep.lambda(),x0.irrep.mu(),S0,rho0,rho};
                                   for(int Ja=std::abs(La-Sa); Ja<=La+Sa; Ja++){
			            for(int Jcd=std::abs(Lcd-Scd); Jcd<=Lcd+Scd; Jcd++){
				     for(int J0=std::abs(L0-S0); J0<=L0+S0; J0++){
                                      if(std::abs(Ja-Jcd)>J0 || J0>Ja+Jcd)continue;
				      double factor7=factor6*(2*Ja+1)*(2*Jcd+1)*(2*J0+1)*sqrt(double((2*Sa+1)*(2*Scd+1)*(2*S0+1)))
					             *wig9jj(2*La,2*Lcd,2*L0,2*Sa,2*Scd,2*S0,2*Ja,2*Jcd,2*J0);
				      if((Ja+Jcd+J0)%2!=0)factor7=-factor7;
                                      for(std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
                                       int SS1p=it2->first[0];
                                       int SS1=it2->first[1];
				       for(int II1ppp=std::abs(2*L1ppp-SS1p); II1ppp<=2*L1ppp+SS1p; II1ppp=II1ppp+2){
                                        for(int II1pp=std::abs(2*L1pp-SS1); II1pp<=2*L1pp+SS1; II1pp=II1pp+2){
					 if(std::abs(II1pp-2*J0)>II1ppp || II1ppp>II1pp+2*J0)continue;
					 double factor8=factor7*(II1ppp+1)*(II1pp+1)*sqrt(double(SS1p+1))
						        *wig9jj(2*L1pp,SS1,II1pp,2*L0,2*S0,2*J0,2*L1ppp,SS1p,II1ppp);
                                         for(int SSp=std::abs(SS1p-1); SSp<=SS1p+1; SSp=SSp+2){
                                          for(int JJ=std::abs(2*Lp-SSp); JJ<=2*Lp+SSp; JJ=JJ+2){
				           for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp=jjp+2){
				            if(std::abs(II1ppp-jjp)>JJ || JJ>II1ppp+jjp)continue;
					    double factor9=factor8*sqrt(double((SSp+1)*(jjp+1)))
						           *wig9jj(2*L1ppp,2*lp,2*Lp,SS1p,1,SSp,II1ppp,jjp,JJ);
					    for(int SS=std::abs(SS1-1); SS<=SS1+1; SS=SS+2){
					     if(std::abs(2*L-SS)>JJ || JJ>2*L+SS)continue;
					     std::array<int,19> key1={lm1p,mu1p,SS1p,Np,xp.irrep.lambda(),xp.irrep.mu(),kappap,Lp,SSp,
                                                                      lm1,mu1,SS1,N,x.irrep.lambda(),x.irrep.mu(),kappa,L,SS,JJ};
                                             for(int jj=std::abs(2*l-1); jj<=2*l+1; jj=jj+2){
					      if(std::abs(II1pp-jj)>JJ || JJ>II1pp+jj || std::abs(2*J0-jjp)>jj || jj>2*J0+jjp)continue;
					      double factor10=factor9*(jj+1)*sqrt(double(SS+1))*wig6jj(II1pp,2*J0,II1ppp,jjp,JJ,jj)
						              *wig9jj(2*L1pp,2*l,2*L,SS1,1,SS,II1pp,jj,JJ);
                                              for(int jja=std::abs(2*la-1); jja<=2*la+1; jja=jja+2){
					       if(std::abs(jja-jj)>2*Ja || 2*Ja>jja+jj || std::abs(jja-2*Jcd)>jj || jj>jja+2*Jcd)continue;
					       double factor11=factor10*sqrt(double(jja+1))*wig6jj(jja,2*Jcd,jjp,2*J0,jj,2*Ja)
						               *wig9jj(2*la,2*l,2*La,1,1,2*Sa,jja,jj,2*Ja);
                                               for(int jjc=std::abs(2*lc-1); jjc<=2*lc+1; jjc=jjc+2){
                                                for(int jjd=std::abs(2*ld-1); jjd<=2*ld+1; jjd=jjd+2){
                                                 if(std::abs(jjc-jjd)>2*Jcd || 2*Jcd>jjc+jjd)continue;
                                                 double factor12=factor11*sqrt(double((jjc+1)*(jjd+1)))
							         *wig9jj(2*lc,2*ld,2*Lcd,1,1,2*Scd,jjc,jjd,2*Jcd);
						 if(((jj-jja-jjc-jjd+II1pp+JJ-jjp)/2)%2!=0)factor12=-factor12;
				                 double vmep=like_int(Na,la,jja,Np,lp,jjp,Nd,ld,jjd,Nc,lc,jjc,Jcd,vp);
				                 double vmen=like_int(Na,la,jja,Np,lp,jjp,Nd,ld,jjd,Nc,lc,jjc,Jcd,vn);
				                 double vmepn=pn_int(Na,la,jja,Np,lp,jjp,Nd,ld,jjd,Nc,lc,jjc,Jcd,vpn);
				                 std::map<std::array<int,19>,std::map<std::array<int,4>,std::array<double,2>>>::iterator
				                   iter=pot_ex.find(key1);
				                 if(iter==pot_ex.end()){
                                                  std::map<std::array<int,4>,std::array<double,2>> mes;
                                                  for(std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>::iterator
                                                      it3=it2->second.begin(); it3!=it2->second.end(); it3++){
                                                   std::array<int,4> key2=it3->first;
				                   std::array<double,2> me;
				                   if(Na<=etamaxp && N<=etamaxp && Nc<=etamaxp && Nd<=etamaxp){
                                                    me[0]=0.5*factor12*vmep*it3->second[keyrho][0];
				                   }else{
                                                    me[0]=0.0;
				                   }
						   if(Na<=etamaxn && N<=etamaxn && Nc<=etamaxn && Nd<=etamaxn){
                                                    me[1]=0.5*factor12*vmen*it3->second[keyrho][1];
                                                   }else{
                                                    me[1]=0.0;
                                                   }
				                   //if(Na<=etamaxn && N<=etamaxp && Nc<=etamaxp && Nd<=etamaxn)me[0]+= in contruction
                                                   if(Na<=etamaxp && N<=etamaxn && Nc<=etamaxn && Nd<=etamaxp)
					             me[1]+=factor12*vmepn*it3->second[keyrho][2];
					             //{std::array<int,15> keyrhopn={N,Na,Nd,Nc,xa.irrep.lambda(),xa.irrep.mu(),Sa,//
                                                      //xcd.irrep.lambda(),xcd.irrep.mu(),Scd,x0.irrep.lambda(),x0.irrep.mu(),S0,rho0,rho};//
						      //double hustota=it3->second[keyrhopn][2];//
						      //if((N+Na+Nd+Nc+xa.irrep.lambda()+xa.irrep.mu()+Sa//
                                                      //    +xcd.irrep.lambda()+xcd.irrep.mu()+Scd)%2!=0)hustota=-hustota;//
						      //me[1]+=factor12*vmepn*hustota;}//
				                   mes[key2]=me;
				                  }
				                  pot_ex[key1]=mes;
				                 }else{
                                                  for(std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>::iterator
                                                      it3=it2->second.begin(); it3!=it2->second.end(); it3++){
                                                   std::array<int,4> key2=it3->first;
                                                   if(Na<=etamaxp && N<=etamaxp && Nc<=etamaxp && Nd<=etamaxp)
						     pot_ex[key1][key2][0]+=0.5*factor12*vmep*it3->second[keyrho][0];
						   if(Na<=etamaxn && N<=etamaxn && Nc<=etamaxn && Nd<=etamaxn)
					             pot_ex[key1][key2][1]+=0.5*factor12*vmen*it3->second[keyrho][1];
						   //if(Na<=etamaxn && N<=etamaxp && Nc<=etamaxp && Nd<=etamaxn)
						   //  pot_ex[key1][key2][0]+= in construction
						   if(Na<=etamaxp && N<=etamaxn && Nc<=etamaxn && Nd<=etamaxp)
					             pot_ex[key1][key2][1]+=factor12*vmepn*it3->second[keyrho][2];
					             //{std::array<int,15> keyrhopn={N,Na,Nd,Nc,xa.irrep.lambda(),xa.irrep.mu(),Sa,//
                                                      //xcd.irrep.lambda(),xcd.irrep.mu(),Scd,x0.irrep.lambda(),x0.irrep.mu(),S0,rho0,rho};//
                                                      //double hustota=it3->second[keyrhopn][2];//
                                                      //if((N+Na+Nd+Nc+xa.irrep.lambda()+xa.irrep.mu()+Sa//
                                                      //    +xcd.irrep.lambda()+xcd.irrep.mu()+Scd)%2!=0)hustota=-hustota;//
                                                      //pot_ex[key1][key2][1]+=factor12*vmepn*hustota;}//
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
  }
 }
}

std::ofstream ViAPkAfile("ViAPkA.dat");
if(!ViAPkAfile){
  std::cout<<"Could not open ViAPkA file"<<std::endl;
  return EXIT_FAILURE;
} 
//ViAPkAfile<<"lm1' mu1' kappa1' L1' SS1' N' lm' mu' kappa' L' SS' lm1 mu1 kappa1 L1 SS1 N lm mu kappa L SS JJ prot_proj neut_proj"<<std::endl;
for(std::map<std::array<int,19>,std::map<std::array<int,4>,std::array<double,2>>>::iterator it1=pot_ex.begin(); it1!=pot_ex.end(); it1++){
  int lm1p=it1->first[0];
  int mu1p=it1->first[1];
  int SS1p=it1->first[2];
  int Np=it1->first[3];
  int lmp=it1->first[4];
  int mup=it1->first[5];
  int kappap=it1->first[6];
  int Lp=it1->first[7];
  int SSp=it1->first[8];
  int lm1=it1->first[9];
  int mu1=it1->first[10];
  int SS1=it1->first[11];
  int N=it1->first[12];
  int lm=it1->first[13];
  int mu=it1->first[14];
  int kappa=it1->first[15];
  int L=it1->first[16];
  int SS=it1->first[17];
  int JJ=it1->first[18];
  for(std::map<std::array<int,4>,std::array<double,2>>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
    ViAPkAfile<<lm1p<<" "<<mu1p<<" "<<it2->first[0]<<" "<<it2->first[1]<<" "<<SS1p<<" "<<Np<<" "
	      <<lmp<<" "<<mup<<" "<<kappap<<" "<<Lp<<" "<<SSp<<" "
	      <<lm1<<" "<<mu1<<" "<<it2->first[2]<<" "<<it2->first[3]<<" "<<SS1<<" "<<N<<" "
	      <<lm<<" "<<mu<<" "<<kappa<<" "<<L<<" "<<SS<<" "<<JJ<<" "<<it2->second[0]<<" "<<it2->second[1]<<std::endl;
  }
}
ViAPkAfile.close();
*/
//*********************************************************
// 2nd term of Hermitized exchange part of potential kernel
//*********************************************************

std::map<std::array<int,19>,std::map<std::array<int,4>,std::array<double,2>>> pot_exh;
// First key is {lm1',mu1',SS1',N',lm',mu',kappa',L',SS',lm1,mu1,SS1,N,lm,mu,kappa,L,SS,JJ}
// Second key is {kappa1',L1',kappa1,L1}
// Value is {proton projectile, neutron projectile}

for(int Np=0; Np<=std::min(Nclustermax,etamax); Np++){
 for(int Nd=0; Nd<=std::min(etamax,N0+Nmax-Np); Nd++){
  for(int Na=0; Na<=etamax; Na++){
   int Nad=Na-Np-Nd;
   int Nbmin=std::max(0,-Nad-Nmax);
   if((Nad+Nbmin)%2!=0)Nbmin++; // so that the TBD operator doesn't change parity
   for(int Nb=Nbmin; Nb<=std::min(std::min(etamax,N0+Nmax-Na),Nmax-Nad); Nb=Nb+2){
    for(auto xab : u3::KroneckerProduct(u3::SU3(Na,0),u3::SU3(Nb,0))){
     for(auto xd : u3::KroneckerProduct(u3::SU3(0,Np),u3::SU3(0,Nd))){
      for(auto x0 : u3::KroneckerProduct(xab.irrep,xd.irrep)){
       for(auto Labtagged : u3::BranchingSO3(xab.irrep)){
	int Lab=int(Labtagged.irrep);
        for(auto Ldtagged : u3::BranchingSO3(xd.irrep)){
         int Ld=int(Ldtagged.irrep);
         for(auto L0tagged : u3::BranchingSO3(x0.irrep)){
          int L0=int(L0tagged.irrep);
          if(std::abs(Lab-Ld)>L0 || L0>Lab+Ld)continue;
	  for(int rho0=1; rho0<=x0.tag; rho0++){
           for(int kappaab=1; kappaab<=Labtagged.tag; kappaab++){
            for(int kappad=1; kappad<=Ldtagged.tag; kappad++){
             for(int kappa0=1; kappa0<=L0tagged.tag; kappa0++){
              double factor1=sqrt(double((2*Lab+1)*(2*Ld+1)*(2*L0+1)))
		             *u3::W(xab.irrep,kappaab,Lab,xd.irrep,kappad,Ld,x0.irrep,kappa0,L0,rho0);
	      if((Np+Nd)%2!=0)factor1=-factor1;
              for(std::map<std::array<int,4>,std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>>::iterator it1=tbrho_by_labels.begin(); it1!=tbrho_by_labels.end(); it1++){
               int lm1p=it1->first[0]; 
               int mu1p=it1->first[1];
               int lm1=it1->first[2];
               int mu1=it1->first[3];
               u3::SU3 x1p(lm1p,mu1p);
               u3::SU3 x1(lm1,mu1);
	       for(int rho=1; rho<=u3::OuterMultiplicity(x1,x0.irrep,x1p); rho++){
                for(auto L1ppptagged : u3::BranchingSO3(x1p)){
                 int L1ppp=int(L1ppptagged.irrep);
		 for(auto L1pptagged : u3::BranchingSO3(x1)){
                  int L1pp=int(L1pptagged.irrep);
		  if(std::abs(L1pp-L0)>L1ppp || L1ppp>L1pp+L0)continue;
		  for(int kappa1ppp=1; kappa1ppp<=L1ppptagged.tag; kappa1ppp++){
                   for(int kappa1pp=1; kappa1pp<=L1pptagged.tag; kappa1pp++){
                    double factor2=factor1*sqrt(double(2*L1ppp+1))*u3::W(x1,kappa1pp,L1pp,x0.irrep,kappa0,L0,x1p,kappa1ppp,L1ppp,rho);
		    for(int N=0; N<=Nclustermax; N++){
		     for(auto x : u3::KroneckerProduct(x1,u3::SU3(N,0))){
		      for(auto Ltagged : u3::BranchingSO3(x.irrep)){
		       int L=int(Ltagged.irrep);
		       for(int kappa=1; kappa<=Ltagged.tag; kappa++){
		        for(auto ltagged : u3::BranchingSO3(u3::SU3(N,0))){
		         int l=int(ltagged.irrep);
			 if(std::abs(L1pp-l)>L || L>L1pp+l)continue;
			 double factor3=factor2*sqrt(double(2*L+1))*u3::W(x1,kappa1pp,L1pp,u3::SU3(N,0),1,l,x.irrep,kappa,L,1);
			 for(auto xp : u3::KroneckerProduct(x1p,u3::SU3(Np,0))){
                          for(auto Lptagged : u3::BranchingSO3(xp.irrep)){
                           int Lp=int(Lptagged.irrep);
                           for(int kappap=1; kappap<=Lptagged.tag; kappap++){
                            for(auto lptagged : u3::BranchingSO3(u3::SU3(Np,0))){
                             int lp=int(lptagged.irrep);
                             if(std::abs(L1ppp-lp)>Lp || Lp>L1ppp+lp)continue;
                             double factor4=factor3*sqrt(double(2*Lp+1))*u3::W(x1p,kappa1ppp,L1ppp,u3::SU3(Np,0),1,lp,xp.irrep,kappap,Lp,1);
			     for(auto ldtagged : u3::BranchingSO3(u3::SU3(0,Nd))){
			      int ld=int(ldtagged.irrep);
			      if(std::abs(lp-ld)>Ld || Ld>lp+ld)continue;
			      double factor5=factor4*u3::W(u3::SU3(0,Np),1,lp,u3::SU3(0,Nd),1,ld,xd.irrep,kappad,Ld,1);
                              for(auto latagged : u3::BranchingSO3(u3::SU3(Na,0))){
                               int la=int(latagged.irrep);
			       for(auto lbtagged : u3::BranchingSO3(u3::SU3(Nb,0))){
                                int lb=int(lbtagged.irrep);
				if(std::abs(la-lb)>Lab || Lab>la+lb)continue;
				double factor6=factor5*u3::W(u3::SU3(Na,0),1,la,u3::SU3(Nb,0),1,lb,xab.irrep,kappaab,Lab,1);
				for(int Sab=0; Sab<=1; Sab++){
                                 for(int Sd=0; Sd<=1; Sd++){
                                  for(int S0=std::abs(Sab-Sd); S0<=Sab+Sd; S0++){
			           std::array<int,15> keyrho={Na,Nb,Np,Nd,xab.irrep.lambda(),xab.irrep.mu(),Sab,
					              xd.irrep.lambda(),xd.irrep.mu(),Sd,x0.irrep.lambda(),x0.irrep.mu(),S0,rho0,rho};
                                   for(int Jab=std::abs(Lab-Sab); Jab<=Lab+Sab; Jab++){
			            for(int Jd=std::abs(Ld-Sd); Jd<=Ld+Sd; Jd++){
				     for(int J0=std::abs(L0-S0); J0<=L0+S0; J0++){
                                      if(std::abs(Jab-Jd)>J0 || J0>Jab+Jd)continue;
				      double factor7=factor6*(2*Jab+1)*(2*Jd+1)*(2*J0+1)*sqrt(double((2*Sab+1)*(2*Sd+1)*(2*S0+1)))
					             *wig9jj(2*Lab,2*Ld,2*L0,2*Sab,2*Sd,2*S0,2*Jab,2*Jd,2*J0);
				      if((Jab+Jd+J0)%2!=0)factor7=-factor7;
                                      for(std::map<std::array<int,2>,std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
                                       int SS1p=it2->first[0];
                                       int SS1=it2->first[1];
				       for(int II1ppp=std::abs(2*L1ppp-SS1p); II1ppp<=2*L1ppp+SS1p; II1ppp=II1ppp+2){
                                        for(int II1pp=std::abs(2*L1pp-SS1); II1pp<=2*L1pp+SS1; II1pp=II1pp+2){
					 if(std::abs(II1pp-2*J0)>II1ppp || II1ppp>II1pp+2*J0)continue;
					 double factor8=factor7*(II1ppp+1)*(II1pp+1)*sqrt(double(SS1p+1))
						        *wig9jj(2*L1pp,SS1,II1pp,2*L0,2*S0,2*J0,2*L1ppp,SS1p,II1ppp);
                                         for(int SSp=std::abs(SS1p-1); SSp<=SS1p+1; SSp=SSp+2){
                                          for(int JJ=std::abs(2*Lp-SSp); JJ<=2*Lp+SSp; JJ=JJ+2){
				           for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp=jjp+2){
				            if(std::abs(II1ppp-jjp)>JJ || JJ>II1ppp+jjp)continue;
					    double factor9=factor8*(jjp+1)*sqrt(double(SSp+1))
						           *wig9jj(2*L1ppp,2*lp,2*Lp,SS1p,1,SSp,II1ppp,jjp,JJ);
					    for(int SS=std::abs(SS1-1); SS<=SS1+1; SS=SS+2){
					     if(std::abs(2*L-SS)>JJ || JJ>2*L+SS)continue;
					     std::array<int,19> key1={lm1p,mu1p,SS1p,Np,xp.irrep.lambda(),xp.irrep.mu(),kappap,Lp,SSp,
                                                                      lm1,mu1,SS1,N,x.irrep.lambda(),x.irrep.mu(),kappa,L,SS,JJ};
                                             for(int jj=std::abs(2*l-1); jj<=2*l+1; jj=jj+2){
					      if(std::abs(II1pp-jj)>JJ || JJ>II1pp+jj || std::abs(2*J0-jjp)>jj || jj>2*J0+jjp)continue;
					      double factor10=factor9*sqrt(double((SS+1)*(jj+1)))*wig6jj(II1pp,2*J0,II1ppp,jjp,JJ,jj)
						              *wig9jj(2*L1pp,2*l,2*L,SS1,1,SS,II1pp,jj,JJ);
					      if(((jj+II1pp+JJ)/2)%2!=0)factor10=-factor10;
                                              for(int jjd=std::abs(2*ld-1); jjd<=2*ld+1; jjd=jjd+2){
					       if(std::abs(jjp-jjd)>2*Jd || 2*Jd>jjp+jjd || std::abs(jjd-2*Jab)>jj || jj>jjd+2*Jab)continue;
					       double factor11=factor10*sqrt(double(jjd+1))*wig6jj(jjd,2*Jd,jjp,2*J0,jj,2*Jab)
						               *wig9jj(2*lp,2*ld,2*Ld,1,1,2*Sd,jjp,jjd,2*Jd);
                                               for(int jja=std::abs(2*la-1); jja<=2*la+1; jja=jja+2){
                                                for(int jjb=std::abs(2*lb-1); jjb<=2*lb+1; jjb=jjb+2){
                                                 if(std::abs(jja-jjb)>2*Jab || 2*Jab>jja+jjb)continue;
                                                 double factor12=factor11*sqrt(double((jja+1)*(jjb+1)))
							         *wig9jj(2*la,2*lb,2*Lab,1,1,2*Sab,jja,jjb,2*Jab);
				                 double vmep=like_int(Na,la,jja,Nb,lb,jjb,Nd,ld,jjd,N,l,jj,Jab,vp);
				                 double vmen=like_int(Na,la,jja,Nb,lb,jjb,Nd,ld,jjd,N,l,jj,Jab,vn);
				                 double vmepn=pn_int(Na,la,jja,Nb,lb,jjb,Nd,ld,jjd,N,l,jj,Jab,vpn);
				                 std::map<std::array<int,19>,std::map<std::array<int,4>,std::array<double,2>>>::iterator
				                   iter=pot_exh.find(key1);
				                 if(iter==pot_exh.end()){
                                                  std::map<std::array<int,4>,std::array<double,2>> mes;
                                                  for(std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>::iterator
                                                      it3=it2->second.begin(); it3!=it2->second.end(); it3++){
                                                   std::array<int,4> key2=it3->first;
				                   std::array<double,2> me;
				                   if(Na<=etamaxp && Nb<=etamaxp && Np<=etamaxp && Nd<=etamaxp){
                                                    me[0]=0.5*factor12*vmep*it3->second[keyrho][0];
				                   }else{
                                                    me[0]=0.0;
				                   }
						   if(Na<=etamaxn && Nb<=etamaxn && Np<=etamaxn && Nd<=etamaxn){
                                                    me[1]=0.5*factor12*vmen*it3->second[keyrho][1];
                                                   }else{
                                                    me[1]=0.0;
                                                   }
				                   //if(Na<=etamaxn && Nb<=etamaxp && Np<=etamaxp && Nd<=etamaxn)me[0]+= in contruction
                                                   if(Na<=etamaxp && Nb<=etamaxn && Np<=etamaxn && Nd<=etamaxp)
					             me[1]+=factor12*vmepn*it3->second[keyrho][2];
					             /*{std::array<int,15> keyrhopn={Nb,Na,Nd,Np,xab.irrep.lambda(),xab.irrep.mu(),Sab,//
                                                      xd.irrep.lambda(),xd.irrep.mu(),Sd,x0.irrep.lambda(),x0.irrep.mu(),S0,rho0,rho};//
						      double hustota=it3->second[keyrhopn][2];//
						      if((Nb+Na+Nd+Np+xab.irrep.lambda()+xab.irrep.mu()+Sab//
                                                          +xd.irrep.lambda()+xd.irrep.mu()+Sd)%2!=0)hustota=-hustota;//
                                                      me[1]+=factor12*vmepn*hustota;}*/
				                   mes[key2]=me;
				                  }
				                  pot_exh[key1]=mes;
				                 }else{
                                                  for(std::map<std::array<int,4>,std::map<std::array<int,15>,std::array<double,3>>>::iterator
                                                      it3=it2->second.begin(); it3!=it2->second.end(); it3++){
                                                   std::array<int,4> key2=it3->first;
                                                   if(Na<=etamaxp && Nb<=etamaxp && Np<=etamaxp && Nd<=etamaxp)
						     pot_exh[key1][key2][0]+=0.5*factor12*vmep*it3->second[keyrho][0];
						   if(Na<=etamaxn && Nb<=etamaxn && Np<=etamaxn && Nd<=etamaxn)
					             pot_exh[key1][key2][1]+=0.5*factor12*vmen*it3->second[keyrho][1];
						   //if(Na<=etamaxn && Nb<=etamaxp && Np<=etamaxp && Nd<=etamaxn)
						   //  pot_exh[key1][key2][0]+= in construction
						   if(Na<=etamaxp && Nb<=etamaxn && Np<=etamaxn && Nd<=etamaxp)
					             pot_exh[key1][key2][1]+=factor12*vmepn*it3->second[keyrho][2];
						     /*{std::array<int,15> keyrhopn={Nb,Na,Nd,Np,xab.irrep.lambda(),xab.irrep.mu(),Sab,//
                                                      xd.irrep.lambda(),xd.irrep.mu(),Sd,x0.irrep.lambda(),x0.irrep.mu(),S0,rho0,rho};//
                                                      double hustota=it3->second[keyrhopn][2];//
                                                      if((Nb+Na+Nd+Np+xab.irrep.lambda()+xab.irrep.mu()+Sab//
                                                          +xd.irrep.lambda()+xd.irrep.mu()+Sd)%2!=0)hustota=-hustota;//
                                                      pot_exh[key1][key2][1]+=factor12*vmepn*hustota;}*/
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
  }
 }
}

std::ofstream PkAViAfile("PkAViA.dat");
if(!PkAViAfile){
  std::cout<<"Could not open PkAViA file"<<std::endl;
  return EXIT_FAILURE;
} 
//PkAViAfile<<"lm1' mu1' kappa1' L1' SS1' N' lm' mu' kappa' L' SS' lm1 mu1 kappa1 L1 SS1 N lm mu kappa L SS JJ prot_proj neut_proj"<<std::endl;
for(std::map<std::array<int,19>,std::map<std::array<int,4>,std::array<double,2>>>::iterator it1=pot_exh.begin(); it1!=pot_exh.end(); it1++){
  int lm1p=it1->first[0];
  int mu1p=it1->first[1];
  int SS1p=it1->first[2];
  int Np=it1->first[3];
  int lmp=it1->first[4];
  int mup=it1->first[5];
  int kappap=it1->first[6];
  int Lp=it1->first[7];
  int SSp=it1->first[8];
  int lm1=it1->first[9];
  int mu1=it1->first[10];
  int SS1=it1->first[11];
  int N=it1->first[12];
  int lm=it1->first[13];
  int mu=it1->first[14];
  int kappa=it1->first[15];
  int L=it1->first[16];
  int SS=it1->first[17];
  int JJ=it1->first[18];
  for(std::map<std::array<int,4>,std::array<double,2>>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
    PkAViAfile<<lm1p<<" "<<mu1p<<" "<<it2->first[0]<<" "<<it2->first[1]<<" "<<SS1p<<" "<<Np<<" "
	      <<lmp<<" "<<mup<<" "<<kappap<<" "<<Lp<<" "<<SSp<<" "
	      <<lm1<<" "<<mu1<<" "<<it2->first[2]<<" "<<it2->first[3]<<" "<<SS1<<" "<<N<<" "
	      <<lm<<" "<<mu<<" "<<kappa<<" "<<L<<" "<<SS<<" "<<JJ<<" "<<it2->second[0]<<" "<<it2->second[1]<<std::endl;
  }
}
PkAViAfile.close();

return EXIT_SUCCESS;
}

#include <cstdio>
#include <fstream>
#include <istream>
#include <iostream>

#include "spncci/computation_control.h"
#include "wigxjpf.h"

int main(int argc, char **argv){

int II1p=2;
int pi1p=1;
int II1=2;
int pi1=1;

if(argc!=2){
  std::cout<<"Usage: "<<argv[0]<<" Nclustermax"<<std::endl;
  return EXIT_FAILURE;
}

const int Nmax=std::stoi(argv[1]);
std::cout<<"Nclustermax = "<<Nmax<<std::endl;
int JJmax=II1+1+2*Nmax;

//***************************************************
// Read MEs of P_iA
//***************************************************

std::map<std::array<int,10>,std::map<std::array<int,5>,std::array<double,2>>> PkA;
// First key is {lm1',mu1',kappa1',L1',2S1',lm1,mu1,kappa1,L1,2S1}
// Second key is {N',N,lm,mu,2S}
// Values are proton and neutron

std::ifstream PkAfile("PiA.dat");
if(!PkAfile){
  std::cout<<"Could not open PiA.dat file!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int lm1p,mu1p,kappa1p,L1p,SS1p,Np,lm1,mu1,kappa1,L1,SS1,N,lm,mu,SS;
  PkAfile>>lm1p>>mu1p>>kappa1p>>L1p>>SS1p>>Np>>lm1>>mu1>>kappa1>>L1>>SS1>>N>>lm>>mu>>SS;
  if(PkAfile){
    std::array<double,2> mes;
    PkAfile>>mes[0]>>mes[1];
    std::array<int,10> key1={lm1p,mu1p,kappa1p,L1p,SS1p,lm1,mu1,kappa1,L1,SS1};
    std::array<int,5> key2={Np,N,lm,mu,SS};
    std::map<std::array<int,10>,std::map<std::array<int,5>,std::array<double,2>>>::iterator it=PkA.find(key1);
    if(it==PkA.end()){
      std::map<std::array<int,5>,std::array<double,2>> map;
      map[key2]=mes;
      PkA[key1]=map;
    }else{
      PkA[key1][key2]=mes;
    }
  }else{
    break;
  }
}
PkAfile.close();

//**************************************************
// Read MEs of V_iA(1-P_ia) + V_jA
//**************************************************

std::map<std::array<int,8>,std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>> ViAPiA;
// First key is {lm1p,mu1p,kappa1p,L1p,lmp,mup,kappap,Lp}
// Second key is {lm1,mu1,kappa1,L1,lm,mu,kappa,L}
// Third key is {SS1p,SSp}
// Fourth key is {SS1,SS}
// Fifth key is {Np,N,JJ}
// Value is {proton projectile, neutron projectile}

std::ifstream ViAPiAfile("ViAPiA.dat");
if(!ViAPiAfile){
  std::cout<<"Could not open ViAPiA.dat file!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int lm1p,mu1p,kappa1p,L1p,SS1p,Np,lmp,mup,kappap,Lp,SSp,lm1,mu1,kappa1,L1,SS1,N,lm,mu,kappa,L,SS,JJ;
  ViAPiAfile>>lm1p>>mu1p>>kappa1p>>L1p>>SS1p>>Np>>lmp>>mup>>kappap>>Lp>>SSp>>lm1>>mu1>>kappa1>>L1>>SS1>>N>>lm>>mu>>kappa>>L>>SS>>JJ;
  if(ViAPiAfile){
    std::array<double,2> value;
    ViAPiAfile>>value[0]>>value[1];
    std::array<int,8> key1={lm1p,mu1p,kappa1p,L1p,lmp,mup,kappap,Lp};
    std::array<int,8> key2={lm1,mu1,kappa1,L1,lm,mu,kappa,L};
    std::array<int,2> key3={SS1p,SSp};
    std::array<int,2> key4={SS1,SS};
    std::array<int,3> key5={Np,N,JJ};
    std::map<std::array<int,8>,std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>>::iterator it1=ViAPiA.find(key1);
    if(it1==ViAPiA.end()){
      std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>> value1;
      std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>> value2;
      std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>> value3;
      std::map<std::array<int,3>,std::array<double,2>> value4;
      value4[key5]=value;
      value3[key4]=value4;
      value2[key3]=value3;
      value1[key2]=value2;
      ViAPiA[key1]=value1;
    }else{
      std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>::iterator it2=ViAPiA[key1].find(key2);
      if(it2==ViAPiA[key1].end()){
        std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>> value2;
        std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>> value3;
        std::map<std::array<int,3>,std::array<double,2>> value4;
        value4[key5]=value;
        value3[key4]=value4;
        value2[key3]=value3;
        ViAPiA[key1][key2]=value2;
      }else{
        std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>::iterator it3=ViAPiA[key1][key2].find(key3);
        if(it3==ViAPiA[key1][key2].end()){
          std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>> value3;
          std::map<std::array<int,3>,std::array<double,2>> value4;
          value4[key5]=value;
          value3[key4]=value4;
          ViAPiA[key1][key2][key3]=value3;
        }else{
          std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>::iterator it4=ViAPiA[key1][key2][key3].find(key4);
          if(it4==ViAPiA[key1][key2][key3].end()){
            std::map<std::array<int,3>,std::array<double,2>> value4;
            value4[key5]=value;
            ViAPiA[key1][key2][key3][key4]=value4;
          }else{
            ViAPiA[key1][key2][key3][key4][key5]=value;
          }
        }
      }
    }
  }else{
    break;
  }
}
ViAPiAfile.close();

//**************************************************
// Read MEs of V_iA*P_ka
//**************************************************

std::map<std::array<int,8>,std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>> ViAPkA;
// First key is {lm1p,mu1p,kappa1p,L1p,lmp,mup,kappap,Lp}
// Second key is {lm1,mu1,kappa1,L1,lm,mu,kappa,L}
// Third key is {SS1p,SSp}
// Fourth key is {SS1,SS}
// Fifth key is {Np,N,JJ}
// Value is {proton projectile, neutron projectile}

std::ifstream ViAPkAfile("ViAPkA.dat");
if(!ViAPkAfile){
  std::cout<<"Could not open ViAPkA.dat file!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int lm1p,mu1p,kappa1p,L1p,SS1p,Np,lmp,mup,kappap,Lp,SSp,lm1,mu1,kappa1,L1,SS1,N,lm,mu,kappa,L,SS,JJ;
  ViAPkAfile>>lm1p>>mu1p>>kappa1p>>L1p>>SS1p>>Np>>lmp>>mup>>kappap>>Lp>>SSp>>lm1>>mu1>>kappa1>>L1>>SS1>>N>>lm>>mu>>kappa>>L>>SS>>JJ;
  if(ViAPkAfile){
    std::array<double,2> value;
    ViAPkAfile>>value[0]>>value[1];
    std::array<int,8> key1={lm1p,mu1p,kappa1p,L1p,lmp,mup,kappap,Lp};
    std::array<int,8> key2={lm1,mu1,kappa1,L1,lm,mu,kappa,L};
    std::array<int,2> key3={SS1p,SSp};
    std::array<int,2> key4={SS1,SS};
    std::array<int,3> key5={Np,N,JJ};
    std::map<std::array<int,8>,std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>>::iterator it1=ViAPkA.find(key1);
    if(it1==ViAPkA.end()){
      std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>> value1;
      std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>> value2;
      std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>> value3;
      std::map<std::array<int,3>,std::array<double,2>> value4;
      value4[key5]=value;
      value3[key4]=value4;
      value2[key3]=value3;
      value1[key2]=value2;
      ViAPkA[key1]=value1;
    }else{
      std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>::iterator it2=ViAPkA[key1].find(key2);
      if(it2==ViAPkA[key1].end()){
        std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>> value2;
        std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>> value3;
        std::map<std::array<int,3>,std::array<double,2>> value4;
        value4[key5]=value;
        value3[key4]=value4;
        value2[key3]=value3;
        ViAPkA[key1][key2]=value2;
      }else{
        std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>::iterator it3=ViAPkA[key1][key2].find(key3);
        if(it3==ViAPkA[key1][key2].end()){
          std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>> value3;
          std::map<std::array<int,3>,std::array<double,2>> value4;
          value4[key5]=value;
          value3[key4]=value4;
          ViAPkA[key1][key2][key3]=value3;
        }else{
          std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>::iterator it4=ViAPkA[key1][key2][key3].find(key4);
          if(it4==ViAPkA[key1][key2][key3].end()){
            std::map<std::array<int,3>,std::array<double,2>> value4;
            value4[key5]=value;
            ViAPkA[key1][key2][key3][key4]=value4;
          }else{
            ViAPkA[key1][key2][key3][key4][key5]=value;
          }
        }
      }
    }
  }else{
    break;
  }
}
ViAPkAfile.close();

//**************************************************
// Read MEs of P_ka*V_iA
//**************************************************

std::map<std::array<int,8>,std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>> PkAViA;
// First key is {lm1p,mu1p,kappa1p,L1p,lmp,mup,kappap,Lp}
// Second key is {lm1,mu1,kappa1,L1,lm,mu,kappa,L}
// Third key is {SS1p,SSp}
// Fourth key is {SS1,SS}
// Fifth key is {Np,N,JJ}
// Value is {proton projectile, neutron projectile}

std::ifstream PkAViAfile("PkAViA.dat");
if(!PkAViAfile){
  std::cout<<"Could not open PkAViA.dat file!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int lm1p,mu1p,kappa1p,L1p,SS1p,Np,lmp,mup,kappap,Lp,SSp,lm1,mu1,kappa1,L1,SS1,N,lm,mu,kappa,L,SS,JJ;
  PkAViAfile>>lm1p>>mu1p>>kappa1p>>L1p>>SS1p>>Np>>lmp>>mup>>kappap>>Lp>>SSp>>lm1>>mu1>>kappa1>>L1>>SS1>>N>>lm>>mu>>kappa>>L>>SS>>JJ;
  if(PkAViAfile){
    std::array<double,2> value;
    PkAViAfile>>value[0]>>value[1];
    std::array<int,8> key1={lm1p,mu1p,kappa1p,L1p,lmp,mup,kappap,Lp};
    std::array<int,8> key2={lm1,mu1,kappa1,L1,lm,mu,kappa,L};
    std::array<int,2> key3={SS1p,SSp};
    std::array<int,2> key4={SS1,SS};
    std::array<int,3> key5={Np,N,JJ};
    std::map<std::array<int,8>,std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>>::iterator it1=PkAViA.find(key1);
    if(it1==PkAViA.end()){
      std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>> value1;
      std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>> value2;
      std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>> value3;
      std::map<std::array<int,3>,std::array<double,2>> value4;
      value4[key5]=value;
      value3[key4]=value4;
      value2[key3]=value3;
      value1[key2]=value2;
      PkAViA[key1]=value1;
    }else{
      std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>::iterator it2=PkAViA[key1].find(key2);
      if(it2==PkAViA[key1].end()){
        std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>> value2;
        std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>> value3;
        std::map<std::array<int,3>,std::array<double,2>> value4;
        value4[key5]=value;
        value3[key4]=value4;
        value2[key3]=value3;
        PkAViA[key1][key2]=value2;
      }else{
        std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>::iterator it3=PkAViA[key1][key2].find(key3);
        if(it3==PkAViA[key1][key2].end()){
          std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>> value3;
          std::map<std::array<int,3>,std::array<double,2>> value4;
          value4[key5]=value;
          value3[key4]=value4;
          PkAViA[key1][key2][key3]=value3;
        }else{
          std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>::iterator it4=PkAViA[key1][key2][key3].find(key4);
          if(it4==PkAViA[key1][key2][key3].end()){
            std::map<std::array<int,3>,std::array<double,2>> value4;
            value4[key5]=value;
            PkAViA[key1][key2][key3][key4]=value4;
          }else{
            PkAViA[key1][key2][key3][key4][key5]=value;
          }
        }
      }
    }
  }else{
    break;
  }
}
PkAViAfile.close();

spncci::InitializeSpNCCI();

//*****************************************
// Calculate kernels
//*****************************************

std::ofstream kernel_file("kernels.dat");
if(!kernel_file){
  std::cout<<"Could not open kernels.dat file"<<std::endl;
  return EXIT_FAILURE;
}

for(int JJ=std::abs(II1%2-1); JJ<=JJmax; JJ=JJ+2){
for(int pi=-1; pi<=1; pi=pi+2){
kernel_file<<JJ<<" "<<pi<<std::endl;
for(int ssp=std::abs(II1p-1); ssp<=II1p+1; ssp=ssp+2){
int lpmin=std::abs(ssp-JJ)/2;
if((pi*pi1p==1 && lpmin%2==1)||(pi*pi1p==-1 && lpmin%2==0))lpmin++;
for(int lp=lpmin; lp<=std::min((ssp+JJ)/2,Nmax); lp=lp+2){
for(int Np=lp; Np<=Nmax; Np=Np+2){
for(int ss=std::abs(II1-1); ss<=II1+1; ss=ss+2){
double factor1=sqrt(double((ssp+1)*(ss+1)*(II1p+1)*(II1+1)));
int lmin=std::abs(ss-JJ)/2;
if((pi*pi1==1 && lmin%2==1)||(pi*pi1==-1 && lmin%2==0))lmin++;
for(int l=lmin; l<=std::min((ss+JJ)/2,Nmax); l=l+2){
for(int N=l; N<=Nmax; N=N+2){
std::array<int,3> key5={Np,N,JJ};
  //*************************************************
  // Norm kernels
  //*************************************************
  double nkerp=0.0;
  double nkern=0.0;
  for(std::map<std::array<int,10>,std::map<std::array<int,5>,std::array<double,2>>>::iterator it=PkA.begin(); it!=PkA.end(); it++){
    u3::SU3 x1p(it->first[0],it->first[1]);
    int kappa1p=it->first[2];
    int L1p=it->first[3];
    int SS1p=it->first[4];
    u3::SU3 x1(it->first[5],it->first[6]);
    int kappa1=it->first[7];
    int L1=it->first[8];
    int SS1=it->first[9];
    for(auto x : u3::KroneckerProduct(x1,u3::SU3(N,0))){
      if(u3::OuterMultiplicity(x1p,u3::SU3(Np,0),x.irrep)==0)continue;
      int lm=x.irrep.lambda();
      int mu=x.irrep.mu();
      for(auto L : u3::BranchingSO3(x.irrep)){
	double sum_kappa_p=0.0;
	double sum_kappa_n=0.0;
        for(int kappa=1; kappa<=L.tag; kappa++){
	  double sum_jp_p=0.0;
	  double sum_jp_n=0.0;
	  for(int jjp=std::abs(1-2*lp); jjp<=1+2*lp; jjp=jjp+2){
            if(JJ<std::abs(jjp-II1p) || JJ>jjp+II1p)continue;
	    double sum_j_p=0.0;
	    double sum_j_n=0.0;
	    for(int jj=std::abs(1-2*l); jj<=1+2*l; jj=jj+2){
              if(JJ<std::abs(jj-II1) || JJ>jj+II1)continue;
	      double sum_S_p=0.0;
	      double sum_S_n=0.0;
	      for(int SS=std::abs(SS1-1); SS<=SS1+1; SS=SS+2){
		if(SS<std::abs(SS1p-1) || SS>SS1p+1)continue;
                double factor6=(SS+1)*wig9jj(2*L1p,2*lp,2*L.irrep,SS1p,1,SS,II1p,jjp,JJ)*wig9jj(2*L1,2*l,2*L.irrep,SS1,1,SS,II1,jj,JJ);
                std::array<int,5> key={Np,N,lm,mu,SS};
		sum_S_p+=factor6*it->second[key][0];
		sum_S_n+=factor6*it->second[key][1];
              }
              double factor5=(jj+1)*wig6jj(II1,1,ss,2*l,JJ,jj);
	      int phase=(II1p+II1+jjp+jj)/2+JJ;
	      if(phase%2!=0)factor5=-factor5;
	      sum_j_p+=factor5*sum_S_p;
	      sum_j_n+=factor5*sum_S_n;
	    }
            double factor4=(jjp+1)*wig6jj(II1p,1,ssp,2*lp,JJ,jjp);
	    sum_jp_p+=factor4*sum_j_p;
	    sum_jp_n+=factor4*sum_j_n;
	  }
	  double factor3=u3::W(x1p,kappa1p,L1p,u3::SU3(Np,0),1,lp,x.irrep,kappa,L.irrep,1)
		        *u3::W(x1,kappa1,L1,u3::SU3(N,0),1,l,x.irrep,kappa,L.irrep,1);
          sum_kappa_p+=factor3*sum_jp_p;
	  sum_kappa_n+=factor3*sum_jp_n;
        }
	int factor2=2*L.irrep+1;
	nkerp+=factor2*sum_kappa_p;
	nkern+=factor2*sum_kappa_n;
      }
    }
  }
  nkerp*=-factor1;
  nkern*=-factor1;
//nkern=0.0;
  //**************************************************
  // Hamiltonian kernels
  //**************************************************
  double pkerdp=0.0;
  double pkerdn=0.0;
  double pkere1p=0.0;
  double pkere1n=0.0;
  double pkere2p=0.0;
  double pkere2n=0.0;
  for(std::map<std::array<int,8>,std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>>::iterator it1=ViAPiA.begin(); it1!=ViAPiA.end(); it1++){
   int L1p=it1->first[3];
   int Lp=it1->first[7];
   if(std::abs(L1p-lp)>Lp || Lp>L1p+lp)continue;
   int lm1p=it1->first[0];
   int mu1p=it1->first[1];
   int lmp=it1->first[4];
   int mup=it1->first[5];
   u3::SU3 x1p(lm1p,mu1p);
   u3::SU3 xp(lmp,mup);
   if(u3::OuterMultiplicity(x1p,u3::SU3(Np,0),xp)==0)continue;
   int kappa1p=it1->first[2];
   int kappap=it1->first[6];
//std::cout<<"lm1p,mu1p,kappa1p,L1p,lmp,mup,kappap,Lp: "<<lm1p<<" "<<mu1p<<" "<<kappa1p<<" "<<L1p<<" "<<lmp<<" "<<mup<<" "<<kappap<<" "<<Lp<<std::endl;
   double factor2=factor1*sqrt(double(2*Lp+1))*u3::W(x1p,kappa1p,L1p,u3::SU3(Np,0),1,lp,xp,kappap,Lp,1);
   for(std::map<std::array<int,8>,std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
    int L1=it2->first[3];
    int L=it2->first[7];
    if(std::abs(L1-l)>L || L>L1+l)continue;
    int lm1=it2->first[0];
    int mu1=it2->first[1];
    int lm=it2->first[4];
    int mu=it2->first[5];
    u3::SU3 x1(lm1,mu1);
    u3::SU3 x(lm,mu);
    if(u3::OuterMultiplicity(x1,u3::SU3(N,0),x)==0)continue;
    int kappa1=it2->first[2];
    int kappa=it2->first[6];
//std::cout<<"lm1,mu1,kappa1,L1,lm,mu,kappa,L: "<<lm1<<" "<<mu1<<" "<<kappa1<<" "<<L1<<" "<<lm<<" "<<mu<<" "<<kappa<<" "<<L<<std::endl;
    double factor3=factor2*sqrt(double(2*L+1))*u3::W(x1,kappa1,L1,u3::SU3(N,0),1,l,x,kappa,L,1);
    for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp=jjp+2){
     if(std::abs(jjp-II1p)>JJ || JJ>jjp+II1p)continue;
//std::cout<<"jjp: "<<jjp<<std::endl;
     double factor4=factor3*(jjp+1)*wig6jj(II1p,1,ssp,2*lp,JJ,jjp);
     if(((II1p+JJ+jjp)/2)%2!=0)factor4=-factor4;
     for(std::map<std::array<int,2>,std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>>::iterator it3=it2->second.begin(); it3!=it2->second.end(); it3++){
      int SS1p=it3->first[0];
      int SSp=it3->first[1];
//std::cout<<"SS1p,SSp: "<<SS1p<<" "<<SSp<<std::endl;
      double factor5=factor4*sqrt(double(SSp+1))*wig9jj(2*L1p,2*lp,2*Lp,SS1p,1,SSp,II1p,jjp,JJ);
      for(int jj=std::abs(2*l-1); jj<=2*l+1; jj=jj+2){
       if(std::abs(jj-II1)>JJ || JJ>jj+II1)continue;
//std::cout<<"jj: "<<jj<<std::endl;
       double factor6=factor5*(jj+1)*wig6jj(II1,1,ss,2*l,JJ,jj);
       if(((II1+JJ+jj)/2)%2!=0)factor6=-factor6;
       for(std::map<std::array<int,2>,std::map<std::array<int,3>,std::array<double,2>>>::iterator it4=it3->second.begin(); it4!=it3->second.end(); it4++){
	std::map<std::array<int,3>,std::array<double,2>>::iterator iter=it4->second.find(key5);
	if(iter==it4->second.end())continue;
        int SS1=it4->first[0];
        int SS=it4->first[1];
//std::cout<<"SS1,SS: "<<SS1<<" "<<SS<<std::endl;
        double factor7=factor6*sqrt(double(SS+1))*wig9jj(2*L1,2*l,2*L,SS1,1,SS,II1,jj,JJ);
        pkerdp+=factor7*it4->second[key5][0];
        pkerdn+=factor7*it4->second[key5][1];
        pkere1p+=factor7*ViAPkA[it1->first][it2->first][it3->first][it4->first][key5][0];
        pkere1n+=factor7*ViAPkA[it1->first][it2->first][it3->first][it4->first][key5][1];
        pkere2p+=factor7*PkAViA[it1->first][it2->first][it3->first][it4->first][key5][0];
        pkere2n+=factor7*PkAViA[it1->first][it2->first][it3->first][it4->first][key5][1];
//if(lmp==lm&&mup==mu&&kappap==kappa&&Lp==L&&SSp==SS)nkern+=factor7*PkA[{lm1p,mu1p,kappa1p,L1p,SS1p,lm1,mu1,kappa1,L1,SS1}][{Np,N,lm,mu,SS}][1];
       }
      }
     }
    }
   }
  }
//  kernel_file<<ssp<<" "<<lp<<" "<<Np<<" "<<ss<<" "<<l<<" "<<N<<" "<<nkern<<" "<<pkerdn<<" "<<-0.5*(pkere1n+pkere2n)<<" "<<pkerdn-0.5*(pkere1n+pkere2n)<<std::endl;
  kernel_file<<ssp<<" "<<lp<<" "<<Np<<" "<<ss<<" "<<l<<" "<<N<<" "<<nkern<<" "<<pkerdn<<" "<<-pkere1n<<"("<<-pkere2n<<")"<<std::endl;
}
}
}
}
}
}
}
}
kernel_file.close();

return EXIT_SUCCESS;
}

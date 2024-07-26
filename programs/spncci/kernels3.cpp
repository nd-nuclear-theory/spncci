#include <cstdio>
#include <fstream>
#include <istream>
#include <iostream>

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

int II1p=2;
int pi1p=1;
int II1=2;
int pi1=1;
int N0=2;

if(argc!=2){
  std::cout<<"Usage: "<<argv[0]<<" Nclustermax"<<std::endl;
  return EXIT_FAILURE;
}

const int Nmax=std::stoi(argv[1]);
std::cout<<"Nclustermax = "<<Nmax<<std::endl;
int JJmax=II1+1+2*Nmax;


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

std::map<std::array<int,7>,std::array<double,2>> obdrmes;
// obdrmes[Np,lp,jjp,N,l,jj,J0][proton,neutron]
std::ifstream obd("spncci_obdrmes.dat");
if(!obd){
  std::cout<<"Could not open spncci_obdrmes.dat file!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int Np,lp,jjp,N,l,jj,J0;
  obd>>Np>>lp>>jjp>>N>>l>>jj>>J0;
  if(obd){
    std::array<double,2> rmes;
    obd>>rmes[0]>>rmes[1];
    obdrmes[{Np,lp,jjp,N,l,jj,J0}]=rmes;
  }else{
    break;
  }
}
obd.close();

std::map<std::array<int,15>,std::array<double,3>> tbdrmes;
// tbdrmes[N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0][proton,neutron,proton-neutron]
std::ifstream tbd("spncci_tbdrmes.dat");
if(!tbd){
  std::cout<<"Could not open spncci_tbdrmes.dat file!"<<std::endl;
  return EXIT_FAILURE;
}
while(true){
  int N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0;
  tbd>>N1>>l1>>jj1>>N2>>l2>>jj2>>N3>>l3>>jj3>>N4>>l4>>jj4>>Jf>>Ji>>J0;
  if(tbd){
    std::array<double,3> rmes;
    tbd>>rmes[0]>>rmes[1]>>rmes[2];
    tbdrmes[{N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0}]=rmes;
  }else{
    break;
  }
}
tbd.close();

spncci::InitializeSpNCCI();

int JJmy=1;
int pimy=1;
int sspmy=1;
int lpmy=0;
int Npmy=0;
int ssmy=1;
int lmy=0;
int Nmy=0;

std::ofstream kernel_file("theirkernels.dat");
if(!kernel_file){
  std::cout<<"Could not open theirkernels.dat file"<<std::endl;
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
int lmin=std::abs(ss-JJ)/2;
if((pi*pi1==1 && lmin%2==1)||(pi*pi1==-1 && lmin%2==0))lmin++;
for(int l=lmin; l<=std::min((ss+JJ)/2,Nmax); l=l+2){ 
for(int N=l; N<=Nmax; N=N+2){
bool write,write2;
if(JJ==JJmy&&pi==pimy&&ssp==sspmy&&lp==lpmy&&Np==Npmy&&ss==ssmy&&l==lmy&&N==Nmy){write=true; write2=true;}else{write=false; write2=false;}

write=false;

 double nker=0.0;
 double factor1=sqrt(double((ssp+1)*(ss+1)));
 for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp=jjp+2){
  if(std::abs(jjp-II1p)>JJ || JJ>jjp+II1p)continue;
  double factor2=factor1*sqrt(double(jjp+1))*wig6jj(II1p,1,ssp,2*lp,JJ,jjp);
  if(((II1p-jjp+JJ)/2)%2!=0)factor2=-factor2;
  for(int jj=std::abs(2*l-1); jj<=2*l+1; jj=jj+2){
   if(std::abs(jj-II1)>JJ || JJ>jj+II1)continue;
   double factor3=factor2*sqrt(double(jj+1))*wig6jj(II1,1,ss,2*l,JJ,jj);
   for(int J0=std::abs(jj-jjp)/2; J0<=(jj+jjp)/2; J0++){
    if(std::abs(II1-2*J0)>II1p || II1p>II1+2*J0)continue;
    double factor4=factor3*(2*J0+1)*wig6jj(II1,2*J0,II1p,jjp,JJ,jj);
    nker+=factor4*obdrmes[{N,l,jj,Np,lp,jjp,J0}][1];
   }
  }
 }

 double pkerd=0.0;
 factor1=-sqrt(double((ssp+1)*(ss+1)));
if(write)std::cout<<"factor1="<<factor1<<std::endl;
 for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp=jjp+2){
  if(std::abs(jjp-II1p)>JJ || JJ>jjp+II1p)continue;
  double factor2=factor1*sqrt(double(jjp+1))*wig6jj(II1p,1,ssp,2*lp,JJ,jjp);
  if(((II1p+jjp+JJ)/2)%2!=0)factor2=-factor2;
if(write){
std::cout<<"jjp: "<<jjp<<std::endl;
std::cout<<"factor2="<<factor2<<std::endl;
}
  for(int jj=std::abs(2*l-1); jj<=2*l+1; jj=jj+2){
   if(std::abs(jj-II1)>JJ || JJ>jj+II1)continue;
   double factor3=factor2*sqrt(double(jj+1))*wig6jj(II1,1,ss,2*l,JJ,jj);
if(write){
std::cout<<" jj: "<<jj<<std::endl; 
std::cout<<" factor3="<<factor3<<std::endl;
}
   for(int Na=0; Na<=Nmax; Na++){ // Nmax is Nclustermax
    for(int Nb=Na%2; Nb<=Nmax; Nb=Nb+2){
     if(Na-Nb<-Nmax+1 || Na-Nb>Nmax-1)continue;
if(write)std::cout<<"  Na Nb: "<<Na<<" "<<Nb<<std::endl;
     for(int la=Na%2; la<=Na; la=la+2){
      for(int lb=Nb%2; lb<=Nb; lb=lb+2){
if(write)std::cout<<"   la lb: "<<la<<" "<<lb<<std::endl;
       for(int jja=std::abs(2*la-1); jja<=2*la+1; jja=jja+2){
        for(int jjb=std::abs(2*lb-1); jjb<=2*lb+1; jjb=jjb+2){
if(write)std::cout<<"    jja jjb: "<<jja<<" "<<jjb<<std::endl;
         for(int J0=std::abs(jja-jjb)/2; J0<=(jja+jjb)/2; J0++){
          if(std::abs(II1-2*J0)>II1p || II1p>II1+2*J0 || std::abs(2*J0-jjp)>jj || jj>2*J0+jjp)continue;
          double factor4=factor3*(2*J0+1)*wig6jj(II1,2*J0,II1p,jjp,JJ,jj);
if(write){
std::cout<<"     J0: "<<J0<<std::endl; 
std::cout<<"     factor4="<<factor4<<std::endl;
}
          for(int Jp=std::abs(jja-jjp)/2; Jp<=(jja+jjp)/2; Jp++){
           if(std::abs(jjb-jj)>2*Jp || 2*Jp>jjb+jj)continue;
	   double factor5=factor4*(2*Jp+1)*wig6jj(jjb,jja,2*J0,jjp,jj,2*Jp);
	   if(((2*Jp-jj-jjb)/2)%2!=0)factor5=-factor5;
if(write){
std::cout<<"      Jp: "<<Jp<<std::endl;
std::cout<<"      factor5="<<factor5<<std::endl;
}
	   double vmen=like_int(Na,la,jja,Np,lp,jjp,Nb,lb,jjb,N,l,jj,Jp,vn);
           double vmepn=pn_int(Na,la,jja,Np,lp,jjp,Nb,lb,jjb,N,l,jj,Jp,vpn);
	   //double vmepn=pn_int(Np,lp,jjp,Na,la,jja,N,l,jj,Nb,lb,jjb,Jp,vpn);
if(write)std::cout<<"pkerd=pkerd+factor5*(vmen*rhon+vmepn*rhop)="<<pkerd<<"+"<<factor5<<"*("<<vmen<<"*"<<obdrmes[{Na,la,jja,Nb,lb,jjb,J0}][1]<<"+"<<vmepn<<"*"<<obdrmes[{Na,la,jja,Nb,lb,jjb,J0}][0]<<")="<<pkerd+factor5*(vmen*obdrmes[{Na,la,jja,Nb,lb,jjb,J0}][1]+vmepn*obdrmes[{Na,la,jja,Nb,lb,jjb,J0}][0])<<std::endl;
	   pkerd+=factor5*(vmen*obdrmes[{Na,la,jja,Nb,lb,jjb,J0}][1]+vmepn*obdrmes[{Na,la,jja,Nb,lb,jjb,J0}][0]);
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

 double pkere=0.0;
 factor1=sqrt(double((ssp+1)*(ss+1)));
if(write2)std::cout<<"factor1="<<factor1<<std::endl;
 for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp=jjp+2){
  if(std::abs(jjp-II1p)>JJ || JJ>jjp+II1p)continue;
  double factor2=factor1*sqrt(double(jjp+1))*wig6jj(II1p,1,ssp,2*lp,JJ,jjp);
if(write2){
std::cout<<"jjp: "<<jjp<<std::endl;
std::cout<<" factor2="<<factor2<<std::endl;
}
  for(int jj=std::abs(2*l-1); jj<=2*l+1; jj=jj+2){
   if(std::abs(jj-II1)>JJ || JJ>jj+II1)continue;
   double factor3=factor2*sqrt(double(jj+1))*wig6jj(II1,1,ss,2*l,JJ,jj);
if(write2){
std::cout<<" jj: "<<jjp<<std::endl;
std::cout<<"  factor3="<<factor3<<std::endl;
}
   for(int Na=0; Na<=std::min(Nmax,N0+Nmax-1-N); Na++){
    for(int Nc=0; Nc<=Nmax; Nc++){
     int Ndmin=std::max(0,Na+N-Nc-Nmax+1);
     if((Na+N-Nc-Ndmin)%2!=0)Ndmin++;
     for(int Nd=Ndmin; Nd<=std::min(std::min(Nmax,N0+Nmax-1-Nc),Na+N-Nc+Nmax-1); Nd=Nd+2){
if(write2)std::cout<<"   Na Nc Nd: "<<Na<<" "<<Nc<<" "<<Nd<<std::endl;
      for(int la=Na%2; la<=Na; la=la+2){
       for(int lc=Nc%2; lc<=Nc; lc=lc+2){
        for(int ld=Nd%2; ld<=Nd; ld=ld+2){
if(write2)std::cout<<"    la lc ld: "<<la<<" "<<lc<<" "<<ld<<std::endl;
         for(int jja=std::abs(2*la-1); jja<=2*la+1; jja=jja+2){
	  for(int jjc=std::abs(2*lc-1); jjc<=2*lc+1; jjc=jjc+2){
           for(int jjd=std::abs(2*ld-1); jjd<=2*ld+1; jjd=jjd+2){
if(write2)std::cout<<"     jja jjc jjd: "<<jja<<" "<<jjc<<" "<<jjd<<std::endl;
	    double d=1.0;
	    if(Na==N && la==l && jja==jj)d=d/sqrt(2.0);
            if(Nc==Nd && lc==ld && jjc==jjd)d=d/sqrt(2.0);
            for(int Ja=std::abs(jja-jj)/2; Ja<=(jja+jj)/2; Ja++){
if(write2)std::cout<<"      Ja: "<<Ja<<std::endl;
	     for(int Jcd=std::abs(jjc-jjd)/2; Jcd<=(jjc+jjd)/2; Jcd++){
              if(std::abs(jja-jjp)>2*Jcd || 2*Jcd>jja+jjp)continue;
if(write2)std::cout<<"       Jcd: "<<Jcd<<std::endl;
	      for(int J0=std::abs(Ja-Jcd); J0<=Ja+Jcd; J0++){
	       if(std::abs(II1-2*J0)>II1p || II1p>II1+2*J0 || std::abs(2*J0-jjp)>jj || jj>2*J0+jjp)continue;
               double factor4=factor3*(2*J0+1)*sqrt(double((2*Ja+1)*(2*Jcd+1)))*wig6jj(II1,2*J0,II1p,jjp,JJ,jj)
		                                                               *wig6jj(2*Ja,2*Jcd,2*J0,jjp,jj,jja);
	       int phase=(II1p+JJ+jjd+jjc-jja)/2+Ja+Jcd+J0;
	       if(phase%2!=0)factor4=-factor4;
	       double vmen=like_int(Na,la,jja,Np,lp,jjp,Nd,ld,jjd,Nc,lc,jjc,Jcd,vn);
	       double vmepn=pn_int(Na,la,jja,Np,lp,jjp,Nd,ld,jjd,Nc,lc,jjc,Jcd,vpn);
	       //double vmepn=pn_int(Np,lp,jjp,Na,la,jja,Nc,lc,jjc,Nd,ld,jjd,Jcd,vpn);
if(write2){
std::cout<<"        J0: "<<J0<<std::endl;
std::cout<<"         factor4="<<factor4<<std::endl;
std::cout<<"         pkere=pkere+factor4*(0.5*vmen*rhon+vmepn*rhopn)="<<pkere<<"+"<<factor4<<"*(0.5*"<<vmen<<"*"<<tbdrmes[{Na,la,jja,N,l,jj,Nc,lc,jjc,Nd,ld,jjd,Ja,Jcd,J0}][1]<<"+"<<vmepn<<"*"<<tbdrmes[{Na,la,jja,N,l,jj,Nc,lc,jjc,Nd,ld,jjd,Ja,Jcd,J0}][2]<<")="<<pkere+factor4*(0.5*vmen*tbdrmes[{Na,la,jja,N,l,jj,Nc,lc,jjc,Nd,ld,jjd,Ja,Jcd,J0}][1]+vmepn*tbdrmes[{Na,la,jja,N,l,jj,Nc,lc,jjc,Nd,ld,jjd,Ja,Jcd,J0}][2])<<std::endl;
}
               pkere+=factor4*(0.5*vmen*tbdrmes[{Na,la,jja,N,l,jj,Nc,lc,jjc,Nd,ld,jjd,Ja,Jcd,J0}][1]
			       +vmepn*tbdrmes[{Na,la,jja,N,l,jj,Nc,lc,jjc,Nd,ld,jjd,Ja,Jcd,J0}][2]);
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

 double pkere2=0.0;
 factor1=sqrt(double((ssp+1)*(ss+1)));
 for(int jjp=std::abs(2*lp-1); jjp<=2*lp+1; jjp=jjp+2){
  if(std::abs(jjp-II1p)>JJ || JJ>jjp+II1p)continue;
  double factor2=factor1*sqrt(double(jjp+1))*wig6jj(II1p,1,ssp,2*lp,JJ,jjp);
  for(int jj=std::abs(2*l-1); jj<=2*l+1; jj=jj+2){
   if(std::abs(jj-II1)>JJ || JJ>jj+II1)continue;
   double factor3=factor2*sqrt(double(jj+1))*wig6jj(II1,1,ss,2*l,JJ,jj);
   for(int Na=0; Na<=Nmax; Na++){
    for(int Nb=0; Nb<=std::min(Nmax,N0+Nmax-1-Na); Nb++){
     int Ndmin=std::max(0,Na+Nb-Np-Nmax+1); // Na+Nb-Np-Nd <= Nmax-1
     if((Na+Nb-Np-Ndmin)%2!=0)Ndmin++;
     for(int Nd=Ndmin; Nd<=std::min(std::min(Nmax,N0+Nmax-1-Np),Na+Nb-Np+Nmax-1); Nd+=2){
      for(int la=Na%2; la<=Na; la+=2){
       for(int lb=Nb%2; lb<=Nb; lb+=2){
        for(int ld=Nd%2; ld<=Nd; ld+=2){
         for(int jja=std::abs(2*la-1); jja<=2*la+1; jja+=2){
	  for(int jjb=std::abs(2*lb-1); jjb<=2*lb+1; jjb+=2){
           for(int jjd=std::abs(2*ld-1); jjd<=2*ld+1; jjd+=2){
            double d=1.0;
            if(Na==Nb && la==lb && jja==jjb)d=d/sqrt(2.0);
            if(Np==Nd && lp==ld && jjp==jjd)d=d/sqrt(2.0);
            for(int Jab=std::abs(jja-jjb)/2; Jab<=(jja+jjb)/2; Jab++){
	     if(std::abs(jjd-jj)>2*Jab || 2*Jab>jjd+jj)continue;
	     for(int Jd=std::abs(jjp-jjd)/2; Jd<=(jjp+jjd)/2; Jd++){
	      for(int J0=std::abs(Jab-Jd); J0<=Jab+Jd; J0++){
	       if(std::abs(II1-2*J0)>II1p || II1p>II1+2*J0 || std::abs(2*J0-jjp)>jj || jj>2*J0+jjp)continue;
               double factor4=factor3*(2*J0+1)*sqrt(double((2*Jab+1)*(2*Jd+1)))*wig6jj(II1,2*J0,II1p,jjp,JJ,jj)
		                                                               *wig6jj(jjd,2*Jd,jjp,2*J0,jj,2*Jab);
	       int phase=(II1p+jjp+JJ)/2+Jd+Jab+J0;
	       if(phase%2!=0)factor4=-factor4;
	       double vmen=like_int(Na,la,jja,Nb,lb,jjb,Nd,ld,jjd,N,l,jj,Jab,vn);
	       double vmepn=pn_int(Na,la,jja,Nb,lb,jjb,Nd,ld,jjd,N,l,jj,Jab,vpn);
               pkere2+=factor4*(0.5*vmen*tbdrmes[{Na,la,jja,Nb,lb,jjb,Np,lp,jjp,Nd,ld,jjd,Jab,Jd,J0}][1]
		 	        +vmepn*tbdrmes[{Na,la,jja,Nb,lb,jjb,Np,lp,jjp,Nd,ld,jjd,Jab,Jd,J0}][2]);
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
 
 kernel_file<<ssp<<" "<<lp<<" "<<Np<<" "<<ss<<" "<<l<<" "<<N<<" "<<nker<<" "<<pkerd<<" "<<pkere<<" "<<pkere2<<std::endl;
}
}
}
}
}
}
}
}

return EXIT_SUCCESS;
}

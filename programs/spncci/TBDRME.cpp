#include <cstdio>
#include <fstream>
#include <istream>
#include <iostream>
#include <omp.h>

#include "spncci/computation_control.h"
#include "wigxjpf.h"

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
std::map<std::array<int,13>,double> amplitude_by_labels;
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
    amplitude_by_labels[{gamma,Nex_sigma,lm_sigma,mu_sigma,SSp,SSn,SS,upsilon,Nex_omega,lm_omega,mu_omega,kappa,L}]=amplitude;
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
std::map<std::array<int,15>,std::array<double,3>> tbdrmes;
// Key is {N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJf,JJi,JJ0}
// Value is {proton, neutron, proton-neutron}
int etamax=std::max(N1vp,N1vn)+Nmax;
int Ntotmax=N0+Nmax;
for(int N1=0; N1<=etamax; N1++){
for(int N2=0; N2<=std::min(etamax,Ntotmax-N1); N2++){
for(int l1=N1%2; l1<=N1; l1+=2){
for(int l2=N2%2; l2<=N2; l2+=2){
for(int jj1=std::abs(2*l1-1); jj1<=2*l1+1; jj1+=2){
for(int jj2=std::abs(2*l2-1); jj2<=2*l2+1; jj2+=2){
for(int N3=0; N3<=etamax; N3++){
for(int l3=N3%2; l3<=N3; l3+=2){
for(int jj3=std::abs(2*l3-1); jj3<=2*l3+1; jj3+=2){
int N4min=std::max(N1+N2-N3-Nmax,0);
if((N1+N2-N3-N4min)%2!=0)N4min++; // so that TBD operator doesn't change parity
for(int N4=N4min; N4<=std::min(std::min(etamax,Ntotmax-N3),N1+N2-N3+Nmax); N4+=2){
for(int l4=N4%2; l4<=N4; l4+=2){
for(int jj4=std::abs(2*l4-1); jj4<=2*l4+1; jj4+=2){
for(int JJf=std::abs(jj1-jj2); JJf<=jj1+jj2; JJf+=2){
for(int JJi=std::abs(jj3-jj4); JJi<=jj3+jj4; JJi+=2){
for(int JJ0=std::max(std::abs(JJf-JJi),std::abs(JJ_ket-JJ_bra)); JJ0<=std::min(JJf+JJi,JJ_ket+JJ_bra); JJ0+=2){
  std::array<double,3> value={0.0,0.0,0.0};
  tbdrmes[{N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJf,JJi,JJ0}]=value;
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
std::cout<<"Number of TBDRMEs: "<<tbdrmes.size()<<std::endl;
std::vector<std::array<int,8>> N1N2l1l2xfLfkf;
for(int N1=0; N1<=etamax; N1++){
for(int N2=0; N2<=std::min(etamax,Ntotmax-N1); N2++){
for(int l1=N1%2; l1<=N1; l1+=2){
for(int l2=N2%2; l2<=N2; l2+=2){
for(auto xf : u3::KroneckerProduct(u3::SU3(N1,0),u3::SU3(N2,0))){
int lmf=xf.irrep.lambda();
int muf=xf.irrep.mu();
for(auto Lftagged : u3::BranchingSO3(xf.irrep)){
int Lf=int(Lftagged.irrep);
if(std::abs(l1-l2)>Lf || Lf>l1+l2)continue;
for(int kappaf=1; kappaf<=Lftagged.tag; kappaf++){
  N1N2l1l2xfLfkf.push_back({N1,N2,l1,l2,lmf,muf,Lf,kappaf});
}
}
}
}
}
}
}
std::cout<<"Number of iterations of outer loop: "<<N1N2l1l2xfLfkf.size()<<std::endl;

#pragma omp parallel
{
//int num_threads=omp_get_num_threads();
//std::cout<<"Number of threads: "<<num_threads<<std::endl;
spncci::InitializeSpNCCI();
u3::PhiCoefCache phi_coef_cache;
#pragma omp for schedule(dynamic)
for(int ind=0; ind<N1N2l1l2xfLfkf.size(); ind++){
//for(int N1=0; N1<=etamax; N1++){
// for(int N2=0; N2<=std::min(etamax,Ntotmax-N1); N2++){
//  for(int l1=N1%2; l1<=N1; l1+=2){
//   for(int l2=N2%2; l2<=N2; l2+=2){
//    for(auto xf : u3::KroneckerProduct(u3::SU3(N1,0),u3::SU3(N2,0))){
//     for(auto Lftagged : u3::BranchingSO3(xf)){
//      int Lf=int(Lftagged.irrep);
//      if(std::abs(l1-l2)>Lf || Lf>l1+l2)continue;
//      for(int kappaf=1; kappaf<=Lftagged.tag; kappaf++){
 int N1=N1N2l1l2xfLfkf[ind][0];
 int N2=N1N2l1l2xfLfkf[ind][1];
 int l1=N1N2l1l2xfLfkf[ind][2];
 int l2=N1N2l1l2xfLfkf[ind][3];
 int lmf=N1N2l1l2xfLfkf[ind][4];
 int muf=N1N2l1l2xfLfkf[ind][5];
 int Lf=N1N2l1l2xfLfkf[ind][6];
 int kappaf=N1N2l1l2xfLfkf[ind][7];
 u3::SU3 xf(lmf,muf);
     int phase1=lmf+muf+N1+N2;
       double factor1=sqrt(double((JJ_ket+1)*(2*Lf+1)))*u3::W(u3::SU3(N1,0),1,l1,u3::SU3(N2,0),1,l2,xf,kappaf,Lf,1);
       for(int N3=0; N3<=etamax; N3++){
        int N4min=std::max(N1+N2-N3-Nmax,0);
        if((N1+N2-N3-N4min)%2!=0)N4min++; // so that TBD operator doesn't change parity
        for(int N4=N4min; N4<=std::min(std::min(etamax,Ntotmax-N3),N1+N2-N3+Nmax); N4+=2){
         for(int l3=N3%2; l3<=N3; l3+=2){
          for(int l4=N4%2; l4<=N4; l4+=2){
           for(auto xi : u3::KroneckerProduct(u3::SU3(0,N3),u3::SU3(0,N4))){
   	    int phase2=phase1+xi.irrep.lambda()+xi.irrep.mu()+N3+N4;
            for(auto Litagged : u3::BranchingSO3(xi.irrep)){
             int Li=int(Litagged.irrep);
             if(std::abs(l3-l4)>Li || Li>l3+l4)continue;
              for(int kappai=1; kappai<=Litagged.tag; kappai++){
              double factor2=factor1*sqrt(double(2*Li+1))*u3::W(u3::SU3(0,N3),1,l3,u3::SU3(0,N4),1,l4,xi.irrep,kappai,Li,1);
              for(auto x0 : u3::KroneckerProduct(xf,xi.irrep)){
	       int phase3=phase2+x0.irrep.lambda()+x0.irrep.mu()+x0.tag;
               for(auto L0tagged : u3::BranchingSO3(x0.irrep)){
	        int L0=int(L0tagged.irrep);
	        if(std::abs(Lf-Li)>L0 || L0>Lf+Li)continue;
	        for(int kappa0=1; kappa0<=L0tagged.tag; kappa0++){
                 for(int rho0=1; rho0<=x0.tag; rho0++){
                  double factor3=factor2*sqrt(double(2*L0+1))*u3::W(xf,kappaf,Lf,xi.irrep,kappai,Li,x0.irrep,kappa0,L0,rho0);
                  for(std::map<std::array<int,13>,double>::iterator
                      itbra=amplitude_by_labels.begin(); itbra!=amplitude_by_labels.end(); itbra++){
                   int gamma_bra=itbra->first[0];
                   int Nex_sigma_bra=itbra->first[1]; 
                   int lm_sigma_bra=itbra->first[2];
                   int mu_sigma_bra=itbra->first[3];
                   int SSp_bra=itbra->first[4];
                   int SSn_bra=itbra->first[5];
                   int SS_bra=itbra->first[6];
                   int upsilon_bra=itbra->first[7];
                   int Nex_omega_bra=itbra->first[8];
                   int lm_omega_bra=itbra->first[9];
                   int mu_omega_bra=itbra->first[10];
                   int kappa_bra=itbra->first[11];
                   int L_bra=itbra->first[12];
	           int Nn_bra=Nex_omega_bra-Nex_sigma_bra;
                   u3::SU3 x_bra(lm_omega_bra,mu_omega_bra);
	           double factor4=factor3*sqrt(double((2*L_bra+1)*(SS_bra+1)))*itbra->second;
                   for(std::map<std::array<int,13>,double>::iterator
                       itket=amplitude_by_labels.begin(); itket!=amplitude_by_labels.end(); itket++){
                    int L_ket=itket->first[12];
                    if(std::abs(L_ket-L0)>L_bra || L_bra>L_ket+L0)continue;
                    int gamma_ket=itket->first[0];
                    int Nex_sigma_ket=itket->first[1];
                    int lm_sigma_ket=itket->first[2];
                    int mu_sigma_ket=itket->first[3];
                    int SSp_ket=itket->first[4];
                    int SSn_ket=itket->first[5];
                    int SS_ket=itket->first[6];
                    int upsilon_ket=itket->first[7];
                    int Nex_omega_ket=itket->first[8];
                    int lm_omega_ket=itket->first[9];
                    int mu_omega_ket=itket->first[10];
                    int kappa_ket=itket->first[11];
		    int Nn_ket=Nex_omega_ket-Nex_sigma_ket;
                    u3::SU3 x_ket(lm_omega_ket,mu_omega_ket);
		    double factor11=sqrt(double((SS_ket+1)*u3::dim(x_ket))/double((SS_bra+1)*u3::dim(x_bra)));
		    int phase=phase3+lm_omega_bra+mu_omega_bra-lm_omega_ket-mu_omega_ket+(SS_ket-SS_bra)/2;
	            double factor5=factor4*itket->second;
                    for(int rho=1; rho<=u3::OuterMultiplicity(x_ket,x0.irrep,x_bra); rho++){
                     double factor6=factor5*u3::W(x_ket,kappa_ket,L_ket,x0.irrep,kappa0,L0,x_bra,kappa_bra,L_bra,rho);
                     for(int jj1=std::abs(2*l1-1); jj1<=2*l1+1; jj1+=2){
                      for(int jj2=std::abs(2*l2-1); jj2<=2*l2+1; jj2+=2){
                       for(int JJf=std::abs(jj1-jj2); JJf<=jj1+jj2; JJf+=2){
	                for(int SSf=0; SSf<=2; SSf+=2){
                         if(std::abs(2*Lf-SSf)>JJf || JJf>2*Lf+SSf)continue;
                         double factor7=factor6*sqrt(double((jj1+1)*(jj2+1)*(JJf+1)*(SSf+1)))*wig9jj(2*l1,2*l2,2*Lf,1,1,SSf,jj1,jj2,JJf);
                         for(int jj3=std::abs(2*l3-1); jj3<=2*l3+1; jj3+=2){
                          for(int jj4=std::abs(2*l4-1); jj4<=2*l4+1; jj4+=2){
                           for(int JJi=std::abs(jj3-jj4); JJi<=jj3+jj4; JJi+=2){
                            for(int SSi=0; SSi<=2; SSi+=2){
                             if(std::abs(2*Li-SSi)>JJi || JJi>2*Li+SSi)continue;
                             double factor8=factor7*sqrt(double((jj3+1)*(jj4+1)*(JJi+1)*(SSi+1)))*wig9jj(2*l3,2*l4,2*Li,1,1,SSi,jj3,jj4,JJi);
                             for(int JJ0=std::max(std::abs(JJf-JJi),std::abs(JJ_ket-JJ_bra)); JJ0<=std::min(JJf+JJi,JJ_ket+JJ_bra); JJ0+=2){
			      std::array<int,15> keyt={N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,JJf,JJi,JJ0};
                              for(int SS0=std::abs(SSf-SSi); SS0<=SSf+SSi; SS0+=2){
	 	               if(std::abs(2*L0-SS0)>JJ0 || JJ0>2*L0+SS0 || std::abs(SS_ket-SS0)>SS_bra || SS_bra>SS_ket+SS0)continue;
                               double factor9=factor8*sqrt(double((JJ0+1)*(SS0+1)))*wig9jj(2*Lf,2*Li,2*L0,SSf,SSi,SS0,JJf,JJi,JJ0)
		   	                      *wig9jj(2*L_ket,SS_ket,JJ_ket,2*L0,SS0,JJ0,2*L_bra,SS_bra,JJ_bra);

                               double TBDRME_p;
 		               double TBDRME_n;
		               double TBDRME_pn;
                               if(Nn_bra>=Nn_ket){
                                std::array<int,37> key={gamma_bra,Nex_sigma_bra,lm_sigma_bra,mu_sigma_bra,SSp_bra,SSn_bra,
                                                        SS_bra,upsilon_bra,Nex_omega_bra,lm_omega_bra,mu_omega_bra,
                                                        gamma_ket,Nex_sigma_ket,lm_sigma_ket,mu_sigma_ket,SSp_ket,SSn_ket,
                                                        SS_ket,upsilon_ket,Nex_omega_ket,lm_omega_ket,mu_omega_ket,
                                                        N1,N2,N3,N4,lmf,muf,SSf,xi.irrep.lambda(),xi.irrep.mu(),SSi,
		                                        x0.irrep.lambda(),x0.irrep.mu(),SS0,rho0,rho};
                                TBDRME_p=proton_TBDRME_by_labels[key];
                                TBDRME_n=neutron_TBDRME_by_labels[key];
                                TBDRME_pn=pn_TBDRME_by_labels[key];
                               }else{
		                TBDRME_p=0.0;
                                TBDRME_n=0.0;
                                TBDRME_pn=0.0;
		                int sign=1;
                                for(int rho0p=1; rho0p<=x0.tag; rho0p++){
		                 sign=-sign;
                                 std::array<int,37> key={gamma_ket,Nex_sigma_ket,lm_sigma_ket,mu_sigma_ket,SSp_ket,SSn_ket,
                                                         SS_ket,upsilon_ket,Nex_omega_ket,lm_omega_ket,mu_omega_ket,
		                                         gamma_bra,Nex_sigma_bra,lm_sigma_bra,mu_sigma_bra,SSp_bra,SSn_bra,
                                                         SS_bra,upsilon_bra,Nex_omega_bra,lm_omega_bra,mu_omega_bra,
                                                         N4,N3,N2,N1,xi.irrep.mu(),xi.irrep.lambda(),SSi,muf,lmf,SSf,
                                                         x0.irrep.mu(),x0.irrep.lambda(),SS0,rho0p,rho};
    	                         double factor0=sign*PhiCached(phi_coef_cache,xf,xi.irrep,x0.irrep,rho0,rho0p);
                                 TBDRME_p+=factor0*proton_TBDRME_by_labels[key];
                                 TBDRME_n+=factor0*neutron_TBDRME_by_labels[key];
	                         TBDRME_pn+=factor0*pn_TBDRME_by_labels[key];
                                }
                                TBDRME_p*=factor11;
		                TBDRME_n*=factor11;
		                TBDRME_pn*=factor11;
		                if(phase%2!=0){
		                 TBDRME_p=-TBDRME_p;
		                 TBDRME_n=-TBDRME_n;
		                 TBDRME_pn=-TBDRME_pn;
		                }
		               }

		               tbdrmes[keyt][0]+=factor9*TBDRME_p;
                               tbdrmes[keyt][1]+=factor9*TBDRME_n;
                               tbdrmes[keyt][2]+=factor9*TBDRME_pn;
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
//      }
//     }
//    }
//   }
//  }
// }
}
phi_coef_cache.clear();
} // End of parallel region

//**********************************************
// Output TBDRMEs
//**********************************************
std::ofstream tbdfile("tbdrmes.dat");
if(!tbdfile){
  std::cout<<"Could not open tbdrmes.dat file"<<std::endl;
  return EXIT_FAILURE;
} 
tbdfile<<"N1 l1 2j1 N2 l2 2j2 N3 l3 2j3 N4 l4 2j4 2Jf 2Ji 2J0 proton neutron proton-neutron"<<std::endl;
for(std::map<std::array<int,15>,std::array<double,3>>::iterator it=tbdrmes.begin(); it!=tbdrmes.end(); it++){
 double trdensfactor=sqrt(double(JJ_bra+1)/double(it->first[14]+1));
 tbdfile<<it->first[0]<<" "<<it->first[1]<<" "<<it->first[2]<<" "<<it->first[3]<<" "<<it->first[4]<<" "<<it->first[5]<<" "
        <<it->first[6]<<" "<<it->first[7]<<" "<<it->first[8]<<" "<<it->first[9]<<" "<<it->first[10]<<" "<<it->first[11]<<" "
        <<it->first[12]<<" "<<it->first[13]<<" "<<it->first[14]<<" "
        <<it->second[0]*trdensfactor<<" "<<it->second[1]*trdensfactor<<" "<<it->second[2]*trdensfactor<<std::endl;
}
tbdfile.close();

return EXIT_SUCCESS;
}


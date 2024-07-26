#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <map>

#include <Eigen/Dense>
#include <su3.h>
#include "wigxjpf.h"

class Spins {
  public:
  // default constructor
  inline Spins() : SSp_(0), SSn_(0), SS_(0) {}
  // constructor
  inline Spins(int SSp, int SSn, int SS) : SSn_(SSn), SSp_(SSp), SS_(SS) {}
  // accessors
  inline int SSp() const {
    return SSp_;
  }
  inline int SSn() const {
    return SSn_;
  }
  inline int SS() const {
    return SS_;
  }
  // key tuple, comparisons
  typedef std::tuple<int,int,int> KeyType;
  inline KeyType Key() const{
    return KeyType(SSp(), SSn(), SS());
  }
  inline friend bool operator == (const Spins& SpSnS1, const Spins& SpSnS2){
    return SpSnS1.Key() == SpSnS2.Key();
  }
  inline friend bool operator < (const Spins& SpSnS1, const Spins& SpSnS2){
    return SpSnS1.Key() < SpSnS2.Key();
  }
  //labels
  private:
  int SSp_, SSn_, SS_;
};

class SU3 {
  public:
  // default constructor
  inline SU3() : lm_(0), mu_(0) {}
  // constructor
  inline SU3(int lm, int mu) : lm_(lm), mu_(mu) {}
  // accessors
  inline int lm() const {
    return lm_;
  }
  inline int mu() const {
    return mu_;
  }
  // key tuple, comparisons
  typedef std::tuple<int,int> KeyType;
  inline KeyType Key() const{
    return KeyType(lm(), mu());
  }
  inline friend bool operator == (const SU3& x1, const SU3& x2){
    return x1.Key() == x2.Key();
  }
  inline friend bool operator < (const SU3& x1, const SU3& x2){
    return x1.Key() < x2.Key();
  }
  //labels
  private:
  int lm_, mu_;
};
/*
class SpinsNex {
  // default constructor
  inline SpinsNex() : Nex_(0) {}
  // constructor
  inline SpinsNex(const Spins& SpSnS, const int& Nex) : SpSnS_(SpSnS), Nex_(Nex) {}
  // accessors
  inline Spins SpSnS() const
  {
    return SpSnS_;
  }
  inline int Nex() const
  {
    return Nex_;
  }
  // key tuple and comparisons
  typedef std::tuple<Spins,int> KeyType;
  inline KeyType Key() const
  {
    return KeyType(SpSnS_,Nex_);
  }
  inline friend bool operator == (const SpinsNex& SpSnSNex1, const SpinsNex& SpSnSNex2)
  {
    return SpSnSNex1.Key() == SpSnSNex2.Key();
  }
  inline friend bool operator < (const SpinsNex& SpSnSNex1, const SpinsNex& SpSnSNex2)
  {
    return SpSnSNex1.Key() < SpSnSNex2.Key();
  }
  // labels
  private:
  Spins SpSnS_;
  int Nex_;
}; 
*/
class Ipini {
  public:
  // default constructor
  inline Ipini() : ip_(0), in_(0), i_(0) {}
  // constructor
  inline Ipini(int ip, int in, int i) : ip_(ip), in_(in), i_(i) {}
  // accessors
  inline int ip() const {
    return ip_;
  }
  inline int in() const {
    return in_;
  }
  inline int i() const {
    return i_;
  }
  // key tuple, comparisons
  typedef std::tuple<int,int,int> KeyType;
  inline KeyType Key() const{
    return KeyType(ip(), in(), i());
  }
  inline friend bool operator == (const Ipini& ipini1, const Ipini& ipini2){
    return ipini1.Key() == ipini2.Key();
  }
  inline friend bool operator < (const Ipini& ipini1, const Ipini& ipini2){
    return ipini1.Key() < ipini2.Key();
  }
  //labels
  private:
  int ip_, in_, i_;
};

class LSU3shellBasisState {
  public:
  // default constructor
  inline LSU3shellBasisState() : Nex_(0) {}
  // constructor
  inline LSU3shellBasisState(const Ipini& ipini, const Spins& SpSnS, const int& Nex, const SU3& x) :
	  ipini_(ipini), SpSnS_(SpSnS), Nex_(Nex), x_(x) {}
  // accessors
  inline Ipini ipini() const
  {
    return ipini_;
  }
  inline Spins SpSnS() const
  {
    return SpSnS_;
  }
  inline int Nex() const
  {
    return Nex_;
  }
  inline SU3 x() const
  {
    return x_;
  }
  // key tuple and comparisons
  typedef std::tuple<Ipini,Spins,int,SU3> KeyType;
  inline KeyType Key() const
  {
    return KeyType(ipini_,SpSnS_,Nex_,x_);
  }
  inline friend bool operator == (const LSU3shellBasisState& bs1, const LSU3shellBasisState& bs2)
  {
    return bs1.Key() == bs2.Key();
  }
  inline friend bool operator < (const LSU3shellBasisState& bs1, const LSU3shellBasisState& bs2)
  {
    return bs1.Key() < bs2.Key();
  }
  // labels
  private:
  Ipini ipini_;
  Spins SpSnS_;
  int Nex_;
  SU3 x_;
};

class ONBasisState {
  public:
  // default constructor
  inline ONBasisState() : gamma_(0), Nex_sigma_(0), upsilon_(0), Nex_omega_(0) {}
  // constructor
  inline ONBasisState(const int& gamma, const int& Nex_sigma, const SU3& x_sigma, const Spins& SpSnS,
		      const int& upsilon, const int& Nex_omega, const SU3& x_omega) :
	  gamma_(gamma), Nex_sigma_(Nex_sigma), x_sigma_(x_sigma), SpSnS_(SpSnS),
	  upsilon_(upsilon), Nex_omega_(Nex_omega), x_omega_(x_omega) {}
  // accessors
  inline int gamma() const
  {
    return gamma_;
  }
  inline int Nex_sigma() const
  {
    return Nex_sigma_;
  }
  inline SU3 x_sigma() const
  {
    return x_sigma_;
  }
  inline Spins SpSnS() const
  {
    return SpSnS_;
  }
  inline int upsilon() const
  {
    return upsilon_;
  }
  inline int Nex_omega() const
  {
    return Nex_omega_;
  }
  inline SU3 x_omega() const
  {
    return x_omega_;
  }
  // key tuple and comparisons
  typedef std::tuple<int,int,SU3,Spins,int,int,SU3> KeyType;
  inline KeyType Key() const
  {
    return KeyType(gamma_,Nex_sigma_,x_sigma_,SpSnS_,upsilon_,Nex_omega_,x_omega_);
  }
  inline friend bool operator == (const ONBasisState& bs1, const ONBasisState& bs2)
  {
    return bs1.Key() == bs2.Key();
  }
  inline friend bool operator < (const ONBasisState& bs1, const ONBasisState& bs2)
  {
    return bs1.Key() < bs2.Key();
  }
  // labels
  private:
  int gamma_,Nex_sigma_,upsilon_,Nex_omega_;
  SU3 x_sigma_,x_omega_;
  Spins SpSnS_;
};

std::array<std::vector<std::array<int,5>>,2> ReadListOfTensors(const std::string& filename){
  std::array<std::vector<std::array<int,5>>,2> tensors;
  std::ifstream file(filename);
  if(!file){
    std::cerr << "Could not open '" << filename << "' input file!" << std::endl;
    su3::finalize();
    exit(EXIT_FAILURE);
  }
  while(true){
    int na, nb, lm0, mu0, ss0;
    file >> na >> nb >> lm0 >> mu0 >> ss0;
    if(file){
      if(((na-nb)/2)*2==na-nb){
        tensors[0].push_back({na, nb, lm0, mu0, ss0});
      }else{
        tensors[1].push_back({na, nb, lm0, mu0, ss0});
      }
    }else{
      break;
    }
  }
  return tensors;
}

void ParseInputFile(const std::stringstream& filename, const std::array<int,5>& tensor,
	   std::map<Spins,std::map<int,std::vector<SU3>>>& SpSnSNexGamma_set,
	   std::map<std::tuple<Spins,int,SU3>,std::vector<Ipini>>& lsu3shell_basis,
	   std::map<std::tuple<LSU3shellBasisState,LSU3shellBasisState,std::array<int,5>>,std::vector<double>>& obdrmes_lsu3shell){
  std::ifstream input_file(filename.str());
  if(!input_file){
    std::cerr << "Could not open file '" << filename.str() << std::endl;
    su3::finalize();
    exit(EXIT_FAILURE);
  }

  while (true) {
    int ip,in,i,SSp_bra,SSn_bra,SS_bra,N_bra,lm_bra,mu_bra,jp,jn,j,SSp_ket,SSn_ket,SS_ket,N_ket,lm_ket,mu_ket;
    input_file>>ip>>in>>i>>SSp_bra>>SSn_bra>>SS_bra>>N_bra>>lm_bra>>mu_bra
              >>jp>>jn>>j>>SSp_ket>>SSn_ket>>SS_ket>>N_ket>>lm_ket>>mu_ket;
    if(input_file){
      Ipini ipini_bra(ip,in,i);
      Spins SpSnS_bra(SSp_bra,SSn_bra,SS_bra);
      SU3 x_bra(lm_bra,mu_bra);
      LSU3shellBasisState bra(ipini_bra,SpSnS_bra,N_bra,x_bra);
      Ipini ipini_ket(jp,jn,j);
      Spins SpSnS_ket(SSp_ket,SSn_ket,SS_ket);
      SU3 x_ket(lm_ket,mu_ket);
      LSU3shellBasisState ket(ipini_ket,SpSnS_ket,N_ket,x_ket);
      int rhomax=su3::mult(lm_ket,mu_ket,tensor[2],tensor[3],lm_bra,mu_bra);
      std::vector<double> rmes(rhomax);
      for(int irho=0; irho<rhomax; irho++){
        input_file>>rmes[irho];
      }
      obdrmes_lsu3shell[{bra,ket,tensor}]=rmes;

      std::map<Spins,std::map<int,std::vector<SU3>>>::iterator it=SpSnSNexGamma_set.find(SpSnS_bra);
      if(it==SpSnSNexGamma_set.end()){
        std::map<int,std::vector<SU3>> x_by_Nex;
        std::vector<SU3> x_set;
	x_set.push_back(x_bra);
	x_by_Nex[N_bra]=x_set;
	SpSnSNexGamma_set[SpSnS_bra]=x_by_Nex;
      }else{
        std::map<int,std::vector<SU3>>::iterator it2=SpSnSNexGamma_set[SpSnS_bra].find(N_bra);
	if(it2==SpSnSNexGamma_set[SpSnS_bra].end()){
          std::vector<SU3> x_set;
          x_set.push_back(x_bra);
	  SpSnSNexGamma_set[SpSnS_bra][N_bra]=x_set;
	}else{
          auto it3=std::find(SpSnSNexGamma_set[SpSnS_bra][N_bra].begin(), SpSnSNexGamma_set[SpSnS_bra][N_bra].end(), x_bra);
          if(it3==SpSnSNexGamma_set[SpSnS_bra][N_bra].end()){
            SpSnSNexGamma_set[SpSnS_bra][N_bra].push_back(x_bra);
          }
	}
      }

      it=SpSnSNexGamma_set.find(SpSnS_ket);
      if(it==SpSnSNexGamma_set.end()){
        std::map<int,std::vector<SU3>> x_by_Nex;
        std::vector<SU3> x_set;
	x_set.push_back(x_ket);
	x_by_Nex[N_ket]=x_set;
	SpSnSNexGamma_set[SpSnS_ket]=x_by_Nex;
      }else{
        std::map<int,std::vector<SU3>>::iterator it2=SpSnSNexGamma_set[SpSnS_ket].find(N_ket);
	if(it2==SpSnSNexGamma_set[SpSnS_ket].end()){
          std::vector<SU3> x_set;
          x_set.push_back(x_ket);
	  SpSnSNexGamma_set[SpSnS_ket][N_ket]=x_set;
	}else{
          auto it3=std::find(SpSnSNexGamma_set[SpSnS_ket][N_ket].begin(), SpSnSNexGamma_set[SpSnS_ket][N_ket].end(), x_ket);
          if(it3==SpSnSNexGamma_set[SpSnS_ket][N_ket].end())SpSnSNexGamma_set[SpSnS_ket][N_ket].push_back(x_ket);
	}
      }

      std::map<std::tuple<Spins,int,SU3>,std::vector<Ipini>>::iterator iter=lsu3shell_basis.find({SpSnS_bra,N_bra,x_bra});
      if(iter==lsu3shell_basis.end()){
        std::vector<Ipini> ipini_set;
	ipini_set.push_back(ipini_bra);
	lsu3shell_basis[{SpSnS_bra,N_bra,x_bra}]=ipini_set;
      }else{
	auto iter2=std::find(lsu3shell_basis[{SpSnS_bra,N_bra,x_bra}].begin(), lsu3shell_basis[{SpSnS_bra,N_bra,x_bra}].end(), ipini_bra);
        if(iter2==lsu3shell_basis[{SpSnS_bra,N_bra,x_bra}].end())lsu3shell_basis[{SpSnS_bra,N_bra,x_bra}].push_back(ipini_bra);
      }

      iter=lsu3shell_basis.find({SpSnS_ket,N_ket,x_ket});
      if(iter==lsu3shell_basis.end()){
        std::vector<Ipini> ipini_set;
        ipini_set.push_back(ipini_ket);
        lsu3shell_basis[{SpSnS_ket,N_ket,x_ket}]=ipini_set;
      }else{
        auto iter2=std::find(lsu3shell_basis[{SpSnS_ket,N_ket,x_ket}].begin(), lsu3shell_basis[{SpSnS_ket,N_ket,x_ket}].end(), ipini_ket);
        if(iter2==lsu3shell_basis[{SpSnS_ket,N_ket,x_ket}].end())lsu3shell_basis[{SpSnS_ket,N_ket,x_ket}].push_back(ipini_ket);
      }

    }else{
      break;
    }

  }
}

std::vector<std::tuple<SU3,int>> KroneckerProduct(const int lm1, const int mu1, const int lm2, const int mu2){
  std::vector<std::tuple<SU3,int>> product;
  for(int lm=0; lm<=lm1+lm2+std::min(mu2,lm1+mu1); lm++){
    for(int mu=0; mu<=mu1+mu2+std::min(lm1,lm2); mu++){
      int rhomax=su3::mult(lm1,mu1,lm2,mu2,lm,mu);
      if(rhomax!=0){
        SU3 x(lm,mu);
        product.push_back({x,rhomax});
      }
    }
  }
  return product;
}

double U(const int lm1, const int mu1, const int lm2, const int mu2, const int lm, const int mu, const int lm3, const int mu3,
         const int lm12, const int mu12, const int rho12, const int rho12_3, const int lm23, const int mu23, const int rho23,
	 const int rho1_23){
  int rho12max = su3::mult(lm1, mu1, lm2, mu2, lm12, mu12);
  int rho12_3max = su3::mult(lm12, mu12, lm3, mu3, lm, mu);
  int rho23max = su3::mult(lm2, mu2, lm3, mu3, lm23, mu23);
  int rho1_23max = su3::mult(lm1, mu1, lm23, mu23, lm, mu);
  if(rho12max==0 || rho12_3max==0 || rho23max==0 || rho1_23max==0)return 0.0;
  std::vector<double> dru3;
  su3::wru3(lm1, mu1, lm2, mu2, lm, mu, lm3, mu3, lm12, mu12, lm23, mu23, rho12max, rho12_3max, rho23max, rho1_23max, dru3);
  int nb = rho12_3max*rho12max;
  int nc = rho23max*nb;
  int i = rho12-1 + (rho12_3-1)*rho12max + (rho23-1)*nb + (rho1_23-1)*nc;
  return dru3[i];
}

double Phi(const int rho, const int rhop, const int lm1, const int mu1, const int lm2, const int mu2, const int lm3, const int mu3){
  int rhomax = su3::mult(lm1, mu1, lm2, mu2, lm3, mu3);
  if(rhomax==1){
    int phase=lm1+mu1+lm2+mu2-lm3-mu3;
    if((phase/2)*2!=phase){
      return -1.0;
    }else{
      return 1.0;
    }
  }else{
    std::vector<double> dzu3;
    su3::wzu3(lm1, mu1, 0, 0, lm3, mu3, lm2, mu2, lm1, mu1, lm2, mu2, 1, rhomax, 1, rhomax, dzu3);
    int i = rho-1 + (rhop-1)*rhomax;
    return dzu3[i];
  }
}

double BosonCreationOperatorME(const int n1p, const int n2p, const int n3p, const int n1, const int n2, const int n3){
  double me;
  if(n1p==n1+2 && n2p==n2 && n3p==n3){
    me=sqrt(double((n1+4)*(n1-n2+2)*(n1-n3+3))/double(2*(n1-n2+3)*(n1-n3+4)));
  }else if(n1p==n1 && n2p==n2+2 && n3p==n3){
    me=sqrt(double((n2+3)*(n1-n2)*(n2-n3+2))/double(2*(n1-n2-1)*(n2-n3+3)));
  }else if(n1p==n1 && n2p==n2 && n3p==n3+2){ 
    me=sqrt(double((n3+2)*(n2-n3)*(n1-n3+1))/double(2*(n1-n3)*(n2-n3-1)));
  }else{
    me=0.0;
  }
  return me;
}

int su3dim(const int lm, const int mu){
  return ((lm+1)*(mu+1)*(lm+mu+2))/2;
}

int su3dim(const SU3& x){
  return ((x.lm()+1)*(x.mu()+1)*(x.lm()+x.mu()+2))/2;
}

double ARME(const int gamma, const int Nex_sigma, const SU3& x_sigma, const Spins& SpSnS, const int upsilonp, const SU3& x_omegap,
	    const int upsilon, const int Nex, const SU3& x_omega,
            std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>& orthonormal_basis,
	    std::map<std::tuple<Spins,int>,std::map<std::array<SU3,2>,Eigen::MatrixXd>>& A){
  return orthonormal_basis[{gamma,Nex_sigma,x_sigma,SpSnS}][Nex+2][x_omegap][upsilonp].transpose()
	 *A[{SpSnS,Nex}][{x_omegap,x_omega}]
	 *orthonormal_basis[{gamma,Nex_sigma,x_sigma,SpSnS}][Nex][x_omega][upsilon];
}

std::vector<double> Recurrence(const ONBasisState& bra, const ONBasisState& ket, const std::array<int, 5>& tensor,
                    std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>& obdrmes,
                    std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>& orthonormal_basis,
                    std::map<std::tuple<Spins,int>,std::map<std::array<SU3,2>,Eigen::MatrixXd>>& A,
                    std::map<std::tuple<int,int,SU3,Spins,int,SU3>,Eigen::MatrixXd>& K,
                    std::map<std::tuple<int,int,SU3,Spins,int,SU3>,Eigen::MatrixXd>& Kinv, int AA, int N1v){
  bool write;
  if(bra.gamma()==1
  && bra.Nex_sigma()==0
  && bra.x_sigma().lm()==2
  && bra.x_sigma().mu()==0
  && bra.SpSnS().SSp()==1
  && bra.SpSnS().SSn()==1
  && bra.SpSnS().SS()==0
  && bra.upsilon()==1
  && bra.Nex_omega()==0
  && bra.x_omega().lm()==2
  && bra.x_omega().mu()==0
  && ket.gamma()==1
  && ket.Nex_sigma()==0
  && ket.x_sigma().lm()==2
  && ket.x_sigma().mu()==0
  && ket.SpSnS().SSp()==1
  && ket.SpSnS().SSn()==1
  && ket.SpSnS().SS()==0
  && ket.upsilon()==1
  && ket.Nex_omega()==6
  && ket.x_omega().lm()==4
  && ket.x_omega().mu()==2
  && tensor[0]==0
  && tensor[1]==6
  && tensor[2]==0
  && tensor[3]==6
  && tensor[4]==0
  ){
    write=true;
  }else{
    write=false;
  }
  write=false;
  if(write)std::cout<<"************************ Recurrence begins ***********************"<<std::endl;

  int rhomax=su3::mult(ket.x_omega().lm(),ket.x_omega().mu(),tensor[2],tensor[3],bra.x_omega().lm(),bra.x_omega().mu());
  std::vector<double> res(rhomax);
  int Nn1=ket.Nex_omega()-ket.Nex_sigma();
  int Nn1p=Nn1-2;
  int Nex_omega1p=ket.Nex_omega()-2;
  int Nex_omegapp=bra.Nex_omega()-2;
  for(int rho=1; rho<=rhomax; rho++){
    res[rho-1]=0.0;
    int ind1=-1;
    // sum over n1
    for(int n11=0; n11<=Nn1; n11=n11+2){
      for(int n12=0; n12<=std::min(n11,Nn1-n11); n12=n12+2){
        int n13=Nn1-n11-n12;
        if(n13>n12)continue;
        int lambda_n1=n11-n12;
        int mu_n1=n12-n13;
        // sum over rho1
        for(int rho1=1; rho1<=su3::mult(ket.x_sigma().lm(),ket.x_sigma().mu(),lambda_n1,mu_n1,ket.x_omega().lm(),ket.x_omega().mu()); rho1++){
          ind1++;
          if(write)std::cout<<"n1, rho1: ["<<n11<<","<<n12<<","<<n13<<"] or "<<Nn1<<"("<<lambda_n1<<","<<mu_n1<<"), "<<rho1<<std::endl;
          // construct set of omega1' labels
          std::vector<SU3> x1p_set;
          for(int n1p1=0; n1p1<=Nn1p; n1p1=n1p1+2){
            for(int n1p2=0; n1p2<=std::min(n1p1,Nn1p-n1p1); n1p2=n1p2+2){
              int n1p3=Nn1p-n1p1-n1p2;
              if(n1p3>n1p2)continue;
              int lambda_n1p=n1p1-n1p2;
              int mu_n1p=n1p2-n1p3;
              for(std::tuple<SU3,int> x1p : KroneckerProduct(ket.x_sigma().lm(),ket.x_sigma().mu(),lambda_n1p,mu_n1p)){
		auto iter=std::find(x1p_set.begin(),x1p_set.end(),std::get<0>(x1p));
                if(iter==x1p_set.end()){
                  x1p_set.push_back(std::get<0>(x1p));
		}
              }
            }
          }
          double sum_omega1p_n1p_rho1p=0.0;
          //sum over omega1'
          for(SU3 x1p : x1p_set){
            double dimx1p=double(su3dim(x1p));
            int ind1p=-1;
            // sum over n1'
            for(int n1p1=0; n1p1<=Nn1p; n1p1=n1p1+2){
              for(int n1p2=0; n1p2<=std::min(n1p1,Nn1p-n1p1); n1p2=n1p2+2){
                int n1p3=Nn1p-n1p1-n1p2;
                if(n1p3>n1p2)continue;
                int lambda_n1p=n1p1-n1p2;
                int mu_n1p=n1p2-n1p3;
                // sum over rho1'
                for(int rho1p=1; rho1p<=su3::mult(ket.x_sigma().lm(),ket.x_sigma().mu(),lambda_n1p,mu_n1p,x1p.lm(),x1p.mu()); rho1p++){
                  ind1p++;
                  if(write)std::cout<<"omega1', n1', rho1': "<<Nex_omega1p<<"("<<x1p.lm()<<","<<x1p.mu()<<"), ["
                          <<n1p1<<","<<n1p2<<","<<n1p3<<"] or "<<Nn1p<<"("<<lambda_n1p<<","<<mu_n1p<<"), "<<rho1p<<std::endl;
                  double sum_upsilon1p=0.0;
                  // sum over upsilon1'
                  for(int upsilon1p=1; upsilon1p<=K[{ket.gamma(),ket.Nex_sigma(),ket.x_sigma(),
                                  ket.SpSnS(),Nex_omega1p,x1p}].cols(); upsilon1p++){
                    if(write)std::cout<<"upsilon1': "<<upsilon1p<<std::endl;
                    ONBasisState ketp(ket.gamma(),ket.Nex_sigma(),ket.x_sigma(),
                                    ket.SpSnS(),upsilon1p,Nex_omega1p,x1p);

                    // first term
                    double term1=0.0;
                    if(bra.Nex_omega()!=bra.Nex_sigma()){
                      // sum over omega''
                      for(std::tuple<SU3,int> xpp : KroneckerProduct(x1p.lm(),x1p.mu(),tensor[2],tensor[3])){
                        if(write)std::cout<<"omega'': "<<Nex_omegapp<<"("<<std::get<0>(xpp).lm()<<","<<std::get<0>(xpp).mu()<<")"<<std::endl;
                        double sum_upsilonpp=0.0;
                        // sum over upsilonpp
                        for(int upsilonpp=1; upsilonpp<=K[{bra.gamma(),bra.Nex_sigma(),bra.x_sigma(),
                             bra.SpSnS(),Nex_omegapp,std::get<0>(xpp)}].cols(); upsilonpp++){
                          if(write)std::cout<<"upsilon'': "<<upsilonpp<<std::endl;
                          ONBasisState brapp(bra.gamma(),bra.Nex_sigma(),bra.x_sigma(),
                                          bra.SpSnS(),upsilonpp,Nex_omegapp,std::get<0>(xpp));
                          double sum_rho3=0.0;
                          std::tuple<ONBasisState,ONBasisState,std::array<int,5>> key={brapp,ketp,tensor};
                          std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,
                                  std::vector<double>>::iterator it=obdrmes.find(key);
                          if(it!=obdrmes.end()){
                            // sum over rho3
                            for(int rho3=1; rho3<=std::get<1>(xpp); rho3++){
                              if(write)std::cout<<"rho3: "<<rho3<<std::endl;
                              double sum_rho4=0.0;
                              // sum over rho4
                              for(int rho4=1; rho4<=std::get<1>(xpp); rho4++){
                                if(write)std::cout<<"rho4: "<<rho4<<std::endl;
                                double sum_rho5=0.0;
                                // sum over rho5
                                for(int rho5=1; rho5<=std::get<1>(xpp); rho5++){
                                  if(write){
                                    std::cout<<"rho5: "<<rho5<<std::endl;
                                    std::cout<<"sum_rho5=sum_rho5+Phi*U="<<sum_rho5<<"+"
                                      <<Phi(rho4,rho5,x1p.lm(),x1p.mu(),std::get<0>(xpp).mu(),std::get<0>(xpp).lm(),tensor[3],tensor[2])
                                      <<"*"<<U(bra.x_omega().lm(),bra.x_omega().mu(),std::get<0>(xpp).mu(),std::get<0>(xpp).lm(),
                                             ket.x_omega().lm(),ket.x_omega().mu(),x1p.lm(),x1p.mu(),2,0,1,1,tensor[3],tensor[2],rho5,rho);
                                  }
                                  sum_rho5+=Phi(rho4,rho5,x1p.lm(),x1p.mu(),std::get<0>(xpp).mu(),std::get<0>(xpp).lm(),tensor[3],tensor[2])
                                          *U(bra.x_omega().lm(),bra.x_omega().mu(),std::get<0>(xpp).mu(),std::get<0>(xpp).lm(),
                                             ket.x_omega().lm(),ket.x_omega().mu(),x1p.lm(),x1p.mu(),2,0,1,1,tensor[3],tensor[2],rho5,rho);
                                  if(write)std::cout<<"="<<sum_rho5<<std::endl;
                                }
                                if(write)std::cout<<"sum_rho4=sum_rho4+Phi*sum_rho5="<<sum_rho4<<"+"
                                  <<Phi(rho3,rho4,std::get<0>(xpp).lm(),std::get<0>(xpp).mu(),tensor[3],tensor[2],x1p.lm(),x1p.mu())
                                  <<"*"<<sum_rho5;
                                sum_rho4+=Phi(rho3,rho4,std::get<0>(xpp).lm(),std::get<0>(xpp).mu(),tensor[3],tensor[2],
                                                x1p.lm(),x1p.mu())*sum_rho5;
                                if(write)std::cout<<"="<<sum_rho4<<std::endl;
                              }
                              if(write)std::cout<<"sum_rho3=sum_rho3+RME*sum_rho4="<<sum_rho3<<"+"<<obdrmes[key][rho3-1]<<"*"<<sum_rho4;
                              sum_rho3+=obdrmes[key][rho3-1]*sum_rho4;
                              if(write)std::cout<<"="<<sum_rho3<<std::endl;
                            }
                          }
                          if(write)std::cout<<"sum_upsilonpp=sum_upsilonpp+ARME*sum_rho3="<<sum_upsilonpp<<"+"
                            <<ARME(bra.gamma(),bra.Nex_sigma(),bra.x_sigma(),bra.SpSnS(),
                              bra.upsilon(),bra.x_omega(),upsilonpp,Nex_omegapp,
                              std::get<0>(xpp),orthonormal_basis,A)<<"*"<<sum_rho3;
                          sum_upsilonpp+=ARME(bra.gamma(),bra.Nex_sigma(),bra.x_sigma(),bra.SpSnS(),
                            bra.upsilon(),bra.x_omega(),upsilonpp,Nex_omegapp,
                            std::get<0>(xpp),orthonormal_basis,A)*sum_rho3;
                          if(write)std::cout<<"="<<sum_upsilonpp<<std::endl;
                        }
                        if(write)std::cout<<"term1=term1+sqrt*sum_upsilonpp="<<term1<<"+"
                          <<sqrt(double(su3dim(std::get<0>(xpp))))<<"*"<<sum_upsilonpp;
                        term1+=sqrt(double(su3dim(std::get<0>(xpp))))*sum_upsilonpp;
                        if(write)std::cout<<"="<<term1<<std::endl;
                      }
                      term1/=sqrt(6.0);
                      int phase=tensor[2]+tensor[3]+ket.x_omega().lm()+ket.x_omega().mu()+bra.x_omega().lm()+bra.x_omega().mu();
                      if((phase/2)*2!=phase)term1=-term1;
                      if(write)std::cout<<"term1=(-1)^phase * term1 / sqrt(6)="<<term1<<std::endl;
                    }
                    if(write)std::cout<<"1st term: "<<term1<<std::endl;

                    // second term
                    double term2=0.0;
                    if(tensor[0]+2<=bra.Nex_omega()+N1v){
                      // sum over Gamma'
                      for(std::tuple<SU3,int> Gammap : KroneckerProduct(tensor[0]+2,0,0,tensor[1])){
                        if(write)std::cout<<"Gamma': ("<<std::get<0>(Gammap).lm()<<","<<std::get<0>(Gammap).mu()<<")"<<std::endl;
                        double factor=sqrt(double(su3dim(std::get<0>(Gammap))))*U(2,0,tensor[0],0,std::get<0>(Gammap).lm(),
                                        std::get<0>(Gammap).mu(),0,tensor[1],tensor[0]+2,0,1,1,tensor[2],tensor[3],1,1);
                        if(((std::get<0>(Gammap).lm()+std::get<0>(Gammap).mu())/2)*2!=std::get<0>(Gammap).lm()+std::get<0>(Gammap).mu())
                                factor=-factor;
                        double sum_rho3=0.0;
                        std::tuple<ONBasisState,ONBasisState,std::array<int,5>>
                                  key={bra,ketp,{tensor[0]+2,tensor[1],std::get<0>(Gammap).lm(),std::get<0>(Gammap).mu(),tensor[4]}};
                        std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,
                                std::vector<double>>::iterator it=obdrmes.find(key);
                        if(it!=obdrmes.end()){
                          // sum over rho3
                          for(int rho3=1; rho3<=su3::mult(x1p.lm(),x1p.mu(),std::get<0>(Gammap).lm(),std::get<0>(Gammap).mu(),
                                                bra.x_omega().lm(),bra.x_omega().mu()); rho3++){
                            if(write){
                              std::cout<<"rho3: "<<rho3<<std::endl;
                              std::cout<<"sum_rho3=sum_rho3+RME*U="<<sum_rho3<<"+"<<obdrmes[key][rho3-1]<<"*"
                                <<U(bra.x_omega().lm(),bra.x_omega().mu(),std::get<0>(Gammap).mu(),std::get<0>(Gammap).lm(),
                                    ket.x_omega().lm(),ket.x_omega().mu(),2,0,x1p.lm(),x1p.mu(),rho3,1,tensor[3],tensor[2],1,rho);
                            }
                            sum_rho3+=obdrmes[key][rho3-1]*U(bra.x_omega().lm(),bra.x_omega().mu(),std::get<0>(Gammap).mu(),
                              std::get<0>(Gammap).lm(),ket.x_omega().lm(),ket.x_omega().mu(),2,0,x1p.lm(),x1p.mu(),rho3,1,
                              tensor[3],tensor[2],1,rho);
                            if(write)std::cout<<"="<<sum_rho3<<std::endl;
                          }
                        }
                        if(write)std::cout<<"term2=term2+factor*sum_rho3="<<term2<<"+"<<factor<<"*"<<sum_rho3;
                        term2+=factor*sum_rho3;
                        if(write)std::cout<<"="<<term2<<std::endl;
                      }
                      term2*=(1.0+1.0/double(AA))*sqrt(double(su3dim(tensor[0],0))/dimx1p);
                      if(((tensor[2]+tensor[3])/2)*2==tensor[2]+tensor[3])term2=-term2;
                      if(write)std::cout<<"term2=factor*term2="<<term2<<std::endl;
                    }
                    if(write)std::cout<<"2nd term: "<<term2<<std::endl;

                    // third term
                    double term3=0.0;
                    if(tensor[1]-2>=0){
                      // sum over Gamma'
                      for(std::tuple<SU3,int> Gammap : KroneckerProduct(tensor[0],0,0,tensor[1]-2)){
                        if(write)std::cout<<"Gamma': ("<<std::get<0>(Gammap).lm()<<","<<std::get<0>(Gammap).mu()<<")"<<std::endl;
                        double factor=sqrt(double(su3dim(std::get<0>(Gammap))))*U(tensor[0],0,0,tensor[1],std::get<0>(Gammap).lm(),
                                        std::get<0>(Gammap).mu(),2,0,tensor[2],tensor[3],1,1,0,tensor[1]-2,1,1);
                        double sum_rho3=0.0;
                        std::tuple<ONBasisState,ONBasisState,std::array<int,5>>
                                  key={bra,ketp,{tensor[0],tensor[1]-2,std::get<0>(Gammap).lm(),std::get<0>(Gammap).mu(),tensor[4]}};
                        std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,
                                std::vector<double>>::iterator it=obdrmes.find(key);
                        if(it!=obdrmes.end()){
                          // sum over rho3
                          for(int rho3=1; rho3<=su3::mult(x1p.lm(),x1p.mu(),std::get<0>(Gammap).lm(),std::get<0>(Gammap).mu(),
                                                bra.x_omega().lm(),bra.x_omega().mu()); rho3++){
                            if(write){
                              std::cout<<"rho3: "<<rho3<<std::endl;
                              std::cout<<"sum_rho3=sum_rho3+RME*U="<<sum_rho3<<"+"<<obdrmes[key][rho3-1]<<"*"
                                <<U(bra.x_omega().lm(),bra.x_omega().mu(),std::get<0>(Gammap).mu(),std::get<0>(Gammap).lm(),
                                    ket.x_omega().lm(),ket.x_omega().mu(),2,0,x1p.lm(),x1p.mu(),rho3,1,tensor[3],tensor[2],1,rho);
                            }
                            sum_rho3+=obdrmes[key][rho3-1]*U(bra.x_omega().lm(),bra.x_omega().mu(),std::get<0>(Gammap).mu(),
                              std::get<0>(Gammap).lm(),ket.x_omega().lm(),ket.x_omega().mu(),2,0,x1p.lm(),x1p.mu(),rho3,1,
                              tensor[3],tensor[2],1,rho);
                            if(write)std::cout<<"="<<sum_rho3<<std::endl;
                          }
                        }
                        if(write)std::cout<<"term3=term3+factor*sum_rho3="<<term3<<"+"<<factor<<"*"<<sum_rho3;
                        term3+=factor*sum_rho3;
                        if(write)std::cout<<"="<<term3<<std::endl;
                      }
                      term3*=(1.0-1.0/double(AA))*sqrt(double(su3dim(tensor[1],0))/dimx1p);
                      if(write)std::cout<<"term3=factor*term3="<<term3<<std::endl;
                    }
                    if(write)std::cout<<"3rd term: "<<term3<<std::endl;

                    // fourth term
                    double term4=0.0;
                    if(tensor[0]+1<=bra.Nex_omega()+N1v && tensor[1]-1>=0){
                      // sum over Gamma'
                      for(std::tuple<SU3,int> Gammap : KroneckerProduct(tensor[2],tensor[3],1,0)){
                        if(su3::mult(tensor[0],0,0,tensor[1]-1,std::get<0>(Gammap).lm(),std::get<0>(Gammap).mu())==0)continue;
                        if(write)std::cout<<"Gamma': ("<<std::get<0>(Gammap).lm()<<","<<std::get<0>(Gammap).mu()<<")"<<std::endl;
                        double sum_Gammapp=0.0;
                        // sum over Gamma''
                        for(std::tuple<SU3,int> Gammapp : KroneckerProduct(tensor[0]+1,0,0,tensor[1]-1)){
                          if(write)std::cout<<"Gamma'': ("<<std::get<0>(Gammapp).lm()<<","<<std::get<0>(Gammapp).mu()<<")"<<std::endl;
                          double factor=sqrt(double(su3dim(std::get<0>(Gammapp))))*U(1,0,tensor[0],0,std::get<0>(Gammapp).lm(),
                            std::get<0>(Gammapp).mu(),0,tensor[1]-1,tensor[0]+1,0,1,1,std::get<0>(Gammap).lm(),std::get<0>(Gammap).mu(),1,1)
                            *U(tensor[2],tensor[3],1,0,std::get<0>(Gammapp).lm(),std::get<0>(Gammapp).mu(),1,0,
                               std::get<0>(Gammap).lm(),std::get<0>(Gammap).mu(),1,1,2,0,1,1);
                          int phase=std::get<0>(Gammap).lm()+std::get<0>(Gammap).mu()+std::get<0>(Gammapp).lm()+std::get<0>(Gammapp).mu();
                          if((phase/2)*2!=phase)factor=-factor;
                          double sum_rho3=0.0;
                          std::tuple<ONBasisState,ONBasisState,std::array<int,5>>
                                  key={bra,ketp,{tensor[0]+1,tensor[1]-1,std::get<0>(Gammapp).lm(),std::get<0>(Gammapp).mu(),tensor[4]}};
                          std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,
                                std::vector<double>>::iterator it=obdrmes.find(key);
                          if(it!=obdrmes.end()){
                            // sum over rho3
                            for(int rho3=1; rho3<=su3::mult(x1p.lm(),x1p.mu(),std::get<0>(Gammapp).lm(),std::get<0>(Gammapp).mu(),
                                                bra.x_omega().lm(),bra.x_omega().mu()); rho3++){
                              if(write){
                                std::cout<<"rho3: "<<rho3<<std::endl;
                                std::cout<<"sum_rho3=sum_rho3+RME*U="<<sum_rho3<<"+"<<obdrmes[key][rho3-1]<<"*"
                                  <<U(bra.x_omega().lm(),bra.x_omega().mu(),std::get<0>(Gammapp).mu(),std::get<0>(Gammapp).lm(),
                                      ket.x_omega().lm(),ket.x_omega().mu(),2,0,x1p.lm(),x1p.mu(),rho3,1,tensor[3],tensor[2],1,rho);
                              }
                              sum_rho3+=obdrmes[key][rho3-1]*U(bra.x_omega().lm(),bra.x_omega().mu(),std::get<0>(Gammapp).mu(),
                                std::get<0>(Gammapp).lm(),ket.x_omega().lm(),ket.x_omega().mu(),2,0,x1p.lm(),x1p.mu(),rho3,1,
                                tensor[3],tensor[2],1,rho);
                              if(write)std::cout<<"="<<sum_rho3<<std::endl;
                            }
                          }
                          if(write)std::cout<<"sum_Gammapp=sum_Gammapp+factor*sum_rho3="<<sum_Gammapp<<"+"<<factor<<"*"<<sum_rho3;
                          sum_Gammapp+=factor*sum_rho3;
                          if(write)std::cout<<"="<<sum_Gammapp<<std::endl;
                        }
                        if(write)std::cout<<"term4=term4+U*sum_Gamapp="<<term4<<"+"
                          <<U(tensor[0],0,0,tensor[1],std::get<0>(Gammap).lm(),std::get<0>(Gammap).mu(),1,0,
                                        tensor[2],tensor[3],1,1,0,tensor[1]-1,1,1)<<"*"<<sum_Gammapp;
                        term4+=U(tensor[0],0,0,tensor[1],std::get<0>(Gammap).lm(),std::get<0>(Gammap).mu(),1,0,
                                        tensor[2],tensor[3],1,1,0,tensor[1]-1,1,1)*sum_Gammapp;
                        if(write)std::cout<<"="<<term4<<std::endl;
                      }
                      term4*=sqrt(double(2*(tensor[1]+2)*(tensor[0]+1))/dimx1p)/double(AA);
                      if(write)std::cout<<"term4=factor*term4="<<term4<<std::endl;
                    }
                    if(write){
                      std::cout<<"4th term: -"<<term4<<std::endl;

                      std::cout<<"sum_upsilon1p=sum_upsilon1p+K*terms="<<sum_upsilon1p<<"+"<<K[{ket.gamma(),ket.Nex_sigma(),
                              ket.x_sigma(),ket.SpSnS(),Nex_omega1p,x1p}]
                                      (ind1p,upsilon1p-1)<<"*"<<(term1+term2+term3-term4);
                    }
                    sum_upsilon1p+=K[{ket.gamma(),ket.Nex_sigma(),ket.x_sigma(),
                                  ket.SpSnS(),Nex_omega1p,x1p}](ind1p,upsilon1p-1)*(term1+term2+term3-term4);
                    if(write)std::cout<<"="<<sum_upsilon1p<<std::endl;
                  }
                  if(write)std::cout<<"sum_omega1p_n1p_rho1p=sum_omega1p_n1p_rho1p+U*X*sum_upsilon1p="<<sum_omega1p_n1p_rho1p<<"+"
                    <<U(2,0,lambda_n1p,mu_n1p,ket.x_omega().lm(),ket.x_omega().mu(),ket.x_sigma().lm(),ket.x_sigma().mu(),
                                  lambda_n1,mu_n1,1,rho1,x1p.lm(),x1p.mu(),rho1p,1)<<"*"
                    <<2.0*BosonCreationOperatorME(n11,n12,n13,n1p1,n1p2,n1p3)/double(Nn1)<<"*"<<sum_upsilon1p;
                  sum_omega1p_n1p_rho1p+=U(2,0,lambda_n1p,mu_n1p,ket.x_omega().lm(),ket.x_omega().mu(),ket.x_sigma().lm(),ket.x_sigma().mu(),
                                  lambda_n1,mu_n1,1,rho1,x1p.lm(),x1p.mu(),rho1p,1)
                          *2.0*BosonCreationOperatorME(n11,n12,n13,n1p1,n1p2,n1p3)*sum_upsilon1p/double(Nn1);
                  if(write)std::cout<<"="<<sum_omega1p_n1p_rho1p<<std::endl;
                }
              }
            }
          }
          if(write)std::cout<<"res=res+Kinv*sum_omega1p_n1p_rho1p="<<res[rho-1]<<"+"<<Kinv[{ket.gamma(),
                  ket.Nex_sigma(),ket.x_sigma(),ket.SpSnS(),
                  ket.Nex_omega(),ket.x_omega()}](ket.upsilon()-1,ind1)<<"*"<<sum_omega1p_n1p_rho1p;
          res[rho-1]+=Kinv[{ket.gamma(),ket.Nex_sigma(),ket.x_sigma(),ket.SpSnS(),
                  ket.Nex_omega(),ket.x_omega()}](ket.upsilon()-1,ind1)*sum_omega1p_n1p_rho1p;
          if(write)std::cout<<"="<<res[rho-1]<<std::endl;
        }
      }
    }
    res[rho-1]*=sqrt(double(su3dim(ket.x_omega()))/double(su3dim(tensor[2],tensor[3])));
    if(write)std::cout<<"res=factor*res="<<res[rho-1]<<std::endl;
  }
  if(write)std::cout<<"************************ Recurrence ends *************************"<<std::endl;
  return res;
}

std::vector<std::array<int,13>> ReadListOfTBTensors(const std::string& filename){
  std::vector<std::array<int,13>> tensors;
  std::ifstream file(filename);
  if(!file){
    std::cerr << "ERROR: Could not open '" << filename << "' input file!" << std::endl;
    su3::finalize();
    exit(EXIT_FAILURE);
  }
  while(true){
    int n1, n2, n3, n4, lmf, muf, ssf, lmi, mui, ssi, lm0, mu0, ss0;
    file >> n1 >> n2 >> n3 >> n4 >> lmf >> muf >> ssf >> lmi >> mui >> ssi >> lm0 >> mu0 >> ss0;
    if(file){
      tensors.push_back({n1, n2, n3, n4, lmf, muf, ssf, lmi, mui, ssi, lm0, mu0, ss0});
    }else{
      break;
    }
  }
  return tensors;
}

void ParseTBInputFile(const std::stringstream& filename, const std::array<int,13>& tensor,
	   std::map<Spins,std::map<int,std::vector<SU3>>>& SpSnSNexGamma_set,
	   std::map<std::tuple<Spins,int,SU3>,std::vector<Ipini>>& lsu3shell_basis, int irho0,
	   std::map<std::tuple<LSU3shellBasisState,LSU3shellBasisState,std::array<int,14>>,std::vector<double>>& tbdrmes_lsu3shell){
  std::ifstream input_file(filename.str());
  if(!input_file){
    std::cerr << "ERROR: Could not open file '" << filename.str() << std::endl;
    su3::finalize();
    exit(EXIT_FAILURE);
  }

  std::array<int,14> tensorirho0={tensor[0],tensor[1],tensor[2],tensor[3],tensor[4],tensor[5],tensor[6],tensor[7],tensor[8],tensor[9],tensor[10],tensor[11],tensor[12],irho0};

  while (true) {
    int ip,in,bra_max,SSp_bra,SSn_bra,SS_bra,N_bra,lm_bra,mu_bra,jp,jn,ket_max,SSp_ket,SSn_ket,SS_ket,N_ket,lm_ket,mu_ket;
    input_file>>ip>>in>>bra_max>>SSp_bra>>SSn_bra>>SS_bra>>N_bra>>lm_bra>>mu_bra
              >>jp>>jn>>ket_max>>SSp_ket>>SSn_ket>>SS_ket>>N_ket>>lm_ket>>mu_ket;
    if(input_file){
      Spins SpSnS_bra(SSp_bra,SSn_bra,SS_bra);
      SU3 x_bra(lm_bra,mu_bra);
      std::map<Spins,std::map<int,std::vector<SU3>>>::iterator it=SpSnSNexGamma_set.find(SpSnS_bra);
      if(it==SpSnSNexGamma_set.end()){
        std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
        su3::finalize();
        exit(EXIT_FAILURE);
      }else{
        std::map<int,std::vector<SU3>>::iterator it2=SpSnSNexGamma_set[SpSnS_bra].find(N_bra);
        if(it2==SpSnSNexGamma_set[SpSnS_bra].end()){
          std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
          su3::finalize();
          exit(EXIT_FAILURE);
        }else{
          auto it3=std::find(SpSnSNexGamma_set[SpSnS_bra][N_bra].begin(), SpSnSNexGamma_set[SpSnS_bra][N_bra].end(), x_bra);
          if(it3==SpSnSNexGamma_set[SpSnS_bra][N_bra].end()){
            std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
            su3::finalize();
            exit(EXIT_FAILURE);
          }
        }
      }
      Spins SpSnS_ket(SSp_ket,SSn_ket,SS_ket);
      SU3 x_ket(lm_ket,mu_ket);
      it=SpSnSNexGamma_set.find(SpSnS_ket);
      if(it==SpSnSNexGamma_set.end()){
        std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
        su3::finalize();
        exit(EXIT_FAILURE);
      }else{
        std::map<int,std::vector<SU3>>::iterator it2=SpSnSNexGamma_set[SpSnS_ket].find(N_ket);
        if(it2==SpSnSNexGamma_set[SpSnS_ket].end()){
          std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
          su3::finalize();
          exit(EXIT_FAILURE);
	}else{
          auto it3=std::find(SpSnSNexGamma_set[SpSnS_ket][N_ket].begin(), SpSnSNexGamma_set[SpSnS_ket][N_ket].end(), x_ket);
          if(it3==SpSnSNexGamma_set[SpSnS_ket][N_ket].end()){
            std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
            su3::finalize();
            exit(EXIT_FAILURE);
          }
        }
      }
      int rhomax=su3::mult(lm_ket,mu_ket,tensor[10],tensor[11],lm_bra,mu_bra);
      for(int i=0; i<bra_max; i++){
	Ipini ipini_bra(ip,in,i);
	LSU3shellBasisState bra(ipini_bra,SpSnS_bra,N_bra,x_bra);
	std::map<std::tuple<Spins,int,SU3>,std::vector<Ipini>>::iterator iter=lsu3shell_basis.find({SpSnS_bra,N_bra,x_bra});
        if(iter==lsu3shell_basis.end()){
          std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
          su3::finalize();
          exit(EXIT_FAILURE);
        }else{
          auto iter2=std::find(lsu3shell_basis[{SpSnS_bra,N_bra,x_bra}].begin(), lsu3shell_basis[{SpSnS_bra,N_bra,x_bra}].end(), ipini_bra);
          if(iter2==lsu3shell_basis[{SpSnS_bra,N_bra,x_bra}].end()){
            std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
            su3::finalize();
            exit(EXIT_FAILURE);
          }
        }
	for(int j=0; j<ket_max; j++){
	  Ipini ipini_ket(jp,jn,j);
	  LSU3shellBasisState ket(ipini_ket,SpSnS_ket,N_ket,x_ket);
	  iter=lsu3shell_basis.find({SpSnS_ket,N_ket,x_ket});
          if(iter==lsu3shell_basis.end()){
            std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
            su3::finalize();
            exit(EXIT_FAILURE);
          }else{
            auto iter2=std::find(lsu3shell_basis[{SpSnS_ket,N_ket,x_ket}].begin(), lsu3shell_basis[{SpSnS_ket,N_ket,x_ket}].end(), ipini_ket);
            if(iter2==lsu3shell_basis[{SpSnS_ket,N_ket,x_ket}].end()){
              std::cout<<"ERROR: New lsu3shell basis state found in TBD RME files"<<std::endl;
              su3::finalize();
              exit(EXIT_FAILURE);
            }
          }
	  std::vector<double> rmes(rhomax);
          for(int irho=0; irho<rhomax; irho++){
            input_file>>rmes[irho];
	  }
	  tbdrmes_lsu3shell[{bra,ket,tensorirho0}]=rmes;
        }
      }

    }else{
      break;
    }

  }
}

int main(int argc, char** argv) {
  if (argc != 7) {
    std::cout << "Usage: " << argv[0] << " <OBD tensor list file> <TBD tensor list file> A Nmax N1vp N1vn" << std::endl;
    return EXIT_FAILURE;
  }
  su3::init(250);
/*
//---------------------------------------------------
int hore=5;
for(int lm1=0; lm1<=hore; lm1++){
  std::cout<<"lm1="<<lm1<<std::endl;
  for(int mu1=0; mu1<=hore; mu1++){
    for(int lm2=0; lm2<=hore; lm2++){
      for(int mu2=0; mu2<=hore; mu2++){
	for(int lm3=0; lm3<=hore; lm3++){
          for(int mu3=0; mu3<=hore; mu3++){
            int rhomax=su3::mult(lm1,mu1,lm2,mu2,lm3,mu3);
	    if(rhomax==0)continue;
	    std::vector<double> su3cgsl,su3cgsr;
	    int reserv=std::max(lm1,mu1)*std::max(lm2,mu2)*std::max(lm3,mu3);
            su3cgsl.reserve(reserv);
	    su3cgsr.reserve(reserv);
	    for(int L1=0; L1<=lm1+mu1; L1++){
              int k1max=su3::kmax(lm1,mu1,L1);
	      if(k1max==0)continue;
	      for(int L2=0; L2<=lm2+mu2; L2++){
                int k2max=su3::kmax(lm2,mu2,L2);
                if(k2max==0)continue;
	        for(int L3=std::abs(L1-L2); L3<=std::min(lm3+mu3,L1+L2); L3++){
                  int k3max=su3::kmax(lm3,mu3,L3);
                  if(k3max==0)continue;
		  double factor=sqrt(double(su3dim(lm3,mu3)*(2*L1+1))/double(su3dim(lm1,mu1)*(2*L3+1)));
		  int phase=lm1+mu1-lm3-mu3+L1+L2-L3;
		  if((phase/2)*2!=phase)factor=-factor;
                  su3::wu3r3w(lm1,mu1,lm2,mu2,lm3,mu3,L1,L2,L3,k3max,k2max,k1max,rhomax,su3cgsl);
		  su3::wu3r3w(lm3,mu3,mu2,lm2,lm1,mu1,L3,L2,L1,k1max,k2max,k3max,rhomax,su3cgsr);
		  int indl=0;
                  for(int ik3=0; ik3<k3max; ik3++){
		    for(int ik2=0; ik2<k2max; ik2++){
		      for(int ik1=0; ik1<k1max; ik1++){
			for(int irho=0; irho<rhomax; irho++){
	                  if(indl!=irho+rhomax*ik1+rhomax*k1max*ik2+rhomax*k1max*k2max*ik3){
                            std::cout<<"ERROR: indexy nesedia"<<std::endl;
			    su3::finalize();
                            return EXIT_SUCCESS;
			  }
			  int indr=irho+rhomax*ik3+rhomax*k3max*ik2+rhomax*k3max*k2max*ik1;
			  double rhs=factor*su3cgsr[indr];
			  if(std::abs(su3cgsl[indl]-rhs)>1.0e-8){
                            std::cout<<su3cgsl[indl]<<"!="<<factor<<"*"<<su3cgsr[indr]<<"="<<rhs<<std::endl;
			    std::cout<<"lm1 mu1 lm2 mu2 lm3 mu3 l1 l2 l3 k1 k2 k3 rho: "<<lm1<<" "<<mu1<<" "<<lm2<<" "<<mu2<<" "
				     <<lm3<<" "<<mu3<<" "<<L1<<" "<<L2<<" "<<L3<<" "<<ik1+1<<" "<<ik2+1<<" "<<ik3+1<<" "<<irho+1<<std::endl;
			    su3::finalize();
                            return EXIT_SUCCESS;
			  }
			  indl++;
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
su3::finalize();
return EXIT_SUCCESS;
//---------------------------------------------------
*/
  const int AA = std::stoi(argv[3]);
  const int Nmax = std::stoi(argv[4]);
  const int N1vp = std::stoi(argv[5]);
  const int N1vn = std::stoi(argv[6]);
  std::cout << "A = " << AA << std::endl;
  std::cout << "Nmax = " << Nmax << std::endl;
  std::cout << "N1vp = " << N1vp << std::endl;
  std::cout << "N1vn = " << N1vn << std::endl;

  std::map<Spins,std::map<int,std::vector<SU3>>> SpSnSNexGamma_set;
  // For each Sp,Sn,S there is map between Nex and vector of SU3 irreps.

  std::map<std::tuple<Spins,int,SU3>,std::vector<Ipini>> lsu3shell_basis;
  // For each Sp,Sn,S,Nex,Gamma there is vector of ip,in,i values.
 
  std::map<std::tuple<LSU3shellBasisState,LSU3shellBasisState,std::array<int,5>>,std::vector<double>>
	  proton_obdrmes_lsu3shell,neutron_obdrmes_lsu3shell;
  std::map<std::tuple<LSU3shellBasisState,LSU3shellBasisState,std::array<int,14>>,std::vector<double>>
          proton_tbdrmes_lsu3shell,neutron_tbdrmes_lsu3shell,pn_tbdrmes_lsu3shell;
  // Third member of the tuple contains N1,N2,N3,N4,lmf,muf,SSf,lmi,mui,SSi,lm0,mu0,SS0,irho0.

  std::array<std::vector<std::array<int,5>>,2> tensors=ReadListOfTensors(argv[1]);
  std::vector<std::array<int,13>> TBtensors=ReadListOfTBTensors(argv[2]);

  // Input and LSU3shell basis construction
  for(int g=0; g<=1; g++){
    for(const std::array<int,5>& tensor : tensors[g]){
      std::stringstream proton_input_filename;
      proton_input_filename << "a+" << tensor[0] << "ta" << tensor[1] << "_" << tensor[2] << "_"
                            << tensor[3] << "_" << tensor[4] << "_proton";
      std::stringstream neutron_input_filename;
      neutron_input_filename << "a+" << tensor[0] << "ta" << tensor[1] << "_" << tensor[2] << "_"
                             << tensor[3] << "_" << tensor[4] << "_neutron";
      ParseInputFile(proton_input_filename,tensor,SpSnSNexGamma_set,lsu3shell_basis,proton_obdrmes_lsu3shell);
      ParseInputFile(neutron_input_filename,tensor,SpSnSNexGamma_set,lsu3shell_basis,neutron_obdrmes_lsu3shell);
    }
  }
  for(const std::array<int,13>& tensor : TBtensors){
    for(int irho0=0; irho0<su3::mult(tensor[4], tensor[5], tensor[7], tensor[8], tensor[10], tensor[11]); irho0++){
      std::stringstream proton_input_filename;
      proton_input_filename << tensor[0] << "_" << tensor[1] << "_" << tensor[2] << "_" << tensor[3] << "--"
  	                    << tensor[4] << "_" << tensor[5] << "_" << tensor[6] << "--"
	  	  	    << tensor[7] << "_" << tensor[8] << "_" << tensor[9] << "--"
			    << irho0 << "_" << tensor[10] << "_" << tensor[11] << "_" << tensor[12] << ".pp";
      std::stringstream neutron_input_filename;
      neutron_input_filename << tensor[0] << "_" << tensor[1] << "_" << tensor[2] << "_" << tensor[3] << "--"
                            << tensor[4] << "_" << tensor[5] << "_" << tensor[6] << "--"
                            << tensor[7] << "_" << tensor[8] << "_" << tensor[9] << "--"
                            << irho0 << "_" << tensor[10] << "_" << tensor[11] << "_" << tensor[12] << ".nn";
      std::stringstream pn_input_filename;
      pn_input_filename << tensor[0] << "_" << tensor[1] << "_" << tensor[2] << "_" << tensor[3] << "--"
                        << tensor[4] << "_" << tensor[5] << "_" << tensor[6] << "--"
                        << tensor[7] << "_" << tensor[8] << "_" << tensor[9] << "--"
                        << irho0 << "_" << tensor[10] << "_" << tensor[11] << "_" << tensor[12] << ".pn";
      ParseTBInputFile(proton_input_filename,tensor,SpSnSNexGamma_set,lsu3shell_basis,irho0,proton_tbdrmes_lsu3shell);
      ParseTBInputFile(neutron_input_filename,tensor,SpSnSNexGamma_set,lsu3shell_basis,irho0,neutron_tbdrmes_lsu3shell);
      ParseTBInputFile(pn_input_filename,tensor,SpSnSNexGamma_set,lsu3shell_basis,irho0,pn_tbdrmes_lsu3shell);
    } 
  }
  std::cout<<"Input and LSU3shell basis construction finished"<<std::endl;

  // Check dimension
  int dim1=0;
  for(std::map<std::tuple<Spins,int,SU3>,std::vector<Ipini>>::iterator it=lsu3shell_basis.begin(); it!=lsu3shell_basis.end(); it++)
     {dim1=dim1+it->second.size();}
  int dim2=0;
  for(std::map<Spins,std::map<int,std::vector<SU3>>>::iterator it=SpSnSNexGamma_set.begin(); it!=SpSnSNexGamma_set.end(); it++){
    for(std::map<int,std::vector<SU3>>::iterator it2=it->second.begin(); it2!=it->second.end(); it2++){
      for(SU3 x : it2->second){dim2=dim2+lsu3shell_basis[{it->first,it2->first,x}].size();}
    }
  }
  if(dim1!=dim2){
    std::cout<<"ERROR: dimensions: "<<dim1<<" "<<dim2<<std::endl;
    su3::finalize();
    return EXIT_FAILURE;
  }
  std::cout<<"Dimensions checked"<<std::endl;

  // Construct B
  std::map<std::tuple<Spins,int>,std::map<std::array<SU3,2>,Eigen::MatrixXd>> B;
  // For each Sp,Sn,S,Nex with Nex<=Nmax-1 there is a set of rectangular matrices - one matrix for each Gamma' within Sp,Sn,S,Nex+1
  // and each Gamma within Sp,Sn,S,Nex.
  for(std::map<Spins,std::map<int,std::vector<SU3>>>::iterator it=SpSnSNexGamma_set.begin(); it!=SpSnSNexGamma_set.end(); it++){
    for(int Nex=0; Nex<Nmax; Nex++){
      std::map<std::array<SU3,2>,Eigen::MatrixXd> B_by_Gammas;
      for(SU3 xp : it->second[Nex+1]){
	for(SU3 x : it->second[Nex]){
          Eigen::MatrixXd Btemp = Eigen::MatrixXd::Zero
		  (lsu3shell_basis[{it->first,Nex+1,xp}].size(),lsu3shell_basis[{it->first,Nex,x}].size());
          for(int i=0; i<lsu3shell_basis[{it->first,Nex+1,xp}].size(); i++){
            LSU3shellBasisState bra(lsu3shell_basis[{it->first,Nex+1,xp}][i],it->first,Nex+1,xp);
	    for(int j=0; j<lsu3shell_basis[{it->first,Nex,x}].size(); j++){
              LSU3shellBasisState ket(lsu3shell_basis[{it->first,Nex,x}][j],it->first,Nex,x);
	      for(int nu=0; nu<=Nmax+std::max(N1vp,N1vn)-1; nu++){
		double obdrmes;
                std::map<std::tuple<LSU3shellBasisState,LSU3shellBasisState,std::array<int,5>>,std::vector<double>>::iterator
			iter=proton_obdrmes_lsu3shell.find({bra,ket,{nu+1,nu,1,0,0}});
		if(iter==proton_obdrmes_lsu3shell.end()){
		  obdrmes=0.0;
		}else{
		  obdrmes=proton_obdrmes_lsu3shell[{bra,ket,{nu+1,nu,1,0,0}}][0];
	        }
		iter=neutron_obdrmes_lsu3shell.find({bra,ket,{nu+1,nu,1,0,0}});
                if(iter!=neutron_obdrmes_lsu3shell.end())obdrmes+=neutron_obdrmes_lsu3shell[{bra,ket,{nu+1,nu,1,0,0}}][0];
		Btemp(i,j)+=sqrt(double((nu+1)*(nu+2)*(nu+3))/3.0)*obdrmes;
              }
            }
          }
	  B_by_Gammas[{xp,x}]=Btemp;
	}
      }
      B[{it->first,Nex}]=B_by_Gammas;
    }
  }
  std::cout<<"B constructed"<<std::endl;
  
  // Construct A
  std::map<std::tuple<Spins,int>,std::map<std::array<SU3,2>,Eigen::MatrixXd>> A;
  // For each Sp,Sn,S,Nex, where Nex<=Nmax-2 is even, there is a set of rectangular matrices - one matrix for each Gamma'
  // within Sp,Sn,S,Nex+2 and each Gamma within Sp,Sn,S,Nex.
  for(std::map<Spins,std::map<int,std::vector<SU3>>>::iterator it=SpSnSNexGamma_set.begin(); it!=SpSnSNexGamma_set.end(); it++){
    for(int Nex=0; Nex<=Nmax-2; Nex=Nex+2){
      std::map<std::array<SU3,2>,Eigen::MatrixXd> A_by_Gammas;
      for(SU3 xp : it->second[Nex+2]){
        for(SU3 x : it->second[Nex]){
          Eigen::MatrixXd Atemp = Eigen::MatrixXd::Zero
                  (lsu3shell_basis[{it->first,Nex+2,xp}].size(),lsu3shell_basis[{it->first,Nex,x}].size());
          for(int i=0; i<lsu3shell_basis[{it->first,Nex+2,xp}].size(); i++){
            LSU3shellBasisState bra(lsu3shell_basis[{it->first,Nex+2,xp}][i],it->first,Nex+2,xp);
            for(int j=0; j<lsu3shell_basis[{it->first,Nex,x}].size(); j++){
              LSU3shellBasisState ket(lsu3shell_basis[{it->first,Nex,x}][j],it->first,Nex,x);
              for(int nu=0; nu<=Nmax+std::max(N1vp,N1vn)-2; nu++){
                double obdrmes;
                std::map<std::tuple<LSU3shellBasisState,LSU3shellBasisState,std::array<int,5>>,std::vector<double>>::iterator
                        iter=proton_obdrmes_lsu3shell.find({bra,ket,{nu+2,nu,2,0,0}});
                if(iter==proton_obdrmes_lsu3shell.end()){
                  obdrmes=0.0;
                }else{
                  obdrmes=proton_obdrmes_lsu3shell[{bra,ket,{nu+2,nu,2,0,0}}][0];
                }
                iter=neutron_obdrmes_lsu3shell.find({bra,ket,{nu+2,nu,2,0,0}});
                if(iter!=neutron_obdrmes_lsu3shell.end())obdrmes+=neutron_obdrmes_lsu3shell[{bra,ket,{nu+2,nu,2,0,0}}][0];
                Atemp(i,j)+=sqrt(double((nu+1)*(nu+2)*(nu+3)*(nu+4))/12.0)*obdrmes;
              }
            }
          }
          for(std::tuple<SU3,int> xpp : KroneckerProduct(x.lm(),x.mu(),1,0)){
            if(su3::mult(std::get<0>(xpp).lm(),std::get<0>(xpp).mu(),1,0,xp.lm(),xp.mu())==0)continue;
            std::map<std::array<SU3,2>,Eigen::MatrixXd>::iterator iter=B[{it->first,Nex+1}].find({xp,std::get<0>(xpp)});
            if(iter!=B[{it->first,Nex+1}].end()){
	      iter=B[{it->first,Nex}].find({std::get<0>(xpp),x});
	      if(iter!=B[{it->first,Nex}].end()){
                Atemp-=(U(x.lm(),x.mu(),1,0,xp.lm(),xp.mu(),1,0,std::get<0>(xpp).lm(),std::get<0>(xpp).mu(),1,1,2,0,1,1)
	               /(sqrt(2.0)*double(AA)))*B[{it->first,Nex+1}][{xp,std::get<0>(xpp)}]*B[{it->first,Nex}][{std::get<0>(xpp),x}];
              }
            }
          }
          A_by_Gammas[{xp,x}]=Atemp;
        }
      }
      A[{it->first,Nex}]=A_by_Gammas;
    }
  }
  std::cout<<"A constructed"<<std::endl;

  // Construct LGIs
  std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>>>> nonorthogonal_basis;
  // Wave function of non-orthogonal basis state is
  // nonorthogonal_basis[{gamma,Nex_sigma,Gamma_sigma,SpSnS}][Nex_omega][Gamma_omega][{Gamma_n,rho}]
  SU3 n00(0,0);
  for(std::map<Spins,std::map<int,std::vector<SU3>>>::iterator it=SpSnSNexGamma_set.begin(); it!=SpSnSNexGamma_set.end(); it++){
    for(SU3 x : it->second[0]){
      int gammamax=lsu3shell_basis[{it->first,0,x}].size();
      for(int gamma=1; gamma<=gammamax; gamma++){
	Eigen::VectorXd wf=Eigen::VectorXd::Zero(gammamax);
	wf(gamma-1)=1.0;
	std::map<std::tuple<SU3,int>,Eigen::VectorXd> wf_by_nrho;
	wf_by_nrho[{n00,1}]=wf;
	std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>> nrho_by_Gamma;
	nrho_by_Gamma[x]=wf_by_nrho;
	std::map<int,std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>>> Gamma_by_Nex;
	Gamma_by_Nex[0]=nrho_by_Gamma;
       	nonorthogonal_basis[{gamma,0,x,it->first}]=Gamma_by_Nex;
      }
    }
  }
  std::cout<<"LGIs constructed"<<std::endl;
/*
  std::cout<<"LGIs:"<<std::endl;
  std::cout<<"gamma N_sigma lm_sigma mu_sigma SSp Sn SS n_lambda n_mu rho N_omega lm_omega mu_omega wave_function"<<std::endl;
  for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>>>>::iterator
		  it=nonorthogonal_basis.begin(); it!=nonorthogonal_basis.end(); it++){
    for(std::map<int,std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>>>::iterator
		    it2=it->second.begin(); it2!=it->second.end(); it2++){
      for(std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>>::iterator it3=it2->second.begin(); it3!=it2->second.end(); it3++){
        for(std::map<std::tuple<SU3,int>,Eigen::VectorXd>::iterator it4=it3->second.begin(); it4!=it3->second.end(); it4++){
          std::cout<<std::get<0>(it->first)<<" "
		   <<std::get<1>(it->first)<<" "
		   <<std::get<2>(it->first).lm()<<" "
		   <<std::get<2>(it->first).mu()<<" "
		   <<std::get<3>(it->first).SSp()<<" "
		   <<std::get<3>(it->first).SSn()<<" "
		   <<std::get<3>(it->first).SS()<<" "
		   <<std::get<0>(it4->first).lm()<<" "
		   <<std::get<0>(it4->first).mu()<<" "
		   <<std::get<1>(it4->first)<<" "
		   <<it2->first<<" "
		   <<it3->first.lm()<<" "
		   <<it3->first.mu()<<" "
		   <<it4->second<<std::endl;
        }
      }
    }
  }
*/
  // Construct non-orthogonal SpNCCI basis
  for(int Nex=2; Nex<=Nmax; Nex=Nex+2){ // Nex=Nn because N_sigma,max=0.
    int Nexp=Nex-2;
    // loop over Sp,Sn,S
    for(std::map<Spins,std::map<int,std::vector<SU3>>>::iterator it=SpSnSNexGamma_set.begin(); it!=SpSnSNexGamma_set.end(); it++){
      // loop over sigma
      for(SU3 x_sigma : it->second[0]){
        int gammamax=lsu3shell_basis[{it->first,0,x_sigma}].size();
	// loop over gamma
        for(int gamma=1; gamma<=gammamax; gamma++){
          std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>> nrho_by_Gamma;
          // construct vector of possible omega values
	  std::vector<SU3> x_set;
	  for(int n1=0; n1<=Nex; n1=n1+2){
            for(int n2=0; n2<=std::min(n1,Nex-n1); n2=n2+2){
              int n3=Nex-n1-n2;
	      if(n3>n2)continue;
	      int lm_n=n1-n2;
	      int mu_n=n2-n3;
	      for(std::tuple<SU3,int> x : KroneckerProduct(x_sigma.lm(),x_sigma.mu(),lm_n,mu_n)){
                auto iter=std::find(x_set.begin(), x_set.end(), std::get<0>(x));
	        if(iter==x_set.end())x_set.push_back(std::get<0>(x));
	      }
            }
	  }
	  // loop over omega
          for(SU3 x : x_set){
            std::map<std::tuple<SU3,int>,Eigen::VectorXd> wf_by_nrho;
	    // loop over n
	    for(int n1=0; n1<=Nex; n1=n1+2){
              for(int n2=0; n2<=std::min(n1,Nex-n1); n2=n2+2){
                int n3=Nex-n1-n2;
                if(n3>n2)continue;
                int lm_n=n1-n2;
                int mu_n=n2-n3;
		SU3 n(lm_n,mu_n);
		// loop over rho
		for(int rho=1; rho<=su3::mult(x_sigma.lm(),x_sigma.mu(),lm_n,mu_n,x.lm(),x.mu()); rho++){
		  // Calculate wave function
                  Eigen::VectorXd wf=Eigen::VectorXd::Zero(lsu3shell_basis[{it->first,Nex,x}].size());
                  // sum over n'
		  for(int np1=0; np1<=Nexp; np1=np1+2){
		    for(int np2=0; np2<=std::min(np1,Nexp-np1); np2=np2+2){
                      int np3=Nexp-np1-np2;
		      if(np3>np2)continue;
		      int lm_np=np1-np2;
		      int mu_np=np2-np3;
		      SU3 np(lm_np,mu_np);
		      double aME=BosonCreationOperatorME(n1,n2,n3,np1,np2,np3);
		      // sum over omega'
		      for(std::tuple<SU3,int> xp : KroneckerProduct(x_sigma.lm(),x_sigma.mu(),lm_np,mu_np)){
                        // sum over rho'
                        for(int rhop=1; rhop<=std::get<1>(xp); rhop++){
                          double factor=aME*U(2,0,lm_np,mu_np,x.lm(),x.mu(),x_sigma.lm(),x_sigma.mu(),lm_n,mu_n,1,rho,
					  std::get<0>(xp).lm(),std::get<0>(xp).mu(),rhop,1);
			  int phase=std::get<0>(xp).lm()+std::get<0>(xp).mu()-x.lm()-x.mu();
			  if((phase/2)*2!=phase)factor=-factor;
			  wf+=factor*A[{it->first,Nexp}][{x,std::get<0>(xp)}]
				  *nonorthogonal_basis[{gamma,0,x_sigma,it->first}][Nexp][std::get<0>(xp)][{np,rhop}];
			}
                      }
		    }
		  }
		  wf_by_nrho[{n,rho}]=2.0*wf/double(Nex);
	        }
              }
            }
            nrho_by_Gamma[x]=wf_by_nrho;
	  }
	  nonorthogonal_basis[{gamma,0,x_sigma,it->first}][Nex]=nrho_by_Gamma;
        }
      }
    }
  }
  std::cout<<"Non-orthogonal SpNCCI basis constructed"<<std::endl;

  // Construct K-matrices and inverse K-matrices
  std::map<std::tuple<int,int,SU3,Spins,int,SU3>,Eigen::MatrixXd> K,Kinv;
  for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>>>>::iterator
		  it1=nonorthogonal_basis.begin(); it1!=nonorthogonal_basis.end(); it1++){
    for(int Nex=0; Nex<=Nmax; Nex=Nex+2){
      for(std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>>::iterator
		      it2=it1->second[Nex].begin(); it2!=it1->second[Nex].end(); it2++){
	int dim=it2->second.size();
        Eigen::MatrixXd KK(dim,dim);
	int bra=-1;
	for(std::map<std::tuple<SU3,int>,Eigen::VectorXd>::iterator itbra=it2->second.begin(); itbra!=it2->second.end(); itbra++){
          bra++;
	  int ket=-1;
          for(std::map<std::tuple<SU3,int>,Eigen::VectorXd>::iterator itket=it2->second.begin(); itket!=it2->second.end(); itket++){
            ket++;
            KK(bra,ket)=itbra->second.transpose()*itket->second;
          }
        }
	std::tuple<int,int,SU3,Spins,int,SU3> key
          ={std::get<0>(it1->first),std::get<1>(it1->first),std::get<2>(it1->first),std::get<3>(it1->first),Nex,it2->first};
        Eigen::MatrixXd Ktmp(dim,dim);
        if(dim==1){
          Ktmp(0,0)=sqrt(KK(0,0));
          K[key]=Ktmp;
          Ktmp(0,0)=1.0/Ktmp(0,0);
          Kinv[key]=Ktmp;
        }else{
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(KK);
          if(eigensolver.info()!=Eigen::Success){
            std::cout<<"ERROR: Eigen failed to diagonalize Hermitian square of K!"<<std::endl;
            su3::finalize();
            return EXIT_FAILURE;
          }
          Eigen::MatrixXd D = Eigen::MatrixXd::Zero(dim,dim);
          for(int ii=0; ii<dim; ii++){
            D(ii,ii)=sqrt(eigensolver.eigenvalues()[ii]);
          }
          Ktmp=eigensolver.eigenvectors()*D;
          K[key]=Ktmp;
          Kinv[key]=Ktmp.inverse();
        }
      }
    }
  }
  std::cout<<"K-matrices and inverse K-matrices constructed"<<std::endl;

  // Print K-matrices
  std::cout<<"K-matrices:"<<std::endl;
  for(std::map<std::tuple<int,int,SU3,Spins,int,SU3>,Eigen::MatrixXd>::iterator it=K.begin(); it!=K.end(); it++){
/*
if(it->second.cols()>1){
Eigen::MatrixXd Ktmp(it->second.rows(),it->second.cols());
Ktmp=it->second.transpose();
Ktmp=-Ktmp;
it->second=Ktmp;
Kinv[it->first]=Ktmp.inverse();
}
*/
    std::cout<<"sigma omega: "<<std::get<1>(it->first)<<"("<<std::get<2>(it->first).lm()<<","<<std::get<2>(it->first).mu()<<") "
	     <<std::get<3>(it->first).SSp()<<" "<<std::get<3>(it->first).SSn()<<" "<<std::get<3>(it->first).SS()<<"     "
	     <<std::get<4>(it->first)<<"("<<std::get<5>(it->first).lm()<<","<<std::get<5>(it->first).mu()<<")"<<std::endl;
    std::cout<<it->second<<std::endl;
//    std::cout<<it->second.inverse()<<std::endl;
  }

  // Construct orthonormal SpNCCI basis
  std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>> orthonormal_basis;
  // Wave function of orthonormal basis state is
  // orthonormal_basis[{gamma,Nex_sigma,Gamma_sigma,SpSnS}][Nex_omega][Gamma_omega][upsilon]
  for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>>>>::iterator
                  it1=nonorthogonal_basis.begin(); it1!=nonorthogonal_basis.end(); it1++){
    std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>> map1;
    for(int Nex=0; Nex<=Nmax; Nex=Nex+2){
      std::map<SU3,std::map<int,Eigen::VectorXd>> map2;
      for(std::map<SU3,std::map<std::tuple<SU3,int>,Eigen::VectorXd>>::iterator
                      it2=it1->second[Nex].begin(); it2!=it1->second[Nex].end(); it2++){
	std::map<int,Eigen::VectorXd> map3;
        int dimK=it2->second.size();
	std::tuple<int,int,SU3,Spins,int,SU3> key={std::get<0>(it1->first),std::get<1>(it1->first),std::get<2>(it1->first),
		std::get<3>(it1->first),Nex,it2->first};
        for(int upsilon=1; upsilon<=dimK; upsilon++){
      	  Eigen::VectorXd wf=Eigen::VectorXd::Zero(lsu3shell_basis[{std::get<3>(it1->first),Nex,it2->first}].size());
	  int ind=-1;
	  for(std::map<std::tuple<SU3,int>,Eigen::VectorXd>::iterator it3=it2->second.begin(); it3!=it2->second.end(); it3++){
            ind++;
            wf+=Kinv[key](upsilon-1,ind)*it3->second;
	  }
	  map3[upsilon]=wf;
        }
	map2[it2->first]=map3;
      }
      map1[Nex]=map2;
    }
    orthonormal_basis[it1->first]=map1;
  }
  std::cout<<"Orthonormal SpNCCI basis constructed"<<std::endl;

  // Check orthonormality of SpNCCI basis
  for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>::iterator
                  it1bra=orthonormal_basis.begin(); it1bra!=orthonormal_basis.end(); it1bra++){
    for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>::iterator
                  it1ket=orthonormal_basis.begin(); it1ket!=orthonormal_basis.end(); it1ket++){
      if(!(std::get<3>(it1bra->first)==std::get<3>(it1ket->first)))continue;
      for(int Nex=0; Nex<=Nmax; Nex=Nex+2){
        for(std::map<SU3,std::map<int,Eigen::VectorXd>>::iterator
		      it2bra=it1bra->second[Nex].begin(); it2bra!=it1bra->second[Nex].end(); it2bra++){
          std::map<SU3,std::map<int,Eigen::VectorXd>>::iterator it2ket=it1ket->second[Nex].find(it2bra->first);
	  if(it2ket!=it1ket->second[Nex].end()){
	    for(std::map<int,Eigen::VectorXd>::iterator it3bra=it2bra->second.begin(); it3bra!=it2bra->second.end(); it3bra++){
   	      for(std::map<int,Eigen::VectorXd>::iterator it3ket=it1ket->second[Nex][it2bra->first].begin();
			      it3ket!=it1ket->second[Nex][it2bra->first].end(); it3ket++){
                double scalprod=it3bra->second.transpose()*it3ket->second;
	        double val;
                if(it1bra->first==it1ket->first && it3bra->first==it3ket->first){
	          val=1.0;
	        }else{
                  val=0.0;
                }
	        if(std::abs(scalprod-val)>1.0e-5){
                  std::cout<<"ERROR: Scalar product "<<scalprod<<" should be "<<val<<std::endl;
                  su3::finalize();
                  return EXIT_FAILURE;
	        }
              }
            }
	  }
        }
      }
    }
  }
  std::cout<<"Orthonormality of SpNCCI basis checked"<<std::endl;

  std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>> proton_tbdrmes,neutron_tbdrmes,pn_tbdrmes;

  // Explicit construction of TBDRMEs in SpNCCI basis
  for(const std::array<int,13>& tens : TBtensors){
    for(int irho0=0; irho0<su3::mult(tens[4], tens[5], tens[7], tens[8], tens[10], tens[11]); irho0++){
      std::array<int,14> tensor={tens[0],tens[1],tens[2],tens[3],tens[4],tens[5],tens[6],tens[7],tens[8],
	                         tens[9],tens[10],tens[11],tens[12],irho0};
      for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>::iterator
                  it1bra=orthonormal_basis.begin(); it1bra!=orthonormal_basis.end(); it1bra++){
        for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>::iterator
                  it1ket=orthonormal_basis.begin(); it1ket!=orthonormal_basis.end(); it1ket++){
	  if(std::abs(std::get<3>(it1ket->first).SS()-tensor[12])>std::get<3>(it1bra->first).SS()
	     || std::get<3>(it1bra->first).SS()>std::get<3>(it1ket->first).SS()+tensor[12])continue;
          for(int Nexbra=0; Nexbra<=Nmax; Nexbra=Nexbra+2){
            int Nexket=Nexbra-tensor[0]-tensor[1]+tensor[2]+tensor[3];
	    if(Nexket<0)continue;
	    if(Nexket>Nmax)break;
            for(std::map<SU3,std::map<int,Eigen::VectorXd>>::iterator
                      it2bra=it1bra->second[Nexbra].begin(); it2bra!=it1bra->second[Nexbra].end(); it2bra++){
	      for(std::map<SU3,std::map<int,Eigen::VectorXd>>::iterator
                      it2ket=it1ket->second[Nexket].begin(); it2ket!=it1ket->second[Nexket].end(); it2ket++){
	        int rhomax=su3::mult(it2ket->first.lm(),it2ket->first.mu(),tensor[10],tensor[11],it2bra->first.lm(),it2bra->first.mu());
	        if(rhomax==0)continue;
                for(std::map<int,Eigen::VectorXd>::iterator it3bra=it2bra->second.begin(); it3bra!=it2bra->second.end(); it3bra++){
		  ONBasisState bra(std::get<0>(it1bra->first),std::get<1>(it1bra->first),std::get<2>(it1bra->first),
		  		   std::get<3>(it1bra->first),it3bra->first,Nexbra,it2bra->first);
		  for(std::map<int,Eigen::VectorXd>::iterator it3ket=it2ket->second.begin(); it3ket!=it2ket->second.end(); it3ket++){
	            ONBasisState ket(std::get<0>(it1ket->first),std::get<1>(it1ket->first),std::get<2>(it1ket->first),
                                     std::get<3>(it1ket->first),it3ket->first,Nexket,it2ket->first);
                    std::vector<double> proton_rmes(rhomax,0.0);
	            std::vector<double> neutron_rmes(rhomax,0.0);
		    std::vector<double> pn_rmes(rhomax,0.0);
		    int ind_bra=-1;
		    for(Ipini ipini_bra : lsu3shell_basis[{std::get<3>(it1bra->first),Nexbra,it2bra->first}]){
		      ind_bra++;
                      LSU3shellBasisState lsu3shell_bra(ipini_bra,std::get<3>(it1bra->first),Nexbra,it2bra->first);
		      int ind_ket=-1;
		      for(Ipini ipini_ket : lsu3shell_basis[{std::get<3>(it1ket->first),Nexket,it2ket->first}]){
		        ind_ket++;
                        LSU3shellBasisState lsu3shell_ket(ipini_ket,std::get<3>(it1ket->first),Nexket,it2ket->first);
		        std::map<std::tuple<LSU3shellBasisState,LSU3shellBasisState,std::array<int,14>>,std::vector<double>>::iterator
			        iter=proton_tbdrmes_lsu3shell.find({lsu3shell_bra,lsu3shell_ket,tensor});
		        if(iter!=proton_tbdrmes_lsu3shell.end()){
			  for(int irho=0; irho<rhomax; irho++){
                            proton_rmes[irho]+=it3bra->second[ind_bra]*proton_tbdrmes_lsu3shell[{lsu3shell_bra,lsu3shell_ket,tensor}][irho]
                                  *it3ket->second[ind_ket];
                          }
		        }
		        iter=neutron_tbdrmes_lsu3shell.find({lsu3shell_bra,lsu3shell_ket,tensor});
		        if(iter!=neutron_tbdrmes_lsu3shell.end()){
			  for(int irho=0; irho<rhomax; irho++){
                            neutron_rmes[irho]+=it3bra->second[ind_bra]*neutron_tbdrmes_lsu3shell[{lsu3shell_bra,lsu3shell_ket,tensor}][irho]
                                *it3ket->second[ind_ket];
                          }
		        }
			iter=pn_tbdrmes_lsu3shell.find({lsu3shell_bra,lsu3shell_ket,tensor});
                        if(iter!=pn_tbdrmes_lsu3shell.end()){
                          for(int irho=0; irho<rhomax; irho++){
                            pn_rmes[irho]+=it3bra->second[ind_bra]*pn_tbdrmes_lsu3shell[{lsu3shell_bra,lsu3shell_ket,tensor}][irho]
                                *it3ket->second[ind_ket];
                          }
                        }
                      }
                    }
                    proton_tbdrmes[{bra,ket,tensor}]=proton_rmes;
		    neutron_tbdrmes[{bra,ket,tensor}]=neutron_rmes;
		    pn_tbdrmes[{bra,ket,tensor}]=pn_rmes;
                  }
                }
	      }
            }
          }
        }
      }
    } 
  }
  std::cout<<"Explicit construction of TBDRMEs in SpNCCI basis finished"<<std::endl;

  // Output explicitly constructed TBDRMEs
  std::ofstream pfile("explicit_tbdrmes_proton.dat");
  if (!pfile) {
    std::cout << "Could not open file explicit_tbdrmes_proton.dat" << std::endl;
    return EXIT_FAILURE;
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                  it=proton_tbdrmes.begin(); it!=proton_tbdrmes.end(); it++){
    pfile<<std::get<0>(it->first).gamma()<<" "
          <<std::get<0>(it->first).Nex_sigma()<<" "<<std::get<0>(it->first).x_sigma().lm()<<" "<<std::get<0>(it->first).x_sigma().mu()<<" "
          <<std::get<0>(it->first).SpSnS().SSp()<<" "<<std::get<0>(it->first).SpSnS().SSn()<<" "<<std::get<0>(it->first).SpSnS().SS()<<" "
          <<std::get<0>(it->first).upsilon()<<" "<<std::get<0>(it->first).Nex_omega()<<" "<<std::get<0>(it->first).x_omega().lm()<<" "
          <<std::get<0>(it->first).x_omega().mu()<<" "<<std::get<1>(it->first).gamma()<<" "
          <<std::get<1>(it->first).Nex_sigma()<<" "<<std::get<1>(it->first).x_sigma().lm()<<" "<<std::get<1>(it->first).x_sigma().mu()<<" "
          <<std::get<1>(it->first).SpSnS().SSp()<<" "<<std::get<1>(it->first).SpSnS().SSn()<<" "<<std::get<1>(it->first).SpSnS().SS()<<" "
          <<std::get<1>(it->first).upsilon()<<" "<<std::get<1>(it->first).Nex_omega()<<" "<<std::get<1>(it->first).x_omega().lm()<<" "
          <<std::get<1>(it->first).x_omega().mu()<<" "
	  <<std::get<2>(it->first)[0]<<" "<<std::get<2>(it->first)[1]<<" "<<std::get<2>(it->first)[2]<<" "
          <<std::get<2>(it->first)[3]<<" "<<std::get<2>(it->first)[4]<<" "<<std::get<2>(it->first)[5]<<" "
	  <<std::get<2>(it->first)[6]<<" "<<std::get<2>(it->first)[7]<<" "<<std::get<2>(it->first)[8]<<" "
	  <<std::get<2>(it->first)[9]<<" "<<std::get<2>(it->first)[10]<<" "<<std::get<2>(it->first)[11]<<" "
	  <<std::get<2>(it->first)[12]<<" "<<std::get<2>(it->first)[13]+1<<" "<<it->second.size();
    for(double rme : it->second){
      pfile<<" "<<rme;
    }
    pfile<<std::endl;
  }
  std::ofstream nfile("explicit_tbdrmes_neutron.dat");
  if (!nfile) {
    std::cout << "Could not open file explicit_tbdrmes_neutron.dat" << std::endl;
    return EXIT_FAILURE;
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                  it=neutron_tbdrmes.begin(); it!=neutron_tbdrmes.end(); it++){
    nfile<<std::get<0>(it->first).gamma()<<" "
          <<std::get<0>(it->first).Nex_sigma()<<" "<<std::get<0>(it->first).x_sigma().lm()<<" "<<std::get<0>(it->first).x_sigma().mu()<<" "
          <<std::get<0>(it->first).SpSnS().SSp()<<" "<<std::get<0>(it->first).SpSnS().SSn()<<" "<<std::get<0>(it->first).SpSnS().SS()<<" "
          <<std::get<0>(it->first).upsilon()<<" "<<std::get<0>(it->first).Nex_omega()<<" "<<std::get<0>(it->first).x_omega().lm()<<" "
          <<std::get<0>(it->first).x_omega().mu()<<" "<<std::get<1>(it->first).gamma()<<" "
          <<std::get<1>(it->first).Nex_sigma()<<" "<<std::get<1>(it->first).x_sigma().lm()<<" "<<std::get<1>(it->first).x_sigma().mu()<<" "
          <<std::get<1>(it->first).SpSnS().SSp()<<" "<<std::get<1>(it->first).SpSnS().SSn()<<" "<<std::get<1>(it->first).SpSnS().SS()<<" "
          <<std::get<1>(it->first).upsilon()<<" "<<std::get<1>(it->first).Nex_omega()<<" "<<std::get<1>(it->first).x_omega().lm()<<" "
          <<std::get<1>(it->first).x_omega().mu()<<" "
          <<std::get<2>(it->first)[0]<<" "<<std::get<2>(it->first)[1]<<" "<<std::get<2>(it->first)[2]<<" "
          <<std::get<2>(it->first)[3]<<" "<<std::get<2>(it->first)[4]<<" "<<std::get<2>(it->first)[5]<<" "
	  <<std::get<2>(it->first)[6]<<" "<<std::get<2>(it->first)[7]<<" "<<std::get<2>(it->first)[8]<<" "
	  <<std::get<2>(it->first)[9]<<" "<<std::get<2>(it->first)[10]<<" "<<std::get<2>(it->first)[11]<<" "
	  <<std::get<2>(it->first)[12]<<" "<<std::get<2>(it->first)[13]+1<<" "<<it->second.size();
    for(double rme : it->second){
      nfile<<" "<<rme;
    }
    nfile<<std::endl;
  }
  std::ofstream pnfile("explicit_tbdrmes_pn.dat");
  if (!pnfile) {
    std::cout << "Could not open file explicit_tbdrmes_pn.dat" << std::endl;
    return EXIT_FAILURE;
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                  it=pn_tbdrmes.begin(); it!=pn_tbdrmes.end(); it++){
    pnfile<<std::get<0>(it->first).gamma()<<" "
          <<std::get<0>(it->first).Nex_sigma()<<" "<<std::get<0>(it->first).x_sigma().lm()<<" "<<std::get<0>(it->first).x_sigma().mu()<<" "
          <<std::get<0>(it->first).SpSnS().SSp()<<" "<<std::get<0>(it->first).SpSnS().SSn()<<" "<<std::get<0>(it->first).SpSnS().SS()<<" "
          <<std::get<0>(it->first).upsilon()<<" "<<std::get<0>(it->first).Nex_omega()<<" "<<std::get<0>(it->first).x_omega().lm()<<" "
          <<std::get<0>(it->first).x_omega().mu()<<" "<<std::get<1>(it->first).gamma()<<" "
          <<std::get<1>(it->first).Nex_sigma()<<" "<<std::get<1>(it->first).x_sigma().lm()<<" "<<std::get<1>(it->first).x_sigma().mu()<<" "
          <<std::get<1>(it->first).SpSnS().SSp()<<" "<<std::get<1>(it->first).SpSnS().SSn()<<" "<<std::get<1>(it->first).SpSnS().SS()<<" "
          <<std::get<1>(it->first).upsilon()<<" "<<std::get<1>(it->first).Nex_omega()<<" "<<std::get<1>(it->first).x_omega().lm()<<" "
          <<std::get<1>(it->first).x_omega().mu()<<" "
          <<std::get<2>(it->first)[0]<<" "<<std::get<2>(it->first)[1]<<" "<<std::get<2>(it->first)[2]<<" "
          <<std::get<2>(it->first)[3]<<" "<<std::get<2>(it->first)[4]<<" "<<std::get<2>(it->first)[5]<<" "
	  <<std::get<2>(it->first)[6]<<" "<<std::get<2>(it->first)[7]<<" "<<std::get<2>(it->first)[8]<<" "
	  <<std::get<2>(it->first)[9]<<" "<<std::get<2>(it->first)[10]<<" "<<std::get<2>(it->first)[11]<<" "
	  <<std::get<2>(it->first)[12]<<" "<<std::get<2>(it->first)[13]+1<<" "<<it->second.size();
    for(double rme : it->second){
      pnfile<<" "<<rme;
    }
    pnfile<<std::endl;
  }
  std::cout<<"Explicitly constructed TBDRMEs in SpNCCI basis written"<<std::endl;

/*
  // Read explicitly constructed TBDRMEs
  std::ifstream ipfile("explicit_tbdrmes_proton.dat");
  if(!ipfile){
    std::cout<<"Could not open explicit_tbdrmes_proton.dat!"<<std::endl;
    su3::finalize();
    exit(EXIT_FAILURE);
  }
  while(true){
    int gammap,Nex_sigmap,lm_sigmap,mu_sigmap,SSp_bra,SSn_bra,SS_bra,upsilonp,Nex_omegap,lm_omegap,mu_omegap,
	gamma,Nex_sigma,lm_sigma,mu_sigma,SSp_ket,SSn_ket,SS_ket,upsilon,Nex_omega,lm_omega,mu_omega,
	N1,N2,N3,N4,lmf,muf,SSf,lmi,mui,SSi,lm0,mu0,SS0,rho0,rhomax;
    ipfile>>gammap>>Nex_sigmap>>lm_sigmap>>mu_sigmap>>SSp_bra>>SSn_bra>>SS_bra>>upsilonp>>Nex_omegap>>lm_omegap>>mu_omegap
          >>gamma>>Nex_sigma>>lm_sigma>>mu_sigma>>SSp_ket>>SSn_ket>>SS_ket>>upsilon>>Nex_omega>>lm_omega>>mu_omega
          >>N1>>N2>>N3>>N4>>lmf>>muf>>SSf>>lmi>>mui>>SSi>>lm0>>mu0>>SS0>>rho0>>rhomax;
    if(ipfile){
      SU3 x_sigmap(lm_sigmap,mu_sigmap);
      Spins SpSnS_bra(SSp_bra,SSn_bra,SS_bra);
      SU3 x_omegap(lm_omegap,mu_omegap);
      SU3 x_sigma(lm_sigma,mu_sigma);
      Spins SpSnS_ket(SSp_ket,SSn_ket,SS_ket);
      SU3 x_omega(lm_omega,mu_omega);
      ONBasisState bra(gammap,Nex_sigmap,x_sigmap,SpSnS_bra,upsilonp,Nex_omegap,x_omegap);
      ONBasisState ket(gamma,Nex_sigma,x_sigma,SpSnS_ket,upsilon,Nex_omega,x_omega);
      std::array<int,14> tensor={N1,N2,N3,N4,lmf,muf,SSf,lmi,mui,SSi,lm0,mu0,SS0,rho0-1};
      std::vector<double> rmes(rhomax);
      for(int irho=0; irho<rhomax; irho++){
        ipfile>>rmes[irho];
      }
      proton_tbdrmes[{bra,ket,tensor}]=rmes;
    }else{
      break;
    }
  }
  std::ifstream infile("explicit_tbdrmes_neutron.dat");
  if(!infile){
    std::cout<<"Could not open explicit_tbdrmes_neutron.dat!"<<std::endl;
    su3::finalize();
    exit(EXIT_FAILURE);
  }
  while(true){
    int gammap,Nex_sigmap,lm_sigmap,mu_sigmap,SSp_bra,SSn_bra,SS_bra,upsilonp,Nex_omegap,lm_omegap,mu_omegap,
        gamma,Nex_sigma,lm_sigma,mu_sigma,SSp_ket,SSn_ket,SS_ket,upsilon,Nex_omega,lm_omega,mu_omega,
        N1,N2,N3,N4,lmf,muf,SSf,lmi,mui,SSi,lm0,mu0,SS0,rho0,rhomax;
    infile>>gammap>>Nex_sigmap>>lm_sigmap>>mu_sigmap>>SSp_bra>>SSn_bra>>SS_bra>>upsilonp>>Nex_omegap>>lm_omegap>>mu_omegap
          >>gamma>>Nex_sigma>>lm_sigma>>mu_sigma>>SSp_ket>>SSn_ket>>SS_ket>>upsilon>>Nex_omega>>lm_omega>>mu_omega
          >>N1>>N2>>N3>>N4>>lmf>>muf>>SSf>>lmi>>mui>>SSi>>lm0>>mu0>>SS0>>rho0>>rhomax;
    if(infile){
      SU3 x_sigmap(lm_sigmap,mu_sigmap);
      Spins SpSnS_bra(SSp_bra,SSn_bra,SS_bra);
      SU3 x_omegap(lm_omegap,mu_omegap);
      SU3 x_sigma(lm_sigma,mu_sigma);
      Spins SpSnS_ket(SSp_ket,SSn_ket,SS_ket);
      SU3 x_omega(lm_omega,mu_omega);
      ONBasisState bra(gammap,Nex_sigmap,x_sigmap,SpSnS_bra,upsilonp,Nex_omegap,x_omegap);
      ONBasisState ket(gamma,Nex_sigma,x_sigma,SpSnS_ket,upsilon,Nex_omega,x_omega);
      std::array<int,14> tensor={N1,N2,N3,N4,lmf,muf,SSf,lmi,mui,SSi,lm0,mu0,SS0,rho0-1};
      std::vector<double> rmes(rhomax);
      for(int irho=0; irho<rhomax; irho++){
        infile>>rmes[irho];
      }
      neutron_tbdrmes[{bra,ket,tensor}]=rmes;
    }else{
      break;
    }
  }
  std::ifstream ipnfile("explicit_tbdrmes_pn.dat");
  if(!ipnfile){
    std::cout<<"Could not open explicit_tbdrmes_pn.dat!"<<std::endl;
    su3::finalize();
    exit(EXIT_FAILURE);
  }
  while(true){
    int gammap,Nex_sigmap,lm_sigmap,mu_sigmap,SSp_bra,SSn_bra,SS_bra,upsilonp,Nex_omegap,lm_omegap,mu_omegap,
        gamma,Nex_sigma,lm_sigma,mu_sigma,SSp_ket,SSn_ket,SS_ket,upsilon,Nex_omega,lm_omega,mu_omega,
        N1,N2,N3,N4,lmf,muf,SSf,lmi,mui,SSi,lm0,mu0,SS0,rho0,rhomax;
    ipnfile>>gammap>>Nex_sigmap>>lm_sigmap>>mu_sigmap>>SSp_bra>>SSn_bra>>SS_bra>>upsilonp>>Nex_omegap>>lm_omegap>>mu_omegap
           >>gamma>>Nex_sigma>>lm_sigma>>mu_sigma>>SSp_ket>>SSn_ket>>SS_ket>>upsilon>>Nex_omega>>lm_omega>>mu_omega
           >>N1>>N2>>N3>>N4>>lmf>>muf>>SSf>>lmi>>mui>>SSi>>lm0>>mu0>>SS0>>rho0>>rhomax;
    if(ipnfile){
      SU3 x_sigmap(lm_sigmap,mu_sigmap);
      Spins SpSnS_bra(SSp_bra,SSn_bra,SS_bra);
      SU3 x_omegap(lm_omegap,mu_omegap);
      SU3 x_sigma(lm_sigma,mu_sigma);
      Spins SpSnS_ket(SSp_ket,SSn_ket,SS_ket);
      SU3 x_omega(lm_omega,mu_omega);
      ONBasisState bra(gammap,Nex_sigmap,x_sigmap,SpSnS_bra,upsilonp,Nex_omegap,x_omegap);
      ONBasisState ket(gamma,Nex_sigma,x_sigma,SpSnS_ket,upsilon,Nex_omega,x_omega);
      std::array<int,14> tensor={N1,N2,N3,N4,lmf,muf,SSf,lmi,mui,SSi,lm0,mu0,SS0,rho0-1};
      std::vector<double> rmes(rhomax);
      for(int irho=0; irho<rhomax; irho++){
        ipnfile>>rmes[irho];
      }
      pn_tbdrmes[{bra,ket,tensor}]=rmes;
    }else{
      break;
    }
  }
  std::cout<<"Explicitly constructed TBDRMEs read"<<std::endl;
*/
  // Check symmetries of explicitly constructed proton TBDRMEs
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
		  it=proton_tbdrmes.begin(); it!=proton_tbdrmes.end(); it++){
    int rho0max=su3::mult(std::get<2>(it->first)[4],std::get<2>(it->first)[5],std::get<2>(it->first)[7],
		          std::get<2>(it->first)[8],std::get<2>(it->first)[10],std::get<2>(it->first)[11]);
    double factor=sqrt(double((std::get<1>(it->first).SpSnS().SS()+1)
		  *su3dim(std::get<1>(it->first).x_omega().lm(),std::get<1>(it->first).x_omega().mu()))
		  /double((std::get<0>(it->first).SpSnS().SS()+1)
		  *su3dim(std::get<0>(it->first).x_omega().lm(),std::get<0>(it->first).x_omega().mu())));
    int phase=std::get<2>(it->first)[10]+std::get<2>(it->first)[11]+std::get<0>(it->first).x_omega().lm()
              +std::get<0>(it->first).x_omega().mu()-std::get<1>(it->first).x_omega().lm()-std::get<1>(it->first).x_omega().mu()
              +(std::get<1>(it->first).SpSnS().SS()-std::get<0>(it->first).SpSnS().SS())/2+rho0max
              +std::get<2>(it->first)[0]+std::get<2>(it->first)[1]+std::get<2>(it->first)[2]+std::get<2>(it->first)[3]
              +std::get<2>(it->first)[4]+std::get<2>(it->first)[5]+std::get<2>(it->first)[7]+std::get<2>(it->first)[8];
    if((phase/2)*2!=phase)factor=-factor;
    int irho=-1;
    for(double rme : it->second){
      irho++;
      double rhs=0.0;
      for(int rho0p=1; rho0p<=rho0max; rho0p++){
	if((rho0p/2)*2==rho0p){
          rhs+=Phi(std::get<2>(it->first)[13]+1,rho0p,std::get<2>(it->first)[4],std::get<2>(it->first)[5],                                                       std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[10],std::get<2>(it->first)[11])
	       *proton_tbdrmes[{std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[3],std::get<2>(it->first)[2],
                 std::get<2>(it->first)[1],std::get<2>(it->first)[0],std::get<2>(it->first)[8],std::get<2>(it->first)[7],
	         std::get<2>(it->first)[9],std::get<2>(it->first)[5],std::get<2>(it->first)[4],std::get<2>(it->first)[6],
	         std::get<2>(it->first)[11],std::get<2>(it->first)[10],std::get<2>(it->first)[12],rho0p-1}}][irho];
        }else{
          rhs-=Phi(std::get<2>(it->first)[13]+1,rho0p,std::get<2>(it->first)[4],std::get<2>(it->first)[5],                                                       std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[10],std::get<2>(it->first)[11])
               *proton_tbdrmes[{std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[3],std::get<2>(it->first)[2],
                 std::get<2>(it->first)[1],std::get<2>(it->first)[0],std::get<2>(it->first)[8],std::get<2>(it->first)[7],
                 std::get<2>(it->first)[9],std::get<2>(it->first)[5],std::get<2>(it->first)[4],std::get<2>(it->first)[6],
                 std::get<2>(it->first)[11],std::get<2>(it->first)[10],std::get<2>(it->first)[12],rho0p-1}}][irho];
	}
      }
      rhs=factor*rhs;
      if(std::abs(rme-rhs)>1.0e-5){
        std::cout<<"ERROR (symmetry of proton TBDRMEs): "<<rme<<"!="<<rhs<<std::endl;
        su3::finalize();
        return EXIT_FAILURE;
      }
      rhs=proton_tbdrmes[{std::get<0>(it->first),std::get<1>(it->first),{std::get<2>(it->first)[1],std::get<2>(it->first)[0],
                   std::get<2>(it->first)[2],std::get<2>(it->first)[3],std::get<2>(it->first)[4],std::get<2>(it->first)[5],
                   std::get<2>(it->first)[6],std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[9],
                   std::get<2>(it->first)[10],std::get<2>(it->first)[11],std::get<2>(it->first)[12],std::get<2>(it->first)[13]}}][irho];
      int faza=std::get<2>(it->first)[0]+std::get<2>(it->first)[1]+std::get<2>(it->first)[4]+std::get<2>(it->first)[5]
               +std::get<2>(it->first)[6]/2;
      if((faza/2)*2!=faza)rhs=-rhs;
      if(std::abs(rme-rhs)>1.0e-5){
        std::cout<<"ERROR (1st little symmetry of proton TBDRMEs): "<<rme<<"!="<<rhs<<std::endl;
        su3::finalize();
        return EXIT_FAILURE;
      }
      rhs=proton_tbdrmes[{std::get<0>(it->first),std::get<1>(it->first),{std::get<2>(it->first)[0],std::get<2>(it->first)[1],
                   std::get<2>(it->first)[3],std::get<2>(it->first)[2],std::get<2>(it->first)[4],std::get<2>(it->first)[5],
                   std::get<2>(it->first)[6],std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[9],
                   std::get<2>(it->first)[10],std::get<2>(it->first)[11],std::get<2>(it->first)[12],std::get<2>(it->first)[13]}}][irho];
      faza=std::get<2>(it->first)[2]+std::get<2>(it->first)[3]+std::get<2>(it->first)[7]+std::get<2>(it->first)[8]
               +std::get<2>(it->first)[9]/2;
      if((faza/2)*2!=faza)rhs=-rhs;
      if(std::abs(rme-rhs)>1.0e-5){
        std::cout<<"ERROR (2nd little symmetry of proton TBDRMEs): "<<rme<<"!="<<rhs<<std::endl;
        su3::finalize();
        return EXIT_FAILURE;
      }
    }
  }
  // Check symmetries of explicitly constructed neutron TBDRMEs
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
		  it=neutron_tbdrmes.begin(); it!=neutron_tbdrmes.end(); it++){
    int rho0max=su3::mult(std::get<2>(it->first)[4],std::get<2>(it->first)[5],std::get<2>(it->first)[7],
		          std::get<2>(it->first)[8],std::get<2>(it->first)[10],std::get<2>(it->first)[11]);
    double factor=sqrt(double((std::get<1>(it->first).SpSnS().SS()+1)
		  *su3dim(std::get<1>(it->first).x_omega().lm(),std::get<1>(it->first).x_omega().mu()))
		  /double((std::get<0>(it->first).SpSnS().SS()+1)
		  *su3dim(std::get<0>(it->first).x_omega().lm(),std::get<0>(it->first).x_omega().mu())));
    int phase=std::get<2>(it->first)[10]+std::get<2>(it->first)[11]+std::get<0>(it->first).x_omega().lm()
              +std::get<0>(it->first).x_omega().mu()-std::get<1>(it->first).x_omega().lm()-std::get<1>(it->first).x_omega().mu()
              +(std::get<1>(it->first).SpSnS().SS()-std::get<0>(it->first).SpSnS().SS())/2+rho0max
              +std::get<2>(it->first)[0]+std::get<2>(it->first)[1]+std::get<2>(it->first)[2]+std::get<2>(it->first)[3]
              +std::get<2>(it->first)[4]+std::get<2>(it->first)[5]+std::get<2>(it->first)[7]+std::get<2>(it->first)[8];
    if((phase/2)*2!=phase)factor=-factor;
    int irho=-1;
    for(double rme : it->second){
      irho++;
      double rhs=0.0;
      for(int rho0p=1; rho0p<=rho0max; rho0p++){
	if((rho0p/2)*2==rho0p){
          rhs+=Phi(std::get<2>(it->first)[13]+1,rho0p,std::get<2>(it->first)[4],std::get<2>(it->first)[5],                                                       std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[10],std::get<2>(it->first)[11])
	       *neutron_tbdrmes[{std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[3],std::get<2>(it->first)[2],
                 std::get<2>(it->first)[1],std::get<2>(it->first)[0],std::get<2>(it->first)[8],std::get<2>(it->first)[7],
	         std::get<2>(it->first)[9],std::get<2>(it->first)[5],std::get<2>(it->first)[4],std::get<2>(it->first)[6],
	         std::get<2>(it->first)[11],std::get<2>(it->first)[10],std::get<2>(it->first)[12],rho0p-1}}][irho];
        }else{
          rhs-=Phi(std::get<2>(it->first)[13]+1,rho0p,std::get<2>(it->first)[4],std::get<2>(it->first)[5],                                                       std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[10],std::get<2>(it->first)[11])
               *neutron_tbdrmes[{std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[3],std::get<2>(it->first)[2],
                 std::get<2>(it->first)[1],std::get<2>(it->first)[0],std::get<2>(it->first)[8],std::get<2>(it->first)[7],
                 std::get<2>(it->first)[9],std::get<2>(it->first)[5],std::get<2>(it->first)[4],std::get<2>(it->first)[6],
                 std::get<2>(it->first)[11],std::get<2>(it->first)[10],std::get<2>(it->first)[12],rho0p-1}}][irho];
	}
      }
      rhs=factor*rhs;
      if(std::abs(rme-rhs)>1.0e-5){
        std::cout<<"ERROR (symmetry of neutron TBDRMEs): "<<rme<<"!="<<rhs<<std::endl;
        su3::finalize();
        return EXIT_FAILURE;
      }
      rhs=neutron_tbdrmes[{std::get<0>(it->first),std::get<1>(it->first),{std::get<2>(it->first)[1],std::get<2>(it->first)[0],
                   std::get<2>(it->first)[2],std::get<2>(it->first)[3],std::get<2>(it->first)[4],std::get<2>(it->first)[5],
                   std::get<2>(it->first)[6],std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[9],
                   std::get<2>(it->first)[10],std::get<2>(it->first)[11],std::get<2>(it->first)[12],std::get<2>(it->first)[13]}}][irho];
      int faza=std::get<2>(it->first)[0]+std::get<2>(it->first)[1]+std::get<2>(it->first)[4]+std::get<2>(it->first)[5]
               +std::get<2>(it->first)[6]/2;
      if((faza/2)*2!=faza)rhs=-rhs;
      if(std::abs(rme-rhs)>1.0e-5){
        std::cout<<"ERROR (1st little symmetry of neutron TBDRMEs): "<<rme<<"!="<<rhs<<std::endl;
        su3::finalize();
        return EXIT_FAILURE;
      }
      rhs=neutron_tbdrmes[{std::get<0>(it->first),std::get<1>(it->first),{std::get<2>(it->first)[0],std::get<2>(it->first)[1],
                   std::get<2>(it->first)[3],std::get<2>(it->first)[2],std::get<2>(it->first)[4],std::get<2>(it->first)[5],
                   std::get<2>(it->first)[6],std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[9],
                   std::get<2>(it->first)[10],std::get<2>(it->first)[11],std::get<2>(it->first)[12],std::get<2>(it->first)[13]}}][irho];
      faza=std::get<2>(it->first)[2]+std::get<2>(it->first)[3]+std::get<2>(it->first)[7]+std::get<2>(it->first)[8]
               +std::get<2>(it->first)[9]/2;
      if((faza/2)*2!=faza)rhs=-rhs;
      if(std::abs(rme-rhs)>1.0e-5){
        std::cout<<"ERROR (2nd little symmetry of neutron TBDRMEs): "<<rme<<"!="<<rhs<<std::endl;
        su3::finalize();
        return EXIT_FAILURE;
      }
    }
  }
  // Check symmetries of explicitly constructed proton-neutron TBDRMEs
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
		  it=pn_tbdrmes.begin(); it!=pn_tbdrmes.end(); it++){
    int rho0max=su3::mult(std::get<2>(it->first)[4],std::get<2>(it->first)[5],std::get<2>(it->first)[7],
		          std::get<2>(it->first)[8],std::get<2>(it->first)[10],std::get<2>(it->first)[11]);
    double factor=sqrt(double((std::get<1>(it->first).SpSnS().SS()+1)
		  *su3dim(std::get<1>(it->first).x_omega().lm(),std::get<1>(it->first).x_omega().mu()))
		  /double((std::get<0>(it->first).SpSnS().SS()+1)
		  *su3dim(std::get<0>(it->first).x_omega().lm(),std::get<0>(it->first).x_omega().mu())));
    int phase=std::get<2>(it->first)[10]+std::get<2>(it->first)[11]+std::get<0>(it->first).x_omega().lm()
	      +std::get<0>(it->first).x_omega().mu()-std::get<1>(it->first).x_omega().lm()-std::get<1>(it->first).x_omega().mu()
	      +(std::get<1>(it->first).SpSnS().SS()-std::get<0>(it->first).SpSnS().SS())/2+rho0max
	      +std::get<2>(it->first)[0]+std::get<2>(it->first)[1]+std::get<2>(it->first)[2]+std::get<2>(it->first)[3]
	      +std::get<2>(it->first)[4]+std::get<2>(it->first)[5]+std::get<2>(it->first)[7]+std::get<2>(it->first)[8];
    if((phase/2)*2!=phase)factor=-factor;
    int irho=-1;
    for(double rme : it->second){
      irho++;
      double rhs=0.0;
      for(int rho0p=1; rho0p<=rho0max; rho0p++){
	if((rho0p/2)*2==rho0p){
          rhs+=Phi(std::get<2>(it->first)[13]+1,rho0p,std::get<2>(it->first)[4],std::get<2>(it->first)[5],                                                       std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[10],std::get<2>(it->first)[11])
	       *pn_tbdrmes[{std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[3],std::get<2>(it->first)[2],
                 std::get<2>(it->first)[1],std::get<2>(it->first)[0],std::get<2>(it->first)[8],std::get<2>(it->first)[7],
	         std::get<2>(it->first)[9],std::get<2>(it->first)[5],std::get<2>(it->first)[4],std::get<2>(it->first)[6],
	         std::get<2>(it->first)[11],std::get<2>(it->first)[10],std::get<2>(it->first)[12],rho0p-1}}][irho];
        }else{
          rhs-=Phi(std::get<2>(it->first)[13]+1,rho0p,std::get<2>(it->first)[4],std::get<2>(it->first)[5],                                                       std::get<2>(it->first)[7],std::get<2>(it->first)[8],std::get<2>(it->first)[10],std::get<2>(it->first)[11])
               *pn_tbdrmes[{std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[3],std::get<2>(it->first)[2],
                 std::get<2>(it->first)[1],std::get<2>(it->first)[0],std::get<2>(it->first)[8],std::get<2>(it->first)[7],
                 std::get<2>(it->first)[9],std::get<2>(it->first)[5],std::get<2>(it->first)[4],std::get<2>(it->first)[6],
                 std::get<2>(it->first)[11],std::get<2>(it->first)[10],std::get<2>(it->first)[12],rho0p-1}}][irho];
	}
      }
      rhs=factor*rhs;
      if(std::abs(rme-rhs)>1.0e-3){
        std::cout<<"ERROR (symmetry of proton-neutron TBDRMEs): "<<rme<<"!="<<rhs<<std::endl;
	std::cout<<"N1 N2 N3 N4 lmf muf SSf lmi mui SSi lm0 mu0 SS0:"<<std::endl;
	std::cout<<std::get<2>(it->first)[0]<<" "<<std::get<2>(it->first)[1]<<" "<<std::get<2>(it->first)[2]<<" "
		 <<std::get<2>(it->first)[3]<<" "<<std::get<2>(it->first)[4]<<" "<<std::get<2>(it->first)[5]<<" "
		 <<std::get<2>(it->first)[6]<<" "<<std::get<2>(it->first)[7]<<" "<<std::get<2>(it->first)[8]<<" "
		 <<std::get<2>(it->first)[9]<<" "<<std::get<2>(it->first)[10]<<" "<<std::get<2>(it->first)[11]<<" "
		 <<std::get<2>(it->first)[12]<<std::endl;
        su3::finalize();
        return EXIT_FAILURE;
      }
    }
  }
  std::cout<<"Symmetries of explicitly constructed TBDRMEs checked"<<std::endl;

  // Read TBDRMEs from spncci
  std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>> proton_tbdrmes_spncci,neutron_tbdrmes_spncci,
                                                                                         pn_tbdrmes_spncci;
  std::ifstream file("spncci_output.dat");
  if(!file){
    std::cerr << "Could not open 'spncci_output.dat' file!" << std::endl;
    su3::finalize();
    exit(EXIT_FAILURE);
  }
  while(true){
    int gammap,Nex_sigmap,lm_sigmap,mu_sigmap,SSp_bra,SSn_bra,SS_bra,upsilonp,Nex_omegap,lm_omegap,mu_omegap,
        gamma,Nex_sigma,lm_sigma,mu_sigma,SSp_ket,SSn_ket,SS_ket,upsilon,Nex_omega,lm_omega,mu_omega,
        N1,N2,N3,N4,lmf,muf,Sf,lmi,mui,Si,lm0,mu0,SS0,rho0,rho,Tz;
    double rme;
    file>>gammap>>Nex_sigmap>>lm_sigmap>>mu_sigmap>>SSp_bra>>SSn_bra>>SS_bra>>upsilonp>>Nex_omegap>>lm_omegap>>mu_omegap
        >>gamma>>Nex_sigma>>lm_sigma>>mu_sigma>>SSp_ket>>SSn_ket>>SS_ket>>upsilon>>Nex_omega>>lm_omega>>mu_omega
        >>N1>>N2>>N3>>N4>>lmf>>muf>>Sf>>lmi>>mui>>Si>>lm0>>mu0>>SS0>>rho0>>rho>>Tz>>rme;
    if(file){
      Spins SpSnS_bra(SSp_bra,SSn_bra,SS_bra);
      Spins SpSnS_ket(SSp_ket,SSn_ket,SS_ket);
      SU3 sigmap(lm_sigmap,mu_sigmap);
      SU3 sigma(lm_sigma,mu_sigma);
      SU3 omegap(lm_omegap,mu_omegap);
      SU3 omega(lm_omega,mu_omega);
      ONBasisState bra(gammap,Nex_sigmap,sigmap,SpSnS_bra,upsilonp,Nex_omegap,omegap);
      ONBasisState ket(gamma,Nex_sigma,sigma,SpSnS_ket,upsilon,Nex_omega,omega);
      std::array<int,14> tensor={N1,N2,N3,N4,lmf,muf,2*Sf,lmi,mui,2*Si,lm0,mu0,SS0,rho0-1};
      std::tuple<ONBasisState,ONBasisState,std::array<int,14>> key={bra,ket,tensor};
      int rhomax=su3::mult(lm_omega,mu_omega,lm0,mu0,lm_omegap,mu_omegap);
      if(Tz==1){
        std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 it=proton_tbdrmes_spncci.find(key);
        if(it==proton_tbdrmes_spncci.end()){
          std::vector<double> rmes(rhomax,0.0);
          proton_tbdrmes_spncci[key]=rmes;
	}
	proton_tbdrmes_spncci[key][rho-1]=rme;
      }else if(Tz==-1){
	std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 it=neutron_tbdrmes_spncci.find(key);
        if(it==neutron_tbdrmes_spncci.end()){
          std::vector<double> rmes(rhomax,0.0);
          neutron_tbdrmes_spncci[key]=rmes;
        }
        neutron_tbdrmes_spncci[key][rho-1]=rme;
      }else{
        std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 it=pn_tbdrmes_spncci.find(key);
        if(it==pn_tbdrmes_spncci.end()){
          std::vector<double> rmes(rhomax,0.0);
          pn_tbdrmes_spncci[key]=rmes;
        }
        pn_tbdrmes_spncci[key][rho-1]=rme;
      }
    }else{
      break;
    }
  }
  std::cout<<"TBDRMEs from spncci read"<<std::endl;

  // Benchmark of explicitly constructed TBDRMEs and those from spncci
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 it=proton_tbdrmes_spncci.begin(); it!=proton_tbdrmes_spncci.end(); it++){
    int gammap=std::get<0>(it->first).gamma();
    int Nex_sigmap=std::get<0>(it->first).Nex_sigma();
    int lm_sigmap=std::get<0>(it->first).x_sigma().lm();
    int mu_sigmap=std::get<0>(it->first).x_sigma().mu();
    int SSp_bra=std::get<0>(it->first).SpSnS().SSp();
    int SSn_bra=std::get<0>(it->first).SpSnS().SSn();
    int SS_bra=std::get<0>(it->first).SpSnS().SS();
    int upsilonp=std::get<0>(it->first).upsilon();
    int Nex_omegap=std::get<0>(it->first).Nex_omega();
    int lm_omegap=std::get<0>(it->first).x_omega().lm();
    int mu_omegap=std::get<0>(it->first).x_omega().mu();
    int gamma=std::get<1>(it->first).gamma();
    int Nex_sigma=std::get<1>(it->first).Nex_sigma();
    int lm_sigma=std::get<1>(it->first).x_sigma().lm();
    int mu_sigma=std::get<1>(it->first).x_sigma().mu();
    int SSp_ket=std::get<1>(it->first).SpSnS().SSp();
    int SSn_ket=std::get<1>(it->first).SpSnS().SSn();
    int SS_ket=std::get<1>(it->first).SpSnS().SS();
    int upsilon=std::get<1>(it->first).upsilon();
    int Nex_omega=std::get<1>(it->first).Nex_omega();
    int lm_omega=std::get<1>(it->first).x_omega().lm();
    int mu_omega=std::get<1>(it->first).x_omega().mu();
    int N1=std::get<2>(it->first)[0];
    int N2=std::get<2>(it->first)[1];
    int N3=std::get<2>(it->first)[2];
    int N4=std::get<2>(it->first)[3];
    int lmf=std::get<2>(it->first)[4];
    int muf=std::get<2>(it->first)[5];
    int Sf=std::get<2>(it->first)[6]/2;
    int lmi=std::get<2>(it->first)[7];
    int mui=std::get<2>(it->first)[8];
    int Si=std::get<2>(it->first)[9]/2;
    int lm0=std::get<2>(it->first)[10];
    int mu0=std::get<2>(it->first)[11];
    int SS0=std::get<2>(it->first)[12];
    int rho0=std::get<2>(it->first)[13]+1;
    std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 iter=proton_tbdrmes.find(it->first);
    if(iter!=proton_tbdrmes.end()){
      int irho=-1;
      for(double rme : it->second){
        irho++;
	if(std::abs(std::abs(rme)-std::abs(proton_tbdrmes[it->first][irho]))>1.0e-5){
          std::cout<<"ERROR (proton): "<<rme<<" "<<proton_tbdrmes[it->first][irho]<<std::endl;
	  std::cout<<gammap<<" "<<Nex_sigmap<<" "<<lm_sigmap<<" "<<mu_sigmap<<" "<<SSp_bra<<" "<<SSn_bra<<" "<<SS_bra
		   <<" "<<upsilonp<<" "<<Nex_omegap<<" "<<lm_omegap<<" "<<mu_omegap
                   <<" "<<gamma<<" "<<Nex_sigma<<" "<<lm_sigma<<" "<<mu_sigma<<" "<<SSp_ket<<" "<<SSn_ket<<" "<<SS_ket
		   <<" "<<upsilon<<" "<<Nex_omega<<" "<<lm_omega<<" "<<mu_omega
                   <<" "<<N1<<" "<<N2<<" "<<N3<<" "<<N4<<" "<<lmf<<" "<<muf<<" "<<Sf<<" "<<lmi<<" "<<mui<<" "<<Si
		   <<" "<<lm0<<" "<<mu0<<" "<<SS0<<" "<<rho0<<" "<<irho+1<<" "<<1<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
	}
      }
    }else{
      std::cout<<"ERROR: didn't find proton RME among those constructed explicitly"<<std::endl;
      std::cout<<gammap<<" "<<Nex_sigmap<<" "<<lm_sigmap<<" "<<mu_sigmap<<" "<<SSp_bra<<" "<<SSn_bra<<" "<<SS_bra
                   <<" "<<upsilonp<<" "<<Nex_omegap<<" "<<lm_omegap<<" "<<mu_omegap
                   <<" "<<gamma<<" "<<Nex_sigma<<" "<<lm_sigma<<" "<<mu_sigma<<" "<<SSp_ket<<" "<<SSn_ket<<" "<<SS_ket
                   <<" "<<upsilon<<" "<<Nex_omega<<" "<<lm_omega<<" "<<mu_omega
                   <<" "<<N1<<" "<<N2<<" "<<N3<<" "<<N4<<" "<<lmf<<" "<<muf<<" "<<Sf<<" "<<lmi<<" "<<mui<<" "<<Si
                   <<" "<<lm0<<" "<<mu0<<" "<<SS0<<" "<<rho0<<" "<<"rho"<<" "<<1<<std::endl;
      su3::finalize();
      return EXIT_FAILURE;
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 it=neutron_tbdrmes_spncci.begin(); it!=neutron_tbdrmes_spncci.end(); it++){
int Nex_bra=std::get<0>(it->first).Nex_omega();
int Nex_ket=std::get<1>(it->first).Nex_omega();
    std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 iter=neutron_tbdrmes.find(it->first);
    if(iter!=neutron_tbdrmes.end()){
      int irho=-1;
      for(double rme : it->second){
        irho++;
        if(std::abs(std::abs(rme)-std::abs(neutron_tbdrmes[it->first][irho]))>1.0e-5){
          std::cout<<"ERROR (neutron): "<<rme<<" "<<neutron_tbdrmes[it->first][irho]<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
        }
      }
    }else{
      std::cout<<"ERROR: didn't find neutron RME among those constructed explicitly"<<std::endl;
      su3::finalize();
      return EXIT_FAILURE;
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 it=pn_tbdrmes_spncci.begin(); it!=pn_tbdrmes_spncci.end(); it++){
    int gammap=std::get<0>(it->first).gamma();
    int Nex_sigmap=std::get<0>(it->first).Nex_sigma();
    int lm_sigmap=std::get<0>(it->first).x_sigma().lm();
    int mu_sigmap=std::get<0>(it->first).x_sigma().mu();
    int SSp_bra=std::get<0>(it->first).SpSnS().SSp();
    int SSn_bra=std::get<0>(it->first).SpSnS().SSn();
    int SS_bra=std::get<0>(it->first).SpSnS().SS();
    int upsilonp=std::get<0>(it->first).upsilon();
    int Nex_omegap=std::get<0>(it->first).Nex_omega();
    int lm_omegap=std::get<0>(it->first).x_omega().lm();
    int mu_omegap=std::get<0>(it->first).x_omega().mu();
    int gamma=std::get<1>(it->first).gamma();
    int Nex_sigma=std::get<1>(it->first).Nex_sigma();
    int lm_sigma=std::get<1>(it->first).x_sigma().lm();
    int mu_sigma=std::get<1>(it->first).x_sigma().mu();
    int SSp_ket=std::get<1>(it->first).SpSnS().SSp();
    int SSn_ket=std::get<1>(it->first).SpSnS().SSn();
    int SS_ket=std::get<1>(it->first).SpSnS().SS();
    int upsilon=std::get<1>(it->first).upsilon();
    int Nex_omega=std::get<1>(it->first).Nex_omega();
    int lm_omega=std::get<1>(it->first).x_omega().lm();
    int mu_omega=std::get<1>(it->first).x_omega().mu();
    int N1=std::get<2>(it->first)[0];
    int N2=std::get<2>(it->first)[1];
    int N3=std::get<2>(it->first)[2];
    int N4=std::get<2>(it->first)[3];
    int lmf=std::get<2>(it->first)[4];
    int muf=std::get<2>(it->first)[5];
    int Sf=std::get<2>(it->first)[6]/2;
    int lmi=std::get<2>(it->first)[7];
    int mui=std::get<2>(it->first)[8];
    int Si=std::get<2>(it->first)[9]/2;
    int lm0=std::get<2>(it->first)[10];
    int mu0=std::get<2>(it->first)[11];
    int SS0=std::get<2>(it->first)[12];
    int rho0=std::get<2>(it->first)[13]+1;
    std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 iter=pn_tbdrmes.find(it->first);
    if(iter!=pn_tbdrmes.end()){
      int irho=-1;
      for(double rme : it->second){
        irho++;
	if(std::abs(std::abs(rme)-std::abs(pn_tbdrmes[it->first][irho]))>1.0e-4){
          std::cout<<"ERROR (proton-neutron): "<<rme<<" "<<pn_tbdrmes[it->first][irho]<<std::endl;
	  std::cout<<gammap<<" "<<Nex_sigmap<<" "<<lm_sigmap<<" "<<mu_sigmap<<" "<<SSp_bra<<" "<<SSn_bra<<" "<<SS_bra
		   <<" "<<upsilonp<<" "<<Nex_omegap<<" "<<lm_omegap<<" "<<mu_omegap
                   <<" "<<gamma<<" "<<Nex_sigma<<" "<<lm_sigma<<" "<<mu_sigma<<" "<<SSp_ket<<" "<<SSn_ket<<" "<<SS_ket
		   <<" "<<upsilon<<" "<<Nex_omega<<" "<<lm_omega<<" "<<mu_omega
                   <<" "<<N1<<" "<<N2<<" "<<N3<<" "<<N4<<" "<<lmf<<" "<<muf<<" "<<Sf<<" "<<lmi<<" "<<mui<<" "<<Si
		   <<" "<<lm0<<" "<<mu0<<" "<<SS0<<" "<<rho0<<" "<<irho+1<<" "<<1<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
	}
      }
    }else{
      std::cout<<"ERROR: didn't find proton-neutron RME among those constructed explicitly"<<std::endl;
      std::cout<<gammap<<" "<<Nex_sigmap<<" "<<lm_sigmap<<" "<<mu_sigmap<<" "<<SSp_bra<<" "<<SSn_bra<<" "<<SS_bra
                   <<" "<<upsilonp<<" "<<Nex_omegap<<" "<<lm_omegap<<" "<<mu_omegap
                   <<" "<<gamma<<" "<<Nex_sigma<<" "<<lm_sigma<<" "<<mu_sigma<<" "<<SSp_ket<<" "<<SSn_ket<<" "<<SS_ket
                   <<" "<<upsilon<<" "<<Nex_omega<<" "<<lm_omega<<" "<<mu_omega
                   <<" "<<N1<<" "<<N2<<" "<<N3<<" "<<N4<<" "<<lmf<<" "<<muf<<" "<<Sf<<" "<<lmi<<" "<<mui<<" "<<Si
                   <<" "<<lm0<<" "<<mu0<<" "<<SS0<<" "<<rho0<<" "<<"rho"<<" "<<1<<std::endl;
      su3::finalize();
      return EXIT_FAILURE;
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 it=proton_tbdrmes.begin(); it!=proton_tbdrmes.end(); it++){
    int gammap=std::get<0>(it->first).gamma();
    int Nex_sigmap=std::get<0>(it->first).Nex_sigma();
    int lm_sigmap=std::get<0>(it->first).x_sigma().lm();
    int mu_sigmap=std::get<0>(it->first).x_sigma().mu();
    int SSp_bra=std::get<0>(it->first).SpSnS().SSp();
    int SSn_bra=std::get<0>(it->first).SpSnS().SSn();
    int SS_bra=std::get<0>(it->first).SpSnS().SS();
    int upsilonp=std::get<0>(it->first).upsilon();
    int Nex_omegap=std::get<0>(it->first).Nex_omega();
    int lm_omegap=std::get<0>(it->first).x_omega().lm();
    int mu_omegap=std::get<0>(it->first).x_omega().mu();
    int gamma=std::get<1>(it->first).gamma();
    int Nex_sigma=std::get<1>(it->first).Nex_sigma();
    int lm_sigma=std::get<1>(it->first).x_sigma().lm();
    int mu_sigma=std::get<1>(it->first).x_sigma().mu();
    int SSp_ket=std::get<1>(it->first).SpSnS().SSp();
    int SSn_ket=std::get<1>(it->first).SpSnS().SSn();
    int SS_ket=std::get<1>(it->first).SpSnS().SS();
    int upsilon=std::get<1>(it->first).upsilon();
    int Nex_omega=std::get<1>(it->first).Nex_omega();
    int lm_omega=std::get<1>(it->first).x_omega().lm();
    int mu_omega=std::get<1>(it->first).x_omega().mu();
    int N1=std::get<2>(it->first)[0];
    int N2=std::get<2>(it->first)[1];
    int N3=std::get<2>(it->first)[2];
    int N4=std::get<2>(it->first)[3];
    int lmf=std::get<2>(it->first)[4];
    int muf=std::get<2>(it->first)[5];
    int Sf=std::get<2>(it->first)[6]/2;
    int lmi=std::get<2>(it->first)[7];
    int mui=std::get<2>(it->first)[8];
    int Si=std::get<2>(it->first)[9]/2;
    int lm0=std::get<2>(it->first)[10];
    int mu0=std::get<2>(it->first)[11];
    int SS0=std::get<2>(it->first)[12];
    int rho0=std::get<2>(it->first)[13]+1;
    if(std::get<0>(it->first).Nex_omega()<std::get<1>(it->first).Nex_omega())continue;
    int irho=-1;
    for(double rme : it->second){
      irho++;
      if(std::abs(rme)>1.0e-5){
	std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
		it2=proton_tbdrmes_spncci.find(it->first);
        if(it2==proton_tbdrmes_spncci.end()){
          std::cout<<"ERROR: didn't find proton RME among those from spncci"<<std::endl;
          std::cout<<gammap<<" "<<Nex_sigmap<<" "<<lm_sigmap<<" "<<mu_sigmap<<" "<<SSp_bra<<" "<<SSn_bra<<" "<<SS_bra
                   <<" "<<upsilonp<<" "<<Nex_omegap<<" "<<lm_omegap<<" "<<mu_omegap
                   <<" "<<gamma<<" "<<Nex_sigma<<" "<<lm_sigma<<" "<<mu_sigma<<" "<<SSp_ket<<" "<<SSn_ket<<" "<<SS_ket
                   <<" "<<upsilon<<" "<<Nex_omega<<" "<<lm_omega<<" "<<mu_omega
                   <<" "<<N1<<" "<<N2<<" "<<N3<<" "<<N4<<" "<<lmf<<" "<<muf<<" "<<Sf<<" "<<lmi<<" "<<mui<<" "<<Si
                   <<" "<<lm0<<" "<<mu0<<" "<<SS0<<" "<<rho0<<" "<<irho+1<<" "<<1<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
	}
	break;
      }
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 it=neutron_tbdrmes.begin(); it!=neutron_tbdrmes.end(); it++){
int Nex_bra=std::get<0>(it->first).Nex_omega();
int Nex_ket=std::get<1>(it->first).Nex_omega();
    if(std::get<0>(it->first).Nex_omega()<std::get<1>(it->first).Nex_omega())continue;
    for(double rme : it->second){
      if(std::abs(rme)>1.0e-5){
        std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                it2=neutron_tbdrmes_spncci.find(it->first);
        if(it2==neutron_tbdrmes_spncci.end()){
          std::cout<<"ERROR: didn't find neutron RME among those from spncci"<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
        }
        break;
      }
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
                 it=pn_tbdrmes.begin(); it!=pn_tbdrmes.end(); it++){
    int gammap=std::get<0>(it->first).gamma();
    int Nex_sigmap=std::get<0>(it->first).Nex_sigma();
    int lm_sigmap=std::get<0>(it->first).x_sigma().lm();
    int mu_sigmap=std::get<0>(it->first).x_sigma().mu();
    int SSp_bra=std::get<0>(it->first).SpSnS().SSp();
    int SSn_bra=std::get<0>(it->first).SpSnS().SSn();
    int SS_bra=std::get<0>(it->first).SpSnS().SS();
    int upsilonp=std::get<0>(it->first).upsilon();
    int Nex_omegap=std::get<0>(it->first).Nex_omega();
    int lm_omegap=std::get<0>(it->first).x_omega().lm();
    int mu_omegap=std::get<0>(it->first).x_omega().mu();
    int gamma=std::get<1>(it->first).gamma();
    int Nex_sigma=std::get<1>(it->first).Nex_sigma();
    int lm_sigma=std::get<1>(it->first).x_sigma().lm();
    int mu_sigma=std::get<1>(it->first).x_sigma().mu();
    int SSp_ket=std::get<1>(it->first).SpSnS().SSp();
    int SSn_ket=std::get<1>(it->first).SpSnS().SSn();
    int SS_ket=std::get<1>(it->first).SpSnS().SS();
    int upsilon=std::get<1>(it->first).upsilon();
    int Nex_omega=std::get<1>(it->first).Nex_omega();
    int lm_omega=std::get<1>(it->first).x_omega().lm();
    int mu_omega=std::get<1>(it->first).x_omega().mu();
    int N1=std::get<2>(it->first)[0];
    int N2=std::get<2>(it->first)[1];
    int N3=std::get<2>(it->first)[2];
    int N4=std::get<2>(it->first)[3];
    int lmf=std::get<2>(it->first)[4];
    int muf=std::get<2>(it->first)[5];
    int Sf=std::get<2>(it->first)[6]/2;
    int lmi=std::get<2>(it->first)[7];
    int mui=std::get<2>(it->first)[8];
    int Si=std::get<2>(it->first)[9]/2;
    int lm0=std::get<2>(it->first)[10];
    int mu0=std::get<2>(it->first)[11];
    int SS0=std::get<2>(it->first)[12];
    int rho0=std::get<2>(it->first)[13]+1;
    if(std::get<0>(it->first).Nex_omega()<std::get<1>(it->first).Nex_omega())continue;
    int irho=-1;
    for(double rme : it->second){
      irho++;
      if(std::abs(rme)>1.0e-5){
	std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,14>>,std::vector<double>>::iterator
		it2=pn_tbdrmes_spncci.find(it->first);
        if(it2==pn_tbdrmes_spncci.end()){
          std::cout<<"ERROR: didn't find proton-neutron RME among those from spncci"<<std::endl;
          std::cout<<gammap<<" "<<Nex_sigmap<<" "<<lm_sigmap<<" "<<mu_sigmap<<" "<<SSp_bra<<" "<<SSn_bra<<" "<<SS_bra
                   <<" "<<upsilonp<<" "<<Nex_omegap<<" "<<lm_omegap<<" "<<mu_omegap
                   <<" "<<gamma<<" "<<Nex_sigma<<" "<<lm_sigma<<" "<<mu_sigma<<" "<<SSp_ket<<" "<<SSn_ket<<" "<<SS_ket
                   <<" "<<upsilon<<" "<<Nex_omega<<" "<<lm_omega<<" "<<mu_omega
                   <<" "<<N1<<" "<<N2<<" "<<N3<<" "<<N4<<" "<<lmf<<" "<<muf<<" "<<Sf<<" "<<lmi<<" "<<mui<<" "<<Si
                   <<" "<<lm0<<" "<<mu0<<" "<<SS0<<" "<<rho0<<" "<<irho+1<<" "<<1<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
	}
	break;
      }
    }
  }
  std::cout<<"Benchmark of explicitly constructed TBDRMEs and those from spncci finished"<<std::endl;

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
/*
  // Prepare proton seeds
  std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>> proton_obdrmes_recursive;
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                  it=proton_obdrmes.begin(); it!=proton_obdrmes.end(); it++){
    if(std::get<0>(it->first).Nex_omega()==std::get<0>(it->first).Nex_sigma()
		    && std::get<1>(it->first).Nex_omega()==std::get<1>(it->first).Nex_sigma()){
      proton_obdrmes_recursive[it->first]=it->second;
    }
  }
  // Prepare neutron seeds
  std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>> neutron_obdrmes_recursive;
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                  it=neutron_obdrmes.begin(); it!=neutron_obdrmes.end(); it++){
    if(std::get<0>(it->first).Nex_omega()==std::get<0>(it->first).Nex_sigma()
                    && std::get<1>(it->first).Nex_omega()==std::get<1>(it->first).Nex_sigma()){
      neutron_obdrmes_recursive[it->first]=it->second;
    }
  }
  std::cout<<"Seeds prepared"<<std::endl;

  // Calculate OBDRMEs with Nn'=0 by recurrence
  for(int Nn=2; Nn<=Nmax; Nn=Nn+2){
    for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>::iterator
		   it1bra=orthonormal_basis.begin(); it1bra!=orthonormal_basis.end(); it1bra++){
      for(std::map<SU3,std::map<int,Eigen::VectorXd>>:: iterator
		      it2bra=it1bra->second[0].begin(); it2bra!=it1bra->second[0].end(); it2bra++){
	for(std::map<int,Eigen::VectorXd>::iterator it3bra=it2bra->second.begin(); it3bra!=it2bra->second.end(); it3bra++){
	  ONBasisState bra(std::get<0>(it1bra->first),std::get<1>(it1bra->first),std::get<2>(it1bra->first),std::get<3>(it1bra->first),
                      it3bra->first,0,it2bra->first);
          for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>::iterator
                   it1ket=orthonormal_basis.begin(); it1ket!=orthonormal_basis.end(); it1ket++){
            for(std::map<SU3,std::map<int,Eigen::VectorXd>>:: iterator
                      it2ket=it1ket->second[Nn].begin(); it2ket!=it1ket->second[Nn].end(); it2ket++){
              for(std::map<int,Eigen::VectorXd>::iterator it3ket=it2ket->second.begin(); it3ket!=it2ket->second.end(); it3ket++){
                ONBasisState ket(std::get<0>(it1ket->first),std::get<1>(it1ket->first),std::get<2>(it1ket->first),std::get<3>(it1ket->first),
                      it3ket->first,Nn,it2ket->first);
                for(const std::array<int,5> tensor : tensors[0]){
		  if(ket.Nex_omega()+tensor[0]-tensor[1]!=bra.Nex_omega())continue;
                  if(std::abs(ket.SpSnS().SS()-tensor[4])>bra.SpSnS().SS() || bra.SpSnS().SS()>ket.SpSnS().SS()+tensor[4])continue;
                  int rhomax=su3::mult(ket.x_omega().lm(),ket.x_omega().mu(),tensor[2],tensor[3],bra.x_omega().lm(),bra.x_omega().mu());
                  if(rhomax==0)continue;
                  std::vector<double> proton_rmes(rhomax);
                  std::vector<double> neutron_rmes(rhomax);
                  proton_rmes=Recurrence(bra,ket,tensor,proton_obdrmes_recursive,orthonormal_basis,A,K,Kinv,AA,N1vp);
                  neutron_rmes=Recurrence(bra,ket,tensor,neutron_obdrmes_recursive,orthonormal_basis,A,K,Kinv,AA,N1vn);
                  proton_obdrmes_recursive[{bra,ket,tensor}]=proton_rmes;
                  neutron_obdrmes_recursive[{bra,ket,tensor}]=neutron_rmes;
                }
	      }
            }
          }
        }
      }
    }
  }
  std::cout<<"OBDRMEs with Nn'=0 calculated"<<std::endl;

  // Calculate proton OBDRMEs with Nn=0 by conjugation
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                  it=proton_obdrmes_recursive.begin(); it!=proton_obdrmes_recursive.end(); it++){
    if(std::get<1>(it->first).Nex_omega()!=std::get<1>(it->first).Nex_sigma()){
      double factor=sqrt(double((std::get<0>(it->first).SpSnS().SS()+1)*su3dim(std::get<0>(it->first).x_omega()))
                        /double((std::get<1>(it->first).SpSnS().SS()+1)*su3dim(std::get<1>(it->first).x_omega())));
      int phase=std::get<2>(it->first)[1]-std::get<2>(it->first)[0]+std::get<1>(it->first).x_omega().lm()
	        +std::get<1>(it->first).x_omega().mu()-std::get<0>(it->first).x_omega().lm()-std::get<0>(it->first).x_omega().mu()
		+(std::get<1>(it->first).SpSnS().SS()-std::get<0>(it->first).SpSnS().SS())/2;
      if((phase/2)*2!=phase)factor=-factor;
      std::tuple<ONBasisState,ONBasisState,std::array<int,5>>
           key={std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[1],std::get<2>(it->first)[0],
                std::get<2>(it->first)[3],std::get<2>(it->first)[2],std::get<2>(it->first)[4]}};
      std::vector<double> rmes(it->second.size());
      for(int irho=0; irho<it->second.size(); irho++){
        rmes[irho]=factor*it->second[irho];
      }
      proton_obdrmes_recursive[key]=rmes;
    }
  }
  // Calculate neutron OBDRMEs with Nn=0 by conjugation
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                  it=neutron_obdrmes_recursive.begin(); it!=neutron_obdrmes_recursive.end(); it++){
    if(std::get<1>(it->first).Nex_omega()!=std::get<1>(it->first).Nex_sigma()){
      double factor=sqrt(double((std::get<0>(it->first).SpSnS().SS()+1)*su3dim(std::get<0>(it->first).x_omega()))
                        /double((std::get<1>(it->first).SpSnS().SS()+1)*su3dim(std::get<1>(it->first).x_omega())));
      int phase=std::get<2>(it->first)[1]-std::get<2>(it->first)[0]+std::get<1>(it->first).x_omega().lm()
                +std::get<1>(it->first).x_omega().mu()-std::get<0>(it->first).x_omega().lm()-std::get<0>(it->first).x_omega().mu()
                +(std::get<1>(it->first).SpSnS().SS()-std::get<0>(it->first).SpSnS().SS())/2;
      if((phase/2)*2!=phase)factor=-factor;
      std::tuple<ONBasisState,ONBasisState,std::array<int,5>>
           key={std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[1],std::get<2>(it->first)[0],
                std::get<2>(it->first)[3],std::get<2>(it->first)[2],std::get<2>(it->first)[4]}};
      std::vector<double> rmes(it->second.size());
      for(int irho=0; irho<it->second.size(); irho++){
        rmes[irho]=factor*it->second[irho];
      }
      neutron_obdrmes_recursive[key]=rmes;
    }
  }
  std::cout<<"OBDRMEs with Nn=0 calculated"<<std::endl;

  // Calculate OBDRMEs with Nn'>=Nn by recurrence
  for(int Nnp=2; Nnp<=Nmax; Nnp=Nnp+2){
    for(int Nn=2; Nn<=Nnp; Nn=Nn+2){
      for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>::iterator
                   it1bra=orthonormal_basis.begin(); it1bra!=orthonormal_basis.end(); it1bra++){
        for(std::map<SU3,std::map<int,Eigen::VectorXd>>:: iterator
                      it2bra=it1bra->second[Nnp].begin(); it2bra!=it1bra->second[Nnp].end(); it2bra++){
          for(std::map<int,Eigen::VectorXd>::iterator it3bra=it2bra->second.begin(); it3bra!=it2bra->second.end(); it3bra++){
            ONBasisState bra(std::get<0>(it1bra->first),std::get<1>(it1bra->first),std::get<2>(it1bra->first),std::get<3>(it1bra->first),
                      it3bra->first,Nnp,it2bra->first);
            for(std::map<std::tuple<int,int,SU3,Spins>,std::map<int,std::map<SU3,std::map<int,Eigen::VectorXd>>>>::iterator
                   it1ket=orthonormal_basis.begin(); it1ket!=orthonormal_basis.end(); it1ket++){
              for(std::map<SU3,std::map<int,Eigen::VectorXd>>:: iterator
                      it2ket=it1ket->second[Nn].begin(); it2ket!=it1ket->second[Nn].end(); it2ket++){
                for(std::map<int,Eigen::VectorXd>::iterator it3ket=it2ket->second.begin(); it3ket!=it2ket->second.end(); it3ket++){
                  ONBasisState ket(std::get<0>(it1ket->first),std::get<1>(it1ket->first),std::get<2>(it1ket->first),
				  std::get<3>(it1ket->first),it3ket->first,Nn,it2ket->first);
                  for(const std::array<int,5> tensor : tensors[0]){
  	            if(ket.Nex_omega()+tensor[0]-tensor[1]!=bra.Nex_omega())continue;
                    if(std::abs(ket.SpSnS().SS()-tensor[4])>bra.SpSnS().SS() || bra.SpSnS().SS()>ket.SpSnS().SS()+tensor[4])continue;
                    int rhomax=su3::mult(ket.x_omega().lm(),ket.x_omega().mu(),tensor[2],tensor[3],bra.x_omega().lm(),bra.x_omega().mu());
                    if(rhomax==0)continue;
                    std::vector<double> proton_rmes(rhomax);
                    std::vector<double> neutron_rmes(rhomax);
                    proton_rmes=Recurrence(bra,ket,tensor,proton_obdrmes_recursive,orthonormal_basis,A,K,Kinv,AA,N1vp);
                    neutron_rmes=Recurrence(bra,ket,tensor,neutron_obdrmes_recursive,orthonormal_basis,A,K,Kinv,AA,N1vn);
                    proton_obdrmes_recursive[{bra,ket,tensor}]=proton_rmes;
                    neutron_obdrmes_recursive[{bra,ket,tensor}]=neutron_rmes;
                  }
	        }
	      }
	    }
 	  }
        }
      }
    }
  }
  std::cout<<"OBDRMEs with Nn'>=Nn calculated"<<std::endl;

  // Calculate proton OBDRMEs with 0<Nn'<Nn by conjugation
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                  it=proton_obdrmes_recursive.begin(); it!=proton_obdrmes_recursive.end(); it++){
    if(std::get<1>(it->first).Nex_omega()!=std::get<1>(it->first).Nex_sigma()
		    && std::get<1>(it->first).Nex_omega()-std::get<1>(it->first).Nex_sigma()
		       <std::get<0>(it->first).Nex_omega()-std::get<0>(it->first).Nex_sigma()){
      double factor=sqrt(double((std::get<0>(it->first).SpSnS().SS()+1)*su3dim(std::get<0>(it->first).x_omega()))
                        /double((std::get<1>(it->first).SpSnS().SS()+1)*su3dim(std::get<1>(it->first).x_omega())));
      int phase=std::get<2>(it->first)[1]-std::get<2>(it->first)[0]+std::get<1>(it->first).x_omega().lm()
                +std::get<1>(it->first).x_omega().mu()-std::get<0>(it->first).x_omega().lm()-std::get<0>(it->first).x_omega().mu()
                +(std::get<1>(it->first).SpSnS().SS()-std::get<0>(it->first).SpSnS().SS())/2;
      if((phase/2)*2!=phase)factor=-factor;
      std::tuple<ONBasisState,ONBasisState,std::array<int,5>>
           key={std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[1],std::get<2>(it->first)[0],
                std::get<2>(it->first)[3],std::get<2>(it->first)[2],std::get<2>(it->first)[4]}};
      std::vector<double> rmes(it->second.size());
      for(int irho=0; irho<it->second.size(); irho++){
        rmes[irho]=factor*it->second[irho];
      }
      proton_obdrmes_recursive[key]=rmes;
    }
  }
  // Calculate neutron OBDRMEs with 0<Nn'<Nn by conjugation
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                  it=neutron_obdrmes_recursive.begin(); it!=neutron_obdrmes_recursive.end(); it++){
    if(std::get<1>(it->first).Nex_omega()!=std::get<1>(it->first).Nex_sigma()
                    && std::get<1>(it->first).Nex_omega()-std::get<1>(it->first).Nex_sigma()
                       <std::get<0>(it->first).Nex_omega()-std::get<0>(it->first).Nex_sigma()){
      double factor=sqrt(double((std::get<0>(it->first).SpSnS().SS()+1)*su3dim(std::get<0>(it->first).x_omega()))
                        /double((std::get<1>(it->first).SpSnS().SS()+1)*su3dim(std::get<1>(it->first).x_omega())));
      int phase=std::get<2>(it->first)[1]-std::get<2>(it->first)[0]+std::get<1>(it->first).x_omega().lm()
                +std::get<1>(it->first).x_omega().mu()-std::get<0>(it->first).x_omega().lm()-std::get<0>(it->first).x_omega().mu()
                +(std::get<1>(it->first).SpSnS().SS()-std::get<0>(it->first).SpSnS().SS())/2;
      if((phase/2)*2!=phase)factor=-factor;
      std::tuple<ONBasisState,ONBasisState,std::array<int,5>>
           key={std::get<1>(it->first),std::get<0>(it->first),{std::get<2>(it->first)[1],std::get<2>(it->first)[0],
                std::get<2>(it->first)[3],std::get<2>(it->first)[2],std::get<2>(it->first)[4]}};
      std::vector<double> rmes(it->second.size());
      for(int irho=0; irho<it->second.size(); irho++){
        rmes[irho]=factor*it->second[irho];
      }
      neutron_obdrmes_recursive[key]=rmes;
    }
  }
  std::cout<<"OBDRMEs with 0<Nn'<Nn calculated"<<std::endl;

  // Benchmark recursively calculated OBDRMEs with those calculated explicitly
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
		 it=proton_obdrmes_recursive.begin(); it!=proton_obdrmes_recursive.end(); it++){
    int irho=-1;
    for(double rme : it->second){
      irho++;
      if(std::abs(rme-proton_obdrmes[it->first][irho])>1.0e-4 ){
		      //(std::get<0>(it->first).Nex_omega()==2 && std::get<1>(it->first).Nex_omega()==4)){
        std::cout<<"ERROR: proton recursive "<<rme<<" != explicit "<<proton_obdrmes[it->first][irho]<<std::endl;
	if(std::get<0>(it->first).Nex_omega()<6 && std::get<1>(it->first).Nex_omega()<6)std::cout<<"Nmax4error"<<std::endl;
	if(std::get<0>(it->first).Nex_omega()==6 && std::get<1>(it->first).Nex_omega()==0)std::cout<<"N60error"<<std::endl;
	if(std::get<0>(it->first).Nex_omega()==6 && std::get<1>(it->first).Nex_omega()==2)std::cout<<"N62error"<<std::endl;
	if(std::get<0>(it->first).Nex_omega()==6 && std::get<1>(it->first).Nex_omega()==4)std::cout<<"N64error"<<std::endl;
	if(std::get<0>(it->first).Nex_omega()==6 && std::get<1>(it->first).Nex_omega()==6)std::cout<<"N66error"<<std::endl;
	if(std::get<0>(it->first).Nex_omega()==0 && std::get<1>(it->first).Nex_omega()==6)std::cout<<"N06error"<<std::endl;
	if(std::get<0>(it->first).Nex_omega()==2 && std::get<1>(it->first).Nex_omega()==6)std::cout<<"N26error"<<std::endl;
	if(std::get<0>(it->first).Nex_omega()==4 && std::get<1>(it->first).Nex_omega()==6)std::cout<<"N46error"<<std::endl;
	std::cout<<"bra: "<<std::get<0>(it->first).gamma()<<" "<<
          std::get<0>(it->first).Nex_sigma()<<" "<<std::get<0>(it->first).x_sigma().lm()<<" "<<std::get<0>(it->first).x_sigma().mu()<<" "<<
          std::get<0>(it->first).SpSnS().SSp()<<" "<<std::get<0>(it->first).SpSnS().SSn()<<" "<<std::get<0>(it->first).SpSnS().SS()<<" "<<
          std::get<0>(it->first).upsilon()<<" "<<std::get<0>(it->first).Nex_omega()<<" "<<std::get<0>(it->first).x_omega().lm()<<" "<<
          std::get<0>(it->first).x_omega().mu()<<std::endl;
	std::cout<<"ket: "<<std::get<1>(it->first).gamma()<<" "<<
          std::get<1>(it->first).Nex_sigma()<<" "<<std::get<1>(it->first).x_sigma().lm()<<" "<<std::get<1>(it->first).x_sigma().mu()<<" "<<
          std::get<1>(it->first).SpSnS().SSp()<<" "<<std::get<1>(it->first).SpSnS().SSn()<<" "<<std::get<1>(it->first).SpSnS().SS()<<" "<<
          std::get<1>(it->first).upsilon()<<" "<<std::get<1>(it->first).Nex_omega()<<" "<<std::get<1>(it->first).x_omega().lm()<<" "<<
          std::get<1>(it->first).x_omega().mu()<<std::endl;
	std::cout<<"tensor: "<<std::get<2>(it->first)[0]<<" "<<std::get<2>(it->first)[1]<<" "<<std::get<2>(it->first)[2]<<" "<<
	  std::get<2>(it->first)[3]<<" "<<std::get<2>(it->first)[4]<<std::endl;
        std::cout<<"rho: "<<irho+1<<std::endl;
//        su3::finalize();
//        return EXIT_FAILURE;
      }
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                 it=proton_obdrmes.begin(); it!=proton_obdrmes.end(); it++){
    std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
	    it2=proton_obdrmes_recursive.find(it->first);
    if(it2==proton_obdrmes_recursive.end()){
      std::cout<<"ERROR: didn't find proton explicit RME among recursive"<<std::endl;
      su3::finalize();
      return EXIT_FAILURE;
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                 it=neutron_obdrmes_recursive.begin(); it!=neutron_obdrmes_recursive.end(); it++){
    int irho=-1;
    for(double rme : it->second){
      irho++;
      if(std::abs(rme-neutron_obdrmes[it->first][irho])>1.0e-4){
        std::cout<<"ERROR: neutron recursive "<<rme<<" != explicit "<<neutron_obdrmes[it->first][irho]<<std::endl;
	if(std::get<0>(it->first).Nex_omega()<6 && std::get<1>(it->first).Nex_omega()<6)std::cout<<"Nmax4error"<<std::endl;
	std::cout<<"bra: "<<std::get<0>(it->first).gamma()<<" "<<
          std::get<0>(it->first).Nex_sigma()<<" "<<std::get<0>(it->first).x_sigma().lm()<<" "<<std::get<0>(it->first).x_sigma().mu()<<" "<<
          std::get<0>(it->first).SpSnS().SSp()<<" "<<std::get<0>(it->first).SpSnS().SSn()<<" "<<std::get<0>(it->first).SpSnS().SS()<<" "<<
          std::get<0>(it->first).upsilon()<<" "<<std::get<0>(it->first).Nex_omega()<<" "<<std::get<0>(it->first).x_omega().lm()<<" "<<
          std::get<0>(it->first).x_omega().mu()<<std::endl;
        std::cout<<"ket: "<<std::get<1>(it->first).gamma()<<" "<<
          std::get<1>(it->first).Nex_sigma()<<" "<<std::get<1>(it->first).x_sigma().lm()<<" "<<std::get<1>(it->first).x_sigma().mu()<<" "<<
          std::get<1>(it->first).SpSnS().SSp()<<" "<<std::get<1>(it->first).SpSnS().SSn()<<" "<<std::get<1>(it->first).SpSnS().SS()<<" "<<
          std::get<1>(it->first).upsilon()<<" "<<std::get<1>(it->first).Nex_omega()<<" "<<std::get<1>(it->first).x_omega().lm()<<" "<<
          std::get<1>(it->first).x_omega().mu()<<std::endl;
        std::cout<<"tensor: "<<std::get<2>(it->first)[0]<<" "<<std::get<2>(it->first)[1]<<" "<<std::get<2>(it->first)[2]<<" "<<
          std::get<2>(it->first)[3]<<" "<<std::get<2>(it->first)[4]<<std::endl;
        std::cout<<"rho: "<<irho+1<<std::endl;
        su3::finalize();
        return EXIT_FAILURE;
      }
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                 it=neutron_obdrmes.begin(); it!=neutron_obdrmes.end(); it++){
    std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
            it2=neutron_obdrmes_recursive.find(it->first);
    if(it2==neutron_obdrmes_recursive.end()){
      std::cout<<"ERROR: didn't find neutron explicit RME among recursive"<<std::endl;
      su3::finalize();
      return EXIT_FAILURE;
    }
  }
  std::cout<<"Benchmark of recursively and explicitly calculated OBDRMEs finished"<<std::endl;

  // Benchmark with results from spncci
  std::map<std::array<int,28>,double> proton_obdrmes_spncci,neutron_obdrmes_spncci;
  std::ifstream file("spncci_output.dat");
  if(!file){
    std::cerr << "Could not open 'spncci_output.dat' file!" << std::endl;
    su3::finalize();
    exit(EXIT_FAILURE);
  }
  while(true){
    int gammap,Nex_sigmap,lm_sigmap,mu_sigmap,SSp_bra,SSn_bra,SS_bra,upsilonp,Nex_omegap,lm_omegap,mu_omegap,
        gamma,Nex_sigma,lm_sigma,mu_sigma,SSp_ket,SSn_ket,SS_ket,upsilon,Nex_omega,lm_omega,mu_omega,
        operator_index,Np,N,lm0,mu0,SS0,rho0;
    double rme;
    file>>gammap>>Nex_sigmap>>lm_sigmap>>mu_sigmap>>SSp_bra>>SSn_bra>>SS_bra>>upsilonp>>Nex_omegap>>lm_omegap>>mu_omegap
        >>gamma>>Nex_sigma>>lm_sigma>>mu_sigma>>SSp_ket>>SSn_ket>>SS_ket>>upsilon>>Nex_omega>>lm_omega>>mu_omega
        >>operator_index>>Np>>N>>lm0>>mu0>>SS0>>rho0>>rme;
    if(file){
      Spins SpSnS_bra(SSp_bra,SSn_bra,SS_bra);
      Spins SpSnS_ket(SSp_ket,SSn_ket,SS_ket);
      SU3 sigmap(lm_sigmap,mu_sigmap);
      SU3 sigma(lm_sigma,mu_sigma);
      SU3 omegap(lm_omegap,mu_omegap);
      SU3 omega(lm_omega,mu_omega);
      ONBasisState bra(gammap,Nex_sigmap,sigmap,SpSnS_bra,upsilonp,Nex_omegap,omegap);
      ONBasisState ket(gamma,Nex_sigma,sigma,SpSnS_ket,upsilon,Nex_omega,omega);
      std::array<int,5> tensor={Np,N,lm0,mu0,SS0};
      std::tuple<ONBasisState,ONBasisState,std::array<int,5>> key={bra,ket,tensor};
      if(operator_index==0){
        std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                 it=proton_obdrmes_recursive.find(key);
        if(it!=proton_obdrmes_recursive.end()){
          double myrme=proton_obdrmes_recursive[key][rho0-1];
          if(std::abs(std::abs(rme)-std::abs(myrme))>1.0e-5){
            std::cout<<"ERROR: "<<rme<<" "<<myrme<<std::endl;
	    std::cout<<"bra: "<<gammap<<" "<<Nex_sigmap<<" "<<lm_sigmap<<" "<<mu_sigmap<<" "<<SSp_bra<<" "<<SSn_bra<<" "<<SS_bra
		     <<" "<<upsilonp<<" "<<Nex_omegap<<" "<<lm_omegap<<" "<<mu_omegap<<std::endl;
	    std::cout<<"ket: "<<gamma<<" "<<Nex_sigma<<" "<<lm_sigma<<" "<<mu_sigma<<" "<<SSp_ket<<" "<<SSn_ket<<" "<<SS_ket
                     <<" "<<upsilon<<" "<<Nex_omega<<" "<<lm_omega<<" "<<mu_omega<<std::endl;
	    std::cout<<"tensor: "<<Np<<" "<<N<<" "<<lm0<<" "<<mu0<<" "<<SS0<<std::endl;
	    std::cout<<"rho: "<<rho0<<std::endl;
            su3::finalize();
            return EXIT_FAILURE;
          }
	}else if(std::abs(rme)>1.0e-5){
          std::cout<<"ERROR: didn't find RME "<<rme<<" among my RMEs"<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
	}
	proton_obdrmes_spncci[{gammap,Nex_sigmap,lm_sigmap,mu_sigmap,SSp_bra,SSn_bra,SS_bra,upsilonp,Nex_omegap,lm_omegap,mu_omegap,
	                       gamma,Nex_sigma,lm_sigma,mu_sigma,SSp_ket,SSn_ket,SS_ket,upsilon,Nex_omega,lm_omega,mu_omega,
	                       Np,N,lm0,mu0,SS0,rho0}]=rme;
      }else{
	std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                 it=neutron_obdrmes_recursive.find(key);
        if(it!=neutron_obdrmes_recursive.end()){
          double myrme=neutron_obdrmes_recursive[key][rho0-1];
          if(std::abs(std::abs(rme)-std::abs(myrme))>1.0e-5){
            std::cout<<"ERROR: "<<rme<<" "<<myrme<<std::endl;
	    std::cout<<"bra: "<<gammap<<" "<<Nex_sigmap<<" "<<lm_sigmap<<" "<<mu_sigmap<<" "<<SSp_bra<<" "<<SSn_bra<<" "<<SS_bra
                     <<" "<<upsilonp<<" "<<Nex_omegap<<" "<<lm_omegap<<" "<<mu_omegap<<std::endl;
            std::cout<<"ket: "<<gamma<<" "<<Nex_sigma<<" "<<lm_sigma<<" "<<mu_sigma<<" "<<SSp_ket<<" "<<SSn_ket<<" "<<SS_ket
                     <<" "<<upsilon<<" "<<Nex_omega<<" "<<lm_omega<<" "<<mu_omega<<std::endl;
            std::cout<<"tensor: "<<Np<<" "<<N<<" "<<lm0<<" "<<mu0<<" "<<SS0<<std::endl;
            std::cout<<"rho: "<<rho0<<std::endl;
            su3::finalize();
            return EXIT_FAILURE;
          }
        }else if(std::abs(rme)>1.0e-5){
          std::cout<<"ERROR: didn't find RME "<<rme<<" among my RMEs"<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
        }
	neutron_obdrmes_spncci[{gammap,Nex_sigmap,lm_sigmap,mu_sigmap,SSp_bra,SSn_bra,SS_bra,upsilonp,Nex_omegap,lm_omegap,mu_omegap,
                               gamma,Nex_sigma,lm_sigma,mu_sigma,SSp_ket,SSn_ket,SS_ket,upsilon,Nex_omega,lm_omega,mu_omega,
                               Np,N,lm0,mu0,SS0,rho0}]=rme;
      }
    }else{
      break;
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                 it=proton_obdrmes_recursive.begin(); it!=proton_obdrmes_recursive.end(); it++){
    if(std::get<0>(it->first).Nex_omega()<std::get<1>(it->first).Nex_omega())continue;
    int rho=0;
    for(double rme : it->second){
      rho++;
      if(std::abs(rme)>1.0e-5){
        std::map<std::array<int,28>,double>::iterator it2=proton_obdrmes_spncci.find({std::get<0>(it->first).gamma(),
	  std::get<0>(it->first).Nex_sigma(),std::get<0>(it->first).x_sigma().lm(),std::get<0>(it->first).x_sigma().mu(),
	  std::get<0>(it->first).SpSnS().SSp(),std::get<0>(it->first).SpSnS().SSn(),std::get<0>(it->first).SpSnS().SS(),
	  std::get<0>(it->first).upsilon(),std::get<0>(it->first).Nex_omega(),std::get<0>(it->first).x_omega().lm(),
	  std::get<0>(it->first).x_omega().mu(),
	  std::get<1>(it->first).gamma(),
          std::get<1>(it->first).Nex_sigma(),std::get<1>(it->first).x_sigma().lm(),std::get<1>(it->first).x_sigma().mu(),
          std::get<1>(it->first).SpSnS().SSp(),std::get<1>(it->first).SpSnS().SSn(),std::get<1>(it->first).SpSnS().SS(),
          std::get<1>(it->first).upsilon(),std::get<1>(it->first).Nex_omega(),std::get<1>(it->first).x_omega().lm(),
          std::get<1>(it->first).x_omega().mu(),
	  std::get<2>(it->first)[0],std::get<2>(it->first)[1],std::get<2>(it->first)[2],std::get<2>(it->first)[3],std::get<2>(it->first)[4],
	  rho});
        if(it2==proton_obdrmes_spncci.end()){
          std::cout<<"ERROR: didn't find proton RME "<<rme<<" among those from spncci"<<std::endl;
	  std::cout<<"bra: "<<std::get<0>(it->first).gamma()<<" "
          <<std::get<0>(it->first).Nex_sigma()<<" "<<std::get<0>(it->first).x_sigma().lm()<<" "<<std::get<0>(it->first).x_sigma().mu()<<" "
          <<std::get<0>(it->first).SpSnS().SSp()<<" "<<std::get<0>(it->first).SpSnS().SSn()<<" "<<std::get<0>(it->first).SpSnS().SS()<<" "
          <<std::get<0>(it->first).upsilon()<<" "<<std::get<0>(it->first).Nex_omega()<<" "<<std::get<0>(it->first).x_omega().lm()<<" "
          <<std::get<0>(it->first).x_omega().mu()<<std::endl;
	  std::cout<<"ket: "<<std::get<1>(it->first).gamma()<<" "
          <<std::get<1>(it->first).Nex_sigma()<<" "<<std::get<1>(it->first).x_sigma().lm()<<" "<<std::get<1>(it->first).x_sigma().mu()<<" "
          <<std::get<1>(it->first).SpSnS().SSp()<<" "<<std::get<1>(it->first).SpSnS().SSn()<<" "<<std::get<1>(it->first).SpSnS().SS()<<" "
          <<std::get<1>(it->first).upsilon()<<" "<<std::get<1>(it->first).Nex_omega()<<" "<<std::get<1>(it->first).x_omega().lm()<<" "
          <<std::get<1>(it->first).x_omega().mu()<<std::endl;
	  std::cout<<"tensor: "<<std::get<2>(it->first)[0]<<" "<<std::get<2>(it->first)[1]<<" "<<std::get<2>(it->first)[2]<<" "
		   <<std::get<2>(it->first)[3]<<" "<<std::get<2>(it->first)[4]<<std::endl;
	  std::cout<<"rho: "<<rho<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
	}
        break;
      }
    }
  }
  for(std::map<std::tuple<ONBasisState,ONBasisState,std::array<int,5>>,std::vector<double>>::iterator
                 it=neutron_obdrmes_recursive.begin(); it!=neutron_obdrmes_recursive.end(); it++){
    if(std::get<0>(it->first).Nex_omega()<std::get<1>(it->first).Nex_omega())continue;
    int rho=0;
    for(double rme : it->second){
      rho++;
      if(std::abs(rme)>1.0e-5){
        std::map<std::array<int,28>,double>::iterator it2=neutron_obdrmes_spncci.find({std::get<0>(it->first).gamma(),
          std::get<0>(it->first).Nex_sigma(),std::get<0>(it->first).x_sigma().lm(),std::get<0>(it->first).x_sigma().mu(),
          std::get<0>(it->first).SpSnS().SSp(),std::get<0>(it->first).SpSnS().SSn(),std::get<0>(it->first).SpSnS().SS(),
          std::get<0>(it->first).upsilon(),std::get<0>(it->first).Nex_omega(),std::get<0>(it->first).x_omega().lm(),
          std::get<0>(it->first).x_omega().mu(),
          std::get<1>(it->first).gamma(),
          std::get<1>(it->first).Nex_sigma(),std::get<1>(it->first).x_sigma().lm(),std::get<1>(it->first).x_sigma().mu(),
          std::get<1>(it->first).SpSnS().SSp(),std::get<1>(it->first).SpSnS().SSn(),std::get<1>(it->first).SpSnS().SS(),
          std::get<1>(it->first).upsilon(),std::get<1>(it->first).Nex_omega(),std::get<1>(it->first).x_omega().lm(),
          std::get<1>(it->first).x_omega().mu(),
          std::get<2>(it->first)[0],std::get<2>(it->first)[1],std::get<2>(it->first)[2],std::get<2>(it->first)[3],std::get<2>(it->first)[4],
          rho});
        if(it2==neutron_obdrmes_spncci.end()){
          std::cout<<"ERROR: didn't find neutron RME "<<rme<<" among those from spncci"<<std::endl;
          std::cout<<"bra: "<<std::get<0>(it->first).gamma()<<" "
          <<std::get<0>(it->first).Nex_sigma()<<" "<<std::get<0>(it->first).x_sigma().lm()<<" "<<std::get<0>(it->first).x_sigma().mu()<<" "
          <<std::get<0>(it->first).SpSnS().SSp()<<" "<<std::get<0>(it->first).SpSnS().SSn()<<" "<<std::get<0>(it->first).SpSnS().SS()<<" "
          <<std::get<0>(it->first).upsilon()<<" "<<std::get<0>(it->first).Nex_omega()<<" "<<std::get<0>(it->first).x_omega().lm()<<" "
          <<std::get<0>(it->first).x_omega().mu()<<std::endl;
          std::cout<<"ket: "<<std::get<1>(it->first).gamma()<<" "
          <<std::get<1>(it->first).Nex_sigma()<<" "<<std::get<1>(it->first).x_sigma().lm()<<" "<<std::get<1>(it->first).x_sigma().mu()<<" "
          <<std::get<1>(it->first).SpSnS().SSp()<<" "<<std::get<1>(it->first).SpSnS().SSn()<<" "<<std::get<1>(it->first).SpSnS().SS()<<" "
          <<std::get<1>(it->first).upsilon()<<" "<<std::get<1>(it->first).Nex_omega()<<" "<<std::get<1>(it->first).x_omega().lm()<<" "
          <<std::get<1>(it->first).x_omega().mu()<<std::endl;
          std::cout<<"tensor: "<<std::get<2>(it->first)[0]<<" "<<std::get<2>(it->first)[1]<<" "<<std::get<2>(it->first)[2]<<" "
                   <<std::get<2>(it->first)[3]<<" "<<std::get<2>(it->first)[4]<<std::endl;
          std::cout<<"rho: "<<rho<<std::endl;
          su3::finalize();
          return EXIT_FAILURE;
        }
        break;
      }
    }
  }
  std::cout<<"Benchmark with spncci finished"<<std::endl;
*/
  su3::finalize();

  return EXIT_SUCCESS;
}

#include <iostream>
#include <map>
#include <iterator>
#include <fstream>
#include "u3shell/u3spn_scheme.h"
#include "sp3rlib/sp3r.h"
#include "lgi/lgi.h"
#include "spncci/spncci_basis.h"
#include "lsu3shell/lsu3shell_basis.h"

// Kronecker product of U(3) x Sp x Sn x S and U(3)
MultiplicityTagged<u3shell::U3SPN>::vector KroneckerProduct(const u3shell::U3SPN& irrep, const u3::U3& omega)
{
  MultiplicityTagged<u3::U3>::vector u3_product=u3::KroneckerProduct(irrep.U3(),omega);
  MultiplicityTagged<u3shell::U3SPN>::vector product;
  product.reserve(u3_product.size());
  for(MultiplicityTagged<u3::U3> omega_tagged : u3_product)
  {    
    MultiplicityTagged<u3shell::U3SPN> irrep_tagged(u3shell::U3SPN(omega_tagged.irrep,irrep.Sp(),irrep.Sn(),irrep.S()),omega_tagged.tag);
    product.push_back(irrep_tagged);
  }
  return product;
}

struct Dimensions
{
  int total,cmf,LGI;
  Dimensions() : total(0), cmf(0), LGI(0) {} // default constructor (if omitted, compiler complains)
  Dimensions(int total_, int cmf_, int LGI_) : total(total_), cmf(cmf_), LGI(LGI_) {}
};

/* struct LS
{
  int L;
  HalfInt S;
  LS() : L(0), S(0) {} // default constructor
  LS(int L_, HalfInt S_) : L(L_), S(S_) {}
}; */

int main(int argc, char **argv)
{

  if(argc!=5){
    std::cout<<"Usage: "<<argv[0]<<" N_sigma_max N_max JJ_min JJ_max"<<std::endl;
    std::cout<<"where N_max is the N_max of the SpNCCI callculation for the target"<<std::endl;
    return EXIT_FAILURE;
  }   

  int N_sigma_max = std::stoi(argv[1]);
  int N_max = std::stoi(argv[2]);
  int JJmin = std::stoi(argv[3]);
  int JJmax = std::stoi(argv[4]);

  int A=6, N_0=2;
  HalfInt Nsigma_0=N_0+HalfInt(3,2)*(A-1);

  std::map<u3shell::U3SPN,Dimensions> dimensions_by_irrep;

/*
  // Digesting LSU3shell output
  std::string filename="Z1-N2-Nmax14_basis_tomas.dat";

  lsu3shell::LSU3ShellBasisTable lsu3_basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;

  lsu3shell::ReadLSU3ShellBasis(
        Nsigma_0, // input
        filename, // input
        lsu3_basis_table, // output
        basis_provenance, // output
        space // output
      );

  for(lsu3shell::LSU3ShellBasisGroupData omega_prime : lsu3_basis_table)
  {
    u3::U3 omega(omega_prime.omegaSPN.N()+HalfInt(3,2),omega_prime.omegaSPN.SU3());
    u3shell::U3SPN irrep(omega,omega_prime.omegaSPN.Sp(),omega_prime.omegaSPN.Sn(),omega_prime.omegaSPN.S());
    std::map<u3shell::U3SPN,Dimensions>::iterator it=dimensions_by_irrep.find(irrep);
    if (it==dimensions_by_irrep.end())
    {
      Dimensions irrep_dimensions(omega_prime.dim,omega_prime.dim,omega_prime.dim);
      dimensions_by_irrep[irrep]=irrep_dimensions;
    }
    else
    {
      dimensions_by_irrep[irrep].total+=omega_prime.dim;
      dimensions_by_irrep[irrep].cmf+=omega_prime.dim;
      dimensions_by_irrep[irrep].LGI+=omega_prime.dim;
    }
  }
*/

  // input from a file
  std::ifstream file ("Z3-N3-Nmax14_u3s-dim.dat");
  while (! file.eof())
  {
    int N_ex,twice_Sp,twice_Sn,twice_S,lambda,mu,dim_tot;
    file >> N_ex >> lambda >> mu >> twice_Sp >> twice_Sn >> twice_S >> dim_tot;
    if(N_ex>N_max)continue;
    HalfInt N(2*(N_ex+N_0)+3*A,2),Sp(twice_Sp,2),Sn(twice_Sn,2),S(twice_S,2); // N=N_ex+N_0+(3/2)*A, i.e. lab frame considered
    u3::SU3 x(lambda,mu);
    u3::U3 omega(N,x);
    u3shell::U3SPN irrep(omega,Sp,Sn,S);
    Dimensions irrep_dimensions(dim_tot,dim_tot,dim_tot);
//    map_irrep_dimensions.insert(std::make_pair(irrep,irrep_dimensions));
    dimensions_by_irrep[irrep]=irrep_dimensions;
  }
  file.close();

  std::ofstream outputfile;
  outputfile.open ("Sp3R_basis_output.out", std::ios::trunc);

  // calculation of CMF dimensions
  for(int N_ex=0; N_ex<=N_max; ++N_ex) // N_ex can start from 1.
  {
    // Calculate CMF dimension at given N_ex by removing spurious contributions coming from all lower subspaces.
    // We already know CMF dimensions for all lower subspaces.

    // Iterate through elements in dimensions_by_irrep that should be sorted by increasing U(1) label
    for(std::map<u3shell::U3SPN,Dimensions>::iterator it = dimensions_by_irrep.begin(); it != dimensions_by_irrep.end(); ++it)
    {
      // it->first is the key of an element of the map, it->second is the value.
      const u3shell::U3SPN& irrep_prime = it->first;
      const Dimensions& irrep_prime_dimensions = it->second;

      // For each irrep' with N_ex'<N_ex quanta look at spurious contribution at N_ex coming from lower CMF irrep'
      int N_ex_prime=(irrep_prime.N().TwiceValue()-2*N_0-3*A)/2;
      if(N_ex_prime>=N_ex)
        break; // Once N_ex_prime=N_ex, we end the loop because elements in dimensions_by_irrep are sorted by increasing U(1) label
      if(irrep_prime_dimensions.cmf==0)
        continue;
      int N_CM=N_ex-N_ex_prime;
      u3::SU3 x_CM(N_CM,0);
      u3::U3 omega_CM(N_CM,x_CM); // automatic type casting of N_CM from int to HalfInt allowed by constructor HalfInt(int)

      for(MultiplicityTagged<u3shell::U3SPN> irrep_tagged : KroneckerProduct(irrep_prime,omega_CM))
      {
        // For each irrep in product irrep' x N_CM(N_CM,0) eliminate spurious contribution to irrep coming from irrep'
        dimensions_by_irrep[irrep_tagged.irrep].cmf-=irrep_prime_dimensions.cmf*irrep_tagged.tag;
      }
    }
  }

  // little sanity check 1
  for(std::map<u3shell::U3SPN,Dimensions>::iterator it = dimensions_by_irrep.begin(); it != dimensions_by_irrep.end(); ++it)
  {
    if(it->second.cmf<0)
    {
       int Nex=(it->first.N().TwiceValue()-2*N_0-3*A)/2;
       outputfile << "little sanity check 1: Dimensions.cmf= " << it->second.cmf << " for Nex,lambda,mu,Sp,Sn,S= " << std::endl;
       outputfile << Nex << "     " << it->first.SU3().lambda() << "     " << it->first.SU3().mu() << "     " << it->first.Sp() << "     " << it->first.Sn() << "     " << it->first.S() << std::endl;
    }
  }

  // initialization of "LGI dimensions/multiplicities"
  for(std::map<u3shell::U3SPN,Dimensions>::iterator it = dimensions_by_irrep.begin(); it != dimensions_by_irrep.end(); ++it)
  {
    it->second.LGI=it->second.cmf;
  }

  std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector> U3SPN_brancing_by_Sp3R_irrep;

  // construction of Sp3R irreps
  for(std::map<u3shell::U3SPN,Dimensions>::iterator it = dimensions_by_irrep.begin(); it != dimensions_by_irrep.end(); ++it)
  {
    const u3shell::U3SPN& sigma = it->first;
    const Dimensions& sigma_dimensions = it->second;
    int N_sigma=(sigma.N().TwiceValue()-2*N_0-3*A)/2;
    if(N_sigma>N_sigma_max)
      break;
    if(sigma_dimensions.LGI==0)
      continue;
    MultiplicityTagged<u3shell::U3SPN> irrep_family(sigma,sigma_dimensions.LGI);
    MultiplicityTagged<u3shell::U3SPN>::vector U3SPN_irreps;
    int N_n_max=N_max-(sigma.N().TwiceValue()-2*N_0-3*A)/2;
    std::vector<u3::U3> raising_polynomial_labels = u3boson::RaisingPolynomialLabels(N_n_max);

    for(u3::U3 n : raising_polynomial_labels)
    {
      for(MultiplicityTagged<u3shell::U3SPN> omega_tagged : KroneckerProduct(sigma,n))
      {
        U3SPN_irreps.push_back(omega_tagged);
        dimensions_by_irrep[omega_tagged.irrep].LGI-=omega_tagged.tag*irrep_family.tag;
      }
    }
    U3SPN_brancing_by_Sp3R_irrep[irrep_family]=U3SPN_irreps;
  }

  // a little sanity check
  for(std::map<u3shell::U3SPN,Dimensions>::iterator it = dimensions_by_irrep.begin(); it != dimensions_by_irrep.end(); ++it)
  {
    if(it->second.LGI<0)
    {
      outputfile << "Dimensions.LGI=" << it->second.LGI << std::endl;
    }
  }

  //Removal of multiple occurencies (DOESN'T WORK BUT DOESN'T BREAK ANYTHING!)
  for(std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector>::iterator
      it = U3SPN_brancing_by_Sp3R_irrep.begin(); it != U3SPN_brancing_by_Sp3R_irrep.end(); ++it)
  {
    std::vector<int> indices;
    indices.clear();
    int i=0;
    for(MultiplicityTagged<u3shell::U3SPN> omega_tagged : it->second)
    {
      ++i;
      MultiplicityTagged<u3shell::U3SPN>::vector::iterator iter;
      iter = find (it->second.begin()+i, it->second.end(), omega_tagged);
      if (iter != it->second.end())
      {
        (*iter).tag+=omega_tagged.tag;
        indices.push_back (i);
      }
    }
    int k=0;
    for(int j : indices)
    {
      ++k;
      it->second.erase (it->second.begin()+j-k);
    }
  }

  //output
/*  std::map<LS,int> dimension_by_LS;
  for(std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector>::iterator
      it = U3SPN_brancing_by_Sp3R_irrep.begin(); it != U3SPN_brancing_by_Sp3R_irrep.end(); ++it)
  {
    int Nsigmaex=(it->first.irrep.N().TwiceValue()-2*N_0-3*A)/2;
    if(Nsigmaex==1||Nsigmaex==3||Nsigmaex==5||Nsigmaex==7||Nsigmaex==9||Nsigmaex==11||Nsigmaex==13||Nsigmaex==15)
      continue;
    for(MultiplicityTagged<u3shell::U3SPN> omega_tagged : it->second)
    {
      for(int L=0; L<=omega_tagged.irrep.SU3().lambda()+omega_tagged.irrep.SU3().mu(); ++L)
      {
        int kappa_max=inner_multiplicity(omega_tagged.irrep.SU3(),L);
        if(kappa_max==0)
          continue;
        LS LS_pair(L,omega_tagged.irrep.S());
        int indicator = 0;
        for(std::map<LS,int>::iterator iter = dimension_by_LS.begin(); iter != dimension_by_LS.end(); ++iter)
        {
          if(iter->first == LS_pair)
          {
            iter->second+=it->first.tag*omega_tagged.tag*kappa_max;
            indicator = 1;
            break;
          }
        }
        if(indicator == 0)
          dimension_by_LS[LS_pair]=it->first.tag*omega_tagged.tag*kappa_max;
      }
    }
  } */

  std::map<HalfInt,int> dimension_by_J;
  for(std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector>::iterator
      it = U3SPN_brancing_by_Sp3R_irrep.begin(); it != U3SPN_brancing_by_Sp3R_irrep.end(); ++it)
  {
    int Nsigmaex=(it->first.irrep.N().TwiceValue()-2*N_0-3*A)/2;
    if(Nsigmaex % 2 != 0) // odd values of Nsigmaex skipped
      continue;
    for(MultiplicityTagged<u3shell::U3SPN> omega_tagged : it->second)
    {
      for(auto L_tagged : u3::BranchingSO3(omega_tagged.irrep.SU3()))
      {
        int L=L_tagged.irrep, kappa_max=L_tagged.tag;
        for(int J_twice=abs(2*L-omega_tagged.irrep.S().TwiceValue()); J_twice<=2*L+omega_tagged.irrep.S().TwiceValue(); J_twice=J_twice+2)
        {
          HalfInt J(J_twice,2);          
          int indicator = 0;
          for(std::map<HalfInt,int>::iterator iter = dimension_by_J.begin(); iter != dimension_by_J.end(); ++iter)
          {
            if(iter->first == J)
            {
              iter->second+=it->first.tag*omega_tagged.tag*kappa_max;
              indicator = 1;
              break;
            }
          }
          if(indicator == 0)
            dimension_by_J[J]=it->first.tag*omega_tagged.tag*kappa_max;
        }
      }
    }
  }

  std::map<u3::U3S,int> dimension_by_U3S_irrep;
  for(std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector>::iterator
      it = U3SPN_brancing_by_Sp3R_irrep.begin(); it != U3SPN_brancing_by_Sp3R_irrep.end(); ++it)
  {
    int Nsigmaex=(it->first.irrep.N().TwiceValue()-2*N_0-3*A)/2;
    if(Nsigmaex % 2 != 0)
      continue;
    for(MultiplicityTagged<u3shell::U3SPN> omega_tagged : it->second)
    {
      u3::U3S omegaU3S(omega_tagged.irrep.U3S());
      int indicator = 0;
      for(std::map<u3::U3S,int>::iterator iter = dimension_by_U3S_irrep.begin(); iter != dimension_by_U3S_irrep.end(); ++iter)
      {
        if(iter->first == omegaU3S)
        {
          iter->second+=it->first.tag*omega_tagged.tag;
          indicator = 1;
          break;
        }
      }
      if(indicator == 0)
        dimension_by_U3S_irrep[omegaU3S]=it->first.tag*omega_tagged.tag;
    }
  }

  std::ofstream myfile;
  myfile.open ("basis_counting.out", std::ios::trunc);

  myfile << "[Sp(3,R) families]" << std::endl;
  myfile << "subspace_index N_ex N lambda mu dim" << std::endl;
  int sigma_index=-1;
  for(std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector>::iterator
      it = U3SPN_brancing_by_Sp3R_irrep.begin(); it != U3SPN_brancing_by_Sp3R_irrep.end(); ++it)
  {
    int Nsigmaex=(it->first.irrep.N().TwiceValue()-2*N_0-3*A)/2;
    if(Nsigmaex % 2 != 0)
      continue;
    ++sigma_index;
    HalfInt Nsigma(it->first.irrep.N().TwiceValue()-3,2);
    myfile << sigma_index << "     " << Nsigmaex << "     " << Nsigma << "     " << it->first.irrep.SU3().lambda() << "     "
           << it->first.irrep.SU3().mu() << "     " << it->first.tag << std::endl;
  }
  myfile << std::endl;

  myfile << "[Sp(3,R) x S_p x S_n x S listing]" << std::endl;
  myfile << "subspace_index N_ex N lambda mu S_p S_n S dim" << std::endl;
  sigma_index=-1;
  for(std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector>::iterator
      it = U3SPN_brancing_by_Sp3R_irrep.begin(); it != U3SPN_brancing_by_Sp3R_irrep.end(); ++it)
  {
    int Nsigmaex=(it->first.irrep.N().TwiceValue()-2*N_0-3*A)/2;
    if(Nsigmaex % 2 != 0)
      continue;
    ++sigma_index;
    HalfInt Nsigma(it->first.irrep.N().TwiceValue()-3,2);
    myfile << sigma_index << "     " << Nsigmaex << "     " << Nsigma << "     " << it->first.irrep.SU3().lambda() << "     "
           << it->first.irrep.SU3().mu() << "     " << it->first.irrep.Sp() << "     " << it->first.irrep.Sn() << "     "
           << it->first.irrep.S() << "     " << it->first.tag << std::endl;
  }
  myfile << std::endl;

  myfile << "[J listing]" << std::endl;
  myfile << "subspace_index J dim" << std::endl;
long int numb=0;
  int J_index=-1;
  for(std::map<HalfInt,int>::iterator it = dimension_by_J.begin(); it != dimension_by_J.end(); ++it)
  {
    ++J_index;
    myfile << J_index << "     " << it->first << "     " << it->second << std::endl;
long int a=it->second;
    numb+=a*a;
  }
std::cout<<"nnz: "<<numb<<std::endl;
  myfile << std::endl;

  myfile << "[U3S listing]" << std::endl;
  myfile << "subspace_index Nex omega_N omega_lambda omega_mu S dim" << std::endl;
  int omega_index=-1;
  int TotalU3Sdim=0;
std::vector<std::tuple<u3::SU3,int>> SU3SSirreps;
  for(std::map<u3::U3S,int>::iterator it = dimension_by_U3S_irrep.begin(); it != dimension_by_U3S_irrep.end(); ++it)
  {
    ++omega_index;
    int Nex=(it->first.U3().N().TwiceValue()-2*N_0-3*A)/2;
    HalfInt N(it->first.U3().N().TwiceValue()-3,2);
    myfile << omega_index << "     " << Nex << "     " << N << "     " << it->first.SU3().lambda() << "     "
           << it->first.SU3().mu() << "     " << it->first.S() << "     " << it->second << std::endl;
    TotalU3Sdim=TotalU3Sdim+it->second;
std::tuple<u3::SU3,int> ser={it->first.SU3(),it->first.S().TwiceValue()};
std::vector<std::tuple<u3::SU3,int>>::iterator iter=std::find(SU3SSirreps.begin(),SU3SSirreps.end(),ser);
if(iter==SU3SSirreps.end())SU3SSirreps.push_back(ser);
  }
  myfile << std::endl;
/*
long int numb_rgm=0;
long int numb_obrho=0;
long int numb_tbrho=0;
for(std::tuple<u3::SU3,int> xSS1 : SU3SSirreps){
// u3::SU3 x1=std::get<0>(xSS1);
 int SS1=std::get<1>(xSS1);
 for(auto L1tagged : u3::BranchingSO3(std::get<0>(xSS1))){
  int L1=int(L1tagged.irrep);
  int kappa1max=L1tagged.tag;
  if(std::abs(2*L1-SS1)>2 || 2>2*L1+SS1)continue;
  for(int N=0; N<=N_max+1; N++){
   for(auto xtagged : u3::KroneckerProduct(std::get<0>(xSS1),u3::SU3(N,0))){
    for(auto Ltagged : u3::BranchingSO3(xtagged.irrep)){
     int L=int(Ltagged.irrep);
     int kappamax=Ltagged.tag;
     for(int SS=std::abs(SS1-1); SS<=SS1+2; SS+=2){
      for(int JJ=std::max(JJmin,std::abs(2*L-SS)); JJ<=std::min(JJmax,2*L+SS); JJ+=2){
       numb_rgm+=kappa1max*kappamax;
      }
     }
    }
   }
  }
  for(std::tuple<u3::SU3,int> xSS1p : SU3SSirreps){
   int SS1p=std::get<1>(xSS1p);
   for(auto L1ptagged : u3::BranchingSO3(std::get<0>(xSS1p))){
    int L1p=int(L1ptagged.irrep);
    int kappa1pmax=L1ptagged.tag;
    if(std::abs(2*L1p-SS1p)>2 || 2>2*L1p+SS1p)continue;
    for(int N=0; N<=N_max+1; N++){
     int Npmin=std::max(0,N-N_max);
     if((N-Npmin)%2!=0)Npmin++; // so that the OBD operator doesn't change parity
     for(int Np=Npmin; Np<=std::min(N_max+1,N+N_max); Np+=2){
      for(auto x0 : u3::KroneckerProduct(u3::SU3(N,0),u3::SU3(0,Np))){
       int rhomax=u3::OuterMultiplicity(std::get<0>(xSS1),x0.irrep,std::get<0>(xSS1p));
       if(rhomax==0)continue;
       for(int S0=0; S0<=1; S0++){
        if(std::abs(SS1-2*S0)>SS1p || SS1p>SS1+2*S0)continue;
         numb_obrho+=kappa1max*kappa1pmax*rhomax;
       }
      }
     }
    }
    for(int N1=0; N1<=N_max+1; N1++){
     for(int N2=0; N2<=std::min(N_max+1,2+N_max-N1); N2++){
      for(int N3=0; N3<=N_max; N3++){
       int N123=N1+N2-N3;
       int N4min=std::max(0,N123-N_max);
       if((N123-N4min)%2!=0)N4min++; // so that the TBD operator doesn't change parity
       for(int N4=N4min; N4<=std::min(std::min(N_max+1,2+N_max-N3),N123+N_max); N4+=2){
        for(auto xf : u3::KroneckerProduct(u3::SU3(N1,0),u3::SU3(N2,0))){
         for(int Sf=0; Sf<=1; Sf++){
          for(auto xi : u3::KroneckerProduct(u3::SU3(0,N3),u3::SU3(0,N4))){
           for(int Si=0; Si<=1; Si++){
            for(auto x0 : u3::KroneckerProduct(xf.irrep,xi.irrep)){
             int rho0max=x0.tag;
             int rhomax=u3::OuterMultiplicity(std::get<0>(xSS1),x0.irrep,std::get<0>(xSS1p));
             if(rhomax==0)continue;
             for(int S0=std::abs(Sf-Si); S0<=Sf+Si; S0++){
              if(std::abs(SS1-2*S0)>SS1p || SS1p>SS1+2*S0)continue;
              numb_tbrho+=kappa1max*kappa1pmax*rhomax*rho0max;
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
std::cout<<"Number of SA-RGM basis states: "<<numb_rgm<<std::endl;
std::cout<<"Number of OB rhos: "<<numb_obrho<<std::endl;
std::cout<<"Number of TB rhos: "<<numb_tbrho<<std::endl;
*/
/*
long int numb_obrho=0;
long int numb_tbrho=0;
for(std::tuple<u3::SU3,int> xSS1 : SU3SSirreps){
 int SS1=std::get<1>(xSS1);
 for(auto L1tagged : u3::BranchingSO3(std::get<0>(xSS1))){
  int L1=int(L1tagged.irrep);
  int kappa1max=L1tagged.tag;
  if(std::abs(2*L1-SS1)>2 || 2>2*L1+SS1)continue;
  for(std::tuple<u3::SU3,int> xSS1p : SU3SSirreps){
   int SS1p=std::get<1>(xSS1p);
   for(auto L1ptagged : u3::BranchingSO3(std::get<0>(xSS1p))){
    int L1p=int(L1ptagged.irrep);
    int kappa1pmax=L1ptagged.tag;
    if(std::abs(2*L1p-SS1p)>2 || 2>2*L1p+SS1p)continue;
    for(int N=0; N<=N_max+1; N++){
     int Npmin=std::max(0,N-N_max);
     if((N-Npmin)%2!=0)Npmin++; // so that the OBD operator doesn't change parity
     for(int Np=Npmin; Np<=std::min(N_max+1,N+N_max); Np+=2){
      for(auto x0 : u3::KroneckerProduct(u3::SU3(N,0),u3::SU3(0,Np))){
       int rhomax=u3::OuterMultiplicity(std::get<0>(xSS1),x0.irrep,std::get<0>(xSS1p));
       if(rhomax==0)continue;
       for(int S0=0; S0<=1; S0++){
        if(std::abs(SS1-2*S0)>SS1p || SS1p>SS1+2*S0)continue;
         numb_obrho+=kappa1max*kappa1pmax*rhomax;
       }
      }
     }
    }
    for(int N1=0; N1<=N_max+1; N1++){
     for(int N2=0; N2<=std::min(N_max+1,2+N_max-N1); N2++){
      for(int N3=0; N3<=N_max; N3++){
       int N123=N1+N2-N3;
       int N4min=std::max(0,N123-N_max);
       if((N123-N4min)%2!=0)N4min++; // so that the TBD operator doesn't change parity
       for(int N4=N4min; N4<=std::min(std::min(N_max+1,2+N_max-N3),N123+N_max); N4+=2){
        for(auto xf : u3::KroneckerProduct(u3::SU3(N1,0),u3::SU3(N2,0))){
         for(int Sf=0; Sf<=1; Sf++){
          for(auto xi : u3::KroneckerProduct(u3::SU3(0,N3),u3::SU3(0,N4))){
           for(int Si=0; Si<=1; Si++){
            for(auto x0 : u3::KroneckerProduct(xf.irrep,xi.irrep)){
             int rho0max=x0.tag;
             int rhomax=u3::OuterMultiplicity(std::get<0>(xSS1),x0.irrep,std::get<0>(xSS1p));
             if(rhomax==0)continue;
             for(int S0=std::abs(Sf-Si); S0<=Sf+Si; S0++){
              if(std::abs(SS1-2*S0)>SS1p || SS1p>SS1+2*S0)continue;
              numb_tbrho+=kappa1max*kappa1pmax*rhomax*rho0max;
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
std::cout<<"Number of OB rhos: "<<numb_obrho<<std::endl;
std::cout<<"Number of TB rhos: "<<numb_tbrho<<std::endl;
*/
/*
long int numb_tbrho=0;
for(std::tuple<u3::SU3,int> xSS1 : SU3SSirreps){
 int SS1=std::get<1>(xSS1);
 for(auto L1tagged : u3::BranchingSO3(std::get<0>(xSS1))){
  int L1=int(L1tagged.irrep);
  int kappa1max=L1tagged.tag;
  if(std::abs(2*L1-SS1)>2 || 2>2*L1+SS1)continue;
  for(std::tuple<u3::SU3,int> xSS1p : SU3SSirreps){
   int SS1p=std::get<1>(xSS1p);
   for(auto L1ptagged : u3::BranchingSO3(std::get<0>(xSS1p))){
    int L1p=int(L1ptagged.irrep);
    int kappa1pmax=L1ptagged.tag;
    if(std::abs(2*L1p-SS1p)>2 || 2>2*L1p+SS1p)continue;

    for(int N1=0; N1<=N_max+1; N1++){
     for(int N2=0; N2<=std::min(N_max+1,2+N_max-N1); N2++){
      for(int N3=0; N3<=N_max; N3++){
       int N123=N1+N2-N3;
       int N4min=std::max(0,N123-N_max);
       if((N123-N4min)%2!=0)N4min++; // so that the TBD operator doesn't change parity
       for(int N4=N4min; N4<=std::min(std::min(N_max+1,2+N_max-N3),N123+N_max); N4+=2){
        for(auto xf : u3::KroneckerProduct(u3::SU3(N1,0),u3::SU3(N2,0))){
         for(int Sf=0; Sf<=1; Sf++){
          for(auto xi : u3::KroneckerProduct(u3::SU3(0,N3),u3::SU3(0,N4))){
           for(int Si=0; Si<=1; Si++){
            for(auto x0 : u3::KroneckerProduct(xf.irrep,xi.irrep)){
	     int rho0max=x0.tag;
             int rhomax=u3::OuterMultiplicity(std::get<0>(xSS1),x0.irrep,std::get<0>(xSS1p));
             if(rhomax==0)continue;
             for(int S0=std::abs(Sf-Si); S0<=Sf+Si; S0++){
              if(std::abs(SS1-2*S0)>SS1p || SS1p>SS1+2*S0)continue;
              numb_tbrho+=kappa1max*kappa1pmax*rhomax*rho0max;
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
std::cout<<"Number of TB rhos: "<<numb_tbrho<<std::endl;
*/
  myfile << "[baby SpNCCI listing]" << std::endl;
  myfile << "subspace_index irrep_family_index Nsigmaex sigma.N sigma.lambda sigma.mu Sp Sn S" << std::endl
         << "Nomegaex omega.N omega.lambda omega.mu sigma_multiplicity omega_multiplicity dim" << std::endl;
  sigma_index=-1;
  omega_index=-1;
  for(std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector>::iterator
      it = U3SPN_brancing_by_Sp3R_irrep.begin(); it != U3SPN_brancing_by_Sp3R_irrep.end(); ++it)
  {
    int Nsigmaex=(it->first.irrep.N().TwiceValue()-2*N_0-3*A)/2;
    if(Nsigmaex % 2 != 0)
      continue;
    ++sigma_index;
    HalfInt Nsigma(it->first.irrep.N().TwiceValue()-3,2);
    for(MultiplicityTagged<u3shell::U3SPN> omega_tagged : it->second)
    {
      ++omega_index;
      int Nomegaex=(omega_tagged.irrep.N().TwiceValue()-2*N_0-3*A)/2;
      HalfInt Nomega(omega_tagged.irrep.N().TwiceValue()-3,2);
      myfile << omega_index << "     " << sigma_index << "     " << Nsigmaex << "     " << Nsigma << "     "
             << it->first.irrep.SU3().lambda() << "     " << it->first.irrep.SU3().mu() << "     " << it->first.irrep.Sp() << "     "
             << it->first.irrep.Sn() << "     " << it->first.irrep.S() << "     " << Nomegaex << "     "
             << Nomega << "     " << omega_tagged.irrep.SU3().lambda() << "     "
             << omega_tagged.irrep.SU3().mu() << "     " << it->first.tag << "     " << omega_tagged.tag << "     "
             << omega_tagged.tag*it->first.tag << std::endl;
    }
  }
  myfile.close();

//*****************************************************************************************

  // Construction of MultiplicityTaggedLGIVector
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  for(std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector>::iterator
      it = U3SPN_brancing_by_Sp3R_irrep.begin(); it != U3SPN_brancing_by_Sp3R_irrep.end(); ++it)
  {
//    int Nsigmaex=(it->first.irrep.N().TwiceValue()-2*N_0-3*A)/2;
    HalfInt Nsigma=it->first.irrep.N();
    int Nsigmaex=int(Nsigma-N_0-HalfInt(3,2)*A);
    if(Nsigmaex % 2 != 0)
      continue; // skip odd Nsigmaex
// constuct new U3SPN object sigma with correct (intrinsic) N and pass it into LGI instead of it->first.irrep!
    u3::U3 sigma_U3(Nsigma-HalfInt(3,2),it->first.irrep.SU3());
    u3shell::U3SPN sigma_U3SPN(sigma_U3,it->first.irrep.Sp(),it->first.irrep.Sn(),it->first.irrep.S());
    lgi::LGI sigma(sigma_U3SPN,Nsigmaex);
    MultiplicityTagged<lgi::LGI> sigma_tagged(sigma,it->first.tag);
    lgi_vector.push_back(sigma_tagged);
  }

  // truncator
//  HalfInt Nsigma_0=N_0+HalfInt(3,2)*(A-1); // Nsigma_0=N_0+(3/2)*(A-1)
  spncci::NmaxTruncator truncator(Nsigma_0,N_max);

  // call GenerateSpNCCISpace
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;
  spncci::GenerateSpNCCISpace(
      lgi_vector, // input
      truncator, // input
      spncci_space, // output
      sigma_irrep_map // output
    );

  // print dimensions
  outputfile << "TotalU3Subspaces: " << spncci::TotalU3Subspaces(spncci_space) << std::endl;
  outputfile << "TotalDimensionU3S: " << spncci::TotalDimensionU3S(spncci_space) << std::endl;
  outputfile << "TotalDimensionU3LS: " << spncci::TotalDimensionU3LS(spncci_space) << std::endl;
  HalfInt J(2);
  outputfile << "TotalDimensionU3LSJConstrained (J=2): " << spncci::TotalDimensionU3LSJConstrained(spncci_space,J) << std::endl;
  outputfile << "TotalDimensionU3LSJAll: " << spncci::TotalDimensionU3LSJAll(spncci_space) << std::endl;
  outputfile << std::endl;
  outputfile << "total U3S dimension: " << TotalU3Sdim << std::endl;

  outputfile.close();

  return 0;
}

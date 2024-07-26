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

int main()
{
  int Z, N, N_max;
  std::cout << "Enter the number of protons" << std::endl;
  std::cin >> Z;
  std::cout << "Enter the number of neutrons" << std::endl;
  std::cin >> N;
  std::cout << "Enter N_max (this must be greater than or equal to the maximal N_sigma_ex in the list of LGIs in irrep_families.dat)" << std::endl;
  std::cin >> N_max;
  std::cout << "Number of protons: " << Z << ", Number of neutrons: " << N << ", N_max = " << N_max << std::endl;
  if(Z>40 || N>40)
  {
    std::cout << "ERROR: Code works for nuclei with Z<=40 and N<=40 only.";
  }
  else
  {
    int A=Z+N, N_0=0;
    Z=Z-2;
    if(Z>0)
    {
      N_0=N_0+Z;
    }
    Z=Z-6;
    if(Z>0)
    {
      N_0=N_0+Z;
    }
    Z=Z-12;
    if(Z>0)
    {
      N_0=N_0+Z;
    }
    N=N-2;
    if(N>0)
    { 
      N_0=N_0+N;
    }
    N=N-6;
    if(N>0)
    { 
      N_0=N_0+N;
    }
    N=N-12;
    if(N>0)
    { 
      N_0=N_0+N;
    }
    std::cout << "Number of oscillator quanta in the lowest configuration: N_0 = " << N_0 << std::endl;

/*
    HalfInt Nsigma_0=N_0+HalfInt(3,2)*(A-1); // N=N_ex+N_0+(3/2)*(A-1)
    lgi::MultiplicityTaggedLGIVector lgi_vector;

    // input from a file containing irrep families
    std::ifstream file ("irrep_families.dat");
    int N_ex_max=0;
    while (! file.eof())
    {
      int N_ex,lambda,mu,twice_Sp,twice_Sn,twice_S,dim;
      file >> N_ex >> lambda >> mu >> twice_Sp >> twice_Sn >> twice_S >> dim;
      N_ex_max=std::max(N_ex_max,N_ex);
      HalfInt Nsigma=N_ex+N_0+HalfInt(3,2)*(A-1),Sp(twice_Sp,2),Sn(twice_Sn,2),S(twice_S,2);
      u3::SU3 x(lambda,mu);
      u3::U3 sigma_U3(Nsigma,x);
      u3shell::U3SPN sigma_U3SPN(sigma_U3,Sp,Sn,S);
      lgi::LGI sigma(sigma_U3SPN,N_ex);
      MultiplicityTagged<lgi::LGI> sigma_tagged(sigma,dim);
      lgi_vector.push_back(sigma_tagged);
    }
    file.close();

    std::cout << "Dimensions:" << std::endl;
    std::cout << "Nmax     TotalDimensionU3LSJAll     TotalDimensionU3LSJConstrained(J=" << J << ")" << std::endl;
    for(int Nmax=N_ex_max; Nmax<=20; Nmax=Nmax+2)
    {
      spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
      spncci::SpNCCISpace spncci_space;
      spncci::SigmaIrrepMap sigma_irrep_map;
      spncci::GenerateSpNCCISpace(
        lgi_vector, // input
        truncator, // input
        spncci_space, // output
        sigma_irrep_map // output
      );
      std::cout << Nmax << "     " << spncci::TotalDimensionU3LSJAll(spncci_space) << "     " << spncci::TotalDimensionU3LSJConstrained(spncci_space,J) << std::endl;
    }
*/

    std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector> U3SPN_branching_by_Sp3R_irrep;
    std::ifstream file ("irrep_families.dat");
    while (! file.eof())
    {
      int N_ex,lambda,mu,twice_Sp,twice_Sn,twice_S,dim_sigma;
      file >> N_ex >> lambda >> mu >> twice_Sp >> twice_Sn >> twice_S >> dim_sigma;
   
      HalfInt HO(2*(N_0+N_ex)+3*(A-1),2),Sp(twice_Sp,2),Sn(twice_Sn,2),S(twice_S,2);
      u3::SU3 x(lambda,mu);
      u3::U3 omega(HO,x);
      u3shell::U3SPN sigma(omega,Sp,Sn,S);
      MultiplicityTagged<u3shell::U3SPN> irrep_family(sigma,dim_sigma);

      MultiplicityTagged<u3shell::U3SPN>::vector U3SPN_irreps;
      std::vector<u3::U3> raising_polynomial_labels = sp3r::RaisingPolynomialLabels(N_max-N_ex);
      for(u3::U3 n : raising_polynomial_labels)
      {
        for(MultiplicityTagged<u3shell::U3SPN> omega_tagged : KroneckerProduct(sigma,n))
        {
          U3SPN_irreps.push_back(omega_tagged);
        }
      }

      U3SPN_branching_by_Sp3R_irrep[irrep_family]=U3SPN_irreps;
    }
    file.close();

    std::map<HalfInt,int> dimension_by_J;
    for(std::map<MultiplicityTagged<u3shell::U3SPN>,MultiplicityTagged<u3shell::U3SPN>::vector>::iterator
        it = U3SPN_branching_by_Sp3R_irrep.begin(); it != U3SPN_branching_by_Sp3R_irrep.end(); ++it)
    {
      for(MultiplicityTagged<u3shell::U3SPN> omega_tagged : it->second)
      {
        for(MultiplicityTagged<int> L_tagged : u3::BranchingSO3(omega_tagged.irrep.SU3()))
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

    std::cout << "J     dimension" << std::endl;
    for(std::map<HalfInt,int>::iterator it = dimension_by_J.begin(); it != dimension_by_J.end(); ++it)
    {
      std::cout << it->first << "     " << it->second << std::endl;
    }

  }
  return 0;
}

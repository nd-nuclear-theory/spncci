/****************************************************************
  u4_test.cpp

  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT

  + 01/15/25 (aem): Created.
****************************************************************/

#include "su4lib/u4.h"
#include "su4lib/u4_test_data.h"
int main()
{
  u4::U4 f(8u,3u,3u,0u);
  std::cout<<f.Str()<<std::endl;
  for(int N=7; N<=9; N++)
    std::cout<<"Has conjuate for N="<<N<<": "<<f.HasSnConjugate(N)<<std::endl;

  u4::U4Conjugate conjugate(f,12);

  std::cout<<conjugate.Str()<<std::endl;

  // std::cout<<u4::OuterMultiplicity(f,{1u,1u,1u,1u},f)<<std::endl;
  // std::cout<<u4::OuterMultiplicity(f,{3u,2u,1u,1u},{8u,6u,4u,3u})<<std::endl;
  std::cout<<u4::OuterMultiplicity(f,{3u,2u,1u,1u},{10u,5u,5u,1u})<<std::endl;

  unsigned int A=12;
  if(true)
  {

  auto product = u4::KroneckerProduct(f, {4u,3u,2u,1u});
  for(const auto& [irrep,mult] : product)
    std::cout<<irrep.Str()<<"  "<<mult<<std::endl;
  // auto tableaux = u4::GenerateU4Tableaux(A);
  // for(const auto& tableau : tableaux){
  //   std::cout<<tableau.Str()<<"  "<<tableau.HasSnConjugate(6u)<<std::endl;
  // }
  }

  u4::U4 h(17,{8,  3,  2}); 
  std::cout<<h.Str()<<" [9,4,3,1]"<< std::endl;

  auto m = u4::OuterMultiplicity({8,3,2,0}, {7,4,2,0}, {8,7,6,5});
  std::cout<<"m "<<m<<std::endl;

  if(true)
  {
    auto product_test = u4::KroneckerProduct(u4::test::f, u4::test::g);
    for(const auto& irrep_tag : product_test)
    {
      const bool not_found = u4::test::test_irreps_fg.find(irrep_tag) == u4::test::test_irreps_fg.end();
      if(not_found)
        std::cout<<"irrep not in product: "<<irrep_tag.irrep.Str()<<" "<<irrep_tag.tag<<std::endl;

      if(product_test.size()!=u4::test::test_irreps_fg.size())
        std::cout<<"missing irreps"<<std::endl;
    }
    std::cout<<"validation complete"<<std::endl;
  }

  if(true)
  {
    auto product_test = u4::KroneckerProduct(u4::test::a, u4::test::b);
    for(const auto& irrep_tag : product_test)
    {
      const bool not_found = u4::test::test_irreps_ab.find(irrep_tag) == u4::test::test_irreps_ab.end();
      if(not_found)
        std::cout<<"irrep not in product: "<<irrep_tag.irrep.Str()<<" "<<irrep_tag.tag<<std::endl;

      if(product_test.size()!=u4::test::test_irreps_ab.size())
        std::cout<<"missing irreps"<<std::endl;
    }
    std::cout<<"validation complete"<<std::endl;
  }

}

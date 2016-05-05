/****************************************************************
  u3coef_test.cpp

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created based on prototype u3.py and 
  T. Dytrych CSU3Master.
****************************************************************/

#include "am/halfint.h"
#include "utilities/utilities.h"
#include "utilities/multiplicity_tagged.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include <map>

void basic_test()
{
  u3::SU3 x1(2,0);
  u3::SU3 x2(2,0);
  u3::SU3 x12(4,0);
  u3::SU3 x3(2,0);
  u3::SU3 x4(1,0);
  u3::SU3 x23(4,0);
  u3::SU3 x34(3,0);
  u3::SU3 x13(4,0);
  u3::SU3 x24(3,0);
  u3::SU3 x(6,0);
  u3::SU3 xx(7,0);


  // U(x1,x2,x,x3,x12,r12,r12_3,x23,r23,r1_23);
  std::cout << U(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;
  std::cout << U(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;
  std::cout << Z(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;

  std::cout << W(x1,1,2,x2,1,2,x12,1,4,1) << std::endl;
  std::cout << W(x1,1,2,x2,1,2,x12,1,2,1) << std::endl;
  std::cout << W(x1,1,2,x2,1,2,x12,1,0,1) << std::endl;

  std::cout << Unitary9LambdaMu(
                                x1,  x2,  x12, 1,
                                x3,  x4,  x34, 1,
                                x13, x24, xx,  1,
                                1,    1,   1
                                )<<std::endl;

}

void block_test()
{

  u3::SU3 x1(2,0);
  u3::SU3 x2(2,0);
  u3::SU3 x12(4,0);
  u3::SU3 x3(2,0);
  u3::SU3 x23(4,0);
  u3::SU3 x(6,0);

  // block access
  std::cout << "U block test" << std::endl;
  u3::UCoefLabels labels(x1, x2, x, x3, x12, x23);
  u3::UCoefBlock block(labels);
  int r12_max, r12_3_max, r23_max, r1_23_max;
  std::tie(r12_max, r12_3_max, r23_max, r1_23_max) = block.Key();
  std::cout << "multiplicities " << r12_max << " " << r12_3_max << " " << r23_max << " " << r1_23_max << std::endl;
  std::cout << block.GetCoef(1,1,1,1) << std::endl;
  std::cout << std::endl;
  
}

void iteration_test()
{

  // label set
  std::vector<u3::UCoefLabels> label_set;


  u3::SU3 x1(10,0);
  u3::SU3 x2(8,5);
  u3::SU3 x3(5,8);
  MultiplicityTagged<u3::SU3>::vector x12_vector=KroneckerProduct(x1,x2);
  MultiplicityTagged<u3::SU3>::vector x23_vector=KroneckerProduct(x2,x3);
  for(int i=0; i<x12_vector.size(); i++)
    {
      u3::SU3 x12=x12_vector[i].irrep;
      MultiplicityTagged<u3::SU3>::vector x_vector=KroneckerProduct(x12,x3);
      for(int j=0; j<x23_vector.size(); j++)
        {
          u3::SU3 x23=x23_vector[j].irrep;
          for(int k=0; k<x_vector.size(); k++)
            {
              u3::SU3 x=x_vector[k].irrep;
              if(u3::OuterMultiplicity(x1,x23,x)>0)
                label_set.push_back(u3::UCoefLabels(x1,x2,x,x3,x12,x23));
            }
        }
    }

  int countHash = 0;
  std::map<std::size_t,int> uniqueHash;
  for (int a=0; a<label_set.size(); a++)
    {
      std::size_t newHash = hash_value(label_set[a]);
      std::cout << label_set[a].Str() << " " << newHash << std::endl;
      countHash++;
      uniqueHash[newHash]++;

    }

  int collisionHash = 0;
  for (std::map<std::size_t,int>::const_iterator it = uniqueHash.begin();
       it != uniqueHash.end(); ++it)
    {
      //      if (it->second > 1)
      std::cout << it->first << "\t" << it->second << std::endl;
      //        collisionHash++;
    }

  std::cout << "Total number of hashes made: " << std::to_string(countHash) << std::endl;
  //std::cout << "Number of hash collision: " << std::to_string(collisionHash) << std::endl;
}


void caching_test()
// Test use of caching wrapper for U coefficients
{

  // generate label set for testing
  u3::SU3 x1(4,0);
  u3::SU3 x2(3,2);
  u3::SU3 x3(2,3);
  // u3::SU3 x1(0,0);
  // u3::SU3 x2(0,0);
  // u3::SU3 x3(0,0);
  std::vector<u3::UCoefLabels> label_set;
  MultiplicityTagged<u3::SU3>::vector x12_values=KroneckerProduct(x1,x2);
  MultiplicityTagged<u3::SU3>::vector x23_values=KroneckerProduct(x2,x3);
  for (auto it12 = x12_values.begin(); it12 != x12_values.end(); ++it12)
    for (auto it23 = x23_values.begin(); it23 != x23_values.end(); ++it23)
      {
        u3::SU3 x12 = it12->irrep;
        u3::SU3 x23 = it23->irrep;
        MultiplicityTagged<u3::SU3>::vector x_values=KroneckerProduct(x12,x3);
        for(auto it = x_values.begin(); it != x_values.end(); ++it)
          {
            u3::SU3 x = it->irrep;
            if(u3::OuterMultiplicity(x1,x23,x)>0)
              {
                u3::UCoefLabels labels(x1,x2,x,x3,x12,x23);

                label_set.push_back(labels);
              }
          }
      }

  // cache coefficients
  std::cout << "Caching coefficients" << std::endl;
  u3::UCoefCache u_coef_cache;
  for (const u3::UCoefLabels& labels : label_set)
    {
      //if (u_coef_cache.count(labels)>0)
      //  std::cout << "  duplicate " << labels.Str() << std::endl;
      if ((u_coef_cache.size()%100)==0)
        std::cout << "  cache size " << u_coef_cache.size() << "..." << std::endl;
      u_coef_cache[labels] = u3::UCoefBlock(labels);
    }
  std::cout << "  cached " << u_coef_cache.size() << std::endl;

  // retrieve labels and compare with on-the-fly values
  std::cout << "Checking cached values" << std::endl;
  for (const u3::UCoefLabels& labels : label_set)
    {
      // retrieve labels
      u3::SU3 x1, x2, x, x3, x12, x23;
      std::tie(x1,x2,x,x3,x12,x23) = labels.Key();

      // retrieve coefficient block
      const u3::UCoefBlock& block = u_coef_cache[labels];

      // retrieve multiplicities
      int r12_max, r12_3_max, r23_max, r1_23_max;
      std::tie(r12_max, r12_3_max, r23_max, r1_23_max) = block.Key();

      // loop over multiplicity indices
      // std::cout << "   " << labels.Str() << std::endl;
      for (int r12 = 1; r12 <= r12_max; ++r12)
        for (int r12_3 = 1; r12_3 <= r12_3_max; ++r12_3)
          for (int r23 = 1; r23 <= r23_max; ++r23)
            for (int r1_23 = 1; r1_23 <= r1_23_max; ++r1_23)
              {
                double coef_direct = u3::U(x1,x2,x,x3,x12,r12,r12_3,x23,r23,r1_23);
                double coef_cached = u3::UCached(u_coef_cache,x1,x2,x,x3,x12,r12,r12_3,x23,r23,r1_23);
                bool compare_ok = (coef_direct == coef_cached);
                if (!compare_ok)
                  std::cout << " " << coef_direct << " " << coef_cached << " " << compare_ok << std::endl;
              }
    }
  std::cout << "Done." << std::endl;

}


int main(int argc, char **argv)
{
  // initialize su3lib
  u3::U3CoefInit();

  // basic tests
  // basic_test();

  // test of block storage
  // block_test();
  
  // test of hash key
  // iteration_test();

  // test cache storage and retrieval
  caching_test();  

}

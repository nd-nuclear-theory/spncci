/****************************************************************
  u3_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created.
  3/8/16 (aem,mac): Add tests for U3ST.

****************************************************************/

#include "sp3rlib/u3.h"

#include <vector>
#include <algorithm>

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // construction and string conversion
  ////////////////////////////////////////////////////////////////

  std::cout << u3::SU3(1,2).Str() << std::endl;
  std::cout 
    << u3::U3(0,0,0).Str()
    << u3::U3(HalfInt(1,2),HalfInt(1,2),HalfInt(1,2)).Str()
    << u3::U3(HalfInt(5,2),u3::SU3(1,0)).Str()
    << std::endl;

  ////////////////////////////////////////////////////////////////
  // comparison operators
  ////////////////////////////////////////////////////////////////

  // SU(3)
  std::vector<u3::SU3> lm_vector;
  lm_vector.push_back(u3::SU3(1,2));
  lm_vector.push_back(u3::SU3(4,3));
  lm_vector.push_back(u3::SU3(1,1));
  lm_vector.push_back(u3::SU3(4,5));
  sort(lm_vector.begin(),lm_vector.end());
  for (int i=0; i<lm_vector.size(); ++i)
    {
      std::cout << lm_vector[i].Str();
    }
  std::cout << std::endl;

  // U(3)
  std::vector<u3::U3> w_vector;
  w_vector.push_back(u3::U3(1,1,1));
  w_vector.push_back(u3::U3(3,1,0));
  w_vector.push_back(u3::U3(3,2,1));
  w_vector.push_back(u3::U3(2,1,1));
  for (int i=0; i<w_vector.size(); ++i)
    {
      std::cout << w_vector[i].Str();
    }
  std::cout << std::endl;
  sort(w_vector.begin(),w_vector.end());
  for (int i=0; i<w_vector.size(); ++i)
    {
      std::cout << w_vector[i].Str();
    }
  std::cout << std::endl;

  // U(3)xSU(2)
  std::vector<u3::U3S> wS_vector;
  wS_vector.push_back(u3::U3S(u3::U3(1,1,1),3));
  wS_vector.push_back(u3::U3S(u3::U3(1,1,1),1));
  wS_vector.push_back(u3::U3S(u3::U3(3,1,0),HalfInt(3,2)));
  wS_vector.push_back(u3::U3S(u3::U3(3,1,0),HalfInt(1,2)));
  for (int i=0; i<wS_vector.size(); ++i)
    {
      std::cout << wS_vector[i].Str();
    }
  std::cout << std::endl;
  sort(wS_vector.begin(),wS_vector.end());
  for (int i=0; i<wS_vector.size(); ++i)
    {
      std::cout << wS_vector[i].Str();
    }
  std::cout << std::endl;

  // U(3)xSU(2)xSU(2)
  std::vector<u3::U3ST> wST_vector;
  wST_vector.push_back(u3::U3ST(u3::U3(1,1,1),3,2));
  wST_vector.push_back(u3::U3ST(u3::U3(1,1,1),1,1));
  wST_vector.push_back(u3::U3ST(u3::U3(3,1,0),HalfInt(3,2),HalfInt(1,2)));
  wST_vector.push_back(u3::U3ST(u3::U3(3,1,0),HalfInt(1,2),1));
  for (int i=0; i<wST_vector.size(); ++i)
    {
      std::cout << wST_vector[i].Str()<<"  ";
    }
  std::cout << std::endl;
  sort(wST_vector.begin(),wST_vector.end());
  for (int i=0; i<wST_vector.size(); ++i)
    {
      std::cout << wST_vector[i].Str()<<"  ";
    }
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // dimension, conjugation, and validation tests
  ////////////////////////////////////////////////////////////////

  for (int i=0; i<lm_vector.size(); ++i)
    {
      std::cout 
	<< lm_vector[i].Str() 
	<< " " << dim(lm_vector[i])
	<< " " << Conjugate(lm_vector[i]).Str() 
	<< " " << ConjugationGrade(lm_vector[i])
	<< std::endl;
    }
  std::cout << std::endl;

  for (int i=0; i<w_vector.size(); ++i)
    {
      std::cout 
	<< w_vector[i].Str() 
	<< " " << dim(w_vector[i])
	<< " " << Conjugate(w_vector[i]).Str() 
	<< " " << ConjugationGrade(w_vector[i])
	<< " " << w_vector[i].Valid()
	<< std::endl;
    }
  std::cout << std::endl;
  std::cout 
    << u3::U3(1,3,2).Str() 
    << " " << u3::U3(1,3,2).Valid() 
    << std::endl;

  for (int i=0; i<wS_vector.size(); ++i)
    {
      std::cout 
	<< wS_vector[i].Str() 
	<< " " << dim(wS_vector[i])
	<< std::endl;
    }
  std::cout << std::endl;

}

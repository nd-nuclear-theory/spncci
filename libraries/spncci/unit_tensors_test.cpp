/****************************************************************
  unit_tensor_test.cpp

  Unit tensor algorithms
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/18/16 (aem,mac): Created.
****************************************************************/
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include <map>
#include <sstream>
#include <vector>

#include "am/am.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "spncci/sp_basis.h"
#include "spncci/unit_tensors.h"

spncci::LGIVectorType lgi_vector;
std::map< u3::U3,std::map<u3::U3,Eigen::MatrixXd> > K_matrix_map;

int main(int argc, char **argv)
{
	u3::U3CoefInit();
	// For generating the lgi_vector, using Li-6 as example;
 	HalfInt Nsigma_0 = HalfInt(11,1);
 	int N1b=2;
	std::string filename = "libraries/spncci/lgi-3-3-2-fql-mini-mini.dat";

	
	spncci::GenerateLGIVector(lgi_vector,filename,Nsigma_0);

  std::vector< std::pair<int,int> > lgi_pair_vector=spncci::LGIPairGenerator(lgi_vector);

  int Nmax=4;

  std::map< int, std::vector<u3::UnitTensor> > unit_sym_map;
  u3::UnitSymmetryGenerator(Nmax,unit_sym_map);

	std:: map< std::pair<int,int>, std::map< u3::UnitTensorRME,Eigen::MatrixXd> > lgi_unit_tensor_rme_map;

  for (int i=0; i<lgi_pair_vector.size(); i++)
  	{
  		std::pair<int,int> lgi_pair=lgi_pair_vector[i];
  		u3::U3 sigmap=lgi_vector[lgi_pair.first].sigma;
  		// operator boson number between to lgi's
  		u3::U3 sigma=lgi_vector[lgi_pair.second].sigma;
  		int N0=int(sigmap.N()-sigma.N());
  		
  		std::map <u3::UnitTensorRME, Eigen::MatrixXd> temp_unit_map;
  		for (int j=0; j<unit_sym_map[N0].size(); j++)
  			{
  				Eigen::MatrixXd temp_matrix(1,1);
  				temp_matrix(0,0)=1;
					
					u3::UnitTensor unit_tensor=unit_sym_map[N0][j];
					int rho0_max=u3::OuterMultiplicity(sigma.SU3(),unit_tensor.omega0.SU3(), sigmap.SU3());
					for (int rho0=1; rho0<=rho0_max; rho0++)
						{
  						if (unit_tensor.rp<=(N1b+(sigmap.N()-Nsigma_0)) && unit_tensor.r<=(N1b+(sigma.N()-Nsigma_0)))
  						{
	  						//std::cout<<unit_tensor.Str()<<std::endl;
	  						temp_unit_map[u3::UnitTensorRME(sigmap,sigma,unit_tensor,rho0)]=temp_matrix;
  							
  						}
  					}
  			}

  		lgi_unit_tensor_rme_map[lgi_pair]=temp_unit_map;

  		// for (auto it=lgi_unit_tensor_rme_map.begin(); it !=lgi_unit_tensor_rme_map.end(); ++it)
  		// 	for (auto i=lgi_unit_tensor_rme_map[it->first].begin(); i !=lgi_unit_tensor_rme_map[it->first].end(); i++)
		  // 		{
		  // 			std::cout <<i->first.tensor.Str()<<"  "<<i->second<<std::endl;
		  // 		}
			//u3::TestFunction(Nmax, lgi_pair, unit_sym_map, lgi_unit_tensor_rme_map);
  		u3::UnitTensorMatrixGenerator(N1b, Nmax, lgi_pair, unit_sym_map,lgi_unit_tensor_rme_map[lgi_pair] );

  	}





}// end main 

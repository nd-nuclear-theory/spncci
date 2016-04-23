#include <cstdio>
#include "spncci/unit_tensors.h"
#include <map>

spncci::LGIVectorType lgi_vector;
std::map< u3::U3,std::map<u3::U3,Eigen::MatrixXd> > K_matrix_map;

void iteration_test()
{
  // tensor set
  std::vector< :: > tensor_set;

  // create testing tensor set and tensorU3Sector set here
  // GenerateUnitTensor(Nmax, map of tensor labels)

  int countHash = 0;
  std::map<std::size_t,int> uniqueHash;
  for (int a=0; a<tensor_set.size(); a++)
    {
      std::size_t newHash = hash_value(tensor_set[a]);
      std::cout << tensor_set[a].Str() << " " << newHash << std::endl;
      countHash++;
      uniqueHash[newHash]++;

    }

  int collisionHash = 0;
  for (std::map<std::size_t,int>::const_iterator it = uniqueHash.begin();
        it != uniqueHash.end(); ++it)
    {
      std::cout << it->first << "\t" << it->second << std::endl;
    }

std::cout << "Total number of hashes made: " << std::to_string(countHash) << std::endl;
//std::cout << "Number of hash collision: " << std::to_string(collisionHash) << std::endl;
}

int Nmax;
int main(int argc, char **argv)
{
  if(argc>1)
      Nmax=std::stoi(argv[1]); 
  else
    Nmax=2;


	u3::U3CoefInit();
	// For generating the lgi_vector, using Li-6 as example;
 	HalfInt Nsigma_0 = HalfInt(11,1);
 	int N1b=2;
  // input file containing LGI's
	std::string filename = "libraries/spncci/lgi-3-3-2-fql-mini-mini.dat";

	// Generate vector of LGI's from input file 
	spncci::GenerateLGIVector(lgi_vector,filename,Nsigma_0);

  spncci::SigmaIrrepMapType sigma_irrep_map;
  spncci::NmaxTruncator truncator(Nsigma_0,Nmax);
  spncci::GenerateSp3RIrreps(lgi_vector,sigma_irrep_map,truncator);

  // Generate list of LGI's for which two-body operators will have non-zero matrix elements 
  std::vector< std::pair<int,int> > lgi_pair_vector=spncci::GenerateLGIPairs(lgi_vector);

  // generate map that stores unit tensor labels keyed by N0
  std::map< int, std::vector<spncci::UnitTensor> > unit_sym_map;
  spncci::GenerateUnitTensors(Nmax,unit_sym_map);

  //initializing map that will store map containing unit tensors.  Outer map is keyed LGI pair. 
  // inner map is keyed by unit tensor matrix element labels of type UnitTensorRME
  // LGI pair -> UnitTensorRME -> Matrix of reduced unit tensor matrix elements for v'v subsector
  std:: map< 
            std::pair<int,int>, 
            std::map<
              std::pair<int,int>,
              std::map< spncci::UnitTensorU3Sector,Eigen::MatrixXd> 
              >
          > lgi_unit_tensor_rme_map;


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Filling out lgi_unit_tensor_rme_map 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  for (int i=0; i<lgi_pair_vector.size(); i++)
  	{
      std::pair<int,int> lgi_pair=lgi_pair_vector[i];
      spncci::LGI lgip=lgi_vector[lgi_pair.first];
      spncci::LGI lgi=lgi_vector[lgi_pair.second];
  		u3::U3 sigmap=lgip.sigma;
  		
      // operator boson number between to lgi's
  		u3::U3 sigma=lgi.sigma;
  		int N0=int(sigmap.N()-sigma.N());
  		
  		std::map <spncci::UnitTensorU3Sector, Eigen::MatrixXd> temp_unit_map;
      int rp,r;
      HalfInt S0;
      u3::U3 omega0;
      //////////////////////////////////////////////////////////////////////////////////////////////
      // Initializing the unit_tensor_rme_map with LGI rm's 
  		//////////////////////////////////////////////////////////////////////////////////////////////
      std::pair<int,int> N0_pair(lgip.Nex,lgi.Nex);
      for (int j=0; j<unit_sym_map[N0].size(); j++)
  			{
  				Eigen::MatrixXd temp_matrix(1,1);
  				temp_matrix(0,0)=1;
					
					spncci::UnitTensor unit_tensor=unit_sym_map[N0][j];
          std::tie (omega0, S0, std::ignore, rp, std::ignore, std::ignore, r, std::ignore,std::ignore)=unit_tensor.Key();
					int rho0_max=u3::OuterMultiplicity(sigma.SU3(),omega0.SU3(), sigmap.SU3());
					for (int rho0=1; rho0<=rho0_max; rho0++)
						{
  						if (
                rp<=(N1b+lgip.Nex) 
                && r<=(N1b+lgi.Nex)
                && abs(lgi.S+S0)>=lgip.S
                )
    						{
  	  						//std::cout<<unit_tensor.Str()<<std::endl;
  	  						temp_unit_map[spncci::UnitTensorU3Sector(sigmap,sigma,unit_tensor,rho0)]=temp_matrix;	
    						}
  					}
  			}
  		lgi_unit_tensor_rme_map[lgi_pair][N0_pair]=temp_unit_map;
      //////////////////////////////////////////////////////////////////////////////////////////////
      // Generating the rme's of the unit tensor for each LGI
      spncci::GenerateUnitTensorMatrix(N1b, Nmax, lgi_pair, unit_sym_map,lgi_unit_tensor_rme_map[lgi_pair] );

  		// for (auto it=lgi_unit_tensor_rme_map.begin(); it !=lgi_unit_tensor_rme_map.end(); ++it)
  		// 	for (auto i=lgi_unit_tensor_rme_map[it->first].begin(); i !=lgi_unit_tensor_rme_map[it->first].end(); i++)
		  		
    //       {
		  // 			std::cout <<i->first.tensor.Str()<<"  "<<i->second<<std::endl;
		  // 		}
			//u3::TestFunction(Nmax, lgi_pair, unit_sym_map, lgi_unit_tensor_rme_map);

  	}
  // test the hashing over an iteration of tensors and u3sectors
  iteration_tset();
}// end main 

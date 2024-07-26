#include <iostream>
#include <vector>
#include <array>
#include <fstream>

#include <su3.h>

std::vector<std::array<int,2>> KroneckerProduct(int lambda1, int mu1, int lambda2, int mu2){
  std::vector<std::array<int,2>> product;
  for(int lambda=0; lambda<=lambda1+lambda2+std::min(mu2,lambda1+mu1); lambda++){
    for(int mu=0; mu<=mu1+mu2+std::min(lambda1,lambda2); mu++){
      int rhomax=su3::mult(lambda1,mu1,lambda2,mu2,lambda,mu);
      if(rhomax!=0){
        product.push_back({lambda,mu});
      }
    }
  }
  return product;
}

int main(int argc, char **argv){
  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " Nmax N1v" << std::endl;
    return EXIT_FAILURE;
  }

  int Nmax = std::stoi(argv[1]);
  int N1v = std::stoi(argv[2]);
  
  std::ofstream file("obd.dat");
  if (!file) {
    std::cout << "Could not open file obd.dat" << std::endl;
    return EXIT_FAILURE;
  }

  int eta_max=Nmax+N1v;

  for(int N0=-Nmax; N0<=Nmax; N0++){       
    for(int etap=0; etap<=eta_max; etap++){ 
      int eta=etap-N0;
      if((eta<0)||(eta>eta_max))continue;
      std::vector<std::array<int,2>> x0_set=KroneckerProduct(etap,0,0,eta);
      for (int SS0=0; SS0<=2; SS0=SS0+2){
        for(int w=0; w<x0_set.size(); w++){
          file<<etap<<" "<<eta<<" "<<x0_set[w][0]<<" "<<x0_set[w][1]<<" "<<SS0<<std::endl;
        }
      }
    }
  }

  return EXIT_SUCCESS;
}

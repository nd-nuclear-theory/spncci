/****************************************************************
  generate_raising_polynomials.cpp
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/16/16 (aem): Created.
****************************************************************/
#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_operator.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  if(argc<6)
    {
      std::cout<<"Syntax: Protons Neutrons Nmax Nstep N1B"<<std::endl;
      std::exit(1);
    }
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  // will be either 1 or 2; 
  int Nstep=std::stoi(argv[4]);
  int N1B=std::stoi(argv[5]);
  assert((Nstep==2)||(Nstep==1));
  int Nmin=Nmax%Nstep;
  int A=N+Z;
  int T0=0;
  int J0=0;
  bool un_u3_restrict=false;
  if( (N==0) || (Z==0) )
    un_u3_restrict=true;
  // Set up unit tensor model space space
  std::string model_space=fmt::format("model_space_{}_{}_Nmax{:02d}.dat",Z,N,Nmax);
  std::ofstream model_stream(model_space);
  model_stream<<Z<<"  "<<N<<"  "<<-1<<std::endl;
  for(int Nex=Nmin; Nex<=Nmax; Nex+=2)
    model_stream<<Nex<<"  "<<-1<<std::endl;
  model_stream.close();

  //begin control file
  // first give specifications for unit tensors, then Brel and Nrel
  std::ofstream control_stream("raising_polynomials.dat");

  std::vector<u3::U3>raising_polynomials=sp3r::RaisingPolynomialLabels(Nmax);

  for(int i=0; i<raising_polynomials.size(); ++i)
  {
    u3::U3 n0(raising_polynomials[i]);
    //Generate Nintr operator up to Nmax cutoff
    std::string Prel_file=fmt::format("Prel_{:02d}_Nmax{:02d}_{:06d}.recoupler",N+Z,Nmax,i);
    u3shell::RelativeUnitTensorCoefficientsU3ST Prel_operator;
    u3shell::RaisingPolynomialRelativeUnitTensorExpansion(n0,0, Nmax, Prel_operator, N+Z);
    lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, Prel_operator, Prel_file, un_u3_restrict);
  }
  for(int i=0; i<raising_polynomials.size(); ++i)
    control_stream<<fmt::format("Prel_{:02d}_Nmax{:02d}_{:06d}.recoupler",N+Z,Nmax,i)<<std::endl;

}
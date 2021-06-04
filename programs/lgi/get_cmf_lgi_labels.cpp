#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "lgi/lgi.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "sp3rlib/sp3r.h"
#include "spncci/spncci_basis.h"
#include "u3shell/u3spn_scheme.h"


typedef std::tuple<int,int,int> SpinTuple;
typedef std::tuple<SpinTuple,SpinTuple,SpinTuple,SpinTuple> CompoundSpinTuple;

// Kronecker product of U(3) x Sp x Sn x S and U(3)
MultiplicityTagged<u3shell::U3SPN>::vector KroneckerProduct(const u3shell::U3SPN& irrep, const u3::U3& omega)
{
  MultiplicityTagged<u3::U3>::vector u3_product = u3::KroneckerProduct(irrep.U3(), omega);
  MultiplicityTagged<u3shell::U3SPN>::vector product;
  product.reserve(u3_product.size());
  for (MultiplicityTagged<u3::U3> omega_tagged : u3_product)
  {
    MultiplicityTagged<u3shell::U3SPN> irrep_tagged(
        u3shell::U3SPN(omega_tagged.irrep, irrep.Sp(), irrep.Sn(), irrep.S()), omega_tagged.tag);
    product.push_back(irrep_tagged);
  }
  return product;
}

// inner multiplicity of L
// int inner_multiplicity(const u3::SU3& x, const int& L)
//{
//  int kappa_max=std::max(0,(x.lambda()+x.mu()+2-L)/2)-std::max(0,(x.lambda()+1-L)/2)-std::max(0,(x.mu()+1-L)/2);
//  return kappa_max;
//}

struct Dimensions
{
  int total, cmf, LGI;
  Dimensions()
      : total(0), cmf(0), LGI(0) {}  // default constructor (if omitted, compiler complains)
  Dimensions(int total_, int cmf_, int LGI_) : total(total_), cmf(cmf_), LGI(LGI_) {}
};

/* struct LS
{
  int L;
  HalfInt S;
  LS() : L(0), S(0) {} // default constructor
  LS(int L_, HalfInt S_) : L(L_), S(S_) {}
}; */

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string input_filename;
  std::string output_filename;
  // mode
  lgi::NuclideType nuclide;
  int Nmax; 
};

RunParameters ProcessArguments(int argc, char **argv)
{
  RunParameters run_parameters;
  // usage message
  if (argc-1 < 5)
    {
      std::cout << "Syntax: get_cmf_lgi_labels Z N Nmax, input_filename output_filename" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  // nuclide
  int Z, N, Nmax;
  {
    std::istringstream parameter_stream(argv[1]);
    parameter_stream >> Z;
    if (!parameter_stream)
    {
      std::cerr << "Expecting numeric value for Z" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  {
    std::istringstream parameter_stream(argv[2]);
    parameter_stream >> N;
    if (!parameter_stream)
    {
      std::cerr << "Expecting numeric value for Z" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  {
    std::istringstream parameter_stream(argv[3]);
    parameter_stream >> Nmax;
    if (!parameter_stream)
    {
      std::cerr << "Expecting numeric value for Nmax" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  }
  run_parameters.nuclide = {Z,N};
  run_parameters.Nmax=Nmax;
  // input filename
  run_parameters.input_filename = argv[4];

  // output filename
  run_parameters.output_filename = argv[5];

  return run_parameters;
}


void read_lsu3shell_basis_dimensions(
  const std::string& input_filename, 
  const int N0,  const int A, const int Nmax,
  std::map<u3shell::U3SPN, Dimensions>& dimensions_by_irrep
  )
  {

    // input from a file
    std::ifstream file(input_filename);
    //  int N_ex_max=0;
    std::string line;
    int line_count = 0;
    int input_Nmax = 0;
    while (mcutils::GetLine(file, line, line_count))
      {
        std::istringstream line_stream(line);
        int N_ex, twice_Sp, twice_Sn, twice_S, lambda, mu, dim_tot;
        line_stream >> N_ex >> lambda >> mu >> twice_Sp >> twice_Sn >> twice_S >> dim_tot;
        mcutils::ParsingCheck(line_stream,line_count,line);

        input_Nmax = std::max(input_Nmax, N_ex);
        
        if (N_ex>Nmax)
          continue;

        HalfInt N(2 * (N_ex + N0) + 3 * (A-1), 2), Sp(twice_Sp, 2), Sn(twice_Sn, 2), S(twice_S, 2);
        u3::SU3 x(lambda, mu);
        u3::U3 omega(N, x);
        u3shell::U3SPN irrep(omega, Sp, Sn, S);
        Dimensions irrep_dimensions(dim_tot, dim_tot, dim_tot);
        dimensions_by_irrep[irrep] = irrep_dimensions;
      }
    if (input_Nmax<Nmax)
      std::cout<<fmt::format("Warning:\n Nmax of input: {}\n Nmax requested: {}",input_Nmax,Nmax)<<std::endl;
    file.close();


  // Sanity check
  for (std::map<u3shell::U3SPN, Dimensions>::iterator it = dimensions_by_irrep.begin();
       it != dimensions_by_irrep.end(); ++it)
    {
      bool valid_dimensions=true;
      const Dimensions& dimensions=it->second; 
      valid_dimensions &= (dimensions.total > 0);
      valid_dimensions &= (dimensions.cmf > 0);
      valid_dimensions &= (dimensions.LGI > 0);

      if(not valid_dimensions)
        std::cout<<fmt::format("Sanity check failed for {}\n Dimensions are: {}  {}  {}",it->first.Str(),dimensions.total,dimensions.cmf,dimensions.LGI)<<std::endl;
      assert(valid_dimensions);
    }
  }


void generate_cmf_lgi(int Nmax,HalfInt Nsigma0,std::map<u3shell::U3SPN, Dimensions>& dimensions_by_irrep,lgi::MultiplicityTaggedLGIVector& lgi_vector)
  {
    // Iterate through the lsu3shell basis and remove CM contaminated states 
    for(int Nex=0; Nex<=Nmax; ++Nex)
      for(const auto& irrep_dimensions : dimensions_by_irrep)
        {         
          const u3shell::U3SPN& irrep=irrep_dimensions.first;
          const auto& dimensions=irrep_dimensions.second;

          int Ncm_max=Nmax-Nex;
          if( (irrep.N()-Nsigma0) == Nex)
            for(int Ncm=1; Ncm<=Ncm_max; Ncm++)
              {
                u3::U3 wcm(Ncm,u3::SU3(Ncm,0));
                MultiplicityTagged<u3shell::U3SPN>::vector cm_irreps=KroneckerProduct(irrep, wcm);
                for(const auto& cm_irrep_tagged : cm_irreps)
                  {
                    assert(cm_irrep_tagged.tag==1);
                    dimensions_by_irrep[cm_irrep_tagged.irrep].cmf -= dimensions.cmf * cm_irrep_tagged.tag;
                  }
              }
        }

    //Copy cmf dimensions to lgi dimensions
    for(auto& irrep_dimension: dimensions_by_irrep)
      irrep_dimension.second.LGI=irrep_dimension.second.cmf;

    //Iterate through basis and identify LGI dimension by substracting 
    //U(3) irreps obtained by laddering from lower grade LGI.
    for(const auto& [lgi,dimensions] : dimensions_by_irrep)
      { 
        HalfInt Sp(lgi.Sp()),Sn(lgi.Sn()), S(lgi.S());
        int Nn_max=Nmax-int(lgi.N()-Nsigma0);
        std::vector<u3::U3> raising_polynomial_labels = sp3r::RaisingPolynomialLabels(Nn_max);

        for(const u3::U3& n : raising_polynomial_labels)
          {
            if (n.N()==0)
              continue;

            MultiplicityTagged<u3::U3>::vector omegas_tagged = u3::KroneckerProduct(lgi.U3(), n);
            for(const auto& [omega,rho_max] : omegas_tagged)
              {
                u3shell::U3SPN omegaSpSnS(omega,Sp,Sn,S);
                dimensions_by_irrep[omegaSpSnS].LGI -= rho_max*dimensions.LGI;
              }
          }
      }

    //Create LGI vector used in SpNCCI basis construction
    for(const auto& [lgi_u3spn,lgi_dims] : dimensions_by_irrep)
      {
        int Nsex=int(lgi_u3spn.N()-Nsigma0);
        lgi::LGI lgi(lgi_u3spn,Nsex);
        lgi_vector.emplace_back(lgi,lgi_dims.LGI);
      }

  }

int main(int argc, char **argv)
{
  auto run_parameters = ProcessArguments(argc, argv);
  //  int Nmax=0, N_sigma_max=0; // N=N_ex+N0+(3/2)*A (lab frame considered)
  HalfInt Nsigma0 = lgi::Nsigma0ForNuclide(run_parameters.nuclide,true);
  int A = run_parameters.nuclide[0] + run_parameters.nuclide[1];
  int N0 = int(Nsigma0 - HalfInt(3, 2) * (A - 1));// CMF N0
  int Nmax=run_parameters.Nmax;
  // std::cerr << A << " " << N0 << " " << Nsigma0 << std::endl;

  std::map<u3shell::U3SPN, Dimensions> dimensions_by_irrep;
  read_lsu3shell_basis_dimensions(run_parameters.input_filename, N0,A, run_parameters.Nmax, dimensions_by_irrep);
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  generate_cmf_lgi(Nmax,Nsigma0,dimensions_by_irrep,lgi_vector);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Get big block dimensions 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // build SpNCCI irrep branchings
  spncci::NmaxTruncator truncator(Nsigma0,Nmax);
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;
  
  bool restrict_sp3r_to_u3_branching=false;
  spncci::GenerateSpNCCISpace(lgi_vector,truncator,spncci_space,sigma_irrep_map,restrict_sp3r_to_u3_branching);

  // For a given sigma, get list of all possible gamma,Sp,Sn,S
  std::map<u3::U3,std::vector<std::tuple<int,HalfInt,HalfInt,HalfInt>>> factored_basis;
  for(const auto& lgi_tagged : lgi_vector)
    {
      const lgi::LGI& lgi = lgi_tagged.irrep;
      int gamma_max=lgi_tagged.tag;
      factored_basis[lgi.U3()].emplace_back(gamma_max,lgi.Sp(),lgi.Sn(),lgi.S());
    }

  // for(const auto& [sigma,upstreams] : factored_basis)
  //     std::cout<<sigma.Str()<<"  "<<upstreams.size()<<std::endl;        

  // Given S1 and S2, what are possible S0 (or T1,T2,T0)
  std::vector<SpinTuple> spin_tuples;
  for(int S1=0; S1<=1; S1++)
    for(int S2=0; S2<=1; S2++)
      for(int S0=abs(S1-S2); S0<=(S1+S2); ++S0)
        {
          spin_tuples.emplace_back(S1,S2,S0);
        }

  int Nv=1;
  for(const auto& [sigma1,upstreams1] : factored_basis)
    for(const auto& [sigma2,upstreams2] : factored_basis)
      {

        // Get max values for N1 and N2 of relative unit tensors
        int N0=int(sigma1.N()-sigma2.N());
        int N1_max=2*Nv+int(sigma1.N()-Nsigma0);
        int N2_max=2*Nv+int(sigma2.N()-Nsigma0); 

        std::map<std::tuple<int,int,u3::SU3>,int> lgi_spatial_labels;
        for (int N2=0; N2<=N2_max; ++N2)  
          {
            int N1=N2+N0; 
            if ((N1>N1_max) or (N1<0))
              continue;

            u3::SU3 x1(N1,0), x2(0,N2);
            MultiplicityTagged<u3::SU3>::vector product_irreps = u3::KroneckerProduct(x1, x2);
            for(int i=0; i<product_irreps.size(); ++i)
              { 
                u3::SU3 x0(product_irreps[i].irrep);
                int rho_max=u3::OuterMultiplicity(sigma2.SU3(), x0,sigma1.SU3());
                if (rho_max>0)
                  {
                    std::tuple<int,int,u3::SU3> labels(N1,N2,x0); 
                    lgi_spatial_labels[labels]=rho_max;
                  }                
              }
          }

        if(lgi_spatial_labels.size()==0)
          continue;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Generating unit tensors for each Nex
        int Nsex1=int(sigma1.N()-Nsigma0);
        int Nsex2=int(sigma2.N()-Nsigma0);

        //Zero initialize with lgi spatial labels 
        std::map<std::pair<int,int>,std::set<std::tuple<int,int,u3::SU3>>> spatial_unit_tensor_labels_by_Nex;
        for(auto& [label,num] : lgi_spatial_labels)
          spatial_unit_tensor_labels_by_Nex[std::pair<int,int>(0,0)].insert(label);

        for(int Nn1=0; Nn1<=(Nmax-Nsex1); Nn1+=2)
          for(int Nn2=0; Nn2<=(Nmax-Nsex2); Nn2+=2)
            {
              std::pair<int,int> Nnpair(Nn1,Nn2);

              // If Nn=0 for both irreps, initialize labels dictionary with lgi spatial labels      
              if(Nn1==0 and Nn2==0)
                for(auto& [label,num] : lgi_spatial_labels)
                  spatial_unit_tensor_labels_by_Nex[Nnpair].insert(label);
              else
                {                  
                  
                  // If both Nn1 and Nn2 are greater than zero, then source unit tensors
                  // also valid for target Nnpair.  (Term 1 of recurrence equation).
                  if (Nn1>0 and Nn2>0)
                    {
                      std::pair<int,int> Nnpair_source(Nn1-2,Nn2-2);
                      const auto& source_labels = spatial_unit_tensor_labels_by_Nex[Nnpair_source];      
                      for(auto& label : source_labels)
                        spatial_unit_tensor_labels_by_Nex[Nnpair].insert(label);
                    }

                  //Checking for all other cases
                  std::vector<std::tuple<int,int,u3::SU3>> test_N1N2x0_values;

                  if(Nn2>0)
                    {                  
                      std::pair<int,int> Nnpair_source(Nn1,Nn2-2);
                      const auto& source_labels = spatial_unit_tensor_labels_by_Nex[Nnpair_source];      
                      for(auto& label : source_labels)
                        {
                          const auto& [N1,N2,x0] = label;
                          if (N1-2>=0)
                            test_N1N2x0_values.emplace_back(N1-2,N2,x0);    
                        
                          if (N2+2<=(N2_max+Nn2))
                            test_N1N2x0_values.emplace_back(N1,N2+2,x0);
                        }
                    }
                  //
                  if(Nn2==0 and Nn1!=0)
                    {
                      std::pair<int,int> Nnpair_source(Nn1-2,Nn2);
                      const auto& source_labels = spatial_unit_tensor_labels_by_Nex[Nnpair_source];      
                      for(auto& label : source_labels)
                        {
                          const auto& [N1,N2,x0] = label;
                          if (N2-2>=0)
                            test_N1N2x0_values.emplace_back(N1,N2-2,x0);
                            
                        
                          if (N1+2<=(N1_max+Nn1))
                            test_N1N2x0_values.emplace_back(N1+2,N2,x0);
                        }
                    }
                  
                  for(auto& [N1_new,N2_new,x0] : test_N1N2x0_values)
                    {
                      MultiplicityTagged<u3::SU3>::vector product_irreps 
                        = u3::KroneckerProduct(u3::SU3(N1_new,0), u3::SU3(0,N2_new));

                      for(const auto& [x0_new,dummy] : product_irreps)
                          if(u3::OuterMultiplicity(x0_new,u3::SU3(2,0),x0))
                            {
                              std::tuple<int,int,u3::SU3> new_label(N1_new,N2_new,x0_new);
                              spatial_unit_tensor_labels_by_Nex[Nnpair].insert(new_label);                            
                            }   
                    }  
                }
            }//end Nn2
          //end Nn1
        // for(const auto& [key,value] : spatial_unit_tensor_labels_by_Nex)
        //   {
        //     const auto& [Nn1,Nn2]=key;
        //     std::cout<<Nn1<<"  "<<Nn2<<"  "<<value.size()<<std::endl;
        //   }                  

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        const auto& irrep1=sigma_irrep_map[sigma1];
        const auto& irrep2=sigma_irrep_map[sigma2];
        // std::vector<int> partition1=PartitionIrrepByNn(irrep1,Nmax);
        // std::vector<int> partition2=PartitionIrrepByNn(irrep2,Nmax);
        std::map<std::tuple<u3::U3,u3::U3,u3::SU3,int,int>,int> spatial_labels;
        for(int index1=0; index1<irrep1.size(); ++index1)
          for(int index2=0; index2<irrep2.size(); ++index2)
            {
              const sp3r::U3Subspace& subspace1 = irrep1.GetSubspace(index1);
              const sp3r::U3Subspace& subspace2 = irrep2.GetSubspace(index2);
              
              const u3::U3& omega1=subspace1.U3();
              const u3::U3& omega2=subspace2.U3();
              int upsilon_max1=subspace1.upsilon_max();
              int upsilon_max2=subspace2.upsilon_max();

              int Nn1=int(omega1.N()-sigma1.N());
              int Nn2=int(omega2.N()-sigma2.N());
              
              std::pair<int,int>Nnpair(Nn1,Nn2);
              const auto& spatial_unit_tensor_labels=spatial_unit_tensor_labels_by_Nex[Nnpair];
              for(const auto& [N1,N2,x0] : spatial_unit_tensor_labels)
                {
                  int rho_max=u3::OuterMultiplicity(omega2.SU3(),x0,omega1.SU3());
                  if(rho_max>0)
                    {
                      std::tuple<u3::U3,u3::U3,u3::SU3,int,int> label(omega1,omega2,x0,N1,N2);
                      spatial_labels[label]=upsilon_max1*upsilon_max2*rho_max;
                    }
                }
            }
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////        
        std::map<CompoundSpinTuple,int> lgi_spin_labels;
        for(const auto& [gamma_max1,Sp1,Sn1,S1] : upstreams1)
          for(const auto& [gamma_max2,Sp2,Sn2,S2] : upstreams2)
            {
              SpinTuple ket_spins(Sp2,Sn2,S2);
              SpinTuple bra_spins(Sp1,Sn1,S1);
              // Selection rules on spin-proton and spin-neutron
              if(abs(Sp1-Sp2)>2)
                continue;

              if(abs(Sn1-Sn2)>2)
                continue;

              for(const auto& spins : spin_tuples)
                {
                  const auto& [Sbar1,Sbar2,S0]=spins;
                  if(not am::AllowedTriangle(S1,S0,S2))
                    continue;

                  for(const auto& isospins : spin_tuples)
                    {
                      CompoundSpinTuple compound_label(bra_spins,ket_spins,spins,isospins);
                      lgi_spin_labels[compound_label]=gamma_max1*gamma_max2;
                    }
                }
            }

        //Sum up spins and upstream quantum numbers
        long int total_lgi_spatial=0;
        for(const auto& it : lgi_spatial_labels)
           total_lgi_spatial+=it.second;


        long int total_spins=0;
        for(const auto& it : lgi_spin_labels)
           total_spins+=it.second;

        long int total_spatial=0;
        for(const auto& it : spatial_labels)
           total_spatial+=it.second;


        std::cout<<fmt::format("{:12}  {:12}  {:10}  {:10}  {:10} ",
          sigma1.Str(),sigma2.Str(),total_lgi_spatial,total_spatial,total_spins
          )<<std::endl;
        
        /// Now do spin 



          
      }

  



  return 0;
}

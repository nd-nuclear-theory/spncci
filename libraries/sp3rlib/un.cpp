/****************************************************************
  un.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  based heavily on lsu3shell CUNMASTER
****************************************************************/
#include "sp3rlib/un.h"

#include "cppformat/format.h"
  

namespace un
{

  typedef std::tuple<int,int,int> U3Weights;

  /*
   * INPUT:
   * (1) weights: 
   * A weight vector of a basis state of U(N) 
   * (2) single_particle_states 
   * a tuple of single particle states of HO in speedometer order for a given shell 
   * A distribution of HO quanta in (z, x, y) directions associated with each of
   * N levels of U(N)  
   *
   * TASK: 
   * calculate labels of U(3) irrep using equation (4) of CPC.
   *
   * OUTPUT: 
   * U(3) labels [N_{z}, N_{x}, N_{y}]
   */

  U3Weights UNWeightToU3Labels(const BasisStateWeightVector& weights, const std::vector<SingleParticleState>& single_particle_states)
  	{
      int nx=0,ny=0,nz=0;
      for(int i=0; i<single_particle_states.size(); i++)
        {
          int spx,spy,spz;
          std::tie(spz,spx,spy)=single_particle_states[i];
          nx+=weights[i]*spx;
          ny+=weights[i]*spy;
          nz+=weights[i]*spz;
        }
      return U3Weights(nz,nx,ny);
  	}


  void GenerateU3Labels(
      const UNLabels& un_labels,
      const std::vector<SingleParticleState>& single_particle_states,
      BasisStateWeightVector& weights,
      std::map<u3::U3,int>& u3_un_multiplicity_map
    )
  {
    int f1,f2,f3;

    int N=un_labels.size()-1;
    
    int sum_galfand_parent_row=0;

    // sum of a parent Gelfand row
    for(int i=0; i<N+1; i++)
      sum_galfand_parent_row+=un_labels[i];

    // evaluate allowed Gelfand patterns based on a parent Gelfand row 
    // and store them in allowed_gelfand_patterns
    std::vector<UNLabels>allowed_gelfand_patterns;
    // number of allowed Gelfand patterns
    int num_allowed_patterns=1;
    //not sure what this is for yet
    std::vector<int>n_labels(N);
    for(int i=0; i<N; i++)
      {
        UNLabels gelfand_pattern;
        int label_min=std::min(un_labels[i], un_labels[i+1]);
        int label_max=std::max(un_labels[i], un_labels[i+1]);
        for(int label=label_min; label<=label_max; gelfand_pattern.push_back(label++));

        allowed_gelfand_patterns.push_back(gelfand_pattern);
        n_labels[i]=gelfand_pattern.size();
        num_allowed_patterns*=gelfand_pattern.size();
      }

    std::vector<unsigned> elements_per_change(N,1);
    for(int i=1; i<N; i++)
      {
        elements_per_change[i]=elements_per_change[i-1]*n_labels[i-1];
      }

    UNLabels gelfand_row(N);
    for(unsigned index=0; index<num_allowed_patterns; index++)
      {
        int sum_gelfand_row=0;
        for(int i=0; i<N; i++)
          {
            int element=(index/elements_per_change[i])%n_labels[i];

            gelfand_row[i]=allowed_gelfand_patterns[i][element];
            sum_gelfand_row+=gelfand_row[i];
          }
        weights[N]=(sum_galfand_parent_row-sum_gelfand_row);
        if(N>1)
          // This condition accomodates the case n=0 case when N=0 and hence one wants to keep
          // weights[0]=sum_gelfand_row 
          GenerateU3Labels(gelfand_row,single_particle_states, weights,u3_un_multiplicity_map);
        else
          {
            if(N==1)
              weights[0]=sum_gelfand_row;
            U3Weights u3_weights=UNWeightToU3Labels(weights,single_particle_states);
            std::tie(f1,f2,f3)=u3_weights;
            u3::U3 w(f1,f2,f3);
            u3_un_multiplicity_map[w]+=1;
          }

      }
  }

  unsigned UNBranchingMultiplicity(const u3::U3 w, const std::map<u3::U3,int>& u3_un_multiplicity_map)
  {
      HalfInt f1(w.f1());
      HalfInt f2(w.f2());
      HalfInt f3(w.f3());

      unsigned multiplicity;
      
      auto it=u3_un_multiplicity_map.find(w);
      auto end=u3_un_multiplicity_map.end();
      if (it == end)
        {
          std::cout <<fmt::format("Error: U3 irrep {} is not present in list of allowed U3 irreps",w.Str())<<std::endl;
          std::cout << std::endl;
          exit(0);
        }
      multiplicity = it->second;

      it=u3_un_multiplicity_map.find(u3::U3(f1+1,f2+1,f3-2));
      multiplicity+=(it==end)? 0: it->second;

      it=u3_un_multiplicity_map.find(u3::U3(f1+2,f2-1,f3-1));
      multiplicity+=(it==end)? 0: it->second;

      it=u3_un_multiplicity_map.find(u3::U3(f1+2,f2,f3-2));
      multiplicity-=(it==end)? 0: it->second;

      it=u3_un_multiplicity_map.find(u3::U3(f1+1,f2-1,f3));
      multiplicity-=(it==end)? 0: it->second;

      it=u3_un_multiplicity_map.find(u3::U3(f1,f2+1,f3-1));
      multiplicity-=(it==end)? 0: it->second;

      return multiplicity;
  }


  void GenerateAllowedSU3xSU2Irreps(
          const unsigned n, const unsigned A, 
          const std::vector<SingleParticleState>& single_particle_states,
          MultiplicityTagged<u3::U3S>::vector& allowed_irreps
        )
  {
    const unsigned N = (n+1)*(n+2)/2;
    if (2*N < A) 
      {
        std::cout << "Error: only " << 2*N << " states available for " << A << " particles!" << std::endl;
        exit(-1);
      }

    HalfInt S_min=!(A%2)?0:HalfInt(1,2);
    HalfInt S_max=HalfInt(A,2);
    unsigned branching_multiplicity;

    for (HalfInt S=S_min; S<=S_max; S++)
      {
        std::vector<int> u2_labels(2);
        u2_labels[0]=(A+TwiceValue(S))/2;
        u2_labels[1]=(A-TwiceValue(S))/2;

        if (u2_labels[0] > N)
          continue;

        UNLabels un_labels(N,0);
        for (size_t i = 0; i < u2_labels[0]; i++)
            un_labels[i] = (i < u2_labels[1]) ? 2 : 1;

        //  Iterate over U(3) labels and calculate multiplicity Mult (corresponds to
        //  alpha in [f] alpha (lm mu)S).  if Mult > 0 => add a given SU(3)xSU(2) irrep into SU3xSU2Labels
        std::map<u3::U3,int>u3_un_multiplicity_map;      
        un::GenerateU3Labels(un_labels, single_particle_states, u3_un_multiplicity_map);
        //  iterate over all [N_{z}, N_{x}, N_{y}] labels spanning U(N) irrep 
        for(auto it=u3_un_multiplicity_map.begin(); it!=u3_un_multiplicity_map.end(); it++)
          {
            u3::U3 w(it->first);
            if (w.Valid()) // this could be U(3) irrep ... let's check it
              {
                branching_multiplicity = UNBranchingMultiplicity(w, u3_un_multiplicity_map);
                //  Calculate multiplicity alpha in irrep of U(N)
                if (branching_multiplicity) 
                  allowed_irreps.push_back(MultiplicityTagged<u3::U3S>(u3::U3S(w,S),branching_multiplicity)); 
                   //SU3MULTList.push_back(UN::SU3(usMult, SU3::LM(U3Labels), SU3::MU(U3Labels)));
              }
          }
      }
  }

  void GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, MultiplicityTagged<u3::U3S>::vector& allowed_irreps)
  {
    // set of sps (Nz, Nx, Ny) for harmonic oscillator shell n
    std::vector<SingleParticleState> single_particle_states;
    int N = (n + 1) * (n + 2) / 2;
    // U3Labels; // [0] -> nz ; [1] -> nx ; [2] -> ny
    
    for (long k = 0; k <= n; k++)
        for (long nx = k; nx >= 0; nx--)
        {
            SingleParticleState sp(n-k,nx,k-nx);
            single_particle_states.push_back(sp);
        }
    GenerateAllowedSU3xSU2Irreps(n, A, single_particle_states, allowed_irreps);
  }

}//end namespace 

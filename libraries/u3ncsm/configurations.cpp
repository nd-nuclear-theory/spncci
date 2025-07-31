/****************************************************************
  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3ncsm/configurations.h"

#include <unordered_set>
#include <cassert>

namespace shell{

  // For a given Nmax, determine maximum number of particles that be in a given shell 
  // in any configurations with Nex<=Nmax
  std::vector<unsigned int> GetShellCapacitiesNmax(const unsigned int Nmax, const shell::Configuration& config0)
  {
    
    int Nshells=config0.size();
    std::vector<unsigned int> shell_capacities_Nmax(Nshells);

    // Get N of valence space 
    int Nv=0; 
    while(
      (config0[Nv+1]>0)&&((Nv+1)<config0.size())
      ){Nv++;}

    // Getting maximum number of particles allowed in a given shell 
    std::vector<int> shell_capacities(Nshells);
    for (int n=0; n<Nshells; n++)
      shell_capacities[n]=shell::Omega(n);

    // For each shell, determine maximum number of particles it will have in an Nmax truncated basis
    for(int Nshell=0; Nshell<Nshells; Nshell++)
      {
        int num_particles = (Nshell>Nv)? std::min(shell_capacities[Nshell],int(Nmax/(Nshell-Nv))): shell_capacities[Nshell];
        // If number of particles is greater than number in valence space, 
        // we need to take into account particles coming from closed shells
        if((num_particles > config0[Nv]) && (Nshell>=Nv))
          {
            num_particles=config0[Nv];
            int Nex=config0[Nv]*(Nshell-Nv);
            for(int n=Nv-1; n>=0; n--)
              {
                for(int i=0; i<config0[n]; ++i)
                  {
                    Nex+=(Nshell-n);
                    
                    if(Nex>Nmax)
                      {
                        n=0;
                        break;
                      }
                      
                    num_particles++;
                  }
              }
          }

        shell_capacities_Nmax[Nshell]=num_particles;
      }
    return shell_capacities_Nmax;
  }


  // Identify the configuration with the particles in the shells with
  // fewest number of oscillator quanta.
  shell::Configuration lowest_pauli_allowed_config(const unsigned int A)
    {    
      int particles_remaining=A;
      shell::Configuration configuration;
      
      int Nshell=0;
      while (particles_remaining>0)
      {
        int omega_shell=shell::Omega(Nshell); 
        configuration.push_back(std::min(omega_shell,particles_remaining));
        particles_remaining+=(-omega_shell);
        Nshell++;
      }
      return configuration;
    }

  // int SumVector(const std::vector<int>& config)
  //   {
  //     int sum=0;
  //     for(auto n : config)
  //       sum+=n;

  //     return sum;
  //   }

  shell::ShellConfigurations get_configurations_Nex(
    const unsigned int A,
    const std::vector<unsigned int>& shell_capacities,
    const shell::ShellConfigurations& seed_configurations
    )
    {
      shell::ShellConfigurations configuration_set;
      for(const auto& seed_config : seed_configurations)
        {
          for(std::size_t i=0; i<seed_config.size()-1; ++i)
            {
              shell::Configuration config=seed_config; 
              --config[i];
              ++config[i+1];
                            
              if (config[i]>=0 && config[i+1]<=shell_capacities[i+1])
              {
                configuration_set.insert(config);
              }
            }
        }
      return configuration_set;
    }

std::vector<shell::ShellConfigurations>
  generate_configurations(const unsigned int A, const unsigned int Nmax, const unsigned int Nshell_max)
  {
    // Initialize container
    std::vector<shell::ShellConfigurations> configurations_by_Nex(Nmax+1);
  
    // Get lowest Pauli allowed configuration
    shell::Configuration configuration0=shell::lowest_pauli_allowed_config(A);
    
    // Determine number of active shells based on number of shells in Nex=0 configuration and Nmax
    // unsigned int num_shells = configuration0.size()+Nmax;
    
    //Resizing lowest pauli allowed configuration to add additional shells
    assert(configuration0.size()<=Nshell_max);
    configuration0.resize(Nshell_max+1,0);

    // Getting maximum number of particles allowed in a given shell 
    std::vector<unsigned int> shell_capacities(Nshell_max+1);
    for (int n=0; n<=Nshell_max; n++)
      shell_capacities[n]=shell::Omega(n);
    
    //Resize container and add Nex=0 configuration 
    configurations_by_Nex[0]={configuration0};

    // Generate remainder of the configuration
    for(int Nex=1; Nex<=Nmax; Nex++)
        configurations_by_Nex[Nex] = shell::get_configurations_Nex(A,shell_capacities,configurations_by_Nex[Nex-1]);         
  
    return configurations_by_Nex;
  }


}//end shell namespces


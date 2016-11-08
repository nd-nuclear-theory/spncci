/****************************************************************
  import_jisp16.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/7/16 (aem,mac): Created to import relative jisp16 files.
****************************************************************/
//#include "u3shell/import_interaction.h"

#include <fstream>

#include "cppformat/format.h"
#include "eigen3/Eigen/Eigen" 
#include "mcutils/parsing.h"

#include "basis/lsjt_scheme.h"
#include "utilities/utilities.h"

Eigen::MatrixXd PopulateMatrix(
    const std::vector<double>& matrix_elements,
    int dim,
    int sector_dim,
    std::string convention
  )
  // Takes a vector of matrix elements of upper triangle, stored in row-major order
  // and stores them in a matrix up to the dimensions given by 
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  //Entering rme's from previous sector into symmetry sector matrix.
  ////////////////////////////////////////////////////////////////////
  // check that previous sector has non-zero entries and if non-zero fill in matrix in vector
  int nmax;
  Eigen::MatrixXd sector;
  bool zeros = std::all_of(matrix_elements.begin(), matrix_elements.end(), [](int i) { return i==0; });
  if (not zeros)
    {
      if (sector_dim)
        {
          nmax=sector_dim;
          sector.resize(nmax+1,nmax+1);
        }
      else
        {
          nmax=dim-1;
          sector.resize(dim,dim);
        }

      int pos=0;
      if(convention=="positive_at_infinity")
        {
          for (int np=0; np<=nmax; np++)
              for(int n=0; n<=np; n++)
              {
                sector(np,n)=matrix_elements[pos]*parity(n+np);
                sector(n,np)=matrix_elements[pos]*parity(n+np);
                ++pos;
              }
        }
      else if(convention=="positive_at_origin")
        {
          for (int np=0; np<=nmax; np++)
              for(int n=0; n<=np; n++)
              {
                sector(np,n)=matrix_elements[pos];
                sector(n,np)=matrix_elements[pos];
                ++pos;
              }
        }
      else
        std::cout<<"invalid convention"<<std::endl;
    }
  return sector;
}


std::vector<Eigen::MatrixXd>
ImportInteraction_JISP(std::string interaction_file, 
  const basis::RelativeSpaceLSJT& space,
  const basis::RelativeSectorsLSJT& sectors,
  int sector_dim,
  std::string convention
  )
{
  /////////////////////////////////////////////////////////////////////////////////////////////
  // If size matrices for each sector are desired, set nmax to a desired dimension. 
  // For example, setting nmax=20 would correspond to 
  // Nmax=40+L.  If all matrix elements are desired, set nmax=0;
  //
  // There are two conventions for rme's "positive_at_origin" and "positive_at_infinity".  
  /////////////////////////////////////////////////////////////////////////////////////////////
	int index;
  std::vector<Eigen::MatrixXd> sector_vector;
  // set size of vector 
  sector_vector.resize(sectors.size());
  // Reading in matrix elements 
  std::ifstream interaction_stream(interaction_file.c_str());
  OpenCheck(bool(interaction_stream),interaction_file);
  std::string line;
  int J, S, Lp, L, dim,ipcut, indent, nmax;
  dim=0;
  double mass, hbar_omega;
  int line_count=0;
  std::vector<double> matrix_elements;
  int num_rows_7, num_rows_1;
  while(std::getline(interaction_stream,line))
    {
      std::istringstream line_stream(line);
      // New sector, obtain symmetry sector labels J, S, L', L
      if(line_count==0)
        {
          /////////////////////////////////////////////////////////////////////////////////////////////////////
          //Entering rme's from previous sector into symmetry sector matrix.
          sector_vector[index]=PopulateMatrix(matrix_elements,dim, sector_dim,convention);

          /////////////////////////////////////////////////////////////////////////////////////////////////////
          // Reading in symmetry sector labels for next sector
          line_stream >> J >> S >> Lp >> L >> ipcut >> dim >> mass >> hbar_omega >> indent;
          //Check for end of file
          if (J==99999)
            continue;
          // Number of matrix elements stored (upper triangle) dim(dim+1)/2
          int num_matrix_elements=dim*(dim+1)/2;
          num_rows_7=num_matrix_elements/7;
          num_rows_1=num_matrix_elements%7;
          // Assuming here that the dim is always 119 or 120.  
          // If 120 then there will be an extra line with a single rme.  
          assert((num_matrix_elements%7==1)||(num_matrix_elements%7==0));
          //////////////////////////////////////////////////////////////////////////////////////////
          // Look up sector index
          // parity is forced by L~g
          int g=L%2;
          int gp=Lp%2;
          // isospin forced by L+S+T requirement
          int T=(L+S+1)%2;

          basis::RelativeSubspaceLSJTLabels ket_labels(L,S,J,T,g);
          basis::RelativeSubspaceLSJTLabels bra_labels(Lp,S,J,T,gp);
          int ket_subspace_index=space.LookUpSubspaceIndex(ket_labels);
          int bra_subspace_index=space.LookUpSubspaceIndex(bra_labels);
          index=sectors.LookUpSectorIndex(bra_subspace_index, ket_subspace_index);
          //////////////////////////////////////////////////////////////////////////////////////////
          matrix_elements.clear();
          line_count+=1;
        }
      // 7 matrix elements per row
      else if(line_count<=num_rows_7)
        {
          double me1, me2,me3,me4,me5,me6,me7;
          line_stream >> me1 >> me2 >> me3 >> me4 >> me5 >> me6 >> me7;
          std::vector<double>me_array={me1,me2,me3,me4,me5,me6,me7};
          matrix_elements.insert(matrix_elements.end(), me_array.begin(), me_array.end());
          line_count+=1;
          if(line_count==num_rows_7+1 && (num_rows_1==0))
            // Finished accumlating rme's into vector, reinitialze to zero
            line_count=0;
        }
      else if(num_rows_1)
        {
          double me;
          line_stream >> me;
          matrix_elements.push_back(me);
          line_count=0;
        }
    }//end while
  return sector_vector;
}//end function

int main(int argc, char **argv)
  {
    int J0=0;
    int T0=0;
    int g0=0;
    int sector_dim=0;
    std::string interaction_file="data/Vrel_JISP16_bare_Jmax4.hw20";

    std::string convention="positive_at_origin";
    // std::string convention="positive_at_infinity";

    // Nmax set to 5 since Lmax=5;
    int Nmax=5;
    int Jmax=4;
    basis::RelativeSpaceLSJT relative_space(Nmax,Jmax);
    basis::RelativeSectorsLSJT relative_sectors(relative_space, J0,T0, g0);


    std::cout<<"importing interaction"<<std::endl;
    std::vector<Eigen::MatrixXd> sector_vector
        =ImportInteraction_JISP(
              interaction_file,relative_space, 
              relative_sectors,
              sector_dim, 
              convention
            );

    for(int i=0; i<sector_vector.size(); ++i)
      {
        std::cout<<sector_vector[i]<<std::endl<<std::endl;
      }
    return 0; 
  }




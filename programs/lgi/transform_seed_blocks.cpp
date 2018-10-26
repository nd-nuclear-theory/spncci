/****************************************************************
transform_seed_blocks.cpp

Code which transforms seed blocks based on different linear combinations
of irreps in the same irrep family 

7/19/18 (aem): Created. 
****************************************************************/
#include <fstream>  

#include "cppformat/format.h"
#include "lgi/lgi_unit_tensors.h"
#include "spncci/io_control.h"
#include "spncci/transform_basis.h"


namespace spncci
{


void GetReorderingMatrix(const std::vector<int>& top_rows,
  const std::vector<int>& bottom_rows,
  spncci::OperatorBlock& reordering_matrix)
{
  int num_rows=top_rows.size()+bottom_rows.size();
  reordering_matrix=spncci::OperatorBlock::Zero(num_rows,num_rows);
  int i=0;
  for(int j : top_rows)
  {
    reordering_matrix(i,j)=1;
    ++i;
  }

  for(int j : bottom_rows)
    {
      reordering_matrix(i,j)=1;
      ++i;
    }
}


void ReorderBlock(
  spncci::OperatorBlock& block,
  spncci::OperatorBlock& reordering_matrix,
  int& num_above_threshold,
  double threshold
)
{
  int gamma_max=block.rows();
  auto norms=block.rowwise().squaredNorm();

  std::vector<int> above_threshold,below_threshold;

  // Determine which irreps are above the threshold 
  for(int g=0; g<gamma_max; ++g)
    {
      if(norms[g]>threshold)
        above_threshold.push_back(g);
      else
        below_threshold.push_back(g);
    }

  spncci::GetReorderingMatrix(above_threshold,below_threshold,reordering_matrix);
  
  //Reorder the block
  block=reordering_matrix*block;

  num_above_threshold=above_threshold.size();
}


// void GetUnitaryTransformationForColumn(
//   spncci::OperatorBlock& block, 
//   spncci::OperatorBlock& transformation_matrix, 
//   int& num_rows,
//   int column
// )
// {
//   int gamma_max=block.rows();

//   //Define pairwise transformation 
//   std::vector<int> top_rows, bottom_rows;
//   spncci::OperatorBlock temp_matrix=spncci::OperatorBlock::Identity(gamma_max,gamma_max);

//   //make num rows even
//   if(num_rows%2)
//     top_rows.push_back(num_rows-1);

//   for(int i=0; i<num_rows-1; i+=2)
//     {
//       double a=block(i,column);
//       double b=block(i+1,column);
//       double beta=sqrt(a*a+b*b);
//       temp_matrix(i,i)=a/beta;
//       temp_matrix(i,i+1)=b/beta;
//       temp_matrix(i+1,i)=-b/beta;
//       temp_matrix(i+1,i+1)=a/beta;
//       top_rows.push_back(i);
//       bottom_rows.push_back(i+1);
//     } 

//   //Finish add rows that were below threshold to start with 
//   for(int i=num_rows; i<gamma_max; ++i)
//     bottom_rows.push_back(i);

//   spncci::OperatorBlock reordering_matrix;
//   spncci::GetReorderingMatrix(top_rows,bottom_rows,reordering_matrix);
//   temp_matrix=reordering_matrix*temp_matrix;

//   block=temp_matrix*block;

//   transformation_matrix=reordering_matrix*temp_matrix*transformation_matrix;
//   num_rows=top_rows.size();
// }


  void ThresholdTruncation(
    spncci::OperatorBlock& block, 
    spncci::OperatorBlock& transformation_matrix,
    double threshold
  )
  {
    int gamma_max=block.rows();

    //Reorder block so irreps with probably above threshold are in the top rows of the block 
    spncci::OperatorBlock reordering_matrix;
    int num_above_threshold;
    spncci::ReorderBlock(block,reordering_matrix,num_above_threshold,threshold);
    transformation_matrix=reordering_matrix.block(0,0,num_above_threshold,gamma_max);
  }


void UnitaryTransformBlock(
  spncci::OperatorBlock& block, 
  spncci::OperatorBlock& transformation_matrix,
  const std::pair<std::string,double>& truncation_mode
  )
{  
  Eigen::JacobiSVD<spncci::OperatorBlock> svd(block,Eigen::ComputeFullU);
  svd.setThreshold(truncation_mode.second);
  transformation_matrix=svd.matrixU().transpose();
  block=transformation_matrix*block;

}



void GetTransformation(
  spncci::OperatorBlock& block, 
  spncci::OperatorBlock& transformation_matrix,
  const std::pair<std::string,double>& truncation_mode
  )
{
  double threshold=truncation_mode.second;
  bool apply_unitary_transformation=false;

  if(truncation_mode.first=="UI" || truncation_mode.first=="IR")
    {

      if(truncation_mode.first=="UI")
        apply_unitary_transformation=true;

      double norm=block.squaredNorm();
      int gamma_max=block.rows();
      int num_cols=block.cols();
      // Maximum probablity of a single irrep
      // std::cout<<block<<std::endl;
      double max_probability=block.rowwise().squaredNorm().maxCoeff();
      std::cout<<"gamma_max:  "<<gamma_max<<"  norm:  "<<norm<<std::endl;
      std::cout<<"initial max probability "<<max_probability<<std::endl;
      std::cout<<block.rowwise().squaredNorm()<<std::endl<<std::endl;    

      // auto temp_block=block;
      
      spncci::OperatorBlock unitary_transformation=spncci::OperatorBlock::Identity(gamma_max,gamma_max);
      if(apply_unitary_transformation==true)
        spncci::UnitaryTransformBlock(block, unitary_transformation, truncation_mode);
      
      // // Eigen::JacobiSVD<spncci::OperatorBlock> svd(block,ComputeThinV);
      // Eigen::JacobiSVD<spncci::OperatorBlock> svd(block,Eigen::ComputeFullU);
      // svd.setThreshold(truncation_mode.second);
      // spncci::OperatorBlock Umatrix=svd.matrixU().transpose();
      // // spncci::OperatorBlock temp_transformation=svd.matrixU().transpose();
      // block=Umatrix*block;
      
      //Reorder block so irreps with probably above threshold are in the top rows of the block 
      spncci::OperatorBlock reordering_matrix;
      int num_above_threshold;
      spncci::ReorderBlock(block,reordering_matrix,num_above_threshold,threshold);

      std::cout<<"populating transformation_matrix"<<std::endl;
      auto temp_transformation=reordering_matrix*unitary_transformation;
      // std::cout<<temp_transformation.rows()<<"  "<<temp_transformation.cols()<<" "<<num_above_threshold<<
      transformation_matrix=temp_transformation.block(0,0,num_above_threshold,gamma_max);
      // std::cout<<transformation_matrix<<std::endl<<std::endl;
      // temp_block=transformation_matrix*temp_block;
      // std::cout<<temp_block.rowwise().squaredNorm()<<std::endl<<std::endl;

      const auto& probabilities=block.rowwise().squaredNorm();
      std::cout<<"final max probabilities 1 "<<probabilities.maxCoeff()<<std::endl;
      std::cout<<probabilities<<std::endl;
      std::cout<<"num above threshold "<<num_above_threshold<<std::endl;
      std::cout<<"---------------------------"<<std::endl<<std::endl;
    }

  else if(truncation_mode.first=="IF")
    {
      double norm=block.squaredNorm();
      std::cout<<"norm is "<<norm<<"  "<<threshold<<std::endl;
      if(norm>threshold)
        {
          int gamma_max=block.rows();
          transformation_matrix=spncci::OperatorBlock::Identity(gamma_max,gamma_max);
          std::cout<<"keep "<<gamma_max<<std::endl;
        }
      else
        {
          transformation_matrix=spncci::OperatorBlock(0,0); //null matrix
          std::cout<<"keep "<<0<<std::endl;
        }
    }
}




  void  DefineIrrepFamilyTransformations2(
  const std::vector<std::pair<int,int>>& Jn_set,
  std::map<int,std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
  std::map<int,std::map<int,int>>& J_index_lookup_table,
  spncci::OperatorBlocks& transformations,
  const std::pair<std::string,double>& truncation_mode,
  double threshold
)
//Set of transformations for a given set of Jn pairs
{
  
  // double threshold=truncation_mode.second;
  
  std::cout<<"defining rotation"<<std::endl;
  int num_irrep_families=irrep_family_blocks.size();
  transformations.resize(num_irrep_families);

  for(auto it=irrep_family_blocks.begin(); it!=irrep_family_blocks.end(); ++it)
  // for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
    {
      int irrep_family_index=it->first;
      std::cout<<"irrep family "<<irrep_family_index<<std::endl;
      const std::vector<spncci::OperatorBlocks>& blocks=it->second;
      // const std::vector<spncci::OperatorBlocks>& blocks=irrep_family_blocks[irrep_family_index];
      // std::cout<<irrep_family_index<<"  "<<blocks.size()<<std::endl;
      if(blocks.size()==0)
        continue;
      
      spncci::OperatorBlock block;
      std::map<int,int>& J_index_table=J_index_lookup_table[irrep_family_index];
      
      std::cout<<"regrouping"<<std::endl;
      spncci::RegroupBlocks(Jn_set, blocks,J_index_table, block);
      // std::cout<<block<<std::endl<<std::endl;  


      int gamma_max=block.rows();
      int irrep_dim=block.cols();
      if(gamma_max==0)
        continue;

      // std::cout<<"get transformation index "<<std::endl; 
      spncci::OperatorBlock& transformation_matrix=transformations[irrep_family_index];

      std::cout<<"transforming"<<std::endl;
      // std::cout<<num_irrep_families<<"  "<<irrep_family_index<<std::endl;
      spncci::GetTransformation(block,transformation_matrix,truncation_mode);
      std::cout<<"----------------------"<<std::endl<<std::endl;
    }

}

}//namespace



int main(int argc, char **argv)
{

  if(argc<4)
  {
    std::cout<<"Syntax: <truncation_file_num> Nsmax Nmax "<<std::endl;
    exit(0);
  }

  int truncation_file_num=std::stoi(argv[1]);
  int Nsmax=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);



  std::array<int,2> nuclide;
  std::vector<std::pair<int,int>> Jn_set;
  int num_Jn;
  double threshold;

  // std::vector<std::pair<int,int>> Jn_set;
  // Jn_set.emplace_back(1,0); // J,n where n is zero based
  // Jn_set.emplace_back(3,0); // J,n where n is zero based

  // double threshold=1e-6*Jn_set.size();

  // std::pair<std::string,double> truncation_mode("None",threshold);
  // std::pair<std::string,double> truncation_mode("Rank",threshold);
  // std::pair<std::string,double> truncation_mode("Threshold",threshold);


  std::ifstream is;
  std::string truncation_file_name=fmt::format("truncation_{:03d}.dat",truncation_file_num);
  is.open(truncation_file_name);

  if(not is)
  {
    std::cout<<"truncation file not found "<<std::endl;
    exit(0);
  }

  // Read in nuclide information
  is>>nuclide[0]>>nuclide[1];

  // Read in truncation information
  std::string truncation_type;
  is>>truncation_type>>threshold>>num_Jn;

  if(truncation_type!="UI" && truncation_type!="IR" && truncation_type!="IF")
    { 
      std::cout<<truncation_type<<" "<< "is and invalid truncation type.  Possible truncation types are "<<std::endl
      <<"unitary_irrep : UI"<<std::endl<<"irrep : IR"<<std::endl<<"irrep_families : IF"<<std::endl;
      std::exit(EXIT_FAILURE);
    }

  std::pair<std::string,double> truncation_mode(truncation_type,threshold);
  
  std::cout<<truncation_type<<"  "<<threshold<<std::endl;
  // Get list of J,n pairs which will contrubute to definition of transformation matrices 
  std::string line;
  int twice_J,n;
  for(int j=0; j<num_Jn; j++)
    {
      is>>twice_J>>n;
      std::cout<<twice_J<<"  "<<n<<std::endl;
      Jn_set.emplace_back(twice_J,n);
    }
  is.close();

  std::cout<<"num Jn included "<<Jn_set.size()<<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////


  // Read in irrep family blocks 
  std::cout<<"read in irrep family blocks"<<std::endl;
  std::map<int,std::vector<spncci::OperatorBlocks>> irrep_family_blocks;
  // std::vector<std::vector<spncci::OperatorBlocks>> irrep_family_blocks;
  std::map<int,std::map<int,int>> J_index_lookup_table;
  std::string irrep_family_blocks_file="irrep_family_blocks_20";
  spncci::ReadIrrepFamilyBlocks(irrep_family_blocks,J_index_lookup_table,irrep_family_blocks_file);


  // Get transformation matrices for each irrep family
  std::cout<<"defined transformations "<<std::endl;
  spncci::OperatorBlocks transformations;
  // std::map<int,spncci::OperatorBlock> transformations;


  spncci::DefineIrrepFamilyTransformations2(
    Jn_set,irrep_family_blocks,J_index_lookup_table,
    transformations,truncation_mode,threshold
  );

  // write truncated lgi families to file 
  std::string truncated_lgi_filename
    =fmt::format("lgi_families_truncated_{:02d}_{:02d}_{:03d}.dat",Nsmax,Nmax,truncation_file_num);
  spncci::WriteTruncatedLGIs(nuclide,transformations,truncated_lgi_filename);

   // write transformations to file
  int Jn_file_num=Jn_set.size();
  std::string transformations_file=fmt::format("transformations_{:02d}_{:02d}_{:03d}.dat",Nsmax,Nmax,truncation_file_num);
  std::cout<<"write transformation matrices "<<std::endl;
  spncci::WriteTransformationMatrices(transformations,transformations_file);

  //Read transformations from file 
  std::cout<<"Read transformation matrices "<<std::endl;
  spncci::OperatorBlocks transformations2;
  spncci::ReadTransformationMatrices(transformations_file,transformations2);

}
/****************************************************************
  vcs.h                       

  Define vector coherent state methods for Sp(3,R).

  Anna E. McCoy
  University of Notre Dame

  Created by Anna E. McCoy on 3/9/16.   
 
  3/9/16 (aem): Created based on prototype vcs.py.

****************************************************************/
#include <Cmath>  
#include "u3.h"
#include <eigen3/Eigen/Eigen>


namespace vcs{
  double BosonCreationRME(const u3::U3& np, const u3::U3& n)
	//  SU(3) Reduced matrix element of a^\dagger boson creation operator
  	{
  		double rme=0;
  		int n1=int(n.f1);
  		int n2=int(n.f2);
  		int n3=int(n.f3);
  		if((np.f1==(np.f1+2))&&(np.f2==n.f2)&&(np.f3==n.f3))
  			rme=sqrt(
  				(n1+4)*(n1-n2+2)*(n1-n3+3)
  				/(2.*(n1-n2+3)*(n1-n3+4))
  				);
  		if ((np.f1==n.f1)&&(np.f2==(n.f2+2))&&(np.f3==n.f3))
  			rme=sqrt(
  				(n2+3)*(n1-n2)*(n2-n3+2)
  				/(2.*(n1-n2-1)*(n2-n3+3))
  				);
  		if ((np.f1==n.f1)&&(np.f2==n.f2)&&(np.f3==(n.f3+2)))
  			rme=sqrt(
  				(n3+2)*(n2-n3)*(n1-n3+1)/(2.*(n1-n3)*(n2-n3-1))
  				);
  		return rme;

  	}

  double SMatrix(u3::U3& s, u3::U3& s, MultiplicityTagged<u3::U3>& nr1, MultiplicityTagged<u3::U3>& nr2)
  {
    double smatrix=0;


    return smatrix;
  }



}  //  namespace 

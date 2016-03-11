/****************************************************************
  vcs.h                       

  Define vector coherent state methods for Sp(3,R).

  Anna E. McCoy
  University of Notre Dame

  Created by Anna E. McCoy on 3/9/16.   
 
  3/9/16 (aem): Created based on prototype vcs.py.

****************************************************************/

#include <cmath>  
//#include <eigen3/Eigen/Eigen>

#include "sp3rlib/u3.h"
#include "sp3rlib/vcs.h"



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

  double SMatrix(u3::U3& s, u3::U3& wp, MultiplicityTagged<u3::U3>& nr1p, MultiplicityTagged<u3::U3>& nr2p)
  {
    double smatrix=0.0;
    u3::U3 n1p=nr1p.irrep;
    u3::U3 n2p=nr2p.irrep;
    std::vector<u3::U3> n1set;
    std::vector<u3::U3> n2set;
    
    if (n1p.f1-2>=n1p.f2)
      n1set.push_back(u3::U3(n1p.f1-2,n1p.f2,n1p.f3));
    if (n1p.f2-2>=n1p.f3)
      n1set.push_back(u3::U3(n1p.f1,n1p.f2-2,n1p.f3));
    if (n1p.f3-2>=0)
      n1set.push_back(u3::U3(n1p.f1,n1p.f2,n1p.f3-2));
    
    if (n2p.f1-2>=n2p.f2)
      n2set.push_back(u3::U3(n2p.f1-2,n2p.f2,n2p.f3));
    if (n2p.f2-2>=n2p.f3)
      n2set.push_back(u3::U3(n2p.f1,n2p.f2-2,n2p.f3));
    if (n2p.f3-2>=0)
      n2set.push_back(u3::U3(n2p.f1,n2p.f2,n2p.f3-2));

    MultiplicityTagged<u3::U3>::vector wset=KroneckerProduct(wp, u3::U3(0,0,-2));
    for (int i=0; i<wset.size(); i++)
      {
        u3::U3 w=wset[i].irrep;
        for(int j=0; j<n1set.size(); j++)
          {
            u3::U3 n1=n1set[j];
            int r1_max=OuterMultiplicity(s.SU3(),n1.SU3(),w.SU3());
            for(int r1=1; r1<=r1_max; r1++)
              {
                MultiplicityTagged<u3::U3> nr1(n1,r1);
                for(int k=0; k<n2set.size(); k++)
                  {
                    u3::U3 n2=n2set[k];
                    int r2_max=OuterMultiplicity(s.SU3(),n2.SU3(),w.SU3());
                    for (int r2=1; r2<=r2_max; r2++)
                      {
                        MultiplicityTagged<u3::U3> nr2(n2,r12);
                        smatrix+=(
                          (Omega(n1p, wp)-Omega(n1,w))
                          *
                          *BosonCreationRME(n1p,n1)
                          )
                      }

                  }
              }
          }

      }





    return smatrix;
  }



}  //  namespace 

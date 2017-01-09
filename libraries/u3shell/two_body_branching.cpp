/****************************************************************
  two_body_branching.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/7/17 (aem): Created.
****************************************************************/
#include "u3shell/two_body_branching.h"
#include "am/wigner_gsl.h"
#include "cppformat/format.h"

extern double zero_threshold;

namespace u3shell
{

	void H2FormatLookUp(int Nmax,std::map<std::tuple<int,int,HalfInt>,int>& h2_lookup)
	  {
	    int index=1;
	    for(int N=0; N<=Nmax; ++N)
	      for(int L=(N%2); L<=N; L+=2)
	        for(HalfInt J=abs(L-HalfInt(1,2)); J<=(L+HalfInt(1,2)); ++J)
	          {
	            std::tuple<int,int,HalfInt> key(N,L,J);
	            h2_lookup[key]=index;
	            std::cout<<fmt::format("{} {} {}  {}",N, L, J, index)<<std::endl;
	            ++index;
	          }
	  }

	void BranchTwoBodyNLST(
	  u3shell::IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
	  std::map<TwoBodyBraketLST,double>& two_body_rmes_lst
	  )
	{
	  for(auto it=indexed_two_body_rmes.begin(); it!=indexed_two_body_rmes.end(); ++it)
	    {
	      u3shell::TwoBodyUnitTensorLabelsU3ST tensor_u3st;
	      int kappa0,L0;

	      u3shell::IndexTwoBodyTensorLabelsU3ST indexed_two_body_tensor=it->first;
	      double rme_u3st=it->second;
	      if(fabs(rme_u3st)<zero_threshold)
	        continue;
	      std::tie(tensor_u3st,kappa0,L0)=indexed_two_body_tensor;
	      u3::SU3 x0=tensor_u3st.x0();
	      HalfInt S0=tensor_u3st.S0();
	      HalfInt T0=tensor_u3st.T0();
	      int rho0=tensor_u3st.rho0();
	      u3shell::TwoBodyStateLabelsU3ST bra=tensor_u3st.bra();
	      u3shell::TwoBodyStateLabelsU3ST ket=tensor_u3st.ket();
	      int N1,N2,N1p,N2p;
	      HalfInt Sp,S,Tp,T;
	      u3::SU3 xp,x;
	      std::tie(N1p,N2p,xp,Sp,Tp)=bra.Key();
	      std::tie(N1,N2,x,S,T)=ket.Key();
	      MultiplicityTagged<int>::vector L_branch=BranchingSO3(x);
	      MultiplicityTagged<int>::vector Lp_branch=BranchingSO3(xp);
	      for(auto lp: Lp_branch)
	        {
	          int Lp=lp.irrep;
	          int kappap_max=lp.tag;
	          for(auto l: L_branch)
	            {
	              int L=l.irrep;
	              int kappa_max=l.tag;
	              // std::cout<<etap%2<<"  "<<etap<<"  "<<eta%2<<eta<<"  "<<std::endl;
	              for(int kappap=1; kappap<=kappap_max; ++kappap)
	                for(int kappa=1; kappa<=kappa_max; ++kappa)
	                  for(int L1p=N1p%2; L1p<=N1p; L1p+=2)
	                    for(int L2p=N2p%2; L2p<=N2p; L2p+=2)
	                      for(int L1=N1%2; L1<=N1; L1+=2)  
	                        for(int L2=N2%2; L2<=N2; L2+=2)
	                          {     
	                            if((abs(L1-L2)>L)||((L1+L2)<L)) //am::triangular
	                              continue;
	                            if((abs(L1p-L2p)>Lp)||((L1p+L2p)<Lp))
	                              continue;
	                            TwoBodyStateLabelsLST bra(N1p,L1p,N2p,L2p,Lp,Sp,Tp);
	                            TwoBodyStateLabelsLST ket(N1, L1,N2,L2,L,S,T);                            
	                            TwoBodyBraketLST braket(L0,S0,T0,bra,ket);
	                            int n1=(N1-L1)/2, n2=(N2-L2)/2,n1p=(N1p-L1p)/2, n2p=(N2p-L2p)/2;
	                            double rme_lst=rme_u3st*u3::W(u3::SU3(N1p,0),1,L1p,u3::SU3(N2p,0),1,L2p,xp,kappap,Lp,1)
	                                            *u3::W(u3::SU3(N1,0),1,L1,u3::SU3(N2,0),1,L2,x,kappa,L,1)
	                                            *u3::W(x,kappa,L,x0,kappa0,L0,xp,kappap,Lp,rho0)
	                                            *parity(n1+n2+n1p+n2p);
	                            two_body_rmes_lst[braket]+=rme_lst;
	                          }
	            }
	        }
	    }
	}
	void BranchTwoBodyLSJT( int Jmax, int J0,
	  const std::map<TwoBodyBraketLST,double>& two_body_rme_lst,
	  std::map<TwoBodyBraketLSJT,double>&two_body_rme_lsjt
	  )
		{
		  for(auto it=two_body_rme_lst.begin(); it!=two_body_rme_lst.end(); it++)
		    {
		      double rme_lst=it->second;
		      if (fabs(rme_lst)<zero_threshold)
		        continue;
		      TwoBodyStateLabelsLST bra_lst,ket_lst;
		      HalfInt S0,T0,S,T,Sp,Tp;
		      int eta1,eta2,eta1p,eta2p,L1,L2,L1p,L2p,L,Lp,L0;
		      std::tie(L0,S0,T0,bra_lst,ket_lst)=it->first;
		      std::tie(eta1,L1,eta2,L2,L,S,T)=ket_lst;
		      std::tie(eta1p,L1p,eta2p,L2p,Lp,Sp,Tp)=bra_lst;

		      int J_max=std::min(Jmax,int(L+S));
		      int Jp_max=std::min(Jmax,int(Lp+Sp));
		      for(HalfInt J=abs(L-S); J<=J_max; ++J)
		          for(HalfInt Jp=abs(Lp-Sp); Jp<=Jp_max; ++Jp)
		            {
		              if((J0<abs(Jp-J)) || (J0>(Jp+J)))
		                continue;

		              double so3coef=am::Unitary9J(L, S, J, L0,S0,J0,Lp,Sp,Jp);
		              TwoBodyStateLabelsLSJT state(eta1,L1,eta2,L2,L,S,J,T);
		              TwoBodyStateLabelsLSJT statep(eta1p,L1p,eta2p,L2p,Lp,Sp,Jp,Tp);
		              TwoBodyBraketLSJT braket(statep,state);
		              two_body_rme_lsjt[braket]+=so3coef*rme_lst;
		              }
		    }
	}

	void branchJJJT(
	  std::map<TwoBodyBraketLSJT,double>&two_body_rme_lsjt,
	  std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt
	  )
	{
	 for(auto it=two_body_rme_lsjt.begin(); it!=two_body_rme_lsjt.end(); ++it)
	    {
	      HalfInt S0,T0,S,T,Sp,Tp,Jp,J,J1p,J2p,J1,J2;
	      int N1,N2,N1p,N2p,L1,L2,L1p,L2p,L,Lp,L0;
	      TwoBodyStateLabelsLSJT bra,ket;
	      std::tie(bra,ket)=it->first;
	      std::tie(N1,L1,N2,L2,L,S,J,T)=ket;
	      std::tie(N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp)=bra;
	      double rme=it->second;
	      // std::cout<<fmt::format("[{} {} {} {}] {} {} {} {} | |[{} {} {} {}] {} {} {} {}     {}",
	      //   N1p,L1p,N2p,L2p,Lp,Sp,Jp,Tp,N1,L1,N2,L2,L,S,J,T,rme)<<std::endl;

	      for(J1p=abs(L1p-HalfInt(1,2)); J1p<=(L1p+HalfInt(1,2)); ++J1p)
	        for(J2p=abs(L2p-HalfInt(1,2)); J2p<=(L2p+HalfInt(1,2)); ++J2p)
	          for(J1=abs(L1-HalfInt(1,2)); J1<=(L1+HalfInt(1,2)); ++J1)
	            for(J2=abs(L2-HalfInt(1,2)); J2<=(L2+HalfInt(1,2)); ++J2)
	              {
	                double norm_factor=1.0;
	                if((N1==N2)&&(L1==L2)&&(J1==J2))
	                  norm_factor*=1/sqrt(2);
	                if((N1p==N2p)&&(L1p==L2p)&&(J1p==J2p))
	                  norm_factor*=1/sqrt(2);

	                if(not am::AllowedTriangle(J1,J2,J))
	                  continue;
	                if(not am::AllowedTriangle(J1p,J2p,Jp))
	                  continue;
	                double coefJp=am::Unitary9J(L1p,HalfInt(1,2),J1p,L2p,HalfInt(1,2),J2p,Lp,Sp,Jp);
	                double coefJ=am::Unitary9J(L1,HalfInt(1,2),J1,L2,HalfInt(1,2),J2,L,S,J);

	                // std::cout<<J1p<<" "<<J2p<<J1<<" "<<J2<<"    "<<rme<<"  "<<norm_factor<<"  "<<coefJp<<"  "<<coefJ<<std::endl;

	                TwoBodyStateLabelsJJJT braJ(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp);
	                TwoBodyStateLabelsJJJT ketJ(N1,L1,J1,N2,L2,J2,J,T);
	                TwoBodyBraketJJJT braket(braJ, ketJ);
	                two_body_rme_jjjt[braket]+=rme*coefJp*coefJ*norm_factor;
	              }
	    }
	}

	void BranchTwoBodyU3STToJJJT(int Jmax, int J0,
	  u3shell::IndexedTwoBodyTensorRMEsU3ST& indexed_two_body_rmes,
	  std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt
	  )
	{
	  std::map<TwoBodyBraketLST,double> two_body_rmes_lst;
	  BranchTwoBodyNLST(indexed_two_body_rmes,two_body_rmes_lst);

	  std::map<TwoBodyBraketLSJT,double> two_body_rme_lsjt;
	  BranchTwoBodyLSJT(Jmax, J0,two_body_rmes_lst,two_body_rme_lsjt);

	  branchJJJT(two_body_rme_lsjt, two_body_rme_jjjt);
	}

	void PrintTwoBodyMatrixElementsJJJT(const std::map<TwoBodyBraketJJJT,double>& two_body_rme_jjjt)
	{
	  for(auto it=two_body_rme_jjjt.begin(); it!=two_body_rme_jjjt.end(); ++it)
	    {
	      TwoBodyStateLabelsJJJT bra,ket;
	      std::tie(bra,ket)=it->first;
	      if(bra>ket)
	        continue;
	      double rme=it->second;
	      int N1p,N2p,N1,N2,L1,L2,L1p,L2p;
	      HalfInt J1,J2,J1p,J2p, Jp,J,Tp,T;
	      std::tie(N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp)=bra;
	      std::tie(N1,L1,J1,N2,L2,J2,J,T)=ket;
	      if(fabs(rme)>zero_threshold)
	        std::cout<<fmt::format("{} {} {} {} {} {} {} {}   {} {} {} {} {} {} {} {}   {}",
	          N1p,L1p,J1p,N2p,L2p,J2p,Jp,Tp,N1,L1,J1,N2,L2,J2,J,T,rme)<<std::endl;
	    }
	}



}
#ifndef CNULLSPACESOLVER_H
#define CNULLSPACESOLVER_H
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/QR>

#include <vector>

//	This class finds solution to a system of (sometimes) undetermined linear
//	equations 
//		Czx * HWS = 0
//		Cxy * HWS = 0
//		Sp  * HWS = 0
//	It uses GSL and Eigen numerical libraries.
//  The methods ordered by their speed:
//	(1) QR decomposition with column pivoting [Eigen]					ColPivotQR_Eigen 
//	(2) QR decomposition with full pivoting [Eigen]						FullPivotQR_Eigen
//	(3) QR decomposition with full pivoting [GSL]						QR_GSL 
//	(4) SVD decomposition (if m>=n) [Eigen] & method (1) in case m<n	SVD_ColPivotQR_Eigen
//	(5) SVD decomposition (if m>=n) [Eigen] & method (2) in case m<n	SVD_FullPivotQR_Eigen
//	(6) SVD decomposition (if m>=n) [GSL] 	& method (3) in case m<n	SVD_QR_GSL
//
//	NOTE:	static DOUBLE PRECISION_LIMIT = std::numeric_limits<float>::epsilon();
//	Warning is issued in case that |max(A * HWS)| > PRECISION_LIMIT
class CNullSpaceSolver {
	public:
	enum SolverType {ColPivotQR_Eigen, FullPivotQR_Eigen, QR_GSL, SVD_ColPivotQR_Eigen, SVD_FullPivotQR_Eigen, SVD_QR_GSL};
	typedef double DOUBLE;
	typedef std::vector<std::pair<std::pair<size_t, size_t>, DOUBLE> > MATRIX;
	typedef std::vector<std::vector<DOUBLE> > HWS;
	private:
	static DOUBLE PRECISION_LIMIT;
	private:
	struct CmpMatrixElems { // This structure enables to sort matrix elements of MATRIX
		public:
		bool operator() (const std::pair<std::pair<size_t, size_t>, DOUBLE>& s1, const std::pair<std::pair<size_t, size_t>, DOUBLE>& s2)
		{
			return (s1.first.second < s2.first.second); // order by row number only, columns won't be sorted
		}
	};
	private:
	SolverType m_Solver;
	MATRIX m_Matrix;
	HWS m_HWSVectors;
	public:
	CNullSpaceSolver(SolverType Solver): m_Solver(Solver) {};
	public:
//	KEY METHOD #1	
	void CalculateNullSpace(const size_t nrows, const size_t ncols, const size_t Mult);
//	KEY METHOD #2	
	const HWS& GetNullSpace() const {return m_HWSVectors;};
	public:
//	These two methods were used during debugging CHwsGen and will be removed in
//	the future.
	void ShowMatrix_matlab(const size_t ncols, const size_t nrows);
	void ShowMatrix_C(const size_t ncols, const size_t nrows);
	public:
	void SaveMatrix(const unsigned n, const unsigned A, const unsigned S2, const unsigned Mult, const unsigned lm, const unsigned mu);
	public:
	inline void SetMatrix(const size_t irow, const size_t icol, const DOUBLE dCoeff) 
	{m_Matrix.push_back(std::make_pair(std::make_pair(icol, irow), dCoeff));};
	inline size_t MatrixSize() const {return m_Matrix.size();};
	inline void Reset() {m_Matrix.resize(0); m_HWSVectors.clear();}; // Question: Do I have to "clear()" m_HWSVectors ?
	private:
//	Drivers for numerical subroutines provided by Eigen library	
	void CalculateNullSpaceSVD_Eigen(const size_t nrows, const size_t ncols, const size_t Mult);
	void CalculateNullSpaceQR_Eigen(const size_t nrows, const size_t ncols, const size_t Mult);

//	Drivers for numerical subroutines provided by GSL library	
	void CalculateNullSpaceSVD_GSL(const size_t nrows, const size_t ncols, const size_t Mult);
	void CalculateNullSpaceQR_GSL(const size_t nrows, const size_t ncols, const size_t Mult);
};

#endif

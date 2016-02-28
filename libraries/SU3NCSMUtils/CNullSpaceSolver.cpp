#include <SU3NCSMUtils/CNullSpaceSolver.h>
#include <iostream>
#include <limits>
#include <sstream>
#include <fstream>

//	Warning is issued in case that |max(A*HWS)| > PRECISION_LIMIT
//	PRECISION_LIMIT = std::numeric_limits<float>::epsilon()
CNullSpaceSolver::DOUBLE CNullSpaceSolver::PRECISION_LIMIT = std::numeric_limits<float>::epsilon();

// Store matrix in text file with simple "row \t column \t coeff" structure.
// Obviously, this is not very "size" efficient format.
// TODO: implement storing in binary file with more efficient structure.
// Such as, e.g., 
// [row1, {column1 coeff1} {column2, coeff2} {column3, coeff3}] 
// [row2, {column1 coeff1} {column4, coeff4}] 
// [row3, {column5 coeff5} {column7, coeff7} {column13, coeff13} {column14, coeff14}]; 
void CNullSpaceSolver::SaveMatrix(	const unsigned n, 
									const unsigned A, 
									const unsigned S2, 
									const unsigned Mult, 
									const unsigned lm, 
									const unsigned mu)
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	// Generate name of the file where hws will be stored. It has this
	// structure: "n#A#S#Mult#lm#mu#.mat". Note that number following S corresponds to the
	// value of spin multiplied by two, i.e. 2*S 
	std::ostringstream sn;
	std::ostringstream sA;
	std::ostringstream sS2;
	std::ostringstream sMult;
	std::ostringstream slm;
	std::ostringstream smu;
	std::string sFileName;
	sn  << n;
	sA  << A;
	sS2 << S2;
	sMult << Mult;
	slm << lm;
	smu << mu;
	sFileName = 'n' + sn.str() + 'A' + sA.str() + 'S' + sS2.str() + '_' + sMult.str() + "lm" + slm.str() + "mu" + smu.str() + ".mat";
	std::ofstream ofs(sFileName.c_str());
	if (!ofs) {
		std::cerr << "Not able to open file " << sFileName << "!" << std::endl;
		exit(EXIT_FAILURE);
	}
	ofs.precision(std::numeric_limits<DOUBLE>::digits10);
	ofs.setf(std::ios::scientific);
	///////////////////////////////////////////////////////////////////////////////////////////////////////
	
	size_t icol, irow;
	DOUBLE coeff;
	size_t nElems = m_Matrix.size();
	for (size_t i = 0; i < nElems; ++i)
	{
		icol = m_Matrix[i].first.first;
		irow = m_Matrix[i].first.second;
		coeff = m_Matrix[i].second;
		ofs << irow << "\t" << icol << "\t" << coeff << std::endl;
	}
}

void CNullSpaceSolver::CalculateNullSpaceSVD_Eigen(const size_t nrows, const size_t ncols, const size_t Mult)
{
// NOTE: SVD interface may change in future as it belongs to the experimental part of eigen
// Eigen does not implement SVD for nrows < ncols. 
	assert(nrows >= ncols); 

	size_t i, j;
	size_t icol, irow, iMult;
	DOUBLE dmax;
	size_t nElems = m_Matrix.size();

//	Create matrix A 
	Eigen::MatrixXd A((int)nrows, (int)ncols);
//	and fill it with zeros
	A.setZero();
	
// 	A <--- m_Matrix
	for (size_t i = 0; i < nElems; ++i)
	{
		icol = m_Matrix[i].first.first;
		irow = m_Matrix[i].first.second;
		A(irow, icol) = m_Matrix[i].second;
	}
//	perform SVD 
	Eigen::JacobiSVD<Eigen::MatrixXd> svdA(A,Eigen::ComputeFullV);
//  last Mult columns of V contains SU(3) HWS	
	const Eigen::MatrixXd& V = svdA.matrixV();

	Eigen::VectorXd null_vector((int)nrows); 
	for (j = (ncols - Mult), iMult = 0; j < ncols; ++j, ++iMult)
	{
		Eigen::VectorXd HWS((int)ncols); 
		for (i = 0; i < ncols; ++i) 
		{
			HWS(i) = V(i, j);
			m_HWSVectors[iMult][i] = HWS(i);
		}
		null_vector = A*HWS;
		dmax = null_vector.maxCoeff();
		std::cout << "max(A * HWS) = " << dmax << std::endl;
		if (fabs(dmax) > PRECISION_LIMIT) 
		{
			std::cerr << "SU(3)xSU(2) HWS are not annihilated by SU(3) and SU(2) raising operators!!!" << std::endl;
		} 
	}
}

void CNullSpaceSolver::CalculateNullSpaceQR_Eigen(const size_t nrows, const size_t ncols, const size_t Mult)
{
	size_t i, j;
	size_t icol, irow, iMult;
	DOUBLE dmax;
	size_t nElems = m_Matrix.size();
	size_t nrows_t = ncols;
	size_t ncols_t = nrows;

//	Create matrix A 
	Eigen::MatrixXd At((int) nrows_t, (int)ncols_t);
	At.setZero();
	
// 	A <--- m_Matrix
	for (size_t i = 0; i < nElems; ++i)
	{
		icol = m_Matrix[i].first.first;
		irow = m_Matrix[i].first.second;
		At(icol, irow) = m_Matrix[i].second; // Notice icol <--> irow
	}

	// Eigen 2 devel: 
	//   Eigen::MatrixXd Q((m_Solver == FullPivotQR_Eigen) ? (At.fullPivotingHouseholderQr().matrixQ()) : (At.colPivotingHouseholderQr().matrixQ()));
	// Eigen3.0 beta4 mixed method names!
//	Eigen::MatrixXd Q((m_Solver == FullPivotQR_Eigen) ? (At.fullPivHouseholderQr().matrixQ()) : (At.colPivHouseholderQr().householderQ()));
	Eigen::MatrixXd Q;
	if (m_Solver == FullPivotQR_Eigen)
	{
		Q = At.fullPivHouseholderQr().matrixQ();
	}
	else
	{
		Q = At.colPivHouseholderQr().householderQ();
	}


	Eigen::VectorXd null_vector((int)ncols_t); 
	Eigen::VectorXd HWS((int)nrows_t); 
	At.transposeInPlace(); // in order to be able to calculate A*HWS we need to Transpose(At)
	for (j = (nrows_t - Mult), iMult = 0; j < nrows_t; ++j, ++iMult)
	{
		for (i = 0; i < nrows_t; ++i) {
			HWS(i) = Q(i, j);
			m_HWSVectors[iMult][i] = HWS(i);
		}
		null_vector = At*HWS; // This is Ok, since At now contains A
		dmax = null_vector.maxCoeff();
		std::cout << "max(A * HWS) = " << dmax << std::endl;
		if (fabs(dmax) > PRECISION_LIMIT) 
		{
			std::cerr << "SU(3)xSU(2) HWS are not annihilated by SU(3) and SU(2) raising operators!!!" << std::endl;
		} 
	}
}

void CNullSpaceSolver::CalculateNullSpaceSVD_GSL(const size_t nrows, const size_t ncols, const size_t Mult)
{
// GSL does not implement SVD for nrows < ncols. In theory, one could use
// SVD on transpose matrix and obtain null space basis in matrix U, but GSL
// performs "thin" SVD and hence columns in U that corresponds are not
// provided.
	assert(nrows >= ncols); 

	size_t i, j;
	size_t icol, irow, iMult;
	size_t nElems = m_Matrix.size();
	DOUBLE dmax;

	gsl_matrix *A = gsl_matrix_calloc (nrows, ncols); 
	gsl_matrix *V = gsl_matrix_calloc (ncols, ncols); 
	gsl_matrix *X = gsl_matrix_calloc (ncols, ncols); 
	gsl_vector *S = gsl_vector_calloc(ncols);
	gsl_vector *work = gsl_vector_calloc(ncols);
	// A = m_Matrix
	for (i = 0; i < nElems; ++i)
	{
		icol = m_Matrix[i].first.first;
		irow = m_Matrix[i].first.second;
		gsl_matrix_set(A, irow, icol, m_Matrix[i].second);
	}

//	gsl_linalg_SV_decomp (A, V, S, work);
	gsl_linalg_SV_decomp_mod (A, X, V, S, work);
		
	// Matrix A is being rewritten by matrix U and hence in order to check
	// SU(3)xSU(2) HWS (basis of null space of A) we need to construct it again
	// so we could evaluate A*HWS and check whether it is null vector
	gsl_matrix_set_zero(A);
	// A = m_Matrix
	for (i = 0; i < nElems; ++i)
	{
		icol = m_Matrix[i].first.first;
		irow = m_Matrix[i].first.second;
		gsl_matrix_set(A, irow, icol, m_Matrix[i].second);
	}
		
	std::cout << "The " << Mult << " lowest singular values: ";
	for (i = (ncols - Mult); i < ncols; ++i) {
		std::cout << gsl_vector_get(S, i) << "\t";
	}
	std::cout << std::endl;

	// Check whether A*HWS = null_vector
	gsl_vector *null_vector = gsl_vector_calloc(nrows);
	for (j = (ncols - Mult), iMult = 0; j < ncols; ++j, ++iMult) 
	{
		gsl_vector_const_view SU3Hws = gsl_matrix_const_column (V, j);
		for (i = 0; i < ncols; ++i)
		{
			m_HWSVectors[iMult][i] = gsl_vector_get(&SU3Hws.vector, i);
		}
		gsl_blas_dgemv(CblasNoTrans, 1.0, A, &SU3Hws.vector, 0.0, null_vector); // null_vector = A*SU3Hws
		dmax = gsl_vector_max(null_vector);
		std::cout << "max(A * HWS) = " << dmax << std::endl;
		if (fabs(dmax) > PRECISION_LIMIT) {
			std::cerr << "SU(3)xSU(2) HWS are not annihilated by SU(3) and SU(2) raising operators!!!" << std::endl;
		} 
	}
	gsl_vector_free (null_vector); 

	gsl_matrix_free (A); 
	gsl_matrix_free (V); 
	gsl_matrix_free (X); 
	gsl_vector_free (S); 
	gsl_vector_free (work); 
}

void CNullSpaceSolver::CalculateNullSpaceQR_GSL(const size_t nrows, const size_t ncols, const size_t Mult)
{
	size_t i, j;
	size_t icol, irow, iMult;
	size_t nElems = m_Matrix.size();
	DOUBLE dmax;
	size_t nrows_t = ncols;
	size_t ncols_t = nrows;

	gsl_matrix * At = gsl_matrix_calloc (nrows_t, ncols_t); 
	// At = Transpose(m_Matrix)
	for (i = 0; i < nElems; ++i)
	{
		icol = m_Matrix[i].first.first;
		irow = m_Matrix[i].first.second;
		gsl_matrix_set(At, icol, irow, m_Matrix[i].second); // Notice icol <--> irow
	}

	gsl_matrix *Q = gsl_matrix_calloc (nrows_t, nrows_t); 
	gsl_matrix *R = gsl_matrix_calloc (nrows_t, ncols_t);
	gsl_vector * tau = gsl_vector_calloc(std::min(ncols_t, nrows_t)); // min(nrowsT, ncolsT) 
	gsl_permutation *p = gsl_permutation_calloc(ncols_t); 
	int signum;
	gsl_vector *norm = gsl_vector_calloc(ncols_t);

	gsl_linalg_QRPT_decomp2 (At, Q, R, tau, p, &signum, norm);

	// Check whether A*HWS = null_vector
	gsl_vector *null_vector = gsl_vector_calloc(ncols_t);
	for (j = nrows_t-Mult, iMult = 0; j < nrows_t; ++j, ++iMult)
	{
		gsl_vector_const_view SU3Hws = gsl_matrix_const_column (Q, j);
		for (i = 0; i < nrows_t; ++i)
		{
			m_HWSVectors[iMult][i] = gsl_vector_get(&SU3Hws.vector, i);
		}
		// A = transpose(At)
		gsl_blas_dgemv(CblasTrans, 1.0, At, &SU3Hws.vector, 0.0, null_vector); // null_vector = transpose(At)*SU3Hws
		dmax = gsl_vector_max(null_vector);

		std::cout << "max(A * HWS) = " << dmax << std::endl;
		if (fabs(dmax) > PRECISION_LIMIT) 
		{
			std::cerr << "SU(3)xSU(2) HWS are not annihilated by SU(3) and SU(2) raising operators!!!" << std::endl;
		} 
	}
	gsl_vector_free(null_vector);

	gsl_matrix_free (At); 
	gsl_matrix_free (Q); 
	gsl_matrix_free (R); 
	gsl_vector_free(tau);
	gsl_vector_free(norm);
	gsl_permutation_free(p);
}

void CNullSpaceSolver::CalculateNullSpace(const size_t nrows, const size_t ncols, const size_t Mult)
{
	// Prepare memory for resulting vectors
	m_HWSVectors.resize(Mult);
	for(size_t i = 0; i < Mult; ++i) {
		m_HWSVectors[i].resize(ncols);
	}

	if (m_Solver <= FullPivotQR_Eigen) {
		CalculateNullSpaceQR_Eigen(nrows, ncols, Mult);
	}
	else if (m_Solver == QR_GSL) {
		CalculateNullSpaceQR_GSL(nrows, ncols, Mult);
	}
	else if (m_Solver == SVD_FullPivotQR_Eigen || m_Solver == SVD_ColPivotQR_Eigen) {
		if (nrows >= ncols) {
			CalculateNullSpaceSVD_Eigen(nrows, ncols, Mult);
		} else {
			CalculateNullSpaceQR_Eigen(nrows, ncols, Mult);
		}
	}
	else if (m_Solver == SVD_QR_GSL)
	{
		if (nrows >= ncols) {
			CalculateNullSpaceSVD_GSL(nrows, ncols, Mult);
		} else {
			CalculateNullSpaceQR_GSL(nrows, ncols, Mult);
		}
	}
}

void CNullSpaceSolver::ShowMatrix_matlab(const size_t ncols, const size_t nrows) 
{
	size_t nElems = m_Matrix.size();
	size_t irow, icol;

	std::sort(m_Matrix.begin(), m_Matrix.end(), CmpMatrixElems());
//	return;
/*
	if (!nElems || nElems > 1000) {
		return;
	}
*/	
	std::cout << "A = [" << std::endl;
	std::vector<DOUBLE> Row(ncols, 0.0);

	irow = m_Matrix[0].first.second;
	for (size_t i = 0; i < nElems; ++i)
	{
		icol = m_Matrix[i].first.first;
		Row[icol] = m_Matrix[i].second;

		if (i < (nElems - 1) && irow == m_Matrix[i+1].first.second) {
			continue;
		}
		else if (i < (nElems - 1) && irow < m_Matrix[i+1].first.second)
		{
			for (size_t j = 0; j < ncols; ++j) 
			{
				std::cout << Row[j] << " ";
			}
			std::cout << ";" << std::endl;
			std::fill_n(Row.begin(), ncols, 0.0);
			irow++;
		} 
		else 
		{
			for (size_t j = 0; j < ncols; ++j) 
			{
				std::cout << Row[j] << " ";
			}
			std::cout << "]" << std::endl << std::endl;
		}
	}
}

void CNullSpaceSolver::ShowMatrix_C(const size_t ncols, const size_t nrows)
{
	size_t nElems = m_Matrix.size();
	size_t irow, icol;

	std::sort(m_Matrix.begin(), m_Matrix.end(), CmpMatrixElems());
//	return;
/*
	if (!nElems || nElems > 1000) {
		return;
	}
*/	
	std::cout << "const int nrows = " << nrows << ";" << std::endl << "const int ncols = " << ncols << ";" << std::endl;
	std::cout << "double Data[" << nrows << "][" << ncols << "] = {" << std::endl;
	std::vector<DOUBLE> Row(ncols, 0.0);

	irow = m_Matrix[0].first.second;
	for (size_t i = 0; i < nElems; ++i)
	{
		icol = m_Matrix[i].first.first;
		Row[icol] = m_Matrix[i].second;

		if (i < (nElems - 1) && irow == m_Matrix[i+1].first.second) {
			continue;
		}
		else if (i < (nElems - 1) && irow < m_Matrix[i+1].first.second)
		{

			std::cout << "{";
			for (size_t j = 0; j < ncols-1; ++j) 
			{
				std::cout << Row[j] << ", ";
			}
			std::cout << Row[ncols-1] << "}," << std::endl;
			std::fill_n(Row.begin(), ncols, 0.0);
			irow++;
		} 
		else 
		{
			std::cout << "{";
			for (size_t j = 0; j < ncols-1; ++j) 
			{
				std::cout << Row[j] << ",";
			}
			std::cout << Row[ncols-1] << "}};" << std::endl << std::endl;
		}
	}				
}

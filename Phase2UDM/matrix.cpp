//----------------------------------------------------------------------
//
// Copyright by Rocscience.  All rights reserved.
//
//----------------------------------------------------------------------

#include "matrix.h"
#include <math.h> //for sqrt()

///#include "vs2008vc6_std_defines.h"
//#include "stdafx.h"

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))

Matrix::Matrix()
{
	v_alloc(0, 0);
}
void Matrix::sizeError(char* where, const Matrix& m, size_t a, size_t b)
{
	//	_STD_CERR << where << ": "
	//		<< m.nRow() << "x" << m.nCol() << " Matrix expected, "
	//		<< a << "x" << b << " found" << _STD_ENDL;
	abort();
}

void Matrix::v_alloc(size_t a, size_t b)
{
	//    _STD_CERR << a << b << _STD_ENDL;
	nrow = a; ncol = b;
	if (nrow > MAX_SIZE || ncol > MAX_SIZE) {
		_v = new double* [a];    // check alloc succeeded !!
		_v[0] = new double[a * b];
		for (size_t i = 1; i < a; ++i)
			_v[i] = _v[i - 1] + b;
	}
	else {
		_v = &vPtr[0];
		_v[0] = &vData[0];
		for (size_t i = 1; i < a; ++i)
			_v[i] = _v[i - 1] + b;
	}
}

void Matrix::initialize(size_t nr, size_t nc, double* f)
{
	v_alloc(nr, nc);
	size_t i, j;
	if (f == 0) { // initialize 0 matrix
		memset(_v[0], 0, sizeof(double) * nrow * ncol);
	}
	else if (nrow == 1) { // initialize row
		for (j = 0; j < ncol; j++) _v[0][j] = f[j];
	}
	else if (ncol == 1) { // initialize column
		for (i = 0; i < nrow; i++) _v[i][0] = f[i];
	}
	else { // initialize matrix
		for (i = 0; i < nrow; i++)
			for (j = 0; j < ncol; j++)
				_v[i][j] = f[i * ncol + j];
	}
}

Matrix::Matrix(int nr, int nc, double* f)
{
	// Check that nr, nc >0
	// allow zero matrix, nr=nc=0?
///	if ( nr<0 || nc<0 )
//		_STD_CERR << "trying to initialize a Matrix of dimensions (" 
//		<< nr << ", " << nc <<")!" << _STD_ENDL;

	size_t nr_safe = size_t(nr);
	size_t nc_safe = size_t(nc);

	initialize(nr_safe, nc_safe, f);
}
/*
Matrix::Matrix(const CMatrix&CM)
{
	v_alloc(CM.GetSizeM(), CM.GetSizeN());
	for (size_t i=0; i<nrow; i++)
		for (size_t j=0; j<ncol; j++)
			_v[i][j] = CM.GetElement(i,j);
}
*/
Matrix::~Matrix()
{
	if (nrow > MAX_SIZE || ncol > MAX_SIZE) {
		delete[] _v[0];
		delete[] _v;
	}
}

void Matrix::resize(size_t irow, size_t icol)
{
	initialize(irow, icol);
	for (size_t i = 0; i < nrow; i++)
		for (size_t j = 0; j < ncol; j++)
			_v[i][j] = 0.0;
}

void Matrix::operator=(const Matrix& m)
{
	if (!sameSize(m.nrow, m.ncol))
		Matrix::sizeError("operator=", *this, m.nrow, m.ncol);
	memcpy(_v[0], m(0), sizeof(double) * ncol * nrow);
}

void Matrix::operator=(double val)
{
	for (size_t i = 0; i < nrow; i++)
		for (size_t j = 0; j < ncol; j++)
			_v[i][j] = val;
}

int Matrix::operator==(const Matrix& m) const
{
	if (!sameSize(m.nrow, m.ncol)) return 0;
	for (size_t i = 0; i < nrow; i++)
		for (size_t j = 0; j < ncol; j++)
			if (_v[i][j] != m(i, j)) return 0;
	return 1;
}

void Matrix::operator+=(const Matrix& m)
{
	if (!sameSize(m.nrow, m.ncol))
		Matrix::sizeError("operator+=", *this, m.nRow(), m.nCol());
	for (size_t i = 0; i < m.nRow(); i++)
		for (size_t j = 0; j < m.nCol(); j++)
			_v[i][j] += m(i, j);
}

Matrix operator+(const Matrix& m, const Matrix& n)
{
	if (!m.sameSize(n.nRow(), n.nCol()))
		Matrix::sizeError("operator+", m, n.nRow(), n.nCol());
	// C++2.0 bug
	//    Matrix rm(m.nRow(),m.nCol());
	size_t nr = m.nRow();
	size_t nc = m.nCol();
	Matrix rm((int)nr, (int)nc);
	for (size_t i = 0; i < m.nRow(); i++)
		for (size_t j = 0; j < m.nCol(); j++)
			rm(i, j) = m(i, j) + n(i, j);
	return rm;
}

Matrix operator-(const Matrix& m, const Matrix& n)
{
	if (!m.sameSize(n.nRow(), n.nCol()))
		Matrix::sizeError("operator-", m, n.nRow(), n.nCol());
	// C++2.0 bug
	//    Matrix rm(m.nRow(),m.nCol());
	size_t nr = m.nRow();
	size_t nc = m.nCol();
	Matrix rm((int)nr, (int)nc);
	for (size_t i = 0; i < m.nRow(); i++)
		for (size_t j = 0; j < m.nCol(); j++)
			rm(i, j) = m(i, j) - n(i, j);
	return rm;
}

Matrix operator-(const Matrix& m)
{
	// C++2.0 bug
	//    Matrix rm(m.nRow(),m.nCol());
	size_t nr = m.nRow();
	size_t nc = m.nCol();
	Matrix rm((int)nr, (int)nc);
	for (size_t i = 0; i < m.nRow(); i++)
		for (size_t j = 0; j < m.nCol(); j++)
			rm(i, j) = -m(i, j);
	return rm;
}

double det(const Matrix& m)
{
	if (m.nRow() != m.nCol()) {
		///		_STD_CERR << "det: not a square matrix" << _STD_ENDL;;
		abort();
	}
	if (m.nRow() == size_t(1)) return m(0, 0);
	if (m.nRow() == size_t(2))
		return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
	if (m.nRow() == size_t(3))
		return (m(0, 0) * m(1, 1) * m(2, 2)
			- m(0, 1) * m(1, 2) * m(2, 0)
			+ m(0, 2) * m(2, 1) * m(0, 1));
	if (m.isUpperTriangle()) {
		double val = 1;
		for (size_t i = 0; i < m.nRow(); i++) val *= m(i, i);
		return val;
	}
	double val = 0;
	int sign = 1;
	for (size_t j = 0; j < m.nCol(); j++) {
		val += sign * m(0, j) * det(m.coFactor(0, j));
		sign *= -1;
	}
	return val;
}

double norm(const Matrix& m)
{
	double val = 0;
	for (size_t i = 0; i < m.nRow(); i++)
		for (size_t j = 0; j < m.nCol(); j++)
			if (m(i, j) > val) val = m(i, j);
	return val;
}

double norm2(const Matrix& m)
{
	double val = 0;
	for (size_t i = 0; i < m.nRow(); i++)
		for (size_t j = 0; j < m.nCol(); j++)
			val += m(i, j) * m(i, j);
	return (val > 1e-32 ? sqrt(val) : 0.0);
}

Matrix operator*(const Matrix& m, const Matrix& n)
{
	if (m.nCol() != n.nRow()) {
		//		_STD_CERR << "operator*: " << m.nCol() <<
		//			"x* Matrix expected "
		//			<< n.nRow() << "x* found." << _STD_ENDL;
		abort();
	}
	// C++2.0 bug
	//    Matrix rm(m.nRow(),n.nCol());
	size_t nr = m.nRow();
	size_t nc = n.nCol();
	Matrix rm((int)nr, (int)nc);
	for (size_t i = 0; i < rm.nRow(); i++) {
		for (size_t j = 0; j < rm.nCol(); j++) {
			double dSubt = 0.0;
			for (size_t k = 0; k < m.nCol(); k++)
				dSubt += m(i, k) * n(k, j);
			rm(i, j) = dSubt;
		}
	}
	return rm;
}

void outer_product(const Matrix& A, const Matrix& B, Matrix& C)
{
	int i, j, k, l, m, n;

	for (i = 0; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					m = 3 * i + k;
					n = 3 * j + l;
					C(m, n) = A(i, j) * B(k, l);
				}

}

Matrix Matrix::t() const
{
	Matrix rm((int)ncol, (int)nrow);
	for (size_t i = 0; i < rm.nrow; i++)
		for (size_t j = 0; j < rm.ncol; j++)
			rm(i, j) = _v[j][i];   // hm ... literal transpose
	return rm;
}

Matrix operator*(double f, const Matrix& m)
{
	// C++2.0 bug
	//    Matrix rm(m.nRow(),m.nCol());
	size_t nr = m.nRow();
	size_t nc = m.nCol();
	Matrix rm((int)nr, (int)nc);
	for (size_t i = 0; i < rm.nRow(); i++)
		for (size_t j = 0; j < rm.nCol(); j++)
			rm(i, j) = f * m(i, j);
	return rm;
}

void Matrix::operator*=(double f)
{
	for (size_t i = 0; i < nrow; i++)
		for (size_t j = 0; j < ncol; j++)
			_v[i][j] *= f;
}

int Matrix::isUpperTriangle() const
{
	for (size_t j = 0; j < ncol; j++)
		for (size_t i = j + 1; i < nrow; i++)
			if (_v[i][j] != 0)
				return 0;
	return 1;
}


Matrix Matrix::coFactor(size_t irow, size_t jcol) const
{
	if (nrow == 1 || ncol == 1) {
		//		_STD_CERR << "coFactor: can't coFactor row or column matrix" <<
		//			_STD_ENDL;
		abort();
	}
	Matrix val((int)(nrow - 1), (int)(ncol - 1));
	size_t getcol, getrow = 0;
	for (size_t i = 0; i < val.nRow(); i++) {
		if (getrow == irow) ++getrow;
		if (getrow == nrow) break;
		getcol = 0;
		for (size_t j = 0; j < val.nCol(); j++) {
			if (getcol == jcol) ++getcol;
			if (getcol == ncol) continue;
			val(i, j) = _v[getrow][getcol];
			++getcol;
		}
		++getrow;
	}
	return val;
}

Matrix Matrix::operator~() const
{
	if (nrow != ncol)
	{
		//		_STD_CERR << "operator~: can't invert a non-square matrix" <<	_STD_ENDL;
		abort();
	}
	if (nrow == 1)
	{
		Matrix T(1, 1);
		T(0, 0) = size_t(1) / _v[0][0];
		return T;
	}
	if (nrow == 2)
	{
		double D = det(*this);
		Matrix T(2, 2);
		T(0, 0) = _v[1][1] / D;
		T(1, 1) = _v[0][0] / D;
		T(0, 1) = -_v[0][1] / D;
		T(1, 0) = -_v[1][0] / D;
		return T;
	}
	if (nrow > 2)
	{
		//		_STD_CERR << "operator~: can't invert a matrix > 2x2" << _STD_ENDL;
		abort();
	}

	return Matrix(0, 0);

}

Matrix Matrix::getInverse() const
{
	if (nrow != ncol) {
		//	_STD_CERR << "operator~: can't invert a non-square matrix" << _STD_ENDL;
		abort();
	}

	// Set result to the matrix so we don't overwrite the matrix
	Matrix result(*this);

	// Gauss elimination from Numerical Recipes in C++, Press et al., 2003
	int icol = -1;
	int irow = -1;
	size_t i, j, k, ll, l;
	double big, dum, pivinv;

	size_t n = nrow;

	// integer arrays used for bookkeeping
	size_t* indxc = new size_t[n];
	size_t* indxr = new size_t[n];
	size_t* ipiv = new size_t[n];
	for (i = 0; i < n; i++)   // should size_t i = 0 or size_t i = size_t(0)
		ipiv[i] = 0;

	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++) {
			if (ipiv[j] != 1) {
				for (k = 0; k < n; k++) {
					if (ipiv[k] == 0) {
						if (fabs(result(j, k)) >= big) {
							big = fabs(result(j, k));
							irow = (int)j;
							icol = (int)k;
						}
					}
				}
			}
		}
		++(ipiv[icol]);  // can this be unsafe

		// we now have pivot element, so interchange rows
		// columns are not physically interchanged, only relabeled using integer arrays
		// for bookkeeping.  Inverse matrix will be scrambled
		if (irow != icol) {
			for (l = 0; l < n; l++) {
				//SWAP(result[irow][l],result[icol][l]);
				double dum = result(irow, l);
				result(irow, l) = result(icol, l);
				result(icol, l) = dum;
			}

		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (result(icol, icol) == 0.0) {
			//			_STD_CERR << "Matrix::getBigInverse() - Trying to invert singular matrix! " << _STD_ENDL;
			abort();
		}
		pivinv = 1.0 / result(icol, icol);
		result(icol, icol) = 1.0;
		for (l = 0; l < n; l++)
			result(icol, l) *= pivinv;
		for (ll = 0; ll < n; ll++) {
			if (static_cast<int>(ll) != icol) {
				dum = result(ll, icol);
				result(ll, icol) = 0.0;
				for (l = 0; l < n; l++)
					result(ll, l) -= result(icol, l) * dum;

			}
		}
	}

	// now unscramble matrix
	int sl;
	for (sl = n - 1; sl >= 0; sl--) {
		if (indxr[sl] != indxc[sl])
			for (k = 0; k < n; k++) {
				//SWAP(result[k][indxr[l]],result[k][indxc[l]]);
				double dum = result(k, indxr[sl]);
				result(k, indxr[sl]) = result(k, indxc[sl]);
				result(k, indxc[sl]) = dum;
			}
	}

	delete[]indxc;
	delete[]indxr;
	delete[]ipiv;

	return result;
}


double* Matrix::GetDataArray()
{
	return &vData[0];
}

const void Matrix::writeMatrix(std::string fileName, int i, int j, int k)
{
	/*
	char buffer[3];
	_itoa_s(i, buffer, 3, 10);
	fileName.append(buffer);
	_itoa_s(j, buffer, 3, 10);
	fileName.append(buffer);
	_itoa_s(k, buffer, 3, 10);
	fileName.append(buffer);
	std::ofstream lf(fileName.c_str());


	for (size_t i = 0; i<nrow; ++i){
		for (size_t j = 0; j<ncol; ++j){
			assert(abs(_v[i][j])<1e10);
			if (abs(_v[i][j])<1e-12) lf << " " << 0;
			else lf << " " << _v[i][j];
		}
		lf << std::endl;
	}
	lf.close();
	*/
}





void Matrix::IdentityMatrix(size_t iRow, size_t iColumn, double dElement)
{

}

// MINE
void Vector::sizeError(char* where, const Vector& v, size_t a)
{
	//	_STD_CERR << where << ": "
	//		<< v.size() << " Vector expected, "
	//		<< a << " found" << _STD_ENDL;
	abort();
}

void Vector::v_alloc(size_t a)
{
	theSize = a;
	if (theSize > MAX_SIZE)
	{
		_v = new double[a];
	}
	else
	{
		_v = &vData[0];
	}
}

Vector::Vector(size_t n)
{
	v_alloc(n);
	memset(_v, 0, sizeof(double) * n);
}

Vector::~Vector()
{
	if (theSize > MAX_SIZE)
		delete[] _v;
}

void Vector::resize(size_t newSize)
{
	if (!sameSize(newSize))
	{
		if (theSize > MAX_SIZE)
			delete[] _v;
		v_alloc(newSize);
		memset(_v, 0, sizeof(double) * newSize);
	}
}

void Vector::resize(size_t newSize, double val)
{
	if (!sameSize(newSize))
	{
		if (theSize > MAX_SIZE)
			delete[] _v;
		v_alloc(newSize);
		memset(_v, 0, sizeof(double) * newSize);
	}
	for (size_t i = 0; i < theSize; i++)
		_v[i] = val;
}


void Vector::operator=(const Vector& v)
{
	if (!sameSize(v.theSize))
		resize(v.theSize);
	memcpy(_v, v._v, sizeof(double) * theSize);
}

void Vector::operator+=(const Vector& v)
{
	if (!sameSize(v.theSize))
		Vector::sizeError("operator=", *this, v.theSize);
	for (size_t i = 0; i < theSize; i++)
		_v[i] += v._v[i];
}

void Vector::operator-=(const Vector& v)
{
	if (!sameSize(v.theSize))
		Vector::sizeError("operator=", *this, v.theSize);
	for (size_t i = 0; i < theSize; i++)
		_v[i] -= v._v[i];
}

void Vector::operator=(double val)
{
	for (size_t i = 0; i < theSize; i++)
		_v[i] = val;
}

int Vector::operator==(const Vector& v) const
{
	if (!sameSize(v.theSize)) return 0;
	for (size_t i = 0; i < theSize; i++)
		if (_v[i] != v._v[i]) return 0;
	return 1;
}

Vector operator+(const Vector& m, const Vector& n)
{
	//K.D Oct 18, 2013
	//if (! m.sameSize(n.size()))
	//Vector::sizeError("operator+", m, n.size());
	size_t ns = min(m.size(), n.size());
	Vector rm(ns);
	for (size_t i = 0; i < ns; i++)
		rm[i] = m[i] + n[i];
	return rm;
}

Vector operator-(const Vector& m, const Vector& n)
{
	//K.D Oct 19,2013
	//if (! m.sameSize(n.size()))
	//Vector::sizeError("operator-", m, n.size());

	size_t ns = min(m.size(), n.size());
	Vector rm(ns);
	for (size_t i = 0; i < ns; i++)
		rm[i] = m[i] - n[i];
	return rm;
}


Vector operator/(const Vector& m, const Vector& n)
{
	if (!m.sameSize(n.size()))
		Vector::sizeError("operator-", m, n.size());
	size_t ns = m.size();
	Vector rm(ns);
	for (size_t i = 0; i < ns; i++)
		rm[i] = m[i] / n[i];
	return rm;
}


Vector operator-(const Vector& m)
{
	size_t ns = m.size();
	Vector rm(ns);
	for (size_t i = 0; i < ns; i++)
		rm[i] = -m[i];
	return rm;
}

double norm2(const Vector& m)
{
	double val = 0;
	size_t sz = m.size();
	for (size_t i = 0; i < sz; i++)
		val += m[i] * m[i];
	return (val > 1e-32 ? sqrt(val) : 0.0);
}

Vector square(const Vector& v)
{
	size_t ns = v.size();
	Vector rv(ns);
	for (size_t i = 0; i < ns; i++)
		rv[i] = v[i] * v[i];
	return rv;
}

double Vector::operator*(const Vector& v)
{
	if (!sameSize(v.size()))
		Vector::sizeError("operator+", *this, v.size());
	double val = 0.0;
	for (size_t i = 0; i < size(); i++)
		val += _v[i] * v._v[i];
	return val;
}

Vector Vector::operator*(double val)
{
	Vector newVec(*this);
	for (size_t i = 0; i < size(); i++)
		newVec[i] *= val;

	return newVec;
}

Vector operator*(double f, const Vector& m)
{
	size_t ns = m.size();
	Vector rm(ns);
	for (size_t i = 0; i < ns; i++)
		rm[i] = f * m[i];
	return rm;
}

Vector operator*(const Vector& v, const Matrix& m)
{
	if (v.size() != m.nRow()) {
		//		_STD_CERR << "operator*: " << v.size() <<
		//			"x* Vector expected "
		//			<< m.nRow() << "x* found." << _STD_ENDL;
		abort();
	}
	size_t nr = m.nRow(), nc = m.nCol();
	Vector rm(nc);
	for (size_t j = 0; j < nc; ++j)
		for (size_t i = 0; i < nr; i++)
			rm[j] += v[i] * m(i, j);
	return rm;
}

Vector operator*(const Matrix& m, const Vector& v)
{
	if (v.size() != m.nCol()) {
		//		_STD_CERR << "operator*: " << v.size() <<
		//			"x* Vector expected "
		//			<< m.nCol() << "x* found." << _STD_ENDL;
		abort();
	}
	size_t nr = m.nRow(), nc = m.nCol();
	Vector rm(nr);
	for (size_t j = 0; j < nc; ++j)
		for (size_t i = 0; i < nr; i++)
			rm[i] += v[j] * m(i, j);
	return rm;
}

void Vector::operator*=(double f)
{
	for (size_t i = 0; i < theSize; i++)
		_v[i] *= f;
}


double Vector::dotProduct(const Vector& v)
{
	if (!(v.theSize == theSize))
		return 0.0;
	double sum = 0.0;
	for (size_t i = 0; i < theSize; i++)
		sum += _v[i] * v(i);

	return sum;
}

void Vector::writeVector(std::string fileName, int i, int j, int k)
{
	/*
	char buffer[3];
	_itoa_s(i, buffer, 3, 10);
	fileName.append(buffer);
	_itoa_s(j, buffer, 3, 10);
	fileName.append(buffer);
	_itoa_s(k, buffer, 3, 10);
	fileName.append(buffer);
	std::ofstream lf(fileName.c_str());
	lf.precision(16);
	for (size_t i = 0; i<theSize; ++i){
		if (abs(_v[i])<1e-12) lf << " " << 0 << std::endl;
		else lf << " " << _v[i] << std::endl;
	}
	lf.close();
	*/
}

Vector Vector::OuterProduct33(const Vector& v2)
{
	Vector V(3);
	V(0) = _v[1] * v2(2) - _v[2] * v2(1);
	V(1) = _v[2] * v2(0) - _v[0] * v2(2);
	V(2) = _v[0] * v2(1) - _v[1] * v2(0);

	return V;
}


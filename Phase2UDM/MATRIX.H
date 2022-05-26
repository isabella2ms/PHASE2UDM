//----------------------------------------------------------------------
//
// Copyright by Rocscience.  All rights reserved.
//
//----------------------------------------------------------------------

#ifndef _MATRIX_H
#define _MATRIX_H

#include "REALTYPE.H"
#include <iostream>
#include <cassert> // for debugging
#include <vector>

class Vector;

class Matrix {
public:

	//! Default constructor
	Matrix();
	Matrix(int nr, int nc, double* = 0);
	Matrix(const Matrix&);
	//Matrix(const CMatrix&);
	~Matrix();

	size_t nRow() const { return nrow; }
	size_t nCol() const { return ncol; }
	void resize(size_t irow, size_t icol);

	double& operator()(size_t irow, size_t icol) const;
	double* operator()(size_t irow) const;
	void operator=(const Matrix&);
	void operator=(double val); // MINE
	void operator+=(const Matrix& m); // MINE
	int operator==(const Matrix&) const;
	void operator*=(double);
	Matrix operator~() const; //  inverse of 0,1,2-matrix only

	int isUpperTriangle() const;
	Matrix coFactor(size_t, size_t) const;
	Matrix t() const;  // transpose

	// get inverse for matrix larger than 2x2, JFH
	Matrix getInverse() const;
	double* GetDataArray();

	int sameSize(size_t a, size_t b) const { return (nrow == a && ncol == b); }
	static void sizeError(char* where, const Matrix&, size_t, size_t);

	const void writeMatrix(std::string fileName, int i, int j, int k);

	void SetElement(size_t iRow, size_t iColumn, double dElement) { _v[iRow][iColumn] = dElement; }
	void SetAllValues(double v);
	// this is only for extracting the 2D Stiffness MAtrix from the 3D one 
	Matrix EMatix3D_to_2D(void);

protected:

	void v_alloc(size_t, size_t);
	void v_dealloc();
	void initialize(size_t, size_t, double* = 0);

	static const size_t MAX_SIZE = 16;
	size_t nrow, ncol;

	double* vPtr[MAX_SIZE];
	double vData[MAX_SIZE * MAX_SIZE];
	double** _v;



};

extern double det(const Matrix&); //determinant
extern double norm(const Matrix&); // 
extern double norm2(const Matrix&); // MINE
extern void outer_product(const Matrix& A, const Matrix& B, Matrix& C);

extern Matrix operator+(const Matrix&, const Matrix&);
extern Matrix operator-(const Matrix&, const Matrix&);
extern Matrix operator*(const Matrix&, const Matrix&);
extern Matrix operator-(const Matrix&);
extern Matrix operator*(double, const Matrix&); // scalar multiply

class Vector {
public:
	void writeVector(std::string fileName, int i, int j, int k);
	Vector(size_t n = 0);
	Vector(const Vector&);
	~Vector();

	size_t size() const { return theSize; }
	void resize(size_t newSize);
	void resize(size_t newSize, double val);
	void SetElement(size_t i, double dElement) { _v[i] = dElement; }

	double& operator()(size_t ndx) const;
	double& operator[](size_t p) const;
	double& operator[](size_t p);

	int operator==(const Vector&) const;
	void operator=(const Vector&);
	void operator=(const std::vector<Real>&);
	void operator=(double val);
	void operator*=(double);
	void operator+=(const Vector& v);
	void operator-=(const Vector& v);
	double operator*(const Vector& v);
	Vector operator*(double val);
	int sameSize(size_t a) const { return (theSize == a); }
	static void sizeError(char* where, const Vector&, size_t);
	double dotProduct(const Vector& v);
	Vector OuterProduct33(const Vector& v2);

	friend class CMatrixSparseCRS;
	friend class CMatrixSparseCRSSym;
	static const int MAX_SIZE = 16;
	double vData[MAX_SIZE];
protected:

	void v_alloc(size_t);



	size_t theSize;
	double* _v; // made private JFH 31 July 2008
};

Vector operator+(const Vector&, const Vector&);
Vector operator-(const Vector&, const Vector&);
Vector operator/(const Vector&, const Vector&);
Vector operator-(const Vector&);
Vector operator*(double, const Vector&);
Vector operator*(const Vector&, const Matrix&);
Vector operator*(const Matrix&, const Vector&);

double norm2(const Vector& v);
Vector square(const Vector& v); // returns a vector with elements squared

// inlined functions

inline Matrix::Matrix(const Matrix& m)
{
	v_alloc(m.nrow, m.ncol);
	*this = m;
}

inline double& Matrix::operator()(size_t irow, size_t icol) const
{

	assert(irow < nrow && icol < ncol);
	return _v[irow][icol];
}

inline double* Matrix::operator()(size_t irow) const
{

	assert(irow < nrow);
	return _v[irow];
}

inline Vector::Vector(const Vector& v)
{
	v_alloc(v.theSize);
	*this = v;
}

inline double& Vector::operator()(size_t i) const
{
	assert(i < theSize);
	return _v[i];
}

inline double& Vector::operator[] (size_t p)
{
	assert(p < theSize);
	return *(_v + p);
}

inline double& Vector::operator[] (size_t p) const
{
	assert(p < theSize);
	return *(_v + p);
}

#endif

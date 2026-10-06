/**
 * @file matrix.cxx
 *
 * @brief implementation of a matrix type
 *
 * copied from
https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Source-File/
 * debugged and extended
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
 @endverbatim
 *
 */

#ifndef __matrix_implementation
#define __matrix_implementation

#include "matrix.h"
#include <iostream>

using namespace std;

// Parameter Constructor
template <typename T>
matrix<T>::matrix(unsigned _rows, unsigned _cols, const T &_initial) {
  mat.resize(_rows);
  for (unsigned i = 0; i < mat.size(); i++) {
    mat[i].resize(_cols, _initial);
  }
  rows = _rows;
  cols = _cols;
}

// Copy Constructor
template <typename T> matrix<T>::matrix(const matrix<T> &rhs) {
  mat = rhs.mat;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

// (Virtual) Destructor
template <typename T> matrix<T>::~matrix() {}

// Assignment Operator
template <typename T> matrix<T> &matrix<T>::operator=(const matrix<T> &rhs) {
  if (&rhs == this)
    return *this;

  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();

  mat.resize(new_rows);
  for (unsigned i = 0; i < mat.size(); i++) {
    mat[i].resize(new_cols);
  }

  for (unsigned i = 0; i < new_rows; i++) {
    for (unsigned j = 0; j < new_cols; j++) {
      mat[i][j] = rhs(i, j);
    }
  }
  rows = new_rows;
  cols = new_cols;

  return *this;
}

// Addition of two matrices
template <typename T> matrix<T> matrix<T>::operator+(const matrix<T> &rhs) const {
  matrix result(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i, j) = this->mat[i][j] + rhs(i, j);
    }
  }

  return result;
}

// Cumulative addition of this matrix and another
template <typename T> matrix<T> &matrix<T>::operator+=(const matrix<T> &rhs) {
  unsigned nrows = rhs.get_rows();
  unsigned ncols = rhs.get_cols();

  for (unsigned i = 0; i < nrows; i++) {
    for (unsigned j = 0; j < ncols; j++) {
      this->mat[i][j] += rhs(i, j);
    }
  }

  return *this;
}

// Subtraction of this matrix and another
template <typename T> matrix<T> matrix<T>::operator-(const matrix<T> &rhs) const {
  unsigned nrows = rhs.get_rows();
  unsigned ncols = rhs.get_cols();
  matrix result(nrows, ncols, 0.0);

  for (unsigned i = 0; i < nrows; i++) {
    for (unsigned j = 0; j < ncols; j++) {
      result(i, j) = this->mat[i][j] - rhs(i, j);
    }
  }

  return result;
}

// Cumulative subtraction of this matrix and another
template <typename T> matrix<T> &matrix<T>::operator-=(const matrix<T> &rhs) {
  unsigned nrows = rhs.get_rows();
  unsigned ncols = rhs.get_cols();

  for (unsigned i = 0; i < nrows; i++) {
    for (unsigned j = 0; j < ncols; j++) {
      this->mat[i][j] -= rhs(i, j);
    }
  }

  return *this;
}

// Left multiplication of this matrix and another
template <typename T> matrix<T> matrix<T>::operator*(const matrix<T> &rhs) const {
  unsigned nrows = this->get_rows();
  unsigned ncols = rhs.get_cols();
  unsigned inner = this->get_cols(); // shared dimension (== rhs.get_rows())
  matrix result(nrows, ncols, 0.0);

  for (unsigned i = 0; i < nrows; i++) {
    for (unsigned j = 0; j < ncols; j++) {
      for (unsigned k = 0; k < inner; k++) {
        result(i, j) += this->mat[i][k] * rhs(k, j);
      }
    }
  }

  return result;
}

// Cumulative left multiplication of this matrix and another
template <typename T> matrix<T> &matrix<T>::operator*=(const matrix<T> &rhs) {
  matrix result = (*this) * rhs;
  (*this) = result;
  return *this;
}

// Calculate a transpose of this matrix
template <typename T> matrix<T> matrix<T>::transpose() const {
  matrix result(cols, rows, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(j, i) = this->mat[i][j];
    }
  }

  return result;
}

// Matrix/scalar addition
template <typename T> matrix<T> matrix<T>::operator+(const T &rhs) const {
  matrix result(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i, j) = this->mat[i][j] + rhs;
    }
  }

  return result;
}

// Matrix/scalar subtraction
template <typename T> matrix<T> matrix<T>::operator-(const T &rhs) const {
  matrix result(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i, j) = this->mat[i][j] - rhs;
    }
  }

  return result;
}

// Matrix/scalar multiplication
template <typename T> matrix<T> matrix<T>::operator*(const T &rhs) const {
  matrix result(rows, cols, 0.0);

  for (unsigned i = 0; i < this->rows; i++) {
    for (unsigned j = 0; j < this->cols; j++) {
      result(i, j) = this->mat[i][j] * rhs;
    }
  }

  return result;
}

// Matrix/scalar division
template <typename T> matrix<T> matrix<T>::operator/(const T &rhs) const {
  matrix result(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i, j) = this->mat[i][j] / rhs;
    }
  }
  return result;
}

// Multiply a matrix with a vector
template <typename T> vector<T> matrix<T>::operator*(const vector<T> &rhs) {
  vector<T> result(rows, 0.0);
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result[i] += this->mat[i][j] * rhs[j];
    }
  }
  return result;
}

// Obtain a vector of the diagonal elements
template <typename T> vector<T> matrix<T>::diag_vec() {
  vector<T> result(rows, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    result[i] = this->mat[i][i];
  }
  return result;
}

// Access the individual elements
template <typename T>
T &matrix<T>::operator()(const unsigned &row, const unsigned &col) {
  return this->mat[row][col];
}

// Access the individual elements (const)
template <typename T>
const T &matrix<T>::operator()(const unsigned &row, const unsigned &col) const {
  return this->mat[row][col];
}

// Get the number of rows of the matrix
template <typename T> unsigned matrix<T>::get_rows() const {
  return this->rows;
}

// Get the number of columns of the matrix
template <typename T> unsigned matrix<T>::get_cols() const {
  return this->cols;
}

#endif

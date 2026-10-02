/**
 * @file matrix.h
 *
 * @brief definition of a matrix type
 *
 * @author Stefan Schippers
 * @verbatim
 $Id: matrix.h 2039 2026-07-20 07:57:32Z iamp $
// SPDX-License-Identifier: MIT
 @endverbatim
 *
 */

#pragma once

#include <vector>

template <typename T> class matrix {

public:
  matrix(unsigned _rows, unsigned _cols, const T &_initial);
  matrix(const matrix<T> &rhs);
  virtual ~matrix();

  // Operator overloading, for "standard" mathematical matrix operations
  matrix<T> &operator=(const matrix<T> &rhs);

  // Matrix mathematical operations
  matrix<T> operator+(const matrix<T> &rhs) const;
  matrix<T> &operator+=(const matrix<T> &rhs);
  matrix<T> operator-(const matrix<T> &rhs) const;
  matrix<T> &operator-=(const matrix<T> &rhs);
  matrix<T> operator*(const matrix<T> &rhs) const;
  matrix<T> &operator*=(const matrix<T> &rhs);
  matrix<T> transpose() const;

  // Matrix/scalar operations
  matrix<T> operator+(const T &rhs) const;
  matrix<T> operator-(const T &rhs) const;
  matrix<T> operator*(const T &rhs) const;
  matrix<T> operator/(const T &rhs) const;

  // Matrix/vector operations
  std::vector<T> operator*(const std::vector<T> &rhs);
  std::vector<T> diag_vec();

  // Access the individual elements
  T &operator()(const unsigned &row, const unsigned &col);
  const T &operator()(const unsigned &row, const unsigned &col) const;

  // Access the row and column sizes
  unsigned get_rows() const;
  unsigned get_cols() const;

private:
  std::vector<std::vector<T>> mat;
  unsigned rows;
  unsigned cols;
};

#include "matrix.cxx"

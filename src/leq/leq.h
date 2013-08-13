#pragma once
#include <Eigen/Dense>

using namespace Eigen;

template <typename DerivedL, typename Derivedb, typename Derivedx>
void forw_sub(const MatrixBase<DerivedL>& L, const MatrixBase<Derivedb>& b,
	          MatrixBase<Derivedx> const & x);

template <typename DerivedU, typename Derivedb, typename Derivedx>
void back_sub(const MatrixBase<DerivedU>& U, const MatrixBase<Derivedb>& b,
	          MatrixBase<Derivedx> const & x);

template <typename Derived>
void lu_fac(const MatrixBase<Derived>& X,
	        MatrixBase<Derived> const & L, MatrixBase<Derived> const & U);

template <typename DerivedA, typename Derivedb, typename Derivedx>
void gauss_elim(const MatrixBase<DerivedA>& A, const MatrixBase<Derivedb>& b,
	            MatrixBase<Derivedx> const & x);

template <typename Derived>
void inv(const MatrixBase<Derived>& X,
	     MatrixBase<Derived> const & Xinv);

#include "leq.hpp"
#pragma once
#include <Eigen/Dense>

using namespace Eigen;

template <class DerivedM, class DerivedV>
void forw_sub(const MatrixBase<DerivedM>& L, const MatrixBase<DerivedV>& b,
	          MatrixBase<DerivedV> const & x);

template <typename DerivedM, typename DerivedV>
void back_sub(const MatrixBase<DerivedM>& U, const MatrixBase<DerivedV>& b,
	          MatrixBase<DerivedV> const & x);

template <typename Derived>
void lu_fac(const MatrixBase<Derived>& X,
	        MatrixBase<Derived> const & L, MatrixBase<Derived> const & U);

template <typename DerivedM, typename DerivedV>
void gauss_elim(const MatrixBase<DerivedM>& A, const MatrixBase<DerivedV>& b,
	            MatrixBase<DerivedV> const & x);

template <typename Derived>
void inv(const MatrixBase<Derived>& X,
	     MatrixBase<Derived> const & Xinv);

#include "leq.hpp"
#pragma once
#include <iostream>
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace std;

template <typename DerivedL, typename DerivedR>
void ASSERT_MATRIX_EQ(const Eigen::DenseBase<DerivedL>& A, const Eigen::DenseBase<DerivedR>& B, double tol=1e-7) {
	ASSERT_EQ(A.rows(), B.rows());
	ASSERT_EQ(A.cols(), B.cols());
	for(int i=0; i < A.rows(); i++) {
		for(int j=0; j < A.cols(); j++) {
			ASSERT_NEAR(A.derived().coeff(i,j), B(i,j), tol);
		}
	}
}
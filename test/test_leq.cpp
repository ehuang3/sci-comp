#include <iostream>
#include <gtest/gtest.h>
#include "test_helper.h"

#include <leq/leq.h>

using namespace std;
using namespace Eigen;

/* ********************************************************************************************* */
TEST(LEQ, FORW_SUB) {
	MatrixXd L(3,3);
	L.setRandom();
	for(int i=0; i < L.rows(); i++) {
		for(int j=0; j < L.cols(); j++) {
			if(i < j)
				L(i,j) = 0;
		}
	}
	VectorXd b = Vector3d::Ones();
	VectorXd x = Vector3d::Zero();
	forw_sub(L,b,x);
	VectorXd r = L*x;
	ASSERT_MATRIX_EQ(r, b);
}
/* ********************************************************************************************* */
TEST(LEQ, BACK_SUB) {
	MatrixXf U(3,3);
	U.setRandom();
	for(int i=0; i < U.rows(); i++) {
		for(int j=0; j < U.cols(); j++) {
			if(i > j)
				U(i,j) = 0;
		}
	}
	VectorXf b = Vector3f::Random();
	VectorXf x = Vector3f::Zero();
	back_sub(U,b,x);
	VectorXf r = U*x;
	ASSERT_MATRIX_EQ(r, b, 1e-6);
}
/* ********************************************************************************************* */
TEST(LEQ, LU_FAC) {
	MatrixXd X(3,3);
	// I.setIdentity();
	X.setRandom();
	MatrixXd L, U;
	lu_fac(X, L, U);
	MatrixXd R = L*U;
	ASSERT_MATRIX_EQ(R, X);

	Matrix3d I3d, L3, U3;
	I3d.setIdentity();
	lu_fac(I3d, L3, U3);
}
/* ********************************************************************************************* */
TEST(LEQ, GAUSS_ELIM) {
	MatrixXd A(10,10);
	VectorXd x(10), b(10);
	A.setRandom();
	x.setOnes();
	b = A*x;
	VectorXd r = x;
	gauss_elim(A,b,r);
	ASSERT_MATRIX_EQ(x,r);
}
/* ********************************************************************************************* */
TEST(LEQ, INV) {
	MatrixXd A(10,10), Ainv(10,10), AAinv;
	A.setRandom();
	inv(A, Ainv);
	AAinv = A*Ainv;
	MatrixXd I = A;
	I.setIdentity();
	ASSERT_MATRIX_EQ(AAinv,I);
}
/* ********************************************************************************************* */
int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
/* ********************************************************************************************* */

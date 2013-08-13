#pragma once

template <typename DerivedL, typename Derivedb, typename Derivedx>
void forw_sub(const MatrixBase<DerivedL>& L, const MatrixBase<Derivedb>& b,
	          MatrixBase<Derivedx> const & x) {
	typedef typename DerivedL::Scalar Scalar;

	MatrixBase<Derivedx>& x_ = const_cast< MatrixBase<Derivedx>& >(x);

	for(int i=0; i < L.rows(); i++) {
		Scalar s = b(i);
		for(int j=0; j < i; j++) {
			s -= L(i,j)*x_(j);
		}
		x_(i) = s/L(i,i);
	}
}

template <typename DerivedU, typename Derivedb, typename Derivedx>
void back_sub(const MatrixBase<DerivedU>& U, const MatrixBase<Derivedb>& b,
	          MatrixBase<Derivedx> const & x) {
	typedef typename DerivedU::Scalar Scalar;

	MatrixBase<Derivedx>& x_ = const_cast< MatrixBase<Derivedx>& >(x);

	for(int i=U.rows()-1; i >= 0; i--) {
		Scalar s = b(i);
		for(int j=U.rows()-1; j > i; j--) {
			s -= U(i,j)*x_(j);
		}
		x_(i) = s/U(i,i);
	}
}

template <typename Derived>
void lu_fac(const MatrixBase<Derived>& X,
	        MatrixBase<Derived> const & L, MatrixBase<Derived> const & U) {

	MatrixBase<Derived>& L_ = const_cast< MatrixBase<Derived>& >(L);
	MatrixBase<Derived>& U_ = const_cast< MatrixBase<Derived>& >(U);

	U_ = X;
	int r = U_.rows();
	for(int i=0; i < U_.rows()-1; i++) {
		U_.block(i+1,i,r-i-1,1) /= U_(i,i);
		for(int j=i+1; j < U_.cols(); j++) {
			U_.block(i+1,j,r-i-1,1) -= U_.block(i+1,i,r-i-1,1)*U_(i,j);
		}
	}

	// L_ = Derived::Identity();
	L_ = X;
	L_.setIdentity();
	for(int i=0; i < U_.rows()-1; i++) {
		L_.block(i+1,i,r-i-1,1) = U_.block(i+1,i,r-i-1,1);
		U_.block(i+1,i,r-i-1,1).setZero();
	}
}

template <typename DerivedA, typename Derivedb, typename Derivedx>
void gauss_elim(const MatrixBase<DerivedA>& A, const MatrixBase<Derivedb>& b,
	            MatrixBase<Derivedx> const & x) {

	MatrixBase<Derivedx>& x_ = const_cast< MatrixBase<Derivedx>& >(x);
	Derivedb y_(b.derived()); // If Derivedb is VectorXd

	DerivedA L_, U_;
	lu_fac(A, L_, U_);

	forw_sub(L_, b, y_);
	back_sub(U_, y_, x_);
}

template <typename Derived>
void inv(const MatrixBase<Derived>& X,
	     MatrixBase<Derived> const & Xinv) {
	typedef typename internal::plain_row_type<Derived>::type RowVectorType;	

	MatrixBase<Derived> I_ = X;
	I_.setIdentity();

	int r = X.rows();
	for(int i=0; i < X.cols(); i++) {
		// gauss_elim(X, I_.block(0,i,r,1), Xinv.block(0,i,r,1));
	}
}
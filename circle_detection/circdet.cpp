/*
 * circdet.cpp
 *
 *  Created on: Dec 12, 2011
 *      Author: ferraraf
 */


#include "circdet.hpp"


namespace circdet
{
// Detect a circle out of data points
void detect_circle( const MatrixXd & x, \
					Vector2d & x_c, \
					double & r )
{
	// Matrices for the Least-Squares equation
	Matrix3d A;
	Vector3d D;

	// Fill up the matrices for the Least-Squares equation
	fillmat( x, A, D );

	// Solve the Least-Squares Equation
	VectorXd X = solveLS( A, D );

	// Calculate the center of the circle
	x_c << X( 0 ), X( 1 );

	// Calculate the radius of the circle
	r = sqrt( X( 0 ) * X( 0 ) + X( 1 ) * X( 1 ) - X( 2 ) );
}

// Fill up the matrices for the Least-Squares equation
void fillmat( const MatrixXd & x, \
			  Matrix3d & A, \
			  Vector3d & D )
{
	// Get the length of the data
	double n = x.rows();

	// Initialize matrix coefficients
	// Matrix A
	double A11 = 0.0;
	double A12 = 0.0;
	double A13 = 0.0;
	double A22 = 0.0;
	double A23 = 0.0;

	// Matrix D
	double D1 = 0.0;
	double D2 = 0.0;
	double D3 = 0.0;

	// Calculate the coefficients
	for ( int i = 0; i < n; i++ )
	{
		// Temporary values
		double temp1 = x( i, 0 ) * x( i, 0 );
		double temp2 = x( i, 1 ) * x( i, 1 );

		// Matrix A
		A11 = A11 + temp1;
		A12 = A12 + x( i, 0 ) * x( i, 1 );
		A13 = A13 - x( i, 0 );
		A22 = A22 + temp2;
		A23 = A23 - x( i, 1 );

		// Matrix D
		D1 = D1 + temp1 * x( i, 0 ) + x( i, 0 ) * temp2;
		D2 = D2 + temp1 * x( i, 1 ) + x( i, 1 ) * temp2;
		D3 = D3 + temp1 + temp2;
	}

	// Matrix A
	A11 = 2 * A11;
	A12 = 2 * A12;
	double A21 = A12;
	A22 = 2 * A22;
	double A31 = -2 * A13;
	double A32 = -2 * A23;
	double A33 = -n;

	A << A11, A12, A13, \
		 A21, A22, A23, \
		 A31, A32, A33;

	// Matrix D
	D << D1, D2, D3;
}

// Solve the Least-Squares Equation
VectorXd solveLS( const Matrix3d & A, \
				  const Vector3d & D, \
				  const double eps )
{
	// Do a singular value decomposition of the matrix
	JacobiSVD< Matrix3d, FullPivHouseholderQRPreconditioner > svd( A, ComputeFullU | ComputeFullV );

	// Set the tolerance
	double tol = eps * max( ( int )A.rows (), ( int )A.cols ()) * svd.singularValues().maxCoeff();

	// Put the singular values into a matrix
	Matrix3d sv = svd.singularValues().asDiagonal();

	// Calculate the singular values according to the tolerance
	for ( int i = 0; i < ( int )sv.cols(); i++ )
	{
		// If the singular values greater than the tolerance
		if ( sv( i, i ) > tol )
		{
			// Take the inverse of the singular values and put it back into the diagonal
			sv(i, i) = 1.0 / sv(i, i);
		}
		else // else
		{
			// Set them to 0
			sv(i, i) = 0.0;
		}
	}

	// Calculate the pseudo inverse
	Matrix3d A_inv = svd.matrixV() * sv.transpose() * svd.matrixU().adjoint();

	// Calculate X
	Vector3d X = A_inv * D;

	return X;
}
}

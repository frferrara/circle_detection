/*
 * circdet.hpp
 *
 *  Created on: Dec 12, 2011
 *      Author: ferraraf
 */

#ifndef CIRCDET_HPP_
#define CIRCDET_HPP_


#include <iostream>
#include <Eigen/Eigen>


using namespace std;
using namespace Eigen;


namespace circdet
{
// Detect a circle out of data points
void detect_circle( const MatrixXd & x, \
					Vector2d & x_c, \
					double & r );

// Fill up the matrices for the Least-Squares equation
void fillmat( const MatrixXd & x, \
			  Matrix3d & A, \
			  Vector3d & D );

// Solve the Least-Squares Equation
VectorXd solveLS( const Matrix3d & A, \
				  const Vector3d & D, \
				  const double eps = numeric_limits<double>::epsilon() );
}


#endif /* CIRCDET_HPP_ */

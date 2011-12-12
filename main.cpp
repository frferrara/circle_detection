/*
 * main.cpp
 *
 *  Created on: Dec 12, 2011
 *      Author: ferraraf
 */


#include "circdet.hpp"
#include <cmath>


int main( int argc, char * argv[] )
{
	cout << "Test of circle detection!" << endl;

	// Create a set of circle points
	// Circle properties
	double x0 = 2.0;
	double y0 = 4.0;
	double r = 2.0;

	// Generation of circle points
	// Point distance
	double dphi = 0.1;

	// Preallocation
	double n = round( 2 * M_PI / dphi ) + 1;
	MatrixXd x;
	x.resize( n, 2 );
	double phi = 0.0;

	// Generation
	for ( int i = 0; i < n; i++ )
	{
		// Calculate the points
		x( i, 0 ) = x0 + r * cos( phi );
		x( i, 1 ) = y0 + r * sin( phi );

		// Augment the angle
		phi = phi + dphi;
	}

	// Circle detection
	Vector2d x_c;
	circdet::detect_circle( x, x_c, r );

	cout << "\nx_c: " << x_c.transpose() << "\nr: " << r << endl;

	return 0;
}

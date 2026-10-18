//$ nocpp

/**
 * @file nmsopt.h
 *
 * @version 2024.6
 *
 * @brief The inclusion file for the CNMSeqOpt class.
 *
 * @section license License
 * 
 * Copyright (c) 2016-2024 Aleksey Vaneev
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef NMSOPT_INCLUDED
#define NMSOPT_INCLUDED

#include "biteaux.h"

/**
 * Sequential Nelder-Mead simplex method. Features custom coefficients tuned
 * to provide better convergence at higher dimensions. Also implements various
 * algorithmic optimizations.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CNMSeqOpt : public CBiteOptBase< double >
{
public:
	CNMSeqOpt()
		: y( NULL )
		, x2( NULL )
	{
	}

	virtual ~CNMSeqOpt()
	{
		delete[] y;
		delete[] x2;
	}

	/**
	 * Function updates dimensionality of *this object.
	 *
	 * @param aParamCount The number of parameters being optimized.
	 * @param PopSize0 The number of elements in population to use. If set to
	 * 0 or negative, the default formula will be used.
	 */

	void updateDims( const int aParamCount, const int PopSize0 = 0 )
	{
		const int aPopSize = ( PopSize0 > 0 ? PopSize0 :
			( aParamCount + 1 ) * 4 );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		N = aParamCount;
		M = aPopSize;
		M1 = aPopSize - 1;
		M1i = 1.0 / M1;

		initBuffers( N, M );
	}

	/**
	 * Function initializes *this optimizer. Performs N=PopSize objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams If not NULL, initial parameter vector, also used as
	 * centroid.
	 * @param InitRadius Initial radius, relative to the default value.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		initCommonVars( rnd );

		double* const xx = x[ 0 ];
		int i;
		int j;

		if( InitParams != NULL )
		{
			copyParams( xx, InitParams );
		}
		else
		{
			for( i = 0; i < N; i++ )
			{
				xx[ i ] = MinValues[ i ] + DiffValues[ i ] * 0.5;
			}
		}

		xlo = 0;

		const double sd = 0.25 * InitRadius;

		for( j = 1; j < M; j++ )
		{
			double* const xj = x[ j ];

			for( i = 0; i < N; i++ )
			{
				xj[ i ] = xx[ i ] + DiffValues[ i ] * rnd.getGaussian() * sd;
			}
		}

		State = stReflection;
	}

	/**
	 * Function performs the parameter optimization iteration that involves 1
	 * objective function evaluation.
	 *
	 * @param rnd Random number generator.
	 * @return The number of non-improving iterations so far.
	 */

	int optimize( CBiteRnd& rnd )
	{
		int i;

		if( DoInitEvals )
		{
			y[ CurPopPos ] = eval( rnd, x[ CurPopPos ]);

			if( y[ CurPopPos ] < y[ xlo ])
			{
				xlo = CurPopPos;
			}

			CurPopPos++;

			if( CurPopPos == M )
			{
				DoInitEvals = false;
				calccent();
			}

			return( 0 );
		}

		StallCount++;

		const double sn = 0.5 * sqrt( ParamCountI );

		static const double alpha = 1.0; // Reflection coeff.
		const double gamma = 1.5 + sn; // Expansion coeff.
		const double rho = -0.75 + sn; // Contraction coeff.
		const double sigma = 1.0 - sn; // Reduction coeff.

		double* const xH = x[ xhi ]; // Highest cost parameter vector.

		switch( State )
		{
			case stReflection:
			{
				for( i = 0; i < N; i++ )
				{
					x1[ i ] = x0[ i ] + alpha * ( x0[ i ] - xH[ i ]);
				}

				y1 = eval( rnd, x1 );

				if( y1 > y[ xlo ] && y1 < y[ xhi2 ])
				{
					copy( x1, y1 );
				}
				else
				{
					if( y1 < y[ xlo ])
					{
						State = stExpansion;
						StallCount--;
					}
					else
					{
						State = stContraction;
					}
				}

				break;
			}

			case stExpansion:
			{
				for( i = 0; i < N; i++ )
				{
					x2[ i ] = x0[ i ] + gamma * ( x0[ i ] - xH[ i ]);
				}

				const double y2 = eval( rnd, x2 );
				xlo = xhi;

				if( y2 < y1 )
				{
					copy( x2, y2 );
				}
				else
				{
					copy( x1, y1 );
				}

				State = stReflection;
				break;
			}

			case stContraction:
			{
				for( i = 0; i < N; i++ )
				{
					x2[ i ] = x0[ i ] + rho * ( x0[ i ] - xH[ i ]);
				}

				const double y2 = eval( rnd, x2 );

				if( y2 < y[ xhi ])
				{
					if( y2 < y[ xlo ])
					{
						xlo = xhi;
					}

					copy( x2, y2 );
					State = stReflection;
				}
				else
				{
					rx = x[ xlo ];
					rj = 0;
					State = stReduction;
				}

				break;
			}

			case stReduction:
			{
				double* xx = x[ rj ];

				if( xx == rx )
				{
					rj++;
					xx = x[ rj ];
				}

				for( i = 0; i < N; i++ )
				{
					xx[ i ] = rx[ i ] + sigma * ( xx[ i ] - rx[ i ]);
				}

				y[ rj ] = eval( rnd, xx );

				if( y[ rj ] < y[ xlo ])
				{
					xlo = rj;
					StallCount = 0;
				}

				rj++;

				if( rj == M || ( rj == M1 && x[ rj ] == rx ))
				{
					calccent();
					State = stReflection;
				}

				break;
			}
		}

		return( StallCount );
	}

private:
	int N; ///< The total number of internal parameter values in use.
	int M; ///< The number of points in a simplex.
	int M1; ///< = M - 1.
	double M1i; ///< = 1.0 / ( M - 1 ).
	int xlo; ///< Current lowest-cost parameter vector.
	int xhi; ///< Current highest-cost parameter vector.
	int xhi2; ///< Current second-highest-cost parameter vector.
	double** x; ///< Parameter vectors for all points.
	double* y; ///< Parameter vector costs.
	double* x0; // Centroid parameter vector.
	double* x1; ///< Temporary parameter vector 1. Passed to stExpansion.
	double y1; ///< Cost of temporary parameter vector 1. Passed to
		///< stExpansion.
	double* x2; ///< Temporary parameter vector 2.
	double* rx; ///< Lowest-cost parameter vector used during reduction.
	int rj; ///< Current vector index during reduction.

	/**
	 * Algorithm's state automata states.
	 */

	enum EState
	{
		stReflection, // Reflection.
		stExpansion, // Expansion.
		stContraction, // Contraction.
		stReduction // Reduction.
	};

	EState State; ///< Current optimization state.

	virtual void initBuffers( const int aParamCount, const int aPopSize,
		const int aCnsCount = 0, const int aObjCount = 1 )
	{
		CBiteOptBase :: initBuffers( aParamCount, aPopSize, aCnsCount,
			aObjCount );

		x = PopParams;
		y = new double[ M ];
		x0 = CentParams;
		x1 = TmpParams;
		x2 = new double[ N ];
	}

	virtual void deleteBuffers()
	{
		CBiteOptBase :: deleteBuffers();

		delete[] y;
		delete[] x2;
	}

	/**
	 * Function finds the highest-cost vector.
	 */

	void findhi()
	{
		xhi2 = ( y[ 0 ] > y[ 1 ]);
		xhi = 1 - xhi2;

		int j;

		for( j = 2; j < M; j++ )
		{
			if( y[ j ] > y[ xhi ])
			{
				xhi2 = xhi;
				xhi = j;
			}
			else
			if( y[ j ] > y[ xhi2 ])
			{
				xhi2 = j;
			}
		}
	}

	/**
	 * Function calculates the centroid vector that excludes a highest-cost
	 * parameter vector.
	 */

	void calccent()
	{
		findhi();

		double* const xc = x0;
		int j;
		int i;

		if( xhi == 0 )
		{
			copyParams( xc, x[ 1 ]);
			j = 1;
		}
		else
		{
			copyParams( xc, x[ 0 ]);

			for( j = 1; j < xhi; j++ )
			{
				const double* const xx = x[ j ];

				for( i = 0; i < N; i++ )
				{
					xc[ i ] += xx[ i ];
				}
			}
		}

		for( j = j + 1; j < M; j++ )
		{
			const double* const xx = x[ j ];

			for( i = 0; i < N; i++ )
			{
				xc[ i ] += xx[ i ];
			}
		}

		for( i = 0; i < N; i++ )
		{
			xc[ i ] *= M1i;
		}
	}

	/**
	 * Function replaces the highest-cost vector with a new vector. Function
	 * also resets the StallCount to 0.
	 *
	 * @param ip Input vector.
	 * @param Cost Input vector's cost.
	 */

	void copy( const double* const ip, const double Cost )
	{
		y[ xhi ] = Cost;
		double* const xH = x[ xhi ];

		copyParams( xH, ip );
		findhi();

		const double* const nxH = x[ xhi ];

		if( xH != nxH )
		{
			int i;

			for( i = 0; i < N; i++ )
			{
				x0[ i ] += ( xH[ i ] - nxH[ i ]) * M1i;
			}
		}

		StallCount = 0;
	}

	/**
	 * Function evaluates parameter vector and applies value range wrapping,
	 * also records a new best solution.
	 *
	 * @param rnd Random number generator.
	 * @param Params Parameter vector to evaluate.
	 */

	double eval( CBiteRnd& rnd, const double* const Params )
	{
		int i;

		for( i = 0; i < N; i++ )
		{
			NewValues[ i ] = wrapParamReal( rnd, Params[ i ], i );
		}

		const double NewCost = fixCostNaN( optcost( NewValues ));
		NewCosts[ 0 ] = NewCost;

		updateBestCost( NewCost, NewValues );

		return( NewCost );
	}
};

#endif // NMSOPT_INCLUDED

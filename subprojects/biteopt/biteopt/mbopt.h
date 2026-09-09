//$ nocpp

/**
 * @file mbopt.h
 *
 * @version 2024.6
 *
 * @brief The inclusion file for the CMiniBiteOpt class.
 *
 * Description is available at https://github.com/avaneev/biteopt
 *
 * E-mail: aleksey.vaneev@gmail.com or info@voxengo.com
 *
 * @section license License
 *
 * Copyright (c) 2024 Aleksey Vaneev
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

#ifndef MBOPT_INCLUDED
#define MBOPT_INCLUDED

#include "nmsopt.h"

/**
 * Minimal BiteOpt optimization class, in real parameter space.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CMiniBiteOpt : public CBiteOptBase< double >
{
public:
	CMiniBiteOpt()
		: ParOpt( this )
	{
		addSel( MethodSel, "MethodSel" );
	}

	/**
	 * Function updates dimensionality of *this object. Function does nothing
	 * if dimensionality has not changed since the last call. This function
	 * should be called at least once before calling the init() function.
	 *
	 * @param aParamCount The number of parameters being optimized.
	 * @param PopSize0 The number of elements in population to use. If set to
	 * 0 or negative, the default formula will be used.
	 */

	void updateDims( const int aParamCount, const int PopSize0 = 0 )
	{
		const int aPopSize = ( PopSize0 > 0 ? PopSize0 :
			calcPopSizeBiteOpt( aParamCount ));

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		initBuffers( aParamCount, aPopSize );

		ParOpt.updateDims( aParamCount, aPopSize * 4 / 3 );
	}

	/**
	 * Function initializes *this optimizer. Does not perform objective
	 * function evaluations.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams If not NULL, initial parameter vector, also used as
	 * centroid for initial population vectors.
	 * @param InitRadius Initial radius, multiplier relative to the default
	 * sigma value.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		initCommonVars( rnd );

		StartSD = 0.25 * InitRadius;

		if( InitParams != NULL )
		{
			copyParams( StartParams, InitParams );
			UseStartParams = true;
		}

		ParOpt.init( rnd, InitParams, InitRadius );
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
		double* Params;
		int i;

		if( DoInitEvals )
		{
			Params = getCurParams();
			LastValues = Params;

			genInitParamsReal( rnd, Params );

			NewCosts[ 0 ] = fixCostNaN( optcost( Params ));
			updateBestCost( NewCosts[ 0 ], Params,
				updatePop( NewCosts[ 0 ], Params ));

			if( CurPopPos == PopSize )
			{
				updateCentroid();
				DoInitEvals = false;
			}

			return( 0 );
		}

		Params = TmpParams;
		LastCosts = NewCosts;
		LastValues = TmpParams;

		bool DoEval = true;
		const int SelMethod = select( MethodSel, rnd );

		if( SelMethod == 0 )
		{
			// A variant of BiteOpt generator 1.

			copyParams( Params, getParamsOrdered( rnd.getPowInt( 3.0, 4 )));

			// Convert real value to normalized value, apply the "bitmask
			// inversion" operation, and convert back.

			i = rnd.getInt( ParamCount );
			Params[ i ] = ( (int64_t) (( Params[ i ] - MinValues[ i ]) *
				DiffValuesI[ i ] * MantMult ) ^
				( IntMantMask >> rnd.getPowInt( 4.0, 32 ))) * MantMultI *
				DiffValues[ i ] + MinValues[ i ];

			const double* const rp2 = getParamsOrdered(
				rnd.getSqrInt( CurPopSize ));

			const double v = rp2[ i ];
			Params[ i ] += ( v - Params[ i ]) * rnd.getTPDF();
			Params[ i ] += ( v - Params[ i ]) * rnd.getTPDF();
			Params[ i ] += ( v - Params[ i ]) * rnd.getTPDF();
			Params[ i ] += ( v - Params[ i ]) * rnd.getTPDF();
		}
		else
		if( SelMethod == 1 )
		{
			// A variant of BiteOpt generator 2.

			const int si1 = rnd.getPowInt( 3.0, 4 );
			const double* const rp1 = getParamsOrdered( si1 );
			const double* const rp3 = getParamsOrdered( CurPopSize1 - si1 );

			const int si2 = 1 + rnd.getInt( CurPopSize1 );
			const double* const rp2 = getParamsOrdered( si2 );

			const int si4 = rnd.getSqrInt( CurPopSize );
			const double* const rp4 = getParamsOrdered( si4 );
			const double* const rp5 = getParamsOrdered( CurPopSize1 - si4 );

			if( rnd.getBit() )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = rp1[ i ] + (( rp2[ i ] - rp3[ i ]) +
						( rp4[ i ] - rp5[ i ])) * 0.5;
				}
			}
			else
			{
				const double* const rp1b = getParamsOrdered(
					rnd.getSqrInt( CurPopSize ));

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = (( rp1[ i ] + rp1b[ i ]) +
						( rp2[ i ] - rp3[ i ]) + ( rp4[ i ] - rp5[ i ])) * 0.5;
				}
			}
		}
		else
		if( SelMethod == 2 )
		{
			if( rnd.getBit() && rnd.getBit() )
			{
				// A variant of BiteOpt generator 3.

				const double* const rp1 = getParamsOrdered(
					rnd.getPowInt( 3.0, 4 ));

				const double* const rp2 = getParamsOrdered(
					rnd.getLogInt( CurPopSize ));

				const double* const cp = getCentroid();

				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = ( rnd.getBit() ? cp[ i ] :
						rp1[ i ] + ( rp1[ i ] - rp2[ i ]));
				}
			}
			else
			{
				// Parallel optimizer.

				const int sc = ParOpt.optimize( rnd );

				LastCosts = ParOpt.getLastCosts();
				LastValues = ParOpt.getLastValues();
				DoEval = false;

				if( sc > ParamCount * 16 )
				{
					ParOpt.init( rnd, getBestParams(), StartSD * 2.0 );
				}
			}
		}

		if( DoEval )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = wrapParamReal( rnd, Params[ i ], i );
			}

			NewCosts[ 0 ] = fixCostNaN( optcost( Params ));
		}

		const int p = updatePop( LastCosts[ 0 ], LastValues, true, 3 );

		if( p > CurPopSize1 )
		{
			applySelsDecr( rnd );

			StallCount++;
		}
		else
		{
			updateBestCost( LastCosts[ 0 ], LastValues, p );
			applySelsIncr( rnd, 1.0 - p * CurPopSizeI );

			StallCount = 0;
		}

		return( StallCount );
	}

protected:
	CBiteSel< 3 > MethodSel; ///< Population generator selector.
	CBiteOptOwned< CNMSeqOpt > ParOpt; ///< Parallel optimizer.
};

#endif // MBOPT_INCLUDED

//$ nocpp

/**
 * @file spheropt.h
 *
 * @version 2024.6
 *
 * @brief The inclusion file for the CSpherOpt class.
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

#ifndef SPHEROPT_INCLUDED
#define SPHEROPT_INCLUDED

#include "biteaux.h"

/**
 * "Converging hyper-spheroid" optimizer class. Simple, converges quite fast.
 *
 * Description is available at https://github.com/avaneev/biteopt
 */

class CSpherOpt : public CBiteOptBase< double >
{
public:
	CSpherOpt()
		: WPopCent( NULL )
		, WPopRad( NULL )
	{
		addSel( CentPowSel, "CentPowSel" );
		addSel( RadPowSel, "RadPowSel" );
		addSel( EvalFacSel, "EvalFacSel" );
	}

	virtual ~CSpherOpt()
	{
		delete[] WPopCent;
		delete[] WPopRad;
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
		const int aPopSize = ( PopSize0 > 0 ? PopSize0 : 14 + aParamCount );

		if( aParamCount == ParamCount && aPopSize == PopSize )
		{
			return;
		}

		initBuffers( aParamCount, aPopSize );

		JitMult = 5.0 * ParamCountI;
		JitOffs = 1.0 - JitMult * 0.5;
	}

	/**
	 * Function initializes *this optimizer.
	 *
	 * @param rnd Random number generator.
	 * @param InitParams If not NULL, initial parameter vector, also used as
	 * centroid.
	 * @param InitRadius Initial radius, multiplier relative to the default
	 * sigma value.
	 */

	void init( CBiteRnd& rnd, const double* const InitParams = NULL,
		const double InitRadius = 1.0 )
	{
		initCommonVars( rnd );

		Radius = 0.5 * InitRadius;
		EvalFac = 2.0;
		cure = 0;
		curem = (int) ceil( CurPopSize * EvalFac );

		// Provide initial centroid and sigma.

		int i;

		if( InitParams == NULL )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				CentParams[ i ] = 0.5;
			}

			DoInitEvals = false;
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				CentParams[ i ] = wrapParam( rnd,
					( InitParams[ i ] - MinValues[ i ]) / DiffValues[ i ]);
			}
		}
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
		double* const Params = getCurParams();
		int i;

		if( DoInitEvals )
		{
			DoInitEvals = false;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = CentParams[ i ];
				NewValues[ i ] = getRealValue( CentParams, i );
			}
		}
		else
		{
			double s2 = 1e-300;

			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = rnd.get() - 0.5;
				s2 += Params[ i ] * Params[ i ];
			}

			const double d = Radius / sqrt( s2 );

			if( ParamCount > 4 )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = wrapParam( rnd,
						CentParams[ i ] + Params[ i ] * d );

					NewValues[ i ] = getRealValue( Params, i );
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					const double m = JitOffs + rnd.get() * JitMult;

					Params[ i ] = wrapParam( rnd,
						CentParams[ i ] + Params[ i ] * d * m );

					NewValues[ i ] = getRealValue( Params, i );
				}
			}
		}

		const double NewCost = fixCostNaN( optcost( NewValues ));
		NewCosts[ 0 ] = NewCost;

		updatePop( NewCost, Params );
		updateBestCost( NewCost, NewValues );

		AvgCost += NewCost;
		cure++;

		if( cure >= curem )
		{
			AvgCost /= cure;

			if( AvgCost < HiBound )
			{
				HiBound = AvgCost;

				applySelsIncr( rnd );
			}
			else
			{
				applySelsDecr( rnd );
			}

			resetCurPopPos();
			AvgCost = 0.0;
			cure = 0;

			update( rnd );

			curem = (int) ceil( CurPopSize * EvalFac );
		}

		StallCount = ( NewCost < HiBound ? 0 : StallCount + 1 );

		return( StallCount );
	}

protected:
	double* WPopCent; ///< Weighting coefficients for centroid.
	double* WPopRad; ///< Weighting coefficients for radius.
	double JitMult; ///< Jitter multiplier.
	double JitOffs; ///< Jitter multiplier offset.
	double Radius; ///< Current radius.
	double EvalFac; ///< Evaluations factor.
	int cure; ///< Current evaluation index.
	int curem; ///< "cure" value threshold.
	CBiteSel< 4 > CentPowSel; ///< Centroid power factor selector.
	CBiteSel< 4 > RadPowSel; ///< Radius power factor selector.
	CBiteSel< 3 > EvalFacSel; ///< EvalFac selector.

	virtual void initBuffers( const int aParamCount, const int aPopSize,
		const int aCnsCount = 0, const int aObjCount = 1 )
	{
		CBiteOptBase< double > :: initBuffers( aParamCount, aPopSize,
			aCnsCount, aObjCount );

		WPopCent = new double[ aPopSize ];
		WPopRad = new double[ aPopSize ];
	}

	virtual void deleteBuffers()
	{
		CBiteOptBase< double > :: deleteBuffers();

		delete[] WPopCent;
		delete[] WPopRad;
	}

	/**
	 * Function updates centroid and radius.
	 *
	 * @param rnd PRNG object.
	 */

	void update( CBiteRnd& rnd )
	{
		static const double WCent[ 4 ] = { 4.5, 6.0, 7.5, 10.0 };
		static const double WRad[ 4 ] = { 14.0, 16.0, 18.0, 20.0 };
		static const double EvalFacs[ 3 ] = { 2.1, 2.0, 1.9 };

		const double CentFac = WCent[ select( CentPowSel, rnd )];
		const double RadFac = WRad[ select( RadPowSel, rnd )];
		EvalFac = EvalFacs[ select( EvalFacSel, rnd )];

		const double lm = 1.0 / curem;
		double s1 = 0.0;
		double s2 = 0.0;
		int i;

		for( i = 0; i < CurPopSize; i++ )
		{
			const double l = 1.0 - i * lm;

			const double v1 = pow( l, CentFac );
			WPopCent[ i ] = v1;
			s1 += v1;

			const double v2 = pow( l, RadFac );
			WPopRad[ i ] = v2;
			s2 += v2;
		}

		s1 = 1.0 / s1;
		s2 = 1.0 / s2;

		const double* ip = PopParams[ 0 ];
		double* const cp = CentParams;
		const double* const wc = WPopCent;
		double w = wc[ 0 ] * s1;

		for( i = 0; i < ParamCount; i++ )
		{
			cp[ i ] = ip[ i ] * w;
		}

		int j;

		for( j = 1; j < CurPopSize; j++ )
		{
			ip = PopParams[ j ];
			w = wc[ j ] * s1;

			for( i = 0; i < ParamCount; i++ )
			{
				cp[ i ] += ip[ i ] * w;
			}
		}

		const double* const rc = WPopRad;
		Radius = 0.0;

		for( j = 0; j < CurPopSize; j++ )
		{
			ip = PopParams[ j ];
			double s = 0.0;

			for( i = 0; i < ParamCount; i++ )
			{
				const double d = ip[ i ] - cp[ i ];
				s += d * d;
			}

			Radius += s * rc[ j ];
		}

		Radius = sqrt( Radius * s2 );
	}
};

#endif // SPHEROPT_INCLUDED

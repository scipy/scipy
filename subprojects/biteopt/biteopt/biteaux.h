//$ nocpp

/**
 * @file biteaux.h
 *
 * @version 2024.6
 *
 * @brief The inclusion file for the CBiteRnd, CBitePop, CBiteParPops,
 * CBiteOptInterface, CBiteOptBase, and CBiteOptOwned classes.
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

#ifndef BITEAUX_INCLUDED
#define BITEAUX_INCLUDED

#include <stdint.h>
#include <math.h>
#include <string.h>

/**
 * Type for an externally-provided random number generator, to be used instead
 * of the default PRNG. Note that if the external produces 64-bit random
 * values, they can be safely truncated/typecasted to the "uint32_t" type. If
 * the external PRNG produces floating-point values, they should be scaled to
 * the 32-bit unsigned integer range. 32-bit PRNG output is required for
 * compatibility with older 32-bit PRNGs (but CBiteRnd is 64-bit PRNG).
 */

typedef uint32_t( *biteopt_rng )( void* rng_data );

/**
 * Class that implements a pseudo-random number generator (PRNG). The default
 * implementation includes a fast high-quality PRNG (2^159 period). See
 * https://github.com/avaneev/prvhash for more details.
 */

class CBiteRnd
{
public:
	/**
	 * Default constructor, calls the init() function.
	 */

	CBiteRnd()
	{
		init( 1 );
	}

	/**
	 * Constructor, calls the init() function.
	 *
	 * @param NewSeed Random seed value.
	 */

	CBiteRnd( const int NewSeed )
	{
		init( NewSeed );
	}

	/**
	 * Function initializes *this PRNG object.
	 *
	 * @param NewSeed New random seed value. Ignored, if "arf" is non-NULL.
	 * @param arf External random number generator to use; NULL: use the
	 * default PRNG. Note that the external RNG should be seeded externally.
	 * @param ardata Data pointer to pass to the "arf" function.
	 */

	void init( const int NewSeed, biteopt_rng const arf = NULL,
		void* const ardata = NULL )
	{
		rf = arf;
		rdata = ardata;

		BitsLeft = 0;
		Seed = (uint64_t) NewSeed;
		lcg = 0;
		Hash = 0;

		// Skip first values to make PRNG "settle down".

		int i;

		for( i = 0; i < 5; i++ )
		{
			advance();
		}
	}

	/**
	 * @return Random number in the range [0; 1).
	 */

	double get()
	{
		return(( advance() >> ( 64 - 53 )) * 0x1p-53 );
	}

	/**
	 * @param N1 Integer value range.
	 * @return Random integer number in the range [0; N1). N1 denotes the
	 * number of bins, not the maximal returned value.
	 */

	int getInt( const int N1 )
	{
		return( (int) ( get() * N1 ));
	}

	/**
	 * @return Random number in the range [0; 1). Beta distribution with
	 * Alpha=0.5, Beta=1 (squared).
	 */

	double getSqr()
	{
		const double v = ( advance() >> ( 64 - 53 )) * 0x1p-53;

		return( v * v );
	}

	/**
	 * @return Random number in the range [0; 1), raised to the specified
	 * power. The function has branching for optimization.
	 */

	double getPow( const double p )
	{
		const double v = ( advance() >> ( 64 - 53 )) * 0x1p-53;

		if( p < 2.0 )
		{
			if( p < 1.0 )
			{
				if( p == 0.5 )
				{
					return( sqrt( v ));
				}

				if( p == 0.25 )
				{
					return( sqrt( sqrt( v )));
				}
			}
			else
			{
				if( p == 1.5 )
				{
					return( v * sqrt( v ));
				}

				if( p == 1.75 )
				{
					const double sv = sqrt( v );
					return( v * sv * sqrt( sv ));
				}

				if( p == 1.0 )
				{
					return( v );
				}
			}
		}
		else
		{
			if( p == 4.0 )
			{
				const double v2 = v * v;
				return( v2 * v2 );
			}

			if( p == 3.0 )
			{
				return( v * v * v );
			}

			if( p == 2.0 )
			{
				return( v * v );
			}
		}

		return( pow( v, p ));
	}

	/**
	 * @return Random number in the range (-1; 1) with approximately
	 * logarithmic PDF, two-lobe with peak at 0.
	 */

	double getLog()
	{
		return( get() * sin( get() * 6.28318530717958648 ));
	}

	/**
	 * @param N1 Integer value range.
	 * @return Random integer number in the range [0; N1). Beta distribution
	 * with Alpha=0.5, Beta=1 (squared). N1 denotes the number of bins, not
	 * the maximal returned value.
	 */

	int getSqrInt( const int N1 )
	{
		return( (int) ( getSqr() * N1 ));
	}

	/**
	 * @param N1 Integer value range.
	 * @return Random integer number in the reversed range (N1; 0]. Beta
	 * distribution with Alpha=0.5, Beta=1 (squared). N1 denotes the number of
	 * bins, not the maximal returned value.
	 */

	int getSqrIntInv( const int N1 )
	{
		return( N1 - (int) ( getSqr() * N1 ) - 1 );
	}

	/**
	 * @param N1 Integer value range.
	 * @return Random integer number in the range [0; N1), raised to the
	 * specified power. N1 denotes the number of bins, not the maximal
	 * returned value.
	 */

	int getPowInt( const double p, const int N1 )
	{
		return( (int) ( getPow( p ) * N1 ));
	}

	/**
	 * @param N1 Integer value range.
	 * @return Random integer number in the range [0; N1), with approximately
	 * logarithmic PDF, peak at 0.
	 */

	int getLogInt( const int N1 )
	{
		return( (int) ( fabs( getLog() ) * N1 ));
	}

	/**
	 * @return Uniformly-distributed random number in the "raw" scale.
	 */

	uint64_t getRaw()
	{
		return( advance() );
	}

	/**
	 * @return TPDF random number in the range (-1; 1).
	 */

	double getTPDF()
	{
		const int64_t v1 = (int64_t) ( advance() >> ( 64 - 53 ));
		const int64_t v2 = (int64_t) ( advance() >> ( 64 - 53 ));

		return(( v1 - v2 ) * 0x1p-53 );
	}

	/**
	 * Function generates a Gaussian-distributed pseudo-random number with
	 * mean=0 and std.dev=1.
	 *
	 * Algorithm is adopted from "Leva, J. L. 1992. "A Fast Normal Random
	 * Number Generator", ACM Transactions on Mathematical Software, vol. 18,
	 * no. 4, pp. 449-453".
	 */

	double getGaussian()
	{
		double q, u, v;

		do
		{
			u = get();
			v = get();

			if( u == 0.0 || v == 0.0 )
			{
				u = 1.0;
				v = 1.0;
			}

			v = 1.7156 * ( v - 0.5 );
			const double x = u - 0.449871;
			const double y = fabs( v ) + 0.386595;
			q = x * x + y * ( 0.19600 * y - 0.25472 * x );

			if( q < 0.27597 )
			{
				break;
			}

		} while( q > 0.27846 || v * v > -4.0 * log( u ) * u * u );

		return( v / u );
	}

	/**
	 * Function returns the next random bit, usually used for 50% probability
	 * evaluations efficiently.
	 */

	int getBit()
	{
		if( BitsLeft == 0 )
		{
			BitPool = advance();

			const int b = (int) ( BitPool & 1 );

			BitsLeft = 63;
			BitPool >>= 1;

			return( b );
		}

		const int b = (int) ( BitPool & 1 );

		BitsLeft--;
		BitPool >>= 1;

		return( b );
	}

protected:
	biteopt_rng rf; ///< External random number generator to use; NULL: use
		///< the default PRNG.
	void* rdata; ///< Data pointer to pass to the "rf" function.
	uint64_t Seed, lcg, Hash; ///< PRNG state variables.
	uint64_t BitPool; ///< Bit pool.
	int BitsLeft; ///< The number of bits left in the bit pool.

	/**
	 * Function advances the PRNG and returns the next PRNG value.
	 */

	uint64_t advance()
	{
		if( rf != NULL )
		{
			uint64_t r = ( *rf )( rdata );
			r |= (uint64_t) ( *rf )( rdata ) << 32;

			return( r );
		}

		Seed *= lcg * 2 + 1;
		const uint64_t rs = Seed >> 32 | Seed << 32;
		Hash += rs + 0xAAAAAAAAAAAAAAAA;
		lcg += Seed + 0x5555555555555555;
		Seed ^= Hash;

		return( lcg ^ rs );
	}
};

/**
 * Probabilistic selector class. Used to keep track of success of various
 * choices. Updates probabilities of future choices based on the selection
 * outcome.
 *
 * The current implementation uses bubble-sort-alike method to update a sparse
 * vector of possible choices. The selection is made as a weighted-random draw
 * of a value from this vector. Previously, the class used simple statistical
 * accumulation of optimization outcomes, to derive probabilites. The current
 * approach is superior in that it has no "memory effects" associated with
 * statistical accumulation.
 *
 * The purpose of the class is to increase a chance of generating an
 * acceptable solution. In practice, this class provides 10-15% more "good"
 * solutions compared to uniformly-random choice selection. This, in turn,
 * improves convergence smoothness and produces more diversity in outcomes in
 * multiple solution attempts of complex multi-modal objective functions.
 */

class CBiteSelBase
{
public:
	/**
	 * Constructor.
	 *
	 * @param Count The number of possible choices, greater than 1.
	 * @param spwr100 Selector's power factor multiplied by 100. At 100,
	 * all choices are equally probable; at 150 the best choice may be
	 * selected up to two times more often than any other choice.
	 */

	CBiteSelBase( const int aCount, const int spwr100 )
		: Count( aCount )
		, SelPower( spwr100 * 0.01 )
		, SelBuf( NULL )
		, SelBufCapacity( 0 )
	{
	}

	~CBiteSelBase()
	{
		delete[] SelBuf;
	}

	/**
	 * This function resets *this object, should be called before calling
	 * other functions, including after object's construction.
	 *
	 * @param rnd PRNG object.
	 * @param ParamCount The number of dimensions being optimized.
	 */

	void reset( CBiteRnd& rnd, const int ParamCount )
	{
		SparseMul = 5;
		CountSp = Count * SparseMul;
		CountSp1 = CountSp - 1;

		const int NewCapacity = SlotCount * CountSp;

		if( NewCapacity > SelBufCapacity )
		{
			delete[] SelBuf;
			SelBufCapacity = NewCapacity;
			SelBuf = new int[ NewCapacity ];
		}

		int j;

		for( j = 0; j < SlotCount; j++ )
		{
			// Fill slot vector with replicas of possible choices.

			int* const sp = SelBuf + j * CountSp;
			Sels[ j ] = sp;
			int i;

			for( i = 0; i < Count; i++ )
			{
				int* const spo = sp + i * SparseMul;
				int k;

				for( k = 0; k < SparseMul; k++ )
				{
					spo[ k ] = i;
				}
			}

			// Randomized swap-mixing.

			for( i = 0; i < CountSp * 5; i++ )
			{
				const int i1 = rnd.getInt( CountSp );
				const int i2 = rnd.getInt( CountSp );

				const int t = sp[ i1 ];
				sp[ i1 ] = sp[ i2 ];
				sp[ i2 ] = t;
			}
		}

		select( rnd );
		IsSelected = false;
	}

	/**
	 * An auxiliary function that returns choice count.
	 */

	int getChoiceCount() const
	{
		return( Count );
	}

	/**
	 * This function should be called when a certain choice is successful.
	 * This function should only be called after a prior select() calls.
	 *
	 * @param rnd PRNG object. May not be used.
	 * @param v Selection increment value (success score), [0; 1].
	 */

	void incr( CBiteRnd& rnd, const double v = 1.0 )
	{
		const int dp = (int) ( -Selp * v * v );

		if( dp < 0 ) // Boost an efficient choice.
		{
			if( dp == -1 )
			{
				Sels[ Slot ][ Selp ] = Sels[ Slot ][ Selp - 1 ];
				Sels[ Slot ][ Selp - 1 ] = Sel;
			}
			else
			{
				const int np = Selp + dp;
				memmove( Sels[ Slot ] + np + 1, Sels[ Slot ] + np,
					-dp * sizeof( Sels[ Slot ][ 0 ]));

				Sels[ Slot ][ np ] = Sel;
			}
		}

		if( Slot > 0 )
		{
			int* const t = Sels[ Slot ];
			Sels[ Slot ] = Sels[ Slot - 1 ];
			Sels[ Slot - 1 ] = t;
		}
	}

	/**
	 * This function should be called when a certain choice is a failure.
	 * This function should only be called after a prior select() calls.
	 *
	 * @param rnd PRNG object. May not be used.
	 */

	void decr( CBiteRnd& rnd )
	{
		if( Selp < CountSp1 ) // Demote an inefficient choice.
		{
			Sels[ Slot ][ Selp ] = Sels[ Slot ][ Selp + 1 ];
			Sels[ Slot ][ Selp + 1 ] = Sel;
		}

		if( Slot < SlotCount - 1 )
		{
			int* const t = Sels[ Slot ];
			Sels[ Slot ] = Sels[ Slot + 1 ];
			Sels[ Slot + 1 ] = t;
		}
	}

	/**
	 * Function produces a random choice index based on the current *this
	 * object's state. Note that "select" functions can only be called once
	 * for a given *this object during the optimize() function call.
	 *
	 * @param rnd PRNG object.
	 */

	int select( CBiteRnd& rnd )
	{
		Slot = rnd.getPowInt( 1.5, SlotCount );
		Selp = rnd.getPowInt( SelPower, CountSp );

		Sel = Sels[ Slot ][ Selp ];
		IsSelected = true;

		return( Sel );
	}

	/**
	 * Function returns the latest made choice index.
	 */

	int getSel() const
	{
		return( Sel );
	}

	/**
	 * Function sets IsSelected value to "false".
	 */

	void unsetIsSelected()
	{
		IsSelected = false;
	}

	/**
	 * Function returns "true" if the selection was made since the recent
	 * increment or decrement.
	 */

	bool getIsSelected() const
	{
		return( IsSelected );
	}

protected:
	static const int SlotCount = 5; ///< The number of choice vectors in use.
	int Count; ///< The number of choices in use.
	int SparseMul; ///< Multiplier used to obtain an actual length of the
		///< choice vector. This multiplier replicates choices in the vector,
		///< increasing precision of the resulting PDF and its stability.
	int CountSp; ///< = Count * SparseMul. The actual length of the choice
		///< vector.
	int CountSp1; ///< = CountSp - 1.
	double SelPower; ///< Selector's power factor.
	int* Sels[ SlotCount ]; ///< Choice vectors.
	int* SelBuf; ///< A singular buffer for Sels vectors.
	int SelBufCapacity; ///< Capacity of SelBuf.
	int Sel; ///< The latest selected choice. Available only after the
		///< select() function calls.
	int Selp; ///< The index of the choice in the Sels vector.
	int Slot; ///< The current Sels vector, depending on incr/decr.
	bool IsSelected; ///< "True" if selection was recently made.
};

/**
 * Templated selector class, for convenient constructor's Count parameter
 * specification.
 *
 * @tparam tCount The number of possible choices, greater than 1.
 * @tparam spwr100 Selector's power factor multiplied by 100. At 100, all
 * choices are equally probable; at 150 the best choice may be selected up to
 * two times more often than any other choice.
 */

template< int tCount, int spwr100 = 150 >
class CBiteSel : public CBiteSelBase
{
public:
	CBiteSel()
		: CBiteSelBase( tCount, spwr100 )
	{
	}
};

/**
 * Class implements storage of population parameter vectors, costs, centroid,
 * and ordering. Note that the derived classes should delete allocated buffers
 * in their destructors, and not relying on the automatic deleteBuffers()
 * function call.
 *
 * @tparam ptype Parameter value storage type.
 */

template< typename ptype >
class CBitePop
{
public:
	CBitePop()
		: MantMult( IntMantMult )
		, MantMultI( 1.0 / IntMantMult )
		, ParamCount( 0 )
		, PopSize( 0 )
		, CnsCount( 0 )
		, ObjCount( 0 )
		, PopParamsBuf( NULL )
		, PopParams( NULL )
		, CentParams( NULL )
	{
	}

	CBitePop( const CBitePop& s )
		: PopParamsBuf( NULL )
		, PopParams( NULL )
		, CentParams( NULL )
	{
		initBuffers( s.ParamCount, s.PopSize, s.CnsCount, s.ObjCount );
		copy( s );
	}

	virtual ~CBitePop()
	{
		delete[] PopParamsBuf;
		delete[] PopParams;
		delete[] CentParams;
	}

	CBitePop& operator = ( const CBitePop& s )
	{
		copy( s );
		return( *this );
	}

	/**
	 * Function initializes all common buffers, and "PopSize" variables. This
	 * function should be called when population's dimensions were changed.
	 * This function calls the deleteBuffers() function to release any
	 * derived classes' allocated buffers. Allocates an additional vector for
	 * temporary use, which is at the same the last vector in the PopParams
	 * array. Derived classes should call this function of the base class.
	 *
	 * @param aParamCount New parameter count.
	 * @param aPopSize New population size. If <= 0, population buffers will
	 * not be allocated.
	 * @param aCnsCount New constraint value count.
	 * @param aObjCount New objective value count. If equals 0, a rank element
	 * will not be auto-added.
	 */

	virtual void initBuffers( const int aParamCount, const int aPopSize,
		const int aCnsCount = 0, const int aObjCount = 1 )
	{
		deleteBuffers();

		ParamCount = aParamCount;
		ParamCountI = 1.0 / ParamCount;
		PopSize = aPopSize;
		PopSize1 = aPopSize - 1;
		CnsCount = aCnsCount;
		ObjCount = aObjCount;
		NeedCentUpdate = false;
		CentLPC = calcLP1Coeff( PopSize );

		PopCnsOffs = aParamCount * sizeof( ptype );
		PopObjOffs = PopCnsOffs + aCnsCount * sizeof( double );
		PopRankOffs = PopObjOffs + aObjCount * sizeof( double );
		PopItemSize = PopRankOffs + ( aObjCount > 0 ? sizeof( double ) : 0 );

		PopParamsBuf = new uint8_t[( aPopSize + 1 ) * PopItemSize ];
		PopParams = new ptype*[ aPopSize + 1 ]; // Last element is temporary.
		CentParams = new ptype[ aParamCount ];

		int i;

		for( i = 0; i <= aPopSize; i++ )
		{
			PopParams[ i ] = (ptype*) ( PopParamsBuf + i * PopItemSize );
		}

		TmpParams = PopParams[ aPopSize ];
	}

	/**
	 * Function copies population from the specified source population. If
	 * *this population has a different size, or is uninitialized, it will
	 * be initialized to source's population size.
	 *
	 * @param s Source population to copy. Should be initalized.
	 */

	void copy( const CBitePop& s )
	{
		if( ParamCount != s.ParamCount || PopSize != s.PopSize ||
			CnsCount != s.CnsCount || ObjCount != s.ObjCount )
		{
			initBuffers( s.ParamCount, s.PopSize, s.CnsCount, s.ObjCount );
		}

		CurPopSize = s.CurPopSize;
		CurPopSizeI = s.CurPopSizeI;
		CurPopSize1 = s.CurPopSize1;
		CurPopPos = s.CurPopPos;
		NeedCentUpdate = s.NeedCentUpdate;
		CentLPC = s.CentLPC;

		int i;

		for( i = 0; i < PopSize; i++ )
		{
			memcpy( PopParams[ i ], s.PopParams[ i ], PopItemSize );
		}

		if( !NeedCentUpdate )
		{
			copyParams( CentParams, s.CentParams );
		}
	}

	/**
	 * Function recalculates centroid based on the current population size.
	 * The NeedCentUpdate variable can be checked if centroid update is
	 * needed. This function resets the NeedCentUpdate to "false". This
	 * function should only be called after the population is filled.
	 */

	void updateCentroid()
	{
		NeedCentUpdate = false;

		const int BatchCount = ( 1 << IntOverBits ) - 1;
		ptype* const cp = CentParams;
		const double cm = 1.0 / PopSize;
		ptype* const tp = TmpParams;
		int i;
		int j;

		if( PopSize <= BatchCount )
		{
			copyParams( tp, PopParams[ 0 ]);

			for( j = 1; j < PopSize; j++ )
			{
				const ptype* const p = PopParams[ j ];

				for( i = 0; i < ParamCount; i++ )
				{
					tp[ i ] += p[ i ];
				}
			}

			for( i = 0; i < ParamCount; i++ )
			{
				cp[ i ] = (ptype) ( tp[ i ] * cm );
			}
		}
		else
		{
			// Batched centroid calculation, for more precision and no integer
			// overflows.

			int pl = PopSize;
			j = 0;
			bool DoCopy = true;

			while( pl > 0 )
			{
				int c = ( pl > BatchCount ? BatchCount : pl );
				pl -= c;
				c--;

				copyParams( tp, PopParams[ j ]);

				while( c > 0 )
				{
					j++;
					const ptype* const p = PopParams[ j ];

					for( i = 0; i < ParamCount; i++ )
					{
						tp[ i ] += p[ i ];
					}

					c--;
				}

				if( DoCopy )
				{
					DoCopy = false;

					for( i = 0; i < ParamCount; i++ )
					{
						cp[ i ] = (ptype) ( tp[ i ] * cm );
					}
				}
				else
				{
					for( i = 0; i < ParamCount; i++ )
					{
						cp[ i ] += (ptype) ( tp[ i ] * cm );
					}
				}
			}
		}
	}

	/**
	 * Function returns pointer to the centroid vector. The NeedUpdateCent
	 * may needed to be checked, and if it is equal to "true", the
	 * updateCentroid() function should be called in advance.
	 */

	const ptype* getCentroid() const
	{
		return( CentParams );
	}

	/**
	 * Function returns population's parameter vector by the specified index
	 * from the ordered list.
	 *
	 * @param i Parameter vector index.
	 */

	ptype* getParamsOrdered( const int i ) const
	{
		return( PopParams[ i ]);
	}

	/**
	 * Function returns pointer to constraint values sub-array within the
	 * specified population vector.
	 *
	 * @param pp Population vector pointer (usually within PopParams).
	 */

	double* getCnsPtr( ptype* const pp ) const
	{
		return( (double*) ( (uintptr_t) pp + (uintptr_t) PopCnsOffs ));
	}

	/**
	 * Function returns pointer to objective values sub-array within the
	 * specified population vector.
	 *
	 * @param pp Population vector pointer (usually within PopParams).
	 */

	double* getObjPtr( ptype* const pp ) const
	{
		return( (double*) ( (uintptr_t) pp + (uintptr_t) PopObjOffs ));
	}

	/**
	 * Function returns pointer to the rank value sub-array within the
	 * specified population vector.
	 *
	 * @param pp Population vector pointer (usually within PopParams).
	 */

	double* getRankPtr( ptype* const pp ) const
	{
		return( (double*) ( (uintptr_t) pp + (uintptr_t) PopRankOffs ));
	}

	/**
	 * Function returns a pointer to array of population vector pointers,
	 * which are sorted in the ascending cost order.
	 */

	ptype** getPopParams() const
	{
		return( PopParams );
	}

	/**
	 * Function returns pointer to the next available parameter vector, at the
	 * initialization stage. When the population was filled, the function
	 * returns pointer to a temporary vector.
	 */

	ptype* getCurParams() const
	{
		return( PopParams[ CurPopPos ]);
	}

	/**
	 * Function returns current population size.
	 */

	int getCurPopSize() const
	{
		return( CurPopSize );
	}

	/**
	 * Function returns current population position.
	 */

	int getCurPopPos() const
	{
		return( CurPopPos );
	}

	/**
	 * Function resets the current population position to zero, and sets
	 * CurPopSize to PopSize. This function is usually called when the
	 * population needs to be completely changed. This function should be
	 * called before any updates to *this population (usually during
	 * optimizer's initialization).
	 */

	void resetCurPopPos()
	{
		CurPopSize = PopSize;
		CurPopSizeI = 1.0 / PopSize;
		CurPopSize1 = PopSize1;

		CurPopPos = 0;
		NeedCentUpdate = false;
		CentLPC = calcLP1Coeff( CurPopSize );
	}

	/**
	 * Function increases current population size, and updates the required
	 * variables. This function can only be called if CurPopSize is less than
	 * PopSize, and previously the whole population was filled. Note that
	 * if centroid update is needed, the NeedCentUpdate should be explicitly
	 * set to "true".
	 */

	void incrCurPopSize()
	{
		CurPopSize++;
		CurPopSizeI = 1.0 / CurPopSize;
		CurPopSize1++;
		CentLPC = calcLP1Coeff( CurPopSize );
	}

	/**
	 * Function decreases current population size, and updates the required
	 * variables. This function can only be called if CurPopSize is greater
	 * than 1, and the whole population was filled. Note that if centroid
	 * update is needed, the NeedCentUpdate should be explicitly set to
	 * "true".
	 */

	void decrCurPopSize()
	{
		CurPopSize--;
		CurPopSizeI = 1.0 / CurPopSize;
		CurPopSize1--;
		CentLPC = calcLP1Coeff( CurPopSize );
	}

	/**
	 * Function removes an item from the population.
	 *
	 * @param p Population's solution index, should be within
	 * the [0; CurPopPos - 1] range.
	 */

	void removeSol( const int p )
	{
		if( CurPopPos == 0 )
		{
			return;
		}

		const int ri = CurPopPos - 1;

		if( p < ri )
		{
			ptype** const pp = PopParams + p;
			ptype* const rp = *pp;
			memmove( pp, pp + 1, ( ri - p ) * sizeof( pp[ 0 ]));
			PopParams[ ri ] = rp;
		}

		CurPopPos--;
	}

	/**
	 * Function replaces the highest-cost previous solution, updates centroid.
	 * This function considers the value of the CurPopPos variable - if it is
	 * smaller than the PopSize, the new solution will be added to
	 * population without any checks.
	 *
	 * @param UpdCost Cost (rank) of the new solution. This value should be
	 * checked for NaN and fixed if needed.
	 * @param UpdParams New parameter values.
	 * @param DoUpdateCentroid "True" if centroid should be updated using
	 * running sum. This update is done for parallel populations.
	 * @param ReplaceThrN8 Solution's index threshold as
	 * CurPopSize mul ReplaceThrN8 div 8. If solution's index is below this
	 * value, and new cost is same as existing cost, a new solution will
	 * replace an existing solution in the population. Such replacing reduces
	 * diversity of competing same-cost best solutions, and usually improves
	 * convergence.
	 * @return Insertion position - greater or equal to PopSize, if the cost
	 * constraint was not met.
	 */

	int updatePop( double UpdCost, const ptype* const UpdParams,
		const bool DoUpdateCentroid = false, const int ReplaceThrN8 = 0 )
	{
		int ri; // Index of population vector to be replaced.

		if( CurPopPos < PopSize )
		{
			ri = CurPopPos;
		}
		else
		{
			ri = PopSize1;

			if( UpdCost > *getRankPtr( PopParams[ ri ]))
			{
				return( PopSize );
			}
		}

		// Perform binary search of a cost within the population.

		int p = 0;
		int i = ri;

		while( p < i )
		{
			const int mid = ( p + i ) >> 1;

			if( *getRankPtr( PopParams[ mid ]) >= UpdCost )
			{
				i = mid;
			}
			else
			{
				p = mid + 1;
			}
		}

		bool DoReplace = false;
		bool IsEqualCost = false;

		if( CurPopPos < PopSize )
		{
			CurPopPos++;
		}
		else
		{
			if( isEqual( UpdCost, *getRankPtr( PopParams[ p ])))
			{
				IsEqualCost = true;

				if( p != 0 && p < CurPopSize * ReplaceThrN8 / 8 &&
					isParams1FartherThan2( PopParams[ p ], UpdParams,
					PopParams[ 0 ]))
				{
					DoReplace = true;
				}
			}
		}

		ptype* rp;

		if( DoReplace )
		{
			rp = PopParams[ p ];
		}
		else
		{
			rp = PopParams[ ri ];

			const int mc = ri - p;
			ptype** const pp = PopParams + p;
			memmove( pp + 1, pp, mc * sizeof( pp[ 0 ]));

			*pp = rp;
		}

		*getObjPtr( rp ) = UpdCost;
		*getRankPtr( rp ) = UpdCost;

		if( rp != UpdParams )
		{
			if( DoUpdateCentroid && !NeedCentUpdate )
			{
				ptype* const cp = CentParams;
				const double lpc = CentLPC;

				for( i = 0; i < ParamCount; i++ )
				{
					cp[ i ] += (ptype) (( UpdParams[ i ] - cp[ i ]) * lpc );
					rp[ i ] = UpdParams[ i ];
				}
			}
			else
			{
				copyParams( rp, UpdParams );
				NeedCentUpdate = true;
			}
		}
		else
		{
			if( DoUpdateCentroid && !NeedCentUpdate )
			{
				ptype* const cp = CentParams;
				const double lpc = CentLPC;

				for( i = 0; i < ParamCount; i++ )
				{
					cp[ i ] += (ptype) (( UpdParams[ i ] - cp[ i ]) * lpc );
				}
			}
			else
			{
				NeedCentUpdate = true;
			}
		}

		return( IsEqualCost ? PopSize : p );
	}

protected:
	static const int IntOverBits = ( sizeof( ptype ) > 4 ? 5 : 3 ); ///< The
		///< number of bits of precision required for integer centroid
		///< calculation and overflows.
	static const int IntMantBits = sizeof( ptype ) * 8 - 1 - IntOverBits; ///<
		///< Mantissa size of the integer parameter values (higher by 1 bit in
		///< practice for real value 1.0). Accounts for a sign bit, and
		///< possible accumulation overflows.
	static const int64_t IntMantMult = 1LL << IntMantBits; ///< Mantissa
		///< multiplier.
	static const int64_t IntMantMultM = -IntMantMult; ///< Negative
		///< IntMantMult.
	static const int64_t IntMantMult2 = ( IntMantMult << 1 ); ///< =
		///< IntMantMult * 2.
	static const int64_t IntMantMask = IntMantMult - 1; ///< Mask that
		///< corresponds to mantissa.
	double MantMult; ///< = IntMantMult.
	double MantMultI; ///< = 1.0 / IntMantMult.

	int ParamCount; ///< The total number of internal parameter values in use.
	double ParamCountI; ///< = 1.0 / ParamCount.
	int PopSize; ///< The size of population in use (maximal).
	int PopSize1; ///< = PopSize - 1.
	int CurPopSize; ///< Current population size.
	int CurPopSize1; ///< = CurPopSize - 1.
	double CurPopSizeI; ///< = 1.0 / CurPopSize.
	int CurPopPos; ///< Current population position, for initial population
		///< update. This variable should be initialized by the optimizer.
	int CnsCount; ///< The number of constraints per solution.
	int ObjCount; ///< The number of objectives per solution.
	uint8_t* PopParamsBuf; ///< Buffer for all PopParams vectors.
	ptype** PopParams; ///< Population parameter vectors. Always kept sorted
		///< in ascending cost order. Each vector represents a complex item,
		///< with additional data stored after parameter values (see Offs
		///< constants).
	size_t PopCnsOffs; ///< Byte offset to the constraint values within
		///< a population item.
	size_t PopObjOffs; ///< Byte offset to the objective values within a
		///< population item.
	size_t PopRankOffs; ///< Byte offset to the rank value within a population
		///< item.
	size_t PopItemSize; ///< Size in bytes of population item.
	ptype* CentParams; ///< Centroid of the parameter vectors.
	bool NeedCentUpdate; ///< "True" if centroid update is needed.
	double CentLPC; ///< Centroid averaging filter coefficient.
	ptype* TmpParams; ///< Temporary parameter vector, points to the last
		///< element of the PopParams array.

	/**
	 * Function deletes buffers previously allocated via the initBuffers()
	 * function. Derived classes should call this function of the base class.
	 */

	virtual void deleteBuffers()
	{
		delete[] PopParamsBuf;
		delete[] PopParams;
		delete[] CentParams;
	}

	/**
	 * An aux function that resets a parameter vector to zero values.
	 *
	 * @param[out] dst Destination vector.
	 */

	void zeroParams( ptype* const dst ) const
	{
		memset( dst, 0, ParamCount * sizeof( dst[ 0 ]));
	}

	/**
	 * An aux function that copies a parameter vector.
	 *
	 * @param[out] dst Destination vector.
	 * @param[in] src Source vector.
	 */

	void copyParams( ptype* const dst, const ptype* const src ) const
	{
		memcpy( dst, src, ParamCount * sizeof( dst[ 0 ]));
	}

	/**
	 * An aux function that copies a real solution vector.
	 *
	 * @param[out] dst Destination vector.
	 * @param[in] src Source vector.
	 */

	void copyValues( double* const dst, const double* const src ) const
	{
		memcpy( dst, src, ParamCount * sizeof( dst[ 0 ]));
	}

	/**
	 * Function compares two floating-point values for equality.
	 *
	 * @param a Value 1.
	 * @param b Value 2.
	 * @param etol Value equality tolerance.
	 * @return "True" if values are equal.
	 */

	static bool isEqual( const double a, const double b,
		const double etol = 0x1p-52 )
	{
		const double d = fabs( b - a );

		if( d == 0.0 )
		{
			return( true );
		}

		const double s = fabs( b ) + fabs( a );

		if( d < s * etol )
		{
			return( true );
		}

		return( false );
	}

	/**
	 * Function compares distances of parameter vectors 1 and 2 to a reference
	 * vector, and returns "true" if distance 1 is greater than 2.
	 *
	 * @param p1 Parameter vector 1.
	 * @param p2 Parameter vector 2.
	 * @param RefParams Reference parameter vector.
	 */

	bool isParams1FartherThan2( const ptype* const p1, const ptype* const p2,
		const ptype* const RefParams ) const
	{
		double s0 = 0.0;
		double s1 = 0.0;
		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			const ptype v = RefParams[ i ];
			const double d0 = (double) ( p1[ i ] - v );
			const double d1 = (double) ( p2[ i ] - v );
			s0 += d0 * d0;
			s1 += d1 * d1;
		}

		return( s0 > s1 );
	}

	/**
	 * Function wraps the specified parameter value so that it stays in the
	 * [0.0; 1.0] range (including in integer range), by wrapping it over the
	 * boundaries using random operator. This operation improves convergence
	 * in comparison to clamping.
	 *
	 * @param rnd PRNG object.
	 * @param v Parameter value to wrap.
	 * @return Wrapped parameter value.
	 */

	static ptype wrapParam( CBiteRnd& rnd, const ptype v )
	{
		if( (ptype) 0.25 == 0 )
		{
			if( v < 0 )
			{
				if( v > IntMantMultM )
				{
					return( (ptype) ( rnd.get() * -v ));
				}

				return( (ptype) ( rnd.getRaw() & IntMantMask ));
			}

			if( v > IntMantMult )
			{
				if( v < IntMantMult2 )
				{
					return( (ptype) ( IntMantMult -
						rnd.get() * ( v - IntMantMult )));
				}

				return( (ptype) ( rnd.getRaw() & IntMantMask ));
			}

			return( v );
		}
		else
		{
			if( v < 0.0 )
			{
				if( v > -1.0 )
				{
					return( (ptype) ( rnd.get() * -v ));
				}

				return( (ptype) rnd.get() );
			}

			if( v > 1.0 )
			{
				if( v < 2.0 )
				{
					return( (ptype) ( 1.0 - rnd.get() * ( v - 1.0 )));
				}

				return( (ptype) rnd.get() );
			}

			return( v );
		}
	}

	/**
	 * Function generates a Gaussian-distributed pseudo-random number, in
	 * integer scale, with the specified mean and std.dev.
	 *
	 * @param rnd PRNG object.
	 * @param sd Standard deviation multiplier.
	 * @param meanInt Mean value, in integer scale.
	 */

	static ptype getGaussianInt( CBiteRnd& rnd, const double sd,
		const ptype meanInt )
	{
		while( true )
		{
			const double r = rnd.getGaussian() * sd;

			if( r > -8.0 && r < 8.0 )
			{
				return( (ptype) ( r * IntMantMult + meanInt ));
			}
		}
	}

	/**
	 * Function calculates averaging coefficient of a "leaky integrator"
	 * 1st order low-pass filter.
	 *
	 * @param Count An approximate number of samples to average.
	 */

	static double calcLP1Coeff( const double Count )
	{
		const double theta = 2.8 / Count;
		const double costheta2 = 2.0 - cos( theta );

		return( 1.0 - ( costheta2 - sqrt( costheta2 * costheta2 - 1.0 )));
	}
};

/**
 * Population class that embeds a dynamically-allocated parallel population
 * objects.
 *
 * @tparam ptype Parameter value storage type.
 */

template< typename ptype >
class CBiteParPops : virtual public CBitePop< ptype >
{
public:
	CBiteParPops()
		: ParPops( NULL )
		, ParPopCount( 0 )
		, ParValues( NULL )
	{
	}

	virtual ~CBiteParPops()
	{
		int i;

		for( i = 0; i < ParPopCount; i++ )
		{
			delete ParPops[ i ];
		}

		delete[] ParPops;
		delete[] ParValues;
	}

protected:
	using CBitePop< ptype > :: ParamCount;
	using CBitePop< ptype > :: PopSize;

	CBitePop< ptype >** ParPops; ///< Parallel population "orbiting" `this`
		///< population.
	int ParPopCount; ///< Parallel population count. This variable should only
		///< be changed via the setParPopCount() function.
	double* ParValues; ///< Temporary value buffer, length equals ParPopCount.

	/**
	 * Function changes the parallel population count, and reallocates
	 * buffers.
	 *
	 * @param NewCount New parallel population count, >= 0.
	 */

	void setParPopCount( const int NewCount )
	{
		while( ParPopCount > NewCount )
		{
			ParPopCount--;
			delete ParPops[ ParPopCount ];
		}

		if( NewCount > ParPopCount )
		{
			CBitePop< ptype >** NewParPops =
				new CBitePop< ptype >*[ NewCount ];

			if( ParPopCount != 0 )
			{
				memcpy( NewParPops, ParPops, sizeof( ParPops[ 0 ]) *
					ParPopCount );
			}

			delete[] ParPops;
			ParPops = NewParPops;

			while( ParPopCount < NewCount )
			{
				ParPops[ ParPopCount ] = new CBitePop< ptype >();
				ParPopCount++;
			}

			delete[] ParValues;
			ParValues = new double[ ParPopCount ];
		}
	}

	/**
	 * Function calculates distances of the specified parameter vector to
	 * parallel populations' centroids.
	 *
	 * @param Params Parameter vector.
	 * @param[out] s Resulting squared distances vector.
	 */

	void calcCentroidDists( const ptype* const Params, double* s ) const
	{
		int k = 0;
		int i;

		while( k < ParPopCount )
		{
			int pc = ParPopCount - k;

			if( pc > 4 )
			{
				pc = 4;
			}

			if( pc == 4 )
			{
				const ptype* const c0 = ParPops[ k ] -> getCentroid();
				const ptype* const c1 = ParPops[ k + 1 ] -> getCentroid();
				const ptype* const c2 = ParPops[ k + 2 ] -> getCentroid();
				const ptype* const c3 = ParPops[ k + 3 ] -> getCentroid();
				double s0 = 0.0;
				double s1 = 0.0;
				double s2 = 0.0;
				double s3 = 0.0;

				for( i = 0; i < ParamCount; i++ )
				{
					const ptype v = Params[ i ];
					const double d0 = (double) ( c0[ i ] - v );
					const double d1 = (double) ( c1[ i ] - v );
					s0 += d0 * d0;
					s1 += d1 * d1;

					const double d2 = (double) ( c2[ i ] - v );
					const double d3 = (double) ( c3[ i ] - v );
					s2 += d2 * d2;
					s3 += d3 * d3;
				}

				s[ 0 ] = s0;
				s[ 1 ] = s1;
				s[ 2 ] = s2;
				s[ 3 ] = s3;
			}
			else
			if( pc == 3 )
			{
				const ptype* const c0 = ParPops[ k ] -> getCentroid();
				const ptype* const c1 = ParPops[ k + 1 ] -> getCentroid();
				const ptype* const c2 = ParPops[ k + 2 ] -> getCentroid();
				double s0 = 0.0;
				double s1 = 0.0;
				double s2 = 0.0;

				for( i = 0; i < ParamCount; i++ )
				{
					const ptype v = Params[ i ];
					const double d0 = (double) ( c0[ i ] - v );
					const double d1 = (double) ( c1[ i ] - v );
					s0 += d0 * d0;
					s1 += d1 * d1;

					const double d2 = (double) ( c2[ i ] - v );
					s2 += d2 * d2;
				}

				s[ 0 ] = s0;
				s[ 1 ] = s1;
				s[ 2 ] = s2;
			}
			else
			if( pc == 2 )
			{
				const ptype* const c0 = ParPops[ k ] -> getCentroid();
				const ptype* const c1 = ParPops[ k + 1 ] -> getCentroid();
				double s0 = 0.0;
				double s1 = 0.0;

				for( i = 0; i < ParamCount; i++ )
				{
					const ptype v = Params[ i ];
					const double d0 = (double) ( c0[ i ] - v );
					const double d1 = (double) ( c1[ i ] - v );
					s0 += d0 * d0;
					s1 += d1 * d1;
				}

				s[ 0 ] = s0;
				s[ 1 ] = s1;
			}
			else
			{
				const ptype* const c0 = ParPops[ k ] -> getCentroid();
				double s0 = 0.0;

				for( i = 0; i < ParamCount; i++ )
				{
					const double d0 = (double) ( c0[ i ] - Params[ i ] );
					s0 += d0 * d0;
				}

				s[ 0 ] = s0;
			}

			k += pc;
			s += pc;
		}
	}

	/**
	 * Function returns index of the parallel population that is most close
	 * to the specified parameter vector. Function returns -1 if the cost
	 * constraint is not met in all parallel populations.
	 *
	 * @param Cost Cost of parameter vector, used to filter considered
	 * parallel population pool.
	 * @param Params Parameter vector.
	 */

	int getMinDistParPop( const double Cost, const ptype* const Params ) const
	{
		double* const s = ParValues;
		calcCentroidDists( Params, ParValues );

		int pp = 0;
		double d = s[ pp ];
		int i;

		for( i = 1; i < ParPopCount; i++ )
		{
			if( s[ i ] <= d )
			{
				pp = i;
				d = s[ i ];
			}
		}

		return( pp );
	}
};

/**
 * Base virtual abstract class that defines common optimizer interfacing
 * functions.
 */

class CBiteOptInterface
{
public:
	CBiteOptInterface()
	{
	}

	virtual ~CBiteOptInterface()
	{
	}

	/**
	 * @return Best parameter vector.
	 */

	virtual const double* getBestParams() const = 0;

	/**
	 * @return Cost of the best parameter vector.
	 */

	virtual double getBestCost() const = 0;

	/**
	 * @return The pointer to array with cost(s) (objectives) obtained at the
	 * latest optcost() function call. This pointer should be valid until
	 * updateDims() or initBuffers() function calls.
	 */

	virtual const double* getLastCosts() const = 0;

	/**
	 * @return The pointer to array with parameter values obtained at the
	 * latest `optcost()` function call. This pointer should be valid until
	 * updateDims() or initBuffers() function calls.
	 */

	virtual const double* getLastValues() const = 0;

	/**
	 * Virtual function that should fill minimal parameter value vector.
	 *
	 * @param[out] p Minimal value vector.
	 */

	virtual void getMinValues( double* const p ) const = 0;

	/**
	 * Virtual function that should fill maximal parameter value vector.
	 *
	 * @param[out] p Maximal value vector.
	 */

	virtual void getMaxValues( double* const p ) const = 0;

	/**
	 * Virtual function (objective function) that should calculate parameter
	 * vector's optimization cost.
	 *
	 * @param p Parameter vector to evaluate.
	 * @return Optimized cost.
	 */

	virtual double optcost( const double* const p ) = 0;
};

/**
 * The base class for optimizers of the "biteopt" project.
 *
 * @tparam ptype Parameter value storage type.
 */

template< typename ptype >
class CBiteOptBase : public CBiteOptInterface,
	virtual protected CBiteParPops< ptype >
{
private:
	CBiteOptBase( const CBiteOptBase& )
	{
		// Copy-construction unsupported.
	}

	CBiteOptBase& operator = ( const CBiteOptBase& )
	{
		// Copying unsupported.
		return( *this );
	}

public:
	static const int MaxSelCount = 64; ///< The maximal number of selectors
		///< that can be added to *this object (for static arrays).

	CBiteOptBase()
		: MinValues( NULL )
		, MaxValues( NULL )
		, DiffValues( NULL )
		, DiffValuesI( NULL )
		, StartParams( NULL )
		, BestValues( NULL )
		, NewValues( NULL )
		, SelCount( 0 )
	{
	}

	virtual ~CBiteOptBase()
	{
		delete[] MinValues;
		delete[] MaxValues;
		delete[] DiffValues;
		delete[] DiffValuesI;
		delete[] StartParams;
		delete[] BestValues;
		delete[] NewValues;
	}

	virtual const double* getBestParams() const
	{
		return( BestValues );
	}

	virtual double getBestCost() const
	{
		return( BestCost );
	}

	virtual const double* getLastCosts() const
	{
		return( LastCosts );
	}

	virtual const double* getLastValues() const
	{
		return( LastValues );
	}

	/**
	 * Function returns a pointer to an array of selectors in use.
	 */

	CBiteSelBase** getSels()
	{
		return( Sels );
	}

	/**
	 * Function returns a pointer to an array of selector names.
	 */

	const char** getSelNames() const
	{
		return( (const char**) SelNames );
	}

	/**
	 * Function returns the number of selectors in use.
	 */

	int getSelCount() const
	{
		return( SelCount );
	}

	/**
	 * Returns the number of iterations without improvement.
	 */

	int getStallCount() const
	{
		return( StallCount );
	}

protected:
	using CBiteParPops< ptype > :: IntMantMult;
	using CBiteParPops< ptype > :: MantMult;
	using CBiteParPops< ptype > :: MantMultI;
	using CBiteParPops< ptype > :: ParamCount;
	using CBiteParPops< ptype > :: CurPopPos;
	using CBiteParPops< ptype > :: resetCurPopPos;
	using CBiteParPops< ptype > :: copyValues;
	using CBiteParPops< ptype > :: wrapParam;
	using CBiteParPops< ptype > :: getGaussianInt;

	double* MinValues; ///< Minimal parameter values.
	double* MaxValues; ///< Maximal parameter values.
	double* DiffValues; ///< Difference between maximal and minimal parameter
		///< values.
	double* DiffValuesI; ///< Inverse DiffValues.
	ptype* StartParams; ///< Starting parameter values.
	bool UseStartParams; ///< "True" if StartParams are available.
	double StartSD; ///< Starting standard deviation.
	double* BestValues; ///< Best parameter vector.
	double BestCost; ///< Cost of the best parameter vector.
	double NewCosts[ 1 ]; ///< Temporary buffer to contain objective function
		///< value(s). Pointer to this buffer is retained in init() call.
	double* NewValues; ///< New parameter values buffer, with real values,
		///< can be also used as a temporary buffer. Pointer to this buffer
		///< is retained in init() call.
	const double* LastCosts; ///< Cost(s) of the latest optcost() call. Points
		///< to NewCosts by default.
	const double* LastValues; ///< Parameter values of the latest optcost()
		///< call. Points to NewValues by default.
	bool DoInitEvals; ///< "True" if initial evaluations should be performed.
	int StallCount; ///< The number of iterations without improvement.
	double HiBound; ///< Higher cost bound, for StallCount estimation. May not
		///< be used by the optimizer.
	double AvgCost; ///< Average cost in the latest batch. May not be used by
		///< the optimizer.
	CBiteSelBase* Sels[ MaxSelCount ]; ///< Pointers to selector objects, for
		///< indexed access in some cases.
	const char* SelNames[ MaxSelCount ]; ///< Selector names.
	int SelCount; ///< The number of selectors in use.
	static const int MaxApplySels = 32; ///< The maximal number of selections
		///< that can be used during a single optimize() function call.
	CBiteSelBase* ApplySels[ MaxApplySels ]; ///< Selectors that were used in
		///< select() function calls during the optimize() function call.
	int ApplySelsCount; ///< The number of select() calls performed during the
		///< optimize() function call.

	virtual void initBuffers( const int aParamCount, const int aPopSize,
		const int aCnsCount = 0, const int aObjCount = 1 )
	{
		CBiteParPops< ptype > :: initBuffers( aParamCount, aPopSize,
			aCnsCount, aObjCount );

		MinValues = new double[ ParamCount ];
		MaxValues = new double[ ParamCount ];
		DiffValues = new double[ ParamCount ];
		DiffValuesI = new double[ ParamCount ];
		StartParams = new ptype[ ParamCount ];
		BestValues = new double[ ParamCount ];
		NewValues = new double[ ParamCount ];
	}

	virtual void deleteBuffers()
	{
		CBiteParPops< ptype > :: deleteBuffers();

		delete[] MinValues;
		delete[] MaxValues;
		delete[] DiffValues;
		delete[] DiffValuesI;
		delete[] StartParams;
		delete[] BestValues;
		delete[] NewValues;
	}

	/**
	 * Function initializes or resets common variables used by optimizers, to
	 * their default values, including value bounds, registered selectors,
	 * calls the resetCurPopPos() and updateDiffValues() functions. This
	 * function is usually called in the init() function of an optimizer.
	 * Sets UseStartParams to "false".
	 */

	void initCommonVars( CBiteRnd& rnd )
	{
		getMinValues( MinValues );
		getMaxValues( MaxValues );
		updateDiffValues();

		resetCurPopPos();

		UseStartParams = false;
		StartSD = 0.25;
		BestCost = 1e300;
		LastCosts = NewCosts;
		LastValues = NewValues;
		DoInitEvals = true;
		StallCount = 0;
		HiBound = 1e300;
		AvgCost = 0.0;
		ApplySelsCount = 0;

		int i;

		for( i = 0; i < SelCount; i++ )
		{
			Sels[ i ] -> reset( rnd, ParamCount );
		}
	}

	/**
	 * Function updates values in the DiffValues array, based on values in the
	 * MinValues and MaxValues arrays.
	 */

	void updateDiffValues()
	{
		int i;

		if( (ptype) 0.25 == 0 )
		{
			for( i = 0; i < ParamCount; i++ )
			{
				const double d = MaxValues[ i ] - MinValues[ i ];

				DiffValues[ i ] = d * MantMultI;
				DiffValuesI[ i ] = MantMult / d;
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				DiffValues[ i ] = MaxValues[ i ] - MinValues[ i ];
				DiffValuesI[ i ] = 1.0 / DiffValues[ i ];
			}
		}
	}

	/**
	 * Function updates BestCost value and BestValues array, if the specified
	 * NewCost is better.
	 *
	 * @param UpdCost New solution's cost.
	 * @param UpdValues New solution's values. The values should be in the
	 * "real" value range.
	 * @param p New solution's position within population. If <0, position
	 * unknown, and the UpdCost should be evaluated.
	 */

	void updateBestCost( const double UpdCost, const double* const UpdValues,
		int p = -1 )
	{
		if( p < 0 )
		{
			if( UpdCost <= BestCost )
			{
				p = 0;
			}
		}

		if( p == 0 )
		{
			BestCost = UpdCost;

			copyValues( BestValues, UpdValues );
		}
	}

	/**
	 * Function checks the supplied cost value for NaN and returns a "fixed"
	 * value.
	 *
	 * @param v Cost value to check.
	 */

	static double fixCostNaN( const double v )
	{
		return( v != v ? 1e300 : v );
	}

	/**
	 * Function calculates population size based on the specified parameter
	 * count, for use with BiteOpt optimizers.
	 *
	 * @param aParamCount The number of parameters in use.
	 */

	static int calcPopSizeBiteOpt( const int aParamCount )
	{
		const double cx = tanh( 0.008 * aParamCount );
		const double psl = 10 + aParamCount * 3;
		const double psh = 11.0 * sqrt( (double) aParamCount );

		return( (int) ( psl * ( 1.0 - cx ) + psh * cx + 0.5 ));
	}

	/**
	 * Function returns specified parameter's value taking into account
	 * minimal and maximal value range.
	 *
	 * @param NormParams Parameter vector of interest, in normalized scale.
	 * @param i Parameter index.
	 */

	double getRealValue( const ptype* const NormParams, const int i ) const
	{
		return( MinValues[ i ] + DiffValues[ i ] * NormParams[ i ]);
	}

	/**
	 * Function wraps the specified parameter value so that it stays in the
	 * [MinValue; MaxValue] real range, by wrapping it over the boundaries
	 * using random operator. This operation improves convergence in
	 * comparison to clamping.
	 *
	 * @param rnd PRNG object.
	 * @param v Parameter value to wrap.
	 * @param i Parameter index.
	 * @return Wrapped parameter value.
	 */

	double wrapParamReal( CBiteRnd& rnd, const double v, const int i ) const
	{
		const double minv = MinValues[ i ];

		if( v < minv )
		{
			const double dv = DiffValues[ i ];

			if( v > minv - dv )
			{
				return( minv + rnd.get() * ( minv - v ));
			}

			return( minv + rnd.get() * dv );
		}

		const double maxv = MaxValues[ i ];

		if( v > maxv )
		{
			const double dv = DiffValues[ i ];

			if( v < maxv + dv )
			{
				return( maxv - rnd.get() * ( v - maxv ));
			}

			return( maxv - rnd.get() * dv );
		}

		return( v );
	}

	/**
	 * Function sets up starting parameters for normalized space.
	 *
	 * @param InitParams Starting parameters in real space. If NULL, the
	 * function does nothing.
	 */

	void setStartParams( const double* const InitParams )
	{
		if( InitParams == NULL )
		{
			return;
		}

		int i;

		for( i = 0; i < ParamCount; i++ )
		{
			StartParams[ i ] = (ptype) (( InitParams[ i ] - MinValues[ i ]) /
				DiffValues[ i ]);
		}

		UseStartParams = true;
	}

	/**
	 * Function generates a random initial solution with Gaussian sampling
	 * using `StartSD` variance, in normalized parameter space. The solution
	 * will be centered around the `StartParams` point if the `UseStartParams`
	 * is `true`. The function also fills the `NewValues` array.
	 *
	 * @param rnd PRNG object.
	 * @param[out] Params Resulting parameter vector.
	 */

	void genInitParams( CBiteRnd& rnd, ptype* const Params ) const
	{
		int i;

		if( UseStartParams )
		{
			if( CurPopPos == 0 )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = wrapParam( rnd, StartParams[ i ]);
					NewValues[ i ] = getRealValue( Params, i );
				}
			}
			else
			{
				if( (ptype) 0.25 == 0 )
				{
					for( i = 0; i < ParamCount; i++ )
					{
						Params[ i ] = wrapParam( rnd, getGaussianInt( rnd,
							StartSD, StartParams[ i ]));

						NewValues[ i ] = getRealValue( Params, i );
					}
				}
				else
				{
					for( i = 0; i < ParamCount; i++ )
					{
						Params[ i ] = wrapParam( rnd, rnd.getGaussian() *
							StartSD + StartParams[ i ]);

						NewValues[ i ] = getRealValue( Params, i );
					}
				}
			}
		}
		else
		{
			if( (ptype) 0.25 == 0 )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = wrapParam( rnd, getGaussianInt( rnd,
						StartSD, IntMantMult >> 1 ));

					NewValues[ i ] = getRealValue( Params, i );
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = wrapParam( rnd, rnd.getGaussian() *
						StartSD + 0.5 );

					NewValues[ i ] = getRealValue( Params, i );
				}
			}
		}
	}

	/**
	 * Function generates a random initial solution with Gaussian sampling
	 * using `StartSD` variance, in real parameter space. The solution will be
	 * centered around the `StartParams` point if the `UseStartParams` is
	 * `true`. The function does not fill the `NewValues` array.
	 *
	 * @param rnd PRNG object.
	 * @param[out] Params Resulting parameter vector.
	 */

	void genInitParamsReal( CBiteRnd& rnd, double* const Params ) const
	{
		int i;

		if( UseStartParams )
		{
			if( CurPopPos == 0 )
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = wrapParamReal( rnd, StartParams[ i ], i );
				}
			}
			else
			{
				for( i = 0; i < ParamCount; i++ )
				{
					Params[ i ] = wrapParamReal( rnd, rnd.getGaussian() *
						StartSD * DiffValues[ i ] + StartParams[ i ], i );
				}
			}
		}
		else
		{
			for( i = 0; i < ParamCount; i++ )
			{
				Params[ i ] = wrapParamReal( rnd, rnd.getGaussian() *
					StartSD * DiffValues[ i ] +
					MinValues[ i ] + DiffValues[ i ] * 0.5, i );
			}
		}
	}

	/**
	 * Function adds a selector to the Sels list.
	 *
	 * @param s Selector object to add.
	 * @param sname Selector's name, should be a static constant.
	 */

	void addSel( CBiteSelBase& s, const char* const sname )
	{
		Sels[ SelCount ] = &s;
		SelNames[ SelCount ] = sname;
		SelCount++;
	}

	/**
	 * Function performs choice selection based on the specified selector, and
	 * adds the selector to apply list.
	 *
	 * @param Sel Selector object.
	 * @param rnd PRNG object.
	 */

	template< class T >
	int select( T& Sel, CBiteRnd& rnd )
	{
		ApplySels[ ApplySelsCount ] = &Sel;
		ApplySelsCount++;

		return( Sel.select( rnd ));
	}

	/**
	 * Function applies selector increments on optimization success.
	 *
	 * @param rnd PRNG object.
	 * @param v Increment value, [0; 1].
	 */

	void applySelsIncr( CBiteRnd& rnd, const double v = 1.0 )
	{
		const int c = ApplySelsCount;
		ApplySelsCount = 0;

		int i;

		for( i = 0; i < c; i++ )
		{
			ApplySels[ i ] -> incr( rnd, v );
		}
	}

	/**
	 * Function applies selector decrements on optimization fail.
	 *
	 * @param rnd PRNG object.
	 */

	void applySelsDecr( CBiteRnd& rnd )
	{
		const int c = ApplySelsCount;
		ApplySelsCount = 0;

		int i;

		for( i = 0; i < c; i++ )
		{
			ApplySels[ i ] -> decr( rnd );
		}
	}
};

/**
 * Wrapping class for parallel optimizers, calls owner's cost function.
 *
 * @tparam The actual optimizer class.
 */

template< class T >
class CBiteOptOwned : public T
{
public:
	CBiteOptOwned( CBiteOptInterface* const aOwner )
		: Owner( aOwner )
	{
	}

	virtual void getMinValues( double* const p ) const
	{
		Owner -> getMinValues( p );
	}

	virtual void getMaxValues( double* const p ) const
	{
		Owner -> getMaxValues( p );
	}

	virtual double optcost( const double* const p )
	{
		return( Owner -> optcost( p ));
	}

protected:
	CBiteOptInterface* Owner; ///< Owner object.
};

#endif // BITEAUX_INCLUDED

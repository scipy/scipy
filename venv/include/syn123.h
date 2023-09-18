/*
	syn123: some audio signal synthesis and format conversion

	copyright 2017-2020 by the mpg123 project,
	free software under the terms of the LGPL 2.1
	see COPYING and AUTHORS files in distribution or http://mpg123.org

	initially written by Thomas Orgis
*/

#ifndef SYN123_H
#define SYN123_H

/** \file syn123.h The header file for the libsyn123 library. */

/* Common audio encoding specification. */
#include <fmt123.h>

/** A macro to check at compile time which set of API functions to expect.
 * This should be incremented at least each time a new symbol is added
 * to the header.
 */
#ifndef SYN123_API_VERSION
#define SYN123_API_VERSION 1
#endif

#ifndef MPG123_EXPORT
/** Defines needed for MS Visual Studio(tm) DLL builds.
 * Every public function must be prefixed with MPG123_EXPORT. When building 
 * the DLL ensure to define BUILD_MPG123_DLL. This makes the function accessible
 * for clients and includes it in the import library which is created together
 * with the DLL. When consuming the DLL ensure to define LINK_MPG123_DLL which
 * imports the functions from the DLL. 
 */
#ifdef BUILD_MPG123_DLL
/* The dll exports. */
#define MPG123_EXPORT __declspec(dllexport)
#else
#ifdef LINK_MPG123_DLL
/* The exe imports. */
#define MPG123_EXPORT __declspec(dllimport)
#else
/* Nothing on normal/UNIX builds */
#define MPG123_EXPORT
#endif
#endif
#endif

/** Support the restrict keyword for handed-in pointers. Defined to
    'restrict' if available.
 */
#ifndef MPG123_RESTRICT
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define MPG123_RESTRICT restrict
#else
#define MPG123_RESTRICT
#endif
#endif

/* Enable use of this file without configure. */
/* Also make sure to define _FILE_OFFSET_BITS, too. */
#ifndef MPG123_NO_CONFIGURE
#include <stdlib.h>
#include <sys/types.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** \defgroup syn123_api syn123 library API
 *  I wanted to create some test signals in different encodings for testing
 *  libout123. Also, they shall serve to verify the decoding of libmpg123
 *  in automated testing. Then, the awful drop-sample resampling inside the
 *  latter needed replacement. As digital filtering is part of the resampler
 *  dance, a chain of those can also be configured and applied. I hope I can
 *  avoid adding yet more functionality. It's a little utility library to
 *  accompany libmpg123 and libout123, OK?
 *
 *  This is what libsyn123 offers:
 *
 *  - signal generation (mix of differing wave shapes, noise, a simulated
 *    Geiger counter for fun)
 *  - format conversion and channel mixing and amplification, with hard and
 *    soft clipping, also optional dithering
 *  - near-zero-latency good-enough quite-fast resampling
 *  - applying digital filters with user-provided coefficients
 *
 *  The usage model for signal generation is this:
 *
 *  1. Create handle with desired output format.
 *  2. Set up synthesis mode with parameters.
 *  3. Repeatedly extract buffers with PCM samples.
 *
 *  If your hardware is slow on floating point operations, you may benefit
 *  from the period buffer in the handle that only needs actual computation
 *  once in the setup function. The frequencies of a wave mix may be fudged
 *  a bit to make for a configured period size.
 *
 *  The usage model for resampling is this:
 *
 *  1. Create handle with any format or re-use a signal generation handle.
 *  2. Set up resampling (rate ratio, channels, dirty mode).
 *  3. Use predictor functions to work out matching sizes of input and
 *     output buffers.
 *  4. Call the resampler on chunks of input, immediately producing matching
 *     chunks of output.
 *
 *  The resampler works on 32 bit float, exclusively. This is what is required
 *  and appropriate for the provided quality.
 *
 *  The usage model for filtering is this:
 *
 *  1. Create or re-use a handle.
 *  2. Set up the filter chain for single or double precision floating point.
 *  3. Apply the filter to chunks of input in succession.
 *
 *  Computations are generally done in floating point, either 32 bit or 64
 *  bit (single or double precision). The functions for encoding conversion,
 *  (de-)interleaving, and interleaved mixing can work without a handle and
 *  only use the buffers you hand in. Some support using a syn123 handle to
 *  provide the temporary buffers for encoding conversion on the fly to apply
 *  flotating-point processing to integer encodings.
 *
 *  Only the functions that are able to return a success code do check
 *  arguments for obvious trouble like NULL pointers. You are supposed to
 *  act responsibly when calling.
 *
 *  The size of buffers is either counted in bytes or samples, which,
 *  depending on the context, refer to individual PCM samples or to what is
 *  more strictly called a frame (one sample for each channel counted
 *  together).
 @{
 */

/** Opaque structure for the libsyn123 handle.
 *
 *  Simple context-free API functions do not need a handle, while
 *  others require it. Those that require it want it as first argument.
 *  Functions taking a handle as last argument after others make optional
 *  use of it (if non-NULL) to enable advanced functionality like
 *  on-the-fly encoding conversion that needs temporary storage.
 */
struct syn123_struct;
/** Typedef shortcut as preferrend name for the handle type. */
typedef struct syn123_struct syn123_handle;

/** Functions that return an integer success code either return
 *  SYN123_OK if everything went fine, or one of the other detailed
 *  error codes.
 */
enum syn123_error
{
	SYN123_OK = 0   /**< no error */
,	SYN123_BAD_HANDLE /**< bad handle given (NULL) */
,	SYN123_BAD_FMT  /**< bad format (rate/channels) */
,	SYN123_BAD_ENC  /**< bad encoding given */
,	SYN123_BAD_CONV /**< unsupported conversion */
,	SYN123_BAD_SIZE /**< buffer size bad (too small) */
,	SYN123_BAD_BUF  /**< bad buffer pointer (NULL) */
,	SYN123_BAD_CHOP /**< byte buffer not cut at sample boundaries */
,	SYN123_DOOM     /**< Disaster, Out Of Memory. */
,	SYN123_WEIRD    /**< An internal error that should never occur. */
,	SYN123_BAD_FREQ /**< Invalid wave frequency given. */
,	SYN123_BAD_SWEEP /**< Invalid sweep curve given. */
,	SYN123_OVERFLOW  /**< Some fatal (integer) overflow that prevents proper operation. */
,	SYN123_NO_DATA   /**< Not enough data to do something. */
,	SYN123_BAD_DATA  /**< Data present, but not usable. */
};

/** Give a short phrase explaining an error code.
 *  \param errcode the code returned by an API function
 *  \return error phrase
 */
const char* syn123_strerror(int errcode);

/** Create new handle with specified output format.
 *  \param rate sampling rate
 *  \param channels channel count (duplicated mono channels)
 *  \param encoding sample encoding (see enum mpg123_enc_enum)
 *  \param maxbuf maximum buffer size in bytes to allow for caching periodic
 *    signals. When this is given, it is attempted to fit any configured signal
 *    into this buffer in the hope it will work being played periodically. Maybe
 *    there will be tricks with non-periodic signals.
 *    A buffer size of zero turns off buffering and signals are always generated
 *    during extraction.
 *  \param err address to store error code, non-NULL if you care
 *  \return Pointer to allocated handle structure or NULL on error.
 */
MPG123_EXPORT
syn123_handle* syn123_new( long rate, int channels, int encoding
,	size_t maxbuf, int *err );

/** Delete a handle.
 *  \param sh the handle to delete
 */
MPG123_EXPORT
void syn123_del(syn123_handle *sh);

/** Enable/disable dithering for conversions.
 *
 *  The default is no dither for conversions to integer encodings. You can
 *  enable dithering after handle creation, or disable it, as you wish.
 *  Enabling the dither resets the random number generator seed to the provided
 *  value or an internal default if the provided value is zero. The dither noise
 *  is unfiltered with triangular distribution (TPDF), as a sensible common
 *  choice. No further filtering of questionable benefit.
 *
 *  \param sh handle to work on
 *  \param dither Disable (0) or enable (1, actually nonzero) dithering.
 *    Positive values > 1 may trigger differing dither modes in future, but not
 *    now.
 *  \param seed optional address to read the initial seed value from (if non-zero)
 *    and to write the current value to
 */
MPG123_EXPORT
int syn123_dither(syn123_handle *sh, int dither, unsigned long *seed);

/** Extract desired amount of data from the generator.
 *  \param sh handle
 *  \param dst destination buffer
 *  \param dst_bytes number of bytes to extract
 *  \return actual number of extracted bytes
 *  (might differ if dst_bytes is no multiple of the PCM frame size)
 */
MPG123_EXPORT
size_t syn123_read(syn123_handle *sh, void *dst, size_t dst_bytes);

/** Wave types */
enum syn123_wave_id
{
	SYN123_WAVE_INVALID = -1 /**< invalid wave pattern */
,	SYN123_WAVE_FLAT = 0 /**< flat line, silence */
,	SYN123_WAVE_SINE     /**< sinusodial wave*/
,	SYN123_WAVE_SQUARE   /**< square wave */
,	SYN123_WAVE_TRIANGLE /**< triangle wave */
,	SYN123_WAVE_SAWTOOTH /**< sawtooth wave */
,	SYN123_WAVE_GAUSS    /**< Gaussian bell shape */
,	SYN123_WAVE_PULSE    /**< pulse shape, x^2 exp(-A x^2)/S */
,	SYN123_WAVE_SHOT     /**< shot (sharper pulse), x^2 exp(-A x)/S
	* (different values for A and S) */
,	SYN123_WAVE_LIMIT /**< valid IDs below that. A newer release of
	* the library might support more. */
};

/** Setup periodic wave generator.
 *  This sets up a series of oscillators with differing wave shapes.
 *  They are multiplied/scaled with each other instead of mixed
 *  in a sum of signals. It's more fun this way to generate interesting
 *  sounds. If you want to mix differing streams with differing volumes,
 *  channel balance and phase shifts, just create multiple single-channel
 *  generators with a convenient format (float encoding comes to mind)
 *  and mix to your heart's desire. You can then still use this library
 *  to get your channel buffers interleaved and converted to something
 *  your output device likes.
 *
 *  You can ensure strict periodicity without possible shifts in phases
 *  due to floating point rounding errors with the buffered variant.
 *  That may adjust your chosen frequencies to be able to keep the limit
 *  on buffer size, but the resulting data is then strictly periodic
 *  without any further computations that may introduce timing errors.
 *  Apart from possibly saving computing time via the precomputed table,
 *  this is the reason to pre-mix multiple waves into a common buffer at
 *  all.
 *
 *  The adjustments of the wave frequencies also include limiting them
 *  between some minimal value and the Nyquist frequency. Without the
 *  buffer, you can happily choose waves that are not resolved at all
 *  by the sampling rate and get the nasty results. Things get nasty
 *  inside the buffer, too, when you approach the Nyquist limit, but that
 *  is life (and mathematics).
 *
 *  The default wave is a 440 Hz sine without phase offset. If any setting
 *  is missing, the default is taken from that.
 *
 *  \param sh handle
 *  \param count number of waves (if zero, one default wave is configured)
 *  \param id array of wave IDs (enum syn123_wave_id), may be NULL
 *  \param freq array of wave frequencies, may be NULL
 *     Your provided frequencies are overwritten with the actual
 *     values if the periodic buffer is chosen.
 *  \param phase array of wave phases, may be NULL
 *  \param backwards array of true (non-zero) or false (zero), indicating
      whether the wave is being inverted in time, may be NULL
 *  \param period address to store the size of the period buffer
 *    in samples (zero if not using the buffer), ignored if NULL
 *  \return success code
 */
MPG123_EXPORT
int syn123_setup_waves( syn123_handle* sh, size_t count
,	int *id, double *freq, double *phase, int* backwards
,	size_t *period );

/** Query current wave generator setup and state.
 *
 *  This lets you extract the setup of the wave generator and
 *  the current phases to be able to re-create it and continue
 *  seamlessly, or maybe intentionally in some tweaked form with
 *  slight phase or frequency shifts.
 *
 *  You only need to set target pointers to non-NULL where you
 *  want the corresponding values extracted. You need to have
 *  the storage prepared for the correct wave count. A common
 *  mode of usage might be to make a first call to only query
 *  the count, allocate storage, then a do a second call for the
 *  wave data.
 *
 *  \param sh handle
 *  \param count address to store number of waves
 *  \param id storage for array of wave IDs
 *  \param freq storage for array of wave frequencies
 *  \param phase storage for array of wave phases
 *  \param backwards storage for array of true (non-zero) or false (zero),
 *    indicating whether the wave is being inverted in time
 *  \param period address to store the size of the period buffer
 *    in samples (zero if not using the buffer)
 *  \return success code
 */
MPG123_EXPORT
int syn123_query_waves( syn123_handle* sh, size_t *count
,	int *id, double *freq, double *phase, int* backwards
,	size_t *period );

/** Return the name of the indicated wave pattern.
 *  \param id The numerical ID of the wave pattern
 *    (out of enum syn123_wave_id).
 *  \return The name string, guaranteed to be non-NULL.
 *    Invalid codes yield the string "???".
 */
MPG123_EXPORT
const char* syn123_wave_name(int id);

/** Return the wave pattern id given a name string.
 *  \param name The name string.
 *  \return The numerical id (out of enum syn123_wave_id).
 *  
 */
MPG123_EXPORT
int syn123_wave_id(const char *name);

/** Types of frequency sweeps.
 *  There are no functions mapping those to/from strings,
 *  as this list is supposed to be fixed and small.
 *  There are only so many types of sweeps that make sense.
 */
enum syn123_sweep_id
{
	SYN123_SWEEP_LIN = 0 /**< linear frequency change */
,	SYN123_SWEEP_QUAD    /**< quadratic frequency change */
,	SYN123_SWEEP_EXP     /**< exponential (octave per time unit) */
,	SYN123_SWEEP_LIMIT   /**< valid IDs less than that */
};

/** Frequency sweep generator.
 *  This generates a sweep from one frequency to another with one
 *  of the available wave shapes over a given time.
 *  While you can just extract your single sweep by exactly reading
 *  the requestet duration, the generator is set up to run the sweep
 *  a bit longer until the beginning phase is reached again
 *  (one sample before that, of course). That way, reasonably smooth 
 *  periodic playback from the buffer is possible without that nasty jump.
 *  Still, a large freqency difference will result in an audible pop,
 *  but at least that is not whole spectrum due to a phase jump.
 *  \param sh handle
 *  \param wave_id wave ID (enum syn123_wave_id)
 *  \param f1 pointer to beginning frequency in Hz (>= 1e-4, please,
 *    a value <= 0 being replaced by the standard frequency and stored for you)
 *  \param f2 ending frequency in Hz (>= 1e-4, please,
 *    in case of exponential sweep: f2-f1 >= 1e-4, too,
 *    a value <= 0 being replaced by the standard frequency and stored for you)
 *  \param sweep_id choice of sweep curve (enum syn123_sweep_id)
 *  \param smooth enable the periodic smoothing, if sensible
 *    (extending signal beyond the sweep to avoid phase jumps, a continuing
 *    starting at given start phase, not the returned endphase)
 *  \param duration duration of sweep in samples (> 1, please)
 *    This theoretically should be an off_t relating to the size of a
 *    file you produce, but off_t in API headers causes headaches.
 *    A 32 bit size_t still gives you over 24 hours of sweeping with 44100 kHz
 *    rate. On top of that, you can use the returned endphase to chop
 *    your monster sweep into pieces.
 *  \param phase initial phase of the sweep
 *  \param backwards invert the waveform in time if true
 *  \param endphase address to store the normal phase that would
 *    smoothly continue the signal without the period correction
 *    (You can create a following sweep that continues smoothly to a new
 *    target frequency by handing in this endphase as initial phase. Combine
 *    that with phases of constant tone and you could simulate a Theremin
 *    player by approximating the reaction to hand movements via sweeps.)
 *  \param period address to store the periodic sample count, usually
 *    being a bit bigger than the duration for getting the phase
 *    back down; does not imply use of the internal period buffer
 *  \param buffer_period address to store the period count only if the
 *    period buffer is actually in use
 */
MPG123_EXPORT
int syn123_setup_sweep( syn123_handle* sh
,	int wave_id, double phase, int backwards
,	int sweep_id, double *f1, double *f2, int smooth, size_t duration
,	double *endphase, size_t *period, size_t *buffer_period );

/** Set up pink noise generator.
 *  This employs the Gardner/McCartney method to the approximate
 *  the real thing. The number of rows in the tree pattern is tunable.
 *  The result is pink noise with around 2.5 dB/octave. Do not expect
 *  more than 32 bits of randomness (or anything;-).
 *  \param sh handle
 *  \param rows rows for the generator algorithm
 *    It maxes out at 30 rows. Below 1 chooses a default.
 *  \param seed a 32 bit seed value for the pseudo-random number generator
 *  \param period optional address to store the size of the enforced period
 *    (zero for endlessly freshly generated signal)
 *  \return success code
 */
MPG123_EXPORT
int syn123_setup_pink( syn123_handle *sh, int rows, unsigned long seed
,	size_t *period );

/** Set up white noise generator.
 *  A simple white noise source using some cheap pseudo RNG. Do not
 *  expect more than 32 bits of randomness (or anything;-).
 *  \param sh handle
 *  \param seed a 32 bit seed value for the pseudo-random number generator
 *  \param period optional address to store the size of the
 *    enforced period (zero for endlessly freshly generated signal)
 */
MPG123_EXPORT
int syn123_setup_white(syn123_handle *sh, unsigned long seed, size_t *period);

/** Set up Geiger counter simulator.
 *  This models a speaker that is triggered by the pulses from
 *  the Geiger-Mueller counter. That creepy ticking sound.
 *  \param sh handle
 *  \param activity average events per second
 *  \param seed a 32 bit seed value for the pseudo-random number generator
 *  \param period optional address to store the size of the enforced period
 *    (zero for endlessly freshly generated signal)
 *  \return success code
 */
MPG123_EXPORT
int syn123_setup_geiger( syn123_handle *sh, double activity
,	unsigned long seed, size_t *period );

/** Set up silence.
 *  This goes back to the vanilla state.
 *  \return success code
 */
MPG123_EXPORT
int syn123_setup_silence(syn123_handle *sh);

/** Convert between supported encodings.
 *  The buffers must not overlap.
 *  Note that syn123 converts -1.0 to -127 for 8 bit signed (and analogous for
 *  other encodings), but still clips asymmetrically at -128 as that is the
 *  range of the type. If you do explicit clipping using syn123_clip(),
 *  the resulting range will be symmetrical inside [-1.0:1.0] and hence
 *  only down to -127 for 8 bit signed encoding.
 *  The conversions only work directly either from anything to double/float or
 *  from double/float to anything. This process is wrapped in the routine if
 *  a handle is provided, using the fixed mixing buffer in that. Clipping
 *  is only handled for floating point to integer conversions. Also, NaN
 *  is set to zero on conversion and counted as clipped.
 *  The ulaw and alaw conversions use Sun's reference implementation of the
 *  G711 standard (differing from libmpg123's big lookup table).
 *
 *  \param dst destination buffer
 *  \param dst_enc destination encoding (enum mpg123_enc_enum)
 *  \param dst_size size of destination buffer in bytes
 *  \param src source buffer
 *  \param src_enc source encoding
 *  \param src_bytes source buffer size in bytes
 *  \param dst_bytes optional address to store the written byte count to
 *  \param clipped optional address to store number of clipped samples to
 *  \param sh an optional syn123_handle which enables arbitrary encoding
 *    conversions by utilizing the contained buffer as intermediate storage,
 *    can be NULL, disabling any conversion not involving floating point
 *    input or output
 *  \return success code
 */
MPG123_EXPORT
int syn123_conv( void * MPG123_RESTRICT dst, int dst_enc, size_t dst_size
,	void * MPG123_RESTRICT src, int src_enc, size_t src_bytes
,	size_t *dst_bytes, size_t *clipped, syn123_handle * sh );

/** The range of decibel values handled by syn123 goes from
 *  -SYN123_DB_LIMIT to +SYN123_DB_LIMIT
 *  This value ensures that a resulting linear volume can still
 *  be expressed using single-precision float.
 *  The resulting amplitude from -500 dB is still small enough
 *  to drive a 32 bit integer sample value orders of magnitude below
 *  1, so it is effectively a zero. Note that physical volume controls
 *  typically give a range as small as 60 dB. You might want to present
 *  a tighter range to the user than +/- 500 dB!
 */
#define SYN123_DB_LIMIT 500

/** Convert decibels to linear volume (amplitude factor).
 *  This just returns pow(10, db/20) in the supported range.
 *  The dB value is limited according to SYN123_DB_LIMIT, with
 *  NaN being put at the lower end of the range. Better silent
 *  than insanely loud.
 *  \param db relative volume in dB
 *  \return linear volume factor
 */
MPG123_EXPORT double syn123_db2lin(double db);

/** Convert linear volume (amplitude factor) to decibels.
 *  This just returns 20*log10(volume) in the supported range.
 *  The returned value is limited according to SYN123_DB_LIMIT, with
 *  NaN being put at the lower end of the range. Better silent
 *  than insanely loud.
 *  \param volume linear volume factor
 *  \return relative volume in dB
 */
MPG123_EXPORT double syn123_lin2db(double volume);

/** Amplify given buffer.
 *  This multiplies all samples by the given floating point value
 *  (possibly converting to/from floating point on the fly, if a
 *  handle with the included working buffer is given).
 *  Also an offset correction is provided.
 *
 *  \param buf the buffer to work on
 *  \param encoding the sample encoding
 *  \param samples number of samples
 *  \param volume linear volume factor (use syn123_db2lin() for
 *    applying a change in dB)
 *  \param offset offset to add to the sample values before
 *    multiplication
 *  \param clipped optional address to store number of clipped samples to
 *  \param sh optional handle to enable work on non-float
 *    encodings
 *  \return success code (e.g. bad encoding without handle)
 */
MPG123_EXPORT
int syn123_amp( void* buf, int encoding, size_t samples
,	double volume, double offset,	size_t *clipped, syn123_handle *sh );

/** Clip samples in buffer to default range.
 *  This only does anything with floating point encoding, but you can always
 *  call it without damage as a no-op on other encodings. After this, the
 *  samples are guaranteed to be in the range [-1,+1]. NaNs are mapped
 *  to zero (and counted as clipped), so they will still sound bad.
 *  If you want to hard clip to a smaller range, use syn123_soft_clip() with
 *  a width of zero.
 *  \param buf buffer to work on
 *  \param encoding sample encoding
 *  \param samples total number of samples
 *  \return number of clipped samples
 */
MPG123_EXPORT
size_t syn123_clip(void *buf, int encoding, size_t samples);

/** Soft clipping / cheap limiting.
 *  This limits the samples above the threshold of limit-width with a
 *  smooth curve, dampening the high-frequency content of the clipping.
 *  This is no proper frequency filter, but just an independent function on
 *  each sample value, also ignorant of the channel count. This can
 *  directly work on float encodings and does nothing on others unless
 *  a handle is provided for on-line conversion.
 *  \param buf buffer to work on
 *  \param encoding sample encoding
 *  \param samples total number of samples
 *  \param limit the limit to clip to (normally 1 for full scale)
 *  \param width smoothing range
 *  \param sh optional handle to work on non-float encodings
 *  \return number of clipped samples
 */
MPG123_EXPORT
size_t syn123_soft_clip( void *buf, int encoding, size_t samples
,	double limit, double width, syn123_handle *sh );

/** Interleave given number of channels into one stream.
 *  A rather trivial functionality, here for completeness. As the
 *  algorithm is agnostic to what is actually stored as a "sample",
 *  the parameter types are so generic that you could use these
 *  functions to arrange huge structs (hence samplesize as size_t)
 *  or whatever. If that makes sense is up to you.
 *  The buffers shall not overlap!
 *  \param dst destination buffer
 *  \param src source buffer array (one per channel)
 *  \param channels channel count
 *  \param samplesize size of one sample
 *  \param samplecount count of samples per channel
 */
MPG123_EXPORT
void syn123_interleave( void * MPG123_RESTRICT dst, void** MPG123_RESTRICT src
,	int channels, size_t samplesize, size_t samplecount );

/** Deinterleave given number of channels out of one stream.
 *  A rather trivial functionality, here for completeness. As the
 *  algorithm is agnostic to what is actually stored as a "sample",
 *  the parameter types are so generic that you could use these
 *  functions to arrange huge structs (hence samplesize as size_t)
 *  or whatever. If that makes sense is up to you.
 *  The buffers must not overlap!
 *  \param dst destination buffer array (one per channel)
 *  \param src source buffer
 *  \param channels channel count
 *  \param samplesize size of one sample (see MPG123_SAMPLESIZE)
 *  \param samplecount count of samples per channel
 */
MPG123_EXPORT
void syn123_deinterleave( void ** MPG123_RESTRICT dst, void * MPG123_RESTRICT src
,	int channels, size_t samplesize, size_t samplecount );

/** Simply copies mono samples into an interleaved stream.
 *  This might be implemented by a call to syn123_interleave(), it might
 *  be optimized to something different. You could have fun measuring that.
 *  \param dst destination buffer
 *  \param src source buffer
 *  \param channels channel count
 *  \param samplesize size of one sample (see MPG123_SAMPLESIZE)
 *  \param samplecount count of samples per channel
 */
MPG123_EXPORT
void syn123_mono2many( void * MPG123_RESTRICT dst, void * MPG123_RESTRICT src
, int channels, size_t samplesize, size_t samplecount );

/** A little helper/reminder on how interleaved format works:
 *  Produce the offset of the given sample for the given channel.
 */
#define SYN123_IOFF(sample, channel, channels) ((sample)*(channels)+(channel))

/** Specify floating point encoding to use for preserving precision in
 *  intermediate computations for given source and destination encoding.
 *  This should return either MPG123_ENC_FLOAT_32 or MPG123_ENC_FLOAT_64,
 *  unless an uncertain future adds things like 16 bit fp ...
 *  This is what syn123_conv() and syn123_mix() will use internally if
 *  intermediate conversion is necessary.
 *  Note that 64 bit floating point material will be mixed in 32 bit if the
 *  destination encoding does not need more precision.
 *  \param src_enc source/input encoding
 *  \param dst_enc destination/output encoding
 *  \return encoding value, zero if none can be chosen (invalid parameters)
 */
MPG123_EXPORT
int syn123_mixenc(int src_enc, int dst_enc);

/** Mix n input channels on top of m output channels.
 *  This takes an interleaved input stream and mixes its channels
 *  into the output stream given a channel matrix (m,n) where
 *  each of the m rows contains the n volume factors (weights)
 *  to apply when summing the samples from the n input channels.
 *  Sample values are added to what is already present unless
 *  initial silence is explicitly requested.
 *  This works directly with identical floating point encodings. It
 *  may have some optimization to work faster with mono or stereo on
 *  either side and slower generic code for arbitrary channel counts.
 *  You can use syn123_conv() to convert from/to input/output encodings
 *  or provide a syn123_handle to do it on the fly.
 *  There are no optimizations for special cases of mixing factors, so
 *  you should always be able to predict the number of floating point
 *  operations being executed.
 *  For fun, you could give the same problem to a BLAS implementation
 *  of your choice and compare the performance;-)
 *  \param dst destination buffer
 *  \param dst_enc output sample encoding, must be MPG123_ENC_FLOAT_32 or
 *   MPG123_ENC_FLOAT_64 unless a syn123_handle is provided
 *  \param dst_channels destination channel count (m)
 *  \param src source buffer
 *  \param src_enc input sample encoding, must be MPG123_ENC_FLOAT_32 or
 *   MPG123_ENC_FLOAT_64 unless a syn123_handle is provided
 *  \param src_channels source channel count (n)
 *  \param mixmatrix mixing factors ((m,n) matrix), same encoding as
 *    the audio data
 *  \param samples count of samples (PCM frames) to work on
 *  \param silence Set to non-zero value to intialize the output
 *    to a silent signal before adding the input.
 *  \param clipped optional address to store number of clipped samples to
 *    (in case of mixing to an integer encoding)
 *  \param sh an optional syn123_handle which enables work on non-float
 *    encodings by utilizing the contained buffer as intermediate storage,
 *    converting to/from float transparently; Note that this may limit
 *    the amount of channels depending on the fixed internal buffer space.
 *    As long as you have up to 128 channels, you should not worry.
 *  \return success code (e.g. bad encoding, channel counts ...)
 */
MPG123_EXPORT
int syn123_mix( void * MPG123_RESTRICT dst, int dst_enc, int dst_channels
,	void * MPG123_RESTRICT src, int src_enc, int src_channels
,	const double * mixmatrix
,	size_t samples, int silence, size_t *clipped, syn123_handle *sh );

/** Set up a generic digital filter.
 *
 *  This takes a filter order N and coefficient set to prepare
 *  the internal state of a digital filter defined by the transfer
 *  function
 *  \f[
 *             
 *    H(z) = \frac
 *            {b_0 + b_1 z^{-1} + ... + b_N z^{-N}}
 *            {1 + a_1 z^{-1} + ... + a_N z^{-N}}
 *  \f]
 *  It is your task to come up with fun values for the coefficients
 *  b_n and a_n to implement various FIR and IIR filters.
 *
 *  Since it is easy to do and useful, this configures not only a single
 *  filter but a chain of them. If you do configure a chain, the choice
 *  of mixenc and channels must match. No conversion between filters.
 *
 *  \param sh mandatory handle
 *  \param append if true, append a filter to the chain, fresh chain otherwise
 *  \param order filter order N (filter length minus one)
 *  \param b nominator coefficients, starting with b_0 (order+1 elements)
 *  \param a denominator coefficients, starting with a_0=1 (order+1 elements).
 *    It is an error to provide a sequence that does not start with 1.
 *    For a non-recursive (FIR) filter, you can set all following
 *    values from a_1 on to zero or choose to provide a NULL pointer.
 *  \param mixenc either MPG123_ENC_FLOAT_32 or MPG123_ENC_FLOAT_64 for
 *    computation in single or double precision, can be zero to refer to
 *    the value demanded by an already present filter
 *  \param channels number of channels in the audio signal, can be zero
 *    to refer to the value demanded by an already present filter
 *  \param init_firstval If non-zero, initialize the filter history with
 *    a constant stream of the first encountered sample instead of zero.
 *  \return success code
 */
MPG123_EXPORT
int syn123_setup_filter( syn123_handle *sh
,	int append, unsigned int order, double *b, double *a
,	int mixenc, int channels, int init_firstval );

/** init_firstval is the effective setting (extreme coefficients may remove the distinction) */

MPG123_EXPORT
int syn123_query_filter( syn123_handle *sh, size_t position
,	size_t *count, unsigned int *order, double *b, double *a
,	int *mixenc, int *channels, int *init_firstval );

/** drop the n last filters */
MPG123_EXPORT
void syn123_drop_filter(syn123_handle *sh, size_t count);

/** Apply a prepared digital filter.
 * 
 *  This applies the filter prepared by syn123_setup_filter to yur
 *  provided buffer in single or double precision float.
 *  Handing in a non-float encoding is an error. You are supposed
 *  to convert and clip before/after applying the filters.
 *
 *  If you got no filters configured, this always succeeds and
 *  does nothing.
 *
 *  \param sh handle
 *  \param buf audio data to work on (channel count matching what
 *     was given to mpg123_setup_filter())
 *  \param encoding audio encoding
 *  \param samples count of samples (PCM frames) in the buffer
 *  \return success code
 */
MPG123_EXPORT
int syn123_filter( syn123_handle *sh
,	void* buf, int encoding, size_t samples );

/** Set up the resampler.
 *
 *  This works independently of the signal generators. You can combine
 *  syn123_setup_resample() with syn123_setup_geiger(), for example.
 *
 *  People can get worked up a lot about differing algorithms for resampling,
 *  while many folks can actually bear the simple drop/repeat method and most
 *  probably do not bother about the distortions from linear resampling.
 *  A testament to this is that in the 18 years of me maintaining mpg123, I
 *  got bugged about the missing dithering and on subtle bias in shuffling
 *  a playlist, but people seem to insist on using the NtoM resampoler inside
 *  libmpg123, despite me warning about its horrible consequences for audio
 *  quality. It is a plain drop-sample implementation. The only good things to
 *  say about it is that it is cheap and is embedded with the sample-accurate
 *  decoder so that you do not have to worry about offsets in terms of input
 *  and output samples.
 *
 *  Anyhow, this is my take on a reasonably good and efficient resampler that is
 *  neither the best-sounding, nor the fastest in terms of CPU time, but gets
 *  by without significant latency. It needs far less computation than usual
 *  high-quality windowed-sinc resampling (libsamplerate), but cannot beat
 *  libsoxr with its FFT-based approach. The less stringent dirty mode (using
 *  only a 72 dB lowpass filter, in practice still close to CD-DA quality)
 *  comes quite close, though.
 *
 *  The selling point is that it produces output samples as soon as you start
 *  feeding, without any buffering of future samples to fill a window for the
 *  FIR filter or the Fourier transform. It employs IIR filters for low-passing,
 *  possibly in multiple stages for decimation, and optimized interpolation
 *  formulas using up to 6 points. These formulas, based on research by
 *  Olli Niemitalo using using Differential Evolution, are what enables a
 *  dynamic range of 108 dB, well above 16 bit CD-DA quality. Simple
 *  cubic splines after low-passing distort up to around -40 dB in my tests.
 *
 *  There is some effective signal delay well below 10 samples. The impulse
 *  response is about 3 samples late, so this is well inside the realm of
 *  (nonlinear) phase shift. The phase distortion looks bad on paper but does
 *  not matter much in the intended domain of application: the final change in
 *  sampling rate before playback on audio hardware, the last filter that is
 *  applied before the sound hits the speakers (or all the other filters
 *  implemented in your audio harware, that you can choose to be ignorant
 *  about). Use better resamplers for mixing in the studio. Use better
 *  resamplers for converting files on disk. For live playback, consider this
 *  one because it is good enough, fast enough, cheap enough.
 *
 *  Note that if you call this function repeatedly, the internal history
 *  is only cleared if you change anything besides the sampling rates. If
 *  only the rates change, the state of the resampler is kept to enable
 *  you to continue on prior data. This means you can vary the resampling
 *  ratio during operation, somewhat smoothly depending on your buffer size.
 *
 *  Also note that even on identical input and output rates, the resampler
 *  will apply the full filtering as for any ratio close to that, including
 *  oversampling. No special shortcuts.
 *
 *  A returned error guarantees that the internal resampler setup has
 *  been cleared (frees a little bit of memory). You can provide zero
 *  inrate, outrate, and channel count at the same time to that effect
 *  without an error message being produced (but still SYN123_BAD_FMT
 *  being returned).
 *
 *  \param sh mandatory handle
 *  \param inrate input sample rate (nominator of ratio)
 *  \param outrate output sample rate (denominator of ratio)
 *  \param channels number of interleaved channels
 *  \param dirty Enable (!= 0) the dirty mode for even more 'good enough'
 *  resampling with less computing time. Offers -72 dB low pass attentuation,
 *  worst-case distortion around that, too, and 85% worst-case bandwidth.
 *  With this set to zero, the normal mode is used, offering at least 108 dB
 *  dynamic range and worst-case bandwidth above 84%.
 *  \param smooth Enable (!=0) extra code smoothing the resampler response to
 *  on-the-fly changes of sampling rate ratio. This involves keeping
 *  some per-stage history to bootstrap additional decimation filters and the
 *  changed final lowpass/interpolation.
 *  \return success code
 */
MPG123_EXPORT
int syn123_setup_resample( syn123_handle *sh, long inrate, long outrate
,	int channels, int dirty, int smooth );

/** Return the maximum allowed value for sample rates given to the resampler.
 *
 *  Not every possible value of the underlying data type is a valid sample
 *  rate for the resampler. It needs some headroom for computations. This
 *  function returns the maximal rate you can specify. For 32-bit long, this
 *  will be above 1e9, for 64-bit long above 4e18. The minimum is 1, of course.
 *  So, with 32 bit, you can speed up/down by a factor of one million if you
 *  want to keep 0.1% precision in the rates.
 *
 *  \return upper sample rate limit
 */
MPG123_EXPORT
long syn123_resample_maxrate(void);

/** Give upper limit for output sample count from the resampler.
 *
 *  Since there is some rounding involved, the exact number of output samples
 *  from the resampler, being given a certain amount of input samples, can
 *  vary (one more or less than expected). This function is here to give you
 *  a safe output buffer size given a certain input buffer size. If you intend
 *  to vary the output rate for a fixed input rate, you may compute the output
 *  buffer size for the largest intended output rate and use that throughout.
 *  The same applies to the input sample count.
 *  A return value of zero indicates an error (zero, negative, or too large
 *  rate given) unless the given input sample count is also zero.
 *  The resampler only produces output when given new input.
 *  \param inrate input sample rate
 *  \param outrate output sample rate
 *  \param ins input sample count for one buffer
 *  \return number of maximum output samples for one buffer, or zero
 *     if no sensible value exists
 */
MPG123_EXPORT
size_t syn123_resample_count(long inrate, long outrate, size_t ins);

/** Return the amount of input samples needed to recreate resample filter state.
 *
 *  This returns a number of input samples that should fill the internal
 *  filter states good enough for output being close to that produced from
 *  the full input since the beginning. Since recursive filters are employed,
 *  there is no exact number that recreates a state apart from the full sequence
 *  that created it the first time.  This number here shall be more than the
 *  non-recursive history, but not by a huge factor. For extreme cases, this
 *  value may be saturated at SIZE_MAX and thus smaller than what is demanded
 *  by the above definition. It is assumed that you define a maximal practical
 *  size of history to consider for your application, anyway.
 *
 *  \param inrate input sample rate
 *  \param outrate output sample rate
 *  \param dirty switch for dirty resampling mode (see syn123_setup_resample())
 *  \return number of input samples to fill history, zero on error
 */
MPG123_EXPORT
size_t syn123_resample_history(long inrate, long outrate, int dirty);

/** Compute the minimal input sample count needed for given output sample count.
 *
 *  The reverse of syn123_resample_count(), in a way. This gives you the
 *  minimum amount of input samples to guarantee at least the desired amount
 *  of output samples. Once you got that, ensure to call syn123_resample_count()
 *  to get a safe buffer size for that amount of input and prepare accordingly.
 *  With this approach, you can ensure that you get your realtime output device
 *  buffer filled with each loop run fetching a bit of input, at the expense
 *  of handling some additional buffering for the returned sample counts above
 *  the minimum.
 *
 *  \param input_rate input sample rate
 *  \param output_rate output sample rate
 *  \param outs desired minimal output sample count for one input buffer
 *  \return number of minimal input samples in one buffer, or zero if no
 *    sensible value exists (invalid input parameters, or zero outs)
 */
MPG123_EXPORT
size_t syn123_resample_incount(long input_rate, long output_rate, size_t outs);

/** Compute the input sample count needed for close to given output sample
 *  count, but never more.
 *
 *  Call this to get a safe fixed count of input samples to make best use
 *  of a preallocated output buffer, filling it as much as safely possible.
 *  This can also be achieved by calling syn123_resample_incount() and reducing
 *  the value until syn123_resample_count() fits into your buffer.
 *
 *  \param input_rate input sample rate
 *  \param output_rate output sample rate
 *  \param outs desired output sample count for one input buffer
 *  \return number of input samples in one buffer, or zero if no
 *    sensible value exists (invalid input parameters, or zero outs)
 */
MPG123_EXPORT
size_t syn123_resample_fillcount(long input_rate, long output_rate, size_t outs);


/** Compute the maximum number of input samples for a resampling buffer.
 *
 *  Upsampling means that you will get more output samples than input samples.
 *  This larger number still needs to fit into the data type for sample
 *  counts. So please don't feed more than the value returned here.
 *
 *  \param input_rate input sample rate
 *  \param output_rate output sample rate
 *  \return maximum safe input sample count
 */
MPG123_EXPORT
size_t syn123_resample_maxincount(long input_rate, long output_rate);

/** Give exact output sample count for feeding given input now.
 *
 *  This gives you the exact number of samples you will get returned
 *  when feeding the resampler the given additional input samples now,
 *  given the current resampler state contained in the handle.
 *
 *  \param sh syn123 handle
 *  \param ins input sample count
 *  \return output sample count (>=0) or error code
 */
MPG123_EXPORT
ssize_t syn123_resample_expect(syn123_handle *sh, size_t ins);

/** Give minimum input sample count needed now for given output.
 *
 *  This give you the minimal number of input samples needed right
 *  now to yield at least the specified amount of output samples.
 *  Since one input sample can result in several output sampels in one
 *  go, you have to check using syn123_resample_expect() how many
 *  output samples to really expect.
 *
 *  \param sh syn123 handle
 *  \param outs output sample count
 *  \return minimal input sample count (>= 0) or error code
 */
MPG123_EXPORT
ssize_t syn123_resample_inexpect(syn123_handle *sh, size_t outs);


#ifndef SYN123_NO_LARGEFUNC

/* The whole block of off_t-using API is optional to be able to build
   the underlying code without confusing the compiler with prototype
   mismatch. */

/* Lightweight large file hackery to enable worry-reduced use of off_t.
   Depending on the size of off_t in your client build, the corresponding
   library function needs to be chosen. */
#ifndef MPG123_NO_CONFIGURE
#if !defined(MPG123_NO_LARGENAME) && 1
#define MPG123_NO_LARGENAME
#endif
#endif

#if defined(_FILE_OFFSET_BITS) && !defined(MPG123_NO_LARGENAME)
#  if _FILE_OFFSET_BITS+0 == 32
#    define syn123_resample_total   syn123_resample_total_32
#    define syn123_resample_intotal syn123_resample_intotal_32
#  elif _FILE_OFFSET_BITS+0 == 64
#    define syn123_resample_total   syn123_resample_total_64
#    define syn123_resample_intotal syn123_resample_intotal_64
#  else
#    error "Unpredicted _FILE_OFFSET_BITS value."
#  endif
#endif

/** Give exact output sample count for total input sample count.
 *
 *  Use this to determine the total length of your output stream
 *  given the length of the input stream. The computation is exact.
 *  But: It is only valid for a constant sampling rate ratio. If you
 *  play with that during runtime, you need to figure out your output
 *  offset yourself.
 *
 *  \param inrate input sample rate
 *  \param outrate output sample rate
 *  \param ins input sample count for the whole stream
 *  \return number of output samples or -1 if the computation fails
 *    (bad/too large sampling rates, integer overflow)
 */
MPG123_EXPORT
off_t syn123_resample_total(long inrate, long outrate, off_t ins);

/** Give minimum input sample count for total output sample count.
 *
 *  You need to feed at least that amount of input samples to get
 *  the desired amount of output samples from the resampler. Depending
 *  on the resampling ratio, you may in fact get more than the desired
 *  amount (one input sample being worth multiple output samples during
 *  upsampling) so make sure to call syn123_resample_total() to get
 *  the exact number of samples you need to prepare for.
 *  Again, the output offset is only meaninful for a constant sampling
 *  rate ratio.
 *
 *  \param inrate input sample rate
 *  \param outrate output sample rate
 *  \param outs output sample count for the whole stream
 *  \return number of input samples or -1 if the computation fails
 *    (bad/too large sampling rates, integer overflow)
 */
MPG123_EXPORT
off_t syn123_resample_intotal(long inrate, long outrate, off_t outs);

#endif

/** Resample input buffer to output buffer.
 *
 *  This executes the resampling configured by syn123_setup_resample(). The
 *  input and output encoding is fixed at single-precision float
 *  (MPG123_ENC_FLOAT_32) and multiple channels are interleaved. There
 *  is no implicit conversion of other encodings since the fixed internal
 *  buffers for that may not fit your chosen extreme resampling ratios. Also,
 *  dealing with double precision does not make sense with the mathematical
 *  limitations of the employed filters.
 *
 *  You are responsible for having your buffers prepared with the correct sizes.
 *  Use syn123_resample_count() to ensure that you are prepared for the correct
 *  number of output samples given your input sample count.
 *
 *  Also, ensuring a minimal number of input samples using
 *  syn123_resample_incount() helps to identify an error situation where zero
 *  samples are returned (which would be valid for low input sample count).
 *  The only error apart from handing in an invalid/unprepared handle is
 *  a too large number of input samples. Check syn123_resample_maxincount().
 *
 *  \param sh handle with prepared resampling method
 *    If this is NULL or if the resampler has not been initialized before, the
 *    function returns zero instead of crashing randomly.
 *  \param dst destination buffer
 *  \param src source buffer
 *  \param samples input samples (PCM frames) in source buffer
 *  \return number of output samples (PCM frames)
 */
MPG123_EXPORT
size_t syn123_resample( syn123_handle *sh,
	float * MPG123_RESTRICT dst, float * MPG123_RESTRICT src, size_t samples );

/** Swap byte order between little/big endian.
 *  \param buf buffer to work on
 *  \param samplesize size of one sample (see MPG123_SAMPLESIZE)
 *  \param samplecount count of samples
 */
MPG123_EXPORT
void syn123_swap_bytes(void* buf, size_t samplesize, size_t samplecount);

/* Wrappers over the above to convert to/from syn123's native byte order
   from/to little or big endian. */

/** Convert from host order to little endian.
 *  \param buf buffer to work on
 *  \param samplesize size of one sample (see MPG123_SAMPLESIZE)
 *  \param samplecount count of samples
 */
MPG123_EXPORT
void syn123_host2le(void *buf, size_t samplesize, size_t samplecount);

/** Convert from host order to big endian.
 *  \param buf buffer to work on
 *  \param samplesize size of one sample (see MPG123_SAMPLESIZE)
 *  \param samplecount count of samples
 */
MPG123_EXPORT
void syn123_host2be(void *buf, size_t samplesize, size_t samplecount);

/** Convert from little endian to host order.
 *  \param buf buffer to work on
 *  \param samplesize size of one sample (see MPG123_SAMPLESIZE)
 *  \param samplecount count of samples
 */
MPG123_EXPORT
void syn123_le2host(void *buf, size_t samplesize, size_t samplecount);

/** Convert from big endian to host order.
 *  \param buf buffer to work on
 *  \param samplesize size of one sample (see MPG123_SAMPLESIZE)
 *  \param samplecount count of samples
 */
MPG123_EXPORT
void syn123_be2host(void *buf, size_t samplesize, size_t samplecount);

/** @} */

#ifdef __cplusplus
}
#endif

#endif

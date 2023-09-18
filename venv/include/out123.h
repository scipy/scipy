/*
	out123: audio output interface

	copyright 1995-2016 by the mpg123 project,
	free software under the terms of the LGPL 2.1

	see COPYING and AUTHORS files in distribution or http://mpg123.org
	initially written as audio.h by Michael Hipp, reworked into out123 API
	by Thomas Orgis
*/

#ifndef _OUT123_H_
#define _OUT123_H_

/** \file out123.h The header file for the libout123 audio output facility. */

/* We only need size_t definition. */
#include <stddef.h>

/* Common audio encoding specification, including a macro for getting
 *  size of encoded samples in bytes. Said macro is still hardcoded
 *  into out123_encsize(). Relying on this one may help an old program
 *  know sizes of encodings added to fmt123.h later on.
 *  If you don't care, just use the macro.
 */
#include <fmt123.h>

/** A macro to check at compile time which set of API functions to expect.
 * This should be incremented at least each time a new symbol is added
 * to the header.
 */
#ifndef OUT123_API_VERSION
#define OUT123_API_VERSION 4
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

/* Earlier versions of libout123 put enums into public API calls,
 * thich is not exactly safe. There are ABI rules, but you can use
 * compiler switches to change the sizes of enums. It is safer not
 * to have them in API calls. Thus, the default is to remap calls and
 * structs to variants that use plain ints. Define MPG123_ENUM_API to
 * prevent that remapping.
 *
 * You might want to define this to increase the chance of your binary
 * working with an older version of the library. But if that is your goal,
 * you should better build with an older version to begin with.
 */
#ifndef MPG123_ENUM_API

#define out123_param    out123_param2
#define out123_getparam out123_getparam2

#endif

#ifdef __cplusplus
extern "C" {
#endif

/** \defgroup out123_api out123 library API
 *  This is out123, a library focused on continuous playback of audio streams
 *  via various platform-specific output methods. It glosses over details of
 *  the native APIs to give an interface close to simply writing data to a
 *  file. There might be the option to tune details like buffer (period) sizes
 *  and the number of them on the device side in future, but the focus of the
 *  library is to ease the use case of just getting that raw audio data out
 *  there, without interruptions.
 *
 *  The basic idea is to create a handle with out123_new() and open a certain
 *  output device (using a certain driver module, possibly build-time defaults)
 *  with out123_open(). Now, you can query the output device for supported
 *  encodings for given rate and channel count with out123_get_encodings() and
 *  decide what to use for actually starting playback with out123_start().
 *
 *  Then, you just need to provide (interleaved pcm) data for playback with
 *  out123_play(), which will block when the device's buffers are full. You get
 *  your timing from that (instead of callbacks). If your program does the
 *  production of the audio data just a little bit faster than the playback,
 *  causing out123_play() to block ever so briefly, you're fine.
 *
 *  You stop playback with out123_stop(), or just close the device and driver
 *  via out123_close(), or even just decide to drop it all and do out123_del()
 *  right away when you're done.
 *
 *  There are other functions for specific needs, but the basic idea should be
 *  covered by the above.
 *
 *  Note that the driver modules that bind to the operating system API for
 *  output might impose restrictions on what you can safely do regarding your
 *  out123_handle and multiple threads or processes. You should be on the safe
 *  side ensuring that you confine usage of a handle to a single thread instead
 *  of passing it around.
 @{
 */

/** Opaque structure for the libout123 handle. */
struct out123_struct;
/** Typedef shortcut as preferrend name for the handle type. */
typedef struct out123_struct out123_handle;

/** Enumeration of codes for the parameters that it is possible to set/get. */
enum out123_parms
{
	OUT123_FLAGS = 1 /**< integer, various flags, see enum #out123_flags */
,	OUT123_PRELOAD /**< float, fraction of buffer to fill before playback */
,	OUT123_GAIN    /**< integer, output device gain (module-specific) */
,	OUT123_VERBOSE /**< integer, verbosity to stderr, >= 0 */
,	OUT123_DEVICEBUFFER /**<
 *  float, length of device buffer in seconds;
 *  This might be ignored, might have only a loose relation to actual
 *  buffer sizes and latency, depending on output driver. Try to tune
 *  this before opening a device if you want to influcence latency or reduce
 *  dropouts. Value <= 0 uses some default, usually favouring stable playback
 *  over low latency. Values above 0.5 are probably too much.
 */
,	OUT123_PROPFLAGS /**< integer, query driver/device property flags (r/o) */
,	OUT123_NAME /**< string, name of this instance (NULL restores default);
 * The value returned by out123_getparam() might be different if the audio
 * backend changed it (to be unique among clients, p.ex.).
 * TODO: The name provided here is used as prefix in diagnostic messages. */
,	OUT123_BINDIR /**< string, path to a program binary directory to use
 * as starting point in the search for the output module directory
 * (e.g. ../lib/mpg123 or ./plugins). The environment variable MPG123_MODDIR
 * is always tried first and the in-built installation path last.
 */
,	OUT123_ADD_FLAGS /**< enable given flags */
,	OUT123_REMOVE_FLAGS /**< disable diven flags */
};

/** Flags to tune out123 behaviour */
enum out123_flags
{
	OUT123_HEADPHONES       = 0x01 /**< output to headphones (if supported) */
,	OUT123_INTERNAL_SPEAKER = 0x02 /**< output to speaker (if supported) */
,	OUT123_LINE_OUT         = 0x04 /**< output to line out (if supported) */
,	OUT123_QUIET            = 0x08 /**< no printouts to standard error */
,	OUT123_KEEP_PLAYING     = 0x10 /**<
 *  When this is set (default), playback continues in a loop when the device
 *  does not consume all given data at once. This happens when encountering
 *  signals (like SIGSTOP, SIGCONT) that cause interruption of the underlying
 *  functions.
 *  Note that this flag is meaningless when the optional buffer is employed,
 *  There, your program will always block until the buffer completely took
 *  over the data given to it via out123_play(), unless a communication error
 *  arises.
 */
,	OUT123_MUTE             = 0x20 /**< software mute (play silent audio) */
};

/** Read-only output driver/device property flags (OUT123_PROPFLAGS). */
enum out123_propflags
{
	OUT123_PROP_LIVE = 0x01 /**< This is a live output, meaning that
 *  special care might be needed for pauses in playback (p.ex. stream
 *  of silence instead of interruption), as opposed to files on disk.
 */
,	OUT123_PROP_PERSISTENT = 0x02 /**< This (live) output does not need
 *  special care for pauses (continues with silence itself),
 *  out123_pause() does nothing to the device.
 */
};

/** Create a new output handle.
 *  This only allocates and initializes memory, so the only possible
 *  error condition is running out of memory.
 * \return pointer to new handle or NULL on error
 */
MPG123_EXPORT
out123_handle *out123_new(void);

/** Delete output handle.
 *  This implies out123_close().
 */
MPG123_EXPORT
void out123_del(out123_handle *ao);

/** Free plain memory allocated within libout123.
 *  This is for library users that are not sure to use the same underlying
 *  memory allocator as libout123. It is just a wrapper over free() in
 *  the underlying C library.
 */
MPG123_EXPORT void out123_free(void *ptr);

/** Error code enumeration
 * API calls return a useful (positve) value or zero (OUT123_OK) on simple
 * success. A negative value (-1 == OUT123_ERR) usually indicates that some
 * error occured. Which one, that can be queried using out123_errcode()
 * and friends.
 */
enum out123_error
{
	OUT123_ERR = -1 /**< generic alias for verbosity, always == -1 */
,	OUT123_OK  = 0  /**< just a name for zero, not going to change */
,	OUT123_DOOM /**< dazzled, out of memory */
,	OUT123_BAD_DRIVER_NAME /**< bad driver name given */
,	OUT123_BAD_DRIVER /**< unspecified issue loading a driver */
,	OUT123_NO_DRIVER /**< no driver loaded */
,	OUT123_NOT_LIVE /**< no active audio device */
,	OUT123_DEV_PLAY /**< some device playback error */
,	OUT123_DEV_OPEN /**< error opening device */
,	OUT123_BUFFER_ERROR /**<
 * Some (really unexpected) error in buffer infrastructure.
 */
,	OUT123_MODULE_ERROR /**< basic failure in module loading */
,	OUT123_ARG_ERROR /**< some bad function arguments supplied */
,	OUT123_BAD_PARAM /**< unknown parameter code */
,	OUT123_SET_RO_PARAM /**< attempt to set read-only parameter */
,	OUT123_BAD_HANDLE /**< bad handle pointer (NULL, usually) */
,	OUT123_NOT_SUPPORTED /**< some requested operation is not supported (right now) */
,	OUT123_DEV_ENUMERATE /**< device enumeration itself failed */
,	OUT123_ERRCOUNT /**< placeholder for shaping arrays */
};

/** Get string representation of last encountered error in the
 *  context of given handle.
 * \param ao handle
 * \return error string
 */
MPG123_EXPORT
const char* out123_strerror(out123_handle *ao);

/** Get the plain errcode intead of a string.
 * Note that this used to return OUT123_ERR instead of
 * OUT123_BAD_HANDLE in case of ao==NULL before mpg123-1.23.5 .
 * \param ao handle
 * \return error code recorded in handle or OUT123_BAD_HANDLE
 */
MPG123_EXPORT
int out123_errcode(out123_handle *ao);

/** Return the error string for a given error code.
 * \param errcode the integer error code
 * \return error string
 */
MPG123_EXPORT
const char* out123_plain_strerror(int errcode);

/** Set a desired output buffer size.
 *  This starts a separate process that handles the audio output, decoupling
 *  the latter from the main process with a memory buffer and saving you the
 *  burden to ensure sparing CPU cycles for actual playback.
 *  This is for applicatons that prefer continuous playback over small latency.
 *  In other words: The kind of applications that out123 is designed for.
 *  This routine always kills off any currently active audio output module /
 *  device, even if you just disable the buffer when there is no buffer.
 *
 *  Keep this in mind for memory-constrainted systems: Activating the
 *  buffer causes a fork of the calling process, doubling the virtual memory
 *  use. Depending on your operating system kernel's behaviour regarding
 *  memory overcommit, it might be wise to call out123_set_buffer() very
 *  early in your program before allocating lots of memory.
 *
 *  There _might_ be a change to threads in future, but for now this is
 *  classic fork with shared memory, working without any threading library.
 *  If your platform or build does not support that, you will always get an
 *  error on trying to set up a non-zero buffer (but the API call will be
 *  present).
 *
 *  Also, if you do intend to use this from a multithreaded program, think
 *  twice and make sure that your setup is happy with forking full-blown
 *  processes off threaded programs. Probably you are better off spawning a
 *  buffer thread yourself.
 *
 * \param ao handle
 * \param buffer_bytes size (bytes) of a memory buffer for decoded audio,
 *    a value of zero disables the buffer.
 * \return 0 on success, OUT123_ERR on error
 */
MPG123_EXPORT
int out123_set_buffer(out123_handle *ao, size_t buffer_bytes);

#ifdef MPG123_ENUM_API
/** Set a parameter on a out123_handle.
 *
 *  Note that this name is mapped to out123_param2() instead unless
 *  MPG123_ENUM_API is defined.
 *
 *  The parameters usually only change what happens on next out123_open, not
 *  incfluencing running operation. There are macros To ease the API a bit:
 *  You can call out123_param_int(ao, code, value) for integer (long) values,
 *  same with out123_param_float() and out123_param_string().
 *
 * \param ao handle
 * \param code parameter code
 * \param value input value for integer parameters
 * \param fvalue input value for floating point parameters
 * \param svalue input value for string parameters (contens are copied)
 * \return 0 on success, OUT123_ERR on error.
 */
MPG123_EXPORT
int out123_param( out123_handle *ao, enum out123_parms code
,                 long value, double fvalue, const char *svalue );
#endif

/** Set a parameter on a out123_handle. No enum.
 *
 *  This is actually called instead of out123_param()
 *  unless MPG123_ENUM_API is defined.
 *
 *  The parameters usually only change what happens on next out123_open, not
 *  incfluencing running operation. There are macros To ease the API a bit:
 *  You can call out123_param_int(ao, code, value) for integer (long) values,
 *  same with out123_param_float() and out123_param_string().
 *
 * \param ao handle
 * \param code parameter code (from enum #out123_parms)
 * \param value input value for integer parameters
 * \param fvalue input value for floating point parameters
 * \param svalue input value for string parameters (contens are copied)
 * \return 0 on success, OUT123_ERR on error.
 */
MPG123_EXPORT
int out123_param2( out123_handle *ao, int code
,                 long value, double fvalue, const char *svalue );


/** Shortcut for out123_param() to set an integer parameter. */
#define out123_param_int(ao, code, value) \
	out123_param((ao), (code), (value), 0., NULL)
/** Shortcut for out123_param() to set a float parameter. */
#define out123_param_float(ao, code, value) \
	out123_param((ao), (code), 0, (value), NULL)
/** Shortcut for out123_param() to set an string parameter. */
#define out123_param_string(ao, code, value) \
	out123_param((ao), (code), 0, 0., (value))

#ifdef MPG123_ENUM_API
/** Get a parameter from an out123_handle.
 *
 *  Note that this name is mapped to out123_param2() instead unless
 *  MPG123_ENUM_API is defined.
 *
 * \param ao handle
 * \param code parameter code
 * \param ret_value output address for integer parameters
 * \param ret_fvalue output address for floating point parameters
 * \param ret_svalue output address for string parameters (pointer to
 *        internal memory, so no messing around, please)
 * \return 0 on success, OUT123_ERR on error (bad parameter name or bad handle).
 */
MPG123_EXPORT
int out123_getparam( out123_handle *ao, enum out123_parms code
,                    long *ret_value, double *ret_fvalue, char* *ret_svalue );
#endif

/** Get a parameter from an out123_handle. No enum.
 *
 *  This is actually called instead of out123_getparam()
 *  unless MPG123_ENUM_API is defined.
 *
 * \param ao handle
 * \param code parameter code (from enum #out123_parms)
 * \param ret_value output address for integer parameters
 * \param ret_fvalue output address for floating point parameters
 * \param ret_svalue output address for string parameters (pointer to
 *        internal memory, so no messing around, please)
 * \return 0 on success, OUT123_ERR on error (bad parameter name or bad handle).
 */
MPG123_EXPORT
int out123_getparam2( out123_handle *ao, int code
,                    long *ret_value, double *ret_fvalue, char* *ret_svalue );

/** Shortcut for out123_getparam() to get an integer parameter. */
#define out123_getparam_int(ao, code, value) \
	out123_getparam((ao), (code), (value), NULL, NULL)
/** Shortcut for out123_getparam() to get a float parameter. */
#define out123_getparam_float(ao, code, value) \
	out123_getparam((ao), (code), NULL, (value), NULL)
/** Shortcut for out123_getparam() to get a string parameter. */
#define out123_getparam_string(ao, code, value) \
	out123_getparam((ao), (code), NULL, NULL, (value))

/** Copy parameters from another out123_handle.
 * \param ao handle
 * \param from_ao the handle to copy parameters from
 * \return 0 in success, -1 on error
 */
MPG123_EXPORT
int out123_param_from(out123_handle *ao, out123_handle* from_ao);

/** Get list of driver modules reachable in system in C argv-style format.
 *
 *  The client is responsible for freeing the memory of both the individual
 *  strings and the lists themselves.  There is out123_stringlists_free()
 *  to assist.
 *
 *  A module that is not loadable because of missing libraries is simply
 *  skipped. You will get stderr messages about that unless OUT123_QUIET was
 *  was set, though. Failure to open the module directory is a serious error,
 *  resulting in negative return value.
 *
 * \param ao handle
 * \param names address for storing list of names
 * \param descr address for storing list of descriptions
 * \return number of drivers found, -1 on error
 */
MPG123_EXPORT
int out123_drivers(out123_handle *ao, char ***names, char ***descr);

/** Get a list of available output devices for a given driver.
 *
 *  If the driver supports enumeration, you can get a listing of possible
 *  output devices. If this list is exhaustive, depends on the driver.
 *  Note that this implies out123_close(). When you have a device already
 *  open, you don't need to look for one anymore. If you really do, just
 *  create another handle.
 *
 *  Your provided pointers are only used for non-negative return values.
 *  In this case, you are responsible for freeing the associated memory of
 *  the strings and the lists themselves. The format of the lists is an
 *  array of char pointers, with the returned count just like the usual
 *  C argv and argc. There is out123_stringlists_free() to assist.
 *
 *  Note: Calling this on a handle with a configured buffer process will
 *  yield #OUT123_NOT_SUPPORTED.
 *
 * \param ao handle
 * \param driver driver name or comma-separated list of names
 *   to try, just like for out123_open(), possibly NULL for some default
 * \param names address for storing list of names
 * \param descr address for storing list of descriptions
 * \param active_driver address for storing a copy of the actually active
 *    driver name (in case you gave a list or NULL as driver), can be NULL
 *    if not interesting
 * \return count of devices or #OUT123_ERR if some error was encountered,
 *   possibly just #OUT123_NOT_SUPPORTED if the driver lacks enumeration support
 */
MPG123_EXPORT
int out123_devices( out123_handle *ao, const char *driver
,	char ***names, char ***descr, char **active_driver );

/** Helper to free string list memory.
 *
 *  This aids in freeing the memory allocated by out123_devices() and
 *  out123_drivers().
 *
 *  Any of the given lists can be NULL and nothing will happen to it.
 *
 * \param name first string list
 * \param descr second string list
 * \param count count of strings
 */
MPG123_EXPORT
void out123_stringlists_free(char **name, char **descr, int count);


/** Open an output device with a certain driver
 *  Note: Opening means that the driver code is loaded and the desired
 *  device name recorded, possibly tested for availability or tentatively
 *  opened. After out123_open(), you can ask for supported encodings
 *  and then really open the device for playback with out123_start().
 * \param ao handle
 * \param driver (comma-separated list of) output driver name(s to try),
 *               NULL for default
 * \param device device name to open, NULL for default 
 *               (stdout for file-based drivers)
 * \return 0 on success, -1 on error.
 */
MPG123_EXPORT
int out123_open(out123_handle *ao, const char* driver, const char* device);

/** Give info about currently loaded driver and device
 *  Any of the return addresses can be NULL if you are not interested in
 *  everything. You get pointers to internal storage. They are valid
 *  as long as the driver/device combination is opened.
 *  The device may be NULL indicating some unnamed default.
 *  TODO: Make the driver modules return names for such defaults.
 * \param ao handle
 * \param driver return address for driver name
 * \param device return address for device name
 * \return 0 on success, -1 on error (i.e. no driver loaded)
 */
MPG123_EXPORT
int out123_driver_info(out123_handle *ao, char **driver, char **device);

/** Close the current output device and driver.
 *  This implies out123_drain() to ensure no data is lost.
 *  With a buffer, that might cause considerable delay during
 *  which your main application is blocked waiting.
 *  Call out123_drop() beforehand if you want to end things
 *  quickly.
 * \param ao handle
 */
MPG123_EXPORT
void out123_close(out123_handle *ao);

/** Get supported audio encodings for given rate and channel count,
 *  for the currently openend audio device.
 *  Usually, a wider range of rates is supported, but the number
 *  of sample encodings is limited, as is the number of channels.
 *  So you can call this with some standard rate and hope that the
 *  returned encodings work also for others, with the tested channel
 *  count.
 *  The return value of -1 on some encountered error conveniently also
 *  does not match any defined format (only 15 bits used for encodings,
 *  so this would even work with 16 bit integers).
 *  This implies out123_stop() to enter query mode.
 * \param ao handle
 * \param rate sampling rate
 * \param channels number of channels
 * \return supported encodings combined with bitwise or, to be checked
 *         against your favourite bitmask, -1 on error
 */
MPG123_EXPORT
int out123_encodings(out123_handle *ao, long rate, int channels);

/** Return the size (in bytes) of one mono sample of the named encoding.
 * \param encoding The encoding value to analyze.
 * \return positive size of encoding in bytes, 0 on invalid encoding. */
MPG123_EXPORT int out123_encsize(int encoding);

/** Get list of supported formats for currently opened audio device.
 *  Given a list of sampling rates and minimal/maximal channel count,
 *  this quickly checks what formats are supported with these
 *  constraints. The first entry is always reserved for a default
 *  format for the output device. If there is no such default,
 *  all values of the format are -1.
 *  For each requested combination of rate and channels, a format entry is
 *  created, possible with encoding value 0 to indicate that this combination
 *  has been tested and rejected. So, when there is no basic error, the
 *  number of returned format entries should be
 *     (ratecount*(maxchannels-minchannels+1)+1)
 *  . But instead of forcing you to guess, this will be allocated by
 *  successful run.
 *  For the first entry, the encoding member is supposed to be a definite
 *  encoding, for the others it is a bitwise combination of all possible
 *  encodings.
 *  This function is more efficient than many calls to out123_encodings().
 * \param ao handle
 * \param rates pointer to an array of sampling rates, may be NULL for none
 * \param ratecount number of provided sampling rates
 * \param minchannels minimal channel count
 * \param maxchannels maximal channel count
 * \param fmtlist return address for array of supported formats
 *        the encoding field of each entry is a combination of all
 *        supported encodings at this rate and channel count;
 *        Memory shall be freed by user.
 * \return number of returned format enries, -1 on error
 */
MPG123_EXPORT
int out123_formats( out123_handle *ao, const long *rates, int ratecount
                  , int minchannels, int maxchannels
                  , struct mpg123_fmt **fmtlist );

/** Get list of encodings known to the library.
 *  You are responsible for freeing the allocated array.
 * \param enclist return address for allocated array of encoding codes
 * \return number of encodings, -1 on error
 */
MPG123_EXPORT
int out123_enc_list(int **enclist);

/** Find encoding code by name.
 * \param name short or long name to find encoding code for
 * \return encoding if found (enum #mpg123_enc_enum), else 0
 */
MPG123_EXPORT
int out123_enc_byname(const char *name);

/** Get name of encoding.
 * \param encoding code (enum #mpg123_enc_enum)
 * \return short name for valid encodings, NULL otherwise
 */
MPG123_EXPORT
const char* out123_enc_name(int encoding);

/** Get long name of encoding.
 * \param encoding code (enum #mpg123_enc_enum)
 * \return long name for valid encodings, NULL otherwise
 */
MPG123_EXPORT
const char* out123_enc_longname(int encoding);

/** Start playback with a certain output format
 *  It might be a good idea to have audio data handy to feed after this
 *  returns with success.
 *  Rationale for not taking a pointer to struct mpg123_fmt: This would
 *  always force you to deal with that type and needlessly enlarge the
 *  shortest possible program.
 * \param ao handle
 * \param encoding sample encoding (values matching libmpg123 API)
 * \param channels number of channels (1 or 2, usually)
 * \param rate sampling rate
 * \return 0 on success, negative on error (bad format, usually)
 */
MPG123_EXPORT
int out123_start( out123_handle *ao
,                 long rate, int channels, int encoding );

/** Pause playback
 *  Interrupt playback, holding any data in the optional buffer.
 *
 *  This closes the audio device if it is a live sink, ready to be re-opened
 *  by out123_continue() or out123_play() with the existing parameters.
 * \param ao handle
 */
MPG123_EXPORT
void out123_pause(out123_handle *ao);

/** Continue playback
 *  The counterpart to out123_pause(). Announce to the driver that playback
 *  shall continue.
 *
 *  Playback might not resume immediately if the optional buffer is configured
 *  to wait for a minimum fill and close to being empty. You can force playback
 *  of the last scrap with out123_drain(), or just by feeding more data with
 *  out123_play(), which will trigger out123_continue() for you, too.
 * \param ao handle
 */
MPG123_EXPORT
void out123_continue(out123_handle *ao);

/** Stop playback.
 *  This waits for pending audio data to drain to the speakers.
 *  You might want to call out123_drop() before stopping if you want
 *  to end things right away.
 * \param ao handle
 */
MPG123_EXPORT
void out123_stop(out123_handle *ao);

/** Hand over data for playback and wait in case audio device is busy.
 *  This survives non-fatal signals like SIGSTOP/SIGCONT and keeps on
 *  playing until the buffer is done with if the flag
 *  OUT123_KEEP_PLAYING ist set (default). So, per default, if
 *  you provided a byte count divisible by the PCM frame size, it is an
 *  error when less bytes than given are played.
 *  To be sure if an error occured, check out123_errcode().
 *  Also note that it is no accident that the buffer parameter is not marked
 *  as constant. Some output drivers might need to do things like swap
 *  byte order. This is done in-place instead of wasting memory on yet
 *  another copy. Software muting also overwrites the data.
 * \param ao handle
 * \param buffer pointer to raw audio data to be played
 * \param bytes number of bytes to read from the buffer
 * \return number of bytes played (might be less than given, even zero)
 */
MPG123_EXPORT
size_t out123_play( out123_handle *ao
                  , void *buffer, size_t bytes );

/** Drop any buffered data, making next provided data play right away.
 *  This does not imply an actual pause in playback.
 *  You are expected to play something, unless you called out123_pause().
 *  Feel free to call out123_stop() afterwards instead for a quicker
 *  exit than the implied out123_drain().
 *  For live sinks, this may include dropping data from their buffers.
 *  For others (files), this only concerns data in the optional buffer.
 * \param ao handle
 */
MPG123_EXPORT
void out123_drop(out123_handle *ao);

/** Drain the output, waiting until all data went to the hardware.
 * This does imply out123_continue() before and out123_pause()
 * after draining.
 * This might involve only the optional buffer process, or the
 * buffers on the audio driver side, too.
 * \param ao handle
 */
MPG123_EXPORT
void out123_drain(out123_handle *ao);

/** Drain the output, but only partially up to the given number of
 *  bytes. This gives you the opportunity to do something while
 *  the optional buffer is writing remaining data instead of having
 *  one atomic API call for it all.
 *
 *  It is wholly expected that the return value of out123_buffered()
 *  before and after calling this has a bigger difference than the
 *  provided limit, as the buffer is writing all the time in the
 *  background.
 *
 *  This is just a plain out123_drain() if the optional buffer is not
 *  in use. Also triggers out123_continue(), but only out123_pause()
 *  if there is no buffered data anymore.
 * \param ao handle
 * \param bytes limit of buffered bytes to drain
 * \return number of bytes drained from buffer
 */
MPG123_EXPORT
void out123_ndrain(out123_handle *ao, size_t bytes);

/** Get an indication of how many bytes reside in the optional buffer.
 * This might get extended to tell the number of bytes queued up in the
 * audio backend, too.
 * \param ao handle
 * \return number of bytes in out123 library buffer
 */
MPG123_EXPORT
size_t out123_buffered(out123_handle *ao);

/** Extract currently used audio format from handle.
 *  matching mpg123_getformat().
 *  Given return addresses may be NULL to indicate no interest.
 * \param ao handle
 * \param rate address for sample rate
 * \param channels address for channel count
 * \param encoding address for encoding
 * \param framesize size of a full PCM frame (for convenience)
 * \return 0 on success, -1 on error
 */
MPG123_EXPORT
int out123_getformat( out123_handle *ao
,	long *rate, int *channels, int *encoding, int *framesize );

/** @} */

#ifdef __cplusplus
}
#endif

#endif


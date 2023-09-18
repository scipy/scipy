---
layout: page
title: libsndfile : the sf_command function.
---

# sf_command

```c
int sf_command (SNDFILE *sndfile, int cmd, void *data, int datasize) ;
```

This function allows the caller to retrieve information from or change aspects
of the library behaviour. Examples include retrieving a string containing the
library version or changing the scaling applied to floating point sample data
during read and write. Most of these operations are performed on a per-file
basis.

The cmd parameter is an integer identifier which is defined in *sndfile.h*. All
of the valid command identifiers have names beginning with "SFC_". Data is
passed to and returned from the library by use of a void pointer. The library
will not read or write more than datasize bytes from the void pointer. For some
calls no data is required in which case data should be NULL and datasize may be
used for some other purpose.

The available commands are as follows:

| Name                                                              | Description                                             |
|:------------------------------------------------------------------|:--------------------------------------------------------|
| [SFC_GET_LIB_VERSION](#sfc_get_lib_version)                       | Retrieve the version of the library as a string.        |
| [SFC_GET_LOG_INFO](#sfc_get_log_info)                             | Retrieve the internal per-file operation log.           |
| [SFC_GET_CURRENT_SF_INFO](#sfc_get_current_sf_info)               | Retrieve `SF_INFO` struct of opened file.               |
| [SFC_CALC_SIGNAL_MAX](#sfc_calc_signal_max)                       | Calculate the measured maximum signal value.            |
| [SFC_CALC_NORM_SIGNAL_MAX](#sfc_calc_norm_signal_max)             | Calculate the measured normalised maximum signal value. |
| [SFC_CALC_MAX_ALL_CHANNELS](#sfc_calc_max_all_channels)           | Calculate the peak value for each channel.              |
| [SFC_CALC_NORM_MAX_ALL_CHANNELS](#sfc_calc_norm_max_all_channels) | Calculate the normalised peak for each channel.         |
| [SFC_GET_SIGNAL_MAX](#sfc_get_signal_max)                         | Retrieve the peak value for the file.                   |
| [SFC_GET_MAX_ALL_CHANNELS](#sfc_get_max_all_channels)             | Retrieve the peak value for each channel.               |
| [SFC_SET_NORM_FLOAT](#sfc_set_norm_float)                         | Set float normalisation behaviour.                      |
| [SFC_SET_NORM_DOUBLE](#sfc_set_norm_double)                       | Set double normalisation behaviour.                     |
| [SFC_GET_NORM_FLOAT](#sfc_get_norm_float)                         | Get float normalisation behaviour.                      |
| [SFC_GET_NORM_DOUBLE](#sfc_get_norm_double)                       | Get double normalisation behaviour.                     |
| [SFC_SET_SCALE_FLOAT_INT_READ](#sfc_set_scale_float_int_read)     | Control scale factor on read.                           |
| [SFC_SET_SCALE_INT_FLOAT_WRITE](#sfc_set_scale_int_float_write)   | Control scale factor on write.                          |
| [SFC_GET_SIMPLE_FORMAT_COUNT](#sfc_get_simple_format_count)       | Get simple formats count.                               |
| [SFC_GET_SIMPLE_FORMAT](#sfc_get_simple_format)                   | Get information about a simple format.                  |
| [SFC_GET_FORMAT_INFO](#sfc_get_format_info)                       | Get information about a major or subtype format.        |
| [SFC_GET_FORMAT_MAJOR_COUNT](#sfc_get_format_major_count)         | Get the number of major formats.                        |
| [SFC_GET_FORMAT_MAJOR](#sfc_get_format_major)                     | Get information about a major format type.              |
| [SFC_GET_FORMAT_SUBTYPE_COUNT](#sfc_get_format_subtype_count)     | Get the number of subformats.                           |
| [SFC_GET_FORMAT_SUBTYPE](#sfc_get_format_subtype)                 | Get information about a subformat.                      |
| [SFC_SET_ADD_PEAK_CHUNK](#sfc_set_add_peak_chunk)                 | Control PEAK chunk write to WAV and AIFF.               |
| [SFC_UPDATE_HEADER_NOW](#sfc_update_header_now)                   | Update the file header in write mode on demand.         |
| [SFC_SET_UPDATE_HEADER_AUTO](#sfc_set_update_header_auto)         | Update the file header on each write.                   |
| [SFC_FILE_TRUNCATE](#sfc_file_truncate)                           | Truncate a file open for write or for read/write.       |
| [SFC_SET_RAW_START_OFFSET](#sfc_set_raw_start_offset)             | Change the data start offset for raw files.             |
| SFC_SET_DITHER_ON_WRITE                                           | Not implemented.                                        |
| SFC_SET_DITHER_ON_READ                                            | Not implemented.                                        |
| SFC_GET_DITHER_INFO_COUNT                                         | Not implemented.                                        |
| SFC_GET_DITHER_INFO                                               | Not implemented.                                        |
| [SFC_SET_CLIPPING](#sfc_set_clipping)                             | Control automatic clipping behaviour.                   |
| [SFC_GET_CLIPPING](#sfc_get_clipping)                             | Get current clipping setting.                           |
| [SFC_GET_EMBED_FILE_INFO](#sfc_get_embed_file_info)               | Get information about embedded audio files.             |
| [SFC_WAVEX_GET_AMBISONIC](#sfc_wavex_get_ambisonic)               | Test a WAVEX file for Ambisonic format.                 |
| [SFC_WAVEX_SET_AMBISONIC](#sfc_wavex_set_ambisonic)               | Modify a WAVEX header for Ambisonic format.             |
| [SFC_SET_VBR_ENCODING_QUALITY](#sfc_set_vbr_encoding_quality)     | Set the Variable Bit Rate encoding quality.             |
| [SFC_SET_OGG_PAGE_LATENCY_MS](#sfc_set_ogg_page_latency_ms)       | Set Ogg page latency for Opus file.                     |
| [SFC_GET_OGG_STREAM_SERIALNO](#sfc_get_ogg_stream_serialno)       | Get Ogg stream serial number.                           |
| [SFC_SET_COMPRESSION_LEVEL](#sfc_set_compression_level)           | Set the compression level.                              |
| [SFC_RAW_DATA_NEEDS_ENDSWAP](#sfc_raw_data_needs_endswap)         | Determine if raw data needs endswapping.                |
| [SFC_GET_BROADCAST_INFO](#sfc_get_broadcast_info)                 | Get the Broadcast Chunk info.                           |
| [SFC_SET_BROADCAST_INFO](#sfc_set_broadcast_info)                 | Set the Broadcast Chunk info.                           |
| [SFC_GET_CHANNEL_MAP_INFO](#sfc_get_channel_map_info)             | Get the channel map info.                               |
| [SFC_SET_CHANNEL_MAP_INFO](#sfc_set_channel_map_info)             | Set the channel map info.                               |
| [SFC_SET_CART_INFO](#sfc_set_cart_info)                           | Set the Cart Chunk info.                                |
| [SFC_GET_CART_INFO](#sfc_get_cart_info)                           | Get the Cart Chunk info.                                |
| [SFC_GET_LOOP_INFO](#sfc_get_loop_info)                           | Get loop info.                                          |
| [SFC_GET_INSTRUMENT](#sfc_get_instrument)                         | Get instrument info.                                    |
| [SFC_SET_INSTRUMENT](#sfc_set_instrument)                         | Set instrument info.                                    |
| [SFC_GET_CUE_COUNT](#sfc_get_cue_count)                           | Get the cue marker count.                               |
| [SFC_GET_CUE](#sfc_get_cue)                                       | Get cue marker info.                                    |
| [SFC_SET_CUE](#sfc_set_cue)                                       | Set cue marker info.                                    |
| [SFC_RF64_AUTO_DOWNGRADE](#sfc_rf64_auto_downgrade)               | Set auto downgrade from RF64 to WAV.                    |
| [SFC_GET_ORIGINAL_SAMPLERATE](#sfc_get_original_samplerate)       | Get original samplerate metadata.                       |
| [SFC_SET_ORIGINAL_SAMPLERATE](#sfc_set_original_samplerate)       | Set original samplerate metadata.                       |
| [SFC_GET_BITRATE_MODE](#sfc_get_bitrate_mode)                     | Get bitrate mode.                                       |
| [SFC_SET_BITRATE_MODE](#sfc_set_bitrate_mode)                     | Set bitrate mode.                                       |

---

## SFC_GET_LIB_VERSION

Retrieve the version of the library as a string.

### Parameters

sndfile
: Not used

cmd
: SFC_GET_LIB_VERSION

data
: A pointer to a char buffer

datasize
: The size of the buffer

### Examples

```c
char  buffer [128] ;
sf_command (NULL, SFC_GET_LIB_VERSION, buffer, sizeof (buffer)) ;
```

### Return value

This call will return the length of the retrieved version string.

### Notes

The string returned in the buffer passed to this function will not overflow the
buffer and will always be null terminated .

## SFC_GET_LOG_INFO

Retrieve the internal per-file operation log.

This log buffer can often contain a good reason for why libsndfile failed to
open a particular file.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_LOG_INFO

data
: A pointer to a char buffer

datasize
: The size of the buffer

Example:

```c
char  buffer [2048] ;
sf_command (sndfile, SFC_GET_LOG_INFO, buffer, sizeof (buffer)) ;
```

### Return value

This call will return the length of the retrieved version string.

### Notes

The string returned in the buffer passed to this function will not overflow the
buffer and will always be null terminated.

## SFC_GET_CURRENT_SF_INFO

Retrieve `SF_INFO` struct of opened file.

`SFC_GET_CURRENT_SF_INFO` command copies `SF_INFO` struct of `sndfile` object to
provided buffer.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_CURRENT_SF_INFO

data
: A pointer to a valid SF_INFO* pointer

datasize
: sizeof (SF_INFO)

### Examples

```c
SF_INFO  sfinfo ;
sf_command (sndfile, SFC_GET_CURRENT_SF_INFO, sfinfo, sizeof (SF_INFO)) ;
```

### Return value

Zero on success, non-zero otherwise.

## SFC_CALC_SIGNAL_MAX

Retrieve the measured maximum signal value. This involves reading through the
whole file which can be slow on large files.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_CALC_SIGNAL_MAX

data
: A pointer to a double

datasize
: sizeof (double)

### Examples

```c
double   max_val ;
sf_command (sndfile, SFC_CALC_SIGNAL_MAX, &max_val, sizeof (max_val)) ;
```

### Return value

Zero on success, non-zero otherwise.

## SFC_CALC_NORM_SIGNAL_MAX

Retrieve the measured normalised maximum signal value. This involves reading
through the whole file which can be slow on large files.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_CALC_NORM_SIGNAL_MAX

data
: A pointer to a double

datasize
: sizeof (double)

### Examples

```c
double   max_val ;
sf_command (sndfile, SFC_CALC_NORM_SIGNAL_MAX, &max_val, sizeof (max_val)) ;
```

### Return value

Zero on success, non-zero otherwise.

## SFC_CALC_MAX_ALL_CHANNELS

Calculate the peak value (ie a single number) for each channel. This involves
reading through the whole file which can be slow on large files.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_CALC_MAX_ALL_CHANNELS

data
: A pointer to a double

datasize
: sizeof (double) * number_of_channels

### Examples

```c
double   peaks [number_of_channels] ;
sf_command (sndfile, SFC_CALC_MAX_ALL_CHANNELS, peaks, sizeof (peaks)) ;
```

### Return value

Zero if peaks have been calculated successfully and non-zero otherwise.

## SFC_CALC_NORM_MAX_ALL_CHANNELS

Calculate the normalised peak for each channel. This involves reading through
the whole file which can be slow on large files.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_CALC_NORM_MAX_ALL_CHANNELS

data
: A pointer to a double

datasize
: sizeof (double) * number_of_channels

### Examples

```c
double   peaks [number_of_channels] ;
sf_command (sndfile, SFC_CALC_NORM_MAX_ALL_CHANNELS, peaks, sizeof (peaks)) ;
```

### Return value

Zero if peaks have been calculated successfully and non-zero otherwise.

## SFC_GET_SIGNAL_MAX

Retrieve the peak value for the file as stored in the file header.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_SIGNAL_MAX

data
: A pointer to a double

datasize
: sizeof (double)

### Examples

```c
double   max_peak ;
sf_command (sndfile, SFC_GET_SIGNAL_MAX, &max_peak, sizeof (max_peak)) ;
```

### Return value

SF_TRUE if the file header contained the peak value. SF_FALSE
otherwise.

## SFC_GET_MAX_ALL_CHANNELS

Retrieve the peak value for the file as stored in the file header.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_SIGNAL_MAX

data
: A pointer to an array of doubles

datasize
: sizeof (double) * number_of_channels

### Example

```c
double   peaks [number_of_channels] ;
sf_command (sndfile, SFC_GET_MAX_ALL_CHANNELS, peaks, sizeof (peaks)) ;
```

### Return value

SF_TRUE if the file header contains per channel peak values for the file,
SF_FALSE otherwise.

## SFC_SET_NORM_FLOAT

This command only affects data read from or written to using the
floating point
functions:

```c
size_t    sf_read_float    (SNDFILE *sndfile, float *ptr, size_t items) ;
size_t    sf_readf_float   (SNDFILE *sndfile, float *ptr, size_t frames) ;

size_t    sf_write_float   (SNDFILE *sndfile, float *ptr, size_t items) ;
size_t    sf_writef_float  (SNDFILE *sndfile, float *ptr, size_t frames) ;
```

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_NORM_FLOAT

data
: NULL

datasize
: SF_TRUE or SF_FALSE

For read operations setting normalisation to SF_TRUE means that the data from
all subsequent reads will be be normalised to the range [-1.0, 1.0].

For write operations, setting normalisation to SF_TRUE means than all data
supplied to the float write functions should be in the range [-1.0, 1.0] and
will be scaled for the file format as necessary.

For both cases, setting normalisation to SF_FALSE means that no scaling will
take place.

### Examples

```c
sf_command (sndfile, SFC_SET_NORM_FLOAT, NULL, SF_TRUE) ;

sf_command (sndfile, SFC_SET_NORM_FLOAT, NULL, SF_FALSE) ;
```

### Return value

Returns the previous float normalisation mode.

## SFC_SET_NORM_DOUBLE

This command only affects data read from or written to using the double
precision floating point
functions:

```c
size_t sf_read_double  (SNDFILE *sndfile, double *ptr, size_t items) ;
size_t sf_readf_double (SNDFILE *sndfile, double *ptr, size_t frames) ;

size_t sf_write_double  (SNDFILE *sndfile, double *ptr, size_t items) ;
size_t sf_writef_double (SNDFILE *sndfile, double *ptr, size_t frames) ;
```

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_NORM_DOUBLE

data
: NULL

datasize
: SF_TRUE or SF_FALSE

For read operations setting normalisation to SF_TRUE means that the data from
all subsequent reads will be be normalised to the range [-1.0, 1.0].

For write operations, setting normalisation to SF_TRUE means than all data
supplied to the double write functions should be in the range [-1.0, 1.0] and
will be scaled for the file format as necessary.

For both cases, setting normalisation to SF_FALSE means that no scaling will
take place.

### Examples

```c
sf_command (sndfile, SFC_SET_NORM_DOUBLE, NULL, SF_TRUE) ;

sf_command (sndfile, SFC_SET_NORM_DOUBLE, NULL, SF_FALSE) ;
```

### Return value

Returns the previous double normalisation mode.

## SFC_GET_NORM_FLOAT

Retrieve the current float normalisation mode.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_NORM_FLOAT

data
: NULL

datasize
: anything

### Examples

```c
normalisation = sf_command (sndfile, SFC_GET_NORM_FLOAT, NULL, 0) ;
```

### Return value

Returns TRUE if normalisation is on and FALSE otherwise.

## SFC_GET_NORM_DOUBLE

Retrieve the current float normalisation mode.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_NORM_DOUBLE

data
: NULL

datasize
: anything

Example:

```c
normalisation = sf_command (sndfile, SFC_GET_NORM_DOUBLE, NULL, 0) ;
```

### Return value

Returns TRUE if normalisation is on and FALSE otherwise.

## SFC_SET_SCALE_FLOAT_INT_READ

Set/clear the scale factor when integer (short/int) data is read from a file
containing floating point data.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd:
SFC_SET_SCALE_FLOAT_INT_READ

data
: NULL

datasize
: TRUE or FALSE

Example:

```c
sf_command (sndfile, SFC_SET_SCALE_FLOAT_INT_READ, NULL, SF_TRUE) ;
```

### Return value

Returns the previous `SFC_SET_SCALE_FLOAT_INT_READ` setting for this file.

## SFC_SET_SCALE_INT_FLOAT_WRITE

Set/clear the scale factor when integer (short/int) data is written to a file as
floating point data.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_SCALE_INT_FLOAT_WRITE

data
: NULL

datasize
: TRUE or FALSE

### Examples

```c
sf_command (sndfile, SFC_SET_SCALE_INT_FLOAT_WRITE, NULL, SF_TRUE) ;
```

### Return value

Returns the previous `SFC_SET_SCALE_INT_FLOAT_WRITE` setting for this file.

## SFC_GET_SIMPLE_FORMAT_COUNT

Retrieve the number of simple formats supported by libsndfile.

### Parameters

sndfile
: Not used.

cmd
: SFC_GET_SIMPLE_FORMAT_COUNT

data
: a pointer to an int

datasize
: sizeof (int)

### Examples

```c
int  count ;
sf_command (sndfile, SFC_GET_SIMPLE_FORMAT_COUNT, &count, sizeof (int)) ;
```

### Return value

`0`.

## SFC_GET_SIMPLE_FORMAT

Retrieve information about a simple format.

### Parameters

sndfile
: Not used.

cmd
: SFC_GET_SIMPLE_FORMAT

data
: a pointer to an  SF_FORMAT_INFO struct

datasize
: sizeof (SF_FORMAT_INFO)

The SF_FORMAT_INFO struct is defined in *sndfile.h* as:

```c
typedef struct
{   int         format ;
    const char  *name ;
    const char  *extension ;
} SF_FORMAT_INFO ;
```

When `sf_command()` is called with `SF_GET_SIMPLE_FORMAT`, the value of the
format field should be the format number (ie 0 \<= format \<= count value
obtained using `SF_GET_SIMPLE_FORMAT_COUNT).

### Examples

```c
SF_FORMAT_INFO  format_info ;
int             k, count ;

sf_command (sndfile, SFC_GET_SIMPLE_FORMAT_COUNT, &count, sizeof (int)) ;

for (k = 0 ; k < count ; k++)
{   format_info.format = k ;
    sf_command (sndfile, SFC_GET_SIMPLE_FORMAT, &format_info, sizeof (format_info)) ;
    printf ("%08x  %s %s\n", format_info.format, format_info.name, format_info.extension) ;
    } ;
```

### Return value

0 on success and non-zero otherwise.

The value of the format field of the `SF_FORMAT_INFO` struct will be a value
which can be placed in the format field of an `SF_INFO` struct when a file is to
be opened for write. The name field will contain a char\* pointer to the name of
the string, eg. "WAV (Microsoft 16 bit PCM)". The extension field will contain
the most commonly used file extension for that file type.

## SFC_GET_FORMAT_INFO

Retrieve information about a major or subtype format.

### Parameters

sndfile
: Not used.

cmd
: SFC_GET_FORMAT_INFO

data
: a pointer to an SF_FORMAT_INFO struct

datasize
: sizeof (SF_FORMAT_INFO)

The `SF_FORMAT_INFO` struct is defined in \<sndfile.h\> as:

```c
typedef struct
{   int         format ;
    const char  *name ;
    const char  *extension ;
} SF_FORMAT_INFO ;
```

When `sf_command()` is called with `SF_GET_FORMAT_INFO`, the format field is
examined and if (format & `SF_FORMAT_TYPEMASK`) is a valid format then the
struct is filled in with information about the given major type. If (format &
`SF_FORMAT_TYPEMASK`) is FALSE and (format & `SF_FORMAT_SUBMASK`) is a valid
subtype format then the struct is filled in with information about the given
subtype.

### Examples

```c
SF_FORMAT_INFO  format_info ;

format_info.format = SF_FORMAT_WAV ;
sf_command (sndfile, SFC_GET_FORMAT_INFO, &format_info, sizeof (format_info)) ;
printf ("%08x  %s %s\n", format_info.format, format_info.name, format_info.extension) ;

format_info.format = SF_FORMAT_ULAW ;
sf_command (sndfile, SFC_GET_FORMAT_INFO, &format_info, sizeof (format_info)) ;
printf ("%08x  %s\n", format_info.format, format_info.name) ;
```

### Return value

0 on success and non-zero otherwise.

## SFC_GET_FORMAT_MAJOR_COUNT

Retrieve the number of major formats.

### Parameters

sndfile
: Not used.

cmd
: SFC_GET_FORMAT_MAJOR_COUNT

data
: a pointer to an int

datasize
: sizeof (int)

### Examples

```c
int  count ;
sf_command (sndfile, SFC_GET_FORMAT_MAJOR_COUNT, &count, sizeof (int)) ;
```

### Return value

0.

## SFC_GET_FORMAT_MAJOR

Retrieve information about a major format type.

### Parameters

sndfile
: Not used.

cmd
: SFC_GET_FORMAT_MAJOR

data
: a pointer to an  SF_FORMAT_INFO struct

datasize
: sizeof (SF_FORMAT_INFO)

### Examples

```c
SF_FORMAT_INFO  format_info ;
int             k, count ;

sf_command (sndfile, SFC_GET_FORMAT_MAJOR_COUNT, &count, sizeof (int)) ;

for (k = 0 ; k < count ; k++)
{   format_info.format = k ;
    sf_command (sndfile, SFC_GET_FORMAT_MAJOR, &format_info, sizeof (format_info)) ;
    printf ("%08x  %s %s\n", format_info.format, format_info.name, format_info.extension) ;
    } ;
```

For a more comprehensive example, see the program `list_formats.c` in the
`examples/` directory of the libsndfile source code distribution.

### Return value

0 on success and non-zero otherwise.

The value of the format field will be one of the major format identifiers such
as `SF_FORMAT_WAV` or `SF_FORMAT`_AIFF. The name field will contain a char\*
pointer to the name of the string, eg. "WAV (Microsoft)". The extension field
will contain the most commonly used file extension for that file type.

## SFC_GET_FORMAT_SUBTYPE_COUNT

Retrieve the number of subformats.

### Parameters

sndfile
: Not used.

cmd
: SFC_GET_FORMAT_SUBTYPE_COUNT

data
: a pointer to an int

datasize
: sizeof (int)

### Examples

```c
int   count ;
sf_command (sndfile, SFC_GET_FORMAT_SUBTYPE_COUNT, &count, sizeof (int)) ;
```

### Return value

Returns zero.

## SFC_GET_FORMAT_SUBTYPE

Enumerate the subtypes (this function does not translate a subtype into a string
describing that subtype). A typical use case might be retrieving a string
description of all subtypes so that a dialog box can be filled in.

### Parameters

sndfile
: Not used.

cmd
: SFC_GET_FORMAT_SUBTYPE

data
: a pointer to an SF_FORMAT_INFO struct

datasize
: sizeof (SF_FORMAT_INFO)

### Examples

Example 1: Retrieve all sybtypes supported by the WAV format.

```c
SF_FORMAT_INFO  format_info ;
int             k, count ;

sf_command (sndfile, SFC_GET_FORMAT_SUBTYPE_COUNT, &count, sizeof (int)) ;

for (k = 0 ; k < count ; k++)
{   format_info.format = k ;
    sf_command (sndfile, SFC_GET_FORMAT_SUBTYPE, &format_info, sizeof (format_info)) ;
    if (! sf_format_check (format_info.format | SF_FORMAT_WAV))
        continue ;
    printf ("%08x  %s\n", format_info.format, format_info.name) ;
    } ;
```

Example 2: Print a string describing the `SF_FORMAT_PCM_16` subtype.

```c
SF_FORMAT_INFO  format_info ;
int             k, count ;

sf_command (sndfile, SFC_GET_FORMAT_SUBTYPE_COUNT, &count, sizeof (int)) ;

for (k = 0 ; k < count ; k++)
{   format_info.format = k ;
    sf_command (sndfile, SFC_GET_FORMAT_SUBTYPE, &format_info, sizeof (format_info)) ;
    if (format_info.format == SF_FORMAT_PCM_16)
    {   printf ("%08x  %s\n", format_info.format, format_info.name) ;
        break ;
        } ;
    } ;
```

For a more comprehensive example, see the program `list_formats.c` in the
`examples/` directory of the libsndfile source code distribution.

### Return value

0 on success and non-zero otherwise.

The value of the format field will be one of the major format identifiers such
as `SF_FORMAT_WAV` or `SF_FORMAT_AIFF`. The name field will contain a char\*
pointer to the name of the string; for instance "WAV (Microsoft)" or "AIFF
(Apple/SGI)". The extension field will be a NULL pointer.

## SFC_SET_ADD_PEAK_CHUNK

By default, WAV and AIFF files which contain floating point data (subtype
`SF_FORMAT_FLOAT` or `SF_FORMAT_DOUBLE`) have a PEAK chunk. By using this
command, the addition of a PEAK chunk can be turned on or off.

**Note**: This call must be made before any data is written to the file.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_ADD_PEAK_CHUNK

data
: Not used (should be NULL)

datasize
: TRUE or FALSE.

### Examples

```c
/* Turn on the PEAK chunk. */
sf_command (sndfile, SFC_SET_ADD_PEAK_CHUNK, NULL, SF_TRUE) ;

/* Turn off the PEAK chunk. */
sf_command (sndfile, SFC_SET_ADD_PEAK_CHUNK, NULL, SF_FALSE) ;
```

### Return value

Returns SF_TRUE if the peak chunk will be written after this call. Returns
SF_FALSE if the peak chunk will not be written after this call.

## SFC_UPDATE_HEADER_NOW

The header of an audio file is normally written by libsndfile when the file is
closed using [**sf_close()**](api.md#file-close-function).

There are however situations where large files are being generated and it would
be nice to have valid data in the header before the file is complete. Using this
command will update the file header to reflect the amount of data written to the
file so far. Other programs opening the file for read (before any more data is
written) will then read a valid sound file header.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_UPDATE_HEADER_NOW

data
: Not used (should be NULL)

datasize
: Not used.

### Examples

```c
/* Update the header now. */
sf_command (sndfile, SFC_UPDATE_HEADER_NOW, NULL, 0) ;
```

### Return value

Returns zero.

## SFC_SET_UPDATE_HEADER_AUTO

Similar to `SFC_UPDATE_HEADER_NOW` but updates the header at the end of every
call to the [sf_write\*](api.md#write) functions.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_UPDATE_HEADER_AUTO

data
: Not used (should be NULL)

datasize
: `SF_TRUE` or `SF_FALSE`

### Examples

```c
/* Turn on auto header update. */
sf_command (sndfile, SFC_SET_UPDATE_HEADER_AUTO, NULL, SF_TRUE) ;

/* Turn off auto header update. */
sf_command (sndfile, SFC_SET_UPDATE_HEADER_AUTO, NULL, SF_FALSE) ;
```

### Return value

TRUE if auto update header is now on; FALSE otherwise.

## SFC_FILE_TRUNCATE

Truncate a file that was opened for write or read/write.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_FILE_TRUNCATE

data
: A pointer to an sf_count_t.

datasize
: sizeof (sf_count_t)

Truncate the file to the number of frames specified by the sf_count_t pointed to
by data. After this command, both the read and the write pointer will be at the
new end of the file. This command will fail (returning non-zero) if the
requested truncate position is beyond the end of the file.

### Examples

```c
/* Truncate the file to a length of 20 frames. */
sf_count_t  frames = 20 ;
sf_command (sndfile, SFC_FILE_TRUNCATE, &frames, sizeof (frames)) ;
```

### Return value

Zero on success, non-zero otherwise.

## SFC_SET_RAW_START_OFFSET

Change the data start offset for files opened up as `SF_FORMAT_RAW`.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_RAW_START_OFFSET

data
: A pointer to an sf_count_t.

datasize
: sizeof (sf_count_t)

For a file opened as format `SF_FORMAT_RAW`, set the data offset to the value
given by `data`.

### Examples

```c
/* Reset the data offset to 5 bytes from the start of the file. */
sf_count_t  offset = 5 ;
sf_command (sndfile, SFC_SET_RAW_START_OFFSET, &offset, sizeof (offset)) ;
```

### Return value

Zero on success, non-zero otherwise.

## SFC_SET_CLIPPING

Turn on/off automatic clipping when doing floating point to integer conversion.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_CLIPPING

data
: NULL

datasize
: SF_TRUE or SF_FALSE.

Turn on (datasize == SF_TRUE) or off (datasize == SF_FALSE) clipping.

### Examples

```c
sf_command (sndfile, SFC_SET_CLIPPING, NULL, SF_TRUE) ;
```

### Return value

Clipping mode (SF_TRUE or SF_FALSE).

## SFC_GET_CLIPPING

Turn on/off automatic clipping when doing floating point to integer conversion.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_CLIPPING

data
: NULL

datasize
: 0

Retrieve the current cliiping setting.

### Examples

```c
sf_command (sndfile, SFC_GET_CLIPPING, NULL, 0) ;
```

### Return value

Clipping mode (SF_TRUE or SF_FALSE).

## SFC_GET_EMBED_FILE_INFO

Get the file offset and file length of a file enbedded within another larger
file.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_EMBED_FILE_INFO

data
: a pointer to an  SF_EMBED_FILE_INFO struct

datasize
: sizeof (SF_EMBED_FILE_INFO)

The `SF_EMBED_FILE_INFO` struct is defined in *sndfile.h* as:

```c
typedef struct
{   sf_count_t  offset ;
    sf_count_t  length ;
} SF_EMBED_FILE_INFO ;
```

### Return value

0 on success and non-zero otherwise.

The value of the offset field of the `SF_EMBED_FILE_INFO` struct will be the
offsets in bytes from the start of the outer file to the start of the audio
file. The value of the offset field of the `SF_EMBED_FILE_INFO` struct will be
the length in bytes of the embedded file.

## SFC_WAVEX_GET_AMBISONIC

Test if the current file has the GUID of a WAVEX file for any of the Ambisonic
formats.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_WAVEX_GET_AMBISONIC

data
: NULL

datasize
: 0

The Ambisonic WAVEX formats are defined here:
<http://dream.cs.bath.ac.uk/researchdev/wave-ex/bformat.html>.

### Return value

`SF_AMBISONIC_NONE(0x40)` or `SF_AMBISONIC_B_FORMAT(0x41)` or zero if the file
format does not support ambisonic formats.

## SFC_WAVEX_SET_AMBISONIC

Set the GUID of a new WAVEX file to indicate an Ambisonics format.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_WAVEX_SET_AMBISONIC

data
: NULL

datasize
: SF_AMBISONIC_NONE or SF_AMBISONIC_B_FORMAT

Turn on (`SF_AMBISONIC_B_FORMAT(0x41)`) or off (`SF_AMBISONIC_NONE(0x40)`)
encoding. This command is currently only supported for files with
`SF_FORMAT_WAVEX` format.

The Ambisonic WAVEX formats are defined here:

<http://dream.cs.bath.ac.uk/researchdev/wave-ex/bformat.html>.

### Return value

Return the ambisonic value that has just been set or zero if the
file format does not support ambisonic encoding.

## SFC_SET_VBR_ENCODING_QUALITY

Set the Variable Bit Rate encoding quality. The encoding quality value
should be between 0.0 (lowest quality) and 1.0 (highest quality).
Currently this command is only implemented for FLAC and Ogg/Vorbis files.
It has no effect on un-compressed file formats.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_VBR_ENCODING_QUALITY

data
: A pointer to a double value

datasize
: sizeof (double)

The command must be sent before any audio data is written to the file.

### Return value

SF_TRUE if VBR encoding quality was set. SF_FALSE otherwise.

## SFC_SET_OGG_PAGE_LATENCY_MS

Set page latency for Ogg Opus file in milliseconds. The value should be between
50.0 and 1600.0. This command is only implemented for Ogg Opus files.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_OGG_PAGE_LATENCY_MS

data
: A pointer to a double value

datasize
: sizeof (double)

### Return value

0 on success and non-zero otherwise.

## SFC_GET_OGG_STREAM_SERIALNO

Get the Ogg stream serial number for files with the Ogg major format. Ogg
stream serail numbers are a randomly chosen 32-bit value, used for
differentiating logical Ogg streams.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_OGG_STREAM_SERIALNO

data
: A pointer to a 32-bit int value

datasize
: sizeof (int32_t) = 4

### Return value

0 on success and non-zero otherwise.

## SFC_SET_COMPRESSION_LEVEL

Set the compression level. The compression level should be between 0.0 (minimum
compression level) and 1.0 (highest compression level). Currently this command is
only implemented for FLAC and Ogg/Vorbis files. It has no effect on
uncompressed file formats.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_COMPRESSION_LEVEL

data
: A pointer to a double value

datasize
: sizeof (double)

The command must be sent before any audio data is written to the file.

### Return value

SF_TRUE if compression level was set. SF_FALSE otherwise.

## SFC_RAW_DATA_NEEDS_ENDSWAP

Determine if raw data read using [sf_read_raw()](api.md#raw) needs to be end
swapped on the host CPU.

For instance, will return SF_TRUE on when reading WAV containing
`SF_FORMAT_PCM_16` data on a big endian machine and `SF_FALSE` on a
little endian machine.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_RAW_DATA_NEEDS_ENDSWAP

data
: NULL

datasize
: 0

### Return value

`SF_TRUE` or `SF_FALSE`.

## SFC_GET_BROADCAST_INFO

Retrieve the Broadcast Extension Chunk from WAV (and related) files.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_BROADCAST_INFO

data
: a pointer to an SF_BROADCAST_INFO struct

datasize
: sizeof (SF_BROADCAST_INFO)

The SF_BROADCAST_INFO struct is defined in *sndfile.h* as:

```c
typedef struct
{   char            description [256] ;
    char            originator [32] ;
    char            originator_reference [32] ;
    char            origination_date [10] ;
    char            origination_time [8] ;
    unsigned int    time_reference_low ;
    unsigned int    time_reference_high ;
    short           version ;
    char            umid [64] ;
    char            reserved [190] ;
    unsigned int    coding_history_size ;
    char            coding_history [256] ;
} SF_BROADCAST_INFO ;
```

### Return value

`SF_TRUE` if the file contained a Broadcast Extension chunk or `SF_FALSE`
otherwise.

## SFC_SET_BROADCAST_INFO

Set the Broadcast Extension Chunk for WAV (and related) files.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_BROADCAST_INFO

data
: a pointer to an SF_BROADCAST_INFO struct

datasize
: sizeof (SF_BROADCAST_INFO)

### Return value

`SF_TRUE` if setting the Broadcast Extension chunk was successful and `SF_FALSE`
otherwise.

## SFC_GET_CHANNEL_MAP_INFO

Retrieve the channel map contained in an AIFF or CAF Channel Layout chunk.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_CHANNEL_MAP_INFO

data
: a pointer to an array of int, the same size as the number of channels in the
file

datasize
: number of channels * sizeof (int)

Channel map positions are defined in an enum in *sndfile.h*:

| Name                                 | Value | Description                                                       |
|:-------------------------------------|:------|:------------------------------------------------------------------|
| SF_CHANNEL_MAP_INVALID               | 0     |                                                                   |
| SF_CHANNEL_MAP_MONO                  | 1     |                                                                   |
| SF_CHANNEL_MAP_LEFT                  | 2     | Apple calls this 'Left'                                           |
| SF_CHANNEL_MAP_RIGHT                 | 3     | Apple calls this 'Right'                                          |
| SF_CHANNEL_MAP_CENTER                | 4     | Apple calls this 'Center'                                         |
| SF_CHANNEL_MAP_FRONT_LEFT            | 5     |                                                                   |
| SF_CHANNEL_MAP_FRONT_RIGHT           | 6     |                                                                   |
| SF_CHANNEL_MAP_FRONT_CENTER          | 7     |                                                                   |
| SF_CHANNEL_MAP_REAR_CENTER           | 8     | Apple calls this 'Center Surround', Msft calls this 'Back Center' |
| SF_CHANNEL_MAP_REAR_LEFT             | 9     | Apple calls this 'Left Surround', Msft calls this 'Back Left'     |
| SF_CHANNEL_MAP_REAR_RIGHT            | 10    | Apple calls this 'Right Surround', Msft calls this 'Back Right'   |
| SF_CHANNEL_MAP_LFE                   | 11    | Apple calls this 'LFEScreen', Msft calls this 'Low Frequency'     |
| SF_CHANNEL_MAP_FRONT_LEFT_OF_CENTER  | 12    | Apple calls this 'Left Center'                                    |
| SF_CHANNEL_MAP_FRONT_RIGHT_OF_CENTER | 13    | Apple calls this 'Right Center'                                   |
| SF_CHANNEL_MAP_SIDE_LEFT             | 14    | Apple calls this 'Left Surround Direct'                           |
| SF_CHANNEL_MAP_SIDE_RIGHT            | 15    | Apple calls this 'Right Surround Direct'                          |
| SF_CHANNEL_MAP_TOP_CENTER            | 16    | Apple calls this 'Top Center Surround'                            |
| SF_CHANNEL_MAP_TOP_FRONT_LEFT        | 17    | Apple calls this 'Vertical Height Left'                           |
| SF_CHANNEL_MAP_TOP_FRONT_RIGHT       | 18    | Apple calls this 'Vertical Height Right'                          |
| SF_CHANNEL_MAP_TOP_FRONT_CENTER      | 19    | Apple calls this 'Vertical Height Center'                         |
| SF_CHANNEL_MAP_TOP_REAR_LEFT         | 20    | Apple and MS call this 'Top Back Left'                            |
| SF_CHANNEL_MAP_TOP_REAR_RIGHT        | 21    | Apple and MS call this 'Top Back Right'                           |
| SF_CHANNEL_MAP_TOP_REAR_CENTER       | 22    | Apple and MS call this 'Top Back Center'                          |
| SF_CHANNEL_MAP_AMBISONIC_B_W         | 23    |                                                                   |
| SF_CHANNEL_MAP_AMBISONIC_B_X         | 24    |                                                                   | 
| SF_CHANNEL_MAP_AMBISONIC_B_Y         | 25    |                                                                   |
| SF_CHANNEL_MAP_AMBISONIC_B_Z         | 26    |                                                                   |
| SF_CHANNEL_MAP_MAX                   | 27    |                                                                   |

### Return value

`SF_TRUE` if the file contained a Channel Layout chunk or `SF_FALSE` otherwise.

## SFC_SET_CHANNEL_MAP_INFO

Set the channel map contained in an AIFF or CAF Channel Layout chunk.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_CHANNEL_MAP_INFO

data
: a pointer to an array of int, the same size as the number of channels in the
file

datasize
: number of channels * sizeof (int)

### Return value

`SF_TRUE` if setting the Channel Layout chunk was successful and `SF_FALSE`
otherwise.

## SFC_GET_CART_INFO

Retrieve the Cart Chunk from WAV (and related) files. Based on AES46 standard
for CartChunk (see [CartChunk.org](http://www.cartchunk.org/) for more
information.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_CART_INFO

data
: a pointer to an SF_CART_INFO struct

datasize
: sizeof (SF_CART_INFO)

The SF_CART_INFO struct is defined in *sndfile.h* as:

```c
#define SF_CART_INFO_VAR(p_tag_text_size) \
    struct
    {       char            version [4] ;
            char            title [64] ;
            char            artist [64] ;
            char            cut_id [64] ;
            char            client_id [64] ;
            char            category [64] ;
            char            classification [64] ;
            char            out_cue [64] ;
            char            start_date [10] ;
            char            start_time [8] ;
            char            end_date [10] ;
            char            end_time [8] ;
            char            producer_app_id [64] ;
            char            producer_app_version [64] ;
            char            user_def [64] ;
            long    level_reference ;
            SF_CART_TIMER   post_timers [8] ;
            char            reserved [276] ;
            char            url [1024] ;
            unsigned int    tag_text_size ;
            char            tag_text[p_tag_text_size] ;
    }
```

### Return value

`SF_TRUE` if the file contained a Cart chunk or `SF_FALSE` otherwise.

## SFC_SET_CART_INFO

Set the Cart Chunk for WAV (and related) files.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_CART_INFO

data
: a pointer to an SF_CART_INFO struct

datasize
: sizeof (SF_CART_INFO)

### Return value

SF_TRUE if setting the Cart chunk was successful and SF_FALSE otherwise.

## SFC_GET_LOOP_INFO

Retrieve loop information for file including time signature, length in beats and
original MIDI base note

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_LOOP_INFO

data
: a pointer to an SF_LOOP_INFO struct

datasize
: sizeof (SF_LOOP_INFO)

The SF_LOOP_INFO struct is defined in *sndfile.h* as:

```c
typedef struct
{   short    time_sig_num ;   /* any positive integer    > 0  */
    short    time_sig_den ;   /* any positive power of 2 > 0  */
    int        loop_mode ;    /* see SF_LOOP enum             */

    int        num_beats ;    /* this is NOT the amount of quarter notes !!!*/
                              /* a full bar of 4/4 is 4 beats */
                              /* a full bar of 7/8 is 7 beats */

    float    bpm ;            /* suggestion, as it can be calculated using other fields:*/
                              /* file's length, file's sampleRate and our time_sig_den*/
                              /* -> bpms are always the amount of _quarter notes_ per minute */

    int    root_key ;         /* MIDI note, or -1 for None */
    int future [6] ;
} SF_LOOP_INFO ;
```

### Examples

```c
SF_LOOP_INFO loop;
sf_command (sndfile, SFC_GET_LOOP_INFO, &loop, sizeof (loop)) ;
```

### Return value

`SF_TRUE` if the file header contains loop information for the file, `SF_FALSE`
otherwise.

## SFC_GET_INSTRUMENT

Retrieve instrument information from file including MIDI base note, keyboard
mapping and looping information (start/stop and mode).

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_INSTRUMENT

data
: a pointer to an SF_INSTRUMENT struct

datasize
: sizeof (SF_INSTRUMENT)

The `SF_INSTRUMENT` struct is defined in *sndfile.h* as:

```c
typedef struct
{   int gain ;
    char basenote, detune ;
    char velocity_lo, velocity_hi ;
    char key_lo, key_hi ;
    int loop_count ;

    struct
    {   int mode ;
        unsigned int start ;
        unsigned int end ;
        unsigned int count ;
    } loops [16] ; /* make variable in a sensible way */
} SF_INSTRUMENT ;
```

`mode` values are defined as:

| Name                | Value | Description |
|:--------------------|:------|:------------|
| SF_LOOP_NONE        | 800   |             |
| SF_LOOP_FORWARD     | 801   |             |
| SF_LOOP_BACKWARD    | 802   |             |
| SF_LOOP_ALTERNATING | 803   |             |

### Examples

```c
SF_INSTRUMENT inst ;
sf_command (sndfile, SFC_GET_INSTRUMENT, &inst, sizeof (inst)) ;
```

### Return value

`SF_TRUE` if the file header contains instrument information for the file,
`SF_FALSE` otherwise.

## SFC_SET_INSTRUMENT

Set the instrument information for the file.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_INSTRUMENT

data
: a pointer to an SF_INSTRUMENT struct

datasize
: sizeof (SF_INSTRUMENT)

### Examples

```c
SF_INSTRUMENT inst ;
sf_command (sndfile, SFC_SET_INSTRUMENT, &inst, sizeof (inst)) ;
```

### Return value

`SF_TRUE` if the file header contains instrument information for the file,
`SF_FALSE` otherwise.

## SFC_GET_CUE_COUNT

Retrieve the number of cue markers available for retrieval using the
[SFC_GET_CUE](#sfc_get_cue) command.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_CUE_COUNT

data
: a pointer to a uint32_t

datasize
: sizeof (uint32_t)

### Examples

```c
uint32_t cue_count ;
sf_command (sndfile, SFC_GET_CUE_COUNT, &cue_count, sizeof (cue_count)) ;
```

### Return value

`SF_TRUE` if the file header contains cue marker information for the file,
`SF_FALSE` otherwise.

## SFC_GET_CUE

Retrieve cue marker information from file.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_CUE

data
: a pointer to an SF_CUES struct

datasize
: sizeof (SF_CUES)

The SF_CUES struct is defined in *sndfile.h* as:

```c
typedef struct
{   int cue_count ;

    struct
    {   int32_t   indx ;
        uint32_t  position ;
        int32_t   fcc_chunk ;
        int32_t   chunk_start ;
        int32_t   block_start ;
        uint32_t  sample_offset ;
        char name [256] ;
    } cue_points [100] ;
} SF_CUES ;
```

There is also an SF_CUES_VAR \#define that allows reading/writing more than 100
cue markers.

### Examples

```c
SF_CUES cues ;
sf_command (sndfile, SFC_GET_CUE, &cues, sizeof (cues)) ;
```

### Return value

`SF_TRUE` if the file header contains cue marker information for the file,
`SF_FALSE` otherwise.

## SFC_SET_CUE

Set the cue marker information for the file.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_CUE

data
: a pointer to an SF_CUES struct

datasize
: sizeof (SF_CUES)

### Examples

```c
SF_CUES cues ;
sf_command (sndfile, SFC_SET_CUE, &cues, sizeof (cues)) ;
```

### Return value

`SF_TRUE` if the file header contains cue marker information for the file,
`SF_FALSE` otherwise.

## SFC_RF64_AUTO_DOWNGRADE

Enable auto downgrade from RF64 to WAV.

The EBU recommendation is that when writing RF64 files and the resulting file is
less than 4Gig in size, it should be downgraded to a WAV file (WAV files have a
maximum size of 4Gig). libsndfile doesn't follow the EBU recommendations
exactly, mainly because the test suite needs to be able test reading/writing
RF64 files without having to generate files larger than 4 gigabytes.

Note: This command should be issued before the first bit of audio data has been
written to the file. Calling this command after audio data has been written will
return the current value of this setting, but will not allow it to be changed.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_RF64_AUTO_DOWNGRADE

data
: NULL

datasize
: SF_TRUE or SF_FALSE

### Examples

```c
/* Enable auto downgrade on file close. */
sf_command (sndfile, SFC_RF64_AUTO_DOWNGRADE, NULL, SF_TRUE) ;
```

### Return value

Returns `SF_TRUE` if `SFC_RF64_AUTO_DOWNGRADE` is set and `SF_FALSE` otherwise.

## SFC_GET_ORIGINAL_SAMPLERATE

Get original samplerate metadata.

The Opus audio codec stores audio data independent of samplerate, but only
supports encoding or decoding at 8000Hz, 12000Hz, 16000Hz, 24000Hz or 48000Hz.
Opus includes a header field to record the original source input samplerate, and
a samplerate converter may be used if needed.

This command gets the original samplerate header field. It does not enable any
(non-existent) samplerate conversion, nor change the current decoder samplerate.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_ORIGINAL_SAMPLERATE

data
: pointer to an integer

datasize
: sizeof (int)

### Examples

```c
/* Get the original sample rate */
int original_samplerate ;
sf_command (sndfile, SFC_GET_ORIGINAL_SAMPLERATE, &original_samplerate, sizeof (original_samplerate)) ;
```

### Return value

Returns `SF_TRUE` on success, `SF_FALSE` otherwise.

The passed integer is set to the value of the original samplerate.

## SFC_SET_ORIGINAL_SAMPLERATE

Set original samplerate metadata.

The Opus audio codec stores audio data independent of samplerate, but only
supports encoding or decoding at 8000Hz, 12000Hz, 16000Hz, 24000Hz or 48000Hz.
Opus includes a header field to record the original source input samplerate, and
a samplerate converter may be used if needed.

When writing an Opus file this command sets the original samplerate header field
to the provided value, which is then stored in the file. This has no effect on
the current encoder samplerate.

When reading an Opus file this command overrides the original samplerate value
as read from the file. libsndfile uses this value to choose what samplerate to
decode at, rounding up to the nearest valid Opus samplerate. After a successful
call, the file samplerate and frames count may have changed.

Note: This command should be issued before the first bit of audio data has been
read from or written to the file.

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_ORIGINAL_SAMPLERATE

data
: pointer to an integer

datasize
: sizeof (int)

### Examples

```c
/* Store the original sample rate as 44100 */
int original_samplerate 44100;
sf_command (sndfile, SFC_SET_ORIGINAL_SAMPLERATE, &original_samplerate, sizeof (input_samplerate)) ;
```

### Return value

Returns SF_TRUE on success, SF_FALSE otherwise.

On write, can only succeed if no data has been written. On read, if successful,
[SFC_GET_CURRENT_SF_INFO](#sfc_get_current_sf_info) should be called to
determine the new frames count and samplerate

## SFC_GET_BITRATE_MODE

Get bitrate mode.

The bitrate mode is one of:

| Name                     | Value | Description       |
|:-------------------------|:------|:------------------|
| SF_BITRATE_MODE_CONSTANT | 800   | Constant bitrate. |
| SF_BITRATE_MODE_AVERAGE  | 801   | Average bitrate.  |
| SF_BITRATE_MODE_VARIABLE | 802   | Variable bitrate. |

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_GET_BITRATE_MODE

data
: NULL

datasize
: anything

### Return value

Returns one of `SF_BITRATE_MODE_XXX` on success, `-1` otherwise.

## SFC_SET_BITRATE_MODE

Set bitrate mode.

The bitrate mode is one of:

| Name                     | Value | Description       |
|:-------------------------|:------|:------------------|
| SF_BITRATE_MODE_CONSTANT | 800   | Constant bitrate. |
| SF_BITRATE_MODE_AVERAGE  | 801   | Average bitrate.  |
| SF_BITRATE_MODE_VARIABLE | 802   | Variable bitrate. |

### Parameters

sndfile
: A valid SNDFILE* pointer

cmd
: SFC_SET_BITRATE_MODE

data
: pointer to an integer

datasize
: sizeof (int)

### Return value

Returns `SF_TRUE` on success, `SF_FALSE` otherwise.

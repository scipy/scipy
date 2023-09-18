---
layout: page
title: The libsndfile API
---

# libsndfile

Libsndfile is a library designed to allow the reading and writing of many different sampled sound file formats (such as
MS Windows WAV and the Apple/SGI AIFF format) through one standard library interface.

During read and write operations, formats are seamlessly converted between the format the application program has
requested or supplied and the file's data format. The application programmer can remain blissfully unaware of issues
such as file endian-ness and data format. See [Note 1](#note-1) and [Note 2](#note-2).

Every effort is made to keep these documents up-to-date, error free and unambiguous. However, since maintaining the
documentation is the least fun part of working on libsndfile, these docs can and do fall behind the behaviour of the
library. If any errors, omissions or ambiguities are found, please notify me (erikd) at mega-nerd dot com.

To supplement this reference documentation, there are simple example programs included in the source code tarball. The
test suite which is also part of the source code tarball is also a good place to look for the correct usage of the
library functions.

**Finally, if you think there is some feature missing from libsndfile, check that it isn't already implemented (and
documented) [here](command.md).**

## Synopsis

```c
#include <stdio.h>;
#include <sndfile.h>;
```

| Name                                                                                                        | Description                                    |
|:------------------------------------------------------------------------------------------------------------|:---------------------------------------        |
| [sf_open, sf_wchar_open](#open)                                                                             | File open functions.                           |
| [sf_open_fd](#open_fd)                                                                                      | Open sound file using file descriptor.         |
| [sf_open_virtual](#open_virtual)                                                                            | Open sound file using virtual API.             |
| [sf_format_check](#check)                                                                                   | Validate sound file info.                      |
| [sf_seek](#seek)                                                                                            | Seek position in sound file.                   |
| [sf_command](command.md)                                                                                    | Command interface.                             |
| [sf_error, sf_strerror, sf_error_number, sf_perror, sf_error_str](#error)                                   | Error functions.                               |
| [sf_close](#close)                                                                                          | File close function.                           |
| [sf_write_sync](#write_sync)                                                                                | Write sync function.                           |
| [sf_read_short, sf_read_int, sf_read_float, sf_read_double](#read)                                          | File items read functions.                     |
| [sf_readf_short, sf_readf_int, sf_readf_float, sf_readf_double](#readf)                                     | File frames read functions.                    |
| [sf_write_short, sf_write_int, sf_write_float, sf_write_double](#write)                                     | File items write functions.                    |
| [sf_writef_short, sf_writef_int, sf_writef_float, sf_writef_double](#writef)                                | File frames write functions.                   |
| [sf_read_raw, sf_write_raw](#raw)                                                                           | Raw read/write functions.                      |
| [sf_get_string, sf_set_string](#string)                                                                     | Functions for reading and writing string data. |
| [sf_version_string](#version_string)                                                                        | Retrieve library version string.                |
| [sf_current_byterate](#current_byterate)                                                                    | Retrieve current byterate.                     |
| [sf_set_chunk, sf_get_chunk_iterator, sf_next_chunk_iterator, sf_get_chunk_size, sf_get_chunk_data](#chunk) | RIFF chunks API.                               |

SNDFILE* is an anonymous pointer to data which is private to the library.

## File Open Function {#open}

```c
SNDFILE*  sf_open    (const char *path, int mode, SF_INFO *sfinfo) ;
```

The sf_open() function opens the sound file at the specified path. The filename is byte encoded, but may be utf-8 on
Linux, while on Mac OS X it will use the filesystem character set. On Windows, there is also a Windows specific
sf_wchar_open() that takes a UTF16_BE encoded filename.

```c
SNDFILE*  sf_wchar_open (LPCWSTR wpath, int mode, SF_INFO *sfinfo) ;
```

The SF_INFO structure is for passing data between the calling function and the library when opening a file for reading
or writing. It is defined in sndfile.h as follows:

```c
typedef struct
{    sf_count_t  frames ;     /* Used to be called samples. */
        int         samplerate ;
        int         channels ;
        int         format ;
        int         sections ;
        int         seekable ;
    } SF_INFO ;
```

The mode parameter for this function can be any one of the following three values:

SFM_READ
: read only mode

SFM_WRITE
: write only mode

SFM_RDWR
: read/write mode

When opening a file for read, the **format** field should be set to zero before
calling **sf_open**(). The only exception to this is the case of RAW files where
the caller has to set the **samplerate**, **channels** and **format** fields to
valid values. All other fields of the structure are filled in by the library.

**Note:** The libsndfile library will reject values ​​for field **channels** that
are greater than `1024`. These value ​​represent the maximum theoretical limit
and may be less for specific formats.

When opening a file for write, the caller must fill in structure members
**samplerate**, **channels**, and **format**.

The **format** field in the above **SF_INFO** structure is made up of the
bit-wise OR of a major format type (values between 0x10000 and 0x08000000), a
minor format type (with values less than 0x10000) and an optional endian-ness
value. The currently understood formats are listed in *sndfile.h* as follows and
also include bitmasks for separating major and minor file types. Not all
combinations of endian-ness and major and minor file types are valid.

| Name                     | Value      | Description                                |
|:-------------------------|:-----------|:-------------------------------------------|
| **Major formats.**                                                                 |
| SF_FORMAT_WAV            | 0x010000   | Microsoft WAV format (little endian).      |
| SF_FORMAT_AIFF           | 0x020000   | Apple/SGI AIFF format (big endian).        |
| SF_FORMAT_AU             | 0x030000   | Sun/NeXT AU format (big endian).           |
| SF_FORMAT_RAW            | 0x040000   | RAW PCM data.                              |
| SF_FORMAT_PAF            | 0x050000   | Ensoniq PARIS file format.                 |
| SF_FORMAT_SVX            | 0x060000   | Amiga IFF / SVX8 / SV16 format.            |
| SF_FORMAT_NIST           | 0x070000   | Sphere NIST format.                        |
| SF_FORMAT_VOC            | 0x080000   | VOC files.                                 |
| SF_FORMAT_IRCAM          | 0x0A0000   | Berkeley/IRCAM/CARL                        |
| SF_FORMAT_W64            | 0x0B0000   | Sonic Foundry's 64 bit RIFF/WAV            |
| SF_FORMAT_MAT4           | 0x0C0000   | Matlab (tm) V4.2 / GNU Octave 2.0          |
| SF_FORMAT_MAT5           | 0x0D0000   | Matlab (tm) V5.0 / GNU Octave 2.1          |
| SF_FORMAT_PVF            | 0x0E0000   | Portable Voice Format                      |
| SF_FORMAT_XI             | 0x0F0000   | Fasttracker 2 Extended Instrument          |
| SF_FORMAT_HTK            | 0x100000   | HMM Tool Kit format                        |
| SF_FORMAT_SDS            | 0x110000   | Midi Sample Dump Standard                  |
| SF_FORMAT_AVR            | 0x120000   | Audio Visual Research                      |
| SF_FORMAT_WAVEX          | 0x130000   | MS WAVE with WAVEFORMATEX                  |
| SF_FORMAT_SD2            | 0x160000   | Sound Designer 2                           |
| SF_FORMAT_FLAC           | 0x170000   | FLAC lossless file format                  |
| SF_FORMAT_CAF            | 0x180000   | Core Audio File format                     |
| SF_FORMAT_WVE            | 0x190000   | Psion WVE format                           |
| SF_FORMAT_OGG            | 0x200000   | Xiph OGG container                         |
| SF_FORMAT_MPC2K          | 0x210000   | Akai MPC 2000 sampler                      |
| SF_FORMAT_RF64           | 0x220000   | RF64 WAV file                              |
| SF_FORMAT_MPEG           | 0x230000   | MPEG-1/2 audio stream                      |
| **Subtypes.**                                                                      |
| SF_FORMAT_PCM_S8         | 0x0001     | Signed 8 bit data                          |
| SF_FORMAT_PCM_16         | 0x0002     | Signed 16 bit data                         |
| SF_FORMAT_PCM_24         | 0x0003     | Signed 24 bit data                         |
| SF_FORMAT_PCM_32         | 0x0004     | Signed 32 bit data                         |
| SF_FORMAT_PCM_U8         | 0x0005     | Unsigned 8 bit data (WAV and RAW only)     |
| SF_FORMAT_FLOAT          | 0x0006     | 32 bit float data                          |
| SF_FORMAT_DOUBLE         | 0x0007     | 64 bit float data                          |
| SF_FORMAT_ULAW           | 0x0010     | U-Law encoded.                             |
| SF_FORMAT_ALAW           | 0x0011     | A-Law encoded.                             |
| SF_FORMAT_IMA_ADPCM      | 0x0012     | IMA ADPCM.                                 |
| SF_FORMAT_MS_ADPCM       | 0x0013     | Microsoft ADPCM.                           |
| SF_FORMAT_GSM610         | 0x0020     | GSM 6.10 encoding.                         |
| SF_FORMAT_VOX_ADPCM      | 0x0021     | OKI / Dialogix ADPCM                       |
| SF_FORMAT_NMS_ADPCM_16   | 0x0022     | 16kbs NMS G721-variant encoding.           |
| SF_FORMAT_NMS_ADPCM_24   | 0x0023     | 24kbs NMS G721-variant encoding.           |
| SF_FORMAT_NMS_ADPCM_32   | 0x0024     | 32kbs NMS G721-variant encoding.           |
| SF_FORMAT_G721_32        | 0x0030     | 32kbs G721 ADPCM encoding.                 |
| SF_FORMAT_G723_24        | 0x0031     | 24kbs G723 ADPCM encoding.                 |
| SF_FORMAT_G723_40        | 0x0032     | 40kbs G723 ADPCM encoding.                 |
| SF_FORMAT_DWVW_12        | 0x0040     | 12 bit Delta Width Variable Word encoding. |
| SF_FORMAT_DWVW_16        | 0x0041     | 16 bit Delta Width Variable Word encoding. |
| SF_FORMAT_DWVW_24        | 0x0042     | 24 bit Delta Width Variable Word encoding. |
| SF_FORMAT_DWVW_N         | 0x0043     | N bit Delta Width Variable Word encoding.  |
| SF_FORMAT_DPCM_8         | 0x0050     | 8 bit differential PCM (XI only)           |
| SF_FORMAT_DPCM_16        | 0x0051     | 16 bit differential PCM (XI only)          |
| SF_FORMAT_VORBIS         | 0x0060     | Xiph Vorbis encoding.                      |
| SF_FORMAT_OPUS           | 0x0064     | Xiph/Skype Opus encoding.                  |
| SF_FORMAT_ALAC_16        | 0x0070     | Apple Lossless Audio Codec (16 bit).       |
| SF_FORMAT_ALAC_20        | 0x0071     | Apple Lossless Audio Codec (20 bit).       |
| SF_FORMAT_ALAC_24        | 0x0072     | Apple Lossless Audio Codec (24 bit).       |
| SF_FORMAT_ALAC_32        | 0x0073     | Apple Lossless Audio Codec (32 bit).       |
| SF_FORMAT_MPEG_LAYER_I   | 0x0080     | MPEG-1 Audio Layer I.                      |
| SF_FORMAT_MPEG_LAYER_II  | 0x0081     | MPEG-1 Audio Layer II.                     |
| SF_FORMAT_MPEG_LAYER_III | 0x0082     | MPEG-2 Audio Layer III.                    |
| **Endian-ness options.**                                                           |
| SF_ENDIAN_FILE           | 0x00000000 | Default file endian-ness.                  |
| SF_ENDIAN_LITTLE         | 0x10000000 | Force little endian-ness.                  |
| SF_ENDIAN_BIG            | 0x20000000 | Force big endian-ness.                     |
| SF_ENDIAN_CPU            | 0x30000000 | Force CPU endian-ness.                     |
| SF_FORMAT_SUBMASK        | 0x0000FFFF |                                            |
| SF_FORMAT_TYPEMASK       | 0x0FFF0000 |                                            |
| SF_FORMAT_ENDMASK        | 0x30000000 |                                            |

Every call to **sf_open**() should be matched with a call to
[**sf_close**()](#close) to free up memory allocated during the call to **sf_open**().

On success, the sf_open function returns a non-NULL pointer which should be passed as the first parameter to all
subsequent libsndfile calls dealing with that audio file. On fail, the sf_open function returns a NULL pointer. An
explanation of the error can obtained by passing NULL to [**sf_strerror**()](#error).

### File Descriptor Open {#open_fd}

```c
SNDFILE*  sf_open_fd (int fd, int mode, SF_INFO *sfinfo, int close_desc) ;
```

**Note:** On Microsoft Windows, this function does not work if the application
and the libsndfile DLL are linked to different versions of the Microsoft C
runtime DLL.

The second open function takes a file descriptor of a file that has already been
opened. Care should be taken to ensure that the mode of the file represented by
the descriptor matches the mode argument. This function is useful in the
following circumstances:

* Opening temporary files securely (ie use the **tmpfile**() to return a FILE*
  pointer and then using fileno() to retrieve the file descriptor which is then
  passed to libsndfile).
* Opening files with file names using OS specific character encodings and then
  passing the file descriptor to **sf_open_fd**().
* Opening sound files embedded within larger files. [More info](embedded_files.md).

Every call to `sf_open_fd`() should be matched with a call to sf_close() to free
up memory allocated during the call to sf_open_fd().

When sf_close() is called, the file descriptor is only closed if the
**close_desc** parameter was TRUE when the sf_open_fd() function was called.

On success, the sf_open_fd() function returns a non-NULL pointer which should be
passed as the first parameter to all subsequent libsndfile calls dealing with
that audio file. On fail, the sf_open_fd() function returns a NULL pointer.

### Virtual File Open Function {#open_virtual}

```c
SNDFILE*    sf_open_virtual (SF_VIRTUAL_IO *sfvirtual, int mode, SF_INFO *sfinfo, void *user_data) ;
```

Opens a soundfile from a virtual file I/O context which is provided by the
caller. This is usually used to interface libsndfile to write/read from memory
with a stream or buffer based system. Apart from the sfvirtual and the user_data
parameters this function behaves like [sf_open()](#open).

```c
    typedef struct
    {    sf_vio_get_filelen  get_filelen ;
          sf_vio_seek         seek ;
          sf_vio_read         read ;
          sf_vio_write        write ;
          sf_vio_tell         tell ;
    } SF_VIRTUAL_IO ;
```

Libsndfile calls the callbacks provided by the SF_VIRTUAL_IO structure when
opening, reading and writing to the virtual file context. The user_data pointer
is a user defined context which will be available in the callbacks.

```c
typedef sf_count_t  (*sf_vio_get_filelen) (void *user_data) ;
typedef sf_count_t  (*sf_vio_seek)        (sf_count_t offset, int whence, void *user_data) ;
typedef sf_count_t  (*sf_vio_read)        (void *ptr, sf_count_t count, void *user_data) ;
typedef sf_count_t  (*sf_vio_write)       (const void *ptr, sf_count_t count, void *user_data) ;
typedef sf_count_t  (*sf_vio_tell)        (void *user_data) ;
```

#### sf_vio_get_filelen

```c
typedef sf_count_t  (*sf_vio_get_filelen) (void *user_data) ;
```

The virtual file context must return the length of the virtual file in bytes.

#### sf_vio_seek

```c
typedef sf_count_t  (*sf_vio_seek)        (sf_count_t offset, int whence, void *user_data) ;
```

The virtual file context must seek to offset using the seek mode provided by
whence which is one of SEEK_CUR, SEEK_SET, SEEK_END.

The return value must contain the new offset in the file.

#### sf_vio_read

```c
typedef sf_count_t  (*sf_vio_read)        (void *ptr, sf_count_t count, void *user_data) ;
```

The virtual file context must copy ("read") "count" bytes into the buffer
provided by ptr and return the count of actually copied bytes.

#### sf_vio_write

```c
typedef sf_count_t  (*sf_vio_write)       (const void *ptr, sf_count_t count, void *user_data) ;
```

The virtual file context must process "count" bytes stored in the buffer passed
with ptr and return the count of actually processed bytes.

#### sf_vio_tell

```c
typedef sf_count_t  (*sf_vio_tell)        (void *user_data) ;
```

Return the current position of the virtual file context.

## Format Check Function {#chek}

```c
int  sf_format_check (const SF_INFO *info) ;
```

This function allows the caller to check if a set of parameters in the SF_INFO
struct is valid before calling [sf_open](#open) (SFM_WRITE).

sf_format_check() returns TRUE if the parameters are valid and FALSE otherwise.

## File Seek Functions

```c
sf_count_t  sf_seek  (SNDFILE *sndfile, sf_count_t frames, int whence) ;
```

The file seek functions work much like lseek in unistd.h with the exception that
the non-audio data is ignored and the seek only moves within the audio data
section of the file. In addition, seeks are defined in number of (multichannel)
frames. Therefore, a seek in a stereo file from the current position forward
with an offset of 1 would skip forward by one sample of both channels.

like lseek(), the whence parameter can be any one of the following three values:

SEEK_SET
: The offset is set to the start of the audio data plus offset (multichannel)
frames.

SEEK_CUR
: The offset is set to its current location plus offset (multichannel) frames.

SEEK_END
: The offset is set to the end of the data plus offset (multichannel) frames.

Internally, libsndfile keeps track of the read and write locations using
separate read and write pointers. If a file has been opened with a mode of
SFM_RDWR, bitwise OR-ing the standard whence values above with either SFM_READ
or SFM_WRITE allows the read and write pointers to be modified separately.
If the SEEK_* values are used on their own, the read and write pointers are
both modified.

Note that the frames offset can be negative and in fact should be when SEEK_END
is used for the whence parameter.

sf_seek will return the offset in (multichannel) frames from the start of the
audio data or -1 if an error occurred (ie an attempt is made to seek beyond the
start or end of the file).

## Error Reporting Functions {#error}

```c
int sf_error (SNDFILE *sndfile) ;
```

This function returns the current error number for the given SNDFILE.

The error number may be one of the following:

| Name                        | Value |
|:----------------------------|:------|
| SF_ERR_NO_ERROR             | 0     |
| SF_ERR_UNRECOGNISED_FORMAT  | 1     |
| SF_ERR_SYSTEM               | 2     |
| SF_ERR_MALFORMED_FILE       | 3     |
| SF_ERR_UNSUPPORTED_ENCODING | 4     |

or any one of many other internal error values.
Applications should only test the return value against error values defined in
\<sndfile.h\>; as the internal error values are subject to change at any time.
For errors not in the above list, the function sf_error_number() can be used to
convert it to an error string.

```c
const char* sf_strerror     (SNDFILE *sndfile) ;
const char* sf_error_number (int errnum) ;
```

The error functions sf_strerror () and sf_error_number () convert the library's
internal error enumerations into text strings.

```c
int sf_perror    (SNDFILE *sndfile) ;
int sf_error_str (SNDFILE *sndfile, char* str, size_t len) ;
```

The functions sf_perror() and sf_error_str() are deprecated and will be dropped
from the library at some later date.

## File Close Function {#close}

```c
int sf_close (SNDFILE *sndfile) ;
```

The close function closes the file, deallocates its internal buffers and returns
0 on success or an error value otherwise.

## Write Sync Function {#write_sync}

```c
void  sf_write_sync  (SNDFILE *sndfile) ;
```

If the file is opened SFM_WRITE or SFM_RDWR, call the operating system's
function to force the writing of all file cache buffers to disk. If the file is
opened SFM_READ no action is taken.

## File Read Functions {#read}

```c
sf_count_t sf_read_short  (SNDFILE *sndfile, short *ptr, sf_count_t items) ;
sf_count_t sf_read_int    (SNDFILE *sndfile, int *ptr, sf_count_t items) ;
sf_count_t sf_read_float  (SNDFILE *sndfile, float *ptr, sf_count_t items) ;
sf_count_t sf_read_double (SNDFILE *sndfile, double *ptr, sf_count_t items) ;
```

{: #readf}
```c
sf_count_t sf_readf_short  (SNDFILE *sndfile, short *ptr, sf_count_t frames) ;
sf_count_t sf_readf_int    (SNDFILE *sndfile, int *ptr, sf_count_t frames) ;
sf_count_t sf_readf_float  (SNDFILE *sndfile, float *ptr, sf_count_t frames) ;
sf_count_t sf_readf_double (SNDFILE *sndfile, double *ptr, sf_count_t frames) ;
```

The file read functions fill the array pointed to by ptr with the requested
number of items or frames.

For the frames-count functions, the frames parameter specifies the number of
frames. A frame is just a block of samples, one for each channel.

**Care must be taken to ensure that there is enough space in the array pointed
to by ptr, to take (frames \* channels) number of items (shorts, ints, floats or
doubles).**

For the items-count functions, the items parameter must be an integer product
of the number of channels or an error will occur. Here, an item is just a
sample.

Note: The only difference between the "items" and "frames" versions of each read
function is the units in which the object count is specified - calling
sf_readf_short() with a count argument of N, on a SNDFILE with C channels, is
the same as calling sf_read_short with a count argument of N\*C. The buffer
pointed to by "ptr" should be the same number of bytes in each case.

Note: The data type used by the calling program and the data format of the file
do not need to be the same. For instance, it is possible to open a 16 bit PCM
encoded WAV file and read the data using sf_read_float(). The library seamlessly
converts between the two formats on-the-fly. See [Note 1](#note-1).

The sf_read_XXXX and sf_readf_XXXX functions return the number of items or
frames read, respectively. Unless the end of the file was reached during the
read, the return value should equal the number of objects requested. Attempts to
read beyond the end of the file will not result in an error but will cause the
read functions to return less than the number of objects requested or 0 if
already at the end of the file. When the buffer is not is not completely filled,
unused buffer space is filled by zeroes.

## File Write Functions {#write}

```c
sf_count_t sf_write_short  (SNDFILE *sndfile, short *ptr, sf_count_t items) ;
sf_count_t sf_write_int    (SNDFILE *sndfile, int *ptr, sf_count_t items) ;
sf_count_t sf_write_float  (SNDFILE *sndfile, float *ptr, sf_count_t items) ;
sf_count_t sf_write_double (SNDFILE *sndfile, double *ptr, sf_count_t items) ;
```

{: #writef}
```c
sf_count_t sf_writef_short  (SNDFILE *sndfile, short *ptr, sf_count_t frames) ;
sf_count_t sf_writef_int    (SNDFILE *sndfile, int *ptr, sf_count_t frames) ;
sf_count_t sf_writef_float  (SNDFILE *sndfile, float *ptr, sf_count_t frames) ;
sf_count_t sf_writef_double (SNDFILE *sndfile, double *ptr, sf_count_t frames) ;
```

The file write functions write the data in the array pointed to by ptr to the
file.

For items-count functions, the items parameter specifies the size of the array
and must be an integer product of the number of channels or an error will occur.

For the frames-count functions, the array is expected to be large enough to hold
a number of items equal to the product of frames and the number of channels.

As with the read functions [above](#read), the only difference in the items and
frames version of each write function is the units in which the buffer size is
specified. Again, the data type used by the calling program and the data format
of the file do not need to be the same ([Note 1](#note-1)).

The sf_write_XXXX and sf_writef_XXXX functions respectively return the number of
items or frames written (which should be the same as the items or frames
parameter).

## Raw File Read and Write Functions {#raw}

```c
sf_count_t sf_read_raw  (SNDFILE *sndfile, void *ptr, sf_count_t bytes) ;
sf_count_t sf_write_raw (SNDFILE *sndfile, void *ptr, sf_count_t bytes) ;
```

**Note:** Unless you are writing an external decoder/encode that uses libsndfile
to handle the file headers, you should not be using these functions.

The raw read and write functions read raw audio data from the audio file (not to
be confused with reading RAW header-less PCM files). The number of bytes read or
written must always be an integer multiple of the number of channels multiplied
by the number of bytes required to represent one sample from one channel.

The raw read and write functions return the number of bytes read or written
(which should be the same as the bytes parameter).

**Note : The result of using of both regular reads/writes and raw reads/writes
on compressed file formats other than SF_FORMAT_ALAW and SF_FORMAT_ULAW is
undefined.**

See also : [SFC_RAW_NEEDS_ENDSWAP](command.md#sfc_raw_needs_endswap).

## Functions for Reading and Writing String Data {#string}

```c
const char* sf_get_string (SNDFILE *sndfile, int str_type) ;
int         sf_set_string (SNDFILE *sndfile, int str_type, const char* str) ;
```

These functions allow strings to be set on files opened for write and to be
retrieved from files opened for read where supported by the given file type. The
**str_type** parameter can be any one of the following string types:

| Name               | Value | Description   |
|:-------------------|:------|:--------------|
| SF_STR_TITLE       | 0x01  | Title.        |
| SF_STR_COPYRIGHT   | 0x02  | Copyright.    |
| SF_STR_SOFTWARE    | 0x03  | Software.     |
| SF_STR_ARTIST      | 0x04  | Artist.       |
| SF_STR_COMMENT     | 0x05  | Comment.      |
| SF_STR_DATE        | 0x06  | Date.         |
| SF_STR_ALBUM       | 0x07  | Album.        |
| SF_STR_LICENSE     | 0x08  | License.      |
| SF_STR_TRACKNUMBER | 0x09  | Track number. |
| SF_STR_GENRE       | 0x10  | Genre.        |

The sf_get_string() function returns the specified string if it exists and a
NULL pointer otherwise. In addition to the string ids above, SF_STR_FIRST (==
SF_STR_TITLE) and SF_STR_LAST (always the same as the highest numbers string id)
are also available to allow iteration over all the available string ids.

The sf_set_string() function sets the string data. It returns zero on success
and non-zero on error.The error code can be converted to a string using
sf_error_number().

Strings passed to and retrieved from these two functions are assumed to be
utf-8. However, while formats like Ogg/Vorbis and FLAC fully support utf-8,
others like WAV and AIFF officially only support ASCII. Writing utf-8 strings to
WAV and AIF files with libsndfile will work when read back with libsndfile, but
may not work with other programs.

The suggested method of dealing with tags retrieved using sf_get_string() is to
assume they are utf-8. Similarly if you have a string in some exotic format like
utf-16, it should be encoded to utf-8 before being written using libsndfile.

## Function for retrieving library version {#version_string}

```c
const char *sf_version_string (void) ;
```

Return the library version string.

## Function for retrieving current byterate {#current_byterate}

```c
int sf_current_byterate (SNDFILE *sndfile) ;
```

Return the current byterate at this point in the file. The byte rate in this
case is the number of bytes per second of audio data. For instance, for a
stereo, 18 bit PCM encoded file with an 16kHz sample rate, the byte rate
would be 2 (stereo) \* 2 (two bytes per sample) * 16000 => 64000 bytes/sec.

For some file formats the returned value will be accurate and exact, for some
it will be a close approximation, for some it will be the average bitrate for
the whole file and for some it will be a time varying value that was accurate
when the file was most recently read or written.

To get the bitrate, multiple this value by 8.

`sf_current_byterate` returns byte per second or -1 if byterate is
unknown.

## Functions to get and set chunks from within a sound file

These functions allow the getting and setting of chunks within a sound file (for
those formats which allow it).

These functions fail safely. Specifically, they will not allow you to overwrite
existing chunks or add extra versions of format specific reserved chunks but
should allow you to retrieve any and all chunks (may not be implemented for all
chunks or all file formats).

### sf_set_chunk

```c
int sf_set_chunk (SNDFILE *sndfile, const SF_CHUNK_INFO *chunk_info) ;
```

Set the specified chunk info (must be done before any audio data is written to
the file). This will fail for format specific reserved chunks. The
`chunk_info->data` pointer must be valid until the file is closed.

The `SF_CHUNK_INFO` struct is documented as follows:

```c
struct SF_CHUNK_INFO
{     char        id [64] ;   /* The chunk identifier. */
    unsigned      id_size ;   /* The size of the chunk identifier. */
    unsigned      datalen ;   /* The size of that data. */
    void          *data ;     /* Pointer to the data. */
} ;
    typedef struct SF_CHUNK_INFO SF_CHUNK_INFO ;
```

`sf_set_chunk` returns `SF_ERR_NO_ERROR` on success or non-zero on failure.

### sf_get_chunk_iterator

```c
SF_CHUNK_ITERATOR *
sf_get_chunk_iterator (SNDFILE *sndfile, const SF_CHUNK_INFO *chunk_info) ;
```

Get an iterator for all chunks matching `chunk_info`.

`SF_CHUNK_ITERATOR` is an opaque structure to an iterator over the all chunks of
a given id and defined as follows:

```c
typedef	struct SF_CHUNK_ITERATOR SF_CHUNK_ITERATOR ;
```

The iterator will point to the first chunk matching `chunk_info`. Chunks are
matching, if (`chunk_info->id`) matches the first (`chunk_info->id_size`) bytes
of a chunk found in the `SNDFILE*` handle. If `chunk_info` is `NULL`, an
iterator to all chunks in the `SNDFILE*` handle is returned. The values of
`chunk_info->datalen` and `chunk_info->data` are ignored. If no matching chunks
are found in the sndfile, `NULL` is returned.

The returned iterator will stay valid until one of the following occurs:

* The sndfile is closed.
* A new chunk is added using [`sf_set_chunk()`](#sf_set_chunk).
* Another chunk iterator function is called on the same `SNDFILE*`
  handle that causes the iterator to be modified.

The memory for the iterator belongs to the SNDFILE* handle and is freed when
[sf_close](#close) is called.

### sf_next_chunk_iterator

```c
sf_next_chunk_iterator (SF_CHUNK_ITERATOR * iterator) ;
```

Iterate through chunks by incrementing the iterator.

Increments the iterator and returns a handle to the new one. After this call,
iterator will no longer be valid, and you must use the newly returned handle
from now on. The returned handle can be used to access the next chunk matching
the criteria as defined in [sf_get_chunk_iterator](#sf_get_chunk_iterator).
If iterator points to the last chunk, this will free all resources associated
with iterator and return `NULL`. The returned iterator will stay valid until
`sf_get_next_chunk_iterator` is called again, the sndfile is closed or a new
chunk us added.

### sf_get_chunk_size

```c
int
sf_get_chunk_size (const SF_CHUNK_ITERATOR * it, SF_CHUNK_INFO * chunk_info) ;
```

Get the size of the specified chunk.

If the specified chunk exists, the size will be returned in the `datalen` field
of the `SF_CHUNK_INFO` struct. Additionally, the id of the chunk will be copied
to the `id` field of the `SF_CHUNK_INFO` struct and it's `id_size` field will be
updated accordingly.

If the chunk doesn't exist `chunk_info->datalen` will be zero, and the `id` and
`id_size` fields will be undefined.

The function will return `SF_ERR_NO_ERROR` on success or non-zero on failure.

### sf_get_chunk_data

```c
int
sf_get_chunk_data (const SF_CHUNK_ITERATOR *it, SF_CHUNK_INFO *chunk_info) ;
```

Get the specified chunk data.

If the specified chunk exists, up to `chunk_info->datalen` bytes of the chunk
data will be copied into the `chunk_info->data` buffer (allocated by the caller)
and the `chunk_info->datalen` field updated to reflect the size of the data. The
`id` and `id_size` field will be updated according to the retrieved chunk. If
the chunk doesn't exist `chunk_info->datalen` will be zero, and the `id` and
`id_size` fields will be undefined.

The function will return `SF_ERR_NO_ERROR` on success or non-zero on failure.

## Note 1

When converting between integer PCM formats of differing size (e.g. using
sf_read_int() to read a 16 bit PCM encoded WAV file) libsndfile obeys one simple
rule:

Whenever integer data is moved from one sized container to another sized
container, the most significant bit in the source container will become the most
significant bit in the destination container.

When converting between integer data and floating point data, different rules
apply. The default behaviour when reading floating point data (sf_read_float()
or sf_read_double ()) from a file with integer data is normalisation. Regardless
of whether data in the file is 8, 16, 24 or 32 bit wide, the data will be read
as floating point data in the range [-1.0, 1.0]. Similarly, data in the range
[-1.0, 1.0] will be written to an integer PCM file so that a data value of 1.0
will be the largest allowable integer for the given bit width. This
normalisation can be turned on or off using the [sf_command](command.md)
interface.

## Note 2

Reading a file containing floating point data (allowable with WAV, AIFF, AU and
other file formats) using integer read methods (sf_read_short() or
sf_read_int()) can produce unexpected results. For instance the data in the file
may have a maximum absolute value &lt; 1.0 which would mean that all sample
values read from the file will be zero. In order to read these files correctly
using integer read methods, it is recommended that you use the
[sf_command](command.md) interface, a command of
[SFC_SET_SCALE_FLOAT_INT_READ](command.md#sfc_set_scale_float_int_read) and a
parameter of SF_TRUE to force correct scaling.

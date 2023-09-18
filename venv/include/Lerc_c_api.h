/*
Copyright 2016 - 2022 Esri

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

A local copy of the license and additional notices are located with the
source distribution at:

http://github.com/Esri/lerc/

Contributors:  Thomas Maurer
*/

#ifndef LERC_API_INCLUDE_GUARD
#define LERC_API_INCLUDE_GUARD

//#define USE_EMSCRIPTEN    // to build a wasm Lerc decoder, install emscripten first

#ifdef USE_EMSCRIPTEN
  #include <emscripten/emscripten.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* LERC version numbers and related macros added in 3.0.0 */

#define LERC_VERSION_MAJOR 4
#define LERC_VERSION_MINOR 0
#define LERC_VERSION_PATCH 0

/* Macro to compute a LERC version number from its components */
#define LERC_COMPUTE_VERSION(maj,min,patch) ((maj)*10000+(min)*100+(patch))

/* Current LERC version from the above version numbers */
#define LERC_VERSION_NUMBER                 \
    LERC_COMPUTE_VERSION(LERC_VERSION_MAJOR, LERC_VERSION_MINOR, LERC_VERSION_PATCH)

/* Macro that returns true if the current LERC version is at least the version specified by (maj,min,patch) */
#define LERC_AT_LEAST_VERSION(maj,min,patch) \
    (LERC_VERSION_NUMBER >= LERC_COMPUTE_VERSION(maj,min,patch))

#if defined _WIN32 || defined __CYGWIN__
#  if defined(LERC_STATIC)
#    define LERCDLL_API
#  elif defined(LERC_EXPORTS)
#    define LERCDLL_API __declspec(dllexport)
#  else
#    define LERCDLL_API __declspec(dllimport)
#  endif
#elif __GNUC__ >= 4
  #define LERCDLL_API __attribute__((visibility("default")))
#else
  #define LERCDLL_API
#endif

  //! C-API for LERC library


  //! Added in version 4.0: 
  //!
  //! - 1) better support 3D and 4D data, allow for lossy encoding even if a noData value is used
  //! - 2) better lossless compression for float and double (pass maxZError = 0)
  //! - 3) allow to pass integer > 32 bit as double (Lerc detects it is all integer and uses that)
  //! - 4) renamed nDim to nDepth (without changing the function signatures)
  //!
  //! More on 1). In version 3.0, for 2D images, the 2D valid / invalid byte masks represent invalid pixels. 
  //! For more than 1 band, different masks per band can be used. No change to that. 
  //! For nDepth > 1, or an array of values per pixel, there is the special case of a mix of valid and invalid values 
  //! at the same pixel. The 2D mask cannot cover this case. 
  //! We have added 4 new functions to version 4.0 to cover this case, see below. If you don't encounter this
  //! "mixed case", you can continue using the same API functions as in version 3.0. 
  //! If you should encounter a Lerc blob that has this mix, both the regular lerc_decode() and 
  //! lerc_getDataRanges() functions will fail with "ErrCode::HasNoData". 
  //! In that case, you need to call the new lerc_decode_4D() function. 
  //! 
  //! More on 2). Better lossless compression for float and double is enabled for all API functions. 
  //! For 1) and 3) you have to call the new "..._4D()" functions, see further below. 


  typedef unsigned int lerc_status;

  //! All output buffers must have been allocated by the caller.

  //! Compute the buffer size in bytes required to hold the compressed input tile. Optional.
  //! You can call lerc_encode(...) directly as long as the output buffer is big enough.

  //! Order of raw input data is top left corner to lower right corner, row by row. This for each band.
  //! Data type is  { char = 0, uchar = 1, short = 2, ushort = 3, int = 4, uint = 5, float = 6, double = 7 }, see Lerc_types.h .
  //! maxZErr is the max compression error per pixel allowed.

  //! The image or mask of valid pixels is optional. Nullptr means all pixels are valid.
  //! If not all pixels are valid, set invalid pixel bytes to 0, valid pixel bytes to 1.
  //! Size of the valid / invalid pixel image is (nCols * nRows * nMasks).

  LERCDLL_API
    lerc_status lerc_computeCompressedSize(
      const void* pData,                 // raw image data, row by row, band by band
      unsigned int dataType,             // char = 0, uchar = 1, short = 2, ushort = 3, int = 4, uint = 5, float = 6, double = 7
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      int nMasks,                        // 0 - all valid, 1 - same mask for all bands, nBands - masks can differ between bands
      const unsigned char* pValidBytes,  // nullptr if all pixels are valid; otherwise 1 byte per pixel (1 = valid, 0 = invalid)
      double maxZErr,                    // max coding error per pixel, defines the precision
      unsigned int* numBytes);           // size of outgoing Lerc blob
      
  //! Encode the input data into a compressed Lerc blob.

  LERCDLL_API
    lerc_status lerc_encode(
      const void* pData,                 // raw image data, row by row, band by band
      unsigned int dataType,             // char = 0, uchar = 1, short = 2, ushort = 3, int = 4, uint = 5, float = 6, double = 7
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      int nMasks,                        // 0 - all valid, 1 - same mask for all bands, nBands - masks can differ between bands
      const unsigned char* pValidBytes,  // nullptr if all pixels are valid; otherwise 1 byte per pixel (1 = valid, 0 = invalid)
      double maxZErr,                    // max coding error per pixel, defines the precision
      unsigned char* pOutBuffer,         // buffer to write to, function fails if buffer too small
      unsigned int outBufferSize,        // size of output buffer
      unsigned int* nBytesWritten);      // number of bytes written to output buffer


  //! Use the 2 functions below to encode to an older codec version

  LERCDLL_API
    lerc_status lerc_computeCompressedSizeForVersion(
      const void* pData,                 // raw image data, row by row, band by band
      int codecVersion,                  // [2 .. 6] for [v2.2 .. v2.6], or -1 for latest codec v2.6
      unsigned int dataType,             // char = 0, uchar = 1, short = 2, ushort = 3, int = 4, uint = 5, float = 6, double = 7
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      int nMasks,                        // 0 - all valid, 1 - same mask for all bands, nBands - masks can differ between bands
      const unsigned char* pValidBytes,  // nullptr if all pixels are valid; otherwise 1 byte per pixel (1 = valid, 0 = invalid)
      double maxZErr,                    // max coding error per pixel, defines the precision
      unsigned int* numBytes);           // size of outgoing Lerc blob
     
  LERCDLL_API
    lerc_status lerc_encodeForVersion(
      const void* pData,                 // raw image data, row by row, band by band
      int codecVersion,                  // [2 .. 6] for [v2.2 .. v2.6], or -1 for latest codec v2.6
      unsigned int dataType,             // char = 0, uchar = 1, short = 2, ushort = 3, int = 4, uint = 5, float = 6, double = 7
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      int nMasks,                        // 0 - all valid, 1 - same mask for all bands, nBands - masks can differ between bands
      const unsigned char* pValidBytes,  // nullptr if all pixels are valid; otherwise 1 byte per pixel (1 = valid, 0 = invalid)
      double maxZErr,                    // max coding error per pixel, defines the precision
      unsigned char* pOutBuffer,         // buffer to write to, function fails if buffer too small
      unsigned int outBufferSize,        // size of output buffer
      unsigned int* nBytesWritten);      // number of bytes written to output buffer


  //! Call this to get info about the compressed Lerc blob. Optional.
  //! Info returned in infoArray is
  //! { version, dataType, nDepth, nCols, nRows, nBands, nValidPixels, blobSize, nMasks, nDepth, nUsesNoDataValue }, see Lerc_types.h .
  //! Info returned in dataRangeArray is { zMin, zMax, maxZErrorUsed }, see Lerc_types.h .
  //! If nDepth > 1 or nBands > 1 the data range [zMin, zMax] is over all values.

  // Remark on function signature. The arrays to be filled may grow in future versions. In order not to break
  // existing code, the function fills these arrays only up to their allocated size.

  // Remark on param blobSize. Usually it is known, either the file size of the blob written to disk,
  // or the size of the blob transmitted. It should be passed accurately for 2 reasons:
  // _ function finds out how many single band Lerc blobs are concatenated, if any
  // _ function checks for truncated file or blob
  // It is OK to pass blobSize too large as long as there is no other (valid) Lerc blob following next.
  // If in doubt, check the code in Lerc::GetLercInfo(...) for the exact logic.

  LERCDLL_API
#ifdef USE_EMSCRIPTEN
    EMSCRIPTEN_KEEPALIVE
#endif
    lerc_status lerc_getBlobInfo(
      const unsigned char* pLercBlob,    // Lerc blob to decode
      unsigned int blobSize,             // blob size in bytes
      unsigned int* infoArray,           // info array with all info needed to allocate the outgoing arrays for calling decode
      double* dataRangeArray,            // quick access to overall data range [zMin, zMax] without having to decode the data
      int infoArraySize,                 // number of elements of infoArray
      int dataRangeArraySize);           // number of elements of dataRangeArray

  //! Call this to quickly get the data ranges [min, max] per dimension and band without having to decode the pixels. Optional.
  //! The 2 output data arrays must have been allocated to the same size (nDepth * nBands).
  //! The output data array's layout is an image with nDepth columns and nBands rows.

  LERCDLL_API
#ifdef USE_EMSCRIPTEN
    EMSCRIPTEN_KEEPALIVE
#endif
    lerc_status lerc_getDataRanges(
      const unsigned char* pLercBlob,    // Lerc blob to decode
      unsigned int blobSize,             // blob size in bytes
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      double* pMins,                     // outgoing minima per dimension and band
      double* pMaxs);                    // outgoing maxima per dimension and band

  //! Decode the compressed Lerc blob into a raw data array.
  //! The data array must have been allocated to size (nDepth * nCols * nRows * nBands * sizeof(dataType)).
  //! The valid pixels array, if not all pixels valid, must have been allocated to size (nCols * nRows * nMasks).

  LERCDLL_API
#ifdef USE_EMSCRIPTEN
    EMSCRIPTEN_KEEPALIVE
#endif
    lerc_status lerc_decode(
      const unsigned char* pLercBlob,    // Lerc blob to decode
      unsigned int blobSize,             // blob size in bytes
      int nMasks,                        // 0, 1, or nBands; return as many masks in the next array
      unsigned char* pValidBytes,        // gets filled if not nullptr, even if all valid
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      unsigned int dataType,             // char = 0, uchar = 1, short = 2, ushort = 3, int = 4, uint = 5, float = 6, double = 7
      void* pData);                      // outgoing data array

  //! Same as above, but decode into double array independent of compressed data type.
  //! Wasteful in memory, but convenient if a caller from Python or C# does not want to deal with
  //! data type conversion, templating, or casting.
  //! Should this api be extended to new data types that don't fit into a double such as int64,
  //! then this function will fail for such compressed data types.

  LERCDLL_API
    lerc_status lerc_decodeToDouble(
      const unsigned char* pLercBlob,    // Lerc blob to decode
      unsigned int blobSize,             // blob size in bytes
      int nMasks,                        // 0, 1, or nBands; return as many masks in the next array
      unsigned char* pValidBytes,        // gets filled if not nullptr, even if all valid
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      double* pData);                    // outgoing data array


  //! Added in version 4.0:
  //!
  //! The 4 functions below are new. The main purpose (and difference to the functions above) is to support, for 3D and 4D data, 
  //! the special case of a mix of valid and invalid values at the same pixel. 
  //! 
  //! Main idea: Lerc has the property that for each 8x8 pixel block the minimum value is always encoded lossless in the block header. 
  //! To enable lossy encoding in the presence of noData values, the original noData value is mapped below the range of the valid values,
  //! if possible. If not possible, it switches to lossless. On decode, that temporary noData value gets mapped back to the original
  //! noData value. 
  //! 
  //! To minimize the occurence of noData values (and for better compression), Lerc tries to move noData values to the byte mask 
  //! wherever possible (e.g., all values at some pixel are invalid). So for a given band the noData values may disappear and get
  //! all moved to the byte mask. Decode only returns a noData value if it is really used. In that case the caller needs to filter
  //! the decoded arrays using both the byte mask returned and the noData value returned. 
  //! 
  //! In addition to the noData support, the new functions can also take integer values > 32 bit (but < 53 bit) as a double array, 
  //! and if all integer, use that for compression. 
  //!
  //! If floating point data contains NaN, Lerc tries to move it to the byte mask or replace it by a passed noData value. 
  //! Note, if not all NaN values can be moved to the mask (mixed case), and no noData value was passed, Lerc will fail. 
  //! It would be wrong to invent a noData value on the tile level. 


  //! Encode functions:
  //!
  //! If you don't use a noData value, are fine with the byte masks, just pass nullptr for the last 2 arguments. 
  //! 
  //! If you do have noData values at pixels that are marked as valid pixels by the byte mask, 
  //! pass 2 arrays of size nBands each, one value per band. 
  //! In pUsesNoData array, for each band, pass 1 for noData value is used, 0 if not.
  //! In noDataValues array, for each band, pass the noData value if there is one. 

  LERCDLL_API
    lerc_status lerc_computeCompressedSize_4D(
      const void* pData,                 // raw image data, row by row, band by band
      unsigned int dataType,             // char = 0, uchar = 1, short = 2, ushort = 3, int = 4, uint = 5, float = 6, double = 7
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      int nMasks,                        // 0 - all valid, 1 - same mask for all bands, nBands - masks can differ between bands
      const unsigned char* pValidBytes,  // nullptr if all pixels are valid; otherwise 1 byte per pixel (1 = valid, 0 = invalid)
      double maxZErr,                    // max coding error per pixel, defines the precision
      unsigned int* numBytes,            // size of outgoing Lerc blob
      const unsigned char* pUsesNoData,  // if there are invalid values not marked by the mask, pass an array of size nBands, 1 - uses noData, 0 - not
      const double* noDataValues);       // same, pass an array of size nBands with noData value per band, or pass nullptr

  LERCDLL_API
    lerc_status lerc_encode_4D(
      const void* pData,                 // raw image data, row by row, band by band
      unsigned int dataType,             // char = 0, uchar = 1, short = 2, ushort = 3, int = 4, uint = 5, float = 6, double = 7
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      int nMasks,                        // 0 - all valid, 1 - same mask for all bands, nBands - masks can differ between bands
      const unsigned char* pValidBytes,  // nullptr if all pixels are valid; otherwise 1 byte per pixel (1 = valid, 0 = invalid)
      double maxZErr,                    // max coding error per pixel, defines the precision
      unsigned char* pOutBuffer,         // buffer to write to, function fails if buffer too small
      unsigned int outBufferSize,        // size of output buffer
      unsigned int* nBytesWritten,       // number of bytes written to output buffer
      const unsigned char* pUsesNoData,  // if there are invalid values not marked by the mask, pass an array of size nBands, 1 - uses noData, 0 - not
      const double* noDataValues);       // same, pass an array of size nBands with noData value per band, or pass nullptr


  //! Decode functions:
  //!
  //! Same as for regular decode, first call lerc_getBlobInfo() to get all info needed from the blob header. 
  //! Check the property (InfoArray::nUsesNoDataValue) to check if there is any noData value used. 
  //!
  //! If not, just pass nullptr for the last 2 arguments. 
  //! 
  //! If yes, pass 2 arrays of size nBands each, one value per band. 
  //! In pUsesNoData array, for each band, 1 means a noData value is used, 0 means not.
  //! In noDataValues array, for each band, it has the noData value if there is one. 
  //! This is the same noData value as passed for encode. 

  LERCDLL_API
#ifdef USE_EMSCRIPTEN
    EMSCRIPTEN_KEEPALIVE
#endif
    lerc_status lerc_decode_4D(
      const unsigned char* pLercBlob,    // Lerc blob to decode
      unsigned int blobSize,             // blob size in bytes
      int nMasks,                        // 0, 1, or nBands; return as many masks in the next array
      unsigned char* pValidBytes,        // gets filled if not nullptr, even if all valid
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      unsigned int dataType,             // char = 0, uchar = 1, short = 2, ushort = 3, int = 4, uint = 5, float = 6, double = 7
      void* pData,                       // outgoing data array
      unsigned char* pUsesNoData,        // pass an array of size nBands, 1 - band uses noData, 0 - not
      double* noDataValues);             // same, pass an array of size nBands to get the noData value per band, if any

  LERCDLL_API
    lerc_status lerc_decodeToDouble_4D(
      const unsigned char* pLercBlob,    // Lerc blob to decode
      unsigned int blobSize,             // blob size in bytes
      int nMasks,                        // 0, 1, or nBands; return as many masks in the next array
      unsigned char* pValidBytes,        // gets filled if not nullptr, even if all valid
      int nDepth,                        // number of values per pixel (e.g., 3 for RGB, data is stored as [RGB, RGB, ...])
      int nCols,                         // number of columns
      int nRows,                         // number of rows
      int nBands,                        // number of bands (e.g., 3 for [RRRR ..., GGGG ..., BBBB ...])
      double* pData,                     // outgoing data array
      unsigned char* pUsesNoData,        // pass an array of size nBands, 1 - band uses noData, 0 - not
      double* noDataValues);             // same, pass an array of size nBands to get the noData value per band, if any


#ifdef __cplusplus
}
#endif

#endif  // LERC_API_INCLUDE_GUARD

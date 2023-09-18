
#pragma once

// You can include this file (if you work in C++) but you don't have to. 
// If you call this api from another language (Python, C#), you see integers. 
// This header file tells you what these integers mean. 
// These enum's may grow in the future. More values can be added. 

namespace LercNS
{
  enum class ErrCode : int
  {
    Ok = 0,
    Failed,
    WrongParam,
    BufferTooSmall,
    NaN,
    HasNoData
  };

  enum class DataType : int
  {
    dt_char = 0,
    dt_uchar,
    dt_short,
    dt_ushort,
    dt_int,
    dt_uint,
    dt_float,
    dt_double
  };

  enum class InfoArrOrder : int
  {
    version = 0,    // codec version
    dataType,
    nDim,    // = nDepth (we renamed nDim to nDepth but don't want to break anything)
    nCols,
    nRows,
    nBands,
    nValidPixels,  // for 1st band
    blobSize,
    nMasks,  // 0 - all valid, 1 - same mask for all bands, nBands - masks can differ between bands
    nDepth,  // = nDim (we renamed nDim to nDepth but don't want to break anything)
    nUsesNoDataValue,  // 0 - no noData value used, nBands - noData value used in 1 or more bands (only possible for nDepth > 1)
    _last
  };

  enum class DataRangeArrOrder : int
  {
    zMin = 0,
    zMax,
    maxZErrUsed,
    _last
  };

}    // namespace


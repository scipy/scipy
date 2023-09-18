/* Generated automatically from machmode.def and config/i386/i386-modes.def
   by genmodes.  */

#ifndef GCC_INSN_MODES_INLINE_H
#define GCC_INSN_MODES_INLINE_H

#if !defined (USED_FOR_TARGET) && GCC_VERSION >= 4001

#ifdef __cplusplus
inline __attribute__((__always_inline__))
#else
extern __inline__ __attribute__((__always_inline__, __gnu_inline__))
#endif
poly_uint16
mode_size_inline (machine_mode mode)
{
  extern poly_uint16_pod mode_size[NUM_MACHINE_MODES];
  gcc_assert (mode >= 0 && mode < NUM_MACHINE_MODES);
  switch (mode)
    {
    case E_VOIDmode: return 0;
    case E_BLKmode: return 0;
    case E_CCmode: return 4;
    case E_CCGCmode: return 4;
    case E_CCGOCmode: return 4;
    case E_CCNOmode: return 4;
    case E_CCGZmode: return 4;
    case E_CCAmode: return 4;
    case E_CCCmode: return 4;
    case E_CCOmode: return 4;
    case E_CCPmode: return 4;
    case E_CCSmode: return 4;
    case E_CCZmode: return 4;
    case E_CCFPmode: return 4;
    case E_BImode: return 1;
    case E_QImode: return 1;
    case E_HImode: return 2;
    case E_SImode: return 4;
    case E_DImode: return 8;
    case E_TImode: return 16;
    case E_OImode: return 32;
    case E_XImode: return 64;
    case E_P2QImode: return 2;
    case E_P2HImode: return 4;
    case E_POImode: return 32;
    case E_QQmode: return 1;
    case E_HQmode: return 2;
    case E_SQmode: return 4;
    case E_DQmode: return 8;
    case E_TQmode: return 16;
    case E_UQQmode: return 1;
    case E_UHQmode: return 2;
    case E_USQmode: return 4;
    case E_UDQmode: return 8;
    case E_UTQmode: return 16;
    case E_HAmode: return 2;
    case E_SAmode: return 4;
    case E_DAmode: return 8;
    case E_TAmode: return 16;
    case E_UHAmode: return 2;
    case E_USAmode: return 4;
    case E_UDAmode: return 8;
    case E_UTAmode: return 16;
    case E_HFmode: return 2;
    case E_SFmode: return 4;
    case E_DFmode: return 8;
    case E_TFmode: return 16;
    case E_SDmode: return 4;
    case E_DDmode: return 8;
    case E_TDmode: return 16;
    case E_CQImode: return 2;
    case E_CP2QImode: return 4;
    case E_CHImode: return 4;
    case E_CP2HImode: return 8;
    case E_CSImode: return 8;
    case E_CDImode: return 16;
    case E_CTImode: return 32;
    case E_CPOImode: return 64;
    case E_COImode: return 64;
    case E_CXImode: return 128;
    case E_HCmode: return 4;
    case E_SCmode: return 8;
    case E_DCmode: return 16;
    case E_TCmode: return 32;
    case E_V2QImode: return 2;
    case E_V4QImode: return 4;
    case E_V2HImode: return 4;
    case E_V1SImode: return 4;
    case E_V8QImode: return 8;
    case E_V4HImode: return 8;
    case E_V2SImode: return 8;
    case E_V1DImode: return 8;
    case E_V12QImode: return 12;
    case E_V6HImode: return 12;
    case E_V14QImode: return 14;
    case E_V16QImode: return 16;
    case E_V8HImode: return 16;
    case E_V4SImode: return 16;
    case E_V2DImode: return 16;
    case E_V1TImode: return 16;
    case E_V32QImode: return 32;
    case E_V16HImode: return 32;
    case E_V8SImode: return 32;
    case E_V4DImode: return 32;
    case E_V2TImode: return 32;
    case E_V64QImode: return 64;
    case E_V32HImode: return 64;
    case E_V16SImode: return 64;
    case E_V8DImode: return 64;
    case E_V4TImode: return 64;
    case E_V128QImode: return 128;
    case E_V64HImode: return 128;
    case E_V32SImode: return 128;
    case E_V16DImode: return 128;
    case E_V8TImode: return 128;
    case E_V64SImode: return 256;
    case E_V2HFmode: return 4;
    case E_V4HFmode: return 8;
    case E_V2SFmode: return 8;
    case E_V6HFmode: return 12;
    case E_V8HFmode: return 16;
    case E_V4SFmode: return 16;
    case E_V2DFmode: return 16;
    case E_V16HFmode: return 32;
    case E_V8SFmode: return 32;
    case E_V4DFmode: return 32;
    case E_V2TFmode: return 32;
    case E_V32HFmode: return 64;
    case E_V16SFmode: return 64;
    case E_V8DFmode: return 64;
    case E_V4TFmode: return 64;
    case E_V64HFmode: return 128;
    case E_V32SFmode: return 128;
    case E_V16DFmode: return 128;
    case E_V8TFmode: return 128;
    case E_V128HFmode: return 256;
    case E_V64SFmode: return 256;
    case E_V32DFmode: return 256;
    case E_V16TFmode: return 256;
    default: return mode_size[mode];
    }
}

#ifdef __cplusplus
inline __attribute__((__always_inline__))
#else
extern __inline__ __attribute__((__always_inline__, __gnu_inline__))
#endif
poly_uint16
mode_nunits_inline (machine_mode mode)
{
  extern const poly_uint16_pod mode_nunits[NUM_MACHINE_MODES];
  switch (mode)
    {
    case E_VOIDmode: return 0;
    case E_BLKmode: return 0;
    case E_CCmode: return 1;
    case E_CCGCmode: return 1;
    case E_CCGOCmode: return 1;
    case E_CCNOmode: return 1;
    case E_CCGZmode: return 1;
    case E_CCAmode: return 1;
    case E_CCCmode: return 1;
    case E_CCOmode: return 1;
    case E_CCPmode: return 1;
    case E_CCSmode: return 1;
    case E_CCZmode: return 1;
    case E_CCFPmode: return 1;
    case E_BImode: return 1;
    case E_QImode: return 1;
    case E_HImode: return 1;
    case E_SImode: return 1;
    case E_DImode: return 1;
    case E_TImode: return 1;
    case E_OImode: return 1;
    case E_XImode: return 1;
    case E_P2QImode: return 1;
    case E_P2HImode: return 1;
    case E_POImode: return 1;
    case E_QQmode: return 1;
    case E_HQmode: return 1;
    case E_SQmode: return 1;
    case E_DQmode: return 1;
    case E_TQmode: return 1;
    case E_UQQmode: return 1;
    case E_UHQmode: return 1;
    case E_USQmode: return 1;
    case E_UDQmode: return 1;
    case E_UTQmode: return 1;
    case E_HAmode: return 1;
    case E_SAmode: return 1;
    case E_DAmode: return 1;
    case E_TAmode: return 1;
    case E_UHAmode: return 1;
    case E_USAmode: return 1;
    case E_UDAmode: return 1;
    case E_UTAmode: return 1;
    case E_HFmode: return 1;
    case E_SFmode: return 1;
    case E_DFmode: return 1;
    case E_XFmode: return 1;
    case E_TFmode: return 1;
    case E_SDmode: return 1;
    case E_DDmode: return 1;
    case E_TDmode: return 1;
    case E_CQImode: return 2;
    case E_CP2QImode: return 2;
    case E_CHImode: return 2;
    case E_CP2HImode: return 2;
    case E_CSImode: return 2;
    case E_CDImode: return 2;
    case E_CTImode: return 2;
    case E_CPOImode: return 2;
    case E_COImode: return 2;
    case E_CXImode: return 2;
    case E_HCmode: return 2;
    case E_SCmode: return 2;
    case E_DCmode: return 2;
    case E_XCmode: return 2;
    case E_TCmode: return 2;
    case E_V2QImode: return 2;
    case E_V4QImode: return 4;
    case E_V2HImode: return 2;
    case E_V1SImode: return 1;
    case E_V8QImode: return 8;
    case E_V4HImode: return 4;
    case E_V2SImode: return 2;
    case E_V1DImode: return 1;
    case E_V12QImode: return 12;
    case E_V6HImode: return 6;
    case E_V14QImode: return 14;
    case E_V16QImode: return 16;
    case E_V8HImode: return 8;
    case E_V4SImode: return 4;
    case E_V2DImode: return 2;
    case E_V1TImode: return 1;
    case E_V32QImode: return 32;
    case E_V16HImode: return 16;
    case E_V8SImode: return 8;
    case E_V4DImode: return 4;
    case E_V2TImode: return 2;
    case E_V64QImode: return 64;
    case E_V32HImode: return 32;
    case E_V16SImode: return 16;
    case E_V8DImode: return 8;
    case E_V4TImode: return 4;
    case E_V128QImode: return 128;
    case E_V64HImode: return 64;
    case E_V32SImode: return 32;
    case E_V16DImode: return 16;
    case E_V8TImode: return 8;
    case E_V64SImode: return 64;
    case E_V2HFmode: return 2;
    case E_V4HFmode: return 4;
    case E_V2SFmode: return 2;
    case E_V6HFmode: return 6;
    case E_V8HFmode: return 8;
    case E_V4SFmode: return 4;
    case E_V2DFmode: return 2;
    case E_V16HFmode: return 16;
    case E_V8SFmode: return 8;
    case E_V4DFmode: return 4;
    case E_V2TFmode: return 2;
    case E_V32HFmode: return 32;
    case E_V16SFmode: return 16;
    case E_V8DFmode: return 8;
    case E_V4TFmode: return 4;
    case E_V64HFmode: return 64;
    case E_V32SFmode: return 32;
    case E_V16DFmode: return 16;
    case E_V8TFmode: return 8;
    case E_V128HFmode: return 128;
    case E_V64SFmode: return 64;
    case E_V32DFmode: return 32;
    case E_V16TFmode: return 16;
    default: return mode_nunits[mode];
    }
}

#ifdef __cplusplus
inline __attribute__((__always_inline__))
#else
extern __inline__ __attribute__((__always_inline__, __gnu_inline__))
#endif
unsigned char
mode_inner_inline (machine_mode mode)
{
  extern const unsigned char mode_inner[NUM_MACHINE_MODES];
  gcc_assert (mode >= 0 && mode < NUM_MACHINE_MODES);
  switch (mode)
    {
    case E_VOIDmode: return E_VOIDmode;
    case E_BLKmode: return E_BLKmode;
    case E_CCmode: return E_CCmode;
    case E_CCGCmode: return E_CCGCmode;
    case E_CCGOCmode: return E_CCGOCmode;
    case E_CCNOmode: return E_CCNOmode;
    case E_CCGZmode: return E_CCGZmode;
    case E_CCAmode: return E_CCAmode;
    case E_CCCmode: return E_CCCmode;
    case E_CCOmode: return E_CCOmode;
    case E_CCPmode: return E_CCPmode;
    case E_CCSmode: return E_CCSmode;
    case E_CCZmode: return E_CCZmode;
    case E_CCFPmode: return E_CCFPmode;
    case E_BImode: return E_BImode;
    case E_QImode: return E_QImode;
    case E_HImode: return E_HImode;
    case E_SImode: return E_SImode;
    case E_DImode: return E_DImode;
    case E_TImode: return E_TImode;
    case E_OImode: return E_OImode;
    case E_XImode: return E_XImode;
    case E_P2QImode: return E_P2QImode;
    case E_P2HImode: return E_P2HImode;
    case E_POImode: return E_POImode;
    case E_QQmode: return E_QQmode;
    case E_HQmode: return E_HQmode;
    case E_SQmode: return E_SQmode;
    case E_DQmode: return E_DQmode;
    case E_TQmode: return E_TQmode;
    case E_UQQmode: return E_UQQmode;
    case E_UHQmode: return E_UHQmode;
    case E_USQmode: return E_USQmode;
    case E_UDQmode: return E_UDQmode;
    case E_UTQmode: return E_UTQmode;
    case E_HAmode: return E_HAmode;
    case E_SAmode: return E_SAmode;
    case E_DAmode: return E_DAmode;
    case E_TAmode: return E_TAmode;
    case E_UHAmode: return E_UHAmode;
    case E_USAmode: return E_USAmode;
    case E_UDAmode: return E_UDAmode;
    case E_UTAmode: return E_UTAmode;
    case E_HFmode: return E_HFmode;
    case E_SFmode: return E_SFmode;
    case E_DFmode: return E_DFmode;
    case E_XFmode: return E_XFmode;
    case E_TFmode: return E_TFmode;
    case E_SDmode: return E_SDmode;
    case E_DDmode: return E_DDmode;
    case E_TDmode: return E_TDmode;
    case E_CQImode: return E_QImode;
    case E_CP2QImode: return E_P2QImode;
    case E_CHImode: return E_HImode;
    case E_CP2HImode: return E_P2HImode;
    case E_CSImode: return E_SImode;
    case E_CDImode: return E_DImode;
    case E_CTImode: return E_TImode;
    case E_CPOImode: return E_POImode;
    case E_COImode: return E_OImode;
    case E_CXImode: return E_XImode;
    case E_HCmode: return E_HFmode;
    case E_SCmode: return E_SFmode;
    case E_DCmode: return E_DFmode;
    case E_XCmode: return E_XFmode;
    case E_TCmode: return E_TFmode;
    case E_V2QImode: return E_QImode;
    case E_V4QImode: return E_QImode;
    case E_V2HImode: return E_HImode;
    case E_V1SImode: return E_SImode;
    case E_V8QImode: return E_QImode;
    case E_V4HImode: return E_HImode;
    case E_V2SImode: return E_SImode;
    case E_V1DImode: return E_DImode;
    case E_V12QImode: return E_QImode;
    case E_V6HImode: return E_HImode;
    case E_V14QImode: return E_QImode;
    case E_V16QImode: return E_QImode;
    case E_V8HImode: return E_HImode;
    case E_V4SImode: return E_SImode;
    case E_V2DImode: return E_DImode;
    case E_V1TImode: return E_TImode;
    case E_V32QImode: return E_QImode;
    case E_V16HImode: return E_HImode;
    case E_V8SImode: return E_SImode;
    case E_V4DImode: return E_DImode;
    case E_V2TImode: return E_TImode;
    case E_V64QImode: return E_QImode;
    case E_V32HImode: return E_HImode;
    case E_V16SImode: return E_SImode;
    case E_V8DImode: return E_DImode;
    case E_V4TImode: return E_TImode;
    case E_V128QImode: return E_QImode;
    case E_V64HImode: return E_HImode;
    case E_V32SImode: return E_SImode;
    case E_V16DImode: return E_DImode;
    case E_V8TImode: return E_TImode;
    case E_V64SImode: return E_SImode;
    case E_V2HFmode: return E_HFmode;
    case E_V4HFmode: return E_HFmode;
    case E_V2SFmode: return E_SFmode;
    case E_V6HFmode: return E_HFmode;
    case E_V8HFmode: return E_HFmode;
    case E_V4SFmode: return E_SFmode;
    case E_V2DFmode: return E_DFmode;
    case E_V16HFmode: return E_HFmode;
    case E_V8SFmode: return E_SFmode;
    case E_V4DFmode: return E_DFmode;
    case E_V2TFmode: return E_TFmode;
    case E_V32HFmode: return E_HFmode;
    case E_V16SFmode: return E_SFmode;
    case E_V8DFmode: return E_DFmode;
    case E_V4TFmode: return E_TFmode;
    case E_V64HFmode: return E_HFmode;
    case E_V32SFmode: return E_SFmode;
    case E_V16DFmode: return E_DFmode;
    case E_V8TFmode: return E_TFmode;
    case E_V128HFmode: return E_HFmode;
    case E_V64SFmode: return E_SFmode;
    case E_V32DFmode: return E_DFmode;
    case E_V16TFmode: return E_TFmode;
    default: return mode_inner[mode];
    }
}

#ifdef __cplusplus
inline __attribute__((__always_inline__))
#else
extern __inline__ __attribute__((__always_inline__, __gnu_inline__))
#endif
unsigned char
mode_unit_size_inline (machine_mode mode)
{
  extern CONST_MODE_UNIT_SIZE unsigned char mode_unit_size[NUM_MACHINE_MODES];
  gcc_assert (mode >= 0 && mode < NUM_MACHINE_MODES);
  switch (mode)
    {
    case E_VOIDmode: return 0;
    case E_BLKmode: return 0;
    case E_CCmode: return 4;
    case E_CCGCmode: return 4;
    case E_CCGOCmode: return 4;
    case E_CCNOmode: return 4;
    case E_CCGZmode: return 4;
    case E_CCAmode: return 4;
    case E_CCCmode: return 4;
    case E_CCOmode: return 4;
    case E_CCPmode: return 4;
    case E_CCSmode: return 4;
    case E_CCZmode: return 4;
    case E_CCFPmode: return 4;
    case E_BImode: return 1;
    case E_QImode: return 1;
    case E_HImode: return 2;
    case E_SImode: return 4;
    case E_DImode: return 8;
    case E_TImode: return 16;
    case E_OImode: return 32;
    case E_XImode: return 64;
    case E_P2QImode: return 2;
    case E_P2HImode: return 4;
    case E_POImode: return 32;
    case E_QQmode: return 1;
    case E_HQmode: return 2;
    case E_SQmode: return 4;
    case E_DQmode: return 8;
    case E_TQmode: return 16;
    case E_UQQmode: return 1;
    case E_UHQmode: return 2;
    case E_USQmode: return 4;
    case E_UDQmode: return 8;
    case E_UTQmode: return 16;
    case E_HAmode: return 2;
    case E_SAmode: return 4;
    case E_DAmode: return 8;
    case E_TAmode: return 16;
    case E_UHAmode: return 2;
    case E_USAmode: return 4;
    case E_UDAmode: return 8;
    case E_UTAmode: return 16;
    case E_HFmode: return 2;
    case E_SFmode: return 4;
    case E_DFmode: return 8;
    case E_TFmode: return 16;
    case E_SDmode: return 4;
    case E_DDmode: return 8;
    case E_TDmode: return 16;
    case E_CQImode: return 1;
    case E_CP2QImode: return 2;
    case E_CHImode: return 2;
    case E_CP2HImode: return 4;
    case E_CSImode: return 4;
    case E_CDImode: return 8;
    case E_CTImode: return 16;
    case E_CPOImode: return 32;
    case E_COImode: return 32;
    case E_CXImode: return 64;
    case E_HCmode: return 2;
    case E_SCmode: return 4;
    case E_DCmode: return 8;
    case E_TCmode: return 16;
    case E_V2QImode: return 1;
    case E_V4QImode: return 1;
    case E_V2HImode: return 2;
    case E_V1SImode: return 4;
    case E_V8QImode: return 1;
    case E_V4HImode: return 2;
    case E_V2SImode: return 4;
    case E_V1DImode: return 8;
    case E_V12QImode: return 1;
    case E_V6HImode: return 2;
    case E_V14QImode: return 1;
    case E_V16QImode: return 1;
    case E_V8HImode: return 2;
    case E_V4SImode: return 4;
    case E_V2DImode: return 8;
    case E_V1TImode: return 16;
    case E_V32QImode: return 1;
    case E_V16HImode: return 2;
    case E_V8SImode: return 4;
    case E_V4DImode: return 8;
    case E_V2TImode: return 16;
    case E_V64QImode: return 1;
    case E_V32HImode: return 2;
    case E_V16SImode: return 4;
    case E_V8DImode: return 8;
    case E_V4TImode: return 16;
    case E_V128QImode: return 1;
    case E_V64HImode: return 2;
    case E_V32SImode: return 4;
    case E_V16DImode: return 8;
    case E_V8TImode: return 16;
    case E_V64SImode: return 4;
    case E_V2HFmode: return 2;
    case E_V4HFmode: return 2;
    case E_V2SFmode: return 4;
    case E_V6HFmode: return 2;
    case E_V8HFmode: return 2;
    case E_V4SFmode: return 4;
    case E_V2DFmode: return 8;
    case E_V16HFmode: return 2;
    case E_V8SFmode: return 4;
    case E_V4DFmode: return 8;
    case E_V2TFmode: return 16;
    case E_V32HFmode: return 2;
    case E_V16SFmode: return 4;
    case E_V8DFmode: return 8;
    case E_V4TFmode: return 16;
    case E_V64HFmode: return 2;
    case E_V32SFmode: return 4;
    case E_V16DFmode: return 8;
    case E_V8TFmode: return 16;
    case E_V128HFmode: return 2;
    case E_V64SFmode: return 4;
    case E_V32DFmode: return 8;
    case E_V16TFmode: return 16;
    default: return mode_unit_size[mode];
    }
}

#ifdef __cplusplus
inline __attribute__((__always_inline__))
#else
extern __inline__ __attribute__((__always_inline__, __gnu_inline__))
#endif
unsigned short
mode_unit_precision_inline (machine_mode mode)
{
  extern const unsigned short mode_unit_precision[NUM_MACHINE_MODES];
  gcc_assert (mode >= 0 && mode < NUM_MACHINE_MODES);
  switch (mode)
    {
    case E_VOIDmode: return 0;
    case E_BLKmode: return 0;
    case E_CCmode: return 4*BITS_PER_UNIT;
    case E_CCGCmode: return 4*BITS_PER_UNIT;
    case E_CCGOCmode: return 4*BITS_PER_UNIT;
    case E_CCNOmode: return 4*BITS_PER_UNIT;
    case E_CCGZmode: return 4*BITS_PER_UNIT;
    case E_CCAmode: return 4*BITS_PER_UNIT;
    case E_CCCmode: return 4*BITS_PER_UNIT;
    case E_CCOmode: return 4*BITS_PER_UNIT;
    case E_CCPmode: return 4*BITS_PER_UNIT;
    case E_CCSmode: return 4*BITS_PER_UNIT;
    case E_CCZmode: return 4*BITS_PER_UNIT;
    case E_CCFPmode: return 4*BITS_PER_UNIT;
    case E_BImode: return 1;
    case E_QImode: return 1*BITS_PER_UNIT;
    case E_HImode: return 2*BITS_PER_UNIT;
    case E_SImode: return 4*BITS_PER_UNIT;
    case E_DImode: return 8*BITS_PER_UNIT;
    case E_TImode: return 16*BITS_PER_UNIT;
    case E_OImode: return 32*BITS_PER_UNIT;
    case E_XImode: return 64*BITS_PER_UNIT;
    case E_P2QImode: return 16;
    case E_P2HImode: return 32;
    case E_POImode: return 160;
    case E_QQmode: return 1*BITS_PER_UNIT;
    case E_HQmode: return 2*BITS_PER_UNIT;
    case E_SQmode: return 4*BITS_PER_UNIT;
    case E_DQmode: return 8*BITS_PER_UNIT;
    case E_TQmode: return 16*BITS_PER_UNIT;
    case E_UQQmode: return 1*BITS_PER_UNIT;
    case E_UHQmode: return 2*BITS_PER_UNIT;
    case E_USQmode: return 4*BITS_PER_UNIT;
    case E_UDQmode: return 8*BITS_PER_UNIT;
    case E_UTQmode: return 16*BITS_PER_UNIT;
    case E_HAmode: return 2*BITS_PER_UNIT;
    case E_SAmode: return 4*BITS_PER_UNIT;
    case E_DAmode: return 8*BITS_PER_UNIT;
    case E_TAmode: return 16*BITS_PER_UNIT;
    case E_UHAmode: return 2*BITS_PER_UNIT;
    case E_USAmode: return 4*BITS_PER_UNIT;
    case E_UDAmode: return 8*BITS_PER_UNIT;
    case E_UTAmode: return 16*BITS_PER_UNIT;
    case E_HFmode: return 2*BITS_PER_UNIT;
    case E_SFmode: return 4*BITS_PER_UNIT;
    case E_DFmode: return 8*BITS_PER_UNIT;
    case E_XFmode: return 80;
    case E_TFmode: return 16*BITS_PER_UNIT;
    case E_SDmode: return 4*BITS_PER_UNIT;
    case E_DDmode: return 8*BITS_PER_UNIT;
    case E_TDmode: return 16*BITS_PER_UNIT;
    case E_CQImode: return 1*BITS_PER_UNIT;
    case E_CP2QImode: return 16;
    case E_CHImode: return 2*BITS_PER_UNIT;
    case E_CP2HImode: return 32;
    case E_CSImode: return 4*BITS_PER_UNIT;
    case E_CDImode: return 8*BITS_PER_UNIT;
    case E_CTImode: return 16*BITS_PER_UNIT;
    case E_CPOImode: return 160;
    case E_COImode: return 32*BITS_PER_UNIT;
    case E_CXImode: return 64*BITS_PER_UNIT;
    case E_HCmode: return 2*BITS_PER_UNIT;
    case E_SCmode: return 4*BITS_PER_UNIT;
    case E_DCmode: return 8*BITS_PER_UNIT;
    case E_XCmode: return 80;
    case E_TCmode: return 16*BITS_PER_UNIT;
    case E_V2QImode: return 1*BITS_PER_UNIT;
    case E_V4QImode: return 1*BITS_PER_UNIT;
    case E_V2HImode: return 2*BITS_PER_UNIT;
    case E_V1SImode: return 4*BITS_PER_UNIT;
    case E_V8QImode: return 1*BITS_PER_UNIT;
    case E_V4HImode: return 2*BITS_PER_UNIT;
    case E_V2SImode: return 4*BITS_PER_UNIT;
    case E_V1DImode: return 8*BITS_PER_UNIT;
    case E_V12QImode: return 1*BITS_PER_UNIT;
    case E_V6HImode: return 2*BITS_PER_UNIT;
    case E_V14QImode: return 1*BITS_PER_UNIT;
    case E_V16QImode: return 1*BITS_PER_UNIT;
    case E_V8HImode: return 2*BITS_PER_UNIT;
    case E_V4SImode: return 4*BITS_PER_UNIT;
    case E_V2DImode: return 8*BITS_PER_UNIT;
    case E_V1TImode: return 16*BITS_PER_UNIT;
    case E_V32QImode: return 1*BITS_PER_UNIT;
    case E_V16HImode: return 2*BITS_PER_UNIT;
    case E_V8SImode: return 4*BITS_PER_UNIT;
    case E_V4DImode: return 8*BITS_PER_UNIT;
    case E_V2TImode: return 16*BITS_PER_UNIT;
    case E_V64QImode: return 1*BITS_PER_UNIT;
    case E_V32HImode: return 2*BITS_PER_UNIT;
    case E_V16SImode: return 4*BITS_PER_UNIT;
    case E_V8DImode: return 8*BITS_PER_UNIT;
    case E_V4TImode: return 16*BITS_PER_UNIT;
    case E_V128QImode: return 1*BITS_PER_UNIT;
    case E_V64HImode: return 2*BITS_PER_UNIT;
    case E_V32SImode: return 4*BITS_PER_UNIT;
    case E_V16DImode: return 8*BITS_PER_UNIT;
    case E_V8TImode: return 16*BITS_PER_UNIT;
    case E_V64SImode: return 4*BITS_PER_UNIT;
    case E_V2HFmode: return 2*BITS_PER_UNIT;
    case E_V4HFmode: return 2*BITS_PER_UNIT;
    case E_V2SFmode: return 4*BITS_PER_UNIT;
    case E_V6HFmode: return 2*BITS_PER_UNIT;
    case E_V8HFmode: return 2*BITS_PER_UNIT;
    case E_V4SFmode: return 4*BITS_PER_UNIT;
    case E_V2DFmode: return 8*BITS_PER_UNIT;
    case E_V16HFmode: return 2*BITS_PER_UNIT;
    case E_V8SFmode: return 4*BITS_PER_UNIT;
    case E_V4DFmode: return 8*BITS_PER_UNIT;
    case E_V2TFmode: return 16*BITS_PER_UNIT;
    case E_V32HFmode: return 2*BITS_PER_UNIT;
    case E_V16SFmode: return 4*BITS_PER_UNIT;
    case E_V8DFmode: return 8*BITS_PER_UNIT;
    case E_V4TFmode: return 16*BITS_PER_UNIT;
    case E_V64HFmode: return 2*BITS_PER_UNIT;
    case E_V32SFmode: return 4*BITS_PER_UNIT;
    case E_V16DFmode: return 8*BITS_PER_UNIT;
    case E_V8TFmode: return 16*BITS_PER_UNIT;
    case E_V128HFmode: return 2*BITS_PER_UNIT;
    case E_V64SFmode: return 4*BITS_PER_UNIT;
    case E_V32DFmode: return 8*BITS_PER_UNIT;
    case E_V16TFmode: return 16*BITS_PER_UNIT;
    default: return mode_unit_precision[mode];
    }
}

#endif /* GCC_VERSION >= 4001 */

#endif /* insn-modes-inline.h */

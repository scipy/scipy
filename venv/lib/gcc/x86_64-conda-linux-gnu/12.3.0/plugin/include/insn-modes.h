/* Generated automatically from machmode.def and config/i386/i386-modes.def
   by genmodes.  */

#ifndef GCC_INSN_MODES_H
#define GCC_INSN_MODES_H

enum machine_mode
{
  E_VOIDmode,              /* machmode.def:193 */
#define HAVE_VOIDmode
#ifdef USE_ENUM_MODES
#define VOIDmode E_VOIDmode
#else
#define VOIDmode ((void) 0, E_VOIDmode)
#endif
  E_BLKmode,               /* machmode.def:197 */
#define HAVE_BLKmode
#ifdef USE_ENUM_MODES
#define BLKmode E_BLKmode
#else
#define BLKmode ((void) 0, E_BLKmode)
#endif
  E_CCmode,                /* machmode.def:235 */
#define HAVE_CCmode
#ifdef USE_ENUM_MODES
#define CCmode E_CCmode
#else
#define CCmode ((void) 0, E_CCmode)
#endif
  E_CCGCmode,              /* config/i386/i386-modes.def:66 */
#define HAVE_CCGCmode
#ifdef USE_ENUM_MODES
#define CCGCmode E_CCGCmode
#else
#define CCGCmode ((void) 0, E_CCGCmode)
#endif
  E_CCGOCmode,             /* config/i386/i386-modes.def:67 */
#define HAVE_CCGOCmode
#ifdef USE_ENUM_MODES
#define CCGOCmode E_CCGOCmode
#else
#define CCGOCmode ((void) 0, E_CCGOCmode)
#endif
  E_CCNOmode,              /* config/i386/i386-modes.def:68 */
#define HAVE_CCNOmode
#ifdef USE_ENUM_MODES
#define CCNOmode E_CCNOmode
#else
#define CCNOmode ((void) 0, E_CCNOmode)
#endif
  E_CCGZmode,              /* config/i386/i386-modes.def:69 */
#define HAVE_CCGZmode
#ifdef USE_ENUM_MODES
#define CCGZmode E_CCGZmode
#else
#define CCGZmode ((void) 0, E_CCGZmode)
#endif
  E_CCAmode,               /* config/i386/i386-modes.def:70 */
#define HAVE_CCAmode
#ifdef USE_ENUM_MODES
#define CCAmode E_CCAmode
#else
#define CCAmode ((void) 0, E_CCAmode)
#endif
  E_CCCmode,               /* config/i386/i386-modes.def:71 */
#define HAVE_CCCmode
#ifdef USE_ENUM_MODES
#define CCCmode E_CCCmode
#else
#define CCCmode ((void) 0, E_CCCmode)
#endif
  E_CCOmode,               /* config/i386/i386-modes.def:72 */
#define HAVE_CCOmode
#ifdef USE_ENUM_MODES
#define CCOmode E_CCOmode
#else
#define CCOmode ((void) 0, E_CCOmode)
#endif
  E_CCPmode,               /* config/i386/i386-modes.def:73 */
#define HAVE_CCPmode
#ifdef USE_ENUM_MODES
#define CCPmode E_CCPmode
#else
#define CCPmode ((void) 0, E_CCPmode)
#endif
  E_CCSmode,               /* config/i386/i386-modes.def:74 */
#define HAVE_CCSmode
#ifdef USE_ENUM_MODES
#define CCSmode E_CCSmode
#else
#define CCSmode ((void) 0, E_CCSmode)
#endif
  E_CCZmode,               /* config/i386/i386-modes.def:75 */
#define HAVE_CCZmode
#ifdef USE_ENUM_MODES
#define CCZmode E_CCZmode
#else
#define CCZmode ((void) 0, E_CCZmode)
#endif
  E_CCFPmode,              /* config/i386/i386-modes.def:77 */
#define HAVE_CCFPmode
#ifdef USE_ENUM_MODES
#define CCFPmode E_CCFPmode
#else
#define CCFPmode ((void) 0, E_CCFPmode)
#endif
  E_BImode,                /* machmode.def:200 */
#define HAVE_BImode
#ifdef USE_ENUM_MODES
#define BImode E_BImode
#else
#define BImode (scalar_int_mode ((scalar_int_mode::from_int) E_BImode))
#endif
  E_QImode,                /* machmode.def:208 */
#define HAVE_QImode
#ifdef USE_ENUM_MODES
#define QImode E_QImode
#else
#define QImode (scalar_int_mode ((scalar_int_mode::from_int) E_QImode))
#endif
  E_HImode,                /* machmode.def:209 */
#define HAVE_HImode
#ifdef USE_ENUM_MODES
#define HImode E_HImode
#else
#define HImode (scalar_int_mode ((scalar_int_mode::from_int) E_HImode))
#endif
  E_SImode,                /* machmode.def:210 */
#define HAVE_SImode
#ifdef USE_ENUM_MODES
#define SImode E_SImode
#else
#define SImode (scalar_int_mode ((scalar_int_mode::from_int) E_SImode))
#endif
  E_DImode,                /* machmode.def:211 */
#define HAVE_DImode
#ifdef USE_ENUM_MODES
#define DImode E_DImode
#else
#define DImode (scalar_int_mode ((scalar_int_mode::from_int) E_DImode))
#endif
  E_TImode,                /* machmode.def:212 */
#define HAVE_TImode
#ifdef USE_ENUM_MODES
#define TImode E_TImode
#else
#define TImode (scalar_int_mode ((scalar_int_mode::from_int) E_TImode))
#endif
  E_OImode,                /* config/i386/i386-modes.def:104 */
#define HAVE_OImode
#ifdef USE_ENUM_MODES
#define OImode E_OImode
#else
#define OImode (scalar_int_mode ((scalar_int_mode::from_int) E_OImode))
#endif
  E_XImode,                /* config/i386/i386-modes.def:105 */
#define HAVE_XImode
#ifdef USE_ENUM_MODES
#define XImode E_XImode
#else
#define XImode (scalar_int_mode ((scalar_int_mode::from_int) E_XImode))
#endif
  E_P2QImode,              /* config/i386/i386-modes.def:110 */
#define HAVE_P2QImode
#ifdef USE_ENUM_MODES
#define P2QImode E_P2QImode
#else
#define P2QImode (scalar_int_mode ((scalar_int_mode::from_int) E_P2QImode))
#endif
  E_P2HImode,              /* config/i386/i386-modes.def:111 */
#define HAVE_P2HImode
#ifdef USE_ENUM_MODES
#define P2HImode E_P2HImode
#else
#define P2HImode (scalar_int_mode ((scalar_int_mode::from_int) E_P2HImode))
#endif
  E_POImode,               /* config/i386/i386-modes.def:116 */
#define HAVE_POImode
#ifdef USE_ENUM_MODES
#define POImode E_POImode
#else
#define POImode (scalar_int_mode ((scalar_int_mode::from_int) E_POImode))
#endif
  E_QQmode,                /* machmode.def:238 */
#define HAVE_QQmode
#ifdef USE_ENUM_MODES
#define QQmode E_QQmode
#else
#define QQmode (scalar_mode ((scalar_mode::from_int) E_QQmode))
#endif
  E_HQmode,                /* machmode.def:239 */
#define HAVE_HQmode
#ifdef USE_ENUM_MODES
#define HQmode E_HQmode
#else
#define HQmode (scalar_mode ((scalar_mode::from_int) E_HQmode))
#endif
  E_SQmode,                /* machmode.def:240 */
#define HAVE_SQmode
#ifdef USE_ENUM_MODES
#define SQmode E_SQmode
#else
#define SQmode (scalar_mode ((scalar_mode::from_int) E_SQmode))
#endif
  E_DQmode,                /* machmode.def:241 */
#define HAVE_DQmode
#ifdef USE_ENUM_MODES
#define DQmode E_DQmode
#else
#define DQmode (scalar_mode ((scalar_mode::from_int) E_DQmode))
#endif
  E_TQmode,                /* machmode.def:242 */
#define HAVE_TQmode
#ifdef USE_ENUM_MODES
#define TQmode E_TQmode
#else
#define TQmode (scalar_mode ((scalar_mode::from_int) E_TQmode))
#endif
  E_UQQmode,               /* machmode.def:244 */
#define HAVE_UQQmode
#ifdef USE_ENUM_MODES
#define UQQmode E_UQQmode
#else
#define UQQmode (scalar_mode ((scalar_mode::from_int) E_UQQmode))
#endif
  E_UHQmode,               /* machmode.def:245 */
#define HAVE_UHQmode
#ifdef USE_ENUM_MODES
#define UHQmode E_UHQmode
#else
#define UHQmode (scalar_mode ((scalar_mode::from_int) E_UHQmode))
#endif
  E_USQmode,               /* machmode.def:246 */
#define HAVE_USQmode
#ifdef USE_ENUM_MODES
#define USQmode E_USQmode
#else
#define USQmode (scalar_mode ((scalar_mode::from_int) E_USQmode))
#endif
  E_UDQmode,               /* machmode.def:247 */
#define HAVE_UDQmode
#ifdef USE_ENUM_MODES
#define UDQmode E_UDQmode
#else
#define UDQmode (scalar_mode ((scalar_mode::from_int) E_UDQmode))
#endif
  E_UTQmode,               /* machmode.def:248 */
#define HAVE_UTQmode
#ifdef USE_ENUM_MODES
#define UTQmode E_UTQmode
#else
#define UTQmode (scalar_mode ((scalar_mode::from_int) E_UTQmode))
#endif
  E_HAmode,                /* machmode.def:250 */
#define HAVE_HAmode
#ifdef USE_ENUM_MODES
#define HAmode E_HAmode
#else
#define HAmode (scalar_mode ((scalar_mode::from_int) E_HAmode))
#endif
  E_SAmode,                /* machmode.def:251 */
#define HAVE_SAmode
#ifdef USE_ENUM_MODES
#define SAmode E_SAmode
#else
#define SAmode (scalar_mode ((scalar_mode::from_int) E_SAmode))
#endif
  E_DAmode,                /* machmode.def:252 */
#define HAVE_DAmode
#ifdef USE_ENUM_MODES
#define DAmode E_DAmode
#else
#define DAmode (scalar_mode ((scalar_mode::from_int) E_DAmode))
#endif
  E_TAmode,                /* machmode.def:253 */
#define HAVE_TAmode
#ifdef USE_ENUM_MODES
#define TAmode E_TAmode
#else
#define TAmode (scalar_mode ((scalar_mode::from_int) E_TAmode))
#endif
  E_UHAmode,               /* machmode.def:255 */
#define HAVE_UHAmode
#ifdef USE_ENUM_MODES
#define UHAmode E_UHAmode
#else
#define UHAmode (scalar_mode ((scalar_mode::from_int) E_UHAmode))
#endif
  E_USAmode,               /* machmode.def:256 */
#define HAVE_USAmode
#ifdef USE_ENUM_MODES
#define USAmode E_USAmode
#else
#define USAmode (scalar_mode ((scalar_mode::from_int) E_USAmode))
#endif
  E_UDAmode,               /* machmode.def:257 */
#define HAVE_UDAmode
#ifdef USE_ENUM_MODES
#define UDAmode E_UDAmode
#else
#define UDAmode (scalar_mode ((scalar_mode::from_int) E_UDAmode))
#endif
  E_UTAmode,               /* machmode.def:258 */
#define HAVE_UTAmode
#ifdef USE_ENUM_MODES
#define UTAmode E_UTAmode
#else
#define UTAmode (scalar_mode ((scalar_mode::from_int) E_UTAmode))
#endif
  E_HFmode,                /* config/i386/i386-modes.def:26 */
#define HAVE_HFmode
#ifdef USE_ENUM_MODES
#define HFmode E_HFmode
#else
#define HFmode (scalar_float_mode ((scalar_float_mode::from_int) E_HFmode))
#endif
  E_SFmode,                /* machmode.def:230 */
#define HAVE_SFmode
#ifdef USE_ENUM_MODES
#define SFmode E_SFmode
#else
#define SFmode (scalar_float_mode ((scalar_float_mode::from_int) E_SFmode))
#endif
  E_DFmode,                /* machmode.def:231 */
#define HAVE_DFmode
#ifdef USE_ENUM_MODES
#define DFmode E_DFmode
#else
#define DFmode (scalar_float_mode ((scalar_float_mode::from_int) E_DFmode))
#endif
  E_XFmode,                /* config/i386/i386-modes.def:24 */
#define HAVE_XFmode
#ifdef USE_ENUM_MODES
#define XFmode E_XFmode
#else
#define XFmode (scalar_float_mode ((scalar_float_mode::from_int) E_XFmode))
#endif
  E_TFmode,                /* config/i386/i386-modes.def:25 */
#define HAVE_TFmode
#ifdef USE_ENUM_MODES
#define TFmode E_TFmode
#else
#define TFmode (scalar_float_mode ((scalar_float_mode::from_int) E_TFmode))
#endif
  E_SDmode,                /* machmode.def:271 */
#define HAVE_SDmode
#ifdef USE_ENUM_MODES
#define SDmode E_SDmode
#else
#define SDmode (scalar_float_mode ((scalar_float_mode::from_int) E_SDmode))
#endif
  E_DDmode,                /* machmode.def:272 */
#define HAVE_DDmode
#ifdef USE_ENUM_MODES
#define DDmode E_DDmode
#else
#define DDmode (scalar_float_mode ((scalar_float_mode::from_int) E_DDmode))
#endif
  E_TDmode,                /* machmode.def:273 */
#define HAVE_TDmode
#ifdef USE_ENUM_MODES
#define TDmode E_TDmode
#else
#define TDmode (scalar_float_mode ((scalar_float_mode::from_int) E_TDmode))
#endif
  E_CQImode,               /* machmode.def:266 */
#define HAVE_CQImode
#ifdef USE_ENUM_MODES
#define CQImode E_CQImode
#else
#define CQImode (complex_mode ((complex_mode::from_int) E_CQImode))
#endif
  E_CP2QImode,             /* machmode.def:267 */
#define HAVE_CP2QImode
#ifdef USE_ENUM_MODES
#define CP2QImode E_CP2QImode
#else
#define CP2QImode (complex_mode ((complex_mode::from_int) E_CP2QImode))
#endif
  E_CHImode,               /* machmode.def:266 */
#define HAVE_CHImode
#ifdef USE_ENUM_MODES
#define CHImode E_CHImode
#else
#define CHImode (complex_mode ((complex_mode::from_int) E_CHImode))
#endif
  E_CP2HImode,             /* machmode.def:267 */
#define HAVE_CP2HImode
#ifdef USE_ENUM_MODES
#define CP2HImode E_CP2HImode
#else
#define CP2HImode (complex_mode ((complex_mode::from_int) E_CP2HImode))
#endif
  E_CSImode,               /* machmode.def:266 */
#define HAVE_CSImode
#ifdef USE_ENUM_MODES
#define CSImode E_CSImode
#else
#define CSImode (complex_mode ((complex_mode::from_int) E_CSImode))
#endif
  E_CDImode,               /* machmode.def:266 */
#define HAVE_CDImode
#ifdef USE_ENUM_MODES
#define CDImode E_CDImode
#else
#define CDImode (complex_mode ((complex_mode::from_int) E_CDImode))
#endif
  E_CTImode,               /* machmode.def:266 */
#define HAVE_CTImode
#ifdef USE_ENUM_MODES
#define CTImode E_CTImode
#else
#define CTImode (complex_mode ((complex_mode::from_int) E_CTImode))
#endif
  E_CPOImode,              /* machmode.def:267 */
#define HAVE_CPOImode
#ifdef USE_ENUM_MODES
#define CPOImode E_CPOImode
#else
#define CPOImode (complex_mode ((complex_mode::from_int) E_CPOImode))
#endif
  E_COImode,               /* machmode.def:266 */
#define HAVE_COImode
#ifdef USE_ENUM_MODES
#define COImode E_COImode
#else
#define COImode (complex_mode ((complex_mode::from_int) E_COImode))
#endif
  E_CXImode,               /* machmode.def:266 */
#define HAVE_CXImode
#ifdef USE_ENUM_MODES
#define CXImode E_CXImode
#else
#define CXImode (complex_mode ((complex_mode::from_int) E_CXImode))
#endif
  E_HCmode,                /* machmode.def:268 */
#define HAVE_HCmode
#ifdef USE_ENUM_MODES
#define HCmode E_HCmode
#else
#define HCmode (complex_mode ((complex_mode::from_int) E_HCmode))
#endif
  E_SCmode,                /* machmode.def:268 */
#define HAVE_SCmode
#ifdef USE_ENUM_MODES
#define SCmode E_SCmode
#else
#define SCmode (complex_mode ((complex_mode::from_int) E_SCmode))
#endif
  E_DCmode,                /* machmode.def:268 */
#define HAVE_DCmode
#ifdef USE_ENUM_MODES
#define DCmode E_DCmode
#else
#define DCmode (complex_mode ((complex_mode::from_int) E_DCmode))
#endif
  E_XCmode,                /* machmode.def:268 */
#define HAVE_XCmode
#ifdef USE_ENUM_MODES
#define XCmode E_XCmode
#else
#define XCmode (complex_mode ((complex_mode::from_int) E_XCmode))
#endif
  E_TCmode,                /* machmode.def:268 */
#define HAVE_TCmode
#ifdef USE_ENUM_MODES
#define TCmode E_TCmode
#else
#define TCmode (complex_mode ((complex_mode::from_int) E_TCmode))
#endif
  E_V2QImode,              /* config/i386/i386-modes.def:98 */
#define HAVE_V2QImode
#ifdef USE_ENUM_MODES
#define V2QImode E_V2QImode
#else
#define V2QImode ((void) 0, E_V2QImode)
#endif
  E_V4QImode,              /* config/i386/i386-modes.def:81 */
#define HAVE_V4QImode
#ifdef USE_ENUM_MODES
#define V4QImode E_V4QImode
#else
#define V4QImode ((void) 0, E_V4QImode)
#endif
  E_V2HImode,              /* config/i386/i386-modes.def:81 */
#define HAVE_V2HImode
#ifdef USE_ENUM_MODES
#define V2HImode E_V2HImode
#else
#define V2HImode ((void) 0, E_V2HImode)
#endif
  E_V1SImode,              /* config/i386/i386-modes.def:97 */
#define HAVE_V1SImode
#ifdef USE_ENUM_MODES
#define V1SImode E_V1SImode
#else
#define V1SImode ((void) 0, E_V1SImode)
#endif
  E_V8QImode,              /* config/i386/i386-modes.def:82 */
#define HAVE_V8QImode
#ifdef USE_ENUM_MODES
#define V8QImode E_V8QImode
#else
#define V8QImode ((void) 0, E_V8QImode)
#endif
  E_V4HImode,              /* config/i386/i386-modes.def:82 */
#define HAVE_V4HImode
#ifdef USE_ENUM_MODES
#define V4HImode E_V4HImode
#else
#define V4HImode ((void) 0, E_V4HImode)
#endif
  E_V2SImode,              /* config/i386/i386-modes.def:82 */
#define HAVE_V2SImode
#ifdef USE_ENUM_MODES
#define V2SImode E_V2SImode
#else
#define V2SImode ((void) 0, E_V2SImode)
#endif
  E_V1DImode,              /* config/i386/i386-modes.def:96 */
#define HAVE_V1DImode
#ifdef USE_ENUM_MODES
#define V1DImode E_V1DImode
#else
#define V1DImode ((void) 0, E_V1DImode)
#endif
  E_V12QImode,             /* config/i386/i386-modes.def:99 */
#define HAVE_V12QImode
#ifdef USE_ENUM_MODES
#define V12QImode E_V12QImode
#else
#define V12QImode ((void) 0, E_V12QImode)
#endif
  E_V6HImode,              /* config/i386/i386-modes.def:101 */
#define HAVE_V6HImode
#ifdef USE_ENUM_MODES
#define V6HImode E_V6HImode
#else
#define V6HImode ((void) 0, E_V6HImode)
#endif
  E_V14QImode,             /* config/i386/i386-modes.def:100 */
#define HAVE_V14QImode
#ifdef USE_ENUM_MODES
#define V14QImode E_V14QImode
#else
#define V14QImode ((void) 0, E_V14QImode)
#endif
  E_V16QImode,             /* config/i386/i386-modes.def:83 */
#define HAVE_V16QImode
#ifdef USE_ENUM_MODES
#define V16QImode E_V16QImode
#else
#define V16QImode ((void) 0, E_V16QImode)
#endif
  E_V8HImode,              /* config/i386/i386-modes.def:83 */
#define HAVE_V8HImode
#ifdef USE_ENUM_MODES
#define V8HImode E_V8HImode
#else
#define V8HImode ((void) 0, E_V8HImode)
#endif
  E_V4SImode,              /* config/i386/i386-modes.def:83 */
#define HAVE_V4SImode
#ifdef USE_ENUM_MODES
#define V4SImode E_V4SImode
#else
#define V4SImode ((void) 0, E_V4SImode)
#endif
  E_V2DImode,              /* config/i386/i386-modes.def:83 */
#define HAVE_V2DImode
#ifdef USE_ENUM_MODES
#define V2DImode E_V2DImode
#else
#define V2DImode ((void) 0, E_V2DImode)
#endif
  E_V1TImode,              /* config/i386/i386-modes.def:95 */
#define HAVE_V1TImode
#ifdef USE_ENUM_MODES
#define V1TImode E_V1TImode
#else
#define V1TImode ((void) 0, E_V1TImode)
#endif
  E_V32QImode,             /* config/i386/i386-modes.def:84 */
#define HAVE_V32QImode
#ifdef USE_ENUM_MODES
#define V32QImode E_V32QImode
#else
#define V32QImode ((void) 0, E_V32QImode)
#endif
  E_V16HImode,             /* config/i386/i386-modes.def:84 */
#define HAVE_V16HImode
#ifdef USE_ENUM_MODES
#define V16HImode E_V16HImode
#else
#define V16HImode ((void) 0, E_V16HImode)
#endif
  E_V8SImode,              /* config/i386/i386-modes.def:84 */
#define HAVE_V8SImode
#ifdef USE_ENUM_MODES
#define V8SImode E_V8SImode
#else
#define V8SImode ((void) 0, E_V8SImode)
#endif
  E_V4DImode,              /* config/i386/i386-modes.def:84 */
#define HAVE_V4DImode
#ifdef USE_ENUM_MODES
#define V4DImode E_V4DImode
#else
#define V4DImode ((void) 0, E_V4DImode)
#endif
  E_V2TImode,              /* config/i386/i386-modes.def:84 */
#define HAVE_V2TImode
#ifdef USE_ENUM_MODES
#define V2TImode E_V2TImode
#else
#define V2TImode ((void) 0, E_V2TImode)
#endif
  E_V64QImode,             /* config/i386/i386-modes.def:85 */
#define HAVE_V64QImode
#ifdef USE_ENUM_MODES
#define V64QImode E_V64QImode
#else
#define V64QImode ((void) 0, E_V64QImode)
#endif
  E_V32HImode,             /* config/i386/i386-modes.def:85 */
#define HAVE_V32HImode
#ifdef USE_ENUM_MODES
#define V32HImode E_V32HImode
#else
#define V32HImode ((void) 0, E_V32HImode)
#endif
  E_V16SImode,             /* config/i386/i386-modes.def:85 */
#define HAVE_V16SImode
#ifdef USE_ENUM_MODES
#define V16SImode E_V16SImode
#else
#define V16SImode ((void) 0, E_V16SImode)
#endif
  E_V8DImode,              /* config/i386/i386-modes.def:85 */
#define HAVE_V8DImode
#ifdef USE_ENUM_MODES
#define V8DImode E_V8DImode
#else
#define V8DImode ((void) 0, E_V8DImode)
#endif
  E_V4TImode,              /* config/i386/i386-modes.def:85 */
#define HAVE_V4TImode
#ifdef USE_ENUM_MODES
#define V4TImode E_V4TImode
#else
#define V4TImode ((void) 0, E_V4TImode)
#endif
  E_V128QImode,            /* config/i386/i386-modes.def:86 */
#define HAVE_V128QImode
#ifdef USE_ENUM_MODES
#define V128QImode E_V128QImode
#else
#define V128QImode ((void) 0, E_V128QImode)
#endif
  E_V64HImode,             /* config/i386/i386-modes.def:86 */
#define HAVE_V64HImode
#ifdef USE_ENUM_MODES
#define V64HImode E_V64HImode
#else
#define V64HImode ((void) 0, E_V64HImode)
#endif
  E_V32SImode,             /* config/i386/i386-modes.def:86 */
#define HAVE_V32SImode
#ifdef USE_ENUM_MODES
#define V32SImode E_V32SImode
#else
#define V32SImode ((void) 0, E_V32SImode)
#endif
  E_V16DImode,             /* config/i386/i386-modes.def:86 */
#define HAVE_V16DImode
#ifdef USE_ENUM_MODES
#define V16DImode E_V16DImode
#else
#define V16DImode ((void) 0, E_V16DImode)
#endif
  E_V8TImode,              /* config/i386/i386-modes.def:86 */
#define HAVE_V8TImode
#ifdef USE_ENUM_MODES
#define V8TImode E_V8TImode
#else
#define V8TImode ((void) 0, E_V8TImode)
#endif
  E_V64SImode,             /* config/i386/i386-modes.def:102 */
#define HAVE_V64SImode
#ifdef USE_ENUM_MODES
#define V64SImode E_V64SImode
#else
#define V64SImode ((void) 0, E_V64SImode)
#endif
  E_V2HFmode,              /* config/i386/i386-modes.def:93 */
#define HAVE_V2HFmode
#ifdef USE_ENUM_MODES
#define V2HFmode E_V2HFmode
#else
#define V2HFmode ((void) 0, E_V2HFmode)
#endif
  E_V4HFmode,              /* config/i386/i386-modes.def:87 */
#define HAVE_V4HFmode
#ifdef USE_ENUM_MODES
#define V4HFmode E_V4HFmode
#else
#define V4HFmode ((void) 0, E_V4HFmode)
#endif
  E_V2SFmode,              /* config/i386/i386-modes.def:87 */
#define HAVE_V2SFmode
#ifdef USE_ENUM_MODES
#define V2SFmode E_V2SFmode
#else
#define V2SFmode ((void) 0, E_V2SFmode)
#endif
  E_V6HFmode,              /* config/i386/i386-modes.def:94 */
#define HAVE_V6HFmode
#ifdef USE_ENUM_MODES
#define V6HFmode E_V6HFmode
#else
#define V6HFmode ((void) 0, E_V6HFmode)
#endif
  E_V8HFmode,              /* config/i386/i386-modes.def:88 */
#define HAVE_V8HFmode
#ifdef USE_ENUM_MODES
#define V8HFmode E_V8HFmode
#else
#define V8HFmode ((void) 0, E_V8HFmode)
#endif
  E_V4SFmode,              /* config/i386/i386-modes.def:88 */
#define HAVE_V4SFmode
#ifdef USE_ENUM_MODES
#define V4SFmode E_V4SFmode
#else
#define V4SFmode ((void) 0, E_V4SFmode)
#endif
  E_V2DFmode,              /* config/i386/i386-modes.def:88 */
#define HAVE_V2DFmode
#ifdef USE_ENUM_MODES
#define V2DFmode E_V2DFmode
#else
#define V2DFmode ((void) 0, E_V2DFmode)
#endif
  E_V16HFmode,             /* config/i386/i386-modes.def:89 */
#define HAVE_V16HFmode
#ifdef USE_ENUM_MODES
#define V16HFmode E_V16HFmode
#else
#define V16HFmode ((void) 0, E_V16HFmode)
#endif
  E_V8SFmode,              /* config/i386/i386-modes.def:89 */
#define HAVE_V8SFmode
#ifdef USE_ENUM_MODES
#define V8SFmode E_V8SFmode
#else
#define V8SFmode ((void) 0, E_V8SFmode)
#endif
  E_V4DFmode,              /* config/i386/i386-modes.def:89 */
#define HAVE_V4DFmode
#ifdef USE_ENUM_MODES
#define V4DFmode E_V4DFmode
#else
#define V4DFmode ((void) 0, E_V4DFmode)
#endif
  E_V2TFmode,              /* config/i386/i386-modes.def:89 */
#define HAVE_V2TFmode
#ifdef USE_ENUM_MODES
#define V2TFmode E_V2TFmode
#else
#define V2TFmode ((void) 0, E_V2TFmode)
#endif
  E_V32HFmode,             /* config/i386/i386-modes.def:90 */
#define HAVE_V32HFmode
#ifdef USE_ENUM_MODES
#define V32HFmode E_V32HFmode
#else
#define V32HFmode ((void) 0, E_V32HFmode)
#endif
  E_V16SFmode,             /* config/i386/i386-modes.def:90 */
#define HAVE_V16SFmode
#ifdef USE_ENUM_MODES
#define V16SFmode E_V16SFmode
#else
#define V16SFmode ((void) 0, E_V16SFmode)
#endif
  E_V8DFmode,              /* config/i386/i386-modes.def:90 */
#define HAVE_V8DFmode
#ifdef USE_ENUM_MODES
#define V8DFmode E_V8DFmode
#else
#define V8DFmode ((void) 0, E_V8DFmode)
#endif
  E_V4TFmode,              /* config/i386/i386-modes.def:90 */
#define HAVE_V4TFmode
#ifdef USE_ENUM_MODES
#define V4TFmode E_V4TFmode
#else
#define V4TFmode ((void) 0, E_V4TFmode)
#endif
  E_V64HFmode,             /* config/i386/i386-modes.def:91 */
#define HAVE_V64HFmode
#ifdef USE_ENUM_MODES
#define V64HFmode E_V64HFmode
#else
#define V64HFmode ((void) 0, E_V64HFmode)
#endif
  E_V32SFmode,             /* config/i386/i386-modes.def:91 */
#define HAVE_V32SFmode
#ifdef USE_ENUM_MODES
#define V32SFmode E_V32SFmode
#else
#define V32SFmode ((void) 0, E_V32SFmode)
#endif
  E_V16DFmode,             /* config/i386/i386-modes.def:91 */
#define HAVE_V16DFmode
#ifdef USE_ENUM_MODES
#define V16DFmode E_V16DFmode
#else
#define V16DFmode ((void) 0, E_V16DFmode)
#endif
  E_V8TFmode,              /* config/i386/i386-modes.def:91 */
#define HAVE_V8TFmode
#ifdef USE_ENUM_MODES
#define V8TFmode E_V8TFmode
#else
#define V8TFmode ((void) 0, E_V8TFmode)
#endif
  E_V128HFmode,            /* config/i386/i386-modes.def:92 */
#define HAVE_V128HFmode
#ifdef USE_ENUM_MODES
#define V128HFmode E_V128HFmode
#else
#define V128HFmode ((void) 0, E_V128HFmode)
#endif
  E_V64SFmode,             /* config/i386/i386-modes.def:92 */
#define HAVE_V64SFmode
#ifdef USE_ENUM_MODES
#define V64SFmode E_V64SFmode
#else
#define V64SFmode ((void) 0, E_V64SFmode)
#endif
  E_V32DFmode,             /* config/i386/i386-modes.def:92 */
#define HAVE_V32DFmode
#ifdef USE_ENUM_MODES
#define V32DFmode E_V32DFmode
#else
#define V32DFmode ((void) 0, E_V32DFmode)
#endif
  E_V16TFmode,             /* config/i386/i386-modes.def:92 */
#define HAVE_V16TFmode
#ifdef USE_ENUM_MODES
#define V16TFmode E_V16TFmode
#else
#define V16TFmode ((void) 0, E_V16TFmode)
#endif
  MAX_MACHINE_MODE,

  MIN_MODE_RANDOM = E_VOIDmode,
  MAX_MODE_RANDOM = E_BLKmode,

  MIN_MODE_CC = E_CCmode,
  MAX_MODE_CC = E_CCFPmode,

  MIN_MODE_BOOL = E_BImode,
  MAX_MODE_BOOL = E_BImode,

  MIN_MODE_INT = E_QImode,
  MAX_MODE_INT = E_XImode,

  MIN_MODE_PARTIAL_INT = E_P2QImode,
  MAX_MODE_PARTIAL_INT = E_POImode,

  MIN_MODE_FRACT = E_QQmode,
  MAX_MODE_FRACT = E_TQmode,

  MIN_MODE_UFRACT = E_UQQmode,
  MAX_MODE_UFRACT = E_UTQmode,

  MIN_MODE_ACCUM = E_HAmode,
  MAX_MODE_ACCUM = E_TAmode,

  MIN_MODE_UACCUM = E_UHAmode,
  MAX_MODE_UACCUM = E_UTAmode,

  MIN_MODE_FLOAT = E_HFmode,
  MAX_MODE_FLOAT = E_TFmode,

  MIN_MODE_DECIMAL_FLOAT = E_SDmode,
  MAX_MODE_DECIMAL_FLOAT = E_TDmode,

  MIN_MODE_COMPLEX_INT = E_CQImode,
  MAX_MODE_COMPLEX_INT = E_CXImode,

  MIN_MODE_COMPLEX_FLOAT = E_HCmode,
  MAX_MODE_COMPLEX_FLOAT = E_TCmode,

  MIN_MODE_VECTOR_BOOL = E_VOIDmode,
  MAX_MODE_VECTOR_BOOL = E_VOIDmode,

  MIN_MODE_VECTOR_INT = E_V2QImode,
  MAX_MODE_VECTOR_INT = E_V64SImode,

  MIN_MODE_VECTOR_FRACT = E_VOIDmode,
  MAX_MODE_VECTOR_FRACT = E_VOIDmode,

  MIN_MODE_VECTOR_UFRACT = E_VOIDmode,
  MAX_MODE_VECTOR_UFRACT = E_VOIDmode,

  MIN_MODE_VECTOR_ACCUM = E_VOIDmode,
  MAX_MODE_VECTOR_ACCUM = E_VOIDmode,

  MIN_MODE_VECTOR_UACCUM = E_VOIDmode,
  MAX_MODE_VECTOR_UACCUM = E_VOIDmode,

  MIN_MODE_VECTOR_FLOAT = E_V2HFmode,
  MAX_MODE_VECTOR_FLOAT = E_V16TFmode,

  MIN_MODE_OPAQUE = E_VOIDmode,
  MAX_MODE_OPAQUE = E_VOIDmode,

  NUM_MACHINE_MODES = MAX_MACHINE_MODE
};

#define NUM_MODE_RANDOM (MAX_MODE_RANDOM - MIN_MODE_RANDOM + 1)
#define NUM_MODE_CC (MAX_MODE_CC - MIN_MODE_CC + 1)
#define NUM_MODE_INT (MAX_MODE_INT - MIN_MODE_INT + 1)
#define NUM_MODE_PARTIAL_INT (MAX_MODE_PARTIAL_INT - MIN_MODE_PARTIAL_INT + 1)
#define NUM_MODE_FRACT (MAX_MODE_FRACT - MIN_MODE_FRACT + 1)
#define NUM_MODE_UFRACT (MAX_MODE_UFRACT - MIN_MODE_UFRACT + 1)
#define NUM_MODE_ACCUM (MAX_MODE_ACCUM - MIN_MODE_ACCUM + 1)
#define NUM_MODE_UACCUM (MAX_MODE_UACCUM - MIN_MODE_UACCUM + 1)
#define NUM_MODE_FLOAT (MAX_MODE_FLOAT - MIN_MODE_FLOAT + 1)
#define NUM_MODE_DECIMAL_FLOAT (MAX_MODE_DECIMAL_FLOAT - MIN_MODE_DECIMAL_FLOAT + 1)
#define NUM_MODE_COMPLEX_INT (MAX_MODE_COMPLEX_INT - MIN_MODE_COMPLEX_INT + 1)
#define NUM_MODE_COMPLEX_FLOAT (MAX_MODE_COMPLEX_FLOAT - MIN_MODE_COMPLEX_FLOAT + 1)
#define NUM_MODE_VECTOR_BOOL 0
#define NUM_MODE_VECTOR_INT (MAX_MODE_VECTOR_INT - MIN_MODE_VECTOR_INT + 1)
#define NUM_MODE_VECTOR_FRACT 0
#define NUM_MODE_VECTOR_UFRACT 0
#define NUM_MODE_VECTOR_ACCUM 0
#define NUM_MODE_VECTOR_UACCUM 0
#define NUM_MODE_VECTOR_FLOAT (MAX_MODE_VECTOR_FLOAT - MIN_MODE_VECTOR_FLOAT + 1)
#define NUM_MODE_OPAQUE 0

#define CONST_MODE_NUNITS const
#define CONST_MODE_PRECISION const
#define CONST_MODE_SIZE
#define CONST_MODE_UNIT_SIZE
#define CONST_MODE_BASE_ALIGN
#define CONST_MODE_IBIT const
#define CONST_MODE_FBIT const
#define CONST_MODE_MASK const

#define BITS_PER_UNIT (8)
#define MAX_BITSIZE_MODE_ANY_INT (64*BITS_PER_UNIT)
#define MAX_BITSIZE_MODE_ANY_MODE (256*BITS_PER_UNIT)
#define NUM_INT_N_ENTS 1
#define NUM_POLY_INT_COEFFS 1

#endif /* insn-modes.h */

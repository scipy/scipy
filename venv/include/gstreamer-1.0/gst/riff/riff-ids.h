/* GStreamer RIFF I/O
 * Copyright (C) 2003 Ronald Bultje <rbultje@ronald.bitfreak.net>
 *
 * riff-ids.h: RIFF IDs and structs
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __GST_RIFF_IDS_H__
#define __GST_RIFF_IDS_H__

#include <gst/gst.h>
#include "riff-prelude.h"

G_BEGIN_DECLS

/* RIFF types */
#define GST_RIFF_RIFF_WAVE GST_MAKE_FOURCC ('W','A','V','E')
#define GST_RIFF_RIFF_AVI  GST_MAKE_FOURCC ('A','V','I',' ')
#define GST_RIFF_RIFF_CDXA GST_MAKE_FOURCC ('C','D','X','A')

/* tags */
#define GST_RIFF_TAG_RIFF GST_MAKE_FOURCC ('R','I','F','F')
#define GST_RIFF_TAG_AVF0 GST_MAKE_FOURCC ('A','V','F','0')
#define GST_RIFF_TAG_RF64 GST_MAKE_FOURCC ('R','F','6','4')
#define GST_RIFF_TAG_RIFX GST_MAKE_FOURCC ('R','I','F','X')
#define GST_RIFF_TAG_LIST GST_MAKE_FOURCC ('L','I','S','T')
#define GST_RIFF_TAG_avih GST_MAKE_FOURCC ('a','v','i','h')
#define GST_RIFF_TAG_strd GST_MAKE_FOURCC ('s','t','r','d')
#define GST_RIFF_TAG_strn GST_MAKE_FOURCC ('s','t','r','n')
#define GST_RIFF_TAG_strh GST_MAKE_FOURCC ('s','t','r','h')
#define GST_RIFF_TAG_strf GST_MAKE_FOURCC ('s','t','r','f')
#define GST_RIFF_TAG_vedt GST_MAKE_FOURCC ('v','e','d','t')
#define GST_RIFF_TAG_JUNK GST_MAKE_FOURCC ('J','U','N','K')
#define GST_RIFF_TAG_JUNQ GST_MAKE_FOURCC ('J','U','N','Q')
#define GST_RIFF_TAG_idx1 GST_MAKE_FOURCC ('i','d','x','1')
#define GST_RIFF_TAG_dmlh GST_MAKE_FOURCC ('d','m','l','h')
#define GST_RIFF_TAG_ID32 GST_MAKE_FOURCC ('I','D','3','2')
#define GST_RIFF_TAG_id3  GST_MAKE_FOURCC ('i','d','3',' ')
#define GST_RIFF_TAG_IDVX GST_MAKE_FOURCC ('I','D','V','X')
/* WAV stuff */
#define GST_RIFF_TAG_fmt  GST_MAKE_FOURCC ('f','m','t',' ')
#define GST_RIFF_TAG_data GST_MAKE_FOURCC ('d','a','t','a')
#define GST_RIFF_TAG_plst GST_MAKE_FOURCC ('p','l','s','t')
#define GST_RIFF_TAG_cue  GST_MAKE_FOURCC ('c','u','e',' ')
#define GST_RIFF_TAG_bext GST_MAKE_FOURCC ('b','e','x','t')
#define GST_RIFF_TAG_BEXT GST_MAKE_FOURCC ('B','E','X','T')
#define GST_RIFF_TAG_fact GST_MAKE_FOURCC ('f','a','c','t')
#define GST_RIFF_TAG_acid GST_MAKE_FOURCC ('a','c','i','d')
#define GST_RIFF_TAG_labl GST_MAKE_FOURCC ('l','a','b','l')
#define GST_RIFF_TAG_note GST_MAKE_FOURCC ('n','o','t','e')
#define GST_RIFF_TAG_smpl GST_MAKE_FOURCC ('s','m','p','l')
#define GST_RIFF_TAG_inst GST_MAKE_FOURCC ('i','n','s','t')
#define GST_RIFF_TAG_FLLR GST_MAKE_FOURCC ('F','L','L','R')

/* LIST types */
#define GST_RIFF_LIST_movi GST_MAKE_FOURCC ('m','o','v','i')
#define GST_RIFF_LIST_hdrl GST_MAKE_FOURCC ('h','d','r','l')
#define GST_RIFF_LIST_odml GST_MAKE_FOURCC ('o','d','m','l')
#define GST_RIFF_LIST_strl GST_MAKE_FOURCC ('s','t','r','l')
#define GST_RIFF_LIST_INFO GST_MAKE_FOURCC ('I','N','F','O')
#define GST_RIFF_LIST_AVIX GST_MAKE_FOURCC ('A','V','I','X')
#define GST_RIFF_LIST_adtl GST_MAKE_FOURCC ('a','d','t','l')
#define GST_RIFF_LIST_ncdt GST_MAKE_FOURCC ('n','c','d','t')

/* fcc types */
#define GST_RIFF_FCC_vids GST_MAKE_FOURCC ('v','i','d','s')
#define GST_RIFF_FCC_auds GST_MAKE_FOURCC ('a','u','d','s')
#define GST_RIFF_FCC_pads GST_MAKE_FOURCC ('p','a','d','s')
#define GST_RIFF_FCC_txts GST_MAKE_FOURCC ('t','x','t','s')
#define GST_RIFF_FCC_vidc GST_MAKE_FOURCC ('v','i','d','c')
#define GST_RIFF_FCC_iavs GST_MAKE_FOURCC ('i','a','v','s')
/* fcc handlers */
#define GST_RIFF_FCCH_RLE  GST_MAKE_FOURCC ('R','L','E',' ')
#define GST_RIFF_FCCH_msvc GST_MAKE_FOURCC ('m','s','v','c')
#define GST_RIFF_FCCH_MSVC GST_MAKE_FOURCC ('M','S','V','C')

/* INFO types - see http://www.saettler.com/RIFFMCI/riffmci.html */
#define GST_RIFF_INFO_IARL GST_MAKE_FOURCC ('I','A','R','L') /* location */
#define GST_RIFF_INFO_IART GST_MAKE_FOURCC ('I','A','R','T') /* artist */
#define GST_RIFF_INFO_ICMS GST_MAKE_FOURCC ('I','C','M','S') /* commissioned */
#define GST_RIFF_INFO_ICMT GST_MAKE_FOURCC ('I','C','M','T') /* comment */
#define GST_RIFF_INFO_ICOP GST_MAKE_FOURCC ('I','C','O','P') /* copyright */
#define GST_RIFF_INFO_ICRD GST_MAKE_FOURCC ('I','C','R','D') /* creation date */
#define GST_RIFF_INFO_ICRP GST_MAKE_FOURCC ('I','C','R','P') /* cropped */
#define GST_RIFF_INFO_IDIM GST_MAKE_FOURCC ('I','D','I','M') /* dimensions */
#define GST_RIFF_INFO_IDPI GST_MAKE_FOURCC ('I','D','P','I') /* dots-per-inch */
#define GST_RIFF_INFO_IENG GST_MAKE_FOURCC ('I','E','N','G') /* engineer(s) */
#define GST_RIFF_INFO_IGNR GST_MAKE_FOURCC ('I','G','N','R') /* genre */
#define GST_RIFF_INFO_IKEY GST_MAKE_FOURCC ('I','K','E','Y') /* keywords */
#define GST_RIFF_INFO_ILGT GST_MAKE_FOURCC ('I','L','G','T') /* lightness */
#define GST_RIFF_INFO_IMED GST_MAKE_FOURCC ('I','M','E','D') /* medium */
#define GST_RIFF_INFO_INAM GST_MAKE_FOURCC ('I','N','A','M') /* name */
#define GST_RIFF_INFO_IPLT GST_MAKE_FOURCC ('I','P','L','T') /* palette setting */
#define GST_RIFF_INFO_IPRD GST_MAKE_FOURCC ('I','P','R','D') /* product (album) */
#define GST_RIFF_INFO_ISBJ GST_MAKE_FOURCC ('I','S','B','J') /* subject */
#define GST_RIFF_INFO_ISFT GST_MAKE_FOURCC ('I','S','F','T') /* software */
#define GST_RIFF_INFO_ISHP GST_MAKE_FOURCC ('I','S','H','P') /* sharpness */
#define GST_RIFF_INFO_ISRC GST_MAKE_FOURCC ('I','S','R','C') /* source */
#define GST_RIFF_INFO_ISRF GST_MAKE_FOURCC ('I','S','R','F') /* source form */
#define GST_RIFF_INFO_ITCH GST_MAKE_FOURCC ('I','T','C','H') /* technician(s) */

#define GST_RIFF_INFO_IAAR GST_MAKE_FOURCC ('I','A','A','R') /* album artist */
#define GST_RIFF_INFO_ITRK GST_MAKE_FOURCC ('I','T','R','K') /* track number */

/* ncdt types - see http://www.sno.phy.queensu.ca/~phil/exiftool/TagNames/Nikon.html#NCDT */
#define GST_RIFF_LIST_nctg GST_MAKE_FOURCC ('n','c','t','g')

/*********Chunk Names***************/
#define GST_RIFF_FF00 GST_MAKE_FOURCC (0xFF,0xFF,0x00,0x00)
#define GST_RIFF_00   GST_MAKE_FOURCC ('0', '0',0x00,0x00)
#define GST_RIFF_01   GST_MAKE_FOURCC ('0', '1',0x00,0x00)
#define GST_RIFF_02   GST_MAKE_FOURCC ('0', '2',0x00,0x00)
#define GST_RIFF_03   GST_MAKE_FOURCC ('0', '3',0x00,0x00)
#define GST_RIFF_04   GST_MAKE_FOURCC ('0', '4',0x00,0x00)
#define GST_RIFF_05   GST_MAKE_FOURCC ('0', '5',0x00,0x00)
#define GST_RIFF_06   GST_MAKE_FOURCC ('0', '6',0x00,0x00)
#define GST_RIFF_07   GST_MAKE_FOURCC ('0', '7',0x00,0x00)
#define GST_RIFF_00pc GST_MAKE_FOURCC ('0', '0', 'p', 'c')
#define GST_RIFF_01pc GST_MAKE_FOURCC ('0', '1', 'p', 'c')
#define GST_RIFF_00dc GST_MAKE_FOURCC ('0', '0', 'd', 'c')
#define GST_RIFF_00dx GST_MAKE_FOURCC ('0', '0', 'd', 'x')
#define GST_RIFF_00db GST_MAKE_FOURCC ('0', '0', 'd', 'b')
#define GST_RIFF_00xx GST_MAKE_FOURCC ('0', '0', 'x', 'x')
#define GST_RIFF_00id GST_MAKE_FOURCC ('0', '0', 'i', 'd')
#define GST_RIFF_00rt GST_MAKE_FOURCC ('0', '0', 'r', 't')
#define GST_RIFF_0021 GST_MAKE_FOURCC ('0', '0', '2', '1')
#define GST_RIFF_00iv GST_MAKE_FOURCC ('0', '0', 'i', 'v')
#define GST_RIFF_0031 GST_MAKE_FOURCC ('0', '0', '3', '1')
#define GST_RIFF_0032 GST_MAKE_FOURCC ('0', '0', '3', '2')
#define GST_RIFF_00vc GST_MAKE_FOURCC ('0', '0', 'v', 'c')
#define GST_RIFF_00xm GST_MAKE_FOURCC ('0', '0', 'x', 'm')
#define GST_RIFF_01wb GST_MAKE_FOURCC ('0', '1', 'w', 'b')
#define GST_RIFF_01dc GST_MAKE_FOURCC ('0', '1', 'd', 'c')
#define GST_RIFF_00__ GST_MAKE_FOURCC ('0', '0', '_', '_')

/*********VIDEO CODECS**************/
#define GST_RIFF_cram GST_MAKE_FOURCC ('c', 'r', 'a', 'm')
#define GST_RIFF_CRAM GST_MAKE_FOURCC ('C', 'R', 'A', 'M')
#define GST_RIFF_wham GST_MAKE_FOURCC ('w', 'h', 'a', 'm')
#define GST_RIFF_WHAM GST_MAKE_FOURCC ('W', 'H', 'A', 'M')
#define GST_RIFF_rgb  GST_MAKE_FOURCC (0x00,0x00,0x00,0x00)
#define GST_RIFF_RGB  GST_MAKE_FOURCC ('R', 'G', 'B', ' ')
#define GST_RIFF_RAW  GST_MAKE_FOURCC ('R', 'A', 'W', ' ')
#define GST_RIFF_DIB  GST_MAKE_FOURCC ('D', 'I', 'B', ' ')
#define GST_RIFF_rle8 GST_MAKE_FOURCC (0x01,0x00,0x00,0x00)
#define GST_RIFF_RLE8 GST_MAKE_FOURCC ('R', 'L', 'E', '8')
#define GST_RIFF_rle4 GST_MAKE_FOURCC (0x02,0x00,0x00,0x00)
#define GST_RIFF_RLE4 GST_MAKE_FOURCC ('R', 'L', 'E', '4')
#define GST_RIFF_none GST_MAKE_FOURCC (0x00,0x00,0xFF,0xFF)
#define GST_RIFF_NONE GST_MAKE_FOURCC ('N', 'O', 'N', 'E')
#define GST_RIFF_pack GST_MAKE_FOURCC (0x01,0x00,0xFF,0xFF)
#define GST_RIFF_PACK GST_MAKE_FOURCC ('P', 'A', 'C', 'K')
#define GST_RIFF_tran GST_MAKE_FOURCC (0x02,0x00,0xFF,0xFF)
#define GST_RIFF_TRAN GST_MAKE_FOURCC ('T', 'R', 'A', 'N')
#define GST_RIFF_ccc  GST_MAKE_FOURCC (0x03,0x00,0xFF,0xFF)
#define GST_RIFF_CCC  GST_MAKE_FOURCC ('C', 'C', 'C', ' ')
#define GST_RIFF_cyuv GST_MAKE_FOURCC ('c', 'y', 'u', 'v')
#define GST_RIFF_CYUV GST_MAKE_FOURCC ('C', 'Y', 'U', 'V')
#define GST_RIFF_jpeg GST_MAKE_FOURCC (0x04,0x00,0xFF,0xFF)
#define GST_RIFF_JPEG GST_MAKE_FOURCC ('J', 'P', 'E', 'G')
#define GST_RIFF_MJPG GST_MAKE_FOURCC ('M', 'J', 'P', 'G')
#define GST_RIFF_mJPG GST_MAKE_FOURCC ('m', 'J', 'P', 'G')
#define GST_RIFF_IJPG GST_MAKE_FOURCC ('I', 'J', 'P', 'G')
#define GST_RIFF_rt21 GST_MAKE_FOURCC ('r', 't', '2', '1')
#define GST_RIFF_RT21 GST_MAKE_FOURCC ('R', 'T', '2', '1')
#define GST_RIFF_iv31 GST_MAKE_FOURCC ('i', 'v', '3', '1')
#define GST_RIFF_IV31 GST_MAKE_FOURCC ('I', 'V', '3', '1')
#define GST_RIFF_iv32 GST_MAKE_FOURCC ('i', 'v', '3', '2')
#define GST_RIFF_IV32 GST_MAKE_FOURCC ('I', 'V', '3', '2')
#define GST_RIFF_iv41 GST_MAKE_FOURCC ('i', 'v', '4', '1')
#define GST_RIFF_IV41 GST_MAKE_FOURCC ('I', 'V', '4', '1')
#define GST_RIFF_iv50 GST_MAKE_FOURCC ('i', 'v', '5', '0')
#define GST_RIFF_IV50 GST_MAKE_FOURCC ('I', 'V', '5', '0')
#define GST_RIFF_cvid GST_MAKE_FOURCC ('c', 'v', 'i', 'd')
#define GST_RIFF_CVID GST_MAKE_FOURCC ('C', 'V', 'I', 'D')
#define GST_RIFF_ULTI GST_MAKE_FOURCC ('U', 'L', 'T', 'I')
#define GST_RIFF_ulti GST_MAKE_FOURCC ('u', 'l', 't', 'i')
#define GST_RIFF_YUV9 GST_MAKE_FOURCC ('Y', 'U', 'V', '9')
#define GST_RIFF_YVU9 GST_MAKE_FOURCC ('Y', 'V', 'U', '9')
#define GST_RIFF_XMPG GST_MAKE_FOURCC ('X', 'M', 'P', 'G')
#define GST_RIFF_xmpg GST_MAKE_FOURCC ('x', 'm', 'p', 'g')
#define GST_RIFF_VDOW GST_MAKE_FOURCC ('V', 'D', 'O', 'W')
#define GST_RIFF_MVI1 GST_MAKE_FOURCC ('M', 'V', 'I', '1')
#define GST_RIFF_v422 GST_MAKE_FOURCC ('v', '4', '2', '2')
#define GST_RIFF_V422 GST_MAKE_FOURCC ('V', '4', '2', '2')
#define GST_RIFF_mvi1 GST_MAKE_FOURCC ('m', 'v', 'i', '1')
#define GST_RIFF_MPIX GST_MAKE_FOURCC (0x04,0x00, 'i', '1')     /* MotionPixels munged their id */
#define GST_RIFF_AURA GST_MAKE_FOURCC ('A', 'U', 'R', 'A')
#define GST_RIFF_DMB1 GST_MAKE_FOURCC ('D', 'M', 'B', '1')
#define GST_RIFF_dmb1 GST_MAKE_FOURCC ('d', 'm', 'b', '1')

#define GST_RIFF_BW10 GST_MAKE_FOURCC ('B', 'W', '1', '0')
#define GST_RIFF_bw10 GST_MAKE_FOURCC ('b', 'w', '1', '0')

#define GST_RIFF_yuy2 GST_MAKE_FOURCC ('y', 'u', 'y', '2')
#define GST_RIFF_YUY2 GST_MAKE_FOURCC ('Y', 'U', 'Y', '2')
#define GST_RIFF_YUV8 GST_MAKE_FOURCC ('Y', 'U', 'V', '8')
#define GST_RIFF_WINX GST_MAKE_FOURCC ('W', 'I', 'N', 'X')
#define GST_RIFF_WPY2 GST_MAKE_FOURCC ('W', 'P', 'Y', '2')
#define GST_RIFF_m263 GST_MAKE_FOURCC ('m', '2', '6', '3')
#define GST_RIFF_M263 GST_MAKE_FOURCC ('M', '2', '6', '3')
#define GST_RIFF_H263 GST_MAKE_FOURCC ('H', '2', '6', '3')
#define GST_RIFF_h263 GST_MAKE_FOURCC ('h', '2', '6', '3')
#define GST_RIFF_i263 GST_MAKE_FOURCC ('i', '2', '6', '3')
#define GST_RIFF_L263 GST_MAKE_FOURCC ('L', '2', '6', '3')
#define GST_RIFF_x263 GST_MAKE_FOURCC ('x', '2', '6', '3')
#define GST_RIFF_VSSH GST_MAKE_FOURCC ( 'V', 'S', 'S', 'H') /* H2.64 */

#define GST_RIFF_Q1_0 GST_MAKE_FOURCC ('Q', '1',0x2e, '0')
#define GST_RIFF_SFMC GST_MAKE_FOURCC ('S', 'F', 'M', 'C')

#define GST_RIFF_y41p GST_MAKE_FOURCC ('y', '4', '1', 'p')
#define GST_RIFF_Y41P GST_MAKE_FOURCC ('Y', '4', '1', 'P')
#define GST_RIFF_yv12 GST_MAKE_FOURCC ('y', 'v', '1', '2')
#define GST_RIFF_YV12 GST_MAKE_FOURCC ('Y', 'V', '1', '2')
#define GST_RIFF_vixl GST_MAKE_FOURCC ('v', 'i', 'x', 'l')
#define GST_RIFF_VIXL GST_MAKE_FOURCC ('V', 'I', 'X', 'L')
#define GST_RIFF_iyuv GST_MAKE_FOURCC ('i', 'y', 'u', 'v')
#define GST_RIFF_IYUV GST_MAKE_FOURCC ('I', 'Y', 'U', 'V')
#define GST_RIFF_i420 GST_MAKE_FOURCC ('i', '4', '2', '0')
#define GST_RIFF_I420 GST_MAKE_FOURCC ('I', '4', '2', '0')
#define GST_RIFF_vyuy GST_MAKE_FOURCC ('v', 'y', 'u', 'y')
#define GST_RIFF_VYUY GST_MAKE_FOURCC ('V', 'Y', 'U', 'Y')

#define GST_RIFF_DIV3 GST_MAKE_FOURCC ('D', 'I', 'V', '3')

#define GST_RIFF_rpza GST_MAKE_FOURCC ('r', 'p', 'z', 'a')
/* And this here's the mistakes that need to be supported */
#define GST_RIFF_azpr GST_MAKE_FOURCC ('a', 'z', 'p', 'r')  /* recognize Apple's rpza mangled? */

/*********** FND in MJPG **********/
#define GST_RIFF_ISFT GST_MAKE_FOURCC ('I', 'S', 'F', 'T')
#define GST_RIFF_IDIT GST_MAKE_FOURCC ('I', 'D', 'I', 'T')

#define GST_RIFF_00AM GST_MAKE_FOURCC ('0', '0', 'A', 'M')
#define GST_RIFF_DISP GST_MAKE_FOURCC ('D', 'I', 'S', 'P')
#define GST_RIFF_ISBJ GST_MAKE_FOURCC ('I', 'S', 'B', 'J')

#define GST_RIFF_rec  GST_MAKE_FOURCC ('r', 'e', 'c', ' ')

/* common data structures */
typedef struct _gst_riff_strh {
  guint32 type;             /* stream type */
  guint32 fcc_handler;       /* fcc_handler */
  guint32 flags;
/* flags values */
#define GST_RIFF_STRH_DISABLED        0x000000001
#define GST_RIFF_STRH_VIDEOPALCHANGES 0x000010000
  guint32 priority;
  guint32 init_frames;       /* initial frames (???) */
  guint32 scale;
  guint32 rate;
  guint32 start;
  guint32 length;
  guint32 bufsize;           /* suggested buffer size */
  guint32 quality;
  guint32 samplesize;
  /* rcFrame, RECT structure(struct of 4 shorts)
  gint32  left;
  gint32  top;
  gint32  right;
  gint32  bottom;
  */
} gst_riff_strh;

typedef struct _gst_riff_strf_vids {       /* == BitMapInfoHeader */
  guint32 size;
  guint32 width;
  guint32 height;
  guint16 planes;
  guint16 bit_cnt;
  guint32 compression;
  guint32 image_size;
  guint32 xpels_meter;
  guint32 ypels_meter;
  guint32 num_colors;        /* used colors */
  guint32 imp_colors;        /* important colors */
  /* may be more for some codecs */
} gst_riff_strf_vids;


typedef struct _gst_riff_strf_auds {       /* == WaveHeader (?) */
  guint16 format;
/**** from public Microsoft RIFF docs ******/
#define GST_RIFF_WAVE_FORMAT_UNKNOWN        (0x0000)
#define GST_RIFF_WAVE_FORMAT_PCM            (0x0001)
#define GST_RIFF_WAVE_FORMAT_ADPCM          (0x0002)
#define GST_RIFF_WAVE_FORMAT_IEEE_FLOAT     (0x0003)
#define GST_RIFF_WAVE_FORMAT_VSELP          (0x0004)
#define GST_RIFF_WAVE_FORMAT_IBM_CVSD       (0x0005)
#define GST_RIFF_WAVE_FORMAT_ALAW           (0x0006)
#define GST_RIFF_WAVE_FORMAT_MULAW          (0x0007)
#define GST_RIFF_WAVE_FORMAT_WMS            (0x000a) /* WMS Windows Media Audio Speech */
#define GST_RIFF_WAVE_FORMAT_OKI_ADPCM      (0x0010)
#define GST_RIFF_WAVE_FORMAT_DVI_ADPCM      (0x0011)
#define GST_RIFF_WAVE_FORMAT_MEDIASPACE_ADPCM (0x0012)
#define GST_RIFF_WAVE_FORMAT_SIERRA_ADPCM   (0x0013)
#define GST_RIFF_WAVE_FORMAT_G723_ADPCM     (0x0014)
#define GST_RIFF_WAVE_FORMAT_DIGISTD        (0x0015)
#define GST_RIFF_WAVE_FORMAT_DIGIFIX        (0x0016)
#define GST_RIFF_WAVE_FORMAT_DIALOGIC_OKI_ADPCM (0x0017)
#define GST_RIFF_WAVE_FORMAT_MEDIAVISION_ADPCM  (0x0018)
#define GST_RIFF_WAVE_FORMAT_CU_CODEC       (0x0019)
#define GST_RIFF_WAVE_FORMAT_YAMAHA_ADPCM   (0x0020)
#define GST_RIFF_WAVE_FORMAT_SONARC         (0x0021)
#define GST_RIFF_WAVE_FORMAT_DSP_TRUESPEECH (0x0022)
#define GST_RIFF_WAVE_FORMAT_ECHOSC1        (0x0023)
#define GST_RIFF_WAVE_FORMAT_AUDIOFILE_AF36 (0x0024)
#define GST_RIFF_WAVE_FORMAT_APTX           (0x0025)
#define GST_RIFF_WAVE_FORMAT_AUDIOFILE_AF10 (0x0026)
#define GST_RIFF_WAVE_FORMAT_PROSODY_1612   (0x0027)
#define GST_RIFF_WAVE_FORMAT_LRC            (0x0028)
#define GST_RIFF_WAVE_FORMAT_DOLBY_AC2      (0x0030)
#define GST_RIFF_WAVE_FORMAT_GSM610         (0x0031)
#define GST_RIFF_WAVE_FORMAT_MSN            (0x0032)
#define GST_RIFF_WAVE_FORMAT_ANTEX_ADPCME   (0x0033)
#define GST_RIFF_WAVE_FORMAT_CONTROL_RES_VQLPC (0x0034)
#define GST_RIFF_WAVE_FORMAT_DIGIREAL       (0x0035)
#define GST_RIFF_WAVE_FORMAT_DIGIADPCM      (0x0036)
#define GST_RIFF_WAVE_FORMAT_CONTROL_RES_CR10  (0x0037)
#define GST_RIFF_WAVE_FORMAT_NMS_VBXADPCM   (0x0038)
#define GST_RIFF_WAVE_FORMAT_CS_IMAADPCM    (0x0039)
#define GST_RIFF_WAVE_FORMAT_ECHOSC3        (0x003A)
#define GST_RIFF_WAVE_FORMAT_ROCKWELL_ADPCM (0x003B)
#define GST_RIFF_WAVE_FORMAT_ROCKWELL_DIGITALK (0x003C)
#define GST_RIFF_WAVE_FORMAT_XEBEC          (0x003D)
#define GST_RIFF_WAVE_FORMAT_ITU_G721_ADPCM (0x0040)
#define GST_RIFF_WAVE_FORMAT_G728_CELP      (0x0041)
#define GST_RIFF_WAVE_FORMAT_MSG723         (0x0042)
#define GST_RIFF_WAVE_FORMAT_ITU_G726_ADPCM (0x0045)
#define GST_RIFF_WAVE_FORMAT_MPEGL12        (0x0050)
#define GST_RIFF_WAVE_FORMAT_RT24           (0x0052)
#define GST_RIFF_WAVE_FORMAT_PAC            (0x0053)
#define GST_RIFF_WAVE_FORMAT_MPEGL3         (0x0055)
#define GST_RIFF_WAVE_FORMAT_AMR_NB         (0x0057)
#define GST_RIFF_WAVE_FORMAT_AMR_WB         (0x0058)
#define GST_RIFF_WAVE_FORMAT_LUCENT_G723    (0x0059)
#define GST_RIFF_WAVE_FORMAT_CIRRUS         (0x0060)
#define GST_RIFF_WAVE_FORMAT_ADPCM_IMA_DK4  (0x0061)  /* not official */
#define GST_RIFF_WAVE_FORMAT_ADPCM_IMA_DK3  (0x0062)  /* not official */
/* FIXME: where are these from? are they used at all? */
#if 0
#define GST_RIFF_WAVE_FORMAT_ESPCM          (0x0061)
#define GST_RIFF_WAVE_FORMAT_VOXWARE        (0x0062)
#endif
#define GST_RIFF_WAVE_FORMAT_CANOPUS_ATRAC  (0x0063)
#define GST_RIFF_WAVE_FORMAT_G726_ADPCM     (0x0064)
#define GST_RIFF_WAVE_FORMAT_G722_ADPCM     (0x0065)
#define GST_RIFF_WAVE_FORMAT_ADPCM_G722     (0x028F)
#define GST_RIFF_WAVE_FORMAT_DSAT_DISPLAY   (0x0067)
#define GST_RIFF_WAVE_FORMAT_ADPCM_IMA_WAV (0x0069)
/* FIXME: where are these from? are they used at all? */
#if 0
#define GST_RIFF_WAVE_FORMAT_VOXWARE_BYTE_ALIGNED (0x0069)
#endif
#define GST_RIFF_WAVE_FORMAT_VOXWARE_AC8    (0x0070)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_AC10   (0x0071)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_AC16   (0x0072)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_AC20   (0x0073)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_METAVOICE (0x0074)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_METASOUND (0x0075)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_RT29HW (0x0076)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_VR12   (0x0077)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_VR18   (0x0078)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_TQ40   (0x0079)
#define GST_RIFF_WAVE_FORMAT_SOFTSOUND      (0x0080)
#define GST_RIFF_WAVE_FORMAT_VOXWARE_TQ60   (0x0081)
#define GST_RIFF_WAVE_FORMAT_MSRT24         (0x0082)
#define GST_RIFF_WAVE_FORMAT_G729A          (0x0083)
#define GST_RIFF_WAVE_FORMAT_MVI_MVI2       (0x0084)
#define GST_RIFF_WAVE_FORMAT_DF_G726        (0x0085)
#define GST_RIFF_WAVE_FORMAT_DF_GSM610      (0x0086)
#define GST_RIFF_WAVE_FORMAT_ISIAUDIO       (0x0088)
#define GST_RIFF_WAVE_FORMAT_ONLIVE         (0x0089)
#define GST_RIFF_WAVE_FORMAT_SBC24          (0x0091)
#define GST_RIFF_WAVE_FORMAT_DOLBY_AC3_SPDIF  (0x0092)
#define GST_RIFF_WAVE_FORMAT_MEDIASONIC_G723  (0x0093)
#define GST_RIFF_WAVE_FORMAT_PROSODY_8KBPS  (0x0094)
#define GST_RIFF_WAVE_FORMAT_ZYXEL_ADPCM    (0x0097)
#define GST_RIFF_WAVE_FORMAT_PHILIPS_LPCBB  (0x0098)
#define GST_RIFF_WAVE_FORMAT_PACKED         (0x0099)
#define GST_RIFF_WAVE_FORMAT_MALDEN_PHONYTALK (0x00A0)
#define GST_RIFF_WAVE_FORMAT_AAC            (0x00ff)
#define GST_RIFF_WAVE_FORMAT_RHETOREX_ADPCM (0x0100)
#define GST_RIFF_IBM_FORMAT_MULAW           (0x0101)
#define GST_RIFF_IBM_FORMAT_ALAW            (0x0102)
#define GST_RIFF_IBM_FORMAT_ADPCM           (0x0103)
#define GST_RIFF_WAVE_FORMAT_VIVO_G723      (0x0111)
#define GST_RIFF_WAVE_FORMAT_VIVO_SIREN     (0x0112)
#define GST_RIFF_WAVE_FORMAT_DIGITAL_G723   (0x0123)
#define GST_RIFF_WAVE_FORMAT_SANYO_LD_ADPCM (0x0125)
#define GST_RIFF_WAVE_FORMAT_SIPROLAB_ACEPLNET  (0x0130)
#define GST_RIFF_WAVE_FORMAT_SIPROLAB_ACELP4800 (0x0131)
#define GST_RIFF_WAVE_FORMAT_SIPROLAB_ACELP8V3  (0x0132)
#define GST_RIFF_WAVE_FORMAT_SIPROLAB_G729   (0x0133)
#define GST_RIFF_WAVE_FORMAT_SIPROLAB_G729A  (0x0134)
#define GST_RIFF_WAVE_FORMAT_SIPROLAB_KELVIN (0x0135)
#define GST_RIFF_WAVE_FORMAT_G726ADPCM       (0x0140)
#define GST_RIFF_WAVE_FORMAT_QUALCOMM_PUREVOICE (0x0150)
#define GST_RIFF_WAVE_FORMAT_QUALCOMM_HALFRATE  (0x0151)
#define GST_RIFF_WAVE_FORMAT_TUBGSM             (0x0155)
#define GST_RIFF_WAVE_FORMAT_WMAV1          (0x0160)
#define GST_RIFF_WAVE_FORMAT_WMAV2          (0x0161)
#define GST_RIFF_WAVE_FORMAT_WMAV3          (0x0162)
#define GST_RIFF_WAVE_FORMAT_WMAV3_L        (0x0163)
#define GST_RIFF_WAVE_FORMAT_CREATIVE_ADPCM (0x0200)
#define GST_RIFF_WAVE_FORMAT_CREATIVE_FASTSPEECH8  (0x0202)
#define GST_RIFF_WAVE_FORMAT_CREATIVE_FASTSPEECH10 (0x0203)
#define GST_RIFF_WAVE_FORMAT_UHER_ADPCM     (0x0210)
#define GST_RIFF_WAVE_FORMAT_QUARTERDECK    (0x0220)
#define GST_RIFF_WAVE_FORMAT_ILINK_VC       (0x0230)
#define GST_RIFF_WAVE_FORMAT_RAW_SPORT      (0x0240)
#define GST_RIFF_WAVE_FORMAT_IPI_HSX        (0x0250)
#define GST_RIFF_WAVE_FORMAT_IPI_RPELP      (0x0251)
#define GST_RIFF_WAVE_FORMAT_CS2            (0x0260)
#define GST_RIFF_WAVE_FORMAT_SONY_ATRAC3    (0x0270)
#define GST_RIFF_WAVE_FORMAT_SIREN          (0x028E)
#define GST_RIFF_WAVE_FORMAT_FM_TOWNS_SND   (0x0300)
#define GST_RIFF_WAVE_FORMAT_BTV_DIGITAL    (0x0400)
#define GST_RIFF_WAVE_FORMAT_IMC            (0x0401)
#define GST_RIFF_WAVE_FORMAT_QDESIGN_MUSIC  (0x0450)
#define GST_RIFF_WAVE_FORMAT_VME_VMPCM      (0x0680)
#define GST_RIFF_WAVE_FORMAT_TPC            (0x0681)
#define GST_RIFF_WAVE_FORMAT_OLIGSM         (0x1000)
#define GST_RIFF_WAVE_FORMAT_OLIADPCM       (0x1001)
#define GST_RIFF_WAVE_FORMAT_OLICELP        (0x1002)
#define GST_RIFF_WAVE_FORMAT_OLISBC         (0x1003)
#define GST_RIFF_WAVE_FORMAT_OLIOPR         (0x1004)
#define GST_RIFF_WAVE_FORMAT_LH_CODEC       (0x1100)
#define GST_RIFF_WAVE_FORMAT_NORRIS         (0x1400)
#define GST_RIFF_WAVE_FORMAT_SOUNDSPACE_MUSICOMPRESS (0x1500)
#define GST_RIFF_WAVE_FORMAT_A52            (0x2000)
#define GST_RIFF_WAVE_FORMAT_DTS            (0x2001)
#define GST_RIFF_WAVE_FORMAT_SONIC          (0x2048)
#define GST_RIFF_WAVE_FORMAT_SONIC_LS       (0x2048)
#define GST_RIFF_WAVE_FORMAT_AAC_AC         (0x4143)
#define GST_RIFF_WAVE_FORMAT_VORBIS1        (0x674f)
#define GST_RIFF_WAVE_FORMAT_VORBIS2        (0x6750)
#define GST_RIFF_WAVE_FORMAT_VORBIS3        (0x6751)
#define GST_RIFF_WAVE_FORMAT_VORBIS1PLUS    (0x676f)
#define GST_RIFF_WAVE_FORMAT_VORBIS2PLUS    (0x6770)
#define GST_RIFF_WAVE_FORMAT_VORBIS3PLUS    (0x6771)
#define GST_RIFF_WAVE_FORMAT_AAC_pm         (0x706d)
#define GST_RIFF_WAVE_FORMAT_GSM_AMR_CBR    (0x7A21)
#define GST_RIFF_WAVE_FORMAT_GSM_AMR_VBR    (0x7A22)
#define GST_RIFF_WAVE_FORMAT_FLAC           (0xF1AC)
#define GST_RIFF_WAVE_FORMAT_EXTENSIBLE     (0xFFFE)
  guint16 channels;
  guint32 rate;
  guint32 av_bps;
  guint16 blockalign;
  guint16 bits_per_sample;
#if 0
  /* missing field */
  guint16 extra_size;
#endif
} gst_riff_strf_auds;

typedef struct _gst_riff_strf_iavs {
  guint32 DVAAuxSrc;
  guint32 DVAAuxCtl;
  guint32 DVAAuxSrc1;
  guint32 DVAAuxCtl1;
  guint32 DVVAuxSrc;
  guint32 DVVAuxCtl;
  guint32 DVReserved1;
  guint32 DVReserved2;
} gst_riff_strf_iavs;

typedef struct _gst_riff_index_entry {
  guint32 id;
  guint32 flags;
#define GST_RIFF_IF_LIST                (0x00000001L)
#define GST_RIFF_IF_KEYFRAME            (0x00000010L)
#define GST_RIFF_IF_NO_TIME             (0x00000100L)
#define GST_RIFF_IF_COMPUSE             (0x0FFF0000L)
  guint32 offset;
  guint32 size;
} gst_riff_index_entry;

typedef struct _gst_riff_dmlh {
  guint32 totalframes;
} gst_riff_dmlh;

/* taken from libsndfile/wav.c (LGPL) */
typedef struct _gst_riff_acid {
  /* 4 bytes (int)     type of file:
   *  this appears to be a bit mask,however some combinations
   *  are probably impossible and/or qualified as "errors"
   *
   *  0x01 On: One Shot         Off: Loop
   *  0x02 On: Root note is Set Off: No root
   *  0x04 On: Stretch is On,   Off: Stretch is OFF
   *  0x08 On: Disk Based       Off: Ram based
   *  0x10 On: ??????????       Off: ????????? (Acidizer puts that ON)
   */
  guint32 loop_type;
  /* 2 bytes (short)      root note
   *  if type 0x10 is OFF : [C,C#,(...),B] -> [0x30 to 0x3B]
   *  if type 0x10 is ON  : [C,C#,(...),B] -> [0x3C to 0x47]
   *  (both types fit on same MIDI pitch albeit different octaves, so who cares)
   */
  guint16 root_note;
  /* 2 bytes (short)      ??? always set to 0x8000
   * 4 bytes (float)      ??? seems to be always 0
   */
  guint16 unknown1;
  gfloat unknown2;
  /* 4 bytes (int)        number of beats
   * 2 bytes (short)      meter denominator   //always 4 in SF/ACID
   * 2 bytes (short)      meter numerator     //always 4 in SF/ACID
   *                      //are we sure about the order?? usually its num/denom
   * 4 bytes (float)      tempo
   */
  guint32 number_of_beats;
  guint16 meter_d, meter_n;
  gfloat tempo;
} gst_riff_acid;

G_END_DECLS

#endif /* __GST_RIFF_IDS_H__ */

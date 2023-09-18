/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QVIDEOFRAMECONVERSIONHELPER_P_H
#define QVIDEOFRAMECONVERSIONHELPER_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <qvideoframe.h>
#include <private/qsimd_p.h>

typedef void (QT_FASTCALL *VideoFrameConvertFunc)(const QVideoFrame &frame, uchar *output);

inline quint32 qConvertBGRA32ToARGB32(quint32 bgra)
{
    return (((bgra & 0xFF000000) >> 24)
            | ((bgra & 0x00FF0000) >> 8)
            | ((bgra & 0x0000FF00) << 8)
            | ((bgra & 0x000000FF) << 24));
}

inline quint32 qConvertBGR24ToARGB32(const uchar *bgr)
{
    return 0xFF000000 | bgr[0] | bgr[1] << 8 | bgr[2] << 16;
}

inline quint32 qConvertBGR565ToARGB32(quint16 bgr)
{
    return 0xff000000
            | ((((bgr) >> 8) & 0xf8) | (((bgr) >> 13) & 0x7))
            | ((((bgr) << 5) & 0xfc00) | (((bgr) >> 1) & 0x300))
            | ((((bgr) << 19) & 0xf80000) | (((bgr) << 14) & 0x70000));
}

inline quint32 qConvertBGR555ToARGB32(quint16 bgr)
{
    return 0xff000000
            | ((((bgr) >> 7) & 0xf8) | (((bgr) >> 12) & 0x7))
            | ((((bgr) << 6) & 0xf800) | (((bgr) << 1) & 0x700))
            | ((((bgr) << 19) & 0xf80000) | (((bgr) << 11) & 0x70000));
}

#define FETCH_INFO_PACKED(frame) \
    const uchar *src = frame.bits(); \
    int stride = frame.bytesPerLine(); \
    int width = frame.width(); \
    int height = frame.height();

#define FETCH_INFO_BIPLANAR(frame) \
    const uchar *plane1 = frame.bits(0); \
    const uchar *plane2 = frame.bits(1); \
    int plane1Stride = frame.bytesPerLine(0); \
    int plane2Stride = frame.bytesPerLine(1); \
    int width = frame.width(); \
    int height = frame.height();

#define FETCH_INFO_TRIPLANAR(frame) \
    const uchar *plane1 = frame.bits(0); \
    const uchar *plane2 = frame.bits(1); \
    const uchar *plane3 = frame.bits(2); \
    int plane1Stride = frame.bytesPerLine(0); \
    int plane2Stride = frame.bytesPerLine(1); \
    int plane3Stride = frame.bytesPerLine(2); \
    int width = frame.width(); \
    int height = frame.height(); \

#define MERGE_LOOPS(width, height, stride, bpp) \
    if (stride == width * bpp) { \
        width *= height; \
        height = 1; \
        stride = 0; \
    }

#define ALIGN(boundary, ptr, x, length) \
    for (; ((reinterpret_cast<qintptr>(ptr) & (boundary - 1)) != 0) && x < length; ++x)

#endif // QVIDEOFRAMECONVERSIONHELPER_P_H


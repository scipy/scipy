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

#ifndef QVIDEOFRAME_H
#define QVIDEOFRAME_H

#include <QtCore/qmetatype.h>
#include <QtCore/qshareddata.h>
#include <QtGui/qimage.h>
#include <QtMultimedia/qabstractvideobuffer.h>
#include <QtCore/qvariant.h>

QT_BEGIN_NAMESPACE

class QSize;

class QVideoFramePrivate;

class Q_MULTIMEDIA_EXPORT QVideoFrame
{
public:
    enum FieldType
    {
        ProgressiveFrame,
        TopField,
        BottomField,
        InterlacedFrame
    };

    enum PixelFormat
    {
        Format_Invalid,
        Format_ARGB32,
        Format_ARGB32_Premultiplied,
        Format_RGB32,
        Format_RGB24,
        Format_RGB565,
        Format_RGB555,
        Format_ARGB8565_Premultiplied,
        Format_BGRA32,
        Format_BGRA32_Premultiplied,
        Format_BGR32,
        Format_BGR24,
        Format_BGR565,
        Format_BGR555,
        Format_BGRA5658_Premultiplied,

        Format_AYUV444,
        Format_AYUV444_Premultiplied,
        Format_YUV444,
        Format_YUV420P,
        Format_YV12,
        Format_UYVY,
        Format_YUYV,
        Format_NV12,
        Format_NV21,
        Format_IMC1,
        Format_IMC2,
        Format_IMC3,
        Format_IMC4,
        Format_Y8,
        Format_Y16,

        Format_Jpeg,

        Format_CameraRaw,
        Format_AdobeDng,
        Format_ABGR32, // ### Qt 6: reorder
        Format_YUV422P,

#ifndef Q_QDOC
        NPixelFormats,
#endif
        Format_User = 1000
    };

    QVideoFrame();
    QVideoFrame(QAbstractVideoBuffer *buffer, const QSize &size, PixelFormat format);
    QVideoFrame(int bytes, const QSize &size, int bytesPerLine, PixelFormat format);
    QVideoFrame(const QImage &image);
    QVideoFrame(const QVideoFrame &other);
    ~QVideoFrame();

    QVideoFrame &operator =(const QVideoFrame &other);
    bool operator==(const QVideoFrame &other) const;
    bool operator!=(const QVideoFrame &other) const;

    QAbstractVideoBuffer *buffer() const;
    bool isValid() const;

    PixelFormat pixelFormat() const;

    QAbstractVideoBuffer::HandleType handleType() const;

    QSize size() const;
    int width() const;
    int height() const;

    FieldType fieldType() const;
    void setFieldType(FieldType);

    bool isMapped() const;
    bool isReadable() const;
    bool isWritable() const;

    QAbstractVideoBuffer::MapMode mapMode() const;

    bool map(QAbstractVideoBuffer::MapMode mode);
    void unmap();

    int bytesPerLine() const;
    int bytesPerLine(int plane) const;

    uchar *bits();
    uchar *bits(int plane);
    const uchar *bits() const;
    const uchar *bits(int plane) const;
    int mappedBytes() const;
    int planeCount() const;

    QVariant handle() const;

    qint64 startTime() const;
    void setStartTime(qint64 time);

    qint64 endTime() const;
    void setEndTime(qint64 time);

    QVariantMap availableMetaData() const;
    QVariant metaData(const QString &key) const;
    void setMetaData(const QString &key, const QVariant &value);

    QImage image() const;

    static PixelFormat pixelFormatFromImageFormat(QImage::Format format);
    static QImage::Format imageFormatFromPixelFormat(PixelFormat format);

private:
    QExplicitlySharedDataPointer<QVideoFramePrivate> d;
};

#ifndef QT_NO_DEBUG_STREAM
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug, const QVideoFrame&);
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug, QVideoFrame::FieldType);
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug, QVideoFrame::PixelFormat);
#endif

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QVideoFrame)
Q_DECLARE_METATYPE(QVideoFrame::FieldType)
Q_DECLARE_METATYPE(QVideoFrame::PixelFormat)

#endif


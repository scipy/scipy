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

#ifndef QVIDEOSURFACEFORMAT_H
#define QVIDEOSURFACEFORMAT_H

#include <QtCore/qlist.h>
#include <QtCore/qpair.h>
#include <QtCore/qshareddata.h>
#include <QtCore/qsize.h>
#include <QtGui/qimage.h>
#include <QtMultimedia/qvideoframe.h>

QT_BEGIN_NAMESPACE


class QDebug;

class QVideoSurfaceFormatPrivate;

class Q_MULTIMEDIA_EXPORT QVideoSurfaceFormat
{
public:
    enum Direction
    {
        TopToBottom,
        BottomToTop
    };

    enum YCbCrColorSpace
    {
        YCbCr_Undefined,
        YCbCr_BT601,
        YCbCr_BT709,
        YCbCr_xvYCC601,
        YCbCr_xvYCC709,
        YCbCr_JPEG,
#ifndef Q_QDOC
        YCbCr_CustomMatrix
#endif
    };

    QVideoSurfaceFormat();
    QVideoSurfaceFormat(
            const QSize &size,
            QVideoFrame::PixelFormat pixelFormat,
            QAbstractVideoBuffer::HandleType handleType = QAbstractVideoBuffer::NoHandle);
    QVideoSurfaceFormat(const QVideoSurfaceFormat &format);
    ~QVideoSurfaceFormat();

    QVideoSurfaceFormat &operator =(const QVideoSurfaceFormat &format);

    bool operator ==(const QVideoSurfaceFormat &format) const;
    bool operator !=(const QVideoSurfaceFormat &format) const;

    bool isValid() const;

    QVideoFrame::PixelFormat pixelFormat() const;
    QAbstractVideoBuffer::HandleType handleType() const;

    QSize frameSize() const;
    void setFrameSize(const QSize &size);
    void setFrameSize(int width, int height);

    int frameWidth() const;
    int frameHeight() const;

    QRect viewport() const;
    void setViewport(const QRect &viewport);

    Direction scanLineDirection() const;
    void setScanLineDirection(Direction direction);

    qreal frameRate() const;
    void setFrameRate(qreal rate);

    QSize pixelAspectRatio() const;
    void setPixelAspectRatio(const QSize &ratio);
    void setPixelAspectRatio(int width, int height);

    YCbCrColorSpace yCbCrColorSpace() const;
    void setYCbCrColorSpace(YCbCrColorSpace colorSpace);

    bool isMirrored() const;
    void setMirrored(bool mirrored);

    QSize sizeHint() const;

    QList<QByteArray> propertyNames() const;
    QVariant property(const char *name) const;
    void setProperty(const char *name, const QVariant &value);

private:
    QSharedDataPointer<QVideoSurfaceFormatPrivate> d;
};

#ifndef QT_NO_DEBUG_STREAM
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug, const QVideoSurfaceFormat &);
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug, QVideoSurfaceFormat::Direction);
Q_MULTIMEDIA_EXPORT QDebug operator<<(QDebug, QVideoSurfaceFormat::YCbCrColorSpace);
#endif

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QVideoSurfaceFormat)
Q_DECLARE_METATYPE(QVideoSurfaceFormat::Direction)
Q_DECLARE_METATYPE(QVideoSurfaceFormat::YCbCrColorSpace)

#endif


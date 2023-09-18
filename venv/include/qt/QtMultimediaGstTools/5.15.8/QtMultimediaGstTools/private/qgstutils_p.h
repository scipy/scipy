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

#ifndef QGSTUTILS_P_H
#define QGSTUTILS_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API. It exists purely as an
// implementation detail. This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <private/qgsttools_global_p.h>
#include <QtCore/qmap.h>
#include <QtCore/qset.h>
#include <QtCore/qvector.h>
#include <gst/gst.h>
#include <gst/video/video.h>
#include <qaudioformat.h>
#include <qcamera.h>
#include <qabstractvideobuffer.h>
#include <qvideoframe.h>
#include <QDebug>

#if GST_CHECK_VERSION(1,0,0)
# define QT_GSTREAMER_PLAYBIN_ELEMENT_NAME "playbin"
# define QT_GSTREAMER_CAMERABIN_ELEMENT_NAME "camerabin"
# define QT_GSTREAMER_COLORCONVERSION_ELEMENT_NAME "videoconvert"
# define QT_GSTREAMER_RAW_AUDIO_MIME "audio/x-raw"
# define QT_GSTREAMER_VIDEOOVERLAY_INTERFACE_NAME "GstVideoOverlay"
#else
# define QT_GSTREAMER_PLAYBIN_ELEMENT_NAME "playbin2"
# define QT_GSTREAMER_CAMERABIN_ELEMENT_NAME "camerabin2"
# define QT_GSTREAMER_COLORCONVERSION_ELEMENT_NAME "ffmpegcolorspace"
# define QT_GSTREAMER_RAW_AUDIO_MIME "audio/x-raw-int"
# define QT_GSTREAMER_VIDEOOVERLAY_INTERFACE_NAME "GstXOverlay"
#endif

QT_BEGIN_NAMESPACE

class QSize;
class QVariant;
class QByteArray;
class QImage;
class QVideoSurfaceFormat;

namespace QGstUtils {
    struct Q_GSTTOOLS_EXPORT CameraInfo
    {
        QString name;
        QString description;
        int orientation;
        QCamera::Position position;
        QByteArray driver;
    };

    Q_GSTTOOLS_EXPORT QMap<QByteArray, QVariant> gstTagListToMap(const GstTagList *list);

    Q_GSTTOOLS_EXPORT QSize capsResolution(const GstCaps *caps);
    Q_GSTTOOLS_EXPORT QSize capsCorrectedResolution(const GstCaps *caps);
    Q_GSTTOOLS_EXPORT QAudioFormat audioFormatForCaps(const GstCaps *caps);
#if GST_CHECK_VERSION(1,0,0)
    Q_GSTTOOLS_EXPORT QAudioFormat audioFormatForSample(GstSample *sample);
#else
    Q_GSTTOOLS_EXPORT QAudioFormat audioFormatForBuffer(GstBuffer *buffer);
#endif
    Q_GSTTOOLS_EXPORT GstCaps *capsForAudioFormat(const QAudioFormat &format);
    Q_GSTTOOLS_EXPORT void initializeGst();
    Q_GSTTOOLS_EXPORT QMultimedia::SupportEstimate hasSupport(const QString &mimeType,
                                             const QStringList &codecs,
                                             const QSet<QString> &supportedMimeTypeSet);

    Q_GSTTOOLS_EXPORT QVector<CameraInfo> enumerateCameras(GstElementFactory *factory = 0);
    Q_GSTTOOLS_EXPORT QList<QByteArray> cameraDevices(GstElementFactory * factory = 0);
    Q_GSTTOOLS_EXPORT QString cameraDescription(const QString &device, GstElementFactory * factory = 0);
    Q_GSTTOOLS_EXPORT QCamera::Position cameraPosition(const QString &device, GstElementFactory * factory = 0);
    Q_GSTTOOLS_EXPORT int cameraOrientation(const QString &device, GstElementFactory * factory = 0);
    Q_GSTTOOLS_EXPORT QByteArray cameraDriver(const QString &device, GstElementFactory * factory = 0);

    Q_GSTTOOLS_EXPORT QSet<QString> supportedMimeTypes(bool (*isValidFactory)(GstElementFactory *factory));

#if GST_CHECK_VERSION(1,0,0)
    Q_GSTTOOLS_EXPORT QImage bufferToImage(GstBuffer *buffer, const GstVideoInfo &info);
    Q_GSTTOOLS_EXPORT QVideoSurfaceFormat formatForCaps(
            GstCaps *caps,
            GstVideoInfo *info = 0,
            QAbstractVideoBuffer::HandleType handleType = QAbstractVideoBuffer::NoHandle);
#else
    Q_GSTTOOLS_EXPORT QImage bufferToImage(GstBuffer *buffer);
    Q_GSTTOOLS_EXPORT QVideoSurfaceFormat formatForCaps(
            GstCaps *caps,
            int *bytesPerLine = 0,
            QAbstractVideoBuffer::HandleType handleType = QAbstractVideoBuffer::NoHandle);
#endif

    Q_GSTTOOLS_EXPORT GstCaps *capsForFormats(const QList<QVideoFrame::PixelFormat> &formats);
    void setFrameTimeStamps(QVideoFrame *frame, GstBuffer *buffer);

    Q_GSTTOOLS_EXPORT void setMetaData(GstElement *element, const QMap<QByteArray, QVariant> &data);
    Q_GSTTOOLS_EXPORT void setMetaData(GstBin *bin, const QMap<QByteArray, QVariant> &data);

    Q_GSTTOOLS_EXPORT GstCaps *videoFilterCaps();

    Q_GSTTOOLS_EXPORT QSize structureResolution(const GstStructure *s);
    Q_GSTTOOLS_EXPORT QVideoFrame::PixelFormat structurePixelFormat(const GstStructure *s, int *bpp = 0);
    Q_GSTTOOLS_EXPORT QSize structurePixelAspectRatio(const GstStructure *s);
    Q_GSTTOOLS_EXPORT QPair<qreal, qreal> structureFrameRateRange(const GstStructure *s);

    Q_GSTTOOLS_EXPORT QString fileExtensionForMimeType(const QString &mimeType);

#if GST_CHECK_VERSION(0,10,30)
    Q_GSTTOOLS_EXPORT QVariant fromGStreamerOrientation(const QVariant &value);
    Q_GSTTOOLS_EXPORT QVariant toGStreamerOrientation(const QVariant &value);
#endif

    Q_GSTTOOLS_EXPORT bool useOpenGL();
}

Q_GSTTOOLS_EXPORT void qt_gst_object_ref_sink(gpointer object);
Q_GSTTOOLS_EXPORT GstCaps *qt_gst_pad_get_current_caps(GstPad *pad);
Q_GSTTOOLS_EXPORT GstCaps *qt_gst_pad_get_caps(GstPad *pad);
Q_GSTTOOLS_EXPORT GstStructure *qt_gst_structure_new_empty(const char *name);
Q_GSTTOOLS_EXPORT gboolean qt_gst_element_query_position(GstElement *element, GstFormat format, gint64 *cur);
Q_GSTTOOLS_EXPORT gboolean qt_gst_element_query_duration(GstElement *element, GstFormat format, gint64 *cur);
Q_GSTTOOLS_EXPORT GstCaps *qt_gst_caps_normalize(GstCaps *caps);
Q_GSTTOOLS_EXPORT const gchar *qt_gst_element_get_factory_name(GstElement *element);
Q_GSTTOOLS_EXPORT gboolean qt_gst_caps_can_intersect(const GstCaps * caps1, const GstCaps * caps2);
Q_GSTTOOLS_EXPORT GList *qt_gst_video_sinks();
Q_GSTTOOLS_EXPORT void qt_gst_util_double_to_fraction(gdouble src, gint *dest_n, gint *dest_d);

Q_GSTTOOLS_EXPORT QDebug operator <<(QDebug debug, GstCaps *caps);

QT_END_NAMESPACE

#endif

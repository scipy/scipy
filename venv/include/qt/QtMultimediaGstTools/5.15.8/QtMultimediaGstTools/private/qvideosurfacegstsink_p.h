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

#ifndef VIDEOSURFACEGSTSINK_P_H
#define VIDEOSURFACEGSTSINK_P_H

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

#include <gst/gst.h>

#if GST_CHECK_VERSION(1,0,0)

#include "qgstvideorenderersink_p.h"

QT_BEGIN_NAMESPACE
typedef QGstVideoRendererSink QVideoSurfaceGstSink;
QT_END_NAMESPACE

#else

#include <gst/video/gstvideosink.h>

#include <QtCore/qlist.h>
#include <QtCore/qmutex.h>
#include <QtCore/qqueue.h>
#include <QtCore/qpointer.h>
#include <QtCore/qwaitcondition.h>
#include <qvideosurfaceformat.h>
#include <qvideoframe.h>
#include <qabstractvideobuffer.h>

#include "qgstbufferpoolinterface_p.h"

QT_BEGIN_NAMESPACE
class QAbstractVideoSurface;

class QVideoSurfaceGstDelegate : public QObject
{
    Q_OBJECT
public:
    QVideoSurfaceGstDelegate(QAbstractVideoSurface *surface);
    ~QVideoSurfaceGstDelegate();

    QList<QVideoFrame::PixelFormat> supportedPixelFormats(
            QAbstractVideoBuffer::HandleType handleType = QAbstractVideoBuffer::NoHandle) const;

    QVideoSurfaceFormat surfaceFormat() const;

    bool start(const QVideoSurfaceFormat &format, int bytesPerLine);
    void stop();

    void unlock();

    bool isActive();

    QGstBufferPoolInterface *pool() { return m_pool; }
    QMutex *poolMutex() { return &m_poolMutex; }
    void clearPoolBuffers();

    void flush();

    GstFlowReturn render(GstBuffer *buffer);

private slots:
    void queuedStart();
    void queuedStop();
    void queuedFlush();
    void queuedRender();

    void updateSupportedFormats();

private:
    QPointer<QAbstractVideoSurface> m_surface;
    QList<QVideoFrame::PixelFormat> m_supportedPixelFormats;
    //pixel formats of buffers pool native type
    QList<QVideoFrame::PixelFormat> m_supportedPoolPixelFormats;
    QGstBufferPoolInterface *m_pool = nullptr;
    QList<QGstBufferPoolInterface *> m_pools;
    QMutex m_poolMutex;
    QMutex m_mutex;
    QWaitCondition m_setupCondition;
    QWaitCondition m_renderCondition;
    QVideoSurfaceFormat m_format;
    QVideoFrame m_frame;
    GstFlowReturn m_renderReturn = GST_FLOW_ERROR;
    int m_bytesPerLine = 0;
    bool m_started = false;
    bool m_startCanceled = false;
};

class QVideoSurfaceGstSink
{
public:
    GstVideoSink parent;

    static QVideoSurfaceGstSink *createSink(QAbstractVideoSurface *surface);
    static void setSurface(QAbstractVideoSurface *surface) { Q_UNUSED(surface); }

private:
    static GType get_type();
    static void class_init(gpointer g_class, gpointer class_data);
    static void base_init(gpointer g_class);
    static void instance_init(GTypeInstance *instance, gpointer g_class);

    static void finalize(GObject *object);

    static void handleShowPrerollChange(GObject *o, GParamSpec *p, gpointer d);

    static GstStateChangeReturn change_state(GstElement *element, GstStateChange transition);

    static GstCaps *get_caps(GstBaseSink *sink);
    static gboolean set_caps(GstBaseSink *sink, GstCaps *caps);

    static GstFlowReturn buffer_alloc(
            GstBaseSink *sink, guint64 offset, guint size, GstCaps *caps, GstBuffer **buffer);

    static gboolean start(GstBaseSink *sink);
    static gboolean stop(GstBaseSink *sink);

    static gboolean unlock(GstBaseSink *sink);

#if GST_CHECK_VERSION(0, 10, 25)
    static GstFlowReturn show_frame(GstVideoSink *sink, GstBuffer *buffer);
#else
    static GstFlowReturn preroll(GstBaseSink *sink, GstBuffer *buffer);
    static GstFlowReturn render(GstBaseSink *sink, GstBuffer *buffer);
#endif

private:
    QVideoSurfaceGstDelegate *delegate = nullptr;

    GstCaps *lastRequestedCaps = nullptr;
    GstCaps *lastBufferCaps = nullptr;
    QVideoSurfaceFormat *lastSurfaceFormat = nullptr;
};

class QVideoSurfaceGstSinkClass
{
public:
    GstVideoSinkClass parent_class;
};

QT_END_NAMESPACE

#endif

#endif

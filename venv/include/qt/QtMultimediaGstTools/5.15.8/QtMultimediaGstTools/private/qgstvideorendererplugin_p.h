/****************************************************************************
**
** Copyright (C) 2016 Jolla Ltd.
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

#ifndef QGSTVIDEORENDERERPLUGIN_P_H
#define QGSTVIDEORENDERERPLUGIN_P_H

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
#include <qabstractvideobuffer.h>
#include <qvideosurfaceformat.h>
#include <QtCore/qobject.h>
#include <QtCore/qplugin.h>

#include <gst/gst.h>

QT_BEGIN_NAMESPACE

class QAbstractVideoSurface;

const QLatin1String QGstVideoRendererPluginKey("gstvideorenderer");

class Q_GSTTOOLS_EXPORT QGstVideoRenderer
{
public:
    virtual ~QGstVideoRenderer() {}

    virtual GstCaps *getCaps(QAbstractVideoSurface *surface) = 0;
    virtual bool start(QAbstractVideoSurface *surface, GstCaps *caps) = 0;
    virtual void stop(QAbstractVideoSurface *surface) = 0;  // surface may be null if unexpectedly deleted.
    virtual bool proposeAllocation(GstQuery *query) = 0;    // may be called from a thread.

    virtual bool present(QAbstractVideoSurface *surface, GstBuffer *buffer) = 0;
    virtual void flush(QAbstractVideoSurface *surface) = 0; // surface may be null if unexpectedly deleted.
};

/*
    Abstract interface for video buffers allocation.
*/
class Q_GSTTOOLS_EXPORT QGstVideoRendererInterface
{
public:
    virtual ~QGstVideoRendererInterface() {}

    virtual QGstVideoRenderer *createRenderer() = 0;
};

#define QGstVideoRendererInterface_iid "org.qt-project.qt.gstvideorenderer/5.4"
Q_DECLARE_INTERFACE(QGstVideoRendererInterface, QGstVideoRendererInterface_iid)

class Q_GSTTOOLS_EXPORT QGstVideoRendererPlugin : public QObject, public QGstVideoRendererInterface
{
    Q_OBJECT
    Q_INTERFACES(QGstVideoRendererInterface)
public:
    explicit QGstVideoRendererPlugin(QObject *parent = 0);
    virtual ~QGstVideoRendererPlugin() {}

    QGstVideoRenderer *createRenderer() override = 0;

};

QT_END_NAMESPACE

#endif

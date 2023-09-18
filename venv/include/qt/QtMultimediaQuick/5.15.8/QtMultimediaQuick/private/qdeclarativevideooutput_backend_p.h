/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Copyright (C) 2016 Research In Motion
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

#ifndef QDECLARATIVEVIDEOOUTPUT_BACKEND_P_H
#define QDECLARATIVEVIDEOOUTPUT_BACKEND_P_H

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

#include <QtCore/qpointer.h>
#include <QtCore/qsize.h>
#include <QtQuick/qquickitem.h>
#include <QtQuick/qsgnode.h>
#include <private/qtmultimediaquickdefs_p.h>

QT_BEGIN_NAMESPACE

class QAbstractVideoSurface;
class QDeclarativeVideoOutput;
class QMediaService;
class QAbstractVideoFilter;

class Q_MULTIMEDIAQUICK_EXPORT QDeclarativeVideoBackend
{
public:
    explicit QDeclarativeVideoBackend(QDeclarativeVideoOutput *parent)
        : q(parent)
    {}

    virtual ~QDeclarativeVideoBackend()
    {}

    virtual bool init(QMediaService *service) = 0;
    virtual void releaseSource() = 0;
    virtual void releaseControl() = 0;
    virtual void itemChange(QQuickItem::ItemChange change,
                            const QQuickItem::ItemChangeData &changeData) = 0;
    virtual QSize nativeSize() const = 0;
    virtual void updateGeometry() = 0;
    virtual QSGNode *updatePaintNode(QSGNode *oldNode, QQuickItem::UpdatePaintNodeData *data) = 0;
    virtual QAbstractVideoSurface *videoSurface() const = 0;

    // The viewport, adjusted for the pixel aspect ratio
    virtual QRectF adjustedViewport() const = 0;

    virtual void appendFilter(QAbstractVideoFilter *filter) { Q_UNUSED(filter); }
    virtual void clearFilters() { }

    virtual void releaseResources() { }
    virtual void invalidateSceneGraph() { }

protected:
    QDeclarativeVideoOutput *q;
    QPointer<QMediaService> m_service;
};

class QDeclarativeVideoBackendFactoryInterface
{
public:
    virtual QDeclarativeVideoBackend *create(QDeclarativeVideoOutput *parent) = 0;
};

#define QDeclarativeVideoBackendFactoryInterface_iid "org.qt-project.qt.declarativevideobackendfactory/5.2"
Q_DECLARE_INTERFACE(QDeclarativeVideoBackendFactoryInterface, QDeclarativeVideoBackendFactoryInterface_iid)

/*
 * Helper - returns true if the given orientation has the same aspect as the default (e.g. 180*n)
 */
namespace {

inline bool qIsDefaultAspect(int o)
{
    return (o % 180) == 0;
}

/*
 * Return the orientation normalized to 0-359
 */
inline int qNormalizedOrientation(int o)
{
    // Negative orientations give negative results
    int o2 = o % 360;
    if (o2 < 0)
        o2 += 360;
    return o2;
}

}

QT_END_NAMESPACE

#endif

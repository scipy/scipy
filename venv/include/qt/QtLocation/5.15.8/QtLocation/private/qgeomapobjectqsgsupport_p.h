/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtLocation module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QGEOMAPOBJECTQSGSUPPORT_P_H
#define QGEOMAPOBJECTQSGSUPPORT_P_H

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

#include <QtLocation/private/qlocationglobal_p.h>
#include <QtLocation/private/qgeomapobject_p.h>
#include <QtLocation/private/qgeomapobject_p_p.h>
#include <QtLocation/private/qmappolylineobjectqsg_p_p.h>
#include <QtLocation/private/qmappolygonobjectqsg_p_p.h>
#include <QtLocation/private/qmapcircleobjectqsg_p_p.h>
#include <QtLocation/private/qmaprouteobjectqsg_p_p.h>
#include <QtLocation/private/qmapiconobjectqsg_p_p.h>
#include <QtLocation/private/qdeclarativepolylinemapitem_p.h>
#include <QtLocation/private/qdeclarativepolygonmapitem_p_p.h>
#include <QtCore/qpointer.h>
#include <memory>

QT_BEGIN_NAMESPACE
struct Q_LOCATION_PRIVATE_EXPORT MapObject {
    MapObject(QPointer<QGeoMapObject> &o, QQSGMapObject *sgo)
        : object(o), sgObject(sgo) {}

    QPointer<QGeoMapObject> object;
    QQSGMapObject *sgObject = nullptr; // this is a QMap*ObjectPrivateQSG. it becomes invalid when the pimpl is destroyed
    VisibleNode *visibleNode = nullptr; // This is a Map*Node (like a MapPolygonNode) that is a QSGNode. This doesn't disappear by itself
    QSGNode *qsgNode = nullptr;
};

class Q_LOCATION_PRIVATE_EXPORT QGeoMapObjectQSGSupport
{
public:
    bool createMapObjectImplementation(QGeoMapObject *obj, QGeoMapPrivate *d);
    QGeoMapObjectPrivate *createMapObjectImplementationPrivate(QGeoMapObject *obj);
    QList<QGeoMapObject *> mapObjects() const;
    void removeMapObject(QGeoMapObject *obj);
    void updateMapObjects(QSGNode *root, QQuickWindow *window);
    void updateObjectsGeometry();

    QList<MapObject> m_mapObjects;
    QList<MapObject> m_pendingMapObjects;
    QList<MapObject> m_removedMapObjects;
    QGeoMap *m_map = nullptr;
    std::unique_ptr<QDeclarativePolygonMapItemPrivateOpenGL::RootNode> m_mapObjectsRootNode;
};

QT_END_NAMESPACE

#endif // QGEOMAPOBJECTQSGSUPPORT_P_H

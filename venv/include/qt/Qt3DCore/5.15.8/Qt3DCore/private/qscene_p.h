/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DCORE_QSCENE_P_H
#define QT3DCORE_QSCENE_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtCore/QScopedPointer>

#include <Qt3DCore/qnode.h>
#include <QtCore/qscopedpointer.h>
#include <Qt3DCore/private/qobservableinterface_p.h>
#include <Qt3DCore/private/qt3dcore_global_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QScenePrivate;
class QAspectEngine;
class NodePostConstructorInit;
class QAbstractArbiter;

typedef QList<QObservableInterface *> QObservableList;

class Q_3DCORE_PRIVATE_EXPORT QScene
{
public:
    QScene(QAspectEngine *engine = nullptr);
    ~QScene();

    QAspectEngine *engine() const;

    void addObservable(QObservableInterface *observable, QNodeId id);
    void addObservable(QNode *observable);
    void removeObservable(QObservableInterface *observable, QNodeId id);
    void removeObservable(QNode *observable);
    QObservableList lookupObservables(QNodeId id) const;

    QNode *lookupNode(QNodeId id) const;
    QVector<QNode *> lookupNodes(const QVector<QNodeId> &ids) const;
    QNodeId nodeIdFromObservable(QObservableInterface *observable) const;

    QNode *rootNode() const;

    void setArbiter(QAbstractArbiter *arbiter);
    QAbstractArbiter *arbiter() const;

    // Component -> Entities
    QVector<QNodeId> entitiesForComponent(QNodeId id) const;
    void addEntityForComponent(QNodeId componentUuid, QNodeId entityUuid);
    void removeEntityForComponent(QNodeId componentUuid, QNodeId entityUuid);
    bool hasEntityForComponent(QNodeId componentUuid, QNodeId entityUuid);

    // Node -> Property Update Data
    struct NodePropertyTrackData
    {
        QNode::PropertyTrackingMode defaultTrackMode = QNode::TrackFinalValues;
        QHash<QString, QNode::PropertyTrackingMode> trackedPropertiesOverrides;
    };
    NodePropertyTrackData lookupNodePropertyTrackData(QNodeId id) const;
    void setPropertyTrackDataForNode(QNodeId id, const NodePropertyTrackData &data);
    void removePropertyTrackDataForNode(QNodeId id);

    NodePostConstructorInit* postConstructorInit() const;

private:
    Q_DECLARE_PRIVATE(QScene)
    QScopedPointer<QScenePrivate> d_ptr;

    void setRootNode(QNode *root);
    friend class QAspectEnginePrivate;
};

} // Qt3D

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DCore::QScene*) // LCOV_EXCL_LINE

#endif // QT3DCORE_QSCENE_P_H

/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QSSGSCENEMANAGER_P_H
#define QSSGSCENEMANAGER_P_H

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

#include <QtCore/QObject>
#include <QtCore/QSet>

#include <QtQuick3D/private/qtquick3dglobal_p.h>

#include "qquick3dobject_p.h"
#include "qquick3dnode_p.h"

QT_BEGIN_NAMESPACE

class QSGDynamicTexture;
class QQuickWindow;
class QSSGBufferManager;

class Q_QUICK3D_PRIVATE_EXPORT QQuick3DSceneManager : public QObject
{
    Q_OBJECT
public:
    explicit QQuick3DSceneManager(QObject *parent = nullptr);
    ~QQuick3DSceneManager() override;

    void setWindow(QQuickWindow *window);
    QQuickWindow *window();

    void dirtyItem(QQuick3DObject *item);
    void cleanup(QSSGRenderGraphObject *item);

    void polishItems();
    void forcePolish();
    void sync();

    void updateDirtyNodes();
    void updateDirtyNode(QQuick3DObject *object);
    void updateDirtyResource(QQuick3DObject *resourceObject);
    void updateDirtySpatialNode(QQuick3DNode *spatialNode);
    void updateBoundingBoxes(const QSSGRef<QSSGBufferManager> &mgr);

    QQuick3DObject *lookUpNode(const QSSGRenderGraphObject *node) const;

    void cleanupNodes();

    QQuick3DObject *dirtySpatialNodeList;
    QQuick3DObject *dirtyResourceList;
    QQuick3DObject *dirtyImageList;
    QList<QQuick3DObject *> dirtyLightList;
    QList<QQuick3DObject *> dirtyBoundingBoxList;
    QList<QSSGRenderGraphObject *> cleanupNodeList;
    QSet<QQuick3DObject *> parentlessItems;
    QVector<QSGDynamicTexture *> qsgDynamicTextures;
    QHash<const QSSGRenderGraphObject *, QQuick3DObject *> m_nodeMap;
    QQuickWindow *m_window = nullptr;
    friend QQuick3DObject;

Q_SIGNALS:
    void needsUpdate();
    void windowChanged();
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuick3DSceneManager)

#endif // QSSGSCENEMANAGER_P_H

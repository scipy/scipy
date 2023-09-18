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

#ifndef QQUICK3DNODE_P_P_H
#define QQUICK3DNODE_P_P_H

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


#include <QtQuick3D/private/qtquick3dglobal_p.h>

#include "qquick3dobject_p.h"
#include "qquick3dnode_p.h"

#include <QtGui/QVector3D>
#include <QtGui/QQuaternion>
#include <QtGui/QMatrix4x4>

QT_BEGIN_NAMESPACE

class QQuick3DNode;

class Q_QUICK3D_PRIVATE_EXPORT QQuick3DNodePrivate : public QQuick3DObjectPrivate
{
    Q_DECLARE_PUBLIC(QQuick3DNode)

public:

    explicit QQuick3DNodePrivate(QQuick3DObjectPrivate::Type t);
    ~QQuick3DNodePrivate();
    void init();

    QMatrix4x4 calculateLocalTransform();
    void calculateGlobalVariables();
    void markSceneTransformDirty();

    inline QMatrix4x4 localRotationMatrix() const;
    inline QMatrix4x4 sceneRotationMatrix() const;

    void emitChangesToSceneTransform();
    bool isSceneTransformRelatedSignal(const QMetaMethod &signal) const;

    void setIsHiddenInEditor(bool isHidden);

    static inline QQuick3DNodePrivate *get(QQuick3DNode *node) { return node->d_func(); }

    QQuaternion m_rotation;
    QVector3D m_eulerRotationAngles;
    QVector3D m_position;
    QVector3D m_scale{ 1.0f, 1.0f, 1.0f };
    QVector3D m_pivot;
    float m_opacity = 1.0f;
    int m_staticFlags = 0;
    bool m_visible = true;
    QMatrix4x4 m_sceneTransform; // Right handed
    bool m_sceneTransformDirty = true;
    int m_sceneTransformConnectionCount = 0;
    bool m_isHiddenInEditor = false;
    bool m_hasInheritedUniformScale = true;
    bool m_eulerRotationDirty = false;
};


QT_END_NAMESPACE

#endif // QQUICK3DNODE_P_P_H



/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
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

#ifndef QSSGMESHBVH_H
#define QSSGMESHBVH_H

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

#include <QtQuick3DUtils/private/qtquick3dutilsglobal_p.h>
#include <QtQuick3DUtils/private/qssgbounds3_p.h>

#include <QtGui/QVector2D>
#include <QtCore/QVector>

QT_BEGIN_NAMESPACE

struct Q_QUICK3DUTILS_EXPORT QSSGMeshBVHNode {
    ~QSSGMeshBVHNode() {
        delete left;
        delete right;
    }

    // Internal
    QSSGMeshBVHNode *left = nullptr;
    QSSGMeshBVHNode *right = nullptr;
    QSSGBounds3 boundingData;
    //splitAxis

    // Leaf
    int offset = 0;
    int count = 0;
};

struct Q_QUICK3DUTILS_EXPORT QSSGMeshBVHTriangle {
    QSSGBounds3 bounds;
    QVector3D vertex1;
    QVector3D vertex2;
    QVector3D vertex3;
    QVector2D uvCoord1;
    QVector2D uvCoord2;
    QVector2D uvCoord3;
};

struct Q_QUICK3DUTILS_EXPORT QSSGMeshBVH
{
    QSSGMeshBVH(const QVector<QSSGMeshBVHNode *> &bvhRoots,
                const QVector<QSSGMeshBVHTriangle *> &bvhTriangles)
        : roots(bvhRoots)
        , triangles(bvhTriangles)
    {}
    ~QSSGMeshBVH();

    QVector<QSSGMeshBVHNode *> roots;
    QVector<QSSGMeshBVHTriangle *> triangles;
};

QT_END_NAMESPACE

#endif // QSSGMESHBVH_H

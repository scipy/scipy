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

#ifndef QSSGMESHBVHBUILDER_H
#define QSSGMESHBVHBUILDER_H

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
#include <QtQuick3DAssetImport/private/qtquick3dassetimportglobal_p.h>
#include <QtQuick3DUtils/private/qssgmeshbvh_p.h>
#include <QtQuick3DAssetImport/private/qssgmeshutilities_p.h>

QT_BEGIN_NAMESPACE

class Q_QUICK3DASSETIMPORT_EXPORT QSSGMeshBVHBuilder
{
public:
    QSSGMeshBVHBuilder(QSSGMeshUtilities::Mesh *mesh);

    QSSGMeshBVH* buildTree();

private:
    enum class Axis
    {
        None = -1,
        X = 0,
        Y = 1,
        Z = 2
    };
    struct Split {
        Axis axis;
        float pos;
    };

    QVector<QSSGMeshBVHTriangle*> calculateTriangleBounds(quint32 indexOffset, quint32 indexCount) const;
    quint32 getIndexBufferValue(quint32 index) const;
    QVector3D getVertexBufferValuePosition(quint32 index) const;
    QVector2D getVertexBufferValueUV0(quint32 index) const;

    QSSGMeshBVHNode *splitNode(QSSGMeshBVHNode *node, quint32 offset, quint32 count, quint32 depth = 0);
    QSSGBounds3 getBounds(quint32 offset, quint32 count) const;
    Split getOptimalSplit(const QSSGBounds3 &nodeBounds, quint32 offset, quint32 count) const;
    static Axis getLongestDimension(const QSSGBounds3 &nodeBounds);
    float getAverageValue(quint32 offset, quint32 count, Axis axis) const;
    quint32 partition(quint32 offset, quint32 count, const Split &split);

    QSSGMeshUtilities::Mesh *m_mesh;
    quint8 *m_baseAddress;
    QSSGRenderComponentType m_indexBufferComponentType;
    QSSGByteView m_indexBufferData;
    QSSGByteView m_vertexBufferData;
    quint32 m_vertexStride;
    bool m_hasPositionData = false;
    quint32 m_vertexPosOffset;
    bool m_hasUVData = false;
    quint32 m_vertexUV0Offset;


    QVector<QSSGMeshBVHTriangle *> m_triangleBounds;
    QVector<QSSGMeshBVHNode *> m_roots;
    quint32 m_maxTreeDepth = 40;
    quint32 m_maxLeafTriangles = 10;
};

QT_END_NAMESPACE

#endif // QSSGMESHBVHBUILDER_H

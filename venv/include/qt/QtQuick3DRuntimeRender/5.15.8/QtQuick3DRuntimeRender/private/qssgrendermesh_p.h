/****************************************************************************
**
** Copyright (C) 2008-2012 NVIDIA Corporation.
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

#ifndef QSSG_RENDER_MESH_H
#define QSSG_RENDER_MESH_H

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

#include <QtQuick3DRender/private/qssgrendervertexbuffer_p.h>
#include <QtQuick3DRender/private/qssgrenderindexbuffer_p.h>
#include <QtQuick3DRender/private/qssgrenderinputassembler_p.h>

#include <QtQuick3DUtils/private/qssgbounds3_p.h>
#include <QtQuick3DUtils/private/qssgmeshbvh_p.h>

QT_BEGIN_NAMESPACE

struct QSSGRenderSubsetBase
{
    quint32 count;
    quint32 offset;
    QSSGBounds3 bounds; // Vertex buffer bounds
    QSSGMeshBVHNode *bvhRoot = nullptr;
    QSSGRenderSubsetBase() = default;
    QSSGRenderSubsetBase(const QSSGRenderSubsetBase &inOther)
        : count(inOther.count)
        , offset(inOther.offset)
        , bounds(inOther.bounds)
        , bvhRoot(inOther.bvhRoot)
    {
    }

    QSSGRenderSubsetBase &operator=(const QSSGRenderSubsetBase &inOther)
    {
        count = inOther.count;
        offset = inOther.offset;
        bounds = inOther.bounds;
        bvhRoot = inOther.bvhRoot;
        return *this;
    }
};

struct QSSGRenderJoint
{
    qint32 jointID;
    qint32 parentID;
    float invBindPose[16];
    float localToGlobalBoneSpace[16];
};

struct QSSGRenderSubset : public QSSGRenderSubsetBase
{
    QSSGRef<QSSGRenderInputAssembler> inputAssembler;
    QSSGRef<QSSGRenderInputAssembler> inputAssemblerDepth;
    QSSGRef<QSSGRenderInputAssembler> inputAssemblerPoints; ///< similar to depth but ignores index buffer.
    QSSGRef<QSSGRenderVertexBuffer> vertexBuffer;
    QSSGRef<QSSGRenderVertexBuffer> posVertexBuffer; ///< separate position buffer for fast depth path rendering
    QSSGRef<QSSGRenderIndexBuffer> indexBuffer;
    QSSGRenderDrawMode primitiveType; ///< primitive type used for drawing
    float edgeTessFactor = 1.0f; ///< edge tessellation amount used for tessellation shaders
    float innerTessFactor = 1.0f; ///< inner tessellation amount used for tessellation shaders
    bool wireframeMode; ///< true if we should draw the object as wireframe ( currently ony if
    /// tessellation is enabled )
    QVector<QSSGRenderJoint> joints;
    QString name;
    QVector<QSSGRenderSubsetBase> subSubsets;

    QSSGRenderSubset() = default;
    QSSGRenderSubset(const QSSGRenderSubset &inOther)
        : QSSGRenderSubsetBase(inOther)
        , inputAssembler(inOther.inputAssembler)
        , inputAssemblerDepth(inOther.inputAssemblerDepth)
        , inputAssemblerPoints(inOther.inputAssemblerPoints)
        , vertexBuffer(inOther.vertexBuffer)
        , posVertexBuffer(inOther.posVertexBuffer)
        , indexBuffer(inOther.indexBuffer)
        , primitiveType(inOther.primitiveType)
        , edgeTessFactor(inOther.edgeTessFactor)
        , innerTessFactor(inOther.innerTessFactor)
        , wireframeMode(inOther.wireframeMode)
        , joints(inOther.joints)
        , name(inOther.name)
        , subSubsets(inOther.subSubsets)
    {
    }
    // Note that subSubsets is *not* copied.
    QSSGRenderSubset(const QSSGRenderSubset &inOther, const QSSGRenderSubsetBase &inBase)
        : QSSGRenderSubsetBase(inBase)
        , inputAssembler(inOther.inputAssembler)
        , inputAssemblerDepth(inOther.inputAssemblerDepth)
        , inputAssemblerPoints(inOther.inputAssemblerPoints)
        , vertexBuffer(inOther.vertexBuffer)
        , posVertexBuffer(inOther.posVertexBuffer)
        , indexBuffer(inOther.indexBuffer)
        , primitiveType(inOther.primitiveType)
        , edgeTessFactor(inOther.edgeTessFactor)
        , innerTessFactor(inOther.innerTessFactor)
        , wireframeMode(inOther.wireframeMode)
        , name(inOther.name)
    {
    }

    QSSGRenderSubset &operator=(const QSSGRenderSubset &inOther)
    {
        if (this != &inOther) {
            QSSGRenderSubsetBase::operator=(inOther);
            inputAssembler = inOther.inputAssembler;
            inputAssemblerDepth = inOther.inputAssemblerDepth;
            vertexBuffer = inOther.vertexBuffer;
            posVertexBuffer = inOther.posVertexBuffer;
            indexBuffer = inOther.indexBuffer;
            primitiveType = inOther.primitiveType;
            edgeTessFactor = inOther.edgeTessFactor;
            innerTessFactor = inOther.innerTessFactor;
            wireframeMode = inOther.wireframeMode;
            joints = inOther.joints;
            name = inOther.name;
            subSubsets = inOther.subSubsets;
        }
        return *this;
    }
};

struct QSSGRenderMeshPath
{
    QString path;
    uint key = 0;

    inline bool isNull() const { return path.isNull(); }

    static QSSGRenderMeshPath create(const QString &path)
    {
        QSSGRenderMeshPath p;
        p.path = path;
        p.key = qHash(path);
        return p;
    }
};

inline bool operator==(const QSSGRenderMeshPath &p1, const QSSGRenderMeshPath &p2)
{
    return (p1.path == p2.path);
}

inline uint qHash(const QSSGRenderMeshPath &path, uint seed) Q_DECL_NOTHROW
{
    return (path.key) ? path.key : qHash(path.path, seed);
}

struct QSSGRenderMesh
{
    Q_DISABLE_COPY(QSSGRenderMesh)

    QVector<QSSGRenderSubset> subsets;
    QVector<QSSGRenderJoint> joints;
    QSSGRenderDrawMode drawMode;
    QSSGRenderWinding winding; // counterclockwise
    quint32 meshId; // Id from the file of this mesh.
    QSSGMeshBVH *bvh = nullptr;
    QVector<QByteArray> inputLayoutInputNames;

    QSSGRenderMesh(QSSGRenderDrawMode inDrawMode, QSSGRenderWinding inWinding, quint32 inMeshId)
        : drawMode(inDrawMode), winding(inWinding), meshId(inMeshId)
    {
    }

    ~QSSGRenderMesh()
    {
        delete bvh;
    }
};
QT_END_NAMESPACE

#endif

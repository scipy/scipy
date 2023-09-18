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

#ifndef QSSG_RENDER_NODE_H
#define QSSG_RENDER_NODE_H

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

#include <QtQuick3DRuntimeRender/private/qssgrendergraphobject_p.h>

#include <QtQuick3DUtils/private/qssgbounds3_p.h>

#include <QtGui/QMatrix4x4>
#include <QtGui/QVector3D>

QT_BEGIN_NAMESPACE

struct QSSGRenderModel;
struct QSSGRenderLight;
struct QSSGRenderCamera;
struct QSSGRenderText;
struct QSSGRenderNode;
class QSSGBufferManager;

class QSSGRenderNodeFilterInterface;

struct Q_QUICK3DRUNTIMERENDER_EXPORT QSSGRenderNode : public QSSGRenderGraphObject
{
    enum class Flag
    {
        Dirty = 1,
        TransformDirty = 1 << 1,
        Active = 1 << 2, ///< Is this exact object active
        Orthographic = 1 << 3,
        PointLight = 1 << 4,
        GloballyActive = 1 << 5, ///< set based in Active and if a parent is active.
        TextDirty = 1 << 6,
        LocallyPickable = 1 << 7,
        GloballyPickable = 1 << 8,
        LayerEnableDepthTest = 1 << 9,
        LayerRenderToTarget = 1 << 10, ///< Does this layer render to the normal render target,
        /// or is it offscreen-only
        ForceLayerOffscreen = 1 << 11, ///< Forces a layer to always use the offscreen rendering
        /// mechanism.  This can be usefulf or caching purposes.
        IgnoreParentTransform = 1 << 12,
        LayerEnableDepthPrePass = 1 << 13, ///< True when we render a depth pass before
        CameraDirty = 1 << 14, ///< True when the camera inheriting from this is dirty
        CameraFrustumProjection = 1 << 15,
        CameraCustomProjection = 1 << 16
    };
    Q_DECLARE_FLAGS(Flags, Flag)

    enum class TransformDirtyFlag : quint8
    {
        TransformNotDirty,
        TransformIsDirty,
    };

    // changing any one of these means you have to
    // set this object dirty
    QQuaternion rotation;
    QVector3D position { 0.0f, 0.0f, 0.0f };
    QVector3D scale { 1.0f, 1.0f, 1.0f };
    QVector3D pivot { 0.0f, 0.0f, 0.0f };
    int staticFlags = 0;

    // This only sets dirty, not transform dirty
    // Opacity of 1 means opaque, opacity of zero means transparent.
    float localOpacity = 1.0f;

    // results of clearing dirty.
    Flags flags {
        Flag::Dirty,
        Flag::TransformDirty,
        Flag::Active,
    };
    // These end up right handed
    QMatrix4x4 localTransform;
    QMatrix4x4 globalTransform;
    float globalOpacity = 1.0f;
    qint32 skeletonId = -1;

    // node graph members.
    QSSGRenderNode *parent = nullptr;
    QSSGRenderNode *nextSibling = nullptr;
    QSSGRenderNode *previousSibling = nullptr;
    QSSGRenderNode *firstChild = nullptr;
    // Property maintained solely by the render system.
    // Depth-first-search index assigned and maintained by render system.
    quint32 dfsIndex = 0;

    QSSGRenderNode();
    QSSGRenderNode(Type type);
    QSSGRenderNode(const QSSGRenderNode &inCloningObject);
    ~QSSGRenderNode() {}

    // Sets this object dirty and walks down the graph setting all
    // children who are not dirty to be dirty.
    void markDirty(TransformDirtyFlag inTransformDirty = TransformDirtyFlag::TransformNotDirty);

    void addChild(QSSGRenderNode &inChild);
    // Specital function for importScene
    void addChildrenToLayer(QSSGRenderNode &inChildren);
    void removeChild(QSSGRenderNode &inChild);
    QSSGRenderNode *getLastChild();

    // Remove this node from the graph.
    // It is no longer the the parent's child lists
    // and all of its children no longer have a parent
    // finally they are no longer siblings of each other.
    void removeFromGraph();

    // Calculate global transform and opacity
    // Walks up the graph ensure all parents are not dirty so they have
    // valid global transforms.
    bool calculateGlobalVariables();

    // Given our rotation order and handedness, calculate the final rotation matrix
    // Only the upper 3x3 of this matrix is filled in.
    // If this object is left handed, then you need to call FlipCoordinateSystem
    // to get a result identical to the result produced in CalculateLocalTransform
    void calculateRotationMatrix(QMatrix4x4 &outMatrix) const;

    // Force the calculation of the local transform
    void calculateLocalTransform();

    /**
     * @brief setup local tranform from a matrix.
     *		  This function decomposes a SRT matrix.
     *		  This will fail if this matrix contains non-affine transformations
     *
     * @param inTransform[in]	input transformation
     *
     * @return true backend type
     */
    void setLocalTransformFromMatrix(QMatrix4x4 &inTransform);

    // Get the bounds of us and our children in our local space.
    QSSGBounds3 getBounds(const QSSGRef<QSSGBufferManager> &inManager,
                            bool inIncludeChildren = true,
                            QSSGRenderNodeFilterInterface *inChildFilter = nullptr) const;
    QSSGBounds3 getChildBounds(const QSSGRef<QSSGBufferManager> &inManager,
                                 QSSGRenderNodeFilterInterface *inChildFilter = nullptr) const;
    // Assumes CalculateGlobalVariables has already been called.
    QVector3D getGlobalPos() const;
    QVector3D getGlobalPivot() const;
    // Pulls the 3rd column out of the global transform.
    QVector3D getDirection() const;
    // Multiplies (0,0,-1) by the inverse transpose of the upper 3x3 of the global transform.
    // This is correct w/r/t to scaling and which the above getDirection is not.
    QVector3D getScalingCorrectDirection() const;

    // outMVP and outNormalMatrix are returned ready to upload to openGL, meaning they are
    // row-major.
    void calculateMVPAndNormalMatrix(const QMatrix4x4 &inViewProjection, QMatrix4x4 &outMVP, QMatrix3x3 &outNormalMatrix) const;

    // This should be in a utility file somewhere
    QMatrix3x3 calculateNormalMatrix() const;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QSSGRenderNode::Flags)

QT_END_NAMESPACE

#endif

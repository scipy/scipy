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

#ifndef QSSG_RENDERER_H
#define QSSG_RENDERER_H

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

#include <QtQuick3DUtils/private/qssgdataref_p.h>
#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendergraphobjectpickquery_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendercamera_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderray_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendernode_p.h>

#include <QtGui/QVector2D>

QT_BEGIN_NAMESPACE

typedef void *QSSGRenderInstanceId;

class QSSGRenderNodeFilterInterface
{
protected:
    virtual ~QSSGRenderNodeFilterInterface() {}

public:
    virtual bool includeNode(const QSSGRenderNode &inNode) = 0;
};
struct QSSGLayerPickSetup
{
    QMatrix4x4 projectionPreMultiply;
    QMatrix4x4 viewProjection;
    QRect scissorRect;
    QSSGLayerPickSetup(const QMatrix4x4 &inProjPreMult, const QMatrix4x4 &inVP, const QRect &inScissor)
        : projectionPreMultiply(inProjPreMult), viewProjection(inVP), scissorRect(inScissor)
    {
    }
    QSSGLayerPickSetup() {}
};

struct QSSGScaleAndPosition
{
    QVector3D position;
    float scale;
    QSSGScaleAndPosition(const QVector3D &inPos, float inScale) : position(inPos), scale(inScale) {}
    QSSGScaleAndPosition() = default;
};

struct QSSGRenderLayer;
class QSSGRendererImpl;
class QSSGRenderContextInterface;

class Q_QUICK3DRUNTIMERENDER_EXPORT QSSGRendererInterface
{
public:
    QAtomicInt ref;
    virtual ~QSSGRendererInterface() {}
    virtual void enableLayerGpuProfiling(bool inEnabled) = 0;
    virtual bool isLayerGpuProfilingEnabled() const = 0;

    // Get the camera that rendered this node last render
    virtual QSSGRenderCamera *cameraForNode(const QSSGRenderNode &inNode) const = 0;
    virtual QSSGOption<QSSGCuboidRect> cameraBounds(const QSSGRenderGraphObject &inObject) = 0;
    // Called when you have changed the number or order of children of a given node.
    virtual void childrenUpdated(QSSGRenderNode &inParent) = 0;

    // The QSSGRenderContextInterface calls these, clients should not.
    virtual void beginFrame() = 0;
    virtual void endFrame() = 0;

    // Setup the vertex and index buffers (but not shader state)
    // and render the quad.  The quad is setup so that its edges
    // go from -1,1 in x,y and its UV coordinates will map naturally
    // to an image.
    virtual void renderQuad() = 0;

    // Render a given texture as flipped to the scene using a given transform.
    virtual void renderFlippedQuad(const QVector2D &inDimensions,
                                   const QMatrix4x4 &inMVP,
                                   QSSGRenderTexture2D &inQuadTexture,
                                   float opacity) = 0;

    // Returns true if this layer or a sibling was dirty.
    virtual bool prepareLayerForRender(QSSGRenderLayer &inLayer, const QSize &surfaceSize) = 0;
    virtual void renderLayer(QSSGRenderLayer &inLayer,
                             const QSize &surfaceSize,
                             bool clear,
                             const QColor &clearColor) = 0;

    // Studio option to disable picking against sub renderers.  This allows better interaction
    // in studio.
    // In pick siblings measn pick the layer siblings; this is the normal behavior.
    // InPickEverything means ignore the node's pick flags; this allows us to only pick things
    // that have handlers
    // in some cases and just pick everything in other things.
    virtual void pickRenderPlugins(bool inPick) = 0;
    virtual QSSGRenderPickResult pick(QSSGRenderLayer &inLayer,
                                        const QVector2D &inViewportDimensions,
                                        const QVector2D &inMouseCoords,
                                        bool inPickSiblings = true,
                                        bool inPickEverything = false) = 0;
    virtual QSSGRenderPickResult syncPick(const QSSGRenderLayer &inLayer,
                                          const QSSGRef<QSSGBufferManager> &bufferManager,
                                          const QVector2D &inViewportDimensions,
                                          const QVector2D &inMouseCoords) = 0;

    // Return the relative hit position, in UV space, of a mouse pick against this object.
    // We need the node in order to figure out which layer rendered this object.
    // We need mapper objects if this is a in a subpresentation because we have to know how
    // to map the mouse coordinates into the subpresentation.  So for instance if inNode is in
    // a subpres then we need to know which image is displaying the subpres in order to map
    // the mouse coordinates into the subpres's render space.
    virtual QSSGOption<QVector2D> facePosition(QSSGRenderNode &inNode,
                                                 QSSGBounds3 inBounds,
                                                 const QMatrix4x4 &inGlobalTransform,
                                                 const QVector2D &inViewportDimensions,
                                                 const QVector2D &inMouseCoords,
                                                 QSSGDataView<QSSGRenderGraphObject *> inMapperObjects,
                                                 QSSGRenderBasisPlanes inIsectPlane) = 0;

    virtual QVector3D unprojectToPosition(QSSGRenderNode &inNode, QVector3D &inPosition, const QVector2D &inMouseVec) const = 0;
    virtual QVector3D unprojectWithDepth(QSSGRenderNode &inNode, QVector3D &inPosition, const QVector3D &inMouseVec) const = 0;
    virtual QVector3D projectPosition(QSSGRenderNode &inNode, const QVector3D &inPosition) const = 0;

    // Roughly equivalent of gluPickMatrix, allows users to setup a perspective transform that
    // will draw some sub component
    // of the layer.  Used in combination with an expected viewport of 0,0,width,height the
    // viewproj matrix returned will center
    // around the center of the viewport and render just the part of the layer around this area.
    // The return value is optional because if the mouse point is completely outside the layer
    // obviously this method is irrelevant.
    virtual QSSGOption<QSSGLayerPickSetup> getLayerPickSetup(QSSGRenderLayer &inLayer,
                                                                 const QVector2D &inMouseCoords,
                                                                 const QSize &inPickDims) = 0;

    // Return the layer's viewport rect after the layer's member variables have been applied.
    // Uses the last rendered viewport rect.
    virtual QSSGOption<QRectF> layerRect(QSSGRenderLayer &inLayer) = 0;
    // Testing function to allow clients to render a layer using a custom view project instead
    // of the one that would be setup
    // using the layer's camera in conjunction with the layer's position,scale.
    virtual void runLayerRender(QSSGRenderLayer &inLayer, const QMatrix4x4 &inViewProjection) = 0;

    // Render the layer's rect onscreen.  Will only render one frame, you need to call this
    // every frame
    // for this to work and be persistent.
    virtual void renderLayerRect(QSSGRenderLayer &inLayer, const QVector3D &inColor) = 0;

    // Called before a layer goes completely out of scope to release any rendering resources
    // related to the layer.
    virtual void releaseLayerRenderResources(QSSGRenderLayer &inLayer) = 0;

    // render Gpu profiler values
    virtual void dumpGpuProfilerStats() = 0;

    // Get the mouse coordinates as they relate to a given layer
    virtual QSSGOption<QVector2D> getLayerMouseCoords(QSSGRenderLayer &inLayer,
                                                        const QVector2D &inMouseCoords,
                                                        const QVector2D &inViewportDimensions,
                                                        bool forceImageIntersect = false) const = 0;

    // Returns true if the renderer expects new frame to be rendered
    // Happens when progressive AA is enabled
    virtual bool rendererRequestsFrames() const = 0;

    static bool isGlEsContext(const QSSGRenderContextType &inContextType);
    static bool isGlEs3Context(const QSSGRenderContextType &inContextType);
    static bool isGl2Context(const QSSGRenderContextType &inContextType);
    static const char *getGlslVesionString(QSSGRenderContextType inContextType);

    static QSSGRef<QSSGRendererInterface> createRenderer(QSSGRenderContextInterface *inContext);
};
QT_END_NAMESPACE

#endif

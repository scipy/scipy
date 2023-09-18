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

#ifndef QSSG_RENDER_LAYER_HELPER_IMPL_H
#define QSSG_RENDER_LAYER_HELPER_IMPL_H

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

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendercamera_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendercontextcore_p.h>

QT_BEGIN_NAMESPACE

/**	An independent, testable entity to encapsulate taking at least:
 *  layer, current viewport rect, current scissor rect, presentation design dimensions
 *	and producing a set of rectangles:
 *	layer viewport rect (inside viewport rect and calculated using outer viewport rect info)
 *	layer scissor rect (inside current scissor rect)
 *	layer camera rect (may be the viewport rect)
 *
 *  In the case where we have to render offscreen for this layer then we need to handle produce
 *	a set of texture dimensions and the layer camera rect ends up being same size but with no
 *offsets.
 *
 *  This object should handle part of any translation from screenspace to global space.
 *	I am using language level access control on this object because it needs specific
 *	interface design that will enable future modifications.
 */
struct Q_AUTOTEST_EXPORT QSSGLayerRenderHelper
{
private:
    QSSGRenderLayer *m_layer = nullptr;
    QSSGRenderCamera *m_camera = nullptr;

    QRectF m_viewport;
    QRectF m_scissor;

public:
    QSSGLayerRenderHelper() = default;

    QSSGLayerRenderHelper(const QRectF &inViewport,
                            const QRectF &inScissor,
                            QSSGRenderLayer &inLayer);

    QSSGRenderLayer *layer() const { return m_layer; }
    QSSGRenderCamera *camera() const { return m_camera; }

    // Does not differ whether offscreen or not, simply states how this layer maps to the
    // presentation
    QRectF viewport() const { return m_viewport; }
    // Does not differ whether offscreen or not, scissor rect of how this layer maps to
    // presentation.
    QRectF scissor() const { return m_scissor; }

    QSize textureDimensions() const;

    QSSGCameraGlobalCalculationResult setupCameraForRender(QSSGRenderCamera &inCamera);

    static QSSGOption<QVector2D> layerMouseCoords(const QRectF &viewport, const QVector2D &inMouseCoords, const QVector2D &inWindowDimensions, bool inForceIntersect);
    static QSSGOption<QSSGRenderRay> pickRay(const QSSGRenderCamera &camera, const QRectF &viewport, const QVector2D &inMouseCoords, const QVector2D &inWindowDimensions, bool inForceIntersect);

    // Checks the various viewports and determines if the layer is visible or not.
    bool isLayerVisible() const;

private:
    // Viewport used when actually rendering.  In the case where this is an offscreen item then
    // it may be different than the layer to presentation viewport.
    QRectF layerRenderViewport() const;
};
QT_END_NAMESPACE

#endif // QSSG_RENDER_LAYER_HELPER_IMPL_H

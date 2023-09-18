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

#ifndef QSSG_RENDER_LAYER_H
#define QSSG_RENDER_LAYER_H

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

#include <QtQuick3DRuntimeRender/private/qssgrendernode_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderer_p.h>

QT_BEGIN_NAMESPACE
class QSSGRenderContextInterface;
struct QSSGRenderPresentation;
struct QSSGRenderEffect;
struct QSSGRenderImage;

// A layer is a special node.  It *always* presents its global transform
// to children as the identity.  It also can optionally have a width or height
// different than the overlying context.  You can think of layers as the transformation
// between a 3d scene graph and a 2D texture.
struct Q_QUICK3DRUNTIMERENDER_EXPORT QSSGRenderLayer : public QSSGRenderNode
{
    enum class AAMode : quint8
    {
        NoAA = 0,
        SSAA,
        MSAA,
        ProgressiveAA
    };

    enum class AAQuality : quint8
    {
        Normal = 2,
        High = 4,
        VeryHigh = 8
    };

    enum class HorizontalField : quint8
    {
        LeftWidth = 0,
        LeftRight,
        WidthRight
    };

    enum class VerticalField : quint8
    {
        TopHeight = 0,
        TopBottom,
        HeightBottom
    };

    enum class UnitType : quint8
    {
        Percent = 0,
        Pixels
    };

    enum class Background : quint8
    {
        Transparent = 0,
        Unspecified,
        Color,
        SkyBox
    };

    enum class BlendMode : quint8
    {
        SourceOver = 0,
        Screen,
        Multiply,
        Add,
        Subtract,
        Overlay,
        ColorBurn,
        ColorDodge
    };

    // First effect in a list of effects.
    QSSGRenderEffect *firstEffect;

    // If a layer has a valid texture path (one that resolves to either a
    // an on-disk image or a offscreen renderer), then it does not render its
    // own source path.  Instead, it renders the offscreen renderer.  Used in this manner,
    // offscreen renderer's also have the option (if they support it) to render directly to the
    // render target given a specific viewport (that is also scissored if necessary).
    QString texturePath;

    //SRenderPlugin *renderPlugin; // Overrides texture path if available.

    QSSGRenderLayer::AAMode antialiasingMode;
    QSSGRenderLayer::AAQuality antialiasingQuality;

    QSSGRenderLayer::Background background;
    QVector3D clearColor;

    BlendMode blendType;

    // TODO: pack
    HorizontalField horizontalFieldValues;
    float m_left;
    UnitType leftUnits;
    float m_width;
    UnitType widthUnits;
    float m_right;
    UnitType rightUnits;

    VerticalField verticalFieldValues;
    float m_top;
    UnitType topUnits;
    float m_height;
    UnitType heightUnits;
    float m_bottom;
    UnitType bottomUnits;

    // Ambient occlusion
    float aoStrength;
    float aoDistance;
    float aoSoftness;
    float aoBias;
    qint32 aoSamplerate;
    bool aoDither;

    // Direct occlusion
    float shadowStrength;
    float shadowDist;
    float shadowSoftness;
    float shadowBias;

    // IBL
    QSSGRenderImage *lightProbe;
    float probeBright;
    bool fastIbl;
    float probeHorizon;
    float probeFov;
    QSSGRenderImage *lightProbe2;
    float probe2Fade;
    float probe2Window;
    float probe2Pos;

    bool temporalAAEnabled;
    float temporalAAStrength;
    float ssaaMultiplier;

    QSSGRenderCamera *activeCamera;
    // It is the used camera for the scene.
    // If activeCamera is not GloballyActive,
    // the first GloballyActive one will be used for render.
    QSSGRenderCamera *renderedCamera;

    QSSGRenderLayer();

    void addEffect(QSSGRenderEffect &inEffect);

    QSSGRenderEffect *getLastEffect();

    QSSGRenderLayer::BlendMode getLayerBlend() { return blendType; }
};
QT_END_NAMESPACE

#endif

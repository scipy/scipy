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

#ifndef QSSG_RENDERER_IMPL_SHADERS_H
#define QSSG_RENDERER_IMPL_SHADERS_H

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

#include <QtQuick3DRender/private/qssgrendershaderprogram_p.h>
#include <QtQuick3DRender/private/qssgrenderprogrampipeline_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderableobjects_p.h>

QT_BEGIN_NAMESPACE

/**
 *	Cached tessellation property lookups this is on a per mesh base
 */
struct QSSGShaderTessellationProperties
{
    QSSGRenderCachedShaderProperty<float> edgeTessLevel; ///< tesselation value for the edges
    QSSGRenderCachedShaderProperty<float> insideTessLevel; ///< tesselation value for the inside
    QSSGRenderCachedShaderProperty<float> phongBlend; ///< blending between linear and phong component
    QSSGRenderCachedShaderProperty<QVector2D> distanceRange; ///< distance range for min and max tess level
    QSSGRenderCachedShaderProperty<float> disableCulling; ///< if set to 1.0 this disables
    /// backface culling optimization in
    /// the tess shader

    QSSGShaderTessellationProperties() = default;
    QSSGShaderTessellationProperties(const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : edgeTessLevel("tessLevelOuter", inShader)
        , insideTessLevel("tessLevelInner", inShader)
        , phongBlend("phongBlend", inShader)
        , distanceRange("distanceRange", inShader)
        , disableCulling("disableCulling", inShader)
    {
    }
};

/**
 *	The results of generating a shader.  Caches all possible variable names into
 *	typesafe objects.
 */
struct QSSGShaderGeneratorGeneratedShader
{
    QAtomicInt ref;
    quint32 layerSetIndex;
    QByteArray queryString;
    QSSGRef<QSSGRenderShaderProgram> shader;
    QSSGRenderCachedShaderProperty<QMatrix4x4> viewportMatrix;
    QSSGShaderTessellationProperties tessellation;

    QSSGShaderGeneratorGeneratedShader(const QByteArray &inQueryString, const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : layerSetIndex(std::numeric_limits<quint32>::max())
        , queryString(inQueryString)
        , shader(inShader)
        , viewportMatrix("viewportMatrix", inShader)
        , tessellation(inShader)
    {
    }

    static quint32 getLayerIndex(const QSSGShaderGeneratorGeneratedShader &inShader)
    {
        return inShader.layerSetIndex;
    }
    static void setLayerIndex(QSSGShaderGeneratorGeneratedShader &inShader, quint32 idx)
    {
        inShader.layerSetIndex = idx;
    }
};

struct QSSGDefaultMaterialRenderableDepthShader
{
    QAtomicInt ref;

    QSSGRef<QSSGRenderShaderProgram> shader;
    QSSGRenderCachedShaderProperty<QMatrix4x4> mvp;

    QSSGDefaultMaterialRenderableDepthShader(const QSSGRef<QSSGRenderShaderProgram> &inShader, QSSGRenderContext &inContext)
        : shader(inShader), mvp("modelViewProjection", inShader)
    {
        // TODO:
        Q_UNUSED(inContext)
    }
};

/**
 *	Cached texture property lookups, used one per texture so a shader generator for N
 *	textures will have an array of N of these lookup objects.
 */
struct QSSGShaderTextureProperties
{
    QSSGRenderCachedShaderProperty<QSSGRenderTexture2D *> sampler;
    QSSGRenderCachedShaderProperty<QVector3D> offsets;
    QSSGRenderCachedShaderProperty<QVector4D> rotations;
    QSSGRenderCachedShaderProperty<QVector2D> size;
    QSSGShaderTextureProperties(const QSSGRef<QSSGRenderShaderProgram> &inShader,
                                  const QByteArray &sampName,
                                  const QByteArray &offName,
                                  const QByteArray &rotName,
                                  const QByteArray &sizeName = QByteArray())
        : sampler(sampName, inShader),
          offsets(offName, inShader),
          rotations(rotName, inShader),
          size(sizeName, inShader)
    {
    }
    QSSGShaderTextureProperties() = default;
};

struct QSSGRenderableDepthPrepassShader
{
    QAtomicInt ref;
    QSSGRef<QSSGRenderShaderProgram> shader;
    QSSGRenderCachedShaderProperty<QMatrix4x4> mvp;
    QSSGRenderCachedShaderProperty<QMatrix4x4> globalTransform;
    QSSGRenderCachedShaderProperty<QMatrix4x4> projection;
    QSSGRenderCachedShaderProperty<QVector3D> cameraPosition;
    QSSGRenderCachedShaderProperty<float> displaceAmount;
    QSSGShaderTextureProperties displacementProps;
    QSSGRenderCachedShaderProperty<QVector2D> cameraProperties;
    QSSGRenderCachedShaderProperty<QVector3D> cameraDirection;
    // QSSGRenderCachedShaderProperty<QMatrix4x4> shadowMv[6];

    // Cache the tessellation property name lookups
    QSSGShaderTessellationProperties tessellation;

    QSSGRenderableDepthPrepassShader(const QSSGRef<QSSGRenderShaderProgram> &inShader, const QSSGRef<QSSGRenderContext> &inContext)
        : shader(inShader)
        , mvp("modelViewProjection", inShader)
        , globalTransform("modelMatrix", inShader)
        , projection("projection", inShader)
        , cameraPosition("cameraPosition", inShader)
        , displaceAmount("displaceAmount", inShader)
        , displacementProps(inShader, "displacementSampler", "displacementMap_offset", "displacementMap_rot")
        , cameraProperties("cameraProperties", inShader)
        , cameraDirection("cameraDirection", inShader)
        , tessellation(inShader)
    {
        Q_UNUSED(inContext)
        /*
            shadowMv[0].m_Shader = &inShader;
            shadowMv[0].m_Constant = inShader.GetShaderConstant( "shadow_mv0" );
            shadowMv[1].m_Shader = &inShader;
            shadowMv[1].m_Constant = inShader.GetShaderConstant( "shadow_mv1" );
            shadowMv[2].m_Shader = &inShader;
            shadowMv[2].m_Constant = inShader.GetShaderConstant( "shadow_mv2" );
            shadowMv[3].m_Shader = &inShader;
            shadowMv[3].m_Constant = inShader.GetShaderConstant( "shadow_mv3" );
            shadowMv[4].m_Shader = &inShader;
            shadowMv[4].m_Constant = inShader.GetShaderConstant( "shadow_mv4" );
            shadowMv[5].m_Shader = &inShader;
            shadowMv[5].m_Constant = inShader.GetShaderConstant( "shadow_mv5" );
            */
    }

    ~QSSGRenderableDepthPrepassShader() {}
};

struct QSSGSkyBoxShader
{
    QAtomicInt ref;
    QSSGRef<QSSGRenderShaderProgram> shader;
    QSSGRenderCachedShaderProperty<QMatrix4x4> viewMatrix;
    QSSGRenderCachedShaderProperty<QMatrix4x4> inverseProjection;
    QSSGRenderCachedShaderProperty<QSSGRenderTexture2D *> skyboxTexture;

    QSSGSkyBoxShader(const QSSGRef<QSSGRenderShaderProgram> &inShader, const QSSGRef<QSSGRenderContext> &inContext)
        : shader(inShader)
        , viewMatrix("viewMatrix", inShader)
        , inverseProjection("inverseProjection", inShader)
        , skyboxTexture("skybox_image", inShader)
    {
        Q_UNUSED(inContext)
    }
    ~QSSGSkyBoxShader() = default;
};

struct QSSGDefaultAoPassShader
{
    QAtomicInt ref;
    QSSGRef<QSSGRenderShaderProgram> shader;
    QSSGRenderCachedShaderProperty<QMatrix4x4> viewMatrix;
    QSSGRenderCachedShaderProperty<QVector2D> cameraProperties;
    QSSGRenderCachedShaderProperty<QVector3D> cameraDirection;
    QSSGRenderCachedShaderProperty<QSSGRenderTexture2D *> depthTexture;
    QSSGRenderCachedShaderProperty<QSSGRenderTextureCube *> cubeTexture;
    QSSGRenderCachedShaderProperty<QVector2D> depthTextureSize;

    QSSGRenderCachedShaderBuffer<QSSGRenderShaderConstantBuffer> aoShadowParams;

    QSSGDefaultAoPassShader(const QSSGRef<QSSGRenderShaderProgram> &inShader, const QSSGRef<QSSGRenderContext> &inContext)
        : shader(inShader)
        , viewMatrix("viewMatrix", inShader)
        , cameraProperties("cameraProperties", inShader)
        , cameraDirection("cameraDirection", inShader)
        , depthTexture("depthTexture", inShader)
        , cubeTexture("depthCube", inShader)
        , depthTextureSize("depthTextureSize", inShader)
        , aoShadowParams("aoShadow", inShader)
    {
        Q_UNUSED(inContext)
    }
    ~QSSGDefaultAoPassShader() {}
};

struct QSSGLayerProgAABlendShader
{
    QAtomicInt ref;
    QSSGRef<QSSGRenderShaderProgram> shader;
    QSSGRenderCachedShaderProperty<QSSGRenderTexture2D *> accumSampler;
    QSSGRenderCachedShaderProperty<QSSGRenderTexture2D *> lastFrame;
    QSSGRenderCachedShaderProperty<QVector2D> blendFactors;
    QSSGLayerProgAABlendShader(const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : shader(inShader), accumSampler("accumulator", inShader), lastFrame("last_frame", inShader), blendFactors("blend_factors", inShader)
    {
    }
};

struct QSSGCompositShader
{
    QAtomicInt ref;
    QSSGRef<QSSGRenderShaderProgram> shader;
    QSSGRenderCachedShaderProperty<QSSGRenderTexture2D *> lastFrame;
    QSSGCompositShader(const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : shader(inShader), lastFrame("last_frame", inShader)
    {
    }
};

struct QSSGLayerLastFrameBlendShader
{
    QAtomicInt ref;
    QSSGRef<QSSGRenderShaderProgram> shader;
    QSSGRenderCachedShaderProperty<QSSGRenderTexture2D *> lastFrame;
    QSSGRenderCachedShaderProperty<float> blendFactor;
    QSSGLayerLastFrameBlendShader(const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : shader(inShader), lastFrame("last_frame", inShader), blendFactor("blend_factor", inShader)
    {
    }
};

struct QSSGFlippedQuadShader
{
    QAtomicInt ref;
    QSSGRef<QSSGRenderShaderProgram> shader;

    QSSGRenderCachedShaderProperty<QMatrix4x4> mvp;
    // Dimensions and offsetting of the image.
    QSSGRenderCachedShaderProperty<QVector2D> dimensions;
    // The fourth member of text color is the opacity
    QSSGRenderCachedShaderProperty<QSSGRenderTexture2D *> sampler;
    // Opacity to use for rendering
    QSSGRenderCachedShaderProperty<float> opacity;

    QSSGFlippedQuadShader(const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : shader(inShader)
        , mvp("modelViewProjection", inShader)
        , dimensions("layer_dimensions", inShader)
        , sampler("layer_image", inShader)
        , opacity("opacity", inShader)
    {
    }
    ~QSSGFlippedQuadShader() {}
};

struct QSSGShadowmapPreblurShader
{
    QAtomicInt ref;
    QSSGRef<QSSGRenderShaderProgram> shader;
    QSSGRenderCachedShaderProperty<QVector2D> cameraProperties;
    QSSGRenderCachedShaderProperty<QSSGRenderTextureCube *> depthCube;
    QSSGRenderCachedShaderProperty<QSSGRenderTexture2D *> depthMap;

    QSSGShadowmapPreblurShader(const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : shader(inShader), cameraProperties("cameraProperties", inShader), depthCube("depthCube", inShader), depthMap("depthSrc", inShader)
    {
    }
    ~QSSGShadowmapPreblurShader() {}
};

struct QSSGGGSGet
{
    quint32 operator()(const QSSGShaderGeneratorGeneratedShader &inShader) { return inShader.layerSetIndex; }
};
struct QSSGGGSSet
{
    void operator()(QSSGShaderGeneratorGeneratedShader &inShader, quint32 idx) { inShader.layerSetIndex = idx; }
};
QT_END_NAMESPACE
#endif

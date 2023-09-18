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

#ifndef QSSG_RENDER_DEFAULT_MATERIAL_H
#define QSSG_RENDER_DEFAULT_MATERIAL_H

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
#include <QtQuick3DRuntimeRender/private/qssgrendermaterialdirty_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderlightmaps_p.h>

#include <QtGui/QVector3D>

QT_BEGIN_NAMESPACE

struct QSSGRenderImage;

struct Q_QUICK3DRUNTIMERENDER_EXPORT QSSGRenderDefaultMaterial : QSSGRenderGraphObject
{
    enum class MaterialLighting : quint8
    {
        NoLighting = 0,
        FragmentLighting
    };
    enum class MaterialBlendMode : quint8
    {
        SourceOver = 0,
        Screen,
        Multiply,
        Overlay,
        ColorBurn,
        ColorDodge
    };
    enum class MaterialSpecularModel : quint8
    {
        Default = 0,
        KGGX,
        KWard
    };
    enum MaterialAlphaMode : quint8
    {
        Opaque = 0,
        Mask,
        Blend,
        Default
    };
    enum TextureChannelMapping : quint8
    {
        R = 0,
        G,
        B,
        A
    };

    // Materials are stored as a linked list on models.
    QSSGRenderGraphObject *nextSibling = nullptr;
    QSSGRenderModel *parent = nullptr;
    QSSGRenderImage *colorMap = nullptr;
    // material section
    QSSGRenderImage *iblProbe = nullptr;
    QSSGRenderImage *emissiveMap = nullptr;
    QSSGRenderImage *specularReflection = nullptr;
    QSSGRenderImage *specularMap = nullptr;
    QSSGRenderImage *roughnessMap = nullptr;
    QSSGRenderImage *opacityMap = nullptr;
    QSSGRenderImage *bumpMap = nullptr;
    QSSGRenderImage *normalMap = nullptr;
    QSSGRenderImage *displacementMap = nullptr;
    QSSGRenderImage *translucencyMap = nullptr;
    QSSGRenderImage *metalnessMap = nullptr;
    QSSGRenderImage *occlusionMap = nullptr;
    // lightmap section
    QSSGRenderLightmaps lightmaps;

    QVector3D specularTint{ 1.0f, 1.0f, 1.0f };
    float ior = 0.2f;
    QVector3D emissiveColor = { 1.0f, 1.0f, 1.0f };
    QVector4D color{ 1.0f, 1.0f, 1.0f, 1.0f }; // colors are 0-1 normalized
    float diffuseLightWrap = 0.0f; // 0 - 1
    float fresnelPower = 0.0f;
    float specularAmount = 0.0f; // 0-??, defaults to 0
    float specularRoughness = 50.0f; // 0-??, defaults to 50
    float metalnessAmount = 0.0f;
    float opacity = 1.0f; // 0-1
    float bumpAmount = 0.0f; // 0-??
    float displaceAmount = 0.0f; // 0-??
    float translucentFalloff = 0.0f; // 0 - ??
    float occlusionAmount = 1.0f; // 0 - 1
    float alphaCutoff = 0.5f; // 0 - 1

    QSSGMaterialDirty dirty;
    MaterialLighting lighting = MaterialLighting::FragmentLighting;
    QSSGRenderDefaultMaterial::MaterialBlendMode blendMode = QSSGRenderDefaultMaterial::MaterialBlendMode::SourceOver;
    QSSGRenderDefaultMaterial::MaterialSpecularModel specularModel = QSSGRenderDefaultMaterial::MaterialSpecularModel::Default;
    QSSGRenderDefaultMaterial::MaterialAlphaMode alphaMode = QSSGRenderDefaultMaterial::Default;
    QSSGCullFaceMode cullMode = QSSGCullFaceMode::Back;
    bool vertexColorsEnabled = false;
    TextureChannelMapping roughnessChannel = TextureChannelMapping::R;
    TextureChannelMapping opacityChannel = TextureChannelMapping::A;
    TextureChannelMapping translucencyChannel = TextureChannelMapping::A;
    TextureChannelMapping metalnessChannel = TextureChannelMapping::R;
    TextureChannelMapping occlusionChannel = TextureChannelMapping::R;

    QSSGRenderDefaultMaterial(Type type = Type::DefaultMaterial);

    bool isSpecularEnabled() const { return specularAmount > .01f; }
    bool isMetalnessEnabled() const { return metalnessAmount > 0.01f; }
    bool isFresnelEnabled() const { return fresnelPower > 0.0f; }
    bool isVertexColorsEnabled() const { return vertexColorsEnabled; }
    bool hasLighting() const { return lighting != MaterialLighting::NoLighting; }
};

QT_END_NAMESPACE

#endif

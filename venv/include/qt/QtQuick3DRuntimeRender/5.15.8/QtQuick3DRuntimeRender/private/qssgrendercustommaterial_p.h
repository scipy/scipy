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

#ifndef QSSG_RENDER_CUSTOM_MATERIAL_H
#define QSSG_RENDER_CUSTOM_MATERIAL_H

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

#include <QtQuick3DRuntimeRender/private/qssgrenderimage_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderlightmaps_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderdynamicobjectsystemcommands_p.h>

#include <QtCore/qurl.h>
#include <QtCore/qvector.h>

QT_BEGIN_NAMESPACE

struct Q_QUICK3DRUNTIMERENDER_EXPORT QSSGRenderCustomMaterial : public QSSGRenderGraphObject
{
    QSSGRenderCustomMaterial() : QSSGRenderGraphObject(Type::CustomMaterial) {}

    struct TextureProperty
    {
        QSSGRenderImage *texImage = nullptr;
        QByteArray name;
        QSSGRenderShaderDataType shaderDataType;
        QSSGRenderTextureMagnifyingOp magFilterType = QSSGRenderTextureMagnifyingOp::Linear;
        QSSGRenderTextureMinifyingOp minFilterType = QSSGRenderTextureMinifyingOp::Linear;
        QSSGRenderTextureCoordOp clampType = QSSGRenderTextureCoordOp::ClampToEdge;
    };

    QVector<TextureProperty> textureProperties;

    struct Property
    {
        Property() = default;
        Property(const QByteArray &name, const QVariant &value, QSSGRenderShaderDataType shaderDataType, int pid = -1)
            : name(name), value(value), shaderDataType(shaderDataType), pid(pid)
        { }
        QByteArray name;
        mutable QVariant value;
        QSSGRenderShaderDataType shaderDataType;
        int pid;
    };

    QVector<Property> properties;

    struct ShaderInfo
    {
        QByteArray version;
        QByteArray type; // I.e., GLSL
        QByteArray shaderPrefix;
    };

    ShaderInfo shaderInfo;
    QMap<QByteArray, QByteArray> shaders;

    QVector<dynamic::QSSGCommand *> commands;

    // IMPORTANT: These flags matches the key produced by a MDL export file
    enum class MaterialShaderKeyValues
    {
        diffuse = 1 << 0,
        specular = 1 << 1,
        cutout = 1 << 2,
        refraction = 1 << 3,
        transparent = 1 << 4,
        displace = 1 << 5,
        transmissive = 1 << 6,
    };
    Q_DECLARE_FLAGS(MaterialShaderKeyFlags, MaterialShaderKeyValues)

    using Flag = QSSGRenderNode::Flag;
    Q_DECLARE_FLAGS(Flags, Flag)

    const char *className = nullptr;

    // lightmap section
    QSSGRenderLightmaps m_lightmaps;
    // material section
    bool m_hasTransparency = false;
    bool m_hasRefraction = false;
    QSSGRenderImage *m_iblProbe = nullptr;
    QSSGRenderImage *m_emissiveMap = nullptr;
    QSSGRenderImage *m_displacementMap = nullptr;
    float m_displaceAmount = 0.0f; ///< depends on the object size
    QSSGCullFaceMode cullMode = QSSGCullFaceMode::Back;

    QSSGRenderGraphObject *m_nextSibling = nullptr;

    MaterialShaderKeyFlags m_shaderKeyValues; ///< input from MDL files

    Flags flags;
    bool m_alwaysDirty = false;

    bool isDielectric() const { return m_shaderKeyValues & MaterialShaderKeyValues::diffuse; }
    bool isSpecularEnabled() const { return m_shaderKeyValues & MaterialShaderKeyValues::specular; }
    bool isCutOutEnabled() const { return m_shaderKeyValues & MaterialShaderKeyValues::cutout; }
    bool isTransmissive() const { return m_shaderKeyValues & MaterialShaderKeyValues::transmissive; }
    bool hasLighting() const { return true; }

    // Dirty
    bool m_dirtyFlagWithInFrame;
    bool isDirty() const { return flags.testFlag(Flag::Dirty) || m_dirtyFlagWithInFrame || m_alwaysDirty; }
    void updateDirtyForFrame()
    {
        m_dirtyFlagWithInFrame = flags.testFlag(Flag::Dirty);
        flags.setFlag(Flag::Dirty, false);
    }
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QSSGRenderCustomMaterial::MaterialShaderKeyFlags)

QT_END_NAMESPACE

#endif

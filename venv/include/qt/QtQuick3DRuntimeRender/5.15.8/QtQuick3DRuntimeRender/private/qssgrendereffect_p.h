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

#ifndef QSSG_RENDER_EFFECT_H
#define QSSG_RENDER_EFFECT_H

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
#include <QtQuick3DRuntimeRender/private/qssgrendernode_p.h>

#include <QtQuick3DRuntimeRender/private/qssgrenderimage_p.h>

QT_BEGIN_NAMESPACE
struct QSSGRenderLayer;
struct QSSGEffectContext;
class QSSGEffectSystem;

namespace dynamic
{
struct QSSGCommand;
}

// Effects are post-render effect applied to the layer.  There can be more than one of
// them and they have completely variable properties.
struct Q_QUICK3DRUNTIMERENDER_EXPORT QSSGRenderEffect : public QSSGRenderGraphObject
{
    QSSGRenderEffect();
    ~QSSGRenderEffect();

    struct TextureProperty
    {
        QSSGRenderImage *texImage = nullptr;
        QByteArray name;
        QSSGRenderShaderDataType shaderDataType;
        QSSGRenderTextureMagnifyingOp magFilterType = QSSGRenderTextureMagnifyingOp::Linear;
        QSSGRenderTextureMinifyingOp minFilterType = QSSGRenderTextureMinifyingOp::Linear;
        QSSGRenderTextureCoordOp clampType = QSSGRenderTextureCoordOp::ClampToEdge;
        QSSGRenderTextureTypeValue usageType;
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

    QSSGRenderLayer *m_layer;
    QSSGRenderEffect *m_nextEffect;
    // Opaque pointer to context type implemented by the effect system.
    // May be null in which case the effect system will generate a new context
    // the first time it needs to render this effect.
    QSSGEffectContext *m_context = nullptr;

    void initialize();

    // If our active flag value changes, then we ask the effect manager
    // to reset our context.
    void setActive(bool inActive, QSSGEffectSystem &inSystem);

    void reset(QSSGEffectSystem &inSystem);

    using Flag = QSSGRenderNode::Flag;
    Q_DECLARE_FLAGS(Flags, Flag)

    QVector<dynamic::QSSGCommand *> commands;

    Flags flags;
    const char *className = nullptr;
    bool requiresDepthTexture = false;
    bool requiresCompilation = true;
    QSSGRenderTextureFormat::Format outputFormat = QSSGRenderTextureFormat::Unknown;
};

QT_END_NAMESPACE

#endif

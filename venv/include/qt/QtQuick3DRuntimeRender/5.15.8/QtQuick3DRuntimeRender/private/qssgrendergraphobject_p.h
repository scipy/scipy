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

#ifndef QSSG_RENDER_GRAPH_OBJECT_H
#define QSSG_RENDER_GRAPH_OBJECT_H

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

#include <QtQuick3DRuntimeRender/private/qtquick3druntimerenderglobal_p.h>
#include <QtCore/QString>

#define QSSG_DEBUG_ID 0

QT_BEGIN_NAMESPACE

// Types should be setup on construction.  Change the type
// at your own risk as the type is used for RTTI purposes.
struct Q_QUICK3DRUNTIMERENDER_EXPORT QSSGRenderGraphObject
{
    enum class Type : quint8
    {
        Unknown = 0,
        Presentation,
        Scene,
        Node,
        Layer,
        Light,
        Camera,
        Model,
        DefaultMaterial,
        PrincipledMaterial,
        Image,
        Effect,
        CustomMaterial,
        RenderPlugin,
        Lightmaps,
        Geometry,
        Item2D,
        LastKnownGraphObjectType,
    };

    QAtomicInt ref;
    // Id's help debugging the object and are optionally set
#if QSSG_DEBUG_ID
    QByteArray id;
#endif
    // Type is used for RTTI purposes down the road.
    Type type;

    QSSGRenderGraphObject(QSSGRenderGraphObject::Type inType) : type(inType) {}
    virtual ~QSSGRenderGraphObject();

    // If you change any detail of the scene graph, or even *breath* on a
    // scene graph object, you need to bump this binary version so at least
    // we know if we can load a file or not.
    static quint32 getSceneGraphBinaryVersion() { return 1; }

    inline bool isMaterialType() const Q_DECL_NOTHROW
    {
        return (type == Type::CustomMaterial || type == Type::DefaultMaterial || type == Type::PrincipledMaterial);
    }

    inline bool isLightmapType() const Q_DECL_NOTHROW
    {
        return (type == Type::Lightmaps || type == Type::DefaultMaterial || type == Type::PrincipledMaterial);
    }

    inline bool isNodeType() const Q_DECL_NOTHROW
    {
        return (type == Type::Node ||
                type == Type::Layer ||
                type == Type::Light ||
                type == Type::Camera ||
                type == Type::Model);
    }

    inline bool isRenderableType() const Q_DECL_NOTHROW
    {
        return (type == Type::Model ||
                type == Type::Item2D);
    }
};

QT_END_NAMESPACE

#endif

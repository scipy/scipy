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

#ifndef QSSG_RENDER_BACKEND_RENDER_STATE_OBJECTS_GL_H
#define QSSG_RENDER_BACKEND_RENDER_STATE_OBJECTS_GL_H

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

QT_BEGIN_NAMESPACE

///< this class handles the shader input variables
class QSSGRenderBackendDepthStencilStateGL
{
public:
    ///< constructor
    QSSGRenderBackendDepthStencilStateGL(bool enableDepth,
                                           bool depthMask,
                                           QSSGRenderBoolOp depthFunc,
                                           bool enableStencil,
                                           QSSGRenderStencilFunction &stencilFuncFront,
                                           QSSGRenderStencilFunction &stencilFuncBack,
                                           QSSGRenderStencilOperation &depthStencilOpFront,
                                           QSSGRenderStencilOperation &depthStencilOpBack)
        : m_depthEnable(enableDepth)
        , m_depthMask(depthMask)
        , m_depthFunc(depthFunc)
        , m_stencilEnable(enableStencil)
        , m_stencilFuncFront(stencilFuncFront)
        , m_stencilFuncBack(stencilFuncBack)
        , m_depthStencilOpFront(depthStencilOpFront)
        , m_depthStencilOpBack(depthStencilOpBack)
    {
    }

    ///< constructor
    QSSGRenderBackendDepthStencilStateGL()
        : m_depthEnable(true), m_depthMask(true), m_depthFunc(QSSGRenderBoolOp::LessThanOrEqual), m_stencilEnable(false)
    {
    }

    ///< destructor
    ~QSSGRenderBackendDepthStencilStateGL() {}

    ///< assignement
    QSSGRenderBackendDepthStencilStateGL &operator=(const QSSGRenderBackendDepthStencilStateGL &rhs)
    {
        // Check for self-assignment!
        if (this == &rhs)
            return *this;

        m_depthEnable = rhs.m_depthEnable;
        m_depthMask = rhs.m_depthMask;
        m_depthFunc = rhs.m_depthFunc;
        m_stencilEnable = rhs.m_stencilEnable;
        m_stencilFuncFront = rhs.m_stencilFuncFront;
        m_stencilFuncBack = rhs.m_stencilFuncBack;
        m_depthStencilOpFront = rhs.m_depthStencilOpFront;
        m_depthStencilOpBack = rhs.m_depthStencilOpBack;

        return *this;
    }

    bool operator==(const QSSGRenderBackendDepthStencilStateGL &other) const
    {
        return (m_depthEnable == other.m_depthEnable && m_depthMask == other.m_depthMask && m_depthFunc == other.m_depthFunc
                && m_stencilEnable == other.m_stencilEnable && m_stencilFuncFront == other.m_stencilFuncFront
                && m_stencilFuncBack == other.m_stencilFuncBack && m_depthStencilOpFront == other.m_depthStencilOpFront
                && m_depthStencilOpBack == other.m_depthStencilOpBack);
    }

    bool m_depthEnable; ///< depth test enabled
    bool m_depthMask; ///< enable / disable depth writes
    QSSGRenderBoolOp m_depthFunc; ///< depth comparison func
    bool m_stencilEnable; ///< enable disable stencil test
    QSSGRenderStencilFunction m_stencilFuncFront; ///< stencil setup for front faces
    QSSGRenderStencilFunction m_stencilFuncBack; ///< stencil setup for back faces
    QSSGRenderStencilOperation m_depthStencilOpFront; ///< depth stencil operation for front faces
    QSSGRenderStencilOperation m_depthStencilOpBack; ///< depth stencil operation for back faces
};

class QSSGRenderBackendMiscStateGL
{
public:
    ///< constructor
    QSSGRenderBackendMiscStateGL() : m_patchVertexCount(1) {}

    quint32 m_patchVertexCount; ///< vertex count for a single patch primitive
};

class QSSGRenderBackendRasterizerStateGL
{
public:
    ///< constructor
    QSSGRenderBackendRasterizerStateGL(float depthBias, float depthScale)
        : m_depthBias(depthBias), m_depthScale(depthScale)
    {
    }
    ///< constructor
    QSSGRenderBackendRasterizerStateGL() : m_depthBias(0.0), m_depthScale(0.0) {}

    QSSGRenderBackendRasterizerStateGL &operator=(const QSSGRenderBackendRasterizerStateGL &rhs)
    {
        // Check for self-assignment!
        if (this == &rhs)
            return *this;

        m_depthBias = rhs.m_depthBias;
        m_depthScale = rhs.m_depthScale;

        return *this;
    }

    bool operator==(const QSSGRenderBackendRasterizerStateGL &other) const
    {
        // TODO: Added fuzzy compare to hide warning, but we should make sure if we actuall want this behavior
        // and disable the warning instead.
        return (qFuzzyCompare(m_depthBias, other.m_depthBias) && qFuzzyCompare(m_depthScale, other.m_depthScale));
    }

    float m_depthBias; ///< depth bias
    float m_depthScale; ///< mulitply constant
};

QT_END_NAMESPACE

#endif

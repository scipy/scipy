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

#ifndef QSSG_RENDER_DEPTH_STENCIL_STATE_H
#define QSSG_RENDER_DEPTH_STENCIL_STATE_H

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
#include <QtQuick3DRender/private/qssgrenderbackend_p.h>

QT_BEGIN_NAMESPACE

class QSSGRenderContext;

// currently this handles only stencil state
class Q_QUICK3DRENDER_EXPORT QSSGRenderDepthStencilState
{
public:
    QAtomicInt ref;

private:
    QSSGRef<QSSGRenderBackend> m_backend; ///< pointer to backend
    QSSGRenderBackend::QSSGRenderBackendDepthStencilStateObject m_handle; ///< opaque backend handle

public:
    /**
     * @brief constructor
     *
     * @param[in] context				Pointer to context
     * @param[in] fnd					Pointer to foundation
     * @param[in] enableDepth			enable depth test
     * @param[in] depthMask				enable depth writes
     * @param[in] depthFunc				depth compare function
     * @param[in] enableStencil			enable stencil test
     * @param[in] stencilFuncFront		stencil setup front faces
     * @param[in] stencilFuncBack		stencil setup back faces
     * @param[in] depthStencilOpFront	depth/stencil operations front faces
     * @param[in] depthStencilOpBack	depth/stencil operations back faces
     *
     * @return No return.
     */
    QSSGRenderDepthStencilState(const QSSGRef<QSSGRenderContext> &context,
                                  bool enableDepth,
                                  bool depthMask,
                                  QSSGRenderBoolOp depthFunc,
                                  bool enableStencil,
                                  QSSGRenderStencilFunction &stencilFuncFront,
                                  QSSGRenderStencilFunction &stencilFuncBack,
                                  QSSGRenderStencilOperation &depthStencilOpFront,
                                  QSSGRenderStencilOperation &depthStencilOpBack);

    ~QSSGRenderDepthStencilState();

    ///< various get functions
    const QSSGRenderStencilFunction stencilFunction(QSSGCullFaceMode face) const
    {
        return (face == QSSGCullFaceMode::Back) ? m_stencilFuncBack : m_stencilFuncFront;
    }
    const QSSGRenderStencilOperation stencilOperation(QSSGCullFaceMode face) const
    {
        return (face == QSSGCullFaceMode::Back) ? m_depthStencilOpBack : m_depthStencilOpFront;
    }
    QSSGRenderBoolOp depthFunction() const { return m_depthFunc; }
    bool depthEnabled() const { return m_depthEnabled; }
    bool stencilEnabled() const { return m_stencilEnabled; }
    bool depthMask() const { return m_depthMask; }

    /**
     * @brief get the backend object handle
     *
     * @return the backend object handle.
     */
    QSSGRenderBackend::QSSGRenderBackendDepthStencilStateObject handle()
    {
        return m_handle;
    }

private:
    bool m_depthEnabled; ///< depth test enabled
    bool m_depthMask; ///< depth writes enabled
    bool m_stencilEnabled; ///< stencil test enabled
    QSSGRenderBoolOp m_depthFunc; ///< depth comparison func
    QSSGRenderStencilFunction m_stencilFuncFront; ///< stencil setup front faces
    QSSGRenderStencilFunction m_stencilFuncBack; ///< stencil setup back faces
    QSSGRenderStencilOperation m_depthStencilOpFront; ///< depth stencil operation front faces
    QSSGRenderStencilOperation m_depthStencilOpBack; ///< depth stencil operation back faces
};

QT_END_NAMESPACE

#endif

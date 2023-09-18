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

#ifndef QSSG_RENDER_INPUT_ASSEMBLER_H
#define QSSG_RENDER_INPUT_ASSEMBLER_H

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

#include <QtQuick3DRender/private/qtquick3drenderglobal_p.h>
#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DRender/private/qssgrenderbackend_p.h>

QT_BEGIN_NAMESPACE

// forward declarations
class QSSGRenderContext;
class QSSGRenderBackend;
class QSSGRenderAttribLayout;

///< this class handles the vertex attribute layout setup
class Q_QUICK3DRENDER_EXPORT QSSGRenderInputAssembler
{
    Q_DISABLE_COPY(QSSGRenderInputAssembler)
public:
    QAtomicInt ref;

    /**
     * @brief constructor
     *
     *	NOTE: The limit for buffers count is currently 16
     *
     * @param[in] context			Pointer to context
     * @param[in] attribLayout		Pointer to QSSGRenderAttribLayout object
     * @param[in] buffers			list of vertex buffers
     * @param[in] indexBuffer		pointer to index buffer. Can be nullptr
     * @param[in] strides			list of strides of the buffer
     * @param[in] offsets			list of offsets into the buffer
     * @param[in] primType			primitive type used for drawing
     * @param[in] patchVertexCount	if primitive is "Patch" this is the vertex count for a
     *single patch
     *
     * @return No return.
     */
    QSSGRenderInputAssembler(const QSSGRef<QSSGRenderContext> &context,
                               const QSSGRef<QSSGRenderAttribLayout> &attribLayout,
                               QSSGDataView<QSSGRef<QSSGRenderVertexBuffer>> buffers,
                               const QSSGRef<QSSGRenderIndexBuffer> &indexBuffer,
                               QSSGDataView<quint32> strides,
                               QSSGDataView<quint32> offsets,
                               QSSGRenderDrawMode primType = QSSGRenderDrawMode::Triangles,
                               quint32 patchVertexCount = 1);
    ~QSSGRenderInputAssembler();

    /**
     * @brief get the backend object handle
     *
     * @return the backend object handle.
     */
    QSSGRenderBackend::QSSGRenderBackendInputAssemblerObject handle() const
    {
        return m_handle;
    }

    /**
     * @brief get the attached index buffer
     *
     * @return the index buffer
     */
    const QSSGRef<QSSGRenderIndexBuffer> &indexBuffer();

    /**
     * @brief get the index count of the attached index buffer (if any)
     *
     * @return the index buffer count
     */
    quint32 indexCount() const;

    /**
     * @brief get the vertex count of the buffer
     *		  Note this makes only sense if we have a single
     *		  interleaved buffer
     *
     * @return the vertex buffer count
     */
    quint32 vertexCount() const;

    /**
     * @brief get the primitive type used for drawing
     *
     * @return primitive type
     */
    QSSGRenderDrawMode drawMode() const { return m_drawMode; }

    /**
     * @brief set the per vertex patch count
     *
     * @return none
     */
    void setPatchVertexCount(quint32 count)
    {
        if (count != m_patchVertexCount) {
            // clamp to 1;
            m_patchVertexCount = (count == 0) ? 1 : count;
            ;
            m_backend->setPatchVertexCount(m_handle, m_patchVertexCount);
        }
    }

private:
    QSSGRef<QSSGRenderContext> m_context; ///< pointer to context
    QSSGRef<QSSGRenderBackend> m_backend; ///< pointer to backend

    QSSGRef<QSSGRenderAttribLayout> m_attribLayout; ///< pointer to attribute layout
    QVector<QSSGRef<QSSGRenderVertexBuffer>> m_vertexBuffers; ///< vertex buffers
    QSSGRef<QSSGRenderIndexBuffer> m_indexBuffer; ///< index buffer
    QSSGDataView<QSSGRenderBackend::QSSGRenderBackendBufferObject> m_vertexbufferHandles; ///< opaque vertex buffer backend handles

    QSSGRenderBackend::QSSGRenderBackendInputAssemblerObject m_handle; ///< opaque backend handle
    QSSGRenderDrawMode m_drawMode; ///< primitive type used for drawing
    quint32 m_patchVertexCount; ///< vertex count if primitive type is patch
};

QT_END_NAMESPACE

#endif

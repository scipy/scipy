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

#ifndef QSSG_RENDER_BACKEND_INPUT_ASSEMBLER_GL_H
#define QSSG_RENDER_BACKEND_INPUT_ASSEMBLER_GL_H

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

#include <QtQuick3DUtils/private/qssgoption_p.h>
#include <QtQuick3DUtils/private/qssgutils_p.h>

#include <QtQuick3DRender/private/qtquick3drenderglobal_p.h>
#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DRender/private/qssgrenderbackend_p.h>

#include <QtCore/QString>

QT_BEGIN_NAMESPACE

struct QSSGRenderBackendLayoutEntryGL
{
    QByteArray m_attribName; ///< must be the same name as used in the vertex shader
    quint8 m_normalize; ///< normalize parameter
    quint32 m_attribIndex; ///< attribute index
    quint32 m_type; ///< GL vertex format type @sa GL_FLOAT, GL_INT
    quint32 m_numComponents; ///< component count. max 4
    quint32 m_inputSlot; ///< Input slot where to fetch the data from
    quint32 m_offset; ///< offset in byte
};

///< this class handles the vertex attribute layout setup
class QSSGRenderBackendAttributeLayoutGL
{
public:
    ///< constructor
    QSSGRenderBackendAttributeLayoutGL(QSSGDataRef<QSSGRenderBackendLayoutEntryGL> entries, quint32 maxInputSlot)
        : m_layoutAttribEntries(entries), m_maxInputSlot(maxInputSlot)
    {
    }
    ///< destructor
    ~QSSGRenderBackendAttributeLayoutGL()
    {
        // We need to release the attrib name
        for (int idx = 0; idx != m_layoutAttribEntries.size(); ++idx)
            m_layoutAttribEntries[idx] = QSSGRenderBackendLayoutEntryGL();
    }

    QSSGRenderBackendLayoutEntryGL *getEntryByName(const QByteArray &entryName) const
    {
        for (int idx = 0; idx != m_layoutAttribEntries.size(); ++idx) {
            if (m_layoutAttribEntries[idx].m_attribName == entryName)
                return &m_layoutAttribEntries.mData[idx];
        }
        return nullptr;
    }

    QSSGOption<QSSGRenderBackendLayoutEntryGL> getEntryByAttribIndex(quint32 attribIndex) const
    {
        for (int idx = 0; idx != m_layoutAttribEntries.size(); ++idx)
        {
            if (m_layoutAttribEntries[idx].m_attribIndex == attribIndex)
                return m_layoutAttribEntries[idx];
        }
        return QSSGEmpty();
    }

    QSSGDataRef<QSSGRenderBackendLayoutEntryGL> m_layoutAttribEntries; ///< vertex attribute layout entries
    qint32 m_maxInputSlot; ///< max used input slot
};

///< this class handles the input assembler setup
class QSSGRenderBackendInputAssemblerGL
{
public:
    ///< constructor
    QSSGRenderBackendInputAssemblerGL(QSSGRenderBackendAttributeLayoutGL *attribLayout,
                                        QSSGDataView<QSSGRenderBackend::QSSGRenderBackendBufferObject> buffers,
                                        const QSSGRenderBackend::QSSGRenderBackendBufferObject indexBuffer,
                                        QSSGDataView<quint32> strides,
                                        QSSGDataView<quint32> offsets,
                                        quint32 patchVertexCount)
        : m_attribLayout(attribLayout)
        , m_vertexbufferHandles(buffers)
        , m_indexbufferHandle(indexBuffer)
        , m_vaoID(0)
        , m_cachedShaderHandle(0)
        , m_patchVertexCount(patchVertexCount)
        , m_strides(strides.size())
        , m_offsets(offsets.size())
    {
        memcpy(m_strides.data(), strides.begin(), strides.size()*sizeof(quint32));
        memcpy(m_offsets.data(), offsets.begin(), offsets.size()*sizeof(quint32));
    }
    ///< destructor
    ~QSSGRenderBackendInputAssemblerGL()
    {
    }

    QSSGRenderBackendAttributeLayoutGL *m_attribLayout; ///< pointer to attribute layout
    QSSGDataView<QSSGRenderBackend::QSSGRenderBackendBufferObject> m_vertexbufferHandles; ///< opaque vertex buffer backend handles
    QSSGRenderBackend::QSSGRenderBackendBufferObject m_indexbufferHandle; ///< opaque index buffer backend handles
    quint32 m_vaoID; ///< this is only used if GL version is greater or equal 3
    quint32 m_cachedShaderHandle; ///< this is the shader id which was last used with this object
    quint32 m_patchVertexCount; ///< vertex count for a single patch primitive
    QVector<quint32> m_strides; ///< buffer strides
    QVector<quint32> m_offsets; ///< buffer offsets
};

QT_END_NAMESPACE

#endif

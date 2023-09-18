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

#ifndef QSSG_RENDER__RENDER_RENDER_BUFFER_H
#define QSSG_RENDER__RENDER_RENDER_BUFFER_H

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

class Q_QUICK3DRENDER_EXPORT QSSGRenderRenderBuffer
{
    Q_DISABLE_COPY(QSSGRenderRenderBuffer)
public:
    QAtomicInt ref;

private:
    QSSGRef<QSSGRenderContext> m_context; ///< pointer to context
    QSSGRef<QSSGRenderBackend> m_backend; ///< pointer to backend
    qint32 m_width; ///< buffer width
    qint32 m_height; ///< buffer height
    QSSGRenderRenderBufferFormat m_storageFormat; ///< buffer storage format

    QSSGRenderBackend::QSSGRenderBackendRenderbufferObject m_handle; ///< opaque backend handle

public:
    /**
     * @brief constructor
     *
     * @param[in] context		Pointer to context
     * @param[in] fnd			Pointer to foundation
     * @param[in] format		Renderbuffer format
     * @param[in] width			Renderbuffer width
     * @param[in] height		Renderbuffer height
     *
     * @return No return.
     */
    QSSGRenderRenderBuffer(const QSSGRef<QSSGRenderContext> &context,
                             QSSGRenderRenderBufferFormat format,
                             quint32 width,
                             quint32 height);

    /// destructor
    ~QSSGRenderRenderBuffer();

    /**
     * @brief query buffer format
     *
     *
     * @return buffer format
     */
    QSSGRenderRenderBufferFormat storageFormat() const { return m_storageFormat; }

    /**
     * @brief query buffer dimension
     *
     *
     * @return QSSGRenderRenderBufferDimensions object
     */
    QSize size() const
    {
        return QSize(m_width, m_height);
    }

    /**
     * @brief constructor
     *
     * @param[in] inDimensions		A dimension object
     *
     * @return buffer format
     */
    void setSize(const QSize &inDimensions);

    /**
     * @brief get the backend object handle
     *
     * @return the backend object handle.
     */
    QSSGRenderBackend::QSSGRenderBackendRenderbufferObject handle()
    {
        return m_handle;
    }
};

QT_END_NAMESPACE

#endif

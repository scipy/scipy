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

#ifndef QSSG_RENDER_STORAGE_BUFFER_H
#define QSSG_RENDER_STORAGE_BUFFER_H

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

#include <QtQuick3DRender/private/qssgrenderdatabuffer_p.h>

#include <QtCore/QString>

QT_BEGIN_NAMESPACE

// forward declaration
class QSSGRenderContext;
class QSSGRenderVertexBuffer;

///< Constant (uniform) buffer representation
class Q_QUICK3DRENDER_EXPORT QSSGRenderStorageBuffer : public QSSGRenderDataBuffer
{
public:
    /**
     * @brief constructor
     *
     * @param[in] context		Pointer to context
     * @param[in] bufferName	Name of the buffer. Must match the name used in programs
     * @param[in] size			Size of the buffer
     * @param[in] usage			Usage of the buffer (e.g. static, dynamic...)
     * @param[in] data			A pointer to the buffer data that is allocated by the
     * application.
     * @param[in] pBuffer		Pointer to the buffer
     *
     * @return No return.
     */
    QSSGRenderStorageBuffer(const QSSGRef<QSSGRenderContext> &context,
                              const QByteArray &bufferName,
                              QSSGRenderBufferUsageType usageType,
                              QSSGByteView data,
                              QSSGRenderDataBuffer *pBuffer = nullptr);

    ///< destructor
    virtual ~QSSGRenderStorageBuffer();

    /**
     * @brief bind the buffer bypasses the context state
     *
     * @return no return.
     */
    void bind() override;

    /**
     * @brief bind the buffer to a shader program
     *
     * @param[in] index			Index of the constant buffer within the program
     *
     * @return no return.
     */
    void bindToShaderProgram(quint32 index);

    /**
     * @brief update the buffer to hardware
     *
     * @return no return.
     */
    void update();

    /**
     * @brief update a piece of memory directly within the storage buffer
     *
     * Note: When you use this function you should know what you are doing.
     *		 The memory layout within C++ must exactly match the memory layout in the
     *shader.
     *		 We use std140 (430) layout which guarantees a specific layout behavior across
     *all HW vendors.
     *		 How the memory layout is computed can be found in the GL spec.
     *
     * @param[in] offset	offset into storage buffer
     * @param[in] data		pointer to data
     *
     * @return no return
     */
    void updateData(qint32 offset, QSSGByteView data);

    /**
     * @brief get the buffer name
     *
     * @return the buffer name
     */
    QByteArray name() const { return m_name; }

private:
    QByteArray m_name; ///< buffer name
    QSSGRenderDataBuffer *m_wrappedBuffer; ///< pointer to wrapped buffer
    bool m_dirty; ///< true if buffer is dirty
};

QT_END_NAMESPACE

#endif

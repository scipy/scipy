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

#ifndef QSSG_RENDER_ATTRIB_LAYOUT_H
#define QSSG_RENDER_ATTRIB_LAYOUT_H

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

// forward declarations
class QSSGRenderContext;
class QSSGRenderBackend;

///< this class handles the vertex attribute layout setup
class Q_QUICK3DRENDER_EXPORT QSSGRenderAttribLayout
{
    Q_DISABLE_COPY(QSSGRenderAttribLayout)
public:
    QAtomicInt ref;
    /**
     * @brief constructor
     *
     * @param[in] context		Pointer to context
     * @param[in] attribs		Pointer to attribute list
     *
     * @return No return.
     */
    QSSGRenderAttribLayout(const QSSGRef<QSSGRenderContext> &context,
                             QSSGDataView<QSSGRenderVertexBufferEntry> attribs);
    ///< destructor
    ~QSSGRenderAttribLayout();

    /**
     * @brief get the backend object handle
     *
     * @return the backend object handle.
     */
    QSSGRenderBackend::QSSGRenderBackendAttribLayoutObject handle() const
    {
        return m_attribLayoutHandle;
    }

private:
    QSSGRef<QSSGRenderContext> m_context; ///< pointer to context
    QSSGRef<QSSGRenderBackend> m_backend; ///< pointer to backend

    QSSGRenderBackend::QSSGRenderBackendAttribLayoutObject m_attribLayoutHandle; ///< opaque backend handle
};

QT_END_NAMESPACE

#endif

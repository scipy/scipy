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

#ifndef QSSG_RENDER_RESOURCE_BUFFER_OBJECTS_H
#define QSSG_RENDER_RESOURCE_BUFFER_OBJECTS_H

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

#include <QtQuick3DRender/private/qssgrendercontext_p.h>
#include <QtQuick3DRender/private/qssgrenderframebuffer_p.h>
#include <QtQuick3DRender/private/qssgrenderrenderbuffer_p.h>

#include <QtQuick3DRuntimeRender/private/qssgrenderresourcemanager_p.h>

QT_BEGIN_NAMESPACE
class QSSGResourceFrameBuffer
{
protected:
    QSSGRef<QSSGResourceManager> m_resourceManager;
    QSSGRef<QSSGRenderFrameBuffer> m_frameBuffer;

public:
    QSSGResourceFrameBuffer(const QSSGRef<QSSGResourceManager> &mgr);
    ~QSSGResourceFrameBuffer();
    bool ensureFrameBuffer();
    void releaseFrameBuffer();

    const QSSGRef<QSSGResourceManager> &getResourceManager() { return m_resourceManager; }
    const QSSGRef<QSSGRenderFrameBuffer> &getFrameBuffer() { return m_frameBuffer; }
    operator const QSSGRef<QSSGRenderFrameBuffer> &() { return m_frameBuffer; }
    const QSSGRef<QSSGRenderFrameBuffer> &operator->()
    {
        Q_ASSERT(m_frameBuffer);
        return m_frameBuffer;
    }
    QSSGRenderFrameBuffer &operator*()
    {
        Q_ASSERT(m_frameBuffer);
        return *m_frameBuffer;
    }
};

class QSSGResourceRenderBuffer
{
protected:
    QSSGRef<QSSGResourceManager> m_resourceManager;
    QSSGRef<QSSGRenderRenderBuffer> m_renderBuffer;
    QSSGRenderRenderBufferFormat m_storageFormat;
    QSize m_dimensions;

public:
    QSSGResourceRenderBuffer(const QSSGRef<QSSGResourceManager> &mgr);
    ~QSSGResourceRenderBuffer();
    bool ensureRenderBuffer(qint32 width, qint32 height, QSSGRenderRenderBufferFormat storageFormat);
    void releaseRenderBuffer();

    operator const QSSGRef<QSSGRenderRenderBuffer> &() { return m_renderBuffer; }
    const QSSGRef<QSSGRenderRenderBuffer> &operator->()
    {
        Q_ASSERT(m_renderBuffer);
        return m_renderBuffer;
    }
    QSSGRenderRenderBuffer &operator*()
    {
        Q_ASSERT(m_renderBuffer);
        return *m_renderBuffer;
    }
};
QT_END_NAMESPACE
#endif

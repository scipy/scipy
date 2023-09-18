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

#ifndef QSSG_RENDER_RESOURCE_TEXTURE_2D_H
#define QSSG_RENDER_RESOURCE_TEXTURE_2D_H

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
#include <QtQuick3DRender/private/qssgrendertexture2d_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderresourcemanager_p.h>

QT_BEGIN_NAMESPACE
class QSSGResourceTexture2D
{
protected:
    QSSGRef<QSSGResourceManager> m_resourceManager;
    QSSGRef<QSSGRenderTexture2D> m_texture;
    QSSGTextureDetails m_textureDetails;

public:
    QSSGResourceTexture2D(const QSSGRef<QSSGResourceManager> &mgr,
                            const QSSGRef<QSSGRenderTexture2D> &inTexture = nullptr);
    // create and allocate the texture right away.
    QSSGResourceTexture2D(const QSSGRef<QSSGResourceManager> &mgr,
                            quint32 width,
                            quint32 height,
                            QSSGRenderTextureFormat inFormat,
                            quint32 inSamples = 1);
    ~QSSGResourceTexture2D();
    // Returns true if the texture matches the specs, false if the texture needs to be
    // reallocated
    bool textureMatches(qint32 width, qint32 height, QSSGRenderTextureFormat inFormat, qint32 inSamples = 1);

    // Returns true if the texture was allocated, false if nothing changed (no allocation).
    // Note this is the exact opposite of TextureMatches.
    bool ensureTexture(qint32 width, qint32 height, QSSGRenderTextureFormat inFormat, qint32 inSamples = 1);

    // Force release the texture.
    void releaseTexture();
    QSSGRenderTexture2D &operator*()
    {
        Q_ASSERT(m_texture);
        return *m_texture;
    }
    const QSSGRef<QSSGRenderTexture2D> &operator->() const
    {
        Q_ASSERT(m_texture);
        return m_texture;
    }
    operator const QSSGRef<QSSGRenderTexture2D>() & { return m_texture; }
    QSSGRef<QSSGRenderTexture2D> getTexture() const { return m_texture; }
    void forgetTexture();
    // Enforces single ownership rules.
    void stealTexture(QSSGResourceTexture2D &inOther);
    void swapTexture(QSSGResourceTexture2D &inOther);
    bool isNull() const { return m_texture.isNull(); }
};

QT_END_NAMESPACE

#endif

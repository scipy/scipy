/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DRENDER_RENDER_ATTACHMENTPACK_P_H
#define QT3DRENDER_RENDER_ATTACHMENTPACK_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <Qt3DRender/qrendertargetoutput.h>
#include <Qt3DRender/private/qt3drender_global_p.h>
#include <QVector>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
namespace Render {

struct Q_3DRENDERSHARED_PRIVATE_EXPORT Attachment
{
    Attachment()
        : m_mipLevel(0)
        , m_layer(0)
        , m_point(QRenderTargetOutput::Color0)
        , m_face(QAbstractTexture::CubeMapNegativeX)
    {}

    QString m_name;
    int m_mipLevel;
    int m_layer;
    Qt3DCore::QNodeId m_textureUuid;
    QRenderTargetOutput::AttachmentPoint m_point;
    QAbstractTexture::CubeMapFace m_face;
};

class RenderTarget;
class RenderTargetSelector;
class AttachmentManager;

class Q_3DRENDERSHARED_PRIVATE_EXPORT AttachmentPack
{
public:
    AttachmentPack();
    AttachmentPack(const RenderTarget *target,
                   AttachmentManager *attachmentManager,
                   const QVector<QRenderTargetOutput::AttachmentPoint> &drawBuffers = QVector<QRenderTargetOutput::AttachmentPoint>());

    QVector<Attachment> attachments() const { return m_attachments; }
    QVector<int> getGlDrawBuffers() const { return m_drawBuffers; }

    // return index of given attachment within actual draw buffers list
    int getDrawBufferIndex(QRenderTargetOutput::AttachmentPoint attachmentPoint) const;

private:
    QVector<Attachment> m_attachments;
    QVector<int> m_drawBuffers;
};

Q_3DRENDERSHARED_PRIVATE_EXPORT bool operator ==(const Attachment &a, const Attachment &b);
Q_3DRENDERSHARED_PRIVATE_EXPORT bool operator !=(const Attachment &a, const Attachment &b);

Q_3DRENDERSHARED_PRIVATE_EXPORT bool operator ==(const AttachmentPack &packA, const AttachmentPack &packB);
Q_3DRENDERSHARED_PRIVATE_EXPORT bool operator !=(const AttachmentPack &packA, const AttachmentPack &packB);

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_ATTACHMENTPACK_P_H

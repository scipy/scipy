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

#ifndef QT3DRENDER_QRENDERTARGETOUTPUT_H
#define QT3DRENDER_QRENDERTARGETOUTPUT_H

#include <Qt3DCore/qnode.h>
#include <Qt3DRender/qt3drender_global.h>
#include <Qt3DRender/QAbstractTexture>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAbstractTexture;
class QRenderTargetOutputPrivate;

class Q_3DRENDERSHARED_EXPORT QRenderTargetOutput : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(AttachmentPoint attachmentPoint READ attachmentPoint WRITE setAttachmentPoint NOTIFY attachmentPointChanged)
    Q_PROPERTY(QAbstractTexture *texture READ texture WRITE setTexture NOTIFY textureChanged)
    Q_PROPERTY(int mipLevel READ mipLevel WRITE setMipLevel NOTIFY mipLevelChanged)
    Q_PROPERTY(int layer READ layer WRITE setLayer NOTIFY layerChanged)
    Q_PROPERTY(Qt3DRender::QAbstractTexture::CubeMapFace face READ face WRITE setFace NOTIFY faceChanged)

public:
    enum AttachmentPoint {
        Color0 = 0,
        Color1,
        Color2,
        Color3,
        Color4,
        Color5,
        Color6,
        Color7,
        Color8,
        Color9,
        Color10,
        Color11,
        Color12,
        Color13,
        Color14,
        Color15,
        Depth,
        Stencil,
        DepthStencil
    };
    Q_ENUM(AttachmentPoint) // LCOV_EXCL_LINE

    explicit QRenderTargetOutput(Qt3DCore::QNode *parent = nullptr);
    ~QRenderTargetOutput();

    AttachmentPoint attachmentPoint() const;
    QAbstractTexture *texture() const;
    int mipLevel() const;
    int layer() const;
    QAbstractTexture::CubeMapFace face() const;

public Q_SLOTS:
    void setAttachmentPoint(AttachmentPoint attachmentPoint);
    void setTexture(QAbstractTexture *texture);
    void setMipLevel(int level);
    void setLayer(int layer);
    void setFace(QAbstractTexture::CubeMapFace face);

Q_SIGNALS:
    void attachmentPointChanged(AttachmentPoint attachmentPoint);
    void textureChanged(QAbstractTexture *texture);
    void mipLevelChanged(int mipLevel);
    void layerChanged(int layer);
    void faceChanged(QAbstractTexture::CubeMapFace face);

protected:
    explicit QRenderTargetOutput(QRenderTargetOutputPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
    Q_DECLARE_PRIVATE(QRenderTargetOutput)
};

} // namespace Qt3DRender

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DRender::QRenderTargetOutput::AttachmentPoint) // LCOV_EXCL_LINE

#endif // QT3DRENDER_QRENDERTARGETOUTPUT_H

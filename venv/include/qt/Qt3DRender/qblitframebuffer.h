/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QBLITFRAMEBUFFER_H
#define QT3DRENDER_QBLITFRAMEBUFFER_H

#include <Qt3DRender/qframegraphnode.h>
#include <Qt3DRender/qrendertargetoutput.h>
#include <QtCore/QRect>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QBlitFramebufferPrivate;
class QRenderTarget;

class Q_3DRENDERSHARED_EXPORT QBlitFramebuffer : public QFrameGraphNode
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QRenderTarget *source READ source WRITE setSource NOTIFY sourceChanged)
    Q_PROPERTY(Qt3DRender::QRenderTarget *destination READ destination WRITE setDestination NOTIFY destinationChanged)
    Q_PROPERTY(QRectF sourceRect READ sourceRect WRITE setSourceRect NOTIFY sourceRectChanged)
    Q_PROPERTY(QRectF destinationRect READ destinationRect WRITE setDestinationRect NOTIFY destinationRectChanged)
    Q_PROPERTY(Qt3DRender::QRenderTargetOutput::AttachmentPoint sourceAttachmentPoint READ sourceAttachmentPoint WRITE setSourceAttachmentPoint NOTIFY sourceAttachmentPointChanged)
    Q_PROPERTY(Qt3DRender::QRenderTargetOutput::AttachmentPoint destinationAttachmentPoint READ destinationAttachmentPoint WRITE setDestinationAttachmentPoint NOTIFY destinationAttachmentPointChanged)
    Q_PROPERTY(InterpolationMethod interpolationMethod READ interpolationMethod WRITE setInterpolationMethod NOTIFY interpolationMethodChanged)
public:
    enum InterpolationMethod {
        Nearest = 0,
        Linear,
    };
    Q_ENUM(InterpolationMethod) // LCOV_EXCL_LINE

    explicit QBlitFramebuffer(Qt3DCore::QNode *parent = nullptr);
    ~QBlitFramebuffer();

    QRenderTarget *source() const;
    QRenderTarget *destination() const;
    QRectF sourceRect() const;
    QRectF destinationRect() const;
    Qt3DRender::QRenderTargetOutput::AttachmentPoint sourceAttachmentPoint() const;
    Qt3DRender::QRenderTargetOutput::AttachmentPoint destinationAttachmentPoint() const;
    InterpolationMethod interpolationMethod() const;

    void setSource(QRenderTarget *source);
    void setDestination(QRenderTarget *destination);
    void setSourceRect(const QRectF &sourceRect);
    void setDestinationRect(const QRectF &destinationRect);
    void setSourceAttachmentPoint(Qt3DRender::QRenderTargetOutput::AttachmentPoint sourceAttachmentPoint);
    void setDestinationAttachmentPoint(Qt3DRender::QRenderTargetOutput::AttachmentPoint destinationAttachmentPoint);
    void setInterpolationMethod(InterpolationMethod interpolationMethod);

Q_SIGNALS:
    void sourceChanged();
    void destinationChanged();
    void sourceRectChanged();
    void destinationRectChanged();
    void sourceAttachmentPointChanged();
    void destinationAttachmentPointChanged();
    void interpolationMethodChanged();

protected:
    explicit QBlitFramebuffer(QBlitFramebufferPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QBlitFramebuffer)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QBLITFRAMEBUFFER_H

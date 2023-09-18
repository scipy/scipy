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

#ifndef QT3DEXTRAS_QFORWARDRENDERER_H
#define QT3DEXTRAS_QFORWARDRENDERER_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qtechniquefilter.h>
#include <QtCore/QRectF>
#include <QtGui/QColor>
#include <Qt3DRender/QClearBuffers>

QT_BEGIN_NAMESPACE

class QSurface;

namespace Qt3DExtras {

class QForwardRendererPrivate;

class Q_3DEXTRASSHARED_EXPORT QForwardRenderer : public Qt3DRender::QTechniqueFilter
{
    Q_OBJECT
    Q_PROPERTY(QObject *surface READ surface WRITE setSurface NOTIFY surfaceChanged)
    Q_PROPERTY(QObject *window READ surface WRITE setSurface NOTIFY surfaceChanged)
    Q_PROPERTY(QRectF viewportRect READ viewportRect WRITE setViewportRect NOTIFY viewportRectChanged)
    Q_PROPERTY(QColor clearColor READ clearColor WRITE setClearColor NOTIFY clearColorChanged)
    Q_PROPERTY(Qt3DRender::QClearBuffers::BufferType buffersToClear READ buffersToClear WRITE setBuffersToClear NOTIFY buffersToClearChanged REVISION 14)
    Q_PROPERTY(Qt3DCore::QEntity *camera READ camera WRITE setCamera NOTIFY cameraChanged)
    Q_PROPERTY(QSize externalRenderTargetSize READ externalRenderTargetSize WRITE setExternalRenderTargetSize NOTIFY externalRenderTargetSizeChanged)
    Q_PROPERTY(bool frustumCulling READ isFrustumCullingEnabled WRITE setFrustumCullingEnabled NOTIFY frustumCullingEnabledChanged)
    Q_PROPERTY(float gamma READ gamma WRITE setGamma NOTIFY gammaChanged REVISION 9)
    Q_PROPERTY(bool showDebugOverlay READ showDebugOverlay WRITE setShowDebugOverlay NOTIFY showDebugOverlayChanged REVISION 15)
public:
    explicit QForwardRenderer(Qt3DCore::QNode *parent = nullptr);
    ~QForwardRenderer();

    QRectF viewportRect() const;
    QColor clearColor() const;
    Qt3DRender::QClearBuffers::BufferType buffersToClear() const;
    Qt3DCore::QEntity *camera() const;
    QObject *surface() const;
    QSize externalRenderTargetSize() const;
    bool isFrustumCullingEnabled() const;
    float gamma() const;
    bool showDebugOverlay() const;

public Q_SLOTS:
    void setViewportRect(const QRectF &viewportRect);
    void setClearColor(const QColor &clearColor);
    void setBuffersToClear(Qt3DRender::QClearBuffers::BufferType);
    void setCamera(Qt3DCore::QEntity *camera);
    void setSurface(QObject * surface);
    void setExternalRenderTargetSize(const QSize &size);
    void setFrustumCullingEnabled(bool enabled);
    void setGamma(float gamma);
    void setShowDebugOverlay(bool showDebugOverlay);

Q_SIGNALS:
    void viewportRectChanged(const QRectF &viewportRect);
    void clearColorChanged(const QColor &clearColor);
    void buffersToClearChanged(Qt3DRender::QClearBuffers::BufferType);
    void cameraChanged(Qt3DCore::QEntity *camera);
    void surfaceChanged(QObject *surface);
    void externalRenderTargetSizeChanged(const QSize &size);
    void frustumCullingEnabledChanged(bool enabled);
    void gammaChanged(float gamma);
    void showDebugOverlayChanged(bool showDebugOverlay);

private:
    Q_DECLARE_PRIVATE(QForwardRenderer)
};

} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QFORWARDRENDERER_H

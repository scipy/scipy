/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QRENDERSURFACESELECTOR_H
#define QT3DRENDER_QRENDERSURFACESELECTOR_H

#include <Qt3DRender/qframegraphnode.h>
#include <Qt3DRender/qt3drender_global.h>
#include <QtCore/QSize>

QT_BEGIN_NAMESPACE

class QSurface;
class QWindow;

namespace Qt3DRender {

class QRenderSurfaceSelectorPrivate;

class Q_3DRENDERSHARED_EXPORT QRenderSurfaceSelector : public Qt3DRender::QFrameGraphNode
{
    Q_OBJECT
    Q_PROPERTY(QObject *surface READ surface WRITE setSurface NOTIFY surfaceChanged)
    Q_PROPERTY(QSize externalRenderTargetSize READ externalRenderTargetSize WRITE setExternalRenderTargetSize NOTIFY externalRenderTargetSizeChanged)
    Q_PROPERTY(float surfacePixelRatio READ surfacePixelRatio WRITE setSurfacePixelRatio NOTIFY surfacePixelRatioChanged)

public:
    explicit QRenderSurfaceSelector(Qt3DCore::QNode *parent = nullptr);
    ~QRenderSurfaceSelector();

    QObject *surface() const;
    QSize externalRenderTargetSize() const;
    float surfacePixelRatio() const;

public Q_SLOTS:
    void setSurface(QObject *surfaceObject);
    void setSurfacePixelRatio(float ratio);
    void setExternalRenderTargetSize(const QSize &size);

Q_SIGNALS:
    void surfaceChanged(QObject *surface);
    void externalRenderTargetSizeChanged(const QSize &size);
    void surfacePixelRatioChanged(float ratio);

protected:
    Q_DECLARE_PRIVATE(QRenderSurfaceSelector)
    explicit QRenderSurfaceSelector(QRenderSurfaceSelectorPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QRENDERSURFACESELECTOR_H

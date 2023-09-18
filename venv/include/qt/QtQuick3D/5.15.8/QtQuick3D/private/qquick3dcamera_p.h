/****************************************************************************
**
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

#ifndef QSSGCAMERA_H
#define QSSGCAMERA_H

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

#include <QtQuick3D/private/qquick3dnode_p.h>

QT_BEGIN_NAMESPACE

struct QSSGRenderCamera;
class Q_QUICK3D_EXPORT QQuick3DCamera : public QQuick3DNode
{
    Q_OBJECT
    Q_PROPERTY(bool frustumCullingEnabled READ frustumCullingEnabled WRITE setFrustumCullingEnabled NOTIFY frustumCullingEnabledChanged)
public:

    enum FieldOfViewOrientation {
        Vertical,
        Horizontal
    };
    Q_ENUM(FieldOfViewOrientation)

    explicit QQuick3DCamera(QQuick3DNode *parent = nullptr);

    Q_INVOKABLE QVector3D mapToViewport(const QVector3D &scenePos) const;
    Q_INVOKABLE QVector3D mapFromViewport(const QVector3D &viewportPos) const;
    QVector3D mapToViewport(const QVector3D &scenePos,
                            qreal width,
                            qreal height);
    QVector3D mapFromViewport(const QVector3D &viewportPos,
                              qreal width,
                              qreal height);

    Q_REVISION(1) Q_INVOKABLE void lookAt(const QVector3D &scenePos);
    Q_REVISION(1) Q_INVOKABLE void lookAt(QQuick3DNode *node);

    QSSGRenderCamera *cameraNode() const;
    void setCameraNode(QSSGRenderCamera *camera) { m_cameraNode = camera; }

    // It will be used only after the scene was drawn.
    // It means that the spatialNode of this camera already was created.
    void updateGlobalVariables(const QRectF &inViewport);

    bool frustumCullingEnabled() const;

public Q_SLOTS:
    void setFrustumCullingEnabled(bool frustumCullingEnabled);

Q_SIGNALS:
    void frustumCullingEnabledChanged();

protected:
    QSSGRenderGraphObject *updateSpatialNode(QSSGRenderGraphObject *node) override;
    virtual bool checkSpatialNode(QSSGRenderCamera *camera) = 0;

private:
    QSSGRenderCamera *m_cameraNode = nullptr;
    bool m_frustumCullingEnabled = false;
};

QT_END_NAMESPACE

#endif // QSSGCAMERA_H

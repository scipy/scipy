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

#ifndef QT3DRENDER_CAMERALENS_H
#define QT3DRENDER_CAMERALENS_H

#include <Qt3DCore/qcomponent.h>
#include <Qt3DRender/qt3drender_global.h>

#include <QtGui/QMatrix4x4>
#include <QtGui/QQuaternion>
#include <QtGui/QVector3D>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QCameraLensPrivate;

class Q_3DRENDERSHARED_EXPORT QCameraLens : public Qt3DCore::QComponent
{
    Q_OBJECT
    Q_PROPERTY(ProjectionType projectionType READ projectionType WRITE setProjectionType NOTIFY projectionTypeChanged)
    Q_PROPERTY(float nearPlane READ nearPlane WRITE setNearPlane NOTIFY nearPlaneChanged)
    Q_PROPERTY(float farPlane READ farPlane WRITE setFarPlane NOTIFY farPlaneChanged)
    Q_PROPERTY(float fieldOfView READ fieldOfView WRITE setFieldOfView NOTIFY fieldOfViewChanged)
    Q_PROPERTY(float aspectRatio READ aspectRatio WRITE setAspectRatio NOTIFY aspectRatioChanged)
    Q_PROPERTY(float left READ left WRITE setLeft NOTIFY leftChanged)
    Q_PROPERTY(float right READ right WRITE setRight NOTIFY rightChanged)
    Q_PROPERTY(float bottom READ bottom WRITE setBottom NOTIFY bottomChanged)
    Q_PROPERTY(float top READ top WRITE setTop NOTIFY topChanged)
    Q_PROPERTY(QMatrix4x4 projectionMatrix READ projectionMatrix WRITE setProjectionMatrix NOTIFY projectionMatrixChanged)
    Q_PROPERTY(float exposure READ exposure WRITE setExposure NOTIFY exposureChanged REVISION 9)

public:
    explicit QCameraLens(QNode *parent = nullptr);
    ~QCameraLens();

    enum ProjectionType {
        OrthographicProjection,
        PerspectiveProjection,
        FrustumProjection,
        CustomProjection
    };
    Q_ENUM(ProjectionType) // LCOV_EXCL_LINE

    ProjectionType projectionType() const;
    float nearPlane() const;
    float farPlane() const;
    float fieldOfView() const;
    float aspectRatio() const;
    float left() const;
    float right() const;
    float bottom() const;
    float top() const;

    QMatrix4x4 projectionMatrix() const;

    void setOrthographicProjection(float left, float right,
                                   float bottom, float top,
                                   float nearPlane, float farPlane);

    void setFrustumProjection(float left, float right,
                              float bottom, float top,
                              float nearPlane, float farPlane);

    void setPerspectiveProjection(float fieldOfView, float aspect,
                                  float nearPlane, float farPlane);

    float exposure() const;

    void viewAll(Qt3DCore::QNodeId cameraId);
    void viewEntity(Qt3DCore::QNodeId entityId, Qt3DCore::QNodeId cameraId);

public Q_SLOTS:
    void setProjectionType(ProjectionType projectionType);
    void setNearPlane(float nearPlane);
    void setFarPlane(float farPlane);
    void setFieldOfView(float fieldOfView);
    void setAspectRatio(float aspectRatio);
    void setLeft(float left);
    void setRight(float right);
    void setBottom(float bottom);
    void setTop(float top);
    void setProjectionMatrix(const QMatrix4x4 &projectionMatrix);
    void setExposure(float exposure);

Q_SIGNALS:
    void projectionTypeChanged(QCameraLens::ProjectionType projectionType);
    void nearPlaneChanged(float nearPlane);
    void farPlaneChanged(float farPlane);
    void fieldOfViewChanged(float fieldOfView);
    void aspectRatioChanged(float aspectRatio);
    void leftChanged(float left);
    void rightChanged(float right);
    void bottomChanged(float bottom);
    void topChanged(float top);
    void projectionMatrixChanged(const QMatrix4x4 &projectionMatrix);
    void exposureChanged(float exposure);
    void viewSphere(const QVector3D &center, float radius);

protected:
    explicit QCameraLens(QCameraLensPrivate &dd, QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QCameraLens)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;
};

} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_CAMERALENS_H

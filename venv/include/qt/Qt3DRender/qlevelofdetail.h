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

#ifndef QT3DRENDER_QLEVELOFDETAIL_H
#define QT3DRENDER_QLEVELOFDETAIL_H

#include <Qt3DCore/qcomponent.h>
#include <Qt3DRender/qt3drender_global.h>
#include <Qt3DRender/qlevelofdetailboundingsphere.h>

#include <QtGui/QVector3D>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QCamera;
class QLevelOfDetailPrivate;

class Q_3DRENDERSHARED_EXPORT QLevelOfDetail : public Qt3DCore::QComponent
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QCamera *camera READ camera WRITE setCamera NOTIFY cameraChanged)
    Q_PROPERTY(int currentIndex READ currentIndex WRITE setCurrentIndex NOTIFY currentIndexChanged)
    Q_PROPERTY(ThresholdType thresholdType READ thresholdType WRITE setThresholdType NOTIFY thresholdTypeChanged)
    Q_PROPERTY(QVector<qreal> thresholds READ thresholds WRITE setThresholds NOTIFY thresholdsChanged)
    Q_PROPERTY(Qt3DRender::QLevelOfDetailBoundingSphere volumeOverride READ volumeOverride WRITE setVolumeOverride NOTIFY volumeOverrideChanged)

public:
    enum ThresholdType {
        DistanceToCameraThreshold,
        ProjectedScreenPixelSizeThreshold,
    };
    Q_ENUM(ThresholdType) // LCOV_EXCL_LINE

    explicit QLevelOfDetail(Qt3DCore::QNode *parent = nullptr);
    ~QLevelOfDetail();

    QCamera *camera() const;
    int currentIndex() const;
    ThresholdType thresholdType() const;
    QVector<qreal> thresholds() const;
    QLevelOfDetailBoundingSphere volumeOverride() const;

    Q_INVOKABLE Qt3DRender::QLevelOfDetailBoundingSphere createBoundingSphere(const QVector3D &center, float radius);

public Q_SLOTS:
    void setCamera(QCamera *camera);
    void setCurrentIndex(int currentIndex);
    void setThresholdType(ThresholdType thresholdType);
    void setThresholds(const QVector<qreal> &thresholds);
    void setVolumeOverride(const QLevelOfDetailBoundingSphere &volumeOverride);

Q_SIGNALS:
    void cameraChanged(QCamera *camera);
    void currentIndexChanged(int currentIndex);
    void thresholdTypeChanged(ThresholdType thresholdType);
    void thresholdsChanged(const QVector<qreal> &thresholds);
    void volumeOverrideChanged(const QLevelOfDetailBoundingSphere &volumeOverride);

protected:
    explicit QLevelOfDetail(QLevelOfDetailPrivate &dd, Qt3DCore::QNode *parent = nullptr);
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
    // TODO Unused remove in Qt6
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

private:
    Q_DECLARE_PRIVATE(QLevelOfDetail)
};

} // namespace Qt3DRender

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DRender::QLevelOfDetailBoundingSphere)

#endif // QT3DRENDER_QLEVELOFDETAIL_H

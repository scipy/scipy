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

#ifndef QT3DEXTRAS_EXTRAS_QUICK_QUICK3DLEVELOFDETAILLOADER_P_H
#define QT3DEXTRAS_EXTRAS_QUICK_QUICK3DLEVELOFDETAILLOADER_P_H

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

#include <Qt3DQuickExtras/private/qt3dquickextras_global_p.h>
#include <Qt3DCore/qentity.h>
#include <Qt3DRender/qlevelofdetail.h>

QT_BEGIN_NAMESPACE

namespace Qt3DExtras {
namespace Extras {
namespace Quick {

class Quick3DLevelOfDetailLoaderPrivate;

class Q_3DQUICKEXTRASSHARED_PRIVATE_EXPORT Quick3DLevelOfDetailLoader : public Qt3DCore::QEntity
{
    Q_OBJECT
    Q_PROPERTY(QVariantList sources READ sources WRITE setSources NOTIFY sourcesChanged)

    Q_PROPERTY(Qt3DRender::QCamera *camera READ camera WRITE setCamera NOTIFY cameraChanged)
    Q_PROPERTY(int currentIndex READ currentIndex WRITE setCurrentIndex NOTIFY currentIndexChanged)
    Q_PROPERTY(Qt3DRender::QLevelOfDetail::ThresholdType thresholdType READ thresholdType WRITE setThresholdType NOTIFY thresholdTypeChanged)
    Q_PROPERTY(QVector<qreal> thresholds READ thresholds WRITE setThresholds NOTIFY thresholdsChanged)
    Q_PROPERTY(Qt3DRender::QLevelOfDetailBoundingSphere volumeOverride READ volumeOverride WRITE setVolumeOverride NOTIFY volumeOverrideChanged)

    Q_PROPERTY(QObject *entity READ entity NOTIFY entityChanged)
    Q_PROPERTY(QUrl source READ source NOTIFY sourceChanged)
public:
    explicit Quick3DLevelOfDetailLoader(QNode *parent = 0);

    QVariantList sources() const;
    void setSources(const QVariantList &sources);

    Qt3DRender::QCamera *camera() const;
    void setCamera(Qt3DRender::QCamera *camera);
    int currentIndex() const;
    void setCurrentIndex(int currentIndex);
    Qt3DRender::QLevelOfDetail::ThresholdType thresholdType() const;
    void setThresholdType(Qt3DRender::QLevelOfDetail::ThresholdType thresholdType);
    QVector<qreal> thresholds() const;
    void setThresholds(const QVector<qreal> &thresholds);
    Qt3DRender::QLevelOfDetailBoundingSphere volumeOverride() const;
    void setVolumeOverride(const Qt3DRender::QLevelOfDetailBoundingSphere &volumeOverride);

    Q_INVOKABLE Qt3DRender::QLevelOfDetailBoundingSphere createBoundingSphere(const QVector3D &center, float radius);

    QObject *entity() const;
    QUrl source() const;

Q_SIGNALS:
    void sourcesChanged();
    void cameraChanged();
    void currentIndexChanged();
    void thresholdTypeChanged();
    void thresholdsChanged();
    void volumeOverrideChanged();
    void entityChanged();
    void sourceChanged();

private:
    Q_DECLARE_PRIVATE(Quick3DLevelOfDetailLoader)
};

} // namespace Quick
} // namespace Extras
} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_EXTRAS_QUICK_QUICK3DLEVELOFDETAILLOADER_P_H

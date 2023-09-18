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

#ifndef QT3DRENDER_QENVIRONMENTLIGHT_H
#define QT3DRENDER_QENVIRONMENTLIGHT_H

#include <Qt3DRender/qt3drender_global.h>
#include <Qt3DCore/qcomponent.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAbstractTexture;
class QEnvironmentLightPrivate;

class Q_3DRENDERSHARED_EXPORT QEnvironmentLight : public Qt3DCore::QComponent
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QAbstractTexture *irradiance READ irradiance WRITE setIrradiance NOTIFY irradianceChanged)
    Q_PROPERTY(Qt3DRender::QAbstractTexture *specular READ specular WRITE setSpecular NOTIFY specularChanged)

public:
    explicit QEnvironmentLight(Qt3DCore::QNode *parent = nullptr);
    ~QEnvironmentLight();

    Qt3DRender::QAbstractTexture *irradiance() const;
    Qt3DRender::QAbstractTexture *specular() const;

public Q_SLOTS:
    void setIrradiance(Qt3DRender::QAbstractTexture *irradiance);
    void setSpecular(Qt3DRender::QAbstractTexture *specular);

protected:
    explicit QEnvironmentLight(QEnvironmentLightPrivate &dd, Qt3DCore::QNode *parent = nullptr);

Q_SIGNALS:
    void irradianceChanged(Qt3DRender::QAbstractTexture *environmentIrradiance);
    void specularChanged(Qt3DRender::QAbstractTexture *environmentSpecular);

private:
    Q_DECLARE_PRIVATE(QEnvironmentLight)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;

    Q_PRIVATE_SLOT(d_func(), void _q_updateEnvMapsSize())
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QENVIRONMENTLIGHT_H

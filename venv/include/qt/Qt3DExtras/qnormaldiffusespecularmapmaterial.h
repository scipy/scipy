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

#ifndef QT3DEXTRAS_QNORMALDIFFUSESPECULARMAPMATERIAL_H
#define QT3DEXTRAS_QNORMALDIFFUSESPECULARMAPMATERIAL_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qmaterial.h>
#include <QtGui/QColor>

QT_BEGIN_NAMESPACE

namespace Qt3DExtras {

class QNormalDiffuseSpecularMapMaterialPrivate;

class Q_3DEXTRASSHARED_EXPORT QNormalDiffuseSpecularMapMaterial : public Qt3DRender::QMaterial
{
    Q_OBJECT
    Q_PROPERTY(QColor ambient READ ambient WRITE setAmbient NOTIFY ambientChanged)
    Q_PROPERTY(Qt3DRender::QAbstractTexture *diffuse READ diffuse WRITE setDiffuse NOTIFY diffuseChanged)
    Q_PROPERTY(Qt3DRender::QAbstractTexture *normal READ normal WRITE setNormal NOTIFY normalChanged)
    Q_PROPERTY(Qt3DRender::QAbstractTexture *specular READ specular WRITE setSpecular NOTIFY specularChanged)
    Q_PROPERTY(float shininess READ shininess WRITE setShininess NOTIFY shininessChanged)
    Q_PROPERTY(float textureScale READ textureScale WRITE setTextureScale NOTIFY textureScaleChanged)

public:
    explicit QNormalDiffuseSpecularMapMaterial(Qt3DCore::QNode *parent = nullptr);
    ~QNormalDiffuseSpecularMapMaterial();

    QColor ambient() const;
    Qt3DRender::QAbstractTexture *diffuse() const;
    Qt3DRender::QAbstractTexture *normal() const;
    Qt3DRender::QAbstractTexture *specular() const;
    float shininess() const;
    float textureScale() const;

public Q_SLOTS:
    void setAmbient(const QColor &ambient);
    void setDiffuse(Qt3DRender::QAbstractTexture *diffuse);
    void setNormal(Qt3DRender::QAbstractTexture *normal);
    void setSpecular(Qt3DRender::QAbstractTexture *specular);
    void setShininess(float shininess);
    void setTextureScale(float textureScale);

Q_SIGNALS:
    void ambientChanged(const QColor &ambient);
    void diffuseChanged(Qt3DRender::QAbstractTexture *diffuse);
    void normalChanged(Qt3DRender::QAbstractTexture *normal);
    void specularChanged(Qt3DRender::QAbstractTexture *specular);
    void shininessChanged(float shininess);
    void textureScaleChanged(float textureScale);

protected:
    QNormalDiffuseSpecularMapMaterial(QNormalDiffuseSpecularMapMaterialPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QNormalDiffuseSpecularMapMaterial)
};

} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QNORMALDIFFUSESPECULARMAPMATERIAL_H

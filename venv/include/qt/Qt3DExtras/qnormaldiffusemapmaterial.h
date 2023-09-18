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

#ifndef QT3DEXTRAS_QNORMALDIFFUSEMAPMATERIAL_H
#define QT3DEXTRAS_QNORMALDIFFUSEMAPMATERIAL_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qmaterial.h>
#include <QtGui/QColor>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAbstractTexture;

} // namespace Qt3DRender

namespace Qt3DExtras {

class QNormalDiffuseMapMaterialPrivate;

class Q_3DEXTRASSHARED_EXPORT QNormalDiffuseMapMaterial : public Qt3DRender::QMaterial
{
    Q_OBJECT
    Q_PROPERTY(QColor ambient READ ambient WRITE setAmbient NOTIFY ambientChanged)
    Q_PROPERTY(QColor specular READ specular WRITE setSpecular NOTIFY specularChanged)
    Q_PROPERTY(Qt3DRender::QAbstractTexture *diffuse READ diffuse WRITE setDiffuse NOTIFY diffuseChanged)
    Q_PROPERTY(Qt3DRender::QAbstractTexture *normal READ normal WRITE setNormal NOTIFY normalChanged)
    Q_PROPERTY(float shininess READ shininess WRITE setShininess NOTIFY shininessChanged)
    Q_PROPERTY(float textureScale READ textureScale WRITE setTextureScale NOTIFY textureScaleChanged)

public:
    explicit QNormalDiffuseMapMaterial(Qt3DCore::QNode *parent = nullptr);
    ~QNormalDiffuseMapMaterial();

    QColor ambient() const;
    QColor specular() const;
    Qt3DRender::QAbstractTexture *diffuse() const;
    Qt3DRender::QAbstractTexture *normal() const;
    float shininess() const;
    float textureScale() const;

public Q_SLOTS:
    void setAmbient(const QColor &ambient);
    void setSpecular(const QColor &specular);
    void setDiffuse(Qt3DRender::QAbstractTexture *diffuse);
    void setNormal(Qt3DRender::QAbstractTexture *normal);
    void setShininess(float shininess);
    void setTextureScale(float textureScale);

Q_SIGNALS:
    void ambientChanged(const QColor &ambient);
    void diffuseChanged(Qt3DRender::QAbstractTexture *diffuse);
    void normalChanged(Qt3DRender::QAbstractTexture *normal);
    void specularChanged(const QColor &specular);
    void shininessChanged(float shininess);
    void textureScaleChanged(float textureScale);

protected:
    QNormalDiffuseMapMaterial(QNormalDiffuseMapMaterialPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QNormalDiffuseMapMaterial)
};

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QNORMALDIFFUSEMAPMATERIAL_H

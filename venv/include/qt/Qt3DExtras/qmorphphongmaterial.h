/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DEXTRAS_QMORPHPHONGMATERIAL_H
#define QT3DEXTRAS_QMORPHPHONGMATERIAL_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qmaterial.h>
#include <QtGui/QColor>

QT_BEGIN_NAMESPACE

namespace Qt3DExtras {

class QMorphPhongMaterialPrivate;

class Q_3DEXTRASSHARED_EXPORT QMorphPhongMaterial : public Qt3DRender::QMaterial
{
    Q_OBJECT
    Q_PROPERTY(QColor ambient READ ambient WRITE setAmbient NOTIFY ambientChanged)
    Q_PROPERTY(QColor diffuse READ diffuse WRITE setDiffuse NOTIFY diffuseChanged)
    Q_PROPERTY(QColor specular READ specular WRITE setSpecular NOTIFY specularChanged)
    Q_PROPERTY(float shininess READ shininess WRITE setShininess NOTIFY shininessChanged)
    Q_PROPERTY(float interpolator READ interpolator WRITE setInterpolator NOTIFY interpolatorChanged)

public:
    explicit QMorphPhongMaterial(Qt3DCore::QNode *parent = nullptr);
    ~QMorphPhongMaterial();

    QColor ambient() const;
    QColor diffuse() const;
    QColor specular() const;
    float shininess() const;
    float interpolator() const;

public Q_SLOTS:
    void setAmbient(const QColor &ambient);
    void setDiffuse(const QColor &diffuse);
    void setSpecular(const QColor &specular);
    void setShininess(float shininess);
    void setInterpolator(float interpolator);

Q_SIGNALS:
    void ambientChanged(const QColor &ambient);
    void diffuseChanged(const QColor &diffuse);
    void specularChanged(const QColor &specular);
    void shininessChanged(float shininess);
    void interpolatorChanged(float interpolator);

private:
    Q_DECLARE_PRIVATE(QMorphPhongMaterial)
};

} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QMORPHPHONGMATERIAL_H


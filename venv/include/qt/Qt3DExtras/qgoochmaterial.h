/****************************************************************************
**
** Copyright (C) 2015 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DEXTRAS_QGOOCHMATERIAL_H
#define QT3DEXTRAS_QGOOCHMATERIAL_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qmaterial.h>
#include <QtGui/QColor>

QT_BEGIN_NAMESPACE

namespace Qt3DExtras {

class QGoochMaterialPrivate;

class Q_3DEXTRASSHARED_EXPORT QGoochMaterial : public Qt3DRender::QMaterial
{
    Q_OBJECT
    Q_PROPERTY(QColor diffuse READ diffuse WRITE setDiffuse NOTIFY diffuseChanged)
    Q_PROPERTY(QColor specular READ specular WRITE setSpecular NOTIFY specularChanged)
    Q_PROPERTY(QColor cool READ cool WRITE setCool NOTIFY coolChanged)
    Q_PROPERTY(QColor warm READ warm WRITE setWarm NOTIFY warmChanged)
    Q_PROPERTY(float alpha READ alpha WRITE setAlpha NOTIFY alphaChanged)
    Q_PROPERTY(float beta READ beta WRITE setBeta NOTIFY betaChanged)
    Q_PROPERTY(float shininess READ shininess WRITE setShininess NOTIFY shininessChanged)

public:
    explicit QGoochMaterial(Qt3DCore::QNode *parent = nullptr);
    ~QGoochMaterial();

    QColor diffuse() const;
    QColor specular() const;
    QColor cool() const;
    QColor warm() const;
    float alpha() const;
    float beta() const;
    float shininess() const;

public Q_SLOTS:
    void setDiffuse(const QColor &diffuse);
    void setSpecular(const QColor &specular);
    void setCool(const QColor &cool);
    void setWarm(const QColor &warm);
    void setAlpha(float alpha);
    void setBeta(float beta);
    void setShininess(float shininess);

Q_SIGNALS:
    void diffuseChanged(const QColor &diffuse);
    void specularChanged(const QColor &specular);
    void coolChanged(const QColor &cool);
    void warmChanged(const QColor &warm);
    void alphaChanged(float alpha);
    void betaChanged(float beta);
    void shininessChanged(float shininess);

protected:
    explicit QGoochMaterial(QGoochMaterialPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QGoochMaterial)
};

} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QGOOCHMATERIAL_H

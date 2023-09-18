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

#ifndef QT3DEXTRAS_QPHONGALPHAMATERIAL_H
#define QT3DEXTRAS_QPHONGALPHAMATERIAL_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qblendequation.h>
#include <Qt3DRender/qblendequationarguments.h>
#include <Qt3DRender/qmaterial.h>
#include <QtGui/QColor>

QT_BEGIN_NAMESPACE

namespace Qt3DExtras {

class QPhongAlphaMaterialPrivate;

class Q_3DEXTRASSHARED_EXPORT QPhongAlphaMaterial : public Qt3DRender::QMaterial
{
    Q_OBJECT
    Q_PROPERTY(QColor ambient READ ambient WRITE setAmbient NOTIFY ambientChanged)
    Q_PROPERTY(QColor diffuse READ diffuse WRITE setDiffuse NOTIFY diffuseChanged)
    Q_PROPERTY(QColor specular READ specular WRITE setSpecular NOTIFY specularChanged)
    Q_PROPERTY(float shininess READ shininess WRITE setShininess NOTIFY shininessChanged)
    Q_PROPERTY(float alpha READ alpha WRITE setAlpha NOTIFY alphaChanged)
    Q_PROPERTY(Qt3DRender::QBlendEquationArguments::Blending sourceRgbArg READ sourceRgbArg WRITE setSourceRgbArg NOTIFY sourceRgbArgChanged)
    Q_PROPERTY(Qt3DRender::QBlendEquationArguments::Blending destinationRgbArg READ destinationRgbArg WRITE setDestinationRgbArg NOTIFY destinationRgbArgChanged)
    Q_PROPERTY(Qt3DRender::QBlendEquationArguments::Blending sourceAlphaArg READ sourceAlphaArg WRITE setSourceAlphaArg NOTIFY sourceAlphaArgChanged)
    Q_PROPERTY(Qt3DRender::QBlendEquationArguments::Blending destinationAlphaArg READ destinationAlphaArg WRITE setDestinationAlphaArg NOTIFY destinationAlphaArgChanged)
    Q_PROPERTY(Qt3DRender::QBlendEquation::BlendFunction blendFunctionArg READ blendFunctionArg WRITE setBlendFunctionArg NOTIFY blendFunctionArgChanged)

public:
    explicit QPhongAlphaMaterial(Qt3DCore::QNode *parent = nullptr);
    ~QPhongAlphaMaterial();

    QColor ambient() const;
    QColor diffuse() const;
    QColor specular() const;
    float shininess() const;
    float alpha() const;

    Qt3DRender::QBlendEquationArguments::Blending sourceRgbArg() const;
    Qt3DRender::QBlendEquationArguments::Blending destinationRgbArg() const;
    Qt3DRender::QBlendEquationArguments::Blending sourceAlphaArg() const;
    Qt3DRender::QBlendEquationArguments::Blending destinationAlphaArg() const;
    Qt3DRender::QBlendEquation::BlendFunction blendFunctionArg() const;

public Q_SLOTS:
    void setAmbient(const QColor &ambient);
    void setDiffuse(const QColor &diffuse);
    void setSpecular(const QColor &specular);
    void setShininess(float shininess);
    void setAlpha(float alpha);
    void setSourceRgbArg(Qt3DRender::QBlendEquationArguments::Blending sourceRgbArg);
    void setDestinationRgbArg(Qt3DRender::QBlendEquationArguments::Blending destinationRgbArg);
    void setSourceAlphaArg(Qt3DRender::QBlendEquationArguments::Blending sourceAlphaArg);
    void setDestinationAlphaArg(Qt3DRender::QBlendEquationArguments::Blending destinationAlphaArg);
    void setBlendFunctionArg(Qt3DRender::QBlendEquation::BlendFunction blendFunctionArg);

Q_SIGNALS:
    void ambientChanged(const QColor &ambient);
    void diffuseChanged(const QColor &diffuse);
    void specularChanged(const QColor &specular);
    void shininessChanged(float shininess);
    void alphaChanged(float alpha);
    void sourceRgbArgChanged(Qt3DRender::QBlendEquationArguments::Blending sourceRgbArg);
    void destinationRgbArgChanged(Qt3DRender::QBlendEquationArguments::Blending destinationRgbArg);
    void sourceAlphaArgChanged(Qt3DRender::QBlendEquationArguments::Blending sourceAlphaArg);
    void destinationAlphaArgChanged(Qt3DRender::QBlendEquationArguments::Blending destinationAlphaArg);
    void blendFunctionArgChanged(Qt3DRender::QBlendEquation::BlendFunction blendFunctionArg);

private:
    Q_DECLARE_PRIVATE(QPhongAlphaMaterial)
};

} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QPHONGALPHAMATERIAL_H

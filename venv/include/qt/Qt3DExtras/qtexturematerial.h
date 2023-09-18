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

#ifndef QT3DEXTRAS_QTEXTUREMATERIAL_H
#define QT3DEXTRAS_QTEXTUREMATERIAL_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qmaterial.h>
#include <QtGui/QVector2D>
#include <QtGui/QMatrix3x3>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAbstractTexture;

} // namespace Qt3DRender

namespace Qt3DExtras {

class QTextureMaterialPrivate;

class Q_3DEXTRASSHARED_EXPORT QTextureMaterial : public Qt3DRender::QMaterial
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QAbstractTexture *texture READ texture WRITE setTexture NOTIFY textureChanged)
    Q_PROPERTY(QVector2D textureOffset READ textureOffset WRITE setTextureOffset NOTIFY textureOffsetChanged)
    Q_PROPERTY(QMatrix3x3 textureTransform READ textureTransform WRITE setTextureTransform NOTIFY textureTransformChanged REVISION 10)
    Q_PROPERTY(bool alphaBlending READ isAlphaBlendingEnabled WRITE setAlphaBlendingEnabled NOTIFY alphaBlendingEnabledChanged REVISION 11)

public:
    explicit QTextureMaterial(Qt3DCore::QNode *parent = nullptr);
    ~QTextureMaterial();

    Qt3DRender::QAbstractTexture *texture() const;
    QVector2D textureOffset() const;
    QMatrix3x3 textureTransform() const;
    bool isAlphaBlendingEnabled() const;

public Q_SLOTS:
    void setTexture(Qt3DRender::QAbstractTexture *texture);
    void setTextureOffset(QVector2D textureOffset);
    void setTextureTransform(const QMatrix3x3 &matrix);
    void setAlphaBlendingEnabled(bool enabled);

Q_SIGNALS:
    void textureChanged(Qt3DRender::QAbstractTexture *texture);
    void textureOffsetChanged(QVector2D textureOffset);
    void textureTransformChanged(const QMatrix3x3 &textureTransform);
    void alphaBlendingEnabledChanged(bool enabled);

private:
    Q_DECLARE_PRIVATE(QTextureMaterial)
};

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QTEXTUREMATERIAL_H

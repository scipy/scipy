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

#ifndef QT3DEXTRAS_QABSTRACTSPRITESHEET_H
#define QT3DEXTRAS_QABSTRACTSPRITESHEET_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DCore/QNode>
#include <QtGui/QVector2D>
#include <QtGui/QMatrix3x3>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAbstractTexture;

} // namespace Qt3DRender

namespace Qt3DExtras {

class QAbstractSpriteSheetPrivate;

class Q_3DEXTRASSHARED_EXPORT QAbstractSpriteSheet : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QAbstractTexture *texture READ texture WRITE setTexture NOTIFY textureChanged)
    Q_PROPERTY(QMatrix3x3 textureTransform READ textureTransform NOTIFY textureTransformChanged)
    Q_PROPERTY(int currentIndex READ currentIndex WRITE setCurrentIndex NOTIFY currentIndexChanged)
public:
    ~QAbstractSpriteSheet();

    Qt3DRender::QAbstractTexture *texture() const;
    QMatrix3x3 textureTransform() const;
    int currentIndex() const;

public Q_SLOTS:
    void setTexture(Qt3DRender::QAbstractTexture *texture);
    void setCurrentIndex(int currentIndex);

Q_SIGNALS:
    void textureChanged(Qt3DRender::QAbstractTexture *texture);
    void textureTransformChanged(const QMatrix3x3 &textureTransform);
    void currentIndexChanged(int currentIndex);

protected:
    explicit QAbstractSpriteSheet(QAbstractSpriteSheetPrivate &d, Qt3DCore::QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QAbstractSpriteSheet)
};

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QABSTRACTSPRITESHEET_H

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

#ifndef QT3DCORE_QJOINT_H
#define QT3DCORE_QJOINT_H

#include <Qt3DCore/qnode.h>
#include <Qt3DCore/qt3dcore_global.h>

#include <QtGui/qmatrix4x4.h>
#include <QtGui/qquaternion.h>
#include <QtGui/qvector3d.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QJointPrivate;

class Q_3DCORESHARED_EXPORT QJoint : public QNode
{
    Q_OBJECT
    Q_PROPERTY(QVector3D scale READ scale WRITE setScale NOTIFY scaleChanged)
    Q_PROPERTY(QQuaternion rotation READ rotation WRITE setRotation NOTIFY rotationChanged)
    Q_PROPERTY(QVector3D translation READ translation WRITE setTranslation NOTIFY translationChanged)
    Q_PROPERTY(QMatrix4x4 inverseBindMatrix READ inverseBindMatrix WRITE setInverseBindMatrix NOTIFY inverseBindMatrixChanged)
    Q_PROPERTY(float rotationX READ rotationX WRITE setRotationX NOTIFY rotationXChanged)
    Q_PROPERTY(float rotationY READ rotationY WRITE setRotationY NOTIFY rotationYChanged)
    Q_PROPERTY(float rotationZ READ rotationZ WRITE setRotationZ NOTIFY rotationZChanged)
    Q_PROPERTY(QString name READ name WRITE setName NOTIFY nameChanged)

public:
    explicit QJoint(Qt3DCore::QNode *parent = nullptr);
    ~QJoint();

    QVector3D scale() const;
    QQuaternion rotation() const;
    QVector3D translation() const;
    QMatrix4x4 inverseBindMatrix() const;
    float rotationX() const;
    float rotationY() const;
    float rotationZ() const;
    QString name() const;

    void addChildJoint(QJoint *joint);
    void removeChildJoint(QJoint *joint);
    QVector<QJoint *> childJoints() const;

public Q_SLOTS:
    void setScale(const QVector3D &scale);
    void setRotation(const QQuaternion &rotation);
    void setTranslation(const QVector3D &translation);
    void setInverseBindMatrix(const QMatrix4x4 &inverseBindMatrix);
    void setRotationX(float rotationX);
    void setRotationY(float rotationY);
    void setRotationZ(float rotationZ);
    void setName(const QString &name);
    void setToIdentity();

Q_SIGNALS:
    void scaleChanged(const QVector3D &scale);
    void rotationChanged(const QQuaternion &rotation);
    void translationChanged(const QVector3D &translation);
    void inverseBindMatrixChanged(const QMatrix4x4 &inverseBindMatrix);
    void rotationXChanged(float rotationX);
    void rotationYChanged(float rotationY);
    void rotationZChanged(float rotationZ);
    void nameChanged(const QString &name);

private:
    Q_DECLARE_PRIVATE(QJoint)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QJOINT_H

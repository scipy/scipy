/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QQUICK3DQUATERNIONUTILS_H
#define QQUICK3DQUATERNIONUTILS_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtQuick3D/private/qtquick3dglobal_p.h>

#include <QtCore/QObject>
#include <QtGui/QQuaternion>
#include <QtGui/QVector3D>

QT_BEGIN_NAMESPACE

class Q_QUICK3D_EXPORT QQuick3DQuaternionUtils : public QObject
{
    Q_OBJECT
public:
    explicit QQuick3DQuaternionUtils(QObject *parent = nullptr);

    Q_INVOKABLE static QQuaternion fromAxesAndAngles(const QVector3D &axis1,
                                                     float angle1,
                                                     const QVector3D &axis2,
                                                     float angle2,
                                                     const QVector3D &axis3,
                                                     float angle3);
    Q_INVOKABLE static QQuaternion fromAxesAndAngles(const QVector3D &axis1,
                                                     float angle1,
                                                     const QVector3D &axis2,
                                                     float angle2);
    Q_INVOKABLE static QQuaternion fromAxisAndAngle(float x, float y, float z, float angle);
    Q_INVOKABLE static QQuaternion fromAxisAndAngle(const QVector3D &axis, float angle);
    Q_INVOKABLE static QQuaternion fromEulerAngles(float x, float y, float z);
    Q_INVOKABLE static QQuaternion fromEulerAngles(const QVector3D &eulerAngles);

    Q_REVISION(1) Q_INVOKABLE static QQuaternion lookAt(const QVector3D &sourcePosition,
                                                        const QVector3D &targetPosition,
                                                        const QVector3D &forwardDirection
                                                        = QVector3D(0, 0, -1),
                                                        const QVector3D &upDirection
                                                        = QVector3D(0, 1, 0));

};

QT_END_NAMESPACE

#endif // QQUICK3DQUATERNIONUTILS_H

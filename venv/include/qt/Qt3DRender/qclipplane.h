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

#ifndef QT3DRENDER_QCLIPPLANE_H
#define QT3DRENDER_QCLIPPLANE_H

#include <Qt3DRender/qrenderstate.h>
#include <QtGui/qvector3d.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QClipPlanePrivate;

class Q_3DRENDERSHARED_EXPORT QClipPlane : public QRenderState
{
    Q_OBJECT
    Q_PROPERTY(int planeIndex READ planeIndex WRITE setPlaneIndex NOTIFY planeIndexChanged)
    Q_PROPERTY(QVector3D normal READ normal WRITE setNormal NOTIFY normalChanged)
    Q_PROPERTY(float distance READ distance WRITE setDistance NOTIFY distanceChanged)
public:
    explicit QClipPlane(Qt3DCore::QNode *parent = nullptr);
    ~QClipPlane();

    int planeIndex() const;
    QVector3D normal() const;
    float distance() const;

public Q_SLOTS:
    void setPlaneIndex(int);
    void setNormal(QVector3D);
    void setDistance(float);

Q_SIGNALS:
    void planeIndexChanged(int planeIndex);
    void normalChanged(QVector3D normal);
    void distanceChanged(float distance);

private:
    Q_DECLARE_PRIVATE(QClipPlane)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QCLIPPLANE_H

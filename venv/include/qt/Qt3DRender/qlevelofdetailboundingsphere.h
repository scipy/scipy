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

#ifndef QT3DRENDER_QLEVELOFDETAILBOUNDINGSPHERE_H
#define QT3DRENDER_QLEVELOFDETAILBOUNDINGSPHERE_H

#include <Qt3DCore/qcomponent.h>
#include <Qt3DRender/qt3drender_global.h>

#include <QtGui/QVector3D>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QLevelOfDetailBoundingSpherePrivate;

class Q_3DRENDERSHARED_EXPORT QLevelOfDetailBoundingSphere
{
    Q_GADGET
    Q_PROPERTY(QVector3D center READ center CONSTANT)
    Q_PROPERTY(float radius READ radius CONSTANT)
public:
    explicit QLevelOfDetailBoundingSphere(QVector3D center = QVector3D(), float radius = 1.0f);
    QLevelOfDetailBoundingSphere(const QLevelOfDetailBoundingSphere &other);
    ~QLevelOfDetailBoundingSphere();

    QLevelOfDetailBoundingSphere &operator =(const QLevelOfDetailBoundingSphere &other);

    QVector3D center() const;
    float radius() const;

    bool isEmpty() const;
    bool operator ==(const QLevelOfDetailBoundingSphere &other) const;
    bool operator !=(const QLevelOfDetailBoundingSphere &other) const;

private:
    QSharedDataPointer<QLevelOfDetailBoundingSpherePrivate> d_ptr;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QLEVELOFDETAILBOUNDINGSPHERE_H

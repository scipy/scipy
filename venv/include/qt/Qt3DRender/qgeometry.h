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

#ifndef QT3DRENDER_QGEOMETRY_H
#define QT3DRENDER_QGEOMETRY_H

#include <Qt3DCore/qnode.h>
#include <Qt3DRender/qt3drender_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAttribute;
class QGeometryPrivate;

class Q_3DRENDERSHARED_EXPORT QGeometry : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QAttribute *boundingVolumePositionAttribute READ boundingVolumePositionAttribute WRITE setBoundingVolumePositionAttribute NOTIFY boundingVolumePositionAttributeChanged)
    Q_PROPERTY(QVector3D minExtent READ minExtent NOTIFY minExtentChanged REVISION 13)
    Q_PROPERTY(QVector3D maxExtent READ maxExtent NOTIFY maxExtentChanged REVISION 13)
public:
    explicit QGeometry(Qt3DCore::QNode *parent = nullptr);
    ~QGeometry();

    QVector<QAttribute *> attributes() const;
    Q_INVOKABLE void addAttribute(Qt3DRender::QAttribute *attribute);
    Q_INVOKABLE void removeAttribute(Qt3DRender::QAttribute *attribute);

    QAttribute *boundingVolumePositionAttribute() const;
    QVector3D minExtent() const;
    QVector3D maxExtent() const;

public Q_SLOTS:
    void setBoundingVolumePositionAttribute(QAttribute *boundingVolumePositionAttribute);

Q_SIGNALS:
    void boundingVolumePositionAttributeChanged(QAttribute *boundingVolumePositionAttribute);
    Q_REVISION(13) void minExtentChanged(const QVector3D &minExtent);
    Q_REVISION(13) void maxExtentChanged(const QVector3D &maxExtent);
protected:
    explicit QGeometry(QGeometryPrivate &dd, Qt3DCore::QNode *parent = nullptr);
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

private:
    Q_DECLARE_PRIVATE(QGeometry)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QGEOMETRY_H

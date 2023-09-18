/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
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

#ifndef Q_QUICK3D_GEOMETRY_P_H
#define Q_QUICK3D_GEOMETRY_P_H

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

#include <QtQuick3D/qquick3dgeometry.h>
#include <QtQuick3D/private/qquick3dobject_p.h>

#include <QtGui/qvector3d.h>

QT_BEGIN_NAMESPACE

class QQuick3DGeometryPrivate : public QQuick3DObjectPrivate
{
public:
    QQuick3DGeometryPrivate();
    static const int MAX_ATTRIBUTE_COUNT = 16;
    QString m_name;
    QByteArray m_vertexBuffer;
    QByteArray m_indexBuffer;
    QQuick3DGeometry::Attribute m_attributes[MAX_ATTRIBUTE_COUNT];
    int m_attributeCount = 0;
    QQuick3DGeometry::PrimitiveType m_primitiveType = QQuick3DGeometry::PrimitiveType::Unknown;
    QVector3D m_min;
    QVector3D m_max;
    int m_stride = 0;
    bool m_nameChanged = true;
    bool m_geometryChanged = true;
    bool m_geometryBoundsChanged = true;
};

QT_END_NAMESPACE

#endif

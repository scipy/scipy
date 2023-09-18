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

#ifndef QT3DEXTRAS_QCYLINDERMESH_H
#define QT3DEXTRAS_QCYLINDERMESH_H

#include <Qt3DExtras/qt3dextras_global.h>
#include <Qt3DRender/qgeometryrenderer.h>

QT_BEGIN_NAMESPACE


namespace Qt3DExtras {

class Q_3DEXTRASSHARED_EXPORT QCylinderMesh : public Qt3DRender::QGeometryRenderer
{
    Q_OBJECT
    Q_PROPERTY(int rings READ rings WRITE setRings NOTIFY ringsChanged)
    Q_PROPERTY(int slices READ slices WRITE setSlices NOTIFY slicesChanged)
    Q_PROPERTY(float radius READ radius WRITE setRadius NOTIFY radiusChanged)
    Q_PROPERTY(float length READ length WRITE setLength NOTIFY lengthChanged)
public:
    explicit QCylinderMesh(Qt3DCore::QNode *parent = nullptr);
    ~QCylinderMesh();

    int rings() const;
    int slices() const;
    float radius() const;
    float length() const;

public Q_SLOTS:
    void setRings(int rings);
    void setSlices(int slices);
    void setRadius(float radius);
    void setLength(float length);

Q_SIGNALS:
    void radiusChanged(float radius);
    void ringsChanged(int rings);
    void slicesChanged(int slices);
    void lengthChanged(float length);

private:
    // As this is a default provided geometry renderer, no one should be able
    // to modify the QGeometryRenderer's properties

    void setInstanceCount(int instanceCount);
    void setVertexCount(int vertexCount);
    void setIndexOffset(int indexOffset);
    void setFirstInstance(int firstInstance);
    void setRestartIndexValue(int index);
    void setPrimitiveRestartEnabled(bool enabled);
    void setGeometry(Qt3DRender::QGeometry *geometry);
    void setPrimitiveType(PrimitiveType primitiveType);
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QCYLINDERMESH_H

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

#ifndef Q_QUICK3D_GEOMETRY_H
#define Q_QUICK3D_GEOMETRY_H

#include <QtQuick3D/qquick3dobject.h>

QT_BEGIN_NAMESPACE

class QQuick3DGeometryPrivate;

class Q_QUICK3D_EXPORT QQuick3DGeometry : public QQuick3DObject
{
    Q_OBJECT
    Q_PROPERTY(QString name READ name WRITE setName NOTIFY nameChanged)
    Q_DECLARE_PRIVATE(QQuick3DGeometry)

public:
    explicit QQuick3DGeometry(QQuick3DObject *parent = nullptr);
    ~QQuick3DGeometry() override;

    enum class PrimitiveType {
        Unknown = 0,
        Points,
        LineStrip,
        Lines = 4,
        TriangleStrip,
        TriangleFan,
        Triangles // Default primitive type
    };

    struct Attribute {
        enum Semantic {
            UnknownSemantic = 0,
            IndexSemantic,
            PositionSemantic, // attr_pos
            NormalSemantic,   // attr_norm
            TexCoordSemantic, // attr_uv0
            TangentSemantic,  // attr_textan
            BinormalSemantic  // attr_binormal
        };
        enum ComponentType {
            DefaultType = 0,
            U16Type = 3,
            U32Type = 5, // Default for IndexSemantic
            F32Type = 10 // Default for other semantics
        };
        Semantic semantic = PositionSemantic;
        int offset = -1;
        ComponentType componentType = DefaultType;
    };

    QString name() const;
    QByteArray vertexBuffer() const;
    QByteArray indexBuffer() const;
    int attributeCount() const;
    Attribute attribute(int index) const;
    PrimitiveType primitiveType() const;
    QVector3D boundsMin() const;
    QVector3D boundsMax() const;
    int stride() const;

    void setVertexData(const QByteArray &data);
    void setIndexData(const QByteArray &data);
    void setStride(int stride);
    void setBounds(const QVector3D &min, const QVector3D &max);
    void setPrimitiveType(PrimitiveType type);

    void addAttribute(Attribute::Semantic semantic, int offset,
                      Attribute::ComponentType componentType);
    void addAttribute(const Attribute &att);

    void clear();

public Q_SLOTS:
    void setName(const QString &name);

Q_SIGNALS:
    void nameChanged();
    void geometryNodeDirty();

protected:
    QSSGRenderGraphObject *updateSpatialNode(QSSGRenderGraphObject *node) override;
    void markAllDirty() override;
};

QT_END_NAMESPACE

#endif // Q_QUICK3D_GEOMETRY_H

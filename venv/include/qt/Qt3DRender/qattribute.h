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

#ifndef QT3DRENDER_QATTRIBUTE_H
#define QT3DRENDER_QATTRIBUTE_H

#include <Qt3DRender/qt3drender_global.h>
#include <Qt3DCore/qnode.h>
#include <QtCore/QSharedPointer>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QBuffer;
class QAttributePrivate;

typedef QSharedPointer<QBuffer> QBufferPtr;

class Q_3DRENDERSHARED_EXPORT QAttribute : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QBuffer *buffer READ buffer WRITE setBuffer NOTIFY bufferChanged)
    Q_PROPERTY(QString name READ name WRITE setName NOTIFY nameChanged)
    Q_PROPERTY(VertexBaseType vertexBaseType READ vertexBaseType WRITE setVertexBaseType NOTIFY vertexBaseTypeChanged)
    Q_PROPERTY(uint vertexSize READ vertexSize WRITE setVertexSize NOTIFY vertexSizeChanged)
    Q_PROPERTY(uint count READ count WRITE setCount NOTIFY countChanged)
    Q_PROPERTY(uint byteStride READ byteStride WRITE setByteStride NOTIFY byteStrideChanged)
    Q_PROPERTY(uint byteOffset READ byteOffset WRITE setByteOffset NOTIFY byteOffsetChanged)
    Q_PROPERTY(uint divisor READ divisor WRITE setDivisor NOTIFY divisorChanged)
    Q_PROPERTY(AttributeType attributeType READ attributeType WRITE setAttributeType NOTIFY attributeTypeChanged)
    Q_PROPERTY(QString defaultPositionAttributeName READ defaultPositionAttributeName CONSTANT)
    Q_PROPERTY(QString defaultNormalAttributeName READ defaultNormalAttributeName CONSTANT)
    Q_PROPERTY(QString defaultColorAttributeName READ defaultColorAttributeName CONSTANT)
    Q_PROPERTY(QString defaultTextureCoordinateAttributeName READ defaultTextureCoordinateAttributeName CONSTANT)
    Q_PROPERTY(QString defaultTextureCoordinate1AttributeName READ defaultTextureCoordinate1AttributeName CONSTANT REVISION 11)
    Q_PROPERTY(QString defaultTextureCoordinate2AttributeName READ defaultTextureCoordinate2AttributeName CONSTANT REVISION 11)
    Q_PROPERTY(QString defaultTangentAttributeName READ defaultTangentAttributeName CONSTANT)
    Q_PROPERTY(QString defaultJointIndicesAttributeName READ defaultJointIndicesAttributeName CONSTANT REVISION 10)
    Q_PROPERTY(QString defaultJointWeightsAttributeName READ defaultJointWeightsAttributeName CONSTANT REVISION 10)

public:
    enum AttributeType {
        VertexAttribute,
        IndexAttribute,
        DrawIndirectAttribute
    };

    Q_ENUM(AttributeType) // LCOV_EXCL_LINE

    enum VertexBaseType {
        Byte = 0,
        UnsignedByte,
        Short,
        UnsignedShort,
        Int,
        UnsignedInt,
        HalfFloat,
        Float,
        Double
    };
    Q_ENUM(VertexBaseType) // LCOV_EXCL_LINE

    explicit QAttribute(QNode *parent = nullptr);
    explicit QAttribute(QBuffer *buf, VertexBaseType vertexBaseType, uint vertexSize, uint count, uint offset = 0, uint stride = 0, QNode *parent = nullptr);
    explicit QAttribute(QBuffer *buf, const QString &name, VertexBaseType vertexBaseType, uint vertexSize, uint count, uint offset = 0, uint stride = 0, QNode *parent = nullptr);
    ~QAttribute();

    QBuffer *buffer() const;
    QString name() const;
    VertexBaseType vertexBaseType() const;
    uint vertexSize() const;
    uint count() const;
    uint byteStride() const;
    uint byteOffset() const;
    uint divisor() const;
    AttributeType attributeType() const;

    Q_INVOKABLE static QString defaultPositionAttributeName();
    Q_INVOKABLE static QString defaultNormalAttributeName();
    Q_INVOKABLE static QString defaultColorAttributeName();
    Q_INVOKABLE static QString defaultTextureCoordinateAttributeName();
    Q_INVOKABLE static QString defaultTangentAttributeName();
    static QString defaultJointIndicesAttributeName();
    static QString defaultJointWeightsAttributeName();
    static QString defaultTextureCoordinate1AttributeName();
    static QString defaultTextureCoordinate2AttributeName();

public Q_SLOTS:
    void setBuffer(QBuffer *buffer);
    void setName(const QString &name);
    void setVertexBaseType(VertexBaseType type);
    void setVertexSize(uint size);
    QT_DEPRECATED void setDataType(VertexBaseType type);
    QT_DEPRECATED void setDataSize(uint size);
    void setCount(uint count);
    void setByteStride(uint byteStride);
    void setByteOffset(uint byteOffset);
    void setDivisor(uint divisor);
    void setAttributeType(AttributeType attributeType);

Q_SIGNALS:
    void bufferChanged(QBuffer *buffer);
    void nameChanged(const QString &name);
    void vertexBaseTypeChanged(VertexBaseType vertexBaseType);
    void vertexSizeChanged(uint vertexSize);
    void dataTypeChanged(VertexBaseType vertexBaseType);
    void dataSizeChanged(uint vertexSize);
    void countChanged(uint count);
    void byteStrideChanged(uint byteStride);
    void byteOffsetChanged(uint byteOffset);
    void divisorChanged(uint divisor);
    void attributeTypeChanged(AttributeType attributeType);

private:
    Q_DECLARE_PRIVATE(QAttribute)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QATTRIBUTE_H

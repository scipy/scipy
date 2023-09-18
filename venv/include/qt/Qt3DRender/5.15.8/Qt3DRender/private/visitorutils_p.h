/****************************************************************************
**
** Copyright (C) 2015 Paul Lemire paul.lemire350@gmail.com
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

#ifndef QT3DRENDER_RENDER_VISITORUTILS_P_H
#define QT3DRENDER_RENDER_VISITORUTILS_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtGlobal>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

namespace Render {

namespace Visitor {

template <QAttribute::VertexBaseType> struct EnumToType;
template <> struct EnumToType<QAttribute::Byte> { typedef const char type; };
template <> struct EnumToType<QAttribute::UnsignedByte> { typedef const uchar type; };
template <> struct EnumToType<QAttribute::Short> { typedef const short type; };
template <> struct EnumToType<QAttribute::UnsignedShort> { typedef const ushort type; };
template <> struct EnumToType<QAttribute::Int> { typedef const int type; };
template <> struct EnumToType<QAttribute::UnsignedInt> { typedef const uint type; };
template <> struct EnumToType<QAttribute::Float> { typedef const float type; };
template <> struct EnumToType<QAttribute::Double> { typedef const double type; };

template<QAttribute::VertexBaseType v>
inline typename EnumToType<v>::type *castToType(const QByteArray &u, uint byteOffset)
{
    return reinterpret_cast< typename EnumToType<v>::type *>(u.constData() + byteOffset);
}

template<typename Func>
void processBuffer(const BufferInfo &info, Func &f)
{
    switch (info.type) {
    case QAttribute::Byte: f(info, castToType<QAttribute::Byte>(info.data, info.byteOffset));
        return;
    case QAttribute::UnsignedByte: f(info, castToType<QAttribute::UnsignedByte>(info.data, info.byteOffset));
        return;
    case QAttribute::Short: f(info, castToType<QAttribute::Short>(info.data, info.byteOffset));
        return;
    case QAttribute::UnsignedShort: f(info, castToType<QAttribute::UnsignedShort>(info.data, info.byteOffset));
        return;
    case QAttribute::Int: f(info, castToType<QAttribute::Int>(info.data, info.byteOffset));
        return;
    case QAttribute::UnsignedInt: f(info, castToType<QAttribute::UnsignedInt>(info.data, info.byteOffset));
        return;
    case QAttribute::Float: f(info, castToType<QAttribute::Float>(info.data, info.byteOffset));
        return;
    case QAttribute::Double: f(info, castToType<QAttribute::Double>(info.data, info.byteOffset));
        return;
    default:
        return;
    }
}

template<typename VertexExecutor, typename IndexExecutor, typename Visitor>
void visitPrimitives(NodeManagers *manager, const GeometryRenderer *renderer, Visitor* visitor)
{
    Geometry *geom = manager->lookupResource<Geometry, GeometryManager>(renderer->geometryId());
    Attribute *positionAttribute = nullptr;
    Attribute *indexAttribute = nullptr;
    Buffer *positionBuffer = nullptr;
    Buffer *indexBuffer = nullptr;

    auto updateStride = [](BufferInfo &info, int stride) {
        if (stride) {
            info.byteStride = stride;
            return;
        }
        switch (info.type) {
        case QAttribute::Byte: info.byteStride = sizeof(qint8) * info.dataSize; return;
        case QAttribute::UnsignedByte: info.byteStride = sizeof(quint8) * info.dataSize; return;
        case QAttribute::Short: info.byteStride = sizeof(qint16) * info.dataSize; return;
        case QAttribute::UnsignedShort: info.byteStride = sizeof(quint16) * info.dataSize; return;
        case QAttribute::Int: info.byteStride = sizeof(qint32) * info.dataSize; return;
        case QAttribute::UnsignedInt: info.byteStride = sizeof(quint32) * info.dataSize; return;
        case QAttribute::Float: info.byteStride = sizeof(float) * info.dataSize; return;
        case QAttribute::Double: info.byteStride = sizeof(double) * info.dataSize; return;
        default: return;
        }
    };

    if (geom) {
        Qt3DRender::Render::Attribute *attribute = nullptr;
        const auto attrIds = geom->attributes();
        for (const Qt3DCore::QNodeId attrId : attrIds) {
            attribute = manager->lookupResource<Attribute, AttributeManager>(attrId);
            if (attribute){
                if (!positionAttribute && attribute->name() == QAttribute::defaultPositionAttributeName())
                    positionAttribute = attribute;
                else if (attribute->attributeType() == QAttribute::IndexAttribute)
                    indexAttribute = attribute;
            }
        }

        if (positionAttribute)
            positionBuffer = manager->lookupResource<Buffer, BufferManager>(positionAttribute->bufferId());
        if (indexAttribute)
            indexBuffer = manager->lookupResource<Buffer, BufferManager>(indexAttribute->bufferId());

        if (positionBuffer) {

            BufferInfo vertexBufferInfo;
            vertexBufferInfo.data = positionBuffer->data();
            vertexBufferInfo.type = positionAttribute->vertexBaseType();
            vertexBufferInfo.byteOffset = positionAttribute->byteOffset();
            vertexBufferInfo.dataSize = positionAttribute->vertexSize();
            vertexBufferInfo.count = positionAttribute->count();
            updateStride(vertexBufferInfo, positionAttribute->byteStride());

            if (indexBuffer) { // Indexed

                BufferInfo indexBufferInfo;
                indexBufferInfo.data = indexBuffer->data();
                indexBufferInfo.type = indexAttribute->vertexBaseType();
                indexBufferInfo.byteOffset = indexAttribute->byteOffset();
                indexBufferInfo.count = indexAttribute->count();
                indexBufferInfo.restartEnabled = renderer->primitiveRestartEnabled();
                indexBufferInfo.restartIndexValue = renderer->restartIndexValue();
                updateStride(indexBufferInfo, indexAttribute->byteStride());

                IndexExecutor executor;
                executor.m_vertexBufferInfo = vertexBufferInfo;
                executor.m_primitiveType = renderer->primitiveType();
                executor.m_visitor = visitor;

                return processBuffer(indexBufferInfo, executor);

            } else { // Non Indexed

                // Check into which type the buffer needs to be casted
                VertexExecutor executor;
                executor.m_primitiveType = renderer->primitiveType();
                executor.m_visitor = visitor;

                return processBuffer(vertexBufferInfo, executor);
            }
        }
    }
}

} // namespace Visitor

} // namespace Render

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_VISITORUTILS_P_H

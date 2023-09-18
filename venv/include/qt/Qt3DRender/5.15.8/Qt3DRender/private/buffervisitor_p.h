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

#ifndef QT3DRENDER_RENDER_BUFFERVISITOR_P_H
#define QT3DRENDER_RENDER_BUFFERVISITOR_P_H

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

#include <Qt3DCore/qnodeid.h>
#include <Qt3DRender/QAttribute>
#include <Qt3DRender/private/trianglesvisitor_p.h>
#include <Qt3DRender/private/attribute_p.h>
#include <Qt3DRender/private/buffer_p.h>
#include <Qt3DRender/private/bufferutils_p.h>
#include <Qt3DRender/private/geometryrenderer_p.h>
#include <Qt3DRender/private/geometry_p.h>
#include <Qt3DRender/private/managers_p.h>
#include <Qt3DRender/private/nodemanagers_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {
class QEntity;
}

namespace Qt3DRender {
namespace Render {


template <typename ValueType, QAttribute::VertexBaseType VertexBaseType, uint dataSize>
class Q_AUTOTEST_EXPORT BufferVisitor
{
public:
    explicit BufferVisitor(NodeManagers *manager)
        : m_manager(manager)
    {
    }
    virtual ~BufferVisitor() { }

    virtual void visit(uint ndx, ValueType x) {
        Q_UNUSED(ndx); Q_UNUSED(x);
    }
    virtual void visit(uint ndx, ValueType x, ValueType y) {
        Q_UNUSED(ndx); Q_UNUSED(x); Q_UNUSED(y);
    }
    virtual void visit(uint ndx, ValueType x, ValueType y, ValueType z) {
        Q_UNUSED(ndx); Q_UNUSED(x); Q_UNUSED(y); Q_UNUSED(z);
    }
    virtual void visit(uint ndx, ValueType x, ValueType y, ValueType z, ValueType w) {
        Q_UNUSED(ndx); Q_UNUSED(x); Q_UNUSED(y); Q_UNUSED(z); Q_UNUSED(w);
    }

    template<typename VertexBufferType, typename IndexBufferType>
    void traverseCoordinateIndexed(VertexBufferType *vertexBuffer,
                                   IndexBufferType *indexBuffer,
                                   int vertexByteStride,
                                   int drawVertexCount,
                                   bool primitiveRestartEnabled,
                                   int primitiveRestartIndex)
    {
        switch (dataSize) {
        case 1: traverseCoordinates1Indexed(vertexBuffer, vertexByteStride, indexBuffer, drawVertexCount,
                                            primitiveRestartEnabled, primitiveRestartIndex);
            break;
        case 2: traverseCoordinates2Indexed(vertexBuffer, vertexByteStride, indexBuffer, drawVertexCount,
                                            primitiveRestartEnabled, primitiveRestartIndex);
            break;
        case 3: traverseCoordinates3Indexed(vertexBuffer, vertexByteStride, indexBuffer, drawVertexCount,
                                            primitiveRestartEnabled, primitiveRestartIndex);
            break;
        case 4: traverseCoordinates4Indexed(vertexBuffer, vertexByteStride, indexBuffer, drawVertexCount,
                                            primitiveRestartEnabled, primitiveRestartIndex);
            break;
        default: Q_UNREACHABLE();
        }
    }

    bool apply(Qt3DRender::Render::Attribute *attribute,
               Qt3DRender::Render::Attribute *indexAttribute,
               int drawVertexCount,
               bool primitiveRestartEnabled,
               int primitiveRestartIndex)
    {
        if (attribute->vertexBaseType() != VertexBaseType)
            return false;
        if (attribute->vertexSize() < dataSize)
            return false;

        auto data = m_manager->lookupResource<Buffer, BufferManager>(attribute->bufferId())->data();
        auto vertexBuffer = BufferTypeInfo::castToType<VertexBaseType>(data, attribute->byteOffset());

        if (indexAttribute) {
            auto indexData = m_manager->lookupResource<Buffer, BufferManager>(indexAttribute->bufferId())->data();
            switch (indexAttribute->vertexBaseType()) {
            case QAttribute::UnsignedShort: {
                auto indexBuffer = BufferTypeInfo::castToType<QAttribute::UnsignedShort>(indexData, indexAttribute->byteOffset());
                traverseCoordinateIndexed(vertexBuffer, indexBuffer, attribute->byteStride(), drawVertexCount,
                                          primitiveRestartEnabled, primitiveRestartIndex);
                break;
            }
            case QAttribute::UnsignedInt: {
                auto indexBuffer = BufferTypeInfo::castToType<QAttribute::UnsignedInt>(indexData, indexAttribute->byteOffset());
                traverseCoordinateIndexed(vertexBuffer, indexBuffer, attribute->byteStride(), drawVertexCount,
                                          primitiveRestartEnabled, primitiveRestartIndex);
                break;
            }
            case QAttribute::UnsignedByte: {
                auto indexBuffer = BufferTypeInfo::castToType<QAttribute::UnsignedByte>(indexData, indexAttribute->byteOffset());
                traverseCoordinateIndexed(vertexBuffer, indexBuffer, attribute->byteStride(), drawVertexCount,
                                          primitiveRestartEnabled, primitiveRestartIndex);
                break;
            }
            default: Q_UNREACHABLE();
            }
        } else {
            switch (dataSize) {
            case 1: traverseCoordinates1(vertexBuffer, attribute->byteStride(), drawVertexCount); break;
            case 2: traverseCoordinates2(vertexBuffer, attribute->byteStride(), drawVertexCount); break;
            case 3: traverseCoordinates3(vertexBuffer, attribute->byteStride(), drawVertexCount); break;
            case 4: traverseCoordinates4(vertexBuffer, attribute->byteStride(), drawVertexCount); break;
            default: Q_UNREACHABLE();
            }
        }

        return true;
    }

protected:

    template <typename Coordinate>
    void traverseCoordinates1(Coordinate *coordinates,
                              const uint byteStride,
                              const uint count)
    {
        const uint stride = byteStride / sizeof(Coordinate);
        for (uint ndx = 0; ndx < count; ++ndx) {
            visit(ndx, coordinates[0]);
            coordinates += stride;
        }
    }

    template <typename Coordinate, typename IndexElem>
    void traverseCoordinates1Indexed(Coordinate *coordinates,
                                     const uint byteStride,
                                     IndexElem *indices,
                                     const uint count,
                                     bool primitiveRestartEnabled,
                                     int primitiveRestartIndex)
    {
        const uint stride = byteStride / sizeof(Coordinate);
        for (uint i = 0; i < count; ++i) {
            if (!primitiveRestartEnabled || (int) indices[i] != primitiveRestartIndex) {
                const uint n = stride * indices[i];
                visit(i, coordinates[n]);
            }
        }
    }

    template <typename Coordinate>
    void traverseCoordinates2(Coordinate *coordinates,
                              const uint byteStride,
                              const uint count)
    {
        const uint stride = byteStride ? byteStride / sizeof(Coordinate) : 2;
        for (uint ndx = 0; ndx < count; ++ndx) {
            visit(ndx, coordinates[0], coordinates[1]);
            coordinates += stride;
        }
    }


    template <typename Coordinate, typename IndexElem>
    void traverseCoordinates2Indexed(Coordinate *coordinates,
                                     const uint byteStride,
                                     IndexElem *indices,
                                     const uint count,
                                     bool primitiveRestartEnabled,
                                     int primitiveRestartIndex)
    {
        const uint stride = byteStride ? byteStride / sizeof(Coordinate) : 2;
        for (uint i = 0; i < count; ++i) {
            if (!primitiveRestartEnabled || (int) indices[i] != primitiveRestartIndex) {
                const uint n = stride * indices[i];
                visit(i, coordinates[n], coordinates[n + 1]);
            }
        }
    }

    template <typename Coordinate>
    void traverseCoordinates3(Coordinate *coordinates,
                              const uint byteStride,
                              const uint count)
    {
        const uint stride = byteStride ? byteStride / sizeof(Coordinate) : 3;
        for (uint ndx = 0; ndx < count; ++ndx) {
            visit(ndx, coordinates[0], coordinates[1], coordinates[2]);
            coordinates += stride;
        }
    }

    template <typename Coordinate, typename IndexElem>
    void traverseCoordinates3Indexed(Coordinate *coordinates,
                                     const uint byteStride,
                                     IndexElem *indices,
                                     const uint count,
                                     bool primitiveRestartEnabled,
                                     int primitiveRestartIndex)
    {
        const uint stride = byteStride ? byteStride / sizeof(Coordinate) : 3;
        for (uint i = 0; i < count; ++i) {
            if (!primitiveRestartEnabled || (int) indices[i] != primitiveRestartIndex) {
                const uint n = stride * indices[i];
                visit(i, coordinates[n], coordinates[n + 1], coordinates[n + 2]);
            }
        }
    }

    template <typename Coordinate>
    void traverseCoordinates4(Coordinate *coordinates,
                              const uint byteStride,
                              const uint count)
    {
        const uint stride = byteStride ? byteStride / sizeof(Coordinate) : 4;
        for (uint ndx = 0; ndx < count; ++ndx) {
            visit(ndx, coordinates[0], coordinates[1], coordinates[2], coordinates[3]);
            coordinates += stride;
        }
    }

    template <typename Coordinate, typename IndexElem>
    void traverseCoordinates4Indexed(Coordinate *coordinates,
                                     const uint byteStride,
                                     IndexElem *indices,
                                     const uint count,
                                     bool primitiveRestartEnabled,
                                     int primitiveRestartIndex)
    {
        const uint stride = byteStride ? byteStride / sizeof(Coordinate) : 4;
        for (uint i = 0; i < count; ++i) {
            if (!primitiveRestartEnabled || (int) indices[i] != primitiveRestartIndex) {
                const uint n = stride * indices[i];
                visit(i, coordinates[n], coordinates[n + 1], coordinates[n + 2], coordinates[n + 3]);
            }
        }
    }

    NodeManagers *m_manager;
};

typedef BufferVisitor<float, QAttribute::Float, 3> Buffer3fVisitor;

} // namespace Render

} // namespace Qt3DRender

QT_END_NAMESPACE


#endif // QT3DRENDER_RENDER_BUFFERVISITOR_P_H

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

#ifndef QT3DRENDER_RENDER_BUFFERUTILS_P_H
#define QT3DRENDER_RENDER_BUFFERUTILS_P_H

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

#include <Qt3DRender/QAttribute>
#include <QByteArray>

QT_BEGIN_NAMESPACE


namespace Qt3DRender {

namespace Render {

class GeometryRenderer;
class NodeManagers;
class Attribute;
class Buffer;

struct BufferInfo
{
    BufferInfo()
        : type(QAttribute::VertexBaseType::Float)
        , dataSize(0)
        , count(0)
        , byteStride(0)
        , byteOffset(0)
        , restartEnabled(false)
        , restartIndexValue(-1)
    {}

    QByteArray data;
    QAttribute::VertexBaseType type;
    uint dataSize;
    uint count;
    uint byteStride;
    uint byteOffset;
    bool restartEnabled;
    int restartIndexValue;
};


namespace BufferTypeInfo {

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
    typename EnumToType<v>::type *castToType(const QByteArray &u, uint byteOffset)
    {
        return reinterpret_cast< typename EnumToType<v>::type *>(u.constData() + byteOffset);
    }

} // namespace BufferTypeInfo

} // namespace Render

} // namespace Qt3DRender

QT_END_NAMESPACE


#endif // QT3DRENDER_RENDER_BUFFERUTILS_P_H

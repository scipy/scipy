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

#ifndef QT3DRENDER_QBUFFER_H
#define QT3DRENDER_QBUFFER_H

#include <Qt3DCore/qnode.h>
#include <Qt3DRender/qt3drender_global.h>
#include <QtCore/QSharedPointer>


QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QBufferPrivate;
class QBufferDataGenerator;
typedef QSharedPointer<QBufferDataGenerator> QBufferDataGeneratorPtr;

class Q_3DRENDERSHARED_EXPORT QBuffer : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(BufferType type READ type WRITE setType NOTIFY typeChanged)
    Q_PROPERTY(UsageType usage READ usage WRITE setUsage NOTIFY usageChanged)
    Q_PROPERTY(bool syncData READ isSyncData WRITE setSyncData NOTIFY syncDataChanged)
    Q_PROPERTY(AccessType accessType READ accessType WRITE setAccessType NOTIFY accessTypeChanged REVISION 9)

public:
    enum BufferType
    {
        VertexBuffer        = 0x8892, // GL_ARRAY_BUFFER
        IndexBuffer         = 0x8893, // GL_ELEMENT_ARRAY_BUFFER
        PixelPackBuffer     = 0x88EB, // GL_PIXEL_PACK_BUFFER
        PixelUnpackBuffer   = 0x88EC, // GL_PIXEL_UNPACK_BUFFER
        UniformBuffer       = 0x8A11, // GL_UNIFORM_BUFFER
        ShaderStorageBuffer = 0x90D2, // GL_SHADER_STORAGE_BUFFER
        DrawIndirectBuffer  = 0x8F3F  // GL_DRAW_INDIRECT_BUFFER
    };
    Q_ENUM(BufferType) // LCOV_EXCL_LINE

    enum UsageType
    {
        StreamDraw          = 0x88E0, // GL_STREAM_DRAW
        StreamRead          = 0x88E1, // GL_STREAM_READ
        StreamCopy          = 0x88E2, // GL_STREAM_COPY
        StaticDraw          = 0x88E4, // GL_STATIC_DRAW
        StaticRead          = 0x88E5, // GL_STATIC_READ
        StaticCopy          = 0x88E6, // GL_STATIC_COPY
        DynamicDraw         = 0x88E8, // GL_DYNAMIC_DRAW
        DynamicRead         = 0x88E9, // GL_DYNAMIC_READ
        DynamicCopy         = 0x88EA  // GL_DYNAMIC_COPY
    };
    Q_ENUM(UsageType) // LCOV_EXCL_LINE

    enum AccessType {
        Write = 0x1,
        Read = 0x2,
        ReadWrite = Write|Read
    };
    Q_ENUM(AccessType) // LCOV_EXCL_LINE

    explicit QBuffer(Qt3DCore::QNode *parent = nullptr);
    QT_DEPRECATED explicit QBuffer(BufferType ty, Qt3DCore::QNode *parent = nullptr);
    ~QBuffer();

    UsageType usage() const;
    QT_DEPRECATED BufferType type() const;
    bool isSyncData() const;
    AccessType accessType() const;

    void setData(const QByteArray &bytes);
    QByteArray data() const;

    Q3D_DECL_DEPRECATED void setDataGenerator(const QBufferDataGeneratorPtr &functor);
    Q3D_DECL_DEPRECATED QBufferDataGeneratorPtr dataGenerator() const;

    Q_INVOKABLE void updateData(int offset, const QByteArray &bytes);

public Q_SLOTS:
    QT_DEPRECATED void setType(BufferType type);
    void setUsage(UsageType usage);
    void setSyncData(bool syncData);
    void setAccessType(AccessType access);

Q_SIGNALS:
    void dataChanged(const QByteArray &bytes);
    void typeChanged(BufferType type);
    void usageChanged(UsageType usage);
    void syncDataChanged(bool syncData);
    void accessTypeChanged(AccessType access);
    void dataAvailable();

protected:
    // TODO Unused remove in Qt6
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

private:
    Q_DECLARE_PRIVATE(QBuffer)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QBUFFER_H

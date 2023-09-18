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

#ifndef QT3DRENDER_QTEXTURE_H
#define QT3DRENDER_QTEXTURE_H

#include <Qt3DRender/qtextureimagedata.h>
#include <Qt3DRender/qtexturewrapmode.h>
#include <Qt3DRender/qtextureimage.h>
#include <Qt3DRender/qabstracttexture.h>

QT_BEGIN_NAMESPACE

#if 0
// Create a forwarding header
#pragma qt_class(QTexture)
#pragma qt_sync_stop_processing
#endif

namespace Qt3DRender {

class Q_3DRENDERSHARED_EXPORT QTexture1D : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTexture1D(Qt3DCore::QNode *parent = nullptr);
    ~QTexture1D();
};

class Q_3DRENDERSHARED_EXPORT QTexture1DArray : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTexture1DArray(Qt3DCore::QNode *parent = nullptr);
    ~QTexture1DArray();
};

class Q_3DRENDERSHARED_EXPORT QTexture2D : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTexture2D(Qt3DCore::QNode *parent = nullptr);
    ~QTexture2D();
};

class Q_3DRENDERSHARED_EXPORT QTexture2DArray : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTexture2DArray(Qt3DCore::QNode *parent = nullptr);
    ~QTexture2DArray();
};

class Q_3DRENDERSHARED_EXPORT QTexture3D : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTexture3D(Qt3DCore::QNode *parent = nullptr);
    ~QTexture3D();
};

class Q_3DRENDERSHARED_EXPORT QTextureCubeMap : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTextureCubeMap(Qt3DCore::QNode *parent = nullptr);
    ~QTextureCubeMap();
};

class Q_3DRENDERSHARED_EXPORT QTextureCubeMapArray : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTextureCubeMapArray(Qt3DCore::QNode *parent = nullptr);
    ~QTextureCubeMapArray();
};

class Q_3DRENDERSHARED_EXPORT QTexture2DMultisample : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTexture2DMultisample(Qt3DCore::QNode *parent = nullptr);
    ~QTexture2DMultisample();
};

class Q_3DRENDERSHARED_EXPORT QTexture2DMultisampleArray : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTexture2DMultisampleArray(Qt3DCore::QNode *parent = nullptr);
    ~QTexture2DMultisampleArray();
};

class Q_3DRENDERSHARED_EXPORT QTextureRectangle : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTextureRectangle(Qt3DCore::QNode *parent = nullptr);
    ~QTextureRectangle();
};

class Q_3DRENDERSHARED_EXPORT QTextureBuffer : public QAbstractTexture
{
    Q_OBJECT
public:
    explicit QTextureBuffer(Qt3DCore::QNode *parent = nullptr);
    ~QTextureBuffer();
};

class QTextureLoaderPrivate;

class Q_3DRENDERSHARED_EXPORT QTextureLoader : public QAbstractTexture
{
    Q_OBJECT
    Q_PROPERTY(QUrl source READ source WRITE setSource NOTIFY sourceChanged)
    Q_PROPERTY(bool mirrored READ isMirrored WRITE setMirrored NOTIFY mirroredChanged)
public:
    explicit QTextureLoader(Qt3DCore::QNode *parent = nullptr);
    ~QTextureLoader();

    QUrl source() const;
    bool isMirrored() const;

public Q_SLOTS:
    void setSource(const QUrl &source);
    void setMirrored(bool mirrored);

Q_SIGNALS:
    void sourceChanged(const QUrl &source);
    void mirroredChanged(bool mirrored);

private:
    Q_DECLARE_PRIVATE(QTextureLoader)
};

class Q_3DRENDERSHARED_EXPORT QSharedGLTexture : public QAbstractTexture
{
    Q_OBJECT
    Q_PROPERTY(int textureId READ textureId WRITE setTextureId NOTIFY textureIdChanged)
public:
    explicit QSharedGLTexture(Qt3DCore::QNode *parent = nullptr);
    ~QSharedGLTexture();

    int textureId() const;

public Q_SLOTS:
    void setTextureId(int id);

Q_SIGNALS:
    void textureIdChanged(int textureId);
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QTEXTURE_H

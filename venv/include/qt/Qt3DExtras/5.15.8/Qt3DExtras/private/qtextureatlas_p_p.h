/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DEXTRAS_QTEXTUREATLAS_P_P_H
#define QT3DEXTRAS_QTEXTUREATLAS_P_P_H

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

#include <QtCore/qscopedpointer.h>
#include <Qt3DRender/private/qabstracttexture_p.h>
#include <Qt3DRender/qtexturegenerator.h>
#include <Qt3DExtras/private/areaallocator_p.h>
#include <Qt3DExtras/private/qtextureatlas_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DExtras {

// Used to store texture info within atlas
struct AtlasTexture
{
    QRect position;
    int padding = 0;
};

// data shared between QTextureAtlasPrivate and the QTextureGenerators
// we use this extra indirection so we can lazily copy the sub-images
// into the actual texture image in the backend texture loader thread.
class QTextureAtlasData
{
public:
    QTextureAtlasData(int w, int h, QImage::Format fmt);
    ~QTextureAtlasData();

    int width() const { return m_image.width(); }
    int height() const { return m_image.height(); }

    void addImage(const AtlasTexture &texture, const QImage &image);
    QByteArray createUpdatedImageData();

private:
    struct Update {
        AtlasTexture textureInfo;
        QImage image;
    };

    QMutex m_mutex;
    QImage m_image;
    QVector<Update> m_updates;
};

typedef QSharedPointer<QTextureAtlasData> QTextureAtlasDataPtr;

class QTextureAtlasPrivate : public Qt3DRender::QAbstractTexturePrivate
{
public:
    QTextureAtlasPrivate();
    ~QTextureAtlasPrivate();

    Q_DECLARE_PUBLIC(QTextureAtlas)

    QTextureAtlas::TextureId m_currId = 1;  // IDs for new sub-textures
    int m_currGen = 0;

    QTextureAtlasDataPtr m_data;
    QScopedPointer<AreaAllocator> m_allocator;
    QOpenGLTexture::PixelFormat m_pixelFormat;
    QHash<QTextureAtlas::TextureId, AtlasTexture> m_textures;
};

class QTextureAtlasGenerator : public Qt3DRender::QTextureGenerator
{
public:
    QTextureAtlasGenerator(const QTextureAtlasPrivate *texAtlas);
    ~QTextureAtlasGenerator();
    Qt3DRender::QTextureDataPtr operator()() override;
    bool operator==(const QTextureGenerator &other) const override;

    QT3D_FUNCTOR(QTextureAtlasGenerator)

private:
    QTextureAtlasDataPtr m_data;
    Qt3DRender::QAbstractTexture::TextureFormat m_format;
    QOpenGLTexture::PixelFormat m_pixelFormat;
    int m_generation;
    Qt3DCore::QNodeId m_atlasId;
};
typedef QSharedPointer<QTextureAtlasGenerator> QTextureAtlasGeneratorPtr;

} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QTEXTUREATLAS_P_P_H

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

#ifndef QT3DRENDER_TEXTUREDATA_P_H
#define QT3DRENDER_TEXTUREDATA_P_H

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

#include "qtextureimagedata.h"
#include <Qt3DRender/private/qt3drender_global_p.h>
#include <functional>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class Q_3DRENDERSHARED_PRIVATE_EXPORT QTextureImageDataPrivate
{
public:
    QTextureImageDataPrivate();

    void setData(const QByteArray &data, int blockSize, bool isCompressed);
    void setData(const QByteArray &data,
                 std::function<QByteArray(QByteArray, int, int, int)> dataExtractor,
                 bool isCompressed);

    bool setCompressedFile(const QString &source);

    QByteArray data(int layer, int face, int mipmapLevel) const;

    int m_width;
    int m_height;
    int m_depth;
    int m_layers;
    int m_faces;
    int m_mipLevels;
    int m_blockSize;
    int m_alignment;

    QOpenGLTexture::Target m_target;
    QOpenGLTexture::TextureFormat m_format;
    QOpenGLTexture::PixelFormat m_pixelFormat;
    QOpenGLTexture::PixelType m_pixelType;

    bool m_isCompressed;
    // ### Qt 6
    // QTextureImageData was originally written with assumptions around the internal data format
    // matching dds layout. This is an ugly, but easy, way to add basic ktx support without any
    // public API changes. Consider https://codereview.qt-project.org/#/c/178474/ for Qt 6.
    bool m_isKtx;
    QByteArray m_data;
    std::function<QByteArray(QByteArray rawData, int layer, int face, int mipmapLevel)> m_dataExtractor;

    static QTextureImageDataPrivate *get(QTextureImageData *imageData);

private:
    int ddsLayerSize() const;
    int ddsFaceSize() const;
    int mipmapLevelSize(int level) const;
    QByteArray ktxData(int layer, int face, int mipmapLevel) const;
};

} // namespace Qt3DRender


QT_END_NAMESPACE

#endif // QT3DRENDER_TEXTUREDATA_P_H

/****************************************************************************
**
** Copyright (C) 2019 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QSHADERIMAGE_H
#define QT3DRENDER_QSHADERIMAGE_H

#include <Qt3DCore/qnode.h>
#include <Qt3DRender/qt3drender_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAbstractTexture;
class QShaderImagePrivate;

class Q_3DRENDERSHARED_EXPORT QShaderImage : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(Qt3DRender::QAbstractTexture *texture READ texture WRITE setTexture NOTIFY textureChanged)
    Q_PROPERTY(bool layered READ layered WRITE setLayered NOTIFY layeredChanged)
    Q_PROPERTY(int mipLevel READ mipLevel WRITE setMipLevel NOTIFY mipLevelChanged)
    Q_PROPERTY(int layer READ layer WRITE setLayer NOTIFY layerChanged)
    Q_PROPERTY(Access access READ access WRITE setAccess NOTIFY accessChanged)
    Q_PROPERTY(ImageFormat format READ format WRITE setFormat NOTIFY formatChanged)

public:
    enum Access {
        ReadOnly = 0,
        WriteOnly,
        ReadWrite
    };
    Q_ENUM(Access)

    enum ImageFormat {
        NoFormat               = 0,         // GL_NONE
        Automatic              = 1,         // The Qt3D engine automatically determines the best format

        // Unsigned normalized formats
        R8_UNorm               = 0x8229,    // GL_R8
        RG8_UNorm              = 0x822B,    // GL_RG8
        RGBA8_UNorm            = 0x8058,    // GL_RGBA8

        R16_UNorm              = 0x822A,    // GL_R16
        RG16_UNorm             = 0x822C,    // GL_RG16
        RGBA16_UNorm           = 0x805B,    // GL_RGBA16

        // Signed normalized formats
        R8_SNorm               = 0x8F94,    // GL_R8_SNORM
        RG8_SNorm              = 0x8F95,    // GL_RG8_SNORM
        RGBA8_SNorm            = 0x8F97,    // GL_RGBA8_SNORM

        R16_SNorm              = 0x8F98,    // GL_R16_SNORM
        RG16_SNorm             = 0x8F99,    // GL_RG16_SNORM
        RGBA16_SNorm           = 0x8F9B,    // GL_RGBA16_SNORM

        // Unsigned integer formats
        R8U                    = 0x8232,    // GL_R8UI
        RG8U                   = 0x8238,    // GL_RG8UI
        RGBA8U                 = 0x8D7C,    // GL_RGBA8UI

        R16U                   = 0x8234,    // GL_R16UI
        RG16U                  = 0x823A,    // GL_RG16UI
        RGBA16U                = 0x8D76,    // GL_RGBA16UI

        R32U                   = 0x8236,    // GL_R32UI
        RG32U                  = 0x823C,    // GL_RG32UI
        RGBA32U                = 0x8D70,    // GL_RGBA32UI

        // Signed integer formats
        R8I                    = 0x8231,    // GL_R8I
        RG8I                   = 0x8237,    // GL_RG8I
        RGBA8I                 = 0x8D8E,    // GL_RGBA8I

        R16I                   = 0x8233,    // GL_R16I
        RG16I                  = 0x8239,    // GL_RG16I
        RGBA16I                = 0x8D88,    // GL_RGBA16I

        R32I                   = 0x8235,    // GL_R32I
        RG32I                  = 0x823B,    // GL_RG32I
        RGBA32I                = 0x8D82,    // GL_RGBA32I

        // Floating point formats
        R16F                   = 0x822D,    // GL_R16F
        RG16F                  = 0x822F,    // GL_RG16F
        RGBA16F                = 0x881A,    // GL_RGBA16F

        R32F                   = 0x822E,    // GL_R32F
        RG32F                  = 0x8230,    // GL_RG32F
        RGBA32F                = 0x8814,    // GL_RGBA32F

        // Packed formats
        RG11B10F               = 0x8C3A,    // GL_R11F_G11F_B10F
        RGB10A2                = 0x8059,    // GL_RGB10_A2
        RGB10A2U               = 0x906F,    // GL_RGB10_A2_UI
    };
    Q_ENUM(ImageFormat)

    explicit QShaderImage(Qt3DCore::QNode *parent = nullptr);
    ~QShaderImage();

    Qt3DRender::QAbstractTexture *texture() const;
    bool layered() const;
    int mipLevel() const;
    int layer() const;
    Access access() const;
    ImageFormat format() const;

public Q_SLOTS:
    void setTexture(Qt3DRender::QAbstractTexture *texture);
    void setLayered(bool layered);
    void setMipLevel(int mipLevel);
    void setLayer(int layer);
    void setAccess(Access access);
    void setFormat(ImageFormat format);

Q_SIGNALS:
    void textureChanged(Qt3DRender::QAbstractTexture *texture);
    void layeredChanged(bool layered);
    void mipLevelChanged(int mipLevel);
    void layerChanged(int layer);
    void accessChanged(Access access);
    void formatChanged(ImageFormat format);

private:
    Q_DECLARE_PRIVATE(QShaderImage)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QSHADERIMAGE_H

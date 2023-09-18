/****************************************************************************
**
** Copyright (C) 2015 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QABSTRACTTEXTURE_H
#define QT3DRENDER_QABSTRACTTEXTURE_H

#include <Qt3DRender/qtextureimagedata.h>
#include <Qt3DRender/qt3drender_global.h>
#include <Qt3DCore/qnode.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QAbstractTexturePrivate;
class QTextureWrapMode;
class QAbstractTextureImage;
class QTextureGenerator;
class QTextureDataUpdate;
typedef QSharedPointer<QTextureGenerator> QTextureGeneratorPtr;

class Q_3DRENDERSHARED_EXPORT QAbstractTexture : public Qt3DCore::QNode
{
    Q_OBJECT
    Q_PROPERTY(Target target READ target CONSTANT)
    Q_PROPERTY(TextureFormat format READ format WRITE setFormat NOTIFY formatChanged)
    Q_PROPERTY(bool generateMipMaps READ generateMipMaps WRITE setGenerateMipMaps NOTIFY generateMipMapsChanged)
    Q_PROPERTY(Qt3DRender::QTextureWrapMode *wrapMode READ wrapMode CONSTANT)
    Q_PROPERTY(Status status READ status NOTIFY statusChanged)
    Q_PROPERTY(int width READ width WRITE setWidth NOTIFY widthChanged)
    Q_PROPERTY(int height READ height WRITE setHeight NOTIFY heightChanged)
    Q_PROPERTY(int depth READ depth WRITE setDepth NOTIFY depthChanged)
    Q_PROPERTY(Filter magnificationFilter READ magnificationFilter WRITE setMagnificationFilter NOTIFY magnificationFilterChanged)
    Q_PROPERTY(Filter minificationFilter READ minificationFilter WRITE setMinificationFilter NOTIFY minificationFilterChanged)
    Q_PROPERTY(float maximumAnisotropy READ maximumAnisotropy WRITE setMaximumAnisotropy NOTIFY maximumAnisotropyChanged)
    Q_PROPERTY(ComparisonFunction comparisonFunction READ comparisonFunction WRITE setComparisonFunction NOTIFY comparisonFunctionChanged)
    Q_PROPERTY(ComparisonMode comparisonMode READ comparisonMode WRITE setComparisonMode NOTIFY comparisonModeChanged)
    Q_PROPERTY(int layers READ layers WRITE setLayers NOTIFY layersChanged)
    Q_PROPERTY(int samples READ samples WRITE setSamples NOTIFY samplesChanged)
    Q_PROPERTY(HandleType handleType READ handleType NOTIFY handleTypeChanged REVISION 13)
    Q_PROPERTY(QVariant handle READ handle NOTIFY handleChanged REVISION 13)

public:

    enum Status {
        None = 0,
        Loading,
        Ready,
        Error
    };
    Q_ENUM(Status) // LCOV_EXCL_LINE

    enum Target {
        TargetAutomatic            = 0,         // Target will be determined by the Qt3D engine
        Target1D                   = 0x0DE0,    // GL_TEXTURE_1D
        Target1DArray              = 0x8C18,    // GL_TEXTURE_1D_ARRAY
        Target2D                   = 0x0DE1,    // GL_TEXTURE_2D
        Target2DArray              = 0x8C1A,    // GL_TEXTURE_2D_ARRAY
        Target3D                   = 0x806F,    // GL_TEXTURE_3D
        TargetCubeMap              = 0x8513,    // GL_TEXTURE_CUBE_MAP
        TargetCubeMapArray         = 0x9009,    // GL_TEXTURE_CUBE_MAP_ARRAY
        Target2DMultisample        = 0x9100,    // GL_TEXTURE_2D_MULTISAMPLE
        Target2DMultisampleArray   = 0x9102,    // GL_TEXTURE_2D_MULTISAMPLE_ARRAY
        TargetRectangle            = 0x84F5,    // GL_TEXTURE_RECTANGLE
        TargetBuffer               = 0x8C2A     // GL_TEXTURE_BUFFER
    };
    Q_ENUM(Target) // LCOV_EXCL_LINE

    enum TextureFormat {
        NoFormat               = 0,         // GL_NONE
        Automatic              = 1,         // The Qt3D engine automatically determines the best format

        // Unsigned normalized formats
        R8_UNorm               = 0x8229,    // GL_R8
        RG8_UNorm              = 0x822B,    // GL_RG8
        RGB8_UNorm             = 0x8051,    // GL_RGB8
        RGBA8_UNorm            = 0x8058,    // GL_RGBA8

        R16_UNorm              = 0x822A,    // GL_R16
        RG16_UNorm             = 0x822C,    // GL_RG16
        RGB16_UNorm            = 0x8054,    // GL_RGB16
        RGBA16_UNorm           = 0x805B,    // GL_RGBA16

        // Signed normalized formats
        R8_SNorm               = 0x8F94,    // GL_R8_SNORM
        RG8_SNorm              = 0x8F95,    // GL_RG8_SNORM
        RGB8_SNorm             = 0x8F96,    // GL_RGB8_SNORM
        RGBA8_SNorm            = 0x8F97,    // GL_RGBA8_SNORM

        R16_SNorm              = 0x8F98,    // GL_R16_SNORM
        RG16_SNorm             = 0x8F99,    // GL_RG16_SNORM
        RGB16_SNorm            = 0x8F9A,    // GL_RGB16_SNORM
        RGBA16_SNorm           = 0x8F9B,    // GL_RGBA16_SNORM

        // Unsigned integer formats
        R8U                    = 0x8232,    // GL_R8UI
        RG8U                   = 0x8238,    // GL_RG8UI
        RGB8U                  = 0x8D7D,    // GL_RGB8UI
        RGBA8U                 = 0x8D7C,    // GL_RGBA8UI

        R16U                   = 0x8234,    // GL_R16UI
        RG16U                  = 0x823A,    // GL_RG16UI
        RGB16U                 = 0x8D77,    // GL_RGB16UI
        RGBA16U                = 0x8D76,    // GL_RGBA16UI

        R32U                   = 0x8236,    // GL_R32UI
        RG32U                  = 0x823C,    // GL_RG32UI
        RGB32U                 = 0x8D71,    // GL_RGB32UI
        RGBA32U                = 0x8D70,    // GL_RGBA32UI

        // Signed integer formats
        R8I                    = 0x8231,    // GL_R8I
        RG8I                   = 0x8237,    // GL_RG8I
        RGB8I                  = 0x8D8F,    // GL_RGB8I
        RGBA8I                 = 0x8D8E,    // GL_RGBA8I

        R16I                   = 0x8233,    // GL_R16I
        RG16I                  = 0x8239,    // GL_RG16I
        RGB16I                 = 0x8D89,    // GL_RGB16I
        RGBA16I                = 0x8D88,    // GL_RGBA16I

        R32I                   = 0x8235,    // GL_R32I
        RG32I                  = 0x823B,    // GL_RG32I
        RGB32I                 = 0x8D83,    // GL_RGB32I
        RGBA32I                = 0x8D82,    // GL_RGBA32I

        // Floating point formats
        R16F                   = 0x822D,    // GL_R16F
        RG16F                  = 0x822F,    // GL_RG16F
        RGB16F                 = 0x881B,    // GL_RGB16F
        RGBA16F                = 0x881A,    // GL_RGBA16F

        R32F                   = 0x822E,    // GL_R32F
        RG32F                  = 0x8230,    // GL_RG32F
        RGB32F                 = 0x8815,    // GL_RGB32F
        RGBA32F                = 0x8814,    // GL_RGBA32F

        // Packed formats
        RGB9E5                 = 0x8C3D,    // GL_RGB9_E5
        RG11B10F               = 0x8C3A,    // GL_R11F_G11F_B10F
        RG3B2                  = 0x2A10,    // GL_R3_G3_B2
        R5G6B5                 = 0x8D62,    // GL_RGB565
        RGB5A1                 = 0x8057,    // GL_RGB5_A1
        RGBA4                  = 0x8056,    // GL_RGBA4
        RGB10A2                = 0x8059,    // GL_RGB10_A2
        RGB10A2U               = 0x906F,    // GL_RGB10_A2UI

        // Depth formats
        D16                    = 0x81A5,    // GL_DEPTH_COMPONENT16
        D24                    = 0x81A6,    // GL_DEPTH_COMPONENT24
        D24S8                  = 0x88F0,    // GL_DEPTH24_STENCIL8
        D32                    = 0x81A7,    // GL_DEPTH_COMPONENT32
        D32F                   = 0x8CAC,    // GL_DEPTH_COMPONENT32F
        D32FS8X24              = 0x8CAD,    // GL_DEPTH32F_STENCIL8

        // Compressed formats
        RGB_DXT1               = 0x83F0,    // GL_COMPRESSED_RGB_S3TC_DXT1_EXT
        RGBA_DXT1              = 0x83F1,    // GL_COMPRESSED_RGBA_S3TC_DXT1_EXT
        RGBA_DXT3              = 0x83F2,    // GL_COMPRESSED_RGBA_S3TC_DXT3_EXT
        RGBA_DXT5              = 0x83F3,    // GL_COMPRESSED_RGBA_S3TC_DXT5_EXT
        R_ATI1N_UNorm          = 0x8DBB,    // GL_COMPRESSED_RED_RGTC1
        R_ATI1N_SNorm          = 0x8DBC,    // GL_COMPRESSED_SIGNED_RED_RGTC1
        RG_ATI2N_UNorm         = 0x8DBD,    // GL_COMPRESSED_RG_RGTC2
        RG_ATI2N_SNorm         = 0x8DBE,    // GL_COMPRESSED_SIGNED_RG_RGTC2
        RGB_BP_UNSIGNED_FLOAT  = 0x8E8F,    // GL_COMPRESSED_RGB_BPTC_UNSIGNED_FLOAT_ARB
        RGB_BP_SIGNED_FLOAT    = 0x8E8E,    // GL_COMPRESSED_RGB_BPTC_SIGNED_FLOAT_ARB
        RGB_BP_UNorm           = 0x8E8C,    // GL_COMPRESSED_RGBA_BPTC_UNORM_ARB
        R11_EAC_UNorm          = 0x9270,    // GL_COMPRESSED_R11_EAC
        R11_EAC_SNorm          = 0x9271,    // GL_COMPRESSED_SIGNED_R11_EAC
        RG11_EAC_UNorm         = 0x9272,    // GL_COMPRESSED_RG11_EAC
        RG11_EAC_SNorm         = 0x9273,    // GL_COMPRESSED_SIGNED_RG11_EAC
        RGB8_ETC2              = 0x9274,    // GL_COMPRESSED_RGB8_ETC2
        SRGB8_ETC2             = 0x9275,    // GL_COMPRESSED_SRGB8_ETC2
        RGB8_PunchThrough_Alpha1_ETC2 = 0x9276, // GL_COMPRESSED_RGB8_PUNCHTHROUGH_ALPHA1_ETC2
        SRGB8_PunchThrough_Alpha1_ETC2 = 0x9277, // GL_COMPRESSED_SRGB8_PUNCHTHROUGH_ALPHA1_ETC2
        RGBA8_ETC2_EAC         = 0x9278,    // GL_COMPRESSED_RGBA8_ETC2_EAC
        SRGB8_Alpha8_ETC2_EAC  = 0x9279,    // GL_COMPRESSED_SRGB8_ALPHA8_ETC2_EAC
        RGB8_ETC1              = 0x8D64,    // GL_ETC1_RGB8_OES

        // sRGB formats
        SRGB8                  = 0x8C41,    // GL_SRGB8
        SRGB8_Alpha8           = 0x8C43,    // GL_SRGB8_ALPHA8
        SRGB_DXT1              = 0x8C4C,    // GL_COMPRESSED_SRGB_S3TC_DXT1_EXT
        SRGB_Alpha_DXT1        = 0x8C4D,    // GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT1_EXT
        SRGB_Alpha_DXT3        = 0x8C4E,    // GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT3_EXT
        SRGB_Alpha_DXT5        = 0x8C4F,    // GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT
        SRGB_BP_UNorm          = 0x8E8D,    // GL_COMPRESSED_SRGB_ALPHA_BPTC_UNORM_ARB

        // ES 2 formats
        DepthFormat            = 0x1902,    // GL_DEPTH_COMPONENT
        AlphaFormat            = 0x1906,    // GL_ALPHA
        RGBFormat              = 0x1907,    // GL_RGB
        RGBAFormat             = 0x1908,    // GL_RGBA
        LuminanceFormat        = 0x1909,    // GL_LUMINANCE
        LuminanceAlphaFormat   = 0x190A
    };
    Q_ENUM(TextureFormat) // LCOV_EXCL_LINE

    enum Filter {
        Nearest                 = 0x2600,   // GL_NEAREST
        Linear                  = 0x2601,   // GL_LINEAR
        NearestMipMapNearest    = 0x2700,   // GL_NEAREST_MIPMAP_NEAREST
        NearestMipMapLinear     = 0x2702,   // GL_NEAREST_MIPMAP_LINEAR
        LinearMipMapNearest     = 0x2701,   // GL_LINEAR_MIPMAP_NEAREST
        LinearMipMapLinear      = 0x2703    // GL_LINEAR_MIPMAP_LINEAR
    };
    Q_ENUM(Filter) // LCOV_EXCL_LINE

    enum CubeMapFace {
        CubeMapPositiveX = 0x8515,  // GL_TEXTURE_CUBE_MAP_POSITIVE_X
        CubeMapNegativeX = 0x8516,  // GL_TEXTURE_CUBE_MAP_NEGATIVE_X
        CubeMapPositiveY = 0x8517,  // GL_TEXTURE_CUBE_MAP_POSITIVE_Y
        CubeMapNegativeY = 0x8518,  // GL_TEXTURE_CUBE_MAP_NEGATIVE_Y
        CubeMapPositiveZ = 0x8519,  // GL_TEXTURE_CUBE_MAP_POSITIVE_Z
        CubeMapNegativeZ = 0x851A,  // GL_TEXTURE_CUBE_MAP_NEGATIVE_Z
        AllFaces
    };
    Q_ENUM(CubeMapFace) // LCOV_EXCL_LINE

    enum ComparisonFunction {
        CompareLessEqual    = 0x0203,   // GL_LEQUAL
        CompareGreaterEqual = 0x0206,   // GL_GEQUAL
        CompareLess         = 0x0201,   // GL_LESS
        CompareGreater      = 0x0204,   // GL_GREATER
        CompareEqual        = 0x0202,   // GL_EQUAL
        CommpareNotEqual    = 0x0205,   // GL_NOTEQUAL
        CompareAlways       = 0x0207,   // GL_ALWAYS
        CompareNever        = 0x0200    // GL_NEVER
    };
    Q_ENUM(ComparisonFunction) // LCOV_EXCL_LINE

    enum ComparisonMode {
        CompareRefToTexture = 0x884E,   // GL_COMPARE_REF_TO_TEXTURE
        CompareNone         = 0x0000    // GL_NONE
    };
    Q_ENUM(ComparisonMode) // LCOV_EXCL_LINE

    enum HandleType {
        NoHandle,
        OpenGLTextureId
    };
    Q_ENUM(HandleType) // LCOV_EXCL_LINE

    ~QAbstractTexture();

    Target target() const;

    TextureFormat format() const;
    bool generateMipMaps() const;

    Status status() const;

    void addTextureImage(QAbstractTextureImage *textureImage);
    void removeTextureImage(QAbstractTextureImage *textureImage);
    QVector<QAbstractTextureImage *> textureImages() const;

    // sampler data - in the future proxy to a Sampler helper
    void setWrapMode(const QTextureWrapMode &wrapMode);
    QTextureWrapMode *wrapMode();

    void setSize(int width, int height=1, int depth=1);

    Filter minificationFilter() const;
    Filter magnificationFilter() const;
    float maximumAnisotropy() const;
    ComparisonFunction comparisonFunction() const;
    ComparisonMode comparisonMode() const;
    int width() const;
    int height() const;
    int depth() const;
    int layers() const;
    int samples() const;
    Q3D_DECL_DEPRECATED QTextureGeneratorPtr dataGenerator() const;
    HandleType handleType() const;
    QVariant handle() const;

    Q_INVOKABLE void updateData(const QTextureDataUpdate &update);


public Q_SLOTS:
    void setFormat(TextureFormat format);
    void setGenerateMipMaps(bool gen);
    void setWidth(int width);
    void setHeight(int height);
    void setDepth(int depth);
    void setMinificationFilter(Filter f);
    void setMagnificationFilter(Filter f);
    void setMaximumAnisotropy(float anisotropy);
    void setComparisonFunction(ComparisonFunction function);
    void setComparisonMode(ComparisonMode mode);
    void setLayers(int layers);
    void setSamples(int samples);

Q_SIGNALS:
    void formatChanged(TextureFormat format);
    void statusChanged(Status status);
    void generateMipMapsChanged(bool generateMipMaps);
    void widthChanged(int width);
    void heightChanged(int height);
    void depthChanged(int depth);
    void magnificationFilterChanged(Filter magnificationFilter);
    void minificationFilterChanged(Filter minificationFilter);
    void maximumAnisotropyChanged(float maximumAnisotropy);
    void comparisonFunctionChanged(ComparisonFunction comparisonFunction);
    void comparisonModeChanged(ComparisonMode comparisonMode);
    void layersChanged(int layers);
    void samplesChanged(int samples);
    Q_REVISION(13) void handleTypeChanged(HandleType handleType);
    Q_REVISION(13) void handleChanged(QVariant handle);

protected:
    explicit QAbstractTexture(Qt3DCore::QNode *parent = nullptr);
    explicit QAbstractTexture(Target target, Qt3DCore::QNode *parent = nullptr);
    explicit QAbstractTexture(QAbstractTexturePrivate &dd, Qt3DCore::QNode *parent = nullptr);
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;

    // TO DO Qt6, should be on private class
    void setStatus(Status status);
    void setHandle(const QVariant &handle);
    void setHandleType(HandleType type);

private:
    Q_DECLARE_PRIVATE(QAbstractTexture)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DRender::QAbstractTexture *) // LCOV_EXCL_LINE

#endif // QT3DRENDER_QABSTRACTTEXTURE_H

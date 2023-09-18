/****************************************************************************
**
** Copyright (C) 2008-2012 NVIDIA Corporation.
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QSSG_RENDER_SHADER_KEY_H
#define QSSG_RENDER_SHADER_KEY_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtQuick3DUtils/private/qssgdataref_p.h>

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>

#include <QtQuick3DRuntimeRender/private/qssgrenderdefaultmaterial_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendertessmodevalues_p.h>

QT_BEGIN_NAMESPACE
// We have an ever expanding set of properties we like to hash into one or more 32 bit
// quantities.
// Furthermore we would like this set of properties to be convertable to string
// So the shader cache file itself is somewhat human readable/diagnosable.
// To do this we create a set of objects that act as properties to the master shader key.
// These objects are tallied in order to figure out their actual offset into the shader key's
// data store.  They are also run through in order to create the string shader cache key.

struct QSSGShaderKeyPropertyBase
{
    const char *name;
    quint32 offset;
    QSSGShaderKeyPropertyBase(const char *inName = "") : name(inName), offset(0) {}
    quint32 getOffset() const { return offset; }
    void setOffset(quint32 of) { offset = of; }

    template<quint32 TBitWidth>
    quint32 getMaskTemplate() const
    {
        quint32 bit = offset % 32;
        quint32 startValue = (1 << TBitWidth) - 1;
        quint32 mask = startValue << bit;
        return mask;
    }

    quint32 getIdx() const { return offset / 32; }

protected:
    void internalToString(QString &ioStr, const char *inBuffer) const
    {
        ioStr.append(QString::fromLocal8Bit(name));
        ioStr.append(QStringLiteral("="));
        ioStr.append(QString::fromLocal8Bit(inBuffer));
    }

    static void internalToString(QString &ioStr, const char *name, bool inValue)
    {
        if (inValue) {
            ioStr.append(QString::fromLocal8Bit(name));
            ioStr.append(QStringLiteral("="));
            ioStr.append(inValue ? QStringLiteral("true") : QStringLiteral("false"));
        }
    }
};

struct QSSGShaderKeyBoolean : public QSSGShaderKeyPropertyBase
{
    enum {
        BitWidth = 1,
    };

    QSSGShaderKeyBoolean(const char *inName = "") : QSSGShaderKeyPropertyBase(inName) {}

    quint32 getMask() const { return getMaskTemplate<BitWidth>(); }
    void setValue(QSSGDataRef<quint32> inDataStore, bool inValue) const
    {
        const qint32 idx = qint32(getIdx());
        Q_ASSERT(idx >= 0 && idx <= INT32_MAX);
        Q_ASSERT(inDataStore.size() > idx);
        quint32 mask = getMask();
        quint32 &target = inDataStore[idx];
        if (inValue) {
            target = target | mask;
        } else {
            mask = ~mask;
            target = target & mask;
        }
    }

    bool getValue(QSSGDataView<quint32> inDataStore) const
    {
        quint32 idx = getIdx();
        quint32 mask = getMask();
        const quint32 &target = inDataStore[idx];
        return (target & mask) ? true : false;
    }

    void toString(QString &ioStr, QSSGDataView<quint32> inKeySet) const
    {
        bool isHigh = getValue(inKeySet);
        internalToString(ioStr, name, isHigh);
    }
};

template<quint32 TBitWidth>
struct QSSGShaderKeyUnsigned : public QSSGShaderKeyPropertyBase
{
    enum {
        BitWidth = TBitWidth,
    };
    QSSGShaderKeyUnsigned(const char *inName = "") : QSSGShaderKeyPropertyBase(inName) {}
    quint32 getMask() const { return getMaskTemplate<BitWidth>(); }
    void setValue(QSSGDataRef<quint32> inDataStore, quint32 inValue) const
    {
        quint32 startValue = (1 << TBitWidth) - 1;
        // Ensure inValue is within range of bit width.
        inValue = inValue & startValue;
        quint32 bit = offset % 32;
        quint32 mask = getMask();
        quint32 idx = getIdx();
        inValue = inValue << bit;
        quint32 &target = inDataStore[idx];
        // Get rid of existing value
        quint32 inverseMask = ~mask;
        target = target & inverseMask;
        target = target | inValue;
    }

    quint32 getValue(QSSGDataView<quint32> inDataStore) const
    {
        quint32 idx = getIdx();
        quint32 bit = offset % 32;
        quint32 mask = getMask();
        const quint32 &target = inDataStore[idx];

        quint32 retval = target & mask;
        retval = retval >> bit;
        return retval;
    }

    void toString(QString &ioStr, QSSGDataView<quint32> inKeySet) const
    {
        quint32 value = getValue(inKeySet);
        char buf[64];
        toStr(value, toDataRef(buf, 64));
        internalToString(ioStr, buf);
    }

private:
    static quint32 toStr(quint32 item, QSSGDataRef<char> buffer)
    {
        // hope the buffer is big enough...
        return static_cast<quint32>(::snprintf(buffer.begin(), buffer.size(), "%u", item));
    }
};

struct QSSGShaderKeyTessellation : public QSSGShaderKeyUnsigned<4>
{
    enum TessellationBits {
        noTessellation = 1 << 0,
        linearTessellation = 1 << 1,
        phongTessellation = 1 << 2,
        npatchTessellation = 1 << 3
    };

    QSSGShaderKeyTessellation(const char *inName = "") : QSSGShaderKeyUnsigned<4>(inName) {}

    bool getBitValue(TessellationBits swizzleBit, QSSGDataView<quint32> inKeySet) const
    {
        return (getValue(inKeySet) & swizzleBit) ? true : false;
    }

    void setBitValue(TessellationBits swizzleBit, bool inValue, QSSGDataRef<quint32> inKeySet)
    {
        quint32 theValue = getValue(inKeySet);
        quint32 mask = swizzleBit;
        if (inValue) {
            theValue = theValue | mask;
        } else {
            mask = ~mask;
            theValue = theValue & mask;
        }
        setValue(inKeySet, theValue);
    }

    void setTessellationMode(QSSGDataRef<quint32> inKeySet, TessellationModeValues tessellationMode, bool val)
    {
        switch (tessellationMode) {
        case TessellationModeValues::NoTessellation:
            setBitValue(noTessellation, val, inKeySet);
            break;
        case TessellationModeValues::Linear:
            setBitValue(linearTessellation, val, inKeySet);
            break;
        case TessellationModeValues::NPatch:
            setBitValue(npatchTessellation, val, inKeySet);
            break;
        case TessellationModeValues::Phong:
            setBitValue(phongTessellation, val, inKeySet);
            break;
        }
    }

    bool isNoTessellation(QSSGDataView<quint32> inKeySet) const { return getBitValue(noTessellation, inKeySet); }
    void setNoTessellation(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(noTessellation, val, inKeySet); }

    bool isLinearTessellation(QSSGDataView<quint32> inKeySet) const
    {
        return getBitValue(linearTessellation, inKeySet);
    }
    void setLinearTessellation(QSSGDataRef<quint32> inKeySet, bool val)
    {
        setBitValue(linearTessellation, val, inKeySet);
    }

    bool isNPatchTessellation(QSSGDataView<quint32> inKeySet) const
    {
        return getBitValue(npatchTessellation, inKeySet);
    }
    void setNPatchTessellation(QSSGDataRef<quint32> inKeySet, bool val)
    {
        setBitValue(npatchTessellation, val, inKeySet);
    }

    bool isPhongTessellation(QSSGDataView<quint32> inKeySet) const
    {
        return getBitValue(phongTessellation, inKeySet);
    }
    void setPhongTessellation(QSSGDataRef<quint32> inKeySet, bool val)
    {
        setBitValue(phongTessellation, val, inKeySet);
    }

    void toString(QString &ioStr, QSSGDataView<quint32> inKeySet) const
    {
        ioStr.append(QString::fromLocal8Bit(name));
        ioStr.append(QStringLiteral("={"));
        internalToString(ioStr, "noTessellation", isNoTessellation(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "linearTessellation", isLinearTessellation(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "npatchTessellation", isNPatchTessellation(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "phongTessellation", isPhongTessellation(inKeySet));
        ioStr.append(QStringLiteral("}"));
    }
};

struct QSSGShaderKeyTextureChannel : public QSSGShaderKeyUnsigned<2>
{
    enum TexturChannelBits {
        R = 0,
        G = 1,
        B = 2,
        A = 3,
    };
    QSSGShaderKeyTextureChannel(const char *inName = "") : QSSGShaderKeyUnsigned<2>(inName) {}

    TexturChannelBits getTextureChannel(QSSGDataView<quint32> inKeySet) const
    {
        return TexturChannelBits(getValue(inKeySet));
    }

    void setTextureChannel(TexturChannelBits channel, QSSGDataRef<quint32> inKeySet)
    {
        setValue(inKeySet, quint32(channel));
    }
    const QString enumToStr[4] = {
        QStringLiteral("R"),
        QStringLiteral("G"),
        QStringLiteral("B"),
        QStringLiteral("A"),
    };
    void toString(QString &ioStr, QSSGDataView<quint32> inKeySet) const
    {
        ioStr.append(QString::fromLocal8Bit(name));
        ioStr.append(QStringLiteral("="));
        ioStr.append(enumToStr[getTextureChannel(inKeySet)]);
    }
};

struct QSSGShaderKeyTextureSwizzle : public QSSGShaderKeyUnsigned<5>
{
    enum TextureSwizzleBits {
        noSwizzle = 1 << 0,
        L8toR8 = 1 << 1,
        A8toR8 = 1 << 2,
        L8A8toRG8 = 1 << 3,
        L16toR16 = 1 << 4
    };

    QSSGShaderKeyTextureSwizzle(const char *inName = "") : QSSGShaderKeyUnsigned<5>(inName) {}

    bool getBitValue(TextureSwizzleBits swizzleBit, QSSGDataView<quint32> inKeySet) const
    {
        return (getValue(inKeySet) & swizzleBit) ? true : false;
    }

    void setBitValue(TextureSwizzleBits swizzleBit, bool inValue, QSSGDataRef<quint32> inKeySet)
    {
        quint32 theValue = getValue(inKeySet);
        quint32 mask = swizzleBit;
        if (inValue) {
            theValue = theValue | mask;
        } else {
            mask = ~mask;
            theValue = theValue & mask;
        }
        setValue(inKeySet, theValue);
    }

    void setSwizzleMode(QSSGDataRef<quint32> inKeySet, QSSGRenderTextureSwizzleMode swizzleMode, bool val)
    {
        switch (swizzleMode) {
        case QSSGRenderTextureSwizzleMode::NoSwizzle:
            setBitValue(noSwizzle, val, inKeySet);
            break;
        case QSSGRenderTextureSwizzleMode::L8toR8:
            setBitValue(L8toR8, val, inKeySet);
            break;
        case QSSGRenderTextureSwizzleMode::A8toR8:
            setBitValue(A8toR8, val, inKeySet);
            break;
        case QSSGRenderTextureSwizzleMode::L8A8toRG8:
            setBitValue(L8A8toRG8, val, inKeySet);
            break;
        case QSSGRenderTextureSwizzleMode::L16toR16:
            setBitValue(L16toR16, val, inKeySet);
            break;
        }
    }

    bool isNoSwizzled(QSSGDataView<quint32> inKeySet) const { return getBitValue(noSwizzle, inKeySet); }
    void setNoSwizzled(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(noSwizzle, val, inKeySet); }

    bool isL8Swizzled(QSSGDataView<quint32> inKeySet) const { return getBitValue(L8toR8, inKeySet); }
    void setL8Swizzled(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(L8toR8, val, inKeySet); }

    bool isA8Swizzled(QSSGDataView<quint32> inKeySet) const { return getBitValue(A8toR8, inKeySet); }
    void setA8Swizzled(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(A8toR8, val, inKeySet); }

    bool isL8A8Swizzled(QSSGDataView<quint32> inKeySet) const { return getBitValue(L8A8toRG8, inKeySet); }
    void setL8A8Swizzled(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(L8A8toRG8, val, inKeySet); }

    bool isL16Swizzled(QSSGDataView<quint32> inKeySet) const { return getBitValue(L16toR16, inKeySet); }
    void setL16Swizzled(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(L16toR16, val, inKeySet); }

    void toString(QString &ioStr, QSSGDataView<quint32> inKeySet) const
    {
        ioStr.append(QString::fromLocal8Bit(name));
        ioStr.append(QStringLiteral("={"));
        internalToString(ioStr, "noswizzle", isNoSwizzled(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "l8swizzle", isL8Swizzled(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "a8swizzle", isA8Swizzled(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "l8a8swizzle", isL8A8Swizzled(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "l16swizzle", isL16Swizzled(inKeySet));
        ioStr.append(QStringLiteral("}"));
    }
};

struct QSSGShaderKeyImageMap : public QSSGShaderKeyUnsigned<6>
{
    enum ImageMapBits {
        Enabled = 1 << 0,
        EnvMap = 1 << 1,
        LightProbe = 1 << 2,
        InvertUV = 1 << 3,
        Premultiplied = 1 << 4,
        Identity = 1 << 5
    };

    QSSGShaderKeyImageMap(const char *inName = "") : QSSGShaderKeyUnsigned<6>(inName) {}

    bool getBitValue(ImageMapBits imageBit, QSSGDataView<quint32> inKeySet) const
    {
        return (getValue(inKeySet) & imageBit) ? true : false;
    }

    void setBitValue(ImageMapBits imageBit, bool inValue, QSSGDataRef<quint32> inKeySet)
    {
        quint32 theValue = getValue(inKeySet);
        quint32 mask = imageBit;
        if (inValue) {
            theValue = theValue | mask;
        } else {
            mask = ~mask;
            theValue = theValue & mask;
        }
        setValue(inKeySet, theValue);
    }

    bool isEnabled(QSSGDataView<quint32> inKeySet) const { return getBitValue(Enabled, inKeySet); }
    void setEnabled(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(Enabled, val, inKeySet); }

    bool isEnvMap(QSSGDataView<quint32> inKeySet) const { return getBitValue(EnvMap, inKeySet); }
    void setEnvMap(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(EnvMap, val, inKeySet); }

    bool isLightProbe(QSSGDataView<quint32> inKeySet) const { return getBitValue(LightProbe, inKeySet); }
    void setLightProbe(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(LightProbe, val, inKeySet); }

    bool isInvertUVMap(QSSGDataView<quint32> inKeySet) const { return getBitValue(InvertUV, inKeySet); }
    void setInvertUVMap(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(InvertUV, val, inKeySet); }

    bool isPremultiplied(QSSGDataView<quint32> inKeySet) const { return getBitValue(Premultiplied, inKeySet); }
    void setPremultiplied(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(Premultiplied, val, inKeySet); }

    bool isIdentityTransform(QSSGDataView<quint32> inKeySet) const { return getBitValue(Identity, inKeySet); }
    void setIdentityTransform(QSSGDataRef<quint32> inKeySet, bool val) { setBitValue(Identity, val, inKeySet); }

    void toString(QString &ioStr, QSSGDataView<quint32> inKeySet) const
    {
        ioStr.append(QString::fromLocal8Bit(name));
        ioStr.append(QStringLiteral("={"));
        internalToString(ioStr, "enabled", isEnabled(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "envMap", isEnvMap(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "lightProbe", isLightProbe(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "invertUV", isInvertUVMap(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "premultiplied", isPremultiplied(inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "identity", isIdentityTransform(inKeySet));
        ioStr.append(QStringLiteral("}"));
    }
};

struct QSSGShaderKeySpecularModel : QSSGShaderKeyUnsigned<2>
{
    QSSGShaderKeySpecularModel(const char *inName = "") : QSSGShaderKeyUnsigned<2>(inName) {}

    void setSpecularModel(QSSGDataRef<quint32> inKeySet, QSSGRenderDefaultMaterial::MaterialSpecularModel inModel)
    {
        setValue(inKeySet, quint32(inModel));
    }

    QSSGRenderDefaultMaterial::MaterialSpecularModel getSpecularModel(QSSGDataView<quint32> inKeySet) const
    {
        return static_cast<QSSGRenderDefaultMaterial::MaterialSpecularModel>(getValue(inKeySet));
    }

    void toString(QString &ioStr, QSSGDataView<quint32> inKeySet) const
    {
        ioStr.append(QString::fromLocal8Bit(name));
        ioStr.append(QStringLiteral("="));
        switch (getSpecularModel(inKeySet)) {
        case QSSGRenderDefaultMaterial::MaterialSpecularModel::KGGX:
            ioStr.append(QStringLiteral("KGGX"));
            break;
        case QSSGRenderDefaultMaterial::MaterialSpecularModel::KWard:
            ioStr.append(QStringLiteral("KWard"));
            break;
        case QSSGRenderDefaultMaterial::MaterialSpecularModel::Default:
            ioStr.append(QStringLiteral("Default"));
            break;
        }
        ioStr.append(QStringLiteral(";"));
    }
};

struct QSSGShaderKeyAlphaMode : QSSGShaderKeyUnsigned<2>
{
    QSSGShaderKeyAlphaMode(const char *inName = "") : QSSGShaderKeyUnsigned<2>(inName) {}

    void setAlphaMode(QSSGDataRef<quint32> inKeySet, QSSGRenderDefaultMaterial::MaterialAlphaMode inMode)
    {
        setValue(inKeySet, quint32(inMode));
    }

    QSSGRenderDefaultMaterial::MaterialAlphaMode getAlphaMode(QSSGDataView<quint32> inKeySet) const
    {
        return static_cast<QSSGRenderDefaultMaterial::MaterialAlphaMode>(getValue(inKeySet));
    }

    void toString(QString &ioStr, QSSGDataView<quint32> inKeySet) const
    {
        ioStr.append(QString::fromLocal8Bit(name));
        ioStr.append(QStringLiteral("="));
        switch (getAlphaMode(inKeySet)) {
        case QSSGRenderDefaultMaterial::MaterialAlphaMode::Opaque:
            ioStr.append(QStringLiteral("Opaque"));
            break;
        case QSSGRenderDefaultMaterial::MaterialAlphaMode::Mask:
            ioStr.append(QStringLiteral("Mask"));
            break;
        case QSSGRenderDefaultMaterial::MaterialAlphaMode::Blend:
            ioStr.append(QStringLiteral("Blend"));
            break;
        case QSSGRenderDefaultMaterial::MaterialAlphaMode::Default:
            ioStr.append(QStringLiteral("Default"));
            break;
        }
        ioStr.append(QStringLiteral(";"));
    }
};

struct QSSGShaderKeyVertexAttribute : public QSSGShaderKeyUnsigned<7>
{
    enum VertexAttributeBits {
        Position = 1 << 0,
        Normal = 1 << 1,
        TexCoord0 = 1 << 2,
        TexCoord1 = 1 << 3,
        Tangent = 1 << 4,
        Binormal = 1 << 5,
        Color = 1 << 6
    };

    QSSGShaderKeyVertexAttribute(const char *inName = "") : QSSGShaderKeyUnsigned<7>(inName) {}

    bool getBitValue(VertexAttributeBits bit, QSSGDataView<quint32> inKeySet) const
    {
        return (getValue(inKeySet) & bit) ? true : false;
    }

    void toString(QString &ioStr, QSSGDataView<quint32> inKeySet) const
    {
        ioStr.append(QString::fromLocal8Bit(name));
        ioStr.append(QStringLiteral("={"));
        internalToString(ioStr, "position", getBitValue(Position, inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "normal", getBitValue(Normal, inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "texcoord0", getBitValue(TexCoord0, inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "texcoord1", getBitValue(TexCoord1, inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "tangent", getBitValue(Tangent, inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "binormal", getBitValue(Binormal, inKeySet));
        ioStr.append(QStringLiteral(";"));
        internalToString(ioStr, "color", getBitValue(Color, inKeySet));
        ioStr.append(QStringLiteral("}"));
    }
};

struct QSSGShaderDefaultMaterialKeyProperties
{
    enum {
        LightCount = 7,
    };
    enum {
        SingleChannelImageCount = 5,
    };
    enum ImageMapNames {
        DiffuseMap = 0,
        EmissiveMap,
        SpecularMap,
        BaseColorMap,
        BumpMap,
        SpecularAmountMap,
        NormalMap,
        DisplacementMap,
        LightmapIndirect,
        LightmapRadiosity,
        LightmapShadow,
        // single channel images
        OpacityMap,
        RoughnessMap,
        MetalnessMap,
        OcclusionMap,
        TranslucencyMap,

        ImageMapCount,
        SingleChannelImagesFirst = OpacityMap
    };
    enum ImageChannelNames {
        OpacityChannel = 0,
        RoughnessChannel,
        MetalnessChannel,
        OcclusionChannel,
        TranslucencyChannel,
    };

    QSSGShaderKeyBoolean m_hasLighting;
    QSSGShaderKeyBoolean m_hasIbl;
    QSSGShaderKeyUnsigned<3> m_lightCount;
    QSSGShaderKeyBoolean m_lightFlags[LightCount];
    QSSGShaderKeyBoolean m_lightSpotFlags[LightCount];
    QSSGShaderKeyBoolean m_lightAreaFlags[LightCount];
    QSSGShaderKeyBoolean m_lightShadowFlags[LightCount];
    QSSGShaderKeyBoolean m_specularEnabled;
    QSSGShaderKeyBoolean m_fresnelEnabled;
    QSSGShaderKeyBoolean m_vertexColorsEnabled;
    QSSGShaderKeySpecularModel m_specularModel;
    QSSGShaderKeyImageMap m_imageMaps[ImageMapCount];
    QSSGShaderKeyTextureSwizzle m_textureSwizzle[ImageMapCount];
    QSSGShaderKeyTextureChannel m_textureChannels[SingleChannelImageCount];
    QSSGShaderKeyTessellation m_tessellationMode;
    QSSGShaderKeyBoolean m_hasSkinning;
    QSSGShaderKeyBoolean m_wireframeMode;
    QSSGShaderKeyBoolean m_isDoubleSided;
    QSSGShaderKeyAlphaMode m_alphaMode;
    QSSGShaderKeyVertexAttribute m_vertexAttributes;

    QSSGShaderDefaultMaterialKeyProperties()
        : m_hasLighting("hasLighting")
        , m_hasIbl("hasIbl")
        , m_lightCount("lightCount")
        , m_specularEnabled("specularEnabled")
        , m_fresnelEnabled("fresnelEnabled")
        , m_vertexColorsEnabled("vertexColorsEnabled")
        , m_specularModel("specularModel")
        , m_tessellationMode("tessellationMode")
        , m_hasSkinning("hasSkinning")
        , m_wireframeMode("wireframeMode")
        , m_isDoubleSided("isDoubleSided")
        , m_alphaMode("alphaMode")
        , m_vertexAttributes("vertexAttributes")
    {
        m_lightFlags[0].name = "light0HasPosition";
        m_lightFlags[1].name = "light1HasPosition";
        m_lightFlags[2].name = "light2HasPosition";
        m_lightFlags[3].name = "light3HasPosition";
        m_lightFlags[4].name = "light4HasPosition";
        m_lightFlags[5].name = "light5HasPosition";
        m_lightFlags[6].name = "light6HasPosition";
        m_lightSpotFlags[0].name = "light0HasSpot";
        m_lightSpotFlags[1].name = "light1HasSpot";
        m_lightSpotFlags[2].name = "light2HasSpot";
        m_lightSpotFlags[3].name = "light3HasSpot";
        m_lightSpotFlags[4].name = "light4HasSpot";
        m_lightSpotFlags[5].name = "light5HasSpot";
        m_lightSpotFlags[6].name = "light6HasSpot";
        m_lightAreaFlags[0].name = "light0HasArea";
        m_lightAreaFlags[1].name = "light1HasArea";
        m_lightAreaFlags[2].name = "light2HasArea";
        m_lightAreaFlags[3].name = "light3HasArea";
        m_lightAreaFlags[4].name = "light4HasArea";
        m_lightAreaFlags[5].name = "light5HasArea";
        m_lightAreaFlags[6].name = "light6HasArea";
        m_lightShadowFlags[0].name = "light0HasShadow";
        m_lightShadowFlags[1].name = "light1HasShadow";
        m_lightShadowFlags[2].name = "light2HasShadow";
        m_lightShadowFlags[3].name = "light3HasShadow";
        m_lightShadowFlags[4].name = "light4HasShadow";
        m_lightShadowFlags[5].name = "light5HasShadow";
        m_lightShadowFlags[6].name = "light6HasShadow";

        m_imageMaps[0].name = "diffuseMap";
        m_imageMaps[1].name = "emissiveMap";
        m_imageMaps[2].name = "specularMap";
        m_imageMaps[3].name = "baseColorMap";
        m_imageMaps[4].name = "bumpMap";
        m_imageMaps[5].name = "specularAmountMap";
        m_imageMaps[6].name = "normalMap";
        m_imageMaps[7].name = "displacementMap";
        m_imageMaps[8].name = "lightmapIndirect";
        m_imageMaps[9].name = "lightmapRadiosity";
        m_imageMaps[10].name = "lightmapShadow";
        m_imageMaps[11].name = "opacityMap";
        m_imageMaps[12].name = "roughnessMap";
        m_imageMaps[13].name = "metalnessMap";
        m_imageMaps[14].name = "occlusionMap";
        m_imageMaps[15].name = "translucencyMap";

        m_textureSwizzle[0].name = "diffuseMap_swizzle";
        m_textureSwizzle[1].name = "emissiveMap_swizzle";
        m_textureSwizzle[2].name = "specularMap_swizzle";
        m_textureSwizzle[3].name = "baseColorMap_swizzle";
        m_textureSwizzle[4].name = "bumpMap_swizzle";
        m_textureSwizzle[5].name = "specularAmountMap_swizzle";
        m_textureSwizzle[6].name = "normalMap_swizzle";
        m_textureSwizzle[7].name = "displacementMap_swizzle";
        m_textureSwizzle[8].name = "lightmapIndirect_swizzle";
        m_textureSwizzle[9].name = "lightmapRadiosity_swizzle";
        m_textureSwizzle[10].name = "lightmapShadow_swizzle";
        m_textureSwizzle[11].name = "opacityMap_swizzle";
        m_textureSwizzle[12].name = "roughnessMap_swizzle";
        m_textureSwizzle[13].name = "metalnessMap_swizzle";
        m_textureSwizzle[14].name = "occlusionMap_swizzle";
        m_textureSwizzle[15].name = "translucencyMap_swizzle";

        m_textureChannels[0].name = "opacityMap_channel";
        m_textureChannels[1].name = "roughnessMap_channel";
        m_textureChannels[2].name = "metalnessMap_channel";
        m_textureChannels[3].name = "occlusionMap_channel";
        m_textureChannels[4].name = "translucencyMap_channel";

        setPropertyOffsets();
    }

    template<typename TVisitor>
    void visitProperties(TVisitor &inVisitor)
    {
        inVisitor.visit(m_hasLighting);
        inVisitor.visit(m_hasIbl);
        inVisitor.visit(m_lightCount);

        for (quint32 idx = 0, end = LightCount; idx < end; ++idx) {
            inVisitor.visit(m_lightFlags[idx]);
        }

        for (quint32 idx = 0, end = LightCount; idx < end; ++idx) {
            inVisitor.visit(m_lightSpotFlags[idx]);
        }

        for (quint32 idx = 0, end = LightCount; idx < end; ++idx) {
            inVisitor.visit(m_lightAreaFlags[idx]);
        }

        for (quint32 idx = 0, end = LightCount; idx < end; ++idx) {
            inVisitor.visit(m_lightShadowFlags[idx]);
        }

        inVisitor.visit(m_specularEnabled);
        inVisitor.visit(m_fresnelEnabled);
        inVisitor.visit(m_vertexColorsEnabled);
        inVisitor.visit(m_specularModel);

        for (quint32 idx = 0, end = ImageMapCount; idx < end; ++idx) {
            inVisitor.visit(m_imageMaps[idx]);
            inVisitor.visit(m_textureSwizzle[idx]);
        }
        for (quint32 idx = 0, end = SingleChannelImageCount; idx < end; ++idx)
            inVisitor.visit(m_textureChannels[idx]);

        inVisitor.visit(m_tessellationMode);
        inVisitor.visit(m_hasSkinning);
        inVisitor.visit(m_wireframeMode);
        inVisitor.visit(m_isDoubleSided);
        inVisitor.visit(m_alphaMode);
        inVisitor.visit(m_vertexAttributes);
    }

    struct OffsetVisitor
    {
        quint32 m_offset;
        OffsetVisitor() : m_offset(0) {}
        template<typename TPropType>
        void visit(TPropType &inProp)
        {
            // if we cross the 32 bit border we just move
            // to the next dword.
            // This cost a few extra bits but prevents tedious errors like
            // loosing shader key bits because they got moved beyond the 32 border
            quint32 bit = m_offset % 32;
            if (bit + TPropType::BitWidth > 31) {
                m_offset += 32 - bit;
            }

            inProp.setOffset(m_offset);
            m_offset += TPropType::BitWidth;
        }
    };

    void setPropertyOffsets()
    {
        OffsetVisitor visitor;
        visitProperties(visitor);
        // If this assert fires, then the default material key needs more bits.
        Q_ASSERT(visitor.m_offset < 320);
    }
};

struct QSSGShaderDefaultMaterialKey
{
    enum {
        DataBufferSize = 10,
    };
    quint32 m_dataBuffer[DataBufferSize];
    uint m_featureSetHash;

    QSSGShaderDefaultMaterialKey(uint inFeatureSetHash) : m_featureSetHash(inFeatureSetHash)
    {
        for (size_t idx = 0; idx < DataBufferSize; ++idx)
            m_dataBuffer[idx] = 0;
    }

    QSSGShaderDefaultMaterialKey() : m_featureSetHash(0)
    {
        for (size_t idx = 0; idx < DataBufferSize; ++idx)
            m_dataBuffer[idx] = 0;
    }

    uint hash() const
    {
        uint retval = 0;
        for (size_t idx = 0; idx < DataBufferSize; ++idx)
            retval = retval ^ qHash(m_dataBuffer[idx]);
        return retval ^ m_featureSetHash;
    }

    bool operator==(const QSSGShaderDefaultMaterialKey &other) const
    {
        bool retval = true;
        for (size_t idx = 0; idx < DataBufferSize && retval; ++idx)
            retval = m_dataBuffer[idx] == other.m_dataBuffer[idx];
        return retval && m_featureSetHash == other.m_featureSetHash;
    }

    // Cast operators to make getting properties easier.
    operator QSSGDataRef<quint32>() { return toDataRef(m_dataBuffer, DataBufferSize); }
    operator QSSGDataView<quint32>() const { return toDataView(m_dataBuffer, DataBufferSize); }

    struct StringVisitor
    {
        QByteArray &m_str;
        QSSGDataView<quint32> m_keyStore;
        StringVisitor(QByteArray &s, QSSGDataView<quint32> ks) : m_str(s), m_keyStore(ks) {}
        template<typename TPropType>
        void visit(const TPropType &prop)
        {
            quint32 originalSize = m_str.size();
            if (m_str.size())
                m_str.append(";");
            QString str = QString::fromUtf8(m_str);
            prop.toString(str, m_keyStore);
            // if the only thing we added was the semicolon
            // then nuke the semicolon
            m_str = str.toLocal8Bit();
            if (originalSize && m_str.size() == int(originalSize + 1))
                m_str.resize(int(originalSize));
        }
    };

    void toString(QByteArray &ioString, QSSGShaderDefaultMaterialKeyProperties &inProperties) const
    {
        StringVisitor theVisitor(ioString, *this);
        inProperties.visitProperties(theVisitor);
    }
};

Q_STATIC_ASSERT(std::is_trivially_destructible<QSSGShaderDefaultMaterialKey>::value);


inline uint qHash(const QSSGShaderDefaultMaterialKey &key)
{
    return key.hash();
}

QT_END_NAMESPACE

#endif

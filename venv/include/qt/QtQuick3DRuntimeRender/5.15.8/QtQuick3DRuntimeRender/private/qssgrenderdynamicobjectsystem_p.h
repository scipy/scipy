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

#ifndef QSSG_RENDER_DYNAMIC_OBJECT_SYSTEM_H
#define QSSG_RENDER_DYNAMIC_OBJECT_SYSTEM_H

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

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>

#include <QtQuick3DRuntimeRender/private/qssgrendershadercache_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendertessmodevalues_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendergraphobject_p.h>

#include <QtGui/QVector2D>

#include <QtCore/QString>
#include <QtCore/qmutex.h>

QT_BEGIN_NAMESPACE
struct QSSGRenderDynamicGraphObject;
// struct SWriteBuffer;
// struct SStrRemapMap;
class QSSGRenderContextInterface;
struct QSSGDynamicObjectClass;

typedef QPair<QByteArray, QByteArray> TStrStrPair;

namespace dynamic {

struct QSSGDynamicShaderMapKey
{
    TStrStrPair m_name;
    ShaderFeatureSetList m_features;
    TessellationModeValues m_tessMode;
    bool m_wireframeMode;
    uint m_hashCode;
    QSSGDynamicShaderMapKey(const TStrStrPair &inName, const ShaderFeatureSetList &inFeatures, TessellationModeValues inTessMode, bool inWireframeMode)
        : m_name(inName), m_tessMode(inTessMode), m_wireframeMode(inWireframeMode)
    {
        for (int i = 0; i < inFeatures.size(); ++i) {
            m_features.append(inFeatures[i]);
        }

        m_hashCode = qHash(m_name) ^ hashShaderFeatureSet(m_features) ^ qHash(m_tessMode) ^ qHash(m_wireframeMode);
    }
    bool operator==(const QSSGDynamicShaderMapKey &inKey) const
    {
        return m_name == inKey.m_name && m_features == inKey.m_features && m_tessMode == inKey.m_tessMode
                && m_wireframeMode == inKey.m_wireframeMode;
    }
};

struct QSSGCommand;

struct QSSGDynamicShaderProgramFlags : public QSSGShaderCacheProgramFlags
{
    TessellationModeValues tessMode = TessellationModeValues::NoTessellation;
    bool wireframeMode = false;

    QSSGDynamicShaderProgramFlags() = default;
    QSSGDynamicShaderProgramFlags(TessellationModeValues inTessMode, bool inWireframeMode)
        : tessMode(inTessMode), wireframeMode(inWireframeMode)
    {
    }

    static const char *wireframeToString(bool inEnable)
    {
        return inEnable ? "wireframeMode:true" : "wireframeMode:false";
    }
};
}

struct QSSGDynamicObjectShaderInfo
{
    QByteArray m_type; ///< shader type (GLSL or HLSL)
    QByteArray m_version; ///< shader version (e.g. 330 vor GLSL)
    bool m_hasGeomShader;
    bool m_isComputeShader;

    QSSGDynamicObjectShaderInfo() : m_hasGeomShader(false), m_isComputeShader(false) {}
    QSSGDynamicObjectShaderInfo(const QByteArray &inType, const QByteArray &inVersion, bool inHasGeomShader, bool inIsComputeShader)
        : m_type(inType), m_version(inVersion), m_hasGeomShader(inHasGeomShader), m_isComputeShader(inIsComputeShader)
    {
    }
};

typedef QPair<QSSGRef<QSSGRenderShaderProgram>, dynamic::QSSGDynamicShaderProgramFlags> TShaderAndFlags;

struct QSSGDynamicObjectSystem
{
    typedef QHash<QByteArray, QByteArray> TPathDataMap;
    typedef QHash<QByteArray, QSSGDynamicObjectShaderInfo> TShaderInfoMap;
    typedef QSet<QString> TPathSet;
    typedef QHash<dynamic::QSSGDynamicShaderMapKey, TShaderAndFlags> TShaderMap;

    QSSGRenderContextInterface *m_context;
    TPathDataMap m_expandedFiles;
    TShaderMap m_shaderMap;
    TShaderInfoMap m_shaderInfoMap;
    QByteArray m_vertShader;
    QByteArray m_fragShader;
    QByteArray m_geometryShader;
    QByteArray m_shaderLibraryVersion;
    QString m_shaderLibraryPlatformDirectory;
    mutable QMutex m_propertyLoadMutex;
    QAtomicInt ref;

    static QString getShaderCodeLibraryDirectory();

    QSSGDynamicObjectSystem(QSSGRenderContextInterface *ctx);

    ~QSSGDynamicObjectSystem();

    void setShaderData(const QByteArray &inPath,
                       const QByteArray &inData,
                       const QByteArray &inShaderType,
                       const QByteArray &inShaderVersion,
                       bool inHasGeomShader,
                       bool inIsComputeShader);

    QByteArray getShaderCacheKey(const QByteArray &inId, const QByteArray &inProgramMacro, const dynamic::QSSGDynamicShaderProgramFlags &inFlags);

    void insertShaderHeaderInformation(QByteArray &theReadBuffer, const QByteArray &inPathToEffect);

    void doInsertShaderHeaderInformation(QByteArray &theReadBuffer, const QByteArray &inPathToEffect);

    QByteArray doLoadShader(const QByteArray &inPathToEffect);

    QStringList getParameters(const QString &str, int begin, int end);

    void insertSnapperDirectives(QString &str);

    QSSGRef<QSSGRenderShaderProgram> compileShader(const QByteArray &inId,
                                                       const QByteArray &inProgramSource,
                                                       const QByteArray &inGeomSource,
                                                       const QByteArray &inProgramMacroName,
                                                       const ShaderFeatureSetList &inFeatureSet,
                                                       const dynamic::QSSGDynamicShaderProgramFlags &inFlags,
                                                       bool inForceCompilation = false);

    // This just returns the custom material shader source without compiling
    QByteArray getShaderSource(const QByteArray &inPath);

    TShaderAndFlags getShaderProgram(const QByteArray &inPath,
                                     const QByteArray &inProgramMacro,
                                     const ShaderFeatureSetList &inFeatureSet,
                                     const dynamic::QSSGDynamicShaderProgramFlags &inFlags,
                                     bool inForceCompilation);

    TShaderAndFlags getDepthPrepassShader(const QByteArray &inPath, const QByteArray &inPMacro, const ShaderFeatureSetList &inFeatureSet);

    void setShaderCodeLibraryVersion(const QByteArray &version);

    QByteArray shaderCodeLibraryVersion();

    void setShaderCodeLibraryPlatformDirectory(const QString &directory);

    QString shaderCodeLibraryPlatformDirectory();
};

QT_END_NAMESPACE

#endif

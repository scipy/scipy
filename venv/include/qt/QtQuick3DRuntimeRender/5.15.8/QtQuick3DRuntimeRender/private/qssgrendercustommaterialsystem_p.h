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

#ifndef QSSG_RENDER_CUSTOM_MATERIAL_SYSTEM_H
#define QSSG_RENDER_CUSTOM_MATERIAL_SYSTEM_H

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

#include <QtQuick3DRuntimeRender/private/qtquick3druntimerenderglobal_p.h>

#include <QtQuick3DRuntimeRender/private/qssgrenderdynamicobjectsystem_p.h>
#include <QtQuick3DRuntimeRender/private/qssgvertexpipelineimpl_p.h>

#include <QtCore/qhash.h>

#include <QtQuick3DRuntimeRender/private/qssgrendercustommaterial_p.h> // Make it possible to forward declare the nested TextureProperty

QT_BEGIN_NAMESPACE

namespace dynamic {
struct QSSGCommand;
}

struct QSSGCustomMaterialRenderContext;
struct QSSGRenderCustomMaterial;
class QSSGMaterialSystem;
struct QSSGRenderSubset;
struct QSSGShaderMapKey;
struct QSSGRenderCustomMaterialShader;
struct QSSGMaterialClass;
struct QSSGCustomMaterialTextureData;
struct QSSGRenderCustomMaterialBuffer;
struct QSSGMaterialOrComputeShader;
namespace dynamic {
struct QSSGBindShader;
struct QSSGApplyInstanceValue;
struct QSSGApplyBlending;
struct QSSGAllocateBuffer;
struct QSSGApplyBufferValue;
struct QSSGBindBuffer;
struct QSSGApplyBlitFramebuffer;
struct QSSGApplyRenderState;
}

// How to handle blend modes?
struct QSSGRenderModel;

class Q_QUICK3DRUNTIMERENDER_EXPORT QSSGMaterialSystem
{
    Q_DISABLE_COPY(QSSGMaterialSystem)
public:
    QAtomicInt ref;

private:
    typedef QHash<QSSGShaderMapKey, QSSGRef<QSSGRenderCustomMaterialShader>> ShaderMap;
    typedef QPair<QByteArray, QByteArray> TStrStrPair;
    typedef QPair<QByteArray, QSSGRef<QSSGCustomMaterialTextureData>> CustomMaterialTextureEntry;

    QSSGRenderContextInterface *context = nullptr;
    ShaderMap shaderMap;
    QVector<CustomMaterialTextureEntry> textureEntries;
    QVector<QSSGRenderCustomMaterialBuffer> allocatedBuffers;
    bool useFastBlits = true;
    QString shaderNameBuilder;
#ifdef QQ3D_UNUSED_TIMER
    QElapsedTimer lastFrameTime;
    float msSinceLastFrame = 0;
#endif // QQ3D_UNUSED_TIMER

    void releaseBuffer(qint32 inIdx);

    QSSGRef<QSSGRenderShaderProgram> getShader(QSSGCustomMaterialRenderContext &inRenderContext,
                                                   const QSSGRenderCustomMaterial &inMaterial,
                                                   const dynamic::QSSGBindShader &inCommand,
                                                   const ShaderFeatureSetList &inFeatureSet,
                                                   const dynamic::QSSGDynamicShaderProgramFlags &inFlags);

    QSSGMaterialOrComputeShader bindShader(QSSGCustomMaterialRenderContext &inRenderContext,
                                             const QSSGRenderCustomMaterial &inMaterial,
                                             const dynamic::QSSGBindShader &inCommand,
                                             const ShaderFeatureSetList &inFeatureSet);

    void doApplyInstanceValue(QSSGRenderCustomMaterial &inMaterial,
                              const QByteArray &propertyName,
                              const QVariant &propertyValue,
                              QSSGRenderShaderDataType inPropertyType,
                              const QSSGRef<QSSGRenderShaderProgram> &inShader);

    void applyInstanceValue(QSSGRenderCustomMaterial &inMaterial,
                            const QSSGRef<QSSGRenderShaderProgram> &inShader,
                            const dynamic::QSSGApplyInstanceValue &inCommand);

    void applyBlending(const dynamic::QSSGApplyBlending &inCommand);
    void applyCullMode(const dynamic::QSSGApplyCullMode &inCommand);

    void applyRenderStateValue(const dynamic::QSSGApplyRenderState &inCommand);

    // we currently only bind a source texture
    QSSGRef<QSSGRenderTexture2D> applyBufferValue(const QSSGRenderCustomMaterial &inMaterial,
                                                            const QSSGRef<QSSGRenderShaderProgram> &inShader,
                                                            const dynamic::QSSGApplyBufferValue &inCommand,
                                                            const QSSGRef<QSSGRenderTexture2D> &inSourceTexture);

    void allocateBuffer(const dynamic::QSSGAllocateBuffer &inCommand, const QSSGRef<QSSGRenderFrameBuffer> &inTarget);

    QSSGRef<QSSGRenderFrameBuffer> bindBuffer(const QSSGRenderCustomMaterial &inMaterial,
                                                  const dynamic::QSSGBindBuffer &inCommand,
                                                  bool &outClearTarget,
                                                  QVector2D &outDestSize);
    void computeScreenCoverage(QSSGCustomMaterialRenderContext &inRenderContext, qint32 *xMin, qint32 *yMin, qint32 *xMax, qint32 *yMax);
    void blitFramebuffer(QSSGCustomMaterialRenderContext &inRenderContext,
                         const dynamic::QSSGApplyBlitFramebuffer &inCommand,
                         const QSSGRef<QSSGRenderFrameBuffer> &inTarget);
    QSSGLayerGlobalRenderProperties getLayerGlobalRenderProperties(QSSGCustomMaterialRenderContext &inRenderContext);
    void renderPass(QSSGCustomMaterialRenderContext &inRenderContext,
                    const QSSGRef<QSSGRenderCustomMaterialShader> &inShader,
                    const QSSGRef<QSSGRenderTexture2D> & /* inSourceTexture */,
                    const QSSGRef<QSSGRenderFrameBuffer> &inFrameBuffer,
                    bool inRenderTargetNeedsClear,
                    const QSSGRef<QSSGRenderInputAssembler> &inAssembler,
                    quint32 inCount,
                    quint32 inOffset,
                    bool applyCullMode);
    void doRenderCustomMaterial(QSSGCustomMaterialRenderContext &inRenderContext,
                                const QSSGRenderCustomMaterial &inMaterial,
                                const ShaderFeatureSetList &inFeatureSet);
    void prepareDisplacementForRender(QSSGRenderCustomMaterial &inMaterial);
    void prepareMaterialForRender(QSSGRenderCustomMaterial &inMaterial);

    qint32 findBuffer(const QByteArray &inName) const;
    bool textureNeedsMips(const QSSGRenderCustomMaterial::TextureProperty *inPropDec, QSSGRenderTexture2D *inTexture);
    void setTexture(const QSSGRef<QSSGRenderShaderProgram> &inShader,
                    const QByteArray &inPropName,
                    const QSSGRef<QSSGRenderTexture2D> &inTexture,
                    const QSSGRenderCustomMaterial::TextureProperty *inPropDec = nullptr,
                    bool needMips = false);

public:
    QSSGMaterialSystem(QSSGRenderContextInterface *ct);

    ~QSSGMaterialSystem();

    void setMaterialClassShader(const QByteArray &inName,
                                const QByteArray &inShaderType,
                                const QByteArray &inShaderVersion,
                                const QByteArray &inShaderData,
                                bool inHasGeomShader,
                                bool inIsComputeShader);

    void setRenderContextInterface(QSSGRenderContextInterface *inContext);

    // Returns true if the material is dirty and thus will produce a different render result
    // than previously.  This effects things like progressive AA.
    bool prepareForRender(const QSSGRenderModel &inModel,
                          const QSSGRenderSubset &inSubset,
                          QSSGRenderCustomMaterial &inMaterial);

    bool renderDepthPrepass(const QMatrix4x4 &inMVP, const QSSGRenderCustomMaterial &inMaterial, const QSSGRenderSubset &inSubset);
    void renderSubset(QSSGCustomMaterialRenderContext &inRenderContext, const ShaderFeatureSetList &inFeatureSet);

    // get shader name
    QByteArray getShaderName(const QSSGRenderCustomMaterial &inMaterial);
    // apply property values
    void applyShaderPropertyValues(const QSSGRenderCustomMaterial &inMaterial, const QSSGRef<QSSGRenderShaderProgram> &inProgram);
    // Called by the uiccontext so this system can clear any per-frame render information.
    void endFrame();
};

struct Q_QUICK3DRUNTIMERENDER_EXPORT QSSGCustomMaterialVertexPipeline : public QSSGVertexPipelineImpl
{
    QSSGRenderContextInterface *m_context;
    TessellationModeValues m_tessMode;

    QSSGCustomMaterialVertexPipeline(QSSGRenderContextInterface *inContext, TessellationModeValues inTessMode);
    void initializeTessControlShader();
    void initializeTessEvaluationShader();
    void finalizeTessControlShader();
    void finalizeTessEvaluationShader();

    // Responsible for beginning all vertex and fragment generation (void main() { etc).
    virtual void beginVertexGeneration(const QSSGShaderDefaultMaterialKey &inKey, quint32 displacementImageIdx, QSSGRenderableImage *displacementImage) override;
    // The fragment shader expects a floating point constant, objectOpacity to be defined
    // post this method.
    virtual void beginFragmentGeneration() override;
    // Output variables may be mangled in some circumstances so the shader generation
    // system needs an abstraction mechanism around this.
    virtual void assignOutput(const QByteArray &inVarName, const QByteArray &inVarValue) override;
    virtual void generateEnvMapReflection(const QSSGShaderDefaultMaterialKey &inKey) override { Q_UNUSED(inKey); }
    virtual void generateViewVector() override {}
    virtual void generateUVCoords(const QSSGShaderDefaultMaterialKey &inKey, quint32 inUVSet) override;
    virtual void generateWorldNormal(const QSSGShaderDefaultMaterialKey &inKey) override;
    virtual void generateObjectNormal() override;
    virtual void generateVarTangentAndBinormal(const QSSGShaderDefaultMaterialKey &inKey) override;
    virtual void generateWorldPosition() override;
    // responsible for closing all vertex and fragment generation
    virtual void endVertexGeneration(bool customShader) override;
    virtual void endFragmentGeneration(bool customShader) override;
    virtual QSSGShaderStageGeneratorInterface &activeStage() override;
    virtual void addInterpolationParameter(const QByteArray &inName, const QByteArray &inType) override;
    virtual void doGenerateUVCoords(quint32 inUVSet) override;
    virtual void doGenerateWorldNormal() override;
    virtual void doGenerateObjectNormal() override;
    virtual void doGenerateWorldPosition() override;
    virtual void doGenerateVarTangent() override;
    virtual void doGenerateVarBinormal() override;
    virtual void doGenerateVertexColor(const QSSGShaderDefaultMaterialKey &inKey) override;
};
QT_END_NAMESPACE
#endif

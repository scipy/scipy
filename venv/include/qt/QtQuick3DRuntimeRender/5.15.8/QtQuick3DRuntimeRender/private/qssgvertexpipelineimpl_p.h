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

#ifndef QSSG_VERTEX_PIPELINE_IMPL_H
#define QSSG_VERTEX_PIPELINE_IMPL_H

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

#include <QtQuick3DRuntimeRender/private/qssgrenderdefaultmaterialshadergenerator_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendertessmodevalues_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendershaderkeys_p.h>

#include <QtCore/QSharedPointer>

QT_BEGIN_NAMESPACE
// Baseclass for the vertex pipelines to be sure we have consistent implementations.
struct QSSGVertexPipelineImpl : public QSSGDefaultMaterialVertexPipelineInterface
{
    enum class GenerationFlag
    {
        UVCoords = 1,
        EnvMapReflection = 1 << 1,
        ViewVector = 1 << 2,
        WorldNormal = 1 << 3,
        ObjectNormal = 1 << 4,
        WorldPosition = 1 << 5,
        TangentBinormal = 1 << 6,
        UVCoords1 = 1 << 7,
        VertexColor = 1 << 8,
    };

    typedef TStrTableStrMap::const_iterator TParamIter;
    typedef QFlags<GenerationFlag> GenerationFlags;

    QSSGRef<QSSGMaterialShaderGeneratorInterface> m_materialGenerator;
    QSSGRef<QSSGShaderProgramGeneratorInterface> m_programGenerator;
    QString m_tempString;

    GenerationFlags m_generationFlags;
    bool m_wireframe;
    TStrTableStrMap m_interpolationParameters;
    quint32 m_displacementIdx;
    QSSGRenderableImage *m_displacementImage;
    QList<QByteArray> m_addedFunctions;

    QSSGVertexPipelineImpl(const QSSGRef<QSSGMaterialShaderGeneratorInterface> &inMaterial,
                           const QSSGRef<QSSGShaderProgramGeneratorInterface> &inProgram,
                           bool inWireframe /* only works if tessellation is true */)

        : m_materialGenerator(inMaterial)
        , m_programGenerator(inProgram)
        , m_wireframe(inWireframe)
        , m_displacementIdx(0)
        , m_displacementImage(nullptr)
    {
    }

    // Trues true if the code was *not* set.
    bool setCode(GenerationFlag inCode)
    {
        if (m_generationFlags & inCode)
            return true;
        m_generationFlags |= inCode;
        return false;
    }
    bool hasCode(GenerationFlag inCode) { return (m_generationFlags & inCode); }
    QSSGRef<QSSGShaderProgramGeneratorInterface> programGenerator() { return m_programGenerator; }

    QSSGShaderStageGeneratorInterface &vertex()
    {
        return *programGenerator()->getStage(QSSGShaderGeneratorStage::Vertex);
    }
    QSSGShaderStageGeneratorInterface &tessControl()
    {
        return *programGenerator()->getStage(QSSGShaderGeneratorStage::TessControl);
    }
    QSSGShaderStageGeneratorInterface &tessEval()
    {
        return *programGenerator()->getStage(QSSGShaderGeneratorStage::TessEval);
    }
    QSSGShaderStageGeneratorInterface &geometry()
    {
        return *programGenerator()->getStage(QSSGShaderGeneratorStage::Geometry);
    }
    QSSGShaderStageGeneratorInterface &fragment()
    {
        return *programGenerator()->getStage(QSSGShaderGeneratorStage::Fragment);
    }
    QSSGRef<QSSGMaterialShaderGeneratorInterface> materialGenerator() { return m_materialGenerator; }

    void setupDisplacement(quint32 displacementImageIdx, QSSGRenderableImage *displacementImage)
    {
        m_displacementIdx = displacementImageIdx;
        m_displacementImage = displacementImage;
    }

    bool hasTessellation() const { return m_programGenerator->getEnabledStages() & QSSGShaderGeneratorStage::TessEval; }
    bool hasGeometryStage() const { return m_programGenerator->getEnabledStages() & QSSGShaderGeneratorStage::Geometry; }
    bool hasDisplacment() const { return m_displacementImage != nullptr; }

    void initializeWireframeGeometryShader()
    {
        if (m_wireframe && programGenerator()->getStage(QSSGShaderGeneratorStage::Geometry)
            && programGenerator()->getStage(QSSGShaderGeneratorStage::TessEval)) {
            QSSGShaderStageGeneratorInterface &geometryShader(*programGenerator()->getStage(QSSGShaderGeneratorStage::Geometry));
            // currently geometry shader is only used for drawing wireframe
            if (m_wireframe) {
                geometryShader.addUniform("viewportMatrix", "mat4");
                geometryShader.addOutgoing("varEdgeDistance", "vec3");
                geometryShader.append("layout (triangles) in;");
                geometryShader.append("layout (triangle_strip, max_vertices = 3) out;");
                geometryShader.append("void main() {");

                // how this all work see
                // http://developer.download.nvidia.com/SDK/10.5/direct3d/Source/SolidWireframe/Doc/SolidWireframe.pdf

                geometryShader.append("// project points to screen space\n"
                                      "    vec3 p0 = vec3(viewportMatrix * (gl_in[0].gl_Position / "
                                      "gl_in[0].gl_Position.w));\n"
                                      "    vec3 p1 = vec3(viewportMatrix * (gl_in[1].gl_Position / "
                                      "gl_in[1].gl_Position.w));\n"
                                      "    vec3 p2 = vec3(viewportMatrix * (gl_in[2].gl_Position / "
                                      "gl_in[2].gl_Position.w));\n"
                                      "// compute triangle heights\n"
                                      "    float e1 = length(p1 - p2);\n"
                                      "    float e2 = length(p2 - p0);\n"
                                      "    float e3 = length(p1 - p0);\n"
                                      "    float alpha = acos( (e2*e2 + e3*e3 - e1*e1) / (2.0*e2*e3) );\n"
                                      "    float beta = acos( (e1*e1 + e3*e3 - e2*e2) / (2.0*e1*e3) );\n"
                                      "    float ha = abs( e3 * sin( beta ) );\n"
                                      "    float hb = abs( e3 * sin( alpha ) );\n"
                                      "    float hc = abs( e2 * sin( alpha ) );\n");
            }
        }
    }

    void finalizeWireframeGeometryShader()
    {
        QSSGShaderStageGeneratorInterface &geometryShader(*programGenerator()->getStage(QSSGShaderGeneratorStage::Geometry));

        if (m_wireframe == true && programGenerator()->getStage(QSSGShaderGeneratorStage::Geometry)
            && programGenerator()->getStage(QSSGShaderGeneratorStage::TessEval)) {
            const char *theExtension("TE[");
            // we always assume triangles
            for (int i = 0; i < 3; i++) {
                char buf[10];
                sprintf(buf, "%d", i);
                for (TStrTableStrMap::iterator iter = m_interpolationParameters.begin(), end = m_interpolationParameters.end();
                     iter != end;
                     ++iter) {
                    geometryShader << "    " << iter.key() << " = " << iter.key() << theExtension << buf << "];\n";
                }

                geometryShader << "    gl_Position = gl_in[" << buf << "].gl_Position;\n";
                // the triangle distance is interpolated through the shader stage
                if (i == 0) {
                    geometryShader << "\n    varEdgeDistance = vec3(ha*"
                                   << "gl_in[" << buf << "].gl_Position.w, 0.0, 0.0);\n";
                } else if (i == 1) {
                    geometryShader << "\n    varEdgeDistance = vec3(0.0, hb*"
                                   << "gl_in[" << buf << "].gl_Position.w, 0.0);\n";
                } else if (i == 2) {
                    geometryShader << "\n    varEdgeDistance = vec3(0.0, 0.0, hc*"
                                   << "gl_in[" << buf << "].gl_Position.w);\n";
                }

                // submit vertex
                geometryShader << "    EmitVertex();\n";
            }
            // end primitive
            geometryShader << "    EndPrimitive();\n";
        }
    }

    virtual void setupTessIncludes(QSSGShaderGeneratorStage inStage, TessellationModeValues inTessMode)
    {
        QSSGShaderStageGeneratorInterface &tessShader(*programGenerator()->getStage(inStage));

        // depending on the selected tessellation mode chose program
        switch (inTessMode) {
        case TessellationModeValues::Phong:
            tessShader.addInclude("tessellationPhong.glsllib");
            break;
        case TessellationModeValues::NPatch:
            tessShader.addInclude("tessellationNPatch.glsllib");
            break;
        default:
            Q_ASSERT(false); // fallthrough intentional
        case TessellationModeValues::Linear:
            tessShader.addInclude("tessellationLinear.glsllib");
            break;
        }
    }

    void generateUVCoords(const QSSGShaderDefaultMaterialKey &inKey, quint32 inUVSet = 0) override
    {
        if (inUVSet == 0 && setCode(GenerationFlag::UVCoords))
            return;
        if (inUVSet == 1 && setCode(GenerationFlag::UVCoords1))
            return;

        Q_ASSERT(inUVSet == 0 || inUVSet == 1);

        if (inUVSet == 0) {
            if (hasAttributeInKey(QSSGShaderKeyVertexAttribute::TexCoord0, inKey)) {
                addInterpolationParameter("varTexCoord0", "vec2");
                doGenerateUVCoords(inUVSet);
            } else {
                fragment() << "    vec2 varTexCoord0 = vec2(0.0);\n";
            }
        } else if (inUVSet == 1) {
            if (hasAttributeInKey(QSSGShaderKeyVertexAttribute::TexCoord1, inKey)) {
                addInterpolationParameter("varTexCoord1", "vec2");
                doGenerateUVCoords(inUVSet);
            } else {
                fragment() << "    vec2 varTexCoord1 = vec2(0.0);\n";
            }
        }
    }
    void generateEnvMapReflection(const QSSGShaderDefaultMaterialKey &inKey) override
    {
        if (setCode(GenerationFlag::EnvMapReflection))
            return;

        generateWorldPosition();
        generateWorldNormal(inKey);
        QSSGShaderStageGeneratorInterface &activeGenerator(activeStage());
        activeGenerator.addInclude("viewProperties.glsllib");
        addInterpolationParameter("var_object_to_camera", "vec3");

        activeGenerator.append("    var_object_to_camera = normalize( local_model_world_position "
                               "- cameraPosition );");

        // World normal cannot be relied upon in the vertex shader because of bump maps.
        fragment().append("    vec3 environment_map_reflection = reflect( "
                          "normalize(var_object_to_camera), world_normal.xyz );");
        fragment().append("    environment_map_reflection *= vec3( 0.5, 0.5, 0 );");
        fragment().append("    environment_map_reflection += vec3( 0.5, 0.5, 1.0 );");
    }
    void generateViewVector() override
    {
        if (setCode(GenerationFlag::ViewVector))
            return;
        generateWorldPosition();
        QSSGShaderStageGeneratorInterface &activeGenerator(activeStage());
        activeGenerator.addInclude("viewProperties.glsllib");
        addInterpolationParameter("varViewVector", "vec3");

        activeGenerator.append("    vec3 local_view_vector = normalize(cameraPosition - "
                               "local_model_world_position);");
        assignOutput("varViewVector", "local_view_vector");
        fragment() << "    vec3 view_vector = normalize(varViewVector);\n";
    }

    // fragment shader expects varying vertex normal
    // lighting in vertex pipeline expects world_normal
    void generateWorldNormal(const QSSGShaderDefaultMaterialKey &inKey) override
    {
        if (setCode(GenerationFlag::WorldNormal))
            return;

        if (hasAttributeInKey(QSSGShaderKeyVertexAttribute::Normal, inKey)) {
            addInterpolationParameter("varNormal", "vec3");
            doGenerateWorldNormal();
        } else {
            generateWorldPosition();
            fragment().append("    vec3 varNormal = cross(dFdx(varWorldPos), dFdy(varWorldPos));");
        }
        fragment().append("    vec3 world_normal = normalize( varNormal );");
    }
    void generateObjectNormal() override
    {
        if (setCode(GenerationFlag::ObjectNormal))
            return;
        doGenerateObjectNormal();
        fragment().append("    vec3 object_normal = normalize(varObjectNormal);");
    }
    void generateWorldPosition() override
    {
        if (setCode(GenerationFlag::WorldPosition))
            return;

        activeStage().addUniform("modelMatrix", "mat4");
        addInterpolationParameter("varWorldPos", "vec3");
        doGenerateWorldPosition();
    }
    void generateVarTangentAndBinormal(const QSSGShaderDefaultMaterialKey &inKey) override
    {
        if (setCode(GenerationFlag::TangentBinormal))
            return;

        // I assume that there is no mesh having only binormal without tangent
        // since it is an abnormal case
        if (hasAttributeInKey(QSSGShaderKeyVertexAttribute::Tangent, inKey)) {
            const bool hasBinormal = hasAttributeInKey(QSSGShaderKeyVertexAttribute::Binormal, inKey);
            addInterpolationParameter("varTangent", "vec3");
            doGenerateVarTangent();
            fragment() << "    vec3 tangent = normalize(varTangent);\n";

            if (hasBinormal) {
                addInterpolationParameter("varBinormal", "vec3");
                doGenerateVarBinormal();
                fragment() << "    vec3 binormal = normalize(varBinormal);\n";
            } else {
                fragment() << "    vec3 binormal = vec3(0.0);\n";
            }
        } else {
            fragment() << "    vec3 tangent = vec3(0.0);\n"
                       << "    vec3 binormal = vec3(0.0);\n";
        }
    }
    void generateVertexColor(const QSSGShaderDefaultMaterialKey &inKey) override
    {
        if (setCode(GenerationFlag::VertexColor))
            return;
        addInterpolationParameter("varColor", "vec4");
        doGenerateVertexColor(inKey);
        fragment().append("    vec4 vertColor = varColor;");
    }

    bool hasActiveWireframe() override { return m_wireframe; }

    void addIncoming(const QByteArray &name, const QByteArray &type) override { activeStage().addIncoming(name, type); }

    void addOutgoing(const QByteArray &name, const QByteArray &type) override { addInterpolationParameter(name, type); }

    void addUniform(const QByteArray &name, const QByteArray &type) override { activeStage().addUniform(name, type); }

    void addInclude(const QByteArray &name) override { activeStage().addInclude(name); }

    void addFunction(const QByteArray &functionName) override
    {
        if (!m_addedFunctions.contains(functionName)) {
            m_addedFunctions.push_back(functionName);
            QByteArray includeName = "func" + functionName + ".glsllib";
            addInclude(includeName);
        }
    }

    void addConstantBuffer(const QByteArray &name, const QByteArray &layout) override
    {
        activeStage().addConstantBuffer(name, layout);
    }
    void addConstantBufferParam(const QByteArray &cbName, const QByteArray &paramName, const QByteArray &type) override
    {
        activeStage().addConstantBufferParam(cbName, paramName, type);
    }

    QSSGShaderStageGeneratorInterface &operator<<(const QByteArray &data) override
    {
        activeStage() << data;
        return *this;
    }

    void append(const QByteArray &data) override { activeStage().append(data); }

    QSSGShaderGeneratorStage stage() const override
    {
        return const_cast<QSSGVertexPipelineImpl *>(this)->activeStage().stage();
    }

    void beginVertexGeneration(const QSSGShaderDefaultMaterialKey &inKey, quint32 displacementImageIdx, QSSGRenderableImage *displacementImage) override = 0;
    void assignOutput(const QByteArray &inVarName, const QByteArray &inVarValueExpr) override = 0;
    void endVertexGeneration(bool customShader) override = 0;

    void beginFragmentGeneration() override = 0;
    void endFragmentGeneration(bool customShader) override = 0;

    virtual QSSGShaderStageGeneratorInterface &activeStage() = 0;
    virtual void addInterpolationParameter(const QByteArray &inParamName, const QByteArray &inParamType) = 0;

    virtual void doGenerateUVCoords(quint32 inUVSet) = 0;
    virtual void doGenerateWorldNormal() = 0;
    virtual void doGenerateObjectNormal() = 0;
    virtual void doGenerateWorldPosition() = 0;
    virtual void doGenerateVarTangent() = 0;
    virtual void doGenerateVarBinormal() = 0;
    virtual void doGenerateVertexColor(const QSSGShaderDefaultMaterialKey &inKey) = 0;
    virtual bool hasAttributeInKey(QSSGShaderKeyVertexAttribute::VertexAttributeBits inAttr, const QSSGShaderDefaultMaterialKey &inKey) {
        // it returns true by default
        Q_UNUSED(inAttr)
        Q_UNUSED(inKey)
        return true;
    }
};
QT_END_NAMESPACE

#endif

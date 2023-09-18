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

#ifndef QSSG_RENDER_SHADER_CODE_GENERATOR_V2_H
#define QSSG_RENDER_SHADER_CODE_GENERATOR_V2_H

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
#include <QtQuick3DRuntimeRender/private/qssgrendershadercodegenerator_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendershadercache_p.h>

#include <QtCore/qstring.h>

QT_BEGIN_NAMESPACE
// So far the generator is only useful for graphics stages,
// it doesn't seem useful for compute stages.
enum class QSSGShaderGeneratorStage
{
    None = 0,
    Vertex = 1,
    TessControl = 1 << 1,
    TessEval = 1 << 2,
    Geometry = 1 << 3,
    Fragment = 1 << 4,
    StageCount = 5,
};

Q_DECLARE_FLAGS(QSSGShaderGeneratorStageFlags, QSSGShaderGeneratorStage)
Q_DECLARE_OPERATORS_FOR_FLAGS(QSSGShaderGeneratorStageFlags)

class Q_QUICK3DRUNTIMERENDER_EXPORT QSSGShaderStageGeneratorInterface
{
protected:
    virtual ~QSSGShaderStageGeneratorInterface();

public:
    virtual void addIncoming(const QByteArray &name, const QByteArray &type) = 0;

    virtual void addOutgoing(const QByteArray &name, const QByteArray &type) = 0;

    virtual void addUniform(const QByteArray &name, const QByteArray &type) = 0;

    virtual void addInclude(const QByteArray &name) = 0;

    virtual void addFunction(const QByteArray &functionName) = 0;

    virtual void addConstantBuffer(const QByteArray &name, const QByteArray &layout) = 0;
    virtual void addConstantBufferParam(const QByteArray &cbName, const QByteArray &paramName, const QByteArray &type) = 0;

    virtual QSSGShaderStageGeneratorInterface &operator<<(const QByteArray &data) = 0;
    virtual void append(const QByteArray &data) = 0;

    virtual QSSGShaderGeneratorStage stage() const = 0;
};

class QSSGRenderContextInterface;

class Q_QUICK3DRUNTIMERENDER_EXPORT QSSGShaderProgramGeneratorInterface
{
public:
    QAtomicInt ref;

    virtual ~QSSGShaderProgramGeneratorInterface() {}
    static QSSGShaderGeneratorStageFlags defaultFlags()
    {
        return QSSGShaderGeneratorStageFlags(QSSGShaderGeneratorStage::Vertex | QSSGShaderGeneratorStage::Fragment);
    }
    virtual void beginProgram(QSSGShaderGeneratorStageFlags inEnabledStages = defaultFlags()) = 0;

    virtual QSSGShaderGeneratorStageFlags getEnabledStages() const = 0;

    // get the stage or nullptr if it has not been created.
    virtual QSSGShaderStageGeneratorInterface *getStage(QSSGShaderGeneratorStage inStage) = 0;

    // Implicit call to end program.

    virtual QSSGRef<QSSGRenderShaderProgram> compileGeneratedShader(const QByteArray &inShaderName,
                                                                        const QSSGShaderCacheProgramFlags &inFlags,
                                                                        const ShaderFeatureSetList &inFeatureSet,
                                                                        bool separableProgram = false) = 0;

    QSSGRef<QSSGRenderShaderProgram> compileGeneratedShader(const QByteArray &inShaderName, bool separableProgram = false);

    static QSSGRef<QSSGShaderProgramGeneratorInterface> createProgramGenerator(QSSGRenderContextInterface *inContext);

    static void outputParaboloidDepthVertex(QSSGShaderStageGeneratorInterface &inGenerator);
    // By convention, the local space result of the TE is stored in vec4 pos local variable.
    // This function expects such state.
    static void outputParaboloidDepthTessEval(QSSGShaderStageGeneratorInterface &inGenerator);
    // Utilities shared among the various different systems.
    static void outputParaboloidDepthFragment(QSSGShaderStageGeneratorInterface &inGenerator);

    static void outputCubeFaceDepthVertex(QSSGShaderStageGeneratorInterface &inGenerator);
    static void outputCubeFaceDepthGeometry(QSSGShaderStageGeneratorInterface &inGenerator);
    static void outputCubeFaceDepthFragment(QSSGShaderStageGeneratorInterface &inGenerator);
};
QT_END_NAMESPACE
#endif

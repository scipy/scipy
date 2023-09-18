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

#ifndef QT3DEXTRAS_QNORMALDIFFUSEMAPMATERIAL_P_H
#define QT3DEXTRAS_QNORMALDIFFUSEMAPMATERIAL_P_H

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

#include <private/qmaterial_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QFilterKey;
class QEffect;
class QAbstractTexture;
class QTechnique;
class QParameter;
class QShaderProgram;
class QShaderProgramBuilder;
class QRenderPass;

} // namespace Qt3DRender

namespace Qt3DExtras {

class QNormalDiffuseMapMaterial;

class QNormalDiffuseMapMaterialPrivate: public Qt3DRender::QMaterialPrivate
{
public:
    QNormalDiffuseMapMaterialPrivate();

    virtual void init();

    void handleAmbientChanged(const QVariant &var);
    void handleDiffuseChanged(const QVariant &var);
    void handleNormalChanged(const QVariant &var);
    void handleSpecularChanged(const QVariant &var);
    void handleShininessChanged(const QVariant &var);
    void handleTextureScaleChanged(const QVariant &var);

    Qt3DRender::QEffect *m_normalDiffuseEffect;
    Qt3DRender::QAbstractTexture *m_diffuseTexture;
    Qt3DRender::QAbstractTexture *m_normalTexture;
    Qt3DRender::QParameter *m_ambientParameter;
    Qt3DRender::QParameter *m_diffuseParameter;
    Qt3DRender::QParameter *m_normalParameter;
    Qt3DRender::QParameter *m_specularParameter;
    Qt3DRender::QParameter *m_shininessParameter;
    Qt3DRender::QParameter *m_textureScaleParameter;
    Qt3DRender::QTechnique *m_normalDiffuseGL3Technique;
    Qt3DRender::QTechnique *m_normalDiffuseGL2Technique;
    Qt3DRender::QTechnique *m_normalDiffuseES2Technique;
    Qt3DRender::QTechnique *m_normalDiffuseRHITechnique;
    Qt3DRender::QRenderPass *m_normalDiffuseGL3RenderPass;
    Qt3DRender::QRenderPass *m_normalDiffuseGL2RenderPass;
    Qt3DRender::QRenderPass *m_normalDiffuseES2RenderPass;
    Qt3DRender::QRenderPass *m_normalDiffuseRHIRenderPass;
    Qt3DRender::QShaderProgram *m_normalDiffuseGL3Shader;
    Qt3DRender::QShaderProgramBuilder *m_normalDiffuseGL3ShaderBuilder;
    Qt3DRender::QShaderProgram *m_normalDiffuseGL2ES2Shader;
    Qt3DRender::QShaderProgramBuilder *m_normalDiffuseGL2ES2ShaderBuilder;
    Qt3DRender::QShaderProgram *m_normalDiffuseRHIShader;
    Qt3DRender::QShaderProgramBuilder *m_normalDiffuseRHIShaderBuilder;
    Qt3DRender::QFilterKey *m_filterKey;

    Q_DECLARE_PUBLIC(QNormalDiffuseMapMaterial)
};

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QNORMALDIFFUSEMAPMATERIAL_P_H


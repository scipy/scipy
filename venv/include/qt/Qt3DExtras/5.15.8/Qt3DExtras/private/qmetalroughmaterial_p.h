/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DEXTRAS_QMETALROUGHMATERIAL_P_H
#define QT3DEXTRAS_QMETALROUGHMATERIAL_P_H

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

#include <Qt3DRender/private/qmaterial_p.h>

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

class QMetalRoughMaterial;

class QMetalRoughMaterialPrivate : public Qt3DRender::QMaterialPrivate
{
public:
    QMetalRoughMaterialPrivate();

    void init();

    void handleTextureScaleChanged(const QVariant &var);

    Qt3DRender::QParameter *m_baseColorParameter;
    Qt3DRender::QParameter *m_metalnessParameter;
    Qt3DRender::QParameter *m_roughnessParameter;
    Qt3DRender::QParameter *m_baseColorMapParameter;
    Qt3DRender::QParameter *m_metalnessMapParameter;
    Qt3DRender::QParameter *m_roughnessMapParameter;
    Qt3DRender::QParameter *m_ambientOcclusionMapParameter;
    Qt3DRender::QParameter *m_normalMapParameter;
    Qt3DRender::QParameter *m_textureScaleParameter;
    Qt3DRender::QEffect *m_metalRoughEffect;
    Qt3DRender::QTechnique *m_metalRoughGL3Technique;
    Qt3DRender::QRenderPass *m_metalRoughGL3RenderPass;
    Qt3DRender::QShaderProgram *m_metalRoughGL3Shader;
    Qt3DRender::QShaderProgramBuilder *m_metalRoughGL3ShaderBuilder;
    Qt3DRender::QTechnique *m_metalRoughES3Technique;
    Qt3DRender::QRenderPass *m_metalRoughES3RenderPass;
    Qt3DRender::QShaderProgram *m_metalRoughES3Shader;
    Qt3DRender::QShaderProgramBuilder *m_metalRoughES3ShaderBuilder;
    Qt3DRender::QTechnique *m_metalRoughRHITechnique;
    Qt3DRender::QRenderPass *m_metalRoughRHIRenderPass;
    Qt3DRender::QShaderProgram *m_metalRoughRHIShader;
    Qt3DRender::QShaderProgramBuilder *m_metalRoughRHIShaderBuilder;
    Qt3DRender::QFilterKey *m_filterKey;

    Q_DECLARE_PUBLIC(QMetalRoughMaterial)
};

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QMETALROUGHMATERIAL_P_H


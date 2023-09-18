/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DEXTRAS_QMORPHPHONGMATERIAL_P_H
#define QT3DEXTRAS_QMORPHPHONGMATERIAL_P_H

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
class QTechnique;
class QParameter;
class QShaderProgram;
class QRenderPass;
class QShaderProgramBuilder;

} // namespace Qt3DRender

namespace Qt3DExtras {

class QMorphPhongMaterial;

class QMorphPhongMaterialPrivate : public Qt3DRender::QMaterialPrivate
{
public:
    QMorphPhongMaterialPrivate();

    void init();

    void handleAmbientChanged(const QVariant &var);
    void handleDiffuseChanged(const QVariant &var);
    void handleSpecularChanged(const QVariant &var);
    void handleShininessChanged(const QVariant &var);
    void handleInterpolatorChanged(const QVariant &var);

    Qt3DRender::QEffect *m_phongEffect;
    Qt3DRender::QParameter *m_ambientParameter;
    Qt3DRender::QParameter *m_diffuseParameter;
    Qt3DRender::QParameter *m_specularParameter;
    Qt3DRender::QParameter *m_shininessParameter;
    Qt3DRender::QParameter *m_interpolatorParameter;
    Qt3DRender::QTechnique *m_phongGL3Technique;
    Qt3DRender::QTechnique *m_phongGL2Technique;
    Qt3DRender::QTechnique *m_phongES2Technique;
    Qt3DRender::QTechnique *m_phongRHITechnique;
    Qt3DRender::QRenderPass *m_phongGL3RenderPass;
    Qt3DRender::QRenderPass *m_phongGL2RenderPass;
    Qt3DRender::QRenderPass *m_phongES2RenderPass;
    Qt3DRender::QRenderPass *m_phongRHIRenderPass;
    Qt3DRender::QShaderProgram *m_phongGL3Shader;
    Qt3DRender::QShaderProgram *m_phongGL2ES2Shader;
    Qt3DRender::QShaderProgram *m_phongRHIShader;
    Qt3DRender::QShaderProgramBuilder *m_phongGL3ShaderBuilder;
    Qt3DRender::QShaderProgramBuilder *m_phongGL2ES2ShaderBuilder;
    Qt3DRender::QShaderProgramBuilder *m_phongRHIShaderBuilder;
    Qt3DRender::QFilterKey *m_filterKey;

    Q_DECLARE_PUBLIC(QMorphPhongMaterial)
};

} // Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QPHONGMATERIAL_P_H



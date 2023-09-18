/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the QtDataVisualization API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef SHADERHELPER_P_H
#define SHADERHELPER_P_H

#include "datavisualizationglobal_p.h"

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class ShaderHelper
{
    public:
    ShaderHelper(QObject *parent,
                 const QString &vertexShader = QString(),
                 const QString &fragmentShader = QString(),
                 const QString &texture = QString(),
                 const QString &depthTexture = QString());
    ~ShaderHelper();

    void setShaders(const QString &vertexShader, const QString &fragmentShader);
    void setTextures(const QString &texture, const QString &depthTexture);

    void initialize();
    bool testCompile();
    void bind();
    void release();
    void setUniformValue(GLint uniform, const QVector2D &value);
    void setUniformValue(GLint uniform, const QVector3D &value);
    void setUniformValue(GLint uniform, const QVector4D &value);
    void setUniformValue(GLint uniform, const QMatrix4x4 &value);
    void setUniformValue(GLint uniform, GLfloat value);
    void setUniformValue(GLint uniform, GLint value);
    void setUniformValueArray(GLint uniform,  const QVector4D *values, int count);

    GLint MVP();
    GLint view();
    GLint model();
    GLint nModel();
    GLint depth();
    GLint lightP();
    GLint lightS();
    GLint ambientS();
    GLint shadowQ();
    GLint color();
    GLint texture();
    GLint shadow();
    GLint gradientMin();
    GLint gradientHeight();
    GLint lightColor();
    GLint volumeSliceIndices();
    GLint colorIndex();
    GLint cameraPositionRelativeToModel();
    GLint color8Bit();
    GLint textureDimensions();
    GLint sampleCount();
    GLint alphaMultiplier();
    GLint preserveOpacity();
    GLint maxBounds();
    GLint minBounds();
    GLint sliceFrameWidth();

    GLint posAtt();
    GLint uvAtt();
    GLint normalAtt();

    private:
    QObject *m_caller;
    QOpenGLShaderProgram *m_program;

    QString m_vertexShaderFile;
    QString m_fragmentShaderFile;

    QString m_textureFile;
    QString m_depthTextureFile;

    GLint m_positionAttr;
    GLint m_uvAttr;
    GLint m_normalAttr;

    GLint m_colorUniform;
    GLint m_viewMatrixUniform;
    GLint m_modelMatrixUniform;
    GLint m_invTransModelMatrixUniform;
    GLint m_depthMatrixUniform;
    GLint m_mvpMatrixUniform;
    GLint m_lightPositionUniform;
    GLint m_lightStrengthUniform;
    GLint m_ambientStrengthUniform;
    GLint m_shadowQualityUniform;
    GLint m_textureUniform;
    GLint m_shadowUniform;
    GLint m_gradientMinUniform;
    GLint m_gradientHeightUniform;
    GLint m_lightColorUniform;
    GLint m_volumeSliceIndicesUniform;
    GLint m_colorIndexUniform;
    GLint m_cameraPositionRelativeToModelUniform;
    GLint m_color8BitUniform;
    GLint m_textureDimensionsUniform;
    GLint m_sampleCountUniform;
    GLint m_alphaMultiplierUniform;
    GLint m_preserveOpacityUniform;
    GLint m_minBoundsUniform;
    GLint m_maxBoundsUniform;
    GLint m_sliceFrameWidthUniform;

    GLboolean m_initialized;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif

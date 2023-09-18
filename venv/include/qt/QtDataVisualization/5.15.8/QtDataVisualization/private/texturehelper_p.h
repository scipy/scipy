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

#ifndef TEXTUREHELPER_P_H
#define TEXTUREHELPER_P_H

#include "datavisualizationglobal_p.h"
#include <QtGui/QRgb>
#include <QtGui/QLinearGradient>
#if !defined(QT_OPENGL_ES_2)
// 3D Textures are not supported by ES set
#  include <QtGui/QOpenGLFunctions_2_1>
#endif

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class TextureHelper : protected QOpenGLFunctions
{
    public:
    TextureHelper();
    ~TextureHelper();

    // Ownership of created texture is transferred to caller
    GLuint create2DTexture(const QImage &image, bool useTrilinearFiltering = false,
                           bool convert = true, bool smoothScale = true, bool clampY = false);
    GLuint create3DTexture(const QVector<uchar> *data, int width, int height, int depth,
                           QImage::Format dataFormat);
    GLuint createCubeMapTexture(const QImage &image, bool useTrilinearFiltering = false);
    // Returns selection texture and inserts generated framebuffers to framebuffer parameters
    GLuint createSelectionTexture(const QSize &size, GLuint &frameBuffer, GLuint &depthBuffer);
    GLuint createCursorPositionTexture(const QSize &size, GLuint &frameBuffer);
    GLuint createUniformTexture(const QColor &color);
    GLuint createGradientTexture(const QLinearGradient &gradient);
    GLuint createDepthTexture(const QSize &size, GLuint textureSize);
    // Returns depth texture and inserts generated framebuffer to parameter
    GLuint createDepthTextureFrameBuffer(const QSize &size, GLuint &frameBuffer, GLuint textureSize);
    void deleteTexture(GLuint *texture);

    private:
    QImage convertToGLFormat(const QImage &srcImage);
    void convertToGLFormatHelper(QImage &dstImage, const QImage &srcImage, GLenum texture_format);
    QRgb qt_gl_convertToGLFormatHelper(QRgb src_pixel, GLenum texture_format);

#if !defined(QT_OPENGL_ES_2)
    QOpenGLFunctions_2_1 *m_openGlFunctions_2_1; // Not owned
#endif
    friend class Bars3DRenderer;
    friend class Surface3DRenderer;
    friend class Scatter3DRenderer;
    friend class Abstract3DRenderer;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif

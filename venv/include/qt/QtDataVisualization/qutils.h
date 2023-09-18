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

#ifndef QUTILS_H
#define QUTILS_H

#include <QtGui/QSurfaceFormat>
#include <QtGui/QOpenGLContext>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QOffscreenSurface>
#include <QtCore/QCoreApplication>

namespace QtDataVisualization {

#ifndef Q_QDOC
static inline QSurfaceFormat qDefaultSurfaceFormat(bool antialias = true) Q_DECL_UNUSED;
#endif
static inline QSurfaceFormat qDefaultSurfaceFormat(bool antialias)
{
    bool isES = false;

    QSurfaceFormat surfaceFormat;

    // Common attributes
    surfaceFormat.setDepthBufferSize(24);
    surfaceFormat.setStencilBufferSize(8);
    surfaceFormat.setSwapBehavior(QSurfaceFormat::DoubleBuffer);
    surfaceFormat.setRenderableType(QSurfaceFormat::DefaultRenderableType);

    QOpenGLContext *ctx = QOpenGLContext::currentContext();
    QOffscreenSurface *dummySurface = nullptr;
    if (!ctx) {
        dummySurface = new QOffscreenSurface();
        dummySurface->setFormat(surfaceFormat);
        dummySurface->create();
        ctx = new QOpenGLContext;
        ctx->setFormat(surfaceFormat);
        ctx->create();
        ctx->makeCurrent(dummySurface);
    }

#if defined(QT_OPENGL_ES_2)
    isES = true;
#elif (QT_VERSION < QT_VERSION_CHECK(5, 3, 0))
    isES = false;
#else
    isES = ctx->isOpenGLES();
#endif

#if (QT_VERSION >= QT_VERSION_CHECK(5, 4, 0))
    // We support only ES2 emulation with software renderer for now
    QString versionStr;
#ifdef Q_OS_WIN
    const GLubyte *openGLVersion = ctx->functions()->glGetString(GL_VERSION);
    versionStr = QString::fromLatin1(reinterpret_cast<const char *>(openGLVersion)).toLower();
#endif
    if (versionStr.contains(QStringLiteral("mesa"))
            || QCoreApplication::testAttribute(Qt::AA_UseSoftwareOpenGL)) {
        qWarning("Only OpenGL ES2 emulation is available for software rendering.");
        isES = true;
    }
#endif

    if (dummySurface) {
        ctx->doneCurrent();
        delete ctx;
        delete dummySurface;
    }

    if (isES) {
        // For ES2 only attributes
        surfaceFormat.setRedBufferSize(8);
        surfaceFormat.setBlueBufferSize(8);
        surfaceFormat.setGreenBufferSize(8);
    } else {
        surfaceFormat.setVersion(2, 1);
        surfaceFormat.setProfile(QSurfaceFormat::CoreProfile);
        // For OpenGL only attributes
        if (antialias)
            surfaceFormat.setSamples(8);
        else
            surfaceFormat.setSamples(0);
    }

    return surfaceFormat;
}

}

#endif

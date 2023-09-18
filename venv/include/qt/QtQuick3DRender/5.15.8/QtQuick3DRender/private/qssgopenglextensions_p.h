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


#ifndef QSSGOPENGLEXTENSIONS_H
#define QSSGOPENGLEXTENSIONS_H

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

#include <QtOpenGLExtensions/QtOpenGLExtensions>

QT_BEGIN_NAMESPACE

/* Some OpenGL extensions that are not (yet) found in Qt's OpenGL extensions.
 * These should be auto-generated and added to QtOpenGLExtensions module */
class QSSGOpenGLExtensionsPrivate : public QAbstractOpenGLExtensionPrivate
{
public:
    void(QOPENGLF_APIENTRYP BlendBarrierNV)();

#if defined(QT_OPENGL_ES) || defined(QT_OPENGL_ES_2_ANGLE)
    void(QOPENGLF_APIENTRYP PatchParameteriEXT)(GLenum, GLint);
    void(QOPENGLF_APIENTRYP QueryCounterEXT)(GLuint, GLenum);
    void(QOPENGLF_APIENTRYP GetQueryObjectui64vEXT)(GLuint, GLenum, GLuint64 *);
    void(QOPENGLF_APIENTRYP BindVertexArrayOES)(GLuint array);
    void(QOPENGLF_APIENTRYP DeleteVertexArraysOES)(GLsizei n, const GLuint *arrays);
    void(QOPENGLF_APIENTRYP GenVertexArraysOES)(GLsizei n, GLuint *arrays);
    GLboolean(QOPENGLF_APIENTRYP IsVertexArrayOES)(GLuint array);
#endif
};

class QSSGOpenGLExtensions : public QAbstractOpenGLExtension
{
public:
    QSSGOpenGLExtensions();

    bool initializeOpenGLFunctions() override;

    void glBlendBarrierNV();

protected:
    Q_DECLARE_PRIVATE(QSSGOpenGLExtensions)
};

inline void QSSGOpenGLExtensions::glBlendBarrierNV()
{
    Q_D(QSSGOpenGLExtensions);
    d->BlendBarrierNV();
}

#if defined(QT_OPENGL_ES) || defined(QT_OPENGL_ES_2_ANGLE)
class QSSGOpenGLES2Extensions : public QSSGOpenGLExtensions
{
public:
    QSSGOpenGLES2Extensions();

    // tesselation shader
    void glPatchParameteriEXT(GLenum pname, GLint value);

    // timer
    void glQueryCounterEXT(GLuint id, GLenum target);
    void glGetQueryObjectui64vEXT(GLuint id, GLenum pname, GLuint64 *params);

    void glBindVertexArrayOES(GLuint array);
    void glDeleteVertexArraysOES(GLsizei n, const GLuint *arrays);
    void glGenVertexArraysOES(GLsizei n, GLuint *arrays);
    GLboolean glIsVertexArrayOES(GLuint array);

    bool initializeOpenGLFunctions() Q_DECL_FINAL;
};

inline void QSSGOpenGLES2Extensions::glPatchParameteriEXT(GLenum pname, GLint value)
{
    Q_D(QSSGOpenGLExtensions);
    d->PatchParameteriEXT(pname, value);
}

inline void QSSGOpenGLES2Extensions::glQueryCounterEXT(GLuint id, GLenum target)
{
    Q_D(QSSGOpenGLExtensions);
    d->QueryCounterEXT(id, target);
}

inline void QSSGOpenGLES2Extensions::glGetQueryObjectui64vEXT(GLuint id, GLenum pname, GLuint64 *params)
{
    Q_D(QSSGOpenGLExtensions);
    d->GetQueryObjectui64vEXT(id, pname, params);
}

inline void QSSGOpenGLES2Extensions::glBindVertexArrayOES(GLuint array)
{
    Q_D(QSSGOpenGLExtensions);
    d->BindVertexArrayOES(array);
}

inline void QSSGOpenGLES2Extensions::glDeleteVertexArraysOES(GLsizei n, const GLuint *arrays)
{
    Q_D(QSSGOpenGLExtensions);
    d->DeleteVertexArraysOES(n, arrays);
}

inline void QSSGOpenGLES2Extensions::glGenVertexArraysOES(GLsizei n, GLuint *arrays)
{
    Q_D(QSSGOpenGLExtensions);
    d->GenVertexArraysOES(n, arrays);
}

inline GLboolean QSSGOpenGLES2Extensions::glIsVertexArrayOES(GLuint array)
{
    Q_D(QSSGOpenGLExtensions);
    return d->IsVertexArrayOES(array);
}

#endif // QT_OPENGL_ES

QT_END_NAMESPACE

#endif // QSSGOPENGLEXTENSIONS_H

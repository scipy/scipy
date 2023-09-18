/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
** Contact: http://www.qt-project.org/legal
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

#ifndef QT3DRENDER_RENDER_SHADERCACHE_H
#define QT3DRENDER_RENDER_SHADERCACHE_H

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

#include <Qt3DRender/private/qt3drender_global_p.h>
#include <Qt3DRender/private/shader_p.h>

#include <QtCore/qhash.h>
#include <QtCore/qmutex.h>
#include <QtCore/qvector.h>

QT_BEGIN_NAMESPACE

class QOpenGLShaderProgram;

namespace Qt3DRender {
namespace Render {

#if defined(QT_BUILD_INTERNAL)
class tst_ShaderCache;
#endif

class Q_3DRENDERSHARED_PRIVATE_EXPORT ShaderCache
{
public:
    ~ShaderCache();

    QOpenGLShaderProgram *getShaderProgramAndAddRef(ProgramDNA dna, Qt3DCore::QNodeId shaderPeerId, bool *wasPresent = nullptr);
    void insert(ProgramDNA dna, Qt3DCore::QNodeId shaderPeerId, QOpenGLShaderProgram *program);
    void removeRef(ProgramDNA dna, Qt3DCore::QNodeId shaderPeerId);
    void purge();
    void clear();

    // Only ever used from the OpenGL submission thread
    QOpenGLShaderProgram *getShaderProgramForDNA(ProgramDNA dna) const;
    QVector<Qt3DCore::QNodeId> shaderIdsForProgram(ProgramDNA dna) const;

private:
    // Only ever used from the OpenGL submission thread
    QHash<ProgramDNA, QOpenGLShaderProgram *> m_programHash;

    // Accessed from both the OpenGL submission thread and the aspect thread
    QHash<ProgramDNA, QVector<Qt3DCore::QNodeId>> m_programRefs;
    QVector<ProgramDNA> m_pendingRemoval;
    QMutex m_refsMutex;

#if defined(QT_BUILD_INTERNAL)
    friend class tst_ShaderCache;
#endif
};

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_SHADERCACHE_H

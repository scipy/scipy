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

#ifndef QSSG_RENDER_BACKEND_SHADER_PROGRAM_GL_H
#define QSSG_RENDER_BACKEND_SHADER_PROGRAM_GL_H

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

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DUtils/private/qssgutils_p.h>
#include <QtQuick3DUtils/private/qssgoption_p.h>
#include <QtCore/QString>

QT_BEGIN_NAMESPACE

struct QSSGRenderBackendShaderInputEntryGL
{
    QByteArray m_attribName; ///< must be the same name as used in the vertex shader
    quint32 m_attribLocation; ///< attribute index
    quint32 m_type; ///< GL vertex format type @sa GL_FLOAT, GL_INT
    quint32 m_numComponents; ///< component count. max 4
};

///< this class handles the shader input variables
class QSSGRenderBackendShaderInputGL
{
public:
    ///< constructor
    QSSGRenderBackendShaderInputGL(QSSGDataRef<QSSGRenderBackendShaderInputEntryGL> entries)
        : m_shaderInputEntries(entries)
    {
    }
    ///< destructor
    ~QSSGRenderBackendShaderInputGL()
    {
        // We need to release the attrib name
        for (int idx = 0; idx != m_shaderInputEntries.size(); ++idx)
            m_shaderInputEntries[idx] = QSSGRenderBackendShaderInputEntryGL();
    }

    QSSGRenderBackendShaderInputEntryGL *getEntryByName(const QByteArray &entryName) const
    {
        for (int idx = 0; idx != m_shaderInputEntries.size(); ++idx) {
            if (m_shaderInputEntries[idx].m_attribName == entryName)
                return &m_shaderInputEntries.mData[idx];
        }
        return nullptr;
    }

    QSSGOption<QSSGRenderBackendShaderInputEntryGL> getEntryByAttribLocation(quint32 attribLocation) const
    {
        for (int idx = 0; idx != m_shaderInputEntries.size(); ++idx) {
            if (m_shaderInputEntries[idx].m_attribLocation == attribLocation)
                return m_shaderInputEntries[idx];
        }
        return QSSGEmpty();
    }

    QSSGDataRef<QSSGRenderBackendShaderInputEntryGL> m_shaderInputEntries; ///< shader input entries
};

///< this class represents the internals of a GL program
class QSSGRenderBackendShaderProgramGL
{
public:
    ///< constructor
    QSSGRenderBackendShaderProgramGL(quint32 programID) : m_programID(programID), m_shaderInput(nullptr) {}

    ///< destructor
    ~QSSGRenderBackendShaderProgramGL()
    {
        if (m_shaderInput) { // Created with malloc, so release with free!
            m_shaderInput->~QSSGRenderBackendShaderInputGL();
            ::free(m_shaderInput);
        }
    }

    quint32 m_programID; ///< this is the OpenGL object ID
    QSSGRenderBackendShaderInputGL *m_shaderInput; ///< pointer to shader input object
};

QT_END_NAMESPACE

#endif

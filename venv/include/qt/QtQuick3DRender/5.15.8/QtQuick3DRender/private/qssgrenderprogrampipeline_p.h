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

#ifndef QSSG_RENDER_PROGRAM_PIPLINE_H
#define QSSG_RENDER_PROGRAM_PIPLINE_H

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

#include <QtQuick3DRender/private/qssgrendercontext_p.h>
#include <QtQuick3DRender/private/qssgrenderbackend_p.h>

QT_BEGIN_NAMESPACE

class QSSGRenderContext;
class QSSGRenderShaderProgram;

///< A program pipeline is a collection of a multiple programs (vertex, fragment, geometry,....)
class Q_QUICK3DRENDER_EXPORT QSSGRenderProgramPipeline
{
    Q_DISABLE_COPY(QSSGRenderProgramPipeline)
public:
    QAtomicInt ref;

protected:
    QSSGRef<QSSGRenderContext> m_context; ///< pointer to context
    QSSGRef<QSSGRenderBackend> m_backend; ///< pointer to backend

public:
    /**
     * @brief constructor
     *
     * @param[in] context		Pointer to render context
     * @param[in] fnd			Pointer to foundation
     *
     * @return No return.
     */
    QSSGRenderProgramPipeline(const QSSGRef<QSSGRenderContext> &context);

    /// @brief destructor
    ~QSSGRenderProgramPipeline();

    /**
     * @brief Query if pipeline is valid
     *
     * @return True if valid.
     */
    bool isValid();

    /**
     * @brief enable / disable a program stage in the pipeline
     *
     * @param[in] pProgram	Pointer to program. If nullptr stage will be disabled
     * @param[in] flags		Flags to which stage this program is bound to. Can more than
     * one stage
     *
     * @return no return.
     */
    void setProgramStages(const QSSGRef<QSSGRenderShaderProgram> &pProgram, QSSGRenderShaderTypeFlags flags);

    /**
     * @brief Make the program pipeline active
     *
     * @return True if valid.
     */
    void bind();

    /**
     * @brief get the backend object handle
     *
     * @return the backend object handle.
     */
    QSSGRenderBackend::QSSGRenderBackendProgramPipeline handle() { return m_handle; }

    /**
     * @brief get the vertex stage program
     *
     * @return the backend object handle.
     */
    QSSGRef<QSSGRenderShaderProgram> vertexStage();

private:
    QSSGRenderBackend::QSSGRenderBackendProgramPipeline m_handle; ///< opaque backend handle

    QSSGRef<QSSGRenderShaderProgram> m_program; ///< for non separable programs this contains the entire program
    QSSGRef<QSSGRenderShaderProgram> m_vertexProgram; ///< for separable programs this contains the vertex program
    QSSGRef<QSSGRenderShaderProgram> m_fragmentProgram; ///< for separable programs this contains the fragment program
    QSSGRef<QSSGRenderShaderProgram> m_tessControlProgram; ///< for separable programs this contains the
    /// tessellation control program
    QSSGRef<QSSGRenderShaderProgram> m_tessEvalProgram; ///< for separable programs this contains the
    /// tessellation evaluation program
    QSSGRef<QSSGRenderShaderProgram> m_geometryProgram; ///< for separable programs this contains the geometry program
    QSSGRef<QSSGRenderShaderProgram> m_computProgram; ///< for separable programs this contains the compute program
};

QT_END_NAMESPACE

#endif

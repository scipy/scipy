/****************************************************************************
**
** Copyright (C) 2016 Paul Lemire
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

#ifndef QT3DRENDER_RENDER_GENERICLAMBDAJOB_H
#define QT3DRENDER_RENDER_GENERICLAMBDAJOB_H

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

#include <Qt3DCore/qaspectjob.h>
#include <Qt3DRender/private/job_common_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

namespace Render {

template<typename T>
class GenericLambdaJob : public Qt3DCore::QAspectJob
{
public:
    explicit GenericLambdaJob(T callable, JobTypes::JobType type = JobTypes::GenericLambda, const char *name = "GenericLambda", int instance = 0)
        : Qt3DCore::QAspectJob()
        , m_callable(callable)
    {
        SET_JOB_RUN_STAT_TYPE_AND_NAME(this, type, name, instance)
    }

    // QAspectJob interface
    void run() final
    {
        m_callable();
    }

private:
    T m_callable;
};

template<typename T>
using GenericLambdaJobPtr = QSharedPointer<GenericLambdaJob<T>>;

template<typename T, typename U>
class GenericLambdaJobAndPostFramePrivate : public Qt3DCore::QAspectJobPrivate
{
public:
    explicit GenericLambdaJobAndPostFramePrivate(U postFrameCallable)
        : m_postFrameCallable(postFrameCallable)
    {}

    ~GenericLambdaJobAndPostFramePrivate() override {}

    void postFrame(Qt3DCore::QAspectManager *manager) override
    {
        m_postFrameCallable(manager);
    }

private:
    U m_postFrameCallable;
};

template<typename T, typename U>
class GenericLambdaJobAndPostFrame : public Qt3DCore::QAspectJob
{
public:
    explicit GenericLambdaJobAndPostFrame(T runCallable, U postFrameCallable, JobTypes::JobType type = JobTypes::GenericLambda, const char *name = "GenericLambda")
        : Qt3DCore::QAspectJob(*new GenericLambdaJobAndPostFramePrivate<T, U>(postFrameCallable))
        , m_runCallable(runCallable)
    {
        SET_JOB_RUN_STAT_TYPE_AND_NAME(this, type, name, 0)
    }

    // QAspectJob interface
    void run() final
    {
        m_runCallable();
    }

private:
    T m_runCallable;
};

template<typename T, typename U>
using GenericLambdaJobAndPostFramePtr = QSharedPointer<GenericLambdaJobAndPostFrame<T, U>>;

} // Render

} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_GENERICLAMBDAJOB_H

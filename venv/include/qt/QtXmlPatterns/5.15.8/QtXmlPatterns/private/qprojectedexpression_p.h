/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtXmlPatterns module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef Patternist_ProjectedExpression_H
#define Patternist_ProjectedExpression_H

#include <private/qitem_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    class ProjectedExpression
    {
    public:
        typedef ProjectedExpression * Ptr;
        typedef QVector<ProjectedExpression::Ptr> Vector;
        virtual ~ProjectedExpression()
        {
        }

        enum Action
        {
            Move = 0,
            Skip = 1,
            Keep = 2,
            KeepSubtree = 4 | Keep
        };

        virtual Action actionForElement(const QXmlName name,
                                        ProjectedExpression::Ptr &next) const
        {
            Q_UNUSED(name);
            Q_UNUSED(next);
            return Skip;
        }

    };

    class ProjectedNodeTest
    {
    public:
        typedef ProjectedNodeTest * Ptr;
        virtual ~ProjectedNodeTest()
        {
        }

        virtual bool isMatch(const QXmlNodeModelIndex::NodeKind kind) const
        {
            Q_UNUSED(kind);
            return false;
        }
    };

    class ProjectedStep : public ProjectedExpression
    {
    public:
        ProjectedStep(const ProjectedNodeTest::Ptr test,
                      const QXmlNodeModelIndex::Axis axis)
        {
            Q_ASSERT(test);
            Q_UNUSED(test);
            Q_UNUSED(axis);
        }

        virtual Action actionForElement(const QXmlName name,
                                        ProjectedExpression::Ptr &next) const
        {
            Q_UNUSED(name);
            Q_UNUSED(next);
            // TODO
            return Skip;
        }

    private:
    };

    class ProjectedPath : public ProjectedExpression
    {
    public:
        ProjectedPath(const ProjectedExpression::Ptr left,
                      const ProjectedExpression::Ptr right) : m_left(left)
        {
            Q_ASSERT(m_left);
            Q_ASSERT(right);
            Q_UNUSED(right);
        }

        virtual Action actionForElement(const QXmlName name,
                                        ProjectedExpression::Ptr &next) const
        {
            ProjectedExpression::Ptr &candidateNext = next;
            const Action a = m_left->actionForElement(name, candidateNext);

            if(a != Skip)
            {
                /* The test accepted it, so let's replace us with the new step. */
                next = candidateNext;
            }

            return a;
        }

    private:
        const ProjectedExpression::Ptr  m_left;
    };
}

QT_END_NAMESPACE

#endif

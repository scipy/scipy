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

#ifndef Patternist_TemplatePattern_H
#define Patternist_TemplatePattern_H

#include <private/qtemplate_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Houses the data necessary for a template pattern.
     *
     * A template pattern is the match pattern, but have had each operand to @c
     * | separated out into a separate TemplatePattern. For instance, the
     * pattern <tt>a | b | c</tt>, becomes three separate TemplatePattern
     * instances.
     *
     * @see TemplateMode
     * @see Template
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     * @since 4.5
     */
    class TemplatePattern : public QSharedData
    {
    public:
        typedef QExplicitlySharedDataPointer<TemplatePattern> Ptr;
        typedef QVector<Ptr> Vector;
        typedef int ID;

        inline TemplatePattern(const Expression::Ptr &matchPattern,
                               const PatternPriority pri,
                               const ID id,
                               const Template::Ptr templ);

        inline PatternPriority priority() const;
        inline const Expression::Ptr &matchPattern() const;
        inline void setMatchPattern(const Expression::Ptr &pattern);
        inline const Template::Ptr &templateTarget() const;
        inline ID id() const;

        /**
         * This ID is used to ensure that, as 6.4 Conflict Resolution for
         * Template Rules reads:
         *
         * "If the pattern contains multiple alternatives separated by |, then
         * the template rule is treated equivalently to a set of template
         * rules, one for each alternative. However, it is not an error if a
         * node matches more than one of the alternatives."
         *
         * For patterns separated by @c |, we have one Template instance for
         * each alternative, but they all have the same ID, hence if several
         * alternatives match, we don't flag it as an error if they have the
         * same ID.
         */
    private:
        Expression::Ptr m_matchPattern;
        PatternPriority m_priority;
        ID              m_id;
        Template::Ptr   m_templateTarget;
        Q_DISABLE_COPY(TemplatePattern)
    };

    TemplatePattern::TemplatePattern(const Expression::Ptr &matchPattern,
                                     const PatternPriority pri,
                                     const ID id,
                                     const Template::Ptr templ) : m_matchPattern(matchPattern)
                                                                , m_priority(pri)
                                                                , m_id(id)
                                                                , m_templateTarget(templ)

    {
        Q_ASSERT(m_matchPattern);
        Q_ASSERT(m_templateTarget);
    }

    const Expression::Ptr &TemplatePattern::matchPattern() const
    {
        return m_matchPattern;
    }

    void TemplatePattern::setMatchPattern(const Expression::Ptr &pattern)
    {
        m_matchPattern = pattern;
    }

    PatternPriority TemplatePattern::priority() const
    {
        return m_priority;
    }

    TemplatePattern::ID TemplatePattern::id() const
    {
        return m_id;
    }

    const Template::Ptr &TemplatePattern::templateTarget() const
    {
        return m_templateTarget;
    }

}

QT_END_NAMESPACE

#endif


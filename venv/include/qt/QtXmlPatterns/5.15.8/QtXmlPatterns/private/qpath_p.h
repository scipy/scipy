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

#ifndef Patternist_Path_H
#define Patternist_Path_H

#include <private/qpaircontainer_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{

    /**
     * @short Implements the path expression, containing two steps, such as in <tt>html/body</tt>.
     *
     * @see <a href="http://www.w3.org/TR/xquery/#id-path-expressions">XQuery 1.0: An
     * XML Query Language, 3.2 Path Expressions</a>
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class Q_AUTOTEST_EXPORT Path : public PairContainer
    {
    public:
        enum Kind
        {
            /**
             * This Path is a plain old path expression as found in XPath.
             * Sorting is performed, and atomics are disallowed as left
             * operand.
             */
            RegularPath = 1,

            /**
             * This Path emulates an @c xsl:for-each instruction. This means no
             * sorting of result, and atomics are allowed as left operand.
             */
            XSLTForEach,

            /**
             * This Path performs the iteration in an @c xsl:apply-templates
             * instruction. This means sorting, and atomics are disallowed
             * as left operand.
             */
            ForApplyTemplate
        };

        Path(const Expression::Ptr &operand1,
             const Expression::Ptr &operand2,
             const Kind kind = RegularPath);
        ~Path();

        virtual Item::Iterator::Ptr evaluateSequence(const DynamicContext::Ptr &context) const;
        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
        virtual void evaluateToSequenceReceiver(const DynamicContext::Ptr &context) const;
        inline Item::Iterator::Ptr mapToSequence(const Item &item,
                                                 const DynamicContext::Ptr &context) const;

        virtual Expression::Ptr typeCheck(const StaticContext::Ptr &context,
                                          const SequenceType::Ptr &reqType);

        virtual SequenceType::List expectedOperandTypes() const;

        /**
         * @returns the static type of the last step where the cardinality is multiplied with
         * the cardinality of the first step's cardinality.
         */
        virtual SequenceType::Ptr staticType() const;

        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;

        virtual Properties properties() const;

        virtual Expression::Ptr compress(const StaticContext::Ptr &context);

        /**
         * @returns the item type of the last step's static type.
         */
        virtual ItemType::Ptr newFocusType() const;

        virtual ID id() const;

        inline void setLast();

        inline Kind kind() const
        {
            return m_kind;
        }

    private:
        typedef QExplicitlySharedDataPointer<const Path> ConstPtr;

        /**
         * One might think this block exists for preventing multiple
         * NodeSortExpressions to be created. However, that is not an issue,
         * since NodeSortExpression optimizes this away anyway.
         *
         * The real reason is to avoid infinite recursion. When our typeCheck()
         * forwards on the type check to the just created
         * NodeSortExpression, it in turn calls typeCheck() on its child, which
         * is this Path. Rince and repeat.
         *
         * We only create node sorts when we're a regular path expression, and
         * not when standing in as a generic map expression. */
        bool        m_hasCreatedSorter;

        /**
         * Whether this path is the step. For instance, in <tt>a/b/c</tt>, the
         * last path has @c c as the right operand.
         */
        bool        m_isLast;

        bool        m_checkXPTY0018;
        const Kind  m_kind;
    };

    void Path::setLast()
    {
        m_isLast = true;
    }
}

QT_END_NAMESPACE

#endif

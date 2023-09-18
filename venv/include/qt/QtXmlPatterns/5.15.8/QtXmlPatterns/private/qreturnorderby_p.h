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

#ifndef Patternist_ReturnOrderBy_H
#define Patternist_ReturnOrderBy_H

#include <private/qorderby_p.h>
#include <private/qunlimitedcontainer_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Together with OrderBy, it implements XQuery 1.0's <tt>order by</tt> expression.
     *
     * ReturnOrderBy evaluates the sort keys and values, and hands it over to
     * OrderBy, which is an AST ancestor, using SortTuples.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     * @ingroup Patternist_expressions
     */
    class ReturnOrderBy : public UnlimitedContainer
    {
    public:
        /**
         * In @p operands the first item is the return expression, and the
         * rest, which is at least one, are the sort keys.
         */
        ReturnOrderBy(const OrderBy::Stability stability,
                      const OrderBy::OrderSpec::Vector &oSpecs,
                      const Expression::List &operands);

        virtual Item evaluateSingleton(const DynamicContext::Ptr &context) const;
        virtual bool evaluateEBV(const DynamicContext::Ptr &context) const;
        virtual SequenceType::Ptr staticType() const;
        virtual SequenceType::List expectedOperandTypes() const;
        virtual Expression::Ptr compress(const StaticContext::Ptr &context);
        virtual ExpressionVisitorResult::Ptr accept(const ExpressionVisitor::Ptr &visitor) const;
        virtual ID id() const;

        inline OrderBy::OrderSpec::Vector orderSpecs() const
        {
            return m_orderSpecs;
        }

        inline OrderBy::Stability stability() const
        {
            return m_stability;
        }

        /**
         * In the case of that we don't have a for-expression beloning us, but
         * only a let clause, this ReturnOrderBy breaks if it stays in the AST.
         * So, by default we assume that we should write ourselves away, unless
         * this function is called. The associated ForClause will call it
         * during typeCheck(), if it exists.
         */
        inline void setStay(const bool a)
        {
            m_flyAway = !a;
        }

        virtual Properties properties() const;
    private:
        /**
         * This variable is unfortunately only used at compile time. However,
         * it's tricky to get rid of it due to how QueryTransformParser would
         * have to be adapted.
         */
        const OrderBy::Stability    m_stability;

        OrderBy::OrderSpec::Vector  m_orderSpecs;

        /**
         * @see stay()
         */
        bool                        m_flyAway;
    };
}

QT_END_NAMESPACE

#endif
